/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/


/*!
  \file chem/pointgroupsymmetry.cc
  \brief implements point group operations

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/chem/pointgroupsymmetry.cc>here</a>.

*/

#include<madness/chem/pointgroupsymmetry.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/functypedefs.h>
#include <madness/fortran_ctypes.h>

namespace madness {


/// (re-)project the argument on the given irreps
template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > projector_irrep::project_on_irreps(
		const std::vector<Function<T,NDIM> >& vrhs,
		const std::vector<std::string>& irreps) const {

	World& world=vrhs.begin()->world();

	// fast return
	if (vrhs.size()==0) return std::vector<Function<T,NDIM> >();
	if (table_.schoenflies_=="C1") {return copy(world,vrhs);}

	// loop over all symmetry operations of this point group and apply them
	// to the input vector; dimension: (nop, nmo)
	std::vector<std::vector<Function<T,NDIM> > > opvrhs(table_.operators_.size());
	for (size_t i=0; i<table_.operators_.size(); ++i) {
		opvrhs[i]=table_.operators_[i](vrhs,false);
	}
	world.gop.fence();
	for (size_t i=0; i<table_.operators_.size(); ++i) compress(world,opvrhs[i],false);
	world.gop.fence();

	std::vector<Function<T,NDIM> > result=zero_functions_compressed<T,NDIM>(world,vrhs.size());

	// loop over all vector elements
	for (size_t ivec=0; ivec<vrhs.size(); ++ivec) {

		// get the characters for this vector element
		const charactertable::characterlineT cline=table_.irreps_.find(irreps[ivec])->second;
		for (size_t iop=0; iop<table_.operators_.size(); ++iop) {
			double character=double(cline[iop]);

			// apply the projection operator P^Gamma v_i = \sum_op \chi_op^Gamma op(v_i)
			result[ivec].gaxpy(1.0,opvrhs[iop][ivec],character/table_.order_,false);
		}
	}
	world.gop.fence();

	truncate(world,result);
	return result;
}


/// symmetrize a vector of functions

/// project a number of functions on a given number of irreps
/// linear dependencies are removed, the functions are orthonormalized
/// vanishing input functions are mapped onto vanishing output functions with irrep "null"
/// the number of input and output function is expected to be the same!
/// @param[in]	vrhs	vector of functions to be projected on the irrep
/// @param[in]	vbra	bra of vrhs if applicable (if bra /= ket), may be empty
/// @param[out]	sirrep	vector with the irrep names corresponding to the result
template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > projector_irrep::apply_symmetry_operators(
		const std::vector<Function<T,NDIM> >& vrhs,
        Function<typename Tensor<T>::scalar_type,NDIM> metric,
		std::vector<std::string>& sirreps) const {

	World& world=vrhs.begin()->world();

	// fast return
	if (vrhs.size()==0) return std::vector<Function<T,NDIM> >();
	if (table_.schoenflies_=="C1") {
		sirreps=std::vector<std::string>(vrhs.size(),"a");
		return copy(world,vrhs);
	}

	// loop over all symmetry operations of this point group and apply them
	// to the input vector; dimension: (nop, nmo)
	std::vector<std::vector<Function<T,NDIM> > > opvrhs(table_.operators_.size());
	for (size_t i=0; i<table_.operators_.size(); ++i) {
		opvrhs[i]=table_.operators_[i](vrhs,false);
	}

	// get all irreps to work on
	std::vector<std::string> all_irreps;
	if (irrep_=="all") {
		all_irreps=table_.mullikan_;
	} else {
		all_irreps.push_back(irrep_);
	}

	// accumulate the linearly combined results here; dimension: (nirrep, nmo)
	std::vector<std::vector<Function<T,NDIM> > > lc_op_vrhs(all_irreps.size());
	for (auto& v : lc_op_vrhs) {
		v=zero_functions_compressed<T,NDIM>(world,vrhs.size(),false);
	}

	// need to fence here before the linear combination with the character weights
	world.gop.fence();

	// loop over all requested irreps
	for (auto& o : opvrhs) change_tree_state(o,compressed);
	for (size_t j=0; j<all_irreps.size(); ++j) {
		std::string& irrep=all_irreps[j];

		// get the characters of this irrep
		const charactertable::characterlineT cline=table_.irreps_.find(irrep)->second;
		MADNESS_ASSERT(table_.operators_.size()==cline.size());


		// do the linear combination; loop over all operators
		// lc_op_vhrs[jrrep] = \sum_iop (character[iop]/table_.order_)*opvhrs[iop];
		for (size_t iop=0; iop<cline.size(); ++iop) {
			double character=double(cline[iop]);
			gaxpy(world,1.0,lc_op_vrhs[j],character/table_.order_,opvrhs[iop],false);
		}
	}

	// need to fence here before the orthogonalization step
	world.gop.fence();
	for (auto& v : lc_op_vrhs) truncate(world,v,0.0,false);
	world.gop.fence();

	// now we have a linearly dependent set of symmetrized functions
	std::vector<Function<T,NDIM> > result1;
	std::vector<std::string> sirreps1;
	std::vector<int> ipiv1;

	// loop over all irreps
	//int iresult=0;	// counter for the result vector ordering
	for (size_t j=0; j<all_irreps.size(); ++j) {

		if (verbosity_>1) print("working on irrep ",all_irreps[j]);
		// compute overlap, include metric if necessary
		Tensor<T> ovlp_t;
		if (metric.is_initialized()) ovlp_t=matrix_inner(world,lc_op_vrhs[j],metric*lc_op_vrhs[j]);
		else ovlp_t=matrix_inner(world,lc_op_vrhs[j],lc_op_vrhs[j]);
		// overlap matrix should be real
		Tensor<double> ovlp=real(ovlp_t);
        Tensor<double> imag_ovlp=imag(ovlp_t);
        MADNESS_CHECK(imag_ovlp.normf()<1.e-10);

//			for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)+=1.e-6;

		if (verbosity_>2) {
			print("ovlp");
			print(ovlp);
		}

		// check for zero rank
		Tensor<double> evec,eval;
		syev(ovlp,evec,eval);

		if (eval.sum()>lindep_) {

			if (verbosity_>1) print("eigenvalues",eval);

			// compute rank-revealing Cholesky, discard linear dependent vectors
			Tensor<integer> piv;
			int rank;
			rr_cholesky(ovlp,lindep_,piv,rank);
			Tensor<double> shrunk=copy(ovlp(Slice(0,rank-1,1),Slice(0,rank-1,1)));
			Tensor<double> inv=inverse(shrunk);

			if (verbosity_>2) {
				print("cholesky");
				print(ovlp);
				print("inverse(cholesky)");
				print(inv);
				print("orthogonality");
				print(inner(shrunk,shrunk,0,0));
			}
			if (verbosity_>1) {
				print("pivot elements");
				print(piv);
			}

			// collect non-vanishing elements
			std::vector<Function<T,NDIM> > a1shrunk;
			for (int i=0; i<rank; ++i) a1shrunk.push_back(lc_op_vrhs[j][piv(i)]);

			// orthonormalize vectors within this irrep
			if (orthonormalize_irreps_) a1shrunk=transform(world,a1shrunk,inv);

			// collect all result functions and their irreps
			for (size_t i=0; i<a1shrunk.size(); ++i) {
				result1.push_back(a1shrunk[i]);
				sirreps1.push_back(all_irreps[j]);
				ipiv1.push_back(piv[i]);
			}
		}
	}

	// check closure
	if (vrhs.size()>1 and (sirreps1.size()>vrhs.size())) {

		print("vrhs.size()            ",vrhs.size());
		print("raw sirreps            ",sirreps1);
		print("raw mapping in ipiv1   ",ipiv1);

		MADNESS_EXCEPTION("\n\nfunction arguments in apply_symmetry_operators are not closed\n\n",1);
	}

	if (verbosity_>1) {
		print("raw sirreps            ",sirreps1);
		print("raw mapping in ipiv1   ",ipiv1);
	}

	// resorting section, does not apply if number of function increases

	// if linear dependencies occur, the mapping of input and output functions
	// might not be unique, i.e. the mapping might look like ipiv={0,1,1,3,4},
	// with double entries (1) and missing entries (2).
	// Imagine a He dimer with the bonding and antibonding orbitals both being mapped
	// to the left atomic orbital, just by numerical noise.
	// Find double and missing entries and replace the double with the missing ones

	std::vector<Function<T,NDIM> > result=zero_functions<T,NDIM>(world,vrhs.size());
	sirreps=std::vector<std::string>(vrhs.size(),"null");

	if (vrhs.size()==1) {
		result=result1;
		sirreps=sirreps1;

	} else {
		// find all double elements
		std::vector<int> ipiv2(ipiv1),double_elements;
		std::sort(ipiv2.begin(), ipiv2.end());
		for (size_t i=1; i<ipiv2.size(); ++i) {
			if (ipiv2[i]==ipiv2[i-1]) double_elements.push_back(ipiv2[i]);
		}

		if (verbosity_>1) print("double elements  ", double_elements);

		// if there are double elements we have to recompute the overlap matrix
		// to safely determine the mapping
		if (double_elements.size()>0) {
			Tensor<T> ovlp2=matrix_inner(world,result1,vrhs);
            Tensor<double> ovlp1=real(ovlp2);
            Tensor<double> imag_ovlp=imag(ovlp2);
            MADNESS_CHECK(imag_ovlp.normf()<1.e-10);


			long index[2];
			// find the largest overlap matrix element in each column
			for (int i=0; i<ovlp1.dim(1); ++i) {

				ovlp1(_,i).absmax(index);
				if (verbosity_>2) {
					print("result overlap matrix");
					print(ovlp1);
					print("index",index[0]);
				}
				ipiv1[i]=index[0];
				// remove this row from the tensor
				ovlp1(index[0],_)=0.0;
			}
		}
		if (verbosity_>1) print("final ipiv1",ipiv1);

		for (size_t i=0; i<ipiv1.size(); ++i) {
			result[ipiv1[i]]=result1[i];
			sirreps[ipiv1[i]]=sirreps1[i];
		}
	}	// end of resorting/mapping section

	if (not keep_ordering_) result=sort_to_irreps(result,sirreps);
	if (verbosity_>0) print("final function irreps: ",sirreps);

	return result;

}



std::vector<std::string> projector_irrep::reduce(const std::vector<std::string> irreps) const {

	bool verbose=false;
	for (const std::string& r : irreps) {
		if (r=="null") return std::vector<std::string> (1,"null");
	}

	// make sure all reps in the input vector exist in the table
	for (auto r : irreps) {
		if (table_.irreps_.find(r)==table_.irreps_.end()) {
			print_character_table();
			print("Irrep", r);
			MADNESS_EXCEPTION("Could not find irrep in table",1);
		}
	}
	if (verbose) {
		print("input irreps in reduce");
		for (auto& i : irreps) print(i);
	}

	// the characters of the direct product of the input representation
	// is the product of their characters
	std::vector<int> reducible(table_.irreps_.find(irreps[0])->second);
	for (std::size_t irep=1; irep<irreps.size(); ++irep) {
		std::vector<int> rep_characters(table_.irreps_.find(irreps[irep])->second);
		for (std::size_t ichar=0; ichar<rep_characters.size(); ++ichar) {
			reducible[ichar]*=rep_characters[ichar];
		}
	}
	if (verbose) {
		print("characters of the reducible product representation");
		for (const int i : reducible) print(i);
	}

	// reduce the reducible representation
	std::vector<std::string> result;
	for (const std::string& irrep : get_all_irreps()) {
		int n_irrep=0;
		std::vector<int> ii=table_.irreps_.find(irrep)->second;
		for (size_t i=0; i<reducible.size(); ++i) {	// sum over all classes/operators
			n_irrep+=reducible[i]*ii[i];
		}
		MADNESS_ASSERT(n_irrep%table_.order_==0);
		n_irrep/=table_.order_;
		if (verbose) {
			print("found irrep",irrep, n_irrep," times");
		}
		if (n_irrep>1) MADNESS_EXCEPTION("cannot handle the same irrep multiple times in reduce",1);
		if (n_irrep==1) result.push_back(irrep);
	}
	return result;
}


charactertable projector_irrep::make_c1_table() const {
	charactertable c1;
	c1.schoenflies_="C1";
	c1.order_=1;
	c1.operators_.push_back(pg_identity());
	c1.irreps_["a"]=vector_factory<int>(1);
	c1.mullikan_=vector_factory<std::string>("a");
	return c1;
}

charactertable projector_irrep::make_cs_table() const {
	charactertable cs;
	cs.schoenflies_="Cs";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_sigma_xy());
	cs.irreps_["a'"]=vector_factory<int>(1,1);
	cs.irreps_["a''"]=vector_factory<int>(1,-1);
	cs.mullikan_=vector_factory<std::string>("a'","a''");
	return cs;
}


charactertable projector_irrep::make_ci_table() const {
	charactertable cs;
	cs.schoenflies_="Ci";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_inversion());
	cs.mullikan_=vector_factory<std::string>("ag","au");
	cs.irreps_["ag"]=vector_factory<int>(1,1);
	cs.irreps_["au"]=vector_factory<int>(1,-1);
	return cs;
}

charactertable projector_irrep::make_c2_table() const {
	charactertable cs;
	cs.schoenflies_="C2";
	cs.order_=2;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.mullikan_=vector_factory<std::string>("a","b");
	cs.irreps_["a"]=vector_factory<int>(1,1);
	cs.irreps_["b"]=vector_factory<int>(1,-1);
	return cs;
}

charactertable projector_irrep::make_c2v_table() const {
	charactertable cs;
	cs.schoenflies_="C2v";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_sigma_xz());
	cs.operators_.push_back(pg_sigma_yz());
	cs.mullikan_=vector_factory<std::string>("a1","a2","b1","b2");
	cs.irreps_["a1"]=vector_factory<int>(1,1,1,1);
	cs.irreps_["a2"]=vector_factory<int>(1,1,-1,-1);
	cs.irreps_["b1"]=vector_factory<int>(1,-1,1,-1);
	cs.irreps_["b2"]=vector_factory<int>(1,-1,-1,1);
	return cs;
}


charactertable projector_irrep::make_c2h_table() const {
	charactertable cs;
	cs.schoenflies_="C2h";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_inversion());
	cs.operators_.push_back(pg_sigma_xy());
	cs.mullikan_=vector_factory<std::string>("ag","bg","au","bu");
	cs.irreps_["ag"]=vector_factory<int>(1, 1, 1, 1);
	cs.irreps_["bg"]=vector_factory<int>(1,-1, 1,-1);
	cs.irreps_["au"]=vector_factory<int>(1, 1,-1,-1);
	cs.irreps_["bu"]=vector_factory<int>(1,-1,-1, 1);
	return cs;
}

charactertable projector_irrep::make_d2_table() const {
	charactertable cs;
	cs.schoenflies_="D2";
	cs.order_=4;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_c2y());
	cs.operators_.push_back(pg_c2x());
	cs.mullikan_=vector_factory<std::string>("a1","b1","b2","b3");
	cs.irreps_["a1"]=vector_factory<int>(1, 1, 1, 1);
	cs.irreps_["b1"]=vector_factory<int>(1, 1,-1,-1);
	cs.irreps_["b2"]=vector_factory<int>(1,-1, 1,-1);
	cs.irreps_["b3"]=vector_factory<int>(1,-1,-1, 1);
	return cs;
}

charactertable projector_irrep::make_d2h_table() const {
	charactertable cs;
	cs.schoenflies_="D2h";
	cs.order_=8;
	cs.operators_.push_back(pg_identity());
	cs.operators_.push_back(pg_c2z());
	cs.operators_.push_back(pg_c2y());
	cs.operators_.push_back(pg_c2x());
	cs.operators_.push_back(pg_inversion());
	cs.operators_.push_back(pg_sigma_xy());
	cs.operators_.push_back(pg_sigma_xz());
	cs.operators_.push_back(pg_sigma_yz());
	cs.mullikan_=vector_factory<std::string>("ag","au","b1g","b1u","b2g","b2u","b3g","b3u");
	cs.irreps_["ag"] =vector_factory<int>(1, 1, 1, 1, 1, 1, 1, 1);
	cs.irreps_["b1g"]=vector_factory<int>(1, 1,-1,-1, 1 ,1,-1,-1);
	cs.irreps_["b2g"]=vector_factory<int>(1,-1, 1,-1, 1,-1, 1,-1);
	cs.irreps_["b3g"]=vector_factory<int>(1,-1,-1, 1, 1,-1,-1, 1);
	cs.irreps_["au"] =vector_factory<int>(1, 1, 1, 1,-1,-1,-1,-1);
	cs.irreps_["b1u"]=vector_factory<int>(1, 1,-1,-1,-1,-1, 1, 1);
	cs.irreps_["b2u"]=vector_factory<int>(1,-1, 1,-1,-1, 1,-1, 1);
	cs.irreps_["b3u"]=vector_factory<int>(1,-1,-1, 1,-1, 1, 1,-1);
	return cs;
}


// explicit instantiation
template std::vector<Function<double,1> > projector_irrep::project_on_irreps (
		const std::vector<Function<double,1> >& vhrs, const std::vector<std::string>& irreps) const;
template std::vector<Function<double,2> > projector_irrep::project_on_irreps (
		const std::vector<Function<double,2> >& vhrs, const std::vector<std::string>& irreps) const;
template std::vector<Function<double,3> > projector_irrep::project_on_irreps (
		const std::vector<Function<double,3> >& vhrs, const std::vector<std::string>& irreps) const;


template std::vector<Function<double,1> >
projector_irrep::apply_symmetry_operators(
		const std::vector<Function<double,1> >& vrhs,
		Function<double,1> metric,
		std::vector<std::string>& sirreps) const;
template std::vector<Function<double,2> >
projector_irrep::apply_symmetry_operators(
		const std::vector<Function<double,2> >& vrhs,
		Function<double,2> metric,
		std::vector<std::string>& sirreps) const;
template std::vector<Function<double,3> >
projector_irrep::apply_symmetry_operators(
		const std::vector<Function<double,3> >& vrhs,
		Function<double,3> metric,
		std::vector<std::string>& sirreps) const;

template std::vector<Function<double_complex,3> >
projector_irrep::apply_symmetry_operators(
        const std::vector<Function<double_complex,3> >& vrhs,
        Function<double,3> metric,
        std::vector<std::string>& sirreps) const;

//template void Class::function(int);
//template class Kinetic<double,2>;
//template class Kinetic<double,3>;
//template class Kinetic<double,4>;
//template class Kinetic<double,5>;
//template class Kinetic<double,6>;
//reproject_on_irrep

} /* namespace madness */
