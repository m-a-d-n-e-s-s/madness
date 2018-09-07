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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file chem/pointgroupsymmetry.cc
  \brief implements point group operations

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/chem/pointgroupsymmetry.cc>here</a>.

*/
#ifndef SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_
#define SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/functypedefs.h>
#include <chem/pointgroupoperator.h>

using namespace madness;

namespace madness {


struct charactertable {
	typedef std::vector<int> characterlineT;

	std::string schoenflies_;					///< Schoenflies symbol of the point group
	int order_;									///< order of the point group
	std::map<std::string,characterlineT> irreps_;
	std::vector<pg_operator> operators_;		///< symmetry operators
	std::vector<std::string> mullikan_;			///< Mullikan symbols of the irreps

};

class projector_irrep {

public:

	/// default ctor
	projector_irrep() = default;

	/// ctor takes the point group and the optionally the irrep
	projector_irrep(std::string pointgroup, std::string irrep="all")
		: lindep_(1.e-3), verbosity_(0), keep_ordering_(true) {
        std::transform(pointgroup.begin(), pointgroup.end(), pointgroup.begin(), ::tolower);
		table_=get_table(pointgroup);
		set_irrep(irrep);
	}

	/// set the irrep on which this projector projects
	void set_irrep(std::string irrep) {
        std::transform(irrep.begin(), irrep.end(), irrep.begin(), ::tolower);

        // check consistency of the input
		if (std::find(table_.mullikan_.begin(), table_.mullikan_.end(), irrep) == table_.mullikan_.end()) {
			if (irrep!="all") {
				print("irrep",irrep);
				print("point group", table_.schoenflies_);
				MADNESS_EXCEPTION("no such irrep in this point group",1);
			}
		}

		irrep_=irrep;
	}

	/// set the linear dependency threshold
	void set_lindep(const double ld) {lindep_=ld;}

	/// set the verbosity level
	void set_verbosity(int v) {verbosity_=v;}

	/// get the verbosity level
	int get_verbosity() {return verbosity_;}

	/// get the point group name
	std::string get_pointgroup() const {return table_.schoenflies_;}

	/// return the character table
	charactertable get_table() const {return table_;}

	/// set the ordering after symmetrization: irreps or keep as is
	void set_ordering(std::string o) {
		if (o=="irrep") keep_ordering_=false;
		else if (o=="keep") keep_ordering_=true;
		else {
			print("projector_irrep: ordering parameter unknown: ",o);
			MADNESS_EXCEPTION("projector_irrep: ordering parameter unknown",1);
		}
	}

	std::vector<std::string> get_all_irreps() const {
		return table_.mullikan_;
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(const Function<T,NDIM>& rhs) const {

		std::vector<Function<T,NDIM> > vrhs;
		vrhs.push_back(rhs);
		return this->operator()(vrhs);

	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs) const {

		Function<T,NDIM> metric;
		std::vector<std::string> sirrep;
		return apply_symmetry_operators(vrhs,metric,sirrep,false);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			const Function<T,NDIM>& metric) const {

		std::vector<std::string> sirrep;
		return apply_symmetry_operators(vrhs,metric,sirrep,false);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			const Function<T,NDIM>& metric,
			std::vector<std::string>& sirreps) const {

		return apply_symmetry_operators(vrhs,metric,sirreps,false);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			std::vector<std::string>& sirreps) const {

		Function<T,NDIM> metric;
		return apply_symmetry_operators(vrhs,metric,sirreps,false);
	}

	/// create a symmetry-adapted basis from a single function

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > create_symmetry_adapted_basis(
			const Function<T,NDIM>& rhs,
			const Function<T,NDIM>& metric,
			std::vector<std::string>& sirreps) {

		std::vector<Function<T,NDIM> > vrhs(1,rhs);
		return apply_symmetry_operators(vrhs,metric,sirreps,true);

	}

	/// create a symmetry-adapted basis from a single function

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > create_symmetry_adapted_basis(
			const Function<T,NDIM>& rhs,
			std::vector<std::string>& sirreps) {

		std::vector<Function<T,NDIM> > vrhs(1,rhs);
		Function<T,NDIM> metric;
		return apply_symmetry_operators(vrhs,metric,sirreps,true);

	}



	/// print the character table
	void print_character_table() const {
		print("character table for point group ",table_.schoenflies_);
		printf("       ");
		for (auto& op : table_.operators_) {
			printf(" %5s ",op.symbol().c_str());
		}
		printf("\n");
		for (auto& mullikan : table_.mullikan_) {
			printf(" %5s ",mullikan.c_str());
			std::vector<int> characters=table_.irreps_.find(mullikan)->second;
			for (auto& character : characters) {
				printf(" %5i ",character);
			}
			printf("\n");
		}
	}

	/// reduce a reducible representation

	/// @param[in]	reps	vector or irrep or reduc. reps whose product will be reduced
	/// @return		vector of irreps constituting the input product
	std::vector<std::string> reduce(const std::vector<std::string> reps) const;

	/// reduce a product of two irreps
	std::vector<std::string> reduce(const std::string irrep1, const std::string irrep2) const {
		return reduce(vector_factory<std::string>(irrep1,irrep2));
	}

	/// reduce a product of three irreps
	std::vector<std::string> reduce(const std::string irrep1, const std::string irrep2,
			const std::string irrep3) const {
		return reduce(vector_factory<std::string>(irrep1,irrep2,irrep3));
	}

	/// reduce a product of four irreps
	std::vector<std::string> reduce(const std::string irrep1, const std::string irrep2,
			const std::string irrep3, const std::string irrep4) const {
		return reduce(vector_factory<std::string>(irrep1,irrep2,irrep3,irrep4));
	}

private:

	charactertable table_;

	/// choose one of the irreps or "all"
	std::string irrep_;

	/// linear dependency threshold for the orthogonalization
	double lindep_;

	/// verbosity level
	int verbosity_;

	/// after projection: result functions being ordered according to irreps or ordering unchanged
	bool keep_ordering_;

	charactertable make_c1_table() const;

	charactertable make_cs_table() const;

	charactertable make_ci_table() const;

	charactertable make_c2_table() const;

	charactertable make_c2v_table() const;

	charactertable make_c2h_table() const;

	charactertable make_d2_table() const;

	charactertable make_d2h_table() const;


	/// return the character table according to the requested point group
	const charactertable get_table(std::string pointgroup) {

        std::transform(pointgroup.begin(), pointgroup.end(), pointgroup.begin(), ::tolower);

        charactertable table;
        if (pointgroup=="c1") table=make_c1_table();
        else if (pointgroup=="ci") table=make_ci_table();
        else if (pointgroup=="cs") table=make_cs_table();
        else if (pointgroup=="c2") table=make_c2_table();
        else if (pointgroup=="c2v") table=make_c2v_table();
        else if (pointgroup=="c2h") table=make_c2h_table();
        else if (pointgroup=="d2") table=make_d2_table();
        else if (pointgroup=="d2h") table=make_d2h_table();
        else {
        	MADNESS_EXCEPTION("unknown group table",1);
        }
		return table;
	}

	/// symmetrize a vector of functions

	/// project a number of functions on a given number of irreps
	/// linear dependencies are removed, the functions are orthonormalized
	/// the number of input and output function is expected to be the same!
	/// @param[in]	vrhs	vector of functions to be projected on the irrep
	/// @param[in]	vbra	bra of vrhs if applicable (if bra /= ket), may be empty
	/// @param[in]	allow_increase	allow increased number of output functions comp. to input
	/// @param[out]	sirrep	vector with the irrep names corresponding to the result
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > apply_symmetry_operators(
			const std::vector<Function<T,NDIM> >& vrhs,
			Function<T,NDIM> metric,
			std::vector<std::string>& sirreps,
			bool allow_increase,
			const bool fence=true) const {

		World& world=vrhs.begin()->world();

		// fast return
		if (vrhs.size()==0) return std::vector<Function<T,NDIM> >();
		if (table_.schoenflies_=="C1") {return copy(world,vrhs);}


		// loop over all symmetry operations of this point group and apply them
		// to the input vector; dimension: (nop, nmo)
		std::vector<std::vector<Function<T,NDIM> > > opvrhs(table_.operators_.size());
		for (int i=0; i<table_.operators_.size(); ++i) {
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
		for (int j=0; j<all_irreps.size(); ++j) {
			std::string& irrep=all_irreps[j];

			// get the characters of this irrep
			const charactertable::characterlineT cline=table_.irreps_.find(irrep)->second;
			MADNESS_ASSERT(table_.operators_.size()==cline.size());


			// do the linear combination; loop over all operators
			// lc_op_vhrs[jrrep] = \sum_iop (character[iop]/table_.order_)*opvhrs[iop];
			for (int iop=0; iop<cline.size(); ++iop) {
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
		int iresult=0;	// counter for the result vector ordering
		for (int j=0; j<all_irreps.size(); ++j) {

			if (verbosity_>1) print("working on irrep ",all_irreps[j]);
			// compute overlap, include metric if necessary
			Tensor<double> ovlp;
			if (metric.is_initialized()) ovlp=matrix_inner(world,lc_op_vrhs[j],metric*lc_op_vrhs[j]);
			else ovlp=matrix_inner(world,lc_op_vrhs[j],lc_op_vrhs[j]);

			for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)+=1.e-6;

			if (verbosity_>1) {
				print("ovlp");
				print(ovlp);
			}

			// check for zero rank
			Tensor<double> evec,eval;
			syev(ovlp,evec,eval);

			if (eval.sum()>1.e-3) {

				if (verbosity_>1) print("eigenvalues",eval);

				// compute rank-revealing Cholesky, discard linear dependent vectors
				Tensor<int> piv;
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

				// transform linearly independent orbitals only
				std::vector<Function<T,NDIM> > a1shrunk;
				for (int i=0; i<rank; ++i) a1shrunk.push_back(lc_op_vrhs[j][piv(i)]);
				std::vector<Function<T,NDIM> > tmp=transform(world,a1shrunk,inv);

				// collect all result functions and their irreps
				for (int i=0; i<a1shrunk.size(); ++i) {
					result1.push_back(tmp[i]);
					sirreps1.push_back(all_irreps[j]);
					ipiv1.push_back(piv[i]);
				}
			}
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

		std::vector<Function<T,NDIM> > result(result1.size());
		sirreps.resize(result1.size());

		if (vrhs.size()==1) {
			result=result1;
			sirreps=sirreps1;

		} else {
			// find all double elements
			std::vector<int> ipiv2(ipiv1),double_elements;
			std::sort(ipiv2.begin(), ipiv2.end());
			for (int i=1; i<ipiv2.size(); ++i) {
				if (ipiv2[i]==ipiv2[i-1]) double_elements.push_back(ipiv2[i]);
			}

			if (verbosity_>1) print("double elements  ", double_elements);

			// if there are double elements we have to recompute the overlap matrix
			// to safely determine the mapping
			if (double_elements.size()>0) {
				Tensor<double> ovlp1=matrix_inner(world,result1,vrhs);

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

			for (int i=0; i<ipiv1.size(); ++i) {
				result[ipiv1[i]]=result1[i];
				sirreps[ipiv1[i]]=sirreps1[i];
			}
		}	// end of resorting/mapping section

		if (not keep_ordering_) result=sort_to_irreps(result,sirreps);
		if (verbosity_>0) print("final function irreps: ",sirreps);

		return result;

	}

	/// given a vector of functions and their irreps, extract only the ones of the requested irrep
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > sort_to_irreps(std::vector<Function<T,NDIM> >& vrhs,
			std::vector<std::string>& sirreps) const {

		std::map<std::string,std::vector<Function<T,NDIM> > > m;
		for (int i=0; i<vrhs.size(); ++i) {
			m[sirreps[i]].push_back(vrhs[i]);
		}

		std::vector<Function<T,NDIM> > result;
		sirreps.clear();
		for (const std::string& irrep : table_.mullikan_) {
			if (m.find(irrep)!=m.end()) {
				// concatenate the functions
				std::vector<Function<T,NDIM> > tmp=m.find(irrep)->second;
				result.insert(result.end(),tmp.begin(),tmp.end());
				// concatenate the strings
				std::vector<std::string> s(tmp.size(),irrep);
				sirreps.insert(sirreps.end(),s.begin(),s.end());
			}
		}
		return result;
	}

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_ */
