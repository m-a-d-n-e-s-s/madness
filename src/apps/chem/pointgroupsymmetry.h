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
	/// ctor takes the point group and the optionally the irrep
	projector_irrep(std::string pointgroup, std::string irrep="") {
        std::transform(pointgroup.begin(), pointgroup.end(), pointgroup.begin(), ::tolower);
		table_=get_table(pointgroup);
		set_irrep(irrep);
	}

	/// set the irrep on which this projector projects
	void set_irrep(std::string irrep) {
        std::transform(irrep.begin(), irrep.end(), irrep.begin(), ::tolower);

        // check consistency of the input
		if (std::find(table_.mullikan_.begin(), table_.mullikan_.end(), irrep) == table_.mullikan_.end()) {
			if (irrep_!="all") {
				print("irrep",irrep);
				print("point group", table_.schoenflies_);
				MADNESS_EXCEPTION("no such irrep in this point group",1);
			}
		}

		irrep_=irrep;
	}

	std::vector<std::string> get_all_irreps() const {
		return table_.mullikan_;
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const Function<T,NDIM>& rhs) const {

		std::vector<Function<T,NDIM> > vrhs;
		vrhs.push_back(rhs);
		return this->operator()(vrhs);

	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs) const {

		typedef std::vector<std::vector<Function<T,NDIM> > > vecvecfunctionT;
		vecvecfunctionT result=apply_symmetry_operators(vrhs,true);
	}


	/// apply symmetry operators on a vector of input functions

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM> std::vector<std::vector<Function<T,NDIM> > >
	apply_symmetry_operators(const std::vector<Function<T,NDIM> >& vrhs, bool fence=true) const {

		// apply the symmetry operators to rhs
		std::vector<std::vector<Function<T,NDIM> > > result(table_.operators_.size());

		if (vrhs.size()==0) return result;
		World& world=vrhs.begin()->world();

		// get all irreps to work on
		std::vector<std::string> all_irreps;
		if (irrep_=="all") {
			all_irreps=table_.mullikan_;
		} else {
			all_irreps.push_back(irrep_);
		}

		// loop over all irreps
		for (auto& irrep : table_.irreps_) {
			// get the irrep line
			const charactertable::characterlineT cline=irrep.second;
			MADNESS_ASSERT(table_.operators_.size()==cline.size());

			for (int i=0; i<cline.size(); ++i) {
				double character=double(cline[i]);
				const pg_operator& op=table_.operators_[i];
				result[i]=(character/table_.order_)*op(vrhs,false);
			}
		}
		if (fence) world.gop.fence();
		return result;
	}

	/// project a number of functions on a given irrep

	/// linear dependencies are removed, the functions are orthonormalized
	/// @param[in]				vrhs	vector of functions to be projected on the irrep
	/// @param[in, optional]	metric	additional metric for overlap and orthonormalization
	/// @param[in, optional]	lindep	linear dependency threshold
	/// @param[in, optional]	verbose	verbosity level
//	template<typename T, std::size_t NDIM>
//	std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vrhs,
//			const Function<T,NDIM> metric=Function<T,NDIM>(),
//			const double lindep=1.e-3, bool verbose=true) const {
//
//		if (vrhs.size()==0) return std::vector<Function<T,NDIM> >();
//		World& world=vrhs[0].world();
//
//		// project input function onto the irrep
//		vector_real_function_3d a1, result(vrhs.size());
//		for (const real_function_3d& f : vrhs) {
//			a1.push_back(operator()(f));
//		}
//
//		// compute overlap, include metric if necessary
//		Tensor<double> ovlp;
//		if (metric.is_initialized()) ovlp=matrix_inner(world,a1,metric*a1);
//		else ovlp=matrix_inner(world,a1,a1);
//
//		if (verbose) {
//			print("ovlp");
//			print(ovlp);
//		}
//
//		// check for zero rank
//		Tensor<double> evec,eval;
//		syev(ovlp,evec,eval);
//
//		if (eval.sum()>1.e-3) {
//
//			if (verbose) print("eigenvalues",eval);
//			print("eigenvalues",eval);
//
//			// compute rank-revealing Cholesky, discard linear dependent vectors
//			Tensor<int> piv;
//			int rank;
//			rr_cholesky(ovlp,lindep,piv,rank);
//			Tensor<double> shrunk=copy(ovlp(Slice(0,rank-1,1),Slice(0,rank-1,1)));
//			Tensor<double> inv=inverse(shrunk);
//
//			if (verbose) {
//				print("cholesky");
//				print(ovlp);
//				print("inverse(cholesky)");
//				print(inv);
//				print("orthogonality");
//				print(inner(shrunk,shrunk,0,0));
//				print("pivot elements");
//				print(piv);
//			}
//
//			// transform linearly independent orbitals only
//			vector_real_function_3d a1shrunk;
//			for (int i=0; i<rank; ++i) a1shrunk.push_back(a1[piv(i)]);
//			vector_real_function_3d tmp=transform(world,a1shrunk,inv);
//
//			// keep original ordering
//			for (int i=0; i<a1shrunk.size(); ++i) {
//				result[piv(i)]=tmp[i];
//			}
//		}
//
//		return result;
//	}

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

	std::map<std::string, charactertable> charactertables_;
	charactertable table_;

	/// choose one of the irreps or "all"
	std::string irrep_;

	/// after projection: result functions being ordered according to irreps or ordering unchanged
	bool irrep_ordering;

	charactertable make_cs_table() const;

	charactertable make_ci_table() const;

	charactertable make_c2_table() const;

	charactertable make_c2v_table() const;

	/// return the character table according to the requested point group
	const charactertable& get_table(std::string pointgroup) {

        std::transform(pointgroup.begin(), pointgroup.end(), pointgroup.begin(), ::tolower);
		std::map<std::string, charactertable>::const_iterator it=charactertables_.find(pointgroup);

		if (it==charactertables_.end()) {
			if (pointgroup=="ci") charactertables_[pointgroup]=make_ci_table();
			if (pointgroup=="cs") charactertables_[pointgroup]=make_cs_table();
			if (pointgroup=="c2") charactertables_[pointgroup]=make_c2_table();
			if (pointgroup=="c2v") charactertables_[pointgroup]=make_c2v_table();
			it=charactertables_.find(pointgroup);
		}
		MADNESS_ASSERT(it!=charactertables_.end());
		return it->second;
	}
};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_ */
