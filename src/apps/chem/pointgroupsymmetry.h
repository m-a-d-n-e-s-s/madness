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
  /trunk/src/apps/chem/pointgroupsymmetry.h>here</a>.

	\par Use cases
	There are 3 use cases for the projector
	 - create_symmetry_adapted_basis()
	   Given a single function generate a set of symmetry-adapted functions of a given irrep
	 - operator()()
	   Given a set of functions transform them to a set of symmetry-adapted function.
	   This is most likely the case for a localized-to-symmetrized transformation of orbitals.
	   The irreps will be determined by the projection based on a rank-revealing Cholesky
	   decomposition. It is therefore necessary that the set of input functions is normalized,
	   so that small eigenvalues are always due to numerics and not to physics.
	   There will be a rough check on closedness, i.e. that the old and the new set of functions
	   project on the same space.
	 - project_on_irreps()
	   Given a set of functions and a set of corresponding irreps, reproject the functions on
	   their irreps
	 Furthermore a reduction method is provided to do some algebra with the irreps

*/
#ifndef SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_
#define SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_

#include <madness/world/MADworld.h>
#include <chem/pointgroupoperator.h>

using namespace madness;

namespace madness {

	template<typename T, std::size_t NDIM>
	class Function;
}



namespace madness {
struct charactertable {
	typedef std::vector<int> characterlineT;

	std::string schoenflies_;					///< Schoenflies symbol of the point group
	int order_;									///< order of the point group
	std::map<std::string,characterlineT> irreps_;
	std::vector<pg_operator> operators_;		///< symmetry operators
	std::vector<std::string> mullikan_;			///< Mullikan symbols of the irreps
	bool is_abelian=true;						///< currently always true

};

/*!
  \file projector_irrep.h
  \brief Defines and implements the symmetry projector class
  \ingroup chem
  \addtogroup chem

  \par Introduction

*/

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

	/// print the parameters of this projector
	void print_info(World& world) const {
		if (world.rank()==0) {
			print("symmetry projector for point group ",table_.schoenflies_);
			print("");
			print("             working on irrep:", irrep_);
			print("  linear dependency threshold:", lindep_);
			print("                keep ordering:", keep_ordering_);
			print("        orthonormalize irreps:", orthonormalize_irreps_);
			print("              verbosity level:", verbosity_);
			print("");
			print_character_table();
		}
	}

	/// set the irrep on which this projector projects
	projector_irrep& set_irrep(std::string irrep) {
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
		return *this;
	}

	/// set the linear dependency threshold
	projector_irrep& set_lindep(const double ld) {
		lindep_=ld;
		return *this;
	}

	/// set the verbosity level
	projector_irrep& set_verbosity(int v) {
		verbosity_=v;
		return *this;
	}

	/// get the verbosity level
	int get_verbosity() const {return verbosity_;}

	/// get the verbosity level
	projector_irrep& set_orthonormalize_irreps(bool flag) {
		orthonormalize_irreps_=flag;
		return *this;
	}

	/// get the verbosity level
	bool get_orthonormalize_irreps() const {return orthonormalize_irreps_;}

	/// get the point group name
	std::string get_pointgroup() const {return table_.schoenflies_;}

	/// return the character table
	charactertable get_table() const {return table_;}

	/// set the ordering after symmetrization: irreps or keep as is
	projector_irrep& set_ordering(const std::string o) {
		if (o=="irrep") keep_ordering_=false;
		else if (o=="keep") keep_ordering_=true;
		else {
			print("projector_irrep: ordering parameter unknown: ",o);
			MADNESS_EXCEPTION("projector_irrep: ordering parameter unknown",1);
		}
		return *this;
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

        Function<typename Tensor<T>::scalar_type,NDIM> metric;
		std::vector<std::string> sirrep;
		return apply_symmetry_operators(vrhs,metric,sirrep);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			const Function<typename Tensor<T>::scalar_type,NDIM>& metric) const {

		std::vector<std::string> sirrep;
		return apply_symmetry_operators(vrhs,metric,sirrep);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			const Function<typename Tensor<T>::scalar_type,NDIM>& metric,
			std::vector<std::string>& sirreps) const {

		return apply_symmetry_operators(vrhs,metric,sirreps);
	}

	/// projector on a given irrep

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > operator()(
			const std::vector<Function<T,NDIM> >& vrhs,
			std::vector<std::string>& sirreps) const {

        Function<typename Tensor<T>::scalar_type,NDIM> metric;
		return apply_symmetry_operators(vrhs,metric,sirreps);
	}

	/// create a symmetry-adapted basis from a single function

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > create_symmetry_adapted_basis(
			const Function<T,NDIM>& rhs,
			const Function<typename Tensor<T>::scalar_type,NDIM>& metric,
			std::vector<std::string>& sirreps) {

		std::vector<Function<T,NDIM> > vrhs(1,rhs);
		return apply_symmetry_operators(vrhs,metric,sirreps);

	}

	/// create a symmetry-adapted basis from a single function

	/// @return	a vector[irreps] of a vector of functions
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > create_symmetry_adapted_basis(
			const Function<T,NDIM>& rhs,
			std::vector<std::string>& sirreps) {

		std::vector<Function<T,NDIM> > vrhs(1,rhs);
        Function<typename Tensor<T>::scalar_type,NDIM> metric;
		return apply_symmetry_operators(vrhs,metric,sirreps);

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

	/// (re-)project the argument on the given irreps

	/// @param[in]	vrhs	the vector of functions to be projected on the irreps given by irreps
	/// @param[in]	irreps	the irreps in order of the functions
	/// @return		vrhs projected on the irreps
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > project_on_irreps(const std::vector<Function<T,NDIM> >& vhrs,
			const std::vector<std::string>& irreps) const;

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

	/// get a mapping canonical to irrep-sorting
	std::vector<int> get_canonical_to_irrep_map(std::vector<std::string> sirreps) const {

		// distribute into bins
		std::map<std::string,std::vector<int> > m;
		for (size_t i=0; i<sirreps.size(); ++i) m[sirreps[i]].push_back(i);

		// concatenate bins
		std::vector<int> map;
		for (auto& irrep : table_.mullikan_) {
			map.insert(end(map), begin(m[irrep]), end(m[irrep]));
		}
		return map;
	}

	/// given a mapping resort
	template<typename R>
	static std::vector<R> resort(const std::vector<int> map,
			const std::vector<R>& vrhs) {
		MADNESS_ASSERT(map.size()==vrhs.size());
		std::vector<R> result(vrhs.size());
		for (size_t i=0; i<map.size(); ++i) result[i]=vrhs[map[i]];
		return result;
	}

	/// given a mapping resort
	template<typename R>
	static std::vector<R> reverse_resort(const std::vector<int> map,
			const std::vector<R>& vrhs) {
		MADNESS_ASSERT(map.size()==vrhs.size());
		std::vector<R> result(vrhs.size());
		for (int i=0; i<map.size(); ++i) result[map[i]]=vrhs[i];
		return result;
	}

private:

	charactertable table_;

	/// choose one of the irreps or "all"
	std::string irrep_="all";

	/// linear dependency threshold for the orthogonalization
	double lindep_=1.e-3;

	/// verbosity level
	int verbosity_=1;

	/// after projection: result functions being ordered according to irreps or ordering unchanged
	bool keep_ordering_=true;

	/// orthonormalize within the irreps or simply discard linear dependent vectors
	bool orthonormalize_irreps_=true;

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
	/// vanishing input functions are mapped onto vanishing output functions with irrep "null"
	/// the number of input and output function is expected to be the same!
	/// @param[in]	vrhs	vector of functions to be projected on the irrep
	/// @param[in]	vbra	bra of vrhs if applicable (if bra /= ket), may be empty
	/// @param[out]	sirrep	vector with the irrep names corresponding to the result
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > apply_symmetry_operators(
			const std::vector<Function<T,NDIM> >& vrhs,
            Function<typename Tensor<T>::scalar_type,NDIM> metric,
			std::vector<std::string>& sirreps) const;

	/// sort the functions according to their irreps
	template<typename T, std::size_t NDIM>
	std::vector<Function<T,NDIM> > sort_to_irreps(std::vector<Function<T,NDIM> >& vrhs,
			std::vector<std::string>& sirreps) const {
		std::vector<int> map=get_canonical_to_irrep_map(sirreps);
		return resort(map,vrhs);
	}

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_POINTGROUPSYMMETRY_H_ */
