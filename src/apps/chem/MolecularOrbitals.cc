/*
 * MolecularOrbitals.cpp
 *
 *  Created on: 11 Jul 2019
 *      Author: fbischoff
 */

#include <chem/MolecularOrbitals.h>
#include <chem/molecule.h>
#include <madness/world/world.h>
#include <madness/world/parallel_archive.h>
#include <chem/molecularbasis.h>

#include<madness/mra/mra.h>

namespace madness {


template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::write_to(std::vector<Function<T,NDIM> >& mo_out, Tensor<double>& eps_out,
		std::vector<std::string>& irrep_out, Tensor<double>& occ_out, std::vector<int>& set_out) const {
	if (mo.size()>0) {
		World& world=mo.front().world();
		mo_out=copy(world,mo);
		eps_out=eps;
		irrep_out=irreps;
		occ_out=occ;
		set_out=localize_sets;
	} else {
		mo_out.clear();
		eps_out.clear();
		irrep_out.clear();
		occ_out.clear();
		set_out.clear();
	}
}

template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::post_process_mos(World& world, const double thresh, const int k) {
	if (mo[0].k() != k) {
		reconstruct(world, mo);
		for (unsigned int i = 0; i < mo.size(); ++i)
			mo[i] = madness::project(mo[i], k, thresh, false);
		world.gop.fence();
	}
	set_thresh(world,mo,thresh);
}


template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::project_ao(World& world, const Tensor<T>& Saomo,
		const std::vector<real_function_3d>& aos) {

	Tensor<double> Saoao=matrix_inner(world, aos, aos);

	Tensor<T> c;
	gesvp(world, convert<T>(Saoao), Saomo, c);
	mo = transform(world, aos, c, 0.0, true);
	truncate(world, mo);
	mo=orthonormalize_symmetric(mo);

}

template class MolecularOrbitals<double,3>;
template class MolecularOrbitals<double_complex,3>;



} /* namespace madness */
