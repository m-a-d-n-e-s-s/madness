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

#include<madness/mra/mra.h>

namespace madness {

template class MolecularOrbitals<double,3>;

template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::write_to(std::vector<Function<T,NDIM> >& mo_out, std::vector<double>& eps_out,
		std::vector<std::string>& irrep_out, std::vector<double>& occ_out, std::vector<int>& set_out) const {
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


} /* namespace madness */
