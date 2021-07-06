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
#include <chem/SCF.h>
#include <chem/pointgroupsymmetry.h>

#include<madness/mra/mra.h>

namespace madness {

template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::recompute_irreps(const std::string pointgroup,
                                                 const Function<typename Tensor<T>::scalar_type,NDIM>& metric) {
    projector_irrep symmetry_projector(pointgroup);
    mo=symmetry_projector(mo,metric,irreps);
    print("in recompute_irreps",irreps);
}


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



/// save MOs in the AO projection for geometry restart

/// compute the aos in MRA projection as:
/// std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, calc.aobasis, calc.molecule);
template<typename T, std::size_t NDIM>
void MolecularOrbitals<T,NDIM>::save_restartaodata(World& world, const Molecule& molecule,
		const MolecularOrbitals<T,NDIM>& amo, const MolecularOrbitals<T,NDIM>& bmo,
		const AtomicBasisSet& aobasis) {

	std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, aobasis,molecule);

	Tensor<T> Saoamo = matrix_inner(world, aos, amo.get_mos());
	Tensor<T> Saobmo = (bmo.get_mos().size()>0) ? matrix_inner(world, aos, bmo.get_mos()) : Tensor<T>();
	if (world.rank() == 0) {
		archive::BinaryFstreamOutputArchive arao("restartaodata");
		arao << aobasis << Saoamo << amo.get_eps() << amo.get_occ() << amo.get_localize_sets();
		if (Saobmo.size()>0) arao << Saobmo << bmo.get_eps() << bmo.get_occ() << bmo.get_localize_sets();
	}
}


/// uses AO-projection as a restart guess

/// @return amo and bmo
template<typename T, std::size_t NDIM>
std::pair<MolecularOrbitals<T,NDIM>, MolecularOrbitals<T,NDIM> >
MolecularOrbitals<T,NDIM>::read_restartaodata(World& world,
		const Molecule& molecule, const bool have_beta) {

	archive::BinaryFstreamInputArchive arao("restartaodata");

	Tensor<T> Saomo;
	Tensor<double> eps, occ;
	std::vector<int> localize_set;
	std::vector<std::string> irrep;
	std::vector<Function<T,NDIM> > dummy_mo;
	AtomicBasisSet aobasis;

	MolecularOrbitals<T,NDIM> amo, bmo;

	try {
		// read alpha orbitals
		arao >> aobasis >> Saomo >> eps >> occ >> localize_set;
		print("reading from aobasis",aobasis.get_name());
		std::vector<Function<double,3> > aos1=SCF::project_ao_basis_only(world, aobasis,molecule);

		amo=MolecularOrbitals<T,NDIM>(dummy_mo,eps,irrep,occ,localize_set);
		amo.project_ao(world,Saomo,aos1);

		if (have_beta) {
			arao >> Saomo >> eps >> occ >> localize_set;
			bmo=MolecularOrbitals<T,NDIM>(dummy_mo,eps,irrep,occ,localize_set);
			bmo.project_ao(world,Saomo,aos1);
		}

	} catch (...) {
		throw std::runtime_error("failed to read restartdata");
	}
	return std::make_pair(amo,bmo);
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
