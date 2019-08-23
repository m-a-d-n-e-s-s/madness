/*
 * MolecularOrbitals.h
 *
 *  Created on: 11 Jul 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_MOLECULARORBITALS_H_
#define SRC_APPS_CHEM_MOLECULARORBITALS_H_

#include<vector>
#include<madness/tensor/tensor.h>
#include <madness/world/parallel_archive.h>
#include <chem/molecule.h>

namespace madness {

template<typename T, std::size_t NDIM> class Function;
class World;
class AtomicBasisSet;

template<typename T, std::size_t NDIM>
class MolecularOrbitals {
public:

	MolecularOrbitals() {}

	MolecularOrbitals(const std::vector<Function<T,NDIM> >& mo, const Tensor<double>& eps,
			const std::vector<std::string>& irrep, const Tensor<double>& occ,
			const std::vector<int>& set)
			: mo(mo), eps(eps), irreps(irrep), occ(occ), localize_sets(set) {
	}

	std::vector<Function<T,NDIM> > get_mos() const {
		return mo;
	}

	Tensor<double> get_eps() const {
		return eps;
	}

	std::vector<std::string> get_irreps() const {
		return irreps;
	}

	std::vector<int> get_localize_sets() const {
		return localize_sets;
	}

	Tensor<double> get_occ() const {
		return occ;
	}

	/// setters will always invalidate all other member variables
	void set_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		invalidate_all();
		mo=mo_new;
	}

	/// updates will keep other member variables
	void update_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		mo=mo_new;
	}

	void update_occ(const Tensor<double>& occ_new) {
		occ=occ_new;
	}

	/// updates will keep other member variables
	void update_mos_and_eps(const std::vector<Function<T,NDIM> >& mo_new,
			const Tensor<double>& eps_new) {
		mo=mo_new;
		eps=copy(eps_new);
	}

	void recompute_irreps() {
	}

	void recompute_localize_sets() {
	}

	void invalidate_all() {
		invalidate_mos();
		invalidate_eps();
		invalidate_irreps();
		invalidate_occ();
		invalidate_localize_sets();
	}

	void invalidate_mos() {
		mo.clear();
	}

	void invalidate_eps() {
		eps.clear();
	}

	void invalidate_irreps() {
		irreps.clear();
	}

	void invalidate_occ() {
		occ.clear();
	}

	void invalidate_localize_sets() {
		localize_sets.clear();
	}

	template <typename Archive>
	void serialize (Archive& ar) {
		ar & mo & eps & irreps & occ & localize_sets;
	}

	friend bool similar(const MolecularOrbitals& mo1, const MolecularOrbitals& mo2, const double thresh=1.e-6) {

		if (mo1.mo.size()!=mo2.mo.size()) return false;
		if (mo1.mo.size()==0) return true;

		World& world=mo1.mo.front().world();
		bool similar=((mo1.eps-mo2.eps).normf()<thresh);
		similar=similar and (norm2(world,mo1.mo-mo2.mo)<thresh);
		similar=similar and (mo1.irreps==mo2.irreps);
		similar=similar and (mo1.localize_sets==mo2.localize_sets);
		return similar;
	}

	void write_to(std::vector<Function<T,NDIM> >& mo_out, Tensor<double>& eps_out,
			std::vector<std::string>& irrep_out, Tensor<double>& occ_out, std::vector<int>& set_out) const;

	/// reads amo and bmo from the restartdata file

	/// @return amo and bmo
	std::pair<MolecularOrbitals<T,NDIM>, MolecularOrbitals<T,NDIM> >
	static read_restartdata(World& world, const Molecule& molecule, const std::size_t nmo_alpha,
			const std::size_t nmo_beta) {
		bool spinrestricted = false;
		double current_energy;
		archive::ParallelInputArchive ar(world, "restartdata");
		ar & current_energy & spinrestricted;

		MolecularOrbitals<T,NDIM> amo, bmo;
		amo.load_mos(ar, molecule, nmo_alpha);
		bool have_beta=(not spinrestricted) and (nmo_beta>0);
		if (have_beta) {
			bmo.load_mos(ar,molecule,nmo_beta);
		}
		return std::make_pair(amo,bmo);
	}

	/// legacy code
	void load_mos(archive::ParallelInputArchive& ar, const Molecule& molecule, const std::size_t nmo_from_input) {

		unsigned int nmo = 0;

		ar & nmo;
		MADNESS_ASSERT(nmo >= nmo_from_input);
		ar & eps & occ & localize_sets;
		mo.resize(nmo);
		for (unsigned int i = 0; i < mo.size(); ++i)
			ar & mo[i];
		unsigned int n_core = molecule.n_core_orb_all();
		if (nmo > nmo_from_input) {
			localize_sets = vector<int>(localize_sets.begin() + n_core,
					localize_sets.begin() + n_core + nmo_from_input);
			mo = std::vector<Function<T,NDIM> >(mo.begin() + n_core,
					mo.begin() + n_core + nmo_from_input);
			eps = copy(eps(Slice(n_core, n_core + nmo_from_input - 1)));
			occ = copy(occ(Slice(n_core, n_core + nmo_from_input - 1)));
		}
	}

	void post_process_mos(World& world, const double thresh, const int k);


	/// save MOs in the AO projection for geometry restart

	/// compute the aos in MRA projection as:
	/// std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, calc.aobasis, calc.molecule);
	static void save_restartaodata(World& world, const std::vector<Function<double,3> > aos,
			const MolecularOrbitals<T,NDIM>& amo, const MolecularOrbitals<T,NDIM>& bmo) {

		Tensor<T> Saoamo = matrix_inner(world, aos, amo.get_mos());
		Tensor<T> Saobmo = (bmo.get_mos().size()>0) ? matrix_inner(world, aos, bmo.get_mos()) : Tensor<T>();
		if (world.rank() == 0) {
			archive::BinaryFstreamOutputArchive arao("restartaodata");
			arao << Saoamo << amo.get_eps() << amo.get_occ() << amo.get_localize_sets();
			if (Saobmo.size()>0) arao << Saobmo << bmo.get_eps() << bmo.get_occ() << bmo.get_localize_sets();
		}
	}


	/// uses AO-projection as a restart guess

	/// @return amo and bmo
	std::pair<MolecularOrbitals<T,NDIM>, MolecularOrbitals<T,NDIM> >
	static read_restartaodata(World& world, const std::vector<Function<double,3> > aos,
			const Molecule& molecule, const bool have_beta) {

		// project ao basis into MRA
		// compute the aos as
//		std::vector<Function<double,3> > aos=SCF::project_ao_basis_only(world, aobasis,molecule);

		archive::BinaryFstreamInputArchive arao("restartaodata");

		Tensor<T> Saomo;
		Tensor<double> eps, occ;
		std::vector<int> localize_set;
		std::vector<std::string> irrep;
		std::vector<Function<T,NDIM> > dummy_mo;

		MolecularOrbitals<T,NDIM> amo, bmo;

		try {
			// read alpha orbitals
			arao >> Saomo >> eps >> occ >> localize_set;
			amo=MolecularOrbitals<T,NDIM>(dummy_mo,eps,irrep,occ,localize_set);
			amo.project_ao(world,Saomo,aos);

			if (have_beta) {
				arao >> Saomo >> eps >> occ >> localize_set;
				bmo=MolecularOrbitals<T,NDIM>(dummy_mo,eps,irrep,occ,localize_set);
				bmo.project_ao(world,Saomo,aos);
			}

		} catch (...) {
			throw std::runtime_error("failed to read restartdata");
		}
		return std::make_pair(amo,bmo);
	}

	void project_ao(World& world, const Tensor<T>& Saomo, const std::vector<Function<double,3> >& aos);

private:
	std::vector<Function<T,NDIM> > mo;
	Tensor<double> eps;
	std::vector<std::string> irreps;
	Tensor<double> occ;
	std::vector<int> localize_sets;

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_MOLECULARORBITALS_H_ */
