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
//class World;
class AtomicBasisSet;

template<typename T, std::size_t NDIM>
class MolecularOrbitals : public archive::ParallelSerializableObject {
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

	void recompute_irreps(const std::string pointgroup,
                       const Function<typename Tensor<T>::scalar_type,NDIM>& metric);

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

	void pretty_print(std::string message) {
	    print(message);
	    print("orbital #   irrep   energy    occupation  localize_set");
        for (int i=mo.size()-1; i>=0; --i) {
//            double n=get_mos()[i].norm2();
            printf("%5d %10s %12.8f  %6.2f  %8d\n", i, get_irreps()[i].c_str(),get_eps()[i],
                   get_occ()[i],get_localize_sets()[i]);
	    }
	}

	template <typename Archive>
	void serialize (Archive& ar) {
		std::size_t nmo=mo.size();
		ar & nmo;
		if (nmo!=mo.size()) mo.resize(nmo);
		for (auto& m : mo) ar & m;
		ar & eps & irreps & occ & localize_sets;
		if (ar.is_input_archive) {
		    if (irreps.size()==0) irreps=std::vector<std::string>(nmo,"unknown");
            if (localize_sets.size()==0) localize_sets=std::vector<int>(nmo,0);
            if (occ.size()==0) occ=Tensor<double>(nmo);
            if (eps.size()==0) eps=Tensor<double>(nmo);
		}
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
	static read_restartdata(World& world, const std::string filename, const Molecule& molecule,
			const std::size_t nmo_alpha, const std::size_t nmo_beta) {
		bool spinrestricted = false;
		double current_energy;
		archive::ParallelInputArchive ar(world, filename.c_str());
		ar & current_energy & spinrestricted;

		MolecularOrbitals<T,NDIM> amo, bmo;
		amo.load_mos(ar, molecule, nmo_alpha);
		bool have_beta=(not spinrestricted) and (nmo_beta>0);
		if (have_beta) {
			bmo.load_mos(ar,molecule,nmo_beta);
		}
		return std::make_pair(amo,bmo);
	}

	/// reads amo and bmo from the restartdata file

	/// @return amo and bmo
	void static save_restartdata(World& world, const std::string filename, const Molecule& molecule,
			const MolecularOrbitals<T,NDIM>& amo, const MolecularOrbitals<T,NDIM>& bmo) {
		bool spinrestricted = false;
		double current_energy=0.0;
		archive::ParallelOutputArchive ar(world, filename.c_str());
		ar & current_energy & spinrestricted;

		amo.save_mos(ar,molecule);
		bmo.save_mos(ar,molecule);
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

	/// legacy code
	void save_mos(archive::ParallelOutputArchive& ar, const Molecule& molecule) const {

		unsigned int nmo=mo.size();

		ar & nmo;
		ar & eps & occ & localize_sets;
		for (unsigned int i = 0; i < mo.size(); ++i)
			ar & mo[i];
		//unsigned int n_core = molecule.n_core_orb_all();
	}

	void post_process_mos(World& world, const double thresh, const int k);


	/// save MOs in the AO projection for geometry restart

	/// compute the aos in MRA projection as:
	static void save_restartaodata(World& world, const Molecule& molecule,
			const MolecularOrbitals<T,NDIM>& amo, const MolecularOrbitals<T,NDIM>& bmo,
			const AtomicBasisSet& aobasis);

	/// uses AO-projection as a restart guess

	/// @return amo and bmo
	std::pair<MolecularOrbitals<T,NDIM>, MolecularOrbitals<T,NDIM> >
	static read_restartaodata(World& world,
			const Molecule& molecule, const bool have_beta);

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
