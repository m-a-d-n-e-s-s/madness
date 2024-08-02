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
#include<madness/chem/molecule.h>
#include<stdio.h>

namespace madness {

template<typename T, std::size_t NDIM> class Function;
//class World;
class AtomicBasisSet;

template<typename T, std::size_t NDIM>
class MolecularOrbitals : public archive::ParallelSerializableObject {
public:

    MolecularOrbitals(const MolecularOrbitals<T,NDIM>& other) = default;

	MolecularOrbitals() = default;

    MolecularOrbitals(const std::vector<Function<T,NDIM> >& mo)
            : mo(mo), eps(), irreps(), occ(), localize_sets() {
    }

    MolecularOrbitals(const std::vector<Function<T,NDIM> >& mo, const Tensor<double>& eps)
            : mo(mo), eps(eps), irreps(), occ(), localize_sets() {
    }

	MolecularOrbitals(const std::vector<Function<T,NDIM> >& mo, const Tensor<double>& eps,
			const std::vector<std::string>& irrep, const Tensor<double>& occ,
			const std::vector<int>& set)
			: mo(mo), eps(eps), irreps(irrep), occ(occ), localize_sets(set) {
	}

    MolecularOrbitals get_subset(const int iset) const {
        auto slices = convert_set_to_slice(localize_sets);
        auto s=slices[iset];
        MolecularOrbitals result;
        MADNESS_CHECK(mo.size()>=size_t(s.end+1));
//        result.mo.assign(mo.begin()+s.start,mo.begin()+s.end+1);
        for (int i=s.start; i<s.end+1; ++i) result.mo.push_back(copy(mo[i]));
        if (eps.size()>0) result.eps=copy(eps(s));
        if (irreps.size()>0) result.irreps.assign(irreps.begin()+s.start,irreps.begin()+s.end+1);
        if (occ.size()>0) result.occ=copy(occ(s));
        if (localize_sets.size()>0) result.localize_sets.assign(localize_sets.begin()+s.start,localize_sets.begin()+s.end+1);

        return result;
    }

    std::size_t size() const {
        return mo.size();
    }

	std::vector<Function<T,NDIM> > get_mos() const {
		return mo;
	}

	[[nodiscard]] Tensor<double> get_eps() const {
		return eps;
	}

	[[nodiscard]] std::vector<std::string> get_irreps() const {
		return irreps;
	}

	std::vector<int> get_localize_sets() const {
		return localize_sets;
	}

	Tensor<double> get_occ() const {
		return occ;
	}

	/// setters will always invalidate all other member variables
	MolecularOrbitals& set_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		invalidate_all();
		mo=mo_new;
        return *this;
	}

	/// updates will keep other member variables
	MolecularOrbitals& update_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		mo=mo_new;
        return *this;
	}

	MolecularOrbitals& update_occ(const Tensor<double>& occ_new) {
		occ=occ_new;
        return *this;
	}

    /// updates will keep other member variables
    MolecularOrbitals& update_localize_set(const std::vector<int>& set) {
        localize_sets=set;
        return *this;
    }

	/// updates will keep other member variables
	MolecularOrbitals& update_mos_and_eps(const std::vector<Function<T,NDIM> >& mo_new,
			const Tensor<double>& eps_new) {
		mo=mo_new;
		eps=copy(eps_new);
        return *this;
	}

	MolecularOrbitals& recompute_irreps(const std::string pointgroup,
                       const Function<typename Tensor<T>::scalar_type,NDIM>& metric);

    /// group orbitals into sets of similar orbital energies for localization
	MolecularOrbitals& recompute_localize_sets(const double bandwidth=1.5) {
        set_all_orbitals_occupied();
        std::size_t nmo = mo.size();
        std::vector<int> set = std::vector<int>(static_cast<size_t>(nmo), 0);
        for (size_t i = 1; i < nmo; ++i) {
            set[i] = set[i - 1];
            // Only the new/boys localizers can tolerate not separating out the core orbitals
            if (eps(i) - eps(i - 1) > bandwidth || get_occ()(i) != 1.0) ++(set[i]);
        }
        update_localize_set(set);
        return *this;
	}

    static std::vector<Slice> convert_set_to_slice(const std::vector<int>& localized_set) {
        std::vector<Slice> blocks;
        long ilo=0;
        for (size_t i=1; i<localized_set.size(); ++i) {
            if (not (localized_set[i]==localized_set[i-1])) {
                blocks.push_back(Slice(ilo, i-1));
                ilo=i;
            }
        }
        // add final block
        blocks.push_back(Slice(ilo,localized_set.size()-1));
        return blocks;
    }


    MolecularOrbitals& set_all_orbitals_occupied() {
        occ=Tensor<double>(mo.size());
        occ=1.0;
        return *this;
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

	void pretty_print(std::string message, std::vector<std::string> flags=std::vector<std::string>()) const {
	    print(message);
        if (flags.size()==0) flags.resize(mo.size());
        std::vector<std::string> irreps=get_irreps();
        if (irreps.size()==0) irreps=std::vector<std::string>(mo.size(),"unknown");
	    print("orbital #   irrep   energy    occupation  localize_set");
        for (int i=mo.size()-1; i>=0; --i) {
//            double n=get_mos()[i].norm2();
            constexpr std::size_t bufsize=1024;
            char buf[bufsize];
            snprintf(buf,bufsize,"%5d %10s %12.8f  %6.2f  %8d  %15s", i, irreps[i].c_str(),get_eps()[i],
                   get_occ()[i],get_localize_sets()[i], flags[i].c_str());
            cout << std::string(buf) <<endl;
	    }
	}


    void print_frozen_orbitals(const long freeze) const {

        World& world=mo.front().world();
        MolecularOrbitals<T, 3> dummy_mo(*this);
        dummy_mo.recompute_localize_sets();
        if (world.rank() == 0) {
            auto flags=std::vector<std::string>(dummy_mo.get_mos().size(),"active");
            for (int i=0; i<freeze; ++i) flags[i]="frozen";
            dummy_mo.pretty_print("diagonal Fock matrix elements with core/valence separation for freezing",flags);
            print("\nfreezing orbitals: ", freeze,"\n");
        }

    }


    /// @param[in] cubefile_header  header of the cube file, from molecule::cubefile_header()
    void print_cubefiles(const std::string name, const std::vector<std::string> cubefile_header) const;

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
        /*
         * save and load from SCF.cc now contain the following extra variables
         * To read Molecular Orbitals properly we need to take care of this data
         * even if it is unused.
         *
         *
            unsigned int version;
            double L;
            int k;
            Molecule molecule;
            std::string xc;

*/
		bool spinrestricted = false;
        double L=0;
        int k1=0;                    // Ignored for restarting, used in response only
        unsigned int version = 4;  // UPDATE THIS IF YOU CHANGE ANYTHING
	//        unsigned int archive_version;
		double current_energy=0, converged_to_thresh=0;
        std::string xc, localize_method;


		archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, filename.c_str());
        ar & version;
		ar & current_energy & spinrestricted;
        ar & L& k1& molecule& xc & localize_method & converged_to_thresh;

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
        // TODO there is a good chance that this needs to be modified if it is intended to be read by SCF save/load
		bool spinrestricted = false;
        unsigned int version=4;
		double current_energy=0.0, converged_to_thresh=1.e10;
        double L=0;
        std::string xc, localize_method;
        int k1=0;                    // Ignored for restarting, used in response only
		archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, filename.c_str());
        ar & version;
        ar & current_energy & spinrestricted;
        ar & L& k1& molecule& xc & localize_method & converged_to_thresh;

		amo.save_mos(ar,molecule);
		bmo.save_mos(ar,molecule);
	}

	/// legacy code
    void load_mos(archive::ParallelInputArchive<>& ar, const Molecule& molecule, const std::size_t nmo_from_input) {

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
        void save_mos(archive::ParallelOutputArchive<>& ar, const Molecule& molecule) const {

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

    std::vector<Vector<typename Tensor<T>::scalar_type,3>> compute_center(
            const Function<typename Tensor<T>::scalar_type,NDIM> metric2=Function<typename Tensor<T>::scalar_type,NDIM>()) const;

private:
	std::vector<Function<T,NDIM> > mo;
	Tensor<double> eps;
	std::vector<std::string> irreps;
	Tensor<double> occ;
	std::vector<int> localize_sets;

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_MOLECULARORBITALS_H_ */
