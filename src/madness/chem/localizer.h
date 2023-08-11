//
// Created by Florian Bischoff on 11/1/21.
//

#ifndef MADNESS_LOCALIZER_H
#define MADNESS_LOCALIZER_H

#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/molecule.h>
#include <madness/tensor/distributed_matrix.h>
#include<madness/mra/mra.h>
#include<madness/chem/distpm.h>


namespace madness {

class SCF;

//template<typename T, std::size_t NDIM>
class Localizer {
public:

    Localizer() = default;

    Localizer(World& world, const AtomicBasisSet& aobasis, const Molecule& molecule,
              const std::vector<Function<double, 3>>& ao);

    Localizer& set_method(const std::string method1) {
        method=method1;
        return *this;
    }

    std::string get_method() {
        return method;
    }

    Localizer& set_metric(const Function<double,3>& R) {
        metric=copy(R);
        return *this;
    }

    Localizer& set_enforce_core_valence_separation(const bool value) {
        enforce_core_valence_separation=value;
        return *this;
    }

    AtomicBasisSet get_aobasis() const {
        return aobasis;
    }

    void print_info() const {
        print("Localizer info");
        print("method  ",method);
        print("aobasis ",aobasis.get_name());
        print("metric  ", metric.is_initialized());
        print("core-valence separation ",enforce_core_valence_separation);
        print("thresh_degenerate ",thresh_degenerate);
    }

    /// localize the orbitals
    template<typename T, std::size_t NDIM>
    MolecularOrbitals<T, NDIM> localize(const MolecularOrbitals<T, NDIM>& mo_in, bool randomize) const;

    /// localize the orbitals, possibly enforce core-valence separation
    template<typename T, std::size_t NDIM>
    MolecularOrbitals<T, NDIM> localize(const MolecularOrbitals<T, NDIM>& mo_in, const Tensor<T>& Fock,
                                        bool randomize) const;

    template<typename T, std::size_t NDIM>
    MolecularOrbitals<T, NDIM> separate_core_valence(const MolecularOrbitals<T, NDIM>& mo_in, const Tensor<T>& Fock) const;

    template<typename T, std::size_t NDIM>
    Tensor<T> compute_localization_matrix(World& world, const MolecularOrbitals<T, NDIM>& mo_in, bool randomize) const;

    /// localize orbitals while enforcing core-valence separation

    /// @param[in]  World   the world
    /// @param[in]  mo_in   the input orbitals
    /// @param[in]  Fock    the Fock matrix for canonicalizing the orbitals first
    /// @param[in]  method  the localization method
    /// @param[in]  tolloc  localization tolerance
    /// @param[in]  randomize   initially randomize the localization procedure
    template<typename T, std::size_t NDIM>
    Tensor<T> compute_core_valence_separation_transformation_matrix(World& world,
                                        const MolecularOrbitals<T, NDIM>& mo_in, const Tensor<T>& Fock) const;

    template<typename T>
    static bool check_core_valence_separation(const Tensor<T>& Fock, const std::vector<int>& localized_set,
                                              const bool silent=false);

    template<typename T>
    static std::size_t determine_frozen_orbitals(const Tensor<T> fmat);

    static bool check_frozen_consistency(const long nfrozen, const std::vector<int>& localize_sets);

        /// given a unitary transformation matrix undo mere reordering
    template<typename T>
    static void undo_reordering(Tensor<T>& U, const Tensor<double>& occ) {
        Tensor<double> eval(U.dim(0)); // dummy tensor
        undo_reordering(U,occ,eval);
    }

    /// given a unitary transformation matrix undo mere reordering
    template<typename T>
    static void undo_reordering(Tensor<T>& U, const Tensor<double>& occ, Tensor<double>& eval);

    /// given a unitary transformation matrix undo rotations between degenerate columns
    template<typename T>
    static void undo_degenerate_rotations(Tensor<T>& U, const Tensor<double>& eval, const double thresh_degenerate);

    /// given a unitary transformation matrix undo rotations within blocks of localized orbitals
    template<typename T>
    static void undo_rotations_within_sets(Tensor<T>& U, const std::vector<int>& localized_set);

    /// find sets of degenerate states/orbitals
    static std::vector<Slice> find_degenerate_blocks(const Tensor<double>& eval, const double thresh_degenerate);

    /// given a unitary transformation matrix undo the rotations within the blocks
    template<typename T>
    static Tensor<T> undo_rotation(const Tensor<T>& U_in, const std::vector<Slice>& blocks);

private:

    template<typename T, std::size_t NDIM>
    DistributedMatrix<T>
    localize_PM(World& world, const std::vector<Function<T, NDIM>>& mo, const std::vector<int>& set,
                const double thresh = 1e-9, const bool randomize = true, const bool doprint = false) const;

    template<typename T, std::size_t NDIM>
    DistributedMatrix<T> localize_boys(World& world,
                                       const std::vector<Function<T, NDIM>>& mo,
                                       const std::vector<int>& set,
                                       const double thresh = 1e-9,
                                       const bool randomize = true,
                                       const bool doprint = false) const;

    template<typename T, std::size_t NDIM>
    DistributedMatrix<T> localize_new(World& world,
                                      const std::vector<Function<T, NDIM>>& mo,
                                      const std::vector<int>& set,
                                      const double thresh = 1e-9,
                                      const bool randomize = true,
                                      const bool doprint = false) const;

    template<typename T>
    inline double DIP(const Tensor<T>& dip, int i, int j, int k, int l) const {
        return dip(i, j, 0) * dip(k, l, 0) + dip(i, j, 1) * dip(k, l, 1) + dip(i, j, 2) * dip(k, l, 2);
    }

    template<typename T>
    Tensor<T> matrix_exponential(const Tensor<T>& A) const;


    std::vector<int> at_to_bf, at_nbf;  /// map atoms to basis functions in the "new" algorithm
    AtomicBasisSet aobasis;             ///
    Molecule molecule;
    std::vector<Function<double, 3>> ao;
    Function<double,3> metric;       /// =R for computing matrix elements of operators
    double thetamax=0.1;                /// maximum rotation(?)
    const double tolloc = 1e-6; // was std::min(1e-6,0.01*dconv) but now trying to avoid unnecessary change
    double thresh_degenerate;           /// when are orbitals degenerate
    bool enforce_core_valence_separation=false;  /// no rotations between core and valence orbitals (distinguished by 'set')
    std::string method="new";           /// localization method

};

}

#endif //MADNESS_LOCALIZER_H
