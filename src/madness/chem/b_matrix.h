// b_matrix.h — B matrix in F12 theory
//
// B_{ij,kl} = <ij| f12 Q12 (Ue - g12 - [K12,f12]) |kl> + (E_kl - E_ij) * X_{ij,kl}
//
// where Ue = [T12,f12] + g12  (electronic regularization potential, regular part),
// and X_{ij,kl} = <ij| f12 Q12 f12 |kl> is the geminal metric matrix.
//
// Derivation (see b_matrix_plan.md for full details):
//   [F12, Q12] = 0  and  Q12^2 = Q12  collapse the two SO projectors.
//   Commuting (F12 - Eij) through f12 produces [F12,f12] + (Ekl-Eij) X.
//   [F12,f12] = [T12,f12] - [K12,f12] = (Ue - g12) - [K12,f12].

#ifndef MADNESS_CHEM_B_MATRIX_H
#define MADNESS_CHEM_B_MATRIX_H

#include <madness/chem/CCStructures.h>
#include <madness/chem/ccpairfunction.h>
#include <madness/chem/CCPotentials.h>
#include <madness/tensor/tensor.h>

namespace madness {

/// Computes the F12/R12 B matrix and its constituent parts.
///
/// All static methods take the Info struct that is already populated by
/// CCPotentials::update_info().  Indexing convention for the returned
/// tensors: the first two indices run over active bra pairs (i,j) and
/// the last two over active ket pairs (k,l), with the freeze offset
/// already subtracted:  B(i-freeze, j-freeze, k-freeze, l-freeze).
struct BMatrix {

    /// Full B matrix: B = B^Ue + B^g + B^KffK + B^X
    static Tensor<double> compute(World& world, const Info& info);

    /// B^Ue_{ij,kl} = +<ij| f12 Q12 Ue |kl>
    static Tensor<double> compute_B_Ue(World& world, const Info& info);

    /// B^g_{ij,kl} = -<ij| f12 Q12 g12 |kl> = -V_{ij,kl}
    static Tensor<double> compute_B_g(World& world, const Info& info);

    /// B^KffK_{ij,kl} = -<ij| f12 Q12 [K12,f12] |kl>
    static Tensor<double> compute_B_KffK(World& world, const Info& info);

    /// B^X_{ij,kl} = (E_kl - E_ij) * X_{ij,kl}
    static Tensor<double> compute_B_X(World& world, const Info& info);

    /// X_{ij,kl} = <ij| f12 Q12 f12 |kl>  (pure 3D, no 6D materialization)
    static Tensor<double> compute_X_matrix(World& world, const Info& info);
};

} // namespace madness
#endif
