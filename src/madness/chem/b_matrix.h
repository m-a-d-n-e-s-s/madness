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
#include <madness/chem/SCFOperators.h>
#include <madness/chem/lowrankfunction.h>
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

    /// Full B matrix: B = B^Ue + B^g + B^KffK + B^X  (LRF commutator, may stall)
    static Tensor<double> compute(World& world, const Info& info);

    /// Full B matrix using exact 6D [K,f12] (no LRF stalling, controlled by thresh_6D)
    static Tensor<double> compute_6d(World& world, const Info& info);

    /// B^Ue_{ij,kl} = +<ij| f12 Q12 Ue |kl>
    static Tensor<double> compute_B_Ue(World& world, const Info& info);

    /// B^g_{ij,kl} = -<ij| f12 Q12 g12 |kl> = -V_{ij,kl}
    static Tensor<double> compute_B_g(World& world, const Info& info);

    /// B^KffK_{ij,kl} = -<ij| f12 Q12 [K12,f12] |kl>  (LRF approximation, may stall)
    static Tensor<double> compute_B_KffK(World& world, const Info& info);

    /// B^KffK_{ij,kl} via exact 6D apply_KffK_6d (no LRF, controlled by thresh_6D)
    static Tensor<double> compute_B_KffK_6d(World& world, const Info& info);

    /// B^X_{ij,kl} = (E_kl - E_ij) * X_{ij,kl}
    static Tensor<double> compute_B_X(World& world, const Info& info);

    /// X_{ij,kl} = <ij| f12 Q12 f12 |kl>  (pure 3D, no 6D materialization)
    static Tensor<double> compute_X_matrix(World& world, const Info& info);

    /// B matrix via LRF factorization of f12:
    ///   B_{ij,kl} = sum_{pq} [ Fhat_pq^{ik} That_pq^{jl} + Shat_pq^{ik} Ghat_pq^{jl} ]
    ///             - (eps_i+eps_j+eps_k+eps_l)/2 * X_{ij,kl}
    /// where Shat, That are projected overlaps and Fhat, Ghat are projected Fock elements
    /// of the LRF orbital products (1-O1)(a_p phi_i).  No 6D LRF of [K,f12] is needed.
    /// fock_ptr: Fock operator from nemo->make_fock_operator(); T will be removed internally.
    static Tensor<double> compute_via_lrf(World& world, const Info& info,
                                           const LowRankFunction<double,6>& lrf_f12,
                                           std::shared_ptr<Fock<double,3>> fock_ptr);

    /// Diagonal B element B_{ij,ij} via an orbital-weighted LRF of f12*phi_i*phi_j.
    ///
    /// lrf_f12_ij must approximate f12(r1,r2)*phi_i(r1)*phi_j(r2).  Its get_g()/get_h()
    /// factors already carry the orbital pre-factor, so no orbital multiplication is
    /// performed here.  The formula is
    ///   B_{ij,ij} = sum_{pq} [ Mhat_pq * Tflat_pq + Sflat_pq * Nhat_pq ]
    ///             - (eps_i + eps_j) * X_lrf
    /// where alpha~_p = (1-O1) g_p, beta~_p = (1-O2) h_p and
    ///   Sflat_pq = <R^2 alpha~_p | alpha~_q>,   Mhat_pq = <R^2 alpha~_p | F1 | alpha~_q>
    ///   Tflat_pq = <R^2 beta~_p  | beta~_q>,    Nhat_pq = <R^2 beta~_p  | F2 | beta~_q>
    ///   X_lrf    = sum_{pq} Sflat_pq * Tflat_pq
    /// i, j are absolute orbital indices (freeze already applied externally if needed).
    static double compute_via_lrf_pair(World& world, const Info& info,
                                        size_t i, size_t j,
                                        const LowRankFunction<double,6>& lrf_f12_ij,
                                        std::shared_ptr<Fock<double,3>> fock_ptr);

    /// Hylleraas pair correlation energy for pair (i,j) at convergence:
    ///   e^{ij} = <ij|g12|u_ij> + hc  +  <ij|g12 Q f12|ij>  -  <u_ij|c_ij>  -  <ij|f12 Q|c_ij>
    /// where  c_ij = (Ue - [K,f12])|ij>  and hc = <u_ij|g12|ij> (= forward term for real orbitals).
    /// i,j are absolute orbital indices (freeze offset NOT subtracted).
    /// u_ij is the converged 6D MP2 pair correction function.
    static double compute_hylleraas_pair(World& world, const Info& info,
                                          int i, int j,
                                          const real_function_6d& u_ij);
};

} // namespace madness
#endif
