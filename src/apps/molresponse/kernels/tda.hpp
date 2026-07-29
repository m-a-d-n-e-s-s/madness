#ifndef MOLRESPONSE_V3_KERNELS_TDA_HPP
#define MOLRESPONSE_V3_KERNELS_TDA_HPP

// =========================================================================
// Kernels<TDA, ClosedShell>  and  Kernels<TDA, OpenShell>
//
// Per-state building blocks for Tamm-Dancoff response. Naming follows
// thesis / molresponse convention:
//
//   Lambda = T·x + V0·x - E0_full·x + gamma     (has kinetic, used to
//                                                 build subspace A)
//   Theta  = V0·x - E0_nodi·x + gamma            (no kinetic, used as
//                                                 BSH input; in FD an
//                                                 extra V_perturbation
//                                                 source term is added)
//
// where E0_full = fock_full and E0_nodi = fock_no_diag. In a canonical
// orbital basis fock_full is diagonal (= diag(eps)) and fock_no_diag is
// zero. The choice is dictated by what BSH absorbs: BSH(omega) inverts
// (T - eps - omega), so theta carries the OFF-diagonal Fock only; the
// subspace matrix A_ij = <x_i|Lambda x_j> uses the full Fock so the
// eigenvalues of (A, S) are the omegas.
//
// Per-state operations (V0x, T0x, E0x_full, E0x_nodi, gamma) are
// individually accessible so:
//   - ES iteration:    build V0x/T0x/E0x_full/gamma, assemble Lambda,
//                      diagonalize, rotate {V0x, E0x_nodi, gamma} by U,
//                      assemble Theta from rotated pieces, BSH apply.
//   - FD iteration:    build V0x/E0x_nodi/gamma, assemble Theta + V_p,
//                      BSH apply.
// Both share the same V0x / E0x_nodi / gamma kernel bodies.
// =========================================================================

#include "../solvers/response_state.hpp"
#include "common_ops.hpp"   // poperatorT + common_ops::{apply_kinetic,
                            // bsh_shift, make_bsh_operators, apply_exchange}
#include "tags.hpp"
#include "two_electron.hpp"  // two_electron::{ExchangePair, apply_gamma}

#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/tensor/tensor.h>

namespace molresponse_v3 {


/// Zeroth-order data that response kernels read against. Field names
/// match `SCF::amo / bmo / aeps / beps / focka / fockb / Qa /
/// Qb` so a v3 reader who's familiar with `SCF.cc` finds the same
/// names. Closed-shell runs leave the `*_b` fields empty.
///
/// TODO(DFT): when c_xc != 1.0 (hybrid / pure DFT), gamma needs the
/// W = xc.apply_xc_kernel(rho) · phi term — see
/// molresponse_legacy/iterate_gamma.cc:114-123. Requires capturing an
/// XCOperator handle here.
///
/// TODO(load-balancing): per-response orbital load balance once we
/// fan out across MPI ranks. Exchange ops in compute_gamma / V0x are
/// O(N²) in n_occupied.
struct ResponseGroundState {
  std::vector<madness::real_function_3d> amo;      // alpha occupied orbitals
  madness::Tensor<double>                aeps;     // alpha orbital energies
  std::vector<madness::real_function_3d> bmo;      // beta occupied orbitals (empty if CS)
  madness::Tensor<double>                beps;

  madness::real_function_3d V_local_alpha;
  madness::real_function_3d V_local_beta;

  madness::Tensor<double>   focka;             // full alpha Fock (for Lambda)
  madness::Tensor<double>   fockb;
  madness::Tensor<double>   focka_no_diag;     // off-diag alpha Fock (for Theta)
  madness::Tensor<double>   fockb_no_diag;

  poperatorT                              coulop;
  madness::QProjector<double, 3>          Qa;    // alpha virtual-space projector
  madness::QProjector<double, 3>          Qb;
  double                                  c_xc = 1.0;
  double                                  lo   = 1.0e-10;

  int                                     n_roots = 0;

  // Cached ground-state HF exchange operators K0 = K[phi_occ, phi_occ], built
  // ONCE per protocol (build_response_ground_state_*). compute_V0x applies these
  // to the response orbitals every iteration; rebuilding the operator each call
  // deep-copies all occupied orbitals (Exchange::set_bra_and_ket → copy(world,
  // bra/ket)) — O(n_occ) churn per apply. K0_beta stays null for closed shell.
  // Built/applied at the protocol's k/thresh; rebuilt when the GS is re-prepared.
  std::shared_ptr<madness::Exchange<double, 3>> K0_alpha;
  std::shared_ptr<madness::Exchange<double, 3>> K0_beta;

  // Cached φ-only Coulomb-convolution tensor g0[i*n+k] = Poisson(φ_i·φ_k)
  // (n_occ² functions), built ONCE per protocol for the FD exchange-tensor path
  // (--fd-tensor; exch::build_g0) and read by assemble_theta_tensor. Empty unless
  // that gate is on — gate 0 / ES / VBC never populate it. Rebuilt per protocol
  // (depends on k/thresh). NOTE: n² RESIDENT functions — a memory-for-wall-time
  // trade (per-iter g0 rebuild → once/protocol); heavy at large n_occ (Inc-3 tiles).
  std::vector<madness::real_function_3d>        g0_alpha;
};


// =========================================================================
// Kernels<TDA, ClosedShell>
// =========================================================================
template <>
struct Kernels<TDA, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<ClosedShell>;

  // ---- per-state operator applications (the shared building blocks) ----

  /// Two-state perturbed density (general form), factor 2 (spin):
  ///   rho = 2 * Sum_p S1.x[p] * S2.x[p].  Linear response is S2 = Phi (ground).
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho1 = dot(world, S1.x_alpha, S2.x_alpha);
    rho1.scale(2.0);
    rho1.truncate();
    return rho1;
  }

  /// Convenience: density of `state` against the ground state {phi}.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    State Phi; Phi.x_alpha = g0.amo;
    return compute_density(world, g0, state, Phi);
  }

  /// Unified two-electron apply. TDA keeps ONE exchange term (drops Y):
  ///   out.x = Q( J*S3.x - c_xc*K(S2.x, S1.x)(S3.x) ).
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  /// gamma = Q( J[rho1]*amo - c_xc * K[amo, x_alpha](amo) ).
  ///
  /// Returns Storage so the (Type, Shell) uniform shape holds: ES Full
  /// will return gamma with both x_alpha and y_alpha components; OpenShell
  /// will also fill the *_beta components. TDA-ClosedShell carries
  /// x_alpha only.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    return out;
  }

  /// V0·x = V_local * x - c_xc * K[phi0, phi0](x).
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    if (g0.c_xc > 0.0) {
      auto k0x = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo,
                                                   state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
    }
    return State{std::move(Vx)};
  }

  /// T·x — kinetic energy applied to the response orbitals.
  static State
  compute_T0x(madness::World &world,
              const ResponseGroundState &/*target*/,
              const State    &state) {
    return State{common_ops::apply_kinetic(world, state.x_alpha)};
  }

  /// E0·x with the FULL Fock matrix (canonical: diag eps; localized:
  /// diag + off-diag). Used to assemble Lambda for the subspace matrix.
  static State
  compute_E0x_full(madness::World &world,
                   const ResponseGroundState &g0,
                   const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           g0.focka, vtol, true)};
  }

  /// E0·x with the OFF-DIAGONAL Fock only. Used to assemble Theta for
  /// the BSH driver, where the diagonal eps is absorbed into the BSH
  /// Green's function. In a canonical basis this is the zero vector.
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           g0.focka_no_diag, vtol, true)};
  }

  // Λ / θ assembly is shell-agnostic: see kernels/assembly.hpp for the
  // State-generic `assemble_lambda` / `assemble_theta` free functions.

  // ---- BSH and residual ----

  /// new_x = Q( BSH(omega) * (-2 * (theta + shift * x)) ).
  /// The shift matches the one inside make_bsh_operators(omega).
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    const double shift = common_ops::bsh_shift(g0.aeps, omega);

    auto rhs = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs, shift, state.x_alpha);
    scale(world, rhs, -2.0);
    truncate(world, rhs);

    auto bsh_ops = common_ops::make_bsh_operators(world, g0.aeps,
                                                  omega, g0.lo);
    auto new_x   = apply(world, bsh_ops, rhs);
    new_x        = g0.Qa(new_x);
    truncate(world, new_x);
    return State{std::move(new_x)};
  }

  /// || x_old - x_new ||.
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto diff = sub(world, state_old.x_alpha, state_new.x_alpha);
    return madness::norm2(world, diff);
  }
};

// =========================================================================
// Kernels<TDA, OpenShell>
//
// UHF TDA — same shape as Static OpenShell but with ONE exchange term
// per spin block (no Y coupling), plus T0x / E0x_full / Lambda used by
// the ES subspace eigenvalue problem.
//
// Density convention: rho = Σ x_α·φ_α + Σ x_β·φ_β  (NO extra factor —
// TDA drops Y, so we don't get the Y=X doubling that Static carries).
// At the restricted limit (α=β orbitals/X equal), this reduces to
// 2·Σ x·φ = Kernels<TDA, ClosedShell>'s density. ✓
//
// V_local is spin-blind (matches HF; spin-resolved Vxc DFT is a TODO).
// =========================================================================
template <>
struct Kernels<TDA, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<OpenShell>;

  /// Two-state perturbed density (general form), factor 1 (TDA drops Y):
  ///   rho = Sum_p ( S1.xa*S2.xa + S1.xb*S2.xb ).  Linear response is S2 = Phi.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho1 = dot(world, S1.x_alpha, S2.x_alpha);
    rho1 += dot(world, S1.x_beta, S2.x_beta);
    rho1.truncate();
    return rho1;
  }

  /// Convenience: density of `state` against the ground state {phi_a, phi_b}.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    State Phi; Phi.x_alpha = g0.amo; Phi.x_beta = g0.bmo;
    return compute_density(world, g0, state, Phi);
  }

  /// Unified two-electron apply (open shell TDA): one exchange term per spin.
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    out.x_beta  = two_electron::apply_gamma(world, J, S3.x_beta,
        {{S2.x_beta,  S1.x_beta }}, g0.Qb, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  /// gamma_α = Q_α( J[rho]*φ_α - c_xc K[φ_α, x_α](φ_α) )
  /// gamma_β = Q_β( J[rho]*φ_β - c_xc K[φ_β, x_β](φ_β) )
  /// (Single exchange term per spin — that's what distinguishes TDA
  /// from Static, which has K[φ,x](φ) + K[x,φ](φ) per spin.)
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    out.x_beta  = two_electron::apply_gamma(world, J_rho, g0.bmo,
        {{g0.bmo, state.x_beta }}, g0.Qb, g0.c_xc, g0.lo);
    return out;
  }

  /// V0·x = V_local · x − c_xc · K[φ, φ](x), per spin.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;

    auto Vx_a = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_ax = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx_a, -g0.c_xc, k0_ax);
    }
    auto Vx_b = mul_sparse(world, g0.V_local_beta, state.x_beta, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_bx = common_ops::apply_ground_exchange(world, g0.K0_beta, g0.bmo, state.x_beta, g0.lo);
      gaxpy(world, 1.0, Vx_b, -g0.c_xc, k0_bx);
    }
    return State{std::move(Vx_a), std::move(Vx_b)};
  }

  /// T·x per spin (kinetic energy on response orbitals — diagonal in spin).
  static State
  compute_T0x(madness::World &world,
              const ResponseGroundState &/*target*/,
              const State    &state) {
    auto Ta = common_ops::apply_kinetic(world, state.x_alpha);
    auto Tb = common_ops::apply_kinetic(world, state.x_beta);
    return State{std::move(Ta), std::move(Tb)};
  }

  /// E0·x with the FULL Fock per spin (Lambda assembly: T0 + V0 − E0_full + γ).
  static State
  compute_E0x_full(madness::World &world,
                   const ResponseGroundState &g0,
                   const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ea = transform(world, state.x_alpha,
                        g0.focka, vtol, true);
    auto Eb = transform(world, state.x_beta,
                        g0.fockb,  vtol, true);
    return State{std::move(Ea), std::move(Eb)};
  }

  /// E0·x with the OFF-DIAGONAL Fock per spin (Theta assembly).
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ea = transform(world, state.x_alpha,
                        g0.focka_no_diag, vtol, true);
    auto Eb = transform(world, state.x_beta,
                        g0.fockb_no_diag,  vtol, true);
    return State{std::move(Ea), std::move(Eb)};
  }

  // Λ / θ assembly is shell-agnostic: see kernels/assembly.hpp.

  /// Per-spin BSH apply, each side with its own ε and shift.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    // alpha
    const double shift_a = common_ops::bsh_shift(g0.aeps, omega);
    auto rhs_a = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_a, shift_a, state.x_alpha);
    scale(world, rhs_a, -2.0);
    truncate(world, rhs_a);
    auto bsh_a = common_ops::make_bsh_operators(world, g0.aeps,
                                                omega, g0.lo);
    auto new_xa = apply(world, bsh_a, rhs_a);
    new_xa      = g0.Qa(new_xa);
    truncate(world, new_xa);

    // beta
    const double shift_b = common_ops::bsh_shift(g0.beps, omega);
    auto rhs_b = madness::copy(world, theta.x_beta);
    gaxpy(world, 1.0, rhs_b, shift_b, state.x_beta);
    scale(world, rhs_b, -2.0);
    truncate(world, rhs_b);
    auto bsh_b = common_ops::make_bsh_operators(world, g0.beps,
                                                omega, g0.lo);
    auto new_xb = apply(world, bsh_b, rhs_b);
    new_xb      = g0.Qb(new_xb);
    truncate(world, new_xb);

    return State{std::move(new_xa), std::move(new_xb)};
  }

  /// Combined-norm residual: sqrt(||Δx_α||² + ||Δx_β||²).
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto da = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto db = sub(world, state_old.x_beta,  state_new.x_beta);
    const double na = madness::norm2(world, da);
    const double nb = madness::norm2(world, db);
    return std::sqrt(na * na + nb * nb);
  }
};

} // namespace molresponse_v3

// Interface enforcement — TDA is ES-only, so it must satisfy the
// extended ESKernel contract (FDKernel + Lambda machinery).
#include "kernel_interface.hpp"
MV3_ASSERT_ES_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::TDA,
                                                ::molresponse_v3::ClosedShell>);
MV3_ASSERT_ES_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::TDA,
                                                ::molresponse_v3::OpenShell>);

#endif // MOLRESPONSE_V3_KERNELS_TDA_HPP
