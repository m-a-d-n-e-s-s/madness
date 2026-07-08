#ifndef MOLRESPONSE_V3_KERNELS_STATIC_HPP
#define MOLRESPONSE_V3_KERNELS_STATIC_HPP

// =========================================================================
// Kernels<Static, ClosedShell>
//
// Frequency-dependent response at omega = 0 (the "static" case). For
// closed shells, Y == X by construction, so Storage = ResponseStateX
// (X only). The per-root building blocks (compute_density, gamma, V0x,
// E0x, bsh_apply, compute_residual_norm) are byte-identical to
// Kernels<TDA, ClosedShell>'s — what's absent here is everything used
// to build the subspace eigenproblem (T0x, E0x_full). The FD solver
// doesn't diagonalize, so omitting those keeps the signature surface
// honest: a Static-FD step cannot accidentally call a routine it has
// no business calling. θ assembly itself is shell-agnostic and lives
// in kernels/assembly.hpp.
//
// The bsh_apply signature takes omega like TDA's so FDSolver can pass
// channel.omega uniformly (== 0.0 here).
// =========================================================================

#include "tags.hpp"
#include "tda.hpp"          // ResponseGroundState, common_ops::*
#include "two_electron.hpp"  // two_electron::{ExchangePair, apply_gamma}
#include "../solvers/response_state.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/world/worldprofile.h>  // PROFILE_BLOCK (perf-model meters; no-op unless WORLD_PROFILE_ENABLE)

namespace molresponse_v3 {

template <>
struct Kernels<Static, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<ClosedShell>;

  /// Singlet closed-shell static response density:
  ///   rho1 = 4 * sum_p amo[p] * x_alpha[p]
  ///        = (spin × 2) · (Y=X × 2) · Σ_p φ_p · x_p
  ///
  /// Factor 4 matches the Coulomb coefficient in the singlet (A+B)
  /// matrix element 4(pa|qb)(X+Y) at Y=X, and is the convention used
  /// by legacy v3 ResponseKernel.hpp::compute_response_density for the
  /// Static case. Self-consistent with Full's density formula at Y=X:
  ///   2·Σφ·(x+y) |_{y=x} = 4·Σφ·x.
  /// Two-state perturbed density (general form), factor 4 (spin x Y=X doubling):
  ///   rho = 4 * Sum_p S1.x[p] * S2.x[p].  Linear response is S2 = Phi (ground).
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho1 = dot(world, S1.x_alpha, S2.x_alpha);
    rho1.scale(4.0);
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

  /// Unified two-electron apply. Static keeps BOTH exchange terms (Y=X limit of
  /// Full): out.x = Q( J*S3.x - c_xc(K(S2.x,S1.x)+K(S1.x,S2.x))(S3.x) ).
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}, {S1.x_alpha, S2.x_alpha}},
        g0.Qa, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  /// gamma = Q( J[rho1]*phi0 - c_xc * (K[phi0, x](phi0) + K[x, phi0](phi0)) ).
  ///
  /// Static = the Y=X limit of Full, so BOTH exchange terms remain.
  /// (TDA drops the second one by neglecting Y; that's the only
  /// difference between Static and TDA gammas at the kernel level.)
  /// Mirrors molresponse_v3/ResponseKernel.hpp::compute_gamma Static case.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    // K[phi0,x](phi0) + K[x,phi0](phi0) -- the second term distinguishes Static from TDA.
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}, {state.x_alpha, g0.amo}},
        g0.Qa, g0.c_xc, g0.lo);
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
      auto k0x = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
    }
    return State{std::move(Vx)};
  }

  /// E0·x with the OFF-DIAGONAL Fock only (diag eps absorbed into BSH).
  /// In a canonical basis this is zero; kept for localized-orbital paths.
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    return State{transform(world, state.x_alpha,
                           g0.focka_no_diag, vtol, true)};
  }

  // θ assembly is shell-agnostic: see kernels/assembly.hpp.

  /// new_x = Q( BSH(omega) * (-2 * (theta + shift * x)) ).
  /// `omega` is passed for shape uniformity with TDA / Full kernels;
  /// FDSolver<Static,*> always passes 0.0.
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
    { PROFILE_BLOCK(rs_projection); new_x = g0.Qa(new_x); }  // FD Q·v projection meter
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
// Kernels<Static, OpenShell>
//
// UHF static response. Storage = ResponseStateX<OpenShell> ({x_alpha,
// x_beta}). The Coulomb response density is spin-summed (single
// real_function_3d carrying total δρ), while exchange and BSH are
// per-spin because exchange is same-spin and eps_α ≠ eps_β.
//
// Density convention (matches legacy v3 ResponseKernel.hpp::
// compute_response_density for unrestricted Static):
//
//     rho = 2 * (Σ x_α · φ0_α + Σ x_β · φ0_β)
//         = (spin_factor=1 for unrestricted)
//           × (y_factor=2 for Static)
//           × Σ_spin Σ_p x_p · φ_p
//
// V_local is shared between spins (it's the spin-blind ground-state
// local potential V_nuc + V_J + V_xc^total). For HF (c_xc=1, no XC
// kernel) this is exact; for hybrid DFT with spin-resolved Vxc the
// per-spin V_local fields would need wiring (open TODO).
// =========================================================================
template <>
struct Kernels<Static, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateX<OpenShell>;

  /// Two-state perturbed density (general form), factor 2 (Y=X doubling):
  ///   rho = 2 * Sum_p ( S1.xa*S2.xa + S1.xb*S2.xb ).  Linear: S2 = Phi.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho = dot(world, S1.x_alpha, S2.x_alpha);
    rho += dot(world, S1.x_beta, S2.x_beta);
    rho.scale(2.0);
    rho.truncate();
    return rho;
  }

  /// Convenience: density of `state` against the ground state {phi_a, phi_b}.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    State Phi; Phi.x_alpha = g0.amo; Phi.x_beta = g0.bmo;
    return compute_density(world, g0, state, Phi);
  }

  /// Unified two-electron apply (open shell static): both exchange terms per spin.
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}, {S1.x_alpha, S2.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    out.x_beta  = two_electron::apply_gamma(world, J, S3.x_beta,
        {{S2.x_beta,  S1.x_beta }, {S1.x_beta,  S2.x_beta }}, g0.Qb, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  /// gamma_α = Q_α( J[rho]*φ_α - c_xc(K[φ_α,x_α](φ_α) + K[x_α,φ_α](φ_α)) )
  /// gamma_β = Q_β( J[rho]*φ_β - c_xc(K[φ_β,x_β](φ_β) + K[x_β,φ_β](φ_β)) )
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}, {state.x_alpha, g0.amo}}, g0.Qa, g0.c_xc, g0.lo);
    out.x_beta  = two_electron::apply_gamma(world, J_rho, g0.bmo,
        {{g0.bmo, state.x_beta }, {state.x_beta,  g0.bmo}}, g0.Qb, g0.c_xc, g0.lo);
    return out;
  }

  /// V0·x for each spin. V_local is shared; K0 is same-spin.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    // Alpha
    auto Vx_a = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_ax = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.x_alpha, g0.lo);
      gaxpy(world, 1.0, Vx_a, -g0.c_xc, k0_ax);
    }
    // Beta
    auto Vx_b = mul_sparse(world, g0.V_local_beta, state.x_beta, vtol);
    if (g0.c_xc > 0.0) {
      auto k0_bx = common_ops::apply_ground_exchange(world, g0.K0_beta, g0.bmo, state.x_beta, g0.lo);
      gaxpy(world, 1.0, Vx_b, -g0.c_xc, k0_bx);
    }
    return State{std::move(Vx_a), std::move(Vx_b)};
  }

  /// E0·x using the per-spin off-diagonal Fock.
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

  // θ assembly is shell-agnostic: see kernels/assembly.hpp.

  /// Per-spin BSH apply. Each spin gets its own shift, BSH ops, and Q.
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
    { PROFILE_BLOCK(rs_projection); new_xa = g0.Qa(new_xa); }  // FD Q·v projection meter
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
    { PROFILE_BLOCK(rs_projection); new_xb = g0.Qb(new_xb); }  // FD Q·v projection meter
    truncate(world, new_xb);

    return State{std::move(new_xa), std::move(new_xb)};
  }

  /// Combined-norm residual: || (x_α, x_β)_old − (x_α, x_β)_new ||.
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

// Interface enforcement — both specializations must satisfy FDKernel.
// Static is FD-only; no Lambda machinery needed.
#include "kernel_interface.hpp"
MV3_ASSERT_FD_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Static,
                                                ::molresponse_v3::ClosedShell>);
MV3_ASSERT_FD_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Static,
                                                ::molresponse_v3::OpenShell>);

#endif // MOLRESPONSE_V3_KERNELS_STATIC_HPP
