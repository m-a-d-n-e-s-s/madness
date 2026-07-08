#ifndef MOLRESPONSE_V3_KERNELS_FULL_HPP
#define MOLRESPONSE_V3_KERNELS_FULL_HPP

// =========================================================================
// Kernels<Full, ClosedShell>
//
// Frequency-dependent response with X and Y amplitudes (omega != 0).
// Closed-shell, so Storage = ResponseStateXY<ClosedShell> ({x_alpha,
// y_alpha}). The response density couples X and Y:
//
//     rho1 = sum_p phi0[p] * (x_alpha[p] + y_alpha[p])
//
// (closed shell, so factor of 1 — the legacy spin-summed factor of 2
// is absorbed into the source term and into the BSH apply on -2·theta;
// see compute_density below for the exact convention used here, which
// matches molresponse_v3/ResponseKernel.hpp::compute_rho1).
//
// All per-Storage methods return Storage. The BSH apply is paired:
// theta.x_alpha is convolved with BSH(+omega) onto new x_alpha;
// theta.y_alpha is convolved with BSH(-omega) onto new y_alpha. Both
// share the same orbital-energy shift logic.
//
// Open-shell Full will reuse the same algorithmic shape with the
// alpha/beta blocks duplicated — defer until Phase 3 of the roadmap.
// =========================================================================

#include "assembly.hpp"     // assemble_lambda (used by apply_A±B helpers)
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
struct Kernels<Full, ClosedShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateXY<ClosedShell>;

  /// Singlet closed-shell Full response density:
  ///   rho1 = 2 * sum_p phi0[p] * (x_alpha[p] + y_alpha[p])
  ///        = spin × (X + Y contributions) × Σφ·(x+y)
  ///
  /// Factor 2 (spin sum) matches the Coulomb coefficient 4(pa|qb) on
  /// the X+Y combination in the singlet (A+B)·(X+Y) sub-block: in
  /// real space, J[2·Σφ(x+y)]·φ = 2·J[Σφ(x+y)]·φ pairs with the
  /// matrix-element 4(pa|qb)(X+Y).
  ///
  /// Self-consistency at Y=X: rho1|Y=X = 2·Σφ·(x+x) = 4·Σφ·x =
  /// `Kernels<Static, ClosedShell>::compute_density`.
  /// Two-state perturbed density (general form): the SAME operation the VBC
  /// quadratic source needs.  rho = 2 * Sum_p ( S1.x[p]*S2.x[p] + S1.y[p]*S2.y[p] ).
  /// Linear response is the case S2 = Phi (ground), giving 2*Sum phi*(x+y).
  /// Factor 2 (spin sum) matches the Coulomb coefficient 4(pa|qb) on X+Y in the
  /// singlet (A+B)(X+Y) sub-block.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho1 = dot(world, S1.x_alpha, S2.x_alpha);
    rho1 += dot(world, S1.y_alpha, S2.y_alpha);
    rho1.scale(2.0);
    rho1.truncate();
    return rho1;
  }

  /// Convenience: perturbed density of `state` against the ground state {phi,phi}.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    State Phi; Phi.x_alpha = g0.amo; Phi.y_alpha = g0.amo;
    return compute_density(world, g0, state, Phi);
  }

  /// Unified two-electron apply: out = Q( J[rho]*S3 - c_xc*K(S1,S2)(S3) ), with
  /// rho built from S1*S2 (returned for the linear iteration's Delta-rho check).
  /// compute_gamma(chi, rho) is the special case S1=chi, S2=S3=Phi.
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}, {S1.y_alpha, S2.y_alpha}},
        g0.Qa, g0.c_xc, g0.lo);
    out.y_alpha = two_electron::apply_gamma(world, J, S3.y_alpha,
        {{S2.y_alpha, S1.y_alpha}, {S1.x_alpha, S2.x_alpha}},
        g0.Qa, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  /// gamma_X = Q( J[rho1]*phi0 - c_xc * (K[phi0, x](phi0) + K[y, phi0](phi0)) )
  /// gamma_Y = Q( J[rho1]*phi0 - c_xc * (K[phi0, y](phi0) + K[x, phi0](phi0)) )
  ///
  /// Cross-channel exchange is what couples X and Y in the RPA equations;
  /// without it the result reduces to TDA on each side. This is apply_g with
  /// S1=state, S2=S3=Phi (ground) but the CALLER's rho1 (callers pass the
  /// folded ρ = 2·Σφ·(x+y)). The Coulomb piece is the SAME on X and Y.
  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    // X: K[phi0, x](phi0) + K[y, phi0](phi0)
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}, {state.y_alpha, g0.amo}},
        g0.Qa, g0.c_xc, g0.lo);
    // Y: K[phi0, y](phi0) + K[x, phi0](phi0)
    out.y_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.y_alpha}, {state.x_alpha, g0.amo}},
        g0.Qa, g0.c_xc, g0.lo);
    return out;
  }

  /// V0·x and V0·y. V_local and K0 each act on x and y independently.
  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vx = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    auto Vy = mul_sparse(world, g0.V_local_alpha, state.y_alpha, vtol);

    if (g0.c_xc > 0.0) {
      auto k0x = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.x_alpha, g0.lo);
      auto k0y = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.y_alpha, g0.lo);
      gaxpy(world, 1.0, Vx, -g0.c_xc, k0x);
      gaxpy(world, 1.0, Vy, -g0.c_xc, k0y);
    }
    return State{std::move(Vx), std::move(Vy)};
  }

  /// E0·x and E0·y with the off-diagonal Fock matrix (diag eps absorbed
  /// into BSH).
  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ex = transform(world, state.x_alpha,
                        g0.focka_no_diag, vtol, true);
    auto Ey = transform(world, state.y_alpha,
                        g0.focka_no_diag, vtol, true);
    return State{std::move(Ex), std::move(Ey)};
  }

  /// T·x and T·y — kinetic energy applied to both response blocks. The
  /// X and Y blocks are independent under T0; pairing happens only at
  /// the BSH step via the paired (+ω, −ω) Green's functions.
  static State
  compute_T0x(madness::World &world,
              const ResponseGroundState &/*g0*/,
              const State    &state) {
    return State{common_ops::apply_kinetic(world, state.x_alpha),
                 common_ops::apply_kinetic(world, state.y_alpha)};
  }

  /// E0·x and E0·y with the FULL Fock matrix (canonical: diag eps;
  /// localized: diag + off-diag). Used to assemble Λ for the subspace
  /// matrix in ESSolver<Full, ClosedShell>.
  static State
  compute_E0x_full(madness::World &world,
                   const ResponseGroundState &g0,
                   const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Ex = transform(world, state.x_alpha, g0.focka, vtol, true);
    auto Ey = transform(world, state.y_alpha, g0.focka, vtol, true);
    return State{std::move(Ex), std::move(Ey)};
  }

  // θ / Λ assembly is shell-agnostic: see kernels/assembly.hpp.

  /// Paired BSH apply.
  ///   new_x = Q( BSH(+omega) * (-2 * (theta_x + shift_+ * x_alpha)) )
  ///   new_y = Q( BSH(-omega) * (-2 * (theta_y + shift_- * y_alpha)) )
  /// Each side computes its own shift independently because the LUMO
  /// + omega vs LUMO - omega check produces different shifts.
  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    // --- X side: BSH at +omega ---------------------------------------
    const double shift_p = common_ops::bsh_shift(g0.aeps, omega);
    auto rhs_x = madness::copy(world, theta.x_alpha);
    gaxpy(world, 1.0, rhs_x, shift_p, state.x_alpha);
    scale(world, rhs_x, -2.0);
    truncate(world, rhs_x);
    auto bsh_p = common_ops::make_bsh_operators(world, g0.aeps,
                                                omega, g0.lo);
    auto new_x = apply(world, bsh_p, rhs_x);
    { PROFILE_BLOCK(rs_projection); new_x = g0.Qa(new_x); }  // FD Q·v projection meter
    truncate(world, new_x);

    // --- Y side: BSH at -omega ---------------------------------------
    const double shift_m = common_ops::bsh_shift(g0.aeps, -omega);
    auto rhs_y = madness::copy(world, theta.y_alpha);
    gaxpy(world, 1.0, rhs_y, shift_m, state.y_alpha);
    scale(world, rhs_y, -2.0);
    truncate(world, rhs_y);
    auto bsh_m = common_ops::make_bsh_operators(world, g0.aeps,
                                                -omega, g0.lo);
    auto new_y = apply(world, bsh_m, rhs_y);
    { PROFILE_BLOCK(rs_projection); new_y = g0.Qa(new_y); }  // FD Q·v projection meter
    truncate(world, new_y);

    return State{std::move(new_x), std::move(new_y)};
  }

  /// || (x_old, y_old) - (x_new, y_new) ||, combined norm.
  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto dx = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto dy = sub(world, state_old.y_alpha, state_new.y_alpha);
    const double nx = madness::norm2(world, dx);
    const double ny = madness::norm2(world, dy);
    return std::sqrt(nx * nx + ny * ny);
  }
};

// =========================================================================
// Kernels<Full, OpenShell>
//
// UHF dynamic response. Storage = ResponseStateXY<OpenShell>
// ({x_α, y_α, x_β, y_β}). Shared Coulomb across spins, per-spin
// exchange + BSH. Cross-channel exchange couples X with Y same as
// the closed-shell Full kernel.
//
// Density convention (legacy v3 unrestricted Full):
//     rho = 1 · 1 · (Σ (x_α + y_α)·φ_α + Σ (x_β + y_β)·φ_β)
//         = Σ_spin Σ_p (x_p + y_p) · φ_p
//
// (No spin or X+Y prefactor — for unrestricted, spin_factor=1; and
// for Full y_factor=1.)
// =========================================================================
template <>
struct Kernels<Full, OpenShell> {

  using vecfuncT = std::vector<madness::real_function_3d>;
  using State    = ResponseStateXY<OpenShell>;

  /// Two-state perturbed density (general form): factor 1, all four channels.
  ///   rho = Sum_p ( xa1*xa2 + ya1*ya2 + xb1*xb2 + yb1*yb2 ).
  /// Linear response is S2 = Phi (ground), giving Sum_spin phi*(x+y).
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &/*g0*/,
                  const State    &S1, const State &S2) {
    auto rho = dot(world, S1.x_alpha, S2.x_alpha);
    rho += dot(world, S1.y_alpha, S2.y_alpha);
    rho += dot(world, S1.x_beta,  S2.x_beta);
    rho += dot(world, S1.y_beta,  S2.y_beta);
    rho.truncate();
    return rho;
  }

  /// Convenience: density of `state` against the ground state {phi_a,phi_a,phi_b,phi_b}.
  static madness::real_function_3d
  compute_density(madness::World &world,
                  const ResponseGroundState &g0,
                  const State    &state) {
    State Phi;
    Phi.x_alpha = g0.amo; Phi.y_alpha = g0.amo;
    Phi.x_beta  = g0.bmo; Phi.y_beta  = g0.bmo;
    return compute_density(world, g0, state, Phi);
  }

  /// Unified two-electron apply (open shell): per-spin Coulomb + cross-channel
  /// exchange, projected by Qa/Qb. compute_gamma(state, rho) is S1=state, S2=S3=Phi.
  static std::pair<State, madness::real_function_3d>
  apply_g(madness::World &world,
          const ResponseGroundState &g0,
          const State &S1, const State &S2, const State &S3) {
    auto rho = compute_density(world, g0, S1, S2);
    auto J   = apply(*g0.coulop, rho);
    State out;
    out.x_alpha = two_electron::apply_gamma(world, J, S3.x_alpha,
        {{S2.x_alpha, S1.x_alpha}, {S1.y_alpha, S2.y_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    out.y_alpha = two_electron::apply_gamma(world, J, S3.y_alpha,
        {{S2.y_alpha, S1.y_alpha}, {S1.x_alpha, S2.x_alpha}}, g0.Qa, g0.c_xc, g0.lo);
    out.x_beta  = two_electron::apply_gamma(world, J, S3.x_beta,
        {{S2.x_beta,  S1.x_beta }, {S1.y_beta,  S2.y_beta }}, g0.Qb, g0.c_xc, g0.lo);
    out.y_beta  = two_electron::apply_gamma(world, J, S3.y_beta,
        {{S2.y_beta,  S1.y_beta }, {S1.x_beta,  S2.x_beta }}, g0.Qb, g0.c_xc, g0.lo);
    return {std::move(out), std::move(rho)};
  }

  static State
  compute_gamma(madness::World &world,
                const ResponseGroundState &g0,
                const State    &state,
                const madness::real_function_3d &rho1) {
    auto J_rho = apply(*g0.coulop, rho1);
    State out;
    // alpha: K[phi,x](phi) + K[y,phi](phi) ; Y swaps x<->y
    out.x_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.x_alpha}, {state.y_alpha, g0.amo}}, g0.Qa, g0.c_xc, g0.lo);
    out.y_alpha = two_electron::apply_gamma(world, J_rho, g0.amo,
        {{g0.amo, state.y_alpha}, {state.x_alpha, g0.amo}}, g0.Qa, g0.c_xc, g0.lo);
    // beta: same shape, beta orbitals, Qb
    out.x_beta  = two_electron::apply_gamma(world, J_rho, g0.bmo,
        {{g0.bmo, state.x_beta }, {state.y_beta,  g0.bmo}}, g0.Qb, g0.c_xc, g0.lo);
    out.y_beta  = two_electron::apply_gamma(world, J_rho, g0.bmo,
        {{g0.bmo, state.y_beta }, {state.x_beta,  g0.bmo}}, g0.Qb, g0.c_xc, g0.lo);
    return out;
  }

  static State
  compute_V0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Vxa = mul_sparse(world, g0.V_local_alpha, state.x_alpha, vtol);
    auto Vya = mul_sparse(world, g0.V_local_alpha, state.y_alpha, vtol);
    auto Vxb = mul_sparse(world, g0.V_local_beta,  state.x_beta,  vtol);
    auto Vyb = mul_sparse(world, g0.V_local_beta,  state.y_beta,  vtol);
    if (g0.c_xc > 0.0) {
      auto k0_ax = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.x_alpha, g0.lo);
      auto k0_ay = common_ops::apply_ground_exchange(world, g0.K0_alpha, g0.amo, state.y_alpha, g0.lo);
      gaxpy(world, 1.0, Vxa, -g0.c_xc, k0_ax);
      gaxpy(world, 1.0, Vya, -g0.c_xc, k0_ay);

      auto k0_bx = common_ops::apply_ground_exchange(world, g0.K0_beta, g0.bmo, state.x_beta, g0.lo);
      auto k0_by = common_ops::apply_ground_exchange(world, g0.K0_beta, g0.bmo, state.y_beta, g0.lo);
      gaxpy(world, 1.0, Vxb, -g0.c_xc, k0_bx);
      gaxpy(world, 1.0, Vyb, -g0.c_xc, k0_by);
    }
    return State{std::move(Vxa), std::move(Vya),
                 std::move(Vxb), std::move(Vyb)};
  }

  static State
  compute_E0x(madness::World &world,
              const ResponseGroundState &g0,
              const State    &state) {
    const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
    auto Exa = transform(world, state.x_alpha,
                         g0.focka_no_diag, vtol, true);
    auto Eya = transform(world, state.y_alpha,
                         g0.focka_no_diag, vtol, true);
    auto Exb = transform(world, state.x_beta,
                         g0.fockb_no_diag,  vtol, true);
    auto Eyb = transform(world, state.y_beta,
                         g0.fockb_no_diag,  vtol, true);
    return State{std::move(Exa), std::move(Eya),
                 std::move(Exb), std::move(Eyb)};
  }

  // θ assembly is shell-agnostic: see kernels/assembly.hpp.

  static State
  bsh_apply(madness::World &world,
            const ResponseGroundState &g0,
            const State    &state,
            const State    &theta,
            double          omega) {
    auto apply_side = [&](const madness::Tensor<double> &eps,
                          const madness::QProjector<double, 3> &Q,
                          const vecfuncT &theta_v,
                          const vecfuncT &x_v,
                          double signed_omega) {
      const double shift = common_ops::bsh_shift(eps, signed_omega);
      auto rhs = madness::copy(world, theta_v);
      gaxpy(world, 1.0, rhs, shift, x_v);
      scale(world, rhs, -2.0);
      truncate(world, rhs);
      auto bsh = common_ops::make_bsh_operators(world, eps, signed_omega,
                                                g0.lo);
      auto out = apply(world, bsh, rhs);
      { PROFILE_BLOCK(rs_projection); out = Q(out); }  // FD Q·v projection meter
      truncate(world, out);
      return out;
    };
    auto new_xa = apply_side(g0.aeps, g0.Qa,
                             theta.x_alpha, state.x_alpha,  omega);
    auto new_ya = apply_side(g0.aeps, g0.Qa,
                             theta.y_alpha, state.y_alpha, -omega);
    auto new_xb = apply_side(g0.beps,  g0.Qb,
                             theta.x_beta,  state.x_beta,   omega);
    auto new_yb = apply_side(g0.beps,  g0.Qb,
                             theta.y_beta,  state.y_beta,  -omega);
    return State{std::move(new_xa), std::move(new_ya),
                 std::move(new_xb), std::move(new_yb)};
  }

  static double
  compute_residual_norm(madness::World &world,
                        const State    &state_old,
                        const State    &state_new) {
    auto dxa = sub(world, state_old.x_alpha, state_new.x_alpha);
    auto dya = sub(world, state_old.y_alpha, state_new.y_alpha);
    auto dxb = sub(world, state_old.x_beta,  state_new.x_beta);
    auto dyb = sub(world, state_old.y_beta,  state_new.y_beta);
    const double na  = madness::norm2(world, dxa);
    const double nya = madness::norm2(world, dya);
    const double nb  = madness::norm2(world, dxb);
    const double nyb = madness::norm2(world, dyb);
    return std::sqrt(na*na + nya*nya + nb*nb + nyb*nyb);
  }
};

} // namespace molresponse_v3

// Interface enforcement.
//   ClosedShell: ESKernel (FD + ES — ESSolver<Full, ClosedShell> uses
//                compute_T0x / compute_E0x_full added above, plus the
//                shell-agnostic assemble_lambda / assemble_theta).
//   OpenShell:   FD-only (open-shell ES is out of scope for the current
//                molresponse_v3 phase).
#include "kernel_interface.hpp"
MV3_ASSERT_ES_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Full,
                                                ::molresponse_v3::ClosedShell>);
MV3_ASSERT_FD_KERNEL(::molresponse_v3::Kernels<::molresponse_v3::Full,
                                                ::molresponse_v3::OpenShell>);

#endif // MOLRESPONSE_V3_KERNELS_FULL_HPP
