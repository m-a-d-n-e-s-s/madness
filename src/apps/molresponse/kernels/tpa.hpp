#ifndef MOLRESPONSE_V3_KERNELS_TPA_HPP
#define MOLRESPONSE_V3_KERNELS_TPA_HPP

// ===========================================================================
// tpa — two-photon-absorption transition moments (single residue of the
// quadratic response), closed-shell / TDA. See TPA_SCOPING.md (§2).
//
// 2PA is "beta with the C-channel replaced by the excited-state eigenvector
// X_f", evaluated at the two-photon resonance omega_1 = omega_2 = omega_f/2.
// The residue of the C-channel linear response at omega -> omega_f IS X_f, a
// HOMOGENEOUS eigenvector: it has no driving operator, so in the C-channel the
// perturbation operator VC_op -> 0. Everything else reuses the Dalton-validated
// beta contraction core verbatim:
//
//   S_ab = beta_abc( A = mu_a-response@(omega_f/2),
//                    VBC = compute_vbc( B = mu_b-response@(omega_f/2), C = X_f,
//                                       VB_op = mu_b, VC_op = 0 ),
//                    B, C = X_f, VA_op = mu_a )
//
// This is a CANDIDATE residue (the "VC_op=0 homogeneous C-channel" reduction);
// the exact surviving-term set + the X_f normalization are pinned numerically
// against the Dalton .QUADRA reference (refs/dalton_tpa.json). Because it reuses
// beta_abc/compute_vbc unchanged, a mismatch is a PHYSICS/normalization signal,
// not a contraction-code bug (alpha/beta/ES all already validate).
//
// The rotationally-invariant strength delta^|| (parallel linear polarization,
// F=G=H=2) is the quantity Dalton reports (two_photon_strengths); S_ab is the
// raw tensor.
// ===========================================================================

#include "beta.hpp"   // beta::beta_abc
#include "vbc.hpp"    // vbc::compute_vbc, ResponseGroundState, ResponseStateXY

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <array>
#include <cmath>

namespace molresponse_v3::tpa {

/// The 3x3 two-photon transition-moment tensor S_ab for one excited state f
/// (closed-shell / TDA). Inputs are already-converged states:
///   Xf        — the ES eigenvector for root f, wrapped as XY (TDA: y_alpha=0).
///   mu_resp   — the dipole response mu_a at omega_f/2, per Cartesian axis a.
///   mu_op     — the raw dipole operators (MomentFunctor), per axis.
/// No electronic solve: pure contraction via beta_abc + compute_vbc.
/// The three residue VBC quadratic sources for one excited state f (C = X_f
/// homogeneous, VC_op = 0), one per B-axis. This is the TPA "second-order
/// perturbation vector" — the analogue of beta's VBC source. assemble_tpa saves
/// these to disk (mirroring beta/Raman vbc_states) for reuse + inspection.
inline std::array<ResponseStateXY<ClosedShell>, 3>
tpa_sources(madness::World &world, const ResponseGroundState &g0,
            const ResponseStateXY<ClosedShell> &Xf,
            const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
            const std::array<madness::real_function_3d, 3> &mu_op) {
  madness::real_function_3d zero_op = madness::copy(mu_op[0]);  // VC_op = 0
  zero_op.scale(0.0);
  std::array<ResponseStateXY<ClosedShell>, 3> vbc_b;
  for (int b = 0; b < 3; ++b)
    vbc_b[b] = vbc::compute_vbc<ClosedShell>(world, g0, mu_resp[b], Xf,
                                             mu_op[b], zero_op);
  return vbc_b;
}

/// Contract PRE-BUILT sources into the 3x3 S tensor (symmetrized). Kept separate
/// from tpa_sources so assemble_tpa can build the sources once, save them, then
/// contract — the same build->save->contract shape as the beta path.
inline madness::Tensor<double>
tpa_moment(madness::World &world, const ResponseGroundState &g0,
           const ResponseStateXY<ClosedShell> &Xf,
           const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
           const std::array<madness::real_function_3d, 3> &mu_op,
           const std::array<ResponseStateXY<ClosedShell>, 3> &vbc_b) {
  using namespace madness;
  Tensor<double> S(3L, 3L);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      S(a, b) = beta::beta_abc<ClosedShell>(world, g0, mu_resp[a], vbc_b[b],
                                            mu_resp[b], Xf, mu_op[a]);

  // Degenerate TPA (two identical photons) => S_ab must be symmetric; the raw
  // beta_abc role split (A vs the B-in-VBC channel) computes ONE ordering, so
  // symmetrize to the physical S_ab = <0|mu_a|p><p|mu_b|f> + (a<->b). Validated
  // against Dalton, whose S is symmetric (upper triangle only).
  Tensor<double> Ssym(3L, 3L);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      Ssym(a, b) = 0.5 * (S(a, b) + S(b, a));
  return Ssym;
}

/// Convenience: build sources + contract in one call (no save).
inline madness::Tensor<double>
tpa_moment(madness::World &world, const ResponseGroundState &g0,
           const ResponseStateXY<ClosedShell> &Xf,
           const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
           const std::array<madness::real_function_3d, 3> &mu_op) {
  return tpa_moment(world, g0, Xf, mu_resp, mu_op,
                    tpa_sources(world, g0, Xf, mu_resp, mu_op));
}

/// SINGLE-RESIDUE 2PA moment (Parker JCTC 2018 / DALTON QRSMO), the corrected
/// channel assignment (TPA_SCOPING §5l-m). The residue of beta is taken in the
/// A channel (omega_sigma = omega_B+omega_C -> omega_f), so:
///   * V^{bc} is built from the TWO PHOTON RESPONSES N^b, N^c (both at
///     omega_f/2) WITH their dipole operators — the exact quadratic RHS beta
///     uses (Eq. 19 / compute_vbc);
///   * the A-channel response is REPLACED by the eigenvector X_f, so the
///     mu_A property terms (beta's b2/b3) drop — the contraction is the bare
///     paired metric inner (beta's b1 only):
///       S_bc = prefactor * ( <x_f|V^{bc}.x> + <y_f|V^{bc}.y> )
///   * normalization: the formula assumes the Parker/TDHF convention
///     X^t X - Y^t Y = 1 over spin-adapted vectors == 0.5 over our spatial
///     functions. Our ES solver normalizes the SPATIAL metric to 1, so its
///     eigenvectors carry sqrt(2) vs that convention -> fold C_N into
///     `prefactor` (pinned empirically via tools/tpa_from_dalton).
/// S is symmetric by construction (compute_vbc symmetrizes the (B,C)+(C,B)
/// halves), so only 6 V^{bc} builds per root.
inline madness::Tensor<double>
tpa_moment_residue(madness::World &world, const ResponseGroundState &g0,
                   const ResponseStateXY<ClosedShell> &Xf,
                   const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
                   const std::array<madness::real_function_3d, 3> &mu_op,
                   double prefactor = 1.0) {
  using namespace madness;
  Tensor<double> S(3L, 3L);
  for (int b = 0; b < 3; ++b) {
    for (int c = b; c < 3; ++c) {
      const double t0 = wall_time();
      if (world.rank() == 0) {
        printf("  [TPA residue] pair %c%c (%d/6): building V^{bc} "
               "(2 orderings x 3 two-electron builds)...\n",
               "xyz"[b], "xyz"[c], b * (5 - b) / 2 + c + 1);
        fflush(stdout);
      }
      auto V = vbc::compute_vbc<ClosedShell>(
          world, g0, mu_resp[static_cast<size_t>(b)],
          mu_resp[static_cast<size_t>(c)], mu_op[static_cast<size_t>(b)],
          mu_op[static_cast<size_t>(c)]);
      const double mx = inner(Xf.x_alpha, V.x_alpha);
      const double my = inner(Xf.y_alpha, V.y_alpha);
      // Diagnostic: both metric pairings. C_N derivation (TPA_SCOPING §5m)
      // forces + if the ES solver's y-sign convention matches the FD one that
      // validated beta's b1; the M− column decides that empirically (the
      // ~6% √2-vs-measured residual is exactly y-pairing sized).
      if (world.rank() == 0)
        printf("  [TPA residue] pair %c%c (%d/6): M+ = %+.6f   M- = %+.6f   "
               "(x-part %+.6f, y-part %+.6f)   [%.1f s]\n",
               "xyz"[b], "xyz"[c], b * (5 - b) / 2 + c + 1, mx + my, mx - my,
               mx, my, wall_time() - t0);
               fflush(stdout);
      S(b, c) = prefactor * (mx + my);
      S(c, b) = S(b, c);
    }
  }
  return S;
}

/// ONE-ELECTRON part of the single-residue contraction: the mu-operator terms
/// of V^{bc} only (compute_vbc_i's fb = -Q(v*c) and fphi = c*<phi|v|phi>,
/// both orderings) — no two-electron builds, so this costs inner products.
/// Validated against DALTON's A2B/A2C/X2 families (<=1.3% on every d-aug-QZ
/// ledger element). Production 2PA is S = C_N*(S_1e + S_E3corr), C_N=sqrt(2)
/// (verified 6/6 vs DALTON 2026-07-22, verify_e3_k6.py).
inline madness::Tensor<double>
tpa_moment_residue_1e(madness::World &world, const ResponseGroundState &g0,
                      const ResponseStateXY<ClosedShell> &Xf,
                      const std::array<ResponseStateXY<ClosedShell>, 3> &mu_resp,
                      const std::array<madness::real_function_3d, 3> &mu_op,
                      double prefactor = 1.0) {
  using namespace madness;
  const vecfuncT &phi = g0.amo;
  std::array<Tensor<double>, 3> mv;   // <phi_i | mu_a | phi_j>
  for (int a = 0; a < 3; ++a)
    mv[a] = matrix_inner(world, phi, mul(world, mu_op[a], phi, true));
  // one ordered half: v = mu_b acting on the C-state (mirrors compute_vbc_i)
  auto half = [&](int b, int c) {
    const auto &cx = mu_resp[static_cast<size_t>(c)].x_alpha;
    const auto &cy = mu_resp[static_cast<size_t>(c)].y_alpha;
    const auto &v = mu_op[static_cast<size_t>(b)];
    vecfuncT vx = g0.Qa(mul(world, v, cx, true)); scale(world, vx, -1.0);
    vecfuncT vy = g0.Qa(mul(world, v, cy, true)); scale(world, vy, -1.0);
    gaxpy(world, 1.0, vx, 1.0, transform(world, cx, mv[b], true));
    gaxpy(world, 1.0, vy, 1.0, transform(world, cy, mv[b], true));
    return inner(Xf.x_alpha, vx) + inner(Xf.y_alpha, vy);
  };
  Tensor<double> S(3L, 3L);
  for (int b = 0; b < 3; ++b)
    for (int c = b; c < 3; ++c)
      S(b, c) = S(c, b) = prefactor * (half(b, c) + half(c, b));
  return S;
}

/// CORRECTED two-electron E[3] part of the TPA single residue (TPA_SCOPING
/// §5q — block-level derivation from DALTON Q3FOCK/C3FCKO with the
/// (x,y)<->(Z,-Y) dictionary and VECB = -swap(N+) [ANTSYM=-1, DIPLEN]).
/// Returns the DALTON-units per-element scalar (== DALTON 'E3 CONTRIBUTION TO
/// SMOM' for this pair; the sqrt(2) is INCLUDED). Add to the UNCHANGED
/// one-electron content. Exactly symmetric under B<->C (analytic proof at the
/// integral level; assert numerically). Three term families:
///   T1: G[D12], D12 = |yb><yf| + |xf><xb| - |zetaD><phi|,
///       zetaD_i = sum_j phi_j (<x^b_j|x^f_i> + <y^f_j|y^b_i>)   [b-f SAME-block]
///   T2: gamma_B (standard channel) against {x^f, y^c, occ-transfers}
///   T3: gamma_F (TRANSPOSED channel — the eigenvector's own kernel, missing
///       from compute_vbc) against {y^b, y^c, occ-transfers}
inline double
tpa_e3_residue(madness::World &world, const ResponseGroundState &g0,
               const ResponseStateXY<ClosedShell> &B,
               const ResponseStateXY<ClosedShell> &C,
               const ResponseStateXY<ClosedShell> &F) {
  using namespace madness;
  const vecfuncT &phi = g0.amo;
  const vecfuncT &xb = B.x_alpha, &yb = B.y_alpha;
  const vecfuncT &xc = C.x_alpha, &yc = C.y_alpha;
  const vecfuncT &xf = F.x_alpha, &yf = F.y_alpha;
  const double cxc = g0.c_xc, lo = g0.lo;

  // occ-space transfer functions (transform(w, v, M)_i = sum_j v_j M(j,i))
  auto zetaD = transform(world, phi,
                         matrix_inner(world, xb, xf) +
                             matrix_inner(world, yf, yb), true);
  vecfuncT nzetaD = copy(world, zetaD);
  scale(world, nzetaD, -1.0);
  auto z_xfxc = transform(world, phi, matrix_inner(world, xf, xc), true);
  auto z_yfyc = transform(world, phi, matrix_inner(world, yf, yc), true);
  auto z_ybxc = transform(world, phi, matrix_inner(world, yb, xc), true);
  auto z_xbyc = transform(world, phi, matrix_inner(world, xb, yc), true);

  // T1: G[D12] = 2J[rho12] - K[D12], D12 = |yb><yf| + |xf><xb| - |zetaD><phi|
  const double t1s = wall_time();
  real_function_3d rho12 = common_ops::dot(world, yb, yf);
  rho12 += common_ops::dot(world, xf, xb);
  rho12 += common_ops::dot(world, phi, nzetaD);
  rho12.scale(2.0);
  rho12.truncate();
  auto J12 = apply(*g0.coulop, rho12);
  auto t1x = two_electron::apply_gamma_raw(
      world, J12, phi, {{yf, yb}, {xb, xf}, {phi, nzetaD}}, cxc, lo);
  auto t1y = two_electron::apply_gamma_raw(
      world, J12, phi, {{yb, yf}, {xf, xb}, {nzetaD, phi}}, cxc, lo);
  const double T1 = inner(xc, t1x) + inner(yc, t1y);

  // T2: gamma_B, standard channel (same pairs as vbc's gbc call)
  const double t2s = wall_time();
  real_function_3d rhoB = common_ops::dot(world, phi, xb);
  rhoB += common_ops::dot(world, phi, yb);
  rhoB.scale(2.0);
  rhoB.truncate();
  auto JB = apply(*g0.coulop, rhoB);
  auto gB = [&](const vecfuncT &t) {
    return two_electron::apply_gamma_raw(world, JB, t, {{xb, phi}, {phi, yb}},
                                         cxc, lo);
  };
  const double T2 = -inner(z_xfxc, gB(phi))    // 2a
                    + inner(xc, gB(xf))        // 2b
                    + inner(yf, gB(yc))        // 2c
                    - inner(phi, gB(z_yfyc));  // 2d

  // T3: gamma_F, TRANSPOSED channel (the eigenvector's kernel)
  real_function_3d rhoF = common_ops::dot(world, phi, xf);
  rhoF += common_ops::dot(world, phi, yf);
  rhoF.scale(2.0);
  rhoF.truncate();
  auto JF = apply(*g0.coulop, rhoF);
  auto gFt = [&](const vecfuncT &t) {
    return two_electron::apply_gamma_raw(world, JF, t, {{phi, xf}, {yf, phi}},
                                         cxc, lo);
  };
  const double t3s = wall_time();
  const double T3 = +inner(xc, gFt(yb))        // 3a
                    - inner(z_ybxc, gFt(phi))  // 3b
                    + inner(xb, gFt(yc))       // 3c
                    - inner(phi, gFt(z_xbyc)); // 3d

  if (world.rank() == 0)
    printf("      [e3corr timing] T1(G[D12]) %.1fs  T2(gammaB) %.1fs  "
           "T3(gammaF) %.1fs   [T1=%+.5f T2=%+.5f T3=%+.5f]\n",
           t2s - t1s, t3s - t2s, wall_time() - t3s, T1, T2, T3);
           fflush(stdout);
  // Global minus: DALTON contracts with VECB = -swap(N^B) for antisymmetric
  // (DIPLEN) operators (ANTSYM=-1, rspprp.F:543); our N^B carries the opposite
  // overall sign. Verified k6 2026-07-22: without it all six ledger elements
  // give E3 ratio -1.00 (spread 0.966..1.003), with it +1.00; 1e unaffected.
  return -std::sqrt(2.0) * (T1 + T2 + T3);
}

/// Two-array variant for convention probes: V^{bc} built from respB[b] and
/// respC[c] (e.g. respC = the (x<->y)-swapped responses, testing the -omega /
/// adjoint form of one E3 channel that DALTON's T3DRV uses — TPA_SCOPING §5n:
/// the one-electron terms match Dalton but E3(N+,N+) does not; the candidate
/// correct object is E3(N+,N-)). Result symmetrized over (b,c).
inline madness::Tensor<double>
tpa_moment_residue_bc(madness::World &world, const ResponseGroundState &g0,
                      const ResponseStateXY<ClosedShell> &Xf,
                      const std::array<ResponseStateXY<ClosedShell>, 3> &respB,
                      const std::array<ResponseStateXY<ClosedShell>, 3> &respC,
                      const std::array<madness::real_function_3d, 3> &opB,
                      const std::array<madness::real_function_3d, 3> &opC,
                      double prefactor = 1.0) {
  using namespace madness;
  Tensor<double> S(3L, 3L);
  for (int b = 0; b < 3; ++b) {
    for (int c = 0; c < 3; ++c) {
      auto V = vbc::compute_vbc<ClosedShell>(
          world, g0, respB[static_cast<size_t>(b)],
          respC[static_cast<size_t>(c)], opB[static_cast<size_t>(b)],
          opC[static_cast<size_t>(c)]);
      const double m =
          inner(Xf.x_alpha, V.x_alpha) + inner(Xf.y_alpha, V.y_alpha);
      S(b, c) = prefactor * m;
    }
  }
  // symmetrize (degenerate photons)
  Tensor<double> Ssym(3L, 3L);
  for (int b = 0; b < 3; ++b)
    for (int c = 0; c < 3; ++c)
      Ssym(b, c) = 0.5 * (S(b, c) + S(c, b));
  return Ssym;
}

/// Rotationally-invariant two-photon strength for parallel linear polarization:
///   delta^|| = Sum_ab ( F S_aa S_bb + G S_ab S_ab + H S_ab S_ba ), F=G=H=2.
/// (Real S here — TDA, undamped, off-resonance.) This is what Dalton's
/// two_photon_strengths reports, so it is the validation target.
inline double
delta_parallel(const madness::Tensor<double> &S) {
  double d = 0.0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      d += 2.0 * S(a, a) * S(b, b) + 2.0 * S(a, b) * S(a, b) +
           2.0 * S(a, b) * S(b, a);
  // 1/30 = the isotropic average over molecular orientations (Monson–McClain,
  // JCP 53:29 1970). Without it delta is ~30x too large. Matches Dalton's
  // two_photon_strengths (D = 2*Df + 4*Dg for the symmetric tensor).
  return d / 30.0;
}

/// Full Monson–McClain two-photon observables from the S tensor, matching
/// Dalton QRSMO (rspvec.F:2803-2821). `omega_ex` = excitation energy (a.u.);
/// both photons are at omega_ex/2. sigma in GM assumes Dalton's 0.1 eV FWHM.
struct Observables {
  double Df, Dg, D_linear, D_circular, R, sigma_linear_gm, sigma_circular_gm;
};

inline Observables
observables(const madness::Tensor<double> &S, double omega_ex) {
  double df = 0.0, dg = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      df += S(i, i) * S(j, j);   // Df = sum_ij S_ii S_jj / 30
      dg += S(i, j) * S(i, j);   // Dg = sum_ij S_ij^2   / 30
    }
  Observables o;
  o.Df = df / 30.0;
  o.Dg = dg / 30.0;
  o.D_linear   =  2.0 * o.Df + 4.0 * o.Dg;
  o.D_circular = -2.0 * o.Df + 6.0 * o.Dg;
  const double den = o.Df + 2.0 * o.Dg;
  o.R = (den != 0.0) ? (-o.Df + 3.0 * o.Dg) / den : 0.0;
  // a.u. -> GM (Dalton AU_TO_GM, rspvec.F:1795-6): 8*pi^2 * alpha * a0[pm]^5 /
  // (c[cm/s] * FWHM[au]); FWHM = 0.1 eV baked in. Numerically ~2.170.
  constexpr double PI      = 3.14159265358979324;
  constexpr double ALPHA   = 7.2973525693e-3;   // fine-structure constant
  constexpr double A0_PM   = 52.9177210903;     // Bohr radius in pm
  constexpr double C_CMS   = 2.99792458e10;     // speed of light, cm/s
  constexpr double FWHM_AU = 0.0036749326;      // 0.1 eV in a.u.
  const double AU_TO_GM = 8.0 * PI * PI * ALPHA * std::pow(A0_PM, 5.0) /
                          (C_CMS * FWHM_AU);
  const double ph2 = (0.5 * omega_ex) * (0.5 * omega_ex);
  o.sigma_linear_gm   = o.D_linear   * ph2 * AU_TO_GM;
  o.sigma_circular_gm = o.D_circular * ph2 * AU_TO_GM;
  return o;
}

} // namespace molresponse_v3::tpa

#endif // MOLRESPONSE_V3_KERNELS_TPA_HPP
