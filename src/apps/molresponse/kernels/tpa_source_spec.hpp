#ifndef MOLRESPONSE_V3_KERNELS_TPA_SOURCE_SPEC_HPP
#define MOLRESPONSE_V3_KERNELS_TPA_SOURCE_SPEC_HPP

// ===========================================================================
// tpa_source_spec — the STATE-FREE two-photon residue source (P,Q) of the
// paper's compact form (eqs. tpa_fbar / tpa_zetabar / tpa_dbc /
// tpa_compact_P,Q), expressed as data over the declarative engine of
// kernels/source_spec.hpp.
//
// This is the promotion (source-unification step 3) of the assembly
// validated in tests/test_tpa_pq_vs_vbc.cpp: (P,Q) is built from the ground
// state and the two photon-leg responses (B,C) ONLY — no excited-state
// vector — and the 2PA moment is then the single contraction
//
//     S_bc ∝ <x^f|P^{BC}> + <y^f|Q^{BC}>            (eq:tpa_plain_master)
//
// The two-electron content reproduces the c-grouped production composition
// -(T1+T2+T3) of tpa::tpa_e3_residue to MRA precision (measured 3e-7,
// gated by tests/test_tpa_pq_vs_vbc.cpp); the optional one-electron content
// reproduces tpa::tpa_moment_residue_1e's per-ordering halves exactly.
//
// Term families, in the (validated) order they are summed, with their
// compact-equation identities [P-channel shown; Q is the conjugate]:
//
//  family B — the property-Fock pair of eq:tpa_compact_P:
//      + Σ_k x^C_k F̄^B_kp        occupied-matrix move, TRANSPOSED block
//      - g'[γ^{B†}] x^C_p         apply move, daggered legs {(φ,x^B),(y^B,φ)}
//    (P carries the kernel-daggered F̄^B of eq:tpa_fbar — "the frequency
//     gate is a single dagger"; Q carries plain F^B.)
//
//  family F — -g'[D^{BC}] φ_p with the pair density D^{BC} of eq:tpa_dbc:
//      response legs  (y^B,x^C) and (y^C,x^B)      [sign -, reoriented vs vbc]
//      relaxation legs (φ, z^{y^B x^C}) and (z^{x^B y^C}, φ)   [sign +,
//        z^{uv}_i = Σ_j φ_j <u_j|v_i> = the occupied-transposed ζ̄ / ζ blocks
//        of eq:tpa_zetabar]
//    Each leg is its own entry (own Coulomb density), matching the validated
//    op-for-op assembly; the ± of eq:tpa_dbc is the entry sign.
//
//  family D — the C-leg kernel of eq:tpa_compact_P. IDENTITY (2026-09-04,
//    term-resolved measurement + adjoint identity g'[γ†]=(g'[γ])†, Ḡ=G^T):
//    family D of ordering (B,C) EQUALS family B of ordering (C,B), entry for
//    entry — it is NOT independent content. Kept here so the per-ordering
//    P^{BC} matches eq:tpa_compact_P as published; the symmetrized builder
//    tpa_pq_spec_sym below exploits the identity (drops D, doubles B).
//      - g'[γ^{C†}] x^B_p         apply move, legs {(φ,x^C),(y^C,φ)}
//      + Σ_k x^B_k Ḡ^C_kp         occupied-matrix move, daggered legs,
//                                 Ḡ^C = <φ|g'[γ^{C†}]|φ> = (G^C)^T
//
//  family 1e (optional, both orderings) — the v^B part of F^B/F̄^B (which
//    the dagger leaves untouched):
//      - Q̂ v^B x^C_p  + Σ_k x^C_k <φ|v^B|φ>_kp     (+ the (C,B) image)
//    exactly tpa::tpa_moment_residue_1e's half(b,c) + half(c,b).
//
// Conventions: closed-shell factor 2 lives in each entry's Coulomb density;
// exchange is scaled by c_xc inside the engine (apply_gamma_raw). The √2
// eigenvector normalization stays at the contraction site, as in production.
// ===========================================================================

#include "source_spec.hpp"
#include "tags.hpp"
#include "tda.hpp"                       // ResponseGroundState
#include "../solvers/response_state.hpp" // ResponseStateXY<ClosedShell>

#include <madness/mra/mra.h>

#include <utility>
#include <vector>

namespace molresponse_v3::tpa {

using vecfuncT = std::vector<madness::real_function_3d>;

namespace pq_detail {

/// Σ_i a_i(r) b_i(r) with per-product truncation — the exact density builder
/// the validated (p,q) assembly used (test_tpa_pq_vs_vbc's vdot; NOT
/// madness::dot, whose truncation pattern differs at rounding level).
inline madness::real_function_3d
vdot(madness::World &world, const vecfuncT &a, const vecfuncT &b) {
  vecfuncT p(a.size());
  for (size_t i = 0; i < a.size(); ++i) p[i] = a[i] * b[i];
  truncate(world, p);
  return sum(world, p);
}

/// Occupied-space transfer block z^{uv}_i = Σ_j φ_j <u_j|v_i>
/// (eq:tpa_zetabar shape; transform convention of kernels/tpa.hpp).
inline vecfuncT
zblk(madness::World &world, const vecfuncT &phi,
     const vecfuncT &u, const vecfuncT &v) {
  return transform(world, phi, matrix_inner(world, u, v), true);
}

} // namespace pq_detail

/// Build the (P,Q) residue source spec for the photon pair (B,C):
/// channels[0] = P (contracts with x^f), channels[1] = Q (with y^f).
/// VB_op / VC_op are the raw one-electron perturbation operators of the two
/// photon legs; pass default-constructed Functions (or omit) to build the
/// TWO-ELECTRON part only (the object gated against tpa::tpa_e3_residue).
inline std::vector<source_spec::SourceSpec>
tpa_pq_spec(madness::World &world, const ResponseGroundState &g0,
            const ResponseStateXY<ClosedShell> &B,
            const ResponseStateXY<ClosedShell> &C,
            const madness::real_function_3d &VB_op = {},
            const madness::real_function_3d &VC_op = {}) {
  using namespace madness;
  using source_spec::apply_entry;
  using source_spec::occupied_matrix_entry;
  using pq_detail::vdot;
  using pq_detail::zblk;

  const vecfuncT &phi = g0.amo;
  const vecfuncT &xb = B.x_alpha, &yb = B.y_alpha;
  const vecfuncT &xc = C.x_alpha, &yc = C.y_alpha;

  // --- pair densities (closed-shell factor 2 IN the density) --------------
  // family B: the photon-B response density γ^B (truncated, as compute_g).
  real_function_3d rho_B = vdot(world, phi, xb);
  rho_B += vdot(world, phi, yb);
  rho_B.scale(2.0);
  rho_B.truncate();
  // family F: the four legs of D^{BC} (eq:tpa_dbc), one density per leg,
  // untruncated — matching the validated per-leg Coulomb builds.
  const vecfuncT z1 = zblk(world, phi, yb, xc);   // z^{y^B x^C}
  const vecfuncT z2 = zblk(world, phi, xb, yc);   // z^{x^B y^C}
  real_function_3d rho_cb = vdot(world, xc, yb);  rho_cb.scale(2.0);
  real_function_3d rho_bc = vdot(world, xb, yc);  rho_bc.scale(2.0);
  real_function_3d rho_z1 = vdot(world, z1, phi); rho_z1.scale(2.0);
  real_function_3d rho_z2 = vdot(world, phi, z2); rho_z2.scale(2.0);
  // family D: the C-leg density γ^C (untruncated, as the validated W build).
  real_function_3d rho_C = vdot(world, xc, phi);
  rho_C += vdot(world, yc, phi);
  rho_C.scale(2.0);

  source_spec::SourceSpec P, Q;

  // --- family B:  Σ_k x^C_k F̄^B_kp  -  g'[γ^{B†}] x^C  (P);  conjugate (Q)
  P.entries.push_back(occupied_matrix_entry(rho_B, {{xb, phi}, {phi, yb}},
                                            phi, xc, +1.0,
                                            /*transpose=*/true));
  P.entries.push_back(apply_entry(rho_B, {{phi, xb}, {yb, phi}},
                                  xc, -1.0, /*Q=*/false));
  Q.entries.push_back(apply_entry(rho_B, {{xb, phi}, {phi, yb}},
                                  yc, -1.0, /*Q=*/false));
  Q.entries.push_back(occupied_matrix_entry(rho_B, {{phi, xb}, {yb, phi}},
                                            phi, yc, +1.0,
                                            /*transpose=*/true));

  // --- family F:  -g'[D^{BC}] φ  (P)  /  -g'[D^{BC†}] φ  (Q), leg by leg
  P.entries.push_back(apply_entry(rho_cb, {{yb, xc}},  phi, -1.0, false));
  P.entries.push_back(apply_entry(rho_z1, {{phi, z1}}, phi, +1.0, false));
  P.entries.push_back(apply_entry(rho_bc, {{yc, xb}},  phi, -1.0, false));
  P.entries.push_back(apply_entry(rho_z2, {{z2, phi}}, phi, +1.0, false));
  Q.entries.push_back(apply_entry(rho_cb, {{xc, yb}},  phi, -1.0, false));
  Q.entries.push_back(apply_entry(rho_z1, {{z1, phi}}, phi, +1.0, false));
  Q.entries.push_back(apply_entry(rho_bc, {{xb, yc}},  phi, -1.0, false));
  Q.entries.push_back(apply_entry(rho_z2, {{phi, z2}}, phi, +1.0, false));

  // --- family D:  -g'[γ^{C†}] x^B + Σ_k x^B_k G^C_kp  (P);  conjugate (Q)
  P.entries.push_back(apply_entry(rho_C, {{phi, xc}, {yc, phi}},
                                  xb, -1.0, false));
  P.entries.push_back(occupied_matrix_entry(rho_C, {{phi, xc}, {yc, phi}},
                                            phi, xb, +1.0,
                                            /*transpose=*/false));
  Q.entries.push_back(apply_entry(rho_C, {{xc, phi}, {phi, yc}},
                                  yb, -1.0, false));
  Q.entries.push_back(occupied_matrix_entry(rho_C, {{xc, phi}, {phi, yc}},
                                            phi, yb, +1.0,
                                            /*transpose=*/false));

  // --- family 1e (optional): -Q̂ v x + Σ_k x_k <φ|v|φ>_kp, both orderings —
  // exactly tpa_moment_residue_1e's half(b,c) + half(c,b).
  if (VB_op.is_initialized()) {
    P.entries.push_back(apply_entry({}, {}, xc, -1.0, /*Q=*/true, VB_op));
    P.entries.push_back(occupied_matrix_entry({}, {}, phi, xc, +1.0,
                                              false, VB_op));
    Q.entries.push_back(apply_entry({}, {}, yc, -1.0, /*Q=*/true, VB_op));
    Q.entries.push_back(occupied_matrix_entry({}, {}, phi, yc, +1.0,
                                              false, VB_op));
  }
  if (VC_op.is_initialized()) {
    P.entries.push_back(apply_entry({}, {}, xb, -1.0, /*Q=*/true, VC_op));
    P.entries.push_back(occupied_matrix_entry({}, {}, phi, xb, +1.0,
                                              false, VC_op));
    Q.entries.push_back(apply_entry({}, {}, yb, -1.0, /*Q=*/true, VC_op));
    Q.entries.push_back(occupied_matrix_entry({}, {}, phi, yb, +1.0,
                                              false, VC_op));
  }

  return {std::move(P), std::move(Q)};
}

/// SYMMETRIZED residue source: P_sym = P^{BC} + P^{CB} (and Q_sym), built
/// via the family-D collapse (see family-D note above): family D of one
/// ordering equals family B of the other, so the two-ordering sum needs only
/// family B (signs doubled) of both orderings plus the D^{BC}/D^{CB} legs —
/// 12 two-electron entries instead of 16 (saves 4 exchange applies per
/// source). Implemented by post-processing tpa_pq_spec so the two builders
/// cannot drift: entries are taken from the validated per-ordering tables.
/// 1e entries are emitted ONCE (a single tpa_pq_spec call already carries
/// both orderings' 1e content). Gate: test_tpa_pq_vs_vbc gate 4 asserts
/// equality with tpa_pq_spec(B,C) + tpa_pq_spec(C,B) (2e part) at
/// summation-order precision.
inline std::vector<source_spec::SourceSpec>
tpa_pq_spec_sym(madness::World &world, const ResponseGroundState &g0,
                const ResponseStateXY<ClosedShell> &B,
                const ResponseStateXY<ClosedShell> &C,
                const madness::real_function_3d &VB_op = {},
                const madness::real_function_3d &VC_op = {}) {
  auto bc = tpa_pq_spec(world, g0, B, C, VB_op, VC_op);
  auto cb = tpa_pq_spec(world, g0, C, B);  // 2e only
  // per-channel entry layout (2e): [0,1]=family B, [2..5]=family F (D legs),
  // [6,7]=family D (dropped — equals the OTHER ordering's family B),
  // [8..]=1e (bc only).
  std::vector<source_spec::SourceSpec> sym(2);
  for (int ch = 0; ch < 2; ++ch) {
    auto &out = sym[static_cast<size_t>(ch)].entries;
    for (int i : {0, 1}) {                       // family B, both orderings, x2
      auto e = bc[static_cast<size_t>(ch)].entries[static_cast<size_t>(i)];
      e.sign *= 2.0;
      out.push_back(std::move(e));
      auto f = cb[static_cast<size_t>(ch)].entries[static_cast<size_t>(i)];
      f.sign *= 2.0;
      out.push_back(std::move(f));
    }
    for (int i : {2, 3, 4, 5}) {                 // family F, both orderings
      out.push_back(bc[static_cast<size_t>(ch)].entries[static_cast<size_t>(i)]);
      out.push_back(cb[static_cast<size_t>(ch)].entries[static_cast<size_t>(i)]);
    }
    for (size_t i = 8; i < bc[static_cast<size_t>(ch)].entries.size(); ++i)
      out.push_back(bc[static_cast<size_t>(ch)].entries[i]);   // 1e, once
  }
  return sym;
}

/// Evaluate the (P,Q) spec into a response-shaped pair (x_alpha = P,
/// y_alpha = Q), truncated. With VB_op/VC_op supplied this is the FULL
/// state-free residue right-hand side of eq:tpa_compact_P/Q (1e + 2e); with
/// them omitted it is the two-electron part only, whose contraction
/// <x^f|P> + <y^f|Q> equals tpa_e3_residue / √2  (== -(T1+T2+T3)).
inline ResponseStateXY<ClosedShell>
assemble_tpa_pq(madness::World &world, const ResponseGroundState &g0,
                const ResponseStateXY<ClosedShell> &B,
                const ResponseStateXY<ClosedShell> &C,
                const madness::real_function_3d &VB_op = {},
                const madness::real_function_3d &VC_op = {}) {
  auto pq = source_spec::assemble_source(
      world, g0, tpa_pq_spec(world, g0, B, C, VB_op, VC_op));
  ResponseStateXY<ClosedShell> out;
  out.x_alpha = std::move(pq[0]);
  out.y_alpha = std::move(pq[1]);
  truncate(world, out.x_alpha);
  truncate(world, out.y_alpha);
  return out;
}

} // namespace molresponse_v3::tpa

#endif // MOLRESPONSE_V3_KERNELS_TPA_SOURCE_SPEC_HPP
