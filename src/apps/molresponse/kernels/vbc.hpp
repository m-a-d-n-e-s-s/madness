#ifndef MOLRESPONSE_V3_KERNELS_VBC_HPP
#define MOLRESPONSE_V3_KERNELS_VBC_HPP

// ===========================================================================
// VBC — the second-order (quadratic) response source V^{BC} for the closed-
// shell 2n+1 contractions: hyperpolarizability beta, Raman (dipole x nuclear
// displacement) and, through the single residue, two-photon absorption. No
// explicit x(2) solve anywhere; the source is only ever contracted.
//
// THE EQUATION (Hurtado/Sekino/Harrison 2026 eq. 19 = Parker/Rappoport/Furche
// 2018 eq. 28 = Salek et al. 2002 eq. 68; derivation and code dictionary in
// madness-workspace/reports/2026-09-09_orientation_derivation). For each
// ordered half (B,C), with F^B = v^B + g'[gamma^B] and
// gamma^B = sum_i |x_i^B><phi_i| + |phi_i><y_i^B|  (the response density at
// +omega_B, the object that drives the x-equation),
//
//   X channel:  V_p = sum_k x_k^C <phi_k|F^B|phi_p>        [M]  occupied matrix
//                   - Q ( F^B x_p^C )                       [A]  apply
//                   - Q ( g'[gamma_L^{BC}] phi_p )          [L]  pair density
//   Y channel:  the same with F^B -> F^B-dagger (= v^B + g'[gamma^B-dagger]),
//               x^C -> y^C, gamma_L -> gamma_L-dagger,
//
// gamma_L^{BC} = |x^B><y^C| (vv block) - |zeta_BC><phi| (oo block, fixed by
// orthonormality / idempotency; make_zeta below). The (C,B) half is added.
// The g'' (second xc kernel) slot [G] is zero for Hartree-Fock and is NOT
// implemented (DFT quadratic response is a separate project).
//
// LEG DICTIONARY (source_spec.hpp): madness::Exchange contracts
// K(bra,ket) f = sum_k ket_k Int bra_k f, so a leg {bra,ket} is the pair
// density |ket><bra|. gamma^B is therefore gamma_legs(x,y,phi) = {phi,x},{y,phi}
// -- the same order the DALTON-validated linear kernel uses
// (kernels/full.hpp compute_gamma) -- and gamma^B-dagger is gamma_dagger_legs.
// Until 2026-09-09 this file wrote the legs in reading order, (x,phi),(phi,y),
// which built every response density transposed (inherited from
// molresponse_v2/VBCMacrotask.hpp). That is invisible at omega=0 (x=y) but
// wrong at finite frequency: H2O SHG omega=0.1 vs DALTON d-aug-cc-pVQZ gave
// beta_zxx +12.9%, beta_xxz -4.0%, while the 2PA source tpa_source_spec.hpp
// (which had the equation's orientation all along) matched DALTON to ~1%.
// The oo (zeta) leg orientation is contraction-irrelevant for HF (the exchange
// of an oo pair density applied to phi is occupied and Q-projected away;
// Coulomb sees only the diagonal) but is written in the equation's orientation
// anyway.
//
// Closed-shell factor 2 lives IN each Coulomb density (v3 convention, matching
// compute_density: rho1 = 2*sum phi0*(x+y)); exchange is scaled by c_xc.
//
// Gates: tests/test_vbc_spec_equivalence (compute_vbc == Q(tpa::quadratic_source)
// on non-Hermitian stand-ins and in the Hermitian limit; static beta + Kleinman),
// tests/test_tpa_pq_vs_vbc ((P,Q) vs V^{BC} two-electron equality).
// ===========================================================================

#include "source_spec.hpp"   // declarative source engine + leg dictionary
#include "tags.hpp"
#include "tda.hpp"   // ResponseGroundState, common_ops::dot
#include "../solvers/response_state.hpp"   // ResponseStateXY<ClosedShell>

#include <madness/mra/mra.h>

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace molresponse_v3::vbc {

using vecfuncT = std::vector<madness::real_function_3d>;

/// Occupied-space relaxation term (idempotency-fixed oo block of gamma^{BC}):
///   zeta_BC[p] = - sum_q phi0[q] <y_B[q] | x_C[p]>,
/// so that -|zeta_BC><phi| has matrix elements K_pq = -<y_p^B|x_q^C>
/// (Parker eq. 24a). Cheap (one matrix_inner + transform).
inline vecfuncT
make_zeta(madness::World &world, const vecfuncT &by, const vecfuncT &cx,
          const vecfuncT &phi0) {
  using namespace madness;
  auto mat = matrix_inner(world, by, cx);
  mat.scale(-1.0);
  return transform(world, phi0, mat, true);
}

/// One ordered (B,C) half of V^{BC} as a two-channel (X,Y) spec, term for
/// term the equation in the header:
///   [L]  -Q g'[gamma_L] phi      legs  |x^B><y^C|, -|zeta_BC><phi|   (X)
///   [A]  -Q F^B x^C              legs  gamma_legs(x^B,y^B,phi), + v x^C
///   [M]  +sum_k x^C_k F^B_kp     same legs, occupied matrix <phi_k|F^B|phi_p>
/// The Y channel carries the daggered densities and the y^C target.
/// rho_B is shared by [A]/[M] and by both channels, so the engine's Coulomb
/// cache convolves it once. `v` is the raw one-electron operator of B.
inline std::vector<source_spec::SourceSpec>
vbc_half_spec(madness::World &world, const ResponseGroundState &g0,
              const vecfuncT &bx, const vecfuncT &by,
              const vecfuncT &cx, const vecfuncT &cy,
              const vecfuncT &zeta_bc,
              const madness::real_function_3d &v) {
  using namespace madness;
  using source_spec::apply_entry;
  using source_spec::gamma_dagger_legs;
  using source_spec::gamma_legs;
  using source_spec::ketbra;
  using source_spec::occupied_matrix_entry;
  const vecfuncT &phi0 = g0.amo;

  // Coulomb densities (diagonals; orientation-blind). zeta_bc carries its -1.
  real_function_3d rho_L = common_ops::dot(world, bx, cy);
  rho_L += common_ops::dot(world, phi0, zeta_bc);
  rho_L.scale(2.0);
  rho_L.truncate();

  real_function_3d rho_B = common_ops::dot(world, bx, phi0);
  rho_B += common_ops::dot(world, phi0, by);
  rho_B.scale(2.0);
  rho_B.truncate();

  source_spec::SourceSpec X, Y;
  // [L]
  X.entries.push_back(apply_entry(rho_L, {ketbra(bx, cy), ketbra(zeta_bc, phi0)},
                                  phi0, -1.0, /*project_Q=*/true));
  // [A]  (one_electron = v folds v x^C into the same entry)
  X.entries.push_back(apply_entry(rho_B, gamma_legs(bx, by, phi0),
                                  cx, -1.0, /*project_Q=*/true, v));
  // [M]
  X.entries.push_back(occupied_matrix_entry(rho_B, gamma_legs(bx, by, phi0),
                                            phi0, cx, +1.0,
                                            /*transpose_matrix=*/false, v));
  // Y channel: daggered densities, y^C target.
  Y.entries.push_back(apply_entry(rho_L, {ketbra(cy, bx), ketbra(phi0, zeta_bc)},
                                  phi0, -1.0, /*project_Q=*/true));
  Y.entries.push_back(apply_entry(rho_B, gamma_dagger_legs(bx, by, phi0),
                                  cy, -1.0, /*project_Q=*/true, v));
  Y.entries.push_back(occupied_matrix_entry(rho_B, gamma_dagger_legs(bx, by, phi0),
                                            phi0, cy, +1.0,
                                            /*transpose_matrix=*/false, v));
  return {std::move(X), std::move(Y)};
}

/// The full V^{BC} for the perturbation pair (B,C) from two CONVERGED
/// first-order states B = (x,y)(omega_B), C = (x,y)(omega_C): both zeta
/// blocks, both ordered halves, evaluated by source_spec::assemble_source.
/// VB_op / VC_op are the raw one-electron perturbation operators of B and C
/// (dipole components, or dV_nuc/dQ for Raman). Closed-shell only.
template <class Shell>
inline ResponseStateXY<Shell>
compute_vbc_spec(madness::World &world, const ResponseGroundState &g0,
                 const ResponseStateXY<Shell> &B,
                 const ResponseStateXY<Shell> &C,
                 const madness::real_function_3d &VB_op,
                 const madness::real_function_3d &VC_op) {
  using namespace madness;
  if constexpr (std::is_same_v<Shell, ClosedShell>) {
    const vecfuncT &phi0 = g0.amo;
    const vecfuncT &bx = B.x_alpha;
    const vecfuncT &by = B.y_alpha;
    const vecfuncT &cx = C.x_alpha;
    const vecfuncT &cy = C.y_alpha;

    auto zeta_bc = make_zeta(world, by, cx, phi0);
    auto zeta_cb = make_zeta(world, cy, bx, phi0);

    auto bc = source_spec::assemble_source(
        world, g0, vbc_half_spec(world, g0, bx, by, cx, cy, zeta_bc, VB_op));
    truncate(world, bc[0]);
    truncate(world, bc[1]);
    auto cb = source_spec::assemble_source(
        world, g0, vbc_half_spec(world, g0, cx, cy, bx, by, zeta_cb, VC_op));
    truncate(world, cb[0]);
    truncate(world, cb[1]);

    ResponseStateXY<ClosedShell> result;
    result.x_alpha = std::move(bc[0]);
    gaxpy(world, 1.0, result.x_alpha, 1.0, cb[0]);
    result.y_alpha = std::move(bc[1]);
    gaxpy(world, 1.0, result.y_alpha, 1.0, cb[1]);
    truncate(world, result.x_alpha);
    truncate(world, result.y_alpha);
    return result;
  } else {
    (void)world; (void)g0; (void)B; (void)C; (void)VB_op; (void)VC_op;
    throw std::runtime_error(
        "compute_vbc_spec: open-shell V^{BC} is not derived (per-spin specs "
        "are future work, step 6b/7).");
  }
}

/// The production entry point (beta, Raman, legacy 2PA arms). Since
/// 2026-09-09 this IS the spec build; the former bespoke builder
/// (compute_g / compute_vbc_i, ported from molresponse_v2) was removed
/// because it wrote the exchange legs transposed — see the header.
template <class Shell>
inline ResponseStateXY<Shell>
compute_vbc(madness::World &world, const ResponseGroundState &g0,
            const ResponseStateXY<Shell> &B,
            const ResponseStateXY<Shell> &C,
            const madness::real_function_3d &VB_op,
            const madness::real_function_3d &VC_op) {
  return compute_vbc_spec<Shell>(world, g0, B, C, VB_op, VC_op);
}

} // namespace molresponse_v3::vbc

#endif // MOLRESPONSE_V3_KERNELS_VBC_HPP
