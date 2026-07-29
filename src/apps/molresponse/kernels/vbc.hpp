#ifndef MOLRESPONSE_V3_KERNELS_VBC_HPP
#define MOLRESPONSE_V3_KERNELS_VBC_HPP

// ===========================================================================
// VBC — the quadratic source for hyperpolarizability beta (2n+1 / VBC
// contraction), closed-shell. Ported from molresponse_v2 SimpleVBCComputer, but
// rebuilt on v3 primitives so the two-electron kernel `compute_g` shares the
// SAME Coulomb / exchange path as the linear kernel
// Kernels<Full, ClosedShell>::compute_gamma (doc 15, beta path). beta is pure
// 2n+1 contraction — no explicit x(2) solve.
//
// compute_g now routes through two_electron::apply_gamma_raw (shared with the
// linear kernels). compute_vbc is templated over <Shell>: closed-shell is the
// v2-ported build (compute_vbc_i, alpha-only); open-shell is guarded (throws)
// until the per-spin make_zeta / compute_vbc_i kernels are derived (step 6b/7).
//
// Convention note (the whole correctness risk): v3 carries the closed-shell
// factor 2 IN the density (compute_density: rho1 = 2*sum phi0*(x+y)), so the
// Coulomb term is J[rho]*phi with no extra factor and exchange is scaled by
// c_xc (= 1 for HF). v2's compute_g instead used 2J - K on an unscaled density
// — identical for HF. We follow the v3 convention so compute_g is consistent
// with the already-Dalton-validated linear kernel.
// ===========================================================================

#include "tags.hpp"
#include "tda.hpp"   // ResponseGroundState, common_ops::{dot, apply_exchange}
#include "two_electron.hpp"  // two_electron::{ExchangePair, apply_gamma_raw}
#include "../solvers/response_state.hpp"   // ResponseStateXY<ClosedShell>

#include <madness/mra/mra.h>

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace molresponse_v3::vbc {

using vecfuncT = std::vector<madness::real_function_3d>;

/// (J[rho] - c_xc*K) applied to (phix, phiy), where the perturbed density
///   rho = 2*( sum_p Aleft[p]*Aright[p] + sum_p Bleft[p]*Bright[p] )
/// carries the closed-shell factor 2 (matching compute_density). J = coulop(rho)
/// thus needs no further factor; exchange is scaled by c_xc. The Y / conjugate
/// channel swaps the exchange bra/ket, exactly as v2 compute_g.
///
/// v2 mapping: result_x = 2*Jx - Kx, result_y = 2*Jy - Ky (HF). Here the 2 is
/// folded into rho and the K coefficient is c_xc (= 1 for HF).
inline std::pair<vecfuncT, vecfuncT>
compute_g(madness::World &world, const ResponseGroundState &g0,
          const vecfuncT &Aleft, const vecfuncT &Aright,
          const vecfuncT &Bleft, const vecfuncT &Bright,
          const vecfuncT &phix,  const vecfuncT &phiy) {
  using namespace madness;

  // rho = 2*(Aleft*Aright + Bleft*Bright)   [closed-shell factor 2 in density]
  real_function_3d rho = common_ops::dot(world, Aleft, Aright);
  rho += common_ops::dot(world, Bleft, Bright);
  rho.scale(2.0);
  rho.truncate();

  // Coulomb + exchange via the shared projection-free core. J = coulop(rho)
  // already carries the factor 2 (in rho); exchange is scaled by c_xc. Same
  // term order as the linear Kernels<Full,ClosedShell> kernel:
  //   X uses (Aleft,Aright) + (Bleft,Bright); Y swaps bra/ket (the conjugate),
  // matching the original v2 compute_g. NOT projected -- the caller (compute_vbc_i)
  // applies Qa where the response-density terms need it and leaves the Fock-matrix
  // term unprojected.
  auto J  = apply(*g0.coulop, rho);
  auto gx = two_electron::apply_gamma_raw(world, J, phix,
      {{Aleft, Aright}, {Bleft, Bright}}, g0.c_xc, g0.lo);
  auto gy = two_electron::apply_gamma_raw(world, J, phiy,
      {{Aright, Aleft}, {Bright, Bleft}}, g0.c_xc, g0.lo);
  return {std::move(gx), std::move(gy)};
}

/// Occupied-space relaxation term:
///   zeta_BC[p] = - sum_q phi0[q] <y_B[q] | x_C[p]>.
/// v2 make_zeta_bc: -1 * transform(phi0, <y_B|x_C>). Cheap (one matrix_inner +
/// transform); negligible vs. compute_g.
inline vecfuncT
make_zeta(madness::World &world, const vecfuncT &by, const vecfuncT &cx,
          const vecfuncT &phi0) {
  using namespace madness;
  auto mat = matrix_inner(world, by, cx);
  mat.scale(-1.0);
  return transform(world, phi0, mat, true);
}

/// One ordered half of the VBC source for the (B,C) pairing — a verbatim port
/// of v2 SimpleVBCComputer::compute_vbc_i, on v3 compute_g + g0.Qa. `v` is the
/// raw (one-electron) perturbation operator for B (e.g. the dipole operator for
/// B's Cartesian direction). Returns (result_x, result_y) for this half; the
/// caller sums the (B,C) and (C,B) halves.
inline std::pair<vecfuncT, vecfuncT>
compute_vbc_i(madness::World &world, const ResponseGroundState &g0,
              const vecfuncT &bx, const vecfuncT &by,
              const vecfuncT &cx, const vecfuncT &cy,
              const vecfuncT &zeta_bc,
              const madness::real_function_3d &v) {
  using namespace madness;
  const vecfuncT &phi0 = g0.amo;

  // gzeta = -Q( g[ bx,cy ; phi0,zeta_bc ](phi0,phi0) )
  auto gzeta = compute_g(world, g0, bx, cy, phi0, zeta_bc, phi0, phi0);
  vecfuncT gzeta_x = g0.Qa(gzeta.first);  scale(world, gzeta_x, -1.0);
  vecfuncT gzeta_y = g0.Qa(gzeta.second); scale(world, gzeta_y, -1.0);

  // g[B,C] response-density two-electron terms, and the phi0-only version.
  auto gbc     = compute_g(world, g0, bx, phi0, phi0, by, cx, cy);
  auto gbc_phi = compute_g(world, g0, bx, phi0, phi0, by, phi0, phi0);

  // fb = -Q( g[B,C] + v*c )
  auto vcx = mul(world, v, cx, true);
  auto vcy = mul(world, v, cy, true);
  vecfuncT fbx = gbc.first;  gaxpy(world, 1.0, fbx, 1.0, vcx); fbx = g0.Qa(fbx); scale(world, fbx, -1.0);
  vecfuncT fby = gbc.second; gaxpy(world, 1.0, fby, 1.0, vcy); fby = g0.Qa(fby); scale(world, fby, -1.0);

  // fphi = transform( c, <phi0 | g[B,C]_phi + v*phi0> )  (Fock-matrix correction)
  auto vb_phi = mul(world, v, phi0, true);
  vecfuncT fb_phi_x = gbc_phi.first;  gaxpy(world, 1.0, fb_phi_x, 1.0, vb_phi);
  vecfuncT fb_phi_y = gbc_phi.second; gaxpy(world, 1.0, fb_phi_y, 1.0, vb_phi);
  auto m_fbx = matrix_inner(world, phi0, fb_phi_x);
  auto m_fby = matrix_inner(world, phi0, fb_phi_y);
  auto fphi_x = transform(world, cx, m_fbx, true);
  auto fphi_y = transform(world, cy, m_fby, true);

  // result = truncate( gzeta + fb + fphi )
  vecfuncT result_x = gzeta_x;
  gaxpy(world, 1.0, result_x, 1.0, fbx);
  gaxpy(world, 1.0, result_x, 1.0, fphi_x);
  vecfuncT result_y = gzeta_y;
  gaxpy(world, 1.0, result_y, 1.0, fby);
  gaxpy(world, 1.0, result_y, 1.0, fphi_y);
  truncate(world, result_x);
  truncate(world, result_y);
  return {std::move(result_x), std::move(result_y)};
}

/// The full VBC quadratic source for the (B, C) perturbation pair from two
/// CONVERGED first-order states B = x(ω_B), C = x(ω_C). Builds both zeta terms
/// and sums the (B,C) and (C,B) halves (compute_vbc_i). VB_op / VC_op are the
/// raw one-electron perturbation operators for B and C. Port of v2
/// SimpleVBCComputer::compute_vbc_response (closed-shell, no x(2) solve).
template <class Shell>
inline ResponseStateXY<Shell>
compute_vbc(madness::World &world, const ResponseGroundState &g0,
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

    auto zeta_bc = make_zeta(world, by, cx, phi0);   // y_B with x_C
    auto zeta_cb = make_zeta(world, cy, bx, phi0);   // y_C with x_B

    auto bc = compute_vbc_i(world, g0, bx, by, cx, cy, zeta_bc, VB_op);
    auto cb = compute_vbc_i(world, g0, cx, cy, bx, by, zeta_cb, VC_op);

    ResponseStateXY<ClosedShell> result;
    result.x_alpha = bc.first;  gaxpy(world, 1.0, result.x_alpha, 1.0, cb.first);
    result.y_alpha = bc.second; gaxpy(world, 1.0, result.y_alpha, 1.0, cb.second);
    truncate(world, result.x_alpha);
    truncate(world, result.y_alpha);
    return result;
  } else {
    (void)world; (void)g0; (void)B; (void)C; (void)VB_op; (void)VC_op;
    throw std::runtime_error(
        "compute_vbc: open-shell VBC quadratic source not yet derived. The "
        "two-electron action is shell-generic (two_electron::apply_gamma_raw "
        "+ Kernels<Full,OpenShell>::compute_density), so the open-shell build is "
        "per-spin make_zeta + compute_vbc_i over alpha AND beta blocks -- future "
        "work (step 6b/7).");
  }
}

} // namespace molresponse_v3::vbc

#endif // MOLRESPONSE_V3_KERNELS_VBC_HPP
