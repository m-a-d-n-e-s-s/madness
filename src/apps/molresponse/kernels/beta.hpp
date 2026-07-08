#ifndef MOLRESPONSE_V3_KERNELS_BETA_HPP
#define MOLRESPONSE_V3_KERNELS_BETA_HPP

// ===========================================================================
// beta — hyperpolarizability assembly (2n+1 contraction), closed-shell. Given
// the converged first-order states and the prebuilt VBC quadratic source, one
// tensor element beta_{A,B,C}. Verbatim port of the contraction core in
// molresponse_v2 PropertyManager::compute_beta (the three inner-product terms +
// the -2 factor), on v3 types + vbc::make_zeta.
//
// Flat structure (v2 DynamicRestrictedResponse = [x_alpha | y_alpha]):
//   beta = -2 * ( b1 + b2 + b3 )
//   b1 = <(-x_A, -y_A) | (VBC.x, VBC.y)>        = -( <x_A|VBC.x> + <y_A|VBC.y> )
//   b2 = < VA * (x_C.y, zeta_bc) | (x_B.x, phi0) >
//      = <VA*x_C.y | x_B.x> + <VA*zeta_bc | phi0>,  zeta_bc = make_zeta(y_B, x_C)
//   b3 = < VA * (x_B.y, zeta_cb) | (x_C.x, phi0) >
//      = <VA*x_B.y | x_C.x> + <VA*zeta_cb | phi0>,  zeta_cb = make_zeta(y_C, x_B)
//
// xA is the A response (X,Y) at omega_sigma = omega_B + omega_C (the response
// being differentiated); for a static A (omega_sigma = 0) the loader supplies
// (x, x). VA is the raw A dipole operator. B, C are the converged first-order
// states at omega_B, omega_C. No explicit x(2) solve.
// ===========================================================================

#include "vbc.hpp"   // vbc::make_zeta, ResponseGroundState, ResponseStateXY

#include <madness/mra/mra.h>

#include <stdexcept>
#include <type_traits>

namespace molresponse_v3::beta {

using vecfuncT = std::vector<madness::real_function_3d>;

/// One beta tensor element, templated over Shell. The CLOSED-SHELL branch is the
/// verbatim v2-ported contraction (byte-identical). OPEN-SHELL is guarded: the
/// open-shell gamma exists (Kernels<Full,OpenShell>::apply_g), but the open-shell
/// make_zeta + per-spin contraction kernels are not yet derived (step 6b/7).
template <class Shell>
inline double
beta_abc(madness::World &world, const ResponseGroundState &g0,
         const ResponseStateXY<Shell> &xA,
         const ResponseStateXY<Shell> &vbc,
         const ResponseStateXY<Shell> &B,
         const ResponseStateXY<Shell> &C,
         const madness::real_function_3d &VA_op) {
  using namespace madness;
  if constexpr (std::is_same_v<Shell, ClosedShell>) {
    const vecfuncT &phi0 = g0.amo;

    // b1 = -( <x_A | VBC.x> + <y_A | VBC.y> )
    const double b1 = -(inner(xA.x_alpha, vbc.x_alpha) +
                        inner(xA.y_alpha, vbc.y_alpha));

    // b2: zeta_bc = make_zeta(y_B, x_C); <VA*x_C.y | x_B.x> + <VA*zeta_bc | phi0>
    auto zeta_bc = vbc::make_zeta(world, B.y_alpha, C.x_alpha, phi0);
    auto va_xcy  = mul(world, VA_op, C.y_alpha, true);
    auto va_zbc  = mul(world, VA_op, zeta_bc,   true);
    const double b2 = inner(va_xcy, B.x_alpha) + inner(va_zbc, phi0);

    // b3: zeta_cb = make_zeta(y_C, x_B); <VA*x_B.y | x_C.x> + <VA*zeta_cb | phi0>
    auto zeta_cb = vbc::make_zeta(world, C.y_alpha, B.x_alpha, phi0);
    auto va_xby  = mul(world, VA_op, B.y_alpha, true);
    auto va_zcb  = mul(world, VA_op, zeta_cb,   true);
    const double b3 = inner(va_xby, C.x_alpha) + inner(va_zcb, phi0);

    return -2.0 * (b1 + b2 + b3);
  } else {
    (void)world; (void)g0; (void)xA; (void)vbc; (void)B; (void)C; (void)VA_op;
    throw std::runtime_error(
        "beta_abc: open-shell beta contraction not yet derived. Open-shell gamma "
        "is available (Kernels<Full,OpenShell>::apply_g); the open-shell make_zeta "
        "+ per-spin contraction kernels are future work (step 6b/7).");
  }
}

} // namespace molresponse_v3::beta

#endif // MOLRESPONSE_V3_KERNELS_BETA_HPP
