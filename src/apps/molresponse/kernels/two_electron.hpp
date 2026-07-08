#ifndef MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
#define MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP

// ===========================================================================
// Shared gamma contraction -- the two-electron part of the perturbed Fock
// operator (doc 15 / kernel unification).
//
// Every Kernels<Type,Shell> assembles the SAME quantity for each response
// component (the X and Y blocks of the response vector, and alpha/beta for
// open shell):
//
//     gamma = Q( J[rho]*apply_to  -  c_xc * Sum_pairs K(bra, ket)(apply_to) )
//
// The ONLY things that differ between (Type, Shell) are (a) the density factor
// + which response components build rho (handled by each kernel's two-state
// `compute_density(S1, S2)`), and (b) the list of exchange (bra, ket) pairs for
// each component (the per-type pairing). This header factors out the common
// boilerplate -- the Coulomb multiply, the c_xc*K accumulation, projection and
// truncation -- so each kernel's `apply_g` / `compute_gamma` is just "build J,
// then name the exchange pairs."
//
// The unified picture: gamma involves THREE response states. S1, S2 build the
// density AND supply the exchange bra/kets; S3 is the state we apply to. Linear
// response is the special case S1 = chi, S2 = S3 = Phi (ground state dressed as
// {phi, ...}); the VBC quadratic source is the case where S1/S2/S3 are the
// B/C-derived response states. compute_gamma is thus apply_g(chi, Phi, Phi).
//
// Exchange bra/ket order is the Dalton-validated convention
// apply_exchange(bra, ket, apply_to, lo) taken verbatim from the original
// per-type compute_gamma bodies (and pinned by tests/test_kernel_equivalence).
// ===========================================================================

#include "common_ops.hpp"   // common_ops::apply_exchange, poperatorT

#include <madness/chem/projector.h>
#include <madness/mra/mra.h>
#include <madness/world/worldprofile.h>  // PROFILE_BLOCK (perf-model meters; no-op unless WORLD_PROFILE_ENABLE)

#include <initializer_list>
#include <vector>

namespace molresponse_v3::two_electron {

using vecfuncT = std::vector<madness::real_function_3d>;

/// One exchange operator K(bra, ket) entering the gamma contraction. `bra` and
/// `ket` are the occupied/response orbital pair that builds the exchange kernel;
/// apply_gamma then acts that operator on `apply_to` (the orbitals of the
/// response component being assembled). The references must outlive the
/// apply_gamma call (they always do -- callers pass kernel-local lvalues inside
/// a braced-init-list argument).
struct ExchangePair {
  const vecfuncT &bra;
  const vecfuncT &ket;
};

/// gamma contraction, projection-free:
///     gamma = J[rho] * apply_to  -  c_xc * Sum_pairs K(bra, ket)(apply_to).
/// This is the two-electron part of the perturbed Fock operator (Coulomb minus
/// scaled exchange) acting on one response component's orbitals -- NO projection,
/// NO truncate. The linear kernels wrap this with Q + truncate (apply_gamma
/// below); VBC/beta use the raw form directly because their Fock-matrix
/// correction term contracts the UNprojected gamma against phi0 (matrix_inner(
/// phi0, gamma) would vanish under Q). `J` is the already-applied Coulomb
/// potential coulop(rho); the caller builds it once and reuses it across the
/// X/Y (and spin) components.
inline vecfuncT
apply_gamma_raw(madness::World &world,
                  const madness::real_function_3d &J,
                  const vecfuncT &apply_to,
                  std::initializer_list<ExchangePair> pairs,
                  double c_xc, double lo) {
  PROFILE_BLOCK(rs_exchange_gamma);  // exchange/γ build (J + K); the meter exchange-thread reports Tx/tile counts into
  using namespace madness;
  auto out = mul(world, J, apply_to, true);
  if (c_xc > 0.0) {
    for (const auto &p : pairs) {
      auto k = common_ops::apply_exchange(world, p.bra, p.ket, apply_to, lo);
      gaxpy(world, 1.0, out, -c_xc, k);
    }
  }
  return out;
}

/// Shared gamma contraction for every Kernels<Type,Shell>::apply_g and
/// compute_gamma. Call once per response component (X, Y; alpha, beta):
///     out = Q( J[rho] * apply_to  -  c_xc * Sum_pairs K(bra, ket)(apply_to) ).
/// Exactly apply_gamma_raw followed by projection onto the virtual space +
/// truncation (so the linear path is byte-identical to the original per-type
/// compute_gamma bodies).
inline vecfuncT
apply_gamma(madness::World &world,
              const madness::real_function_3d &J,
              const vecfuncT &apply_to,
              std::initializer_list<ExchangePair> pairs,
              const madness::QProjector<double, 3> &Q,
              double c_xc, double lo) {
  auto out = apply_gamma_raw(world, J, apply_to, pairs, c_xc, lo);
  out = Q(out);
  madness::truncate(world, out);
  return out;
}

} // namespace molresponse_v3::two_electron

#endif // MOLRESPONSE_V3_KERNELS_TWO_ELECTRON_HPP
