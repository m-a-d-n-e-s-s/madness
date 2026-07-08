#ifndef MOLRESPONSE_V3_KERNELS_ASSEMBLY_HPP
#define MOLRESPONSE_V3_KERNELS_ASSEMBLY_HPP

// =========================================================================
// State-generic assembly of the BSH right-hand sides.
//
//   assemble_theta  :  θ = V0x − E0x + γ       (FD + ES recompute)
//   assemble_lambda :  Λ = T0x + V0x − E0x_full + γ   (ES subspace matrix)
//
// These are pure linear combinations of State-shaped buffers. The pattern
// is identical across every (Type, Shell) specialization — gaxpy each
// component, then truncate — so the work belongs at the State level via
// `axpy` / `truncate_all`, not redundantly inlined in each kernel.
//
// State must provide:
//   - `State State::copy(madness::World&) const`          deep copy
//   - `void  State::axpy(madness::World&, double, const State&)`
//   - `void  State::truncate_all(madness::World&, double)`
//
// ResponseStateX<…> and ResponseStateXY<…> both satisfy this.
// =========================================================================

#include <madness/mra/mra.h>                    // FunctionDefaults

#include "../solvers/response_state.hpp"

namespace molresponse_v3 {

/// Λ = T0x + V0x − E0x_full + γ
///
/// Used by ESSolver as the symmetric subspace-projection assembly. The
/// `E0x_full` term is the diagonal-plus-off-diagonal block; in TDA this
/// is the only `E0x` representation, in Full it carries both X and Y
/// blocks.
template <typename State>
inline State assemble_lambda(madness::World &world,
                             const State &T0x,
                             const State &V0x,
                             const State &E0x_full,
                             const State &gamma) {
  State L = T0x.copy(world);
  L.axpy(world,  1.0, V0x);
  L.axpy(world, -1.0, E0x_full);
  L.axpy(world,  1.0, gamma);
  L.truncate_all(world, madness::FunctionDefaults<3>::get_thresh());
  return L;
}

/// θ = V0x − E0x + γ
///
/// The BSH right-hand side used by both FDSolver (every step) and
/// ESSolver in its recompute path. Streaming theta assembly in those
/// solvers does not call this helper — they accumulate directly into a
/// working buffer to avoid the extra deep copy here. This helper exists
/// for the non-streaming call sites and for the ES recompute path where
/// V0x, E0x, γ are already materialized.
template <typename State>
inline State assemble_theta(madness::World &world,
                            const State &V0x,
                            const State &E0x,
                            const State &gamma) {
  State T = V0x.copy(world);
  T.axpy(world, -1.0, E0x);
  T.axpy(world,  1.0, gamma);
  T.truncate_all(world, madness::FunctionDefaults<3>::get_thresh());
  return T;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_ASSEMBLY_HPP
