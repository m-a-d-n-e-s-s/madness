#ifndef MOLRESPONSE_V3_SOLVERS_FD_TARGET_HPP
#define MOLRESPONSE_V3_SOLVERS_FD_TARGET_HPP

// =========================================================================
// FDProblem — read-only problem definition for the frequency-dependent
// solver. Wraps the ground-state-side `ResponseGroundState` (phi0, eps,
// V_local, fock, coulop, Q, c_xc, lo) that the per-root kernels already
// use, and adds a list of `FDPerturbation` entries — one per
// (perturbation, frequency) combination the solver should iterate.
//
// Symmetric with ESProblem<Type, Shell> in es_solver.hpp.
// =========================================================================

#include "../kernels/tda.hpp"          // ResponseGroundState, vecfuncT, poperatorT
#include "../kernels/tags.hpp"
#include "response_state.hpp"

#include <madness/mra/mra.h>

#include <vector>

namespace molresponse_v3 {

/// One frequency-dependent channel: a fixed omega and the precomputed
/// perturbation source `V_pert · phi0` packed into the same Storage
/// shape the solver iterates against. For closed-shell static/dipole
/// perturbations, source_x == source_y at construction (caller's
/// concern); for Full at omega != 0 they evolve independently inside
/// step() because the BSH operators differ on ± omega.
template <typename Shell>
struct FDPerturbationX {
  double omega = 0.0;
  ResponseStateX<Shell> source;  // V_pert · phi0  (X-only carrier)
};
template <typename Shell>
struct FDPerturbationXY {
  double omega = 0.0;
  ResponseStateXY<Shell> source; // V_pert · phi0  (X+Y carrier)
};

/// Pick the FDPerturbation shape that matches the Storage shape for a
/// given (Type, Shell). FD only supports Static and Full — TDA-FD is
/// not a thing in this code base (TDA is reserved for the excited-state
/// eigenvalue problem). Trying to instantiate FDSolver<TDA, ...> fails
/// at compile time because no FDPerturbationOf<TDA, *> mapping exists.
template <typename Type, typename Shell> struct FDPerturbationOf;
template <typename Shell> struct FDPerturbationOf<Static, Shell> {
  using type = FDPerturbationX<Shell>;
};
template <typename Shell> struct FDPerturbationOf<Full, Shell> {
  using type = FDPerturbationXY<Shell>;
};
template <typename Type, typename Shell>
using FDPerturbationOf_t = typename FDPerturbationOf<Type, Shell>::type;

/// Read-only problem definition for FDSolver<Type, Shell>.
///   gs        — the prepared ground state (zeroth-order data) the
///               kernels read from. Same object ESSolver uses.
///   responses — one entry per (perturbation, frequency) combination
///               the solver should iterate. State::responses[r] gives
///               the current response state for responses[r].
template <typename Type, typename Shell>
struct FDProblem {
  ResponseGroundState gs;
  std::vector<FDPerturbationOf_t<Type, Shell>> responses;

  int n_responses() const { return static_cast<int>(responses.size()); }
};

// TODO(rename): when the FD design stabilises, rename `ResponseGroundState` ->
// `KernelTarget` since it carries the shared ground-state-side kernel
// inputs (not ES-specific). Defer the rename until after Phase 1c so
// it lands as one isolated cleanup commit.

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_TARGET_HPP
