#ifndef MOLRESPONSE_V3_SOLVERS_ITERATE_HPP
#define MOLRESPONSE_V3_SOLVERS_ITERATE_HPP

// =========================================================================
// solvers::iterate — single driver template for every solver type.
//
// Any concrete solver (ClosedShellTdaSolver, ClosedShellFullSolver,
// OpenShellTdaSolver, ClosedShellFrequencyDependentSolver, ...) is
// driven by this one function. The step order is auditable here in
// one place; no solver class can re-order its own pipeline silently.
//
// The Solver concept (duck-typed, no inheritance required):
//
//   typename Solver::State
//
//   State Solver::step(State in);
//   bool  Solver::converged(const State&) const;
//
// Solvers that need optional acceleration (KAIN) or step restriction
// fold those into ::step. The driver itself stays minimal — one loop,
// one convergence check, no policy.
//
// See docs/12_solver_architecture_sketch.md.
// =========================================================================

namespace molresponse_v3::solvers {

struct IteratePolicy {
  int max_iters = 25;
};

template <typename Solver>
auto iterate(Solver &solver, typename Solver::State state,
             const IteratePolicy &policy) -> typename Solver::State {
  bool converged = false;
  for (int k = 0; k < policy.max_iters; ++k) {
    state = solver.step(state);
    if (solver.converged(state)) { converged = true; break; }
  }
  solver.print_final(state, converged);
  return state;
}

} // namespace molresponse_v3::solvers

#endif // MOLRESPONSE_V3_SOLVERS_ITERATE_HPP
