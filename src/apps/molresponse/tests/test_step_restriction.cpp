// ===========================================================================
// test_step_restriction.cpp — pins the two ResponseSubspaceKain step-restriction
// modes (workstream C). KAIN disabled, so apply() runs ONLY step restriction.
//
// Build a 2-function state with per-function update norms ||d0||=1.0, ||d1||=0.1
// and cap maxrotn=0.5:
//   PerOrbital — f0 clamped to 0.5, f1 (0.1 < cap) untouched.
//   PerState   — total ||Δ||=sqrt(1.01)=1.005 > cap, ONE scale 0.5/1.005 applied
//                to BOTH: f0 -> 0.4975, f1 -> 0.04975.
// That f1 (0.1 vs 0.0498) is the discriminator between the two modes.
//
//   test_step_restriction        (MPI=1, fast; PASS/FAIL printed, rc reflects it)
// ===========================================================================

#include "../solvers/convergence_policy.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/response_subspace_kain.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

static double gauss0(const coord_3d &r) {
  const double x = r[0] - 0.5;
  return std::exp(-(x * x + r[1] * r[1] + r[2] * r[2]));
}
static double gauss1(const coord_3d &r) {
  const double x = r[0] + 0.5;
  return std::exp(-(x * x + r[1] * r[1] + r[2] * r[2]));
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);
    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(1e-6);

    using State = ResponseStateX<ClosedShell>;

    // Base functions and unit-norm perturbation directions.
    real_function_3d a = real_factory_3d(world).f(gauss0);
    real_function_3d b = real_factory_3d(world).f(gauss1);
    real_function_3d d0 = real_factory_3d(world).f(gauss1);
    real_function_3d d1 = real_factory_3d(world).f(gauss0);
    d0.scale(1.0 / d0.norm2());          // ||d0|| = 1.0
    d1.scale(0.1 / d1.norm2());          // ||d1|| = 0.1

    auto make_old = [&] { State s; s.x_alpha = {copy(a), copy(b)}; return s; };
    auto make_new = [&] {
      State s; s.x_alpha = {copy(a), copy(b)};
      s.x_alpha[0] += d0; s.x_alpha[1] += d1;          // updates
      return s;
    };

    const double maxrotn = 0.5;
    auto run = [&](ConvergencePolicy::StepRestrictMode mode) {
      ConvergencePolicy pol;
      pol.kain = false;                 // isolate step restriction
      pol.maxrotn = maxrotn;
      pol.step_restrict_mode = mode;
      ResponseSubspaceKain<State> kain(world, pol);
      std::vector<State> in{make_old()};
      std::vector<State> out{make_new()};
      kain.apply(in, out, 0);
      auto o = make_old();
      const double n0 = (out[0].x_alpha[0] - o.x_alpha[0]).norm2();
      const double n1 = (out[0].x_alpha[1] - o.x_alpha[1]).norm2();
      return std::pair<double, double>{n0, n1};
    };

    auto [po0, po1] = run(ConvergencePolicy::StepRestrictMode::PerOrbital);
    auto [ps0, ps1] = run(ConvergencePolicy::StepRestrictMode::PerState);

    const double tol = 5e-3;
    const double ps_scale = maxrotn / std::sqrt(1.0 * 1.0 + 0.1 * 0.1);  // 0.4975
    bool ok_po = std::abs(po0 - 0.5) < tol && std::abs(po1 - 0.1) < tol;
    bool ok_ps = std::abs(ps0 - ps_scale) < tol &&
                 std::abs(ps1 - 0.1 * ps_scale) < tol;

    if (world.rank() == 0) {
      print("\n=== step-restriction modes (maxrotn=", maxrotn, ") ===");
      print("  PerOrbital:  ||Δf0||=", po0, " (exp 0.5)   ||Δf1||=", po1,
            " (exp 0.1)   ->", ok_po ? "OK" : "BAD");
      print("  PerState:    ||Δf0||=", ps0, " (exp", ps_scale, ") ||Δf1||=", ps1,
            " (exp", 0.1 * ps_scale, ") ->", ok_ps ? "OK" : "BAD");
      bool ok = ok_po && ok_ps;
      print("\n", ok ? "PASSED" : "FAILED", " (step-restriction PerOrbital vs PerState)");
      rc = ok ? 0 : 1;
    }
    world.gop.broadcast(rc, 0);
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("EXCEPTION:", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
