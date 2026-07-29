// ===========================================================================
// test_fd_kfallback_load.cpp — pins the k-CONSISTENCY contract of
// try_load_fd_state (review fix, confirmed HIGH).
//
// Scenario (the historical "functions have different k" abort): an FD leg —
// e.g. beta's A@ω_σ=2ω, closest to a pole — exists on disk ONLY at a coarse
// rung (1e-4/k6, statics-only ladder), while property assembly runs at the
// fine protocol (1e-6/k8). The coarser-rung-fallback loader used to return
// the archive's SAVED-k functions; the first mixed-k inner() in beta_abc then
// hard-aborted the whole job at the end of an expensive run.
//
// Contract under test: try_load_fd_state re-projects a non-exact source to
// the ACTIVE (k, thresh) inside the loader, so every consumer receives
// active-key functions, while exact/source_protocol_key stay the honest
// provenance record (row_accuracy reports the coarse source).
//
//   1. save a Full/ClosedShell FD state at the COARSE protocol (1e-4, k6);
//   2. switch to the FINE protocol (1e-6, k8) and try_load_fd_state:
//      -> exact == false, source_protocol_key == the coarse key,
//      -> every returned function is at k8,
//      -> inner() against a fresh k8 function succeeds (the beta contraction);
//   3. save at the fine key too and reload -> exact == true (no reprojection).
//
//   test_fd_kfallback_load        (MPI=1, fast; PASS/FAIL printed, rc reflects it)
// ===========================================================================

#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../solvers/fd_save_load.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>

using namespace madness;
using namespace molresponse_v3;

static int n_fail = 0;
#define EXPECT(cond, msg)                                                      \
  do {                                                                         \
    const bool ok = (cond);                                                    \
    if (world.rank() == 0)                                                     \
      std::printf("  %-60s %s\n", msg, ok ? "PASS" : "FAIL");                  \
    if (!ok) ++n_fail;                                                         \
  } while (0)

static double gauss0(const coord_3d &r) {
  const double x = r[0] - 0.5;
  return std::exp(-(x * x + r[1] * r[1] + r[2] * r[2]));
}
static double gauss1(const coord_3d &r) {
  const double x = r[0] + 0.5;
  return std::exp(-0.5 * (x * x + r[1] * r[1] + r[2] * r[2]));
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);

    const double coarse_t = 1e-4, fine_t = 1e-6;
    const int    coarse_k = 6,    fine_k = 8;
    const std::string coarse_key = protocol_key(coarse_t, coarse_k);
    const std::string fine_key   = protocol_key(fine_t,   fine_k);

    const Perturbation pert = Perturbation::dipole(2);   // dipole_z
    const double freq = 0.0856;                          // ω_σ = 2ω (SHG-like)

    // Rank-consistent scratch dir (every rank passes the same path into the
    // collective archive I/O); wiped up front so a stale run can't pollute.
    namespace fs = std::filesystem;
    const std::string dir =
        (fs::temp_directory_path() / "molresponse_v3_kfallback_test").string();
    if (world.rank() == 0) {
      std::error_code ec;
      fs::remove_all(dir, ec);
    }
    world.gop.fence();

    using Solver = FDSolver<Full, ClosedShell>;

    // ---- 1. save at the COARSE protocol only --------------------------------
    FunctionDefaults<3>::set_k(coarse_k);
    FunctionDefaults<3>::set_thresh(coarse_t);
    {
      Solver::State st;
      st.responses.resize(1);
      st.responses[0].x_alpha = {real_factory_3d(world).f(gauss0)};
      st.responses[0].y_alpha = {real_factory_3d(world).f(gauss1)};
      st.last_bsh_residual = {3.0e-5};
      save_fd_state<Full, ClosedShell>(world, st, dir, pert, freq,
                                       /*converged=*/true);
    }

    // ---- 2. load at the FINE protocol (coarser-rung fallback) ---------------
    FunctionDefaults<3>::set_k(fine_k);
    FunctionDefaults<3>::set_thresh(fine_t);
    if (world.rank() == 0)
      print("\n=== coarse->fine fallback load (", coarse_key, "->", fine_key,
            ") ===");
    {
      auto r = try_load_fd_state<Full, ClosedShell>(world, dir, pert, freq);
      EXPECT(bool(r), "fallback load finds the coarse source");
      if (r) {
        EXPECT(!r->exact, "exact == false (coarser source)");
        EXPECT(r->source_protocol_key == coarse_key,
               "source_protocol_key == coarse key (honest provenance)");
        bool all_k_fine = true;
        for (auto *blk : r->state.responses[0].blocks())
          for (auto &fn : *blk)
            if (fn.k() != fine_k) all_k_fine = false;
        EXPECT(all_k_fine, "loader returned ACTIVE-k functions (k8)");

        // The beta contraction that used to abort: inner() between the loaded
        // state and a function built at the active protocol.
        real_function_3d v = real_factory_3d(world).f(gauss1);
        const double ip = inner(r->state.responses[0].x_alpha[0], v);
        EXPECT(std::isfinite(ip) && std::abs(ip) > 0.0,
               "mixed-provenance inner() succeeds (no k-mismatch abort)");
      }
    }

    // ---- 3. exact load at the fine key is untouched --------------------------
    {
      Solver::State st;
      st.responses.resize(1);
      st.responses[0].x_alpha = {real_factory_3d(world).f(gauss0)};
      st.responses[0].y_alpha = {real_factory_3d(world).f(gauss1)};
      st.last_bsh_residual = {8.0e-7};
      save_fd_state<Full, ClosedShell>(world, st, dir, pert, freq,
                                       /*converged=*/true);
      auto r = try_load_fd_state<Full, ClosedShell>(world, dir, pert, freq);
      EXPECT(bool(r) && r->exact, "exact-key load reports exact == true");
      if (r)
        EXPECT(r->source_protocol_key == fine_key,
               "exact load picks the active key over the coarse one");
    }

    world.gop.fence();
    if (world.rank() == 0) {
      std::error_code ec;
      fs::remove_all(dir, ec);
      print(n_fail == 0 ? "\ntest_fd_kfallback_load: ALL PASS"
                        : "\ntest_fd_kfallback_load: FAILURES");
    }
    rc = (n_fail == 0) ? 0 : 1;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
