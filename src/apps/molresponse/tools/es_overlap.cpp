// es_overlap — response-metric overlap of a SEED excited-state vector with each
// root of a converged ES bundle. Identifies which converged MADNESS root
// corresponds to the operator-selected seed by CHARACTER (max |<seed|root>|),
// not by energy — robust when the target sits in a dense/Rydberg manifold where
// energy-ordering is ambiguous. This is the "(b) overlap-based pinning" step of
// the dalton-cis-seed targeting thread: pick the state that best matches the
// cheap-calc (operator-selected) eigenvector.
//
// Overlap uses the response metric rs::metric_inner (TDA: sum_i <x_i|x_i>),
// normalized: ov_i = |<seed|root_i>| / sqrt(<seed|seed> <root_i|root_i>).
//
// Both bundles are loaded via load_es_roots; run at the write-NP, or write both
// bundles with the HDF5 backend (np-portable) and run this at any NP.
//
// Usage:
//   es_overlap --seed-dir=<es__KEY seed bundle> --bundle-dir=<converged es bundle>

#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"                 // TDA, ClosedShell
#include "../kernels/response_space_ops.hpp"   // rs::metric_inner
#include "../solvers/response_state.hpp"       // ResponseStateX
#include "../solvers/es_solver.hpp"            // ESSolver<TDA,ClosedShell>::State
#include "../solvers/es_save_load.hpp"         // load_es_roots

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <string>

using namespace madness;
using namespace molresponse_v3;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  startup(world, argc, argv, true);
  commandlineparser parser(argc, argv);

  if (!parser.key_exists("seed-dir") || !parser.key_exists("bundle-dir")) {
    if (world.rank() == 0)
      print("Usage: es_overlap --seed-dir=<es__KEY seed bundle> "
            "--bundle-dir=<converged es bundle>");
    finalize();
    return 2;
  }
  const std::string seed_dir = parser.value_raw("seed-dir");
  const std::string bundle_dir = parser.value_raw("bundle-dir");

  {
    auto seed = load_es_roots<TDA, ClosedShell>(world, seed_dir);
    auto conv = load_es_roots<TDA, ClosedShell>(world, bundle_dir);
    MADNESS_CHECK(!seed.roots.empty());

    const auto &s = seed.roots[0];
    const double snn = rs::metric_inner(s, s);
    const double seed_omega =
        (seed.omega.size() > 0) ? seed.omega(0L) : 0.0;

    int best = -1;
    double best_ov = -1.0;
    if (world.rank() == 0) {
      print("\n=== es_overlap ===");
      print("  seed   :", seed_dir, " omega=", seed_omega,
            "au (", 27.2114 * seed_omega, "eV)");
      print("  bundle :", bundle_dir, " n_roots=",
            static_cast<int>(conv.roots.size()));
      print("  root    omega(au)    omega(eV)   |<seed|root>|");
    }
    for (std::size_t i = 0; i < conv.roots.size(); ++i) {
      const auto &r = conv.roots[i];
      const double rnn = rs::metric_inner(r, r);
      const double sr = rs::metric_inner(s, r);
      const double ov =
          (snn > 0.0 && rnn > 0.0) ? std::abs(sr) / std::sqrt(snn * rnn) : 0.0;
      const double om =
          (i < static_cast<std::size_t>(conv.omega.size())) ? conv.omega(i) : 0.0;
      if (world.rank() == 0)
        std::printf("  %3zu   %11.6f  %10.3f     %.4f\n", i, om, 27.2114 * om,
                    ov);
      if (ov > best_ov) { best_ov = ov; best = static_cast<int>(i); }
    }
    if (world.rank() == 0 && best >= 0) {
      const double om = conv.omega(best);
      std::printf("\nBEST MATCH: root %d  omega=%.6f au (%.3f eV)  overlap=%.4f\n",
                  best, om, 27.2114 * om, best_ov);
      std::printf("(seed was %.3f eV; energy-lowest root is #0 at %.3f eV)\n",
                  27.2114 * seed_omega,
                  27.2114 * (conv.omega.size() > 0 ? conv.omega(0L) : 0.0));
    }
  }

  finalize();
  return 0;
}
