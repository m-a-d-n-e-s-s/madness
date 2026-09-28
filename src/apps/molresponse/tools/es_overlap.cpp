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
//   es_overlap --seed-dir=<es__KEY seed bundle> --bundle-dir=<converged es bundle> [--full]
// Prints the full |<seed_r|root_i>| matrix (all seed roots x all converged
// roots) plus the best match per seed root. --full loads Full (X,Y) bundles and
// uses the RPA metric <X|X> - <Y|Y>; default TDA (X only).

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
#include <cstdio>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;


template <typename Type>
static int run_overlap(World &world, const std::string &seed_dir,
                       const std::string &bundle_dir) {
  auto seed = load_es_roots<Type, ClosedShell>(world, seed_dir);
  auto conv = load_es_roots<Type, ClosedShell>(world, bundle_dir);
  MADNESS_CHECK(!seed.roots.empty());
  const std::size_t R = seed.roots.size(), N = conv.roots.size();

  std::vector<double> snn(R), rnn(N);
  for (std::size_t r = 0; r < R; ++r) snn[r] = rs::metric_inner(seed.roots[r], seed.roots[r]);
  for (std::size_t i = 0; i < N; ++i) rnn[i] = rs::metric_inner(conv.roots[i], conv.roots[i]);
  std::vector<std::vector<double>> ov(R, std::vector<double>(N, 0.0));
  for (std::size_t r = 0; r < R; ++r)
    for (std::size_t i = 0; i < N; ++i) {
      const double sr = rs::metric_inner(seed.roots[r], conv.roots[i]);   // collective
      ov[r][i] = (snn[r] > 0.0 && rnn[i] > 0.0) ? std::abs(sr) / std::sqrt(snn[r] * rnn[i]) : 0.0;
    }

  if (world.rank() == 0) {
    print("\n=== es_overlap ===");
    print("  seed   :", seed_dir, " n_roots=", (int)R);
    print("  bundle :", bundle_dir, " n_roots=", (int)N);
    std::printf("  converged roots: ");
    for (std::size_t i = 0; i < N; ++i)
      std::printf("  %2zu:%8.4f", i, i < (std::size_t)conv.omega.size() ? conv.omega(long(i)) : 0.0);
    std::printf("  (au)\n  |<seed_r|root_i>| (rows = seed roots, omega_seed in au):\n");
    for (std::size_t r = 0; r < R; ++r) {
      const double os = r < (std::size_t)seed.omega.size() ? seed.omega(long(r)) : 0.0;
      std::printf("  seed %2zu %8.4f :", r, os);
      std::size_t best = 0; for (std::size_t i = 1; i < N; ++i) if (ov[r][i] > ov[r][best]) best = i;
      double sumsq = 0.0; for (double v : ov[r]) sumsq += v * v;
      for (std::size_t i = 0; i < N; ++i) std::printf(" %7.4f", ov[r][i]);
      std::printf("   best=%zu (%.4f)  ||P_conv seed||=%.4f\n", best, ov[r][best], std::sqrt(sumsq));
    }
    std::printf("  (||P_conv seed|| = norm of the seed's projection onto the span of the "
                "converged roots; 1 = the seed lies in that span)\n");
  }
  return 0;
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  startup(world, argc, argv, true);
  commandlineparser parser(argc, argv);

  if (!parser.key_exists("seed-dir") || !parser.key_exists("bundle-dir")) {
    if (world.rank() == 0)
      print("Usage: es_overlap --seed-dir=<es__KEY seed bundle> "
            "--bundle-dir=<converged es bundle> [--full]");
    finalize();
    return 2;
  }
  const std::string seed_dir = parser.value_raw("seed-dir");
  const std::string bundle_dir = parser.value_raw("bundle-dir");
  if (parser.key_exists("full")) run_overlap<Full>(world, seed_dir, bundle_dir);
  else                           run_overlap<TDA>(world, seed_dir, bundle_dir);

  finalize();
  return 0;
}
