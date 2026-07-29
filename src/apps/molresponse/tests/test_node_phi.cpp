// test_node_phi — Inc 2 of the ES node-φ state-parallel design (docs/21).
//
// Proves the memory claim the whole design rests on: shipping the ground state
// from the universe into NODE-ALIGNED subworlds gives ONE DISTRIBUTED copy per
// node (per-rank φ ≈ |φ|/ranks_per_node), NOT a per-rank replica (which is the
// v2 OOM cause). Also de-risks the cross-world copy outside MacroTaskQ — run
// with a timeout; a hang here would tell us the transfer needs MacroTaskQ's
// fence discipline.
//
// Method (content-independent — synthetic gaussians stand in for occupied
// orbitals): build φ in the universe (distributed over all P ranks), copy each
// into this rank's node-subworld, and compare coefficient counts.
//   - size()        : total coeffs of a function (collective in its world)
//   - size_local()  : THIS rank's local coeffs  -> the distributed/replicated tell
// For a DISTRIBUTED φ in the subworld:  Σ_{node ranks} size_local == size()  (1 copy/node)
// For a REPLICATED φ:                   Σ_{node ranks} size_local == R*size() (R copies)
//
//   mpirun -np <N> ./test_node_phi [n_orbitals]

#include "../solvers/node_subworlds.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/ranks_and_hosts.h>

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {
/// Gaussian at center c, exponent a — a cheap stand-in for an occupied orbital.
class GaussFunctor : public FunctionFunctorInterface<double, 3> {
  coord_3d c_;
  double   a_;
public:
  GaussFunctor(const coord_3d &c, double a) : c_(c), a_(a) {}
  double operator()(const coord_3d &r) const override {
    const double x = r[0] - c_[0], y = r[1] - c_[1], z = r[2] - c_[2];
    return std::exp(-a_ * (x * x + y * y + z * z));
  }
};
} // namespace

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);
  startup(universe, argc, argv, true);

  // CLI:  ./test_node_phi [N] [k] [thresh] [budget_GiB]   (defaults preserve old run)
  const int    N          = (argc > 1) ? std::atoi(argv[1]) : 12;   // ~ occupied count
  const int    k          = (argc > 2) ? std::atoi(argv[2]) : 6;
  const double thresh     = (argc > 3) ? std::atof(argv[3]) : 1e-4;
  const double budget_GiB = (argc > 4) ? std::atof(argv[4]) : 72.0; // per-task budget

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);

  // --- build φ in the UNIVERSE (distributed over all ranks) ---
  std::vector<real_function_3d> phi(N);
  for (int i = 0; i < N; ++i) {
    coord_3d c{{-8.0 + 1.6 * i, 0.5 * ((i % 3) - 1), -0.3 * ((i % 5) - 2)}};
    const double a = 0.6 + 0.05 * i;
    phi[i] = real_factory_3d(universe).functor(
        real_functor_3d(new GaussFunctor(c, a)));
  }
  truncate(universe, phi);
  universe.gop.fence();

  std::size_t uni_total = 0, uni_local = 0;
  for (auto &f : phi) { uni_total += f.size(); uni_local += f.size_local(); }

  // --- node-aligned subworlds + cross-world copy of φ into each node ---
  NodeSubworldInfo info;
  auto subworld = make_node_aligned_subworld(universe, &info);

  std::vector<real_function_3d> phi_sub(N);
  for (int i = 0; i < N; ++i) phi_sub[i] = copy(*subworld, phi[i]);
  subworld->gop.fence();

  std::size_t sub_total = 0, sub_local = 0;
  for (auto &f : phi_sub) { sub_total += f.size(); sub_local += f.size_local(); }

  // Σ local over the node's ranks: == sub_total if distributed (1 copy/node),
  // == R*sub_total if replicated per rank.
  std::size_t node_sum_local = sub_local;
  subworld->gop.sum(node_sum_local);

  const double ratio = (sub_total > 0)
                           ? double(node_sum_local) / double(sub_total)
                           : 0.0;
  const int    R     = info.subworld_size;
  // distributed: ratio ≈ 1; replicated: ratio ≈ R. Pass with generous slack.
  const int ok_local = (ratio < 1.0 + 0.5) ? 1 : 0;
  int ok = ok_local;
  universe.gop.fence();
  universe.gop.sum(&ok, 1);
  const bool all_ok = (ok == universe.size());

  // Per-rank diagnostic line.
  print("  [rank", info.universe_rank, "] host=", info.hostname,
        " uni_local=", (long)uni_local, " sub_local=", (long)sub_local,
        " (node R=", R, ")");
  universe.gop.fence();

  if (universe.rank() == 0) {
    print("\n=== node-φ memory check (", N, "orbitals,",
          info.n_nodes, "nodes ) ===");
    print("  |φ| total coeffs            =", (long)uni_total);
    print("  per-rank coeffs, universe    =", (long)(uni_total / universe.size()),
          " (one copy spread over all", universe.size(), "ranks)");
    print("  per-rank coeffs, node-sub    =", (long)(sub_total / (R > 0 ? R : 1)),
          " (one copy per node, over its", R, "ranks)");
    print("  per-rank coeffs, IF replicated=", (long)sub_total,
          " (v2 per-rank replication — what we AVOID)");
    print("  Σ(local)/size ratio per node =", ratio,
          " (≈1 distributed=1 copy/node ; ≈R replicated)");
    print("  node-φ = ONE distributed copy per node:", all_ok ? "PASS" : "FAIL");

    // Absolute per-rank φ memory at this (k, thresh) — ground-state contribution
    // only (the dominant OOM term; real RSS adds response-vector overhead). The
    // R× reduction (per_rank_node == if_replicated / R) is an EXACT structural fact,
    // independent of how realistic the synthetic gaussians are. Emit as one
    // preformatted string so key=value tokens stay space-free for grep/parsing.
    const double B = 8.0 / (1024.0 * 1024.0 * 1024.0);  // 1 coeff = 1 double -> GiB
    const double per_rank_node_GiB = double(sub_total) / (R > 0 ? R : 1) * B;
    const double if_replicated_GiB = double(sub_total) * B;
    const bool   fits = per_rank_node_GiB <= budget_GiB;
    std::ostringstream os;
    os << "NODE_PHI_MEM N=" << N << " k=" << k << " thresh=" << thresh
       << " nodes=" << info.n_nodes << " R=" << R
       << " per_rank_node_GiB=" << per_rank_node_GiB
       << " if_replicated_GiB=" << if_replicated_GiB
       << " budget_GiB=" << budget_GiB << " fits=" << (fits ? "YES" : "NO");
    print(os.str());

    print("\nNODE_PHI_TEST:", all_ok ? "PASS" : "FAIL");
  }

  // Teardown ordering matters: ALL Functions must be destroyed before
  // finalize() (and subworld Functions + the subworld World before the
  // universe), else destructors run against a torn-down runtime ->
  // "RecursiveMutex::lock() failed" abort.
  phi_sub.clear();
  subworld->gop.fence();
  subworld.reset();
  universe.gop.fence();
  phi.clear();
  universe.gop.fence();

  finalize();
  return all_ok ? 0 : 1;
}
