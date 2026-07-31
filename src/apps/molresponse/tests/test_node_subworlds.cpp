// test_node_subworlds — Inc 1 of the ES node-φ state-parallel design (docs/21).
//
// Builds one subworld per physical node (Split_type SHARED), prints the
// host->ranks map, and asserts the "one subworld per node" invariant
// (subworld size == ranks on that host). No solver, no MRA solve — pure
// communicator-partitioning smoke test. Run multi-rank to exercise it:
//   mpirun -np <N> ./test_node_subworlds
// PASS iff every rank's subworld size matches its host's rank count.

#include "../solvers/node_subworlds.hpp"

#include <madness/world/MADworld.h>

using namespace madness;
using namespace molresponse_v3;

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);

  NodeSubworldInfo info;
  auto subworld = make_node_aligned_subworld(universe, &info);

  // Per-rank diagnostic line (each rank prints its own placement).
  print("  [rank", info.universe_rank, "/", info.universe_size, "] host=",
        info.hostname, " subworld_rank=", info.subworld_rank,
        " subworld_size=", info.subworld_size, " nodes=", info.n_nodes);
  universe.gop.fence();

  const bool ok = verify_one_subworld_per_node(universe, *subworld);

  if (universe.rank() == 0)
    print("\nNODE_SUBWORLD_TEST:", ok ? "PASS" : "FAIL");

  // The subworld World must be torn down before finalize / the universe.
  subworld.reset();
  universe.gop.fence();

  finalize();
  return ok ? 0 : 1;
}
