#ifndef MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP
#define MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP

// =========================================================================
// Node-aligned subworld creation — Inc 1 of the ES node-φ state-parallel
// design (docs/21).
//
// One subworld per PHYSICAL NODE, via MPI_Comm_split_type(SHARED): every rank
// on the same host lands in the same subworld. This is the prerequisite for
// the Stage-2 design — with node-aligned subworlds the ground state can be
// Cloud-copied into each node-subworld ONCE and left DISTRIBUTED over that
// node's ranks, giving exactly one φ copy per node (per-rank φ identical to
// the current single-World cost). That sidesteps the v2 per-RANK replication
// OOM without needing shared-memory Function storage (which MADNESS lacks).
//
// MacroTaskQ's default create_worlds uses color = rank % nsubworld, which
// INTERLEAVES subworlds across nodes (safempi Split). We want the opposite:
// node-CONTIGUOUS subworlds, which Split_type(SHARED) gives directly.
//
// This header is placement/diagnostics only — no solver behavior. The build
// fan-out that uses these subworlds is a later increment (docs/21 §5).
// =========================================================================

#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/world/safempi.h>
#include <madness/world/ranks_and_hosts.h>   // ranks_per_host, get_hostname
#include <madness/world/print.h>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Diagnostics describing how the universe split into node-aligned subworlds.
struct NodeSubworldInfo {
  std::string hostname;        ///< this rank's host
  int         universe_rank = 0;
  int         universe_size = 0;
  int         subworld_rank = 0;   ///< rank within this node's subworld
  int         subworld_size = 0;   ///< ranks in THIS subworld
  int         n_nodes       = 0;   ///< distinct hosts in the universe
  int         node_index    = 0;   ///< position of this host in the sorted host map
  // F2f (doc 32): sub-node granularity. groups_per_node subworlds per physical
  // node (1 = node-aligned). subworld_index = within-node group (0..P-1). gid =
  // GLOBAL subworld id = node_index*groups_per_node + subworld_index — the
  // deterministic round-robin partition key. n_subworlds = n_nodes*groups_per_node.
  int         groups_per_node = 1;
  int         subworld_index  = 0;
  int         gid             = 0;
  int         n_subworlds      = 1;
};

/// Create one subworld per physical node (MPI_COMM_TYPE_SHARED). All ranks on
/// the same host land in the same subworld; the `Key = universe.rank()` keeps
/// the within-node ordering deterministic. Collective on `universe`. Fills
/// `info` if non-null. The returned World must be destroyed (reset) before
/// `finalize()` and before the universe.
/// Create `groups_per_node` subworlds per physical node (doc 32 F2f). Composes:
///   1. node split (MPI_COMM_TYPE_SHARED) — guarantees no subworld straddles a
///      node, so the ground state stays node-local (one φ copy per node maximum).
///   2. within-node CONTIGUOUS split into `groups_per_node` equal rank slices.
/// Under a locality-ordered launch (`--map-by numa` / `--rank-by core`) the
/// contiguous within-node slices land on NUMA domains; otherwise they are just P
/// equal slices. groups_per_node=1 ⇒ node-aligned. Larger P trades memory (φ
/// replicated P× per node) for more state-parallelism — the small-system regime.
/// Assumes uniform ranks/node (the usual `--ntasks-per-node` launch); P is clamped
/// to the within-node rank count. Collective on `universe`. The returned World
/// must be reset() before finalize() and before the universe.
inline std::shared_ptr<madness::World>
make_subworld_pool(madness::World &universe, int groups_per_node,
                   NodeSubworldInfo *info = nullptr) {
  if (groups_per_node < 1) groups_per_node = 1;
  // Node grouping by HOSTNAME, not Split_type(SHARED): under srun --mpi=pmix
  // on Seawulf the SHARED split returns SINGLETON comms (PMIx doesn't publish
  // locality), which silently disabled fan-out — 2 ranks on one node came
  // back wn_size=1 => gpn=1 => n_subworlds=1 (gate 2082328). ranks_per_host
  // is collective and hostname-based, so it sees node topology regardless of
  // launcher; rph is a std::map (sorted keys), so the node index — and thus
  // every color below — is identical on every rank.
  const auto rph = madness::ranks_per_host(universe);   // collective: host -> [ranks]
  const std::string host = madness::get_hostname();
  int nidx = 0;
  for (const auto &kv : rph) { if (kv.first == host) break; ++nidx; }
  SafeMPI::Intracomm node_comm =
      universe.mpi.comm().Split(/*Color=*/nidx, /*Key=*/universe.rank());
  const int wn_rank = node_comm.Get_rank();
  const int wn_size = node_comm.Get_size();
  // The clamp must be UNIVERSE-uniform: with uneven ranks-per-node a per-node
  // min gives different nodes different gpn — and gid = nidx*gpn + color then
  // COLLIDES across nodes (two subworlds solving the same items and writing
  // the same archives/shards) while other items go unowned. Clamp to the
  // smallest node's rank count instead (collective).
  int min_wn = wn_size;
  universe.gop.min(min_wn);
  const int gpn     = std::min(groups_per_node, min_wn);   // uniform on every rank
  const int color   = (wn_rank * gpn) / wn_size;           // 0..gpn-1, contiguous
  SafeMPI::Intracomm sub_comm = node_comm.Split(color, /*Key=*/wn_rank);
  auto subworld = std::make_shared<madness::World>(sub_comm);
  universe.gop.fence();

  if (info) {
    info->hostname        = host;
    info->universe_rank   = universe.rank();
    info->universe_size   = universe.size();
    info->subworld_rank   = subworld->rank();
    info->subworld_size   = subworld->size();
    info->n_nodes         = static_cast<int>(rph.size());
    info->node_index      = nidx;
    info->groups_per_node = gpn;
    info->subworld_index  = color;
    info->gid             = nidx * gpn + color;            // global subworld id
    info->n_subworlds     = static_cast<int>(rph.size()) * gpn;
  }
  return subworld;
}

/// One subworld per physical node — the groups_per_node=1 case of the pool.
inline std::shared_ptr<madness::World>
make_node_aligned_subworld(madness::World &universe,
                           NodeSubworldInfo *info = nullptr) {
  return make_subworld_pool(universe, 1, info);
}

/// Print (rank 0) the host -> ranks map and check the "one subworld per node"
/// invariant: on every rank, the subworld size must equal the number of
/// universe ranks sharing that rank's host. Returns true iff the invariant
/// holds on ALL ranks. Collective on `universe`.
inline bool
verify_one_subworld_per_node(madness::World &universe,
                             madness::World &subworld) {
  auto rph = madness::ranks_per_host(universe);   // full map on every rank
  const std::string host = madness::get_hostname();
  const int expected = static_cast<int>(rph.count(host) ? rph[host].size() : 0);
  const int actual   = subworld.size();
  int ok = (expected == actual && expected > 0) ? 1 : 0;
  universe.gop.fence();
  universe.gop.sum(&ok, 1);                        // total #ranks that passed
  const bool all_ok = (ok == universe.size());

  if (universe.rank() == 0) {
    madness::print("=== node-aligned subworld map ===");
    madness::print("  universe ranks =", universe.size(),
                   "   distinct nodes =", static_cast<int>(rph.size()));
    for (const auto &kv : rph)
      madness::print("  node", kv.first, ":", static_cast<int>(kv.second.size()),
                     "ranks ->", kv.second);
    madness::print("  one-subworld-per-node:", all_ok ? "PASS" : "FAIL",
                   " (", ok, "/", universe.size(), "ranks consistent)");
  }
  return all_ok;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_NODE_SUBWORLDS_HPP
