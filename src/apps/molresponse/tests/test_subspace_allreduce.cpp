// test_subspace_allreduce — keystone of the persistent-subworld recompute design
// (docs/21, docs/22 §8). Proves the user's "collective to construct the matrix":
// the subspace matrix can be assembled from per-NODE partial columns + a universe
// allreduce of the M×M SCALAR tensor — WITHOUT moving the functions to a common
// world — and equals the direct universe matrix.
//
// Model: X (M root-functions) is REPLICATED into each node-subworld. Roots are
// partitioned round-robin by node (root j -> node j % nnodes). Each node computes
// the columns it owns via a subworld-collective matrix_inner over its replicated
// X; only subworld-rank-0 contributes (the result is replicated across the
// subworld's ranks, so others would double-count). Universe-sum the partials ->
// full S. Compare to S computed directly on the universe.
//
//   mpirun -np <N> ./test_subspace_allreduce [M]
// PASS iff max|S_dist - S_direct| < 1e-9.

#include "../solvers/node_subworlds.hpp"

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>            // matrix_inner
#include <madness/world/MADworld.h>
#include <madness/world/ranks_and_hosts.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {
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
/// Node index = position of this rank's host in the (sorted) host map — the
/// same ordering on every rank, so round-robin column ownership is consistent.
int my_node_index(World &universe) {
  auto rph = ranks_per_host(universe);
  const std::string host = get_hostname();
  int idx = 0;
  for (const auto &kv : rph) { if (kv.first == host) return idx; ++idx; }
  return -1;
}
} // namespace

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);
  startup(universe, argc, argv, true);

  // CLI:  ./test_subspace_allreduce [M] [k] [thresh]   (defaults preserve old run)
  const int    M      = (argc > 1) ? std::atoi(argv[1]) : 6;   // # roots
  const int    k      = (argc > 2) ? std::atoi(argv[2]) : 6;
  const double thresh = (argc > 3) ? std::atof(argv[3]) : 1e-5;

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);

  // X: M root-functions in the universe (distributed over all ranks).
  std::vector<real_function_3d> X(M);
  for (int j = 0; j < M; ++j) {
    coord_3d c{{-6.0 + 2.0 * j, 0.4 * ((j % 3) - 1), 0.0}};
    const double a = 0.5 + 0.07 * j;
    X[j] = real_factory_3d(universe).functor(
        real_functor_3d(new GaussFunctor(c, a)));
  }
  truncate(universe, X);
  universe.gop.fence();

  // DIRECT reference: S computed on the universe.
  Tensor<double> S_direct = matrix_inner(universe, X, X);   // M x M

  // Node-subworlds; replicate X into each (one distributed copy per node, Inc 2).
  NodeSubworldInfo info;
  auto sub = make_node_aligned_subworld(universe, &info);
  const int G   = info.n_nodes;
  const int nid = my_node_index(universe);

  std::vector<real_function_3d> Xs(M);
  for (int j = 0; j < M; ++j) Xs[j] = copy(*sub, X[j]);
  sub->gop.fence();

  // Per-rank X-replication footprint = the design's extra cost (X replicated as one
  // distributed copy per node across G subworlds). size() is collective in the
  // subworld -> compute on ALL ranks, gate only the print below.
  std::size_t x_total = 0;
  for (auto &f : Xs) x_total += f.size();

  // Per-node: full S over the node's replicated X (subworld-collective), then
  // keep only the columns this node OWNS (j % G == nid).
  Tensor<double> S_local = matrix_inner(*sub, Xs, Xs);      // M x M, replicated in subworld
  Tensor<double> S_part(M, M);
  S_part.fill(0.0);
  for (int j = 0; j < M; ++j)
    if ((j % G) == nid)
      for (int i = 0; i < M; ++i) S_part(i, j) = S_local(i, j);

  // S_local is replicated across the subworld's ranks -> only rank 0 of each
  // subworld contributes, else the universe-sum double-counts by ranks/node.
  if (sub->rank() != 0) S_part.fill(0.0);
  universe.gop.fence();
  universe.gop.sum(S_part.ptr(), M * M);                   // disjoint columns -> full S
  Tensor<double> S_dist = S_part;

  double maxdiff = 0.0;
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < M; ++j)
      maxdiff = std::max(maxdiff, std::abs(S_dist(i, j) - S_direct(i, j)));
  const bool ok = (maxdiff < 1e-9);

  if (universe.rank() == 0) {
    print("=== distributed subspace-matrix allreduce  (M =", M,
          ", nodes =", G, ") ===");
    print("  built per-node partial columns, universe-summed the M×M scalars;");
    print("  Λ/functions NEVER crossed worlds.");
    print("  max|S_dist - S_direct| =", maxdiff, " (tol 1e-9)");

    const double per_rank_X_GiB = double(x_total) /
        (info.subworld_size > 0 ? info.subworld_size : 1) *
        8.0 / (1024.0 * 1024.0 * 1024.0);
    std::ostringstream os;
    os << "SUBSPACE_MEM M=" << M << " k=" << k << " thresh=" << thresh
       << " nodes=" << G << " per_rank_X_GiB=" << per_rank_X_GiB;
    print(os.str());

    print("\nSUBSPACE_ALLREDUCE_TEST:", ok ? "PASS" : "FAIL");
  }

  // Teardown: all Functions before finalize; subworld before universe.
  Xs.clear();
  sub->gop.fence();
  sub.reset();
  universe.gop.fence();
  X.clear();
  universe.gop.fence();

  finalize();
  return ok ? 0 : 1;
}
