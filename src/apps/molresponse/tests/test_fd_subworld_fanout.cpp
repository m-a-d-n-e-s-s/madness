// test_fd_subworld_fanout — F1 of the FD state-parallel design (docs/31, 32).
//
// PROVES: independent FD (dipole) response states solved each in its OWN
// node-aligned subworld give a converged polarizability tensor α BIT-IDENTICAL
// to the single-World path. FD states are coupling-free, so the distributed
// iteration is pure  fan-out → solve → (filesystem) gather → universe assemble.
// This is the gate before any CalcManager::run() change (F2, doc 32 §5).
//
// ALSO de-risks, for free, the archive np-robustness the real design needs (the
// e544ae063 hazard, doc 32 §3) — GroundState's ctor is private (from_archive is
// the ONLY path), so the subworld GS is loaded with from_archive(*sub, …), and
// each subworld's converged FD state is saved with save_fd_state in the subworld
// then re-read by the universe in assemble_alpha. F1 tests BOTH directions of
// nio/nproc round-tripping.
//
//   mpirun -np <N> ./test_fd_subworld_fanout --archive=<moldft restart> \
//          [--thresh=1e-4] [--maxiter=25] [--calc-dir=<scratch>]
// PASS iff  max|α_sub − α_ref| < 1e-9.  Needs ≥2 nodes for a real partition;
// 1 node is a no-op partition and MUST still match exactly (free regression).
//
// Includes mirror tests/test_calc_manager_run.cpp (the single-World driver this
// A/Bs against) + tests/test_subspace_allreduce.cpp (the node-index helper).

#include "../GroundState.hpp"
#include "../ResponsePropertyPlanner.hpp"      // ResponsePlan, FDRequest
#include "../ResponseProtocol.hpp"             // set_response_protocol, protocol_key
#include "../calc/calc_executor.hpp"           // ExecutorContext, solve_fd_protocol, assemble_alpha
#include "../calc/calc_manager.hpp"            // Perturbation, NodeAction
#include "../solvers/convergence_policy.hpp"
#include "../solvers/node_subworlds.hpp"       // make_node_aligned_subworld
#include "../solvers/response_metadata.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/misc/info.h>                 // commandlineparser
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/ranks_and_hosts.h>     // ranks_per_host, get_hostname

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

// Node index = position of this rank's host in the (sorted) host map — the
// deterministic round-robin key (copied from test_subspace_allreduce.cpp:45).
int my_node_index(World &universe) {
  auto rph = ranks_per_host(universe);
  const std::string host = get_hostname();
  int idx = 0;
  for (const auto &kv : rph) { if (kv.first == host) return idx; ++idx; }
  return 0;
}

// Load the molecule from the moldft calc_info next to the archive (copied from
// test_calc_manager_run.cpp:122-141) — needed by GroundState::from_archive.
Molecule load_molecule(const std::string &archive_path) {
  Molecule molecule;
  auto dir = std::filesystem::path(archive_path).parent_path();
  for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    auto cand = dir / name;
    if (std::filesystem::exists(cand)) {
      std::ifstream ifs(cand);
      nlohmann::json j; ifs >> j;
      nlohmann::json mj;
      if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
        mj = j["tasks"][0]["molecule"];
      else if (j.contains("molecule"))
        mj = j["molecule"];
      if (!mj.is_null()) molecule.from_json(mj);
      break;
    }
  }
  return molecule;
}

std::string fock_json_near(const std::string &archive_path) {
  auto dir = std::filesystem::path(archive_path).parent_path();
  for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
    auto cand = dir / name;
    if (std::filesystem::exists(cand)) return cand.string();
  }
  return {};
}

// Read the 3×3 α tensor that assemble_alpha wrote under properties/alpha/<key>.
Tensor<double> read_alpha(const std::string &calc_dir, const std::string &key) {
  Tensor<double> a(3, 3);
  auto meta = ResponseMetadata::load_or_create(calc_dir + "/response_metadata.json");
  const auto &j = meta.json();
  // Layout (add_property, response_metadata.hpp): properties.alpha[<key>] is an
  // ARRAY of records (push_back per assemble call); each record has .alpha = the
  // 3×3 matrix. Take the most recent record. Return zeros if any level is absent
  // (the vacuous-pass guard in main then flags it as "no usable α").
  if (!j.contains("properties")) return a;
  const auto &props = j["properties"];
  if (!props.contains("alpha") || !props["alpha"].contains(key)) return a;
  const auto &arr = props["alpha"][key];
  if (!arr.is_array() || arr.empty() || !arr.back().contains("alpha")) return a;
  const auto &m = arr.back()["alpha"];
  for (long i = 0; i < 3; ++i)
    for (long jj = 0; jj < 3; ++jj) a(i, jj) = m[i][jj].get<double>();
  return a;
}

ResponsePlan dipole_plan(double thresh) {
  ResponsePlan plan;
  for (int ax = 0; ax < 3; ++ax)
    plan.fd.push_back(FDRequest{Perturbation::dipole(ax), 0.0, {thresh}});
  return plan;
}

double maxabs_diff(const Tensor<double> &a, const Tensor<double> &b) {
  double m = 0.0;
  for (long i = 0; i < a.dim(0); ++i)
    for (long j = 0; j < a.dim(1); ++j)
      m = std::max(m, std::abs(a(i, j) - b(i, j)));
  return m;
}

} // namespace

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);
  startup(universe, argc, argv, true);
  commandlineparser parser(argc, argv);

  if (!parser.key_exists("archive")) {
    if (universe.rank() == 0)
      print("Usage: test_fd_subworld_fanout --archive=<moldft restart> "
            "[--thresh=1e-4] [--maxiter=25] [--calc-dir=<scratch>]");
    finalize();
    return 2;
  }
  const std::string archive   = parser.value_raw("archive");
  const double      thresh    = parser.key_exists("thresh")
                                    ? std::stod(parser.value("thresh")) : 1e-4;
  const int         max_iters = parser.key_exists("maxiter")
                                    ? std::stoi(parser.value("maxiter")) : 25;
  const std::string base_dir  = parser.key_exists("calc-dir")
                                    ? parser.value_raw("calc-dir")
                                    : std::string("fd_subworld_fanout");
  const std::string dir_ref = base_dir + "/ref";
  const std::string dir_sub = base_dir + "/sub";
  if (universe.rank() == 0) {
    std::filesystem::create_directories(dir_ref);
    std::filesystem::create_directories(dir_sub);
  }
  universe.gop.fence();

  // Teardown discipline (matches test_es_build_subworld.cpp; doc 32 §2): scope
  // every World-bound object — GroundState (MO Functions + cached operators +
  // v_local_), the ExecutorContexts, the subworld remnants — inside this block
  // so they ALL destruct BEFORE finalize(). A Function/operator destructing
  // AFTER finalize() tries to lock a torn-down runtime mutex →
  // "RecursiveMutex::lock() failed … not an initialized mutex object"
  // (SIGABRT, rc=134), even though the A/B comparison already PASSed.
  bool ok = false;
  {
  const Molecule  mol  = load_molecule(archive);
  const std::string fj = fock_json_near(archive);

  ConvergencePolicy policy;
  policy.step_restrict_mode = ConvergencePolicy::StepRestrictMode::PerState;

  // -------------------------------------------------------------------------
  // (A) REFERENCE — solve all 3 dipole FD states in the UNIVERSE, assemble α.
  // -------------------------------------------------------------------------
  auto gs_u = GroundState::from_archive(universe, archive, mol);
  const double L = gs_u.L();
  set_response_protocol(universe, L, thresh);
  {
    const double t0 = FunctionDefaults<3>::get_thresh();
    auto cop = poperatorT(CoulombOperatorPtr(universe, gs_u.params().lo(), 0.001 * t0));
    gs_u.prepare(universe, 0.001 * t0, cop, fj);
  }
  ExecutorContext ctx_ref(universe, gs_u, L, fj,
                          ExecutorSettings{policy, PrintLevel::Normal, dir_ref, max_iters});
  for (int ax = 0; ax < 3; ++ax)
    solve_fd_protocol<Static, ClosedShell>(ctx_ref, Perturbation::dipole(ax),
                                           0.0, thresh, NodeAction::Fresh);
  universe.gop.fence();
  assemble_alpha(ctx_ref, dipole_plan(thresh), thresh);   // → dir_ref metadata
  const std::string key = protocol_key();
  const Tensor<double> alpha_ref = read_alpha(dir_ref, key);
  universe.gop.fence();

  // -------------------------------------------------------------------------
  // (B) SUBWORLD — partition the 3 states round-robin by node; each subworld
  //     loads its own GS, solves its owned states, saves to the SHARED dir_sub.
  // -------------------------------------------------------------------------
  NodeSubworldInfo info;
  auto sub = make_node_aligned_subworld(universe, &info);
  const int G   = info.n_nodes;
  const int nid = my_node_index(universe);
  {
    // S1 discipline: point the global default pmap at the subworld so everything
    // built inside (GS MOs via from_archive, sources, response vectors) is
    // subworld-local; restore the universe pmap BEFORE sub.reset().
    FunctionDefaults<3>::set_default_pmap(*sub);
    set_response_protocol(*sub, L, thresh);

    auto gs_s = GroundState::from_archive(*sub, archive, mol);   // ← np-robustness (read)
    {
      const double t0 = FunctionDefaults<3>::get_thresh();
      auto cop = poperatorT(CoulombOperatorPtr(*sub, gs_s.params().lo(), 0.001 * t0));
      gs_s.prepare(*sub, 0.001 * t0, cop, fj);
    }
    ExecutorContext ctx_sub(*sub, gs_s, L, fj,
                            ExecutorSettings{policy, PrintLevel::Normal, dir_sub, max_iters});

    for (int ax = 0; ax < 3; ++ax)
      if (ax % G == nid)                                        // round-robin partition
        solve_fd_protocol<Static, ClosedShell>(ctx_sub, Perturbation::dipole(ax),
                                               0.0, thresh, NodeAction::Fresh);
    sub->gop.fence();
  }  // ctx_sub / gs_s (subworld Functions) destruct here, sub still alive
  sub->gop.fence();
  FunctionDefaults<3>::set_default_pmap(universe);   // restore BEFORE reset
  sub.reset();
  universe.gop.fence();

  // Gather: each subworld wrote its states' archives to the shared dir_sub
  // (distinct files per axis → no collision). Assemble α in the UNIVERSE,
  // which re-reads those archives (← np-robustness, write-then-read).
  set_response_protocol(universe, L, thresh);
  ExecutorContext ctx_gather(universe, gs_u, L, fj,
                             ExecutorSettings{policy, PrintLevel::Normal, dir_sub, max_iters});
  assemble_alpha(ctx_gather, dipole_plan(thresh), thresh);      // → dir_sub metadata
  const Tensor<double> alpha_sub = read_alpha(dir_sub, key);
  universe.gop.fence();

  // -------------------------------------------------------------------------
  // A/B gate. NB: the comparison is only meaningful on a CONVERGED reference —
  // unconverged FD states are path-chaotic (not bit-stable across decomposition)
  // and, if a state never converges, assemble_alpha writes no α so read_alpha
  // returns zeros. Guard against a VACUOUS pass (0 vs 0) by requiring α_ref to
  // be non-trivial: a real h2o α is O(1-10), so a near-zero ref means "didn't
  // produce a usable α" (raise --maxiter / use a converging mol+thresh), FAIL.
  double ref_mag = 0.0;
  for (long i = 0; i < alpha_ref.dim(0); ++i)
    for (long j = 0; j < alpha_ref.dim(1); ++j)
      ref_mag = std::max(ref_mag, std::abs(alpha_ref(i, j)));
  const double diff      = maxabs_diff(alpha_sub, alpha_ref);
  const bool   ref_ok    = (ref_mag > 1e-6);     // reference actually produced an α
  const bool   match     = (diff < 1e-9);
  ok                     = ref_ok && match;
  if (universe.rank() == 0) {
    print("\n=== FD subworld fan-out A/B (nodes=", G, " ranks=", universe.size(),
          " thresh=", thresh, ") ===");
    print("  α_ref =\n", alpha_ref);
    print("  α_sub =\n", alpha_sub);
    print("  max|α_ref| =", ref_mag, " (must be > 1e-6 — else no usable α)");
    print("  max|α_sub − α_ref| =", diff, " (tol 1e-9)");
    if (!ref_ok)
      print("  WARNING: reference α ~ 0 — FD did not converge into a usable α; "
            "raise --maxiter or pick a converging mol/thresh. NOT a valid A/B.");
    print("\nFD_SUBWORLD_FANOUT_TEST:", ok ? "PASS" : "FAIL");
  }
  universe.gop.fence();
  }  // World-bound objects (gs_u, ctx_*, subworld remnants) destruct here, while
     // the runtime is still alive — see the teardown-discipline note above.
  finalize();
  return ok ? 0 : 1;
}
