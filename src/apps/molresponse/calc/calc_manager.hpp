#ifndef MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP
#define MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP

// ===========================================================================
// CalcManager — the scheduling core (doc 15, increment 15a).
//
// This header is the WORLD-FREE, PURE half of the calc manager: it turns a
// ResponsePlan (doc 14 Tier A) into a DAG of CalcNodes, reconciles each node
// against what is already on disk (response_metadata.json), and emits the
// two-phase, perturbation-parallel schedule of (node, protocol) WorkItems.
//
// Nothing here touches madness::World, MPI, MRA functions, or the archive.
// Every function takes the aggregate metadata as a const nlohmann::json& so a
// test can pass a literal blob. The single World-bound seam is the abstract
// ICalcExecutor (declared at the bottom): the real FD/ES executor and
// CalcManager::run() live in calc/calc_executor.hpp (next increment) and are
// validated on the allocation.
//
// Design decisions vs doc 15:
//   * build_dag takes `int n_atoms`, not a full Molecule. The only reason the
//     planner needs the molecule is to expand the nuclear_all sentinel into
//     one node per atom; n_atoms is all that requires. CalcManager::build()
//     (executor side) extracts natom() from the Molecule and forwards it. This
//     keeps the whole header molecule-independent and trivially testable.
//   * reconcile_protocol / prerequisites_converged / schedule are keyed off the protocol step THRESHOLD
//     (a double), not a pre-built protocol_key string. k is derived once via
//     default_k_for_thresh, so callers only thread the ramp of thresholds.
//     (15a assumes the standard thresh->k table; an override_k ramp is a
//     later concern — see protocol_key_at().)
//   * prerequisites_converged resolves a node's prerequisites against the DAG (id -> node) so
//     it can reconstruct each prerequisite's metadata identity. In 15a no
//     concrete schedulable node carries a hard edge (DerivedFD "*" is a
//     symbolic expansion template, excluded from scheduling until ES
//     converges); prerequisites_converged becomes load-bearing for VBC in 15b and is
//     implemented + tested now so 15b inherits a working gate.
//   * The reconcile diverged-split reads `entry.value("diverged", false)`
//     defensively: the FD/ES save paths record `converged` today but not
//     `diverged`. Recording it is a one-line follow-up on the persistence
//     side; until then a missing field reads as not-diverged (-> Resume),
//     which is the safe default.
// ===========================================================================

#include "../Perturbations.hpp"          // Perturbation
#include "../ResponseProtocol.hpp"        // protocol_key, default_k_for_thresh
#include "../ResponsePropertyPlanner.hpp" // ResponsePlan, FDRequest, ESRequest, ...
#include "../solvers/response_metadata.hpp" // ResponseMetadata::freq_key

#include <madness/external/nlohmann_json/json.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Vocabulary (scheduler term -> response physics):
//   CalcNode       a response STATE to compute: x(ω) for a perturbation (FD),
//                  or an excited-state bundle (ES).
//   perturbation   the operator defining a unit of work (dipole_a, nuc_iα); the
//                  unit of parallelism — distinct perturbations are independent.
//                  (An ES bundle is its own unit, keyed by its id.)
//   seed_state     per perturbation, the state solved first (static ω=0 if
//                  requested, else lowest |ω|) that seeds frequency-continuation
//                  for the rest — solve α(0), continue to α(ω).
//   protocol step  one (thresh, k) MRA accuracy level on the coarse->fine ramp.
//   wave           a set of states with no prerequisite between them — safe to
//                  solve concurrently.
//   prerequisite   a state that must be CONVERGED before this one (correctness).
//   seed_from      a soft starting-guess hint (efficiency; falls back to fresh).
// ---------------------------------------------------------------------------

namespace molresponse_v3 {

// ---------------------------------------------------------------------------
// Schedulable unit.
// ---------------------------------------------------------------------------

enum class CalcKind {
  FD,         // concrete frequency-dependent dipole solve
  NuclearFD,  // nuclear-displacement FD (shares the fd_states subtree)
  ES,         // excited-state bundle
  DerivedFD,  // FD whose frequency comes from a converged ES root
  VBC,        // 15b: quadratic source (placeholder; not emitted yet)
};

/// One response state with its own protocol ladder. `id` is a stable string
/// that doubles as the dependency key (implicit-by-identity — see doc 15).
struct CalcNode {
  CalcKind            kind = CalcKind::FD;
  std::string         id;

  // FD / NuclearFD / DerivedFD payload
  Perturbation        pert;
  double              freq = 0.0;   // DerivedFD: 0 until expanded from a root

  // VBC payload: the C-side (B uses pert/freq above).
  Perturbation        pert_c;
  double              freq_c = 0.0;

  // ES payload
  bool                tda = true;
  int                 n_roots = 0;

  std::vector<double> protocols;    // coarse -> fine ladder (from the plan)

  // DerivedFD symbolic reference: "*" (all roots) or "es_root_0003".
  std::string         es_root_id;

  // DerivedFD: FD frequency = es_freq_factor * (ES root energy). 0.5 keeps the
  // FD off the response pole at ω = ωₙ (two-photon resonance, 2ω = ωₙ). 1.0 for
  // ordinary nodes (unused there).
  double              es_freq_factor = 1.0;

  // Hard edges (correctness): node ids that must be converged before this
  // node runs. Resolved by prerequisites_converged against the DAG + metadata.
  std::vector<std::string> prerequisites;

  // Soft edge: a seed HINT (node id), never required for correctness. The
  // executor always picks the nearest converged state via try_load_fd_state's
  // restart precedence; this only biases intra-perturbation ordering.
  std::string         seed_from;

  /// A symbolic node is an expansion template, not directly schedulable:
  /// a DerivedFD that still references "*" (or no concrete frequency yet).
  bool is_symbolic() const {
    return kind == CalcKind::DerivedFD && es_root_id == "*";
  }
};

/// Reconcile verdict for one (node, protocol step). See the doc-15 table.
enum class NodeAction { Skip, Restart, Resume, Fresh };

/// Stringify a NodeAction for logs / the R1c scheduler trace (Output.diagnostics).
inline const char *node_action_name(NodeAction a) {
  switch (a) {
    case NodeAction::Skip:    return "skip";
    case NodeAction::Restart: return "restart";
    case NodeAction::Resume:  return "resume";
    case NodeAction::Fresh:   return "fresh";
  }
  return "?";
}

/// One (node, protocol) step of work.
struct WorkItem {
  const CalcNode *node = nullptr;
  double          thresh = 0.0;   // this protocol step's truncation threshold
  NodeAction      action = NodeAction::Fresh;
};

// ---------------------------------------------------------------------------
// Identity helpers (free functions — also used by the executor later).
// ---------------------------------------------------------------------------

/// Canonical protocol_key for a protocol step, derived from the threshold alone using
/// the standard thresh->k table. (15a: override_k ramps unsupported.)
inline std::string protocol_key_at(double thresh) {
  return protocol_key(thresh, default_k_for_thresh(thresh));
}

inline std::string fd_node_id(const Perturbation &pert, double freq) {
  return "fd:" + pert.description() + "@" + ResponseMetadata::freq_key(freq);
}
inline std::string es_node_id(bool tda, int n_roots) {
  return std::string("es:") + (tda ? "tda" : "full") + "_n"
       + std::to_string(n_roots);
}
inline std::string derived_fd_node_id(const Perturbation &pert,
                                       const std::string &es_root_id) {
  return "dfd:" + pert.description() + "@" + es_root_id;
}
inline std::string vbc_node_id(const Perturbation &pb, const Perturbation &pc,
                               double fb, double fc) {
  return "vbc:" + pb.description() + "__" + pc.description() + "@"
       + ResponseMetadata::freq_key(fb) + "_" + ResponseMetadata::freq_key(fc);
}

/// The perturbation a node belongs to — the unit of parallelism (distinct
/// perturbations are independent and share every wave). FD / NuclearFD /
/// DerivedFD key on their perturbation operator; an ES bundle is its own unit,
/// so it keys on its own id.
inline std::string perturbation_key(const CalcNode &n) {
  if (n.kind == CalcKind::ES || n.kind == CalcKind::VBC) return n.id;
  return n.pert.description();
}

// ---------------------------------------------------------------------------
namespace detail_calc {

/// Sorted (coarse -> fine) deduplicated protocol ramp.
inline std::vector<double> sorted_ramp(std::vector<double> r) {
  std::sort(r.begin(), r.end(), std::greater<double>{});
  r.erase(std::unique(r.begin(), r.end(),
                      [](double a, double b) {
                        return std::abs(a - b) <= 1e-15 * std::max(a, b);
                      }),
          r.end());
  return r;
}

inline bool ramp_contains(const std::vector<double> &ramp, double thresh) {
  for (double t : ramp)
    if (std::abs(t - thresh) <= 1e-12 * std::max(t, thresh)) return true;
  return false;
}

inline bool entry_converged(const nlohmann::json &e) {
  return e.value("converged", false);
}
inline bool entry_diverged(const nlohmann::json &e) {
  return e.value("diverged", false);
}

/// FD status entry at an exact (pert, key, freq) or nullptr.
inline const nlohmann::json *
fd_entry(const nlohmann::json &meta, const std::string &pert,
         const std::string &key, const std::string &fkey) {
  if (!meta.contains("fd_states")) return nullptr;
  const auto &fs = meta["fd_states"];
  if (!fs.contains(pert) || !fs[pert].contains(key) ||
      !fs[pert][key].contains(fkey))
    return nullptr;
  return &fs[pert][key][fkey];
}

// FD coarser-source selection now lives in one shared helper —
// best_usable_fd_source_key (response_metadata.hpp) — used by both the
// reconcile verdict here and the archive load in try_load_fd_state, so the two
// can never disagree. (Replaced the old converged-only has_coarser_converged_fd,
// which missed coarser PARTIALS and didn't exclude diverged seeds.)

/// ES bundle entry at a key, or nullptr.
inline const nlohmann::json *es_entry(const nlohmann::json &meta,
                                      const std::string &key) {
  if (!meta.contains("excited_states") ||
      !meta["excited_states"].contains(key))
    return nullptr;
  return &meta["excited_states"][key];
}

/// VBC status entry at vbc_states/<id>/<key>, or nullptr.
inline const nlohmann::json *
vbc_entry(const nlohmann::json &meta, const std::string &id,
          const std::string &key) {
  if (!meta.contains("vbc_states")) return nullptr;
  const auto &vs = meta["vbc_states"];
  if (!vs.contains(id) || !vs[id].contains(key)) return nullptr;
  return &vs[id][key];
}

inline bool has_coarser_converged_vbc(const nlohmann::json &meta,
                                      const std::string &id,
                                      double target_thresh, int target_k) {
  if (!meta.contains("vbc_states") || !meta["vbc_states"].contains(id))
    return false;
  if (!meta.contains("protocols")) return false;
  const auto &protos = meta["protocols"];
  for (const auto &[k2, ent] : meta["vbc_states"][id].items()) {
    if (!entry_converged(ent)) continue;
    if (!protos.contains(k2)) continue;
    const double t  = protos[k2].value("thresh", 0.0);
    const int    kk = protos[k2].value("k", 0);
    if (t >= target_thresh && kk <= target_k) return true;
  }
  return false;
}

inline bool has_coarser_converged_es(const nlohmann::json &meta,
                                     double target_thresh, int target_k) {
  if (!meta.contains("excited_states") || !meta.contains("protocols"))
    return false;
  const auto &protos = meta["protocols"];
  for (const auto &[key, ent] : meta["excited_states"].items()) {
    if (!entry_converged(ent)) continue;
    if (!protos.contains(key)) continue;
    const double t  = protos[key].value("thresh", 0.0);
    const int    kk = protos[key].value("k", 0);
    if (t >= target_thresh && kk <= target_k) return true;
  }
  return false;
}

} // namespace detail_calc

// ---------------------------------------------------------------------------
// build_dag — ResponsePlan -> nodes. Pure.
// ---------------------------------------------------------------------------

/// Expand a ResponsePlan into schedulable nodes:
///   * fd          -> one FD node each; dynamic ω gets seed_from = the same
///                    perturbation's static (ω=0) node if the plan has one.
///   * es          -> one ES node each.
///   * derived_fd  -> one SYMBOLIC DerivedFD node each (es_root_id kept),
///                    prerequisites = {the plan's ES node id}. Expanded to
///                    concrete per-root nodes by the executor once the ES
///                    bundle converges (run(), next increment).
///   * nuclear_fd  -> all-atoms sentinel (pert.atom < 0) expands to one node
///                    per atom in [0, n_atoms); a concrete atom passes through.
std::vector<CalcNode> build_dag(const ResponsePlan &plan, int n_atoms) {
  std::vector<CalcNode> dag;

  // Pre-index static dipole frequencies per perturbation for seed hints.
  auto has_static_fd = [&](const std::string &pdesc) {
    for (const auto &r : plan.fd)
      if (r.pert.description() == pdesc && std::abs(r.freq) < 1e-12)
        return true;
    return false;
  };

  for (const auto &r : plan.fd) {
    CalcNode n;
    n.kind      = CalcKind::FD;
    n.pert      = r.pert;
    n.freq      = r.freq;
    n.protocols = r.protocols;
    n.id        = fd_node_id(r.pert, r.freq);
    if (std::abs(r.freq) >= 1e-12 && has_static_fd(r.pert.description()))
      n.seed_from = fd_node_id(r.pert, 0.0);
    dag.push_back(std::move(n));
  }

  std::string first_es_id;
  for (const auto &r : plan.es) {
    CalcNode n;
    n.kind      = CalcKind::ES;
    n.tda       = r.tda;
    n.n_roots   = r.n_roots;
    n.protocols = r.protocols;
    n.id        = es_node_id(r.tda, r.n_roots);
    if (first_es_id.empty()) first_es_id = n.id;
    dag.push_back(std::move(n));
  }

  for (const auto &r : plan.derived_fd) {
    CalcNode n;
    n.kind       = CalcKind::DerivedFD;
    n.pert       = r.pert;
    n.es_root_id     = r.es_root_id;
    n.es_freq_factor = r.es_freq_factor;
    n.protocols      = r.protocols;
    n.id             = derived_fd_node_id(r.pert, r.es_root_id);
    if (!first_es_id.empty()) n.prerequisites.push_back(first_es_id);
    dag.push_back(std::move(n));
  }

  for (const auto &r : plan.nuclear_fd) {
    if (r.pert.is_all_atoms()) {
      for (int a = 0; a < n_atoms; ++a) {
        CalcNode n;
        n.kind      = CalcKind::NuclearFD;
        n.pert      = Perturbation::nuclear(a, r.pert.axis);
        n.freq      = r.freq;
        n.protocols = r.protocols;
        n.id        = fd_node_id(n.pert, r.freq);
        dag.push_back(std::move(n));
      }
    } else {
      CalcNode n;
      n.kind      = CalcKind::NuclearFD;
      n.pert      = r.pert;
      n.freq      = r.freq;
      n.protocols = r.protocols;
      n.id        = fd_node_id(r.pert, r.freq);
      dag.push_back(std::move(n));
    }
  }

  for (const auto &r : plan.vbc) {
    CalcNode n;
    n.kind          = CalcKind::VBC;
    n.pert          = r.pert_b;     // B side
    n.freq          = r.freq_b;
    n.pert_c        = r.pert_c;     // C side
    n.freq_c        = r.freq_c;
    n.protocols     = r.protocols;
    n.id            = vbc_node_id(r.pert_b, r.pert_c, r.freq_b, r.freq_c);
    // Hard prerequisites: the two converged first-order states (same protocol).
    n.prerequisites = {fd_node_id(r.pert_b, r.freq_b),
                       fd_node_id(r.pert_c, r.freq_c)};
    dag.push_back(std::move(n));
  }

  return dag;
}

// ---------------------------------------------------------------------------
// reconcile_protocol — what to do with one node at one protocol step. Pure.
// ---------------------------------------------------------------------------

/// Implements the doc-15 reconcile table at threshold `thresh`:
///   converged at this protocol step                       -> Skip
///   else converged at a coarser-or-equal protocol step    -> Restart
///   else present, not converged, not diverged    -> Resume
///   else present, diverged                        -> Fresh
///   else absent                                   -> Fresh
///
/// HONEST-CLIMB refinement (FD family, active when max_iters > 0): a rung whose
/// last attempt already used >= max_iters iterations without diverging is DONE
/// for this rung -> Skip. The ladder is a guess-refinement scheme: the coarse
/// rung's job is to produce the best seed for the next rung (which Restart-seeds
/// from it and projects up), NOT to block the climb. `converged` stays honestly
/// false — the property layer's per-row accuracy carries the verdict, and the
/// consumer (cm_record/cm_check) applies the dconv policy. Re-running with a
/// LARGER --maxiter turns the same entry back into Resume (iter < max_iters), so
/// "continue iterating" is still available. max_iters == 0 = legacy strict
/// behavior (Resume forever; run() then halts on no-progress).
NodeAction reconcile_protocol(const CalcNode &node, const nlohmann::json &meta,
                          double thresh, int max_iters = 0) {
  const int         k   = default_k_for_thresh(thresh);
  const std::string key = protocol_key(thresh, k);

  // A DIVERGED exact-key entry must NOT resume the blown-up state — but it
  // should still re-seed from the nearest COARSER converged neighbor when one
  // exists (doc-15 reconcile table), not fall back to a bare x0=0 Fresh solve
  // and throw away the coarse-to-fine ladder. Restart re-seeds; Fresh only
  // when there is no usable coarser source. (Safe against a persistently-
  // diverging channel: a deterministic re-solve from the same coarser seed
  // saves identical metadata, which run()'s progress-aware no-progress guard
  // then quarantines — no infinite Restart loop.)
  if (node.kind == CalcKind::ES) {
    if (const auto *e = detail_calc::es_entry(meta, key)) {
      if (detail_calc::entry_converged(*e)) return NodeAction::Skip;
      if (detail_calc::entry_diverged(*e))
        return detail_calc::has_coarser_converged_es(meta, thresh, k)
                   ? NodeAction::Restart : NodeAction::Fresh;
      return NodeAction::Resume;
    }
    return detail_calc::has_coarser_converged_es(meta, thresh, k)
               ? NodeAction::Restart
               : NodeAction::Fresh;
  }

  if (node.kind == CalcKind::VBC) {
    if (const auto *e = detail_calc::vbc_entry(meta, node.id, key)) {
      if (detail_calc::entry_converged(*e)) return NodeAction::Skip;
      if (detail_calc::entry_diverged(*e))
        return detail_calc::has_coarser_converged_vbc(meta, node.id, thresh, k)
                   ? NodeAction::Restart : NodeAction::Fresh;
      return NodeAction::Resume;
    }
    return detail_calc::has_coarser_converged_vbc(meta, node.id, thresh, k)
               ? NodeAction::Restart
               : NodeAction::Fresh;
  }

  // FD / NuclearFD / (concrete) DerivedFD all live in fd_states.
  const std::string pert = node.pert.description();
  const std::string fkey = ResponseMetadata::freq_key(node.freq);
  if (const auto *e = detail_calc::fd_entry(meta, pert, key, fkey)) {
    if (detail_calc::entry_converged(*e)) return NodeAction::Skip;
    if (detail_calc::entry_diverged(*e))
      // Re-seed from the nearest coarser NON-diverged source (best_usable
      // excludes the diverged exact entry itself); Fresh only if none exists.
      return best_usable_fd_source_key(meta, pert, fkey, thresh, k, key).empty()
                 ? NodeAction::Fresh : NodeAction::Restart;
    // Honest-climb: budget exhausted at this rung -> done here, climb (see
    // the doc block above). `iter` is the last attempt's count (save-side).
    if (max_iters > 0 && e->value("iter", 0) >= max_iters)
      return NodeAction::Skip;
    return NodeAction::Resume;
  }
  // No exact entry: a coarser-or-equal NON-diverged snapshot (converged OR a
  // partial) is a usable restart seed — a coarse partial need not be converged
  // to seed the finer step. Same helper try_load_fd_state uses, so the verdict
  // and the load can never disagree (doc-15 reconcile refinement #5).
  return !best_usable_fd_source_key(meta, pert, fkey, thresh, k, key).empty()
             ? NodeAction::Restart
             : NodeAction::Fresh;
}

// ---------------------------------------------------------------------------
// prerequisites_converged — readiness gate (implicit-by-identity). Pure.
// ---------------------------------------------------------------------------

/// True iff every hard_dep of `node` is converged at protocol step `thresh`. Dep ids
/// are resolved against `dag`; each dep's convergence is checked at the SAME
/// protocol step key (the doc-15 "reconstruct prerequisites' identity at the same
/// protocol_key" contract). A dep id that resolves to no node, or to a node
/// with no converged entry, gates the node out.
bool prerequisites_converged(const CalcNode &node, const std::vector<CalcNode> &dag,
                const nlohmann::json &meta, double thresh) {
  if (node.prerequisites.empty()) return true;
  const int         k   = default_k_for_thresh(thresh);
  const std::string key = protocol_key(thresh, k);

  auto find_node = [&](const std::string &id) -> const CalcNode * {
    for (const auto &n : dag)
      if (n.id == id) return &n;
    return nullptr;
  };

  for (const auto &dep_id : node.prerequisites) {
    const CalcNode *dep = find_node(dep_id);
    if (!dep) return false;
    bool ok = false;
    if (dep->kind == CalcKind::ES) {
      const auto *e = detail_calc::es_entry(meta, key);
      ok = e && detail_calc::entry_converged(*e);
    } else {
      const auto *e = detail_calc::fd_entry(meta, dep->pert.description(),
                                            key,
                                            ResponseMetadata::freq_key(dep->freq));
      ok = e && detail_calc::entry_converged(*e);
    }
    if (!ok) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// schedule — two-phase, perturbation-parallel waves of WorkItems. Pure.
// ---------------------------------------------------------------------------

/// Emit work in waves (a wave = items with no edges between them, safe to run
/// concurrently). Channels (perturbation / ES bundle) are the unit of
/// parallelism and share every wave. Ordering:
///   Phase A (protocol step = ramp.front()):
///     wave A1 — each perturbation's SEED STATE (ω=0 if present, else lowest |ω|),
///     wave A2 — all remaining nodes at this protocol step.
///   Phase B (protocol step in ramp[1..]): one wave per protocol step of every participating
///     node (embarrassingly parallel).
/// A node participates in a protocol step only if the protocol step is in its protocols ladder
/// and reconcile_protocol != Skip. Symbolic nodes (DerivedFD "*") are excluded —
/// they are expansion templates resolved post-ES. Gated nodes (deps not
/// ready) are dropped from the wave and picked up on a later re-schedule
/// (after run() updates the metadata). Empty waves are omitted.
///
/// CONTRACT — only `waves.front()` has FINAL actions. Every wave is reconciled
/// against the SAME metadata snapshot, so a later wave's NodeActions are
/// provisional: e.g. a Phase-B item may read `Fresh` here because P0 has not
/// run yet, when it should be `Restart` once P0 lands. Membership is correct;
/// actions downstream of the front wave are not. This is safe because
/// `CalcManager::run` consumes ONLY the front wave and then re-schedules
/// against fresh metadata — so each wave's actions are final by the time it
/// becomes the front. Do not consume non-front waves' actions directly. (The
/// full list is returned anyway: it is cheap and makes the schedule legible /
/// unit-testable.)
std::vector<std::vector<WorkItem>>
schedule(const std::vector<CalcNode> &dag,
         const std::vector<double> &global_ramp,
         const nlohmann::json &meta, int max_iters = 0) {
  using detail_calc::ramp_contains;
  std::vector<std::vector<WorkItem>> waves;

  const std::vector<double> ramp = detail_calc::sorted_ramp(global_ramp);
  if (ramp.empty()) return waves;

  // Schedulable (non-symbolic) nodes only.
  std::vector<const CalcNode *> nodes;
  for (const auto &n : dag)
    if (!n.is_symbolic()) nodes.push_back(&n);

  // Does node participate at this protocol step (in ladder, not Skip, deps ready)?
  auto make_item = [&](const CalcNode *n, double thresh,
                       WorkItem &out) -> bool {
    if (!ramp_contains(n->protocols, thresh)) return false;
    if (!prerequisites_converged(*n, dag, meta, thresh))   return false;
    const NodeAction a = reconcile_protocol(*n, meta, thresh, max_iters);
    if (a == NodeAction::Skip) return false;
    out = WorkItem{n, thresh, a};
    return true;
  };

  // ---- Phase A: seed states then the rest, at ramp.front() --------------------
  const double p0 = ramp.front();

  // Per-perturbation seed state = the node with freq closest to 0 (static
  // preferred), among nodes that participate at p0. Precomputed once (O(n));
  // std::map keeps A1 in a deterministic (rank-consistent) perturbation order.
  std::map<std::string, const CalcNode *> seed_by_perturbation;
  for (const CalcNode *n : nodes) {
    if (!ramp_contains(n->protocols, p0)) continue;
    const std::string ch = perturbation_key(*n);
    auto it = seed_by_perturbation.find(ch);
    if (it == seed_by_perturbation.end() ||
        std::abs(n->freq) < std::abs(it->second->freq))
      seed_by_perturbation[ch] = n;
  }

  std::vector<WorkItem> waveA1, waveA2;
  for (const auto &[ch, seed] : seed_by_perturbation) {
    WorkItem it;
    if (make_item(seed, p0, it)) waveA1.push_back(it);
  }
  // A2: every other p0 node that isn't its perturbation's seed state.
  for (const CalcNode *n : nodes) {
    auto a = seed_by_perturbation.find(perturbation_key(*n));
    if (a != seed_by_perturbation.end() && a->second == n) continue;
    WorkItem it;
    if (make_item(n, p0, it)) waveA2.push_back(it);
  }
  if (!waveA1.empty()) waves.push_back(std::move(waveA1));
  if (!waveA2.empty()) waves.push_back(std::move(waveA2));

  // ---- Phase B: one wave per finer protocol step -----------------------------------
  for (size_t i = 1; i < ramp.size(); ++i) {
    std::vector<WorkItem> wave;
    for (const CalcNode *n : nodes) {
      WorkItem it;
      if (make_item(n, ramp[i], it)) wave.push_back(it);
    }
    if (!wave.empty()) waves.push_back(std::move(wave));
  }

  return waves;
}

// ---------------------------------------------------------------------------
// Executor seam (World-bound impl + CalcManager::run land in the next
// increment, calc/calc_executor.hpp). Declared here so the scheduling side
// and its tests agree on the contract.
// ---------------------------------------------------------------------------

struct NodeResult {
  bool                converged = false;
  std::string         reached_protocol_key;
  std::vector<double> es_root_freqs;   // vestigial: DerivedFD expansion is now
                                       // metadata-driven (expand_converged_es
                                       // reads roots from response_metadata.json).
                                       // run()'s loop discards NodeResult, so this
                                       // is no longer populated or read.
};

/// Runs one protocol step of one node: load nearest seed -> iterate at this protocol
/// -> save (+metrics). The only World-bound surface of the calc manager.
struct ICalcExecutor {
  virtual NodeResult run_protocol(const WorkItem &item) = 0;
  virtual ~ICalcExecutor() = default;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_MANAGER_HPP
