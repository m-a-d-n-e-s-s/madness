// ===========================================================================
// CalcManager scheduling-core tests (doc 15, 15a). Pure C++, no MPI / no MRA
// runtime — exercises build_dag, reconcile_protocol, prerequisites_converged, and schedule
// against hand-written response_metadata.json blobs.
// ===========================================================================

#include "../calc/calc_manager.hpp"
#include "../ResponsePropertyPlanner.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <string>

using namespace molresponse_v3;
using nlohmann::json;

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

const CalcNode *find_id(const std::vector<CalcNode> &dag, const std::string &id) {
  for (const auto &n : dag)
    if (n.id == id) return &n;
  return nullptr;
}

bool wave_has(const std::vector<WorkItem> &w, const std::string &id) {
  for (const auto &it : w)
    if (it.node->id == id) return true;
  return false;
}

[[maybe_unused]] NodeAction action_in(const std::vector<WorkItem> &w, const std::string &id) {
  for (const auto &it : w)
    if (it.node->id == id) return it.action;
  return NodeAction::Skip; // sentinel "not present"
}

// Minimal empty aggregate.
json empty_meta() {
  return json{{"schema_version", 1},
              {"protocols", json::object()},
              {"fd_states", json::object()},
              {"excited_states", json::object()},
              {"properties", json::object()}};
}

// Register a protocol key in the registry.
void put_protocol(json &m, double thresh) {
  const int k = default_k_for_thresh(thresh);
  m["protocols"][protocol_key(thresh, k)] = {{"thresh", thresh}, {"k", k},
                                             {"index", -1}};
}

// Write an fd_states entry. `iter` = the last attempt's iteration count
// (save-side state.iter) — drives reconcile's honest-climb rule.
void put_fd(json &m, const std::string &pert, double thresh, double freq,
            bool converged, bool diverged = false, int iter = 0,
            bool stalled = false) {
  const std::string key = protocol_key_at(thresh);
  const std::string fk  = ResponseMetadata::freq_key(freq);
  put_protocol(m, thresh);
  m["fd_states"][pert][key][fk] = {{"freq", freq},
                                   {"converged", converged},
                                   {"diverged", diverged},
                                   {"iter", iter},
                                   {"stalled", stalled}};
}

void put_es(json &m, double thresh, bool converged, bool diverged = false) {
  const std::string key = protocol_key_at(thresh);
  put_protocol(m, thresh);
  m["excited_states"][key] = {{"converged", converged},
                              {"diverged", diverged}};
}

// ES bundle entry with the identity + progress fields save_es_roots writes
// (type "tda"/"full", shell, n_roots, iter, stalled). put_es above writes the
// legacy shape (no identity fields).
void put_es_bundle(json &m, double thresh, bool converged, const std::string &type,
                   int n_roots, bool stalled = false, int iter = 0,
                   bool diverged = false) {
  const std::string key = protocol_key_at(thresh);
  put_protocol(m, thresh);
  m["excited_states"][key] = {{"type", type},         {"shell", "closed"},
                              {"n_roots", n_roots},   {"bundle_dir", "es__" + key},
                              {"converged", converged}, {"diverged", diverged},
                              {"stalled", stalled},   {"iter", iter}};
}

void put_vbc(json &m, const std::string &id, double thresh, bool converged,
             bool diverged = false) {
  const std::string key = protocol_key_at(thresh);
  put_protocol(m, thresh);
  m["vbc_states"][id][key] = {{"converged", converged}, {"diverged", diverged}};
}

} // namespace

int main() {
  const std::vector<double> P = {1e-4, 1e-6};

  // ====== build_dag: polarizability on x,y,z ===============================
  std::printf("=== build_dag: alpha(0.057) x/y/z ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Polarizability;
    r.frequencies = {0.057};
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), /*n_atoms=*/3);
    EXPECT(dag.size() == 3, "3 FD nodes");
    EXPECT(find_id(dag, "fd:dipole_x@f0.05700") != nullptr, "dipole_x node id");
    const CalcNode *zx = find_id(dag, "fd:dipole_z@f0.05700");
    EXPECT(zx && zx->kind == CalcKind::FD && std::abs(zx->freq - 0.057) < 1e-9,
           "dipole_z node fields");
    EXPECT(zx && zx->seed_from.empty(),
           "no static seed hint when no ω=0 in plan");
  }

  // ====== build_dag: SHG seeds dynamics from static of same perturbation ========
  std::printf("=== build_dag: static-seed hint ===\n");
  {
    // Two dipole_x FD: ω=0 and ω=0.1 -> dynamic should hint at the static.
    ResponsePlan plan;
    plan.fd.push_back({Perturbation::dipole(0), 0.0, P});
    plan.fd.push_back({Perturbation::dipole(0), 0.1, P});
    auto dag = build_dag(plan, 0);
    const CalcNode *dyn = find_id(dag, "fd:dipole_x@f0.10000");
    EXPECT(dyn && dyn->seed_from == "fd:dipole_x@f0.00000",
           "dynamic seeds from same-perturbation static");
    const CalcNode *stat = find_id(dag, "fd:dipole_x@f0.00000");
    EXPECT(stat && stat->seed_from.empty(), "static has no seed hint");
  }

  // ====== build_dag: nuclear_all expands to 3N ============================
  std::printf("=== build_dag: nuclear_all expansion ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Nuclear;
    r.frequencies = {0.0};
    r.protocol_thresholds = P;
    auto plan = plan_one(r);
    auto dag = build_dag(plan, /*n_atoms=*/2);  // 2 atoms × 3 axes = 6 nuc nodes
    int nuc = 0;
    for (const auto &n : dag)
      if (n.kind == CalcKind::NuclearFD) ++nuc;
    EXPECT(nuc == 6, "2 atoms × 3 axes = 6 NuclearFD nodes");
    EXPECT(find_id(dag, "fd:nuc_0_x@f0.00000") != nullptr, "nuc atom0 x");
    EXPECT(find_id(dag, "fd:nuc_1_z@f0.00000") != nullptr, "nuc atom1 z");
    for (const auto &n : dag)
      EXPECT(!(n.kind == CalcKind::NuclearFD && n.pert.is_all_atoms()),
             "no unexpanded all-atoms sentinel remains");
  }

  // ====== build_dag: resonance Raman -> ES + symbolic derived FD ===========
  std::printf("=== build_dag: resonant gradient (ES + derived) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 4;
    r.tpa = true;   // this case pins the DERIVED-FD DAG mechanics, so ask for the legs
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);
    const CalcNode *es = find_id(dag, "es:tda_n4");
    EXPECT(es && es->kind == CalcKind::ES && es->n_roots == 4, "ES node");
    const CalcNode *d = find_id(dag, "dfd:dipole_x@*");
    EXPECT(d && d->is_symbolic(), "derived FD is symbolic");
    EXPECT(d && d->prerequisites.size() == 1 && d->prerequisites[0] == "es:tda_n4",
           "derived FD hard-deps the ES bundle");
  }

  // ====== reconcile_protocol: the 5-row table =================================
  std::printf("=== reconcile_protocol table ===\n");
  {
    CalcNode fd;
    fd.kind = CalcKind::FD;
    fd.pert = Perturbation::dipole(0);
    fd.freq = 0.057;
    fd.protocols = P;
    fd.id = fd_node_id(fd.pert, fd.freq);

    EXPECT(reconcile_protocol(fd, empty_meta(), 1e-6) == NodeAction::Fresh,
           "absent -> Fresh");

    json m1 = empty_meta();
    put_fd(m1, "dipole_x", 1e-6, 0.057, /*converged=*/true);
    EXPECT(reconcile_protocol(fd, m1, 1e-6) == NodeAction::Skip,
           "converged at protocol step -> Skip");

    json m2 = empty_meta();
    put_fd(m2, "dipole_x", 1e-4, 0.057, /*converged=*/true);  // coarser only
    EXPECT(reconcile_protocol(fd, m2, 1e-6) == NodeAction::Restart,
           "converged at coarser protocol step -> Restart");

    json m3 = empty_meta();
    put_fd(m3, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/false);
    EXPECT(reconcile_protocol(fd, m3, 1e-6) == NodeAction::Resume,
           "present, not converged, not diverged -> Resume");

    json m4 = empty_meta();
    put_fd(m4, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/true);
    EXPECT(reconcile_protocol(fd, m4, 1e-6) == NodeAction::Fresh,
           "present, diverged -> Fresh");

    // refinement #5: a coarser-or-equal PARTIAL (not converged, not diverged) is
    // a usable restart seed — must be Restart, not Fresh (else the good partial
    // is discarded and the finer step re-seeds from zero).
    json m5 = empty_meta();
    put_fd(m5, "dipole_x", 1e-4, 0.057, /*converged=*/false, /*diverged=*/false);
    EXPECT(reconcile_protocol(fd, m5, 1e-6) == NodeAction::Restart,
           "coarser PARTIAL (not converged, not diverged) -> Restart");

    // but a coarser DIVERGED-only snapshot is never a seed -> Fresh.
    json m6 = empty_meta();
    put_fd(m6, "dipole_x", 1e-4, 0.057, /*converged=*/false, /*diverged=*/true);
    EXPECT(reconcile_protocol(fd, m6, 1e-6) == NodeAction::Fresh,
           "coarser DIVERGED only -> Fresh (never seed a blown-up state)");

    // ---- honest-climb: a budget-exhausted rung is DONE, the ladder climbs ----
    // (`converged` stays false; only the reconcile verdict flips to Skip so the
    // next rung's Restart takes over. max_iters == 0 keeps legacy semantics.)
    json m7 = empty_meta();
    put_fd(m7, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/false,
           /*iter=*/8);
    EXPECT(reconcile_protocol(fd, m7, 1e-6, /*max_iters=*/8) == NodeAction::Skip,
           "budget exhausted (iter >= max_iters), not diverged -> Skip (climb)");
    EXPECT(reconcile_protocol(fd, m7, 1e-6, /*max_iters=*/25) == NodeAction::Resume,
           "budget remains (iter < max_iters) -> Resume (continue iterating)");
    EXPECT(reconcile_protocol(fd, m7, 1e-6) == NodeAction::Resume,
           "max_iters unspecified (0) -> legacy Resume");
    // ---- a STALLED rung is exhausted too, however few iterations it used ----
    // The plateau detector gives up at stall_window (~6) iterations, far short
    // of the shipped decks' max_iters (60), so `iter >= max_iters` never fired
    // for it: reconcile said Resume, the re-solve plateaued in the same place,
    // and reconcile said Resume again — forever. Re-solving a plateaued state
    // is deterministic, so there is nothing to gain from the next attempt.
    json m7s = empty_meta();
    put_fd(m7s, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/false,
           /*iter=*/6, /*stalled=*/true);
    EXPECT(reconcile_protocol(fd, m7s, 1e-6, /*max_iters=*/60) == NodeAction::Skip,
           "stalled well under the budget -> Skip (climb), not Resume-forever");
    // Not stalled at the same iteration count keeps the old answer.
    json m7n = empty_meta();
    put_fd(m7n, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/false,
           /*iter=*/6, /*stalled=*/false);
    EXPECT(reconcile_protocol(fd, m7n, 1e-6, /*max_iters=*/60) == NodeAction::Resume,
           "same iter, not stalled -> Resume (still converging)");

    // A diverged state never climbs — divergence beats the exhausted-budget rule.
    json m8 = empty_meta();
    put_fd(m8, "dipole_x", 1e-6, 0.057, /*converged=*/false, /*diverged=*/true,
           /*iter=*/8);
    EXPECT(reconcile_protocol(fd, m8, 1e-6, /*max_iters=*/8) == NodeAction::Fresh,
           "diverged at budget -> Fresh (divergence beats honest-climb)");
  }

  // ====== reconcile_protocol: ES path =========================================
  std::printf("=== reconcile_protocol ES ===\n");
  {
    CalcNode es;
    es.kind = CalcKind::ES; es.tda = true; es.n_roots = 3; es.protocols = P;
    es.id = es_node_id(true, 3);
    json m = empty_meta();
    put_es(m, 1e-6, /*converged=*/true);
    EXPECT(reconcile_protocol(es, m, 1e-6) == NodeAction::Skip, "ES converged -> Skip");
    json m2 = empty_meta();
    put_es(m2, 1e-4, /*converged=*/true);
    EXPECT(reconcile_protocol(es, m2, 1e-6) == NodeAction::Restart,
           "ES coarser converged -> Restart");
  }

  // ====== prerequisites_converged: VBC-like synthetic dep ==============================
  std::printf("=== prerequisites_converged (synthetic VBC dep) ===\n");
  {
    std::vector<CalcNode> dag;
    CalcNode b; b.kind = CalcKind::FD; b.pert = Perturbation::dipole(0);
    b.freq = 0.05; b.protocols = P; b.id = fd_node_id(b.pert, b.freq);
    CalcNode c; c.kind = CalcKind::FD; c.pert = Perturbation::dipole(1);
    c.freq = 0.07; c.protocols = P; c.id = fd_node_id(c.pert, c.freq);
    CalcNode vbc; vbc.kind = CalcKind::VBC; vbc.id = "vbc:test";
    vbc.protocols = P; vbc.prerequisites = {b.id, c.id};
    dag = {b, c, vbc};

    json m = empty_meta();
    EXPECT(!prerequisites_converged(vbc, dag, m, 1e-6), "deps absent -> not ready");
    put_fd(m, "dipole_x", 1e-6, 0.05, true);
    EXPECT(!prerequisites_converged(vbc, dag, m, 1e-6), "one dep ready -> still not ready");
    put_fd(m, "dipole_y", 1e-6, 0.07, true);
    EXPECT(prerequisites_converged(vbc, dag, m, 1e-6), "both deps converged -> ready");
    // Same deps but at a coarser protocol step must NOT satisfy a finer protocol step.
    json m2 = empty_meta();
    put_fd(m2, "dipole_x", 1e-4, 0.05, true);
    put_fd(m2, "dipole_y", 1e-4, 0.07, true);
    EXPECT(!prerequisites_converged(vbc, dag, m2, 1e-6),
           "deps converged only at coarser protocol step -> not ready at finer");
  }

  // ====== schedule: two-phase, perturbation-parallel ===========================
  std::printf("=== schedule: alpha x/y/z over [1e-4,1e-6] ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Polarizability;
    r.frequencies = {0.0, 0.1};            // static + dynamic per axis
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);  // 6 FD nodes (3 axes × 2 freqs)
    auto waves = schedule(dag, P, empty_meta());

    // Expect: A1 (3 seed states = the ω=0 of each axis), A2 (3 dynamics @1e-4),
    //         then B (all 6 @1e-6).
    EXPECT(waves.size() == 3, "3 waves: A1, A2, B");
    if (waves.size() == 3) {
      EXPECT(waves[0].size() == 3, "A1 = 3 perturbation seed states");
      EXPECT(wave_has(waves[0], "fd:dipole_x@f0.00000") &&
             wave_has(waves[0], "fd:dipole_y@f0.00000") &&
             wave_has(waves[0], "fd:dipole_z@f0.00000"),
             "A1 holds the static (seed state) of each perturbation");
      EXPECT(waves[1].size() == 3, "A2 = 3 dynamics @ P0");
      EXPECT(wave_has(waves[1], "fd:dipole_x@f0.10000"),
             "A2 holds the dynamic frequencies");
      EXPECT(waves[2].size() == 6, "B = all 6 nodes @ 1e-6");
      for (const auto &it : waves[0])
        EXPECT(std::abs(it.thresh - 1e-4) < 1e-12 &&
               it.action == NodeAction::Fresh,
               "A1 items: thresh=1e-4, Fresh");
    }
  }

  // ====== schedule: skips converged, symbolic excluded ====================
  std::printf("=== schedule: skip converged + exclude symbolic ===\n");
  {
    // Resonant + tpa: ES bundle + symbolic derived. The case tests that the
    // scheduler EXCLUDES symbolic nodes, so the plan must contain one —
    // without tpa the roots-only plan has none and the loop passes trivially.
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::PolarizabilityGradient;
    r.gradient_mode = GradientMode::Resonant;
    r.n_roots = 2;
    r.tpa = true;
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);
    auto waves = schedule(dag, P, empty_meta());
    for (const auto &w : waves)
      EXPECT(!wave_has(w, "dfd:dipole_x@*"),
             "symbolic derived FD never scheduled");
    // Only the ES node is schedulable -> at least one wave with es:tda_n2.
    bool es_seen = false;
    for (const auto &w : waves) es_seen |= wave_has(w, "es:tda_n2");
    EXPECT(es_seen, "ES node is scheduled");

    // Now mark the ES bundle converged at both protocol steps -> everything Skips.
    json m = empty_meta();
    put_es(m, 1e-4, true);
    put_es(m, 1e-6, true);
    auto waves2 = schedule(dag, P, m);
    EXPECT(waves2.empty(), "fully-converged plan schedules no work");
  }

  // ====== build_dag: Hyperpolarizability SHG -> FD + VBC nodes =============
  std::printf("=== build_dag: beta SHG (VBC nodes + prerequisites) ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Hyperpolarizability;
    r.beta_process = BetaProcess::SHG;
    r.frequencies = {0.05};
    r.axes = {'x', 'y'};
    r.protocol_thresholds = P;
    auto dag = build_dag(plan_one(r), 0);
    // FD at omega and 2*omega (beta_freqs SHG) for x and y -> 4 FD nodes.
    EXPECT(find_id(dag, "fd:dipole_x@f0.05000") != nullptr, "FD x @ omega");
    EXPECT(find_id(dag, "fd:dipole_x@f0.10000") != nullptr, "FD x @ 2omega");
    // VBC: one per (B,C) axis pair at (omega,omega) -> 4 nodes (xx,xy,yx,yy).
    int nvbc = 0;
    for (const auto &n : dag) if (n.kind == CalcKind::VBC) ++nvbc;
    EXPECT(nvbc == 4, "4 VBC nodes (x,y x x,y)");
    const CalcNode *vbc_xy =
        find_id(dag, "vbc:dipole_x__dipole_y@f0.05000_f0.05000");
    EXPECT(vbc_xy != nullptr, "VBC(x,y) node id");
    EXPECT(vbc_xy && vbc_xy->prerequisites.size() == 2 &&
           vbc_xy->prerequisites[0] == "fd:dipole_x@f0.05000" &&
           vbc_xy->prerequisites[1] == "fd:dipole_y@f0.05000",
           "VBC(x,y) prerequisites = FD(x,omega), FD(y,omega)");
  }

  // ====== prerequisites_converged gate for a VBC node =====================
  std::printf("=== VBC prerequisites_converged gate ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Hyperpolarizability;
    r.beta_process = BetaProcess::SHG;
    r.frequencies = {0.05};
    r.axes = {'x'};
    r.protocol_thresholds = {1e-4};
    auto dag = build_dag(plan_one(r), 0);
    const CalcNode *vbc =
        find_id(dag, "vbc:dipole_x__dipole_x@f0.05000_f0.05000");
    EXPECT(vbc != nullptr, "VBC(x,x) present");
    json m = empty_meta();
    EXPECT(vbc && !prerequisites_converged(*vbc, dag, m, 1e-4),
           "FD prereq absent -> VBC not ready");
    put_fd(m, "dipole_x", 1e-4, 0.05, /*converged=*/true);
    EXPECT(vbc && prerequisites_converged(*vbc, dag, m, 1e-4),
           "FD(x,0.05) converged -> VBC ready");
  }

  // ====== reconcile_protocol VBC branch + schedule gating =================
  std::printf("=== reconcile_protocol VBC + schedule gating ===\n");
  {
    ResponsePropertyRequest r;
    r.kind = ResponsePropertyKind::Hyperpolarizability;
    r.beta_process = BetaProcess::SHG;
    r.frequencies = {0.05};
    r.axes = {'x'};
    r.protocol_thresholds = {1e-4};
    auto dag = build_dag(plan_one(r), 0);
    const CalcNode *vbc =
        find_id(dag, "vbc:dipole_x__dipole_x@f0.05000_f0.05000");
    json m = empty_meta();
    EXPECT(reconcile_protocol(*vbc, m, 1e-4) == NodeAction::Fresh,
           "VBC absent -> Fresh");
    put_vbc(m, vbc->id, 1e-4, /*converged=*/true);
    EXPECT(reconcile_protocol(*vbc, m, 1e-4) == NodeAction::Skip,
           "VBC converged -> Skip");

    // schedule: VBC excluded until its FD prereqs converge.
    json m0 = empty_meta();
    auto w0 = schedule(dag, {1e-4}, m0);
    bool vbc_scheduled = false;
    for (const auto &wave : w0)
      if (wave_has(wave, vbc->id)) vbc_scheduled = true;
    EXPECT(!vbc_scheduled, "VBC not scheduled while FD prereq unconverged");
  }

  // ====== reconcile ES: the bundle on disk must BE the requested bundle ======
  // review/findings C5: the ES verdict keyed on protocol only, so a converged TDA
  // bundle satisfied a later RPA (or larger n_roots) request in the same calc dir.
  std::printf("=== reconcile ES: bundle identity (type, n_roots) ===\n");
  {
    ResponsePlan plan;
    plan.es.push_back({/*tda=*/false, /*n_roots=*/2, {1e-4}});
    auto dag = build_dag(plan, 0);
    const CalcNode *es = find_id(dag, es_node_id(false, 2));
    json m = empty_meta();
    put_es_bundle(m, 1e-4, /*converged=*/true, "tda", 2);
    EXPECT(reconcile_protocol(*es, m, 1e-4) != NodeAction::Skip,
           "converged TDA bundle does not satisfy an RPA request");
    json m2 = empty_meta();
    put_es_bundle(m2, 1e-4, /*converged=*/true, "full", 1);
    EXPECT(reconcile_protocol(*es, m2, 1e-4) != NodeAction::Skip,
           "converged 1-root bundle does not satisfy a 2-root request");
    json m3 = empty_meta();
    put_es_bundle(m3, 1e-4, /*converged=*/true, "full", 2);
    EXPECT(reconcile_protocol(*es, m3, 1e-4) == NodeAction::Skip,
           "matching converged bundle -> Skip");
    json m4 = empty_meta();
    put_es(m4, 1e-4, /*converged=*/true);   // legacy entry: no type / n_roots
    EXPECT(reconcile_protocol(*es, m4, 1e-4) == NodeAction::Skip,
           "legacy entry without identity fields still Skips");

    // The derived-FD prerequisite must see the same identity: a converged TDA
    // bundle does not unblock legs planned for the RPA bundle.
    plan.derived_fd.push_back({Perturbation::dipole(2), "*", {1e-4}});
    auto dag2 = build_dag(plan, 0);
    const CalcNode *dfd = find_id(dag2, derived_fd_node_id(Perturbation::dipole(2), "*"));
    EXPECT(dfd && !prerequisites_converged(*dfd, dag2, m, 1e-4),
           "TDA bundle does not satisfy the RPA bundle prerequisite");
    EXPECT(dfd && prerequisites_converged(*dfd, dag2, m3, 1e-4),
           "matching bundle satisfies the prerequisite");
  }

  // ====== reconcile ES: honest climb (a stalled / exhausted rung is done) =====
  // review/findings C12: ES had no analogue of the FD honest-climb rule, so a
  // stalled bundle was re-queued as Resume on every pass with no bound.
  std::printf("=== reconcile ES: honest climb ===\n");
  {
    ResponsePlan plan;
    plan.es.push_back({/*tda=*/true, /*n_roots=*/1, {1e-4}});
    auto dag = build_dag(plan, 0);
    const CalcNode *es = find_id(dag, es_node_id(true, 1));
    json ms = empty_meta();
    put_es_bundle(ms, 1e-4, /*converged=*/false, "tda", 1, /*stalled=*/true, /*iter=*/7);
    EXPECT(reconcile_protocol(*es, ms, 1e-4, /*max_iters=*/30) == NodeAction::Skip,
           "stalled unconverged ES bundle -> Skip");
    json mx = empty_meta();
    put_es_bundle(mx, 1e-4, /*converged=*/false, "tda", 1, /*stalled=*/false, /*iter=*/30);
    EXPECT(reconcile_protocol(*es, mx, 1e-4, /*max_iters=*/30) == NodeAction::Skip,
           "budget-exhausted ES bundle -> Skip");
    json mp = empty_meta();
    put_es_bundle(mp, 1e-4, /*converged=*/false, "tda", 1, /*stalled=*/false, /*iter=*/5);
    EXPECT(reconcile_protocol(*es, mp, 1e-4, /*max_iters=*/30) == NodeAction::Resume,
           "partial ES bundle with budget left -> Resume");
    EXPECT(reconcile_protocol(*es, ms, 1e-4, /*max_iters=*/0) == NodeAction::Resume,
           "max_iters == 0 keeps the legacy Resume");
  }

  // ====== ES restart source: never a diverged bundle ==========================
  // review/findings C4: try_load_es_bundle preferred the exact-key bundle even
  // when it had diverged, so a Restart re-seeded from the blown-up state. The
  // selection is now a pure helper shared with the loader (mirrors the FD one).
  std::printf("=== best_usable_es_source_key ===\n");
  {
    const std::string k4 = protocol_key_at(1e-4), k6 = protocol_key_at(1e-6);
    json m = empty_meta();
    put_es_bundle(m, 1e-4, /*converged=*/true, "full", 2);
    put_es_bundle(m, 1e-6, /*converged=*/false, "full", 2, false, 3, /*diverged=*/true);
    EXPECT(best_usable_es_source_key(m, "full", "closed", 2, 1e-6, 8, k6) == k4,
           "diverged exact bundle is skipped; the coarser converged one is used");
    json mp = empty_meta();
    put_es_bundle(mp, 1e-4, /*converged=*/true, "full", 2);
    put_es_bundle(mp, 1e-6, /*converged=*/false, "full", 2, false, 3);
    EXPECT(best_usable_es_source_key(mp, "full", "closed", 2, 1e-6, 8, k6) == k6,
           "non-diverged exact partial is preferred");
    EXPECT(best_usable_es_source_key(m, "tda", "closed", 2, 1e-6, 8, k6).empty(),
           "type mismatch -> no source");
    EXPECT(best_usable_es_source_key(m, "full", "closed", 3, 1e-6, 8, k6).empty(),
           "n_roots mismatch -> no source");
    EXPECT(best_usable_es_source_key(empty_meta(), "full", "closed", 2, 1e-6, 8, k6).empty(),
           "empty metadata -> no source");
  }

  std::printf("\n%s  (%d failures)\n", failed ? "FAILED" : "PASSED", failed);
  return failed ? 1 : 0;
}
