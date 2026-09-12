#ifndef MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
#define MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP

// -----------------------------------------------------------------------------
// run_response — the single public entry point (doc 16, L1).
//
// One self-contained call: given a checkpoint + a response plan + settings, load
// the ground state, drive CalcManager, assemble Tier-A properties, and return a
// structured Output. It owns no scheduling, ownership, or file-format logic —
// those live in CalcManager and the persistence layer. main.cpp, the test
// runner, and (later, R3) madqc all build a ResponseWorkflowInput and call this,
// so the engine is testable, scriptable, and Python-bindable through one seam.
//
// R0a scope: wraps EXACTLY today's flow (no behavior change). The timing /
// diagnostics / exports Output slots are reserved here and filled by R1 / R2.
// -----------------------------------------------------------------------------

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"          // set_response_protocol
#include "../ResponsePropertyPlanner.hpp"   // ResponsePlan
#include "../calc/calc_executor.hpp"        // ExecutorSettings/Context, FdResponseExecutor, assemble_*
#include "../calc/calc_manager.hpp"         // CalcManager
#include "../solvers/gs_fingerprint.hpp"    // GS-archive restart-safety gate
#include "../solvers/response_metadata.hpp"

#include <nlohmann/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace molresponse_v3 {

// ---------------------------------------------------------------------------
// Public contract (doc 16). Input is World-free and self-contained; Output is
// pure JSON so madqc/tests/Python consume one object.
// ---------------------------------------------------------------------------

struct ResponseWorkflowInput {
  std::string         archive_file;   // SCF/moldft checkpoint (ground state)
  std::vector<double> protocols;      // coarse->fine truncation-threshold ladder
  // R0a: the Tier-A plan is pre-built by the caller (plan_one/merge_plans, plus
  // any es-full override). FUTURE: carry std::vector<ResponsePropertyRequest>
  // and lower to a plan inside run_response once the TDA/Full choice is a
  // request field (needed by madqc, R3).
  ResponsePlan        plan;
  ExecutorSettings    settings;       // convergence policy + ES/seed/accept knobs (World-free)
  std::optional<madness::Molecule> molecule;  // else built from the archive dir
};

struct ResponseWorkflowOutput {
  nlohmann::json properties;   // Tier-A tensors (mirror of metadata "properties")
  nlohmann::json metadata;     // full response_metadata.json
  nlohmann::json timing;       // R1 — reserved
  nlohmann::json diagnostics;  // R1 — reserved
  nlohmann::json exports;      // R2 — reserved
  nlohmann::json debug_log;    // optional iteration-level
  int            rc = 0;
};

// ---------------------------------------------------------------------------
// Small helpers factored out of the old driver (the archive-adjacent lookups).
// ---------------------------------------------------------------------------

/// Build a Molecule from moldft/mad.calc_info.json next to the archive. Returns
/// an empty Molecule if none is found (harmless for pure-alpha; nuclear/Raman
/// needs it for the per-atom expansion). Mirrors the old test-runner logic.
inline madness::Molecule molecule_from_archive_dir(const std::string &archive_file) {
  madness::Molecule molecule;
  const auto dir = std::filesystem::path(archive_file).parent_path();
  for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    const auto candidate = dir / name;
    if (!std::filesystem::exists(candidate)) continue;
    std::ifstream ifs(candidate);
    nlohmann::json j;
    ifs >> j;
    nlohmann::json mol_json;
    if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
      mol_json = j["tasks"][0]["molecule"];
    else if (j.contains("molecule"))
      mol_json = j["molecule"];
    if (!mol_json.is_null()) molecule.from_json(mol_json);
    break;
  }
  return molecule;
}

/// Resolve the moldft/mad.fock.json path next to the archive ("" if none).
inline std::string fock_json_from_archive_dir(const std::string &archive_file) {
  const auto dir = std::filesystem::path(archive_file).parent_path();
  for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
    const auto candidate = dir / name;
    if (std::filesystem::exists(candidate)) return candidate.string();
  }
  return {};
}

// ---------------------------------------------------------------------------
// Stage timing (R1a, doc 16 L3). Rank-0 wall/cpu deltas — the "shape of the
// run" (is time in solving or assembly? which stage dominates?). wall_time() /
// cpu_time() are LOCAL (not collective); stage boundaries are bracketed by the
// collectives inside each stage, so rank-0's wall is representative. Per-rank-max
// reduction is a later refinement.
// ---------------------------------------------------------------------------
namespace detail_workflow {
struct StageTimer {
  double w0 = madness::wall_time();
  double c0 = madness::cpu_time();
  nlohmann::json lap() const {
    return nlohmann::json{{"wall_s", madness::wall_time() - w0},
                          {"cpu_s", madness::cpu_time() - c0}};
  }
};
} // namespace detail_workflow

// ---------------------------------------------------------------------------
// run_response — load ground -> plan -> CalcManager -> assemble -> Output.
// ---------------------------------------------------------------------------

/// Core (stages 2–4): DAG build → CalcManager solve → Tier-A assembly → collect.
/// `gs` must already be loaded AND prepared at the coarsest protocol; `L` is the
/// cubic-cell half-edge and `fock_json` the moldft Fock path ("" = compute).
/// Fills `Output.timing` for plan_build/solve/assemble (the caller owns load +
/// total). Shared by the CLI entry (run_response) and the madqc adapter, which
/// supplies a `gs` built from an in-memory SCF rather than an archive (R3).
inline ResponseWorkflowOutput
run_response_with_ground(madness::World &world, GroundState &gs, double L,
                         const std::string &fock_json,
                         const ResponseWorkflowInput &in) {
  ResponseWorkflowOutput out;
  MADNESS_CHECK(!in.protocols.empty());
  nlohmann::json timing;

  // 2a. Build the DAG.
  detail_workflow::StageTimer t_build;
  CalcManager::Policy mgr_policy;
  mgr_policy.max_iters_per_step = in.settings.max_iters;
  mgr_policy.fd_subworlds       = in.settings.fd_subworlds;   // F2
  CalcManager mgr(in.plan, in.settings.calc_dir, mgr_policy);
  mgr.build(gs.molecule().natom());
  timing["plan_build"] = t_build.lap();

  // 2a'. GS-archive fingerprint gate (restart safety — see gs_fingerprint.hpp
  // for the failure mode: phase-flipped regenerated orbitals + cached X/Y =>
  // silently wrong properties). Runs BEFORE mgr.run so no restart state is
  // touched on a mismatch. Rank 0 hashes the archive and reads the stamp; the
  // verdict is broadcast so an abort is collective (a rank-0-only throw would
  // hang the other ranks in the solve). MADRESPONSE_ALLOW_GS_MISMATCH=1
  // downgrades the abort to a restamp — the stale states then remain loadable,
  // so it is only sane when the caller KNOWS the states match this archive.
  if (!in.archive_file.empty()) {
    const std::string meta_path =
        in.settings.calc_dir + "/response_metadata.json";
    // 0 = match/proceed, 1 = stamped fresh, 2 = warned+restamped, 3 = abort,
    // 4 = rank-0 error (unreadable archive / corrupt or newer-schema metadata /
    //     ENOSPC on stamp save). Everything inside the rank-0 block can throw,
    //     and a rank-0-only unwind BEFORE the broadcast strands the other ranks
    //     in it — so catch, encode as gate=4, and abort collectively.
    int gate = 0;
    if (world.rank() == 0) try {
      // The gate stamps the metadata BEFORE mgr.run creates anything, so the
      // calc dir may not exist yet (standalone runs point calc_dir at a fresh
      // subdir; the madqc path happens to run in an existing cwd). Create it
      // first — otherwise the stamp save() fails to open its .tmp and the
      // collective error path aborts a perfectly good run.
      std::filesystem::create_directories(in.settings.calc_dir);
      const GsFingerprint fp = gs_archive_fingerprint(in.archive_file);
      auto meta = ResponseMetadata::load_or_create(meta_path);
      const std::string stored = meta.ground_state_fingerprint();
      switch (gs_fingerprint_verdict(meta.json(), fp.hex)) {
        case GsGateVerdict::Match:    gate = 0; break;
        case GsGateVerdict::FreshDir: gate = 1; break;
        case GsGateVerdict::MissingStamp:
          gate = 2;
          madness::print(
              "[GS-FINGERPRINT] WARNING: this calc dir has response states "
              "but no ground-state stamp (metadata predates the gate). "
              "Cannot verify they belong to", in.archive_file,
              "— stamping it now; if the archive was regenerated since those "
              "states were written, start a fresh calc dir.");
          break;
        case GsGateVerdict::Mismatch: {
          const char *env = std::getenv("MADRESPONSE_ALLOW_GS_MISMATCH");
          if (env && env[0] == '1') {
            gate = 2;
            madness::print(
                "[GS-FINGERPRINT] OVERRIDE (MADRESPONSE_ALLOW_GS_MISMATCH=1): "
                "stamp", stored, "!= current", fp.hex,
                "— restamping and continuing. The existing response states "
                "remain loadable; results are only trustworthy if these "
                "archives are phase-identical.");
          } else {
            gate = 3;
            madness::print(
                "[GS-FINGERPRINT] MISMATCH in", meta_path, "\n"
                "  stored :", stored,
                "(stamped when this dir's response states were written)\n"
                "  current:", fp.hex, "(", in.archive_file, ",", fp.bytes,
                "bytes,", fp.nparts, "part(s))\n"
                "  The cached response states belong to a DIFFERENT ground-state "
                "archive. Orbitals of a regenerated ground state can be "
                "phase-flipped (identical physics, opposite signs), and reusing "
                "these states would corrupt properties silently. Point the run "
                "at the original archive, or use a fresh calc dir, or set "
                "MADRESPONSE_ALLOW_GS_MISMATCH=1 if you know better.");
          }
          break;
        }
      }
      if (gate == 1 || gate == 2) {
        meta.set_ground_state(in.archive_file, fp.hex, fp.bytes, fp.nparts);
        meta.save();
      }
    } catch (const std::exception &e) {
      gate = 4;
      madness::print("[GS-FINGERPRINT] ERROR during gate on rank 0:", e.what());
    }
    world.gop.broadcast_serializable(gate, 0);
    if (gate == 3)
      MADNESS_EXCEPTION(
          "GS-archive fingerprint mismatch — refusing to reuse restart state "
          "(details printed by rank 0; see [GS-FINGERPRINT])", 0);
    if (gate == 4)
      MADNESS_EXCEPTION(
          "GS-fingerprint gate failed on rank 0 (unreadable archive, corrupt "
          "or newer-schema metadata, or stamp-save error — see "
          "[GS-FINGERPRINT] ERROR); aborting collectively", 0);
  }

  // 2a''. F2 restart safety: sweep in any metadata shards stranded by an
  // interrupted subworld run, BEFORE reconcile reads the canonical file —
  // otherwise finished states would be invisible and silently re-solved.
  // (The stranded shards belong to the same GS the gate above just verified.)
  // Same collective-error discipline as the gate: a rank-0 throw here
  // (corrupt shard, ENOSPC on the canonical save) must not strand the other
  // ranks at the fence.
  {
    int sweep_bad = 0;
    if (world.rank() == 0) try {
      const int n =
          ResponseMetadata::merge_stale_state_shards(in.settings.calc_dir);
      if (n > 0)
        madness::print("SHARD_SWEEP  merged", n,
                       "stale metadata shard(s) from an interrupted run");
    } catch (const std::exception &e) {
      sweep_bad = 1;
      madness::print("[SHARD-SWEEP-ERROR]", e.what());
    }
    world.gop.max(sweep_bad);
    if (sweep_bad)
      MADNESS_EXCEPTION("stale-shard sweep failed on rank 0 "
                        "(see [SHARD-SWEEP-ERROR]); aborting collectively", 0);
  }
  world.gop.fence();

  // 2b. Drive the calc manager (the solve).
  detail_workflow::StageTimer t_solve;
  ExecutorContext ctx(world, gs, L, fock_json, in.settings);
  FdResponseExecutor exec(ctx);

  // F2 (doc 32 §5): the per-subworld solve closure. Only built when an archive is
  // available — GroundState's only construction path is from_archive (doc 32 §3),
  // so the madqc in-memory-GS path can't fan out (it leaves fan_out null → the
  // single-World path). Each node-subworld loads its own GS and solves its owned
  // FD items, writing to its node_index metadata shard (merged by rank 0). This
  // is exactly the F1-proven block, lifted into the live run().
  CalcManager::SubworldSolve fan_out;
  if (in.settings.fd_subworlds > 0 && in.archive_file.empty() &&
      world.rank() == 0)
    madness::print("F2: fd_subworlds =", in.settings.fd_subworlds,
                   "requested but no archive_file — staying single-World "
                   "(subworlds reload the ground state from the archive)");
  if (in.settings.fd_subworlds > 0 && !in.archive_file.empty()) {
    const std::string  archive = in.archive_file;
    const madness::Molecule mol = gs.molecule();
    const std::string  fock = fock_json;
    const double       L_   = L;
    const ExecutorSettings base = in.settings;
    fan_out = [archive, mol, fock, L_, base](
                  madness::World &sub, const std::vector<WorkItem> &items,
                  double thresh, int gid, const std::string &log_prefix) {
      // F2d-ii: verbose=false suppresses the redundant per-subworld GS/protocol
      // banners (the universe already printed them once).
      set_response_protocol(sub, L_, thresh, /*override_k=*/-1, /*verbose=*/false);
      auto gs_s = GroundState::from_archive(sub, archive, mol);
      const double t0 = madness::FunctionDefaults<3>::get_thresh();
      auto cop = madness::poperatorT(
          madness::CoulombOperatorPtr(sub, gs_s.params().lo(), 0.001 * t0));
      gs_s.prepare(sub, 0.001 * t0, cop, fock, /*verbose=*/false);
      ExecutorSettings s = base;
      s.metadata_shard = std::to_string(gid);  // F2: per-subworld metadata shard
      s.log_prefix     = log_prefix;     // F2d-i: tag this subworld's per-state lines
      s.log_group      = gid;            // F2d-i: MEMORY_HWM group= field
      ExecutorContext ctx_s(sub, gs_s, L_, fock, s);
      FdResponseExecutor exec_s(ctx_s);
      for (const auto &it : items) exec_s.run_protocol(it);
    };
  }
  nlohmann::json sched_diag = mgr.run(world, exec, fan_out);   // R1c scheduler trace
  timing["solve"] = t_solve.lap();

  // 3. Tier-A property assembly (off the solve path). A mixed request may want
  //    several properties at once (e.g. alpha + beta), so assemble each that the
  //    plan supports — not one-XOR-the-other. assemble_alpha self-guards (returns
  //    if no dipole FD; it ignores nuclear FD) so it is safe to always attempt;
  //    beta/raman come from VBC. (ES-derived scalar properties are a later step.)
  detail_workflow::StageTimer t_assemble;
  assemble_alpha(ctx, in.plan, in.protocols.back());
  if (!in.plan.vbc.empty())
    assemble_beta(ctx, in.plan, in.protocols.back());
  // 2PA: derived FD legs are only planned when the request asked for the
  // two-photon contraction (request.tpa), so their presence IS the gate.
  // assemble_tpa self-guards (returns if no ES bundle; ClosedShell only).
  if (!in.plan.derived_fd.empty())
    assemble_tpa(ctx, in.plan, in.protocols.back());
  timing["assemble"] = t_assemble.lap();

  // 4. Collect the Output from the aggregate metadata (rank 0 authoritative).
  if (world.rank() == 0) {
    auto meta = ResponseMetadata::load_or_create(
        in.settings.calc_dir + "/response_metadata.json");
    // Provenance: stamp the effective restart-I/O backend into the metadata
    // (madqc review R3 — an HDF5 run must be distinguishable from a native one).
#ifdef MADNESS_HAS_HDF5
    meta.set_io_info(hdf5_io_enabled() ? "hdf5" : "native", true);
#else
    meta.set_io_info("native", false);
#endif
    meta.save();
    out.metadata = meta.json();
    if (out.metadata.contains("properties"))
      out.properties = out.metadata["properties"];
    out.timing = std::move(timing);
    out.diagnostics = std::move(sched_diag);   // R1c scheduler trace
  }
  return out;
}

inline ResponseWorkflowOutput
run_response(madness::World &world, const ResponseWorkflowInput &in) {
  MADNESS_CHECK(!in.protocols.empty());
  const detail_workflow::StageTimer t_total;
  detail_workflow::StageTimer t_load;

  // 1. Protocol + ground state from the ARCHIVE (CLI/standalone path).
  //    set_response_protocol to the coarsest rung; the executor re-prepares per
  //    protocol during the solve.
  const auto header = GroundState::read_archive_header(world, in.archive_file);
  set_response_protocol(world, header.L, in.protocols.front());

  const madness::Molecule molecule =
      in.molecule ? *in.molecule : molecule_from_archive_dir(in.archive_file);
  auto gs = GroundState::from_archive(world, in.archive_file, molecule);
  if (world.rank() == 0) gs.print_info();

  const std::string fock_json = fock_json_from_archive_dir(in.archive_file);
  const double cur_thresh = madness::FunctionDefaults<3>::get_thresh();
  auto coulop = madness::poperatorT(
      madness::CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
  gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);
  nlohmann::json load_timing = t_load.lap();   // capture before the core

  ResponseWorkflowOutput out =
      run_response_with_ground(world, gs, header.L, fock_json, in);

  if (world.rank() == 0) {
    out.timing["load"]  = std::move(load_timing);
    out.timing["total"] = t_total.lap();
    if (in.settings.print_level >= PrintLevel::Normal)
      madness::print("[TIMING] load_wall_s=", out.timing["load"]["wall_s"],
                     "  solve_wall_s=", out.timing["solve"]["wall_s"],
                     "  assemble_wall_s=", out.timing["assemble"]["wall_s"],
                     "  total_wall_s=", out.timing["total"]["wall_s"]);
  }
  return out;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ORCHESTRATOR_RESPONSE_WORKFLOW_HPP
