#ifndef MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP

// =========================================================================
// FD save/load (Inc 13c-i: save only; load + restart precedence land in
// 13c-ii / 13c-iii).
//
// FD persistence is per-point — one Storage archive per (perturbation,
// protocol, freq) triple. The unit of persistence is therefore simpler than
// the ES bundle: a single ResponseStateX<Shell> or ResponseStateXY<Shell>
// archive plus a metadata entry in the unified response_metadata.json.
//
// Disk layout (per doc 13):
//
//   <calc dir>/
//   ├── response_metadata.json                      ← 13b aggregator
//   ├── dipole_x__1e-06_k8__f0.05700.00000          ← FD point (parallel
//   ├── dipole_x__1e-06_k8__f0.05700.00001            archive, .000N per rank)
//   └── es_bundle__1e-06_k8/                        ← ES bundle
//       ├── roots.json
//       └── root_0.00000 ...
//
// response_filename(pert, protocol_key, freq) is the join key for both the
// archive name and the metadata path (fd_states/<pert>/<key>/<fkey>) so a
// property layer can find FD and ES inputs by string compare on protocol_key.
//
// Collective discipline: state.responses[0].save() is collective on `world`.
// The metadata upsert is rank-0 only (filesystem op) and the function fences
// after both finish. Caller does NOT need to add its own fence.
// =========================================================================

#include "../Perturbations.hpp"          // Perturbation
#include "../ResponseProtocol.hpp"        // protocol_key()
#include "../kernels/tags.hpp"            // Static/Full/TDA, ClosedShell/OpenShell
#include "fd_solver.hpp"
#include "response_metadata.hpp"
#include "response_state.hpp"
#include "state_metrics.hpp"              // measure_state (per-state mem/iters)

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

namespace detail_fd_save_load {

template <typename Type>
inline const char *type_tag() {
  if constexpr (std::is_same_v<Type, Static>) return "static";
  else if constexpr (std::is_same_v<Type, Full>) return "full";
  else if constexpr (std::is_same_v<Type, TDA>)  return "tda";
  else return "unknown";
}

template <typename Shell>
inline const char *shell_tag() {
  if constexpr (std::is_same_v<Shell, ClosedShell>) return "closed_shell";
  else if constexpr (std::is_same_v<Shell, OpenShell>) return "open_shell";
  else return "unknown";
}

/// np-guard for FD/VBC archive loads (the analog of load_es_roots' guard —
/// same nio=1 primitive, same documented heap-corruption reproduction on a
/// np mismatch). Only np-locked backends are guarded: the HDF5 blob is
/// gathered to one client at write time, so it reloads at any np.
/// writer_nproc == 0 = legacy entry with no recorded count — proceed (same-np
/// is the common case; the caller's rank-0 print surfaces the assumption).
/// All inputs are rank-uniform (broadcast by the caller), so the throw is
/// collective and cannot deadlock.
inline void check_writer_nproc(madness::World &world, IoBackend backend,
                               int writer_nproc, const char *who,
                               const std::string &archive) {
  if (io_backend_np_portable(backend)) return;
  if (writer_nproc == 0 || writer_nproc == static_cast<int>(world.size()))
    return;  // no recorded count (legacy) or an exact match — nothing to guard.

  // Cross-process-count load of a native nio=1 parallel archive. This is
  // np-PORTABLE: the archive stores nio in-file and ParallelInputArchive
  // redistributes the tree across the current world on read. Verified
  // bit-identical np=2 -> np=1 on the FD/ES bundles (max reldiff 0). The former
  // hard block predated the atomic-write + collective-throw I/O fixes that
  // resolved the actual crashes (truncated archives, mismatched collectives);
  // with those in, cross-np native restart is safe. Proceed by default so runs
  // are resource-count agnostic. MADRESPONSE_STRICT_NP=1 restores the hard error
  // (e.g. as a backstop if a future nio>1 writer is introduced).
  const char *strict = std::getenv("MADRESPONSE_STRICT_NP");
  if (strict && strict[0] == '1') {
    throw std::runtime_error(
        std::string(who) + ": archive " + archive + " was written with " +
        std::to_string(writer_nproc) + " process(es) but is being loaded with " +
        std::to_string(world.size()) +
        " — MADRESPONSE_STRICT_NP=1 forbids cross-process-count native restart. "
        "Re-run with " + std::to_string(writer_nproc) +
        " rank(s), unset MADRESPONSE_STRICT_NP, or use the HDF5 backend.");
  }
  if (world.rank() == 0)
    std::fprintf(stderr,
        "[NP] %s: %s written with %d rank(s), loading with %d — native nio=1 "
        "archive is np-portable (redistributed on read); proceeding.\n",
        who, archive.c_str(), writer_nproc, static_cast<int>(world.size()));
}

} // namespace detail_fd_save_load

/// Canonical archive name for an FD point. Per doc 13:
///   <pert>__<protocol_key>__<freq_key>      dipole_x__1e-06_k8__f0.05700
inline std::string response_filename(const std::string &pert,
                                     const std::string &protocol_key_str,
                                     double freq) {
  return pert + "__" + protocol_key_str + "__"
       + ResponseMetadata::freq_key(freq);
}

/// Save one FD point. Writes the Storage archive (collective) and upserts
/// the entry into response_metadata.json (rank-0). Assumes the active
/// FunctionDefaults<3> (k, thresh) is the protocol the state was solved at —
/// the caller's protocol step has already set this.
///
/// Invariant: state.responses.size() == 1. The skeleton (and present v3 FD
/// driver) solve one perturbation/freq per State; multi-channel save is a
/// later concern.
template <typename Type, typename Shell>
void save_fd_state(madness::World &world,
                   const typename FDSolver<Type, Shell>::State &state,
                   const std::string &dir,
                   const Perturbation &pert,
                   double freq,
                   bool converged,
                   const std::string &seed = std::string(),
                   bool accepted = false,
                   double wall_s = 0.0,     // R1b: point-solve wall time
                   const std::string &metadata_shard = std::string(),
                   const std::string &log_prefix = std::string(),  // F2d
                   int log_group = -1) {                            // F2d
  // F2 (doc 32 §5.3): when metadata_shard is set the metadata upsert goes to a
  // per-group shard (response_metadata.group<tag>.json) so concurrent node-
  // subworlds never race the canonical file; rank 0 merges them after the fence.
  // "" = write the canonical file (single-World path). Per-state ARCHIVES are
  // distinct files keyed by pert/freq and always land in the shared dir.
  MADNESS_CHECK(state.responses.size() == 1);

  if (world.rank() == 0) {
    std::filesystem::create_directories(dir);
  }
  world.gop.fence();

  const double      thresh = madness::FunctionDefaults<3>::get_thresh();
  const int         k_now  = madness::FunctionDefaults<3>::get_k();
  const std::string key    = protocol_key(thresh, k_now);
  const std::string fkey   = ResponseMetadata::freq_key(freq);
  const std::string pdesc  = pert.description();
  const std::string archive_basename = response_filename(pdesc, key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  // (1) Collective binary save — same per-state primitive ES uses per-root.
  // The returned backend is recorded in the metadata entry so load opens the
  // SAME physical file family (a stale twin from a backend toggle can never
  // shadow this state). save() also removed the other backend's twin.
  const IoBackend backend = state.responses[0].save(world, archive_path);

  // (1b) Collective per-state metrics (coeffs/bytes/rss/iters/wall) — every rank.
  StateMetrics metrics =
      measure_state(world, state.responses[0], state.iter);
  metrics.wall_s = wall_s;   // R1b (measure_state is collective; this is local)

  // (2) Rank-0 metadata upsert.
  if (world.rank() == 0) {
    const std::string meta_path =
        dir + "/response_metadata" +
        (metadata_shard.empty() ? "" : ".group" + metadata_shard) + ".json";
    auto meta = ResponseMetadata::load_or_create(meta_path);

    // Register the protocol if this is the first artifact at this key.
    // `index` is unknown to this writer (it's a property of the caller's
    // ramp); a later orchestrator can overwrite it with the real ordering.
    if (!meta.json()["protocols"].contains(key)) {
      meta.set_protocol(key, thresh, k_now, /*index=*/-1);
    }

    const double bsh_res =
        state.last_bsh_residual.empty()
            ? 0.0
            : *std::max_element(state.last_bsh_residual.begin(),
                                state.last_bsh_residual.end());

    nlohmann::json entry = {
        {"freq",         freq},
        {"type",         detail_fd_save_load::type_tag<Type>()},
        {"shell",        detail_fd_save_load::shell_tag<Shell>()},
        {"converged",    converged},
        // `accepted` = converged was forced true by best-effort acceptance at
        // maxiter (--accept-at-maxiter), NOT a real target hit. The gates
        // (reconcile/prerequisites) read `converged`; this field preserves the
        // honest verdict — inspect bsh_residual to judge the deliverable quality.
        {"accepted",     accepted},
        {"diverged",     state.diverged},
        {"iter",         state.iter},
        {"bsh_residual", bsh_res},
        {"seed",         seed},   // initial-guess origin: source/fd_restart/es_root
        {"archive",      archive_basename},
        // Which physical file family holds the archive ("native"/"hdf5") —
        // loaders open THIS backend; bare .h5-existence detect is legacy-only.
        {"backend",      io_backend_tag(backend)},
        // #processes that wrote it. Native nio=1 parallel archives are
        // np-locked (worlddc assumes #writers == #readers; a mismatch silently
        // corrupts the container) — loaders refuse a native np-mismatch.
        {"writer_nproc", static_cast<int>(world.size())},
        {"metrics",      metrics.to_json()},
    };
    meta.set_fd_state(pdesc, key, fkey, entry);
    meta.save();

    // F2d: prepend the per-subworld tag (empty ⇒ unchanged, G=0 byte-identical).
    auto save_line = [&](auto &&...a) {
      if (log_prefix.empty()) madness::print(a...);
      else                    madness::print(log_prefix, a...);
    };
    save_line("[SAVE] fd_state: pert=", pdesc,
              "  protocol_key=", key,
              "  freq=", freq,
              "  archive=", archive_basename,
              "  bsh_res=", bsh_res,
              "  converged=", converged,
              (accepted ? "  (ACCEPTED best-effort @ maxiter)" : ""));
    // R1b: machine-readable per-state memory high-water mark (worst task, via
    // gop.max in measure_state) at this protocol boundary — feeds the R4
    // memory-scaling model / pre-flight abort (L2). Greppable: ^MEMORY_HWM.
    // F2d: `group=<gid>` field (only in subworld mode) lets the studies/perf-model
    // parsers attribute the line to its subworld instead of treating G copies as one.
    // Two-branch (not trailing optional args) so G=0 stays byte-identical — a
    // trailing empty print() arg would emit a stray separator space.
    if (log_group >= 0)
      madness::print("MEMORY_HWM  kind=fd  protocol=", key,
                     "  rss_gb_max=", metrics.rss_gb, "  coeffs=", metrics.coeffs,
                     "  wall_s=", wall_s, "  id=", pdesc, "@", fkey,
                     "  group=", log_group);
    else
      madness::print("MEMORY_HWM  kind=fd  protocol=", key,
                     "  rss_gb_max=", metrics.rss_gb, "  coeffs=", metrics.coeffs,
                     "  wall_s=", wall_s, "  id=", pdesc, "@", fkey);
  }
  world.gop.fence();
}

/// Exact-match load (Inc 13c-ii). Looks up the archive at
///   dir/response_filename(pert.description(), protocol_key(), freq)
/// using the ACTIVE FunctionDefaults<3> (k, thresh) — caller is expected to
/// have called set_response_protocol() first so the key matches.
///
/// Returns a fresh State with one response slot: responses[0] = loaded
/// Storage; iter=0; diverged=false; last_bsh_residual seeded from the
/// metadata entry's bsh_residual (for diagnostic continuity — the solver
/// recomputes on the first step). Density-residual history is not
/// persisted; the first-iter Δρ guard inside FDSolver::step() handles a
/// missing rho_alpha_prev (same path as ESSolver after es_load).
///
/// Throws if either the archive or the metadata entry is missing.
/// Cross-protocol fallback (load at a lower protocol than the active one)
/// is deferred to try_load_fd_state in 13c-iii.
template <typename Type, typename Shell>
typename FDSolver<Type, Shell>::State
load_fd_state(madness::World &world,
              const std::string &dir,
              const Perturbation &pert,
              double freq) {
  using State   = typename FDSolver<Type, Shell>::State;
  using Storage = typename FDSolver<Type, Shell>::Storage;

  const std::string key   = protocol_key();
  const std::string fkey  = ResponseMetadata::freq_key(freq);
  const std::string pdesc = pert.description();
  const std::string archive_basename = response_filename(pdesc, key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  // Rank-0 reads the metadata entry (required); broadcasts the small
  // numeric diagnostics. Errors are collected into `err` and re-thrown on
  // EVERY rank after the broadcast — a rank-0-only throw before a collective
  // deadlocks any caller that catches. The archive load below is collective.
  double      bsh_res        = 0.0;
  int         iter_at_save   = 0;
  int         converged_int  = 0;   // bool via int for gop.broadcast
  int         writer_nproc   = 0;   // 0 = legacy entry (no recorded count)
  std::string backend_tag;          // "" = legacy entry -> auto-detect
  std::string err;
  if (world.rank() == 0) {
    const std::string meta_path = dir + "/response_metadata.json";
    try {
      if (!std::filesystem::exists(meta_path)) {
        throw std::runtime_error("load_fd_state: missing " + meta_path);
      }
      auto meta = ResponseMetadata::load_or_create(meta_path);
      const auto &j = meta.json();
      if (!j["fd_states"].contains(pdesc) ||
          !j["fd_states"][pdesc].contains(key) ||
          !j["fd_states"][pdesc][key].contains(fkey)) {
        throw std::runtime_error(
            "load_fd_state: no fd_states/" + pdesc + "/" + key + "/" + fkey +
            " in " + meta_path);
      }
      const auto &e = j["fd_states"][pdesc][key][fkey];
      bsh_res       = e.value("bsh_residual", 0.0);
      iter_at_save  = e.value("iter", 0);
      converged_int = e.value("converged", false) ? 1 : 0;
      writer_nproc  = e.value("writer_nproc", 0);
      backend_tag   = e.value("backend", std::string{});
    } catch (const std::exception &ex) {
      err = ex.what();
    }
  }
  world.gop.broadcast_serializable(err, 0);
  if (!err.empty()) throw std::runtime_error(err);  // collective throw
  world.gop.broadcast(bsh_res,       0);
  world.gop.broadcast(iter_at_save,  0);
  world.gop.broadcast(converged_int, 0);
  world.gop.broadcast(writer_nproc,  0);
  world.gop.broadcast_serializable(backend_tag, 0);

  const IoBackend backend = io_backend_from_tag(backend_tag);
  detail_fd_save_load::check_writer_nproc(world, backend, writer_nproc,
                                          "load_fd_state", archive_basename);

  // Collective binary load — same primitive ES uses per-root.
  State s;
  s.responses.resize(1);
  s.responses[0] = Storage::load(world, archive_path, backend);

  s.last_bsh_residual     = {bsh_res};
  s.last_density_residual = {};      // not persisted; first-iter recomputes
  s.rho_alpha_prev        .clear();  // first-iter Δρ guard handles empty
  s.iter                  = 0;
  s.diverged              = false;

  if (world.rank() == 0) {
    madness::print("[LOAD] fd_state: pert=", pdesc,
                   "  protocol_key=", key,
                   "  freq=", freq,
                   "  archive=", archive_basename,
                   "  bsh_res_at_save=", bsh_res,
                   "  iter_at_save=", iter_at_save,
                   "  converged_at_save=", converged_int != 0,
                   "  backend=", io_backend_tag(backend));
  }
  return s;
}

/// Result of try_load_fd_state. The loader itself re-projects a non-exact
/// (coarser-rung) source to the ACTIVE (k, thresh) before returning, so every
/// consumer — restart seeding, the VBC build, alpha/beta property assembly —
/// receives functions at the active key and inner() is always well-defined
/// (mixing k is a hard MADNESS abort). `exact` / `source_protocol_key` remain
/// the honest provenance record: re-projecting a coarser state UP makes it
/// representable at the active k but adds NO accuracy.
template <typename Type, typename Shell>
struct FDRestartResult {
  typename FDSolver<Type, Shell>::State state;
  std::string source_protocol_key;
  bool        exact = false;
};

/// Restart precedence (Inc 13c-iii). Reads response_metadata.json and picks
/// the best on-disk bundle for (pert, freq):
///
///   1. exact match at the active protocol_key   → exact = true
///   2. else any saved (thresh, k) that's COARSER-OR-EQUAL to active
///      (saved_thresh >= active_thresh AND saved_k <= active_k):
///      pick max k, then min thresh — closest to active = best initial guess
///   3. else nullopt — caller seeds fresh
///
/// Returns nullopt if no metadata file or no compatible record. The
/// archive load itself is collective, so this function must be called by
/// every rank.
template <typename Type, typename Shell>
std::optional<FDRestartResult<Type, Shell>>
try_load_fd_state(madness::World &world,
                  const std::string &dir,
                  const Perturbation &pert,
                  double freq) {
  using State   = typename FDSolver<Type, Shell>::State;
  using Storage = typename FDSolver<Type, Shell>::Storage;

  const std::string active_key = protocol_key();
  const double active_thresh   = madness::FunctionDefaults<3>::get_thresh();
  const int    active_k        = madness::FunctionDefaults<3>::get_k();
  const std::string fkey       = ResponseMetadata::freq_key(freq);
  const std::string pdesc      = pert.description();

  // Rank-0 picks the source_key. Empty string ↔ no compatible match. Selection
  // goes through the shared best_usable_fd_source_key helper — the SAME one
  // reconcile_protocol uses — so the load can never pick a source the reconcile
  // verdict didn't sanction. It also excludes `diverged` snapshots (never seed a
  // blown-up state) and accepts coarser PARTIALS, not just converged ones.
  std::string source_key;
  std::string backend_tag;   // chosen entry's recorded backend ("" = legacy)
  int         writer_nproc = 0;
  if (world.rank() == 0) {
    const std::string meta_path = dir + "/response_metadata.json";
    if (std::filesystem::exists(meta_path)) {
      auto meta = ResponseMetadata::load_or_create(meta_path);
      source_key = best_usable_fd_source_key(meta.json(), pdesc, fkey,
                                             active_thresh, active_k, active_key);
      if (!source_key.empty()) {
        const auto &e = meta.json()["fd_states"][pdesc][source_key][fkey];
        backend_tag  = e.value("backend", std::string{});
        writer_nproc = e.value("writer_nproc", 0);
      }
    }
  }
  world.gop.broadcast_serializable(source_key, 0);
  world.gop.broadcast_serializable(backend_tag, 0);
  world.gop.broadcast(writer_nproc, 0);

  if (source_key.empty()) {
    if (world.rank() == 0) {
      madness::print("[LOAD] try_load_fd_state: no compatible bundle for pert=",
                     pdesc, " freq=", freq, " in ", dir,
                     " (active_key=", active_key, ")");
    }
    return std::nullopt;
  }

  // Collective load, at the entry's recorded backend (np-guarded for native).
  const std::string archive_basename = response_filename(pdesc, source_key, freq);
  const std::string archive_path     = dir + "/" + archive_basename;

  const IoBackend backend = io_backend_from_tag(backend_tag);
  detail_fd_save_load::check_writer_nproc(world, backend, writer_nproc,
                                          "try_load_fd_state", archive_basename);

  FDRestartResult<Type, Shell> r;
  r.state.responses.resize(1);
  r.state.responses[0] = Storage::load(world, archive_path, backend);
  r.state.last_bsh_residual    .clear();
  r.state.last_density_residual.clear();
  r.state.rho_alpha_prev       .clear();
  r.state.iter     = 0;
  r.state.diverged = false;
  r.source_protocol_key = source_key;
  r.exact = (source_key == active_key);

  // k-CONSISTENCY (review fix, confirmed HIGH): a coarser-rung archive comes
  // off disk at its SAVED k, not the active one. Re-project HERE — collective,
  // all ranks — so the contract "loader returns active-key functions" holds
  // for every consumer at once (assemble_beta used to contract a coarse
  // fallback straight into inner() and abort with "functions have different k").
  if (!r.exact) {
    for (auto *blk : r.state.responses[0].blocks())
      for (auto &fn : *blk) fn = madness::project(fn, active_k, active_thresh);
    world.gop.fence();
  }

  if (world.rank() == 0) {
    madness::print("[LOAD] try_load_fd_state: pert=", pdesc,
                   "  freq=", freq,
                   "  source_protocol_key=", source_key,
                   "  active=", active_key,
                   "  exact=", r.exact,
                   "  archive=", archive_basename,
                   "  backend=", io_backend_tag(backend));
  }
  return r;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_SAVE_LOAD_HPP
