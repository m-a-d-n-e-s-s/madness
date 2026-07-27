#ifndef MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP

// =========================================================================
// Solve / save / restart for ESSolver<Type, Shell>.
//
// On disk, an "ES roots directory" is laid out as:
//
//   my_dir/
//   ├── roots.json    ← sidecar: type, shell, n_roots, k, thresh, per-root
//   │                   metadata (slot, omega, residuals, file name)
//   ├── root_0        ← MADNESS parallel archive (one per slot)
//   ├── root_1
//   └── ...
//
// Per-root binary content uses the existing `ResponseStateX<Shell>::save`
// / `::load` machinery in response_state.hpp — this header only adds the
// sidecar + a directory-level wrapper.
//
// `roots.json` is the *index*: it tells the loader (a) what type/shell to
// expect, (b) which file holds which slot, (c) the ω / residual recorded
// at save time, (d) the per-root STABLE identity (stable_index / root_id /
// display_name) plus the bundle's slot_permutation (slot -> stable_index).
// ω is recorded for sanity (caller compares against a reference) — load
// does not permute by ω. The slot_permutation is what lets a restart match
// roots by stable identity rather than by raw array slot (doc 03, Inc 9).
//
// On load: `iter` is reset to 0 (fresh count for the resumed solver),
// `rho_alpha_prev` is left empty (ESSolver::step's first-iter check is
// already gated by `iter <= 1`). The KAIN buffer of the resumed solver
// starts empty too — that's the right call for a restart test.
//
// Cross-protocol restart: the bundle's native k/thresh is recorded but
// NOT enforced. The caller's protocol-step `prepare` hook re-projects
// each root to the active k/thresh via `madness::project`, same path the
// per-protocol stepper already uses.
//
// Collective discipline: parallel archive I/O is collective on `world`
// and must be called by every rank. The JSON sidecar read/write is
// rank-0-only (filesystem op); rank 0 broadcasts the parsed metadata.
// =========================================================================

#include "es_root_identity.hpp"
#include "es_solver.hpp"
#include "response_metadata.hpp"     // ResponseMetadata (doc 13 aggregate)
#include "response_state.hpp"
#include "state_metrics.hpp"         // process_rss_gb (bundle RSS)
#include "../ResponseProtocol.hpp"   // protocol_key(thresh, k)
#include "../kernels/tags.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cstdlib>       // std::getenv (MADRESPONSE_STRICT_NP)
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

namespace detail_save_load {

template <typename Type>
inline const char *type_tag() {
  if constexpr (std::is_same_v<Type, TDA>)  return "tda";
  else if constexpr (std::is_same_v<Type, Full>) return "full";
  else return "unknown";
}

template <typename Shell>
inline const char *shell_tag() {
  if constexpr (std::is_same_v<Shell, ClosedShell>) return "closed_shell";
  else if constexpr (std::is_same_v<Shell, OpenShell>) return "open_shell";
  else return "unknown";
}

inline std::string root_file(int slot) {
  return "root_" + std::to_string(slot);
}

} // namespace detail_save_load

/// Save an ESSolver state to a directory. Writes one parallel archive
/// per root + a `roots.json` sidecar with per-slot metadata.
template <typename Type, typename Shell>
void save_es_roots(madness::World &world,
                   const typename ESSolver<Type, Shell>::State &state,
                   const std::string &dir,
                   bool converged,
                   double wall_s = 0.0,            // R1b: bundle-solve wall time
                   bool register_aggregate = true) {
  const int n_roots = static_cast<int>(state.roots.size());
  MADNESS_CHECK(n_roots > 0);

  if (world.rank() == 0) {
    std::filesystem::create_directories(dir);
  }
  world.gop.fence();

  // Per-root binary archives — collective. All roots share one backend (the
  // opt-in is stable across a save); it is recorded in roots.json so the
  // loader opens the same file family (stale-twin discipline, see
  // response_state.hpp).
  IoBackend backend = IoBackend::Native;
  for (int s = 0; s < n_roots; ++s) {
    const std::string path = dir + "/" + detail_save_load::root_file(s);
    backend = state.roots[s].save(world, path);
  }

  // Per-root coefficient counts + worst-task bundle RSS — collective on
  // every rank, batched into a single sum reduce. coeffs[s] feeds the
  // memory-scaling study; rss_gb is the OOM-relevant per-task figure.
  std::vector<std::size_t> root_coeffs(static_cast<size_t>(n_roots), 0);
  for (int s = 0; s < n_roots; ++s)
    for (const auto &f : state.roots[s].flatten())
      root_coeffs[static_cast<size_t>(s)] += f.size_local();
  world.gop.sum(root_coeffs.data(), static_cast<size_t>(n_roots));
  double bundle_rss_gb = process_rss_gb();
  world.gop.max(bundle_rss_gb);

  if (world.rank() == 0) {
    const double thresh = madness::FunctionDefaults<3>::get_thresh();

    // Stable identity: use what the state carries, else assign 0..N-1 so a
    // save always records identity. Display names are grouped at this
    // protocol boundary (tol = 10 * thresh; matches doc 03).
    std::vector<int> stable_index = state.stable_index;
    assign_initial_stable_index(stable_index, n_roots);
    const std::vector<std::string> display_names =
        assign_display_names(state.omega, 10.0 * thresh);

    // R1b: uniform bundle metrics — sum the per-root coeffs, worst-task RSS,
    // iters, wall. The same StateMetrics block FD/VBC write, so every state
    // type reports the same shape (feeds the R4 memory-scaling model / L2).
    StateMetrics bundle_metrics;
    for (int s = 0; s < n_roots; ++s)
      bundle_metrics.coeffs += root_coeffs[static_cast<size_t>(s)];
    bundle_metrics.bytes  = bundle_metrics.coeffs * sizeof(double);
    bundle_metrics.rss_gb = bundle_rss_gb;
    bundle_metrics.iters  = state.iter;
    bundle_metrics.wall_s = wall_s;

    nlohmann::json j;
    j["type"]      = detail_save_load::type_tag<Type>();
    j["shell"]     = detail_save_load::shell_tag<Shell>();
    j["n_roots"]   = n_roots;
    // The per-root archives are nio=1 parallel archives; MADNESS can only reload
    // them with the same #processes that wrote them (it assumes #writers ==
    // #readers — a mismatch silently loads the wrong coefficient sets into a
    // wrong-pmap WorldContainer and corrupts the heap). Record the writer count
    // so load_es_roots can verify it (see the np-guard there). world.size() is a
    // local query, safe on rank 0.
    j["writer_nproc"] = world.size();
    // Which file family holds the per-root archives ("native"/"hdf5"); the
    // loader opens THIS backend (bare .h5 existence detect is legacy-only)
    // and skips the np-guard for hdf5 (single-client blob, np-portable).
    j["backend"] = io_backend_tag(backend);
    const int k_now = madness::FunctionDefaults<3>::get_k();
    j["k"]            = k_now;
    j["thresh"]       = thresh;
    j["protocol_key"] = protocol_key(thresh, k_now);  // doc 13 join key
    j["iter"]         = state.iter;
    j["converged"]    = converged;
    j["diverged"]     = state.diverged;

    // slot_permutation[slot] = stable_index — the cross-protocol root map.
    j["slot_permutation"] = stable_index;

    nlohmann::json roots_arr = nlohmann::json::array();
    const long n_omega = state.omega.size();
    for (int s = 0; s < n_roots; ++s) {
      nlohmann::json entry;
      entry["slot"]         = s;
      entry["stable_index"] = stable_index[s];
      entry["root_id"]      = make_root_id(stable_index[s]);
      entry["display_name"] =
          (static_cast<size_t>(s) < display_names.size()) ? display_names[s]
                                                          : std::string();
      entry["omega"] = (s < n_omega) ? state.omega(s)
                                     : std::numeric_limits<double>::quiet_NaN();
      entry["bsh_residual"] =
          (static_cast<size_t>(s) < state.last_bsh_residual.size())
              ? state.last_bsh_residual[s]
              : std::numeric_limits<double>::quiet_NaN();
      entry["density_residual"] =
          (static_cast<size_t>(s) < state.last_density_residual.size())
              ? state.last_density_residual[s]
              : std::numeric_limits<double>::quiet_NaN();
      entry["file"] = detail_save_load::root_file(s);
      entry["coeffs"] = root_coeffs[static_cast<size_t>(s)];
      entry["bytes"]  = root_coeffs[static_cast<size_t>(s)] * sizeof(double);
      roots_arr.push_back(entry);
    }
    j["roots"]  = roots_arr;
    j["rss_gb"] = bundle_rss_gb;  // worst-task RSS at this protocol (kept for compat)
    j["metrics"] = bundle_metrics.to_json();  // R1b: uniform metrics block

    // Atomic write (tmp + rename) so a crash mid-write can't leave a half-written
    // index that load_es_roots would mis-parse. The per-root archives above
    // already exist, so the index becoming visible last is the safe order.
    const std::string roots_json = dir + "/roots.json";
    const std::string roots_tmp  = roots_json + ".tmp";
    { std::ofstream out(roots_tmp); out << j.dump(2) << "\n"; }
    std::filesystem::rename(roots_tmp, roots_json);

    // 13d: upsert into the calc-level aggregate response_metadata.json — UNLESS
    // this is a side bundle (e.g. a cached warmup guess) that must not appear in
    // the shared excited_states index (it would collide with / be picked up as
    // the main es__<key> restart bundle). roots.json above is always written so
    // load_es_roots(dir) can read the side bundle directly.
    if (register_aggregate) {
      namespace fs = std::filesystem;
      fs::path bundle_path(dir);
      fs::path calc_dir = bundle_path.parent_path();
      if (calc_dir.empty()) calc_dir = ".";
      const std::string aggregate_path =
          (calc_dir / "response_metadata.json").string();

      auto meta = ResponseMetadata::load_or_create(aggregate_path);
      const std::string key = protocol_key(thresh, k_now);
      if (!meta.json()["protocols"].contains(key)) {
        meta.set_protocol(key, thresh, k_now, /*index=*/-1);
      }

      nlohmann::json bundle_entry = {
          {"type",             detail_save_load::type_tag<Type>()},
          {"shell",            detail_save_load::shell_tag<Shell>()},
          {"n_roots",          n_roots},
          {"backend",          io_backend_tag(backend)},
          {"bundle_dir",       bundle_path.filename().string()},
          {"converged",        converged},
          {"diverged",         state.diverged},
          {"slot_permutation", stable_index},
          {"roots",            roots_arr},
          {"iter",             state.iter},
          {"rss_gb",           bundle_rss_gb},
          {"metrics",          bundle_metrics.to_json()},  // R1b: uniform block
      };
      meta.set_es_bundle(key, bundle_entry);
      meta.save();

      madness::print("[SAVE] es_bundle: protocol_key=", key,
                     "  bundle_dir=", bundle_path.filename().string(),
                     "  n_roots=", n_roots,
                     "  aggregate=", aggregate_path);
      // R1b: uniform memory high-water mark (see fd_save_load).
      madness::print("MEMORY_HWM  kind=es  protocol=", key,
                     "  rss_gb_max=", bundle_metrics.rss_gb,
                     "  coeffs=", bundle_metrics.coeffs, "  wall_s=", wall_s,
                     "  id=", detail_save_load::type_tag<Type>(), "_n", n_roots);
    } else {
      madness::print("[SAVE] es side-bundle (no aggregate): dir=", dir,
                     "  n_roots=", n_roots);
    }
  }
  world.gop.fence();
}

/// Load an ESSolver state previously written by `save_es_roots`.
///
/// Validates that the sidecar's `type` / `shell` / `n_roots` match the
/// requested template instantiation. Per-slot ω + residuals are
/// restored; `iter` resets to 0 and `rho_alpha_prev` stays empty so the
/// resumed solver picks up from a clean count.
template <typename Type, typename Shell>
typename ESSolver<Type, Shell>::State
load_es_roots(madness::World &world, const std::string &dir) {
  using State = typename ESSolver<Type, Shell>::State;

  // Rank 0 reads + parses; broadcast the n_roots scalar and the per-slot
  // ω / residuals via Tensor<double> + vector<double> hops. The
  // expected type/shell strings are baked in at compile time, so we
  // only need rank 0 to assert them — no broadcast needed for those.
  int n_roots = 0;
  int iter_at_save = 0;
  int writer_nproc = 0;   // #processes that wrote the per-root archives (0 = legacy)
  std::string backend_tag;  // recorded backend ("" = legacy -> auto-detect)
  madness::Tensor<double> omega_recorded;
  std::vector<double> bsh_residual_recorded;
  std::vector<double> drho_residual_recorded;
  std::vector<int> stable_index_recorded;

  // Rank 0 decides the load status; the status is broadcast BEFORE anyone
  // throws so a bad/missing/mismatched roots.json fails COLLECTIVELY. A
  // rank-0-only throw here would leave every other rank blocked in the
  // broadcasts below (multi-rank hang) — same discipline as the np-guard
  // further down.
  std::string load_error;
  if (world.rank() == 0) {
    try {
      std::ifstream in(dir + "/roots.json");
      if (!in) {
        throw std::runtime_error("load_es_roots: cannot open " + dir +
                                 "/roots.json");
      }
      nlohmann::json j;
      in >> j;

      const std::string saved_type  = j.value("type",  std::string{});
      const std::string saved_shell = j.value("shell", std::string{});
      const std::string want_type   = detail_save_load::type_tag<Type>();
      const std::string want_shell  = detail_save_load::shell_tag<Shell>();
      if (saved_type != want_type || saved_shell != want_shell) {
        throw std::runtime_error(
            "load_es_roots: saved type/shell = " + saved_type + "/" +
            saved_shell + ", requested = " + want_type + "/" + want_shell);
      }

      n_roots      = j.value("n_roots", 0);
      iter_at_save = j.value("iter", 0);
      writer_nproc = j.value("writer_nproc", 0);   // 0 = legacy bundle (no count)
      backend_tag  = j.value("backend", std::string{});
      if (n_roots <= 0) {
        throw std::runtime_error("load_es_roots: n_roots <= 0 in " + dir +
                                 "/roots.json");
      }
      if (!j.contains("roots") || !j["roots"].is_array() ||
          static_cast<int>(j["roots"].size()) != n_roots) {
        throw std::runtime_error(
            "load_es_roots: roots[] array missing or wrong length in " +
            dir + "/roots.json");
      }

      omega_recorded         = madness::Tensor<double>(n_roots);
      bsh_residual_recorded  .assign(n_roots, 0.0);
      drho_residual_recorded .assign(n_roots, 0.0);
      stable_index_recorded  .assign(n_roots, 0);
      for (int s = 0; s < n_roots; ++s) {
        const auto &e = j["roots"][s];
        const int slot = e.value("slot", -1);
        if (slot != s) {
          throw std::runtime_error(
              "load_es_roots: roots[" + std::to_string(s) + "].slot = " +
              std::to_string(slot) + " (expected " + std::to_string(s) +
              ") in " + dir + "/roots.json");
        }
        omega_recorded(s)         = e.value("omega", 0.0);
        bsh_residual_recorded[s]  = e.value("bsh_residual", 0.0);
        drho_residual_recorded[s] = e.value("density_residual", 0.0);
        // Stable identity: legacy files predate it — fall back to slot so
        // the loaded bundle still carries a well-formed identity vector.
        stable_index_recorded[s]  = e.value("stable_index", s);
      }
    } catch (const std::exception &e) {
      load_error = e.what();
    } catch (...) {
      load_error = "load_es_roots: unknown error reading " + dir +
                   "/roots.json";
    }
  }
  world.gop.broadcast_serializable(load_error, 0);
  if (!load_error.empty()) throw std::runtime_error(load_error);  // collective

  // Broadcast the metadata to every rank.
  world.gop.broadcast(n_roots, 0);
  world.gop.broadcast(iter_at_save, 0);
  world.gop.broadcast(writer_nproc, 0);
  world.gop.broadcast_serializable(backend_tag, 0);
  if (world.rank() != 0) {
    omega_recorded         = madness::Tensor<double>(n_roots);
    bsh_residual_recorded  .assign(n_roots, 0.0);
    drho_residual_recorded .assign(n_roots, 0.0);
    stable_index_recorded  .assign(n_roots, 0);
  }
  world.gop.broadcast(omega_recorded.ptr(), n_roots, 0);
  world.gop.broadcast(bsh_residual_recorded.data(),  n_roots, 0);
  world.gop.broadcast(drho_residual_recorded.data(), n_roots, 0);
  world.gop.broadcast(stable_index_recorded.data(),  n_roots, 0);

  // np-guard (all ranks — writer_nproc/world.size() are rank-uniform here, so the
  // throw is collective and cannot deadlock). The nio=1 per-root archives can
  // only be reloaded by the same #processes that wrote them; a mismatch silently
  // corrupts the WorldContainer and aborts at teardown (this is the parked ES
  // heap-OOB, reproduced by --es-analyze-only --es-load-only at a different np).
  // Fail cleanly instead. writer_nproc==0 is a legacy bundle with no recorded
  // count — proceed with a warning (same-np is the common case).
  // Cross-process-count ES restart. The nio=1 per-root archives are np-portable
  // (ParallelInputArchive redistributes on read; verified bit-identical np=2->1
  // on the FD/ES bundles). The former hard block + the "parked ES heap-OOB"
  // predated the atomic-write + collective-throw I/O fixes that resolved the
  // real crashes; with those in, cross-np ES restart is safe. Proceed by default
  // (resource-count agnostic); MADRESPONSE_STRICT_NP=1 restores the hard error.
  const IoBackend backend = io_backend_from_tag(backend_tag);
  if (!io_backend_np_portable(backend) &&
      writer_nproc != 0 && writer_nproc != world.size()) {
    const char *strict = std::getenv("MADRESPONSE_STRICT_NP");
    if (strict && strict[0] == '1') {
      throw std::runtime_error(
          "load_es_roots: ES bundle in " + dir + " was written with " +
          std::to_string(writer_nproc) + " process(es) but loaded with " +
          std::to_string(world.size()) +
          " — MADRESPONSE_STRICT_NP=1 forbids cross-process-count native "
          "restart. Re-run with " + std::to_string(writer_nproc) +
          " rank(s), unset MADRESPONSE_STRICT_NP, or use the HDF5 backend.");
    }
    if (world.rank() == 0)
      madness::print("[LOAD] ES bundle in", dir, "written with", writer_nproc,
                     "rank(s), loading with", world.size(),
                     "— native nio=1 is np-portable; proceeding.");
  }
  if (writer_nproc == 0 && world.rank() == 0) {
    madness::print("[LOAD] WARNING: ES bundle in", dir,
                   "has no recorded writer_nproc (legacy); assuming it was "
                   "written with the current", world.size(),
                   "process(es) — a mismatch may still crash.");
  }

  // Per-root binary archives — collective, at the recorded backend.
  State s;
  s.roots.reserve(n_roots);
  for (int b = 0; b < n_roots; ++b) {
    const std::string path = dir + "/" + detail_save_load::root_file(b);
    s.roots.push_back(
        std::remove_reference_t<decltype(s.roots[0])>::load(world, path,
                                                            backend));
  }

  s.omega                 = std::move(omega_recorded);
  s.last_bsh_residual     = std::move(bsh_residual_recorded);
  s.last_density_residual = std::move(drho_residual_recorded);
  s.stable_index          = std::move(stable_index_recorded);
  s.iter                  = 0;       // fresh count for the resumed solver
  s.diverged              = false;
  // rho_alpha_prev left empty — first-iter density-residual check in
  // ESSolver::step() is already gated by iter <= 1.

  if (world.rank() == 0) {
    madness::print("[LOAD] es_roots: n_roots=", n_roots,
                   " iter_at_save=", iter_at_save,
                   " dir=", dir);
    for (int b = 0; b < n_roots; ++b) {
      madness::print("       slot=", b,
                     "  root_id=", make_root_id(s.stable_index[b]),
                     "  omega=", s.omega(b),
                     "  bsh_res=", s.last_bsh_residual[b],
                     "  drho_res=", s.last_density_residual[b]);
    }
  }

  return s;
}

/// ES restart precedence — the cross-protocol analog of try_load_fd_state.
/// Reads response_metadata.json at `calc_dir`, finds the best
/// `excited_states/<key>` bundle for the requested (Type, Shell), and
/// loads it via the existing load_es_roots (which reads the per-bundle
/// roots.json — loader's truth per the doc 13 authority split).
///
/// Precedence:
///   1. exact match at the active protocol_key       → exact=true
///   2. else coarser-or-equal saved (thresh, k):
///      pick max k, then min thresh — closest = best initial guess
///   3. else nullopt — caller falls back to fresh guess / warmup
///
/// Non-exact loads are auto-reprojected by iterate_protocol's first
/// prepare() — the same path the solver already takes every protocol
/// step, so no extra projection code is needed in the test driver.
///
/// `calc_dir` is the directory containing response_metadata.json (one
/// level up from any es_bundle subdir). The aggregate's `bundle_dir`
/// field gives the basename of the bundle subdir to load.
template <typename Type, typename Shell>
struct ESRestartResult {
  typename ESSolver<Type, Shell>::State state;
  std::string source_protocol_key;
  std::string bundle_dir;            // basename, relative to calc_dir
  bool        exact = false;
};

template <typename Type, typename Shell>
std::optional<ESRestartResult<Type, Shell>>
try_load_es_bundle(madness::World &world, const std::string &calc_dir) {
  const std::string active_key = protocol_key();
  const double active_thresh   = madness::FunctionDefaults<3>::get_thresh();
  const int    active_k        = madness::FunctionDefaults<3>::get_k();
  const std::string want_type  = detail_save_load::type_tag<Type>();
  const std::string want_shell = detail_save_load::shell_tag<Shell>();

  // Rank-0 picks (source_key, bundle_dir). Both empty ↔ no match.
  std::string source_key;
  std::string bundle_dir;
  if (world.rank() == 0) {
    const std::string meta_path = calc_dir + "/response_metadata.json";
    if (std::filesystem::exists(meta_path)) {
      auto meta = ResponseMetadata::load_or_create(meta_path);
      const auto &j = meta.json();
      if (j.contains("excited_states") && j["excited_states"].is_object() &&
          j.contains("protocols")      && j["protocols"].is_object()) {
        struct Cand { std::string key; std::string bdir; double thresh; int k; };
        std::vector<Cand> cands;
        for (const auto &[key, ent] : j["excited_states"].items()) {
          // Type/shell must match the requested instantiation — the loader
          // also validates, but filter here so we don't pick an unloadable
          // candidate over a loadable one.
          if (ent.value("type",  std::string{}) != want_type)  continue;
          if (ent.value("shell", std::string{}) != want_shell) continue;
          if (!j["protocols"].contains(key))                    continue;
          const double t  = j["protocols"][key].value("thresh", 0.0);
          const int    kk = j["protocols"][key].value("k", 0);
          if (t >= active_thresh && kk <= active_k) {
            cands.push_back({
                key,
                ent.value("bundle_dir", std::string{}),
                t, kk});
          }
        }
        // Exact match wins.
        for (const auto &c : cands) {
          if (c.key == active_key) {
            source_key = c.key;
            bundle_dir = c.bdir;
            break;
          }
        }
        // Else closest-to-active: max k, then min thresh.
        if (source_key.empty() && !cands.empty()) {
          auto best = std::max_element(
              cands.begin(), cands.end(),
              [](const Cand &a, const Cand &b) {
                if (a.k != b.k)      return a.k < b.k;        // higher k wins
                return a.thresh > b.thresh;                   // smaller thresh wins
              });
          source_key = best->key;
          bundle_dir = best->bdir;
        }
      }
    }
  }
  world.gop.broadcast_serializable(source_key, 0);
  world.gop.broadcast_serializable(bundle_dir, 0);

  if (source_key.empty() || bundle_dir.empty()) {
    if (world.rank() == 0) {
      madness::print("[LOAD] try_load_es_bundle: no compatible bundle in ",
                     calc_dir, " for type=", want_type,
                     " shell=", want_shell,
                     " active_key=", active_key);
    }
    return std::nullopt;
  }

  // Collective load — delegate to load_es_roots (reads roots.json + per-
  // root archives). Aggregate's `bundle_dir` is a basename, so resolve
  // against calc_dir.
  const std::string bundle_path =
      (std::filesystem::path(calc_dir) / bundle_dir).string();

  ESRestartResult<Type, Shell> r;
  r.state = load_es_roots<Type, Shell>(world, bundle_path);
  r.source_protocol_key = source_key;
  r.bundle_dir          = bundle_dir;
  r.exact               = (source_key == active_key);

  if (world.rank() == 0) {
    madness::print("[LOAD] try_load_es_bundle: source_protocol_key=", source_key,
                   "  active=", active_key,
                   "  exact=", r.exact,
                   "  bundle_dir=", bundle_dir);
  }
  return r;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_ES_SAVE_LOAD_HPP
