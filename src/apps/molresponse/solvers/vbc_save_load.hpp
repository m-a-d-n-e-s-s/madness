#ifndef MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP
#define MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP

// =========================================================================
// VBC persistence (beta-ii-b). The VBC quadratic source is a
// ResponseStateXY<ClosedShell> (same Storage as an FD Full state); we save it
// as a parallel archive under the calc dir plus a vbc_states/<vbc_id>/<key>
// entry in response_metadata.json. VBC is rebuilt from its converged FD inputs
// at each protocol, so there is no try_load_vbc — reconcile's Skip (converged at
// this protocol_key) is what avoids recomputation; restart is handled by the
// FD inputs already being on disk.
// =========================================================================

#include "../ResponseProtocol.hpp"   // protocol_key()
#include "fd_save_load.hpp"          // detail_fd_save_load::check_writer_nproc
#include "response_metadata.hpp"
#include "response_state.hpp"
#include "state_metrics.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <filesystem>
#include <optional>
#include <string>

namespace molresponse_v3 {

/// Filesystem-safe archive basename for a VBC node at a protocol:
///   "vbc:dipole_x__dipole_y@fB_fC"  ->  "vbc_dipole_x__dipole_y_fB_fC__<key>"
inline std::string vbc_archive_basename(const std::string &vbc_id,
                                        const std::string &key) {
  std::string base = vbc_id;
  for (char &ch : base)
    if (ch == ':' || ch == '@') ch = '_';
  return base + "__" + key;
}

/// Save one VBC source: collective archive write + rank-0 metadata upsert.
/// Assumes the active FunctionDefaults<3> (k, thresh) is the build protocol.
template <class Shell>
inline void save_vbc_state(madness::World &world,
                           const ResponseStateXY<Shell> &vbc,
                           const std::string &dir,
                           const std::string &vbc_id,
                           bool converged,
                           double wall_s = 0.0,     // R1b: build wall time
                           const std::string &metadata_shard = std::string(), // F2g
                           const std::string &log_prefix = std::string(),     // F2d
                           int log_group = -1) {                               // F2d
  // F2g: in subworld mode the VBC metadata upsert goes to a per-group shard
  // (response_metadata.group<tag>.json), merged by rank 0 after the fence — same
  // discipline as save_fd_state. "" = canonical file (single-World path).
  if (world.rank() == 0) std::filesystem::create_directories(dir);
  world.gop.fence();

  const double      thresh = madness::FunctionDefaults<3>::get_thresh();
  const int         k_now  = madness::FunctionDefaults<3>::get_k();
  const std::string key    = protocol_key(thresh, k_now);
  const std::string base   = vbc_archive_basename(vbc_id, key);
  const std::string archive_path = dir + "/" + base;

  // Collective; the returned backend is recorded so load_vbc opens the same
  // physical file family (see fd_save_load — same stale-twin discipline).
  const IoBackend backend = vbc.save(world, archive_path);
  StateMetrics metrics = measure_state(world, vbc, /*iter=*/0);
  metrics.wall_s = wall_s;   // R1b (uniform with FD/ES)

  if (world.rank() == 0) {
    const std::string meta_path =
        dir + "/response_metadata" +
        (metadata_shard.empty() ? "" : ".group" + metadata_shard) + ".json";
    auto meta = ResponseMetadata::load_or_create(meta_path);
    if (!meta.json()["protocols"].contains(key))
      meta.set_protocol(key, thresh, k_now, /*index=*/-1);
    nlohmann::json entry = {
        {"converged", converged},
        {"diverged",  false},
        {"archive",   base},
        // Backend + writer count, same contract as fd_states entries: load
        // opens the recorded backend; native archives are np-locked.
        {"backend",      io_backend_tag(backend)},
        {"writer_nproc", static_cast<int>(world.size())},
        {"metrics",   metrics.to_json()},
    };
    meta.set_vbc_state(vbc_id, key, entry);
    meta.save();
    // F2d: prepend the per-subworld tag (empty ⇒ unchanged, G=0 byte-identical).
    if (log_prefix.empty())
      madness::print("[SAVE] vbc_state: id=", vbc_id, "  protocol_key=", key,
                     "  archive=", base, "  converged=", converged);
    else
      madness::print(log_prefix, "[SAVE] vbc_state: id=", vbc_id,
                     "  protocol_key=", key, "  archive=", base,
                     "  converged=", converged);
    // R1b: uniform memory high-water mark (see fd_save_load). F2d group= field.
    if (log_group >= 0)
      madness::print("MEMORY_HWM  kind=vbc  protocol=", key,
                     "  rss_gb_max=", metrics.rss_gb, "  coeffs=", metrics.coeffs,
                     "  wall_s=", wall_s, "  id=", vbc_id, "  group=", log_group);
    else
      madness::print("MEMORY_HWM  kind=vbc  protocol=", key,
                     "  rss_gb_max=", metrics.rss_gb, "  coeffs=", metrics.coeffs,
                     "  wall_s=", wall_s, "  id=", vbc_id);
  }
  world.gop.fence();
}

/// Load a VBC source built at the ACTIVE protocol (exact key) for the contraction.
/// Returns nullopt if no converged vbc_states/<id>/<active-key> entry exists.
/// NOTE: exact-key-only by design — unlike try_load_fd_state there is no
/// coarser-rung fallback, so the returned functions are always at the active
/// (k, thresh) and need no re-projection (k-consistency contract).
template <class Shell>
inline std::optional<ResponseStateXY<Shell>>
load_vbc(madness::World &world, const std::string &dir, const std::string &vbc_id) {
  const std::string key = protocol_key();
  std::string archive;
  std::string backend_tag;   // recorded backend ("" = legacy -> auto-detect)
  int         writer_nproc = 0;
  if (world.rank() == 0) {
    const std::string mp = dir + "/response_metadata.json";
    if (std::filesystem::exists(mp)) {
      auto meta = ResponseMetadata::load_or_create(mp);
      const auto &j = meta.json();
      if (j.contains("vbc_states") && j["vbc_states"].contains(vbc_id) &&
          j["vbc_states"][vbc_id].contains(key)) {
        const auto &e = j["vbc_states"][vbc_id][key];
        if (e.value("converged", false)) {
          archive      = e.value("archive", std::string{});
          backend_tag  = e.value("backend", std::string{});
          writer_nproc = e.value("writer_nproc", 0);
        }
      }
    }
  }
  world.gop.broadcast_serializable(archive, 0);
  world.gop.broadcast_serializable(backend_tag, 0);
  world.gop.broadcast(writer_nproc, 0);
  if (archive.empty()) return std::nullopt;
  const IoBackend backend = io_backend_from_tag(backend_tag);
  detail_fd_save_load::check_writer_nproc(world, backend, writer_nproc,
                                          "load_vbc", archive);
  return ResponseStateXY<Shell>::load(world, dir + "/" + archive, backend);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_VBC_SAVE_LOAD_HPP
