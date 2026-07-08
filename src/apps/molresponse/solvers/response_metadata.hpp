#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP

// =========================================================================
// ResponseMetadata — the unified response_metadata.json writer (doc 13).
//
// Single source of truth for the aggregate metadata file FD and ES
// persistence share. Schema:
//
//   {
//     "schema_version": 1,
//     "protocols":       { "<protocol_key>": {thresh, k, index} },
//     "fd_states":       { "<pert>": { "<protocol_key>": { "<freq_key>":
//                          {freq, type, shell, converged, iter,
//                           bsh_residual, archive} } } },
//     "excited_states":  { "<protocol_key>": {type, shell, n_roots,
//                          bundle_dir, converged, slot_permutation,
//                          roots[]} },
//     "properties":      { "<name>": { "<protocol_key>": [ {...}, ... ] } }
//   }
//
// The writer is pure JSON — no World, no MPI dependency. CALLERS guard with
//   if (world.rank() == 0) { ResponseMetadata::load_or_create(...) ... save(); }
//   world.gop.fence();
// so other ranks observe the file after rank 0 has finished writing.
//
// Atomic-write contract: save() writes to `<path>.tmp` then renames over
// `<path>`. On POSIX the rename is atomic within a filesystem, so a crash
// mid-write leaves either the old file intact or the new one — never a
// half-written aggregate. Two concurrent writers on the same file are not
// supported; serialize at the orchestration layer.
//
// The writer takes payloads as raw nlohmann::json blobs rather than typed
// structs: keeps the FD/ES save paths in charge of what to record, and the
// writer in charge of where it lands. Helper builders (typed -> json) can
// live next to the solver code that owns the types.
// =========================================================================

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

class ResponseMetadata {
public:
  static constexpr int kSchemaVersion = 1;

  /// Load the file if it exists, otherwise seed a fresh schema. Rank-0 only
  /// (filesystem op); the caller is responsible for guarding with rank().
  static ResponseMetadata load_or_create(const std::string &path) {
    ResponseMetadata m;
    m.path_ = path;
    if (std::filesystem::exists(path)) {
      std::ifstream in(path);
      if (!in) {
        throw std::runtime_error("ResponseMetadata: cannot read " + path);
      }
      in >> m.j_;
      const int v = m.j_.value("schema_version", 0);
      if (v != kSchemaVersion) {
        throw std::runtime_error(
            "ResponseMetadata: schema_version " + std::to_string(v) +
            " in " + path + " not recognized (expected " +
            std::to_string(kSchemaVersion) + ")");
      }
      // Ensure required top-level keys exist for upsert paths below.
      for (const char *k : {"protocols", "fd_states",
                             "excited_states", "properties", "vbc_states"}) {
        if (!m.j_.contains(k)) m.j_[k] = nlohmann::json::object();
      }
    } else {
      m.j_ = nlohmann::json{
          {"schema_version", kSchemaVersion},
          {"protocols",      nlohmann::json::object()},
          {"fd_states",      nlohmann::json::object()},
          {"excited_states", nlohmann::json::object()},
          {"properties",     nlohmann::json::object()},
          {"vbc_states",     nlohmann::json::object()},
      };
    }
    return m;
  }

  /// Upsert a protocol registry entry. Overwrites any existing entry at
  /// `key` (the (thresh, k) physical identity is stable — only `index` can
  /// legitimately change across runs with different ramps).
  void set_protocol(const std::string &key, double thresh, int k, int index) {
    j_["protocols"][key] = nlohmann::json{
        {"thresh", thresh}, {"k", k}, {"index", index}};
  }

  /// Upsert one FD point: fd_states/<pert>/<protocol_key>/<freq_key> = entry.
  void set_fd_state(const std::string &pert,
                    const std::string &protocol_key,
                    const std::string &freq_key,
                    const nlohmann::json &entry) {
    j_["fd_states"][pert][protocol_key][freq_key] = entry;
  }

  /// Replace the ES bundle entry at this protocol. ES is bundle-at-a-time —
  /// the whole entry is rewritten on save, not patched root-by-root.
  void set_es_bundle(const std::string &protocol_key,
                     const nlohmann::json &entry) {
    j_["excited_states"][protocol_key] = entry;
  }

  /// Upsert one VBC quadratic source: vbc_states/<vbc_id>/<protocol_key> = entry.
  void set_vbc_state(const std::string &vbc_id,
                     const std::string &protocol_key,
                     const nlohmann::json &entry) {
    j_["vbc_states"][vbc_id][protocol_key] = entry;
  }

  /// Record the effective I/O configuration (provenance): which restart backend
  /// wrote this run's state archives, and whether the binary had HDF5 compiled
  /// in. Without this an HDF5 run and a native run are indistinguishable.
  void set_io_info(const std::string &backend, bool hdf5_compiled) {
    j_["io"] = {{"backend", backend}, {"hdf5_compiled", hdf5_compiled}};
  }

  /// Append a property record to properties/<name>/<protocol_key>[]. The
  /// caller stamps it with whatever provenance is meaningful (es_root_id,
  /// fd_freq, the value, etc. — see doc 13's matching contract).
  void add_property(const std::string &name,
                    const std::string &protocol_key,
                    const nlohmann::json &record) {
    auto &arr = j_["properties"][name][protocol_key];
    if (!arr.is_array()) arr = nlohmann::json::array();
    arr.push_back(record);
  }

  /// Atomic write. Throws on filesystem error.
  void save() const {
    const std::string tmp = path_ + ".tmp";
    {
      std::ofstream out(tmp);
      if (!out) throw std::runtime_error(
          "ResponseMetadata: cannot open for write: " + tmp);
      out << j_.dump(2) << "\n";
    }
    std::filesystem::rename(tmp, path_);
  }

  /// Direct read access — for tests and the (rare) case a caller needs to
  /// query something the typed setters don't cover.
  const nlohmann::json &json() const { return j_; }
  const std::string    &path() const { return path_; }

  /// F2 (doc 32 §5.3): merge the per-group state metadata shards written by the
  /// subworlds (response_metadata.group<gid>.json) into the canonical file. The
  /// FD (and F2g VBC) states are DISJOINT across subworlds — each (pert/freq) FD
  /// point and each VBC id lives in exactly one subworld — so this is a
  /// conflict-free union, done through the typed setters (never a raw write).
  /// Removes the shards after a successful canonical save. Rank-0 only (caller
  /// guards). Idempotent: a missing shard is skipped, so re-running is safe.
  static void merge_state_shards(const std::string &calc_dir, int n_groups) {
    auto canon = load_or_create(calc_dir + "/response_metadata.json");
    for (int g = 0; g < n_groups; ++g) {
      const std::string sp =
          calc_dir + "/response_metadata.group" + std::to_string(g) + ".json";
      if (!std::filesystem::exists(sp)) continue;
      std::ifstream in(sp);
      if (!in) continue;
      nlohmann::json sj;
      in >> sj;
      // NB: iterate the lvalue member (sj["fd_states"]), NOT a .value(...)
      // temporary — .items() on a temporary json dangles.
      if (sj.contains("fd_states") && sj["fd_states"].is_object())
        for (const auto &pert : sj["fd_states"].items())
          for (const auto &pk : pert.value().items())
            for (const auto &fk : pk.value().items())
              canon.set_fd_state(pert.key(), pk.key(), fk.key(), fk.value());
      // F2g: VBC quadratic-source states (vbc_states/<id>/<protocol_key>).
      if (sj.contains("vbc_states") && sj["vbc_states"].is_object())
        for (const auto &id : sj["vbc_states"].items())
          for (const auto &pk : id.value().items())
            canon.set_vbc_state(id.key(), pk.key(), pk.value());
      // Protocol registry is informational; union any keys the canonical lacks.
      if (sj.contains("protocols") && sj["protocols"].is_object())
        for (const auto &p : sj["protocols"].items())
          if (!canon.j_["protocols"].contains(p.key()))
            canon.j_["protocols"][p.key()] = p.value();
    }
    canon.save();
    for (int g = 0; g < n_groups; ++g)
      std::filesystem::remove(
          calc_dir + "/response_metadata.group" + std::to_string(g) + ".json");
  }

  /// Canonical frequency key for fd_states / archive naming. f%.5f matches
  /// the dimensionless precision used in v2 and is enough to distinguish
  /// any physically relevant frequency we'd solve at. Doc 13.
  static std::string freq_key(double freq) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "f%.5f", freq);
    return {buf};
  }

private:
  nlohmann::json j_;
  std::string    path_;
};

/// Pick the best usable FD restart/seed source for (pert, freq_key) at target
/// protocol (target_thresh, target_k), from the aggregate metadata json.
/// "Usable" = present in fd_states, listed in the protocols registry,
/// COARSER-OR-EQUAL to the target (saved_thresh >= target_thresh AND
/// saved_k <= target_k), and NOT `diverged`. Convergence is NOT required:
/// a coarse partial is a perfectly good seed for the finer step. Preference:
/// exact target_key first, else closest-to-target (max k, then min thresh).
/// Returns "" if no usable source exists.
///
/// SINGLE SOURCE OF TRUTH for both the reconcile verdict (calc_manager's
/// reconcile_protocol) and the archive load (fd_save_load's try_load_fd_state),
/// so the two can never disagree — the doc-15 "share one pick-best-usable-source
/// helper" contract. Pure: reads only the json, no World / filesystem.
inline std::string
best_usable_fd_source_key(const nlohmann::json &meta, const std::string &pert,
                          const std::string &freq_key, double target_thresh,
                          int target_k, const std::string &target_key) {
  if (!meta.contains("fd_states") || !meta["fd_states"].contains(pert))
    return {};
  if (!meta.contains("protocols") || !meta["protocols"].is_object()) return {};
  const auto &protos = meta["protocols"];

  struct Cand { std::string key; double thresh; int k; };
  std::vector<Cand> cands;
  for (const auto &[key, ent] : meta["fd_states"][pert].items()) {
    if (!ent.contains(freq_key))                continue;
    if (ent[freq_key].value("diverged", false)) continue;  // never seed a blown-up state
    if (!protos.contains(key))                  continue;
    const double t  = protos[key].value("thresh", 0.0);
    const int    kk = protos[key].value("k", 0);
    if (t >= target_thresh && kk <= target_k) cands.push_back({key, t, kk});
  }
  if (cands.empty()) return {};

  for (const auto &c : cands)                    // exact wins
    if (c.key == target_key) return c.key;
  return std::max_element(                        // else closest-to-target
             cands.begin(), cands.end(),
             [](const Cand &a, const Cand &b) {
               if (a.k != b.k) return a.k < b.k;       // higher k wins
               return a.thresh > b.thresh;             // smaller thresh wins
             })
      ->key;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_METADATA_HPP
