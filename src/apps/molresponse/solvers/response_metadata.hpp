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
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
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

  /// Attach the seeded-solve root-identity guard block to the ES bundle entry
  /// at this protocol: excited_states/<protocol_key>/seed_guard = entry. Kept
  /// separate from set_es_bundle (which full-replaces the bundle entry on every
  /// save) because the guard is evaluated once, AFTER the final save.
  void set_es_seed_guard(const std::string &protocol_key,
                         const nlohmann::json &entry) {
    j_["excited_states"][protocol_key]["seed_guard"] = entry;
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

  /// Record external-seed provenance (the `dalton.dir` import contract,
  /// showcase W3): which DALTON calculation seeded this dir's initial guesses
  /// — {dalton_dir, basis, method, geometry_hash, mra_writer{version,commit}}.
  /// Full-replace upsert (one seed source per calc dir by contract).
  void set_seeded_from(const nlohmann::json &prov) {
    j_["seeded_from"] = prov;
  }

  /// Stamp the ground-state archive identity this calc dir's response states
  /// were built against (restart safety — the GS fingerprint gate). `hex` is
  /// the FNV-1a-64 of the archive part bytes; a later run with a different
  /// archive (regenerated orbitals => possible phase flips) must not reuse
  /// the cached response vectors.
  void set_ground_state(const std::string &archive, const std::string &hex,
                        std::uint64_t bytes, int nparts) {
    j_["ground_state"] = {{"archive", archive},
                          {"fnv1a64", hex},
                          {"bytes", bytes},
                          {"nparts", nparts}};
  }

  /// The stamped GS fingerprint, or "" if this metadata predates the stamp.
  std::string ground_state_fingerprint() const {
    if (j_.contains("ground_state") && j_["ground_state"].is_object())
      return j_["ground_state"].value("fnv1a64", "");
    return {};
  }

  /// Upsert a property record into properties/<name>/<protocol_key>[]. The
  /// caller stamps it with whatever provenance is meaningful (es_root_id,
  /// fd_freq, the value, etc. — see doc 13's matching contract).
  ///
  /// Identity contract (review fix — append-only rows made restarts accumulate
  /// duplicate/stale rows): a record's identity within (name, protocol_key) is
  /// the tuple of its identity fields (see property_identity). A re-run that
  /// re-assembles the same row REPLACES the old one; rows with different
  /// identities append; history across DIFFERENT protocols is preserved by the
  /// protocol_key level. Records carrying none of the identity fields keep the
  /// legacy append behaviour.
  void add_property(const std::string &name,
                    const std::string &protocol_key,
                    const nlohmann::json &record) {
    auto &arr = j_["properties"][name][protocol_key];
    if (!arr.is_array()) arr = nlohmann::json::array();
    const nlohmann::json id = property_identity(record);
    if (!id.empty())
      for (auto &existing : arr)
        if (property_identity(existing) == id) { existing = record; return; }
    arr.push_back(record);
  }

  /// The identity fields of a property record: alpha rows key on
  /// (omega, directions); beta/raman rows on (A, B, C, freq_b, freq_c);
  /// ES-derived rows on (es_root_id, fd_freq). Doubles round-trip exactly
  /// through nlohmann::json, so equality on re-assembled rows is exact.
  static nlohmann::json property_identity(const nlohmann::json &record) {
    static constexpr const char *keys[] = {"omega",  "directions", "A",
                                           "B",      "C",          "freq_b",
                                           "freq_c", "es_root_id", "fd_freq"};
    nlohmann::json id = nlohmann::json::object();
    for (const char *k : keys)
      if (record.contains(k)) id[k] = record[k];
    return id;
  }

  /// Record planned work the run ended WITHOUT executing (review fix: honest-
  /// climb + hard VBC prerequisite gates could silently drop beta/raman work
  /// while stop_reason claimed 'complete'). Full-replace upsert: run() writes
  /// the current list on every exit, so a later run that completes the work
  /// clears it (empty array = nothing dropped).
  void set_dropped_work(const nlohmann::json &items) {
    j_["run_summary"]["dropped_work"] = items;
  }

  /// Process-wide read-only switch. Concurrent per-root verification
  /// processes (--tpa-roots) share one calc dir; with this set every save()
  /// is a no-op so siblings cannot race on the index (in-memory mutations
  /// like the derived-FD expansion still happen and are read back normally).
  static bool &read_only() { static bool ro = false; return ro; }

  /// Atomic write. Throws on filesystem error — INCLUDING a failed write or
  /// close (ENOSPC etc.): an unchecked stream would let the rename install a
  /// truncated tmp over a good index. The bad tmp is removed on failure.
  void save() const {
    if (read_only()) return;
    // Fixed tmp name, deliberately: a stranded tmp from a killed save is
    // CONSUMED by the next save() rather than accumulating (enforced by
    // test_response_metadata "re-save consumed/replaced the stranded tmp").
    // Concurrent writers on one file are not supported by contract — the
    // read_only() switch above is how processes that share a calc dir avoid
    // writing at all.
    const std::string tmp = path_ + ".tmp";
    {
      std::ofstream out(tmp);
      if (!out) throw std::runtime_error(
          "ResponseMetadata: cannot open for write: " + tmp);
      out << j_.dump(2) << "\n";
      out.close();
      if (out.fail()) {
        std::error_code ec;
        std::filesystem::remove(tmp, ec);
        throw std::runtime_error(
            "ResponseMetadata: write to " + tmp +
            " failed (disk full?) — keeping the previous " + path_);
      }
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
    for (int g = 0; g < n_groups; ++g)
      merge_shard_into(canon, calc_dir + "/response_metadata.group" +
                                  std::to_string(g) + ".json");
    canon.save();
    for (int g = 0; g < n_groups; ++g)
      std::filesystem::remove(
          calc_dir + "/response_metadata.group" + std::to_string(g) + ".json");
  }

  /// F2 restart safety: merge any per-group shards STRANDED by an interrupted
  /// run. The normal path (merge_state_shards above) collapses shards right
  /// after each fan-out wave — but a kill between fan-out and merge leaves
  /// response_metadata.group<g>.json files behind, and a restart reading only
  /// the canonical file would silently re-solve states that are already on
  /// disk. Discovers shards by name (the interrupted run's G is unknown and
  /// may differ from this run's), merges through the same typed-setter union,
  /// and removes them. Returns the number of shards merged; 0 = nothing to
  /// do. Rank-0 only (caller guards); idempotent.
  static int merge_stale_state_shards(const std::string &calc_dir) {
    namespace fs = std::filesystem;
    const std::string prefix = "response_metadata.group";
    const std::string suffix = ".json";
    std::vector<std::string> shards;
    std::error_code ec;
    for (const auto &e : fs::directory_iterator(calc_dir, ec)) {
      if (ec || !e.is_regular_file()) continue;
      const std::string name = e.path().filename().string();
      if (name.size() <= prefix.size() + suffix.size()) continue;
      if (name.compare(0, prefix.size(), prefix) != 0) continue;
      if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
        continue;
      const std::string gid =
          name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
      if (gid.empty() ||
          gid.find_first_not_of("0123456789") != std::string::npos)
        continue;  // not a shard (e.g. a stray .json.tmp or foreign file)
      shards.push_back(e.path().string());
    }
    if (shards.empty()) return 0;
    std::sort(shards.begin(), shards.end());
    auto canon = load_or_create(calc_dir + "/response_metadata.json");
    for (const auto &sp : shards) merge_shard_into(canon, sp);
    canon.save();
    for (const auto &sp : shards) fs::remove(sp);
    return static_cast<int>(shards.size());
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
  /// The single shard->canonical union both merge entries share. FD and VBC
  /// states are DISJOINT across subworlds, so this is conflict-free; done
  /// through the typed setters (never a raw write). Missing/unreadable shard
  /// files are skipped (idempotence).
  static void merge_shard_into(ResponseMetadata &canon, const std::string &sp) {
    if (!std::filesystem::exists(sp)) return;
    std::ifstream in(sp);
    if (!in) return;
    nlohmann::json sj;
    // Shards are written atomically (save() tmp+rename), so a parse failure
    // means external damage — quarantine it (keep the evidence, unblock every
    // future startup) rather than throwing the whole run down.
    try {
      in >> sj;
    } catch (const std::exception &) {
      in.close();
      std::error_code ec;
      std::filesystem::rename(sp, sp + ".corrupt", ec);
      std::fprintf(stderr,
                   "ResponseMetadata: shard %s unparsable — quarantined as "
                   "%s.corrupt (its states will be re-solved)\n",
                   sp.c_str(), sp.c_str());
      return;
    }
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
