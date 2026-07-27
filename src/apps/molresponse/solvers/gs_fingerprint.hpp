#ifndef MOLRESPONSE_V3_GS_FINGERPRINT_HPP
#define MOLRESPONSE_V3_GS_FINGERPRINT_HPP

// gs_fingerprint.hpp — ground-state archive identity for restart safety.
//
// WHY: a response calc dir caches converged response vectors (X, Y) keyed by
// protocol, and every one of them is expressed relative to ONE set of
// ground-state orbitals. HF/KS orbitals are only determined up to a
// per-orbital sign (and rotations among degenerate orbitals) — a REGENERATED
// moldft archive can be physically identical yet phase-flipped, and restart
// vectors from the old archive then assemble into silently corrupted
// properties. Observed in the wild: a calc dir holding two GS archives
// produced a converged-looking alpha with relative error exactly 2.0 (pure
// sign flip). Nothing crashes; the number is just wrong.
//
// THE GATE: the first run stamps the archive's byte-level fingerprint into
// response_metadata.json (ground_state/ block, via the metadata layer); every
// later run recomputes it and refuses to proceed on mismatch. Byte-level (not
// physics-level) identity is exactly right here: a byte-identical copy keeps
// its phases, a regenerated archive does not — even when the physics agrees.
//
// This header is pure filesystem + JSON — no World, no MPI (unit-testable in
// the `unit` ctest tier). The collective wrapper (rank-0 compute + broadcast
// + all-ranks abort) lives in orchestrator/response_workflow.hpp.

#include <madness/external/nlohmann_json/json.hpp>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

struct GsFingerprint {
  std::string   hex;        // 16-char lowercase FNV-1a-64 of all part bytes
  std::uint64_t bytes = 0;  // total bytes hashed
  int           nparts = 0; // archive part files found
};

inline std::uint64_t fnv1a64_update(std::uint64_t h, const char *p,
                                    std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<unsigned char>(p[i]);
    h *= 1099511628211ULL;
  }
  return h;
}

inline constexpr std::uint64_t kFnv1a64Basis = 14695981039346656037ULL;

/// The physical files behind a MADNESS archive base name: either the bare
/// file itself, or the ParallelOutputArchive parts "<base>.00000", ".00001",
/// ... in ascending order (the "%s.%5.5d" client naming).
inline std::vector<std::string> gs_archive_parts(const std::string &base) {
  namespace fs = std::filesystem;
  std::vector<std::string> parts;
  if (fs::exists(base) && fs::is_regular_file(base)) {
    parts.push_back(base);
    return parts;
  }
  for (int i = 0;; ++i) {
    char suf[16];
    std::snprintf(suf, sizeof(suf), ".%05d", i);
    const std::string p = base + suf;
    if (!fs::exists(p)) break;
    parts.push_back(p);
  }
  return parts;
}

/// FNV-1a-64 over the archive part bytes, streamed in order. Throws if no
/// part files exist (callers run this after a successful GS load, so absence
/// means the base path is wrong, not that the archive is optional).
inline GsFingerprint gs_archive_fingerprint(const std::string &base) {
  const auto parts = gs_archive_parts(base);
  if (parts.empty())
    throw std::runtime_error(
        "gs_archive_fingerprint: no archive part files found for base: " +
        base);
  GsFingerprint fp;
  fp.nparts = static_cast<int>(parts.size());
  std::uint64_t h = kFnv1a64Basis;
  std::vector<char> buf(1 << 20);
  for (const auto &p : parts) {
    std::ifstream in(p, std::ios::binary);
    if (!in)
      throw std::runtime_error("gs_archive_fingerprint: cannot read: " + p);
    while (in) {
      in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
      const auto got = static_cast<std::size_t>(in.gcount());
      if (got == 0) break;
      h = fnv1a64_update(h, buf.data(), got);
      fp.bytes += got;
    }
  }
  char hex[17];
  std::snprintf(hex, sizeof(hex), "%016llx",
                static_cast<unsigned long long>(h));
  fp.hex = hex;
  return fp;
}

/// True if the metadata already holds any response state entries — i.e.
/// restart data exists that a mismatched ground state could poison.
inline bool metadata_has_response_states(const nlohmann::json &j) {
  for (const char *k : {"fd_states", "excited_states", "vbc_states"})
    if (j.contains(k) && j[k].is_object() && !j[k].empty()) return true;
  return false;
}

enum class GsGateVerdict {
  FreshDir,     // no stamp, no states -> stamp and proceed
  Match,        // stamp equals current archive -> proceed
  MissingStamp, // states exist but predate the stamp -> warn, stamp, proceed
  Mismatch      // stamp differs -> the restart state belongs to another GS
};

inline GsGateVerdict gs_fingerprint_verdict(const nlohmann::json &meta,
                                            const std::string &current_hex) {
  std::string stored;
  if (meta.contains("ground_state") && meta["ground_state"].is_object())
    stored = meta["ground_state"].value("fnv1a64", "");
  if (stored.empty())
    return metadata_has_response_states(meta) ? GsGateVerdict::MissingStamp
                                              : GsGateVerdict::FreshDir;
  return stored == current_hex ? GsGateVerdict::Match : GsGateVerdict::Mismatch;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_GS_FINGERPRINT_HPP
