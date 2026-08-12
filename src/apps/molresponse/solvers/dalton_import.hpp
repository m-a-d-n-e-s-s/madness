#ifndef MOLRESPONSE_V3_SOLVERS_DALTON_IMPORT_HPP
#define MOLRESPONSE_V3_SOLVERS_DALTON_IMPORT_HPP

// =========================================================================
// dalton_import.hpp — the `dalton.dir` import contract (property-showcase W3).
//
// IMPORT-ONLY: madness never invokes DALTON. Given a directory produced by a
// DALTON linear-response run (*QUADRA/QRLRVE — e.g. gecko's alpha/beta decks),
// this header
//
//   1. LOCATES the artifacts by convention (explicit overrides allowed):
//        molden : <dir>/molden.inp, else a unique *.molden.inp / *.molden
//        RSPVEC : <dir>/RSPVEC,     else a unique *.RSPVEC
//        .out   : a unique <dir>/*.out (optional — provenance + alpha check)
//      If molden/RSPVEC are not loose in <dir>, a unique <dir>/*.tar.gz (the
//      DALTON output archive) is extracted — members RSPVEC + molden.inp only,
//      via a rank-0 shell-out to `tar` — into <calc_dir>/dalton_import/.
//      Multiple candidates for any slot are a HARD ERROR naming them (pass the
//      explicit override instead); silent guessing is forbidden.
//
//   2. FINGERPRINTS before any consumption:
//        geometry : the molden [Atoms] block vs the active Molecule — HARD
//                   ERROR beyond tol (both geometries + max deviation printed).
//        method   : from the .out ("Wave function type --- HF ---"); anything
//                   non-HF is a HARD ERROR (DFT seeding is not validated).
//        basis    : from the .out ('Basis set used is "..."'), recorded.
//      Provenance {dalton_dir, basis, method, geometry_hash} plus the MADNESS
//      build that wrote the MRA seed artifacts (version/commit — the archive
//      format is build-locked) goes into response_metadata.json via the
//      metadata layer (ResponseMetadata::set_seeded_from).
//
//   3. SEEDS the FD states: for each planned (dipole perturbation, freq) the
//      matching RSPVEC record is located at the EXACT frequency (1e-9 au) —
//      frequencies are exact 1-to-1 BY CONSTRUCTION (both inputs come from one
//      description); a mismatch means the description was violated, so it is a
//      HARD ERROR listing requested vs available, never a nearest-frequency
//      fallback. The record's (Z, Y) blocks are projected to MRA with the
//      PINNED convention (validated on quad_H2O aug-cc-pVTZ, 7 digits against
//      DALTON's own QRLRVE <<A;A>> at every frequency and axis):
//
//          x_i = -sum_a Z_ai phi_a        y_i = +sum_a Y_ai phi_a
//
//      so that alpha_ab = -2 [<x^a|v_b phi> + <y^a|v_b phi>] (the validated
//      molresponse convention, calc_executor assemble_alpha) reproduces
//      DALTON's printed value. Q-projected against the ACTIVE MADNESS ground
//      state, then written as an FD bundle under the synthetic protocol key
//      "<active_key>_dseed" — registered in the protocols table at the active
//      (thresh, k), so the UNCHANGED restart machinery (best_usable_fd_source_
//      key / try_load_fd_state) picks it up as a coarser-or-equal source:
//      reconcile says Restart, the loader returns exact=false, and the seed is
//      an initial guess to refine — never mistaken for a solved rung.
//
// Collective discipline: run_dalton_import() must be called by every rank.
// Filesystem writes (extraction, metadata) are rank-0; parses are replicated
// (read-only); every throw is broadcast-then-throw so no rank diverges before
// a collective.
// =========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../ResponsePropertyPlanner.hpp"
#include "../tools/dalton_gto.hpp"
#include "../tools/dalton_mra.hpp"
#include "../tools/dalton_rspvec.hpp"
#include "fd_save_load.hpp"        // response_filename
#include "gs_fingerprint.hpp"      // fnv1a64_update / kFnv1a64Basis
#include "response_metadata.hpp"
#include "response_state.hpp"

#include <nlohmann/json.hpp>
#include <madness/chem/molecule.h>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor_lapack.h>   // syev (Loewdin gauge polish)
#include <madness/world/MADworld.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace molresponse_v3 {

// -------------------------------------------------------------------------
// Manifest — everything the import knows about the DALTON directory.
// -------------------------------------------------------------------------
struct DaltonManifest {
  std::string dir;           // the dalton.dir as given
  std::string molden_path;   // resolved molden file
  std::string rspvec_path;   // resolved RSPVEC binary
  std::string out_path;      // resolved .out ("" = none found; provenance-degraded)
  std::string basis  = "unknown";   // from the .out
  std::string method = "unknown";   // "HF" | "unknown" (non-HF is rejected earlier)
  std::string geometry_hash;        // FNV-1a-64 over the canonical geometry string
  std::vector<int> atomic_numbers;                // molden [Atoms]
  std::vector<std::array<double, 3>> coords;      // bohr

  /// The provenance block stamped into response_metadata.json (doc: showcase
  /// interface contract `seeded_from: {dalton_dir, basis, method,
  /// geometry_hash}` + the MADNESS build identity of the MRA artifacts).
  nlohmann::json provenance() const {
    return nlohmann::json{
        {"dalton_dir",    dir},
        {"molden",        molden_path},
        {"rspvec",        rspvec_path},
        {"out",           out_path},
        {"basis",         basis},
        {"method",        method},
        {"geometry_hash", geometry_hash},
        // Which MADNESS build wrote the MRA seed artifacts. Mandatory since
        // the archive format is build-locked ("incompatible MADNESS version"
        // on load) — without this a rejected seed dir is undiagnosable.
        {"mra_writer", nlohmann::json{
            {"version",    madness::info::version()},
            {"git_commit", madness::info::git_commit()},
        }},
    };
  }
};

namespace detail_dalton_import {

namespace fs = std::filesystem;

/// Candidates in `dir` whose filename matches: exact `name`, or suffix `ext`.
inline std::vector<std::string> find_candidates(const std::string &dir,
                                                const std::string &exact,
                                                const std::string &suffix) {
  std::vector<std::string> hits;
  std::error_code ec;
  // Exact conventional name wins outright.
  if (!exact.empty() && fs::exists(dir + "/" + exact, ec))
    return {dir + "/" + exact};
  for (const auto &e : fs::directory_iterator(dir, ec)) {
    if (ec || !e.is_regular_file()) continue;
    const std::string name = e.path().filename().string();
    if (!suffix.empty() && name.size() > suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
      hits.push_back(e.path().string());
  }
  std::sort(hits.begin(), hits.end());
  return hits;
}

inline std::string list_join(const std::vector<std::string> &v) {
  std::string s;
  for (const auto &x : v) { if (!s.empty()) s += ", "; s += x; }
  return s;
}

/// Resolve one artifact slot: explicit override > exact conventional name >
/// unique suffix match. Zero or multiple suffix matches -> "" (caller decides
/// whether that is fatal); an override that does not exist IS fatal here.
inline std::string resolve_slot(const std::string &dir,
                                const std::string &override_path,
                                const std::string &exact,
                                const std::string &suffix,
                                std::string &err, const char *what) {
  if (!override_path.empty()) {
    if (!fs::exists(override_path)) {
      err = std::string("dalton import: explicit ") + what + " override '" +
            override_path + "' does not exist";
    }
    return override_path;
  }
  auto c = find_candidates(dir, exact, suffix);
  if (c.size() == 1) return c[0];
  if (c.size() > 1) {
    err = std::string("dalton import: multiple ") + what + " candidates in " +
          dir + " (" + list_join(c) + ") — pass the explicit override";
  }
  return {};
}

} // namespace detail_dalton_import

// -------------------------------------------------------------------------
// locate — by-convention artifact discovery (+ tarball extraction).
// Pure rank-0 function: performs filesystem reads and (only when molden/
// RSPVEC are not loose and a unique *.tar.gz exists) an extraction into
// `extract_dir` via a shell-out to `tar`. Throws on any ambiguity.
// -------------------------------------------------------------------------
inline DaltonManifest locate_dalton_dir(const std::string &dir,
                                        const std::string &extract_dir,
                                        const std::string &molden_override = {},
                                        const std::string &rspvec_override = {},
                                        const std::string &out_override = {}) {
  namespace fs = std::filesystem;
  using namespace detail_dalton_import;

  if (!fs::exists(dir) || !fs::is_directory(dir))
    throw std::runtime_error("dalton import: dalton.dir '" + dir +
                             "' is not a directory");

  DaltonManifest m;
  m.dir = dir;
  std::string err;

  m.molden_path = resolve_slot(dir, molden_override, "molden.inp",
                               ".molden.inp", err, "molden");
  if (err.empty() && m.molden_path.empty())
    m.molden_path = resolve_slot(dir, "", "", ".molden", err, "molden");
  if (!err.empty()) throw std::runtime_error(err);

  m.rspvec_path = resolve_slot(dir, rspvec_override, "RSPVEC",
                               ".RSPVEC", err, "RSPVEC");
  if (!err.empty()) throw std::runtime_error(err);

  // Not loose -> try the DALTON output tarball (unique *.tar.gz). Extraction
  // is documented behavior: `tar -xzf <tar> -C <extract_dir> RSPVEC molden.inp`
  // (only the two members), so pre-extraction is never required but always
  // honored. tar's exit code is ignored on purpose (one member may be absent);
  // what matters is which files exist afterwards.
  if (m.molden_path.empty() || m.rspvec_path.empty()) {
    auto tars = find_candidates(dir, "", ".tar.gz");
    if (tars.size() > 1)
      throw std::runtime_error(
          "dalton import: molden/RSPVEC not loose in " + dir +
          " and multiple *.tar.gz candidates (" + list_join(tars) +
          ") — extract manually or pass explicit overrides");
    if (tars.size() == 1) {
      fs::create_directories(extract_dir);
      const std::string cmd = "tar -xzf '" + tars[0] + "' -C '" + extract_dir +
                              "' RSPVEC molden.inp 2>/dev/null";
      (void)std::system(cmd.c_str());
      if (m.molden_path.empty() && fs::exists(extract_dir + "/molden.inp"))
        m.molden_path = extract_dir + "/molden.inp";
      if (m.rspvec_path.empty() && fs::exists(extract_dir + "/RSPVEC"))
        m.rspvec_path = extract_dir + "/RSPVEC";
    }
  }
  if (m.molden_path.empty())
    throw std::runtime_error(
        "dalton import: no molden file in " + dir +
        " (looked for molden.inp, *.molden.inp, *.molden, and inside a unique "
        "*.tar.gz) — pass --dalton-molden=PATH");
  if (m.rspvec_path.empty())
    throw std::runtime_error(
        "dalton import: no RSPVEC in " + dir +
        " (looked for RSPVEC, *.RSPVEC, and inside a unique *.tar.gz) — pass "
        "--dalton-rspvec=PATH");

  {
    auto outs = find_candidates(dir, "", ".out");
    if (!out_override.empty()) {
      if (!fs::exists(out_override))
        throw std::runtime_error("dalton import: explicit out override '" +
                                 out_override + "' does not exist");
      m.out_path = out_override;
    } else if (outs.size() == 1) {
      m.out_path = outs[0];
    } else if (outs.size() > 1) {
      throw std::runtime_error(
          "dalton import: multiple *.out candidates in " + dir + " (" +
          list_join(outs) + ") — pass --dalton-out=PATH");
    }
    // zero .out files: allowed (provenance-degraded; caller warns).
  }

  // Geometry from the molden (authoritative: it is what the MOs live on).
  {
    DaltonMoldenResult mol = read_molden(m.molden_path);
    m.atomic_numbers = mol.atomic_numbers;
    m.coords         = mol.coords;
    std::uint64_t h = kFnv1a64Basis;
    char line[128];
    for (size_t i = 0; i < m.coords.size(); ++i) {
      std::snprintf(line, sizeof line, "%d %.6f %.6f %.6f\n",
                    m.atomic_numbers[i], m.coords[i][0], m.coords[i][1],
                    m.coords[i][2]);
      h = fnv1a64_update(h, line, std::strlen(line));
    }
    char hex[17];
    std::snprintf(hex, sizeof hex, "%016llx",
                  static_cast<unsigned long long>(h));
    m.geometry_hash = hex;
  }

  // Method + basis from the .out (when present).
  if (!m.out_path.empty()) {
    std::ifstream in(m.out_path);
    std::string line;
    const std::regex basis_re("Basis set used is \"([^\"]+)\"");
    const std::regex wf_re("Wave function type\\s+-+\\s*(.+?)\\s*-+\\s*$");
    std::smatch sm;
    while (std::getline(in, line)) {
      if (m.basis == "unknown" && std::regex_search(line, sm, basis_re))
        m.basis = sm[1];
      if (m.method == "unknown" && std::regex_search(line, sm, wf_re))
        m.method = sm[1];
      if (m.basis != "unknown" && m.method != "unknown") break;
    }
  }
  return m;
}

// -------------------------------------------------------------------------
// Geometry fingerprint — DALTON geometry vs the active Molecule.
// -------------------------------------------------------------------------
struct DaltonGeometryCheck {
  bool        ok = false;
  double      max_dev = 0.0;   // bohr
  std::string report;          // both geometries side by side + verdict
};

inline DaltonGeometryCheck
fingerprint_dalton_geometry(const DaltonManifest &m,
                            const madness::Molecule &mol,
                            double tol_bohr = 1e-4) {
  DaltonGeometryCheck c;
  std::ostringstream r;
  r << "dalton import geometry fingerprint (bohr; tol=" << tol_bohr << "):\n";
  r << "  DALTON (" << m.molden_path << ")  vs  active Molecule\n";
  const size_t nd = m.coords.size(), nm = mol.natom();
  if (nd != nm) {
    r << "  ATOM-COUNT MISMATCH: dalton natom=" << nd
      << "  molecule natom=" << nm << "\n";
    for (size_t i = 0; i < nd; ++i)
      r << "    dalton  Z=" << m.atomic_numbers[i] << "  " << m.coords[i][0]
        << " " << m.coords[i][1] << " " << m.coords[i][2] << "\n";
    for (size_t i = 0; i < nm; ++i) {
      const auto &a = mol.get_atom(static_cast<unsigned int>(i));
      r << "    active  Z=" << a.atomic_number << "  " << a.x << " " << a.y
        << " " << a.z << "\n";
    }
    c.ok = false;
    c.max_dev = std::numeric_limits<double>::infinity();
    c.report = r.str();
    return c;
  }
  bool z_ok = true;
  for (size_t i = 0; i < nd; ++i) {
    const auto &a = mol.get_atom(static_cast<unsigned int>(i));
    const double dx = m.coords[i][0] - a.x;
    const double dy = m.coords[i][1] - a.y;
    const double dz = m.coords[i][2] - a.z;
    const double dev = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
    c.max_dev = std::max(c.max_dev, dev);
    const bool zmatch =
        (m.atomic_numbers[i] == static_cast<int>(a.atomic_number));
    z_ok = z_ok && zmatch;
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "  atom %2zu  Z %2d|%2u  dalton % .8f % .8f % .8f   active "
                  "% .8f % .8f % .8f   |dev| %.3e%s\n",
                  i, m.atomic_numbers[i], a.atomic_number, m.coords[i][0],
                  m.coords[i][1], m.coords[i][2], a.x, a.y, a.z, dev,
                  zmatch ? "" : "  <-- ELEMENT MISMATCH");
    r << buf;
  }
  c.ok = z_ok && (c.max_dev <= tol_bohr);
  r << "  max deviation = " << c.max_dev << " bohr  ->  "
    << (c.ok ? "MATCH" : "MISMATCH") << "\n";
  c.report = r.str();
  return c;
}

// -------------------------------------------------------------------------
// Exact 1-to-1 frequency matching (pure; unit-testable without a World).
// -------------------------------------------------------------------------

/// |dw| tolerance for "the same frequency". Guards only against decimal->
/// binary representation noise — NEVER a nearest-frequency fallback (the
/// physically closest distinct grid points are >= 1e-5 apart by the freq_key
/// f%.5f contract, 4 orders above this).
inline constexpr double kDaltonFreqTol = 1e-9;

inline const char *dalton_dipole_label(int axis) {
  static const char *labs[3] = {"XDIPLEN", "YDIPLEN", "ZDIPLEN"};
  return labs[axis];
}

/// Find the RSPVEC record for (dipole axis, freq) at the EXACT frequency.
inline const RspVecEntry *
find_dalton_dipole_entry(const std::vector<RspVecEntry> &entries, int axis,
                         double freq) {
  for (const auto &e : entries)
    if (std::string(e.lab1) == dalton_dipole_label(axis) &&
        std::abs(e.freq1 - freq) <= kDaltonFreqTol)
      return &e;
  return nullptr;
}

/// Verify every requested (axis, freq) has an exact RSPVEC record. Returns ""
/// when all are present, else the complete hard-error message listing the
/// missing requests AND every available dipole frequency (per label).
inline std::string
match_dalton_frequencies(const std::vector<RspVecEntry> &entries,
                         const std::vector<std::pair<int, double>> &wanted,
                         const std::string &rspvec_path) {
  std::vector<std::string> missing;
  for (const auto &[axis, freq] : wanted)
    if (!find_dalton_dipole_entry(entries, axis, freq)) {
      char b[64];
      std::snprintf(b, sizeof b, "%s @ %.9f", dalton_dipole_label(axis), freq);
      missing.push_back(b);
    }
  if (missing.empty()) return {};

  std::map<std::string, std::vector<double>> avail;
  for (const auto &e : entries) {
    const std::string l(e.lab1);
    if (l == "XDIPLEN" || l == "YDIPLEN" || l == "ZDIPLEN")
      avail[l].push_back(e.freq1);
  }
  std::ostringstream os;
  os << "dalton import: requested frequencies ABSENT from " << rspvec_path
     << " (exact 1-to-1 contract, no nearest-frequency fallback):\n";
  os << "  missing:\n";
  for (const auto &s : missing) os << "    " << s << "\n";
  os << "  available:\n";
  for (const auto &[lab, fr] : avail) {
    os << "    " << lab << " :";
    char b[32];
    for (double f : fr) { std::snprintf(b, sizeof b, " %.9f", f); os << b; }
    os << "\n";
  }
  os << "  (both inputs are generated from one description — a mismatch "
        "means the description was violated)";
  return os.str();
}

// -------------------------------------------------------------------------
// FD seeding.
// -------------------------------------------------------------------------
struct DaltonSeedReport {
  int n_seeded  = 0;
  int n_skipped = 0;   // usable on-disk state already present
};

/// The synthetic protocol key the seed bundles are written under. Registered
/// in the protocols table at the ACTIVE (thresh, k) so best_usable_fd_source_
/// key treats it as coarser-or-equal to every rung, but never equal to a real
/// rung key -> the loader always returns exact=false ("initial guess, refine
/// me"), and rung entries never overwrite the seed record (provenance stays).
inline std::string dalton_seed_protocol_key() {
  return protocol_key() + "_dseed";
}

/// Seed the planned dipole FD states from the DALTON RSPVEC. Collective.
///
/// PRECONDITIONS: set_response_protocol(...) + gs.prepare(...) have run at the
/// ramp's FIRST rung (the driver does this before planning) — the seed is
/// projected at the active (k, thresh, cell) and Q-projected against gs.
///
/// Frequencies are matched EXACTLY (|dw| <= 1e-9 au). Any planned dipole
/// frequency absent from the RSPVEC is a hard error listing requested vs
/// available. Non-dipole requests (nuclear FD) have no DALTON counterpart and
/// are left to the normal zero-guess path.
inline DaltonSeedReport
seed_fd_from_dalton(madness::World &world, GroundState &gs,
                    const std::vector<FDRequest> &fd_requests,
                    const std::string &calc_dir, const DaltonManifest &m) {
  using namespace madness;
  DaltonSeedReport rep;

  if (!gs.is_spin_restricted())
    throw std::runtime_error(
        "dalton import: FD seeding is ClosedShell-only (RSPVEC import has no "
        "open-shell mapping)");

  // The distinct dipole (axis, freq) targets, plan order preserved.
  struct Item { int axis; double freq; };
  std::vector<Item> items;
  {
    std::set<std::pair<int, std::string>> seen;
    for (const auto &r : fd_requests) {
      if (r.pert.kind != Perturbation::Kind::Dipole) continue;
      auto key = std::make_pair(r.pert.axis, ResponseMetadata::freq_key(r.freq));
      if (seen.insert(key).second) items.push_back({r.pert.axis, r.freq});
    }
  }
  if (items.empty()) {
    if (world.rank() == 0)
      print("[DALTON-SEED] no dipole FD requests in the plan — nothing to seed");
    return rep;
  }

  // Replicated read-only parses (files exist on every rank after the caller's
  // locate + fence; identical bytes -> identical state on every rank).
  auto rsp = read_rspvec(m.rspvec_path);
  const auto &info    = rsp.first;
  const auto &entries = rsp.second;
  DaltonMoldenResult molden = read_molden(m.molden_path);

  const int n_ao  = molden.n_ao;
  const int n_mo  = molden.n_mo;
  const int n_occ = info.nish[0];
  const int n_vir = n_mo - n_occ;
  if (info.nsym != 1)
    throw std::runtime_error(
        "dalton import: RSPVEC has NSYM=" + std::to_string(info.nsym) +
        " — only C1 (Nosymmetry) DALTON runs are supported");
  if (static_cast<size_t>(n_occ) != gs.orbitals_alpha().size())
    throw std::runtime_error(
        "dalton import: DALTON n_occ=" + std::to_string(n_occ) +
        " does not match the ground state's " +
        std::to_string(gs.orbitals_alpha().size()) +
        " occupied orbitals — different molecule/charge?");

  // ---- exact 1-to-1 frequency match (no nearest fallback) -----------------
  {
    std::vector<std::pair<int, double>> wanted;
    for (const auto &it : items) wanted.emplace_back(it.axis, it.freq);
    const std::string err =
        match_dalton_frequencies(entries, wanted, m.rspvec_path);
    if (!err.empty()) throw std::runtime_error(err);
  }

  // ---- skip decisions (rank 0 reads metadata; broadcast) ------------------
  // If a usable on-disk source already exists for (pert, freq) at the active
  // rung, the run would restart from IT — re-projecting a DALTON seed under
  // it would be wasted work, so skip (logged). A diverged-only dir has no
  // usable source and is re-seeded.
  const double active_thresh = FunctionDefaults<3>::get_thresh();
  const int    active_k      = FunctionDefaults<3>::get_k();
  const std::string active_key = protocol_key();
  const std::string seed_key   = dalton_seed_protocol_key();

  std::vector<int> skip(items.size(), 0);
  if (world.rank() == 0) {
    const std::string meta_path = calc_dir + "/response_metadata.json";
    if (std::filesystem::exists(meta_path)) {
      auto meta = ResponseMetadata::load_or_create(meta_path);
      for (size_t i = 0; i < items.size(); ++i) {
        const std::string pdesc =
            Perturbation::dipole(items[i].axis).description();
        const std::string fkey = ResponseMetadata::freq_key(items[i].freq);
        if (!best_usable_fd_source_key(meta.json(), pdesc, fkey, active_thresh,
                                       active_k, active_key)
                 .empty())
          skip[i] = 1;
      }
    }
  }
  world.gop.broadcast_serializable(skip, 0);

  if (world.rank() == 0) {
    print("[DALTON-SEED] rspvec=", m.rspvec_path);
    print("[DALTON-SEED] molden=", m.molden_path, "  n_ao=", n_ao,
          " n_mo=", n_mo, " n_occ=", n_occ, " n_vir=", n_vir);
    print("[DALTON-SEED] basis=", m.basis, "  method=", m.method,
          "  geometry_hash=", m.geometry_hash);
    print("[DALTON-SEED] seed protocol key=", seed_key,
          "  (active=", active_key, ")");
  }

  // ---- occupied-orbital GAUGE ROTATION -------------------------------------
  // DALTON's response vectors are indexed by its CANONICAL occupied MOs; the
  // FD solver's x_i belong to the MADNESS ground-state orbitals — typically
  // LOCALIZED (moldft `localize new`). Response functions transform
  // covariantly with the occupied orbitals, so without this the per-orbital
  // pairing is scrambled and the seed is worthless (observed: seed-implied
  // alpha_zz 3.56 vs 8.39 on the h2o A/B). Build the occupied-occupied
  // overlap M(i,j) = <phi^DAL_i | phi^MAD_j> in MRA space, Loewdin-polish it
  // to an exact unitary U = M (M^T M)^{-1/2}, and rotate every seed block:
  //     x^MAD_j = sum_i x^DAL_i U(i,j)     (same for y).
  // The eigenvalues of M^T M are ALSO the quantitative ground-state import
  // fidelity: ~1 when the two occupied spaces coincide. An eigenvalue below
  // 0.5 means a MADNESS occupied orbital has little support in the DALTON
  // occupied space — a different state entirely — hard error.
  Tensor<double> U;
  {
    bool any = false;
    for (size_t i = 0; i < items.size(); ++i) any = any || !skip[i];
    if (any) {
      vector_real_function_3d phi_dal;
      for (int i = 0; i < n_occ; ++i) {
        std::vector<double> w(
            molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i) * n_ao,
            molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i + 1) * n_ao);
        phi_dal.push_back(project_dalton_weights(world, molden.basis,
                                                 std::move(w), active_thresh));
      }
      truncate(world, phi_dal, active_thresh);
      Tensor<double> M = matrix_inner(world, phi_dal, gs.orbitals_alpha());
      Tensor<double> S = inner(transpose(M), M);   // M^T M, n_occ x n_occ
      Tensor<double> V, ev;
      syev(S, V, ev);
      double ev_min = ev(0L);
      for (long k = 0; k < ev.dim(0); ++k) ev_min = std::min(ev_min, ev(k));
      if (world.rank() == 0)
        print("[DALTON-SEED] occupied-gauge overlap M^T M eigenvalues:", ev,
              "  (import fidelity; 1 = perfect span match)");
      if (ev_min < 0.5) {
        std::ostringstream os;
        os << "dalton import: occupied-space mismatch — min eigenvalue of the "
              "DALTON/MADNESS occupied overlap M^T M is " << ev_min
           << " (want ~1). The DALTON ground state does not span the MADNESS "
              "occupied orbitals (different state/charge/geometry?)";
        throw std::runtime_error(os.str());
      }
      Tensor<double> Sinvhalf(static_cast<long>(n_occ),
                              static_cast<long>(n_occ));
      for (long a = 0; a < n_occ; ++a)
        for (long b = 0; b < n_occ; ++b) {
          double s = 0.0;
          for (long k = 0; k < n_occ; ++k)
            s += V(a, k) * V(b, k) / std::sqrt(ev(k));
          Sinvhalf(a, b) = s;
        }
      U = inner(M, Sinvhalf);   // Loewdin: exactly unitary
    }
  }

  // DALTON QRLRVE <<A;A>>(w) diagonal values from the .out (when present) —
  // diagnostic comparison only, never a gate (the seed's implied alpha uses
  // the MRA ground state, so basis-quality differences are EXPECTED — they
  // are the point of refining the seed).
  std::map<std::pair<int, long>, double> dal_alpha;   // (axis, round(w*1e5))
  if (world.rank() == 0 && !m.out_path.empty()) {
    std::ifstream in(m.out_path);
    std::string line;
    const std::regex re("<< ([XYZ])DIPLEN\\s*;\\s*([XYZ])DIPLEN\\s*>>\\s*\\(\\s*"
                        "([-0-9.]+)\\)\\s*:\\s*([-0-9.eEdD+]+)");
    std::smatch sm;
    while (std::getline(in, line)) {
      if (!std::regex_search(line, sm, re)) continue;
      if (sm[1] != sm[2]) continue;
      const int axis = std::string("XYZ").find(sm[1].str()[0]);
      std::string val = sm[4];
      for (auto &ch : val) if (ch == 'D' || ch == 'd') ch = 'e';
      dal_alpha[{axis, std::lround(std::stod(sm[3]) * 1e5)}] = std::stod(val);
    }
  }

  // ---- project + save ------------------------------------------------------
  // Per-axis Q(mu*phi) for the implied-alpha diagnostic, built lazily.
  std::map<int, vector_real_function_3d> vdiag;
  auto v_of = [&](int axis) -> const vector_real_function_3d & {
    auto it = vdiag.find(axis);
    if (it == vdiag.end())
      it = vdiag.emplace(axis, dipole_perturbation(world, gs, axis)).first;
    return it->second;
  };

  bool wrote_any = false;
  for (size_t i = 0; i < items.size(); ++i) {
    const int    axis = items[i].axis;
    const double freq = items[i].freq;
    const Perturbation pert = Perturbation::dipole(axis);
    const std::string  pdesc = pert.description();
    const std::string  fkey  = ResponseMetadata::freq_key(freq);
    if (skip[i]) {
      ++rep.n_skipped;
      if (world.rank() == 0)
        print("[DALTON-SEED] ", pdesc, "@", freq,
              " — usable on-disk state already present; NOT overwriting");
      continue;
    }
    const RspVecEntry *e = find_dalton_dipole_entry(entries, axis, freq);   // guaranteed above
    auto [Z, Y] = split_ov(e->vec, n_occ, n_vir);

    // Pinned convention: x = -Z·phi_vir, y = +Y·phi_vir (header block).
    auto x = project_dalton_ov_block(world, molden.basis, molden.mo_coeffs,
                                     n_ao, n_mo, n_occ, n_vir, Z,
                                     active_thresh, -1.0);
    vector_real_function_3d y;
    const bool is_static = std::abs(freq) < 1e-12;
    if (!is_static) {
      if (!Y.empty()) {
        y = project_dalton_ov_block(world, molden.basis, molden.mo_coeffs,
                                    n_ao, n_mo, n_occ, n_vir, Y,
                                    active_thresh, +1.0);
      } else {  // TDA-style record (no Y block): y = 0
        y = madness::copy(world, x);
        madness::scale(world, y, 0.0);
      }
    }

    // Rotate into the MADNESS occupied gauge (see the U block above), then
    // Q-project against the ACTIVE MADNESS ground state (removes occupied-
    // space contamination from the finite-basis <-> MRA orbital mismatch).
    x = transform(world, x, U);
    x = gs.Q()(x);
    truncate(world, x, active_thresh);
    if (!is_static) {
      y = transform(world, y, U);
      y = gs.Q()(y);
      truncate(world, y, active_thresh);
    }

    // Implied-alpha diagnostic against DALTON's own QRLRVE value.
    {
      const auto &v = v_of(axis);
      const double ax_val =
          is_static ? -4.0 * inner(x, v) : -2.0 * (inner(x, v) + inner(y, v));
      double dal = std::numeric_limits<double>::quiet_NaN();
      if (world.rank() == 0) {
        auto it = dal_alpha.find({axis, std::lround(freq * 1e5)});
        if (it != dal_alpha.end()) dal = it->second;
        const char ax_c = "xyz"[axis];
        char msg[160];
        if (dal == dal)
          std::snprintf(msg, sizeof msg,
                        "seed-implied alpha_%c%c = %.6f  (DALTON <<A;A>> = %.6f)",
                        ax_c, ax_c, ax_val, dal);
        else
          std::snprintf(msg, sizeof msg,
                        "seed-implied alpha_%c%c = %.6f  (no DALTON .out value)",
                        ax_c, ax_c, ax_val);
        print("[DALTON-SEED] ", pdesc, "@", freq, "  ", msg);
      }
    }

    // Save the bundle under the seed key; metadata via the layer (rank 0).
    const std::string archive_basename = response_filename(pdesc, seed_key, freq);
    const std::string archive_path     = calc_dir + "/" + archive_basename;
    if (world.rank() == 0) std::filesystem::create_directories(calc_dir);
    world.gop.fence();

    IoBackend backend;
    const char *type_tag;
    if (is_static) {
      ResponseStateX<ClosedShell> st;
      st.x_alpha = std::move(x);
      backend = st.save(world, archive_path);
      type_tag = "static";
    } else {
      ResponseStateXY<ClosedShell> st;
      st.x_alpha = std::move(x);
      st.y_alpha = std::move(y);
      backend = st.save(world, archive_path);
      type_tag = "full";
    }

    if (world.rank() == 0) {
      auto meta = ResponseMetadata::load_or_create(
          calc_dir + "/response_metadata.json");
      if (!meta.json()["protocols"].contains(seed_key))
        meta.set_protocol(seed_key, active_thresh, active_k, /*index=*/-1);
      nlohmann::json entry = {
          {"freq",         freq},
          {"type",         type_tag},
          {"shell",        "closed_shell"},
          {"converged",    false},   // a seed is an initial guess, never a result
          {"accepted",     false},
          {"diverged",     false},
          {"iter",         0},
          {"bsh_residual", 1.0},     // unknown; the first solve step measures it
          {"seed",         "dalton_import"},
          {"archive",      archive_basename},
          {"backend",      io_backend_tag(backend)},
          {"writer_nproc", static_cast<int>(world.size())},
          {"dalton_freq",  e->freq1},
      };
      meta.set_fd_state(pdesc, seed_key, fkey, entry);
      meta.set_seeded_from(m.provenance());
      meta.save();
      print("[DALTON-SEED] wrote ", pdesc, "@", freq, "  type=", type_tag,
            "  archive=", archive_basename);
    }
    world.gop.fence();
    ++rep.n_seeded;
    wrote_any = true;
  }
  (void)wrote_any;

  if (world.rank() == 0)
    print("[DALTON-SEED] done: seeded=", rep.n_seeded,
          "  skipped(existing)=", rep.n_skipped);
  world.gop.fence();
  return rep;
}

// -------------------------------------------------------------------------
// run_dalton_import — the one-call driver/adapter entry point. Collective.
// locate (rank 0, broadcast) -> geometry fingerprint (HARD error) -> method
// gate (non-HF = HARD error) -> seed the FD plan.
// -------------------------------------------------------------------------
inline DaltonSeedReport
run_dalton_import(madness::World &world, GroundState &gs,
                  const madness::Molecule &molecule,
                  const ResponsePlan &plan, const std::string &calc_dir,
                  const std::string &dalton_dir,
                  const std::string &molden_override = {},
                  const std::string &rspvec_override = {},
                  const std::string &out_override = {},
                  double geometry_tol_bohr = 1e-4) {
  DaltonManifest m;
  std::string err;

  // Rank 0 locates (may extract) + fingerprints; everything is broadcast so
  // every rank throws the SAME error or proceeds with the SAME manifest.
  std::string report;
  int ok = 0;
  if (world.rank() == 0) {
    try {
      m = locate_dalton_dir(dalton_dir, calc_dir + "/dalton_import",
                            molden_override, rspvec_override, out_override);
      if (molecule.natom() == 0)
        throw std::runtime_error(
            "dalton import: active Molecule is EMPTY — cannot fingerprint the "
            "geometry (moldft calc_info.json missing next to the archive?)");
      auto check = fingerprint_dalton_geometry(m, molecule, geometry_tol_bohr);
      report = check.report;
      if (!check.ok)
        throw std::runtime_error(
            "dalton import: GEOMETRY FINGERPRINT MISMATCH — the DALTON "
            "directory belongs to a different geometry.\n" + check.report);
      if (m.method != "HF") {
        if (m.method == "unknown" && m.out_path.empty()) {
          // no .out at all: method unverifiable — warn, proceed (recorded).
          std::fprintf(stderr,
                       "[DALTON-SEED] WARNING: no .out in %s — method/basis "
                       "unverified (recorded as 'unknown')\n",
                       dalton_dir.c_str());
        } else {
          throw std::runtime_error(
              "dalton import: DALTON wave function type is '" + m.method +
              "' (from " + m.out_path +
              ") — only HF seeds are supported (DFT mapping unvalidated)");
        }
      }
      ok = 1;
    } catch (const std::exception &ex) {
      err = ex.what();
    }
  }
  world.gop.broadcast_serializable(err, 0);
  if (!err.empty()) throw std::runtime_error(err);   // collective throw
  world.gop.broadcast(ok, 0);

  // Broadcast the manifest fields the other ranks need.
  world.gop.broadcast_serializable(m.dir, 0);
  world.gop.broadcast_serializable(m.molden_path, 0);
  world.gop.broadcast_serializable(m.rspvec_path, 0);
  world.gop.broadcast_serializable(m.out_path, 0);
  world.gop.broadcast_serializable(m.basis, 0);
  world.gop.broadcast_serializable(m.method, 0);
  world.gop.broadcast_serializable(m.geometry_hash, 0);
  world.gop.fence();   // extraction visible before the replicated parses

  if (world.rank() == 0) {
    print("[DALTON-SEED] fingerprint OK:");
    print(report);
  }

  return seed_fd_from_dalton(world, gs, plan.fd, calc_dir, m);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_DALTON_IMPORT_HPP
