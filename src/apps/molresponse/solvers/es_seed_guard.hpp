#ifndef MOLRESPONSE_V3_SOLVERS_ES_SEED_GUARD_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SEED_GUARD_HPP

// =========================================================================
// Root-identity guard for SEEDED excited-state solves (W5 finding).
//
// A seeded ES solve TRACKS the seeded states: it refines the eigenvectors
// it was handed to convergence, it does NOT search for the N lowest
// states. With a poor seeding basis this is a silent failure mode — the
// solve converges fast, reports converged=true, and quietly delivers a
// DIFFERENT (higher) state in place of a low one the seed basis could not
// describe (H2O/cc-pVDZ pilot: seeded root 3 converged in 2 iters to
// 0.4626 au while every augmented seed finds the true 4th-lowest at
// 0.4096 au — and the seed's own DALTON ladder AGREES with 0.4626, so the
// refinement was faithful; the seed was simply the wrong state).
//
// The guard therefore provides VISIBILITY plus optional cross-checks, not
// an error: after a seeded solve it reports, per converged root,
//
//   seed_overlap  — normalized response-metric overlap of the converged
//                   root with its identity-matched seed root (matched via
//                   stable_index, which travels through slot sorts),
//   omega_seed    — the seed's recorded ω for that root. For a bundle
//                   written by seed_from_dalton this IS the DALTON ladder
//                   value passed via --omegas (saved into roots.json).
//   omega_shift   — ω_converged − ω_seed,
//
// prints a rank-0 summary table, records the block into
// response_metadata.json (excited_states/<key>/seed_guard — via the
// ResponseMetadata layer, never hand-written), and warns loudly when
//   * a converged root has LOW overlap (< kSeedGuardOverlapWarn) with
//     EVERY seed root — the solve left the seeded basin, or
//   * every root kept high overlap (> kSeedGuardTrackingOverlap) and the
//     whole solve took <= n_roots iterations — pure seed tracking; the
//     N-lowest guarantee does NOT hold (informational note), or
//
// Collectivity: evaluate_es_seed_guard computes MRA inner products and is
// collective on `world` (call from every rank). The report it returns is
// rank-uniform (built from replicated tensors); printing and the metadata
// write are rank-0-only inside the print/record helpers.
// =========================================================================

#include "es_root_identity.hpp"           // make_root_id
#include "response_metadata.hpp"          // ResponseMetadata (aggregate writer)
#include "response_state.hpp"             // ResponseStateX/XY storage
#include "../kernels/response_space_ops.hpp"  // rs::metric (RPA-metric overlap)

#include <nlohmann/json.hpp>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Basin-escape threshold: a converged root whose best overlap with EVERY
/// seed root is below this warns loudly.
inline constexpr double kSeedGuardOverlapWarn = 0.5;
/// Pure-tracking threshold: when every root overlaps its seed above this
/// AND the solve took <= n_roots iterations, the tracking note prints.
inline constexpr double kSeedGuardTrackingOverlap = 0.9;

/// Snapshot of the seed bundle a seeded ES solve started from — deep copies
/// of the root vectors (the solver mutates its own in place) plus the seed's
/// recorded ω ladder and stable identity. Captured at load time by
/// capture_es_seed_reference; consumed once after the solve.
template <typename Storage>
struct EsSeedReference {
  std::vector<Storage>    roots;         // deep copies, at the bundle's (k,thresh)
  madness::Tensor<double> omega;         // recorded seed ω (DALTON ladder for
                                         // seed_from_dalton bundles)
  std::vector<int>        stable_index;  // seed slot -> stable identity
  std::string             bundle_dir;    // provenance (basename under calc_dir)
  std::string             source_protocol_key;
};

/// Deep-copy the just-loaded seed state so the guard can compare against it
/// after the solve. Collective (madness::copy of every function).
template <typename StateT>
auto capture_es_seed_reference(madness::World &world, const StateT &s,
                               const std::string &bundle_dir,
                               const std::string &source_protocol_key) {
  using Storage = typename std::decay_t<decltype(s.roots)>::value_type;
  EsSeedReference<Storage> ref;
  ref.roots.reserve(s.roots.size());
  for (const auto &r : s.roots) ref.roots.push_back(r.copy(world));
  ref.omega = madness::copy(s.omega);
  ref.stable_index = s.stable_index;
  // Legacy bundles may lack identity — fall back to slot order, matching
  // what load_es_roots/assign_initial_stable_index would do.
  if (ref.stable_index.size() != ref.roots.size()) {
    ref.stable_index.resize(ref.roots.size());
    for (std::size_t i = 0; i < ref.roots.size(); ++i)
      ref.stable_index[i] = static_cast<int>(i);
  }
  return ref;
}

/// Per-root guard record (one converged slot).
struct EsSeedGuardRoot {
  int         slot          = -1;   // converged-state slot (ascending ω)
  int         stable_index  = -1;
  std::string root_id;
  double      omega         = 0.0;  // converged ω
  int         seed_slot     = -1;   // identity-matched seed slot (-1 = none)
  double      omega_seed    = 0.0;  // seed's recorded ω (valid if seed_slot >= 0)
  double      omega_shift   = 0.0;  // omega - omega_seed (valid if seed_slot >= 0)
  double      seed_overlap  = 0.0;  // vs the identity-matched seed root
  int         best_seed_slot = -1;  // argmax_j overlap
  double      best_overlap  = 0.0;  // max_j overlap (basin-escape discriminator)
};

/// Full guard report — rank-uniform; to_json() feeds the metadata layer.
struct EsSeedGuardReport {
  std::vector<EsSeedGuardRoot> roots;
  bool   seeded    = false;
  int    iters     = 0;      // total iterations of the seeded solve
  int    n_roots   = 0;
  bool   converged = false;  // bundle-level verdict of the solve
  bool   tracking_note = false;   // all seed_overlap > 0.9 && iters <= n_roots
  std::vector<int> basin_escape_slots;   // best_overlap < 0.5
  std::string seed_bundle_dir;
  std::string seed_protocol_key;

  nlohmann::json to_json() const {
    nlohmann::json j;
    j["seeded"]             = seeded;
    j["iters"]              = iters;
    j["n_roots"]            = n_roots;
    j["converged"]          = converged;
    nlohmann::json arr = nlohmann::json::array();
    for (const auto &r : roots) {
      nlohmann::json e;
      e["slot"]         = r.slot;
      e["stable_index"] = r.stable_index;
      e["root_id"]      = r.root_id;
      e["omega"]        = r.omega;
      if (seeded) {
        if (r.seed_slot >= 0) {
          e["seed_slot"]   = r.seed_slot;
          e["omega_seed"]  = r.omega_seed;
          e["omega_shift"] = r.omega_shift;
        }
        e["seed_overlap"]   = r.seed_overlap;
        e["best_seed_slot"] = r.best_seed_slot;
        e["best_overlap"]   = r.best_overlap;
      }
      arr.push_back(e);
    }
    j["roots"] = arr;
    if (seeded) {
      j["seed_bundle_dir"]    = seed_bundle_dir;
      j["seed_protocol_key"]  = seed_protocol_key;
      j["tracking_note"]      = tracking_note;
      j["overlap_warn_threshold"]     = kSeedGuardOverlapWarn;
      j["tracking_overlap_threshold"] = kSeedGuardTrackingOverlap;
      j["basin_escape_slots"] = basin_escape_slots;
    }
    return j;
  }
};

/// Evaluate the guard: overlap matrix (response metric, normalized) between
/// the converged roots and the seed roots, identity-matched via stable_index,
/// plus the tracking / basin-escape verdicts. Collective.
///
/// The seed copies live at the bundle's native (k, thresh); they are
/// re-projected here to the ACTIVE FunctionDefaults so the inner products
/// are well-defined at the final protocol (same path prepare() uses).
template <typename StateT, typename Storage>
EsSeedGuardReport
evaluate_es_seed_guard(madness::World &world, const StateT &sf, bool converged,
                       const EsSeedReference<Storage> &seed) {
  (void)world;  // collectivity contract: the rs:: inner products below are
                // collective on the functions' world; keep the parameter so
                // call sites read as the collective they are.
  EsSeedGuardReport rep;
  rep.seeded            = true;
  rep.iters             = sf.iter;
  rep.n_roots           = static_cast<int>(sf.roots.size());
  rep.converged         = converged;
  rep.seed_bundle_dir   = seed.bundle_dir;
  rep.seed_protocol_key = seed.source_protocol_key;

  const int M = static_cast<int>(sf.roots.size());
  const int Ns = static_cast<int>(seed.roots.size());
  if (M == 0 || Ns == 0) return rep;

  // Re-project the seed copies to the active (k, thresh). Assigning the
  // projection rebinds the local copy's handles only.
  std::vector<Storage> seedp = seed.roots;
  {
    const int    k = madness::FunctionDefaults<3>::get_k();
    const double t = madness::FunctionDefaults<3>::get_thresh();
    for (auto &r : seedp)
      for (auto *blk : r.blocks())
        for (auto &fn : *blk) fn = madness::project(fn, k, t);
  }

  // Overlap in the response metric (⟨X|X'⟩ for TDA, ⟨X|X'⟩−⟨Y|Y'⟩ for Full —
  // the same metric the solver orthonormalizes in), normalized per pair.
  auto O   = rs::metric(sf.roots, seedp);   // M x Ns
  auto Dcc = rs::metric(sf.roots, sf.roots);
  auto Dss = rs::metric(seedp, seedp);
  auto norm_ov = [&](int i, int j) {
    const double d = std::sqrt(std::abs(Dcc(i, i)) * std::abs(Dss(j, j)));
    return (d > 1.0e-14) ? std::abs(O(i, j)) / d : 0.0;
  };

  // Converged slot -> its seed slot, matched by stable identity (slots may
  // have been permuted by sort_state_by_omega; stable_index travelled along).
  std::vector<int> conv_stable = sf.stable_index;
  if (static_cast<int>(conv_stable.size()) != M) {
    conv_stable.resize(M);
    for (int i = 0; i < M; ++i) conv_stable[i] = i;
  }

  const long n_omega_seed = seed.omega.dim(0);
  for (int i = 0; i < M; ++i) {
    EsSeedGuardRoot r;
    r.slot         = i;
    r.stable_index = conv_stable[i];
    r.root_id      = make_root_id(r.stable_index);
    r.omega        = (i < sf.omega.dim(0)) ? sf.omega(i) : 0.0;

    for (int j = 0; j < Ns; ++j)
      if (seed.stable_index[j] == conv_stable[i]) { r.seed_slot = j; break; }
    if (r.seed_slot >= 0) {
      r.seed_overlap = norm_ov(i, r.seed_slot);
      if (r.seed_slot < n_omega_seed) {
        r.omega_seed  = seed.omega(r.seed_slot);
        r.omega_shift = r.omega - r.omega_seed;
      }
    }
    for (int j = 0; j < Ns; ++j) {
      const double ov = norm_ov(i, j);
      if (ov > r.best_overlap) { r.best_overlap = ov; r.best_seed_slot = j; }
    }
    if (r.best_overlap < kSeedGuardOverlapWarn)
      rep.basin_escape_slots.push_back(i);
    rep.roots.push_back(r);
  }

  // Pure tracking: every root kept high overlap with its own seed AND the
  // solve needed no more iterations than it has roots.
  bool all_high = true;
  for (const auto &r : rep.roots)
    if (r.seed_overlap <= kSeedGuardTrackingOverlap) { all_high = false; break; }
  rep.tracking_note = all_high && rep.iters <= rep.n_roots;
  return rep;
}

/// Rank-0 summary table + the loud warnings. Plain words, greppable tags.
inline void print_es_seed_guard(madness::World &world,
                                const EsSeedGuardReport &rep) {
  if (world.rank() != 0) return;
  {
  std::printf("\n[SEED-GUARD] seeded ES solve vs seed bundle %s "
              "(source key %s): iters=%d n_roots=%d converged=%d\n",
              rep.seed_bundle_dir.c_str(), rep.seed_protocol_key.c_str(),
              rep.iters, rep.n_roots, rep.converged ? 1 : 0);
  std::printf("[SEED-GUARD]  slot  root_id       omega_conv   omega_seed  "
              "omega_shift  seed_overlap  best(seed,ov)\n");
  for (const auto &r : rep.roots) {
    if (r.seed_slot >= 0)
      std::printf("[SEED-GUARD]  %4d  %-12s  %10.6f  %10.6f  %+10.6f  "
                  "%12.6f  (%d, %.6f)\n",
                  r.slot, r.root_id.c_str(), r.omega, r.omega_seed,
                  r.omega_shift, r.seed_overlap, r.best_seed_slot,
                  r.best_overlap);
    else
      std::printf("[SEED-GUARD]  %4d  %-12s  %10.6f  %10s  %10s  %12.6f  "
                  "(%d, %.6f)\n",
                  r.slot, r.root_id.c_str(), r.omega, "-", "-",
                  r.seed_overlap, r.best_seed_slot, r.best_overlap);
  }
  }  // seeded table

  for (int s : rep.basin_escape_slots) {
    const auto &r = rep.roots[static_cast<std::size_t>(s)];
    std::printf(
        "[SEED-GUARD] WARNING: root %d (%s, omega=%.6f) has overlap < %.2f "
        "with EVERY seed root (best %.4f vs seed %d) — the solve LEFT the "
        "seeded basin; this root is not the state the seed described.\n",
        r.slot, r.root_id.c_str(), r.omega, kSeedGuardOverlapWarn,
        r.best_overlap, r.best_seed_slot);
  }

  if (rep.tracking_note) {
    std::printf(
        "[SEED-GUARD] NOTE: every root kept overlap > %.2f with its seed and "
        "the solve converged in %d iteration(s) (<= n_roots=%d) — this was "
        "pure SEED TRACKING. A seeded ES solve refines the states it was "
        "handed; it does NOT search for the %d lowest states. If the seeding "
        "basis missed a low-lying state, that state is silently absent here. "
        "Cross-check against a better-basis ladder if one is available.\n",
        kSeedGuardTrackingOverlap, rep.iters, rep.n_roots, rep.n_roots);
  }

  std::fflush(stdout);
}

/// Record the guard block into the calc-level aggregate at
/// excited_states/<protocol_key>/seed_guard — through the ResponseMetadata
/// layer. Rank-0 write + fence (same discipline as save_es_roots).
inline void record_es_seed_guard(madness::World &world,
                                 const std::string &calc_dir,
                                 const std::string &protocol_key,
                                 const EsSeedGuardReport &rep) {
  if (world.rank() == 0) {
    auto meta =
        ResponseMetadata::load_or_create(calc_dir + "/response_metadata.json");
    meta.set_es_seed_guard(protocol_key, rep.to_json());
    meta.save();
  }
  world.gop.fence();
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_ES_SEED_GUARD_HPP
