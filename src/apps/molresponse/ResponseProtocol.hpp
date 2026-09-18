#ifndef MOLRESPONSE_V3_RESPONSEPROTOCOL_HPP
#define MOLRESPONSE_V3_RESPONSEPROTOCOL_HPP

#include <madness/mra/funcdefaults.h>
#include <madness/mra/operator.h>
#include <madness/world/world.h>

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

/// Map a wavelet truncation threshold to its conventional polynomial order.
/// Mirrors the table in molresponse_v2/ResponseManager.cpp.
inline int default_k_for_thresh(double thresh) {
    if (thresh >= 0.9e-2) return 4;
    if (thresh >= 0.9e-4) return 6;
    if (thresh >= 0.9e-6) return 8;
    if (thresh >= 0.9e-8) return 10;
    return 12;
}

/// Inverse mapping: standard wavelet thresh for a given polynomial order k.
/// The archive's stored `k` is reliable; its `converged_for_thresh` is an
/// SCF tracking variable that is not necessarily the final protocol thresh.
/// So when reconstructing protocol intent from an archive, derive thresh
/// from k via this table rather than trusting the header's thresh field.
inline double default_thresh_for_k(int k) {
    switch (k) {
        case 4:  return 1e-2;
        case 6:  return 1e-4;
        case 8:  return 1e-6;
        case 10: return 1e-8;
        case 12: return 1e-10;
        default: return 1e-6;   // unknown k — safe middle of the table
    }
}

/// Canonical, filename-safe protocol identity shared by FD and ES
/// persistence (doc 13). The join key for property matching is the physical
/// accuracy (thresh, k) — NOT a positional ramp index, which collides across
/// runs with different ramps and can't assert "same accuracy". `%.0e` yields
/// a stable two-digit exponent on glibc, e.g. protocol_key(1e-6, 8) ==
/// "1e-06_k8". Used for both archive-name suffixes and JSON keys.
inline std::string protocol_key(double thresh, int k) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.0e_k%d", thresh, k);
    return {buf};
}

/// protocol_key built from the active FunctionDefaults<3> — the common funnel
/// both solvers configure via set_response_protocol() before solving/saving.
inline std::string protocol_key() {
    return protocol_key(FunctionDefaults<3>::get_thresh(),
                        FunctionDefaults<3>::get_k());
}

/// Configure FunctionDefaults<3> for one response-protocol step.
///
/// This is the v3 analog of `ResponseManager::setProtocol` in molresponse_v2.
/// Both the test harness and the production `molresponse_v3` app need to call
/// this *before* loading orbitals so that re-projection lands at the intended
/// (k, thresh) — otherwise a tightly-converged ground state silently
/// downsamples to MADNESS's defaults (k=6, thresh=1e-4).
///
/// @param world      MPI world
/// @param L          half-edge of the cubic simulation cell (FunctionDefaults
///                   will be set to [-L, L]^3)
/// @param thresh     wavelet truncation threshold for this protocol step
/// @param override_k optional polynomial order; if -1 the standard mapping
///                   from `default_k_for_thresh(thresh)` is used
inline void set_response_protocol(World& world,
                                  double L,
                                  double thresh,
                                  int override_k = -1,
                                  bool verbose = true) {  // F2d-ii: false in subworlds
    int k = (override_k > 0) ? override_k : default_k_for_thresh(thresh);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_project_randomize(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    // Convolution caches are k/thresh-specific; clear so cached operators
    // built at the previous protocol can't be reused incorrectly.
    GaussianConvolution1DCache<double>::map.clear();

    if (verbose && world.rank() == 0) {
        print("PROTOCOL_SET  L=", L, "  thresh=", thresh, "  k=", k,
              (override_k > 0 ? "  (k overridden)" : ""));
    }
    world.gop.fence();
}

/// Build the standard coarse-to-fine protocol ramp ending at `target_thresh`.
///
/// Walks down by factor 100 from `coarse` to `target`, matching MADNESS's
/// usual protocol table:
///   target=1e-4  →  [1e-4]
///   target=1e-6  →  [1e-4, 1e-6]
///   target=1e-8  →  [1e-4, 1e-6, 1e-8]
///   target=1e-10 →  [1e-4, 1e-6, 1e-8, 1e-10]
///
/// If `target` is coarser than `coarse`, returns just `[target]` (no ramp).
inline std::vector<double> build_protocol_ramp(double target_thresh,
                                                double coarse = 1e-4) {
    std::vector<double> ramp;
    if (target_thresh >= coarse) {
        ramp.push_back(target_thresh);
        return ramp;
    }
    // Step from coarse down to target by ×100 increments. Use a small
    // multiplicative tolerance to compare doubles cleanly.
    constexpr double step_factor = 0.01;       // each step tightens by 100x
    constexpr double cmp_eps = 1.001;
    double t = coarse;
    while (t >= target_thresh * cmp_eps) {
        ramp.push_back(t);
        t *= step_factor;
    }
    // Make sure the final entry is exactly the target (handles non-decade
    // targets like 5e-7 cleanly).
    if (ramp.empty() || ramp.back() / cmp_eps > target_thresh) {
        ramp.push_back(target_thresh);
    } else {
        ramp.back() = target_thresh;
    }
    return ramp;
}

/// Protocol-start control for SEEDED response runs (`seed.start_rung` /
/// --seed-start-rung). Shared by the standalone driver and the madqc deck
/// adapter so the two surfaces can never disagree.
///
///   mode == "coarse" (default) — no-op: the full coarse->fine ladder runs,
///     today's behavior unchanged.
///   mode == "fine" AND have_seed — collapse the ladder to its FINEST rung.
///     The caller then prepares the ground state and projects the seed at
///     that rung's (k, thresh), and the manager schedules only that rung.
///
/// Rationale (W6 between-pole finding, note 2026-08-17): between poles the
/// coarse rung (1e-4/k6) hits maxiter unconverged in every arm — cold or
/// seeded — and launders away the seed's head start; by the time the fine
/// rung starts, both arms inherit the same half-converged coarse state. A
/// DALTON seed is already at the physics; refining it at the coarse rung
/// first is a detour. This matches the validated ES warm-start semantics
/// ("a near-converged seed lets the solve start directly at the production
/// resolution and skip the coarse warmup rung entirely" —
/// docs/guides/excited_states.md), applied to the FD ladder.
///
/// mode == "fine" WITHOUT a seed is a no-op (returns false); the caller
/// warns. A COLD fine start needs no flag — it is just a single-rung
/// protocol request. Any other mode string throws (deck typo = loud error).
///
/// Returns true iff the ladder was collapsed.
inline bool apply_seed_start_rung(std::vector<double>& protocol,
                                  const std::string& mode,
                                  bool have_seed) {
    if (mode.empty() || mode == "coarse") return false;
    if (mode != "fine")
        throw std::invalid_argument(
            "seed start rung must be 'coarse' or 'fine', got '" + mode + "'");
    if (!have_seed || protocol.size() <= 1) return false;
    const double fine = *std::min_element(protocol.begin(), protocol.end());
    protocol.assign(1, fine);
    return true;
}

/// Parse a comma-separated protocol list, e.g. "1e-4,1e-6,1e-8".
/// Whitespace around values is allowed. Returns empty vector if input is
/// blank.
inline std::vector<double> parse_protocol_csv(const std::string& csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        // trim whitespace
        size_t a = tok.find_first_not_of(" \t");
        size_t b = tok.find_last_not_of(" \t");
        if (a == std::string::npos) continue;
        out.push_back(std::stod(tok.substr(a, b - a + 1)));
    }
    return out;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSEPROTOCOL_HPP
