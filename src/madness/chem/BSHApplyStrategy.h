/*
 * BSHApplyStrategy.h
 *
 * Backend selection for the compute_residual BSH (Helmholtz) apply.
 *
 * The decision layer of the tile/macrotask dispatch: a pure function of run-level
 * context (no World / FunctionDefaults), so it is unit-testable in isolation and the
 * choice is reproducible. See docs/bsh_apply_macrotask.md "IMPLEMENTATION SPEC".
 *
 * Model (settled by the investigation): macrotask's only tax is the cross-world input
 * gather (scales with protocol); tile pays a roughly-constant compute inefficiency plus
 * inter-node convolution comm (scales with node count) and holds a larger nonstandard
 * working set. Net: tile is the speed choice ONLY in the high-protocol + few-node +
 * fits-in-memory corner; macrotask wins or ties everywhere else, and is the only backend
 * that fits the large/high-precision/multinode target.
 *
 * Dispatch (the multinode choice is constant across the run -> decide-once cumulative win):
 *   - multinode            -> macrotask  (rank-local apply avoids inter-node conv comm,
 *                                          bounds memory; reaps the cumulative low/medium
 *                                          win for large systems)
 *   - single node, fits    -> tile       (faster: no gather, no inter-node comm)
 *   - would not fit         -> macrotask  (OOM-safe; a single-node safety net since tile's
 *                                          per-rank peak is bounded by the tile size)
 * Topology is the primary signal; the byte estimate is the OOM safety net (it cannot be
 * exact across molecule types, so it is deliberately not the load-bearing signal).
 */

#ifndef MADNESS_CHEM_BSHAPPLYSTRATEGY_H
#define MADNESS_CHEM_BSHAPPLYSTRATEGY_H

#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace madness {

/// Which compute_residual BSH-apply backend to run.
enum class BSHBackend { Tile, Macrotask };

/// Run-level inputs to the backend decision. Plain numbers (no World) so the chooser is
/// pure. The "decide once" that captures the multinode cumulative win is carried by
/// `n_nodes` (constant across a run); `k` is the CURRENT protocol's wavelet order, used
/// only for the single-node OOM estimate -- re-evaluating it per protocol is safe (and
/// slightly sharper: it picks tile while tile fits and switches exactly when it would not).
struct BSHApplyContext {
    long   nmo           = 0;    ///< number of orbitals to apply (psi.size())
    long   k             = 0;    ///< wavelet order at the current protocol
    long   nproc         = 1;    ///< total MPI ranks (world.size())
    long   n_nodes       = 1;    ///< distinct compute nodes spanned by the ranks
    double mem_budget_gb = 0.0;  ///< per-RANK memory budget (node RAM / ranks-per-node);
                                 ///< <= 0 disables the OOM gate
    int    ndim          = 3;
};

/// The chosen backend + knob settings, with a human-readable `reason` (logged at
/// print_level>=N so the dispatch is never a black box).
struct BSHApplyPlan {
    BSHBackend  backend  = BSHBackend::Tile;
    long        max_tile = 0;    ///< tile: cap orbitals/tile (0 = auto)
    long        nworld   = 0;    ///< macrotask: #subworlds (0 = auto = nproc)
    long        batch    = 0;    ///< macrotask: orbitals/task (0 = partitioner default);
                                 ///< set per-protocol via choose_bsh_batch()
    bool        redistribute = false;  ///< macrotask: pre-localize the operand to single-owner
                                 ///< batches up front (choose_bsh_redistribute); tight protocol only
    std::string reason;
};

/// Coarse, conservative estimate of the tile backend's peak resident set per rank (GB) at
/// the current protocol. Tiling bounds the in-flight set to `tile_size` orbitals'
/// nonstandard-form coefficients, so the per-rank peak is ~`(tile_size/nproc)` orbitals'
/// worth of `(2k)^ndim` coefficient blocks. The single constant `C_GB` is calibrated so a
/// C40-like point (nmo=161, k=10, nproc=8) lands near its measured ~34.5 GB/rank.
///
/// NOTE: this cannot be exact across molecule types (tree structure varies) and does NOT
/// model the multinode distributed-apply blow-up that actually OOMed w213 -- that case is
/// handled by topology (multinode -> macrotask), not this estimate. Used ONLY as the
/// single-node OOM safety net. Refined by Tier-2 calibration; see docs.
inline double estimate_tile_peak_per_rank_gb(const BSHApplyContext& c) {
    const long min_tile = 10;
    const long nproc = std::max<long>(1, c.nproc);
    const long tile_size = std::min(std::max<long>(0, c.nmo), min_tile * nproc);
    const double per_rank_orbitals = double(tile_size) / double(nproc);
    const double coeff_block = std::pow(double(2 * std::max<long>(1, c.k)), c.ndim);
    // Calibrated so C40 (nmo=161,k=10,nproc=8) -> ~34.4 GB/rank. PROVISIONAL: that
    // reference is from an ABBREVIATED (restart-ramp, 1 iter/protocol) run, so a full SCF
    // may peak higher -> this likely UNDER-estimates. The SAFETY margin below + topology
    // being the primary signal absorb most of it; re-calibrate against a full-run peakRSS
    // in Tier-3 (docs/bsh_apply_macrotask.md test plan).
    const double C_GB = 4.3e-4;
    return per_rank_orbitals * coeff_block * C_GB;
}

/// Decide the BSH-apply backend (decide-once, on peak demand). Topology primary; byte
/// estimate as the single-node OOM safety net. Fills `reason`; leaves the knob fields at
/// auto (0) -- batch is set per-protocol by choose_bsh_batch().
inline BSHApplyPlan choose_bsh_apply_plan(const BSHApplyContext& c) {
    BSHApplyPlan p;
    const double w = estimate_tile_peak_per_rank_gb(c);
    const double SAFETY = 0.8;   ///< headroom; a mid-apply bad_alloc is fatal
    char buf[256];
    if (c.mem_budget_gb > 0.0 && w > SAFETY * c.mem_budget_gb) {
        p.backend = BSHBackend::Macrotask;
        std::snprintf(buf, sizeof(buf),
            "macrotask: est. tile peak %.0f GB/rank exceeds %.0f%% of %.0f GB/rank budget (OOM-safe)",
            w, 100.0 * SAFETY, c.mem_budget_gb);
    } else if (c.n_nodes >= 2) {
        p.backend = BSHBackend::Macrotask;
        std::snprintf(buf, sizeof(buf),
            "macrotask: %ld nodes (multinode) -> rank-local apply avoids inter-node convolution comm + bounds memory",
            c.n_nodes);
    } else {
        p.backend = BSHBackend::Tile;
        std::snprintf(buf, sizeof(buf),
            "tile: single node, est. peak %.0f GB/rank fits -> faster (no gather, no inter-node comm)", w);
    }
    p.reason = buf;
    return p;
}

/// Topology+protocol-aware batch (orbitals per macrotask apply task). Two competing effects:
/// a larger batch amortizes the per-task framework/gather overhead (speed) but enlarges the
/// per-task nonstandard-form working set (memory) and, intra-node, becomes a compute-monolith
/// that blocks peers' gather-serving (the capacity-bound penalty).
///
///   - loose/medium protocol (thresh > 1e-7): the per-orbital working set is tiny, so batch
///     up freely (`b_lo`) -- pure overhead amortization, ~0 memory cost.
///   - tight protocol: default b=1 (protects memory; avoids the monolith penalty -- validated
///     intra-node, e.g. valinomycin 8/node b=4 was +20% SLOWER at 1e-8). Widen ONLY in the
///     inter-node, thinly-packed regime (>=2 nodes AND few ranks/node), where the apply is
///     framework/overhead-bound rather than monolith-bound: many orbitals spread over many
///     one-rank subworlds => many waves whose per-task overhead dominates, so a modest batch
///     amortizes it (w213 60rk/2-per: b~8 ~24% faster apply at 1e-8). The widen is gated on a
///     memory budget (per-task working set ~ batch*(2k)^ndim*C_GB, the same coeff model as the
///     tile estimator on one rank) kept to a small fraction of the per-rank budget, since the
///     apply shares the node with the (unmodeled) resident floor; with no budget we cannot
///     bound it, so stay at b=1.
///
/// Caller must still cap the result at floor(nmo / nsubworld). `thresh` is the CURRENT protocol.
/// PROVISIONAL constants (APPLY_FRACTION, B_HI, the ranks/node boundary) -- calibrated from a
/// single inter-node point (w213); see docs/bsh_apply_macrotask.md "Batch size in the model".
/// Whether to pre-redistribute the macrotask operand onto single-owner batches before the
/// apply (so the framework's auto_copy gather becomes a local, straggler-free fetch). Only
/// meaningful for the Macrotask backend -- the caller gates on that.
///
/// Fires ONLY at tight protocol (thresh <= 1e-7). There the naive one-sided auto_copy gather
/// degenerates into serve-starvation -- a single owner rank throttles all concurrent pulls
/// (measured 114 s at valinomycin 1e-8, k10) -- and the coalesced redistribute + local
/// fast-path collapse it (Apply BSH 265 -> 169 s multinode, -36%; -29% single-node C40 OOM).
/// At loose/medium protocol (thresh > 1e-7) the gather is small and NOT serve-starved, so the
/// redistribute's transport + materialization overhead is a NET LOSS (measured: k8 +65% cross-
/// run) and its ceiling is ~1% of wall behind transient system jitter -- so it is gated OFF.
/// Same 1e-7 boundary as choose_bsh_batch (the loose/tight protocol divide).
inline bool choose_bsh_redistribute(double thresh) {
    return thresh <= 1.0e-7;
}

inline long choose_bsh_batch(const BSHApplyContext& c, double thresh, long b_lo = 4) {
    if (thresh > 1.0e-7) return b_lo;          // loose/medium: working set tiny, amortize freely

    const long rpn = std::max<long>(1, c.nproc / std::max<long>(1, c.n_nodes));
    const bool inter_node_sparse = (c.n_nodes >= 2) && (rpn <= 2);
    if (!inter_node_sparse || c.mem_budget_gb <= 0.0) return 1;   // monolith regime / unbounded -> safe

    const double APPLY_FRACTION = 0.10;        ///< apply working set may use ~10% of per-rank budget
    const long   B_HI           = 8;           ///< cap (measured w213 sweet spot; avoids monolith)
    const double per_orb_ws = std::pow(2.0 * double(std::max<long>(1, c.k)), c.ndim) * 4.3e-4;  // GB
    const long b_mem = std::max<long>(1, long(APPLY_FRACTION * c.mem_budget_gb / std::max(1e-9, per_orb_ws)));
    return std::min<long>(B_HI, b_mem);
}

}  // namespace madness

#endif  // MADNESS_CHEM_BSHAPPLYSTRATEGY_H
