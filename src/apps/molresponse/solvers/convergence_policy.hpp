#ifndef MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
#define MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP

// =========================================================================
// ConvergencePolicy — shared by ESSolver<T,S> and (future) FDSolver<T,S>.
//
// Policy carries USER intent (target dconv, behaviour switches). At
// each protocol level the policy resolves to EFFECTIVE TARGETS using the
// SCF convention (SCF.cc::solve, lines 2143 + 2382):
//
//   dconv          = max(thresh, dconv_user)          // can't beat the protocol
//   density_target = density_residual_factor * dconv  // SCF uses max(5,natom)
//   bsh_target     = bsh_residual_factor     * dconv  // SCF uses 5.0
//
// `dconv = max(thresh, dconv_user)` is the key idea: at a given protocol
// you cannot resolve a residual below ~thresh (past that you're just
// resolving noise until the next, deeper protocol refines the wavelet
// basis), so dconv floors at thresh; and you never ask for tighter than
// the user wants, so it also floors at dconv_user. BOTH gates then ride at
// ~5× dconv (SCF.cc:2382 uses 5*dconv for bsh and dconv*max(5,natom) for
// density). The 5× headroom is what lets a single coarse protocol reach
// "converged" rather than stalling right at the thresh noise floor — at
// thresh == dconv_user the targets are 5*thresh, comfortably above the
// per-protocol floor.
//
// TODO(study): the factor 5.0 and dconv_user defaults are inherited from
// SCF; response may want different values. Increasing thresh (deeper
// protocol) is what actually buys resolution — within a protocol the
// residual floors at noise after a few iters. A sweep over
// bsh_residual_factor / protocol ladder is planned to pick response-
// appropriate values.
// =========================================================================

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace molresponse_v3 {

struct ConvergencePolicy {
  // User-facing density-convergence target. The effective target at each
  // protocol is a multiple of max(thresh, dconv_user) — see
  // effective_for_thresh.
  double dconv_user = 1.0e-4;

  // Both gates ride at a multiple of the (protocol-floored) dconv,
  // matching SCF::solve (SCF.cc:2382):
  //   bsh     < 5.0 * dconv
  //   density < dconv * max(5, natom)   -> 5.0 for small molecules
  // Using 5.0 for both gives each gate ~5x headroom over the per-protocol
  // noise floor, so a single coarse protocol can actually reach
  // "converged" instead of stalling right at thresh. (SCF's density gate
  // loosens further with atom count; we keep a flat factor here since the
  // response per-state density residual isn't naturally natom-scaled —
  // revisit in the planned convergence study.) Exposed so the study can
  // sweep them independently.
  double bsh_residual_factor     = 5.0;
  double density_residual_factor = 5.0;
  // ES convergence: the eigenvalue (omega) is the physical deliverable. For a
  // near-degenerate ES root the raw BSH amplitude residual ||Δx|| jitters in the
  // loosely-constrained degenerate direction long after omega + density are
  // settled, so ES convergence is gated on density + |Δω| (NOT ||Δx||; see
  // ESSolver::es_root_converged). omega target = omega_residual_factor * dconv.
  // Rides at the same ~5× headroom as bsh/density (see the note above): at a
  // rung where thresh == dconv the |Δω| floor is the truncation noise ~thresh,
  // so a factor of 1.0 makes the gate EQUAL the noise floor — it can never
  // pass at the coarse rung (the 1e-4/k6 0-of-4 stall).
  double omega_residual_factor   = 5.0;

  // Cluster-unmix threshold factor in rs::diagonalize. The legacy
  // path uses 100·thresh for TDA (loose) and 10·thresh for Full/RPA
  // (tighter — clusters with smaller separation get polar-decomp
  // unmixed). Tighter unmix is more aggressive at freezing eigenvector
  // mixing within a cluster, which stabilises slot identity across
  // iters when omega's are closely spaced (the Li ω=0.04 regime).
  // Default 100 matches the TDA convention; set to 10 for Full or any
  // run with closely-clustered omegas.
  double cluster_unmix_factor = 100.0;

  // ESSolver memory mode. When true, step() streams Lambda assembly
  // and DROPS V0x/E0x/gamma after the subspace rotation, then
  // RECOMPUTES the rotated pieces to assemble Theta. Costs ~2× the
  // per-root kernel work (V0x/E0x/gamma evaluated twice per iter)
  // but cuts peak memory from ~8M Storage to ~3-4M. Set true when
  // n_occ × n_roots × k³ × leaves is the memory bottleneck.
  // Default false: current "keep pieces, rotate them, assemble Theta
  // from rotated pieces" layout.
  bool stream_theta = false;

  // Gate the FD theta assembly onto the tensor-exchange layer
  // (kernels/exchange_ctx.hpp): false (default) = the per-op REFERENCE path
  // (compute_V0x − compute_E0x + compute_gamma, untouched); true = build_ctx +
  // assemble_theta (ClosedShell Static/Full only; A/B-to-thresh vs reference).
  // CLI: --fd-tensor (doc 28 §4 Inc 1).
  bool exchange_tensor = false;

  // Tile size for the tensor-exchange Tx/Ty/g0 builds over the φ-row index
  // (kernels/exchange_ctx.hpp build_pair_tensors): 0 (default) = no tiling (one
  // fused Poisson wave, peak ~n²); >0 = block the φ rows to bound the peak to
  // ~tile·n at the cost of more waves (memory ↔ fences). Bit-identical to tile=0.
  // Used when exchange_tensor is on (FD) and by ESSolver's --es-tensor g0 build.
  // CLI: --fd-tensor-tile=N (doc 33 Inc-3b).
  int exchange_tile = 0;

  // Diverging-residual bail-out. Triggers only on a runaway residual
  // (BSH residual > guard). The legacy ES guard was 2.0 in normalised
  // units, but FD's iter-1 residual from a "x = perturbation" guess is
  // O(perturbation norm) which can easily exceed 2.0 for closed-shell
  // first-row molecules — so 2.0 false-fires on iter 1 in FD. Bumped
  // to 1e3 (effectively off) until we wire a relative-growth check
  // (compare iter k vs iter k-1, bail if it grew by >2×).
  double explosion_guard = 1.0e3;

  // Minimum iters before we even check convergence (lets KAIN warm up).
  int min_iters_before_conv = 0;

  // ---- Plateau (stall) detector, 2026-09-11 ----
  // Within one protocol the residual bottoms out at a truncation-noise floor
  // (measured: BSH amplitude residual 2-6e-5 at thresh 1e-6/k8, 4-8e-6 at
  // 1e-8/k10, while the density residual keeps falling). A target below that
  // floor is never met and the solve burns every remaining iteration doing
  // nothing (the 1e-8 closeout attempt with dconv = thresh: 60 iterations per
  // leg, flat). The detector watches the normalised gate distance
  //   g = max_c max(bsh_c / bsh_target, drho_c / density_target)   (FD)
  //   g = max_s max(drho_s / density_target, |dw_s| / omega_target) (ES, active roots)
  // and declares a stall when g has not improved by at least `stall_ratio`
  // (relative) over the last `stall_window` iterations while still > 1. A
  // stalled solve exits the loop with converged = false and State::stalled =
  // true; the executor records `stalled` in the metadata and applies the same
  // best-effort acceptance it applies at maxiter (--accept-at-maxiter). So the
  // verdict is unchanged - only the wasted iterations are gone. The history is
  // reset whenever the protocol thresh changes (targets change with it).
  // stall_window <= 0 disables. Deck: response { stall.window 6  stall.ratio 0.1 }.
  int    stall_window = 6;
  double stall_ratio  = 0.10;

  /// Plateau test. `tracks` holds one history per gated quantity (FD: the BSH
  /// residual and the density change; ES: the density change and |dw|), each
  /// normalised by its own target, one entry per iteration, newest last.
  ///
  /// A solve is stalled only when EVERY quantity has stopped improving. Testing
  /// the max over quantities instead was wrong and cost a production run: the
  /// nuclear leg of a Raman solve has a density change that barely moves while
  /// its BSH residual falls by a factor of two every few iterations, and the
  /// max-gate reads that as a plateau and stops a converging solve (h2o oxygen
  /// displacements, 2026-09-16: stopped at iteration 7 with residual 1.2 where
  /// the equivalent hydrogen legs reached 9e-5 by iteration 11).
  ///
  /// A track is "not improving" when its newest entry is above 1 (its target is
  /// not met) and is no smaller than (1 - stall_ratio) x its entry
  /// `stall_window` iterations earlier. Tracks are compared only over entries
  /// that measure the same thing, so a history must be at least window+1 long;
  /// a non-finite entry (a quantity with no measurement yet, e.g. the first
  /// iteration's density change or an ES eigenvalue step) makes its track
  /// inconclusive, and one inconclusive track is enough to withhold a stall.
  bool plateau(const std::vector<std::vector<double>> &tracks) const {
    if (stall_window <= 0 || tracks.empty()) return false;
    const auto w = static_cast<std::size_t>(stall_window);
    for (const auto &h : tracks) {
      if (h.size() <= w) return false;                    // not enough history yet
      const double now = h[h.size() - 1], then = h[h.size() - 1 - w];
      if (!std::isfinite(now) || !std::isfinite(then)) return false;  // inconclusive
      if (!(now > 1.0)) return false;                     // this one has met its target
      if (now <= (1.0 - stall_ratio) * then) return false;  // this one is still improving
    }
    return true;
  }

  /// Per-root plateau test (ES). `tracks` are ONE root's histories (density
  /// change, |dw|), normalised as for plateau(). Unlike plateau(), a track that
  /// has met its target does not block the verdict: the root is stuck when it
  /// has at least one unmet track and every unmet track is flat. plateau() takes
  /// the max over roots, so a root whose |dw| has converged (second order, and
  /// usually first to settle) while its density change sits flat above target
  /// makes the |dw| track read "met", and the solve-level test never fires.
  /// A met track needs only its newest entry; unmet tracks need the same
  /// finite window+1 history as plateau().
  bool root_plateau(const std::vector<std::vector<double>> &tracks) const {
    if (stall_window <= 0 || tracks.empty()) return false;
    const auto w = static_cast<std::size_t>(stall_window);
    bool unmet = false;
    for (const auto &h : tracks) {
      if (h.size() <= w) return false;                    // not enough history yet
      const double now = h[h.size() - 1], then = h[h.size() - 1 - w];
      if (!std::isfinite(now)) return false;              // inconclusive
      if (now <= 1.0) continue;                           // met: need not improve
      if (!std::isfinite(then)) return false;             // inconclusive
      if (now <= (1.0 - stall_ratio) * then) return false;  // still improving
      unmet = true;
    }
    return unmet;
  }

  // Lock debounce (ESSolver full-deflation locking): a root must satisfy the
  // convergence criterion for this many CONSECUTIVE iters before it is locked.
  // Prevents premature locking of an unsettled root (which poisoned the
  // deflation -> spurious roots in the 4-root h2o solid_lock sweep). >=2.
  int lock_min_pass = 2;

  // ---- KAIN acceleration + step restriction ----
  // KAIN is enabled by default. The solver allocates an
  // XNonlinearSolver over the flat state vector; size of the KAIN
  // subspace history is kain_maxsub. Disable by setting kain=false.
  bool kain = true;
  int  kain_maxsub = 5;
  // Blowup safeguard threshold for KAIN coefficient max-abs. When the
  // KAIN solve returns coefficients with |c|max exceeding this value
  // (after rcond escalation has been exhausted), the per-state history
  // is cleared and the iteration falls back to raw BSH for that step.
  //
  // SCF::update_subspace uses 3.0 — appropriate for orbital optimization
  // where large |c| genuinely indicates instability (Brillouin theorem
  // etc.). RESPONSE iterations are linear; large |c| just means strong
  // acceleration toward a polarizable mode (e.g. Li at α≈170). The
  // SCF-style 3.0 trips bailout every iter for such systems and kills
  // KAIN entirely. Bumped to a permissive 100 by default so KAIN can
  // actually accelerate strongly-polarizable response iterations;
  // back-set to 3.0 if you want strict SCF semantics.
  double kain_cmax_cap = 100.0;
  // KAIN hold-off for rough starts (2026-09-10). KAIN is applied in an
  // iteration only when that iteration's raw BSH residual (max over roots)
  // is already below this value; above it the step is the plain BSH update
  // and no history is recorded. Why: a DALTON-seeded ES solve enters with a
  // 20-25 % residual; at its second iteration KAIN extrapolated from a
  // two-vector history with coefficients of 5-8, the step-restriction cap
  // then scaled that garbage direction into the state, and the next Ritz
  // matrix had a negative eigenvalue (the "eps_core ghost", closeout
  // attempts 3-10, lih/h2o/c2h4). The cold path never exposes KAIN to such a
  // start because its TDA warm-up is KAIN-free. 0 disables the hold-off.
  double kain_min_residual = 0.1;
  // Warmup oversampling factor — used by the run_oversampled_tda_warmup
  // helper. The warm-up phase runs with ceil(warmup_oversample_factor *
  // n_roots) trial states; after warmup completes the lowest n_roots
  // are kept for the main iteration. Mirrors legacy iterate_trial's
  // "trial bundle is 2× requested states, then select_functions picks
  // the N lowest" pattern (ExcitedResponse.cpp:83-105 and select_
  // functions). Default 1.0 (no oversampling — back-compat). Set to
  // 2.0 for cold CIS-style guesses; the extra trial vectors absorb
  // higher-root contamination so the kept N converge cleanly.
  // Only takes effect if tda_warmup_iters > 0.
  double warmup_oversample_factor = 1.0;

  // TDA-warmup iters: KAIN is disabled for the first `tda_warmup_iters`
  // step() calls, then turns on. Mirrors legacy iterate_trial's
  // "no-KAIN BSH-only filter" pre-pass (ExcitedResponse.cpp:456-666):
  // a few power-iteration steps clean up a cold guess before KAIN's
  // per-slot history starts recording. Without this, KAIN memorizes
  // the unstable iter-0/iter-1 slot ordering and slot identity
  // flickers across iters when the rotation re-permutes. Default 0
  // (KAIN on from iter 1, current behaviour). Recommended 5-10 for
  // CIS-style guesses, 0 for restart from converged lower-protocol
  // state. (Affects ESSolver only; ignored by FDSolver.)
  int  tda_warmup_iters = 0;
  // Step-restriction: if ||x_new − x_old|| > maxrotn after KAIN,
  // damp the step toward x_old. Set < 0 to disable.
  double maxrotn = 0.5;

  // Step-restriction granularity:
  //   PerOrbital — clamp each response function independently (the legacy
  //                SCF::do_step_restriction behaviour). A single runaway
  //                orbital is damped without touching the well-behaved ones.
  //   PerState   — (DEFAULT, set in ca2c21dfc) measure ||v_new - v_old|| over
  //                the WHOLE flattened state (all functions of a root / channel)
  //                and scale the ENTIRE vector by one factor maxrotn/||diff||.
  //                Damps a coherent rotation of the whole state (useful when a
  //                near-degenerate root rotates as a block rather than one
  //                orbital running away).
  enum class StepRestrictMode { PerOrbital, PerState };
  StepRestrictMode step_restrict_mode = StepRestrictMode::PerState;

  // Full-deflation locking of converged ES roots (workstream B). When true, a
  // root that meets the convergence test is locked: removed from the subspace +
  // rotation, used only to deflate (orthogonalize) the active roots, and skipped
  // by KAIN/step-restriction. Stops the per-iter rotation from re-mixing already-
  // converged roots (the near-degenerate failure mode). Default off. ESSolver only.
  bool lock_converged = false;

  struct Targets {
    double bsh_residual;     // ‖x_old − x_new‖ cap (FD gate; ES sanity only)
    double density_residual; // ‖ρ_new − ρ_old‖ cap
    double omega_residual;   // |ω_new − ω_old| cap (ES eigenvalue gate)
  };

  Targets effective_for_thresh(double thresh) const {
    // SCF convention (SCF.cc:2143): can't resolve below the protocol's
    // thresh, won't over-tighten past the user's request.
    const double dconv = std::max(thresh, dconv_user);
    Targets t;
    t.density_residual = density_residual_factor * dconv;  // SCF.cc:2382 (da < dconv*max(5,natom))
    t.bsh_residual     = bsh_residual_factor * dconv;      // SCF.cc:2382 (bsh < 5*dconv)
    t.omega_residual   = omega_residual_factor * dconv;    // ES eigenvalue gate
    return t;
  }
};

/// ES solve-level verdict from per-root flags, one entry per slot: the solve is
/// stalled when at least one active (unlocked) root has plateaued
/// (ConvergencePolicy::root_plateau) and every other active root is either
/// converged or plateaued too, so no active root is still making progress.
inline bool es_roots_stalled(const std::vector<char> &active,
                             const std::vector<char> &converged,
                             const std::vector<char> &plateaued) {
  bool stuck = false;
  for (std::size_t s = 0; s < active.size(); ++s) {
    if (!active[s]) continue;
    const bool conv = s < converged.size() && converged[s];
    const bool plat = s < plateaued.size() && plateaued[s];
    if (conv) continue;
    if (!plat) return false;  // this root is still iterating
    stuck = true;
  }
  return stuck;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_CONVERGENCE_POLICY_HPP
