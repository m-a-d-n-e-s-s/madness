#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP

// =========================================================================
// ResponseSubspaceKain<Storage> — hand-rolled per-state KAIN +
// per-orbital step restriction, modelled on SCF::update_subspace
// (src/madness/chem/SCF.cc:1916).
//
// Direct port of the SCF algorithm rather than XNonlinearSolver, so we
// can inspect the coefficient vector `c` and recover when KAIN's
// subspace becomes ill-conditioned. Each state (root for ES, channel
// for FD) gets its own independent history and Q matrix — different
// states converge at different rates and a shared history blends
// inconsistent iterates.
//
// Residual sign convention (matches SCF::update_subspace):
//
//     r = v − F(v)         (positive when v hasn't converged)
//     Q[i,j] = <v_i | r_j>
//     new_v = Σ_m c_m · (v_m − r_m) = Σ_m c_m · F(v_m)     ← KAIN(Q, rcond)
//
// where F is the BSH fixed-point map. The XNonlinearSolver-based
// KainAccelerator used the OPPOSITE residual sign (r = F(v) − v); we
// flip here so the SCF formula `c·(v − r)` works verbatim.
//
// Blowup safeguard (simplified vs SCF's `goto restart`):
//
//   while: c = KAIN(Q, rcond)
//          if |c|_max < 3:           accept
//          elif rcond < 0.01:        rcond *= 100, retry
//          else:                     clear history, return raw BSH
//                                    (skip KAIN this iter — simpler than
//                                    SCF's recursive restart)
//
// TDA warmup (policy.tda_warmup_iters):
//   For the first `tda_warmup_iters` apply() calls, KAIN extrapolation
//   is skipped AND no history is recorded. After warmup, history
//   recording begins fresh — KAIN starts learning from the warmed-up
//   bundle, not from the volatile iter-0/iter-1 state. Mirrors legacy
//   iterate_trial which runs pure BSH+GS for `guess_max_iter` iters
//   before handing off to the KAIN-accelerated main iter
//   (ExcitedResponse.cpp:456-666). Step restriction still applies
//   through warmup, so individual orbitals can't run away.
//
// Step restriction (per-orbital, mirrors SCF::do_step_restriction at
// SCF.cc:2058) is applied AFTER KAIN, also AFTER raw BSH when KAIN is
// disabled or skipped. For each function p in the flat state:
//
//     if ||v[p] − new_v[p]|| > maxrotn:
//         new_v[p] := s · new_v[p] + (1−s) · v[p],  s = maxrotn / norm
//
// =========================================================================

#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "response_state.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>   // madness::KAIN

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace molresponse_v3 {

template <typename Storage>
class ResponseSubspaceKain {
public:
  using vecfuncT = std::vector<madness::real_function_3d>;

  ResponseSubspaceKain(madness::World &world, const ConvergencePolicy &policy)
      : world_(world), policy_(policy) {}

  /// Apply per-state KAIN (if enabled) then per-orbital step restriction.
  ///   in[s]  = x_old for state s
  ///   out[s] = raw BSH output for state s; overwritten in place
  ///
  /// Warmup behaviour: skips KAIN extrapolation for the first
  /// `policy.tda_warmup_iters` calls (in addition to the always-skipped
  /// iter-0). During warmup the history is NOT recorded — when KAIN
  /// turns on, it starts fresh on the warmed-up bundle so its
  /// extrapolation isn't poisoned by the rapidly-changing early-iter
  /// state. Step restriction still applies through the warmup.
  /// `diag_level`: 0 = silent; ≥1 = one `[KAIN]` line per state per call
  /// with history length, rcond, max|c|, bailout flag, plus
  /// `[STEP-REST]` when at least one orbital is clamped. Caller derives
  /// this from print_level (uniform across ranks). All collective work
  /// (norm2s, sub, gaxpy, inner) is performed unconditionally; only
  /// the printf is rank-gated.
  void apply(const std::vector<Storage> &in,
             std::vector<Storage>       &out,
             int                         diag_level = 0,
             const std::vector<char>    *locked    = nullptr) {
    if (states_.size() != in.size()) {
      states_.clear();
      states_.resize(in.size());
    }

    const bool in_warmup = (calls_consumed_ < policy_.tda_warmup_iters);
    auto is_locked = [&](std::size_t s) {
      return locked && s < locked->size() && (*locked)[s];
    };

    if (policy_.kain && !in_warmup) {
      if (first_kain_call_) {
        first_kain_call_ = false;
        // First post-warmup call: skip extrapolation; history empty.
        // Step restriction still applies below.
        if (diag_level >= 1 && world_.rank() == 0) {
          printf("[KAIN] iter=%d SKIP (first post-warmup call, history empty)\n",
                 calls_consumed_);
          fflush(stdout);
        }
      } else {
        for (std::size_t s = 0; s < in.size(); ++s) {
          if (is_locked(s)) continue;   // locked root: freeze, keep history
          apply_one_state(s, in[s], out[s], diag_level);
        }
      }
    }

    apply_step_restriction(in, out, diag_level, locked);
    ++calls_consumed_;
  }

  /// Drop all per-state KAIN histories. Call between protocol steps.
  /// Resets the warmup counter too — protocol transitions get a fresh
  /// warmup window (the projected state at the new k/thresh is a cold
  /// guess relative to the post-BSH dynamics at the new resolution).
  void reset() {
    states_.clear();
    first_kain_call_ = true;
    calls_consumed_  = 0;
  }

  void set_policy(const ConvergencePolicy &p) { policy_ = p; }

private:
  struct PerStateHistory {
    std::vector<std::pair<vecfuncT, vecfuncT>> entries; // (v_m, r_m)
    madness::Tensor<double>                    Q;        // m × m
  };

  /// One state's KAIN update. Mirrors SCF::update_subspace lines
  /// 1957-2028 step-for-step. `diag_level` controls per-call printing.
  void apply_one_state(std::size_t s_idx, const Storage &s_in,
                       Storage &s_out, int diag_level) {
    auto v = s_in.flatten();
    auto x_bsh = s_out.flatten();
    // SCF convention: residual r = v − F(v) (NOT F(v) − v)
    auto r = madness::sub(world_, v, x_bsh);

    auto &state = states_[s_idx];
    state.entries.push_back({v, r});
    const long m = static_cast<long>(state.entries.size());

    // ---- build new row + column of Q (SCF:1962-1979) -------------------
    madness::Tensor<double> ms(m), sm(m);
    for (long k = 0; k < m; ++k) {
      ms[k] = madness::inner(v, state.entries[k].second);
      sm[k] = madness::inner(state.entries[k].first, r);
    }
    madness::Tensor<double> newQ(m, m);
    if (m > 1) {
      newQ(madness::Slice(0, -2), madness::Slice(0, -2)) = state.Q;
    }
    newQ(m - 1, madness::_) = ms;
    newQ(madness::_, m - 1) = sm;
    state.Q = newQ;

    // ---- regularized solve with rcond escalation (SCF:1983-2004) -------
    // Coefficient cap (policy.kain_cmax_cap) controls the bailout
    // threshold — SCF uses 3.0 (linear-orbital instability bound),
    // response iterations typically need 50-100 to accelerate strongly
    // polarizable systems where the leading eigenvector dominates the
    // residual and the natural KAIN c are intrinsically large.
    double rcond = 1.0e-12;
    int    n_rcond_escalations = 0;
    const double cmax_cap = policy_.kain_cmax_cap;
    madness::Tensor<double> c;
    while (true) {
      c = madness::KAIN(state.Q, rcond);
      if (c.absmax() < cmax_cap) break;
      if (rcond < 0.01) {
        rcond *= 100.0;
        ++n_rcond_escalations;
        continue;
      }
      // Bail-out: clear history, leave s_out as raw BSH. KAIN will
      // start re-learning from the next iter. (Simpler than SCF's
      // recursive restart; the user opted for "clear and skip".)
      const double cmax_bail = c.absmax();
      state.entries.clear();
      state.Q = madness::Tensor<double>();
      if (diag_level >= 1 && world_.rank() == 0) {
        printf("[KAIN] iter=%d state=%zu m=%ld rcond=%.1e |c|max=%.3f "
               "n_escalations=%d  ★ BAILED (cleared history, fell back "
               "to raw BSH)\n",
               calls_consumed_, s_idx, m, rcond, cmax_bail,
               n_rcond_escalations);
        fflush(stdout);
      }
      return;
    }

    // ---- linear combination (SCF:2013-2028) ----------------------------
    // new_v = Σ_k c_k · (v_k − r_k) = Σ_k c_k · F(v_k)
    auto new_v = madness::zero_functions<double, 3>(world_, v.size());
    for (long k = 0; k < m; ++k) {
      const auto &vm = state.entries[k].first;
      const auto &rm = state.entries[k].second;
      madness::gaxpy(world_, 1.0, new_v,  c(k), vm);
      madness::gaxpy(world_, 1.0, new_v, -c(k), rm);
    }

    // ---- sliding window (SCF:2030-2035) --------------------------------
    if (m >= policy_.kain_maxsub) {
      state.entries.erase(state.entries.begin());
      state.Q = state.Q(madness::Slice(1, -1), madness::Slice(1, -1));
    }

    s_out.from_flat(new_v);

    // ---- diagnostic --------------------------------------------------
    if (diag_level >= 1 && world_.rank() == 0) {
      // Compact summary: identify whether KAIN is dominated by the
      // latest entry (|c[m-1]| ≈ 1, rest ≈ 0 → "use raw BSH only",
      // no acceleration) or actually mixing history (|c| spread).
      double c_latest = std::fabs(c(m - 1));
      double c_max    = c.absmax();
      printf("[KAIN] iter=%d state=%zu m=%ld rcond=%.1e "
             "|c|max=%.3f c[m-1]=%.3f  n_escalations=%d  c=[",
             calls_consumed_, s_idx, m, rcond, c_max, c_latest,
             n_rcond_escalations);
      for (long k = 0; k < m; ++k) printf(" %+.3e", c(k));
      printf(" ]\n");
      fflush(stdout);
    }
  }

  /// Per-orbital step restriction (SCF::do_step_restriction, lines
  /// 2058-2084). Iterates per function in the flat state; each gets
  /// its own scale factor. A single bad orbital doesn't dominate the
  /// per-state norm and over-restrict the well-behaved ones.
  void apply_step_restriction(const std::vector<Storage> &in,
                              std::vector<Storage>       &out,
                              int                         diag_level,
                              const std::vector<char>    *locked = nullptr) {
    const double maxrotn = policy_.maxrotn;
    if (maxrotn <= 0.0) return;

    const bool per_state =
        (policy_.step_restrict_mode ==
         ConvergencePolicy::StepRestrictMode::PerState);

    for (std::size_t s = 0; s < in.size(); ++s) {
      if (locked && s < locked->size() && (*locked)[s]) continue;
      auto v_old = in[s].flatten();
      auto v_new = out[s].flatten();
      // Per-orbital diff norms via batched norm2s (collective-safe).
      auto diff  = madness::sub(world_, v_new, v_old);
      auto norms = madness::norm2s(world_, diff);
      bool changed = false;

      if (per_state) {
        // ---- state-wise: one norm over the whole flattened state, one scale.
        double total2 = 0.0;
        for (double n : norms) total2 += n * n;
        const double total = std::sqrt(total2);
        if (total > maxrotn) {
          const double scale = maxrotn / total;
          for (std::size_t p = 0; p < v_new.size(); ++p)
            v_new[p].gaxpy(scale, v_old[p], 1.0 - scale, false);
          changed = true;
          world_.gop.fence();
          out[s].from_flat(v_new);
          if (diag_level >= 1 && world_.rank() == 0) {
            printf("[STEP-REST] iter=%d state=%zu STATE-WISE  ||Δstate||=%.3e "
                   "scale=%.3f  cap=%.3e\n",
                   calls_consumed_, s, total, scale, maxrotn);
            fflush(stdout);
          }
        }
        continue;
      }

      // ---- per-orbital (default): clamp each function independently.
      int  n_clamped = 0;
      double max_norm = 0.0;
      for (std::size_t p = 0; p < v_new.size(); ++p) {
        max_norm = std::max(max_norm, norms[p]);
        if (norms[p] > maxrotn) {
          const double scale = maxrotn / norms[p];
          // v_new[p] := scale·v_new[p] + (1−scale)·v_old[p]
          v_new[p].gaxpy(scale, v_old[p], 1.0 - scale, false);
          changed = true;
          ++n_clamped;
        }
      }
      if (changed) {
        world_.gop.fence();
        out[s].from_flat(v_new);
      }
      // Diagnostic: print whenever the cap actually fired. Tells us
      // whether step restriction is the binding constraint on the
      // iteration speed.
      if (diag_level >= 1 && n_clamped > 0 && world_.rank() == 0) {
        printf("[STEP-REST] iter=%d state=%zu clamped=%d/%zu  "
               "max_unrestricted=%.3e  cap=%.3e\n",
               calls_consumed_, s, n_clamped, v_new.size(),
               max_norm, maxrotn);
        fflush(stdout);
      }
    }
  }

  madness::World               &world_;
  ConvergencePolicy             policy_;
  std::vector<PerStateHistory>  states_;
  bool                          first_kain_call_ = true;
  int                           calls_consumed_  = 0;
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_SUBSPACE_KAIN_HPP
