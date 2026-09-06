#ifndef MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
#define MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP

// ===========================================================================
// CalcManager executor — the WORLD-BOUND half of the calc manager (doc 15,
// 15a). The scheduling core (calc_manager.hpp) decides what runs and in what
// order; this header actually drives the solves.
//
//   FdResponseExecutor : ICalcExecutor
//       run_protocol(WorkItem) -> one protocol step of an FD/NuclearFD
//       solve. Dispatches Static/Full × Closed/Open at runtime, builds the
//       perturbation source, honors the reconcile action for seed loading,
//       runs solvers::iterate_protocol over the single protocol step threshold (with
//       a reproject-prepare and a save post_step), and reports convergence.
//
//   CalcManager
//       build(n_atoms) -> dag from the plan; run(world, exec) -> the
//       one-wave-per-pass reschedule loop. Each pass reloads the aggregate
//       metadata and re-schedules, so every protocol step's action is computed against
//       up-to-date disk state (the protocol-step-level barrier the design decided on).
//
// Scope (15a): FD (dipole + nuclear) closed-shell, dipole open-shell. ES
// bundles and DerivedFD expansion are recognized but not yet solved here —
// run_protocol returns not-converged with a clear message and run() stops cleanly
// once it can make no further progress. The ES executor + derived expansion
// are the next increment; the expansion hook below is already wired so only
// the ES solve has to land.
//
// This header reuses, unchanged, exactly the machinery the FD skeleton test
// exercises: GroundState::prepare, build_response_ground_state_*,
// dipole/nuclear_perturbation, FDSolver, iterate_protocol, try_load_fd_state,
// save_fd_state.
// ===========================================================================

#include "calc_manager.hpp"

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/full.hpp"
#include "../kernels/static.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/fd_problem.hpp"
#include "../solvers/fd_save_load.hpp"
#include "../kernels/beta.hpp"
#include "../kernels/tpa_source_spec.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/es_analysis.hpp"
#include "../solvers/es_save_load.hpp"
#include "../solvers/es_seed_guard.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/fd_solver.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/node_subworlds.hpp"
#include "../solvers/response_metadata.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/vbc_save_load.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// Plain, World-free executor settings — everything the *caller* specifies
/// (convergence policy, ES/seed/accept knobs). Carried by value in
/// `ResponseWorkflowInput` (R0a) so the public contract holds no World-bound
/// state. `ExecutorContext` inherits these, so every `ctx.<knob>` access in this
/// file resolves to an inherited member and is unchanged by the split.
struct ExecutorSettings {
  ConvergencePolicy policy;
  PrintLevel        print_level = PrintLevel::Normal;
  std::string       calc_dir;      // holds response_metadata.json + archives
  int               max_iters = 25;
  // Cross-type seed toggle (optional): a derived FD starts from its converged
  // ES-root vector rather than the perturbation source. Default OFF — measured
  // neutral-to-slightly-worse at ωₙ/2 (off-resonance, H2); the ES root is the
  // right guess only for a future AT-resonance (factor 1.0) derived FD. FUTURE:
  // make the seed selectable per node, and allow a MIXTURE of roots to target
  // in-between frequencies (not just a single root's vector).
  bool              seed_derived_from_es_root = false;
  // Excited-state (Full / TDA-warmup) solve settings — defaults for the Full
  // closed-shell path (random guess, 10 warmup iters, oversampled warmup, KAIN).
  ESGuessMode       es_guess              = ESGuessMode::SolidHarmonics;  // sweep-validated default
  // AO basis projected for the VirtualAO guess (CLI --es-guess-basis / deck
  // response.excited.guess_basis). Radial rank per l-sector = #shells(l) −
  // #occupied(l), so a larger basis reaches further up the radial ladder
  // (Rydberg series want d-aug-cc-pvtz or similar). Ignored by other modes.
  std::string       es_guess_basis        = "aug-cc-pvdz";
  // --tpa-residue: contract the 2PA moment with the corrected single-residue
  // form (tpa::tpa_moment_residue — X_f in the residue slot against V^{bc}
  // built from the two photon responses) instead of the legacy beta-reuse
  // candidate. --tpa-prefactor scales the residue moment (normalization C_N).
  // DEFAULT TRUE since 2026-09-05: running WITHOUT this flag silently used the
  // legacy beta-reuse contraction, which is the "reuse V^BC and negate" form —
  // measured 3-11% below the validated composition (and DALTON) on H2O
  // per-root moments (reports/2026-09-04_report_pq_vbc_reconciliation). The
  // legacy form is kept behind --tpa-legacy as the measured comparison arm.
  bool              tpa_residue           = true;
  // Contraction composition for the two-electron residue when tpa_residue:
  // false (DEFAULT) = the Parker-style state-free (P,Q) source via the spec
  // engine (kernels/tpa_source_spec.hpp, symmetrized builder) — the form that
  // matches V^BC term-for-term (one dagger + one reorientation). true = the
  // c-grouped composition (tpa_e3_residue) — gate-equal (test_tpa_pq_vs_vbc,
  // 3e-7) and kept as the production comparison arm (--tpa-cgrouped).
  bool              tpa_cgrouped          = false;
  double            tpa_prefactor         = 1.0;
  // --tpa-decompose: also compute the zero-operator (pure two-electron E3)
  // variant of the residue and print the per-element E3/1e split + the
  // phase- and normalization-invariant fraction f = E3/total (compare vs the
  // patched-DALTON 'E3 CONTRIBUTION TO SMOM' ledger, TPA_SCOPING §5n).
  bool              tpa_decompose         = false;
  // --tpa-diag-only: print the cross-code diagnostics (norms, transition
  // dipoles, alpha(w/2)) from the loaded states and SKIP the contraction and
  // all metadata writes — a fast read-only pass over a converged calc dir.
  bool              tpa_diag_only         = false;
  // --tpa-roots=0,2 : restrict the TPA loop to these root indices (0-based).
  // Enables poor-man's macrotasking: N independent single-rank processes, one
  // per root, on one node. When set, the metadata/property write is SKIPPED
  // (read-only pass) so concurrent per-root processes cannot race.
  std::vector<int>  tpa_roots;
  int               es_tda_warmup_iters   = 10;
  // Warmup oversampling: keep the lowest n_roots of ceil(factor*n_roots)
  // partially-converged warmup trials. 3.0 (was 2.0) because on the PRODUCTION
  // ladder (moldft climbing to 1e-6/k8) the extra cut-line margin lets the
  // downselect keep h2o's true 4th state 0.4097 au and match Dalton
  // d-aug-cc-pVQZ to 0.03% (ladder job 2084168), vs a 2.4% missing-state error
  // at 2.0 (job 2083837, which reported the 5th state 0.420 in its place).
  // KNOWN LIMITATION (follow-up): the underlying cause is the TDA warmup
  // MISORDERING the near-degenerate 0.40/0.41 cluster at the coarse rung — a
  // 1e-4/k6-only run still misorders it (the true 4th sits just above the
  // 4-root cut) and does not fully tighten that root in 25 iters. Oversampling
  // only buys margin, it can't fix an ordering inversion at the cut; the real
  // fix is an energy-ordered guess (ESGuessMode::VirtualAO, the NWChem
  // CIS-diagonal path) or more warmup iters. Cost of 3.0 is warmup-only
  // (coarse rung, es_tda_warmup_iters iters).
  double            es_warmup_oversample  = 3.0;
  int               es_kain_maxsub        = 8;
  double            es_maxrotn            = 0.5;
  // Delay KAIN onset in the MAIN ES solve by this many iters (pure BSH +
  // step-restriction first, so roots stabilize before KAIN starts recording
  // history). 0 = KAIN from iter 1 (previous behaviour). Distinct from
  // es_tda_warmup_iters, which is the separate oversampled-warmup pre-pass.
  int               es_main_kain_delay    = 0;
  // Full-deflation locking of converged ES roots (workstream B). ON by default
  // (sweep-validated: required to converge near-degenerate h2o/c2h4 clusters).
  bool              es_lock_converged     = true;
  // Cache the (promoted) TDA warmup guess to a side bundle es_warmup__<pkey> and
  // reuse it on later runs, so iterating on the Full main solve skips the warmup.
  // Opt-in; the main solve never overwrites the cache (stable fixed start point).
  bool              es_warmup_cache       = false;
  // Persistent location for the warmup cache (empty = ctx.calc_dir). Set this to
  // a fixed per-fixture path so the cached guess survives across fresh calc dirs
  // (cm_es makes a new timestamped calc dir each run). Setting it implies caching.
  std::string       es_warmup_cache_dir;
  // Inc-3c tensor-layer ES γ (--es-tensor; TDA/ClosedShell only): the bundle's
  // exchange runs off the per-protocol cached g0 tensor + per-root contractions
  // (ESSolver::set_gamma_tensor / tda_batch::compute_gamma_flat) instead of one
  // Exchange operator per root per iter. A/B-to-floor vs the reference (~1e-3
  // rel in θ → ~1e-6 in a converged ω), NOT bitwise. Off by default.
  bool              es_gamma_tensor   = false;
  // Best-effort acceptance at maxiter (FD nodes). When true, a non-diverged FD
  // solve that exhausts max_iters WITHOUT meeting the strict target is recorded
  // converged (with an `accepted` marker + its real residual). NOTE: ladder
  // CLIMBING no longer needs this — reconcile's honest-climb rule (a budget-
  // exhausted rung Skips forward with `converged` honestly false; see
  // calc_manager.hpp reconcile_protocol) is the default. This flag's remaining
  // purpose is to FORCE `converged=true` so the downstream VBC/property
  // prerequisite gates unblock on a stubborn channel (e.g. the cusped
  // nuclear-displacement Raman FD). Off by default (honest verdicts). The finest
  // protocol in --protocol is the de-facto "final" rung; its recorded
  // bsh_residual is the verdict either way.
  bool              accept_at_maxiter = false;
  // F2 (doc 32 §5): fan independent states out across subworlds. This is the
  // number of subworlds PER PHYSICAL NODE (F2f): 0 = single-World reference
  // (byte-identical); 1 = one subworld per node (node-aligned); >1 = sub-node /
  // NUMA packing (φ replicated P× per node, more state-parallelism — small-system
  // regime). Total subworlds G = n_nodes × this. G≤1 short-circuits to the G=0
  // path. Fans FD / NuclearFD / VBC (F2g); ES stays single-World.
  int               fd_subworlds = 0;
  // F2 (doc 32 §5.3): when non-empty, FD metadata writes go to a per-group shard
  // response_metadata.group<tag>.json instead of the canonical file, so concurrent
  // subworlds never race it; universe rank 0 merges the shards after the fence.
  // "" = write the canonical file directly (the default / single-World path).
  std::string       metadata_shard;
  // F2d (doc 32 §5.6): per-subworld log attribution. Both are the EMPTY/-1
  // sentinel in single-World mode ⇒ output byte-identical to the G=0 reference.
  // Set by the fan-out closure. A non-empty log_prefix is ALSO the "I am in a
  // subworld" signal that suppresses the redundant GS/protocol banners (F2d-ii).
  std::string       log_prefix;     // "[g2/3 xm005] " prepended to per-state human lines
  int               log_group = -1; // node index → MEMORY_HWM `group=` field (-1 = omit)
};

/// Everything an FD protocol step solve needs beyond the node itself. The
/// World-bound part (`world` / `gs` borrowed — the executor lives within the
/// solve scope) plus the plain `ExecutorSettings` (inherited). The orchestrator
/// builds this from a loaded GroundState + the Input's settings.
struct ExecutorContext : ExecutorSettings {
  madness::World   &world;
  GroundState      &gs;            // mutable: re-prepared at each protocol
  double            L = 0.0;       // cubic cell half-edge (from archive header)
  std::string       fock_json;     // moldft.fock.json path ("" -> none)

  ExecutorContext(madness::World &w, GroundState &g, double L_,
                  std::string fock, ExecutorSettings s = {})
      : ExecutorSettings(std::move(s)), world(w), gs(g), L(L_),
        fock_json(std::move(fock)) {}
};

// ---------------------------------------------------------------------------
namespace detail_exec {

inline bool is_static_freq(double freq) { return std::abs(freq) < 1e-12; }

/// Parse a derived-FD provenance label "es_root_NNNN" -> NNNN; -1 otherwise.
inline int es_root_index(const std::string &label) {
  const std::string pre = "es_root_";
  if (label.rfind(pre, 0) != 0) return -1;
  try { return std::stoi(label.substr(pre.size())); } catch (...) { return -1; }
}

/// Load a converged first-order state as an (X,Y) pair for the VBC build.
/// Dynamic (freq != 0): the Full FD state directly. Static (freq == 0): the
/// Static FD state with Y = X (the static limit). nullopt if not on disk.
template <class Shell>
inline std::optional<ResponseStateXY<Shell>>
load_fd_as_xy(madness::World &world, const std::string &calc_dir,
              const Perturbation &pert, double freq) {
  if (is_static_freq(freq)) {
    auto r = try_load_fd_state<Static, Shell>(world, calc_dir, pert, freq);
    if (!r) return std::nullopt;
    ResponseStateXY<Shell> xy;
    xy.x_alpha = madness::copy(world, r->state.responses[0].x_alpha);
    xy.y_alpha = madness::copy(world, r->state.responses[0].x_alpha);  // static: Y = X
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      xy.x_beta = madness::copy(world, r->state.responses[0].x_beta);
      xy.y_beta = madness::copy(world, r->state.responses[0].x_beta);
    }
    return xy;
  }
  auto r = try_load_fd_state<Full, Shell>(world, calc_dir, pert, freq);
  if (!r) return std::nullopt;
  return r->state.responses[0];
}

/// Alpha-spin perturbation source vector for one (pert) at the active
/// protocol. Throws for sources 15a does not build.
inline vector_real_function_3d
pert_source_alpha(madness::World &world, GroundState &gs,
                  const Perturbation &pert) {
  switch (pert.kind) {
    case Perturbation::Kind::Dipole:
      return dipole_perturbation(world, gs, pert.axis);
    case Perturbation::Kind::NuclearDisplacement:
      return nuclear_perturbation(world, gs, pert.atom, pert.axis);
    case Perturbation::Kind::Magnetic:
      throw std::runtime_error(
          "calc executor: magnetic perturbation source not implemented (15a)");
  }
  throw std::runtime_error("calc executor: unknown perturbation kind");
}

/// Beta-spin source (open shell). Only dipole has a beta builder today.
inline vector_real_function_3d
pert_source_beta(madness::World &world, GroundState &gs,
                 const Perturbation &pert) {
  if (pert.kind == Perturbation::Kind::Dipole)
    return dipole_perturbation_beta(world, gs, pert.axis);
  throw std::runtime_error(
      "calc executor: open-shell source for this perturbation not "
      "implemented (15a) — only dipole has a beta builder");
}

/// Fill the FD source carrier (== the solver Storage shape) for (Type, Shell).
template <typename Type, typename Shell, typename Carrier>
void fill_source(madness::World &world, GroundState &gs,
                 const Perturbation &pert, Carrier &src) {
  const auto pa = pert_source_alpha(world, gs, pert);
  src.x_alpha = pa;
  if constexpr (std::is_same_v<Type, Full>) src.y_alpha = pa;
  if constexpr (std::is_same_v<Shell, OpenShell>) {
    const auto pb = pert_source_beta(world, gs, pert);
    src.x_beta = pb;
    if constexpr (std::is_same_v<Type, Full>) src.y_beta = pb;
  }
}

/// Build the (Type, Shell)-correct GS-side kernel inputs.
template <typename Shell>
ResponseGroundState build_gs(madness::World &world, GroundState &gs) {
  if constexpr (std::is_same_v<Shell, ClosedShell>)
    return build_response_ground_state_closed_shell(
        world, gs, gs.hf_exchange_coefficient(), gs.params().lo());
  else
    return build_response_ground_state_open_shell(
        world, gs, gs.hf_exchange_coefficient(), gs.params().lo());
}

/// Re-project every in-flight component to a new wavelet basis. Generic over
/// (Type, Shell) via Storage::blocks().
template <typename State>
void reproject_state(State &st, int k, double thresh) {
  for (auto &ch : st.responses)
    for (auto *blk : ch.blocks())
      for (auto &fn : *blk) fn = madness::project(fn, k, thresh);
  for (auto &rho : st.rho_alpha_prev)
    rho = madness::project(rho, k, thresh);
}

} // namespace detail_exec

// ---------------------------------------------------------------------------
// One protocol step of one FD node.
// ---------------------------------------------------------------------------

/// Solve (pert, freq) at a single protocol `thresh`. `action` is the reconcile
/// verdict: Fresh starts from the perturbation guess (never loads, so a
/// diverged archive can't poison the solve); Restart/Resume load the nearest
/// converged-or-partial seed via try_load_fd_state (which re-projects a
/// coarser source to the active key). Saves the result (+ metrics) through
/// save_fd_state.
template <typename Type, typename Shell>
NodeResult solve_fd_protocol(ExecutorContext &ctx, const Perturbation &pert,
                         double freq, double thresh, NodeAction action,
                         const std::string &es_root_id = std::string()) {
  using namespace madness;
  using Solver = FDSolver<Type, Shell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  const double step_w0 = madness::wall_time();   // R1b: point-solve wall timer

  // Bring FunctionDefaults<3> + the ground state to this protocol so the
  // source builders and try_load (which keys on protocol_key()) are correct.
  const bool gs_verbose = ctx.log_prefix.empty();  // F2d-ii: quiet in subworlds
  set_response_protocol(world, ctx.L, thresh, /*override_k=*/-1, gs_verbose);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json, gs_verbose);
  }

  // Target at this protocol (rebuilt inside prepare for the ramp's single step).
  FDProblem<Type, Shell> tgt;
  tgt.gs = detail_exec::build_gs<Shell>(world, gs);
  tgt.responses.resize(1);
  tgt.responses[0].omega = freq;
  detail_exec::fill_source<Type, Shell>(world, gs, pert,
                                        tgt.responses[0].source);

  // Initial-guess precedence (the source carrier IS a valid Storage):
  //   1. disk FD partial at this (pert, freq) — restart (skipped if Fresh);
  //   2. else the converged ES-root vector if this is a derived FD and the
  //      cross-type seed is enabled (optional, ctx.seed_derived_from_es_root);
  //   3. else the perturbation source.
  typename Solver::State s0;
  s0.responses.resize(1);
  // Fresh guess x0 = 0 (matches molresponse_v2, which zero-allocates the response):
  // the first BSH step then gives G(-2*V_p*phi), the bounded uncoupled first-order
  // response. Seeding x0 = the source directly blows up iter-1 for a large/peaked
  // source (nuclear displacement) — V0*x0 and the BSH shift*x0 are O(||source||),
  // which for nuclear exceeds the explosion guard. Same converged fixed point
  // either way (linear response is unique), so alpha/beta values are unchanged;
  // this only fixes the transient / divergence.
  // DEEP-copy the source's shape then zero it. A shallow copy (operator=) shares
  // the Function handles, so scaling it would also zero the perturbation source
  // itself -> no driving term -> x stays 0 (the "converged in 1 iter, res=0,
  // value=0" bug). madness::copy gives an independent buffer.
  s0.responses[0] = tgt.responses[0].source;    // shape (shallow)
  {
    auto zflat = madness::copy(world, s0.responses[0].flatten());
    madness::scale(world, zflat, 0.0);
    s0.responses[0].from_flat(zflat);           // independent zeros
  }
  const char *seed_kind = "zero";

  bool seeded = false;
  if (action != NodeAction::Fresh) {
    auto loaded =
        try_load_fd_state<Type, Shell>(world, ctx.calc_dir, pert, freq);
    if (loaded) {
      s0 = std::move(loaded->state);
      seeded = true;
      seed_kind = "fd_restart";
    }
  }
  // Cross-type seed: a derived FD (es_root_id set) starts from its converged
  // ES-root vector. The root is ResponseStateX (X only); for FD Full we seed
  // x = y = X. Only Full/ClosedShell derived FDs use it, and only when enabled.
  if constexpr (std::is_same_v<Type, Full> &&
                std::is_same_v<Shell, ClosedShell>) {
    if (!seeded && ctx.seed_derived_from_es_root && !es_root_id.empty()) {
      const int idx = detail_exec::es_root_index(es_root_id);
      if (idx >= 0) {
        auto esb = try_load_es_bundle<TDA, ClosedShell>(world, ctx.calc_dir);
        if (esb && idx < static_cast<int>(esb->state.roots.size())) {
          s0.responses[0].x_alpha =
              madness::copy(world, esb->state.roots[idx].x_alpha);
          s0.responses[0].y_alpha =
              madness::copy(world, esb->state.roots[idx].x_alpha);
          seeded = true;
          seed_kind = "es_root";
        }
      }
    }
  }

  Solver solver(world, tgt, ctx.policy, ctx.print_level, ctx.log_prefix);  // F2d tag

  // The protocol + ground state + target are already set up above for `thresh`.
  // iterate_protocol calls prepare() before the (single) step; re-doing the
  // expensive gs.prepare()/build_gs() there would DOUBLE the per-step prep cost
  // (ground-state projection dominates — ~half of wall time per solve). So only
  // re-prepare if the protocol actually changed; otherwise just re-project the
  // (possibly coarser-loaded) state into the active basis.
  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, typename Solver::State &st) {
    const double cur_t = FunctionDefaults<3>::get_thresh();
    const bool changed =
        std::abs(th - prepared_t) > 1e-15 * std::max(th, prepared_t) ||
        std::abs(cur_t - th)      > 1e-15 * std::max(cur_t, th);
    if (changed) {
      set_response_protocol(world, ctx.L, th, /*override_k=*/-1, gs_verbose);
      const double new_t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
      gs.prepare(world, 0.001 * new_t, coulop, ctx.fock_json, gs_verbose);

      FDProblem<Type, Shell> nt;
      nt.gs = detail_exec::build_gs<Shell>(world, gs);
      nt.responses.resize(1);
      nt.responses[0].omega = freq;
      detail_exec::fill_source<Type, Shell>(world, gs, pert,
                                            nt.responses[0].source);
      solv.set_target(std::move(nt));
      prepared_t = new_t;
    }
    detail_exec::reproject_state(st, FunctionDefaults<3>::get_k(),
                                 FunctionDefaults<3>::get_thresh());
  };

  // Real convergence = not diverged AND residuals within this protocol's
  // targets. A max-iters STALL must read as unconverged so it is neither
  // skipped on restart nor fed to property assembly as if good.
  auto converged_now = [](const typename Solver::State &st, const Solver &sv) {
    if (st.diverged) return false;
    double mb = 0.0, md = 0.0;
    for (double r : st.last_bsh_residual)     mb = std::max(mb, r);
    for (double r : st.last_density_residual) md = std::max(md, r);
    const auto &t = sv.targets();
    return mb <= t.bsh_residual && md <= t.density_residual;
  };

  // Best-effort acceptance (--accept-at-maxiter): a non-diverged solve that
  // exhausts max_iters without meeting the strict target is accepted so the
  // node climbs the ladder / unblocks VBC instead of stalling. The `accepted`
  // flag + recorded residual keep the verdict honest.
  auto accepted_now = [&](const typename Solver::State &st, const Solver &sv) {
    return ctx.accept_at_maxiter && !st.diverged &&
           st.iter >= ctx.max_iters && !converged_now(st, sv);
  };

  auto post_step = [&](double, Solver &solv, typename Solver::State &st) {
    const bool strict   = converged_now(st, solv);
    const bool accepted = accepted_now(st, solv);
    const double wall_s = madness::wall_time() - step_w0;   // R1b
    save_fd_state<Type, Shell>(world, st, ctx.calc_dir, pert, freq,
                               /*converged=*/strict || accepted,
                               /*seed=*/seed_kind, /*accepted=*/accepted,
                               /*wall_s=*/wall_s,
                               /*metadata_shard=*/ctx.metadata_shard,
                               /*log_prefix=*/ctx.log_prefix,
                               /*log_group=*/ctx.log_group);
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  const bool strict   = converged_now(sf, solver);
  const bool accepted = accepted_now(sf, solver);
  r.converged = strict || accepted;
  r.reached_protocol_key = protocol_key();  // active defaults reflect this protocol step
  if (world.rank() == 0) {
    const char *acc = accepted ? " (ACCEPTED best-effort @ maxiter — strict target "
                                 "NOT met; see bsh_residual)" : "";
    // F2d: prepend the subworld tag (empty ⇒ unchanged, G=0 byte-identical).
    if (ctx.log_prefix.empty())
      madness::print("[CALC] fd solve: pert=", pert.description(), " freq=", freq,
                     " thresh=", thresh, " seed=", seed_kind, " iters=", sf.iter,
                     " converged=", r.converged, acc);
    else
      madness::print(ctx.log_prefix,
                     "[CALC] fd solve: pert=", pert.description(), " freq=", freq,
                     " thresh=", thresh, " seed=", seed_kind, " iters=", sf.iter,
                     " converged=", r.converged, acc);
  }
  return r;
}

// ---------------------------------------------------------------------------
// Excited-state bundle solve (15b-i: TDA closed-shell only).
// ---------------------------------------------------------------------------

/// Solve a TDA closed-shell excited-state bundle of `n_roots` at one protocol
/// step. Mirrors solve_fd_protocol: set protocol -> prepare GS -> build problem
/// -> guess (Fresh = fresh solid-harmonics trial; else nearest on-disk bundle
/// via the restart precedence) -> iterate_protocol over the single step (guarded
/// set_gs/reproject prepare + save_es_roots post-step). The bundle saves to a
/// per-protocol subdir `es__<key>` under calc_dir so lower-protocol restart
/// works; the converged excitation energies are persisted to
/// response_metadata.json (excited_states/<key>/roots[].omega), which is what
/// DerivedFD expansion (expand_converged_es) reads — making expansion restart-safe.
inline NodeResult solve_es_tda_closed_shell(ExecutorContext &ctx, int n_roots,
                                            double thresh, NodeAction action) {
  using namespace madness;
  using Solver = ESSolver<TDA, ClosedShell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  const double step_w0 = madness::wall_time();   // R1b: bundle-solve wall timer

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto problem = build_es_problem_tda<ClosedShell>(world, gs, n_roots, c_xc, lo);

  // Warmup policy (KAIN on, kain_maxsub/maxrotn from ctx); the main solve
  // disables the warmup window — the fresh guess is an oversampled TDA warmup
  // (>= 2x roots by default), same recipe as the Full path.
  ConvergencePolicy warm_policy = ctx.policy;
  warm_policy.kain                     = true;
  warm_policy.kain_maxsub              = ctx.es_kain_maxsub;
  warm_policy.maxrotn                  = ctx.es_maxrotn;
  warm_policy.tda_warmup_iters         = ctx.es_tda_warmup_iters;
  warm_policy.warmup_oversample_factor = ctx.es_warmup_oversample;
  ConvergencePolicy main_policy = warm_policy;
  main_policy.tda_warmup_iters = ctx.es_main_kain_delay;
  main_policy.lock_converged   = ctx.es_lock_converged;  // warmup never locks

  Solver::State s0;
  bool seeded = false;
  // Root-identity guard (W5): remember WHERE the seed came from; the guard
  // reloads it from disk at evaluation time (no bundle copy held in memory
  // during the solve — see es_seed_guard.hpp).
  std::optional<EsSeedReference> seed_ref;
  if (action != NodeAction::Fresh) {
    auto loaded = try_load_es_bundle<TDA, ClosedShell>(world, ctx.calc_dir);
    if (loaded) {
      seed_ref = EsSeedReference{ctx.calc_dir + "/" + loaded->bundle_dir,
                                 loaded->bundle_dir,
                                 loaded->source_protocol_key};
      s0 = std::move(loaded->state);
      seeded = true;
    }
  }
  if (!seeded) {
    const long n_warm = std::max<long>(
        n_roots, static_cast<long>(
                     std::ceil(ctx.es_warmup_oversample *
                               static_cast<double>(n_roots))));
    s0 = run_oversampled_tda_warmup<ClosedShell>(
        world, gs, n_roots, n_warm, ctx.es_tda_warmup_iters, warm_policy,
        c_xc, lo, ctx.print_level, ctx.es_guess, ctx.es_guess_basis);
  }

  Solver solver(world, std::move(problem), main_policy, ctx.print_level);
  solver.set_gamma_tensor(ctx.es_gamma_tensor);  // Inc-3c: tensor-layer γ gate
  // Review HIGH: the locked step variant (lock_converged, the production
  // default) uses the per-root REFERENCE γ path and ignores gamma_tensor_, so
  // --es-tensor is a silent no-op unless --no-lock-converged is also set. Warn
  // loudly rather than let a user believe the untiled-crash mitigation is
  // active when it is not. (Teaching the locked variant to use the tensor path
  // is tracked as a follow-up.)
  if (ctx.es_gamma_tensor && ctx.es_lock_converged && world.rank() == 0)
    madness::print(
        "[ES-TENSOR] WARNING: --es-tensor was requested but lock_converged is "
        "on (the default) — the locked step uses the REFERENCE gamma path, so "
        "the tensor (Inc-3c) path is NOT active. Pass --no-lock-converged to "
        "actually use --es-tensor (e.g. to avoid the untiled n_occ>=34 crash).");

  // Single-step guard (as in solve_fd_protocol): re-prepare the ground state
  // only on an actual protocol change; otherwise just re-project the roots.
  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, Solver::State &st) {
    const double cur_t = FunctionDefaults<3>::get_thresh();
    const bool changed =
        std::abs(th - prepared_t) > 1e-15 * std::max(th, prepared_t) ||
        std::abs(cur_t - th)      > 1e-15 * std::max(cur_t, th);
    if (changed) {
      set_response_protocol(world, ctx.L, th);
      const double new_t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
      gs.prepare(world, 0.001 * new_t, coulop, ctx.fock_json);
      solv.set_gs(build_response_ground_state_closed_shell(
          world, gs, gs.hf_exchange_coefficient(), gs.params().lo()));
      prepared_t = new_t;
    }
    const int    k  = FunctionDefaults<3>::get_k();
    const double tt = FunctionDefaults<3>::get_thresh();
    for (auto &root : st.roots)
      for (auto &fn : root.x_alpha) fn = madness::project(fn, k, tt);
    for (auto &rho : st.rho_alpha_prev) rho = madness::project(rho, k, tt);
  };

  // Report/save convergence with the SAME criterion the iteration stops on
  // (ESSolver::converged = energy+density, not the jittering BSH amplitude).
  auto converged_now = [](const Solver::State &st, const Solver &sv) {
    return !st.diverged && sv.converged(st);
  };

  auto post_step = [&](double, Solver &solv, Solver::State &st) {
    const std::string dir = ctx.calc_dir + "/es__" + protocol_key();
    const double wall_s = madness::wall_time() - step_w0;   // R1b
    save_es_roots<TDA, ClosedShell>(world, st, dir,
                                    /*converged=*/converged_now(st, solv),
                                    /*wall_s=*/wall_s);
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();
  // Root-identity guard (W5): visibility on what a SEEDED solve actually did
  // (seed overlap / ω shift per root, basin-escape + pure-tracking warnings)
  // Collective evaluate; rank-0 print + metadata write. Fresh (unseeded)
  // solves skip the guard entirely.
  if (seed_ref) {
    auto guard = evaluate_es_seed_guard<TDA, ClosedShell>(world, sf, r.converged, *seed_ref);
    print_es_seed_guard(world, guard);
    record_es_seed_guard(world, ctx.calc_dir, protocol_key(), guard);
  }
  // Post-convergence transition-property report (legacy TDDFT::analysis +
  // analyze_vectors). Runs on the IN-MEMORY converged state `sf` at the solve's
  // own process count and writes only a rank-0 JSON — it does NOT reload the
  // bundle, so it never hits the cross-np load path that caused the parked ES
  // heap-OOB (now guarded in load_es_roots). Re-enabled with that guard in place.
  // Ungated (was print_level-gated): oscillator strengths / transition moments
  // are part of the machine-readable property set — computing them and throwing
  // them away at low print levels made the ES output incomplete (madqc review).
  if (r.converged)
    report_es_analysis<TDA, ClosedShell>(
        world, gs, sf, ctx.print_level,
        ctx.calc_dir + "/es_analysis__" + protocol_key() + ".json");
  // NB: ES roots are persisted to response_metadata.json via save_es_roots;
  // DerivedFD expansion reads them from there (expand_converged_es), NOT from
  // this return value (run()'s loop discards it). So no in-memory copy here.
  return r;
}

// ---------------------------------------------------------------------------
// Excited-state bundle solve — Full (paired X,Y) closed-shell (15b-ii).
// ---------------------------------------------------------------------------

/// Solve a Full (paired X,Y) closed-shell ES bundle of `n_roots` at one protocol
/// step via the direct ESSolver<Full, ClosedShell>. Fresh guess: a random-guess
/// TDA warmup (es_tda_warmup_iters, oversampled by es_warmup_oversample)
/// promoted to (X, Y=0) — the Full solver couples X and Y from iter 1. Restart:
/// load the (X,Y) bundle directly (no promote). Saves a proper "full" (X,Y)
/// bundle. Convergence is residual-based (same fields as TDA).
inline NodeResult solve_es_full_closed_shell(ExecutorContext &ctx, int n_roots,
                                             double thresh, NodeAction action) {
  using namespace madness;
  using Solver = ESSolver<Full, ClosedShell>;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  const double step_w0 = madness::wall_time();   // R1b: bundle-solve wall timer

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }

  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto problem = build_es_problem_full<ClosedShell>(world, gs, n_roots, c_xc, lo);

  // Warmup policy (KAIN on, kain_maxsub / maxrotn from ctx); the main solve
  // disables the TDA warmup window (it was done explicitly for the fresh guess).
  ConvergencePolicy warm_policy = ctx.policy;
  warm_policy.kain                     = true;
  warm_policy.kain_maxsub              = ctx.es_kain_maxsub;
  warm_policy.maxrotn                  = ctx.es_maxrotn;
  warm_policy.tda_warmup_iters         = ctx.es_tda_warmup_iters;
  warm_policy.warmup_oversample_factor = ctx.es_warmup_oversample;
  ConvergencePolicy main_policy = warm_policy;
  main_policy.tda_warmup_iters = ctx.es_main_kain_delay;
  main_policy.lock_converged   = ctx.es_lock_converged;  // warmup never locks

  Solver::State s0;
  bool seeded = false;
  // Root-identity guard (W5): deep-copy the seed for the post-solve overlap /
  // ω-shift report (see es_seed_guard.hpp). NB: the cached TDA warmup guess
  // below is NOT a seed in this sense (it is this solver's own cold-start
  // artifact), so seed_ref stays empty on that path.
  std::optional<EsSeedReference> seed_ref;
  if (action != NodeAction::Fresh) {
    auto loaded = try_load_es_bundle<Full, ClosedShell>(world, ctx.calc_dir);
    if (loaded) {
      seed_ref = EsSeedReference{ctx.calc_dir + "/" + loaded->bundle_dir,
                                 loaded->bundle_dir,
                                 loaded->source_protocol_key};
      s0 = std::move(loaded->state);
      seeded = true;
    }
  }
  if (!seeded) {
    const std::string cache_base =
        ctx.es_warmup_cache_dir.empty() ? ctx.calc_dir : ctx.es_warmup_cache_dir;
    const std::string warm_cache = cache_base + "/es_warmup__" + protocol_key();
    bool used_cache = false;
    if (ctx.es_warmup_cache) {
      int have = 0;
      if (world.rank() == 0)
        have = std::filesystem::exists(warm_cache + "/roots.json") ? 1 : 0;
      world.gop.broadcast(have, 0);
      if (have) {
        try {
          s0 = load_es_roots<Full, ClosedShell>(world, warm_cache);
          used_cache = true;
          if (world.rank() == 0)
            madness::print("[CALC] solve_es_full: reusing cached TDA warmup guess"
                           " from", warm_cache, "(skipping warmup)");
        } catch (const std::exception &e) {
          if (world.rank() == 0)
            madness::print("[CALC] solve_es_full: cached warmup unusable (",
                           e.what(), ") — running a fresh warmup");
        }
      }
    }
    if (!used_cache) {
      const long n_warm = std::max<long>(
          n_roots, static_cast<long>(
                       std::ceil(ctx.es_warmup_oversample *
                                 static_cast<double>(n_roots))));
      auto tda = run_oversampled_tda_warmup<ClosedShell>(
          world, gs, n_roots, n_warm, ctx.es_tda_warmup_iters, warm_policy,
          c_xc, lo, ctx.print_level, ctx.es_guess, ctx.es_guess_basis);
      s0 = promote_tda_to_full_closed_shell(world, tda);
      if (ctx.es_warmup_cache) {
        save_es_roots<Full, ClosedShell>(world, s0, warm_cache,
                                         /*converged=*/false, /*wall_s=*/0.0,
                                         /*register_aggregate=*/false);
        if (world.rank() == 0)
          madness::print("[CALC] solve_es_full: cached TDA warmup guess ->",
                         warm_cache);
      }
    }
  }

  Solver solver(world, std::move(problem), main_policy, ctx.print_level);

  double prepared_t = t0;
  auto prepare = [&](double th, Solver &solv, Solver::State &st) {
    const double cur_t = FunctionDefaults<3>::get_thresh();
    const bool changed =
        std::abs(th - prepared_t) > 1e-15 * std::max(th, prepared_t) ||
        std::abs(cur_t - th)      > 1e-15 * std::max(cur_t, th);
    if (changed) {
      set_response_protocol(world, ctx.L, th);
      const double new_t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
      gs.prepare(world, 0.001 * new_t, coulop, ctx.fock_json);
      solv.set_gs(build_response_ground_state_closed_shell(world, gs, c_xc, lo));
      prepared_t = new_t;
    }
    const int    k  = FunctionDefaults<3>::get_k();
    const double tt = FunctionDefaults<3>::get_thresh();
    for (auto &root : st.roots)
      for (auto *blk : root.blocks())
        for (auto &fn : *blk) fn = madness::project(fn, k, tt);
    for (auto &rho : st.rho_alpha_prev) rho = madness::project(rho, k, tt);
  };

  // Report/save convergence with the SAME criterion the iteration stops on
  // (ESSolver::converged = energy+density, not the jittering BSH amplitude).
  auto converged_now = [](const Solver::State &st, const Solver &sv) {
    return !st.diverged && sv.converged(st);
  };

  auto post_step = [&](double, Solver &solv, Solver::State &st) {
    const std::string dir = ctx.calc_dir + "/es__" + protocol_key();
    const double wall_s = madness::wall_time() - step_w0;   // R1b
    save_es_roots<Full, ClosedShell>(world, st, dir,
                                     /*converged=*/converged_now(st, solv),
                                     /*wall_s=*/wall_s);
  };

  solvers::IterateProtocolPolicy pp;
  pp.max_iters_per_step = ctx.max_iters;
  const std::vector<double> one_protocol = {thresh};
  auto sf = solvers::iterate_protocol(solver, std::move(s0), one_protocol,
                                      prepare, post_step, pp);

  NodeResult r;
  r.converged = converged_now(sf, solver);
  r.reached_protocol_key = protocol_key();
  // Root-identity guard (W5) — same contract as the TDA path above.
  if (seed_ref) {
    auto guard = evaluate_es_seed_guard<Full, ClosedShell>(world, sf, r.converged, *seed_ref);
    print_es_seed_guard(world, guard);
    record_es_seed_guard(world, ctx.calc_dir, protocol_key(), guard);
  }
  // Post-convergence transition-property report (legacy TDDFT::analysis +
  // analyze_vectors). Runs on the in-memory `sf` (no bundle reload), so it never
  // hits the cross-np load path that caused the parked ES heap-OOB (now guarded
  // in load_es_roots). See the TDA path above for the full note. Re-enabled.
  if (r.converged && ctx.print_level >= PrintLevel::Normal)
    report_es_analysis<Full, ClosedShell>(
        world, gs, sf, ctx.print_level,
        ctx.calc_dir + "/es_analysis__" + protocol_key() + ".json");
  // NB: ES roots are persisted to response_metadata.json via save_es_roots;
  // DerivedFD expansion reads them from there (expand_converged_es), NOT from
  // this return value (run()'s loop discards it). So no in-memory copy here.
  return r;
}

// ---------------------------------------------------------------------------
// VBC quadratic-source build (beta-ii-b, closed-shell). Not iterative: load the
// two converged first-order states at this protocol, build the source, save.
// ---------------------------------------------------------------------------

inline NodeResult
solve_vbc_closed_shell(ExecutorContext &ctx, const CalcNode &node, double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  const bool gs_verbose = ctx.log_prefix.empty();  // F2d-ii: quiet in subworlds
  set_response_protocol(world, ctx.L, thresh, /*override_k=*/-1, gs_verbose);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json, gs_verbose);
  }

  // The two converged first-order states (X,Y) at this protocol. The
  // prerequisites_converged gate guarantees they are converged before we run.
  auto B = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, node.pert,   node.freq);
  auto C = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, node.pert_c, node.freq_c);
  if (!B || !C) {
    if (world.rank() == 0)
      madness::print("[CALC] solve_vbc:", node.id,
                     "— a converged FD prerequisite is missing; skipping.");
    return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
  }

  const double step_w0 = madness::wall_time();   // R1b: VBC build wall timer
  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);

  // Raw one-electron perturbation operators for B and C (dipole only for now).
  // Generic per-kind operator (dipole -> r, nuclear -> dV/dQ) so a VBC pair can
  // mix types (Raman = dipole x nuclear).
  auto VB_op = perturbation_operator(world, gs, node.pert);
  auto VC_op = perturbation_operator(world, gs, node.pert_c);

  auto vbc_src = vbc::compute_vbc<ClosedShell>(world, g0, *B, *C, VB_op, VC_op);
  const double wall_s = madness::wall_time() - step_w0;   // R1b
  save_vbc_state<ClosedShell>(world, vbc_src, ctx.calc_dir, node.id,
                              /*converged=*/true, /*wall_s=*/wall_s,
                              /*metadata_shard=*/ctx.metadata_shard,  // F2g
                              /*log_prefix=*/ctx.log_prefix,          // F2d
                              /*log_group=*/ctx.log_group);           // F2d

  if (world.rank() == 0) {
    if (ctx.log_prefix.empty())
      madness::print("[CALC] solve_vbc: built", node.id, "at", protocol_key());
    else
      madness::print(ctx.log_prefix, "[CALC] solve_vbc: built", node.id, "at",
                     protocol_key());
  }
  return NodeResult{/*converged=*/true, protocol_key(), {}};
}

// ---------------------------------------------------------------------------
// FdResponseExecutor — dispatches the runtime (Type × Shell) for FD nodes.
// ---------------------------------------------------------------------------

class FdResponseExecutor : public ICalcExecutor {
public:
  explicit FdResponseExecutor(ExecutorContext ctx) : ctx_(ctx) {}

  NodeResult run_protocol(const WorkItem &item) override {
    const CalcNode &node = *item.node;

    if (node.kind == CalcKind::ES) {
      const bool restricted = ctx_.gs.is_spin_restricted();
      if (restricted) {
        if (ctx_.world.rank() == 0)
          madness::print("[CALC] run_protocol: ES bundle", node.id,
                         " n_roots=", node.n_roots, " thresh=", item.thresh,
                         (node.tda ? " (TDA closed-shell)"
                                   : " (Full closed-shell)"),
                         " action=", node_action_name(item.action));
        return node.tda
                   ? solve_es_tda_closed_shell(ctx_, node.n_roots, item.thresh,
                                               item.action)
                   : solve_es_full_closed_shell(ctx_, node.n_roots, item.thresh,
                                                item.action);
      }
      // Open-shell excited states are out of scope (closed-shell only; future).
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: ES bundle", node.id,
                       "on an OPEN-SHELL ground state — open-shell excited "
                       "states are OUT OF SCOPE (closed-shell only; future "
                       "research). Skipping.");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }
    if (node.kind == CalcKind::VBC) {
      if (ctx_.gs.is_spin_restricted()) {
        if (ctx_.world.rank() == 0)
          madness::print("[CALC] run_protocol: VBC", node.id,
                         " B=", node.pert.description(), "@", node.freq,
                         " C=", node.pert_c.description(), "@", node.freq_c,
                         " thresh=", item.thresh);
        return solve_vbc_closed_shell(ctx_, node, item.thresh);
      }
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: VBC", node.id,
                       "on an open-shell ground state is out of scope. Skipping.");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }
    if (node.kind == CalcKind::DerivedFD) {
      // DerivedFD should have been promoted to FD on expansion.
      if (ctx_.world.rank() == 0)
        madness::print("[CALC] run_protocol: node", node.id,
                       "kind DerivedFD should be promoted to FD — skipping");
      return NodeResult{/*converged=*/false, /*reached_protocol_key=*/"", {}};
    }

    const bool restricted = ctx_.gs.is_spin_restricted();
    const bool is_static   = detail_exec::is_static_freq(node.freq);

    if (ctx_.world.rank() == 0) {
      madness::print("[CALC] run_protocol: node", node.id,
                     " pert=", node.pert.description(),
                     " freq=", node.freq, " thresh=", item.thresh,
                     " type=", (is_static ? "static" : "full"),
                     " shell=", (restricted ? "closed" : "open"),
                     " action=", node_action_name(item.action));
    }

    if (restricted && is_static)
      return solve_fd_protocol<Static, ClosedShell>(ctx_, node.pert, node.freq,
                                                 item.thresh, item.action, node.es_root_id);
    if (restricted && !is_static)
      return solve_fd_protocol<Full, ClosedShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action, node.es_root_id);
    if (!restricted && is_static)
      return solve_fd_protocol<Static, OpenShell>(ctx_, node.pert, node.freq,
                                               item.thresh, item.action, node.es_root_id);
    return solve_fd_protocol<Full, OpenShell>(ctx_, node.pert, node.freq,
                                          item.thresh, item.action, node.es_root_id);
  }

private:
  ExecutorContext ctx_;
};

// ---------------------------------------------------------------------------
// CalcManager — top-level driver.
// ---------------------------------------------------------------------------

/// Seed-selection strategy (nearest converged; the knob is the metric, not a
/// fixed topology — see doc 15). Hoisted to namespace scope so it is a
/// complete type when used as a `= {}` default argument inside CalcManager.
enum class SeedStrategy { NearestConverged };
struct CalcManagerPolicy {
  SeedStrategy seed = SeedStrategy::NearestConverged;
  int          max_iters_per_step = 25;
  int          fd_subworlds = 0;   // F2f: subworlds PER NODE (0=off, 1=node-aligned)
};

// ===========================================================================
// Execution model (15a — single group). WHY all ranks must agree on the plan.
//
// A MADNESS Function is distributed across ALL ranks of the communicator, so
// solving one response state is a COLLECTIVE operation: every rank holds a
// piece and they collaborate on every step (BSH apply, fences, ...). 15a runs
// single-group — the whole communicator solves ONE state at a time; a "wave"
// of independent states is simply looped over (all ranks per state), NOT fanned
// out across ranks. Distributing a wave's states to rank SUBGROUPS sized to fit
// memory is the 15c design (STATE_PARALLEL / subworlds); the wave structure
// exists precisely so that later layer can partition it.
//
// Consequence — DETERMINISM. Every rank runs schedule() independently, then all
// ranks collectively solve waves.front()[i] together. If two ranks computed
// different waves they would issue mismatched collective calls and DEADLOCK.
// So the schedule must come out byte-identical on every rank: same
// response_metadata.json in (all ranks read the same file after a fence),
// deterministic logic out (sorted ramp + std::map perturbation order + pure
// reconcile). The std::map (not unordered_map) is a correctness requirement,
// not a style choice.
//
// run() in one sentence: forever — re-read what's done from disk, compute the
// next batch of independent work, solve it (all ranks together, one state at a
// time), save each; stop when nothing is left, or when a batch repeats
// unchanged (no progress possible).
// ===========================================================================
class CalcManager {
public:
  using Policy = CalcManagerPolicy;

  CalcManager(ResponsePlan plan, std::string calc_dir, Policy policy = {})
      : plan_(std::move(plan)), calc_dir_(std::move(calc_dir)),
        policy_(policy) {}

  /// Build the node DAG. n_atoms expands the nuclear_all sentinel.
  void build(int n_atoms) { dag_ = build_dag(plan_, n_atoms); }

  const std::vector<CalcNode> &dag() const { return dag_; }
  const std::string &calc_dir() const { return calc_dir_; }

  /// Drive every node to its target protocol against `exec`. One wave per
  /// pass: each pass reloads the aggregate metadata, re-schedules (so actions
  /// reflect current disk), runs the first non-empty wave, fences, then
  /// expands any newly-converged ES bundle into concrete DerivedFD nodes.
  /// Terminates when the schedule is empty (all Skip). A wave that repeats
  /// with no progress is QUARANTINED (its nodes leave the schedule; trace
  /// records a stall_event) and the loop continues with independent work —
  /// stop_reason is 'complete_with_stalled' + diag.stalled_nodes when
  /// anything was quarantined, so a stuck node neither spins the loop nor
  /// starves unrelated perturbation channels (review fix). On exit the plan
  /// is audited against the metadata: planned VBC nodes that never converged
  /// at their top rung (gated out by unconverged prerequisites, or stalled)
  /// are reported via stop_reason 'complete_with_dropped_beta' (or
  /// '..._stalled_and_dropped_beta') + diag.dropped_beta, and recorded under
  /// run_summary/dropped_work in the metadata — a run with silently dropped
  /// beta/raman work is never labelled plain 'complete' (review fix).
  /// Returns the R1c scheduler trace (doc 16 L3) -> Output.diagnostics: the
  /// wave-by-wave reconcile actions (id/thresh/action per item), stop_reason,
  /// and pass count. Built identically on every rank (schedule() is
  /// deterministic), so the value is rank-consistent.
  /// F2 (doc 32 §5): inner per-subworld solve — build the ground state in `sub`
  /// from the archive, solve exactly `items` (its owned FD subset) at `thresh`,
  /// save to disk under the `node_index` metadata shard. Supplied by the
  /// orchestrator (which holds the archive path). Empty ⇒ the single-World path.
  using SubworldSolve =
      std::function<void(madness::World &sub, const std::vector<WorkItem> &items,
                         double thresh, int gid,
                         const std::string &log_prefix)>;

  nlohmann::json run(madness::World &world, ICalcExecutor &exec,
                     SubworldSolve fan_out = {}) {
    const std::vector<double> ramp = global_ramp();
    const std::string meta_path = calc_dir_ + "/response_metadata.json";

    nlohmann::json diag;                       // R1c scheduler trace
    diag["schedule"] = nlohmann::json::array();
    int pass = 0;

    std::string last_sig;
    std::string last_progress;           // fingerprint of the saved solve state
    std::set<std::string> stalled_ids;   // quarantined no-progress nodes
    for (;;) {
      world.gop.fence();
      auto meta = ResponseMetadata::load_or_create(meta_path);
      // Grow the DAG from any ES bundle converged ON DISK (idempotent and
      // restart-safe: reads roots from metadata, not in-memory run state).
      expand_converged_es(world, meta.json());
      // max_iters flows into reconcile so a budget-exhausted rung climbs the
      // ladder (honest-climb) instead of Resume-looping into the no-progress halt.
      auto waves = schedule(dag_, ramp, meta.json(), policy_.max_iters_per_step);
      // Review fix (confirmed HIGH — front-wave starvation): quarantine nodes
      // that stopped making progress instead of halting the WHOLE run. Only
      // the front wave ever executes, so a deterministically-stuck node (e.g.
      // a diverging ES bundle) used to starve every independent item queued
      // in later waves. Stalled nodes are excluded from scheduling, recorded
      // in the trace, and reported honestly in stop_reason at the end.
      for (auto &w : waves) {
        w.erase(std::remove_if(w.begin(), w.end(),
                               [&](const WorkItem &it) {
                                 return stalled_ids.count(it.node->id) > 0;
                               }),
                w.end());
      }
      waves.erase(std::remove_if(waves.begin(), waves.end(),
                                 [](const auto &w) { return w.empty(); }),
                  waves.end());
      if (waves.empty()) {
        // Review fix (confirmed MEDIUM — silent beta drop): an empty schedule
        // does NOT mean every planned node executed. Honest-climb can walk a
        // stubborn FD prerequisite up the whole ladder unconverged; the
        // dependent VBC node is then gated out of every wave (gated nodes are
        // invisible to schedule()), and the beta/raman it feeds is dropped
        // without a trace. Audit the plan against the on-disk metadata: any
        // VBC node without a converged entry at its top rung never delivered.
        nlohmann::json dropped = nlohmann::json::array();
        for (const auto &n : dag_) {
          if (n.kind != CalcKind::VBC || n.protocols.empty()) continue;
          const std::string top = protocol_key_at(n.protocols.back());
          const auto &j = meta.json();
          const bool has_entry = j.contains("vbc_states") &&
                                 j["vbc_states"].contains(n.id) &&
                                 j["vbc_states"][n.id].contains(top);
          if (has_entry && j["vbc_states"][n.id][top].value("converged", false))
            continue;
          // NB: solve_vbc is one-shot and always saves converged=true, so a
          // present VBC entry is caught by the `continue` above — the
          // "built ... not converged" arm is currently unreachable for VBC and
          // kept only as defensive coverage if VBC ever gains partial saves
          // (review LOW).
          const char *reason =
              stalled_ids.count(n.id)
                  ? "stalled (quarantined by the no-progress guard)"
                  : has_entry
                        ? "built at the top rung but not converged"
                        : "prerequisites never converged (gated out of every wave)";
          dropped.push_back({{"id", n.id},
                             {"top_protocol_key", top},
                             {"reason", reason}});
        }
        // Dropped VBC ⇒ the beta/raman rows it feeds cannot be assembled.
        // Record the audit in the metadata (rank 0, through the layer) even
        // when empty, so a later completing run CLEARS a stale drop list.
        if (world.rank() == 0) {
          auto meta_out = ResponseMetadata::load_or_create(meta_path);
          meta_out.set_dropped_work(dropped);
          meta_out.save();
        }

        if (stalled_ids.empty() && dropped.empty()) {
          if (world.rank() == 0)
            madness::print("[CALC] run: nothing left to schedule — done");
          diag["stop_reason"] = "complete";
        } else {
          if (world.rank() == 0 && !stalled_ids.empty()) {
            madness::print("[CALC] run: all remaining work is STALLED —",
                           (int)stalled_ids.size(),
                           "node(s) made no progress and were quarantined:");
            for (const auto &id : stalled_ids) madness::print("    stalled:", id);
          }
          if (world.rank() == 0 && !dropped.empty()) {
            madness::print("[CALC] run: WARNING —", (int)dropped.size(),
                           "planned VBC/beta node(s) were never solved "
                           "(beta/raman rows will be missing):");
            for (const auto &d : dropped)
              madness::print("    dropped:", d.value("id", std::string("?")),
                             "—", d.value("reason", std::string("?")));
          }
          if (dropped.empty())
            diag["stop_reason"] = "complete_with_stalled";
          else if (stalled_ids.empty())
            diag["stop_reason"] = "complete_with_dropped_beta";
          else
            diag["stop_reason"] = "complete_with_stalled_and_dropped_beta";
          if (!stalled_ids.empty()) diag["stalled_nodes"] = stalled_ids;
          if (!dropped.empty())     diag["dropped_beta"]  = dropped;
        }
        break;
      }

      const std::string sig = wave_signature(waves.front());
      // Progress fingerprint: the saved solve state. ONLY the front wave runs
      // each pass, so if the front wave repeats (sig == last_sig) AND the saved
      // state is byte-identical to the previous pass, the last attempt changed
      // NOTHING → genuinely stuck → quarantine. If the state changed (a Resume
      // that advanced iter/residual), it is still converging — give it another
      // pass. Fixes the false-quarantine of a multi-pass ES/VBC Resume (and the
      // legacy max_iters==0 "Resume forever" FD mode), whose reconcile branches
      // don't honest-climb the schedule shape the way the FD/max_iters branch
      // does, so the bare id@protocol signature repeats even mid-convergence.
      const auto &mj = meta.json();
      std::string progress;
      for (const char *k : {"fd_states", "excited_states", "vbc_states"})
        if (mj.contains(k)) progress += mj[k].dump();
      if (sig == last_sig && progress == last_progress) {
        // No progress on this wave: quarantine its nodes and try the rest of
        // the schedule. The run only stops when nothing unquarantined remains.
        if (world.rank() == 0)
          madness::print("[CALC] run: no progress on wave {", sig,
                         "} (saved state unchanged) — quarantining",
                         (int)waves.front().size(),
                         "node(s) and continuing with independent work");
        nlohmann::json srec = {{"pass", pass}, {"wave", sig},
                               {"quarantined", nlohmann::json::array()}};
        for (const auto &it : waves.front()) {
          stalled_ids.insert(it.node->id);
          srec["quarantined"].push_back(it.node->id);
        }
        diag["stall_events"].push_back(std::move(srec));
        last_sig.clear(); last_progress.clear();  // fresh window for the next wave
        continue;
      }
      last_sig = sig;
      last_progress = progress;

      // R1c: record this wave (id/thresh/action per item) + protocol markers.
      // A wave is one protocol level, so wthresh = the front item's threshold.
      const auto &wave = waves.front();
      const double wthresh = wave.front().thresh;
      const int    wk      = default_k_for_thresh(wthresh);
      int pidx = 0;
      for (size_t i = 0; i < ramp.size(); ++i)
        if (std::abs(ramp[i] - wthresh) <= 1e-12 * std::max(ramp[i], wthresh)) {
          pidx = static_cast<int>(i);
          break;
        }
      nlohmann::json wrec = {{"pass", pass}, {"protocol_index", pidx},
                             {"thresh", wthresh}, {"k", wk},
                             {"items", nlohmann::json::array()}};
      for (const auto &it : wave)
        wrec["items"].push_back({{"id", it.node->id},
                                 {"thresh", it.thresh},
                                 {"action", node_action_name(it.action)}});
      diag["schedule"].push_back(std::move(wrec));

      if (world.rank() == 0) {
        madness::print("[CALC] run: wave of", (int)wave.size(),
                       "protocol step(s): {", sig, "}");
        madness::print("PROTOCOL_START  pass=", pass, "  protocol_index=", pidx,
                       "  thresh=", wthresh, "  k=", wk,
                       "  items=", (int)wave.size());
      }
      const bool fan_enabled = fan_out && policy_.fd_subworlds > 0;
      if (!fan_enabled) {
        // Single-World reference path — byte-for-byte the pre-F2 behaviour
        // (every rank solves each wave item together, in wave order).
        for (const auto &item : wave)
          exec.run_protocol(item);
      } else {
        // F2 (doc 32 §5): split the wave into the fan-out-eligible subset (FD /
        // NuclearFD / VBC — all independent per item) and the remainder (ES, which
        // stays single-World). The two sub-phases run SEQUENTIALLY on all universe
        // ranks — never concurrently across communicators (would deadlock). VBC
        // items (F2g) read their converged FD prerequisites from the shared
        // calc_dir; the scheduler's dependency gate guarantees those archives are
        // already on disk from a prior pass, so there is no intra-wave race.
        std::vector<WorkItem> fan_items, rest;
        for (const auto &it : wave) {
          const CalcKind k = it.node->kind;
          ((k == CalcKind::FD || k == CalcKind::NuclearFD || k == CalcKind::VBC)
               ? fan_items
               : rest)
              .push_back(it);
        }

        if (!fan_items.empty()) {
          // F2f: groups_per_node = policy_.fd_subworlds (1 = node-aligned; larger
          // = sub-node / NUMA packing for small systems). G = total subworlds.
          // Guard the pool construction with the same collective-error discipline
          // as the fan_out() call below (review: it was unprotected). This
          // converts a SYMMETRIC/post-collective failure — World ctor, bad_alloc,
          // a Split all ranks fail — into a clean collective abort instead of a
          // hang. (An ASYMMETRIC throw inside make_subworld_pool's own internal
          // collectives — ranks_per_host gather/broadcast, Split, gop.fence — is
          // inherently unrecoverable in MPI: the non-throwing ranks are already
          // blocked inside those collectives and never reach the max() below.)
          NodeSubworldInfo info;
          std::shared_ptr<madness::World> sub;
          std::string pool_err;
          try {
            sub = make_subworld_pool(world, policy_.fd_subworlds, &info);
          } catch (const std::exception &e) { pool_err = e.what(); }
            catch (...) { pool_err = "unknown exception in make_subworld_pool"; }
          int pool_bad = pool_err.empty() ? 0 : 1;
          world.gop.max(pool_bad);
          if (pool_bad) {
            if (!pool_err.empty())
              madness::print("[SUBWORLD-POOL-ERROR] universe rank", world.rank(),
                             ":", pool_err);
            MADNESS_EXCEPTION("subworld pool construction failed (see "
                              "[SUBWORLD-POOL-ERROR]); aborting collectively", 0);
          }
          const int G = info.n_subworlds;
          if (G <= 1) {
            // One subworld total (single node, P=1) ⇒ no real partition (subworld
            // == universe). Literal G=0 path, no redundant GS reload (doc 32 §7).
            if (world.rank() == 0)
              madness::print("SUBWORLD_FANOUT skipped: n_subworlds=", G,
                             " (nodes=", info.n_nodes, " groups_per_node=",
                             info.groups_per_node, ") — universe path");
            sub.reset();
            for (const auto &it : fan_items) exec.run_protocol(it);
          } else {
            std::vector<WorkItem> mine;   // round-robin by global gid (deterministic)
            for (size_t i = 0; i < fan_items.size(); ++i)
              if (static_cast<int>(i) % G == info.gid)
                mine.push_back(fan_items[i]);
            // F2d: per-subworld log tag "[g{gid}/{G} {host}] " (doc 32 §5.6).
            const std::string lp = "[g" + std::to_string(info.gid) + "/" +
                                   std::to_string(G) + " " + info.hostname + "] ";
            if (world.rank() == 0)
              madness::print("SUBWORLD_FANOUT  pass=", pass, "  n_subworlds=", G,
                             "  groups_per_node=", info.groups_per_node,
                             "  fan_items=", (int)fan_items.size());
            // S1 pmap discipline: point the default pmap at the subworld so
            // everything built inside is subworld-local; restore BEFORE reset.
            // Exception safety: a subworld-local throw must NOT skip the pmap
            // restore or the universe fence below — the failing ranks would
            // reach the driver's top-level catch and finalize while every
            // other rank blocks in the fence forever. Catch, restore, fence,
            // then agree collectively and abort together. (A throw that leaves
            // SIBLING ranks of the same subworld inside a subworld collective
            // is inherently unrecoverable in MPI — this handles the
            // synchronized-failure cases: archive open, shard save, guards.)
            std::string fan_err;
            madness::FunctionDefaults<3>::set_default_pmap(*sub);
            try {
              fan_out(*sub, mine, wthresh, info.gid, lp);
              sub->gop.fence();
            } catch (const std::exception &e) {
              fan_err = e.what();
            } catch (...) {
              fan_err = "unknown exception in subworld fan-out";
            }
            madness::FunctionDefaults<3>::set_default_pmap(world);
            sub.reset();
            if (!fan_err.empty())
              madness::print("[FANOUT-ERROR] universe rank", world.rank(),
                             "(gid", info.gid, "):", fan_err);
            int fan_bad = fan_err.empty() ? 0 : 1;
            world.gop.max(fan_bad);
            if (fan_bad) {
              world.gop.fence();
              MADNESS_EXCEPTION(
                  "subworld fan-out failed on at least one subworld "
                  "(see [FANOUT-ERROR]); aborting collectively", 0);
            }
          }
          world.gop.fence();
          // F2b/F2g: collapse the per-group FD+VBC metadata shards into the
          // canonical file (universe rank 0; through the metadata layer). Only
          // when we actually sharded (G > 1). The merge itself can throw
          // (corrupt shard, ENOSPC) — that must be collective too, not a
          // rank-0-only unwind that strands everyone else at the fence.
          int merge_bad = 0;
          if (world.rank() == 0 && G > 1) {
            try {
              ResponseMetadata::merge_state_shards(calc_dir_, G);
            } catch (const std::exception &e) {
              merge_bad = 1;
              madness::print("[SHARD-MERGE-ERROR]", e.what());
            }
          }
          world.gop.max(merge_bad);
          if (merge_bad)
            MADNESS_EXCEPTION("per-wave shard merge failed on rank 0 "
                              "(see [SHARD-MERGE-ERROR])", 0);
          world.gop.fence();
        } else if (policy_.fd_subworlds > 0 && !rest.empty() &&
                   world.rank() == 0) {
          // Review MED: this wave is ES-only (ES is deliberately excluded from
          // fan-out — it stays single-World), so the requested `subworlds` had
          // no effect on it. Say so, or the knob looks silently ignored.
          madness::print("SUBWORLD_FANOUT skipped: wave is ES-only (",
                         (int)rest.size(),
                         "item(s)); ES runs single-World, `subworlds` applies "
                         "only to FD/VBC waves");
        }
        for (const auto &it : rest) exec.run_protocol(it);   // ES: universe
      }
      world.gop.fence();
      if (world.rank() == 0)
        madness::print("PROTOCOL_DONE  pass=", pass, "  protocol_index=", pidx,
                       "  thresh=", wthresh, "  items=", (int)wave.size());
      ++pass;
    }
    diag["passes"] = pass;
    return diag;
  }

private:
  /// Coarse->fine union of every node's protocol ladder.
  std::vector<double> global_ramp() const {
    std::vector<double> all;
    for (const auto &n : dag_)
      for (double t : n.protocols) all.push_back(t);
    return detail_calc::sorted_ramp(std::move(all));
  }

  /// Stable, rank-deterministic signature of a wave (sorted id@protocol list).
  static std::string wave_signature(const std::vector<WorkItem> &wave) {
    std::vector<std::string> parts;
    parts.reserve(wave.size());
    for (const auto &it : wave)
      // protocol_key, NOT freq_key: freq_key's fixed-decimal format collapses
      // 1e-6 and 1e-8 to the same string ("f0.00000"), so two consecutive
      // LADDER rungs got identical signatures and the no-progress guard
      // falsely halted between them (exposed by honest-climb, which makes
      // consecutive front waves be DIFFERENT rungs of the same node).
      parts.push_back(it.node->id + "@" +
                      protocol_key(it.thresh, default_k_for_thresh(it.thresh)));
    std::sort(parts.begin(), parts.end());
    std::string s;
    for (const auto &p : parts) { s += p; s += ';'; }
    return s;
  }

  /// Turn each symbolic DerivedFD "*" whose ES bundle has CONVERGED ON DISK into
  /// one promoted FD node per root (freq = root energy at the ES bundle's top
  /// protocol). Reads roots straight from response_metadata.json so it is
  /// restart-safe (works on a resumed run where the ES converged in a previous
  /// process) and idempotent (new nodes are deduped against existing dag_ ids;
  /// re-running is a no-op). Deterministic across ranks: same metadata + same
  /// dag_ in, same nodes appended in the same order.
  void expand_converged_es(madness::World &world, const nlohmann::json &meta) {
    if (!meta.contains("excited_states")) return;
    std::set<std::string> have;
    for (const auto &n : dag_) have.insert(n.id);

    std::vector<CalcNode> additions;
    for (const auto &sym : dag_) {
      if (!sym.is_symbolic()) continue;
      if (sym.prerequisites.empty()) continue;          // no ES dependency
      const std::string &es_id = sym.prerequisites.front();
      // Resolve the ES node to its top protocol (derived FDs start once the ES
      // bundle has reached its finest protocol).
      const CalcNode *es = nullptr;
      for (const auto &n : dag_) if (n.id == es_id) { es = &n; break; }
      if (!es || es->protocols.empty()) continue;
      const std::string es_key = protocol_key_at(es->protocols.back());
      if (!meta["excited_states"].contains(es_key)) continue;
      const auto &b = meta["excited_states"][es_key];
      if (!b.value("converged", false)) continue;       // ES not converged yet
      if (!b.contains("roots") || !b["roots"].is_array()) continue;

      for (size_t i = 0; i < b["roots"].size(); ++i) {
        const double w = b["roots"][i].value(
            "omega", std::numeric_limits<double>::quiet_NaN());
        if (!(w == w)) continue;                         // skip NaN omega
        const double fd_freq = sym.es_freq_factor * w;   // e.g. ωₙ/2 (off-pole)
        CalcNode c;
        // A derived FD at ω = factor·(root energy) is an ordinary FD point ->
        // kind FD, so the FD executor solves it and skip/restart/seed logic is
        // uniform. factor=0.5 puts it at the two-photon resonance ωₙ/2, off the
        // linear-response pole at ω = ωₙ.
        c.kind       = CalcKind::FD;
        c.pert       = sym.pert;
        c.freq       = fd_freq;
        c.protocols  = sym.protocols;
        c.es_root_id = make_es_root_label(static_cast<int>(i));  // provenance
        c.id         = fd_node_id(c.pert, fd_freq);  // id keys on the node's TRUE freq (ωₙ/2), not ωₙ
        if (have.count(c.id)) continue;                  // dedup -> idempotent
        have.insert(c.id);
        additions.push_back(std::move(c));
      }
    }
    if (!additions.empty()) {
      if (world.rank() == 0)
        madness::print("[CALC] expanded", (int)additions.size(),
                       "derived FD node(s) at ω = es_freq_factor·ωₙ from converged ES roots (metadata)");
      for (auto &c : additions) dag_.push_back(std::move(c));
    }
  }

  static std::string make_es_root_label(int i) {
    char buf[24];
    std::snprintf(buf, sizeof buf, "es_root_%04d", i);
    return buf;
  }

  ResponsePlan          plan_;
  std::string           calc_dir_;
  Policy                policy_;
  std::vector<CalcNode> dag_;
};

// ---------------------------------------------------------------------------
// beta property assembly (Tier-A). A post-run step: contract the converged A
// states against each VBC source (+ zeta terms) into the beta tensor, recorded
// under properties/beta. Runs at one protocol (typically the top of the ramp).
// ---------------------------------------------------------------------------

inline char beta_axis_name(int a) { return a == 0 ? 'x' : a == 1 ? 'y' : 'z'; }

inline void assemble_beta(ExecutorContext &ctx, const ResponsePlan &plan,
                          double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  if (plan.vbc.empty()) return;
  // Same ClosedShell-only reader constraint as assemble_alpha — refuse loudly.
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[ASSEMBLE] open-shell beta/raman assembly is not implemented — "
            "SKIPPED (ClosedShell-only readers; see REVIEW_FINDINGS).");
    return;
  }

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);
  const std::string key = protocol_key();

  if (world.rank() == 0)
    madness::print("\n=== beta assembly  protocol_key=", key, " ===");

  // Accumulate property rows and write them to response_metadata.json in ONE
  // load/add-all/save after the loops (rank 0), instead of reloading + rewriting
  // the whole file per (vbc, axis). `key` is constant across all rows.
  std::vector<std::pair<std::string, nlohmann::json>> rows;  // (property key, row)

  // Rank-0 metadata snapshot for per-row accuracy: each beta/raman row records
  // the honest verdict of the FD states that built it (same discipline as
  // alpha's row_accuracy; madqc review R4). Uses the SAME source-selection
  // helper the loader uses, so reported accuracy == the state actually loaded.
  nlohmann::json acc_meta;
  if (world.rank() == 0)
    acc_meta = ResponseMetadata::load_or_create(
                   ctx.calc_dir + "/response_metadata.json").json();
  const double acc_thr = madness::FunctionDefaults<3>::get_thresh();
  const int    acc_k   = madness::FunctionDefaults<3>::get_k();
  auto fd_acc = [&](const Perturbation &p, double f) -> nlohmann::json {
    nlohmann::json a;   // rank-0 only: rows are assembled on rank 0
    const std::string chan = p.description();
    const std::string fk   = ResponseMetadata::freq_key(f);
    const std::string sk =
        best_usable_fd_source_key(acc_meta, chan, fk, acc_thr, acc_k, key);
    a["input"] = chan;
    a["freq"]  = f;
    a["source_protocol_key"] = sk;
    if (!sk.empty() && acc_meta["fd_states"][chan][sk].contains(fk)) {
      const auto &e = acc_meta["fd_states"][chan][sk][fk];
      // Same semantics as alpha's row_accuracy (review fix): `converged` is the
      // honest STRICT verdict — an accepted-at-maxiter source (metadata
      // converged forced true to unblock the VBC gate) reports converged=false
      // + accepted=true here, never a silently over-stated accuracy claim.
      const bool acc = e.value("accepted", false);
      a["converged"]    = e.value("converged", false) && !acc;
      a["accepted"]     = acc;
      a["bsh_residual"] = e.value("bsh_residual", 0.0);
    }
    return a;
  };

  for (const auto &vr : plan.vbc) {
    const double ws = vr.freq_b + vr.freq_c;   // omega_sigma = omega_B + omega_C
    const std::string vbc_id =
        vbc_node_id(vr.pert_b, vr.pert_c, vr.freq_b, vr.freq_c);

    const bool is_raman =
        (vr.pert_c.kind == Perturbation::Kind::NuclearDisplacement);
    const char *pkey = is_raman ? "raman" : "beta";
    const char *ptag = is_raman ? "[RAMAN]" : "[BETA]";
    auto vbc = load_vbc<ClosedShell>(world, ctx.calc_dir, vbc_id);
    auto B = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, vr.pert_b, vr.freq_b);
    auto C = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, vr.pert_c, vr.freq_c);
    if (!vbc || !B || !C) {
      if (world.rank() == 0)
        madness::print("[BETA] skip", vbc_id, "— missing VBC or FD input");
      continue;
    }

    for (int a = 0; a < 3; ++a) {
      const Perturbation pA = Perturbation::dipole(a);
      auto xA = detail_exec::load_fd_as_xy<ClosedShell>(world, ctx.calc_dir, pA, ws);
      if (!xA) {
        if (world.rank() == 0)
          madness::print("[BETA] skip A=", beta_axis_name(a), "— missing FD",
                         pA.description(), "@", ws);
        continue;
      }
      auto VA_op = dipole_operator(world, a);
      const double b = beta::beta_abc<ClosedShell>(world, g0, *xA, *vbc, *B, *C, VA_op);

      if (world.rank() == 0) {
        madness::print(ptag, " A=", beta_axis_name(a),
                       "  B=", vr.pert_b.description(),
                       "  C=", vr.pert_c.description(),
                       "  fB=", vr.freq_b, "  fC=", vr.freq_c,
                       "  value=", b);
        nlohmann::json row{{"A", std::string(1, beta_axis_name(a))},
                           {"B", vr.pert_b.description()},
                           {"C", vr.pert_c.description()},
                           {"freq_b", vr.freq_b},
                           {"freq_c", vr.freq_c},
                           {"beta", b}};
        // Honest per-row accuracy: the three FD inputs (A@ws, B@fB, C@fC) and
        // the VBC state's own verdict at this protocol.
        row["row_accuracy"] = nlohmann::json::array(
            {fd_acc(pA, ws), fd_acc(vr.pert_b, vr.freq_b),
             fd_acc(vr.pert_c, vr.freq_c)});
        if (acc_meta.contains("vbc_states") &&
            acc_meta["vbc_states"].contains(vbc_id) &&
            acc_meta["vbc_states"][vbc_id].contains(key)) {
          const auto &v = acc_meta["vbc_states"][vbc_id][key];
          row["vbc_accuracy"] = {{"converged", v.value("converged", false)},
                                 {"diverged", v.value("diverged", false)}};
        }
        // Aggregate row verdict (review MED): mirror alpha's row-level
        // converged/max_bsh_residual so the .out summary and the website viewer
        // get ONE honest flag per beta/raman row instead of only the per-leg
        // sub-structures. The row is converged iff all three FD legs met the
        // strict target (VBC itself is non-iterative, so its always-true
        // verdict adds no signal — the FD legs carry the real accuracy).
        {
          bool all_conv = true; double max_res = 0.0;
          for (const auto &leg : row["row_accuracy"]) {
            all_conv = all_conv && leg.value("converged", false);
            max_res  = std::max(max_res, leg.value("bsh_residual", 0.0));
          }
          row["converged"]        = all_conv;
          row["max_bsh_residual"] = max_res;
        }
        rows.emplace_back(pkey, std::move(row));
      }
      world.gop.fence();
    }
  }

  if (world.rank() == 0 && !rows.empty()) {
    auto meta = ResponseMetadata::load_or_create(
        ctx.calc_dir + "/response_metadata.json");
    for (const auto &[pk, row] : rows) meta.add_property(pk, key, row);
    meta.save();
  }
}

// ---------------------------------------------------------------------------
// assemble_tpa — Tier-B two-photon-absorption assembly (post-solve, off the
// critical path; mirrors assemble_beta). Needs a converged ES bundle + the
// derived dipole FD at omega_f/2 (both produced by the resonant/ES plan). For
// each root f: load X_f + the three mu_a responses at omega_f/2, contract via
// tpa::tpa_moment (the beta-residue with a homogeneous C-channel), record the
// S tensor + delta^|| under properties/tpa. ClosedShell/TDA only.
// The residue form is a CANDIDATE validated against refs/dalton_tpa.json.
// ---------------------------------------------------------------------------
/// Load the ES bundle of EsType and return per-root (omega, X_f as XY). Full
/// (RPA) carries the real de-excitation Y; TDA sets Y=0. Templated so the
/// EsType Storage (X-only vs X,Y) is handled at compile time.
template <typename EsType>
inline bool load_tpa_es_xy(madness::World &world, const std::string &calc_dir,
                           std::vector<double> &omega_out,
                           std::vector<ResponseStateXY<ClosedShell>> &Xf_out) {
  auto es = try_load_es_bundle<EsType, ClosedShell>(world, calc_dir);
  if (!es) return false;
  auto &state = es->state;
  const long nr = state.omega.dim(0);
  for (long f = 0; f < nr; ++f) {
    omega_out.push_back(state.omega(f));
    ResponseStateXY<ClosedShell> Xf;
    Xf.x_alpha = state.roots[f].x_alpha;
    if constexpr (std::is_same_v<EsType, Full>) {
      Xf.y_alpha = state.roots[f].y_alpha;               // RPA de-excitation
    } else {
      Xf.y_alpha = madness::copy(world, Xf.x_alpha);      // TDA: y = 0
      madness::scale(world, Xf.y_alpha, 0.0);
    }
    Xf_out.push_back(std::move(Xf));
  }
  return true;
}

inline void assemble_tpa(ExecutorContext &ctx, const ResponsePlan &plan,
                         double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;
  if (plan.es.empty()) return;
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[TPA] open-shell TPA not implemented — SKIPPED (ClosedShell only).");
    return;
  }

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const double c_xc = gs.hf_exchange_coefficient();
  const double lo   = gs.params().lo();
  auto g0 = build_response_ground_state_closed_shell(world, gs, c_xc, lo);
  const std::string key = protocol_key();

  // Dispatch on the ES method recorded in the bundle metadata (tda vs full/RPA).
  // Dalton's TPA reference is full RPA; --es-full gives the matching v3 bundle.
  int es_full = 0;
  if (world.rank() == 0) {
    auto meta = ResponseMetadata::load_or_create(
        ctx.calc_dir + "/response_metadata.json");
    const auto &j = meta.json();
    if (j.contains("excited_states") && j["excited_states"].contains(key))
      es_full = (j["excited_states"][key].value("type", "tda") == "full") ? 1 : 0;
  }
  world.gop.broadcast(es_full, 0);

  std::vector<double> omegas;
  std::vector<ResponseStateXY<ClosedShell>> Xfs;
  const double t_es0 = madness::wall_time();
  const bool loaded =
      es_full ? load_tpa_es_xy<Full>(world, ctx.calc_dir, omegas, Xfs)
              : load_tpa_es_xy<TDA>(world, ctx.calc_dir, omegas, Xfs);
  if (world.rank() == 0) {
    printf("[TPA timing] ES bundle load: %.1f s\n", madness::wall_time() - t_es0);
    fflush(stdout);
  }
  if (!loaded) {
    if (world.rank() == 0)
      print("[TPA] no ES bundle under", ctx.calc_dir, "— SKIPPED");
    return;
  }
  const long nroots = static_cast<long>(omegas.size());
  if (world.rank() == 0)
    print("[TPA] ES method:", es_full ? "Full (RPA)" : "TDA");

  std::array<real_function_3d, 3> mu_op{dipole_operator(world, 0),
                                        dipole_operator(world, 1),
                                        dipole_operator(world, 2)};

  if (world.rank() == 0)
    print("\n=== TPA assembly  protocol_key=", key, "  n_roots=", nroots, " ===");

  std::vector<nlohmann::json> rows;
  // rank-0 table data (S tensor + observables per root) for the Dalton-style print.
  std::vector<int>                     tbl_f;
  std::vector<double>                  tbl_w;
  std::vector<madness::Tensor<double>> tbl_S;
  std::vector<tpa::Observables>        tbl_o;
  for (long f = 0; f < nroots; ++f) {
    if (!ctx.tpa_roots.empty() &&
        std::find(ctx.tpa_roots.begin(), ctx.tpa_roots.end(),
                  static_cast<int>(f)) == ctx.tpa_roots.end())
      continue;   // --tpa-roots filter (parallel per-root verification)
    const double wf = omegas[f];
    const double wf_half = 0.5 * wf;
    const ResponseStateXY<ClosedShell> &Xf = Xfs[f];   // real Y for RPA, 0 for TDA

    std::array<ResponseStateXY<ClosedShell>, 3> mu_resp;
    bool ok = true;
    const double t_fd0 = madness::wall_time();
    for (int a = 0; a < 3; ++a) {
      auto r = detail_exec::load_fd_as_xy<ClosedShell>(
          world, ctx.calc_dir, Perturbation::dipole(a), wf_half);
      if (!r) { ok = false; break; }
      mu_resp[a] = std::move(*r);
    }
    if (world.rank() == 0)
      printf("[TPA timing] root %ld/%ld: 3 FD loads @ w=%.5f: %.1f s\n",
             f + 1, nroots, wf_half, madness::wall_time() - t_fd0);
             fflush(stdout);
    if (!ok) {
      if (world.rank() == 0)
        print("[TPA] root", f, " omega_f=", wf,
              " — missing dipole FD @", wf_half, " — skip");
      continue;
    }

    // ---- cross-code diagnostics (PrintLevel >= Normal or --tpa-diag-only):
    // vector norms, transition dipoles/oscillators, and alpha(omega_f/2) from
    // the SAME loaded states the contraction uses — line-for-line comparable
    // with the patched-DALTON QRSMONORM + QRLRVE output (TPA_SCOPING §5n/§5o).
    if (static_cast<int>(ctx.print_level) >= 1 || ctx.tpa_diag_only) {
      const double xx = inner(Xf.x_alpha, Xf.x_alpha);
      const double yy = inner(Xf.y_alpha, Xf.y_alpha);
      if (world.rank() == 0)
        printf("  [TPA diag] root %ld  omega=%.8f (%.3f eV)\n"
               "    Xf norms: |x|=%.10f  |y|=%.10f  x2-y2=%.10f "
               "(DALTON C(EXCI): Z2-Y2=0.5 => ours/DALTON = sqrt2)\n",
               f, wf, wf * 27.211386245988, std::sqrt(xx), std::sqrt(yy),
               xx - yy);
      std::array<vecfuncT, 3> mu_amo;
      for (int a = 0; a < 3; ++a) mu_amo[a] = mul(world, mu_op[a], g0.amo, true);
      // transition dipole + oscillator (Parker-normalized: M = -sqrt2 * t)
      for (int a = 0; a < 3; ++a) {
        const double t = inner(mu_amo[a], Xf.x_alpha) + inner(mu_amo[a], Xf.y_alpha);
        if (world.rank() == 0 && std::abs(t) > 1e-6)
          printf("    transition dipole %c: t=%+.8f  M=-sqrt2*t=%+.8f  "
                 "osc=(2/3)w*2t^2=%.6f\n",
                 "xyz"[a], t, -std::sqrt(2.0) * t,
                 (2.0 / 3.0) * wf * 2.0 * t * t);
      }
      // alpha(omega_f/2) full 3x3 from the loaded FD states (validated formula)
      for (int a = 0; a < 3; ++a) {
        const double nx = inner(mu_resp[a].x_alpha, mu_resp[a].x_alpha);
        const double ny = inner(mu_resp[a].y_alpha, mu_resp[a].y_alpha);
        if (world.rank() == 0)
          printf("    N^%c(w/2) norms: |x|=%.10f  |y|=%.10f\n", "xyz"[a],
                 std::sqrt(nx), std::sqrt(ny));
      }
      if (world.rank() == 0) printf("    alpha(w=%.6f):", wf_half);
      for (int a = 0; a < 3; ++a) {
        const double al = -2.0 * (inner(mu_resp[a].x_alpha, mu_amo[a]) +
                                  inner(mu_resp[a].y_alpha, mu_amo[a]));
        if (world.rank() == 0) printf("  %c%c=%.6f", "xyz"[a], "xyz"[a], al);
      }
      if (world.rank() == 0) printf("   (DALTON: QRLRVE <<A;A>> at same w)\n");
    }
    if (ctx.tpa_diag_only) continue;   // diagnostics only — no contraction,
                                       // no metadata writes (read-only pass)

    madness::Tensor<double> S;
    std::vector<std::string> src_files;
    if (ctx.tpa_residue) {
      // PRODUCTION composition (TPA_SCOPING §5m + §5q, verified 6/6 vs DALTON
      // d-aug-QZ 2026-07-22 via verify_e3_k6.py):
      //   S = C_N * ( S_1e + S_E3corr ),   C_N = sqrt(2)
      // S_1e: mu-operator terms of V^{bc} only (tpa_moment_residue_1e; cheap).
      // S_E3corr: corrected two-electron composition (tpa_e3_residue; DALTON
      // units WITH the sqrt2 — divided out here so C_N is applied exactly
      // once and ctx.tpa_prefactor stays a pure A/B knob on top).
      // BOTH orderings are computed off-diagonal: the composition is
      // analytically b<->c symmetric, so the asymmetry is a built-in
      // correctness assertion (DALTON's E3 shows the same invariance).
      if (world.rank() == 0) {
        print("[TPA] contraction: single-residue, S = sqrt2*(1e + E3corr),"
              " 2e composition =",
              (ctx.tpa_cgrouped ? "c-grouped (tpa_e3_residue)"
                                : "Parker (P,Q) spec [production]"),
              " prefactor =", ctx.tpa_prefactor);
        fflush(stdout);
      }
      const auto S_1e =
          tpa::tpa_moment_residue_1e(world, g0, Xf, mu_resp, mu_op, 1.0);
      madness::Tensor<double> S_e3c(3L, 3L);
      double max_asym = 0.0;
      for (int b3 = 0; b3 < 3; ++b3) {
        for (int c3 = b3; c3 < 3; ++c3) {
          if (world.rank() == 0) {
            printf("  [TPA e3corr] root %ld pair %c%c%s (%s) ...\n", f,
                   "xyz"[b3], "xyz"[c3],
                   (b3 != c3 ? " (both orderings)" : ""),
                   (ctx.tpa_cgrouped ? "c-grouped" : "Parker spec"));
            fflush(stdout);
          }
          double e = 0.0;
          if (!ctx.tpa_cgrouped) {
            // PRODUCTION (2026-09-05): the Parker-style state-free (P,Q)
            // source, symmetrized builder (family-D collapse). Gate 1 of
            // test_tpa_pq_vs_vbc pins <x^f|P>+<y^f|Q> (one ordering, 2e)
            // == tpa_e3_residue/sqrt2, and the sym builder returns the
            // ordering SUM, so 0.5*<f|P_sym,Q_sym> lands in the same units
            // as the averaged e/sqrt2 below — the outer sqrt2*(1e+E3corr)
            // stays untouched.
            auto pq = source_spec::assemble_source(
                world, g0,
                tpa::tpa_pq_spec_sym(world, g0,
                                     mu_resp[static_cast<size_t>(b3)],
                                     mu_resp[static_cast<size_t>(c3)]));
            const double sp =
                inner(world, Xf.x_alpha, pq[0]).sum() +
                inner(world, Xf.y_alpha, pq[1]).sum();
            e = 0.5 * sp;
            S_e3c(b3, c3) = S_e3c(c3, b3) = e;
          } else {
            const double ebc = tpa::tpa_e3_residue(
                world, g0, mu_resp[static_cast<size_t>(b3)],
                mu_resp[static_cast<size_t>(c3)], Xf);
            e = ebc;
            if (b3 != c3) {
              const double ecb = tpa::tpa_e3_residue(
                  world, g0, mu_resp[static_cast<size_t>(c3)],
                  mu_resp[static_cast<size_t>(b3)], Xf);
              max_asym = std::max(max_asym, std::abs(ebc - ecb));
              e = 0.5 * (ebc + ecb);
            }
            S_e3c(b3, c3) = S_e3c(c3, b3) = e / std::sqrt(2.0);
          }
        }
      }
      S = madness::copy(S_1e);
      S += S_e3c;
      S.scale(std::sqrt(2.0) * ctx.tpa_prefactor);
      if (world.rank() == 0)
        printf("  [TPA] root %ld assembled: S = sqrt2*(1e + E3corr)   "
               "b<->c max asym = %.2e\n", f, max_asym);
      if (ctx.tpa_decompose) {
        // A/B diagnostics vs the LEGACY vbc-based contraction (old E3
        // composition, kept for comparison): table stays in prefactor-1
        // "ours" units so verify_e3_k6.py's x sqrt(2) recovers DALTON.
        auto S_full = tpa::tpa_moment_residue(world, g0, Xf, mu_resp, mu_op,
                                              1.0);
        real_function_3d zop = madness::copy(mu_op[0]); zop.scale(0.0);
        std::array<real_function_3d, 3> zops{zop, madness::copy(zop),
                                             madness::copy(zop)};
        auto S_e3 = tpa::tpa_moment_residue(world, g0, Xf, mu_resp, zops,
                                            1.0);
        if (world.rank() == 0) {
          double max_1e_dev = 0.0;   // direct 1e vs (legacy full - legacy E3)
          for (int b3 = 0; b3 < 3; ++b3)
            for (int c3 = 0; c3 < 3; ++c3)
              max_1e_dev = std::max(
                  max_1e_dev, std::abs(S_1e(b3, c3) - (S_full(b3, c3) -
                                                       S_e3(b3, c3))));
          printf("  [TPA decompose] root %ld (E3old, E3corr, 1e, total)   "
                 "b<->c max asym = %.2e   1e direct-vs-derived max dev = %.2e\n",
                 f, max_asym, max_1e_dev);
          const char *nm[6] = {"xx", "yy", "zz", "xy", "xz", "yz"};
          const int ii[6] = {0, 1, 2, 0, 0, 1}, jj[6] = {0, 1, 2, 1, 2, 2};
          for (int e = 0; e < 6; ++e) {
            const double tot = S_full(ii[e], jj[e]), e3 = S_e3(ii[e], jj[e]);
            printf("    %s: E3=%+.6f  E3corr=%+.6f  1e=%+.6f  tot=%+.6f  "
                   "totcorr=%+.6f\n",
                   nm[e], e3, S_e3c(ii[e], jj[e]), tot - e3, tot,
                   tot - e3 + S_e3c(ii[e], jj[e]));
                   fflush(stdout);
          }
        }
      }
    } else {
      // Legacy candidate: build the 3 axis sources (the TPA analogue of beta's
      // VBC quadratic source), SAVE them to disk (mirrors beta/Raman
      // vbc_states — reusable + inspectable/isosurface-able), then contract.
      auto vbc_b = tpa::tpa_sources(world, g0, Xf, mu_resp, mu_op);
      for (int b = 0; b < 3; ++b) {
        const std::string src = ctx.calc_dir + "/tpa_src__root" +
            std::to_string(f) + "__" + std::string(1, "xyz"[b]) + "__" + key;
        vbc_b[b].save(world, src);
        src_files.push_back(std::filesystem::path(src).filename().string());
      }
      S = tpa::tpa_moment(world, g0, Xf, mu_resp, mu_op, vbc_b);
    }
    const auto obs = tpa::observables(S, wf);

    if (world.rank() == 0) {
      nlohmann::json Sj = nlohmann::json::array();
      for (int a = 0; a < 3; ++a) {
        nlohmann::json ra = nlohmann::json::array();
        for (int b = 0; b < 3; ++b) ra.push_back(S(a, b));
        Sj.push_back(ra);
      }
      rows.push_back({{"es_root_id", static_cast<int>(f)},
                      {"omega", wf},
                      {"omega_ev", wf * 27.211386245988},
                      {"source_files", src_files},
                      {"writer_nproc", static_cast<int>(world.size())},
                      {"S", Sj},
                      {"Df", obs.Df}, {"Dg", obs.Dg},
                      {"D_linear", obs.D_linear}, {"D_circular", obs.D_circular},
                      {"R", obs.R},
                      {"sigma_linear_gm", obs.sigma_linear_gm},
                      {"sigma_circular_gm", obs.sigma_circular_gm},
                      {"delta_parallel", obs.D_linear}});
      tbl_f.push_back(static_cast<int>(f)); tbl_w.push_back(wf);
      tbl_S.push_back(S); tbl_o.push_back(obs);
    }
    world.gop.fence();
  }

  // Dalton-style output (matches rspvec.F QRSMO): tensor S table + the
  // Monson-McClain Df/Dg/D(lin,circ)/sigma/R summary. No point-group symmetry
  // in MRA, so one manifold (no Sym column). Diffable against Dalton .out.
  if (world.rank() == 0 && !tbl_f.empty()) {
    const double H2EV = 27.211386245988;
    printf("\n                  +--------------------------------+\n");
    printf("                  | Two-photon transition tensor S |\n");
    printf("                  +--------------------------------+\n");
    printf("     -----------------------------------------------------------------\n");
    printf("      No  Energy(eV)     Sxx     Syy     Szz     Sxy     Sxz     Syz\n");
    printf("     -----------------------------------------------------------------\n");
    for (size_t r = 0; r < tbl_f.size(); ++r)
      printf("     %3d   %8.2f  %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, tbl_S[r](0, 0), tbl_S[r](1, 1),
             tbl_S[r](2, 2), tbl_S[r](0, 1), tbl_S[r](0, 2), tbl_S[r](1, 2));
    printf("     -----------------------------------------------------------------\n");
    printf("\n     D = 2*Df+4*Dg (Linear);  D = -2*Df+6*Dg (Circular)\n");
    printf("     Df = sum_ij S_ii*S_jj /30;   Dg = sum_ij S_ij^2 /30\n");
    printf("     sigma = D*(E/2)^2*AU_TO_GM (GM, 0.1 eV FWHM);  R = (-Df+3Dg)/(Df+2Dg)\n");
    printf("\n                   +-----------------------------------+\n");
    printf("                   |   Two-photon absorption summary   |\n");
    printf("                   +-----------------------------------+\n");
    printf("      No Energy(eV) Polarization      Df         Dg          D        sigma      R\n");
    printf("     ---------------------------------------------------------------------------------\n");
    for (size_t r = 0; r < tbl_f.size(); ++r) {
      const auto &o = tbl_o[r];
      printf("     %3d %8.2f  Linear    %10.3E %10.3E %10.3E %10.3E %6.2f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, o.Df, o.Dg, o.D_linear, o.sigma_linear_gm, o.R);
      printf("     %3d %8.2f  Circular  %10.3E %10.3E %10.3E %10.3E %6.2f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, o.Df, o.Dg, o.D_circular, o.sigma_circular_gm, o.R);
    }
    printf("     ---------------------------------------------------------------------------------\n");
    fflush(stdout);
  }

  // Root-filtered runs are read-only (concurrent per-root processes must not
  // race on the metadata file); full runs record properties/tpa as before.
  if (world.rank() == 0 && !rows.empty() && ctx.tpa_roots.empty()) {
    auto meta = ResponseMetadata::load_or_create(
        ctx.calc_dir + "/response_metadata.json");
    for (auto &row : rows) meta.add_property("tpa", key, row);
    meta.save();
  }
}

// ---------------------------------------------------------------------------
// Polarizability assembly (Tier-A, closed-shell). After the FD states converge,
// contract each converged response against the dipole perturbation vectors:
//
//   alpha_ij(omega) = -2 * ( <x_i | v_j> + <y_i | v_j> )
//
// where x_i,y_i are the converged FD response (direction i, frequency omega) and
// v_j = Q(mu_j * phi0) is the dipole perturbation source (direction j). For the
// static case load_fd_as_xy supplies y = x, so this reduces to -4<x|v> -- the
// same factors as molresponse_v2 PropertyManager::compute_alpha (-2 dynamic /
// -4 static). Printed as [ALPHA ...] and recorded under properties/alpha.
// ---------------------------------------------------------------------------
inline void assemble_alpha(ExecutorContext &ctx, const ResponsePlan &plan,
                           double thresh) {
  using namespace madness;
  World &world = ctx.world;
  GroundState &gs = ctx.gs;

  // Review finding (confirmed HIGH): the assembly readers below are hardcoded
  // to the ClosedShell archive layout; an open-shell state deserializes
  // MISALIGNED (crash or garbage alpha). The open-shell SOLVES are fine and
  // saved — refuse the assembly loudly instead of fabricating numbers.
  if (!gs.is_spin_restricted()) {
    if (world.rank() == 0)
      print("[ASSEMBLE] open-shell property assembly is not implemented yet — "
            "response states are solved and saved, but alpha assembly is "
            "SKIPPED (ClosedShell-only readers; see REVIEW_FINDINGS).");
    return;
  }

  // Dipole directions + frequencies present in the plan.
  std::set<int>    axes;
  std::set<double> freqs;
  for (const auto &r : plan.fd)
    if (r.pert.kind == Perturbation::Kind::Dipole) {
      axes.insert(r.pert.axis);
      freqs.insert(r.freq);
    }
  if (axes.empty() || freqs.empty()) return;
  const std::vector<int>    ax(axes.begin(),  axes.end());
  const std::vector<double> fr(freqs.begin(), freqs.end());

  set_response_protocol(world, ctx.L, thresh);
  const double t0 = FunctionDefaults<3>::get_thresh();
  {
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
    gs.prepare(world, 0.001 * t0, coulop, ctx.fock_json);
  }
  const std::string key = protocol_key();

  // Dipole perturbation source per direction (built once at this protocol).
  std::map<int, vector_real_function_3d> v;
  for (int a : ax) v[a] = dipole_perturbation(world, gs, a);

  if (world.rank() == 0)
    madness::print("=== polarizability (alpha) assembly  protocol_key=", key,
                   " ===");

  const int k_now = madness::FunctionDefaults<3>::get_k();

  for (double w : fr) {
    Tensor<double> alpha(static_cast<long>(ax.size()),
                         static_cast<long>(ax.size()));
    // Per-ROW accuracy: the honest verdict of the state that built each direction
    // row — tied to the SOURCE protocol it came from, NOT the reprojection k below.
    std::vector<std::string> src_key (ax.size());
    std::vector<int>         row_conv(ax.size(), 0);   // converged AND not accepted@maxiter
    std::vector<int>         row_acc (ax.size(), 0);   // accepted@maxiter (converged forced true)
    std::vector<double>      row_res (ax.size(), 0.0); // bsh residual at the source
    bool present = true;   // does every channel have SOME usable (non-diverged) state?

    for (size_t i = 0; i < ax.size() && present; ++i) {
      const std::string chan = Perturbation::dipole(ax[i]).description();
      const std::string fkey = ResponseMetadata::freq_key(w);
      // Pick the SAME source the loader will use (shared best-usable helper) and
      // read its honest accuracy, so the reported accuracy matches the state that
      // actually builds the row. Rank 0 decides; broadcasts keep the collective
      // load + reproject below in lockstep across ranks.
      std::string sk; int cv = 0, ac = 0; double rs = 0.0;
      if (world.rank() == 0) {
        auto meta = ResponseMetadata::load_or_create(
            ctx.calc_dir + "/response_metadata.json");
        sk = best_usable_fd_source_key(meta.json(), chan, fkey, thresh, k_now, key);
        if (!sk.empty()) {
          const auto &e = meta.json()["fd_states"][chan][sk][fkey];
          ac = e.value("accepted",  false) ? 1 : 0;
          cv = (e.value("converged", false) && !ac) ? 1 : 0;
          rs = e.value("bsh_residual", 0.0);
        }
      }
      world.gop.broadcast_serializable(sk, 0);
      world.gop.broadcast(cv, 0);
      world.gop.broadcast(ac, 0);
      world.gop.broadcast(rs, 0);

      if (sk.empty()) {   // no usable state at all -> cannot form this row
        if (world.rank() == 0)
          madness::print("[ALPHA] omega=", w, " dir=", beta_axis_name(ax[i]),
                         " — no usable FD state; tensor not assembled");
        present = false;
        break;
      }
      src_key[i] = sk; row_conv[i] = cv; row_acc[i] = ac; row_res[i] = rs;

      auto Xi = detail_exec::load_fd_as_xy<ClosedShell>(
          world, ctx.calc_dir, Perturbation::dipole(ax[i]), w);
      if (!Xi) { present = false; break; }   // metadata usable but archive gone

      // k-CONSISTENCY: try_load_fd_state re-projects any coarser source to the
      // active (k, thresh) inside the loader, so Xi is at k_now here — inner()
      // below is well-defined no matter which protocol the source came from.
      // The per-row verdict recorded above stays tied to the source protocol.

      for (size_t j = 0; j < ax.size(); ++j)
        alpha(static_cast<long>(i), static_cast<long>(j)) =
            -2.0 * (inner(Xi->x_alpha, v[ax[j]]) + inner(Xi->y_alpha, v[ax[j]]));
    }
    if (!present) continue;   // a whole direction had NO usable state at all

    bool all_conv = true; double max_res = 0.0;
    for (size_t i = 0; i < ax.size(); ++i) {
      all_conv = all_conv && (row_conv[i] != 0);
      if (row_res[i] > max_res) max_res = row_res[i];
    }

    if (world.rank() == 0) {
      madness::print("[ALPHA] omega=", w, "  directions=",
                     [&]{ std::string d; for (int a : ax) d += beta_axis_name(a); return d; }(),
                     "  converged=", all_conv, "  max_bsh_res=", max_res);
      for (size_t i = 0; i < ax.size(); ++i) {
        madness::print("[ALPHA]  row dir=", beta_axis_name(ax[i]),
                       " source_protocol=", src_key[i],
                       " converged=", (row_conv[i] != 0),
                       (row_acc[i] ? " (ACCEPTED@maxiter)" : ""),
                       " bsh_res=", row_res[i]);
        for (size_t j = 0; j < ax.size(); ++j)
          madness::print("[ALPHA]  alpha_", beta_axis_name(ax[i]),
                         beta_axis_name(ax[j]), " (omega=", w, ") =",
                         alpha(static_cast<long>(i), static_cast<long>(j)));
      }
      madness::print("[ALPHA] tensor(omega=", w, ") =\n", alpha);

      auto meta = ResponseMetadata::load_or_create(
          ctx.calc_dir + "/response_metadata.json");
      nlohmann::json row;
      row["omega"] = w;
      std::string dirs; for (int a : ax) dirs += beta_axis_name(a);
      row["directions"] = dirs;
      nlohmann::json mat = nlohmann::json::array();
      for (size_t i = 0; i < ax.size(); ++i) {
        nlohmann::json r = nlohmann::json::array();
        for (size_t j = 0; j < ax.size(); ++j)
          r.push_back(alpha(static_cast<long>(i), static_cast<long>(j)));
        mat.push_back(r);
      }
      row["alpha"] = mat;
      // Per-row accuracy (honest: reflects each row's SOURCE protocol, not k_now).
      nlohmann::json racc = nlohmann::json::array();
      for (size_t i = 0; i < ax.size(); ++i) {
        std::string dname; dname += beta_axis_name(ax[i]);
        racc.push_back({{"direction", dname},
                        {"source_protocol_key", src_key[i]},
                        {"converged", row_conv[i] != 0},
                        {"accepted",  row_acc[i]  != 0},
                        {"bsh_residual", row_res[i]}});
      }
      row["row_accuracy"]     = racc;
      row["converged"]        = all_conv;    // tensor-level: every row met dconv
      row["max_bsh_residual"] = max_res;
      meta.add_property("alpha", key, row);
      meta.save();
    }
    world.gop.fence();
  }
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_CALC_CALC_EXECUTOR_HPP
