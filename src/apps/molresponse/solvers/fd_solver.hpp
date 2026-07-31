#ifndef MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP
#define MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP

// =========================================================================
// FDSolver<Type, Shell> — frequency-dependent response, batched over
// (perturbation, frequency) responses at a single protocol step.
//
// step() shape vs ESSolver:
//
//   1. per-channel: V0x, E0x, gamma, rho            (no T0x, no E0x_full)
//   2. (skip — no Lambda assembly)
//   3. (skip — omega fixed, no subspace eigenproblem)
//   4. (skip — no rotation)
//   5. Theta = V0x - E0x + gamma + perturbation_source
//   6. BSH apply at fixed omega for this channel → new_x [+ new_y]
//      residual = || (x,y)_old - (x,y)_new ||
//   7. explosion guard + banner
//
// Same `ConvergencePolicy`, same `iterate` / `iterate_protocol`,
// same Storage-typed kernel signatures as ESSolver. The differences
// are concentrated here: no T0x / Lambda / sygv / rotation, BSH ops
// could be cached per response (deferred — built inside K::bsh_apply
// for now, same as ES).
//
// TODO(streaming): step() currently materializes V0x, E0x, gamma as
// separate vectors of Storage and assembles theta from them — peak
// is ~6M Storage instances in flight (in + out + V0x + E0x + gamma
// + theta). Refactor to stream the assembly per response: theta starts
// as V0x, then E0x is computed as a temporary and subtracted in
// place, then gamma is added in place, then perturbation_source is
// added. Peak drops to in + out + theta + ~1 temporary ≈ 2M + 3.
// This also sets up future per-response MacroTaskQ parallelism: each
// worker owns 1 response's buffers, not the whole bundle.
//
// TODO(VBC): for cubic response properties (β, dα/dQ for Raman) we
// need a Kernels<Type, Shell>::compute_second_order_vector kernel
// mirroring molresponse_v2/VBCMacrotask.hpp::compute_vbc_i. Pattern:
// given two first-order response vectors B and C, build ζ_BC =
// Σ y_B·x_C, apply compute_g (a 3-input J+K-like operator), apply
// V_pert·c_x, project onto virtuals, combine. Same shape but on 2×
// the storage. Plan the Kernels API so this plugs in cleanly.
// =========================================================================

#include "../kernels/exchange_ctx.hpp"  // exch::assemble_theta_tensor (--fd-tensor gate)
#include "../kernels/full.hpp"     // Kernels<Full, ClosedShell>
#include "../kernels/static.hpp"   // Kernels<Static, ClosedShell>
#include "../kernels/tags.hpp"
#include "convergence_policy.hpp"
#include "fd_problem.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"

#include <madness/mra/mra.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace molresponse_v3 {

template <typename Type, typename Shell>
class FDSolver {
public:
  using K        = Kernels<Type, Shell>;
  using Storage  = StorageOf_t<Type, Shell>;
  using Channel  = FDPerturbationOf_t<Type, Shell>;
  using vecfuncT = std::vector<madness::real_function_3d>;

  struct State {
    std::vector<Storage>                    responses;
    std::vector<double>                     last_bsh_residual;
    std::vector<madness::real_function_3d>  rho_alpha_prev;
    std::vector<double>                     last_density_residual;
    int                                     iter = 0;
    bool                                    diverged = false;
  };

  FDSolver(madness::World &world, FDProblem<Type, Shell> target,
           ConvergencePolicy policy,
           PrintLevel print_level = PrintLevel::Normal,
           std::string log_prefix = {})
      : world_(world), target_(std::move(target)),
        policy_(policy), kain_(world, policy),
        print_level_(print_level), log_prefix_(std::move(log_prefix)) {
    refresh_convergence_targets();
    print_header();
  }

  // -- accessors / config -------------------------------------------------
  madness::World&     world() const             { return world_; }
  PrintLevel          print_level() const       { return print_level_; }
  void                set_print_level(PrintLevel p) { print_level_ = p; }
  ConvergencePolicy   convergence_policy() const{ return policy_; }
  void                set_convergence_policy(ConvergencePolicy p) {
    policy_ = p;
    kain_.set_policy(p);
    refresh_convergence_targets();
  }
  void                set_target(FDProblem<Type, Shell> t) {
    target_ = std::move(t);
    // KAIN's stored iterates live at the old k/thresh; drop them on
    // protocol-driven target swaps.
    kain_.reset();
  }
  const FDProblem<Type, Shell>& target() const { return target_; }

  /// Enable per-iter convergence logging to a CSV file. One row per
  /// (iter, channel) appended at the end of each step(). Pass empty
  /// string to disable.
  void set_log_path(const std::string &path) {
    log_path_ = path;
    log_header_written_ = false;
  }

  void refresh_convergence_targets() {
    targets_ = policy_.effective_for_thresh(
        madness::FunctionDefaults<3>::get_thresh());
  }
  const ConvergencePolicy::Targets& targets() const { return targets_; }

private:
  // F2d (doc 32 §5.6): gated print. Centralizes the Normal/rank-0 guard and the
  // per-subworld tag. EMPTY prefix ⇒ print(...) verbatim (G=0 byte-identical);
  // non-empty ⇒ the "[g{gid}/{G} {host}] " tag leads each line.
  template <class... Ts> void plog(Ts &&...a) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    if (log_prefix_.empty()) print(std::forward<Ts>(a)...);
    else                     print(log_prefix_, std::forward<Ts>(a)...);
  }

  void print_header() const {
    plog("");
    plog("FDSolver<", type_name(), ", ClosedShell>  n_responses =",
          target_.n_responses(),
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", target_.gs.c_xc);
    plog("  policy: dconv_user =", policy_.dconv_user,
          " bsh_target =", targets_.bsh_residual,
          " density_target =", targets_.density_residual);
    plog("  responses:");
    for (int c = 0; c < target_.n_responses(); ++c) {
      plog("    [", c, "]  omega =", target_.responses[c].omega);
    }
    plog("");
  }

  static constexpr const char* type_name() {
    if constexpr (std::is_same_v<Type, Static>) return "Static";
    if constexpr (std::is_same_v<Type, Full>)   return "Full";
    if constexpr (std::is_same_v<Type, TDA>)    return "TDA";
    return "Unknown";
  }

  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double max_res  = 0.0;
    for (double r : out.last_bsh_residual) max_res = std::max(max_res, r);
    double max_drho = 0.0;
    for (double r : out.last_density_residual) max_drho = std::max(max_drho, r);
    plog("iter", out.iter,
          "  max_res =", max_res, "  max_dρ =", max_drho);
    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
        double dr = (c < out.last_density_residual.size())
                        ? out.last_density_residual[c] : 0.0;
        plog("  ch", c, "  omega =", target_.responses[c].omega,
              "  res =", out.last_bsh_residual[c],
              "  dρ =", dr);
      }
    }
  }

  /// Add the perturbation source to theta in place. Differs across
  /// Storage shapes — Storage = ResponseStateX adds to x_alpha only;
  /// ResponseStateXY adds to both x_alpha and y_alpha. For symmetric
  /// (real, time-independent) dipole perturbations the X and Y source
  /// are equal, so caller fills source.x_alpha == source.y_alpha for
  /// Full; for Static there's only one component anyway.
  static void add_perturbation_source(madness::World &world,
                                       Storage &theta,
                                       const Channel &ch) {
    // X-channel α
    madness::gaxpy(world, 1.0, theta.x_alpha, 1.0, ch.source.x_alpha);
    madness::truncate(world, theta.x_alpha);
    // Y-channel α (Full only)
    if constexpr (std::is_same_v<Type, Full>) {
      madness::gaxpy(world, 1.0, theta.y_alpha, 1.0, ch.source.y_alpha);
      madness::truncate(world, theta.y_alpha);
    }
    // β-spin (OpenShell only): X-channel β, and Y-channel β for Full.
    if constexpr (std::is_same_v<Shell, OpenShell>) {
      madness::gaxpy(world, 1.0, theta.x_beta, 1.0, ch.source.x_beta);
      madness::truncate(world, theta.x_beta);
      if constexpr (std::is_same_v<Type, Full>) {
        madness::gaxpy(world, 1.0, theta.y_beta, 1.0, ch.source.y_beta);
        madness::truncate(world, theta.y_beta);
      }
    }
  }

public:
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    plog("");
    // Diverged-first: converged() returns true on diverged so iterate<>
    // would otherwise mislabel the exit as "Converged".
    if (s.diverged)       plog("Stopped at iter", s.iter,
                                "(diverged — residual exceeded "
                                "explosion guard).");
    else if (converged)   plog("Converged in", s.iter, "iters.");
    else                  plog("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    if (!s.last_bsh_residual.empty()) {
      plog("  residuals   =");
      for (size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
        double dr = (c < s.last_density_residual.size())
                        ? s.last_density_residual[c] : 0.0;
        plog("    ch", c, "  omega =", target_.responses[c].omega,
              "  bsh =", s.last_bsh_residual[c], "  dρ =", dr);
      }
    }
    plog("");
  }

  /// One outer iteration. Streamed-theta layout: for each response
  /// in turn, build theta in place (V0x → += gamma → -= E0x → += source)
  /// then BSH-apply, freeing the intermediates before moving to the next
  /// response. Peak memory ≈ `in + out + 1 theta + 1 temp + 1 new_x`
  /// instead of the materialize-V0x/E0x/gamma-for-all-responses-up-front
  /// layout's `~5M Storage` peak.
  State step(State in) {
    State out;
    out.iter      = in.iter + 1;
    out.responses = in.responses;
    const int M   = target_.n_responses();
    const double thr = madness::FunctionDefaults<3>::get_thresh();

    out.last_density_residual.assign(M, 0.0);
    out.last_bsh_residual.assign(M, 0.0);
    out.rho_alpha_prev.resize(M);

    // Inc-2: build the φ-only g0 exchange tensor ONCE per protocol (cached on the
    // ground state) so assemble_theta_tensor doesn't rebuild it every iteration.
    // Only for the ClosedShell Static/Full tensor path (matches the θ gate below).
    if constexpr (std::is_same_v<typename K::State, ResponseStateX<ClosedShell>> ||
                  std::is_same_v<typename K::State, ResponseStateXY<ClosedShell>>) {
      if (policy_.exchange_tensor && target_.gs.g0_alpha.empty())
        target_.gs.g0_alpha =
            exch::build_g0(world_, target_.gs, thr * 0.1, policy_.exchange_tile);
    }

    for (int r = 0; r < M; ++r) {
      // ρ (one density function per response; kept for next iter's Δρ check)
      auto rho = K::compute_density(world_, target_.gs, in.responses[r]);
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho - in.rho_alpha_prev[r];
        out.last_density_residual[r] = drho.norm2();
      }

      // θ = V0x − E0x + γ. Gate 1 (policy_.exchange_tensor, --fd-tensor): the
      // fused tensor-exchange assembly (exch::, ClosedShell Static/Full only);
      // gate 0 (default) = the per-op reference path, streamed via Storage::axpy,
      // UNTOUCHED. doc 28 §4 Inc 1; A/B-to-thresh between the two gates.
      typename K::State theta;
      bool tensor_done = false;
      if constexpr (std::is_same_v<typename K::State, ResponseStateX<ClosedShell>> ||
                    std::is_same_v<typename K::State, ResponseStateXY<ClosedShell>>) {
        if (policy_.exchange_tensor) {
          theta = exch::assemble_theta_tensor(world_, target_.gs,
                                              in.responses[r], rho,
                                              policy_.exchange_tile);
          tensor_done = true;
        }
      }
      if (!tensor_done) {
        theta = K::compute_V0x(world_, target_.gs, in.responses[r]);
        {
          auto E0x = K::compute_E0x(world_, target_.gs, in.responses[r]);
          theta.axpy(world_, -1.0, E0x);                        // theta -= E0x
        }
        {
          auto gamma = K::compute_gamma(world_, target_.gs,
                                        in.responses[r], rho);
          theta.axpy(world_, +1.0, gamma);                      // theta += γ
        }
      }
      add_perturbation_source(world_, theta, target_.responses[r]); // theta += V_p
      theta.truncate_all(world_, thr);

      // BSH apply
      auto x_new = K::bsh_apply(world_, target_.gs, in.responses[r],
                                theta, target_.responses[r].omega);
      out.last_bsh_residual[r] =
          K::compute_residual_norm(world_, in.responses[r], x_new);

      out.responses[r]      = std::move(x_new);
      out.rho_alpha_prev[r] = std::move(rho);
    }

    // ---- KAIN + step restriction (over flat in/out responses) ------------
    // diag_level pulled from print_level_; gated uniformly across ranks
    // inside ResponseSubspaceKain (collective work is unconditional).
    const int kain_diag =
        (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
    kain_.apply(in.responses, out.responses, kain_diag);

    // Explosion guard
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Converged: per-channel BSH residual < bsh_target AND
  ///            per-channel Δρ < density_target  (iter >= 2)
  ///            OR diverged.
  bool converged(const State &s) const {
    if (s.diverged) return true;
    if (s.iter < policy_.min_iters_before_conv) return false;
    if (s.last_bsh_residual.empty()) return false;

    double mx_bsh = 0.0;
    for (double r : s.last_bsh_residual) mx_bsh = std::max(mx_bsh, r);
    if (mx_bsh >= targets_.bsh_residual) return false;

    if (s.iter <= 1) return true;
    double mx_drho = 0.0;
    for (double r : s.last_density_residual) mx_drho = std::max(mx_drho, r);
    return mx_drho < targets_.density_residual;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int c = 0; c < target_.n_responses(); ++c) {
      s.responses[c].save(world_, path_prefix + ".channel_"
                                   + std::to_string(c));
    }
  }

private:
  /// Append one row per channel at the end of step(). Schema:
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged
  /// `state` here is the channel index; `omega` is the fixed channel
  /// frequency (not changing iter-to-iter like ES). Header is written
  /// once per file. Only rank 0 writes.
  void append_convergence_log(const State &s) {
    if (log_path_.empty() || world_.rank() != 0) return;
    std::ofstream out(log_path_, std::ios::app);
    if (!out) return;
    if (!log_header_written_) {
      out << "iter,protocol_thresh,state,omega,bsh_residual,"
          << "density_residual,diverged\n";
      log_header_written_ = true;
    }
    const double pthr = madness::FunctionDefaults<3>::get_thresh();
    const int M = target_.n_responses();
    out.precision(12);
    for (int c = 0; c < M; ++c) {
      const double bsh  = (static_cast<size_t>(c) < s.last_bsh_residual.size())
                              ? s.last_bsh_residual[c] : 0.0;
      const double drho = (static_cast<size_t>(c) < s.last_density_residual.size())
                              ? s.last_density_residual[c] : 0.0;
      out << s.iter << ',' << pthr << ',' << c << ','
          << target_.responses[c].omega << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << '\n';
    }
  }

  madness::World            &world_;
  FDProblem<Type, Shell>      target_;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  ResponseSubspaceKain<Storage>   kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_prefix_;   // F2d: per-subworld tag ("" = none)
  std::string                log_path_;
  bool                       log_header_written_ = false;
};

// Call-site shape (cf. ESSolver):
//
//   FDProblem<Static, ClosedShell> tgt{...};                // gs + responses
//   FDSolver<Static, ClosedShell> solver(world, tgt,
//                                        ConvergencePolicy{1e-4});
//   typename decltype(solver)::State s0 = ...;             // initial guess
//   auto sf = solvers::iterate_protocol(
//       solver, s0, {1e-4, 1e-6}, prepare_fn,
//       solvers::IterateProtocolPolicy{25});
//   solver.save(sf, "fd_solve");

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_FD_SOLVER_HPP
