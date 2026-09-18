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
#include <limits>
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
    /// Per-channel ||theta||, the BSH source norm at the last iteration. This is
    /// the constant in the error bound the convergence guide states,
    /// |d alpha| <= |c| ||theta|| ||r_B||, which until now nothing in the code
    /// evaluated -- the bound was documented and unevaluable. Recorded per
    /// channel so the budget can be checked at run time.
    std::vector<double>                     last_theta_norm;
    /// Per-channel source overlap sum_i <v phi_i | x_i> (+ y block for Full):
    /// the diagonal property estimate of this channel, tracked per iteration
    /// so the property's own convergence can be read off the log (2026-09-11).
    std::vector<double>                     last_property;
    /// Normalised distance per iteration (newest last), ONE TRACK PER GATED
    /// QUANTITY: [0] = max bsh/bsh_target, [1] = max drho/density_target.
    /// A stall needs every track flat (ConvergencePolicy::plateau); reset when
    /// the thresh changes.
    std::vector<std::vector<double>>        gate_history;
    double                                  gate_thresh = -1.0;
    int                                     iter = 0;
    bool                                    diverged = false;
    /// Set by step() when the residual has plateaued above the targets
    /// (ConvergencePolicy::stall_*). converged() returns true on it so the
    /// loop exits; the executor records it and applies best-effort acceptance.
    bool                                    stalled = false;
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

  /// <z|(A-w)|z> for the operator this solver actually solves, formed directly.
  ///
  /// Read off bsh_apply rather than assumed: the fixed point satisfies
  ///
  ///     (T - eps_i - w) x_i + (V0 x)_i - (E0 x)_i + gamma_i(x) = -v_i
  ///
  /// and theta(z) = V0 z - E0 z + gamma(z) + v collects every term but the
  /// kinetic and diagonal-energy ones, so
  ///
  ///     <z|(A-w)|z> = sum_i [ <z_i|T|z_i> - (eps_i + w) <z_i|z_i> ]
  ///                   + <z|theta(z)> - <z|v>
  ///
  /// No Green's operator and no projector enter, which is the point: the two
  /// earlier attempts at this both went wrong by inferring the quadratic form
  /// from the BSH map (first ignoring its level shift, then assuming Q passed
  /// through the inner product because theta was in the Q space -- compute_V0x
  /// returns V*x unprojected, so it is not). <z_i|T|z_i> comes from
  /// madness::Kinetic, which uses the gradient trick and never applies a
  /// Laplacian to an MRA function.
  ///
  /// Collective. Static only: at finite w the operator is symmetric in the
  /// paired metric, not in this inner product.
  double quadratic_form(const Storage &z, int r) {
    auto rho   = K::compute_density(world_, target_.gs, z);
    auto theta = assemble_theta_for(z, r, rho);
    madness::Kinetic<double, 3> T(world_);
    const auto zb = z.blocks();
    double acc = 0.0;
    for (std::size_t b = 0; b < zb.size(); ++b) {
      const auto &zv = *zb[b];
      if (zv.empty()) continue;
      std::vector<madness::real_function_3d> v(zv.begin(), zv.end());
      madness::Tensor<double> tm = T(v, v);            // collective
      madness::Tensor<double> ov = madness::inner(world_, v, v);
      const double w = target_.responses[r].omega;
      for (std::size_t i = 0; i < v.size(); ++i)
        acc += tm(static_cast<long>(i), static_cast<long>(i))
             - (target_.gs.aeps(static_cast<long>(i)) + w) * ov(static_cast<long>(i));
    }
    return acc + storage_inner(world_, z, theta)
               - source_overlap(world_, z, target_.responses[r]);
  }

  /// Stationary (Hylleraas) estimate of channel `r`'s source overlap <v|x>.
  ///
  /// The equation is (A-w)x = -v, so the functional stationary at the solution
  /// is S[z] = 2<v|z> + <z|(A-w)|z>, with S[x] = <v|x>. The naive estimate
  /// <v|z> is linear in z and therefore first order in the response error; S is
  /// stationary, so the SAME vector gives an error quadratic in it -- exactly
  /// -<dx|(A-w)|dx>, and exactly, because S is quadratic and the expansion
  /// terminates.
  ///
  /// Evaluated at the converged vector itself: no extra BSH cycle, one kernel
  /// build for theta(z) and one gradient pass for the kinetic term.
  double stationary_property(const State &s, int r) {
    if (r < 0 || r >= static_cast<int>(s.responses.size()))
      return std::numeric_limits<double>::quiet_NaN();
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &z = s.responses[r];
      return 2.0 * source_overlap(world_, z, target_.responses[r])
           + quadratic_form(z, r);
    }
  }

  /// Is S actually stationary at the converged vector? Perturb it along a fixed
  /// direction by eps and 2*eps and compare the two deviations: a stationary
  /// functional gives a ratio near 4, a first-order one near 2.
  ///
  /// This is the check the two earlier attempts skipped. It needs no reference
  /// value, no sweep and no curve fit -- it interrogates the functional itself,
  /// on one converged state, so a wrong quadratic form is caught before any
  /// compute is spent on measuring its slope.
  /// NB the size of `eps`. The signal here is eps^2 <d|(A-w)|d>, recovered by
  /// cancellation from terms of order ||theta|| ~ 10; at eps = 1e-3 that signal
  /// is ~1e-6, which is the MRA representation noise at thresh = 1e-6, and the
  /// ratio then measures the noise rather than the functional (the first run of
  /// this check returned ~8, i.e. a cubic law, which a quadratic functional
  /// cannot produce). eps must be large enough that eps^2 clears the floor.
  double stationarity_ratio(const State &s, int r, double eps = 1e-1) {
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &x = s.responses[r];
      const double S0 = stationary_property(s, r);
      // direction: the perturbation source, Q-projected like x and independent
      // of the error already in x
      auto dir = target_.responses[r].source;
      const double dn = storage_norm(world_, dir);
      if (!(dn > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      auto probe = [&](double a) {
        State t = s;
        t.responses[r] = x;
        t.responses[r].axpy(world_, a / dn, dir);
        return stationary_property(t, r);
      };
      const double d1 = std::abs(probe(eps)      - S0);
      const double d2 = std::abs(probe(2.0 * eps) - S0);
      return (d1 > 0.0) ? d2 / d1 : std::numeric_limits<double>::quiet_NaN();
    }
  }

  /// The LINEAR coefficient of S about the converged vector, extracted by a
  /// central difference so no fitting is involved:
  ///
  ///     S(+eps) - S(-eps) = 2 eps alpha,   alpha = 2 <d | (A-w) x + v>
  ///
  /// alpha is twice the projected residual contracted with the probe direction
  /// d (which lies in the Q space, being the perturbation source). If the
  /// quadratic form matched the equation the solver solves, alpha would BE the
  /// residual overlap and would fall with it. Reported beside ||r_B|| so the
  /// two can be compared directly at several convergence levels.
  double linear_coefficient(const State &s, int r, double eps = 1e-1) {
    if constexpr (!std::is_same_v<Type, Static>) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      const Storage &x = s.responses[r];
      auto dir = target_.responses[r].source;
      const double dn = storage_norm(world_, dir);
      if (!(dn > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      auto probe = [&](double a) {
        State t = s;
        t.responses[r] = x;
        t.responses[r].axpy(world_, a / dn, dir);
        return stationary_property(t, r);
      };
      return (probe(eps) - probe(-eps)) / (2.0 * eps);
    }
  }

  /// One parseable line per channel: the naive and stationary estimates, the
  /// residual, ||theta||, the documented error bound evaluated, and the
  /// stationarity self-check. No-op unless the executor asks (--hylleraas).
  void report_property_accuracy(const State &s, const std::string &pert,
                                double freq, double c_conv = -2.0) {
    for (int r = 0; r < target_.n_responses(); ++r) {
      const double naive = (r < static_cast<int>(s.last_property.size()))
                               ? s.last_property[r] : 0.0;
      const double res   = (r < static_cast<int>(s.last_bsh_residual.size()))
                               ? s.last_bsh_residual[r] : 0.0;
      const double tn    = (r < static_cast<int>(s.last_theta_norm.size()))
                               ? s.last_theta_norm[r] : 0.0;
      const double S     = stationary_property(s, r);      // collective
      // two probe sizes: a correct quadratic form gives 4 at BOTH, and the
      // agreement between them is what says the probe is above the noise floor
      const double ratio  = stationarity_ratio(s, r, 1e-1);   // collective
      const double ratio2 = stationarity_ratio(s, r, 5e-2);   // collective
      const double alin   = linear_coefficient(s, r);         // collective
      plog("[PROPACC] pert=", pert, " freq=", freq,
           " thresh=", madness::FunctionDefaults<3>::get_thresh(),
           " k=", madness::FunctionDefaults<3>::get_k(),
           " ch=", r, " iters=", s.iter,
           " res=", res, " theta_norm=", tn,
           " bound=", std::abs(c_conv) * tn * res,
           " naive=", naive, " S=", S, " diff=", S - naive,
           " stat_ratio=", ratio, " stat_ratio_half=", ratio2,
           " lin_coeff=", alin);
    }
  }

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
          "  max_res =", max_res, "  max_dρ =", max_drho,
          "  gate =", out.gate_history.empty() || out.gate_history[0].empty()
                          ? 0.0 : out.gate_history[0].back());
    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
        double dr = (c < out.last_density_residual.size())
                        ? out.last_density_residual[c] : 0.0;
        double pr = (c < out.last_property.size()) ? out.last_property[c] : 0.0;
        plog("  ch", c, "  omega =", target_.responses[c].omega,
              "  res =", out.last_bsh_residual[c],
              "  dρ =", dr, "  <v|x> =", pr);
      }
    }
  }

  /// Source overlap sum over blocks of sum_i <source_i | x_i>: the channel's
  /// diagonal property estimate (alpha_vv up to the convention factor). Cheap
  /// (n_occ inner products per block) and collective.
  static double source_overlap(madness::World &world, const Storage &x,
                               const Channel &ch) {
    double s = 0.0;
    const auto xb = x.blocks();
    const auto vb = ch.source.blocks();
    for (std::size_t b = 0; b < xb.size(); ++b) {
      const auto &xv = *xb[b];
      const auto &vv = *vb[b];
      const std::size_t n = std::min(xv.size(), vv.size());
      if (n == 0) continue;
      madness::Tensor<double> ov = madness::inner(world,
          std::vector<madness::real_function_3d>(xv.begin(), xv.begin() + n),
          std::vector<madness::real_function_3d>(vv.begin(), vv.begin() + n));
      s += ov.sum();
    }
    return s;
  }

  /// theta = V0z - E0z + gamma(z) + V_p for an ARBITRARY z, mirroring the
  /// assembly step() performs inline. The stationary estimate must see the same
  /// operator the iteration solves.
  Storage assemble_theta_for(const Storage &z, int r,
                             const madness::real_function_3d &rho) {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    Storage theta = K::compute_V0x(world_, target_.gs, z);
    {
      auto E0x = K::compute_E0x(world_, target_.gs, z);
      theta.axpy(world_, -1.0, E0x);
    }
    {
      auto gamma = K::compute_gamma(world_, target_.gs, z, rho);
      theta.axpy(world_, +1.0, gamma);
    }
    add_perturbation_source(world_, theta, target_.responses[r]);
    theta.truncate_all(world_, thr);
    return theta;
  }

  /// Inner product of two Storages, block by block. source_overlap is this
  /// against the channel's perturbation; the stationary property estimate and
  /// the error-bound constant need the general form.
  static double storage_inner(madness::World &world, const Storage &a,
                              const Storage &b) {
    double s = 0.0;
    const auto ab = a.blocks();
    const auto bb = b.blocks();
    for (std::size_t k = 0; k < ab.size() && k < bb.size(); ++k) {
      const auto &av = *ab[k];
      const auto &bv = *bb[k];
      const std::size_t n = std::min(av.size(), bv.size());
      if (n == 0) continue;
      madness::Tensor<double> ov = madness::inner(world,
          std::vector<madness::real_function_3d>(av.begin(), av.begin() + n),
          std::vector<madness::real_function_3d>(bv.begin(), bv.begin() + n));
      s += ov.sum();
    }
    return s;
  }

  static double storage_norm(madness::World &world, const Storage &a) {
    return std::sqrt(std::max(0.0, storage_inner(world, a, a)));
  }

  /// Plateau bookkeeping: append this iteration's normalised distance for EACH
  /// gated quantity (bsh, density) to its own track (reset on a thresh change)
  /// and set out.stalled per policy. Tracking them separately is what keeps a
  /// slow-moving density from masking a residual that is still falling.
  void update_stall(State &out, const State &in) const {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    out.gate_history = (in.gate_thresh == thr)
                           ? in.gate_history
                           : std::vector<std::vector<double>>(2);
    out.gate_thresh  = thr;
    if (out.gate_history.size() < 2) out.gate_history.resize(2);
    double g_bsh = 0.0, g_rho = 0.0;
    for (std::size_t c = 0; c < out.last_bsh_residual.size(); ++c) {
      g_bsh = std::max(g_bsh, out.last_bsh_residual[c] / targets_.bsh_residual);
      if (c < out.last_density_residual.size())
        g_rho = std::max(g_rho, out.last_density_residual[c] / targets_.density_residual);
    }
    // iter 1 has no previous density, so drho is 0 and would read as "target
    // met"; mark it as no measurement so the track cannot be compared against it.
    out.gate_history[0].push_back(g_bsh);
    out.gate_history[1].push_back(out.iter > 1 ? g_rho
                                   : std::numeric_limits<double>::infinity());
    out.stalled = !out.diverged && policy_.plateau(out.gate_history);
    if (out.stalled && print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      const std::size_t n = out.gate_history[0].size();
      const std::size_t w = static_cast<std::size_t>(policy_.stall_window);
      plog("[STALL] iter", out.iter, ": every gated quantity flat over", policy_.stall_window,
           "iters (<", 100.0 * policy_.stall_ratio, "% improvement) — bsh",
           out.gate_history[0][n - 1], "vs", out.gate_history[0][n - 1 - w],
           ", density", out.gate_history[1][n - 1], "vs", out.gate_history[1][n - 1 - w],
           "(normalised); stopping.");
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
    else if (s.stalled)   plog("Stopped at iter", s.iter,
                                "(stalled — residual plateau above the "
                                "targets for", policy_.stall_window,
                                "iters; not converged).");
    else if (converged)   plog("Converged in", s.iter, "iters.");
    else                  plog("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    if (!s.last_bsh_residual.empty()) {
      plog("  residuals   =");
      for (size_t c = 0; c < s.last_bsh_residual.size(); ++c) {
        double dr = (c < s.last_density_residual.size())
                        ? s.last_density_residual[c] : 0.0;
        double pr = (c < s.last_property.size()) ? s.last_property[c] : 0.0;
        plog("    ch", c, "  omega =", target_.responses[c].omega,
              "  bsh =", s.last_bsh_residual[c], "  dρ =", dr,
              "  <v|x> =", pr);
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
    out.last_property.assign(M, 0.0);
    out.last_theta_norm.assign(M, 0.0);
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
      // The constant in the error bound docs/guides/convergence.md states,
      // |d alpha| <= |c| ||theta|| ||r_B||. Until now nothing evaluated it, so
      // the bound was documented and uncheckable; theta is already in hand here.
      out.last_theta_norm[r] = storage_norm(world_, theta);

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

    // Property trace (after KAIN: this is the iterate that will be saved).
    for (int r = 0; r < M; ++r)
      out.last_property[r] = source_overlap(world_, out.responses[r],
                                            target_.responses[r]);
    update_stall(out, in);

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Converged: per-channel BSH residual < bsh_target AND
  ///            per-channel Δρ < density_target  (iter >= 2)
  ///            OR diverged.
  bool converged(const State &s) const {
    if (s.diverged) return true;
    if (s.stalled)  return true;   // exit; executor reads State::stalled
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
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged,
  ///   property,gate,stalled
  /// `state` here is the channel index; `omega` is the fixed channel
  /// frequency (not changing iter-to-iter like ES). `property` is the source
  /// overlap <v|x> of the channel, `gate` the normalised gate distance
  /// (plateau detector). Header is written once per file. Only rank 0 writes.
  void append_convergence_log(const State &s) {
    if (log_path_.empty() || world_.rank() != 0) return;
    std::ofstream out(log_path_, std::ios::app);
    if (!out) return;
    if (!log_header_written_) {
      out << "iter,protocol_thresh,state,omega,bsh_residual,"
          << "density_residual,diverged,property,gate,stalled\n";
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
      const double prop = (static_cast<size_t>(c) < s.last_property.size())
                              ? s.last_property[c] : 0.0;
      const double gate = (s.gate_history.empty() || s.gate_history[0].empty())
                              ? 0.0 : s.gate_history[0].back();
      out << s.iter << ',' << pthr << ',' << c << ','
          << target_.responses[c].omega << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << ','
          << prop << ',' << gate << ',' << (s.stalled ? 1 : 0) << '\n';
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
