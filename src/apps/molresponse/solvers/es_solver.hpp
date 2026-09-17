#ifndef MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP

// =========================================================================
// ESSolver<Type, Shell> — one skeleton, six instantiations.
//
// step() body is identical across (Type, Shell). Every per-root
// kernel call has the signature K::foo(world_, gs_, state[,
// intermediate]); the Kernels<Type, Shell> specialization decides
// which fields it reads. Bundle-level ops (subspace matrix,
// diagonalize, rotate) are shared free functions in
// kernels/bundle_helpers.hpp.
//
// FDSolver<Type, Shell> will share the SAME Kernels<T,S> bodies for
// per-root quantities (gamma, V0·x, residual) — the FD/ES axis is in
// the outer iteration shape only.
//
// Convergence policy and per-iter density tracking live here so the
// FD solver can reuse [[ConvergencePolicy]] without duplication.
// =========================================================================

#include "../kernels/assembly.hpp"   // assemble_lambda, assemble_theta
#include "../kernels/full.hpp"  // Kernels<Full, ClosedShell> (ES-capable)
#include "../kernels/response_space_ops.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tda.hpp"   // Kernels<TDA, *>
#include "../kernels/tda_batch.hpp"  // tda_batch:: bundle-batched V0x/T0x/BSH
#include "convergence_policy.hpp"
#include "es_root_identity.hpp"
#include "iterate.hpp"
#include "response_state.hpp"
#include "response_subspace_kain.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <limits>
#include <array>
#include <fstream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

/// CODATA 2018 hartree -> eV; excitation energies are quoted in eV in the
/// literature (e.g. Salek et al.), so every omega print carries both.
inline constexpr double kHartreeToEV = 27.211386245988;


/// Read-only problem definition for ESSolver<Type, Shell>. Symmetric
/// with FDProblem<Type, Shell> in fd_problem.hpp:
///
///   gs       — the prepared ground state (zeroth-order data) the
///              kernels read from. Same object FDSolver uses.
///   n_roots  — solver-level count of excitation roots to track.
///
/// `n_roots` is a solver concern, not a ground-state property, so it
/// lives on the problem alongside gs rather than inside the gs struct.
/// Debug helper: the coordinate function r[d] (for orbital centroids).
struct CoordinateFunctor : public madness::FunctionFunctorInterface<double, 3> {
  int d;
  explicit CoordinateFunctor(int dd) : d(dd) {}
  double operator()(const madness::coord_3d &r) const override { return r[d]; }
};

template <typename Type, typename Shell>
struct ESProblem {
  ResponseGroundState gs;
  int                 n_roots = 0;
};

template <typename Type, typename Shell>
class ESSolver {
public:
  using K       = Kernels<Type, Shell>;
  using Storage = StorageOf_t<Type, Shell>;
  using vecfuncT = std::vector<madness::real_function_3d>;

  struct State {
    std::vector<Storage>                    roots;
    madness::Tensor<double>                 omega;
    /// Stable identity per slot: stable_index[s] is the permanent
    /// identity of the root currently in roots[s]. Empty until assigned
    /// (see ensure_root_identity); travels through sort_state_by_omega
    /// and save/load so a root keeps its id across reorderings/protocols.
    std::vector<int>                        stable_index;
    std::vector<double>                     last_bsh_residual;
    /// Per-root previous-iter density (for Δρ tracking). Empty on
    /// iter 0, populated thereafter; index s tracks the same root
    /// slot as roots[s] before rotation.
    std::vector<madness::real_function_3d>  rho_alpha_prev;
    std::vector<double>                     last_density_residual;
    /// Per-root |ω_new − ω_old| (eigenvalue residual). The ES convergence /
    /// lock criterion is gated on this + density, NOT the BSH amplitude
    /// residual (which jitters for near-degenerate roots). Populated each step.
    std::vector<double>                     last_omega_residual;
    /// Per-root sticky "locked" flag (workstream B, full deflation). When
    /// policy.lock_converged is set, a converged root is locked: excluded from
    /// the subspace + rotation, used only to deflate the active roots, and
    /// skipped by KAIN/step restriction. Empty = none (the default path never
    /// touches it). Cleared at each protocol (wavelet-thresh) change.
    std::vector<char>                       locked;
    /// Normalised distance per iteration (newest last), ONE TRACK PER GATED
    /// QUANTITY: [0] = max drho/density_target, [1] = max |dw|/omega_target,
    /// over the ACTIVE roots. A stall needs both flat; reset on a thresh change.
    std::vector<std::vector<double>>        gate_history;
    double                                  gate_thresh = -1.0;
    int                                     iter = 0;
    /// Set by step() when the explosion guard trips; iterate<>
    /// terminates on this so we don't burn iters on a runaway state.
    bool                                    diverged = false;
    /// Set by step() when the gated residuals (density, |dw|) of the active
    /// roots have plateaued above the targets (ConvergencePolicy::stall_*).
    /// converged() returns true on it so the loop exits; the executor reports
    /// the bundle as NOT converged and records `stalled`.
    bool                                    stalled = false;
  };

  /// Preferred ctor — caller supplies a problem and a policy.
  ESSolver(madness::World &world, ESProblem<Type, Shell> problem,
           ConvergencePolicy policy,
           PrintLevel print_level = PrintLevel::Normal)
      : world_(world), gs_(std::move(problem.gs)),
        n_roots_(problem.n_roots),
        policy_(policy), kain_(world, policy),
        print_level_(print_level) {
    refresh_convergence_targets();
    print_header();
  }

  /// Future variant — low-memory Lambda assembly:
  ///   The current step() materializes T0x, V0x, E0x_full, gamma as
  ///   four separate response_spaces before summing them into Lambda.
  ///   For larger systems where holding three or four response_spaces
  ///   simultaneously is prohibitive, an alternative is to fold each
  ///   term in place: build Lambda by += T0x, += V0x, -= E0x_full, +=
  ///   gamma, freeing each piece as soon as it has been added. That
  ///   would replace `assemble_lambda(...)` with a streaming-accumulator
  ///   form. Leaving the current shape for clarity; revisit when memory
  ///   pressure on a target system demands it.

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

  /// Swap in a fresh ground state — used by iterate_protocol's
  /// prepare callback after the ground-state side has been re-prepped
  /// at a new (k, thresh). Carries everything protocol-sensitive that
  /// isn't part of State (phi0, V_local, coulop, Q, fock). KAIN's
  /// stored iterates are reset because they live at the old wavelet
  /// basis. n_roots is unchanged.
  void set_gs(ResponseGroundState gs) {
    gs_ = std::move(gs);
    kain_.reset();
    lock_pass_.clear();
    lock_protocol_thresh_ = -1.0;  // re-earn locks at the new protocol
  }
  const ResponseGroundState& gs()      const { return gs_; }
  int                        n_roots() const { return n_roots_; }

  /// Stage-1 bundle batching (closed-shell TDA only): compute V0x/T0x and
  /// the BSH apply over the flattened M·n_occ bundle in one pass instead of
  /// M per-root passes. Pure per-function maps ⇒ numerically identical to
  /// the per-root path (see kernels/tda_batch.hpp). Default OFF; --es-batch
  /// flips it. Has effect only for TDA/ClosedShell with M>1.
  void set_batched(bool b) { batched_ = b; }
  bool batched() const     { return batched_; }

  /// Inc-3c tensor-layer γ (closed-shell TDA only): build the φ-only g0
  /// exchange tensor once per protocol (cached on gs_; set_gs invalidates it)
  /// and compute the bundle's γ as one Coulomb wave + per-root g0
  /// contractions instead of M per-root Exchange-operator applies
  /// (kernels/tda_batch.hpp compute_gamma_flat). A/B-to-floor (~1e-3 REL) vs
  /// the reference — NOT bit-identical like set_batched — so it is a
  /// SEPARATE gate (--es-tensor). Default OFF = per-root reference γ.
  void set_gamma_tensor(bool b) { gamma_tensor_ = b; }
  bool gamma_tensor() const     { return gamma_tensor_; }

  /// Per-phase wall-time instrumentation (diagnostic). When on, step_rotate_
  /// pieces fences at phase boundaries and prints one rank-0 line per step:
  ///   ES_TIMING iter=N gamma=.. build=.. subspace=.. rotate=.. bsh=.. total=..
  /// The boundary fences are EXTRA syncs (perturb timing slightly) added only
  /// under this flag — normal runs are untouched. Default OFF; --es-time flips.
  void set_time_phases(bool b) { time_phases_ = b; }
  bool time_phases() const     { return time_phases_; }

  /// Enable per-iter convergence logging to a CSV file. One row per
  /// (iter, state) appended at the end of each step(). Header is
  /// written on the first row of a fresh file. Pass empty string to
  /// disable.
  void set_log_path(const std::string &path) {
    log_path_ = path;
    log_header_written_ = false;
  }

  /// Called by iterate_protocol after FunctionDefaults<3>::thresh has
  /// changed. Re-derives effective BSH / density targets from the
  /// fresh thresh.
  void refresh_convergence_targets() {
    targets_ = policy_.effective_for_thresh(
        madness::FunctionDefaults<3>::get_thresh());
  }
  const ConvergencePolicy::Targets& targets() const { return targets_; }

private:
  void print_header() const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    // Type/Shell are compile-time tags; name them (the banner used to say
    // "TDA, ClosedShell" for every instantiation, which misread Full solves).
    const char *type_name  = std::is_same_v<Type, TDA> ? "TDA" : "Full";
    const char *shell_name = std::is_same_v<Shell, ClosedShell> ? "ClosedShell" : "OpenShell";
    print(std::string("ESSolver<") + type_name + ", " + shell_name + ">  num_roots =", n_roots_,
          " thresh =", madness::FunctionDefaults<3>::get_thresh(),
          " c_xc =", gs_.c_xc);
    print("  policy: dconv_user =", policy_.dconv_user,
          " bsh_target =", targets_.bsh_residual,
          " density_target =", targets_.density_residual,
          " explosion_guard =", policy_.explosion_guard);
    print("");
  }

  void print_iter_banner(const State &out) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    double max_res  = 0.0;
    for (double r : out.last_bsh_residual) max_res = std::max(max_res, r);
    double max_drho = 0.0;
    for (double r : out.last_density_residual) max_drho = std::max(max_drho, r);
    print("iter", out.iter, "  omega =", out.omega,
          "  max_res =", max_res, "  max_dρ =", max_drho,
          "  gate =", out.gate_history.empty() || out.gate_history[0].empty()
                          ? 0.0 : out.gate_history[0].back());
    print("             omega(eV) =", out.omega * kHartreeToEV);

    if (print_level_ >= PrintLevel::Verbose) {
      for (size_t s = 0; s < out.last_bsh_residual.size(); ++s) {
        double dr = (s < out.last_density_residual.size())
                        ? out.last_density_residual[s] : 0.0;
        print("  root", s, "  res =", out.last_bsh_residual[s],
              "  dρ =", dr,
              "  omega =", out.omega(static_cast<long>(s)));
      }
    }
  }

  // Debug-tier per-iter dumps: A, S, U, omega tensors + selected
  // <X|.|X> inner-product matrices.
  void print_debug_iter(const madness::Tensor<double> &A,
                        const madness::Tensor<double> &S_mat,
                        const madness::Tensor<double> &omega_new,
                        const madness::Tensor<double> &U) const {
    if (print_level_ < PrintLevel::Debug || world_.rank() != 0) return;
    print("[DEBUG] A (subspace):");  print(A);
    print("[DEBUG] S (overlap):");    print(S_mat);
    print("[DEBUG] omega:");          print(omega_new);
    print("[DEBUG] U:");              print(U);
  }

  /// Slot-identity diagnostic — emitted per iter at Verbose+. Shows
  /// whether dominance-sort swapped slots (perm != identity) and
  /// whether the post-fixup U(i,i) stayed near ±1 (slot identity
  /// preserved) or dropped (slot picked up character from another).
  /// Pair this with the per-state residual columns to correlate
  /// slot-flip events with KAIN history corruption.
  void print_rot_slots(int iter,
                       const rs::DiagonalizeResult &r) const {
    if (print_level_ < PrintLevel::Verbose || world_.rank() != 0) return;
    const long M = static_cast<long>(r.dominance_perm.size());
    bool ident = true;
    for (long i = 0; i < M; ++i)
      if (r.dominance_perm[i] != i) { ident = false; break; }
    printf("[ROT-SLOTS] iter=%d  perm=[", iter);
    for (long i = 0; i < M; ++i)
      printf(" %ld", r.dominance_perm[i]);
    printf(" ]  U_diag=[");
    for (long i = 0; i < r.U_diag.dim(0); ++i)
      printf(" %+.3f", r.U_diag(i));
    printf(" ]  omega_asc=[");
    for (long i = 0; i < r.omega_ascending.dim(0); ++i)
      printf(" %+.5f", r.omega_ascending(i));
    printf(" ]%s\n", ident ? "  identity" : "  REORDERED");
    fflush(stdout);
  }

  void print_debug_norms(int s, const State &in,
                         const std::vector<madness::real_function_3d> &V0x_s,
                         double t_diag,   // ½Σ|∇x|² (gradient-form kinetic diagonal)
                         const std::vector<madness::real_function_3d> &gamma_s,
                         const std::vector<madness::real_function_3d> &theta_s) const {
    if (print_level_ < PrintLevel::Debug) return;
    // Collective: every rank must call norm2.
    double nx = madness::norm2(world_, in.roots[s].x_alpha);
    double nv = madness::norm2(world_, V0x_s);
    double nt = t_diag;
    double ng = madness::norm2(world_, gamma_s);
    double nh = madness::norm2(world_, theta_s);
    if (world_.rank() != 0) return;
    printf("[NORMS] iter=%d root=%d  |x|=%.3e  |V0x|=%.3e  <x|T|x>=%.3e  "
           "|gamma|=%.3e  |theta|=%.3e\n",
           in.iter, s, nx, nv, nt, ng, nh);
    fflush(stdout);
  }

public:
  /// Called by iterate<> after the loop exits. Reports converged or
  /// max-iter-reached banner.
  void print_final(const State &s, bool converged) const {
    if (print_level_ < PrintLevel::Normal || world_.rank() != 0) return;
    print("");
    // iterate<> sets `converged=true` whenever solver.converged()
    // returns true. converged() also returns true on diverged so the
    // loop exits — so distinguish here by checking s.diverged first.
    if (s.diverged)       print("Stopped at iter", s.iter,
                                "(diverged — residual exceeded "
                                "explosion guard).");
    else if (s.stalled)   print("Stopped at iter", s.iter,
                                "(stalled — density/omega plateau above the "
                                "targets for", policy_.stall_window,
                                "iters; not converged).");
    else if (converged)   print("Converged in", s.iter, "iters.");
    else                  print("Stopped at iter", s.iter,
                                "(max iters reached, not converged).");
    print("  omega_final =", s.omega);
    print("  omega_final(eV) =", s.omega * kHartreeToEV);
    if (!s.last_bsh_residual.empty()) {
      print("  residuals   =");
      for (size_t i = 0; i < s.last_bsh_residual.size(); ++i) {
        double dr = (i < s.last_density_residual.size())
                        ? s.last_density_residual[i] : 0.0;
        print("    root", i, ":  bsh =", s.last_bsh_residual[i],
              "  dρ =", dr);
      }
    }
    print("");
  }

  /// One outer iteration. Dispatches on `policy_.stream_theta`:
  ///   - false (default): step_rotate_pieces() — keeps V0x/E0x/gamma
  ///     for all roots in memory, rotates them by U, assembles Theta
  ///     from the rotated pieces. ~8M Storage peak. Faster (no
  ///     recompute) but heavier.
  ///   - true: step_recompute_pieces() — drops V0x/E0x/gamma after
  ///     building Lambda, rotates only X, then recomputes V0x/E0x/
  ///     gamma against the rotated X to stream Theta in place.
  ///     ~3-4M Storage peak. ~2× the per-iter kernel work but big
  ///     memory savings.
  State step(State in) {
    if (policy_.lock_converged)            // workstream B: full-deflation locking
      return step_rotate_pieces_locked(std::move(in));
    return policy_.stream_theta
        ? step_recompute_pieces(std::move(in))
        : step_rotate_pieces(std::move(in));
  }

private:
  /// Debug: the four blocks that sum to the Ritz matrix, so an indefinite A can
  /// be attributed to a piece. T = ½<∇x|∇x> (gradient form), V = <x|V0x>,
  /// E = <x|E0x_full>, G = <x|γ>; A = T + V − E + G. Collective (inner
  /// products); prints on rank 0 at Debug.
  void print_debug_pieces(const std::vector<Storage> &roots,
                          const madness::Tensor<double> &T,
                          const std::vector<Storage> &V0x,
                          const std::vector<Storage> &E0x_full,
                          const std::vector<Storage> &gamma) const {
    if (print_level_ < PrintLevel::Debug) return;
    auto V = rs::inner(roots, V0x);
    auto E = rs::inner(roots, E0x_full);
    auto G = rs::inner(roots, gamma);
    // Block norms per root (collective): |x|^2 and, for paired states, |y|^2 —
    // the RPA metric normalizes |x|^2 - |y|^2 = 1, so these show how far the
    // vectors have moved from the TDA limit.
    std::vector<double> nx(roots.size(), 0.0), ny(roots.size(), 0.0);
    for (std::size_t i = 0; i < roots.size(); ++i) {
      nx[i] = madness::inner(world_, roots[i].x_alpha, roots[i].x_alpha).sum();
      if constexpr (rs::detail::has_y_alpha<Storage>::value)
        ny[i] = madness::inner(world_, roots[i].y_alpha, roots[i].y_alpha).sum();
    }
    // (a) Occupied overlaps of the (projected) roots: Q x should leave
    //     max_i |<phi_i|x_j>| ~ thresh. Collective.
    std::vector<double> ovx(roots.size(), 0.0), ovy(roots.size(), 0.0);
    for (std::size_t j = 0; j < roots.size(); ++j) {
      auto mx = madness::matrix_inner(world_, gs_.amo, roots[j].x_alpha);
      ovx[j] = mx.absmax();
      if constexpr (rs::detail::has_y_alpha<Storage>::value) {
        auto my = madness::matrix_inner(world_, gs_.amo, roots[j].y_alpha);
        ovy[j] = my.absmax();
      }
    }
    // (b) Self-check of the Fock pieces on the occupied orbitals themselves:
    //     1/2<grad phi|grad phi> + <phi|V_local phi> - <phi|K0 phi> must equal the
    //     Fock diagonal F_ii (g0.focka). A mismatch means the pieces, not the
    //     vectors, are wrong. Collective.
    std::vector<double> fchk_T(gs_.amo.size()), fchk_V(gs_.amo.size()), fchk_K(gs_.amo.size());
    {
      const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
      auto Vphi = mul_sparse(world_, gs_.V_local_alpha, gs_.amo, vtol);
      auto Kphi = common_ops::apply_ground_exchange(world_, gs_.K0_alpha, gs_.amo, gs_.amo, gs_.lo);
      for (int d = 0; d < 3; ++d) {
        madness::real_derivative_3d D(world_, d);
        auto g = apply(world_, D, gs_.amo);
        for (std::size_t i = 0; i < gs_.amo.size(); ++i) fchk_T[i] += 0.5 * madness::inner(g[i], g[i]);
      }
      for (std::size_t i = 0; i < gs_.amo.size(); ++i) {
        fchk_V[i] = madness::inner(gs_.amo[i], Vphi[i]);
        fchk_K[i] = madness::inner(gs_.amo[i], Kphi[i]);
      }
    }
    // (c) Consistency of the ground-state object itself: orbital Gram matrix,
    //     projector applied to its own orbitals, |V_local|, orbital centroids
    //     (x,y,z) and <phi|V_local|phi> through a second code path (operator*).
    madness::Tensor<double> gram = madness::matrix_inner(world_, gs_.amo, gs_.amo);
    double qamo_max = 0.0;
    {
      auto qphi = gs_.Qa(gs_.amo);
      auto ov = madness::matrix_inner(world_, gs_.amo, qphi);
      qamo_max = ov.absmax();
    }
    const double vloc_norm = gs_.V_local_alpha.norm2();
    std::vector<double> vphi2(gs_.amo.size(), 0.0);
    std::vector<std::array<double,3>> cen(gs_.amo.size());
    {
      for (int d = 0; d < 3; ++d) {
        std::shared_ptr<madness::FunctionFunctorInterface<double, 3>> cf(new CoordinateFunctor(d));
        madness::real_function_3d rd = madness::real_factory_3d(world_).functor(cf);
        for (std::size_t i = 0; i < gs_.amo.size(); ++i) {
          auto sq = gs_.amo[i] * gs_.amo[i];
          cen[i][d] = madness::inner(sq, rd);
        }
      }
      for (std::size_t i = 0; i < gs_.amo.size(); ++i) {
        auto vp = gs_.V_local_alpha * gs_.amo[i];
        vphi2[i] = madness::inner(gs_.amo[i], vp);
      }
    }
    // (d) V_local sampled along the z axis (x=y=0): a neutral molecule's
    //     V_nuc + J must be strongly negative near the nuclei and vanish far out.
    std::vector<double> zs = {-3.0, -1.0, -0.3, 0.0, 0.3, 1.0, 1.5, 2.0, 2.7, 3.5, 5.0, 8.0, 15.0, 40.0};
    std::vector<double> vz(zs.size(), 0.0);
    {
      auto vl = madness::copy(gs_.V_local_alpha);
      vl.reconstruct();
      for (std::size_t i = 0; i < zs.size(); ++i) vz[i] = vl(madness::coord_3d{0.0, 0.0, zs[i]});
    }
    if (world_.rank() != 0) return;
    printf("[DEBUG] GS object: |V_local|_2 = %.6e   max|<phi|Q phi>| = %.3e\n", vloc_norm, qamo_max);
    printf("[DEBUG] V_local(0,0,z):");
    for (std::size_t i = 0; i < zs.size(); ++i) printf("  z=%.1f:%.4f", zs[i], vz[i]);
    printf("\n");
    printf("[DEBUG] GS orbital Gram <phi_i|phi_j>:\n"); print(gram);
    printf("[DEBUG] GS orbital centroids (i: x y z) and <phi|V_local|phi> via operator*:\n");
    for (std::size_t i = 0; i < gs_.amo.size(); ++i)
      printf("   %2zu: %9.4f %9.4f %9.4f   Vphi*=%10.5f\n", i, cen[i][0], cen[i][1], cen[i][2], vphi2[i]);
    printf("[DEBUG] Fock self-check on occupied orbitals (i: T V -K | sum | focka_ii):\n");
    for (std::size_t i = 0; i < gs_.amo.size(); ++i)
      printf("   %2zu: %10.5f %10.5f %10.5f | %10.5f | %10.5f\n", i, fchk_T[i], fchk_V[i], -fchk_K[i],
             fchk_T[i] + fchk_V[i] - fchk_K[i], gs_.focka(long(i), long(i)));
    printf("[DEBUG] occupied overlaps after projection (root: max|<phi|x>| max|<phi|y>|):\n");
    for (std::size_t j = 0; j < roots.size(); ++j) printf("   %2zu: %10.3e %10.3e\n", j, ovx[j], ovy[j]);
    printf("[DEBUG] block norms (root: |x|^2 |y|^2):\n");
    for (std::size_t i = 0; i < roots.size(); ++i)
      printf("   %2zu: %10.5f %10.5f\n", i, nx[i], ny[i]);
    print("[DEBUG] T = 1/2<grad x|grad x>:"); print(T);
    print("[DEBUG] V = <x|V0 x>:");           print(V);
    print("[DEBUG] E = <x|E0_full x>:");      print(E);
    print("[DEBUG] G = <x|gamma>:");          print(G);
    const long n = T.dim(0);
    printf("[DEBUG] diag pieces (root: T V E G  | T+V-E+G):\n");
    for (long i = 0; i < n; ++i)
      printf("   %2ld: %10.5f %10.5f %10.5f %10.5f  | %10.5f\n", i, T(i,i), V(i,i), E(i,i), G(i,i),
             T(i,i)+V(i,i)-E(i,i)+G(i,i));
    fflush(stdout);
  }

  /// KAIN hold-off (policy.kain_min_residual): true when THIS iteration's raw
  /// BSH residual is still too large for a subspace extrapolation to be
  /// trustworthy. Rank-uniform (residuals are replicated). Prints at Verbose.
  bool kain_held_off(const State &out) const {
    if (policy_.kain_min_residual <= 0.0 || out.last_bsh_residual.empty()) return false;
    double mx = 0.0;
    for (double r : out.last_bsh_residual) mx = std::max(mx, r);
    const bool hold = mx > policy_.kain_min_residual;
    if (hold && print_level_ >= PrintLevel::Verbose && world_.rank() == 0)
      printf("[KAIN] iter=%d HOLD (raw BSH max_res=%.3e > kain_min_residual=%.3e): plain BSH step, no history\n",
             out.iter, mx, policy_.kain_min_residual);
    return hold;
  }

  /// Final stage shared by both variants — KAIN apply, explosion
  /// guard, banner, log. Mutates out in place.
  void finalize_iter(const State &in, State &out,
                     const std::vector<Storage> &x_pre_bsh) {
    (void)in;
    const int kain_diag =
        (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
    if (!kain_held_off(out)) kain_.apply(x_pre_bsh, out.roots, kain_diag);
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }
    update_stall(out, in);
    print_iter_banner(out);
    append_convergence_log(out);
  }

public:
  /// "Keep pieces, rotate them, assemble Theta from rotated pieces."
  /// Default variant — fastest, highest memory.
  /// Algorithm:
  ///   1. Per root, build V0x, E0x_full, E0x (no-diag), gamma.
  ///   2. Assemble Lambda = V0x - E0x_full + gamma (no kinetic term).
  ///   3. Build A = <X | Lambda> + ½<∇X|∇X> (rs::kinetic_gram), S = <X | X>,
  ///      diagonalize → omega, U.
  ///   4. Rotate {X, V0x, E0x, gamma} by U  (T0x and E0x_full are
  ///      Lambda-only and can be discarded).
  ///   5. Assemble Theta = V0x - E0x + gamma  (from rotated pieces).
  ///   6. BSH apply → new X; residual = || X_old - X_new ||.
  State step_rotate_pieces(State in) {
    State out;
    out.iter  = in.iter + 1;
    out.roots = in.roots;
    // Root identity travels with the state: dropping it here erased every
    // root_id on the first step and made restarts silently relabel roots
    // (the locked variant always carried it — this is the same contract).
    out.stable_index = in.stable_index;
    const int M = n_roots_;

    // ---- per-phase wall timing (diagnostic, --es-time) --------------------
    // MADNESS ops are async (queued, run at the next fence), so attributing
    // time needs a fence at each phase boundary. `lap(acc)` fences, adds the
    // elapsed since the last mark to `acc`, and advances the mark. The fences
    // are added ONLY when time_phases_ is on; otherwise lap() is a no-op.
    double t_gamma = 0, t_build = 0, t_subspace = 0, t_rotate = 0, t_bsh = 0;
    double tmark = 0.0;
    auto lap = [&](double &acc) {
      if (!time_phases_) return;
      world_.gop.fence();
      double now = madness::wall_time();
      acc += now - tmark;
      tmark = now;
    };
    if (time_phases_) { world_.gop.fence(); tmark = madness::wall_time(); }

    // ---- 0. top-of-iter discipline ----------------------------------------
    // Strip ground-orbital contamination that accumulated during the
    // previous iter's BSH + KAIN steps (else BSH amplifies it into the
    // ω = -ε_core ghost eigenvector), then re-orthonormalize so the
    // subspace step sees a clean basis.
    rs::project(out.roots, gs_.Qa, gs_.Qb);
    rs::orthonormalize(world_, out.roots);

    // ---- 1. per-root building blocks --------------------------------------
    // Intermediates are Storage-typed: TDA-ClosedShell carries only
    // x_alpha; Full-ClosedShell will carry x_alpha + y_alpha; OpenShell
    // adds the *_beta members. Kernel signatures don't change between
    // (Type, Shell); the wrapping is invisible to step().
    std::vector<Storage> V0x(M), E0x_full(M), E0x(M), gamma(M);
    std::vector<madness::real_function_3d> rho_alpha(M);
    out.last_density_residual.assign(M, 0.0);

    // Stage-1 batching: build V0x for ALL roots in one bundle pass
    // (pure per-function maps — identical to the per-root calls). gamma /
    // E0x / density stay per-root (see kernels/tda_batch.hpp for why).
    bool did_batch_v0t0 = false;
    if constexpr (std::is_same_v<Type, TDA> &&
                  std::is_same_v<Shell, ClosedShell>) {
      if (batched_ && M > 1) {
        const std::size_t n = out.roots[0].x_alpha.size();
        auto Xf  = tda_batch::flatten(out.roots);
        auto V0f = tda_batch::compute_V0x_flat(world_, gs_, Xf);
        tda_batch::unflatten_into(V0f, n, V0x);
        did_batch_v0t0 = true;
      }
    }
    lap(t_build);  // top-of-iter project/ortho + batched V0x/T0x

    // gamma pass (density + Coulomb-exchange γ) — the suspected hot phase, so
    // timed on its own. Split from the ops pass below ONLY for attribution;
    // the per-root kernel calls are independent, so results are unchanged.
    //
    // Gate 1 (gamma_tensor_, --es-tensor; TDA/ClosedShell only): the Inc-3c
    // tensor-layer γ — densities per root (unchanged), then g0 built once per
    // protocol (cached on gs_; set_gs hands us a fresh GS each protocol so the
    // empty-guard invalidates correctly) + compute_gamma_flat (one Coulomb
    // wave + per-root g0 contractions; NO per-iter exchange convolutions).
    // Gate 0 (default) = the per-root reference loop, untouched.
    bool did_gamma_tensor = false;
    if constexpr (std::is_same_v<Type, TDA> &&
                  std::is_same_v<Shell, ClosedShell>) {
      if (gamma_tensor_) {
        for (int s = 0; s < M; ++s) {
          rho_alpha[s] = K::compute_density(world_, gs_, out.roots[s]);
          if (!in.rho_alpha_prev.empty()) {
            auto drho = rho_alpha[s] - in.rho_alpha_prev[s];
            out.last_density_residual[s] = drho.norm2();
          }
        }
        if (gs_.g0_alpha.empty())
          gs_.g0_alpha = exch::build_g0(
              world_, gs_, madness::FunctionDefaults<3>::get_thresh() * 0.1,
              policy_.exchange_tile);  // untiled n_occ² wave crashes at n_occ≥34
        auto gflat = tda_batch::compute_gamma_flat(world_, gs_, out.roots,
                                                   rho_alpha, gs_.g0_alpha);
        for (int s = 0; s < M; ++s) gamma[s] = std::move(gflat[s]);
        did_gamma_tensor = true;
      }
    }
    if (!did_gamma_tensor) {
      for (int s = 0; s < M; ++s) {
        rho_alpha[s] = K::compute_density(world_, gs_, out.roots[s]);
        gamma[s]     = K::compute_gamma(world_, gs_, out.roots[s], rho_alpha[s]);
        // Δρ vs previous iter (per slot — see note in State).
        if (!in.rho_alpha_prev.empty()) {
          auto drho = rho_alpha[s] - in.rho_alpha_prev[s];
          out.last_density_residual[s] = drho.norm2();
        }
      }
    }
    out.rho_alpha_prev = std::move(rho_alpha);
    lap(t_gamma);

    // ops pass (per-root V0x if not batched, + E0x_full, E0x).
    for (int s = 0; s < M; ++s) {
      if (!did_batch_v0t0) {
        V0x[s]    = K::compute_V0x(world_, gs_, out.roots[s]);
      }
      E0x_full[s] = K::compute_E0x_full(world_, gs_, out.roots[s]);
      E0x[s]      = K::compute_E0x(world_, gs_, out.roots[s]);
    }

    // ---- 2. assemble Lambda per root --------------------------------------
    // Λ here EXCLUDES the kinetic term; it enters A in gradient form below.
    std::vector<Storage> lambda(M);
    for (int s = 0; s < M; ++s) {
      lambda[s] = assemble_lambda(world_, V0x[s], E0x_full[s], gamma[s]);
    }
    lap(t_build);  // ops pass + Lambda assembly

    // ---- 3. subspace A, S, diagonalize ------------------------------------
    // rs::inner / rs::metric / rs::transform take std::vector<State>
    // directly and flatten α (+β for OpenShell) under the hood. The
    // subspace matrix A uses the Euclidean (+) inner:
    //   A_ij = ⟨x_α_i|Λ_α_j⟩ (+ ⟨y_α_i|Λ_y_j⟩ for Full, + β blocks for
    //   OpenShell).
    // The overlap S uses rs::metric — the indefinite RPA metric
    // M = diag(I_X, −I_Y):
    //   S_ij = ⟨X_i|X_j⟩ − ⟨Y_i|Y_j⟩   (Full)
    //   S_ij = ⟨X_i|X_j⟩               (TDA — no Y block, M = I, so
    //                                    rs::metric == rs::inner here).
    // This is the RPA subspace eigenproblem A c = ω S c with the
    // correct symplectic metric (legacy ExcitedResponse.cpp:871). The
    // rotation U mixes ROOTS not spins, so applying it to the flat
    // concat equals applying per spin block.
    // Kinetic block in gradient form (symmetric PSD); see rs::kinetic_gram.
    const auto Tg = rs::kinetic_gram(world_, out.roots);
    auto A     = rs::inner(out.roots, lambda);
    A += Tg;
    print_debug_pieces(out.roots, Tg, V0x, E0x_full, gamma);
    auto S_mat = rs::metric(out.roots, out.roots);
    auto diag_result = rs::diagonalize(A, S_mat,
                                       /*thresh_degenerate=*/-1.0,
                                       policy_.cluster_unmix_factor,
                                       &world_);  // subworld-safe sygvp
    auto &omega_new = diag_result.omega;
    auto &U         = diag_result.U;
    out.omega = omega_new;
    fill_omega_residual(out, in);
    print_debug_iter(A, S_mat, omega_new, U);
    print_rot_slots(out.iter, diag_result);
    lap(t_subspace);

    // ---- 4. rotate the response pieces needed for Theta -------------------
    rs::transform(world_, out.roots, U);
    rs::transform(world_, V0x,        U);
    rs::transform(world_, E0x,        U);
    rs::transform(world_, gamma,      U);

    // ---- 5. assemble Theta from rotated pieces ---------------------------
    std::vector<Storage> theta(M);
    for (int s = 0; s < M; ++s) {
      theta[s] = assemble_theta(world_, V0x[s], E0x[s], gamma[s]);
    }
    lap(t_rotate);  // rotation transforms + Theta assembly

    // ---- 6. BSH apply + residual ------------------------------------------
    // out.roots[s] currently holds the rotated state (pre-BSH). BSH
    // updates each root in place; we capture the raw BSH residual
    // BEFORE KAIN damping (it's what tracks convergence of the
    // underlying dynamics, not the damped step size).
    //
    // Save the rotated pre-BSH state so KAIN's "x_old" is the
    // self-consistent input to the BSH step, not the pre-rotation
    // input. Without this, KAIN's history would mix pre/post rotation
    // and produce garbage iterates.
    std::vector<Storage> x_pre_bsh = out.roots;
    out.last_bsh_residual.assign(M, 0.0);

    // Stage-1 batching: one BSH apply over the flattened M·n_occ bundle
    // (per-root ω ⇒ per-root op set, concatenated). Residual per slot from
    // a single norm2s over the rotated-minus-new difference.
    bool did_batch_bsh = false;
    if constexpr (std::is_same_v<Type, TDA> &&
                  std::is_same_v<Shell, ClosedShell>) {
      if (batched_ && M > 1) {
        const std::size_t n = out.roots[0].x_alpha.size();
        auto Xf  = tda_batch::flatten(out.roots);   // rotated, pre-BSH
        auto Thf = tda_batch::flatten(theta);
        auto NXf = tda_batch::bsh_apply_flat(world_, gs_, Xf, Thf,
                                             omega_new, n);
        auto diff = madness::sub(world_, Xf, NXf);
        auto dn   = madness::norm2s(world_, diff);  // per-function, M·n
        for (int s = 0; s < M; ++s) {
          double acc = 0.0;
          for (std::size_t p = 0; p < n; ++p)
            acc += dn[s * n + p] * dn[s * n + p];
          out.last_bsh_residual[s] = std::sqrt(acc);
        }
        tda_batch::unflatten_into(NXf, n, out.roots);
        for (int s = 0; s < M; ++s)
          print_debug_norms(s, in, V0x[s].x_alpha, Tg(s, s),
                            gamma[s].x_alpha, theta[s].x_alpha);
        did_batch_bsh = true;
      }
    }

    if (!did_batch_bsh) {
      for (int s = 0; s < M; ++s) {
        auto x_new = K::bsh_apply(world_, gs_, out.roots[s],
                                  theta[s], omega_new(s));
        out.last_bsh_residual[s] =
            K::compute_residual_norm(world_, out.roots[s], x_new);
        print_debug_norms(s, in, V0x[s].x_alpha, Tg(s, s),
                          gamma[s].x_alpha, theta[s].x_alpha);
        out.roots[s] = std::move(x_new);
      }
    }

    // ---- 7. KAIN acceleration + step restriction --------------------------
    // Feed (rotated x_pre_bsh) → (raw BSH output in out.roots).
    // policy_.kain=false disables.
    {
      const int kain_diag =
          (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
      if (!kain_held_off(out)) kain_.apply(x_pre_bsh, out.roots, kain_diag);
    }
    lap(t_bsh);  // BSH apply + residual + KAIN

    // Explosion guard — if any per-root BSH residual blows up past
    // policy_.explosion_guard, mark the state diverged and let iterate<>
    // stop. Legacy iterate_excited.cc:164 uses the same heuristic.
    for (double r : out.last_bsh_residual) {
      if (r > policy_.explosion_guard) { out.diverged = true; break; }
    }

    if (time_phases_ && world_.rank() == 0) {
      const double tot = t_gamma + t_build + t_subspace + t_rotate + t_bsh;
      madness::print("ES_TIMING iter=", out.iter,
                     " gamma=", t_gamma, " build=", t_build,
                     " subspace=", t_subspace, " rotate=", t_rotate,
                     " bsh=", t_bsh, " total=", tot, " (s, rank0)");
    }

    update_stall(out, in);

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// "Stream Lambda, drop pieces, recompute V0x/E0x/gamma after
  /// rotation to assemble Theta in place." Memory-conscious variant.
  /// Algorithm:
  ///   0. Top-of-iter Q + GS (shared).
  ///   1+2. Per root, stream Lambda = V0x − E0x_full + gamma (kinetic block
  ///        added to A in gradient form),
  ///        freeing each kernel temporary after it's folded in.
  ///   3. Subspace A = <X|Λ>, S = <X|X>, diagonalize → omega, U.
  ///   4. Rotate X only (drop Lambda).
  ///   5. RECOMPUTE V0x, E0x_nodi, gamma against rotated X, stream
  ///      Theta = V0x − E0x + gamma in place.
  ///   6. BSH apply → new X.
  ///   7. KAIN (shared).
  ///
  /// Peak memory: ~3-4M Storage (vs ~8M for step_rotate_pieces).
  /// Cost: ~2× the per-iter kernel work (each V0x/E0x/gamma computed
  /// twice — once for Lambda, once for Theta).
  State step_recompute_pieces(State in) {
    State out;
    out.iter  = in.iter + 1;
    out.roots = in.roots;
    out.stable_index = in.stable_index;  // root identity travels with the state
    const int M = n_roots_;
    const double thr = madness::FunctionDefaults<3>::get_thresh();

    // ---- 0. top-of-iter discipline ----------------------------------------
    rs::project(out.roots, gs_.Qa, gs_.Qb);
    rs::orthonormalize(world_, out.roots);

    // ---- 1 + 2. per-root Lambda assembly (streamed) -----------------------
    std::vector<Storage> lambda(M);
    std::vector<madness::real_function_3d> rho_alpha(M);
    out.last_density_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      rho_alpha[s] = K::compute_density(world_, gs_, out.roots[s]);
      if (!in.rho_alpha_prev.empty()) {
        auto drho = rho_alpha[s] - in.rho_alpha_prev[s];
        out.last_density_residual[s] = drho.norm2();
      }

      // Stream Lambda = V0x − E0x_full + gamma (kinetic block added to A in
      // gradient form below — rs::kinetic_gram). Each temporary is scoped:
      // built, axpy'd into lambda[s], then freed.
      lambda[s] = K::compute_V0x(world_, gs_, out.roots[s]);
      {
        auto Efull = K::compute_E0x_full(world_, gs_, out.roots[s]);
        lambda[s].axpy(world_, -1.0, Efull);
      }
      {
        auto G = K::compute_gamma(world_, gs_, out.roots[s], rho_alpha[s]);
        lambda[s].axpy(world_, +1.0, G);
      }
      lambda[s].truncate_all(world_, thr);
    }
    out.rho_alpha_prev = std::move(rho_alpha);

    // ---- 3. subspace A, S, diagonalize ------------------------------------
    // See step_rotate_pieces for the State-overloaded rs::* contract:
    // flatten/from_flat happen under the hood, so OpenShell α+β bundle
    // metric and slot-preserving rotation come for free.
    auto A     = rs::inner(out.roots, lambda);
    A += rs::kinetic_gram(world_, out.roots);   // gradient-form kinetic block
    auto S_mat = rs::metric(out.roots, out.roots);
    auto diag_result = rs::diagonalize(A, S_mat,
                                       /*thresh_degenerate=*/-1.0,
                                       policy_.cluster_unmix_factor,
                                       &world_);  // subworld-safe sygvp
    auto &omega_new = diag_result.omega;
    auto &U         = diag_result.U;
    out.omega = omega_new;
    fill_omega_residual(out, in);
    print_debug_iter(A, S_mat, omega_new, U);
    print_rot_slots(out.iter, diag_result);

    // ---- 4. rotate X only; drop Lambda ------------------------------------
    rs::transform(world_, out.roots, U);
    lambda.clear();  // free Lambda — not needed past here

    // ---- 5 + 6. Recompute V0x/E0x/gamma against rotated X, stream Theta
    //              + BSH apply per root --------------------------------------
    // KAIN's "x_old" is the rotated state (pre-BSH), same as
    // step_rotate_pieces. Snapshot before BSH overwrites out.roots.
    std::vector<Storage> x_pre_bsh = out.roots;
    out.last_bsh_residual.assign(M, 0.0);
    for (int s = 0; s < M; ++s) {
      // ρ on rotated X (different from the pre-rotation ρ we used for
      // Lambda — gamma here needs the rotated coupling).
      auto rho_rot = K::compute_density(world_, gs_, out.roots[s]);

      // Stream Theta = V0x − E0x_nodi + gamma, scoped temporaries.
      auto theta = K::compute_V0x(world_, gs_, out.roots[s]);
      {
        auto E = K::compute_E0x(world_, gs_, out.roots[s]);
        theta.axpy(world_, -1.0, E);
      }
      {
        auto G = K::compute_gamma(world_, gs_, out.roots[s], rho_rot);
        theta.axpy(world_, +1.0, G);
      }
      theta.truncate_all(world_, thr);

      // BSH
      auto x_new = K::bsh_apply(world_, gs_, out.roots[s],
                                theta, omega_new(s));
      out.last_bsh_residual[s] =
          K::compute_residual_norm(world_, out.roots[s], x_new);
      out.roots[s] = std::move(x_new);
    }

    // ---- 7. KAIN + guard + banner + log -----------------------------------
    finalize_iter(in, out, x_pre_bsh);
    return out;
  }

  /// Full-deflation locking variant (workstream B; policy.lock_converged).
  /// Converged roots are LOCKED: frozen, removed from the subspace + rotation,
  /// used only to orthogonalize (deflate) the active roots, and skipped by
  /// KAIN/step-restriction (passed a locked mask so slot indices stay stable).
  /// This stops the per-iteration subspace rotation from re-mixing already-
  /// converged roots — the near-degenerate failure mode. Rotate-pieces layout.
  /// Locks are per-protocol (cleared when the wavelet thresh changes). The
  /// deflation overlap uses rs::metric_inner (Euclidean for TDA; the indefinite
  /// RPA metric for Full).
  State step_rotate_pieces_locked(State in) {
    const int M = n_roots_;
    State out;
    out.iter         = in.iter + 1;
    out.roots        = in.roots;
    out.stable_index = in.stable_index;

    // Per-protocol locks: clear when the wavelet thresh changes.
    const double cur_thresh = madness::FunctionDefaults<3>::get_thresh();
    out.locked = in.locked;
    if (cur_thresh != lock_protocol_thresh_) {
      out.locked.assign(M, 0);
      lock_pass_.assign(M, 0);
      lock_protocol_thresh_ = cur_thresh;
    }
    if (static_cast<int>(out.locked.size()) != M) out.locked.assign(M, 0);
    if (static_cast<int>(lock_pass_.size())  != M) lock_pass_.assign(M, 0);

    std::vector<int> act;
    for (int s = 0; s < M; ++s) if (!out.locked[s]) act.push_back(s);
    const int nA = static_cast<int>(act.size());

    // Output bookkeeping sized to M; locked slots carry their prior values.
    out.last_bsh_residual.assign(M, 0.0);
    out.last_density_residual.assign(M, 0.0);
    out.omega = madness::Tensor<double>(M);
    out.rho_alpha_prev = in.rho_alpha_prev;
    if (static_cast<int>(out.rho_alpha_prev.size()) != M)
      out.rho_alpha_prev.resize(M);
    for (int s = 0; s < M; ++s) {
      if (!out.locked[s]) continue;
      if (s < static_cast<int>(in.last_bsh_residual.size()))
        out.last_bsh_residual[s] = in.last_bsh_residual[s];
      if (s < static_cast<int>(in.last_density_residual.size()))
        out.last_density_residual[s] = in.last_density_residual[s];
      if (in.omega.dim(0) == M) out.omega(s) = in.omega(s);
    }

    if (nA == 0) {  // every root converged + locked
      update_stall(out, in);
      print_iter_banner(out);
      append_convergence_log(out);
      return out;
    }

    // ---- top-of-iter: project all; deflate active vs locked; orthonormalize.
    rs::project(out.roots, gs_.Qa, gs_.Qb);
    std::vector<Storage> L_roots;
    for (int s = 0; s < M; ++s) if (out.locked[s]) L_roots.push_back(out.roots[s]);
    std::vector<Storage> A_roots;
    A_roots.reserve(nA);
    for (int i : act) A_roots.push_back(out.roots[i]);
    for (auto &a : A_roots)
      for (auto &L : L_roots) {
        const double ov = rs::metric_inner(L, a);   // deflate: a -= <L|a> L
        a.axpy(world_, -ov, L);
      }
    rs::orthonormalize(world_, A_roots);
    for (int k = 0; k < nA; ++k) out.roots[act[k]] = A_roots[k];

    // ---- per-active-root building blocks ----------------------------------
    std::vector<Storage> V0x(nA), E0x_full(nA), E0x(nA), gamma(nA);
    std::vector<madness::real_function_3d> rho_alpha(nA);
    for (int k = 0; k < nA; ++k) {
      const int s = act[k];
      rho_alpha[k] = K::compute_density(world_, gs_, out.roots[s]);
      gamma[k]     = K::compute_gamma(world_, gs_, out.roots[s], rho_alpha[k]);
      V0x[k]       = K::compute_V0x(world_, gs_, out.roots[s]);
      E0x_full[k]  = K::compute_E0x_full(world_, gs_, out.roots[s]);
      E0x[k]       = K::compute_E0x(world_, gs_, out.roots[s]);
      if (s < static_cast<int>(in.rho_alpha_prev.size()) &&
          in.rho_alpha_prev[s].is_initialized()) {
        auto drho = rho_alpha[k] - in.rho_alpha_prev[s];
        out.last_density_residual[s] = drho.norm2();
      }
    }

    // ---- subspace over the ACTIVE block only ------------------------------
    std::vector<Storage> lambda(nA);   // V0x − E0x_full + γ; kinetic added to A below
    for (int k = 0; k < nA; ++k)
      lambda[k] = assemble_lambda(world_, V0x[k], E0x_full[k], gamma[k]);

    std::vector<Storage> act_roots;
    act_roots.reserve(nA);
    for (int i : act) act_roots.push_back(out.roots[i]);
    auto A     = rs::inner(act_roots, lambda);
    {
      const auto Tg = rs::kinetic_gram(world_, act_roots);   // gradient-form kinetic block
      A += Tg;
      print_debug_pieces(act_roots, Tg, V0x, E0x_full, gamma);
    }
    auto S_mat = rs::metric(act_roots, act_roots);
    auto dr    = rs::diagonalize(A, S_mat, /*thresh_degenerate=*/-1.0,
                                 policy_.cluster_unmix_factor,
                                 &world_);  // subworld-safe sygvp
    auto &omega_act = dr.omega;
    auto &U         = dr.U;
    print_debug_iter(A, S_mat, omega_act, U);   // Debug: subspace A, S, omega, U (active block)
    print_rot_slots(out.iter, dr);

    // ---- rotate active roots + Theta pieces by U --------------------------
    rs::transform(world_, act_roots, U);
    rs::transform(world_, V0x,       U);
    rs::transform(world_, E0x,       U);
    rs::transform(world_, gamma,     U);
    for (int k = 0; k < nA; ++k) out.roots[act[k]] = act_roots[k];

    // ---- Theta + BSH + residual (active) ----------------------------------
    std::vector<Storage> x_pre_bsh = out.roots;  // full M; locked slots frozen
    for (int k = 0; k < nA; ++k) {
      const int s = act[k];
      auto theta = assemble_theta(world_, V0x[k], E0x[k], gamma[k]);
      auto x_new = K::bsh_apply(world_, gs_, out.roots[s], theta, omega_act(k));
      out.last_bsh_residual[s] =
          K::compute_residual_norm(world_, out.roots[s], x_new);
      out.roots[s]          = std::move(x_new);
      out.omega(s)          = omega_act(k);
      out.rho_alpha_prev[s] = std::move(rho_alpha[k]);
    }

    fill_omega_residual(out, in);  // out.omega now fully populated (locked + active)

    // ---- KAIN + step restriction on ACTIVE slots only (masked) ------------
    {
      const int kain_diag = (print_level_ >= PrintLevel::Verbose) ? 1 : 0;
      if (!kain_held_off(out)) kain_.apply(x_pre_bsh, out.roots, kain_diag, &out.locked);
    }

    // ---- explosion guard (active) -----------------------------------------
    for (int i : act)
      if (out.last_bsh_residual[i] > policy_.explosion_guard) {
        out.diverged = true;
        break;
      }

    // ---- lock active roots that pass the (energy+density) criterion for
    //      lock_min_pass CONSECUTIVE iters (debounce against premature locking).
    for (int i : act) {
      if (es_root_converged(out, i)) {
        if (++lock_pass_[i] >= policy_.lock_min_pass) out.locked[i] = 1;
      } else {
        lock_pass_[i] = 0;
      }
    }
    if (print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      int nl = 0;
      for (char c : out.locked) nl += (c ? 1 : 0);
      print("  [lock] active =", nA, " newly+previously locked =", nl, "/", M);
    }

    update_stall(out, in);

    print_iter_banner(out);
    append_convergence_log(out);
    return out;
  }

  /// Fill out.last_omega_residual[s] = |ω_out(s) − ω_in(s)| (eigenvalue step).
  /// Requires out.omega fully populated. Sentinel 1e9 when no prior ω exists
  /// (iter 0 / size mismatch) so a root can't be called converged on it.
  void fill_omega_residual(State &out, const State &in) const {
    const long M = out.omega.dim(0);
    out.last_omega_residual.assign(static_cast<size_t>(M), 1.0e9);
    if (in.omega.dim(0) != M) return;
    for (long s = 0; s < M; ++s)
      out.last_omega_residual[s] = std::abs(out.omega(s) - in.omega(s));
  }

  /// Plateau bookkeeping: one track per gated quantity (density change, |dw|)
  /// over the ACTIVE (unlocked) roots, each normalised by its own target. A
  /// stall needs BOTH flat — testing the max over them stops a solve in which
  /// one quantity is still falling (see ConvergencePolicy::plateau). An
  /// iteration with no measurable |dw| (sentinel) records +inf, which makes
  /// that track inconclusive and withholds the stall.
  void update_stall(State &out, const State &in) const {
    const double thr = madness::FunctionDefaults<3>::get_thresh();
    out.gate_history = (in.gate_thresh == thr)
                           ? in.gate_history
                           : std::vector<std::vector<double>>(2);
    out.gate_thresh  = thr;
    if (out.gate_history.size() < 2) out.gate_history.resize(2);
    const double inf = std::numeric_limits<double>::infinity();
    double g_rho = 0.0, g_om = 0.0;
    bool any = false, om_ok = true;
    const std::size_t M = out.last_density_residual.size();
    for (std::size_t s2 = 0; s2 < M; ++s2) {
      if (s2 < out.locked.size() && out.locked[s2]) continue;
      any = true;
      g_rho = std::max(g_rho, out.last_density_residual[s2] / targets_.density_residual);
      if (s2 >= out.last_omega_residual.size() || out.last_omega_residual[s2] >= 1.0e8)
        om_ok = false;
      else
        g_om = std::max(g_om, out.last_omega_residual[s2] / targets_.omega_residual);
    }
    const bool measured = any && out.iter > 1;
    out.gate_history[0].push_back(measured ? g_rho : inf);
    out.gate_history[1].push_back(measured && om_ok ? g_om : inf);
    out.stalled = !out.diverged && policy_.plateau(out.gate_history);
    if (out.stalled && print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      const std::size_t n = out.gate_history[0].size();
      const std::size_t w = static_cast<std::size_t>(policy_.stall_window);
      print("[STALL] iter", out.iter, ": density and |dw| both flat over",
            policy_.stall_window, "iters (<", 100.0 * policy_.stall_ratio,
            "% improvement) — density", out.gate_history[0][n - 1], "vs",
            out.gate_history[0][n - 1 - w], ", |dw|", out.gate_history[1][n - 1],
            "vs", out.gate_history[1][n - 1 - w], "(normalised); stopping.");
    }
  }

  /// Per-root ES convergence: energy + density (NOT the BSH amplitude, which
  /// jitters for near-degenerate roots after ω/ρ are settled). Needs ≥2 iters
  /// (one to measure Δω). Drives both converged() and the lock trigger so they
  /// share one definition. explosion_guard still catches runaway amplitudes.
  bool es_root_converged(const State &out, int s) const {
    if (out.iter <= std::max(1, policy_.min_iters_before_conv)) return false;
    if (s >= static_cast<int>(out.last_density_residual.size())) return false;
    if (s >= static_cast<int>(out.last_omega_residual.size()))   return false;
    return out.last_density_residual[s] < targets_.density_residual &&
           out.last_omega_residual[s]   < targets_.omega_residual;
  }

  /// Converged: diverged, OR every root passes es_root_converged (energy +
  /// density). With locking on this is equivalent to "all roots locked".
  bool converged(const State &s) const {
    if (s.diverged) return true;  // exit; print_final reports diverged
    if (s.stalled)  return true;  // exit; executor reads State::stalled
    if (s.iter < policy_.min_iters_before_conv) return false;
    const long M = s.omega.dim(0);
    if (M == 0) return false;
    for (long i = 0; i < M; ++i)
      if (!es_root_converged(s, static_cast<int>(i))) return false;
    return true;
  }

  void save(const State &s, const std::string &path_prefix) const {
    for (int b = 0; b < n_roots_; ++b) {
      s.roots[b].save(world_, path_prefix + ".root_" + std::to_string(b));
    }
  }

  /// Permute the state so omega is ascending. The diagonal-dominance
  /// sort inside rs::diagonalize preserves slot identity from the
  /// initial guess (so per-slot densities don't ping-pong across
  /// near-degeneracies), but the initial guess slot order is
  /// arbitrary. Call this once after convergence to get a canonical
  /// ordering for output / archive comparison against legacy.
  void sort_state_by_omega(State &s) const {
    const long M = s.omega.dim(0);
    if (M <= 1) return;

    std::vector<long> perm(M);
    std::iota(perm.begin(), perm.end(), 0L);
    std::sort(perm.begin(), perm.end(),
              [&](long a, long b) { return s.omega(a) < s.omega(b); });

    // No-op fast path.
    bool already_sorted = true;
    for (long i = 0; i < M; ++i) {
      if (perm[i] != i) { already_sorted = false; break; }
    }
    if (already_sorted) return;

    auto omega_new = madness::copy(s.omega);
    std::vector<Storage> roots_new(M);
    std::vector<double>  bsh_new(M, 0.0), drho_new(M, 0.0);
    std::vector<madness::real_function_3d> rho_prev_new(s.rho_alpha_prev.size());
    std::vector<char>    locked_new(s.locked.empty() ? 0 : M, 0);
    for (long i = 0; i < M; ++i) {
      const long src = perm[i];
      omega_new(i)  = s.omega(src);
      roots_new[i]  = std::move(s.roots[src]);
      if (src < static_cast<long>(s.last_bsh_residual.size()))
        bsh_new[i] = s.last_bsh_residual[src];
      if (src < static_cast<long>(s.last_density_residual.size()))
        drho_new[i] = s.last_density_residual[src];
      if (src < static_cast<long>(s.rho_alpha_prev.size()))
        rho_prev_new[i] = s.rho_alpha_prev[src];
      if (!s.locked.empty() && src < static_cast<long>(s.locked.size()))
        locked_new[i] = s.locked[src];
    }
    s.omega                 = omega_new;
    s.roots                 = std::move(roots_new);
    s.last_bsh_residual     = std::move(bsh_new);
    s.last_density_residual = std::move(drho_new);
    s.rho_alpha_prev        = std::move(rho_prev_new);
    if (!s.locked.empty()) s.locked = std::move(locked_new);
    // Keep each root's stable identity attached to its data as slots move.
    permute_stable_index(s.stable_index, perm);

    if (print_level_ >= PrintLevel::Normal && world_.rank() == 0) {
      print("  [sort_state_by_omega] permuted slots to ascending omega:");
      for (long i = 0; i < M; ++i) {
        const int sid = (i < static_cast<long>(s.stable_index.size()))
                            ? s.stable_index[i] : -1;
        print("    new slot", i, "<- old slot", perm[i],
              "  omega =", s.omega(i),
              "  root_id =", (sid >= 0 ? make_root_id(sid) : "(none)"));
      }
    }
  }

  /// Assign initial stable identities (slot i -> stable i) if none exist.
  /// No-op when stable_index is already populated (e.g. loaded from disk),
  /// so a restart keeps the identity it was saved with. Called once before
  /// the protocol ramp via iterate_protocol's SFINAE hook.
  void ensure_root_identity(State &s) const {
    assign_initial_stable_index(s.stable_index,
                                static_cast<int>(s.roots.size()));
  }

private:
  /// Append one row per state at the end of step() to a CSV file at
  /// `log_path_`. Schema:
  ///   iter,protocol_thresh,state,omega,bsh_residual,density_residual,diverged
  /// Header is written exactly once per file (tracked via flag). Only
  /// rank 0 writes; CSV is robust against concurrent reads.
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
    const long M = s.omega.dim(0);
    out.precision(12);
    for (long st = 0; st < M; ++st) {
      const double drho = (static_cast<size_t>(st) < s.last_density_residual.size())
                              ? s.last_density_residual[st] : 0.0;
      const double bsh  = (static_cast<size_t>(st) < s.last_bsh_residual.size())
                              ? s.last_bsh_residual[st] : 0.0;
      out << s.iter << ',' << pthr << ',' << st << ','
          << s.omega(st) << ',' << bsh << ','
          << drho << ',' << (s.diverged ? 1 : 0) << '\n';
    }
  }

  madness::World            &world_;
  ResponseGroundState        gs_;
  int                        n_roots_ = 0;
  ConvergencePolicy          policy_;
  ConvergencePolicy::Targets targets_{};
  // Protocol thresh at which the current locks were earned; locks clear when it
  // changes (a coarse-protocol convergence isn't valid at the finer one).
  double                     lock_protocol_thresh_ = -1.0;
  // Per-slot consecutive-pass counter for the lock debounce (lock_min_pass).
  std::vector<int>           lock_pass_;
  ResponseSubspaceKain<Storage>   kain_;
  PrintLevel                 print_level_ = PrintLevel::Normal;
  std::string                log_path_;
  bool                       log_header_written_ = false;
  // Stage-1 bundle batching toggle (see set_batched). Closed-shell TDA only.
  bool                       batched_ = false;
  // Inc-3c tensor-layer γ toggle (see set_gamma_tensor). Closed-shell TDA only.
  bool                       gamma_tensor_ = false;
  // Per-phase wall-time instrumentation toggle (see set_time_phases).
  bool                       time_phases_ = false;
};

// Build helpers (ResponseGroundState + initial guess) live in
// build_response_ground_state.hpp so this header doesn't depend on
// GroundState or the legacy ESSolverGuess.

} // namespace molresponse_v3

// Call-site shape — identical across (Type, Shell):
//
//   auto problem = build_es_problem_tda<ClosedShell>(world, gs, n_roots);
//   ESSolver<TDA, ClosedShell> solver(world, problem,
//                                     ConvergencePolicy{1e-4});
//   typename ESSolver<TDA, ClosedShell>::State s0 =
//       build_initial_guess(...);
//   auto sf = solvers::iterate(solver, s0, {/*max_iters=*/25});
//   solver.save(sf, save_prefix);
//
// With protocol ladder:
//   sf = solvers::iterate_protocol(
//          solver, s0,
//          /*thresholds=*/{1e-4, 1e-6, 1e-7},
//          /*set_protocol=*/[&](double th){
//             set_response_protocol(world, L, th, override_k);
//          },
//          /*max_iters_per_step=*/25);

#endif // MOLRESPONSE_V3_SOLVERS_ES_SOLVER_HPP
