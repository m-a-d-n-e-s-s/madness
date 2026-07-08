// =========================================================================
// test_v3_fd_skeleton.cpp — validates the FD-side architecture
//
//   kernels/static.hpp              (Kernels<Static, ClosedShell>)
//   kernels/full.hpp                (Kernels<Full,   ClosedShell>)
//   solvers/fd_target.hpp           (FDProblem + FDPerturbation)
//   solvers/fd_solver.hpp           (FDSolver<Type, Shell>)
//   solvers/iterate.hpp             (driver template, shared with ES)
//   solvers/iterate_protocol.hpp    (protocol ladder, shared with ES)
//
// Closed-shell only. Solves one dipole component at one frequency
// (`--axis=z --omega=0.0 --type=static` by default). Reports the
// resulting polarizability component and compares against a reference
// value if one is tabulated for the fixture; otherwise just prints
// "INFO" so first-run-on-a-new-fixture is not a hard failure.
// =========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseKernel.hpp"       // alpha_factor
#include "../ResponseProtocol.hpp"
#include "../kernels/full.hpp"
#include "../kernels/static.hpp"
#include "../kernels/tags.hpp"
#include "../solvers/build_response_ground_state.hpp"  // build_response_ground_state_closed_shell
#include "../solvers/convergence_policy.hpp"
#include "../solvers/fd_save_load.hpp"
#include "../solvers/fd_solver.hpp"
#include "../solvers/fd_problem.hpp"
#include "../solvers/iterate.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/response_state.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

// -------------------------------------------------------------------------
// Reference table — one entry per (fixture, type, omega, axis).
//
// Two references per entry when available:
//
//   alpha_legacy  — MRA legacy fd_solve at the same protocol / basis.
//                   When populated this is the TIGHT pass/fail target.
//   alpha_dalton  — Gaussian HF response (Dalton). MRA vs Gaussian
//                   basis differ slightly so this is INFORMATIONAL,
//                   not pass/fail (~1% agreement is the expectation).
//
// Sourced from
//   /gpfs/projects/rjh/adrian/madness_studies/refs/dalton_polarizability.json
// -------------------------------------------------------------------------
struct FDRef {
  const char* name;
  long n_alpha;
  const char* type;       // "static" or "full"
  double omega;
  int axis;               // 0=x, 1=y, 2=z
  double alpha_legacy;    // NaN = no reference yet
  double alpha_dalton;    // NaN = no reference yet
  const char* dalton_basis;
  const char* note;
};
static const FDRef FD_REFS[] = {
  // All DALTON values from
  //   /gpfs/projects/rjh/adrian/madness_studies/refs/dalton_polarizability.json
  // at aug-cc-pVQZ, HF response. Geometries used by the fixtures live in
  // /gpfs/projects/rjh/adrian/madness/src/apps/molresponse/tests/
  //   fixtures/systems/<mol>_hf/geometry.xyz  (or symlinked equivalents).
  //
  // ---- H2 RHF (k=8 fixture) ----
  // α_xx = α_yy by symmetry; α_zz is the parallel direction.
  {"H2",  1, "static", 0.0000, 0, std::nan(""), 4.615629498931, "aug-cc-pVQZ", ""},
  {"H2",  1, "static", 0.0000, 2, std::nan(""), 6.455562574995, "aug-cc-pVQZ", ""},
  {"H2",  1, "full",   0.0400, 0, std::nan(""), 4.639605330052, "aug-cc-pVQZ", ""},
  {"H2",  1, "full",   0.0400, 2, std::nan(""), 6.495415513427, "aug-cc-pVQZ", ""},
  {"H2",  1, "full",   0.0656, 0, std::nan(""), 4.680754908456, "aug-cc-pVQZ", ""},
  {"H2",  1, "full",   0.0656, 2, std::nan(""), 6.563975177480, "aug-cc-pVQZ", ""},
  // ---- H2O RHF (geometry: O on z, H in y-z plane; C2 axis = z) ----
  {"H2O", 5, "static", 0.0000, 0, std::nan(""), 7.847668958740, "aug-cc-pVQZ", ""},
  {"H2O", 5, "static", 0.0000, 1, std::nan(""), 9.189420520225, "aug-cc-pVQZ", ""},
  {"H2O", 5, "static", 0.0000, 2, std::nan(""), 8.488078839469, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0400, 0, std::nan(""), 7.889312581818, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0400, 1, std::nan(""), 9.223199511206, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0400, 2, std::nan(""), 8.524442446141, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0656, 0, std::nan(""), 7.961721918531, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0656, 1, std::nan(""), 9.281008699488, "aug-cc-pVQZ", ""},
  {"H2O", 5, "full",   0.0656, 2, std::nan(""), 8.587021478808, "aug-cc-pVQZ", ""},
  // ---- Li UHF (atom, spherically symmetric: α_xx = α_yy = α_zz) ----
  // n_alpha = 2 (1s, 2s↑), n_beta = 1 (1s↓). DALTON @ aug-cc-pVQZ.
  // ω=0.0656 deliberately omitted: it's right at the 2s→2p pole
  // (Dalton gives α≈3596 at that frequency — sub-1% basis differences
  // diverge there, so it's not a useful comparison point).
  {"Li",  2, "static", 0.0000, 0, std::nan(""), 170.0017472721, "aug-cc-pVQZ", ""},
  {"Li",  2, "static", 0.0000, 1, std::nan(""), 170.0017472721, "aug-cc-pVQZ", ""},
  {"Li",  2, "static", 0.0000, 2, std::nan(""), 170.0017472721, "aug-cc-pVQZ", ""},
  {"Li",  2, "full",   0.0400, 0, std::nan(""), 262.3065652631, "aug-cc-pVQZ", ""},
  {"Li",  2, "full",   0.0400, 1, std::nan(""), 262.3065652631, "aug-cc-pVQZ", ""},
  {"Li",  2, "full",   0.0400, 2, std::nan(""), 262.3065652631, "aug-cc-pVQZ", ""},
  // ---- C2H4 RHF (geometry: C=C along z, H in y-z plane) ----
  // Strong z-polarization from the π system; x is out-of-plane.
  {"C2H4", 8, "static", 0.0000, 0, std::nan(""), 23.007427109450, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "static", 0.0000, 1, std::nan(""), 24.805237569200, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "static", 0.0000, 2, std::nan(""), 37.026818276860, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0400, 0, std::nan(""), 23.214891500110, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0400, 1, std::nan(""), 24.940308248400, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0400, 2, std::nan(""), 37.497540248800, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0656, 0, std::nan(""), 23.580863199260, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0656, 1, std::nan(""), 25.172962001480, "aug-cc-pVQZ", ""},
  {"C2H4", 8, "full",   0.0656, 2, std::nan(""), 38.335555432940, "aug-cc-pVQZ", ""},
};
static const FDRef* find_ref(long n_alpha, const std::string& type,
                              double omega, int axis) {
  for (const auto& r : FD_REFS) {
    if (r.n_alpha == n_alpha && std::string(r.type) == type &&
        std::abs(r.omega - omega) < 1e-9 && r.axis == axis) {
      return &r;
    }
  }
  return nullptr;
}

static int parse_axis(const std::string& s) {
  if (s.empty()) return 2;
  char c = std::tolower(s[0]);
  return (c == 'x') ? 0 : (c == 'y') ? 1 : 2;
}

// Run the static-closed-shell pathway. Returns final alpha_zz_like.
static double run_static_closed_shell(
    World& world, GroundState& gs, double L, int override_k,
    const std::vector<double>& protocol_thresholds,
    int axis, int max_iters, const ConvergencePolicy& policy, PrintLevel print_level,
    const std::string& fock_json,
    const std::string& log_path = "",
    const std::string& save_dir = "",
    const std::string& load_dir = "") {
  using Solver = FDSolver<Static, ClosedShell>;

  // Build target at the first protocol level (rebuilt inside `prepare`).
  ResponseStateX<ClosedShell> source;
  source.x_alpha = dipole_perturbation(world, gs, axis);

  FDProblem<Static, ClosedShell> tgt;
  tgt.gs = build_response_ground_state_closed_shell(world, gs,
                                             gs.hf_exchange_coefficient(),
                                             gs.params().lo());
  tgt.responses.resize(1);
  tgt.responses[0].omega  = 0.0;
  tgt.responses[0].source = source;

  // Initial guess: use the perturbation itself as starting X (matches
  // legacy fd_solve's pattern).
  Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0].x_alpha = source.x_alpha;

  if (!load_dir.empty()) {
    auto loaded = try_load_fd_state<Static, ClosedShell>(
        world, load_dir, Perturbation::dipole(axis), /*freq=*/0.0);
    if (loaded) {
      // Non-exact loads (coarser saved protocol): the first prepare() in
      // iterate_protocol re-projects to the active (k, thresh).
      s0 = std::move(loaded->state);
    }
  }

  Solver solver(world, tgt, policy, print_level);
  if (!log_path.empty()) solver.set_log_path(log_path);

  auto prepare = [&](double thresh, Solver& solv, Solver::State& st) {
    set_response_protocol(world, L, thresh, override_k);
    const int    new_k = FunctionDefaults<3>::get_k();
    const double new_t = FunctionDefaults<3>::get_thresh();
    if (world.rank() == 0) {
      print("\n--- protocol step: thresh =", new_t,
            "  k =", new_k, " ---");
    }
    auto coulop_new = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
    gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);

    // Rebuild the source (V_pert · phi0) at the new protocol —
    // dipole_perturbation projects against the new phi0 / k.
    // Fresh target — reuse the shared ES builder for the GS side.
    FDProblem<Static, ClosedShell> new_tgt;
    new_tgt.gs = build_response_ground_state_closed_shell(
        world, gs, gs.hf_exchange_coefficient(),
        gs.params().lo());
    new_tgt.responses.resize(1);
    new_tgt.responses[0].omega          = 0.0;
    new_tgt.responses[0].source.x_alpha =
        dipole_perturbation(world, gs, axis);
    solv.set_target(std::move(new_tgt));

    // Reproject in-flight state to the new wavelet basis.
    for (auto& ch : st.responses) {
      for (auto& fn : ch.x_alpha) fn = madness::project(fn, new_k, new_t);
    }
    for (auto& rho : st.rho_alpha_prev)
      rho = madness::project(rho, new_k, new_t);
  };

  solvers::IterateProtocolPolicy proto_policy;
  proto_policy.max_iters_per_step = max_iters;
  // Save a snapshot at each protocol step. The hook fires while the
  // active FunctionDefaults<3> still reflects the just-completed protocol,
  // so save_fd_state picks up the right protocol_key.
  auto post_step = [&](double, Solver& /*solv*/, Solver::State& st) {
    if (!save_dir.empty()) {
      save_fd_state<Static, ClosedShell>(world, st, save_dir,
                                          Perturbation::dipole(axis),
                                          /*freq=*/0.0,
                                          /*converged=*/!st.diverged);
    }
  };
  auto sf = solvers::iterate_protocol(solver, s0, protocol_thresholds,
                                       prepare, post_step, proto_policy);

  // alpha_axis,axis = af · <V_pert · phi0 | x_alpha>
  // af = -4 for restricted static.
  const double af = alpha_factor(ResponseType::Static, true);
  auto src_final = dipole_perturbation(world, gs, axis);
  const double ip = inner(src_final, sf.responses[0].x_alpha);
  return af * ip;
}

// Run the full-closed-shell pathway. Returns final alpha_zz_like.
static double run_full_closed_shell(
    World& world, GroundState& gs, double L, int override_k,
    const std::vector<double>& protocol_thresholds, double omega,
    int axis, int max_iters, const ConvergencePolicy& policy, PrintLevel print_level,
    const std::string& fock_json,
    const std::string& log_path = "",
    const std::string& save_dir = "",
    const std::string& load_dir = "") {
  using Solver = FDSolver<Full, ClosedShell>;

  ResponseStateXY<ClosedShell> source;
  source.x_alpha = dipole_perturbation(world, gs, axis);
  source.y_alpha = dipole_perturbation(world, gs, axis);

  FDProblem<Full, ClosedShell> tgt;
  tgt.gs = build_response_ground_state_closed_shell(world, gs,
                                             gs.hf_exchange_coefficient(),
                                             gs.params().lo());
  tgt.responses.resize(1);
  tgt.responses[0].omega  = omega;
  tgt.responses[0].source = source;

  Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0].x_alpha = source.x_alpha;
  s0.responses[0].y_alpha = source.y_alpha;

  if (!load_dir.empty()) {
    auto loaded = try_load_fd_state<Full, ClosedShell>(
        world, load_dir, Perturbation::dipole(axis), /*freq=*/omega);
    if (loaded) s0 = std::move(loaded->state);
  }

  Solver solver(world, tgt, policy, print_level);
  if (!log_path.empty()) solver.set_log_path(log_path);

  auto prepare = [&](double thresh, Solver& solv, Solver::State& st) {
    set_response_protocol(world, L, thresh, override_k);
    const int    new_k = FunctionDefaults<3>::get_k();
    const double new_t = FunctionDefaults<3>::get_thresh();
    if (world.rank() == 0) {
      print("\n--- protocol step: thresh =", new_t,
            "  k =", new_k, " ---");
    }
    auto coulop_new = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
    gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);

    auto src_new = dipole_perturbation(world, gs, axis);

    FDProblem<Full, ClosedShell> new_tgt;
    new_tgt.gs = build_response_ground_state_closed_shell(
        world, gs, gs.hf_exchange_coefficient(),
        gs.params().lo());
    new_tgt.responses.resize(1);
    new_tgt.responses[0].omega           = omega;
    new_tgt.responses[0].source.x_alpha  = src_new;
    new_tgt.responses[0].source.y_alpha  = src_new;
    solv.set_target(std::move(new_tgt));

    for (auto& ch : st.responses) {
      for (auto& fn : ch.x_alpha) fn = madness::project(fn, new_k, new_t);
      for (auto& fn : ch.y_alpha) fn = madness::project(fn, new_k, new_t);
    }
    for (auto& rho : st.rho_alpha_prev)
      rho = madness::project(rho, new_k, new_t);
  };

  solvers::IterateProtocolPolicy proto_policy;
  proto_policy.max_iters_per_step = max_iters;
  auto post_step = [&](double, Solver& /*solv*/, Solver::State& st) {
    if (!save_dir.empty()) {
      save_fd_state<Full, ClosedShell>(world, st, save_dir,
                                        Perturbation::dipole(axis),
                                        /*freq=*/omega,
                                        /*converged=*/!st.diverged);
    }
  };
  auto sf = solvers::iterate_protocol(solver, s0, protocol_thresholds,
                                       prepare, post_step, proto_policy);

  // alpha_axis,axis = af · (<src | x_alpha> + <src | y_alpha>)
  // af = -2 for restricted Full.
  const double af = alpha_factor(ResponseType::Full, true);
  auto src_final = dipole_perturbation(world, gs, axis);
  const double ip_x = inner(src_final, sf.responses[0].x_alpha);
  const double ip_y = inner(src_final, sf.responses[0].y_alpha);
  return af * (ip_x + ip_y);
}

// =========================================================================
// OpenShell pathways (Li UHF and friends)
// =========================================================================

static double run_static_open_shell(
    World& world, GroundState& gs, double L, int override_k,
    const std::vector<double>& protocol_thresholds,
    int axis, int max_iters, const ConvergencePolicy& policy, PrintLevel print_level,
    const std::string& fock_json,
    const std::string& log_path = "",
    const std::string& save_dir = "",
    const std::string& load_dir = "") {
  using Solver = FDSolver<Static, OpenShell>;

  auto src_alpha = dipole_perturbation(world, gs, axis);
  auto src_beta  = dipole_perturbation_beta(world, gs, axis);

  ResponseStateX<OpenShell> source;
  source.x_alpha = src_alpha;
  source.x_beta  = src_beta;

  FDProblem<Static, OpenShell> tgt;
  tgt.gs = build_response_ground_state_open_shell(world, gs,
                                          gs.hf_exchange_coefficient(),
                                          gs.params().lo());
  tgt.responses.resize(1);
  tgt.responses[0].omega  = 0.0;
  tgt.responses[0].source = source;

  Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0].x_alpha = src_alpha;
  s0.responses[0].x_beta  = src_beta;

  if (!load_dir.empty()) {
    auto loaded = try_load_fd_state<Static, OpenShell>(
        world, load_dir, Perturbation::dipole(axis), /*freq=*/0.0);
    if (loaded) s0 = std::move(loaded->state);
  }

  Solver solver(world, tgt, policy, print_level);
  if (!log_path.empty()) solver.set_log_path(log_path);

  auto prepare = [&](double thresh, Solver& solv, Solver::State& st) {
    set_response_protocol(world, L, thresh, override_k);
    const int    new_k = FunctionDefaults<3>::get_k();
    const double new_t = FunctionDefaults<3>::get_thresh();
    if (world.rank() == 0) {
      print("\n--- protocol step: thresh =", new_t,
            "  k =", new_k, " ---");
    }
    auto coulop_new = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
    gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);

    FDProblem<Static, OpenShell> new_tgt;
    new_tgt.gs = build_response_ground_state_open_shell(
        world, gs, gs.hf_exchange_coefficient(),
        gs.params().lo());
    new_tgt.responses.resize(1);
    new_tgt.responses[0].omega = 0.0;
    new_tgt.responses[0].source.x_alpha = dipole_perturbation(world, gs, axis);
    new_tgt.responses[0].source.x_beta  =
        dipole_perturbation_beta(world, gs, axis);
    solv.set_target(std::move(new_tgt));

    for (auto& ch : st.responses) {
      for (auto& fn : ch.x_alpha) fn = madness::project(fn, new_k, new_t);
      for (auto& fn : ch.x_beta)  fn = madness::project(fn, new_k, new_t);
    }
    for (auto& rho : st.rho_alpha_prev)
      rho = madness::project(rho, new_k, new_t);
  };

  solvers::IterateProtocolPolicy proto_policy;
  proto_policy.max_iters_per_step = max_iters;
  auto post_step = [&](double, Solver& /*solv*/, Solver::State& st) {
    if (!save_dir.empty()) {
      save_fd_state<Static, OpenShell>(world, st, save_dir,
                                        Perturbation::dipole(axis),
                                        /*freq=*/0.0,
                                        /*converged=*/!st.diverged);
    }
  };
  auto sf = solvers::iterate_protocol(solver, s0, protocol_thresholds,
                                       prepare, post_step, proto_policy);

  // alpha = af · (<src_α | x_α> + <src_β | x_β>)
  // af = -2 for unrestricted Static.
  const double af = alpha_factor(ResponseType::Static, false);
  auto sa = dipole_perturbation(world, gs, axis);
  auto sb = dipole_perturbation_beta(world, gs, axis);
  const double ip_a = inner(sa, sf.responses[0].x_alpha);
  const double ip_b = inner(sb, sf.responses[0].x_beta);
  return af * (ip_a + ip_b);
}

static double run_full_open_shell(
    World& world, GroundState& gs, double L, int override_k,
    const std::vector<double>& protocol_thresholds, double omega,
    int axis, int max_iters, const ConvergencePolicy& policy, PrintLevel print_level,
    const std::string& fock_json,
    const std::string& log_path = "",
    const std::string& save_dir = "",
    const std::string& load_dir = "") {
  using Solver = FDSolver<Full, OpenShell>;

  auto src_alpha = dipole_perturbation(world, gs, axis);
  auto src_beta  = dipole_perturbation_beta(world, gs, axis);

  ResponseStateXY<OpenShell> source;
  source.x_alpha = src_alpha;
  source.y_alpha = src_alpha;
  source.x_beta  = src_beta;
  source.y_beta  = src_beta;

  FDProblem<Full, OpenShell> tgt;
  tgt.gs = build_response_ground_state_open_shell(world, gs,
                                          gs.hf_exchange_coefficient(),
                                          gs.params().lo());
  tgt.responses.resize(1);
  tgt.responses[0].omega  = omega;
  tgt.responses[0].source = source;

  Solver::State s0;
  s0.responses.resize(1);
  s0.responses[0].x_alpha = src_alpha;
  s0.responses[0].y_alpha = src_alpha;
  s0.responses[0].x_beta  = src_beta;
  s0.responses[0].y_beta  = src_beta;

  if (!load_dir.empty()) {
    auto loaded = try_load_fd_state<Full, OpenShell>(
        world, load_dir, Perturbation::dipole(axis), /*freq=*/omega);
    if (loaded) s0 = std::move(loaded->state);
  }

  Solver solver(world, tgt, policy, print_level);
  if (!log_path.empty()) solver.set_log_path(log_path);

  auto prepare = [&](double thresh, Solver& solv, Solver::State& st) {
    set_response_protocol(world, L, thresh, override_k);
    const int    new_k = FunctionDefaults<3>::get_k();
    const double new_t = FunctionDefaults<3>::get_thresh();
    if (world.rank() == 0) {
      print("\n--- protocol step: thresh =", new_t,
            "  k =", new_k, " ---");
    }
    auto coulop_new = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
    gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);

    auto sa = dipole_perturbation(world, gs, axis);
    auto sb = dipole_perturbation_beta(world, gs, axis);

    FDProblem<Full, OpenShell> new_tgt;
    new_tgt.gs = build_response_ground_state_open_shell(
        world, gs, gs.hf_exchange_coefficient(),
        gs.params().lo());
    new_tgt.responses.resize(1);
    new_tgt.responses[0].omega           = omega;
    new_tgt.responses[0].source.x_alpha  = sa;
    new_tgt.responses[0].source.y_alpha  = sa;
    new_tgt.responses[0].source.x_beta   = sb;
    new_tgt.responses[0].source.y_beta   = sb;
    solv.set_target(std::move(new_tgt));

    for (auto& ch : st.responses) {
      for (auto& fn : ch.x_alpha) fn = madness::project(fn, new_k, new_t);
      for (auto& fn : ch.y_alpha) fn = madness::project(fn, new_k, new_t);
      for (auto& fn : ch.x_beta)  fn = madness::project(fn, new_k, new_t);
      for (auto& fn : ch.y_beta)  fn = madness::project(fn, new_k, new_t);
    }
    for (auto& rho : st.rho_alpha_prev)
      rho = madness::project(rho, new_k, new_t);
  };

  solvers::IterateProtocolPolicy proto_policy;
  proto_policy.max_iters_per_step = max_iters;
  auto post_step = [&](double, Solver& /*solv*/, Solver::State& st) {
    if (!save_dir.empty()) {
      save_fd_state<Full, OpenShell>(world, st, save_dir,
                                      Perturbation::dipole(axis),
                                      /*freq=*/omega,
                                      /*converged=*/!st.diverged);
    }
  };
  auto sf = solvers::iterate_protocol(solver, s0, protocol_thresholds,
                                       prepare, post_step, proto_policy);

  // alpha = af · ( <src_α | x_α> + <src_α | y_α> + <src_β | x_β> + <src_β | y_β> )
  // af = -1 for unrestricted Full.
  const double af = alpha_factor(ResponseType::Full, false);
  auto sa = dipole_perturbation(world, gs, axis);
  auto sb = dipole_perturbation_beta(world, gs, axis);
  const double ip_xa = inner(sa, sf.responses[0].x_alpha);
  const double ip_ya = inner(sa, sf.responses[0].y_alpha);
  const double ip_xb = inner(sb, sf.responses[0].x_beta);
  const double ip_yb = inner(sb, sf.responses[0].y_beta);
  return af * (ip_xa + ip_ya + ip_xb + ip_yb);
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    if (world.rank() == 0) {
      print("\n=========================================================");
      print(" molresponse_v3 FD Skeleton Test  (Static/Full × ClosedShell)");
      print("=========================================================\n");
    }
    // Inner scope so MADNESS objects destruct before finalize().
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0) {
          print("Usage: test_v3_fd_skeleton --archive=<path> "
                "[--type=static|full] [--omega=X] [--axis=x|y|z] "
                "[--maxiter=N] [--dconv=X] "
                "[--thresh=X | --protocol=th1,th2,...] [--k=N] "
                "[--tol=X] [--print-level=0..3] "
                "[--log-convergence=PATH] [--save=DIR] [--load=DIR] "
                "[--no-kain] [--kain-maxsub=N] [--maxrotn=X] "
                "[--kain-cmax=X]");
        }
        finalize();
        return 1;
      }

      const std::string archive_path = parser.value_raw("archive");
      const std::string type_str = parser.key_exists("type")
                                       ? parser.value("type") : "static";
      const double omega = parser.key_exists("omega")
                               ? std::stod(parser.value("omega")) : 0.0;
      const int axis = parse_axis(parser.key_exists("axis")
                                      ? parser.value("axis") : "z");
      const int max_iters = parser.key_exists("maxiter")
                                ? std::stoi(parser.value("maxiter")) : 25;
      const double dconv = parser.key_exists("dconv")
                               ? std::stod(parser.value("dconv")) : 1e-4;
      const double tol = parser.key_exists("tol")
                             ? std::stod(parser.value("tol")) : 5e-3;
      const int pl_int = parser.key_exists("print-level")
                             ? std::stoi(parser.value("print-level")) : 1;
      // KAIN acceleration knobs — same flag set as the ES driver so
      // FD can be A/B'd kain-on vs kain-off across (Type × Shell).
      const bool   no_kain     = parser.key_exists("no-kain");
      const int    kain_maxsub = parser.key_exists("kain-maxsub")
                                     ? std::stoi(parser.value("kain-maxsub"))
                                     : 10;
      const double maxrotn     = parser.key_exists("maxrotn")
                                     ? std::stod(parser.value("maxrotn"))
                                     : 0.5;
      const double kain_cmax   = parser.key_exists("kain-cmax")
                                     ? std::stod(parser.value("kain-cmax"))
                                     : 100.0;
      const PrintLevel print_level = static_cast<PrintLevel>(
          std::max(0, std::min(3, pl_int)));
      const std::string log_path =
          parser.key_exists("log-convergence")
              ? parser.value_raw("log-convergence")
              : std::string();
      const std::string save_dir =
          parser.key_exists("save") ? parser.value_raw("save")
                                    : std::string();
      const std::string load_dir =
          parser.key_exists("load") ? parser.value_raw("load")
                                    : std::string();

      auto header = GroundState::read_archive_header(world, archive_path);
      const int override_k = parser.key_exists("k")
                                 ? std::stoi(parser.value("k")) : -1;

      std::vector<double> protocol_thresholds;
      if (parser.key_exists("protocol")) {
        protocol_thresholds = parse_protocol_csv(parser.value("protocol"));
      } else {
        double single = parser.key_exists("thresh")
                            ? std::stod(parser.value("thresh"))
                            : default_thresh_for_k(header.k);
        protocol_thresholds.push_back(single);
      }
      set_response_protocol(world, header.L, protocol_thresholds.front(),
                            override_k);

      Molecule molecule;
      auto archive_dir = std::filesystem::path(archive_path).parent_path();
      for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
        auto candidate = archive_dir / name;
        if (std::filesystem::exists(candidate)) {
          std::ifstream ifs(candidate);
          nlohmann::json j; ifs >> j;
          nlohmann::json mol_json;
          if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
            mol_json = j["tasks"][0]["molecule"];
          else if (j.contains("molecule"))
            mol_json = j["molecule"];
          if (!mol_json.is_null()) molecule.from_json(mol_json);
          break;
        }
      }
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      if (world.rank() == 0) gs.print_info();

      std::string fock_json;
      for (const auto& name : {"moldft.fock.json", "mad.fock.json"}) {
        auto candidate = archive_dir / name;
        if (std::filesystem::exists(candidate)) {
          fock_json = candidate.string();
          break;
        }
      }
      const double cur_thresh = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(
          CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
      gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);

      const bool restricted = gs.is_spin_restricted();
      const char* axis_label = (axis == 0) ? "x" : (axis == 1) ? "y" : "z";
      if (world.rank() == 0) {
        print("\nRUN CONFIG:");
        print("  type       =", type_str);
        print("  shell      =", restricted ? "ClosedShell" : "OpenShell");
        print("  omega      =", omega);
        print("  axis       =", axis_label);
        print("  dconv_user =", dconv);
        print("  kain       =", no_kain ? "off" : "on",
              "  kain_maxsub =", kain_maxsub,
              "  maxrotn =", maxrotn,
              "  kain_cmax =", kain_cmax);
        print("  max_iters  =", max_iters, " per protocol step");
        print("  protocol   =", protocol_thresholds);
      }

      // Assemble the full policy once; threaded into every dispatch.
      ConvergencePolicy policy;
      policy.dconv_user  = dconv;
      policy.kain        = !no_kain;
      policy.kain_maxsub = kain_maxsub;
      policy.maxrotn     = maxrotn;
      policy.kain_cmax_cap = kain_cmax;

      double alpha = 0.0;
      if (type_str == "static" && restricted) {
        alpha = run_static_closed_shell(world, gs, header.L, override_k,
                                         protocol_thresholds, axis,
                                         max_iters, policy, print_level,
                                         fock_json, log_path, save_dir, load_dir);
      } else if (type_str == "full" && restricted) {
        alpha = run_full_closed_shell(world, gs, header.L, override_k,
                                       protocol_thresholds, omega, axis,
                                       max_iters, policy, print_level,
                                       fock_json, log_path, save_dir, load_dir);
      } else if (type_str == "static" && !restricted) {
        alpha = run_static_open_shell(world, gs, header.L, override_k,
                                       protocol_thresholds, axis,
                                       max_iters, policy, print_level,
                                       fock_json, log_path, save_dir, load_dir);
      } else if (type_str == "full" && !restricted) {
        alpha = run_full_open_shell(world, gs, header.L, override_k,
                                     protocol_thresholds, omega, axis,
                                     max_iters, policy, print_level,
                                     fock_json, log_path);
      } else {
        if (world.rank() == 0) print("Unknown --type:", type_str);
        finalize();
        return 1;
      }

      if (world.rank() == 0) {
        print("\n=== Final result ===");
        print("  alpha_", axis_label, axis_label, " =", alpha,
              "  (type =", type_str, ", omega =", omega, ")");
      }

      // Reference checks.
      //   alpha_legacy (MRA-on-MRA): tight PASS/FAIL when populated.
      //   alpha_dalton (Gaussian):   informational INFO with ~1% expected
      //                              agreement window (MRA vs Gaussian
      //                              basis-set difference — NOT pass/fail).
      const FDRef* ref = find_ref(gs.num_alpha(), type_str, omega, axis);
      if (!ref && world.rank() == 0) {
        print("\nNo reference tabulated for n_alpha=", gs.num_alpha(),
              " type=", type_str, " omega=", omega,
              " axis=", axis_label,
              " — running without comparison.");
      }
      if (ref) {
        if (!std::isnan(ref->alpha_legacy)) {
          const double err = std::abs(alpha - ref->alpha_legacy);
          const bool ok    = err < tol;
          if (world.rank() == 0) {
            print("\n=== PASS/FAIL  vs legacy fd_solve (", ref->name, ") ===");
            print("  [", (ok ? "PASS" : "FAIL"), "]  alpha  ref=",
                  ref->alpha_legacy, "  actual=", alpha,
                  "  err=", err, "  tol=", tol);
          }
          rc = ok ? 0 : 1;
        }
        if (!std::isnan(ref->alpha_dalton)) {
          const double err     = std::abs(alpha - ref->alpha_dalton);
          const double pct_err = 100.0 * err / std::abs(ref->alpha_dalton);
          if (world.rank() == 0) {
            print("\n=== INFO  vs DALTON HF (", ref->name, ", basis ",
                  ref->dalton_basis,
                  ") — informational, not pass/fail ===");
            print("  alpha  ref=", ref->alpha_dalton,
                  "  actual=", alpha,
                  "  err=", err, " (", pct_err, "%)");
            // Basis-set difference between MRA (k=8, thresh=1e-6) and
            // Gaussian aug-cc-pVQZ for HF response is typically <1%
            // for closed-shell first-row molecules. >5% strongly hints
            // at an implementation bug in v3.
            if (pct_err > 5.0) {
              print("  WARNING: >5% disagreement with Dalton — likely a "
                    "v3 implementation bug, not a basis-set effect.");
            }
          }
        }
      }
      world.gop.fence();
    }  // end inner scope
    world.gop.fence();
    finalize();
    return rc;
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("ERROR:", e.what());
    finalize();
    return 2;
  }
}
