#include "ResponseSolver.hpp"
#include "../molresponse/response_macrotask.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "functypedefs.h"
#include "projector.h"

vector_real_function_3d StaticRestrictedSolver::ComputeRSH(
    World &world, const GroundStateData &gs, const ResponseVector &vecs,
    const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
    const ResponseManager &rm, ResponseDebugLogger &logger) {

  auto &x = std::get<StaticRestrictedResponse>(vecs).x_alpha;
  auto &all_x = std::get<StaticRestrictedResponse>(vecs).flat;

  auto num_orbitals = gs.orbitals.size();

  std::vector<int> state_index;
  std::vector<int> ii;

  int i = 0;
  int si = 0;
  for (int j = 0; j < num_orbitals; j++) {
    state_index.push_back(si);
    ii.push_back(i++);
  }

  ResponseComputeGroundExchange t0;
  MacroTask g0_task(world, t0);
  ResponseComputeGammaX tresponse;
  MacroTask gx_task(world, tresponse);

  auto c_xc = gs.xcf_.hf_exchange_coefficient();
  vector_real_function_3d k0;
  vector_real_function_3d gx;

  DEBUG_TIMED_BLOCK(world, &logger, "g0_task", {
    k0 = g0_task(ii, state_index, all_x, gs.orbitals, true);
  });
  DEBUG_TIMED_BLOCK(world, &logger, "gx_task", {
    gx = gx_task(ii, state_index, all_x, gs.orbitals, true);
  });

  auto v_local = gs.V_local * x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);

  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  rsh = gs.Qhat(rsh); // project out the ground state density from the response
  return rsh;
}

vector_real_function_3d DynamicRestrictedSolver::ComputeRSH(
    World &world, const GroundStateData &gs, const ResponseVector &vecs,
    const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
    const ResponseManager &rm, ResponseDebugLogger &logger) {

  auto &x = std::get<DynamicRestrictedResponse>(vecs).x_alpha;
  auto &y = std::get<DynamicRestrictedResponse>(vecs).y_alpha;
  auto &all_x = std::get<DynamicRestrictedResponse>(vecs).flat;

  auto num_orbitals = gs.orbitals.size();

  std::vector<int> state_index;
  std::vector<int> ii;

  int i = 0;
  for (int j = 0; j < 2 * num_orbitals; j++) {
    state_index.push_back(0);
    ii.push_back(i++);
  }

  ResponseComputeGroundExchange t0;
  MacroTask g0_task(world, t0);
  ResponseComputeGammaX tresponse;
  MacroTask gx_task(world, tresponse);

  vector_real_function_3d k0;
  vector_real_function_3d gx;

  DEBUG_TIMED_BLOCK(world, &logger, "g0_task", {
    k0 = g0_task(ii, state_index, all_x, gs.orbitals, false);
  });
  DEBUG_TIMED_BLOCK(world, &logger, "gx_task", {
    gx = gx_task(ii, state_index, all_x, gs.orbitals, false);
  });

  auto c_xc = gs.xcf_.hf_exchange_coefficient();
  auto v_local = gs.V_local * all_x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
  auto epsilony = transform(world, y, gs.Hamiltonian_no_diag, true);
  epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  truncate(world, thetax);
  rsh = gs.Qhat(rsh);
  return rsh;
}

// (v0 + g0) * x + g[x] * phi0 + vp

bool StaticRestrictedSolver::iterate(World &world, const ResponseManager &rm,
                                     const GroundStateData &gs,
                                     const ResponseState &state,
                                     ResponseVector &response,
                                     ResponseDebugLogger &logger,
                                     size_t max_iter, double conv_thresh) {

  DEBUG_LOG_VALUE(world, &logger, "orbital norms", norm2s(world, gs.orbitals));

  auto &rvec = std::get<StaticRestrictedResponse>(response);
  auto &x = std::get<StaticRestrictedResponse>(response).x_alpha;
  auto &all_x = std::get<StaticRestrictedResponse>(response).flat;

  auto thresh = FunctionDefaults<3>::get_thresh();
  const double dconv = std::max(FunctionDefaults<3>::get_thresh(), conv_thresh);
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), gs.molecule.natom()));
  const double residual_target = density_target * 10.0;

  auto freq = state.current_frequency();
  const size_t n = x.size();

  DEBUG_LOG_VALUE(world, &logger, "x_norm",
                  ResponseSolverUtils::inner(world, x, x));
  DEBUG_LOG_VALUE(world, &logger, "all_x",
                  ResponseSolverUtils::inner(world, all_x, all_x));
  // Set up RHS (perturbation vector)
  vector_real_function_3d Vp = state.perturbation_vector(world, gs);

  DEBUG_LOG_VALUE(world, &logger, "Vp_norm",
                  ResponseSolverUtils::inner(world, Vp, Vp));
  auto &phi0 = gs.orbitals;
  const auto &orbital_energies = gs.getEnergies();
  // set up the bsh operators
  auto bsh_x = std::vector<poperatorT>(n);
  double x_shifts = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -.05 - (freq + orbital_energies[long(n) - 1]);
  }
  DEBUG_TIMED_BLOCK(world, &logger, "make_bsh_operators", {
    bsh_x = ResponseSolverUtils::make_bsh_operators_response(
        world, x_shifts, freq, orbital_energies, rm.params().lo());
  });

  // Create KAIN solver
  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())), false);
  solver.set_maxsub(static_cast<int>(rm.params().maxsub()));

  // Make initial response density
  auto drho = StaticRestrictedSolver::compute_density(world, x, phi0);
  functionT drho_old;
  // Main iteration
  for (size_t iter = 0; iter < max_iter; ++iter) {
    logger.begin_iteration(iter);
    drho_old = copy(drho);

    DEBUG_LOG_VALUE(world, &logger, "<x|x>",
                    ResponseSolverUtils::inner(world, rvec.flat, rvec.flat));

    vector_real_function_3d x_new;
    DEBUG_TIMED_BLOCK(world, &logger, "compute_rsh", {
      x_new = StaticRestrictedSolver::ComputeRSH(world, gs, rvec, Vp, bsh_x, rm,
                                                 logger);
    });

    // 2. Form residual r = x_new - x
    auto rx = x_new - all_x;

    vector_real_function_3d temp_vec(n);
    // 3. DIIS / KAIN extrapolation

    DEBUG_TIMED_BLOCK(world, &logger, "KAIN step",
                      { temp_vec = solver.update(all_x, rx); });
    // 4. Apply step restriction
    double res_norm;
    DEBUG_TIMED_BLOCK(world, &logger, "step_restriction", {
      res_norm = ResponseSolverUtils::do_step_restriction(
          world, all_x, temp_vec, "a", rm.params().maxrotn());
    });
    DEBUG_LOG_VALUE(world, &logger, "res_norm", res_norm);

    rvec.flat = copy(world, temp_vec);
    rvec.sync();
    // 5. Compute updated response density
    drho = StaticRestrictedSolver::compute_density(world, rvec.flat, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);

    auto alpha = -4.0 * ResponseSolverUtils::inner(world, rvec.flat, Vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);

    // 6. Convergence check
    if (world.rank() == 0) {
      print("Iter", iter, "  Δρ =", drho_change, "  xres =", res_norm,
            " alpha = ", alpha);
    }

    if (drho_change < density_target && res_norm < conv_thresh) {
      if (world.rank() == 0)
        print("✓ Converged in", iter, "iterations.");

      // logger.log_value("converged", true);
      logger.end_iteration();
      logger.finalize_state();
      // 🔄 Sync flat vector back into structured view
      auto &rvec = std::get<StaticRestrictedResponse>(response);
      rvec.sync();

      return true;
    }
    logger.end_iteration();
  }

  if (world.rank() == 0)
    print("⚠️  Reached max iterations without convergence.");

  logger.finalize_state();
  return false;
}

bool DynamicRestrictedSolver::iterate(World &world, const ResponseManager &rm,
                                      const GroundStateData &gs,
                                      const ResponseState &state,
                                      ResponseVector &response,
                                      ResponseDebugLogger &logger,
                                      size_t max_iter, double conv_thresh) {

  auto &rvec = std::get<DynamicRestrictedResponse>(response);
  auto &x = rvec.x_alpha;
  auto &y = rvec.y_alpha;
  auto &all_x = rvec.flat;

  const auto thresh = FunctionDefaults<3>::get_thresh();
  const double dconv = std::max(thresh, conv_thresh);
  const double density_target =
      dconv * static_cast<double>(std::max(size_t(5), gs.molecule.natom()));
  const double freq = state.current_frequency();
  const size_t n = x.size();

  // Set up perturbation
  auto Vp = state.perturbation_vector(world, gs);
  auto vp_copy = copy(world, Vp);
  Vp.insert(Vp.end(), vp_copy.begin(), vp_copy.end());

  const auto &phi0 = gs.orbitals;
  const auto &orbital_energies = gs.getEnergies();

  // Construct BSH operators for ω and -ω
  double x_shift = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    print("*** Applying energy shift to BSH operator");
    x_shift = -.05 - (freq + orbital_energies[long(n) - 1]);
  }

  auto bsh_x = ResponseSolverUtils::make_bsh_operators_response(
      world, x_shift, freq, orbital_energies, rm.params().lo());
  auto bsh_y = ResponseSolverUtils::make_bsh_operators_response(
      world, 0.0, -freq, orbital_energies, rm.params().lo());

  bsh_x.insert(bsh_x.end(), bsh_y.begin(), bsh_y.end());

  // Initial density and solver
  auto drho = DynamicRestrictedSolver::compute_density(world, all_x, phi0);
  functionT drho_old;

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())), false);

  // Main iteration loop
  for (size_t iter = 0; iter < max_iter; ++iter) {
    logger.begin_iteration(iter);
    drho_old = copy(drho);

    DEBUG_LOG_VALUE(world, &logger, "<x|x>",
                    ResponseSolverUtils::inner(world, rvec.flat, rvec.flat));
    vector_real_function_3d x_new;
    DEBUG_TIMED_BLOCK(world, &logger, "compute_rsh", {
      x_new = DynamicRestrictedSolver::ComputeRSH(world, gs, rvec, Vp, bsh_x,
                                                  rm, logger);
    });
    // 2. Form residual r = x_new - x
    auto rx = x_new - all_x;

    vector_real_function_3d temp_vec(2 * n);
    // 3. DIIS / KAIN step

    DEBUG_TIMED_BLOCK(world, &logger, "KAIN step",
                      { temp_vec = solver.update(all_x, rx); });

    // 4. Apply step restriction
    double res_norm;
    DEBUG_TIMED_BLOCK(world, &logger, "step_restriction", {
      res_norm = ResponseSolverUtils::do_step_restriction(
          world, all_x, temp_vec, "a", rm.params().maxrotn());
    });
    DEBUG_LOG_VALUE(world, &logger, "res_norm", res_norm);
    // 5. Update response vector
    rvec.flat = copy(world, temp_vec);
    rvec.sync();
    // 6. Compute updated response density
    drho = DynamicRestrictedSolver::compute_density(world, rvec.flat, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);
    auto alpha = -2.0 * ResponseSolverUtils::inner(world, rvec.flat, Vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);
    // 7. Convergence check
    if (world.rank() == 0) {
      print("Iter", iter, "  Δρ =", drho_change, "  xres =", res_norm,
            " alpha = ", alpha);
    }

    if (drho_change < density_target && res_norm < conv_thresh) {
      if (world.rank() == 0)
        print("✓ Converged in", iter, "iterations.");

      logger.end_iteration();
      logger.finalize_state();
      // 🔄 Sync flat vector back into structured view
      auto &rvec = std::get<DynamicRestrictedResponse>(response);
      rvec.sync();

      return true;
    }
    logger.end_iteration();
  }

  if (world.rank() == 0)
    print("⚠️  Reached max iterations without convergence.");
  logger.finalize_state();
  return false;
}
