#pragma once
#include <madness/world/world.h>

#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseIO.hpp"
#include "ResponseInitializer.hpp"
#include "ResponseManager.hpp"
#include "ResponseMetaData.hpp"
#include "ResponseSolver.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseState.hpp"
#define NOT_IMPLEMENTED_THROW \
  throw std::runtime_error("This solver is not yet implemented.");

template <typename ResponseType>
inline bool iterate(World &world, const ResponseManager &rm,
                    const GroundStateData &gs,
                    const LinearResponseDescriptor &state,
                    ResponseType &response, ResponseDebugLogger &logger,
                    size_t max_iter, double conv_thresh) {
  using Policy = ResponseSolverPolicy<ResponseType>;

  auto &rvec = response;
  auto &all_x = rvec.flat;

  // const auto thresh = FunctionDefaults<3>::get_thresh();
  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh() * 10, conv_thresh);
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), gs.molecule.natom()));
  const auto x_residual_target = density_target * 10.0;

  auto vp = perturbation_vector(world, gs, state);

  auto &phi0 = gs.orbitals;
  const auto &orbital_energies = gs.getEnergies();

  // First difference, Make bsh operators is different for each solver
  auto bsh_ops = Policy::make_bsh_operators(
      world, rm, state.current_frequency(), orbital_energies,
      static_cast<int>(gs.orbitals.size()), logger);

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())),
      /*do_printing*/ false);

  auto drho = Policy::compute_density(world, rvec, phi0);
  functionT drho_old;

  for (size_t iter = 0; iter < max_iter; ++iter) {
    logger.begin_iteration(iter);
    drho_old = copy(drho);
    // Inner product of response state
    DEBUG_LOG_VALUE(world, &logger, "<x|x>",
                    ResponseSolverUtils::inner(world, rvec.flat, rvec.flat));
    // 1. Coupled-response equations
    vector_real_function_3d x_new;
    DEBUG_TIMED_BLOCK(world, &logger, "compute_rsh", {
      x_new = Policy::CoupledResponseEquations(world, gs, rvec, vp, bsh_ops, rm,
                                               logger);
    });
    // 2. Form residual r = x_new - x
    auto residuals = x_new - all_x;
    vector_real_function_3d kain_x;
    DEBUG_TIMED_BLOCK(world, &logger, "KAIN step",
                      { kain_x = solver.update(all_x, residuals); });

    // 3. Compute norm of difference;
    double res_norm = norm2(world, sub(world, all_x, x_new));
    DEBUG_LOG_VALUE(world, &logger, "res_norm", res_norm);
    // 4. Do step restriction
    if (res_norm > rm.params().maxrotn()) {
      DEBUG_TIMED_BLOCK(world, &logger, "step_restriction", {
        ResponseSolverUtils::do_step_restriction(world, all_x, kain_x, res_norm,
                                                 "a", rm.params().maxrotn());
      });
    }
    // 5. Update response vector
    rvec.flat = copy(world, kain_x);
    rvec.sync();
    // 6. Compute updated response density
    drho = Policy::compute_density(world, rvec, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);
    auto alpha =
        Policy::alpha_factor * ResponseSolverUtils::inner(world, rvec.flat, vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);
    // 7. Convergence check
    if (world.rank() == 0) {
      ResponseSolverUtils::print_iteration_line(iter, res_norm, drho_change,
                                                alpha, density_target,
                                                x_residual_target);
    }
    if (drho_change < density_target && res_norm < x_residual_target) {
      if (world.rank() == 0) print("✓ Converged in", iter, "iterations.");
      logger.end_iteration();
      rvec.sync();
      logger.finalize_state();
      return true;
    }
    logger.end_iteration();
  }
  if (world.rank() == 0)
    print("⚠️  Reached max iterations without convergence.");
  logger.end_iteration();

  if (world.rank() == 0) {
    madness::print("📊 Iteration summary for", state.description());
    logger.print_timing_table(state.description());
    logger.print_values_table(state.description());
  }
  logger.finalize_state();
  return false;
};

inline bool solve_response_vector(
    World &world, const ResponseManager &rm, const GroundStateData &gs,
    const LinearResponseDescriptor &state,
    ResponseVector &response_variant,  // the std::variant<…>
    ResponseDebugLogger &logger, size_t max_iter = 10,
    double conv_thresh = 1e-4) {
  return std::visit(
      overloaded{[&](StaticRestrictedResponse &r) {
                   return iterate(world, rm, gs, state, r, logger, max_iter,
                                  conv_thresh);
                 },
                 [&](DynamicRestrictedResponse &r) {
                   return iterate(world, rm, gs, state, r, logger, max_iter,
                                  conv_thresh);
                 },
                 [&](StaticUnrestrictedResponse &r) {
                   throw std::runtime_error(
                       "Static unrestricted response not implemented yet");
                   return false;
                   /*return iterate(world, rm, gs, state, r, logger, max_iter,*/
                   /*               conv_thresh);*/
                 },
                 [&](DynamicUnrestrictedResponse &r) {
                   throw std::runtime_error(
                       "Dynamic unrestricted response not implemented yet");
                   return false;
                   /*iterate(world, rm, gs, state, r, logger, max_iter,*/
                   /*                          conv_thresh);*/
                 }},
      response_variant);
}

inline void promote_response_vector(World &world,
                                    const ResponseVector &promote_from,
                                    ResponseVector &promote_to) {
  if (std::holds_alternative<StaticRestrictedResponse>(promote_from)) {
    if (world.rank() == 0)
      madness::print("🔁 Promoting static restricted → dynamic restricted");
    const auto &prev_resp = std::get<StaticRestrictedResponse>(promote_from);

    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.flatten();
    promote_to = current_resp;

  } else if (std::holds_alternative<StaticUnrestrictedResponse>(promote_from)) {
    if (world.rank() == 0)
      madness::print("🔁 Promoting static unrestricted → dynamic unrestricted");
    const auto &prev_resp = std::get<StaticUnrestrictedResponse>(promote_from);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_beta = copy(world, prev_resp.x_beta);
    current_resp.flatten();
    promote_to = current_resp;

  } else if (std::holds_alternative<DynamicRestrictedResponse>(promote_from)) {
    if (world.rank() == 0)
      madness::print("📥 Copying dynamic restricted response");
    const auto &prev_resp = std::get<DynamicRestrictedResponse>(promote_from);

    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.flatten();
    promote_to = current_resp;

  } else if (std::holds_alternative<DynamicUnrestrictedResponse>(
                 promote_from)) {
    if (world.rank() == 0)
      madness::print("📥 Copying dynamic unrestricted response");
    const auto &prev_resp = std::get<DynamicUnrestrictedResponse>(promote_from);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.y_beta = copy(world, prev_resp.y_beta);
    current_resp.flatten();
    promote_to = current_resp;
  } else {
    throw std::runtime_error(
        "Unknown response variant in promote_response_vector");
  }
}

inline void computeFrequencyLoop(World &world, const ResponseManager &rm,
                                 LinearResponseDescriptor &state,
                                 const GroundStateData &ground_state,
                                 ResponseMetadata &metadata,
                                 ResponseDebugLogger &logger) {
  // const auto &frequencies = state.frequencies;

  auto state_id = state.perturbationDescription();
  double protocol = state.current_threshold();
  size_t thresh_index = state.current_thresh_index;

  bool at_final_protocol = state.at_final_threshold();
  bool is_unrestricted = !ground_state.isSpinRestricted();
  auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());

  ResponseVector previous_response =
      make_response_vector(num_orbitals, state.is_static(), is_unrestricted);
  bool have_previous_freq_response = false;
  auto pertDesc = state.perturbationDescription();

  // Frequency loop
  for (size_t i = state.current_frequency_index; i < state.frequencies.size();
       i++) {
    state.set_frequency_index(i);
    bool is_static = state.is_static();
    double freq = state.current_frequency();
    auto freq_index = state.current_frequency_index;
    bool is_saved =
        metadata.is_saved(pertDesc, protocol, freq);  // Check if already saved
    bool should_solve =
        !is_saved ||
        (at_final_protocol && !metadata.is_converged(pertDesc, protocol, freq));
    if (!should_solve) {
      if (world.rank() == 0) {
        print("⚠️  Skipping frequency", freq, "at protocol", protocol,
              "for state:", pertDesc);
      }
      continue;
    }
    world.gop.fence();
    ResponseVector guess =
        make_response_vector(num_orbitals, is_static, is_unrestricted);

    /*if (world.rank() == 0) {*/
    /*  madness::print("🔄 Attempting to load response vector for",
     * state.description());*/
    /*}*/
    // At this point, I know that I'm either not saved or not at the final
    if (is_saved && load_response_vector(world, num_orbitals, state, guess,
                                         thresh_index, freq_index)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from disk.");
      }
    } else if (thresh_index > 0 &&
               load_response_vector(world, num_orbitals, state, guess,
                                    thresh_index - 1, freq_index)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from previous protocol.");
      }
    } else if (!is_static) {
      // Now i try to load the previous frequency response
      // if it's in memory, I can use it
      if (have_previous_freq_response) {
        if (world.rank() == 0) {
          madness::print("📂 Promoting previous in-memory response as guess.");
        }
      } else {
        load_response_vector(world, num_orbitals, state, previous_response,
                             thresh_index, freq_index - 1);
      }

      world.gop.fence();
      promote_response_vector(world, previous_response, guess);
    } else {
      // Static case: just initialize
      guess = initialize_guess_vector(world, ground_state, state);
      if (world.rank() == 0) {
        madness::print("📂 Initialized guess vector.");
      }
    }

    // Run the solver with logging
    logger.start_state(state);

    auto max_iter = rm.params().maxiter();
    auto conv_thresh = rm.params().dconv();

    bool converged = solve_response_vector(
        world, rm, ground_state, state, guess, logger, max_iter, conv_thresh);
    logger.finalize_state();

    world.gop.fence();

    // Always save results (even if not fully converged)
    save_response_vector(world, state, guess);
    metadata.mark_saved(state_id, protocol, freq);
    metadata.mark_converged(state_id, protocol, freq, converged);

    if (world.rank() == 0) {
      madness::print("💾 Saved response vector at protocol:", protocol,
                     state.description());
    }

    previous_response = guess;
    have_previous_freq_response = true;
  }

  // reset current frequency index
  state.set_frequency_index(0);
  // If all frequencies done at this protocol, update status
  state.is_converged =
      state.at_final_threshold() && metadata.final_converged(state_id);
}
