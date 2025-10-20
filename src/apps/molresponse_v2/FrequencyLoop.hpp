#pragma once
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseIO.hpp"
#include "ResponseInitializer.hpp"
#include "ResponseManager.hpp"
#include "ResponseRecord.hpp"
#include "ResponseSolver.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseState.hpp"

#include <madness/world/world.h>
#define NOT_IMPLEMENTED_THROW                                                  \
  throw std::runtime_error("This solver is not yet implemented.");

template <typename ResponseType>
bool iterate(World &world, const ResponseManager &rm, const GroundStateData &gs,
             const LinearResponseDescriptor &state, ResponseType &response,
             ResponseDebugLogger &logger, size_t max_iter, double conv_thresh) {
  // using Policy = ResponseSolverPolicy<ResponseType>;

  auto &rvec = response;
  auto &all_x = rvec.flat;

  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh() * 10, conv_thresh);
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), gs.molecule.natom()));
  const auto x_residual_target = density_target * 10.0;

  auto vp = perturbation_vector(world, gs, state);

  auto &phi0 = gs.orbitals;
  const auto &orbital_energies = gs.getEnergies();

  // First difference, Make bsh operators is different for each solver
  auto bsh_ops =
      make_bsh_operators(world, rm, state.current_frequency(), orbital_energies,
                         static_cast<int>(gs.orbitals.size()), logger, rvec);

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())),
      /*do_printing*/ false);

  auto drho = compute_density(world, rvec, phi0);
  functionT drho_old;

  logger.start_state(state);
  for (size_t iter = 0; iter < max_iter; ++iter) {
    logger.begin_iteration(iter);
    drho_old = copy(drho);
    // Inner product of response state
    DEBUG_LOG_VALUE(world, &logger, "<x|x>",
                    ResponseSolverUtils::inner(world, rvec.flat, rvec.flat));
    // 1. Coupled-response equations
    vector_real_function_3d x_new;
    DEBUG_TIMED_BLOCK(world, &logger, "compute_rsh", {
      x_new =
          CoupledResponseEquations(world, gs, rvec, vp, bsh_ops, rm, logger);
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
    if (res_norm > rm.params().maxrotn() && false) {
      DEBUG_TIMED_BLOCK(world, &logger, "step_restriction", {
        ResponseSolverUtils::do_step_restriction(world, all_x, kain_x, res_norm,
                                                 "a", rm.params().maxrotn());
      });
    }
    // 5. Update response vector
    rvec.flat = copy(world, kain_x);
    rvec.sync();
    // 6. Compute updated response density
    drho = compute_density(world, rvec, phi0);
    double drho_change = (drho - drho_old).norm2();
    DEBUG_LOG_VALUE(world, &logger, "drho_change", drho_change);
    auto alpha =
        alpha_factor(rvec) * ResponseSolverUtils::inner(world, rvec.flat, vp);
    DEBUG_LOG_VALUE(world, &logger, "alpha", alpha);
    // 7. Convergence check
    logger.end_iteration();
    if (world.rank() == 0) {
      ResponseSolverUtils::print_iteration_line(iter, res_norm, drho_change,
                                                alpha, density_target,
                                                x_residual_target);
    }
    if (drho_change < density_target && res_norm < x_residual_target) {
      rvec.sync();
      return true;
    }
  }
  return false;
};

inline bool
solve_response_vector(World &world, const ResponseManager &rm,
                      const GroundStateData &gs,
                      const LinearResponseDescriptor &state,
                      ResponseVector &response_variant, // the std::variant<…>
                      ResponseDebugLogger &logger, size_t max_iter = 10,
                      double conv_thresh = 1e-4) {
  return std::visit(
      overloaded{// Only wire types that are READY:
                 [&](StaticRestrictedResponse &r) {
                   return iterate(world, rm, gs, state, r, logger, max_iter,
                                  conv_thresh);
                 },
                 [&](DynamicRestrictedResponse &r) {
                   return iterate(world, rm, gs, state, r, logger, max_iter,
                                  conv_thresh);
                 },
                 [&](auto &) -> bool {
                   throw std::logic_error(
                       "This response type isn’t implemented yet");
                 }},
      response_variant);
}

// return std::visit(overloaded{[&](StaticRestrictedResponse &r) {
//                                return iterate(world, rm, gs, state, r,
//                                logger, max_iter, conv_thresh);
//                              },
//                              [&](DynamicRestrictedResponse &r) {
//                                return iterate(world, rm, gs, state, r,
//                                logger, max_iter, conv_thresh);
//                              },
//                              [&](StaticUnrestrictedResponse &r) {
//                                throw std::runtime_error("Static unrestricted
//                                response not implemented yet"); return false;
//                                /*return iterate(world, rm, gs, state, r,
//                                logger, max_iter,*/
//                                /*               conv_thresh);*/
//                              },
//                              [&](DynamicUnrestrictedResponse &r) {
//                                throw std::runtime_error("Dynamic unrestricted
//                                response not implemented yet"); return false;
//                                /*iterate(world, rm, gs, state, r, logger,
//                                max_iter,*/
//                                /*                          conv_thresh);*/
//                              }},
//                   response_variant);

inline void promote_response_vector(World &world, const ResponseVector &x_in,
                                    ResponseVector &x_out) {
  if (std::holds_alternative<StaticRestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("🔁 Promoting static restricted → dynamic restricted");
    }
    const auto &prev_resp = std::get<StaticRestrictedResponse>(x_in);
    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<StaticUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {

      madness::print("🔁 Promoting static unrestricted → dynamic unrestricted");
    }
    const auto &prev_resp = std::get<StaticUnrestrictedResponse>(x_in);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_beta = copy(world, prev_resp.x_beta);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<DynamicRestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("📥 Copying dynamic restricted response");
    }
    const auto &prev_resp = std::get<DynamicRestrictedResponse>(x_in);

    DynamicRestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.flatten();
    x_out = current_resp;
  } else if (std::holds_alternative<DynamicUnrestrictedResponse>(x_in)) {
    if (world.rank() == 0) {
      madness::print("📥 Copying dynamic unrestricted response");
    }
    const auto &prev_resp = std::get<DynamicUnrestrictedResponse>(x_in);

    DynamicUnrestrictedResponse current_resp;
    current_resp.x_alpha = copy(world, prev_resp.x_alpha);
    current_resp.x_beta = copy(world, prev_resp.x_beta);
    current_resp.y_alpha = copy(world, prev_resp.y_alpha);
    current_resp.y_beta = copy(world, prev_resp.y_beta);
    current_resp.flatten();
    x_out = current_resp;
  } else {
    throw std::runtime_error(
        "Unknown response variant in promote_response_vector");
  }
}

inline void computeFrequencyLoop(World &world,
                                 const ResponseManager &response_manager,
                                 LinearResponseDescriptor &state_desc,
                                 const GroundStateData &ground_state,
                                 ResponseRecord2 &response_record,
                                 ResponseDebugLogger &logger) {

  bool at_final_protocol = state_desc.at_final_threshold();
  bool is_unrestricted = !ground_state.isSpinRestricted();
  auto num_orbitals = static_cast<int>(ground_state.getNumOrbitals());

  ResponseVector previous_response = make_response_vector(
      num_orbitals, state_desc.is_static(), is_unrestricted);
  bool have_previous_freq_response = false;
  auto pertDesc = state_desc.perturbationDescription();
  auto thresh_index = state_desc.current_thresh_index;
  // if (world.rank() == 0) {
  //   print("Entering frequency loop for state:", pertDesc,
  //         "at protocol:", ResponseRecord2::protocol_key(protocol),
  //         "threshold index:", thresh_index);
  //
  //   print("Frequencies:", state_desc.frequencies);
  // }

  ResponseVector x_0 = make_response_vector(num_orbitals, state_desc.is_static(),
                                           is_unrestricted);
  for (size_t i = state_desc.current_frequency_index;
       i < state_desc.frequencies.size(); i++) {
    state_desc.set_frequency_index(i);
    auto freq_index = state_desc.current_frequency_index;
    bool is_saved =
        response_record.is_saved(state_desc); // Check if already saved
    bool should_solve =
        !is_saved ||
        (at_final_protocol && !response_record.is_converged(state_desc));
    if (!should_solve) {
      // if (world.rank() == 0) {
      //   print("⚠️  Skipping frequency",
      //         ResponseRecord2::freq_key(state_desc.current_frequency()),
      //         "at protocol", ResponseRecord2::protocol_key(protocol),
      //         "for state:", pertDesc);
      // }
      continue;
    }
    world.gop.fence();

    if (is_saved && load_response_vector(world, num_orbitals, state_desc,
                                         thresh_index, freq_index, x_0)) {
      // if (world.rank() == 0) {
      //   madness::print("📂 Loaded response vector from disk.");
      // }
    } else if (thresh_index > 0 &&
               load_response_vector(world, num_orbitals, state_desc,
                                    thresh_index - 1, freq_index, x_0)) {
      // if (world.rank() == 0) {
      //   madness::print("📂 Loaded response vector from previous protocol.");
      // }
    } else if (!state_desc.is_static()) {
      // Now i try to load the previous frequency response
      // if it's in memory, I can use it
      if (have_previous_freq_response) {
        // if (world.rank() == 0) {
        //   madness::print("📂 Promoting previous in-memory response as
        //   guess.");
        // }
      } else {
        load_response_vector(world, num_orbitals, state_desc, thresh_index,
                             freq_index - 1, x_0);
      }
      world.gop.fence();
      if (state_desc.is_static(freq_index - 1)) {
        promote_response_vector(world, x_0, x_0);
      }
    } else {
      // Static case: just initialize
      x_0 = initialize_guess_vector(world, ground_state, state_desc);
      // if (world.rank() == 0) {
      //   madness::print("📂 Initialized guess vector.");
      // }
    }

    // Run the solver with logging
    logger.start_state(state_desc);

    auto max_iter = response_manager.params().maxiter();
    auto conv_thresh = response_manager.params().dconv();

    bool converged =
        solve_response_vector(world, response_manager, ground_state, state_desc,
                              x_0, logger, max_iter, conv_thresh);

    // if (world.rank() == 0) {
    //   // Print final convergence status
    //   print("State:", pertDesc, "Frequency:",
    //         ResponseRecord2::freq_key(state_desc.current_frequency()),
    //         "at protocol:", ResponseRecord2::protocol_key(protocol),
    //         "threshold index:", thresh_index,
    //         converged ? "✅ Converged" : "❌ Not Converged");
    // }
    if (world.rank() == 0) {
      logger.print_timing_table(state_desc);
      logger.print_values_table(state_desc);
    }
    world.gop.fence();
    // save and record the response vector
    save_response_vector(world, state_desc, x_0);
    world.gop.fence();
    response_record.record_status(state_desc, converged);

    previous_response = x_0;
    have_previous_freq_response = true;
  }
  state_desc.set_frequency_index(0);
}
