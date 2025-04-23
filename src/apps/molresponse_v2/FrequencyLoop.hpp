#pragma once
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseIO.hpp"
#include "ResponseInitializer.hpp"
#include "ResponseManager.hpp"
#include "ResponseMetaData.hpp"
#include "ResponseSolver.hpp"
#include "ResponseState.hpp"

#include <madness/world/world.h>
#define NOT_IMPLEMENTED_THROW                                                  \
  throw std::runtime_error("This solver is not yet implemented.");

inline bool
solve_response_vector(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, ResponseDebugLogger &logger,
                      size_t max_iter = 10, double conv_thresh = 1e-4) {
  return std::visit(
      [&](auto &vec) -> bool {
        using T = std::decay_t<decltype(vec)>;

        if constexpr (std::is_same_v<T, StaticRestrictedResponse>) {
          return StaticRestrictedSolver::iterate(world, rm, gs, state, response,
                                                 logger, max_iter, conv_thresh);
        } else if constexpr (std::is_same_v<T, DynamicRestrictedResponse>) {
          return DynamicRestrictedSolver::iterate(
              world, rm, gs, state, response, logger, max_iter, conv_thresh);
        } else if constexpr (std::is_same_v<T, StaticUnrestrictedResponse>) {
          return StaticUnrestrictedSolver::iterate(
              world, rm, gs, state, response, logger, max_iter, conv_thresh);
        } else if constexpr (std::is_same_v<T, DynamicUnrestrictedResponse>) {
          return DynamicUnrestrictedSolver::iterate(
              world, rm, gs, state, response, logger, max_iter, conv_thresh);
        } else {
          throw std::runtime_error("Unknown ResponseVector type "
                                   "in solver.");
        }
      },
      response);
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
                                 ResponseState &state,
                                 const GroundStateData &ground_state,
                                 ResponseMetadata &metadata,
                                 ResponseDebugLogger &logger) {

  const auto &frequencies = state.frequencies;
  const auto state_id = state.perturbationDescription();
  double protocol = state.current_threshold();
  size_t thresh_index = state.current_thresh_index;

  bool at_final_protocol = state.at_final_threshold();
  bool is_unrestricted = !ground_state.isSpinRestricted();
  size_t num_orbitals = ground_state.getNumOrbitals();

  ResponseVector previous_response =
      make_response_vector(num_orbitals, state.is_static(), is_unrestricted);
  bool have_previous_freq_response = false;

  if (world.rank() == 0) {
    madness::print("🚀 Starting calculation at protocol:", protocol,
                   "for state:", state_id);
  }

  // Frequency loop
  for (size_t i = state.current_frequency_index; i < state.frequencies.size();
       i++) {
    state.set_frequency_index(i);
    bool is_static = state.is_static();
    double freq = state.current_frequency();
    auto freq_index = state.current_frequency_index;
    bool is_saved =
        metadata.is_saved(state_id, protocol, freq); // Check if already saved
    bool should_solve =
        !is_saved ||
        (at_final_protocol && !metadata.is_converged(state_id, protocol, freq));
    if (!should_solve) {
      if (world.rank() == 0) {
        print("⚠️  Skipping frequency", freq, "at protocol", protocol,
              "for state:", state_id);
      }
      continue;
    }
    world.gop.fence();
    ResponseVector guess =
        make_response_vector(num_orbitals, is_static, is_unrestricted);

    if (world.rank() == 0) {
      madness::print("🔄 Attempting to load response vector for",
                     state.description());
    }
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
    }

    // Run the solver with logging
    logger.start_state(state);

    auto max_iter = rm.params().maxiter();
    auto conv_thresh = rm.params().dconv();

    bool converged = solve_response_vector(
        world, rm, ground_state, state, guess, logger, max_iter, conv_thresh);
    logger.finalize_state();

    if (world.rank() == 0) {
      madness::print("📊 Iteration summary for", state.description());
      logger.print_timing_table(state.description());
      logger.print_values_table(state.description());
    }
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
