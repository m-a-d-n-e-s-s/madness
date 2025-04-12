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
                      size_t max_iter = 5, double conv_thresh = 1e-6) {
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

inline void computeFrequencyLoop(World &world, const ResponseManager &rm,
                                 ResponseState &state,
                                 const GroundStateData &ground_state,
                                 ResponseMetadata &metadata,
                                 ResponseDebugLogger &logger) {

  const auto &frequencies = state.frequencies;
  const auto state_id = state.perturbationDescription();
  double protocol = state.current_threshold();
  size_t thresh_index = state.current_thresh_index;

  bool is_static = state.is_static();
  bool is_unrestricted = !ground_state.isSpinRestricted();
  size_t num_orbitals = ground_state.getNumOrbitals();

  ResponseVector previous_response;
  bool have_previous_freq_response = false;

  // Early check if already computed and saved at this protocol
  if (metadata.is_saved(state_id, protocol)) {
    if (world.rank() == 0) {
      madness::print("✅ Already saved state at protocol:", protocol, state_id);
    }
    return;
  }

  if (world.rank() == 0) {
    madness::print("🚀 Starting calculation at protocol:", protocol,
                   "for state:", state_id);
  }

  // Frequency loop
  for (; !state.at_final_frequency(); state.advance_frequency()) {
    double freq = state.current_frequency();
    size_t freq_index = state.current_frequency_index;

    if (metadata.is_converged(state_id, freq, protocol)) {
      if (world.rank() == 0) {
        madness::print("✅ Frequency already converged:", state.description());
      }
      continue;
    }

    if (world.rank() == 0) {
      madness::print("\n🔄 Processing:", state.description());
    }

    ResponseVector response =
        make_response_vector(num_orbitals, is_static, is_unrestricted);

    // 1. Attempt to load current frequency/protocol response from disk
    if (load_response_vector(world, state, freq_index, thresh_index,
                             response)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from current protocol.");
      }
    }

    // 2. Try loading from the previous protocol
    else if (thresh_index > 0 &&
             load_response_vector(world, state, freq_index, thresh_index - 1,
                                  response)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from previous protocol.");
      }
    }

    // 3. For dynamic cases: use previous frequency's solution as guess
    else if (!is_static && have_previous_freq_response) {
      if (world.rank() == 0) {
        madness::print(
            "📂 Using previous frequency response as initial guess.");
      }
      response = previous_response;
    }

    // 4. Initialize fresh guess (no available data)
    else {
      if (world.rank() == 0) {
        madness::print(
            "🆕 No prior data found. Initializing new guess vector.");
      }
      response = initialize_guess_vector(world, ground_state, state);
    }

    // Run the solver with logging
    logger.start_state(state);
    bool converged =
        solve_response_vector(world, rm, ground_state, state, response, logger);
    logger.finalize_state();

    if (world.rank() == 0) {
      madness::print("📊 Iteration summary for", state.description());
      logger.print_timing_table(state.description());
      logger.print_values_table(state.description());
    }

    // Always save results (even if not fully converged)
    save_response_vector(world, state, response);
    metadata.mark_saved(state_id, protocol);

    if (world.rank() == 0) {
      madness::print("💾 Saved response vector at protocol:", protocol,
                     state.description());
    }

    // Explicitly mark convergence status if at final protocol
    if (state.at_final_threshold() && converged) {
      metadata.mark_final_converged(state_id, true);
      metadata.mark_converged(state_id, freq, protocol);
      if (world.rank() == 0) {
        madness::print("🎯 Final convergence achieved for",
                       state.description());
      }
    }

    previous_response = response;
    have_previous_freq_response = true;
  }

  // If all frequencies done at this protocol, update status
  state.is_converged =
      state.at_final_threshold() && metadata.final_converged(state_id);
}
