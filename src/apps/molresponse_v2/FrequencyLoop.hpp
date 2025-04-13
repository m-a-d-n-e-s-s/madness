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

inline void promote_static_to_dynamic(const StaticRestrictedResponse &src,
                                      DynamicRestrictedResponse &dst) {
  assert(src.x_alpha.size() == dst.x_alpha.size());
  assert(dst.y_alpha.size() ==
         dst.x_alpha.size()); // Should already be initialized
                              //
  auto &world = src.x_alpha[0].world();
  dst.x_alpha = copy(world, src.x_alpha);
  dst.y_alpha =
      copy(world, src.x_alpha); // Use the same function as initial guess
}

inline void copy_dynamic_guess(const DynamicRestrictedResponse &src,
                               DynamicRestrictedResponse &dst) {

  assert(src.x_alpha.size() == dst.x_alpha.size());
  assert(src.y_alpha.size() == dst.y_alpha.size());

  auto &world = src.x_alpha[0].world();
  dst.x_alpha.resize(src.x_alpha.size());
  dst.y_alpha.resize(src.y_alpha.size());
}

bool promote_from_previous_frequency(World &world,
                                     const ResponseState &current_state,
                                     size_t freq_index, size_t thresh_index,
                                     size_t num_orbitals, bool is_unrestricted,
                                     ResponseVector &response,
                                     ResponseVector &previous_response) {
  size_t prev_index = freq_index - 1;

  // Clone current state and adjust frequency
  ResponseState prev_state = current_state;
  prev_state.set_frequency_index(prev_index);

  // Allocate previous as static
  previous_response = make_response_vector(num_orbitals,
                                           /*is_static=*/true, is_unrestricted);

  if (!load_response_vector(world, prev_state, prev_index, thresh_index,
                            previous_response)) {
    madness::print("❌ Failed to load previous frequency from disk.");
    return false;
  }

  // Promote from static to dynamic
  response = make_response_vector(num_orbitals,
                                  /*is_static=*/false, is_unrestricted);

  const auto &static_guess =
      std::get<StaticRestrictedResponse>(previous_response);
  auto &dyn_guess = std::get<DynamicRestrictedResponse>(response);
  promote_static_to_dynamic(static_guess, dyn_guess);

  return true;
}

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
  bool at_final_protocol = state.at_final_threshold();
  bool is_unrestricted = !ground_state.isSpinRestricted();
  size_t num_orbitals = ground_state.getNumOrbitals();

  ResponseVector previous_response;
  bool have_previous_freq_response = false;

  if (world.rank() == 0) {
    madness::print("🚀 Starting calculation at protocol:", protocol,
                   "for state:", state_id);
  }

  // Frequency loop
  for (size_t i = state.current_frequency_index; i < state.frequencies.size();
       i++) {
    state.set_frequency_index(i);
    double freq = state.current_frequency();
    auto freq_index = state.current_frequency_index;

    // Skip if already converged at final threshold
    if (metadata.is_converged(state_id, freq)) {
      if (world.rank() == 0)
        madness::print("✅ Frequency already converged:", state.description());
      continue;
    }

    bool is_saved =
        metadata.is_saved(state_id, freq, protocol); // Check if already saved
    //
    bool should_solve = !is_saved || at_final_protocol;

    if (!should_solve) {
      if (world.rank() == 0) {
        print("⚠️  Skipping frequency", freq, "at protocol", protocol,
              "for state:", state_id);
      }
      continue;
    }
    ResponseVector response;
    // If we are at the final protocol, we 1. try the current state, 2. try same
    // freq, lower protocol, 3. use previous frquency, same protocol,4.
    // initialize
    if (world.rank() == 0) {
      madness::print("🔄 Attempting to load response vector for",
                     state.description());
    }

    if (is_saved &&
        load_response_vector(world, state, i, thresh_index, response) &&
        !at_final_protocol) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from current protocol.");
      }
    } else if (thresh_index > 0 &&
               load_response_vector(world, state, i, thresh_index - 1,
                                    response)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response vector from previous protocol.");
      }
    } else if (!is_static) {
      if (have_previous_freq_response) {
        if (world.rank() == 0)
          madness::print("📂 Promoting previous in-memory response as guess.");

        if (std::holds_alternative<StaticRestrictedResponse>(
                previous_response)) {
          const auto &static_guess =
              std::get<StaticRestrictedResponse>(previous_response);
          auto &dyn_guess = std::get<DynamicRestrictedResponse>(response);
          promote_static_to_dynamic(static_guess, dyn_guess);
        } else {
          const auto &dyn_prev =
              std::get<DynamicRestrictedResponse>(previous_response);
          auto &dyn_guess = std::get<DynamicRestrictedResponse>(response);
          copy_dynamic_guess(dyn_prev, dyn_guess);
        }

      } else if (freq_index > 0 &&
                 promote_from_previous_frequency(
                     world, state, freq_index, thresh_index, num_orbitals,
                     is_unrestricted, response, previous_response)) {
        have_previous_freq_response = true;

      } else {
        madness::print("❌ Could not promote from previous frequency.");
        response = initialize_guess_vector(world, ground_state, state);
      }
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
    world.gop.fence();

    // Always save results (even if not fully converged)
    save_response_vector(world, state, response);
    metadata.mark_saved(state_id, freq, protocol);

    if (world.rank() == 0) {
      madness::print("💾 Saved response vector at protocol:", protocol,
                     state.description());
    }

    // Explicitly mark convergence status if at final protocol
    if (state.at_final_threshold() && converged) {
      metadata.mark_final_converged(state_id, true);
      metadata.mark_converged(state_id, freq);
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
