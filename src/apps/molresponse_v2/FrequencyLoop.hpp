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
                      size_t max_iter = 15, double conv_thresh = 1e-6) {
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
                                 GlobalMetadataManager &global_metadata,
                                 ResponseDebugLogger &logger) {

  const auto &frequencies = state.frequencies;
  const size_t num_frequencies = frequencies.size();
  const auto state_id = state.perturbationDescription();
  bool all_converged = true;

  double thresh = state.current_threshold();
  size_t freq_index = state.current_frequency_index;
  size_t thresh_index = state.current_thresh_index;
  size_t num_orbitals = ground_state.getNumOrbitals();
  bool is_static = state.is_static();
  bool is_unrestricted = !ground_state.isSpinRestricted();

  ResponseVector previous_response;
  bool have_prev = false;

  if (world.rank() == 0) {
    madness::print("🔄 Frequency loop for", state.description(),
                   "(thresh =", thresh, ")");
  }

  int iteration = 0;
  for (; !state.at_final_frequency(); state.advance_frequency()) {
    double freq = state.current_frequency();
    freq_index = state.current_frequency_index;

    if (global_metadata.is_converged(state_id, freq, thresh)) {
      if (world.rank() == 0) {
        madness::print("✅ Already converged:", state.description());
      }
      continue;
    }

    if (world.rank() == 0) {
      madness::print("🔄 Frequency loop for", state.description(),
                     "(thresh =", thresh, ") - iteration:", iteration);
    }

    ResponseVector response =
        make_response_vector(num_orbitals, is_static, is_unrestricted);

    // ============================
    // 1. Try to load current state
    // ============================
    if (load_response_vector(world, state, freq_index, thresh_index,
                             response)) {
      if (world.rank() == 0) {
        madness::print("📂 Loaded response "
                       "vector from file.");
      }
    }
    // ============================
    // 2. Try previous threshold at same
    // frequency
    // ============================
    else if (thresh_index > 0) {
      size_t prev_thresh_index = thresh_index - 1;
      if (load_response_vector(world, state, freq_index, prev_thresh_index,
                               response)) {
        if (world.rank() == 0) {
          madness::print("📂 Loaded response vector "
                         "from "
                         "previous threshold.");
        }
      } else {
        if (world.rank() == 0) {
          madness::print("❌ Failed to load previous "
                         "threshold response.");
        }
        response =
            make_response_vector(num_orbitals, is_static, is_unrestricted);
      }
    }

    // ============================
    // 3. Use previous frequency result as
    // guess (dynamic only)
    // ============================
    else if (!is_static && freq_index > 0 && have_prev) {
      if (world.rank() == 0) {
        madness::print("📂 Using previous frequency "
                       "result as guess.");
      }
      response = previous_response;
    }

    // ============================
    // 4. No previous data → fresh
    // ============================
    else {
      if (world.rank() == 0) {
        madness::print("❌ No previous data. "
                       "Using fresh guess.");
      }
      // here we have to intialize a guess
      // which is taken as the perturbation
      response = initialize_guess_vector(world, ground_state, state);
    }

    // ============================
    // Solve
    // ============================
    if (world.rank() == 0) {
      madness::print("🔄 Solving for", state.description());
    }

    // TODO: Call actual solver here
    // solve_response_vector(world,
    // ground_state, state, response);
    bool converged =
        solve_response_vector(world, rm, ground_state, state, response, logger);
    logger.finalize_state(state);
    if(world.rank() == 0) {
      logger.pretty_print_summary(state.description());
    }
    iteration++;

    if (converged) {
      // ============================
      // ============================
      save_response_vector(world, state, response);

      global_metadata.mark_converged(state_id, freq, thresh);
      if (world.rank() == 0) {
        madness::print("✅ Saved and marked "
                       "as converged:",
                       state.description());
      }

      previous_response = response;
      have_prev = true;
    }
  }

  state.is_converged = true;
}
