#pragma once
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include <filesystem>
#include <iostream>
#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

template <typename StateType>
void save_response_vector(World &world, const StateType &state,
                          const ResponseVector &response) {

  static_assert(std::is_base_of_v<AbstractResponseDescriptor, StateType>,
                "StateType must be derived from AbstractResponseDescriptor");

  auto filename = state.response_filename();
  if (world.rank() == 0) {
    std::cout << "✅ Saving response vector to: " << filename << std::endl;
  }
  auto ar = archive::ParallelOutputArchive(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();
  ar & current_k;

  std::visit(
      [&](const auto &vec) {
        int i = 0;

        for (const auto &f : vec.flat) {
          if (world.rank() == 0) {
            print("Saving response vector ", i++);
          }
          ar & f;
        }
      },
      response);
}

template <typename StateType>
inline bool
load_response_vector(World &world, const int &num_orbitals, StateType &state,
                     ResponseVector &load_response, size_t thresh_index = 0,
                     size_t freq_index = 0) {

  static_assert(std::is_base_of_v<AbstractResponseDescriptor, StateType>);

  auto filename = state.response_filename(thresh_index, freq_index);
  if (world.rank() == 0) {
    std::cout << "📥 Loading response vector from: " << filename << std::endl;
  }
  auto spin_restricted = state.is_spin_restricted();
  // start with whatever the state thinks, but override if it's second‐order
  bool is_static = state.is_static();

  load_response =
      make_response_vector(num_orbitals, is_static, !spin_restricted);

  if (!fs::exists(filename + ".00000")) {
    if (world.rank() == 0) {
      std::cerr << "❌ Response file not found: " << filename << std::endl;
    }
    return false;
  }
  archive::ParallelInputArchive ar(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();
  int loaded_k;

  ar & loaded_k;

  FunctionDefaults<3>::set_k(loaded_k);

  std::visit(
      [&](auto &v) {
        int i = 0;
        for (auto &f : v.flat) {
          if (world.rank() == 0) {
            print("Loading response vector ", i++);
          }

          ar & f;
        }
        v.sync();
      },
      load_response);

  FunctionDefaults<3>::set_k(current_k);

  std::visit(
      [&](auto &v) {
        auto &flat_vec = v.flat;

        if (current_k != loaded_k) {
          if (world.rank() == 0) {
            print("Reconstructing response vector with k=", current_k);
          }
          reconstruct(world, flat_vec);
          for (auto &v : flat_vec) {
            v = project(v, current_k, FunctionDefaults<3>::get_thresh(), true);
          }

          truncate(world, flat_vec, FunctionDefaults<3>::get_thresh());
        } else {
          truncate(world, flat_vec, FunctionDefaults<3>::get_thresh());
        }
        v.sync();
      },
      load_response);

  return true;
}
