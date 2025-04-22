#pragma once
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include <filesystem>
#include <iostream>
#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

inline void save_response_vector(World &world, const ResponseState &state,
                                 const ResponseVector &response) {
  auto filename = state.response_filename();
  if (world.rank() == 0) {
    std::cout << "✅ Saving response vector to: " << filename << std::endl;
  }

  auto spin_restricted = state.is_spin_restricted();
  auto is_static = state.is_static();
  auto ar = archive::ParallelOutputArchive(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();

  ar & current_k;

  if (spin_restricted && is_static) {
    for (auto &v : std::get<StaticRestrictedResponse>(response).x_alpha)
      ar & v;
  } else if (spin_restricted && !is_static) {
    for (auto &v : std::get<DynamicRestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<DynamicRestrictedResponse>(response).y_alpha)
      ar & v;
  } else if (!spin_restricted && is_static) {
    for (auto &v : std::get<StaticUnrestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<StaticUnrestrictedResponse>(response).x_beta)
      ar & v;
  } else { // unrestricted && dynamic
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).y_alpha)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).x_beta)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).y_beta)
      ar & v;
  }
}

inline bool load_response_vector(World &world, const int &num_orbitals,
                                 const ResponseState &state,
                                 size_t thresh_index, size_t freq_index,
                                 ResponseVector &load_response) {

  auto filename = state.response_filename(thresh_index, freq_index);
  auto spin_restricted = state.is_spin_restricted();
  auto is_static = std::abs(state.frequencies[freq_index]) < 1e-8;
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

  if (spin_restricted && is_static) {
    auto &response = std::get<StaticRestrictedResponse>(load_response);
    for (auto &v : response.x_alpha) {
      ar & v;
    }
    response.flatten();

    auto norm_load = norm2s(world, response.flat);
    if (world.rank() == 0) {
      print("Loaded response vector norm:", norm_load);
    }
  } else if (spin_restricted && !is_static) {
    auto &response = std::get<DynamicRestrictedResponse>(load_response);
    for (auto &v : response.x_alpha) {
      ar & v;
    }
    for (auto &v : response.y_alpha) {
      ar & v;
    }
    response.flatten();
  } else if (!spin_restricted && is_static) {
    auto &response = std::get<StaticUnrestrictedResponse>(load_response);
    for (auto &v : response.x_alpha) {
      ar & v;
    }
    for (auto &v : response.x_beta) {
      ar & v;
    }
    response.flatten();
  } else { // unrestricted && dynamic
    auto &response = std::get<DynamicUnrestrictedResponse>(load_response);
    for (auto &v : response.x_alpha) {
      ar & v;
    }
    for (auto &v : response.y_alpha) {
      ar & v;
    }
    for (auto &v : response.x_beta) {
      ar & v;
    }
    for (auto &v : response.y_beta) {
      ar & v;
    }
    response.flatten();
  }

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

inline void save_vbc_vector(World &world, const SecondOrderResponseState &state,
                            const ResponseVector &response) {
  auto filename = state.response_filename();
  if (world.rank() == 0) {
    std::cout << "✅ Saving VBC vector to: " << filename << std::endl;
  }

  auto ar = archive::ParallelOutputArchive(world, filename.c_str());
  auto current_k = FunctionDefaults<3>::get_k();
  ar & current_k;

  if (state.spin_restricted && state.is_static()) {
    for (auto &v : std::get<StaticRestrictedResponse>(response).x_alpha)
      ar & v;
  } else if (state.spin_restricted && !state.is_static()) {
    for (auto &v : std::get<DynamicRestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<DynamicRestrictedResponse>(response).y_alpha)
      ar & v;
  } else if (!state.spin_restricted && state.is_static()) {
    for (auto &v : std::get<StaticUnrestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<StaticUnrestrictedResponse>(response).x_beta)
      ar & v;
  } else {
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).x_alpha)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).y_alpha)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).x_beta)
      ar & v;
    for (auto &v : std::get<DynamicUnrestrictedResponse>(response).y_beta)
      ar & v;
  }
}
