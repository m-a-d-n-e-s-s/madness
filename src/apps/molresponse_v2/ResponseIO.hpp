#pragma once
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <madness/external/nlohmann_json/json.hpp>
#include <variant>

using json = nlohmann::json;
namespace fs = std::filesystem;

inline void save_response_vector(World &world, const ResponseState &state,
                                 const ResponseVector &response) {
  auto filename = state.response_filename();
  if (world.rank() == 0) {
    std::cout << "✅ Saving response vector to: " << filename << std::endl;
  }

  auto ar = archive::ParallelOutputArchive(world, filename.c_str());

  // What type of response is it?
  std::visit(
      [&](const auto &vec) {
        if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                     StaticRestrictedResponse>) {
          for (auto &v : vec.x_alpha)
            ar & v;
          for (auto &v : vec.x_alpha)
            ar & v;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            DynamicRestrictedResponse>) {
          for (auto &v : vec.x_alpha)
            ar & v;

          for (auto &v : vec.y_alpha)
            ar & v;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            StaticUnrestrictedResponse>) {
          for (auto &v : vec.x_alpha)
            ar & v;
          for (auto &v : vec.x_alpha)
            ar & v;
          for (auto &v : vec.x_beta)
            ar & v;
          for (auto &v : vec.x_beta)
            ar & v;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            DynamicUnrestrictedResponse>) {
          for (auto &v : vec.x_alpha)
            ar & v;
          for (auto &v : vec.y_alpha)
            ar & v;
          for (auto &v : vec.x_beta)
            ar & v;
          for (auto &v : vec.y_beta)
            ar & v;
        }
      },
      response);
}

inline bool load_response_vector(World &world, const ResponseState &state,
                                 size_t freq_index, size_t thresh_index,
                                 ResponseVector &response) {

  std::string filename = state.response_filename();
  if (!fs::exists(filename)) {
    if (world.rank() == 0) {
      std::cerr << "❌ Response file not found: " << filename << std::endl;
    }
    return false;
  }
  archive::ParallelInputArchive ar(world, state.response_filename().c_str());

  std::visit(
      [&](auto &vec) {
        if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                     StaticRestrictedResponse>) {
          for (auto &v : vec.x_alpha) {
            ar & v;
          }
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            DynamicRestrictedResponse>) {
          for (auto &v : vec.x_alpha) {
            ar & v;
          }

          for (auto &v : vec.y_alpha) {
            ar & v;
          }
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            StaticUnrestrictedResponse>) {
          for (auto &v : vec.x_alpha) {
            ar & v;
          }
          for (auto &v : vec.x_beta) {
            ar & v;
          }
        } else if constexpr (std::is_same_v<std::decay_t<decltype(vec)>,
                                            DynamicUnrestrictedResponse>) {
          for (auto &v : vec.x_alpha) {
            ar & v;
          }
          for (auto &v : vec.y_alpha) {
            ar & v;
          }
          for (auto &v : vec.x_beta) {
            ar & v;
          }
          for (auto &v : vec.y_beta) {
            ar & v;
          }
        }
      },
      response);
  return true;
}
