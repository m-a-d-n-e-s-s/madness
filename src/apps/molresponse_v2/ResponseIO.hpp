#pragma once
#include <filesystem>
#include <iostream>
#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseState.hpp"
#include "ResponseVector.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

template <typename Desc>
void save_response_vector(World &world, const Desc &desc,
                          const ResponseVector &response) {
  auto filename = desc.response_filename();
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

template <typename Desc>
inline bool load_response_vector(World &world, const int &num_orbitals,
                                 const Desc &desc,
                                 ResponseVector &load_response,
                                 size_t thresh_index = 0,
                                 size_t freq_index = 0) {
  auto filename = desc.response_filename(thresh_index, freq_index);
  if (!fs::exists(filename + ".00000")) {
    if (world.rank() == 0) {
      std::cerr << "❌ Response file not found: " << filename << std::endl;
    }
    return false;
  }
  if (world.rank() == 0) {
    std::cout << "Loading response vector from: " << filename << std::endl;
  }

  load_response = desc.make_vector(num_orbitals, freq_index);

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
