#pragma once

// Rank-0 reads a JSON file and broadcasts it to every rank (collective).
// v3-native port of apps/molresponse_v2/broadcast_json.hpp (M1 decoupling
// Stage 1): identical semantics, but namespaced and without the v2 header's
// global-scope `using` declarations.

#include "madness/external/nlohmann_json/json.hpp"
#include <fstream>
#include <madness/world/world.h>
#include <sstream>

namespace molresponse_v3 {

inline nlohmann::json broadcast_json_file(madness::World &world,
                                          const std::string &filepath) {
  std::string json_str;

  if (world.rank() == 0) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
      std::cerr << "ERROR JSON_FILE_OPEN_FAILED file=" << filepath
                << std::endl;
      throw std::runtime_error("Failed to open JSON file");
    }
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    json_str = buffer.str();
    ifs.close();
  }

  world.gop.fence();
  world.gop.broadcast_serializable(json_str, 0);
  world.gop.fence();

  return nlohmann::json::parse(json_str);
}

} // namespace molresponse_v3
