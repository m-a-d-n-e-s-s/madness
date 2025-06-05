#pragma once

#include "madness/external/nlohmann_json/json.hpp"
#include <fstream>
#include <madness/world/world.h>
#include <sstream>

using json = nlohmann::json;
using madness::World;

inline json broadcast_json_file(World &world, const std::string &filepath) {
  std::string json_str;

  if (world.rank() == 0) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
      std::cerr << "❌ Error opening file: " << filepath << std::endl;
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

  return json::parse(json_str);
}
