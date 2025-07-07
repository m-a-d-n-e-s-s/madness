
#pragma once
#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;

inline json &global_inner_contributions() {
  static json g_inner_contributions = json::object();
  return g_inner_contributions;
}
