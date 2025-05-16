//
// Created by adrianhurtado on 1/18/22.
//
#include "response_parameters.h"

using json = nlohmann::json;

void madness::from_json(const json& j, ResponseParameters& p) {
  p.from_json(j);
  print("--------JSON-Parameters------------\n", p.print_to_string());
}
bool madness::operator==(const ResponseParameters& p1, const ResponseParameters& p2) { return p1.to_json() == p2.to_json(); }
bool madness::operator!=(const ResponseParameters& p1, const ResponseParameters& p2) { return !(p1 == p2); }
