//
// Created by adrianhurtado on 1/18/22.
//
#include "response_parameters.h"
using json = nlohmann::json;

void ResponseParameters::to_json(json& j) {
  json j_params={};

  // TODO Is there a way to the get member for every parameter even though get is a template function?
  for (auto& p : parameters) {
    auto param_type = p.second.get_type();
    if (param_type == "i") {
      j_params[p.first] = get<int>(p.first);
      // if vector of double
    } else if (param_type == "d") {
      j_params[p.first] = get<double>(p.first);

      // if vector of bool
    } else if (param_type == "b") {
      j_params[p.first] = get<bool>(p.first);

      // if vector of doubles?
    } else if (param_type == "St6vectorIdSaIdEE") {
      j_params[p.first] = get<std::vector<double>>(p.first);
    } else if (param_type == "NSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE") {
      auto sval = get<std::string>(p.first);
      if (!sval.empty()) continue;
      j_params[p.first] = sval;
      // size t
    } else if (p.second.get_type() == "m") {
      j_params[p.first] = get<size_t>(p.first);
    }
  }
  j["response_parameters"]=j_params;
}
