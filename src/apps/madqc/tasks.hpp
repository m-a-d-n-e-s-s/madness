
// Tasks.hpp defines the available tasks in MADQC which are
// energy, optimization, property
//
//
// We will implement energy and optimization tasks using a strategy design
// pattern
//

#ifndef MADQC_TASKS_HPP
#define MADQC_TASKS_HPP

#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>
#include <madness/world/print.h>

namespace madness {

using property_map = std::map<std::string, std::map<std::string, bool>>;

struct TaskParameters {

  std::string driver;
  std::string method;
  property_map properties;

  void print()const  {
    madness::print("------------Task Parameters---------------");
    madness::print("Driver: ", driver);
    madness::print("Method: ", method);
    madness::print("Properties: ");
    for (const auto& [model, props] : properties) {
      madness::print("Model: ", model);
      for (const auto& [prop, value] : props) {
        madness::print(prop, " : ", value);
      }
    }
    madness::print("-------------------------------------------");
  }
};

void to_json(nlohmann::json& j, const TaskParameters& t) {
  j = nlohmann::json{{"driver", t.driver}, {"method", t.method}, {"properties", t.properties}};
}

void from_json(const nlohmann::json& j, TaskParameters& t) {
  j.at("driver").get_to(t.driver);
  j.at("method").get_to(t.method);
  j.at("properties").get_to(t.properties);
}

/*struct TaskParameters : public QCCalculationParametersBase {*/
/*  TaskParameters(const TaskParameters& other) = default;*/
/**/
/*  TaskParameters(World& world, const commandlineparser& parser) : TaskParameters() {*/
/*    read_input_and_commandline_options(world, parser, "task");*/
/*  }*/
/*  TaskParameters() {*/
/*    initialize<std::string>("driver", "property", "Type of calculation that will be run",*/
/*                            {"property", "optimize", "custom"});*/
/*    initialize<std::string>("model", "dft", "Select the application to use for the calculation",*/
/*                            {"dft", "response", "mp2", "cis", "cc2", "mp3"});*/
/*    initialize<property_map>("properties", {}, "Properties to be calculated for each model");*/
/*  }*/
/**/
/*  using QCCalculationParametersBase::read_input_and_commandline_options;*/
/**/
/*  [[nodiscard]] std::string get_driver() const { return get<std::string>("driver"); }*/
/*  [[nodiscard]] std::string get_method() const { return get<std::string>("method"); }*/
/*};*/

}  // namespace madness
#endif  // MADQC_TASKS_HPP
