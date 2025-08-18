#pragma once
#include <madness/chem/CCParameters.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/oep.h>
#include <madness/mra/QCCalculationParametersBase.h>

#include <madness/chem/ResponseParameters.hpp>

using namespace madness;
using path = std::filesystem::path;

struct OptimizationParameters : public QCCalculationParametersBase {
  static constexpr char const *tag = "optimization";
  OptimizationParameters(const OptimizationParameters &other) = default;

  OptimizationParameters(World &world, const commandlineparser &parser)
      : OptimizationParameters() {
    read_input_and_commandline_options(world, parser, tag);
  }
  OptimizationParameters() {
    initialize<int>("maxiter", 20, "optimization maxiter");

    initialize<bool>("initial_hessian", false,
                     "compute inital hessian for optimization");
    initialize<std::string>("algopt", "bfgs", "algorithm used for optimization",
                            {"bfgs", "cg"});
    initialize<double>("value_precision", 1.e-5, "value precision");
    initialize<double>("gradient_precision", 1.e-4, "gradient precision");
    initialize<bool>("geometry_tolerence", false, "geometry tolerance");
  }

  std::string get_tag() const override { return std::string(tag); }

  using QCCalculationParametersBase::read_input_and_commandline_options;

  void print() const {
    madness::print("------------Optimization Parameters---------------");
    madness::print("Maxiter: ", get<int>("maxiter"));
    madness::print("Initial Hessian: ", get<bool>("initial_hessian"));
    madness::print("Algorithm: ", get<std::string>("algopt"));
    madness::print("Value Precision: ", get<double>("value_precision"));
    madness::print("Gradient Precision: ", get<double>("gradient_precision"));
    madness::print("Geometry Tolerance: ", get<bool>("geometry_tolerence"));
    madness::print("-------------------------------------------");
  }

  [[nodiscard]] std::string get_method() const {
    return get<std::string>("method");
  }
  [[nodiscard]] int get_maxiter() const { return get<int>("maxiter"); }
  [[nodiscard]] bool get_initial_hessian() const {
    return get<bool>("initial_hessian");
  }
  [[nodiscard]] std::string get_algopt() const {
    return get<std::string>("algopt");
  }
  [[nodiscard]] double get_value_precision() const {
    return get<double>("value_precision");
  }
  [[nodiscard]] double get_gradient_precision() const {
    return get<double>("gradient_precision");
  }
  [[nodiscard]] bool get_geometry_tolerence() const {
    return get<bool>("geometry_tolerence");
  }
};

template <typename... Groups>
class ParameterManager {
  std::tuple<Groups...> groups_;
  commandlineparser parser_;
  nlohmann::json all_input_json_;
  std::string prefix_;

  World &world_;

  // helper to invoke each group’s JSON export:
  template <typename G>
  void addGroupJson() {
    auto const &g = std::get<G>(groups_);
    auto j = g.to_json_if_precedence("defined");
    if (world_.rank() == 0) {
      madness::print("Group: ", G::tag, " JSON: ", j.dump(4));
    }
    if (!j.is_null()) all_input_json_[G::tag] = j;
  }

 public:
  ParameterManager() : world_(World::get_default()) {}

  /// "Master" ctor: takes any single intput file, JSON or plain-text
  // ParameterManager(World &w, const path &filename) : world_(w) {
  ParameterManager(World &w, const commandlineparser &parser)
      : world_(w), parser_(parser) {
    // parser_.set_keyval("input", filename);
    //
    std::string inputfile = parser.value("input");
    std::string prefix = commandlineparser::remove_extension(
        commandlineparser::base_name(inputfile));
    if (prefix != "input") {
      prefix_ = prefix;
    } else {
      prefix_ = "mad";
    }

    const path &filename = parser_.value("input");

    if (is_json_file(filename)) {
      auto j = read_json_file(filename);
      initFromJson(j);
      // invoke each group’s JSON parser:
    } else {
      // plain-text file
      initFromText(filename);
    }
    set_derived_values();
  }

  /// here comes some logic for the calculation, e.g. the number of electrons
  /// derived from the molecule
  void set_derived_values() {
    this->get<CalculationParameters>().set_derived_values(
        this->get<Molecule>());
  }

  /// dump out the merged JSON
  [[nodiscard]] nlohmann::json const &getAllInputJson() const {
    return all_input_json_;
  }

  /// access a particular group by type:
  template <typename G>
  G const &get() const {
    return std::get<G>(groups_);
  }
  template <typename G>
  G &get() {
    return std::get<G>(groups_);
  }
  template <typename G>
  void set(G const &g) {
    std::get<G>(groups_) = g;
  }
  std::string prefix() const { return prefix_; }

  /// pretty-print everything
  void print_all() const { (print_group_if_defined<Groups>(), ...); }

 private:
  void initFromJson(nlohmann::json const &j) {
    (
        [&] {
          if (j.contains(Groups::tag)) {
            if (world_.rank() == 0) {
              madness::print("Group: ", Groups::tag,
                             " JSON: ", j.at(Groups::tag).dump(4));
            }
            std::get<Groups>(groups_).from_json(j.at(Groups::tag));
          }
        }(),
        ...);

    std::string inputfile = parser_.value("input");
    std::string prefix = commandlineparser::remove_extension(
        commandlineparser::base_name(inputfile));
    if (prefix != "input") {
      prefix_ = prefix;
    } else {
      prefix_ = "mad";
    }

    all_input_json_ = j;
  }
  // 1) read from a plain-text “.inp” file
  void initFromText(const path &filename) {
    // parser_.set_keyval("input", filename);
    // invoke each group’s file+CLI parser:
    ((void)(std::get<Groups>(groups_) = Groups(world_, parser_)), ...);

    // collect JSON for any defined keys:
    ((void)addGroupJson<Groups>(), ...);
  }

  template <typename G>
  void print_group_if_defined() const {
    auto const &g = std::get<G>(groups_);
    // grab only the user-defined values:
    auto j = g.to_json_if_precedence("defined");
    // json.empty() is true if no user-defined values
    if (!j.empty()) {
      g.print();
    }
  }
  static bool is_json_file(const path &f) {
    std::ifstream input_file_stream(f);
    bool is_json = json::accept(input_file_stream);
    input_file_stream.close();
    return is_json;
  }
  static json read_json_file(const path &input_file) {
    std::ifstream input_file_stream(input_file);
    auto j = json::parse(input_file_stream);
    input_file_stream.close();
    return j;
  }
};

// Define a concrete aliased ParameterManager type
using Params =
    ParameterManager<CalculationParameters, ResponseParameters,
                     Nemo::NemoCalculationParameters, OptimizationParameters,
                     OEP_Parameters, TDHFParameters, CCParameters, Molecule>;
