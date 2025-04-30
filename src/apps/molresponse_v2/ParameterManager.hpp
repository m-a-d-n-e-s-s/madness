#include <QCCalculationParametersBase.h>

using namespace madness;

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
class ParameterManager : public madness::QCCalculationParametersBase {
  std::tuple<Groups...> groups_;
  commandlineparser parser_;
  nlohmann::json all_input_json_;

  World &world_;

  // helper to invoke each group’s JSON export:
  template <typename G> void addGroupJson() {
    auto const &g = std::get<G>(groups_);
    auto j = g.to_json_if_precedence("defined");
    if (!j.is_null())
      all_input_json_[G::tag] = j;
  }

public:
  // 1) read from a plain-text “.inp” file
  ParameterManager(World &w, std::string const &filename) : world_(w) {
    parser_.set_keyval("input", filename);
    // invoke each group’s file+CLI parser:
    ((void)(std::get<Groups>(groups_) = Groups(world_, parser_)), ...);

    // collect JSON for any defined keys:
    ((void)addGroupJson<Groups>(), ...);
  }

  // 2) or read from an existing JSON
  ParameterManager(World &w, nlohmann::json const &j) : world_(w) {
    (
        [&] {
          if (j.contains(Groups::tag)) {
            std::get<Groups>(groups_).from_json(j.at(Groups::tag));
          }
        }(),
        ...);
  }

  /// dump out the merged JSON
  [[nodiscard]] nlohmann::json const &getAllInputJson() const {
    return all_input_json_;
  }

  /// access a particular group by type:
  template <typename G> G const &get() const { return std::get<G>(groups_); }


  /// pretty-print everything
  void print_all() const { (std::get<Groups>(groups_).print(), ...); }
};
