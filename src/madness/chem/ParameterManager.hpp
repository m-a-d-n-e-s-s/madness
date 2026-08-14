#pragma once
#include <madness/chem/CCParameters.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/ResponseParameters.hpp>
#include <madness/chem/TDHF.h>
#include <madness/chem/oep.h>
#include <madness/mra/QCCalculationParametersBase.h>

#include <type_traits>

using namespace madness;
using path = std::filesystem::path;

struct OptimizationParameters : public QCCalculationParametersBase {
  static constexpr const char *tag = "optimization";
  OptimizationParameters(const OptimizationParameters &other) = default;

  OptimizationParameters(World &world, const commandlineparser &parser)
      : OptimizationParameters() {
    read_input_and_commandline_options(world, parser, tag);
  }
  OptimizationParameters() {
    // No `method` here: the reference method comes from --wf=<scf|nemo>, and
    // --optimize says to optimize its geometry. Naming it twice invites the two
    // to disagree.
    initialize<int>("maxiter", 20, "maximum number of geometry steps");
    initialize<bool>("initial_hessian", false,
                     "compute inital hessian for optimization");
    initialize<std::string>("algopt", "bfgs", "hessian update used by MolOpt",
                            {"bfgs", "sr1"});
    initialize<double>("maxstep", 0.1,
                       "maximum step in any cartesian coordinate (a.u.)");
    // Convergence thresholds and assumed precisions. The defaults below are
    // only fallbacks: OptimizeDriver derives all five from the wavefunction
    // threshold (protocol().back()) with set_derived_value, which beats a default
    // but yields to anything the deck sets -- so a deck can tighten or loosen any
    // one of them, and otherwise they track the accuracy of the calculation.
    initialize<double>("etol", 1.e-5, "convergence: energy change");
    initialize<double>("gtol", 1.e-4, "convergence: maximum gradient element");
    initialize<double>("xtol", 1.e-4, "convergence: maximum cartesian step");
    initialize<double>("value_precision", 1.e-5, "assumed precision of the energy");
    initialize<double>("gradient_precision", 1.e-4,
                       "assumed precision of the gradient");
  }

  std::string get_tag() const override { return std::string(tag); }

  using QCCalculationParametersBase::read_input_and_commandline_options;
  using QCCalculationParametersBase::print;

  [[nodiscard]] int get_maxiter() const { return get<int>("maxiter"); }
  [[nodiscard]] bool get_initial_hessian() const {
    return get<bool>("initial_hessian");
  }
  [[nodiscard]] std::string get_algopt() const {
    return get<std::string>("algopt");
  }
  [[nodiscard]] double get_maxstep() const { return get<double>("maxstep"); }
  [[nodiscard]] double get_etol() const { return get<double>("etol"); }
  [[nodiscard]] double get_gtol() const { return get<double>("gtol"); }
  [[nodiscard]] double get_xtol() const { return get<double>("xtol"); }
  [[nodiscard]] double get_value_precision() const {
    return get<double>("value_precision");
  }
  [[nodiscard]] double get_gradient_precision() const {
    return get<double>("gradient_precision");
  }
};

/// Deck-level `io` block (roadmap change 5): run-wide I/O configuration.
/// `backend` selects the restart-archive family every task that persists
/// MRA state writes (today: response; moldft & friends as they adopt it) —
/// it generalizes the response-only `response.hdf5` knob, which is kept as
/// a working alias. Future parallel-HDF5 knobs (collective I/O, per-dataset
/// compression) belong in this block, not per-task.
struct IOParameters : public QCCalculationParametersBase {
  static constexpr const char *tag = "io";
  IOParameters(const IOParameters &other) = default;

  IOParameters(World &world, const commandlineparser &parser)
      : IOParameters() {
    read_input_and_commandline_options(world, parser, tag);
  }
  IOParameters() {
    initialize<std::string>(
        "backend", "native",
        "restart-archive backend for tasks that support it: native "
        "(MADNESS parallel archives) or hdf5 (single .h5 blobs; requires a "
        "-DMADNESS_ENABLE_HDF5=ON build)",
        {"native", "hdf5"});
  }

  std::string get_tag() const override { return std::string(tag); }

  using QCCalculationParametersBase::read_input_and_commandline_options;

  void print() const {
    madness::print("------------IO Parameters---------------");
    madness::print("Backend: ", get<std::string>("backend"));
    madness::print("-------------------------------------------");
  }

  [[nodiscard]] std::string backend() const {
    return get<std::string>("backend");
  }
  [[nodiscard]] bool hdf5() const { return backend() == "hdf5"; }
};

template <typename... Groups> class ParameterManager {
  std::tuple<Groups...> groups_;
  commandlineparser parser_;
  nlohmann::json all_input_json_;

  World &world_;

  // helper to invoke each group’s JSON export:
  template <typename G> void addGroupJson() {
    const auto &g = std::get<G>(groups_);
    auto j = g.to_json_if_precedence("defined");
    // if (world_.rank() == 0) {
    //   madness::print("Group: ", G::tag, " JSON: ", j.dump(4));
    // }
    if (!j.is_null())
      all_input_json_[G::tag] = j;
  }

public:
  ParameterManager() : world_(World::get_default()) {}

  /// "Master" ctor: takes any single intput file, JSON or plain-text
  // ParameterManager(World &w, const path &filename) : world_(w) {
  ParameterManager(World &w, const commandlineparser &parser)
      : parser_(parser), world_(w) {
    // parser_.set_keyval("input", filename);
    //
    const bool user_defined_prefix = parser_.key_exists("user_defined_prefix");
    std::string inputfile = parser.value("input");
    std::string prefix_from_input = commandlineparser::remove_extension(
        commandlineparser::base_name(inputfile));

    const path &filename = parser_.value("input");

    if (is_json_file(filename)) {
      auto j = read_json_file(filename);
      initFromJson(j);
      // invoke each group’s JSON parser:
    } else {
      // plain-text file
      initFromText(filename);
    }

    // Prefix precedence, highest first: --prefix on the command line, then a
    // `prefix` in the deck's dft block, then the input file's base name, then
    // the "mad" default. This must be layered in AFTER the deck is parsed:
    // set_derived_value does not override an already-defined value, so a
    // deck-level prefix survives while the file-name fallback still applies
    // when the deck is silent. CalculationParameters::prefix() is the single
    // source of truth (see prefix() below) — computing it up front into a
    // separate member shadowed the deck value.
    auto &cparam = get<CalculationParameters>();
    if (user_defined_prefix) {
      cparam.set_user_defined_value("prefix", parser_.value("prefix"));
    } else if (prefix_from_input != "input") {
      cparam.set_derived_value("prefix", prefix_from_input);
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
  [[nodiscard]] const nlohmann::json &getAllInputJson() const {
    return all_input_json_;
  }

  /// access a particular group by type:
  template <typename G> const G &get() const { return std::get<G>(groups_); }
  template <typename G> G &get() { return std::get<G>(groups_); }
  template <typename G> void set(const G &g) { std::get<G>(groups_) = g; }
  /// The run-level output prefix. Delegates to CalculationParameters so a
  /// `prefix` set in the deck (or later via set_derived_value) is honored by
  /// everything that names files or directories from it.
  std::string prefix() const { return get<CalculationParameters>().prefix(); }

  /// pretty-print everything
  void print_all() const { (print_group_if_defined<Groups>(), ...); }

private:
  void initFromJson(const nlohmann::json &j) {
    (
        [&] {
          if (j.contains(Groups::tag)) {
            std::get<Groups>(groups_).from_json(j.at(Groups::tag));
          }
          // Layer command-line overrides (--group=key=val) ON TOP of the JSON,
          // exactly as initFromText gets them through the file+CLI parser.
          // Without this, JSON-format decks silently ignored CLI overrides
          // (review MED). Applies to defaults too when the tag is absent.
          // Only the QCCalculationParametersBase-derived groups have the
          // override machinery — Molecule (also in the tuple) does not, so
          // guard with if constexpr.
          if constexpr (std::is_base_of_v<QCCalculationParametersBase, Groups>) {
            std::get<Groups>(groups_).read_commandline_options(world_, parser_,
                                                               Groups::tag);
          }
        }(),
        ...);

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

  template <typename G> void print_group_if_defined() const {
    const auto &g = std::get<G>(groups_);
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
                     OEP_Parameters, CCParameters, TDHFParameters, Molecule,
                     IOParameters>;
