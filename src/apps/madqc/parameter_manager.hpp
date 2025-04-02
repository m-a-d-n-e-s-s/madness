#ifndef MADNESS_DFT_RESPONSE_PARAMETER_MANAGER_HPP
#define MADNESS_DFT_RESPONSE_PARAMETER_MANAGER_HPP
#include "tasks.hpp"
#include <apps/molresponse/response_parameters.h>
#include <filesystem>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <utility>

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

using namespace madness;

struct OptimizationParameters : public QCCalculationParametersBase {
  OptimizationParameters(const OptimizationParameters &other) = default;

  OptimizationParameters(World &world, const commandlineparser &parser)
      : OptimizationParameters() {
    read_input_and_commandline_options(world, parser, "optimization");
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
    madness::print("Method: ", get<std::string>("method"));
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
class ParameterManager {
private:
  path input_file_path;
  path input_file_json_path;

  path final_input_file_path = "final_input";
  path final_input_json_path = "final_input.json";

  json all_input_json;
  commandlineparser parser;

  Molecule molecule;

  TaskParameters task_params{};
  CalculationParameters moldft_params{};
  ResponseParameters molresponse_params{};
  OptimizationParameters optimization_params{};
  /**
   * Reads the chemical parameters from standard madness input file
   * and writes the input to json if parameters are defined
   * Available parameters are for
   *  - dft
   *  - response
   *
   *
   * @param world : madness world
   * @param input_file : path to the input file
   */

  void read_input_file(World &world, const path &input_file) {
    parser.set_keyval("input", input_file.string());

    /*task_params = TaskParameters(world, this->parser);*/
    /*if (world.rank() == 0) {*/
    /*  task_params.print();*/
    /*}*/
    /*auto task_json = task_params.to_json_if_precedence("defined");*/
    /*if (world.rank() == 0) {*/
    /*  print("Task json: ", task_json.dump(4));*/
    /*}*/
    /*if (!task_json.is_null()) {*/
    /*  all_input_json["task"] = task_json;*/
    /*}*/

    moldft_params = CalculationParameters(world, this->parser);
    // keeps a json version of the inputs
    // only outputs if the parameters are defined in the input file (aka
    // precedence is defined)
    auto moldft_json = moldft_params.to_json_if_precedence("defined");
    if (!moldft_json.is_null()) {
      all_input_json["dft"] = moldft_params.to_json_if_precedence("defined");
    }
    molresponse_params = ResponseParameters(world, this->parser);
    auto molresponse_json = molresponse_params.to_json_if_precedence("defined");
    if (!molresponse_json.is_null()) {
      all_input_json["response"] = molresponse_json;
    }
  }
  /** Reads the molecule from the json file and orients it
   * if paramater is set to true.  The molecule is then copied
   * back to the json file to ensure that the geometry is
   * consistent
   *
   *
   * @param j
   */
  void read_molecule_from_json_and_orient(json &j) {
    if (j.contains("molecule")) {
      molecule.from_json(j["molecule"]);
      j["molecule"] = molecule.to_json();
    } else {
      throw std::runtime_error("Molecule not found in input file");
    }
  }

  /**
   * Reads the chemical parameters from the json file
   * Available parameters are for
   *  - dft
   *  - response
   *
   *
   * @param j
   */
  void input_from_json(const json &j) {
    if (j.contains("task")) {
      all_input_json["task"] = j["task"];
      json &j = all_input_json["task"];
      task_params = j.get<TaskParameters>();
    }
    if (j.contains("dft")) {
      all_input_json["dft"] = j["dft"];
      moldft_params.from_json(j["dft"]);
    }
    if (j.contains("response")) {
      all_input_json["response"] = j["response"];
      molresponse_params.from_json(j["response"]);
    }
    if (j.contains("optimize")) {
      all_input_json["optimize"] = j["optimize"];
      optimization_params.from_json(j["optimize"]);
    }
  }

public:
  ParameterManager() = default;

  void print_file_paths() const {
    ::print("------------Parameter Manager---------------");
    ::print("Input File Path: ", input_file_path);
    ::print("-------------------------------------------");
  }

  static void help() {
    print_header2("help page for MADNESS DFT and Response Properties Code ");
  }
  void print_params() const {
    ::print("------------Parameter Manager---------------");
    ::print("Molecule: ");
    molecule.print();
    ::print("Task Parameters: ");
    task_params.print();
    ::print("Moldft Parameters: ");
    moldft_params.print();
    ::print("Molresponse Parameters: ");
    molresponse_params.print();
    ::print("-------------------------------------------");
  }

  /**
   * Broadcasts the parameters to all the nodes from 0
   *
   * @param world
   */

  explicit ParameterManager(World &world, const path &input_file) {
    // First read the molecule file because we have one

    std::ifstream input_file_stream(input_file);
    bool is_json = json::accept(input_file_stream);
    input_file_stream.close();
    if (world.rank() == 0) {
      print("Input file path: ", input_file);
    }

    if (is_json) {
      input_file_stream.open(input_file);
      all_input_json = json::parse(input_file_stream);
      input_file_stream.close();
      input_from_json(all_input_json);
      read_molecule_from_json_and_orient(all_input_json);
    } else {
      read_molecule_file(input_file);
      all_input_json["molecule"] = molecule.to_json();
      read_input_file(world, input_file);
    }

    if (world.rank() == 0) {
      print(all_input_json.dump(4));
    }
  }

  void read_molecule_file(const path &mol_file) {
    bool is_json = is_json_file(mol_file);
    if (is_json) {
      print("Molecule file is json");

      std::ifstream mol_stream(mol_file);
      auto mol_json = json::parse(mol_stream);
      molecule.from_json(mol_json["molecule"]);
      mol_stream.close();
    } else {
      std::ifstream mol_stream(mol_file);
      molecule.read(mol_stream);
      mol_stream.close();
    }
    all_input_json["molecule"] = molecule.to_json();
  }

  static bool is_json_file(const path &input_file) {
    std::ifstream input_file_stream(input_file);
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

  // The intention of this second constructor is to allow the user to
  // specify the molecule file and the input file separately.  The specific
  // use cases are when one wants to create a database of molecule
  // calculations all with the same basic input.
  explicit ParameterManager(World &world,
                            const std::pair<path, path> &input_files) {
    auto [input_file, mol_file] = input_files;

    read_molecule_file(mol_file);
    if (world.rank() == 0) {
      print("Input file path: ", input_file);
    }
    auto is_json = is_json_file(input_file);

    if (is_json) {
      if (world.rank() == 0) {
        print(input_file, " is a json file");
      }
      auto input_j = read_json_file(input_file);
      input_from_json(input_j);
    } else {
      if (world.rank() == 0) {
        print(input_file, " is not json file");
      }
      read_input_file(world, input_file);
    }
    if (world.rank() == 0) {
      print("Input files: ", input_file, mol_file);
      print(all_input_json.dump(4));
    }
  }

  [[nodiscard]] json get_input_json() const { return all_input_json; }
  [[nodiscard]] auto get_task_params() const -> const TaskParameters & {
    return task_params;
  }
  [[nodiscard]] auto
  get_optimization_params() const -> const OptimizationParameters & {
    return optimization_params;
  }

  [[nodiscard]] auto
  get_moldft_params() const -> const CalculationParameters & {
    return moldft_params;
  }
  [[nodiscard]] auto
  get_molresponse_params() const -> const ResponseParameters & {
    return molresponse_params;
  }
  [[nodiscard]] auto get_molecule() const -> const Molecule & {
    return molecule;
  }

  void write_moldft_json(std::ostream &os) {
    os << std::setw(4) << all_input_json["dft"];
    os << std::setw(4) << all_input_json["molecule"];
  }

  void write_response_json(std::ostream &os) {
    os << std::setw(4) << all_input_json["response"];
  }
};

// create a helper class for checking equivalence of two parameter class
//
bool operator==(const ParameterManager &lhs, const ParameterManager &rhs) {
  return lhs.get_input_json() == rhs.get_input_json();
}

#endif
