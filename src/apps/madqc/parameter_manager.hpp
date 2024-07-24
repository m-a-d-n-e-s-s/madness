#include <apps/molresponse/response_parameters.h>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <utility>
#include "tasks.hpp"
#include "utils.hpp"

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

using namespace madness;

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

 public:
  ParameterManager() = default;

  void print_file_paths() const {
    ::print("------------Parameter Manager---------------");
    ::print("Input File Path: ", input_file_path);
    ::print("Final Input File Path: ", final_input_file_path);
    ::print("Input File Json Path: ", input_file_json_path);
    ::print("-------------------------------------------");
  }

  static void help() {
    print_header2("help page for MADNESS DFT and Response Properties Code ");
    print(
        "This code is designed to run DFT and Response Property "
        "calculations");
    print(
        "Within the input one defines both the ground and response "
        "calculations in the input file by specifiying the dft and response "
        "blocks");
    print(
        "By defining the quadratic block one can compute quadratic response "
        "properties such as the hyperpolarizability");
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

  /** Reads the molecule from the json file and orients it
   * if paramater is set to true.  The molecule is then copied
   * back to the json file to ensure that the geometry is
   * consistent
   *
   *
   * @param j
   */
  void read_molecule_from_json_and_orient(json& j) {
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
  void input_from_json(const json& j) {
    if (j.contains("task")) {
      all_input_json["task"] = j["task"];
      task_params.from_json(j["task"]);
    }
    if (j.contains("dft")) {
      all_input_json["dft"] = j["dft"];
      moldft_params.from_json(j["dft"]);
    }
    if (j.contains("response")) {
      all_input_json["response"] = j["response"];
      molresponse_params.from_json(j["response"]);
    }
  }

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

  void read_input_file(World& world, const path& input_file) {
    parser.set_keyval("input", input_file.string());

    task_params = TaskParameters(world, this->parser);
    if (world.rank() == 0) {
      task_params.print();
    }
    auto task_json = task_params.to_json_if_precedence("defined");
    if (world.rank() == 0) {
      print("Task json: ", task_json.dump(4));
    }
    if (!task_json.is_null()) {
      all_input_json["task"] = task_json;
    }

    moldft_params = CalculationParameters(world, this->parser);
    // keeps a json version of the inputs
    // only outputs if the parameters are defined in the input file (aka precedence is defined)
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

  /**
   * Broadcasts the parameters to all the nodes from 0
   *
   * @param world
   */

  explicit ParameterManager(World& world, const path& input_file) {
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

  void read_molecule_file(const path& mol_file) {
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

  static bool is_json_file(const path& input_file) {
    std::ifstream input_file_stream(input_file);
    bool is_json = json::accept(input_file_stream);
    input_file_stream.close();
    return is_json;
  }

  static json read_json_file(const path& input_file) {
    std::ifstream input_file_stream(input_file);
    auto j = json::parse(input_file_stream);
    input_file_stream.close();
    return j;
  }

  // The intention of this second constructor is to allow the user to
  // specify the molecule file and the input file separately.  The specific
  // use cases are when one wants to create a database of molecule
  // calculations all with the same basic input.
  explicit ParameterManager(World& world,
                            const std::pair<path, path>& input_files) {
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

  [[nodiscard]] auto get_moldft_params() const -> const CalculationParameters& {
    return moldft_params;
  }
  [[nodiscard]] auto get_molresponse_params() const
      -> const ResponseParameters& {
    return molresponse_params;
  }
  [[nodiscard]] auto get_molecule() const -> const Molecule& {
    return molecule;
  }

  void write_moldft_json(std::ostream& os) {
    os << std::setw(4) << all_input_json["dft"];
    os << std::setw(4) << all_input_json["molecule"];
  }

  void write_response_json(std::ostream& os) {
    os << std::setw(4) << all_input_json["response"];
  }
};

// create a helper class for checking equivalence of two parameter class
//
bool operator==(const ParameterManager& lhs, const ParameterManager& rhs) {
  return lhs.get_input_json() == rhs.get_input_json();
}
