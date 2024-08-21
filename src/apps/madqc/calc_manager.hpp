#ifndef CALC_MANAGER_HPP
#define CALC_MANAGER_HPP

#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <apps/molresponse/FrequencyResponse.hpp>
#include <apps/molresponse/ResponseExceptions.hpp>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <utility>
#include <vector>
#include "tensor_json.hpp"
#include "utils.hpp"

using json = nlohmann::json;
using path = std::filesystem::path;
using ResponseInput = std::tuple<std::string, std::string, std::vector<double>>;

// A word for a class which takes input from another class and
// therefore has to define the input interface
class InputInterface {

 protected:
  vector<std::string> input_names;

 public:
  explicit InputInterface(const std::vector<std::string>& names) : input_names(names) {}
};

class CalculationStrategy {
 public:
  virtual json calcPaths(const path& root) = 0;
  virtual void compute(World& world, const json& calc_paths) = 0;
  virtual ~CalculationStrategy() = default;
  // Virtual function that outputs the properties of the calculation// this requires that each derived class defines the SCF/Calc object as a member variable
  // to get access to the properties of the calculation
  //virtual json output_properites(std::vector<std::string> properties) = 0;
  //
  // Doesn't work because SCF does not have a default constructor... if it did, then I could hold the SCF object as a member variable
  //
  // For now, For each calc_strategy we need to the define the ouput_properties function that takes the properties requested and the object used to compute them...
  //
  //
  //
  // All calculations should have a name where they output their data in the json object
  CalculationStrategy(std::string calc_name, vector<std::string> properties)
      : name(std::move(calc_name)), requested_properties(std::move(properties)){};
  CalculationStrategy() = default;

 protected:
  std::string name;
  std::vector<std::string> requested_properties;
};

class CompositeCalculationStrategy : public CalculationStrategy {
 private:
  std::vector<std::unique_ptr<CalculationStrategy>> strategies;

 public:
  void addStrategy(std::unique_ptr<CalculationStrategy> strategy) { strategies.push_back(std::move(strategy)); }

  json calcPaths(const path& root) override {
    json result = {};
    for (const auto& strategy : strategies) {
      json paths = strategy->calcPaths(root);
      // I know there should only be one key in the json object corresponding the calc_name.
      // Insert json as j[calc_name]={}
      for (auto& [key, value] : paths.items()) {
        result[key] = value;
      }
    }
    print(result.dump(4));

    // for the keys in the json object, create the directories
    // if they don't exist
    //
    for (const auto& [key, value] : result.items()) {
      if (!std::filesystem::exists(key)) {
        std::filesystem::create_directory(key);
      }
    }

    for (const auto& calc_name : result) {
      auto calc_dirs = calc_name["calculation"];
      print(calc_dirs);
      for (const auto& calc_dir : calc_dirs) {
        if (!std::filesystem::exists(calc_dir)) {
          std::filesystem::create_directory(calc_dir);
        }
      }
    }

    return result;
  }

  void compute(World& world, const json& calc_paths) override {
    for (const auto& strategy : strategies) {
      strategy->compute(world, calc_paths);
    }
  }
};

class MoldftCalculationStrategy : public CalculationStrategy {

  CalculationParameters parameters;
  Molecule molecule;
  json paths;
  path output_path;

  std::vector<std::string> available_properties = {"energy", "gradient", "dipole"};

  static json compute_property(World& world, const std::string& property, SCF& calc, double energy) {
    json result;
    if (property == "energy") {
      result["energy"] = energy;
    } else if (property == "gradient") {
      auto gradient = calc.derivatives(world, calc.make_density(world, calc.aocc, calc.amo));
      result["gradient"] = tensor_to_json(gradient);  // here we need to make json representation of a tensor....
    } else if (property == "dipole") {
      auto dipole = calc.dipole(world, calc.make_density(world, calc.aocc, calc.amo));
      result["dipole"] = tensor_to_json(dipole);
    } else {
      throw std::runtime_error("Property not available");
    }
    return result;
  }

  void output_properites(World& world, SCF& calc, double energy) {

    // This is where a property interface would be useful. for properties use the properties interface
    // to get the properties of the calculation
    json results;
    paths[name]["output"]["properties"] = {};
    auto& output = paths[name]["output"]["properties"];

    for (const auto& property : requested_properties) {
      output[property] = compute_property(world, property, calc, energy);
    }
    if (world.rank() == 0) {
      print("output: ", output.dump(4));
      std::ofstream ofs(output_path);
      ofs << output.dump(4);
      ofs.close();
    }
  }

 public:
  // Notice here that I am passing the parameters of the calculation as well as name
  MoldftCalculationStrategy(const CalculationParameters& params, Molecule mol, std::string calc_name = "moldft",
                            std::vector<std::string> properties = {"energy"})
      : parameters(params),
        molecule(std::move(mol)),
        CalculationStrategy(std::move(calc_name), std::move(properties)){};

  json calcPaths(const path& root) override {
    // All calculations need to create a json with the top level key being the name of the calculation
    paths[name] = {};
    auto& moldft = paths[name];

    auto base_root = root / name;

    // these are the possible paths that are needed for the calculation
    // I need a way to inject the ability to compute different properties from the calculation.
    // For example, if we are doing a finite-difference of the dipole moment, we need to inject the ability
    // to compute the dipole moment.  This will normally be done by the parameters.
    // Should the constructor have a paramater that indicates which properties to compute?
    //
    // For example, [energy, gradient,dipole,quadrupole] etc.
    //
    //
    // For an optimization, we need to compute both the energy and the gradient
    //
    // For a finite-difference calculation, we need to compute a specific property
    //
    // Maybe the last stage of the calculation should be a post-processing step that computes the properties
    //
    // And the first step of calculation should be something that checks whether the properties can be computed
    // given a method.
    //
    moldft["calculation"] = base_root;
    moldft["restart"] = base_root / "moldft.restartdata.00000";
    moldft["output"] = {};
    moldft["output"]["calc_info"] = base_root / "moldft.calc_info.json";
    moldft["output"]["scf_info"] = base_root / "moldft.scf_info.json";
    moldft["output"]["properties"] = base_root / "output.json";
    output_path = moldft["output"]["properties"].get<path>();

    return paths;
  }

  void compute(World& world, const json& calc_paths) override {
    paths = calc_paths;
    // the first step is to look for the moldft paths
    json moldft_paths = paths[name];
    if (world.rank() == 0)
      print(moldft_paths.dump(4));

    // Get the paths from the json object
    path moldft_path = moldft_paths["calculation"];
    path moldft_restart = moldft_paths["restart"];
    path moldft_calc_info_path = moldft_paths["output"]["calc_info"];

    // print the relevant paths
    if (world.rank() == 0) {
      std::cout << "Calculation Path: " << moldft_paths.dump(4) << std::endl;
    }

    if (!std::filesystem::exists(moldft_path)) {
      std::filesystem::create_directory(moldft_path);
    }

    std::filesystem::current_path(moldft_path);

    json calcInfo;
    auto param1 = parameters;
    world.gop.broadcast_serializable(param1, 0);

    if (world.rank() == 0) {
      ::print("-------------Running moldft------------");
    }

    if (std::filesystem::exists(moldft_restart) && std::filesystem::exists(moldft_calc_info_path)) {
      // if both exist, read the calc_info json
      std::ifstream ifs(moldft_calc_info_path);

      auto moldft_calc_info = json::parse(ifs);
      if (world.rank() == 0) {
        std::cout << "time: " << moldft_calc_info["time"] << std::endl;
        std::cout << "MOLDFT return energy: " << moldft_calc_info["return_energy"] << std::endl;
      }
    } else {
      // if params are different run and if restart exists and if im asking to
      if (std::filesystem::exists(moldft_restart)) {
        param1.set_user_defined_value<bool>("restart", true);
      }
      world.gop.fence();
      if (world.rank() == 0) {

        json moldft_input_json = {};
        moldft_input_json["dft"] = parameters.to_json_if_precedence("defined");
        moldft_input_json["molecule"] = molecule.to_json();

        std::ofstream ofs("moldft.in");
        write_moldft_input(moldft_input_json, ofs);
        ofs.close();
      }
      world.gop.fence();
      commandlineparser parser;
      parser.set_keyval("input", "moldft.in");

      if (world.rank() == 0)
        ::print("input filename: ", parser.value("input"));

      print_meminfo(world.rank(), "startup");
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

      std::cout.precision(6);
      SCF calc(world, parser);
      if (world.rank() == 0) {
        ::print("\n\n");
        ::print(
            " MADNESS Hartree-Fock and Density Functional Theory "
            "Program");
        ::print(
            " ----------------------------------------------------------"
            "\n");
        calc.param.print("dft");
      }
      if (world.size() > 1) {
        calc.set_protocol<3>(world, 1e-4);
        calc.make_nuclear_potential(world);
        calc.initial_load_bal(world);
      }
      calc.set_protocol<3>(world, calc.param.protocol()[0]);
      MolecularEnergy E(world, calc);                                  // does
      double energy = E.value(calc.molecule.get_all_coords().flat());  // ugh!
      //
      calc.output_calc_info_schema();

      output_properites(world, calc, energy);
    }

    // Add actual calculation logic here
  }
};

/*
class MP2CalculationStrategy : public CalculationStrategy {
 public:
  void runCalculation(World& world, const json& paths) override {
    std::cout << "Running MP2 Calculation\n";
    std::cout << "Calc Directory: " << paths["mp2_calc_dir"] << "\n";
    std::cout << "Restart Path: " << paths["mp2_restart"] << "\n";
    std::cout << "Calc Info Path: " << paths["mp2_calc_info"] << "\n";
    // Add actual calculation logic here
  }
};
*/

class ResponseConfig {

 public:
  std::string calc_name = "response";
  std::string perturbation;
  std::string xc;
  std::vector<double> frequencies;
  static path restart_path(const std::filesystem::path& calc_path) {

    auto save_path = std::filesystem::path(calc_path);
    auto run_name = calc_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;
    save_path += ".00000";

    return save_path;
  }

  /**
     * generates the frequency response path using the format
     * [property]_[xc]_[1-100]
     *
     * where 1-100 corresponds a frequency of 1.100
     *
     * @param moldft_path
     * @param property
     * @param frequency
     * @param xc
     * @return
     */
  [[nodiscard]] auto calc_path(const path& root, const double& frequency) const -> std::filesystem::path {
    std::string s_frequency = std::to_string(frequency);
    auto sp = s_frequency.find('.');
    s_frequency = s_frequency.replace(sp, sp, "-");
    std::string run_name = this->perturbation + "_" + this->xc + "_" + s_frequency;
    return root / std::filesystem::path(run_name);
  }
  explicit ResponseConfig(std::string calc_name, ResponseInput input)
      : calc_name(std::move(calc_name)),
        perturbation(std::get<0>(input)),
        xc(std::get<1>(input)),
        frequencies(std::get<2>(input)) {}
  explicit ResponseConfig(ResponseInput input)
      : perturbation(std::get<0>(input)), xc(std::get<1>(input)), frequencies(std::get<2>(input)) {}
};

class LinearResponseStrategy : public CalculationStrategy, InputInterface {

  ResponseParameters parameters;
  std::string op;
  std::string xc;
  std::vector<double> freqs;
  std::string calc_name = "response";
  ResponseConfig config;
  json paths;
  path output_path;
  std::vector<std::string> available_properties = {"alpha"};

 public:
  explicit LinearResponseStrategy(const ResponseParameters& params, const ResponseInput& r_input,
                                  std::string name = "response",
                                  const std::vector<std::string>& input_names = {"moldft"})
      : parameters(params), calc_name(std::move(name)), config(name, r_input), InputInterface(input_names) {}

  json calcPaths(const path& root) override {

    json paths = {};
    paths[calc_name] = {};
    auto& response = paths[calc_name];
    response["frequencies"] = config.frequencies;
    response["calculation"] = {};
    response["output"] = {};
    response["restart"] = {};

    auto base_path = root / calc_name;

    response["properties"] = {};
    response["properties"]["alpha"] = base_path / "alpha.json";

    for (const auto& frequency : config.frequencies) {
      auto frequency_run_path = config.calc_path(base_path, frequency);
      auto restart_path = ResponseConfig::restart_path(frequency_run_path);
      response["calculation"].push_back(frequency_run_path);
      response["restart"].push_back(restart_path);
      response["output"].push_back(frequency_run_path / "response_base.json");
    }
    return paths;
  }

  static void append_to_alpha_json(const double& omega, const std::vector<std::string>& ij, const Tensor<double>& alpha,
                                   nlohmann::ordered_json& alpha_json) {
    auto num_unique_elements = ij.size();
    for (int i = 0; i < num_unique_elements; i++) {
      alpha_json["omega"].push_back(omega);
      alpha_json["ij"].push_back(ij[i]);
      alpha_json["alpha"].push_back(alpha[i]);
    }
  }
  static void add_alpha_i_to_json(nlohmann::ordered_json& alpha_i, nlohmann::ordered_json& alpha_json) {
    print(alpha_json.dump(4));

    alpha_json["omega"].insert(alpha_json["omega"].end(), alpha_i["omega"].begin(), alpha_i["omega"].end());
    alpha_json["ij"].insert(alpha_json["ij"].end(), alpha_i["ij"].begin(), alpha_i["ij"].end());
    alpha_json["alpha"].insert(alpha_json["alpha"].end(), alpha_i["alpha"].begin(), alpha_i["alpha"].end());
  }

  /**
     *
     * @param world
     * @param filename
     * @param frequency
     * @param property
     * @param xc
     * @param moldft_path
     * @param restart_path
     * @return
     
     */
  bool runFrequency(World& world, ResponseParameters& r_params, double frequency) {

    auto op = r_params.perturbation();

    // Set the response parameters

    if (world.rank() == 0) {

      json input_json = {};
      input_json["response"] = r_params.to_json_if_precedence("defined");
      std::ofstream out("response.in");
      write_json_to_input_file(input_json, {"response"}, out);
    }
    bool converged = false;
    // if rbase exists and converged I just return save path and true
    if (world.rank() == 0) {
      ::print("Checking if response has converged for frequency: ", frequency);

      if (std::filesystem::exists("response_base.json")) {
        {
          std::ifstream ifs("response_base.json");
          json response_base;
          ifs >> response_base;
          if (response_base["converged"] && response_base["precision"]["dconv"] == r_params.dconv())
            converged = true;
          ::print("Response has converged: ", converged);
        }
      }
    }
    world.gop.broadcast(converged, 0);
    // add logic to compute alpha.json if needed

    if (converged) {
      return true;
    } else {

      world.gop.fence();
      if (world.rank() == 0) {
        ::print("Running response calculation for frequency: ", frequency);
      }

      GroundStateCalculation ground_calculation{world, r_params.archive()};
      Molecule molecule = ground_calculation.molecule();
      r_params.set_ground_state_calculation_data(ground_calculation);
      r_params.set_derived_values(world, molecule);
      CalcParams calc_params = {ground_calculation, molecule, r_params};

      RHS_Generator rhs_generator;
      if (op == "dipole") {
        rhs_generator = dipole_generator;
      } else {
        rhs_generator = nuclear_generator;
      }
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
      FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
      if (world.rank() == 0) {
        ::print("\n\n");
        ::print(
            " MADNESS Time-Dependent Density Functional Theory Response "
            "Program");
        ::print(" ----------------------------------------------------------\n");
        ::print("\n");
        calc_params.molecule.print();
        ::print("\n");
        calc_params.response_parameters.print("response");
        // put the response parameters in a j_molrespone json object
        calc_params.response_parameters.to_json(calc.j_molresponse);
      }
      calc.solve(world);
      auto [omega, polar_omega] = calc.get_response_data();
      // flatten polar_omega

      auto alpha = polar_omega.flat();

      std::vector<std::string> ij{"XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"};
      nlohmann::ordered_json alpha_json;
      alpha_json["omega"] = {};
      alpha_json["ij"] = {};
      alpha_json["alpha"] = {};

      append_to_alpha_json(omega, ij, alpha, alpha_json);
      calc.j_molresponse["properties"] = {};
      calc.j_molresponse["properties"]["alpha"] = alpha_json;

      // set protocol to the first
      if (world.rank() == 0) {
        // calc.time_data.to_json(calc.j_molresponse);
        calc.output_json();
      }
      // calc.time_data.print_data();
      return calc.j_molresponse["converged"];
    }
  }

  void compute(World& world, const json& all_paths) override {

    // First example of using input names
    paths = all_paths;
    auto moldft_name = input_names[0];
    print("moldft_name: ", moldft_name);
    print("all paths: ", paths.dump(4));

    auto& moldft_paths = paths[moldft_name];
    auto moldft_restart = moldft_paths["restart"].get<std::string>();
    moldft_restart = std::string(path(moldft_restart).replace_extension(""));

    auto& my_response = paths[calc_name];
    auto& calc_paths = my_response["calculation"];

    auto restart_paths = my_response["restart"].get<std::vector<std::string>>();
    auto output_paths = my_response["output"].get<std::vector<std::string>>();
    auto alpha_path = my_response["properties"]["alpha"].get<std::string>();

    size_t num_freqs = calc_paths.size();
    // I need to analyze alpha.json to see which frequencies are missing.
    // Or I just write them to seperate files and then combine them later?

    auto freqs = config.frequencies;

    bool last_converged = false;
    nlohmann::ordered_json alpha_json;
    alpha_json["omega"] = json::array();
    alpha_json["ij"] = json::array();
    alpha_json["alpha"] = json::array();

    for (size_t i = 0; i < num_freqs; i++) {

      auto freq_i = freqs[i];

      if (!std::filesystem::exists(calc_paths[i])) {
        std::filesystem::create_directory(calc_paths[i]);
      }
      std::filesystem::current_path(calc_paths[i]);

      print("current path: ", std::filesystem::current_path());
      print("calc path: ", calc_paths[i]);
      print("freq: ", freq_i);
      // if the last converged is true, then we can restart from the last save path

      bool restart = true;
      path save_path = restart_paths[i];  // current restart path aka save path
      path restart_path = (i > 0) ? restart_paths[i - 1] : restart_paths[i];
      std::string save_string = save_path.filename().stem();
      if (world.rank() == 0) {
        ::print("-------------Running response------------ at frequency: ", freq_i);
        ::print("moldft restart path", moldft_restart);
        ::print("restart path", restart_path);
        ::print("save path", save_path);
        ::print("save string", save_string);
      }

      ResponseParameters r_params = parameters;

      r_params.set_user_defined_value("omega", freq_i);
      r_params.set_user_defined_value("archive", moldft_restart);
      if (last_converged || i == 0) {

        if (world.rank() == 0) {
          r_params.set_user_defined_value("save", true);
          r_params.set_user_defined_value("save_file", save_string);
          if (restart) {  // if we are trying a restart calculation
            if (std::filesystem::exists(save_path)) {
              // if the save path exists then we know we can
              //  restart from the previous save
              r_params.set_user_defined_value("restart", true);
              r_params.set_user_defined_value("restart_file", save_string);
            } else if (std::filesystem::exists(restart_path)) {
              ::print("restart path exists", restart_path);
              r_params.set_user_defined_value("restart", true);
              ::print(restart_path.parent_path().stem());
              ::print(restart_path.filename().stem());

              // get the directory of the restart path
              auto new_restart_path = path("../") / restart_path.parent_path().stem() / restart_path.filename().stem();

              // format restart path to be ../restart_path/restart_path

              //
              ::print("new restart file: ", restart_path);
              r_params.set_user_defined_value("restart_file", new_restart_path.string());
              // Then we restart from the previous file instead
            } else {
              r_params.set_user_defined_value("restart", false);
            }
            // neither file exists therefore you need to start from fresh
          }
        }
      } else {
        throw Response_Convergence_Error{};
      }
      world.gop.broadcast_serializable(r_params, 0);
      last_converged = runFrequency(world, r_params, freq_i);

      nlohmann::ordered_json alpha_i;

      if (last_converged) {
        // read output_paths[i]
        std::ifstream ifs(output_paths[i]);
        json response_base_i;
        ifs >> response_base_i;
        alpha_i = response_base_i["properties"]["alpha"];
      }
      if (world.rank() == 0) {
        print("last converged: ", last_converged);
        print(alpha_i.dump(4));
      }

      // combine alpha_json_i with alpha_json
      add_alpha_i_to_json(alpha_i, alpha_json);
    }

    std::ofstream out_file(alpha_path);

    if (world.rank() == 0) {
      out_file << alpha_json.dump(4);
    }
  }
};

class ResponseHyper : public CalculationStrategy, InputInterface {

  ResponseParameters parameters;
  std::string op;
  std::string xc;
  std::vector<double> freqs;
  std::string calc_name = "hyper";
  ResponseConfig config;
  json paths;
  path output_path;
  std::vector<std::string> available_properties = {"beta"};

 public:
  explicit ResponseHyper(const ResponseParameters& params, const ResponseInput& r_input, std::string name = "hyper",
                         const std::vector<std::string>& input_names = {"moldft", "response"})
      : parameters(params), calc_name(std::move(name)), config(name, r_input), InputInterface(input_names) {
    print("input names: ", input_names);
    if (input_names.size() != 2) {
      throw std::runtime_error("ResponseHyper requires two input names");
    }
  }

  json calcPaths(const path& root) override {

    json paths = {};
    paths[calc_name] = {};
    auto& response = paths[calc_name];
    response["output"] = {};

    auto base_path = root / calc_name;

    response["calculation"] = base_path;
    response["output"]["properties"] = {};
    response["properties"]["beta"] = base_path / "beta.json";
    return paths;
  }

  void compute(World& world, const json& all_paths) override {

    paths = all_paths;

    auto moldft_name = input_names[0];
    auto response_name = input_names[1];
    if (world.rank() == 0) {
      print("moldft_name: ", moldft_name);
      print("response_name: ", response_name);
    }

    auto& moldft_paths = paths[moldft_name];
    auto moldft_path = moldft_paths["calculation"].get<std::string>();
    auto moldft_restart = moldft_paths["restart"].get<std::string>();
    moldft_restart = std::string(path(moldft_restart).replace_extension(""));

    auto& response_paths = paths[response_name];
    auto& calc_paths = response_paths["calculation"];

    std::vector<path> restart_paths = response_paths["restart"].get<std::vector<path>>();

    auto freqs = response_paths["frequencies"].get<std::vector<double>>();
    auto output_paths = response_paths["output"].get<std::vector<std::string>>();
    print("freqs: ", freqs);
    // Run the calculations
    bool last_converged = false;
    nlohmann::ordered_json beta_json;
    beta_json["Afreq"] = json::array();
    beta_json["Bfreq"] = json::array();
    beta_json["Cfreq"] = json::array();
    beta_json["A"] = json::array();
    beta_json["B"] = json::array();
    beta_json["C"] = json::array();

    beta_json["Beta"] = json::array();

    try {
      auto num_freqs = freqs.size();

      if (world.rank() == 0) {
        ::print("Running quadratic response calculations");
      }
      std::filesystem::current_path(moldft_path);
      RHS_Generator rhs_generator;
      if (config.perturbation == "dipole") {
        rhs_generator = dipole_generator;
      } else {
        rhs_generator = nuclear_generator;
      }
      if (world.rank() == 0) {
        ::print("Set up rhs generator");
      }

      auto set_hyperpolarizability_parameters = [&]() {
        ResponseParameters quad_parameters{};

        auto molresponse_params = parameters;

        quad_parameters.set_user_defined_value("quadratic", true);
        quad_parameters.set_user_defined_value("freq_range", molresponse_params.freq_range());
        quad_parameters.set_user_defined_value("hfexalg", molresponse_params.hfexalg());

        if (config.perturbation == "dipole") {
          quad_parameters.set_user_defined_value("dipole", true);
          quad_parameters.set_derived_value<size_t>("states", 3);
        }

        auto final_protocol = molresponse_params.protocol().back();
        quad_parameters.set_user_defined_value<vector<double>>("protocol", {final_protocol});

        return quad_parameters;
      };

      auto quad_parameters = set_hyperpolarizability_parameters();

      // auto calc_params = initialize_calc_params(world,
      // std::string("quad.in"));
      commandlineparser parser;
      std::string moldft_archive = "moldft.restartdata";
      GroundStateCalculation ground_calculation{world, moldft_archive};
      if (world.rank() == 0) {
        ground_calculation.print_params();
      }
      Molecule molecule = ground_calculation.molecule();
      quad_parameters.set_ground_state_calculation_data(ground_calculation);
      if (world.rank() == 0) {
        quad_parameters.print();
      }

      world.gop.fence();
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

      QuadraticResponse quad_calculation{
          world,
          {ground_calculation, molecule, quad_parameters},
          rhs_generator,
      };

      num_freqs = (freqs.size() / 2) + 1;

      for (int b = 0; b < num_freqs; b++) {
        for (int c = 0; c < num_freqs; c++) {

          ::print(world.rank(), "b = ", b, " c = ", c);

          auto omega_a = freqs[b + c];
          auto omega_b = freqs[b];
          auto omega_c = freqs[c];

          auto restartA = restart_paths[b + c];
          auto restartB = restart_paths[b];
          auto restartC = restart_paths[c];

          std::array<double, 3> omegas{omega_a, omega_b, omega_c};
          std::array<path, 3> restarts{restartA.replace_extension(""), restartB.replace_extension(""),
                                       restartC.replace_extension("")};

          quad_calculation.set_x_data(world, omegas, restarts);
          auto [beta, beta_directions] = quad_calculation.compute_beta_v2(world, omega_b, omega_c);

          if (world.rank() == 0) {
            ::print("Beta values for omega_A", " = -(omega_", b, " + omega_", c, ") = -", omega_a, " = (", omega_b,
                    " + ", omega_c, ")");
            {
              for (int i = 0; i < beta_directions.size(); i++) {
                std::cout << std::fixed << std::setprecision(5) << "i = " << i + 1 << ", beta[" << beta_directions[i]
                          << "]" << " = " << beta[i] << std::endl;
              }
            }
            append_to_beta_json({-1.0 * omega_a, omega_b, omega_c}, beta_directions, beta, beta_json);
            std::ofstream outfile("beta.json");
            if (outfile.is_open()) {
              outfile << beta_json.dump(4);
              outfile.close();
            }
          }
        }
      }
    } catch (Response_Convergence_Error& e) {
      if (world.rank() == 0) {
        ::print(
            "First order response calculations haven't been run and "
            "can't be run");
        ::print("Quadratic response calculations can't be run");
      }
    }
  }
  static void append_to_beta_json(const std::array<double, 3>& omega, const std::vector<std::string>& beta_directions,
                                  const Tensor<double>& beta, nlohmann::ordered_json& beta_json) {
    auto num_unique_elements = beta_directions.size();
    for (int i = 0; i < num_unique_elements; i++) {

      auto ijk = beta_directions[i];
      auto beta_value = beta[i];
      auto A = ijk[0];
      auto B = ijk[1];
      auto C = ijk[2];

      beta_json["Afreq"].push_back(omega[0]);
      beta_json["Bfreq"].push_back(omega[1]);
      beta_json["Cfreq"].push_back(omega[2]);

      beta_json["A"].push_back(std::string(1, A));
      beta_json["B"].push_back(std::string(1, B));
      beta_json["C"].push_back(std::string(1, C));
      beta_json["Beta"].push_back(beta_value);
    }
  }

  void add_beta_i_to_json(nlohmann::ordered_json& beta_i, nlohmann::ordered_json& beta_json) {
    print(beta_json.dump(4));

    for (auto& [key, value] : beta_i.items()) {
      beta_json[key].insert(beta_json[key].end(), value.begin(), value.end());
    }
  }
};

/**/
/*class WriteResponseVTKOutputStrategy : public CalculationStrategy {*/
/**/
/*  ResponseParameters parameters;*/
/*  std::string op;*/
/*  std::string name;*/
/**/
/* public:*/
/*  explicit WriteResponseVTKOutputStrategy(const ResponseParameters& params) : parameters(params){};*/
/*  json calcPaths(const path& root) override { return {}; }*/
/*  void runCalculation(World& world, const json& paths) override {*/
/**/
/*    auto& moldft_paths = paths["moldft"];*/
/*    auto moldft_restart = moldft_paths["restart"].get<std::string>();*/
/*    moldft_restart = std::string(path(moldft_restart).replace_extension(""));*/
/**/
/*    auto& response_paths = paths["response"];*/
/*    auto& calc_paths = response_paths["calculation"];*/
/**/
/*    auto restart_paths = response_paths["restart"].get<std::vector<std::string>>();*/
/*    auto output_paths = response_paths["output"].get<std::vector<std::string>>();*/
/*    auto alpha_path = response_paths["properties"]["alpha"].get<std::string>();*/
/*    auto freqs = response_paths["frequencies"].get<std::vector<double>>();*/
/*    auto num_freqs = freqs.size();*/
/**/
/*    if (world.rank() == 0) {*/
/*      print("Running VTK output for response calculations");*/
/*      print("Number of frequencies: ", num_freqs);*/
/*      print("Frequencies: ", freqs);*/
/*      print("Output paths: ", output_paths);*/
/*      print("Restart paths: ", restart_paths);*/
/*    }*/
/**/
/*    for (size_t i = 0; i < num_freqs; i++) {*/
/*      auto freq_i = freqs[i];*/
/**/
/*      if (!std::filesystem::exists(calc_paths[i])) {*/
/*        std::filesystem::create_directory(calc_paths[i]);*/
/*      }*/
/*      std::filesystem::current_path(calc_paths[i]);*/
/*      print("current path: ", std::filesystem::current_path());*/
/*      print("calc path: ", calc_paths[i]);*/
/*      print("freq: ", freq_i);*/
/**/
/*      path restart_file_i = restart_paths[i];*/
/*      path response_base_i = output_paths[i];*/
/**/
/*      double last_protocol;*/
/*      bool converged = false;*/
/**/
/*      if (world.rank() == 0) {*/
/*        std::ifstream ifs(response_base_i);*/
/*        json response_base;*/
/*        ifs >> response_base;*/
/*        last_protocol = *response_base["parameters"]["protocol"].get<std::vector<double>>().end();*/
/*        converged = response_base["converged"].get<bool>();*/
/*        print("Thresh of restart data: ", last_protocol);*/
/*        print("Converged: ", converged);*/
/*      }*/
/**/
/*      world.gop.broadcast(converged, 0);*/
/*      world.gop.broadcast(last_protocol, 0);*/
/**/
/*      ResponseParameters r_params = parameters;*/
/*      r_params.set_user_defined_value("omega", freq_i);*/
/*      r_params.set_user_defined_value("archive", moldft_restart);*/
/*      r_params.set_user_defined_value("restart", true);*/
/*      std::string restart_file_string = restart_file_i.filename().stem();*/
/*      r_params.set_user_defined_value("restart_file", restart_file_string);*/
/*      if (converged) {*/
/**/
/*        GroundStateCalculation ground_calculation{world, r_params.archive()};*/
/*        Molecule molecule = ground_calculation.molecule();*/
/*        r_params.set_ground_state_calculation_data(ground_calculation);*/
/*        r_params.set_derived_values(world, molecule);*/
/*        CalcParams calc_params = {ground_calculation, molecule, r_params};*/
/**/
/*        RHS_Generator rhs_generator;*/
/*        if (op == "dipole") {*/
/*          rhs_generator = dipole_generator;*/
/*        } else {*/
/*          rhs_generator = nuclear_generator;*/
/*        }*/
/**/
/*        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));*/
/*        FrequencyResponse calc(world, calc_params, freq_i, rhs_generator);*/
/*        if (world.rank() == 0) {*/
/*          print("\n\n");*/
/*          print(*/
/*              " MADNESS Time-Dependent Density Functional Theory Response "*/
/*              "Program");*/
/*          print(" ----------------------------------------------------------\n");*/
/*          // put the response parameters in a j_molrespone json object*/
/*        }*/
/*        std::string plot_name = "response_" + std::to_string(freq_i);*/
/*        calc.write_vtk(world, 100, parameters.L() / 20, plot_name);*/
/**/
/*        world.gop.fence();*/
/*        // Then we restart from the previous file instead*/
/*      } else {*/
/*      }*/
/*    }*/
/*  }*/
/*};*/
/** /*/

// CalcManager class
// Takes path manager and creates json of paths
// Generates the paths for the calculations in a dry run
// Then runs the calculations
class CalcManager {
 private:
  CompositeCalculationStrategy strategies;
  path root;

 public:
  void addStrategy(std::unique_ptr<CalculationStrategy> newStrategy) { strategies.addStrategy(std::move(newStrategy)); }
  json get_path_json() { return strategies.calcPaths(root); }

  json runCalculations(World& world) {
    json paths = strategies.calcPaths(root);
    // output the paths to disk
    if (world.rank() == 0) {
      std::ofstream ofs("paths.json");
      ofs << paths.dump(4);
      ofs.close();
    }

    if (world.rank() == 0) {
      print(paths.dump(4));
    }

    strategies.compute(world, paths);
    return paths;
  }

  explicit CalcManager(path base_root = std::filesystem::current_path()) : root(std::move(base_root)) {
    strategies = CompositeCalculationStrategy{};
  };
};

#endif
