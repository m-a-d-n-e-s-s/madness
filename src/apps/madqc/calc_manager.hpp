#ifndef CALC_MANAGER_HPP
#define CALC_MANAGER_HPP


#include <apps/molresponse/response_parameters.h>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <utility>
#include <vector>
#include "path_manager.hpp"
#include "utils.hpp"

using json = nlohmann::json;
using path = std::filesystem::path;

class CalculationStrategy {
 public:
  virtual void runCalculation(World& world, const json& paths) = 0;
  virtual ~CalculationStrategy() = default;
};

class CompositeCalculationStrategy : public CalculationStrategy {
 private:
  std::vector<std::unique_ptr<CalculationStrategy>> strategies;

 public:
  void addStrategy(std::unique_ptr<CalculationStrategy> strategy) {
    strategies.push_back(std::move(strategy));
  }

  void runCalculation(World& world, const json& paths) override {
    for (const auto& strategy : strategies) {
      strategy->runCalculation(world, paths);
    }
  }
};

class MoldftCalculationStrategy : public CalculationStrategy {

  CalculationParameters parameters;
  json input_json;
  Molecule molecule;

 public:
  void runCalculation(World& world, const json& moldft_paths) override {

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

    if (std::filesystem::exists(moldft_restart) &&
        std::filesystem::exists(moldft_calc_info_path)) {
      // if both exist, read the calc_info json
      std::ifstream ifs(moldft_calc_info_path);

      auto moldft_calc_info = json::parse(ifs);
      if (world.rank() == 0) {
        std::cout << "time: " << moldft_calc_info["time"] << std::endl;
        std::cout << "MOLDFT return energy: "
                  << moldft_calc_info["return_energy"] << std::endl;
      }
    } else {
      // if params are different run and if restart exists and if im asking to
      if (std::filesystem::exists(moldft_restart)) {
        param1.set_user_defined_value<bool>("restart", true);
      }
      world.gop.fence();
      if (world.rank() == 0) {

        auto moldft_input_json = parameters.to_json_if_precedence("defined");
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
      // vama
      calc.set_protocol<3>(world, calc.param.protocol()[0]);
      // calc.set_protocol<3>(world, 1e-4);
      world.gop.fence();
      MolecularEnergy ME(world, calc);
      world.gop.fence();
      // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
      ME.value(calc.molecule.get_all_coords().flat());  // ugh!
      world.gop.fence();
      const real_function_3d rho =
          2.0 * calc.make_density(world, calc.aocc, calc.amo);
      auto dipole_t = calc.dipole(world, rho);
      std::map<std::string, double> results;
      results["scf_energy"] = calc.current_energy;
      world.gop.fence();
      if (world.rank() == 0) {
        calc.output_scf_info_schema(results, dipole_t);
        ME.output_calc_info_schema();
      }
    }

    // Add actual calculation logic here
  }

  MoldftCalculationStrategy(const CalculationParameters& params, Molecule mol)
      : parameters(params), molecule(std::move(mol)){};
};

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

class ResponseCalculationStrategy : public CalculationStrategy {
 public:
  void runCalculation(World& world, const json& paths) override {
    std::cout << "Running Response Calculation\n";
    for (const auto& dir : paths["response_calc_dirs"]) {
      std::cout << "Calc Directory: " << dir << "\n";
    }
    for (const auto& restart : paths["response_restarts"]) {
      std::cout << "Restart Path: " << restart << "\n";
    }
    for (const auto& outfile : paths["response_outfiles"]) {
      std::cout << "Output File: " << outfile << "\n";
    }
    // Add actual calculation logic here
  }
};

class CalcManager {
 private:
  std::unique_ptr<CalculationStrategy> strategy;
  PathManager pathManager;

 public:
  void setCalculationStrategy(
      std::unique_ptr<CalculationStrategy> newStrategy) {
    strategy = std::move(newStrategy);
  }

  void addPathStrategy(std::unique_ptr<PathStrategy> pathStrategy) {
    pathManager.addStrategy(std::move(pathStrategy));
  }

  void runCalculations(World& world, const path& root) {
    json paths = pathManager.generateCalcPaths(root);

    json allPaths;
    allPaths.update(paths);

    strategy->runCalculation(world, allPaths);
  }
};
#endif
