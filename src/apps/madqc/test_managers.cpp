// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include <filesystem>
#include "calc_manager.hpp"
#include "madness/external/catch/catch.hpp"
#include "parameter_manager.hpp"
#include "path_manager.hpp"

using path = std::filesystem::path;

int main(int argc, char* argv[]) {
  int result = 0;
  {
    World& world = madness::initialize(argc, argv);

    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }
    world.gop.fence();

    finalize();
  }
  return result;
}

TEST_CASE("MoldftStrategy", "PathStrategy") {
  World& world = World::get_default();

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>("moldft_1"));
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>("moldft_2"));

  json paths = path_manager.generateCalcPaths("root");
  std::cout << paths.dump(4) << std::endl;
};

TEST_CASE("ResponsePathStrategy", "PathStrategy") {
  World& world = World::get_default();

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());

  std::vector<double> freq_range = {0.0, 0.056, 0.1};
  std::string perturbation = "dipole";
  std::string xc = "lda";
  ResponseInput input = std::make_tuple(perturbation, xc, freq_range);

  path_manager.addStrategy(std::make_unique<ResponsePathStrategy>(input, "response"));
  json paths = path_manager.generateCalcPaths("root");
  // output the paths
  std::cout << paths.dump(4) << std::endl;

  ResponseInput input_2 = std::make_tuple("nuclear", "xc", freq_range);
  PathManager response_manager2;
  response_manager2.addStrategy(std::make_unique<MoldftPathStrategy>());
  response_manager2.addStrategy(std::make_unique<ResponsePathStrategy>(input, "response"));

  json paths2 = response_manager2.generateCalcPaths("");
  // output the paths
  std::cout << paths2.dump(4) << std::endl;
};

TEST_CASE("ExcitedStatePath", "PathStrategy") {
  World& world = World::get_default();

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());

  int nums_states = 5;
  std::string xc = "lda";

  path_manager.addStrategy(std::make_unique<ExcitedStatePathStrategy>("excited_states", xc, nums_states));

  json paths = path_manager.generateCalcPaths("root");
  std::cout << paths.dump(4) << std::endl;
}

TEST_CASE("MP2PathStrategy", "PathStrategy") {
  World& world = World::get_default();

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());
  path_manager.addStrategy(std::make_unique<MP2PathStrategy>());

  json paths = path_manager.generateCalcPaths("root");
  std::cout << paths.dump(4) << std::endl;
}

TEST_CASE("MOLDFT Calculation") {

  World& world = World::get_default();

  ParameterManager param_manager;

  param_manager = ParameterManager(world, {"input.json"});
  auto params = param_manager.get_moldft_params();
  auto molecule = param_manager.get_molecule();

  // this is where I we create our calculation
  CalcManager calc_manager;
  calc_manager.addPathStrategy(std::make_unique<MoldftPathStrategy>());
  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);

  calc_manager.setCalculationStrategy(std::move(moldft_calc));
  // get cwd
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
  //reset the current path to the original path
  std::filesystem::current_path(cwd);
}

TEST_CASE("Response Calculation") {

  World& world = World::get_default();

  ParameterManager param_manager;
  param_manager = ParameterManager(world, {"input.json"});

  std::string perturbation = "dipole";
  std::string xc = "hf";
  auto response_params = param_manager.get_molresponse_params();
  std::vector<double> freq_range = response_params.freq_range();

  //auto freq_range = response_params.freq_range();
  //std::vector<double> freq_range = {0.0, 0.056, 0.1};
  if (world.rank() == 0) {
    print("Running Response Calculation");
    print("Perturbation: ", perturbation);
    print("XC: ", xc);
    print("Frequency Range: ", freq_range);
  }

  auto params = param_manager.get_moldft_params();
  auto molecule = param_manager.get_molecule();

  // this is where I we create our calculation
  CalcManager calc_manager;
  calc_manager.addPathStrategy(std::make_unique<MoldftPathStrategy>());
  ResponseInput r_input = std::make_tuple(perturbation, xc, freq_range);
  calc_manager.addPathStrategy(std::make_unique<ResponsePathStrategy>(r_input, "response"));

  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  auto response_calc = std::make_unique<ResponseCalculationStrategy>(response_params);

  calc_manager.setCalculationStrategy(std::move(moldft_calc));
  calc_manager.setCalculationStrategy(std::move(response_calc));

  // get cwd
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
  std::filesystem::current_path(cwd);
}

TEST_CASE("Hyperpolarizability Calculation") {

  World& world = World::get_default();

  ParameterManager param_manager;
  param_manager = ParameterManager(world, {"input.json"});

  std::string perturbation = "dipole";
  std::string xc = "hf";
  auto response_params = param_manager.get_molresponse_params();
  std::vector<double> freq_range = response_params.freq_range();

  //auto freq_range = response_params.freq_range();
  //std::vector<double> freq_range = {0.0, 0.056, 0.1};
  if (world.rank() == 0) {
    print("Running Response Calculation");
    print("Perturbation: ", perturbation);
    print("XC: ", xc);
    print("Frequency Range: ", freq_range);
  }
  ResponseInput r_input = std::make_tuple(perturbation, xc, freq_range);
  auto params = param_manager.get_moldft_params();
  auto molecule = param_manager.get_molecule();

  auto beta_path = std::make_unique<HyperPolarizabilityPathStrategy>(r_input);

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());
  path_manager.addStrategy(std::make_unique<HyperPolarizabilityPathStrategy>(r_input));
  path_manager.generateCalcPaths("");

  CalcManager calc_manager;
  calc_manager.addPathStrategy(std::make_unique<MoldftPathStrategy>());
  calc_manager.addPathStrategy(std::make_unique<ResponsePathStrategy>(r_input, "response"));
  calc_manager.addPathStrategy(std::make_unique<HyperPolarizabilityPathStrategy>(r_input));

  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  auto response_calc = std::make_unique<ResponseCalculationStrategy>(response_params);
  auto hyper_calc = std::make_unique<HyperPolarizabilityCalcStrategy>(response_params);

  calc_manager.setCalculationStrategy(std::move(moldft_calc));
  calc_manager.setCalculationStrategy(std::move(response_calc));
  calc_manager.setCalculationStrategy(std::move(hyper_calc));

  // get cwd
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
  std::filesystem::current_path(cwd);

  // this is where I we create our calculation
  /*CalcManager calc_manager;*/
  /*calc_manager.addPathStrategy(std::make_unique<MoldftPathStrategy>());*/
  /*calc_manager.addPathStrategy(*/
  /*    st#d::make_unique<ResponsePathStrategy>("response", r_input));*/
  /**/
  /*auto moldft_calc =*/
  /*    std::make_unique<MoldftCalculationStrategy>(params, molecule);*/
  /**/
  /*auto response_calc =*/
  /*    std::make_unique<ResponseCalculationStrategy>(response_params);*/
  /**/
  /*calc_manager.setCalculationStrategy(std::move(moldft_calc));*/
  /*calc_manager.setCalculationStrategy(std::move(response_calc));*/
  /**/
  /*// get cwd*/
  /*path cwd = std::filesystem::current_path();*/
  /*calc_manager.runCalculations(world, cwd);*/
}
