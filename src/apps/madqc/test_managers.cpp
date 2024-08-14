// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include <filesystem>
#include "calc_manager.hpp"
#include "madness/external/catch/catch.hpp"
#include "parameter_manager.hpp"

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

TEST_CASE("MOLDFT Calculation") {

  World& world = World::get_default();

  ParameterManager param_manager;

  param_manager = ParameterManager(world, {"input.json"});
  auto params = param_manager.get_moldft_params();
  auto molecule = param_manager.get_molecule();

  // this is where I we create our calculation
  CalcManager calc_manager;
  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  calc_manager.addStrategy(std::move(moldft_calc));
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
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
  ResponseInput r_input = std::make_tuple(perturbation, xc, freq_range);
  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  auto response_calc = std::make_unique<ResponseCalculationStrategy>(response_params, r_input);

  calc_manager.addStrategy(std::move(moldft_calc));
  calc_manager.addStrategy(std::move(response_calc));

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

  CalcManager calc_manager;

  auto set_freqs = [&]() {
    vector<double> freqs_copy = freq_range;
    auto num_freqs = freq_range.size();
    auto compare_freqs = [](double x, double y) {
      return std::abs(x - y) < 1e-3;
    };

    for (int i = 0; i < num_freqs; i++) {    // for i=0:n-1
      for (int j = i; j < num_freqs; j++) {  // for j = i  omega_3=-(omega_1+omega_2)
        auto omega_1 = freq_range[i];
        auto omega_2 = freq_range[j];
        auto omega_3 = omega_1 + omega_2;

        // if you can find omega_3 in freq_copy skip
        if (omega_2 == 0.0)
          continue;
        if (std::find_if(freqs_copy.begin(), freqs_copy.end(), [&](double x) { return compare_freqs(x, omega_3); }) !=
            freqs_copy.end()) {
          continue;
        } else {
          freqs_copy.push_back(omega_3);
        }
      }
    }
    return freqs_copy;
  };

  freq_range = set_freqs();
  print("Frequency Range: ", freq_range);

  ResponseInput h_input = std::make_tuple(perturbation, xc, freq_range);

  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  auto response_calc = std::make_unique<ResponseCalculationStrategy>(response_params, h_input);
  auto hyper_calc = std::make_unique<HyperPolarizabilityCalcStrategy>(response_params, h_input);

  calc_manager.addStrategy(std::move(moldft_calc));
  calc_manager.addStrategy(std::move(response_calc));
  calc_manager.addStrategy(std::move(hyper_calc));

  // get cwd
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
  std::filesystem::current_path(cwd);
}

TEST_CASE("Output Response VTK") {

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

  CalcManager calc_manager;

  auto moldft_calc = std::make_unique<MoldftCalculationStrategy>(params, molecule);
  auto response_calc = std::make_unique<ResponseCalculationStrategy>(response_params, r_input);
  auto vtk_plot = std::make_unique<WriteResponseVTKOutputStrategy>(response_params);

  calc_manager.addStrategy(std::move(moldft_calc));
  calc_manager.addStrategy(std::move(response_calc));
  calc_manager.addStrategy(std::move(vtk_plot));

  // get cwd
  path cwd = std::filesystem::current_path();
  calc_manager.runCalculations(world, cwd);
  std::filesystem::current_path(cwd);
}
