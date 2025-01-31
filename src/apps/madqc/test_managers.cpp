// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include "calc_factory.hpp"
#include "calc_manager.hpp"
#include "madness/external/catch/catch.hpp"
#include "parameter_manager.hpp"
#include <filesystem>

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
TEST_CASE("Constructing an Optimization Calculation") {

  World& world = World::get_default();

  ParameterManager params;
  // Initialize the necessary components
  path input_file("input.json");
  if (world.rank() == 0) {
    print("Input file found");
    print("Parsing Command Line");
  }
  params = ParameterManager(world, input_file);

  auto task_params = params.get_task_params();

  auto method = task_params.method;
  auto driver = task_params.driver;
  auto properties = task_params.properties;
  auto molecule = params.get_molecule();

  if (world.rank() == 0) {
    task_params.print();
  }
  std::unique_ptr<CalculationDriver> calc_manager;
  calc_manager = createEnergyDriver(world, method, params, properties);

  path cwd = std::filesystem::current_path();
  // This should create a new path_manager with root = optimize1
  calc_manager->setRoot(cwd / "optimize1");

  auto& opt_params = params.get_optimization_params();

  MolOpt opt(opt_params.get_maxiter(),             // geoometry max iter
             0.1,                                  // geometry step size
             opt_params.get_value_precision(),     // value precision
             opt_params.get_geometry_tolerence(),  // geometry tolerance
             1e-3,                                 // XTOL
             1e-5,                                 // EPREC
             opt_params.get_gradient_precision(),  // gradient precision
             (world.rank() == 0) ? 1 : 0,          // print_level
             opt_params.get_algopt());             // algorithm options

  auto new_molecule = opt.optimize(molecule, *calc_manager);
}

TEST_CASE("Copy driver with new name... or nested") {

  World& world = World::get_default();

  ParameterManager params;
  // Initialize the necessary components
  path input_file("input.json");
  if (world.rank() == 0) {
    print("Input file found");
    print("Parsing Command Line");
  }
  params = ParameterManager(world, input_file);

  auto task_params = params.get_task_params();

  auto method = task_params.method;
  auto driver = task_params.driver;
  auto properties = task_params.properties;
  auto molecule = params.get_molecule();

  if (world.rank() == 0) {
    task_params.print();
  }
  auto calc_manager = createEnergyDriver(world, method, params, properties);
  calc_manager->setRoot("optimize1");
  auto new_manager = calc_manager->clone();
  path cwd = std::filesystem::current_path();
  calc_manager->runCalculations(molecule.get_all_coords().flat());

  std::filesystem::current_path(cwd);
  new_manager->setRoot("optimize2");
  new_manager->runCalculations(molecule.get_all_coords().flat());
  std::filesystem::current_path(cwd);

  // create a driver with moldft
}

TEST_CASE("Hyperpolarizability Calculation") {}

/*

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
  calc_manager.runCalculations(world);
  std::filesystem::current_path(cwd);
}
*/
