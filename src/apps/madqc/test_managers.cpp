// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include <filesystem>
#include "madness/external/catch/catch.hpp"
#include "path_manager.hpp"
#include "calc_manager.hpp"

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

  path_manager.addStrategy(
      std::make_unique<ResponsePathStrategy>("response", input));
  json paths = path_manager.generateCalcPaths("root");
  // output the paths
  std::cout << paths.dump(4) << std::endl;

  ResponseInput input_2 = std::make_tuple("nuclear", "xc", freq_range);
  PathManager response_manager2;
  response_manager2.addStrategy(std::make_unique<MoldftPathStrategy>());
  response_manager2.addStrategy(
      std::make_unique<ResponsePathStrategy>("response", input));

  json paths2 = response_manager2.generateCalcPaths("root");
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

  path_manager.addStrategy(std::make_unique<ExcitedStatePathStrategy>(
      "excited_states", xc, nums_states));

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


TEST_CASE("MOLDFT Calculation"){

  World& world = World::get_default();

  CalcManager calc_manager;






}









