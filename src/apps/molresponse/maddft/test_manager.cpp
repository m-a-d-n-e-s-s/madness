// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include <filesystem>
#include "maddft/response_manager.hpp"
#include "madness/external/catch/catch.hpp"
#include "molresponse/FrequencyResponse.hpp"

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

TEST_CASE("Testing Parameters Class", "[Parameters]") {
  World& world = World::get_default();

  path input_file("input.json");

  if (world.rank() == 0) {
    print("Input file found");
    print("Parsing Command Line");
  }

  auto params = ParameterManager(world, input_file);

  if (world.rank() == 0) {
    print("Checking Parameters");
    REQUIRE(params.get_run_moldft() == true);
    REQUIRE(params.get_molresponse_params().quadratic() == true);
    REQUIRE(params.get_molresponse_params().freq_range() ==
            std::vector<double>{0.0, 0.056, 0.1});
  }
  world.gop.fence();
}
