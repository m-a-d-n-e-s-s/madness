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

TEST_CASE("Test Parameters Class", "[Parameters]") {
  World& world = World::get_default();

  auto params_json = ParameterManager(world, {"input.json"});

  if (world.rank() == 0) {
    print("Checking Parameters");
    REQUIRE(params_json.get_run_moldft() == true);
    REQUIRE(params_json.get_molresponse_params().quadratic() == true);
    REQUIRE(params_json.get_molresponse_params().freq_range() ==
            std::vector<double>{0.0, 0.056, 0.1});
  }
  world.gop.fence();

  auto param_input = ParameterManager(world, path{"input"});

  if (world.rank() == 0) {
    print("Checking Parameters");
    REQUIRE(param_input.get_run_moldft() == true);
    REQUIRE(param_input.get_molresponse_params().quadratic() == true);
    REQUIRE(param_input.get_molresponse_params().freq_range() ==
            std::vector<double>{0.0, 0.056, 0.1});
  }
  // check if parameters are the same

  if (world.rank() == 0) {
    bool check_equal = params_json == param_input;
    if (!check_equal) {
      print("Parameters are not equal");
      print("params:");
      print(params_json.get_input_json());
      print("param_input:");
      print(param_input.get_input_json());

      auto diff = json::diff(params_json.get_input_json(),
                             param_input.get_input_json());
      print("diff:");
      print(diff.dump(2));

    } else {
      print("Parameters are equal");
    }
    CHECK(check_equal);
  }
}
