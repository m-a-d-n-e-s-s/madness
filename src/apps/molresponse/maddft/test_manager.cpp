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

TEST_CASE("Basic Input", "[Parameters]") {
  World& world = World::get_default();

  auto check_parameters = [&](auto& params) {
    if (world.rank() == 0) {
      print("Checking Parameters");
    }
    REQUIRE(params.get_run_moldft() == true);
    REQUIRE(params.get_molresponse_params().quadratic() == true);
    REQUIRE(params.get_molresponse_params().freq_range() ==
            std::vector<double>{0.0, 0.056, 0.1});
  };
  auto check_params_equal = [&](auto& params_json, auto& param_input) {
    bool check_equal = params_json == param_input;
    if (!check_equal) {
      auto diff = json::diff(params_json.get_input_json(),
                             param_input.get_input_json());
      print("diff:");
      print(diff.dump(2));
    }
    REQUIRE(check_equal);
  };

  auto params_json = ParameterManager(world, {"input.json"});
  check_parameters(params_json);

  auto param_input = ParameterManager(world, path{"input"});
  check_parameters(param_input);
  check_params_equal(params_json, param_input);
}

TEST_CASE("MOLDFT ONLY", "[Parameters]") {
  World& world = World::get_default();

  auto check_moldft = [&](auto& params) {
    REQUIRE(params.get_run_moldft() == true);
    REQUIRE(params.get_run_response() == false);
  };

  auto check_params_equal = [&](auto& params_json, auto& param_input) {
    bool check_equal = params_json == param_input;
    if (!check_equal) {
      auto diff = json::diff(params_json.get_input_json(),
                             param_input.get_input_json());
      print("diff:");
      print(diff.dump(2));
    }
    REQUIRE(check_equal);
  };

  auto params_json = ParameterManager(world, {"moldft_input.json"});
  check_moldft(params_json);

  auto param_input = ParameterManager(world, path{"moldft_input"});
  check_moldft(param_input);

  // check if parameters are the same
  check_params_equal(params_json, param_input);
}
