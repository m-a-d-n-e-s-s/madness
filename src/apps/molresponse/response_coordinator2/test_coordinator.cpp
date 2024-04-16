//
// Created by adrianhurtado on 1/1/22.
//
#include "CalculationParameters.h"

#define CATCH_CONFIG_RUNNER

#include "FrequencyResponse.hpp"
#include "coordinator.hpp"
#include "madness/external/catch/catch.hpp"
#include "string"
#include <filesystem>


using path = std::filesystem::path;

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

}


std::pair<CalculationParameters, ResponseParameters> get_params(const std::string &json_input_path) {

    std::ifstream input_stream(json_input_path);
    json input_json = json::parse(input_stream);

    json moldft_json = input_json["moldft"];
    json molresponse_json = input_json["molresponse"];

    CalculationParameters moldft_params;
    moldft_params.from_json(moldft_json);
    moldft_params.print();

    ResponseParameters molresponse_params;
    molresponse_params.from_json(molresponse_json);
    molresponse_params.print();
    return std::make_pair(moldft_params, molresponse_params);

}

TEST_CASE("Basic Response Manager Test") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path json_input_path("resources/inputs/input.json");
    path molecule_path("resources/molecules/H2O.mol");
    // print full path
    // Read in json

    auto [moldft_params, molresponse_params] = get_params(json_input_path);
    auto response_manager = ResponseCalcManager(world, molecule_path, moldft_params, molresponse_params);
    response_manager.print();
    // The json is converted into a temporary getKW file which is then read by the parser.
    // Now we need to write a function
}

TEST_CASE("Param Hash Generation Test") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path mol1_path("resources/molecules/H2O.mol");
    path mol2_path("resources/molecules/Be.mol");
    path json_input_path("resources/inputs/input.json");

    auto [moldft_params, molresponse_params] = get_params(json_input_path);


    auto response_manager = ResponseCalcManager(world, mol1_path, moldft_params, molresponse_params);
    auto response_manager2 = ResponseCalcManager(world, mol2_path, moldft_params, molresponse_params);

    // Test that the param directories are the same
    CHECK(response_manager.get_param_path()== response_manager2.get_param_path());
    // Test that the molecule directories are different
    CHECK(response_manager.get_moldft_path()!= response_manager2.get_moldft_path());

}

TEST_CASE("Test Manager Run MOLDFT"){ 
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path mol1_path("resources/molecules/He.mol");
    path json_input_path("resources/inputs/input.json");

    auto [moldft_params, molresponse_params] = get_params(json_input_path);


    auto response_manager = ResponseCalcManager(world, mol1_path, moldft_params, molresponse_params);

    response_manager.run_moldft(world,true);
}
