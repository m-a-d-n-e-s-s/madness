//
// Created by adrianhurtado on 1/1/22.
//
#include "CalculationParameters.h"
#define CATCH_CONFIG_RUNNER

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "coordinator.hpp"
#include "madness/external/catch/catch.hpp"
#include "madness/tensor/tensor_json.hpp"
#include "response_functions.h"
#include "string"
#include "x_space.h"
#include <filesystem>


using path = std::filesystem::path;

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}


TEST_CASE("Hash Generation Test") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    Molecule molecule = Molecule();
    path json_input_path("resources/inputs/input.json");
    print("Full path to json input file: ", json_input_path.string());
    path molecule_path("resources/molecules/H2O.mol");
    // print full path
    // Read in json

    std::ifstream input_stream(json_input_path);
    json input_json = json::parse(input_stream);
    print("Input json read from file: ");
    print(input_json.dump(4));

    json moldft_json = input_json["moldft"];
    json molresponse_json = input_json["molresponse"];
    commandlineparser parser;

    CalculationParameters moldft_params;
    moldft_params.from_json(moldft_json);
    moldft_params.print();

    ResponseParameters molresponse_params;
    molresponse_params.from_json(molresponse_json);
    molresponse_params.print();


    auto response_manager = ResponseCalcManager(world, molecule_path, moldft_params, molresponse_params);
    response_manager.print();


    if (world.rank() == 0) print("input filename: ", parser.value("input"));


    // The json is converted into a temporary getKW file which is then read by the parser.
    // Now we need to write a function
}
