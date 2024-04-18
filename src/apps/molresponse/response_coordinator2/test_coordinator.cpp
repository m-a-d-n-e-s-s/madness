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
#include <utility>


using path = std::filesystem::path;

// Global structure to store command line arguments
struct CommandLineArgs {
    int argc;
    std::vector<std::string> argv;

    CommandLineArgs() : argc(0) {}

    void parse(int _argc, char *_argv[]) {
        argc = _argc;
        argv.clear();
        for (int i = 0; i < argc; ++i) {
            argv.emplace_back(_argv[i]);
        }
    }
};

// Global instance
CommandLineArgs cmdArgs;
// Create function which converts json input to and outputs to a output stream




int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);

    // Parse command line arguments
    cmdArgs.parse(argc, argv);

    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

}


std::pair<CalculationParameters, ResponseParameters> get_params(const std::string &json_input_path) {

    std::ifstream input_stream(json_input_path);
    json input_json = json::parse(input_stream);

    json moldft_json = input_json["dft"];
    json molresponse_json = input_json["response"];
    CalculationParameters moldft_params;


    moldft_params.from_json(moldft_json);
    ResponseParameters molresponse_params;
    molresponse_params.from_json(molresponse_json);
    return std::make_pair(moldft_params, molresponse_params);

}



TEST_CASE("Testing Parameters Class", "[Parameters]") {
    using namespace madness;

    World &world = World::get_default();
    auto testArgs(cmdArgs);

    // commandlineparser parser(2, std::vector{std::string{"binaryname"}, std::string{"input"}});
    commandlineparser parser;
    parser.set_defaults();

    if (world.rank() == 0) print("input filename: ", parser.value("input"));

    ParameterManager params(world, parser);
    ParameterManager params2(world, params.get_input_json());

    CHECK(params.get_input_json() == params2.get_input_json());

}

TEST_CASE("Paramter Reader From input path", "[JSON]") {
    using namespace madness;

    World &world = World::get_default();
    auto testArgs(cmdArgs);


    path path_input("input");
    path path_json("example.json");


    ParameterManager params(world, path_input);
    ParameterManager params2(world, path_json);

    CHECK(params.get_input_json() == params2.get_input_json());

}

TEST_CASE("Define parameters with input and mol file separately") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path input_json("resources/inputs/input.json");
    path mol_input("resources/molecules/H2.mol");


    ParameterManager params(world, input_json, mol_input);
    auto response_manager = ResponseCalcManager(world, params);

    response_manager.run_moldft(world,false);
}

TEST_CASE("Run MOLDFT") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);

    path json_input_path("example.json");

    ParameterManager params(world, json_input_path);
    auto response_manager = ResponseCalcManager(world, params);

    response_manager.run_moldft(world,false);

    // The json is converted into a temporary getKW file which is then read by the parser.
    // Now we need to write a function
}

TEST_CASE("Run MOLDFT + MOLRESPONSE") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();
    std::cout.precision(6);


    path input_json("resources/inputs/freq_input.json");
    path mol_input("resources/molecules/He.mol");


    ParameterManager params(world, input_json, mol_input);
    auto response_manager = ResponseCalcManager(world, params);

    response_manager.run_moldft(world,true);
    response_manager.run_molresponse(world);

    // The json is converted into a temporary getKW file which is then read by the parser.
    // Now we need to write a function
}


