//
// Created by adrianhurtado on 1/1/22.
//
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

    // step 1 is to read molecule from molecule file or
    // to read in the molecule directory from geometry input
    Molecule molecule = Molecule();
    molecule.read(world, "molecule.in");

    std::string filename = "response.in";

    commandlineparser parser(argc, argv);

    if (parser.key_exists("help")) {
        FrequencyResponse::help();
    } else if (parser.key_exists("print_parameters")) {
        FrequencyResponse::print_parameters();
    } else {
        molresponse::start_timer(world);

    }



    const std::string molecule_name{argv[1]};
    const std::string xc{argv[2]};
    const std::string op{argv[3]};
    const std::string precision{argv[4]};
    const std::string static_calc{argv[5]};
    if (precision != "high" && precision != "low" && precision != "super") {
        if (world.rank() == 0) { std::cout << "Set precision to low high super" << std::endl; }
        return 1;
    }
    auto schema = runSchema(world, xc);

    // read in


}

