//
// Created by adrianhurtado on 1/1/22.
//
#define CATCH_CONFIG_RUNNER

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
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

TEST_CASE("Testing New X SPACE") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();

    std::cout.precision(6);
    std::string filename = "response.in";
    std::string load_x = "restart_dipole_hf_0-000000.000000";
    double frequency = 0.0;
    std::string property = "dipole";

    auto calc_params = initialize_calc_params(world, std::string(filename));

    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
    calc.load(world, load_x);
    auto x = calc.get_chi();

    auto m = x.num_states();
    x.active.resize(m);
    int i = 0;
    for (auto &ai: x.active) { ai = i++; }
    print(x.active);


    // A calculation is defined by a molecule, functional, and operator
    // xc inclu
}
