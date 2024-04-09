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
    double frequency = 0.0;
    std::string property = "dipole";

    auto calc_params = initialize_calc_params(world, std::string(filename));
    print(calc_params.response_parameters.restart_file());

    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
    calc.load(world, calc_params.response_parameters.save_file());
    auto x = calc.get_chi();
    print(x.num_states());
    print(x.num_orbitals());
    //calc.check_k(world,1e-4,7);

    auto PQ = calc.generator(world, calc);
    PQ.truncate();


    int i = 0;
    for (auto &ai: x.active) { ai = i++; }
    print(x.active);

    auto pn = PQ.norm2s();
    auto xn = x.norm2s();

    auto add = x + PQ;
    print("x: ", xn);
    print("p: ", pn);
    print("x+p", add.norm2s());

    vector<bool> converged{1, 0, 0};
    int b = 0;
    x.active.remove_if([&](auto x) { return converged[b++]; });
    b = 0;

    print("x active", x.active);
    auto add_1 = x + PQ;
    print("x+p after remove", add_1.norm2s());
    x.reset_active();
    auto add_2 = x + PQ;
    auto norm2 = add_2.norm2s();
    print("x+p after reset", add_2.norm2s());
    CHECK(add_2.norm2s()[0] != add_1.norm2s()[0]);

    b = 0;
    x.active.remove_if([&](auto x) { return converged[b++]; });
    x += PQ + x;
    print("x+=pQ ", x.norm2s());
    x = x - PQ;
    print("x=x-PQ ", x.norm2s());

    auto base = {1e-3, 1e-5, 1e-4};
    vector<std::pair<double, double>> ri(base.size());

    std::transform(base.begin(), base.end(), ri.begin(), [](auto vi) {
        return std::pair<double, double>{vi, vi * 5};
    });

    std::pair<double, double> threshold{1e-3, 5e-4};

    std::transform(ri.begin(), ri.end(), converged.begin(), [&](auto &resi) {
        if (resi.first < threshold.first && resi.second < threshold.second) {
            return true;
        } else {
            return false;
        }
    });
    x.reset_active();
    print(converged);
    b = 0;
    x.active.remove_if([&](auto x) { return converged[b++]; });
    print(x.norm2s());
    print((x + PQ).norm2s());
}

TEST_CASE("Testing xspace zero functions") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();

    std::cout.precision(6);
    std::string filename = "response.in";
    double frequency = 0.0;
    std::string property = "dipole";

    auto calc_params = initialize_calc_params(world, std::string(filename));
    print(calc_params.response_parameters.restart_file());

    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
    calc.load(world, calc_params.response_parameters.save_file());
    auto x = calc.get_chi();


    auto y = x.copy();
    auto z = y;
    auto zeros=X_space::zero_functions(world,x.num_states(),x.num_states());

    print(y.norm2s());
    print(x.norm2s());
    print(z.norm2s());
    print(zeros.norm2s());

    auto xy=calc.response_context.inner(x,y);
    print(xy);



}
