//
// Created by adrianhurtado on 1/1/22.
//


#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "madness/tensor/tensor_json.hpp"
#include "response_functions.h"
#include "response_parameters.h"
#include "runners.hpp"


using path = std::filesystem::path;

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    auto molecule_path = root;
    molecule_path += "/molecules";

    const std::string molecule_name = "Be";
    const std::string xc = "hf";
    const std::string op = "excited-state";
    const std::string precision = "low";
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    // operators include (excited-state)
    auto schema = ResponseCalcManager(world, xc);
    auto mol_path = addPath(schema.molecules, molecule_name);
    auto m_schema = moldftSchema(world, molecule_name, xc, schema);
    auto excited_schema = excitedSchema(schema, m_schema);
    excited_schema.print();
    ResponseParameters r_params{};
    set_excited_parameters(world, r_params, excited_schema.xc,
                           excited_schema.num_states, precision);
    create_excited_paths(excited_schema);
    std::filesystem::current_path(excited_schema.excited_state_run_path);
    set_and_write_restart_excited_parameters(r_params, excited_schema, false);
    auto xc_path = create_xc_path_and_directory(root, xc);

    auto calc_params = initialize_calc_params(world, "response.in");

    ExcitedResponse calc(world, calc_params);
    ExcitedTester tester{world, calc, 1e-4};
    auto guess = tester.test_ao_guess(world, calc);


    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}
