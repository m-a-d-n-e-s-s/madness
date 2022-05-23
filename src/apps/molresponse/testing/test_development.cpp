//
// Created by adrianhurtado on 1/1/22.
//

#include <fstream>

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "response_parameters.h"
#include "runners.hpp"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"


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
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    // operators include (excited-state)
    auto schema = runSchema(xc);
    auto mol_path = addPath(schema.molecule_path, molecule_name);
    auto m_schema = moldftSchema(molecule_name, xc, schema);
    auto excited_schema = excitedSchema(schema, m_schema);
    excited_schema.print();
    ResponseParameters r_params{};
    set_excited_parameters(r_params, excited_schema.xc, excited_schema.num_states);
    create_excited_paths(r_params, excited_schema, false);
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
