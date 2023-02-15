//
// Created by adrianhurtado on 1/1/22.
//
#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "/chem/SCF.h"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "runners.hpp"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static inline int file_exists(const char *input_name) {
    struct stat buffer {};
    size_t rc = stat(input_name, &buffer);
    return (rc == 0);
}

#endif

using path = std::filesystem::path;

using namespace madness;


int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);

    std::cout.precision(6);

    if (argc != 3) {

        std::cout << "Wrong number of inputs" << std::endl;
        return 1;
    }

    const std::string molecule_name{argv[1]};
    const std::string xc{argv[2]};
    const std::string op = "excited-state";


    auto schema = runSchema(xc);
    auto mol_path = addPath(schema.molecule_path, molecule_name);

    auto m_schema = moldftSchema(molecule_name, xc, schema);
    m_schema.print();
    //moldft(world, m_schema, false);
    auto excited_schema = excitedSchema(schema, m_schema);
    excited_schema.print();
    //bool success = runExcited(world, excited_schema, true);

    // Set the response parameters
    ResponseParameters r_params{};

    set_excited_parameters(r_params, excited_schema.xc, excited_schema.num_states);
    create_excited_paths(r_params, excited_schema, true);
    std::filesystem::current_path(excited_schema.excited_state_run_path);
    set_and_write_restart_excited_parameters(r_params, excited_schema, true);

    auto calc_params = initialize_calc_params(world, "response.in");

    auto molecule = calc_params.molecule;
    std::string all_pg[] = {"c1", "cs", "c2", "ci", "c2v", "c2h", "d2", "d2h"};

    auto g_orbitals = calc_params.ground_calculation.orbitals();


    double error = 0.0;

    auto pg = all_pg[4];//"c2v

    print("point group", pg);
    projector_irrep proj(pg);

    for (const std::string &irrep: proj.get_all_irreps()) {// loop over all irreps
        print(" irrep", irrep);
        proj.set_irrep(irrep);
        real_function_3d f1 = proj(g_orbitals[0])[0];// result is the first element of result vector

        charactertable table = proj.get_table();
        for (int i = 0; i < table.order_; ++i) {// loop over all symmetry operations
            const pg_operator syop = table.operators_[i];
            const int character = table.irreps_[irrep][i];
            double n1 = (f1 - character * syop(f1)).norm2();
            print("  operator, character, norm ", syop.name(), character, n1);
            error += n1;
        }
    }

    /*
        ExcitedResponse calc(world, calc_params);
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS Time-Dependent Density Functional Theory Response "
                  "Program");
            print(" ----------------------------------------------------------\n");
            print("\n");
            calc_params.molecule.print();
            print("\n");
            calc_params.response_parameters.print("response");
            // put the response parameters in a j_molrespone json object
            calc_params.response_parameters.to_json(calc.j_molresponse);
        // set protocol to the first
        calc.solve(world);
        calc.output_json();
        return true;
         */

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}