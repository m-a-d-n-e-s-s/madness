//
// Created by adrianhurtado on 1/1/22.
//
#define CATCH_CONFIG_RUNNER

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
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
    struct stat buffer{};
    size_t rc = stat(input_name, &buffer);
    return (rc == 0);
}

#endif

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    std::cout.precision(6);
    std::string filename = "response.in";
    auto calc_params = initialize_calc_params(world, filename);

    try {
        ExcitedResponse calc(world, calc_params);
        // Warm and fuzzy for the user
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
        }
        // Come up with an initial OK data map
        // set protocol to the first
        calc.solve(world);
        calc.output_json();
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}
