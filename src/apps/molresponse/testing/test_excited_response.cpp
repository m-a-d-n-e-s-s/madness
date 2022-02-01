//
// Created by adrianhurtado on 1/1/22.
//
#define CATCH_CONFIG_RUNNER
#include "apps/external_headers/catch.hpp"
#include "apps/chem/SCF.h"
#include "madness/world/worldmem.h"
#include "apps/external_headers/tensor_json.hpp"
#include "response_functions.h"
#include "ResponseExceptions.hpp"
#include "ExcitedResponse.hpp"
#include "string"
#include "timer.h"
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

TEST_CASE("Test Response Functions", "Testing Basic Functionality of Response Functions") {

    int argc = 1;
    char **argv = nullptr;
    initialize(argc, argv);// initializes a world argument with argc and argv
    {

        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);


        const char *input_file = "excited.in";
        if (world.rank() == 0) { print("input filename: ", input_file); }
        if (!file_exists(input_file)) { throw Input_Error{}; }


        auto calc_params = initialize_calc_params(world, std::string(input_file));

        auto &[ground_calculation, molecule, response_parameters] = calc_params;


        vecfuncT ground_orbitals = ground_calculation.orbitals();

        print("Printing Norms of Ground Orbitals",norm2s_T(world, ground_orbitals));

        ExcitedResponse calc(world, calc_params);
        calc.solve(world);


        finalize();
    }
}
int main(int argc, char *argv[]) {
    //World& world=initialize(argc, argv);// initializes a world argument with argc and argv
    // World world(SafeMPI::COMM_WORLD);
    // startup(world, argc, argv, true);
    try {
        int result = Catch::Session().run(argc, argv);
        return result;
    } catch (const SafeMPI::Exception &e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException &e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException &e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char *s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string &s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception &e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}
