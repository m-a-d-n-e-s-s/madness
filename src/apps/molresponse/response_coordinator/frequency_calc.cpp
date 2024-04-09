//
// Created by adrianhurtado on 1/1/22.
//
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "response_functions.h"
#include "timer.h"

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

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    std::cout.precision(6);
    std::string filename = "response.in";

    try {
        auto calc_params = initialize_calc_params(world, filename);
        auto omega = calc_params.response_parameters.omega();
        FrequencyResponse calc(world, calc_params, omega, dipole_generator);
        // Warm and fuzzy for the user
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS Time-Dependent Density Functional Theory Frequency Response "
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
    } catch (const SafeMPI::Exception &e) { print(e); } catch (const madness::MadnessException &e) {
        std::cout << e << std::endl;
    } catch (const madness::TensorException &e) { print(e); } catch (const char *s) {
        print(s);
    } catch (const std::string &s) { print(s); } catch (const std::exception &e) {
        print(e.what());
    } catch (...) { error("caught unhandled exception"); }

    return result;
}
