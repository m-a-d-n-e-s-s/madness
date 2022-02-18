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
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}


bool is_equal(const Tensor<double> &a, const Tensor<double> &b, double thresh) {


    // check if same elements and same number of dimensions
    if (a.size() != b.size() && a.ndim() != b.ndim()) { return false; }

    // now check if the dimensions are the same
    if (std::equal(a.dims(), a.dims() + a.ndim(), b.dims())) {

        auto af = a.flat();
        auto bf = b.flat();
        return std::equal(af.ptr(), af.ptr() + af.size(), bf.ptr(),
                          [&thresh](auto aa, auto bb) { return abs(aa - bb) <= thresh; });
    };
}

void set_and_write(ResponseParameters &r_params, std::string property, std::string xc,
                   double frequency);

static std::string set_frequency_path_and_restart(ResponseParameters parameters, std::string property,
                                           double frequency, std::filesystem::path path1,bool restart);

TEST_CASE("Test Dipole constructor") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);
    try {
        auto moldft_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/orbital_analysis/"
                "10_Be");
        std::filesystem::current_path(moldft_path);
        auto restart_path = moldft_path;
        std::string xc{"hf"};
        std::string property{"dipole"};
        double frequency{0.0};

        ResponseParameters r_params{};
        bool restart{false};
        set_and_write(r_params, property, xc, frequency);
        auto filename = set_frequency_path_and_restart(r_params, property, frequency, moldft_path,restart);


        auto calc_params = initialize_calc_params(world, filename);
        FrequencyResponse calc(world, calc_params,frequency, dipole_generator);
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
}

