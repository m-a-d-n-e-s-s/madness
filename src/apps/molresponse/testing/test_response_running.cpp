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
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"
#include "runners.hpp"

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

TEST_CASE("Test Gamma Functions Response ") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    try {
        auto response_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/orbital_analysis/"
                "11_Ne/excited_state");
        std::filesystem::current_path(response_path);
        std::cout << "before initialize" << std::endl;
        initialize_excited_restart(world, "restart_excited.in", 4, "hf");

        auto calc_params = initialize_calc_params(world, "restart_excited.in");
        auto[ground_calculation, molecule, r_params] = calc_params;
        r_params.set_user_defined_value("print_level", 10);
        vecfuncT ground_orbitals = ground_calculation.orbitals();
        print(norm2s_T(world, ground_orbitals));
        density_vector rho = set_density_type(world, r_params, ground_calculation);

        TDDFT old_calc(world, rho);
        ExcitedResponse new_calc(world, calc_params);

        double thresh = 1e-4;
        ResponseTester tester;
        tester.load_calc(world, &new_calc, thresh);

        TDDFT_Tester tddft;
        tddft.load_calc(world, &old_calc, thresh);

        auto old_chi = old_calc.Chi;
        auto new_chi = new_calc.get_chi();

        SECTION("Testing Gamma") {

            auto new_gamma = tester.compute_gamma_full(world, &new_calc, thresh);
            auto old_gamma = tddft.compute_gamma_full(world, &old_calc, thresh);

            auto xgx_old = inner(old_chi, old_gamma);
            auto xgx_new = inner(new_chi, new_gamma);

            bool gamma_equal = is_equal(xgx_old, xgx_new, thresh);
            print("old", xgx_old);
            print("new", xgx_new);

            REQUIRE(gamma_equal);
        }

        SECTION("Testing V0X and FOX") {

            auto[VOX, FOX] = tester.compute_VFOX(world, &new_calc, true);
            auto[oVOX, oFOX] = tddft.compute_VFOX(world, &old_calc, true);

            auto old_xVx = inner(old_chi, oVOX);
            auto old_xFx = inner(old_chi, oFOX);
            auto new_xVx = inner(new_chi, VOX);
            auto new_xFx = inner(new_chi, FOX);
            print("VOX first");
            print("old\n", old_xVx);
            print("new\n", new_xVx);
            bool V_equal = is_equal(old_xVx, new_xVx, thresh);
            // because I can't get the tensor == to work
            REQUIRE(V_equal);

            print("VOX first");
            print("old\n", old_xFx);
            print("new\n", new_xFx);

            bool F_equal = is_equal(old_xFx, new_xFx, thresh);
            // because I can't get the tensor == to work
            REQUIRE(V_equal);
        }SECTION("Testing Lambda") {

            auto new_lambda = tester.compute_lambda_X(world, &new_calc, thresh);
            auto old_lambda = tddft.compute_lambda_X(world, &old_calc, thresh);

            auto old_xlx = inner(old_chi, old_lambda);
            auto new_xlx = inner(new_chi, new_lambda);
            print("old\n", old_xlx);
            print("new\n", new_xlx);
            bool _equal = is_equal(old_xlx, new_xlx, thresh);
            // because I can't get the tensor == to work


            REQUIRE(_equal);
        }


    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}

TEST_CASE("Run Frequency Response ") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    try {
        auto moldft_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/orbital_analysis/"
                "10_Be");
        std::filesystem::current_path(moldft_path);
        auto restart_path = moldft_path;
        auto[next_restart, success] =
        RunResponse(world, "response.in", 0, "dipole", "hf", moldft_path, restart_path);
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}

TEST_CASE("Run A Few Frequency Response ") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    try {
        auto moldft_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/orbital_analysis/"
                "10_Be");

        std::filesystem::current_path(moldft_path);
        std::vector<double> frequencies = {0, 0.025, 0.050, 0.075, 0.1};
        runFrequencyTests(world, moldft_path, frequencies, "hf", "dipole");
        // add a restart path
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}

