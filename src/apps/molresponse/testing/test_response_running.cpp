//
// Created by adrianhurtado on 1/1/22.
//
#define CATCH_CONFIG_RUNNER
#include "ExcitedResponse.hpp"
#include "ResponseExceptions.hpp"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
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

    // Here we run all the tests
    result = Catch::Session().run(argc, argv);
    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

void set_default_response_parameters(ResponseParameters &r_params) {

    r_params.set_user_defined_value("maxiter", size_t(10));
    r_params.set_user_defined_value("archive", std::string("../restartdata"));
    r_params.set_user_defined_value("kain", true);
    r_params.set_user_defined_value("maxsub", size_t(10));
}

void RunResponse(World &world, std::string filename, double frequency, std::string property,
                 std::string xc, std::filesystem::path moldft_path) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_default_response_parameters(r_params);
    r_params.set_user_defined_value("xc", xc);
    if (property == "dipole") {
        r_params.set_user_defined_value("dipole", true);
    } else if (property == "nuclear") {
        r_params.set_user_defined_value("nuclear", true);
    }

    r_params.set_user_defined_value("omega", frequency);
    r_params.set_user_defined_value("first_order", true);
    r_params.set_user_defined_value("plot_all_orbitals", true);

    r_params.set_user_defined_value("save", true);
    std::string s_frequency = std::to_string(frequency);
    print(s_frequency);
    auto sp = s_frequency.find(".");
    s_frequency = s_frequency.replace(sp, sp, "-");
    print(s_frequency);
    std::string run_name = property + "_" + xc + "_" + s_frequency;
    print(run_name);

    // set r_params to restart true if restart file exist

    auto run_path = moldft_path;
    run_path += "/";
    run_path += std::filesystem::path(run_name);
    std::cout << run_path << endl;

    if (std::filesystem::is_directory(run_path)) {

        cout << "Response directory found " << std::endl;

    } else {// create the file
        std::filesystem::create_directory(run_path);
        cout << "Creating response_path directory" << std::endl;
    }

    auto restart_path = run_path;
    std::filesystem::current_path(run_path);
    std::string restart_string = "restart_" + run_name;
    restart_path += "/";
    restart_path += restart_string;
    r_params.set_user_defined_value("save_file", restart_string);

    if (std::filesystem::exists(restart_path)) {
        r_params.set_user_defined_value("restart", true);
        r_params.set_user_defined_value("restart_file", restart_string);
    } else {
        std::cout << "Restart File Does Not Exist!!!" << endl;
    }

    write_response_input(r_params, filename);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    auto &[ground_calculation, molecule, r_params_copy] = calc_params;
    vecfuncT ground_orbitals = ground_calculation.orbitals();

    print(norm2s_T(world, ground_orbitals));

    /*
    ExcitedResponse calc(world, calc_params);

    calc.solve(world);
    calc.output_json();
     */
}


TEST_CASE("Run Frequency Response ") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    try {
        auto moldft_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/hf");
        std::filesystem::current_path(moldft_path);

        RunResponse(world, "response.in", 0, "dipole", "hf", moldft_path);
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}

TEST_CASE("Run A Few Frequency Response ") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    try {
        auto moldft_path = std::filesystem::path(
                "/home/adrianhurtado/projects/madness-test-suite/tests_response/hf");
        std::filesystem::current_path(moldft_path);
        std::vector<double> frequencies = {0, 0.025, 0.050, 0.075, 0.1};
        for (const auto &freq: frequencies) {
            RunResponse(world, "response.in", freq, "dipole", "hf", moldft_path);
        }
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}
