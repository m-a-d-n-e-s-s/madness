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


TEST_CASE("Run ground and excited-state") {
    // Set up the run directories
    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    const std::string molecule_name = "Be";
    const std::string xc = "hf";
    const std::string op = "excited-state";


    auto schema=runSchema(xc);
    auto mol_path=addPath(schema.molecule_path,molecule_name);

    try {

        auto num_states = set_excited_states(schema.rdb, schema.molecule_path, molecule_name, xc);
        auto moldft_path = run_moldft_path(world, schema.xc_path, xc, mol_path, molecule_name);

        try {
            runExcitedStates(world, moldft_path, num_states, xc);
        } catch (const SafeMPI::Exception &e) {
            print(e);
        } catch (const madness::MadnessException &e) {
            std::cout << e << std::endl;
        } catch (const madness::TensorException &e) {
            print(e);
        } catch (const char *s) { print(s); } catch (const std::string &s) {
            print(s);
        } catch (const std::exception &e) { print(e.what()); } catch (...) {
            error("caught unhandled exception");
        }

    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }

}
