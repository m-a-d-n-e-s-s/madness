//
// Created by adrianhurtado on 2/11/22.
//
#define CATCH_CONFIG_RUNNER

#include <xc.h>

#include "ExcitedResponse.hpp"
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
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

// Run a single chosen molecule and check results
TEST_CASE("Test MOLDFT JSON") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    const std::string xc = "hf";

    auto schema = runSchema(xc);
    auto m_schema = moldftSchema("Be", xc, schema);



    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 2);
    param1.set_user_defined_value<std::string>("xc", xc);
    param1.set_user_defined_value<double>("l", 200);
    param1.set_user_defined_value<vector<double>>("protocol", {1e-4});
    param1.set_user_defined_value<std::string>("localize", "canon");
    // write restart file
    write_test_input test_input(param1, "moldft.in", m_schema.mol_path);// molecule HF commandlineparser parser;


    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());
    SCF calc(world, parser);
    calc.set_protocol<3>(world, 1e-4);
    MolecularEnergy ME(world, calc);
    // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
    ME.value(calc.molecule.get_all_coords().flat());// ugh!
    world.gop.fence();
    ME.output_calc_info_schema();

}