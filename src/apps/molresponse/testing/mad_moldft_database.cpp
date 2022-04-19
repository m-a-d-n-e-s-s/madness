//
// Created by adrianhurtado on 1/1/22.
//
#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "apps/chem/SCF.h"
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

    try {

        auto num_states = set_excited_states(schema.rdb, schema.molecule_path, molecule_name, xc);
        auto m_schema = moldftSchema(molecule_name, xc, schema);
        m_schema.print();
        run_moldft_path(world, m_schema);

    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }


    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}


