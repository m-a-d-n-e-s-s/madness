//
// Created by adrianhurtado on 2/17/22.
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
#include "response_data_base.hpp"

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

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    // first step is to read the molecule directory for molecules... check if it exists else throw error

    std::string property = "dipole";
    auto molecule_path = root;
    molecule_path += "/molecules";
    std::string xc = "hf";
    auto xc_path = root;
    xc_path += "/";
    xc_path += std::filesystem::path(xc);
    if (std::filesystem::is_directory(xc_path)) {
        cout << "XC directory found " << xc << "\n";
    } else {// create the file
        std::filesystem::create_directory(xc_path);
        cout << "Creating XC directory for " << xc << ":\n";
    }
    try {

        if (std::filesystem::exists("molecules/frequency.json")) {
            std::ifstream ifs("molecules/frequency.json");
            std::cout << "Trying to read frequency.json" << std::endl;
            json j_read;
            ifs >> j_read;
            std::cout << "READ IT" << std::endl;

        } else {
            json data = generate_response_data(molecule_path, xc, property, {0});
            std::ofstream ofs("molecules/frequency.json");
            ofs << std::setw(4) << data << std::endl;
            ResponseDataBase response_data_base(data);
            ResponseDataBase rbd;
        }


    } catch (const std::exception &e) {
        cout << e.what() << std::endl;
    }
}