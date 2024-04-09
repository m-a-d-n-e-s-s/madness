//
// Created by adrianhurtado on 2/17/22.
//

#include "ResponseExceptions.hpp"
#include "madness/chem/SCF.h"
#include "madness/tensor/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_data_base.hpp"
#include "response_functions.h"
#include "string"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
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

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    // first step is to read the molecule directory for molecules... check if it exists else throw error

    std::string op = "excited-state";
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
            json excited_data = generate_excited_data(molecule_path, xc, 4);
            std::ofstream ofs("molecules/frequency.json");
            j_read.merge_patch(excited_data);
            ofs << std::setw(4) << j_read << std::endl;

        } else {
            json data = generate_excited_data(molecule_path, xc, 4);
            std::ofstream ofs("molecules/frequency.json");
            ofs << std::setw(4) << data << std::endl;
        }


    } catch (const std::exception &e) { cout << e.what() << std::endl; }
}