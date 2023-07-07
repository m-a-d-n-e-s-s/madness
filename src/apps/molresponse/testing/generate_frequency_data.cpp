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

    std::string op = "nuclear";
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
    auto freq_json_path = "json_data/frequency.json";
    try {
        if (std::filesystem::exists(freq_json_path)) {
            std::ifstream ifs(freq_json_path);
            std::cout << "Trying to read frequency.json" << std::endl;
            json j_read;
            ifs >> j_read;
            std::cout << "READ IT" << std::endl;
            json data = generate_response_data(molecule_path, xc, op, {0});

            print(data);
            j_read.merge_patch(data);
            // make the keyword and add the data
            std::ofstream ofs(freq_json_path);
            ofs << std::setw(4) << j_read << std::endl;

        } else {
            json data;
            json new_data = generate_response_data(molecule_path, xc, op, {0});
            std::ofstream ofs(freq_json_path);
            ofs << std::setw(4) << data << std::endl;
        }


    } catch (const std::exception &e) { cout << e.what() << std::endl; }
}