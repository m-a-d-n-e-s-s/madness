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
    ResponseDataBase response_data_base = ResponseDataBase();

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
        if (std::filesystem::is_directory(molecule_path)) {
            for (const std::filesystem::directory_entry &mol_path:
                    std::filesystem::directory_iterator(molecule_path)) {

                std::filesystem::current_path(xc_path);
                //for each molecule in molecules directory
                bool moldft_results_exists = false;
                json moldft_answers;

                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    std::cout << "\n\n----------------------------------------------------\n";
                    std::cout << "Beginning Tests for Molecule: " << molecule_name << "\n";

                    // We would like to read the moldft results from corresponding json "molecule_name.json"
                    auto response_json_path =
                            generate_response_json_path(molecule_path, molecule_name, xc, property);

                    json response_json;
                    if (std::filesystem::exists(response_json_path)) {
                        std::cout << "response_json exists:" << std::endl;
                        response_json =
                                response_data_base.retrieve_data(molecule_name, xc, property);
                    } else {
                        try {
                            response_data_base.output_data(response_json_path.string(), molecule_name, xc, property);
                        } catch (json::exception &e) {
                            std::cout << e.what() << std::endl;
                        }
                    }
                    std::cout << response_json << std::endl;
                }
            }
            generate_response_data(molecule_path, xc, property, {0});


            // if the results exists then save them into answers
        }
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}