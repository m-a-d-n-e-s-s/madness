//
// Created by adrianhurtado on 2/11/22.
//
#define CATCH_CONFIG_RUNNER

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
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

TEST_CASE("Run MOLDFT/RESPONSE") {


    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);
// INPUTS
// root  is the current path
// molecule_path is the path where molecules are


    auto root = std::filesystem::current_path();//="/"+molecule_name;
    // first step is to read the molecule directory for molecules... check if it exists else throw error
    auto molecule_path = root;
    molecule_path += "/molecules";
    std::string xc = "hf";
    auto xc_path = create_xc_path_and_directory(root, xc);
    std::string property = "dipole";

    ResponseDataBase response_data_base = ResponseDataBase();
    if (std::filesystem::exists("molecules/frequency.json")) {
        std::ifstream ifs("molecules/frequency.json");
        std::cout << "Trying to read frequency.json" << std::endl;
        json j_read;
        ifs >> j_read;
        std::cout << "READ IT" << std::endl;
        response_data_base.set_data(j_read);
    } else {
        json data = generate_excited_data(molecule_path, xc, 4);
        std::ofstream ofs("molecules/frequency.json");
        ofs << std::setw(4) << data << std::endl;
        response_data_base.set_data(data);
    }
    try {
        if (std::filesystem::is_directory(molecule_path)) {
            // for every molecule within the molecule path
            for (const std::filesystem::directory_entry &mol_path:
                    std::filesystem::directory_iterator(molecule_path)) {
                size_t num_states{0};
                std::filesystem::current_path(xc_path);

                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    std::cout << "\n\n----------------------------------------------------\n";
                    std::cout << "Beginning Tests for Molecule: " << molecule_name << "\n";

                    num_states = set_excited_states(response_data_base, molecule_path, molecule_name, xc);
                    auto moldft_path = run_moldft_path(world, xc_path, xc, mol_path, molecule_name);
                    // states.
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
                    // now check if the answers exist.  if the answers do not exist run
                    // response else check the answers
                }
            }
            std::cout << "Please check what happens when I get to this point of the loop"
                      << std::endl;
            // Now check if restart file exists and if calc_info.json exists
        } else {
            std::cout << "did not find molecules" << std::endl;
        }
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}
