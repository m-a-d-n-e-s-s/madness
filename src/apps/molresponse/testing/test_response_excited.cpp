//
// Created by adrianhurtado on 2/11/22.
//

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
    std::cout.precision(6);

    std::string xc = "hf";
    auto schema = runSchema(xc);

    try {
        if (std::filesystem::is_directory(schema.molecule_path)) {
            // for every molecule within the molecule path
            for (const std::filesystem::directory_entry &mol_path:
                    std::filesystem::directory_iterator(schema.molecule_path)) {
                size_t num_states{0};
                std::filesystem::current_path(schema.xc_path);

                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    std::cout << "\n\n----------------------------------------------------\n";
                    std::cout << "Beginning Tests for Molecule: " << molecule_name << "\n";

                    num_states = set_excited_states(schema.rdb, schema.molecule_path, molecule_name, xc);
                    auto m_schema = moldftSchema(molecule_name, xc, schema);
                    m_schema.print();
                    run_moldft_path(world, m_schema);
                    auto excited_schema = excitedSchema(m_schema, num_states);
                    excited_schema.print();


                    // states.
                    try {
                        bool success = runExcited(world, excited_schema, true);
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
