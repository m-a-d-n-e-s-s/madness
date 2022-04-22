//
// Created by adrianhurtado on 2/11/22.
//

#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_data_base.hpp"
#include "response_functions.h"
#include "runners.hpp"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"


using namespace madness;

int main(int argc, char *argv[]) {

    if (argc != 3) {

        std::cout << "Wrong number of inputs" << std::endl;
        return 1;
    }

    World &world = madness::initialize(argc, argv);
    world.gop.fence();
    startup(world, argc, argv);

    std::cout.precision(6);

    const std::string xc{argv[1]};
    const std::string op{argv[2]};
    auto schema = runSchema(xc);

    try {
        if (std::filesystem::is_directory(schema.molecule_path)) {
            for (const std::filesystem::directory_entry &mol_path:
                 std::filesystem::directory_iterator(schema.molecule_path)) {

                std::filesystem::current_path(schema.xc_path);
                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    try {

                        auto m_schema = moldftSchema(molecule_name, xc, schema);
                        m_schema.print();
                        moldft(world, m_schema, true);
                        auto f_schema = frequencySchema(schema, m_schema, op);
                        runFrequencyTests(world, f_schema);

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
    return 0;
}
