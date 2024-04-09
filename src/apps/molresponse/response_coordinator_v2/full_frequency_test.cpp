//
// Created by adrianhurtado on 2/11/22.
//

#include "ResponseExceptions.hpp"
#include "madness/external/nlohmann_json/json.hpp"
#include "madness/tensor/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "runners.hpp"


using namespace madness;

int main(int argc, char *argv[]) {

    if (argc != 4) {

        std::cout << "Wrong number of inputs" << std::endl;
        return 1;
    }

    World &world = madness::initialize(argc, argv);
    world.gop.fence();
    startup(world, argc, argv);

    std::cout.precision(6);

    // set last keyword to high to set high prec
    const std::string xc{argv[1]};
    const std::string op{argv[2]};
    const std::string precision{argv[3]};
    if (precision != "high" && precision != "low" && precision != "super") {
        if (world.rank() == 0) {
            std::cout << "Set precision to low high super" << std::endl;
        }
        return 1;
    }

    auto schema = runSchema(world, xc);

    try {
        if (std::filesystem::is_directory(schema.molecules)) {
            for (const std::filesystem::directory_entry &mol_path:
                 std::filesystem::directory_iterator(schema.molecules)) {

                std::filesystem::current_path(schema.xc_path);
                if (mol_path.path().extension() == ".mol") {
                    print(mol_path.path());
                    auto molecule_name = mol_path.path().stem();
                    try {

                        auto m_schema = moldftSchema(world, molecule_name, xc, schema);
                        moldft(world, m_schema, true, false, precision);
                        auto f_schema = frequencySchema(world, schema, m_schema, op, false);
                        runFrequencyTests(world, f_schema, precision);

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
