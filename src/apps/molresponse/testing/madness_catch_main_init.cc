#define CATCH_CONFIG_RUNNER
#include <chem/write_test_input.h>

#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "madness/world/worldmem.h"


int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    try {

        {
            CalculationParameters param1;
            param1.set_user_defined_value("maxiter", 2);
            param1.set_user_defined_value("protocol", std::vector<double>({1.e-3}));
            // write restart file
            // write restart file
            write_test_input test_input(param1, "hf");// molecule HF
            commandlineparser parser;
            parser.set_keyval("input", test_input.filename());
            SCF calc(world, parser);
            calc.set_protocol<3>(world, 1e-4);
            MolecularEnergy ME(world, calc);
            //double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
            ME.value(calc.molecule.get_all_coords().flat());// ugh!
            ME.output_calc_info_schema();

            world.gop.fence();
        }

        // Here we run all the tests
        result = Catch::Session().run(argc, argv);
        return result;


    } catch (const SafeMPI::Exception &e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException &e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException &e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char *s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string &s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception &e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) { error("caught unhandled exception"); }

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}
