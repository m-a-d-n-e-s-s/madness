/*
 *
 *   Main source file for the Dirac Fock code
 *
 */


#include <madness/misc/info.h>
#include "DF.h"    // All Dirac-Fock functions/objects enter through this
#include <stdlib.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
// static inline int file_exists(const char * inpname)
// {
//      struct stat buffer;
//      int rc = stat(inpname, &buffer);
//      return (rc==0);
// }
#endif

int main(int argc, char **argv) {

    World& world = initialize(argc, argv);
    if (world.rank() == 0) print_header1("DIRAC -- molecular Dirac-Fock code");

    // Load info for MADNESS numerical routines
    startup(world, argc, argv, true);
    if (world.rank() == 0) print(info::print_revision_information());

    commandlineparser parser(argc, argv);
    if (parser.key_exists("help")) {
        DF::help();

    } else if (parser.key_exists("print_parameters")) {
        DF::print_parameters();

    } else {

        try {
            // Create the DF object
            DF my_calc(world, parser.value("input").c_str());

            // Have it iterate to convergence
            my_calc.solve(world);
        } catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        }
        catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        }
        catch (const madness::TensorException& e) {
            print(e);
            error("caught a Tensor exception");
        }
        catch (const char *s) {
            print(s);
            error("caught a string exception");
        }
        catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        }
        catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        }
        catch (...) {
            error("caught unhandled exception");
        }

        // Nearly all memory will be freed at this point
        world.gop.fence();
        world.gop.fence();
        print_stats(world);
    } // world is dead -- ready to finalize
    finalize();
    return 0;
}

