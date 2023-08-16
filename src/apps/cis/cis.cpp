/*
 * cis.cpp
 *
 *  Created on: Aug 10, 2015
 *      Author: Jakob S. Kottmann
 */

#include <madness/chem/TDHF.h>
#include <madness/mra/commandlineparser.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>

using namespace madness;


int main(int argc, char **argv) {
    int error = 0;
    World& world = initialize(argc, argv);
    try {
        if (world.rank() == 0) {
            print_header1("CIS -- compute DFT and Hartree-Fock excited states in CIS/TDA approximation");
        }

        //const double time_start = wall_time();
        std::cout.precision(6);

        startup(world, argc, argv, true);
        if (world.rank()==0) print(info::print_revision_information());

        printf("starting at time %.1f\n", wall_time());
        print_meminfo(world.rank(), "startup");

        if (world.rank() == 0) {
            std::cout << "\n\n";
            std::cout << "-------------------------------------------------------------------------------------\n";
            std::cout << "SOLVING MRA-CIS as described in \n";
            std::cout << "J.S. Kottmann, S. Höfener ,F.A. Bischoff\n";
            std::cout
                    << "Numerically accurate linear response-properties in the configuration-interaction singles (CIS) approximation \n";
            std::cout << "Phys. Chem. Chem. Phys., 2015, 17, 31453-31462\n";
            std::cout << "DOI: 10.1039/C5CP00345H\n";
            std::cout << "-------------------------------------------------------------------------------------\n";
            std::cout << "\n\n";
        }

        commandlineparser parser(argc, argv);
        parser.print_map();

        if (parser.key_exists("help")) {
            TDHF::help();

        } else if (parser.key_exists("print_parameters")) {
            TDHF::print_parameters();

        } else if (parser.key_exists("test")) {
            print("entering test mode");
            error = TDHF::test(world, parser);

        } else {

            TDHF tdhf(world, parser);

            print_header2("input section");
            if (world.rank() == 0) {
                tdhf.get_calcparam().print("dft","end");
                print("");
                tdhf.get_parameters().print("response","end");
                tdhf.get_calc()->molecule.print();
            }

            // solve the CIS equations
            const double time_scf_start = wall_time();
            tdhf.prepare_calculation();
            const double time_scf_end = wall_time();
            if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

            const double time_cis_start = wall_time();
            std::vector<CC_vecfunction> roots = tdhf.solve_cis();
            const double time_cis_end = wall_time();
            if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

            if (world.rank() == 0) {
                std::cout << std::setfill(' ');
                std::cout << "\n\n\n";
                std::cout << "--------------------------------------------------\n";
                std::cout << "MRA-CIS ended \n";
                std::cout << "--------------------------------------------------\n";
                std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
                std::cout << std::setw(25) << "time cis" << " = " << time_cis_end - time_cis_start << "\n";
                std::cout << "--------------------------------------------------\n";
            }
            tdhf.analyze(roots);
        }
        world.gop.fence();
        if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
        print_stats(world);
    } catch (std::exception& e) {
        print("an exception occurred");
        print(e.what());
    } catch (...) {
        print("an unknown exception occurred");
    }
    finalize();
    return error;
}



