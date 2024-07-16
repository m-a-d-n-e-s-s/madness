/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
*/

/// \file moldft.cc
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density functional and Hartree-Fock code

#include<madness/chem/SCF.h>
#include"madness/mra/commandlineparser.h"
#include<madness/chem/molopt.h>
#include <madness/world/worldmem.h>
#include <madness/misc/info.h>


#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// static inline int file_exists(const char *inpname) {
//     struct stat buffer;
//     int rc = stat(inpname, &buffer);
//     return (rc == 0);
// }

#endif

using namespace madness;


static double ttt, sss;

static void START_TIMER(World &world) {
    world.gop.fence();
    ttt = wall_time();
    sss = cpu_time();
}

static void END_TIMER(World &world, const char *msg) {
    ttt = wall_time() - ttt;
    sss = cpu_time() - sss;
    if (world.rank() == 0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

int main(int argc, char **argv) {

    World &world = initialize(argc, argv);
    if (world.rank() == 0) {
        print_header1("MOLDFT -- molecular DFT and Hartree-Fock code");
    }

    { // limit lifetime of world so that finalize() can execute cleanly
        START_TIMER(world);
        try {
            // Load info for MADNESS numerical routines
            startup(world, argc, argv, true);
            if (world.rank() == 0) print(info::print_revision_information());

            commandlineparser parser(argc, argv);

            if (parser.key_exists("help")) {
                SCF::help();

            } else if (parser.key_exists("print_parameters")) {
                SCF::print_parameters();

            } else {
                if (world.rank() == 0) print("input filename: ", parser.value("input"));


                print_meminfo(world.rank(), "startup");
                FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3> >(world)));

                std::cout.precision(6);
                SCF calc(world, parser);

                // Warm and fuzzy for the user
                if (world.rank() == 0) {
                    print("\n\n");
                    print(" MADNESS Hartree-Fock and Density Functional Theory Program");
                    print(" ----------------------------------------------------------\n");
                    print("\n");
                    calc.molecule.print();
                    print("\n");
                    calc.param.print("dft");
                }
                END_TIMER(world, "initialize");
                // Come up with an initial OK data map
                if (world.size() > 1) {
                    calc.set_protocol<3>(world, 1e-4);
                    calc.make_nuclear_potential(world);
                    calc.initial_load_bal(world);
                }
//vama
                calc.set_protocol<3>(world, calc.param.protocol()[0]);


                if (calc.param.gopt()) {
                    // print("\n\n Geometry Optimization                      ");
                    // print(" ----------------------------------------------------------\n");
                    // calc.param.gprint(world);

                    // Tensor<double> geomcoord = calc.molecule.get_all_coords().flat();
                    // QuasiNewton geom(std::shared_ptr<OptimizationTargetInterface>(new MolecularEnergy(world, calc)),
                    //                  calc.param.gmaxiter,
                    //                  calc.param.gtol,  //tol
                    //                  calc.param.gval,  //value prec
                    //                  calc.param.gprec); // grad prec
                    // geom.set_update(calc.param.algopt);
                    // geom.set_test(calc.param.gtest);
                    // long ncoord = calc.molecule.natom()*3;
                    // Tensor<double> h(ncoord,ncoord);
                    // for (int i=0; i<ncoord; ++i) h(i,i) = 0.5;
                    // geom.set_hessian(h);
                    // geom.optimize(geomcoord);

                    MolOpt opt(calc.param.gmaxiter(),
                               0.1,
                               calc.param.gval(),
                               calc.param.gtol(),
                               1e-3, //XTOL
                               1e-5, //EPREC
                               calc.param.gprec(),
                               (world.rank() == 0) ? 1 : 0, //print_level
                               calc.param.algopt());

                    MolecularEnergy target(world, calc);
                    opt.optimize(calc.molecule, target);
                } else {
                    MolecularEnergy E(world, calc);
                    double energy = E.value(calc.molecule.get_all_coords().flat()); // ugh!
                    if ((world.rank() == 0) and (calc.param.print_level() > 0)) {
                        printf("final energy=%16.8f \n", energy);
                        E.output_calc_info_schema();
                    }


                    functionT rho = calc.make_density(world, calc.aocc, calc.amo);
                    functionT brho = rho;
                    if (calc.param.nbeta() != 0 && !calc.param.spin_restricted())
                        brho = calc.make_density(world, calc.bocc, calc.bmo);
                    rho.gaxpy(1.0, brho, 1.0);

                    if (calc.param.derivatives()){
                        auto gradient=calc.derivatives(world, rho);
                        calc.e_data.add_gradient(gradient);
                    }
                    // automatically print dipole moment and output scf info
                    std::map<std::string, double> results;
                    results["scf_energy"] = calc.current_energy;
                    auto dipole_t = calc.dipole(world, rho);
                    calc.output_scf_info_schema(results, dipole_t);
                }

                //        if (calc.param.twoint) {
                //Tensor<double> g = calc.twoint(world,calc.amo);
                //cout << g;
                // }

                calc.do_plots(world);

            }
        }
        catch (const SafeMPI::Exception &e) {
            print(e);
            error("caught an MPI exception");
        }
        catch (const madness::MadnessException &e) {
            print(e);
            error("caught a MADNESS exception");
        }
        catch (const madness::TensorException &e) {
            print(e);
            error("caught a Tensor exception");
        }
        catch (const char *s) {
            print(s);
            error("caught a string exception");
        }
        catch (const std::string &s) {
            print(s);
            error("caught a string (class) exception");
        }
        catch (const std::exception &e) {
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
