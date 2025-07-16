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

/*!
  \file examples/znemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

*/

#if !defined USE_GENTENSOR

#include <madness/chem/znemo.h>
using namespace madness;

int main(int argc, char** argv) {

    World& world=initialize(argc, argv,false);
    if (world.rank() == 0) {
        print_header1(" ZNEMO -- complex Hartree-Fock using numerical exponential molecular orbitals");
        printf("starting at time %.1f\n", wall_time());
    }

    startup(world,argc,argv,true);
    std::cout.precision(6);
    if (world.rank()==0) print(info::print_revision_information());


    commandlineparser parser(argc,argv);
    parser.set_keyval("no_orient","true");  // needs to be passed to molecule *before* construction
    if (parser.key_exists("help")) {
        Znemo::help();

    } else if (parser.key_exists("print_parameters")) {
        Znemo::print_parameters();

    } else {

        try {

            bool analyze_only=parser.key_exists("analyze");

            std::shared_ptr<Znemo> znemo(new Znemo(world,parser));

            // optimize the geometry if requested
            if (znemo->get_calc_param().gopt()) {
                print_header2("Geometry Optimization");

                Tensor<double> geomcoord = znemo->molecule().get_all_coords().flat();
                MolecularOptimizer geom(world,parser,znemo);
                geom.parameters.set_derived_value<std::vector<std::string> >("remove_dof",{});
                geom.parameters.print("geometry optimization parameters","end");
                geom.optimize(geomcoord);
            } else {

                // compute the energy to get converged orbitals
                double energy=0.0;
                if (not analyze_only) {
                    energy=znemo->value();
                } else {
                    auto zmos=znemo->read_reference();
                    for (int i=0; i<znemo->get_calc_param().nalpha(); ++i) znemo->amo.push_back(zmos.first.get_mos()[i]);
                    for (int i=0; i<znemo->get_calc_param().nbeta(); ++i) znemo->bmo.push_back(zmos.second.get_mos()[i]);
                    energy=znemo->analyze();
                }

                if (world.rank()==0) {
                    printf("final energy   %12.8f\n", energy);
                    printf("finished at time %.1f\n", wall_time());
                }
            }

        } catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        } catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        } catch (const madness::TensorException& e) {
            print(e);
            error("caught a Tensor exception");
        } catch (const char* s) {
            print(s);
            error("caught a string exception");
        } catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        } catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        } catch (...) {
            error("caught unhandled exception");
        }
    }

    finalize();
    return 0;
}
#else

#include <iostream>
int main() {
    std::cout << "U need to configure with -D ENABLE_GENTENSOR=0 to enable znemo\n";
    return 0;
}

#endif
