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
  \file examples/nemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

*/

#include <madness/chem/SCFOperators.h>
#include <madness/chem/molecular_optimizer.h>
#include <madness/chem/nemo.h>
#include <madness/chem/projector.h>
#include <madness/misc/info.h>

using namespace madness;


int main(int argc, char** argv) {

    World& world=initialize(argc, argv,false);
    if (world.rank() == 0) {
        print_header1("NEMO -- Hartree-Fock using numerical exponential molecular orbitals");
    	printf("starting at time %.1f\n", wall_time());
    }

    startup(world,argc,argv,true);
    std::cout.precision(6);
    if (world.rank()==0) print(info::print_revision_information());


    commandlineparser parser(argc,argv);
    if (parser.key_exists("help")) {
        Nemo::help();

    } else if (parser.key_exists("print_parameters")) {
        Nemo::print_parameters();

    } else {

        try {

            std::shared_ptr<Nemo> nemo(new Nemo(world, parser));

            print_header2("input section");
            if (world.rank() == 0) {
                nemo->get_calc_param().print("dft", "end");
                nemo->molecule().print();
            }

            // optimize the geometry if requested
            if (nemo->get_calc_param().gopt()) {
                print_header2("Geometry Optimization");
                MolecularOptimizer geom(world, parser, nemo);
                geom.parameters.print("geoopt", "end");

                // compute the energy to get converged orbitals
                print_header2("computing initial wave function");
                nemo->value();

                // reduce print level
                nemo->get_calc()->param.set_derived_value("print_level",2);
                nemo->get_calc()->param.set_derived_value("print_level",2);
                nemo->get_calc()->set_print_timings(false);

                // compute initial hessian
                if (nemo->get_calc_param().ginitial_hessian()) {
                    Tensor<double> hess = nemo->hessian(nemo->get_calc()->molecule.get_all_coords());
                    geom.set_hessian(hess);
                }

                print_header2("Starting geometry optimization");
                Tensor<double> geomcoord = nemo->get_calc()->molecule.get_all_coords().flat();
                geom.optimize(geomcoord);
            }


            double energy=nemo->value();
            if (world.rank() == 0) {
                printf("final energy   %12.8f\n", energy);
                printf("finished at time %.1f\n", wall_time());
            }

            // compute the hessian
            if (nemo->get_nemo_param().hessian()) nemo->hessian(nemo->get_calc()->molecule.get_all_coords());


        } catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        } catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        } catch (const madness::TensorException& e) {
            print(e);
            error("caught a Tensor exception");
        } catch (const char *s) {
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
