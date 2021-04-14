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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/nemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

*/

#include <chem/nemo.h>
#include <chem/molecular_optimizer.h>
#include <chem/SCFOperators.h>
#include <chem/projector.h>
#include <madness/misc/gitinfo.h>

using namespace madness;


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  NEMO -- Hartree-Fock using numerical exponential molecular orbitals \n");
    	printf("starting at time %.1f\n", wall_time());

    }
    startup(world,argc,argv);
    std::cout.precision(6);

    if (world.rank()==0) print(info::print_revision_information());

    try {

        commandlineparser parser(argc,argv);
        std::shared_ptr<Nemo> nemo(new Nemo(world,parser));
        if (world.rank()==0) nemo->get_param().print("dft","end");
        if (world.rank()==0) nemo->get_calc()->param.print("dft","end");

        // optimize the geometry if requested
        if (nemo->get_param().gopt()) {
            print("\n\n Geometry Optimization                      ");
            print(" ----------------------------------------------------------\n");
//            calc->param.gprint(world);

            Tensor<double> geomcoord =nemo->get_calc()->molecule.get_all_coords().flat();
//            MolecularOptimizer geom(std::shared_ptr<MolecularOptimizationTargetInterface>(new Nemo(world, calc)),
            MolecularOptimizer geom(world,nemo);
//            MolecularOptimizer geom(nemo,
//                    calc->param.gmaxiter(),
//                    calc->param.gtol(),  //tol
//                    calc->param.gval(),  //value prec
//                    calc->param.gprec()); // grad prec
//            geom.set_update(calc->param.algopt);
//            geom.set_test(calc->param.gtest);

            // compute initial hessian
            if (nemo->get_param().ginitial_hessian()) {
                nemo->value();
                Tensor<double> hess=nemo->hessian(nemo->get_calc()->molecule.get_all_coords());
                geom.set_hessian(hess);
            }
            geom.optimize(geomcoord);
        } else {

            // compute the energy to get converged orbitals
//            Nemo nemo(world,calc);
            const double energy=nemo->value();
            if (world.rank()==0) {
                printf("final energy   %12.8f\n", energy);
                printf("finished at time %.1f\n", wall_time());
            }

        }

        // compute the hessian
        if (nemo->param.hessian()) nemo->hessian(nemo->get_calc()->molecule.get_all_coords());


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


    finalize();
    return 0;
}
