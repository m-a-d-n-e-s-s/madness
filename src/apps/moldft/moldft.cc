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
/// \defgroup moldft The molecular density funcitonal and Hartree-Fock code


#include <chem/SCF.h>

int main(int argc, char** argv) {
    TAU_START("main()");
    TAU_START("initialize()");
    initialize(argc, argv);
    TAU_STOP("initialize()");

    TAU_START("World lifetime");
    { // limit lifetime of world so that finalize() can execute cleanly
      World world(SafeMPI::COMM_WORLD);

      try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
	print_meminfo(world.rank(), "startup");
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap< Key<3> >(world)));

        std::cout.precision(6);

        // Process 0 reads input information and broadcasts
        const char * inpname = "input";
        for (int i=1; i<argc; i++) {
            if (argv[i][0] != '-') {
                inpname = argv[i];
                break;
            }
        }
        if (world.rank() == 0) print(inpname);
        SCF calc(world, inpname);

        // Warm and fuzzy for the user
        if (world.rank() == 0) {
          print("\n\n");
          print(" MADNESS Hartree-Fock and Density Functional Theory Program");
          print(" ----------------------------------------------------------\n");
          print("\n");
          calc.molecule.print();
          print("\n");
          calc.param.print(world);
        }

        // Come up with an initial OK data map
        if (world.size() > 1) {
          calc.set_protocol<3>(world,calc.param.econv);
          calc.make_nuclear_potential(world);
          calc.initial_load_bal(world);
        }
//vama
        calc.set_protocol<3>(world,calc.param.protocol_data[0]);


        if ( calc.param.gopt) {
          print("\n\n Geometry Optimization                      ");
          print(" ----------------------------------------------------------\n");
          calc.param.gprint(world);

          Tensor<double> geomcoord = calc.molecule.get_all_coords().flat();
          QuasiNewton geom(std::shared_ptr<OptimizationTargetInterface>(new MolecularEnergy(world, calc)),
                           calc.param.gmaxiter,
                           calc.param.gtol,  //tol
                           calc.param.gval,  //value prec
                           calc.param.gprec); // grad prec
          geom.set_update(calc.param.algopt);
          geom.set_test(calc.param.gtest);
          long ncoord = calc.molecule.natom()*3;
          Tensor<double> h(ncoord,ncoord);
          for (int i=0; i<ncoord; ++i) h(i,i) = 0.5;
          geom.set_hessian(h);
          geom.optimize(geomcoord);
        }
        else if (calc.param.tdksprop) {
          print("\n\n Propagation of Kohn-Sham equation                      ");
          print(" ----------------------------------------------------------\n");
//          calc.propagate(world,VextCosFunctor<double>(world,new DipoleFunctor(2),0.1),0);
          calc.propagate(world,0.1,0);
        }
        else {
          MolecularEnergy E(world, calc);
          E.value(calc.molecule.get_all_coords().flat()); // ugh!
          if (calc.param.derivatives) calc.derivatives(world);
          if (calc.param.dipole) calc.dipole(world);
        }

        //        if (calc.param.twoint) {
        //Tensor<double> g = calc.twoint(world,calc.amo);
        //cout << g;
        // }

        calc.do_plots(world);

      }
      catch (const SafeMPI::Exception& e) {
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
      catch (char* s) {
        print(s);
        error("caught a string exception");
      }
      catch (const char* s) {
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
    TAU_STOP("World lifetime");
    finalize();
    TAU_STOP("main()");

    return 0;
}
