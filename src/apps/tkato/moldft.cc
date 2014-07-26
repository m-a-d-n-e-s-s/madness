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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <tensor/solvers.h>
using namespace madness;


#include <tkato/scfparam.h>
#include <tkato/scf.h>
#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>
#include <moldft/corepotential.h>
#include <moldft/xcfunctional.h>

class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;

// Computes molecular energy as a function of the geometry
// This is cludgy ... need better factorization of functionality
// between calculation, main program and this ... or just merge it all.
class MolecularEnergy : public OptimizationTargetInterface {
    World& world;
    SCF& calc;
    mutable double coords_sum;     // sum of square of coords at last solved geometry
    mutable double E; //< Current energy

public:
    MolecularEnergy(World& world, SCF& calc)
        : world(world)
        , calc(calc)
        , coords_sum(-1.0)
    {}

    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        double xsq = x.sumsq();
        if (xsq == coords_sum) {
            return calc.current_energy;
        }
        calc.molecule.set_all_coords(x.reshape(calc.molecule.natom(),3));
        coords_sum = xsq;

        // The below is missing convergence test logic, etc.

        // Make the nuclear potential, initial orbitals, etc.
        for (unsigned int proto=0; proto<calc.param.protocol_data.size(); proto++) {
            calc.set_protocol(world,calc.param.protocol_data[proto]);
            calc.make_nuclear_potential(world);
            calc.project_ao_basis(world);

            if (proto == 0) {
                if (calc.param.restart) {
                    calc.load_mos(world);
                }
                else {
                    calc.initial_guess(world);
                    calc.param.restart = true;
                }
            }
            else {
                calc.project(world);
            }

            // If the basis for the inital guess was not sto-3g
            // switch to sto-3g since this is needed for analysis
            // of the MOs and orbital localization
            if (calc.param.aobasis != "sto-3g") {
                calc.param.aobasis = "sto-3g";
                calc.project_ao_basis(world);
            }

            calc.solve(world);
            calc.save_mos(world);
        }
        return calc.current_energy;
    }

    madness::Tensor<double> gradient(const Tensor<double>& x) {
        value(x); // Ensures DFT equations are solved at this geometry

        return calc.derivatives(world);
    }
};



int main(int argc, char** argv) {
    initialize(argc, argv);

    { // limit lifetime of world so that finalize() can execute cleanly
      World world(SafeMPI::COMM_WORLD);

      try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

        std::cout.precision(6);

        // Process 0 reads input information and broadcasts
        SCF calc(world, "input");

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
          calc.set_protocol(world,1e-4);
          calc.make_nuclear_potential(world);
          calc.initial_load_bal(world);
        }

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
          geom.optimize(geomcoord);
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
    finalize();

    return 0;
}
