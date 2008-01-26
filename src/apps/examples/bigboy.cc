/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

/// \file bigboy.cc
/// \brief artificial benchmark for Cray XT testing

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
using namespace madness;

#include <examples/molecule.h>

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

static const double PI = 3.1415926535897932384;

static inline double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return sqrt(xx*xx + yy*yy + zz*zz);
}

class MolecularDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;

public:
    MolecularDensityFunctor(const Molecule& molecule) 
        : molecule(molecule)
    {}

    double operator()(const coordT& x) const {
        double sum = 0.0;
        for (int i=0; i<molecule.natom(); i++) {
            const Atom& atom = molecule.get_atom(i);
            double r = distance(x[0],x[1],x[2],atom.x,atom.y,atom.z);
            sum += exp(-2*r + 1e-6);
        }
        return sum;
    }
};


void bigboy(World& world, const Molecule& molecule) {
    // Solves Poisson's equation for a molecule-like charge density

    START_TIMER;
    functionT rho = factoryT(world).functor(functorT(new MolecularDensityFunctor(molecule)));
    END_TIMER("project density");
    START_TIMER;
    LoadBalImpl<3> lb(rho);
    world.gop.fence();
    FunctionDefaults<3>::pmap = lb.load_balance();
    rho = copy(rho, FunctionDefaults<3>::pmap);
    END_TIMER("load balancing");
    START_TIMER;
    double rho_norm = rho.norm2();
    END_TIMER("density norm");

    if (world.rank() == 0) print("rho norm", rho_norm);

    SeparatedConvolution<double,3> op = CoulombOperator<double,3>(world, 
                                                                  FunctionDefaults<3>::k, 
                                                                  1e-3, 
                                                                  FunctionDefaults<3>::thresh);
    world.am.print_stats();
    world.taskq.print_stats();

    START_TIMER;
    functionT V = apply(op,rho);
    END_TIMER("Coulomb operator");

    world.am.print_stats();
    world.taskq.print_stats();

    world.gop.fence();
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        
        // Process 0 reads the molecule information and
        // broadcasts to everyone else
        Molecule molecule;
        if (world.rank() == 0) molecule.read_file("bigboy.input");
        world.gop.broadcast_serializable(molecule, 0);
        
        // Use a cell big enough to have exp(-sqrt(2*I)*r) decay to
        // 1e-6 with I=1ev=0.037Eh --> need 50 a.u. either side of the molecule
        double L = molecule.bounding_cube();
        L += 50.0;
        for (int i=0; i<3; i++) {
            FunctionDefaults<3>::cell(i,0) = -L;
            FunctionDefaults<3>::cell(i,1) =  L;
        }
        
        // Setup initial defaults for numerical functions
        FunctionDefaults<3>::k = 12;
        FunctionDefaults<3>::thresh = 1e-6;
        FunctionDefaults<3>::refine = true;
        FunctionDefaults<3>::initial_level = 2;
        FunctionDefaults<3>::truncate_mode = 1;  
        
        // Warm and fuzzy for the user
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS example Hartree-Fock program");
            print(" ------------------------------------\n");
            molecule.print();
            print("\n");
            print("            box size ", L);
            print(" number of processes ", world.size());
            //if (nelec&1) throw "Closed shell only --- number of electrons is odd";
        }
        
        bigboy(world, molecule);

        world.gop.fence();

    } catch (const MPI::Exception& e) {
        //        print(e);
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

    MPI::Finalize();
    
    return 0;
}

