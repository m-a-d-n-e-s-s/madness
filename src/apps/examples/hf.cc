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

  
  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file hartree-fock.cc
/// \brief example solution of the closed-shell HF equations

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

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;

public:
    MolecularPotentialFunctor(const Molecule& molecule) 
        : molecule(molecule)
    {}

    double operator()(const coordT& x) const {
        return molecule.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};


double guess(const Vector<double,3>& r) {
  return exp(-1.5*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])+1e-2);
}

void hf_solve(World& world, const Molecule& molecule, int nmo) {
    // Presently nmo is ignored and we are doing a single electron
    // calculation for the purpose of debugging

    double tol = FunctionDefaults<3>::thresh*0.1; // <<<<<< why?  is not the box size or truncate mode

    START_TIMER;
    functionT nuclear_potential = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(tol); 
    END_TIMER("project potential");
    START_TIMER;
    functionT psi = factoryT(world).f(guess);
    END_TIMER("project guess psi");

    double orbital_energy = -0.5;
    for (int iter=0; iter<1; iter++) {
        START_TIMER;
        psi.truncate().scale(1.0/psi.norm2());
        END_TIMER("truncate psi");

        START_TIMER;
        functionT Vpsi = nuclear_potential*psi;
        END_TIMER("V*psi");

        START_TIMER;
        Vpsi.truncate();
        END_TIMER("truncate Vpsi");

        operatorT op = BSHOperator<double,3>(world, sqrt(-2.0*orbital_energy), 
                                             FunctionDefaults<3>::k, 1e-2, FunctionDefaults<3>::thresh);

        START_TIMER;
        functionT psi_new = apply(op,Vpsi).scale(-2.0);
        END_TIMER("BSH operator * Vpsi");

//         const double* nflop = op.get_nflop();
//         for (int i=0; i<64; i++) print(i, nflop[i]);

        if (world.rank() == 0) print(iter,"delta_psi =",(psi-psi_new).norm2(),psi.norm2(), psi_new.norm2());
        psi = psi_new;
    }
        
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
        if (world.rank() == 0) molecule.read_file("input");
        world.gop.broadcast_serializable(molecule, 0);
        
        int nelec = int(molecule.total_nuclear_charge());
        int nmo = nelec/2;
        
        // Use a cell big enough to have exp(-sqrt(2*I)*r) decay to
        // 1e-6 with I=1ev=0.037Eh --> need 50 a.u. either side of the molecule
        double L = molecule.bounding_cube();
        L += 50.0;
        for (int i=0; i<3; i++) {
            FunctionDefaults<3>::cell(i,0) = -L;
            FunctionDefaults<3>::cell(i,1) =  L;
        }
        
        // Setup initial defaults for numerical functions
        FunctionDefaults<3>::k = 8;
        FunctionDefaults<3>::thresh = 1e-6;
        FunctionDefaults<3>::refine = true;
        FunctionDefaults<3>::initial_level = 2;
        FunctionDefaults<3>::truncate_mode = 0;  
        
        // Warm and fuzzy for the user
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS example Hartree-Fock program");
            print(" ------------------------------------\n");
            molecule.print();
            print("\n");
            print("            box size ", L);
            print(" number of electrons ", nelec);
            print(" number of  orbitals ", nmo);
            //if (nelec&1) throw "Closed shell only --- number of electrons is odd";
        }
        
        hf_solve(world, molecule, nmo);

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

