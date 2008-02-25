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

/// \file hartree-fock.cc
/// \brief example solution of the closed-shell HF equations

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
using namespace madness;

#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr<operatorT> poperatorT;

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

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis) 
        : molecule(molecule), aobasis(aobasis)
    {}

    double operator()(const coordT& x) const {
        return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
    }
};


class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
    const AtomicBasisFunction aofunc;
public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc) : aofunc(aofunc)
    {}
 
    double operator()(const coordT& x) const {
        return aofunc(x[0], x[1], x[2]);
    }
};


vector<functionT> project_ao_basis(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis) {
    vector<functionT> ao(aobasis.nbf(molecule));

    for (int i=0; i<aobasis.nbf(molecule); i++) {
        functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule,i)));
        ao[i] = factoryT(world).functor(aofunc).initial_level(3).nofence();
    }
    world.gop.fence();

    vector<double> norms = norm2(world, ao);
    if (world.rank() == 0) {
        for (int i=0; i<aobasis.nbf(molecule); i++) {
            if (world.rank() == 0) print(i,"ao.norm", norms[i]);
            norms[i] = 1.0/norms[i];
        }
    }

    scale(world, ao, norms);
    
    return ao;
}


Tensor<double> make_kinetic_energy_matrix(World& world, const vector<functionT>& v) {

    reconstruct(world, v);

    int n = v.size();
    Tensor<double> r(n,n);
    for (int axis=0; axis<3; axis++) {
        vector<functionT> dv = diff(world,v,axis);
        r += inner(world, dv, dv);

        dv.clear(); world.gop.fence(); // Allow function memory to be freed
    }
    
    return r.scale(0.5);
}

vector<functionT> apply_potential(World& world, const functionT& vnuc, const vector<functionT>& v) {
    return mul(world, vnuc, v);
}

Tensor<double> make_potential_energy_matrix(World& world, const functionT& vnuc, const vector<functionT>& v) {
    vector<functionT> Vpsi = apply_potential(world,vnuc,v);
    compress(world,Vpsi,false);
    compress(world,v);
    Tensor<double> r = inner(world,v,Vpsi);

    Vpsi.clear(); world.gop.fence(); // Allow function memory to be freed
    return r;
}

struct CalculationParameters {
    int nalpha;
    int nbeta;
    bool spin_restricted;

    void read_file(const std::string& filename) {
        std::ifstream f(filename.c_str());
        position_stream(f, "dft");
        f >> nalpha >> nbeta >> spin_restricted;
    }

    template <typename Archive>
    void serialize(Archive& ar) {ar & nalpha & nbeta & spin_restricted;}
};


Tensor<double> make_fock_matrix(World& world, const functionT& vnuc, const vector<functionT>& mos) {
    Tensor<double> kinetic = make_kinetic_energy_matrix(world, mos);
    Tensor<double> potential = make_potential_energy_matrix(world, vnuc, mos);
    return kinetic + potential;
}

vector<functionT> make_guess_orbitals(World& world, const Molecule& molecule, 
                                      const AtomicBasisSet& aobasis, const CalculationParameters& param,
                                      const functionT& vnuc, Tensor<double>& e)
{
    vector<functionT> ao = project_ao_basis(world, molecule, aobasis);
    Tensor<double> overlap = inner(world,ao,ao);
    Tensor<double> kinetic = make_kinetic_energy_matrix(world, ao);
    Tensor<double> potential = make_potential_energy_matrix(world, vnuc, ao);
    Tensor<double> fock = kinetic + potential;
    Tensor<double> c;
    sygv(fock, overlap, 1, &c, &e);
    if (world.rank() == 0) {
        print("THIS iS THE OVERLAP MATRIX");
        print(overlap);
        print("THIS iS THE KINETIC MATRIX");
        print(kinetic);
        print("THIS iS THE POTENTIAL MATRIX");
        print(potential);
        print("THESE ARE THE EIGENVECTORS");
        print(c);
        print("THESE ARE THE EIGENVALUES");
        print(e);
    }

    int nmo = max(param.nalpha, param.nbeta);

    ao = transform(world, ao, c(_,Slice(0,nmo-1)));

    overlap = inner(world,ao,ao);
    if (world.rank() == 0) {
        print("NEW overlap matrix");
        print(overlap);
    }

    e = e(Slice(0,nmo-1));

    return ao;
}
    
vector<poperatorT> make_bsh_operators(World& world, const Tensor<double>& evals, int k, double lo, double tol) {
    int nmo = evals.dim[0];
    vector<poperatorT> ops(nmo);
    for (int i=0; i<nmo; i++) {
        ops[i] = poperatorT(BSHOperatorPtr<double,3>(world, sqrt(-2.0*evals(i)), k, lo, tol));
    }
    return ops;
}


void normalize(World& world, vector<functionT>& v) {
    vector<double> nn = norm2(world, v);
    for (unsigned int i=0; i<v.size(); i++) v[i].scale(1.0/nn[i]);
}

void hf_solve(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, const CalculationParameters& param) {
    const double lo = molecule.smallest_length_scale();
    if (world.rank() == 0) print("smallest length scale", lo);

    //operatorT coulombOP = CoulombOperator<double, NDIM>(world, FunctionDefaults<3>::k, smallest_lengthscale, 

    double tol = FunctionDefaults<3>::thresh;
    int k = FunctionDefaults<3>::k;
    Tensor<double> evals;

    START_TIMER;
    functionT vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(tol*0.1); 
    END_TIMER("project potential");

//     START_TIMER;
//     functionT guess_density = factoryT(world).functor(functorT(new MolecularGuessDensityFunctor(molecule,aobasis))).thresh(tol); 
//     END_TIMER("project potential");

//     guess_density.compress();
//     print("and the number of electrons is", guess_density.trace());
//     guess_density.reconstruct();
//     print("and the number of electrons is", guess_density.trace());

    START_TIMER;
    vector<functionT> mos = make_guess_orbitals(world, molecule, aobasis, param, vnuc, evals);
    END_TIMER("Make guess mos");

    int nmo = mos.size();

    for (int iter=0; iter<10; iter++) {
        START_TIMER;
        vector<functionT> Vpsi = apply_potential(world,vnuc,mos);
        END_TIMER("Apply potential");

        truncate(world,Vpsi);

        START_TIMER;
        vector<poperatorT> ops = make_bsh_operators(world, evals, k, lo, tol);
        vector<functionT> new_mo = apply(world, ops, Vpsi);
        END_TIMER("Apply BSH operators");

        truncate(world, new_mo);
        for (int i=0; i<nmo; i++) mos.push_back(new_mo[i]);

        ops.clear();
        Vpsi.clear();
        new_mo.clear();
        world.gop.fence();

        START_TIMER;
        Tensor<double> overlap = inner(world,mos,mos);
        Tensor<double> fock_matrix = make_fock_matrix(world, vnuc, mos);
        END_TIMER("Make matrices");
        
        Tensor<double> c, e;
        sygv(fock_matrix, overlap, 1, &c, &e); // In practice better to avoid this entirely!

        START_TIMER;
        mos = transform(world, mos, c(_,Slice(0,nmo-1)));
        truncate(world, mos);
        normalize(world, mos);
        END_TIMER("Transform");
        evals = e(Slice(0,nmo-1));

        if (world.rank() == 0) {
            print(iter,"EVALS");
            print(evals);
        }
    }

//     START_TIMER;
//     functionT psi = factoryT(world).f(guess);
//     END_TIMER("project guess psi");

//     double orbital_energy = -0.5;
//     for (int iter=0; iter<20; iter++) {
//         START_TIMER;
//         psi.truncate().scale(1.0/psi.norm2());
//         END_TIMER("truncate psi");

//         START_TIMER;
//         functionT Vpsi = vnuc*psi;
//         END_TIMER("V*psi");

//         START_TIMER;
//         Vpsi.truncate();
//         END_TIMER("truncate Vpsi");

//         operatorT op = BSHOperator<double,3>(world, sqrt(-2.0*orbital_energy), 
//                                              FunctionDefaults<3>::k, 1e-2, FunctionDefaults<3>::thresh);

//         START_TIMER;
//         functionT psi_new = apply(op,Vpsi).scale(-2.0);
//         END_TIMER("BSH operator * Vpsi");

// //         const double* nflop = op.get_nflop();
// //         for (int i=0; i<64; i++) print(i, nflop[i]);

//         double deltanorm = (psi-psi_new).norm2();
//         double newnorm = psi_new.norm2();
//         if (world.rank() == 0) print(iter,"delta norm =",deltanorm, "     new norm =",newnorm);
//         psi = psi_new;
//     }
        
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        
        // Process 0 reads the molecule information and broadcasts
        Molecule molecule;
        CalculationParameters param;
        if (world.rank() == 0) {
            molecule.read_file("input");
            param.read_file("input");
        }
            
        world.gop.broadcast_serializable(molecule, 0);
        world.gop.broadcast_serializable(param, 0);
        
        // Process 0 reads the LCAO guess information and broadcasts
        AtomicBasisSet aobasis;
        if (world.rank() == 0) aobasis.read_file("sto-3g");
        world.gop.broadcast_serializable(aobasis, 0);

        // Use a cell big enough to have exp(-sqrt(2*I)*r) decay to
        // 1e-6 with I=1ev=0.037Eh --> need 50 a.u. either side of the molecule
        double L = molecule.bounding_cube();
        L += 50.0;
        for (int i=0; i<3; i++) {
            FunctionDefaults<3>::cell(i,0) = -L;
            FunctionDefaults<3>::cell(i,1) =  L;
        }
        
        // Setup initial defaults for numerical functions
        FunctionDefaults<3>::k = 6;
        FunctionDefaults<3>::thresh = 1e-4;
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
            print(" number of electrons ", param.nalpha, param.nbeta);
            print("     spin restricted ", param.spin_restricted);
        }
        
        hf_solve(world, molecule, aobasis, param);

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

