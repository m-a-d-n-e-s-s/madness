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
#include <ctime>
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



struct CalculationParameters {
    // !!! If you add more data don't forget to add them to serialize method.

    // First list input parameters
    double charge;              ///< Total molecular charge
    double smear;               ///< Smearing parameter
    double econv;               ///< Energy convergence
    double dconv;               ///< Density convergence
    double L;                   ///< User coordinates box size
    double maxrotn;             ///< Step restriction used in autoshift algorithm
    int nvalpha;                ///< Number of alpha virtuals to compute
    int nvbeta;                 ///< Number of beta virtuals to compute
    int nopen;                  ///< Number of unpaired electrons = napha-nbeta
    bool spin_restricted;       ///< True if spin restricted
    // Next list inferred parameters
    int nalpha;                 ///< Number of alpha spin electrons
    int nbeta;                  ///< Number of beta  spin electrons
    int nmo_alpha;              ///< Number of alpha spin molecular orbitals
    int nmo_beta;               ///< Number of beta  spin molecular orbitals
    double lo;                  ///< Smallest length scale we need to resolve
    
    CalculationParameters()
        : charge(0.0)
        , smear(0.0)
        , econv(1e-5)
        , dconv(1e-4)
        , L(0.0)
        , maxrotn(0.25)
        , nvalpha(1)
        , nvbeta(1)
        , nopen(0)
        , spin_restricted(true)
        , nalpha(0)
        , nbeta(0)
        , nmo_alpha(0)
        , nmo_beta(0)
        , lo(1e-10)
    {}

    void read_file(const std::string& filename) {
        std::ifstream f(filename.c_str());
        position_stream(f, "dft");
        string s;
        while (f >> s) {
            if (s == "end") {
                break;
            } else if (s == "charge") {
                f >> charge;
            } else if (s == "smear") {
                f >> smear;
            } else if (s == "econv") {
                f >> econv;
            } else if (s == "dconv") {
                f >> dconv;
            } else if (s == "L") {
                f >> L;
            } else if (s == "maxrotn") {
                f >> maxrotn;
            } else if (s == "nvalpha") {
                f >> nvalpha;
            } else if (s == "nvbeta") {
                f >> nvbeta;
            } else if (s == "nopen") {
                f >> nopen;
            } else if (s == "unrestricted") {
                spin_restricted = false;
            } else if (s == "restricted") {
                spin_restricted = true;
            } else {
                std::cout << "moldft: unrecognized input keyword " << s << std::endl;
                MADNESS_EXCEPTION("input error",0);
            }
        }
    }

    void set_molecular_info(const Molecule& molecule, const AtomicBasisSet& aobasis) {
        double z = molecule.total_nuclear_charge();
        int nelec = z - charge;
        if (fabs(nelec+charge-z) > 1e-6) {
            error("non-integer number of electrons?", nelec+charge-z);
        }
        nalpha = (nelec + nopen)/2;
        nbeta  = (nelec - nopen)/2;
        if (nalpha < 0) error("negative number of alpha electrons?", nalpha);
        if (nbeta < 0) error("negative number of beta electrons?", nbeta);
        if ((nalpha+nbeta) != nelec) error("nalpha+nbeta != nelec", nalpha+nbeta);
        nmo_alpha = nalpha + nvalpha;
        nmo_beta = nbeta + nvbeta;
        if (nalpha != nbeta) spin_restricted = true;

        // Ensure we have enough basis functions to guess the requested
        // number of states ... a minimal basis for a closed-shell atom
        // might not have any functions for virtuals.
        int nbf = aobasis.nbf(molecule);
        nmo_alpha = min(nbf,nmo_alpha);  
        nmo_beta = min(nbf,nmo_beta);  
        if (nalpha>nbf || nbeta>nbf) error("too few basis functions?", nbf);
        
        // Unless overridden by the user use a cell big enough to
        // have exp(-sqrt(2*I)*r) decay to 1e-6 with I=1ev=0.037Eh
        // --> need 50 a.u. either side of the molecule
        
        if (L == 0.0) {
            L = molecule.bounding_cube() + 50.0;
        }

        lo = molecule.smallest_length_scale();
    }
    
    void print(World& world) const {
        time_t t = time((time_t *) 0);
        char *tmp = ctime(&t);
        tmp[strlen(tmp)-1] = 0; // lose the trailing newline
        madness::print(" date of calculation ", tmp);
        madness::print(" number of processes ", world.size());
        madness::print("        total charge ", charge);
        madness::print("            smearing ", smear);
        madness::print(" number of electrons ", nalpha, nbeta);
        madness::print("  number of orbitals ", nmo_alpha, nmo_beta);
        madness::print("     spin restricted ", spin_restricted);
        madness::print("  energy convergence ", econv);
        madness::print(" density convergence ", dconv);
        madness::print("    maximum rotation ", maxrotn);
    }
    
    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & L & nvalpha & nvbeta & nopen & spin_restricted;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta;
    }
};

struct Calculation {
    Molecule molecule;            ///< Molecular coordinates, etc.
    CalculationParameters param;  ///< User input data, nalpha, etc.
    AtomicBasisSet aobasis;       ///< Currently always the STO-3G basis
    functionT vnuc;               ///< The effective nuclear potential
    vector<functionT> amo, bmo;   ///< alpha and beta molecular orbitals
    Tensor<double> aocc, bocc;    ///< alpha and beta occupation numbers
    Tensor<double> aeps, beps;    ///< alpha and beta BSH shifts (eigenvalues if canonical)
    poperatorT coulop;            ///< Coulomb Green function
    double vtol;                  ///< Tolerance used for potentials and density

    Calculation(World& world, const char* filename) {
        if (world.rank() == 0) {
            molecule.read_file(filename);
            param.read_file(filename);
            aobasis.read_file("sto-3g");
            molecule.center();
            param.set_molecular_info(molecule,aobasis);
        }
            
        world.gop.broadcast_serializable(molecule, 0);
        world.gop.broadcast_serializable(param, 0);
        world.gop.broadcast_serializable(aobasis, 0);

        for (int i=0; i<3; i++) {
            FunctionDefaults<3>::cell(i,0) = -param.L;
            FunctionDefaults<3>::cell(i,1) =  param.L;
        }
        
        // Setup initial defaults for numerical functions
        FunctionDefaults<3>::k = 6;
        FunctionDefaults<3>::thresh = 1e-4;
        FunctionDefaults<3>::refine = true;
        FunctionDefaults<3>::initial_level = 2;
        FunctionDefaults<3>::truncate_mode = 1;  

        double safety = 0.1;
        vtol = FunctionDefaults<3>::thresh*safety;

        coulop = poperatorT(CoulombOperatorPtr<double, 3>(world, 
                                                          FunctionDefaults<3>::k, 
                                                          param.lo, 
                                                          vtol));
    }

    void make_nuclear_potential(World& world) {
        START_TIMER;
        vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(vtol); 
        vnuc.truncate();
        vnuc.reconstruct();
        END_TIMER("project nuclear potential");
    }

    vector<functionT> project_ao_basis(World& world) {
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

    Tensor<double> kinetic_energy_matrix(World& world, const vector<functionT>& v) {
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

    
    /// Initializes alpha and beta mos, occupation numbers, eigenvalues
    void initial_guess(World& world) {
        vector<functionT> ao   = project_ao_basis(world);
        Tensor<double> overlap = inner(world,ao,ao);
        Tensor<double> kinetic = kinetic_energy_matrix(world, ao);

        vector<functionT> vpsi = mul(world, vnuc, ao);

        Tensor<double> potential = inner(world, vpsi, ao); 
        vpsi.clear(); 

        Tensor<double> fock = kinetic + potential;
        fock = 0.5*(fock + transpose(fock));

        Tensor<double> c,e;
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

        amo = transform(world, ao, c(_,Slice(0,param.nmo_alpha-1)));
        truncate(world, amo);
        normalize(world, amo);

        print("this is the new overlap");
        print(inner(world,amo,amo));

        aeps = e(Slice(0,param.nmo_alpha-1))*0.5; // 0.5 from using bare core Hamiltonian
        aocc = Tensor<double>(param.nalpha);
        for (int i=0; i<param.nalpha; i++) aocc[i] = 1.0;

        if (!param.spin_restricted) {
            bmo = transform(world, ao, c(_,Slice(0,param.nmo_beta-1)));
            truncate(world, bmo);
            normalize(world, bmo);
            beps = e(Slice(0,param.nmo_beta-1));
            bocc = Tensor<double>(param.nbeta);
            for (int i=0; i<param.nbeta; i++) bocc[i] = 1.0;
        }
    }

    functionT make_density(World& world, const Tensor<double>& occ, const vector<functionT>& v) {
        vector<functionT> vsq = square(world, v);
        compress(world,vsq);
        functionT rho = factoryT(world).thresh(vtol);
        rho.compress();
        for (unsigned int i=0; i<vsq.size(); i++) {
            rho.gaxpy(1.0,vsq[i],occ[i],false);
        }
        world.gop.fence();
        vsq.clear();
        world.gop.fence();
        return rho;
    }

    vector<poperatorT> make_bsh_operators(World& world, const Tensor<double>& evals) {
        int nmo = evals.dim[0];
        vector<poperatorT> ops(nmo);
        int k = FunctionDefaults<3>::k;
        double tol = FunctionDefaults<3>::thresh;
        for (int i=0; i<nmo; i++) {
            double eps = evals(i);
            if (eps > 0) eps = -0.05;
            ops[i] = poperatorT(BSHOperatorPtr<double,3>(world, sqrt(-2.0*eps), k, param.lo, tol));
        }
        return ops;
    }

    vector<functionT> apply_hf_exchange(World& world, 
                                        const Tensor<double>& occ,
                                        const vector<functionT>& psi,
                                        const vector<functionT>& f) {
        int nocc = psi.size();
        int nf = f.size();

        // Problem here is balancing memory usage vs. parallel efficiency.
        // Once we start using localized orbitals both the occupied
        // and target functions will have limited support and hence
        // simply parallelizing either the loop over f or over occupied
        // will not generate real concurrency ... need to parallelize
        // them both.  
        //
        // For now just parallelize one loop but will need more
        // intelligence soon.

        vector<functionT> Kf = zero_functions<double,3>(world, nf);

        compress(world,Kf);
        reconstruct(world,psi);
        for (int i=0; i<nocc; i++) {
            if (occ[i] > 0.0) {
                vector<functionT> psif = mul(world, psi[i], f);
                set_thresh(world, psif, vtol);  //<<<<<<<<<<<<<<<<<<<<<<<<< since cannot yet put in apply

                truncate(world,psif);
                psif = apply(world, *coulop, psif);
                truncate(world, psif);

                psif = mul(world, psi[i], psif);

                gaxpy(world, 1.0, Kf, occ[i], psif);
            }
        }
        truncate(world, Kf, vtol);
        return Kf;
    }
        

    void update(World& world, 
                const functionT& vlocal, 
                Tensor<double>& occ, 
                Tensor<double>& eps,
                vector<functionT>& psi) 
    {
        int nmo = psi.size();
        vector<functionT> Vpsi = mul(world, vlocal, psi);
        vector<functionT> Kpsi = apply_hf_exchange(world, occ, psi, psi);
        gaxpy(world, 1.0, Vpsi, -1.0, Kpsi);
        Kpsi.clear();
        Vpsi[0].scale(-2.0);
        truncate(world,Vpsi);
        vector<poperatorT> ops = make_bsh_operators(world, eps);

        set_thresh(world, Vpsi, FunctionDefaults<3>::thresh);  //<<<<< Since cannot set in apply
        vector<functionT> new_psi = apply(world, ops, Vpsi);

        ops.clear();
        truncate(world, new_psi);
        Vpsi[0].scale(-0.5);

        // Approx orthog of new to old (so can level shift and increase sparsity)
        compress(world, psi);
        compress(world, new_psi);
        Tensor<double> oldnew = inner(world, psi, new_psi);
        for (int i=0; i<nmo; i++) {
            for (int j=0; j<nmo; j++) {
                new_psi[i].gaxpy(1.0,psi[j],-oldnew(j,i), false);
            }
        }
        world.gop.fence();
        truncate(world,new_psi);

        vector<functionT> Vnew_psi = mul(world, vlocal, new_psi);
        Kpsi = apply_hf_exchange(world, occ, psi, new_psi);
        gaxpy(world, 1.0, Vnew_psi, -1.0, Kpsi);
        Kpsi.clear();

        for (int i=0; i<nmo; i++) {
            psi.push_back(new_psi[i]);
            Vpsi.push_back(Vnew_psi[i]);
        }
        new_psi.clear();
        Vnew_psi.clear();

        Tensor<double> potential = inner(world, Vpsi, psi);
        Vpsi.clear();
        world.gop.fence(); // Frees the memory

        Tensor<double> overlap = inner(world, psi, psi);
        if (world.rank() == 0) {print("overlap"); print(overlap);}

        Tensor<double> fock = potential + kinetic_energy_matrix(world, psi);
        fock = (fock + transpose(fock))*0.5;
        if (world.rank() == 0) {print("fock"); print(fock);}

        // Examine Fock matrix elements to assess convergence
        double maxocc = -1e300;
        double minvirt = 1e300;
        double maxoffd = 0.0;
        for (int i=0; i<nmo; i++) {
            for (int j=0; j<nmo; j++) {
                maxoffd = max(maxoffd, fabs(fock(i,j+nmo)/sqrt(overlap(i,i)*overlap(j+nmo,j+nmo))));
            }
        }
        for (int i=0; i<2*nmo; i++) {
            if (i<nmo && occ(i)>0.0) maxocc = max(maxocc,fock(i,i));
            else minvirt = min(minvirt,fock(i,i)/overlap(i,i));
        }
        for (int i=0; i<2*nmo; i++) {
            if (occ(i) > 0.0) maxocc = max(maxocc,fock(i,i));
            else minvirt = min(minvirt,fock(i,i));
        }
        double mingap = minvirt - maxocc;

        // Determine automatic level shift by constraining the largest
        // Jacobi rotation mixing occupied and virtual.  The
        // denominator must be positive and we want the overall
        // rotation maxoffd/(mingap+autoshift) to be less than
        // maxrotn=0.25=sin(15 degrees).  Thus, ...
        double autoshift = max(0.0, maxoffd/param.maxrotn - mingap);
        if (world.rank() == 0) {
            print("Max. gradient", maxoffd);
            print("     Min. gap", mingap);
            print("    Autoshift", autoshift);
        }            
                              
        for (int i=0; i<nmo; i++) fock(i,i) -= autoshift;//param.lshift; // Apply level shift
        Tensor<double> c, e;
        sygv(fock, overlap, 1, &c, &e);
        for (int i=0; i<nmo; i++) e(i) += autoshift;//param.lshift; // Undo level shift

        psi = transform(world, psi, c(_,Slice(0,nmo-1)));
        truncate(world, psi);
        normalize(world, psi);
        eps = e(Slice(0,nmo-1));

        print("PSI thresh at end of iter", psi[0].thresh());
        // NEED TO UPDATE OCC HERE IF SMEARING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

    void solve(World& world) {
        for (int iter=0; iter<30; iter++) {
            //aeps[0] = -4.732;
            //aeps[1] = -0.309;
            functionT arho = make_density(world, aocc, amo);
            functionT brho = param.spin_restricted ? arho : make_density(world, bocc, bmo);
            functionT vlocal = vnuc + apply(*coulop, arho+brho);
            vlocal.truncate();
            // For DFT add into vlocal the exchange-correlation potentials

            // Could do a lot of the alpha and beta work in parallel
            // if we need more parallelism ... however it will double
            // the memory too!

            update(world, vlocal, aocc, aeps, amo);
            if (!param.spin_restricted) update(world, vlocal, bocc, beps, bmo);

            if (world.rank() == 0) {
                print(iter,"alpha evals");
                print(aeps);
            }
        }
    }
};



int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        
        // Process 0 reads input information and broadcasts
        Calculation calc(world, "input");
        
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

        // Make the nuclear potential, initial orbitals, etc.
        calc.make_nuclear_potential(world);
        calc.initial_guess(world);
        
        calc.solve(world);

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

