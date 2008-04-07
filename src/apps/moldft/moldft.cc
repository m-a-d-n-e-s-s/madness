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

/// \file moldft.cc
/// \brief Molecular HF and DFT code

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

Tensor<double> sqrt(const Tensor<double>& s, double tol=1e-8) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    Tensor<double> c, e;
    //s.gaxpy(0.5,transpose(s),0.5); // Ensure exact symmetry
    syev(s, &c, &e);
    for (int i=0; i<n; i++) {
        if (e(i) < -tol) {
            MADNESS_EXCEPTION("Matrix square root: negative eigenvalue",i);
        }
        else if (e(i) < tol) { // Ugh ..
            print("Matrix square root: Warning: small eigenvalue ", i, e(i));
            e(i) = tol;
        }
        e(i) = 1.0/sqrt(e(i));
    }
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            c(j,i) *= e(i);
        }
    }
    return c;
}

Tensor<double> energy_weighted_orthog(const Tensor<double>& s, const Tensor<double> eps) {
    int n=s.dim[0], m=s.dim[1];
    MADNESS_ASSERT(n==m);
    Tensor<double> d(n,n);
    for (int i=0; i<n; i++) d(i,i) = eps(i);
    Tensor<double> c, e;
    sygv(d, s, 1, &c, &e);
    return c;
}


template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

struct CalculationParameters {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!                                                                   !!!
    // !!! If you add more data don't forget to add them to serialize method !!!
    // !!!                                                                   !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    int maxiter;                ///< Maximum number of iterations
    bool spin_restricted;       ///< True if spin restricted
    // Next list inferred parameters
    int nalpha;                 ///< Number of alpha spin electrons
    int nbeta;                  ///< Number of beta  spin electrons
    int nmo_alpha;              ///< Number of alpha spin molecular orbitals
    int nmo_beta;               ///< Number of beta  spin molecular orbitals
    double lo;                  ///< Smallest length scale we need to resolve

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & L & nvalpha & nvbeta & nopen & spin_restricted;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
    }

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
        , maxiter(20)
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
            } else if (s == "maxiter") {
                f >> maxiter;
            } else {
                std::cout << "moldft: unrecognized input keyword " << s << std::endl;
                MADNESS_EXCEPTION("input error",0);
            }
            if (nopen != 0) spin_restricted = false;
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
        if (nalpha != nbeta) spin_restricted = false;

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
        set_protocol(world, 1e-4);
    }

    void set_protocol(World& world, double thresh) {
        FunctionDefaults<3>::thresh = thresh;
        if (thresh >= 1e-4) FunctionDefaults<3>::k = 6;
        else if (thresh >= 1e-6) FunctionDefaults<3>::k = 8;
        else if (thresh >= 1e-8) FunctionDefaults<3>::k = 10;
        else FunctionDefaults<3>::k = 12;

        FunctionDefaults<3>::refine = true;
        FunctionDefaults<3>::initial_level = 2;
        FunctionDefaults<3>::truncate_mode = 1;  
        double safety = 0.1;
        vtol = FunctionDefaults<3>::thresh*safety;

        coulop = poperatorT(CoulombOperatorPtr<double, 3>(world, 
                                                          FunctionDefaults<3>::k, 
                                                          param.lo, 
                                                          vtol));
        if (world.rank() == 0) {
            print("\nSolving with thresh",thresh, "and k", FunctionDefaults<3>::k, "\n");
        }
    }

    void project(World& world) {
        reconstruct(world,amo);
        for (unsigned int i=0; i<amo.size(); i++) {
            amo[i] = madness::project(amo[i], FunctionDefaults<3>::k, FunctionDefaults<3>::thresh, false);
        }
        world.gop.fence();
        truncate(world,amo);
        normalize(world,amo);
        if (param.nbeta && !param.spin_restricted) {
            reconstruct(world,bmo);
            for (unsigned int i=0; i<bmo.size(); i++) {
                bmo[i] = madness::project(bmo[i], FunctionDefaults<3>::k, FunctionDefaults<3>::thresh, false);
            }
            world.gop.fence();
            truncate(world,bmo);
            normalize(world,bmo);
        }
    }
        
    void make_nuclear_potential(World& world) {
        START_TIMER;
        vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(vtol); 
        vnuc.truncate();
        vnuc.reconstruct();
        END_TIMER("Project vnuclear");
    }

    vector<functionT> project_ao_basis(World& world) {
        vector<functionT> ao(aobasis.nbf(molecule));

        for (int i=0; i<aobasis.nbf(molecule); i++) {
            functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule,i)));
            ao[i] = factoryT(world).functor(aofunc).initial_level(3).nofence();
        }
        world.gop.fence();

        vector<double> norms = norm2(world, ao);

        for (int i=0; i<aobasis.nbf(molecule); i++) {
            if (world.rank() == 0) print(i,"ao.norm", norms[i]);
            norms[i] = 1.0/norms[i];
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
            r += matrix_inner(world, dv, dv, true);
            dv.clear(); world.gop.fence(); // Allow function memory to be freed
        }
        
        return r.scale(0.5);
    }

    
    struct GuessDensity : public FunctionFunctorInterface<double,3> {
        const Molecule& molecule;
        const AtomicBasisSet& aobasis;
        double operator()(const coordT& x) const {
            return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
        }
        GuessDensity(const Molecule& molecule, const AtomicBasisSet& aobasis)
            : molecule(molecule), aobasis(aobasis) {}
    };


    /// Initializes alpha and beta mos, occupation numbers, eigenvalues
    void initial_guess(World& world) {
        START_TIMER;
        functionT rho = factoryT(world).functor(functorT(new GuessDensity(molecule, aobasis)));
        END_TIMER("guess density");
        double nel = rho.trace();
        if (world.rank() == 0) print("guess dens trace", nel);

        if (world.size() > 1) {
            LoadBalImpl<3> lb(rho, lbcost<double,3>);
            lb.load_balance();
            FunctionDefaults<3>::pmap = lb.load_balance();
            world.gop.fence();
            rho = copy(rho, FunctionDefaults<3>::pmap, false);
            vnuc = copy(vnuc, FunctionDefaults<3>::pmap, false);
        }

        START_TIMER;
        vector<functionT> ao   = project_ao_basis(world);
        END_TIMER("project ao basis");
        START_TIMER;
        Tensor<double> overlap = matrix_inner(world,ao,ao,true);
        END_TIMER("make overlap");
        START_TIMER;
        Tensor<double> kinetic = kinetic_energy_matrix(world, ao);
        END_TIMER("make KE matrix");

        functionT vlocal;
        if (param.nalpha+param.nbeta > 1) {
            START_TIMER;
            vlocal = vnuc + apply(*coulop, rho) * 0.9;            
            END_TIMER("guess Coulomb potn");
            vlocal.truncate();
        }
        else {
            vlocal = vnuc;
        }

        vlocal.reconstruct();
        reconstruct(world, ao);
        START_TIMER;
        //vector<functionT> vpsi = mul(world, vlocal, ao);
        vector<functionT> vpsi = mul_sparse(world, vlocal, ao, vtol);
        END_TIMER("make V*psi");
        truncate(world, vpsi);

        Tensor<double> potential = matrix_inner(world, vpsi, ao, true); 
        END_TIMER("make PE matrix");
        vpsi.clear(); 
        world.gop.fence();

        Tensor<double> fock = kinetic + potential;
        fock = 0.5*(fock + transpose(fock));

        Tensor<double> c,e;
        sygv(fock, overlap, 1, &c, &e);

        if (world.rank() == 0) {
//             print("THIS iS THE OVERLAP MATRIX");
//             print(overlap);
//             print("THIS iS THE KINETIC MATRIX");
//             print(kinetic);
//             print("THIS iS THE POTENTIAL MATRIX");
//             print(potential);
//             print("THESE ARE THE EIGENVECTORS");
//             print(c);
             print("THESE ARE THE INITIAL EIGENVALUES");
             print(e);
         }

        compress(world,ao);
        world.gop.fence();
        START_TIMER;
        amo = transform(world, ao, c(_,Slice(0,param.nmo_alpha-1)));
        world.gop.fence();
        END_TIMER("transform initial aMOs");
        truncate(world, amo);
        normalize(world, amo);

        aeps = e(Slice(0,param.nmo_alpha-1));

        aocc = Tensor<double>(param.nalpha);
        for (int i=0; i<param.nalpha; i++) aocc[i] = 1.0;

        if (param.nbeta && !param.spin_restricted) {
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
                //vector<functionT> psif = mul(world, psi[i], f);
                vector<functionT> psif = mul_sparse(world, psi[i], f, vtol);
                set_thresh(world, psif, vtol);  //<<<<<<<<<<<<<<<<<<<<<<<<< since cannot yet put in apply

                truncate(world,psif);
                psif = apply(world, *coulop, psif);
                truncate(world, psif);

                //psif = mul(world, psi[i], psif);
                psif = mul_sparse(world, psi[i], psif, vtol);

                gaxpy(world, 1.0, Kf, occ[i], psif);
            }
        }
        truncate(world, Kf, vtol);
        return Kf;
    }
        
    vector<functionT> 
    apply_potential(World& world, 
                    const Tensor<double>& occ,
                    const vector<functionT>& psi, 
                    const functionT& vlocal) 
    {
        //vector<functionT> Vpsi = mul(world, vlocal, psi);
        vector<functionT> Vpsi = mul_sparse(world, vlocal, psi, vtol);
        vector<functionT> Kpsi = apply_hf_exchange(world, occ, psi, psi);
        gaxpy(world, 1.0, Vpsi, -1.0, Kpsi);
        Kpsi.clear();
        truncate(world,Vpsi);
        world.gop.fence(); // ensure memory is cleared
        return Vpsi;
    }

    /// Updates the orbitals and eigenvalues of one spin
    void update(World& world, 
                Tensor<double>& occ, 
                Tensor<double>& eps,
                vector<functionT>& psi,
                vector<functionT>& Vpsi) 
    {
        int nmo = psi.size();
        vector<double> fac(nmo,-2.0);
        scale(world, Vpsi, fac);
        
        vector<poperatorT> ops = make_bsh_operators(world, eps);
        set_thresh(world, Vpsi, FunctionDefaults<3>::thresh);  // <<<<< Since cannot set in apply
        
        START_TIMER;
        vector<functionT> new_psi = apply(world, ops, Vpsi);
        END_TIMER("Apply BSH");

        ops.clear();            // free memory
        
        vector<functionT> r = sub(world, psi, new_psi); // residuals
        vector<double> rnorm = norm2(world, r);
        vector<double> new_norm = norm2(world, new_psi);
        Tensor<double> Vpr = inner(world, Vpsi, r); // Numerator of delta_epsilon
        Tensor<double> deps(nmo);
        for (int i=0; i<nmo; i++) deps(i) = 0.5*Vpr(i)/(new_norm[i]*new_norm[i]);

        Vpsi.clear(); 
        normalize(world, new_psi);

        Tensor<double> new_eps(nmo);

        for (int i=0; i<nmo; i++) {
            double step = (rnorm[i] < param.maxrotn) ? 1.0 : rnorm[i]/param.maxrotn;
            if (step!=1.0 && world.rank()==0) {
                print("  restricting step for orbital ", i, step);
            }
            psi[i].gaxpy(1.0-step, new_psi[i], step, false);
            double dd = deps[i]*step;
            if (abs(dd) > abs(0.1*eps[i])) dd = (dd/abs(dd))*0.1*abs(eps[i]);
            new_eps[i] = eps[i] + dd;
            if (new_eps[i] > -0.05) new_eps = -0.05;
        }


        if (world.rank() == 0) {
            print("rnorms");
            print(rnorm);
            //print(deps);
            //print(new_eps);
        }

        eps = new_eps;


        world.gop.fence();
        new_psi.clear(); // free memory

        truncate(world, psi);

        START_TIMER;
        // Orthog the new orbitals using sqrt(overlap).
        // Try instead an energy weighted orthogonalization
        Tensor<double> c = energy_weighted_orthog(matrix_inner(world, psi, psi, true), eps);
        psi = transform(world, psi, c);
        truncate(world, psi);
        normalize(world, psi);
        END_TIMER("Eweight orthog");

        return;
    }

    /// Diagonalizes the fock matrix and also returns energy contribution from this spin

    /// Also computes energy contributions for this spin (consistent with the input 
    /// orbitals not the output)
    void diag_fock_matrix(World& world, 
                          vector<functionT>& psi, 
                          vector<functionT>& Vpsi,
                          Tensor<double>& occ,
                          Tensor<double>& evals,
                          double& ekinetic, double& enuclear, double& etwo)
    {
        // This is unsatisfactory for large molecules, but for now will have to suffice.
        Tensor<double> overlap = matrix_inner(world, psi, psi, true);
        Tensor<double> ke = kinetic_energy_matrix(world, psi);
        Tensor<double> pe = matrix_inner(world, Vpsi, psi, true);

        //vector<functionT> vnucpsi = mul(world, vnuc, psi, vtol);
        vector<functionT> vnucpsi = mul_sparse(world, vnuc, psi, vtol);
        truncate(world, vnucpsi);
        Tensor<double> en = inner(world, vnucpsi, psi);
        vnucpsi.clear();

        int nocc = occ.size;
        ekinetic = enuclear = etwo = 0.0;
        for (int i=0; i<nocc; i++) {
            ekinetic += occ[i]*ke(i,i);
            enuclear += occ[i]*en(i);
            etwo += occ[i]*pe(i,i);
        }
        etwo -= enuclear;
        etwo *= 0.5;

        en = Tensor<double>();

        Tensor<double> fock = ke + pe;
        pe = Tensor<double>();

        fock.gaxpy(0.5,transpose(fock),0.5);

        Tensor<double> c;
        sygv(fock, overlap, 1, &c, &evals);

        Vpsi = transform(world, Vpsi, c);
        psi = transform(world, psi, c);

        truncate(world, psi);
        normalize(world, psi);
    }

    void loadbal(World& world) {
        if (world.size() == 1) return;
        LoadBalImpl<3> lb(vnuc, lbcost<double,3>);
        for (unsigned int i=0; i<amo.size(); i++) {
            lb.add_tree(amo[i],lbcost<double,3>);
        }
        if (param.nbeta && !param.spin_restricted) {
            for (unsigned int i=0; i<bmo.size(); i++) {
                lb.add_tree(bmo[i],lbcost<double,3>);
            }
        }
        lb.load_balance();
        FunctionDefaults<3>::pmap = lb.load_balance();
        world.gop.fence();
        vnuc = copy(vnuc, FunctionDefaults<3>::pmap, false);
        for (unsigned int i=0; i<amo.size(); i++) {
            amo[i] = copy(amo[i], FunctionDefaults<3>::pmap, false);
        }
        if (param.nbeta && !param.spin_restricted) {
            for (unsigned int i=0; i<bmo.size(); i++) {
                bmo[i] = copy(bmo[i], FunctionDefaults<3>::pmap, false);
            }
        }
        world.gop.fence();
    }

    void solve(World& world) {
        functionT arho_old, brho_old;
        for (int iter=0; iter<param.maxiter; iter++) {
            if (world.rank()==0) print("\nIteration", iter,"\n");
            START_TIMER;
            loadbal(world);
            if (iter > 0) {
               arho_old = copy(arho_old, FunctionDefaults<3>::pmap, false);
               if (param.nbeta && !param.spin_restricted) {
                   brho_old = copy(brho_old, FunctionDefaults<3>::pmap, false);
               }
            }
            END_TIMER("Load balancing");

            START_TIMER;
            functionT arho = make_density(world, aocc, amo);
            functionT brho = factoryT(world);
            if (param.nbeta) {
                if (param.spin_restricted) brho = arho;
                else brho = make_density(world, bocc, bmo);
            }
            END_TIMER("Make densities");

            double da, db;
            if (iter > 0) {
                da = (arho - arho_old).norm2();
                db = 0.0;
                if (param.nbeta && !param.spin_restricted) db = (brho - brho_old).norm2();
                if (world.rank()==0) print("delta rho", da, db);
            }
            arho_old = arho;
            brho_old = brho;

            functionT rho = arho+brho;
            rho.truncate();
            START_TIMER;
            functionT vlocal = apply(*coulop, rho) + vnuc;
            END_TIMER("Coulomb");
            rho.clear(false);
            vlocal.truncate(); // For DFT must add exchange-correlation into vlocal

            START_TIMER;
            vector<functionT> Vpsia = apply_potential(world, aocc, amo, vlocal);
            vector<functionT> Vpsib;
            if (param.nbeta && !param.spin_restricted) Vpsib = apply_potential(world, bocc, bmo, vlocal);
            END_TIMER("Apply potential");

            double ekinetic, enuclear, etwo;
            START_TIMER;
            diag_fock_matrix(world, amo, Vpsia, aocc, aeps, ekinetic, enuclear, etwo);
            if (param.spin_restricted) {
                ekinetic*=2;
                enuclear*=2;
                etwo*=2;
            }
            else if (param.nbeta) {
                double ekineticb, enuclearb, etwob;
                diag_fock_matrix(world, bmo, Vpsib, bocc, beps, ekineticb, enuclearb, etwob);
                ekinetic+=ekineticb; enuclear+=enuclearb; etwo+=etwob;
            }
            END_TIMER("Diag and transform");

            if (world.rank() == 0) {
                double enrep = molecule.nuclear_repulsion_energy();
                double etot = ekinetic + enuclear + etwo + enrep;
                print("              kinetic", ekinetic);
                print("   nuclear attraction", enuclear);
                print("         two-electron", etwo);
                print("    nuclear-repulsion", enrep);
                print("                total", etot);
                print(" ");
                print("alpha evals");
                print(aeps);
                if (param.nbeta && !param.spin_restricted) {
                    print("beta evals");
                    print(beps);
                }
            }

            if (iter > 0) {
                double dconv = max(FunctionDefaults<3>::thresh, param.dconv);
                if (da<dconv && db<dconv) {
                    if (world.rank()==0) print("\nConverged!\n");
                    return;
                }
            }


            //update(world, vlocal, aocc, aeps, amo, Vpsia);
            //if (!param.spin_restricted) update(world, vlocal, bocc, beps, bmo, Vpsib);
            
            update(world, aocc, aeps, amo, Vpsia); 
            if (param.nbeta && !param.spin_restricted) update(world, bocc, beps, bmo, Vpsib);

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
        calc.set_protocol(world,1e-4);
        calc.make_nuclear_potential(world);
        calc.initial_guess(world);
        calc.solve(world);

        calc.set_protocol(world,1e-6);
        calc.make_nuclear_potential(world);
        calc.project(world);
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
    } catch (char* s) {
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

