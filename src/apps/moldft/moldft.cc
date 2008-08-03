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
#include <tensor/random.h>
#include <ctime>
using namespace madness;

extern int x_rks_s__(const double *rho, double *f, double *dfdra);
extern int c_rks_vwn5__(const double *rho, double *f, double *dfdra);
// extern int x_uks_s__(const double *rho, double *f, double *dfdra);
// extern int c_uks_vwn5__(const double *rho, double *f, double *dfdra);

#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>


/// Simplified (?) version of BLAS-1 DROT(N, DX, INCX, DY, INCY, DC, DS)
void drot(long n, double* RESTRICT a, double* RESTRICT b, double s, double c, long inc) {
    if (inc == 1) {
        for (long i=0; i<n; i++) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
    else {
        for (long i=0; i<(n*inc); i+=inc) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
}

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

typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef vector<functionT> vecfuncT;
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

class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coordT& x) const {
        return x[axis];
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
    bool lda;                   ///< True if LDA (HF if false)
    int plotlo,plothi;          ///< Range of MOs to print (for both spins if polarized)
    bool plotdens;              ///< If true print the density at convergence
    bool plotcoul;              ///< If true plot the total coulomb potential at convergence
    // Next list inferred parameters
    int nalpha;                 ///< Number of alpha spin electrons
    int nbeta;                  ///< Number of beta  spin electrons
    int nmo_alpha;              ///< Number of alpha spin molecular orbitals
    int nmo_beta;               ///< Number of beta  spin molecular orbitals
    double lo;                  ///< Smallest length scale we need to resolve

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & L & maxrotn & nvalpha & nvbeta & nopen & maxiter & spin_restricted & lda;
        ar & plotlo & plothi & plotdens & plotcoul;
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
        , lda(true)
        , plotlo(0)
        , plothi(-1)
        , plotdens(false)
        , plotcoul(false)
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
            } else if (s == "lda") {
                lda = true;
            } else if (s == "hf") {
                lda = false;
            } else if (s == "plotmos") {
                f >> plotlo >> plothi;
            } else if (s == "plotdens") {
                plotdens = true;
            } else if (s == "plotcoul") {
                plotcoul = true;
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
        nvalpha = nmo_alpha - nalpha;
        nvbeta = nmo_beta - nbeta;
        
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
        const char* calctype[2] = {"Hartree-Fock","LDA"};
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
        madness::print("    calculation type ", calctype[int(lda)]);
    }
    
};

struct Calculation {
    Molecule molecule;            ///< Molecular coordinates, etc.
    CalculationParameters param;  ///< User input data, nalpha, etc.
    AtomicBasisSet aobasis;       ///< Currently always the STO-3G basis
    functionT vnuc;               ///< The effective nuclear potential
    vecfuncT amo, bmo;   ///< alpha and beta molecular orbitals
    Tensor<double> aocc, bocc;    ///< alpha and beta occupation numbers
    Tensor<double> aeps, beps;    ///< alpha and beta energy shifts (eigenvalues if canonical)
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

        FunctionDefaults<3>::set_cubic_cell(-param.L,param.L);
        
        // Setup initial defaults for numerical functions
        set_protocol(world, 1e-4);
    }

    void set_protocol(World& world, double thresh) {
        FunctionDefaults<3>::set_thresh(thresh);
        if (thresh >= 1e-2) FunctionDefaults<3>::set_k(4);
        else if (thresh >= 1e-4) FunctionDefaults<3>::set_k(6);
        else if (thresh >= 1e-6) FunctionDefaults<3>::set_k(8);
        else if (thresh >= 1e-8) FunctionDefaults<3>::set_k(10);
        else FunctionDefaults<3>::set_k(12);

        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_initial_level(2);
        FunctionDefaults<3>::set_truncate_mode(1);  
        FunctionDefaults<3>::set_autorefine(false);  

        double safety = 0.1;
        vtol = FunctionDefaults<3>::get_thresh()*safety;

        coulop = poperatorT(CoulombOperatorPtr<double, 3>(world, 
                                                          FunctionDefaults<3>::get_k(), 
                                                          param.lo, 
                                                          vtol/safety)); // No need for safety here!
        if (world.rank() == 0) {
            print("\nSolving with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), 
                  "   conv", max(thresh, param.dconv), "\n");
        }
    }

    void project(World& world) {
        reconstruct(world,amo);
        for (unsigned int i=0; i<amo.size(); i++) {
            amo[i] = madness::project(amo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
        }
        world.gop.fence();
        truncate(world,amo);
        normalize(world,amo);
        if (param.nbeta && !param.spin_restricted) {
            reconstruct(world,bmo);
            for (unsigned int i=0; i<bmo.size(); i++) {
                bmo[i] = madness::project(bmo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
            }
            world.gop.fence();
            truncate(world,bmo);
            normalize(world,bmo);
        }
    }
        
    void make_nuclear_potential(World& world) {
        START_TIMER;
        vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(vtol).truncate_on_project(); 
        //vnuc.truncate();
        vnuc.reconstruct();
        END_TIMER("Project vnuclear");
    }

    vecfuncT project_ao_basis(World& world) {
        vecfuncT ao(aobasis.nbf(molecule));

        Level initial_level = 3;
        for (int i=0; i<aobasis.nbf(molecule); i++) {
            functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule,i)));
//             if (world.rank() == 0) {
//                 aobasis.get_atomic_basis_function(molecule,i).print_me(cout);
//             }
            ao[i] = factoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
        }
        world.gop.fence();
        //truncate(world, ao);

        vector<double> norms;

        while (1) {
            norms = norm2(world, ao);
            initial_level += 2;
            if (initial_level >= 11) throw "project_ao_basis: projection failed?";
            int nredone = 0;
            for (int i=0; i<aobasis.nbf(molecule); i++) {
                if (norms[i] < 0.99) {
                    nredone++;
                    //if (world.rank() == 0) print("re-projecting ao basis function", i,"at level",initial_level);
                    functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule,i)));
                    ao[i] = factoryT(world).functor(aofunc).initial_level(6).truncate_on_project().nofence();
                }
            }
            world.gop.fence();
            if (nredone == 0) break;
        }

        for (int i=0; i<aobasis.nbf(molecule); i++) {
            if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
            norms[i] = 1.0/norms[i];
        }

        scale(world, ao, norms);
        
        return ao;
    }

    void localize(World& world, vecfuncT& mo) {
        // Form the ao-mo overlap matrix
        
        truncate(world,mo);

        vecfuncT ao = project_ao_basis(world);
        Tensor<double> S = matrix_inner(world, mo, ao);
        Tensor<double> Sorig = copy(S);

        long nmo = mo.size();
        long nao = ao.size();

        Tensor<double> U(nmo,nmo);
        for (long i=0; i<nmo; i++) U(i,i) = 1.0;

        print("S before");
        print(S);

        const double thresh = 1e-9; // Final convergence test
        double tol = 0.1;           // Current test
        long ndone=0;               // Total no. of rotations performed


        for (long iter=0; iter<100; iter++) {
            // Compute the objective function to track convergence
            double sum = 0.0;
            for (long i=0; i<nmo; i++) {
                for (long mu=0; mu<nao; mu++) {
                    double si = S(i,mu);
                    sum += si*si*si*si;
                }
            }

            //tol *= 0.25;
            tol *= 0.5;
            if (tol < thresh) tol = thresh;
            long ndone_iter=0; // No. of rotations done this iteration
            double screen = tol*sum;

            print("iteration", iter, sum, ndone);
            for (long i=0; i<nmo; i++) {
                for (long j=0; j<i; j++) {
                    double g = 0.0;
                    double h = 0.0;
                    double sij = 0.0; // Estimate of the overlap of |i> and |j>
                    for (long mu=0; mu<nao; mu++) {
                        double si = S(i,mu);
                        double sj = S(j,mu);
                        g += si*sj*sj*sj - si*si*si*sj;
                        h += 6*si*si*sj*sj - si*si*si*si - sj*sj*sj*sj;
                        sij += si*si*sj*sj;
                    }

                    double thetamax = 0.5;  // circa 30 degrees
                    if (h == 0.0) {
                        //print("     exactly zero h", i, j, h);
                        h = -1e-20;
                        thetamax *= 0.5;
                    }
                    else if (h > 0.0) {
                        //print("     forcing negative h", i, j, h);
                        h = -h;
                        thetamax *= 0.5;
                    }

                    double theta = -g/h;

                    if (fabs(theta) > thetamax) {
                        //print("  restricting", i, j);
                        if (g < 0) theta = -thetamax;
                        else       theta =  thetamax*0.8; // Breaking symmetry is good
                    }

                    if (iter == 0 && sij > 0.01 && theta < 0.01) {
                        //print("   randomizing", i, j, h, g, sij, theta);
                        theta += (RandomNumber<double>() - 0.5);
                    }
                    
                    if (fabs(theta) > screen) {
                        ndone_iter++;
                        //print("    ", i, j, g, h, theta);
                        
                        double c = cos(theta);
                        double s = sin(theta);
                        drot(nao, &S(i,0), &S(j,0), s, c, 1);
                        drot(nmo, &U(i,0), &U(j,0), s, c, 1);
                    }
                }
            }
            ndone += ndone_iter;
            if (ndone_iter==0 && tol==thresh) {
                print("Converged!", ndone);
                break;
            }
        }
        U = transpose(U);

        print(S);

        print("err in SU", (S-inner(transpose(U),Sorig)).normf());

        vecfuncT monew = transform(world, mo, U);
        truncate(world,monew);

        Tensor<double> Snew = matrix_inner(world, monew, ao);

        print("err in Snew", (S-Snew).normf());
        print(Snew);

        long nold=0, nnew=0;
        for (int i=0; i<nmo; i++) {
            int nn = monew[i].size();
            int no = mo[i].size();
            nnew += nn;
            nold += no;
            print(i, no, nn);
        }
        print("totals", nold, nnew);

        throw "done";
    }

    void analyze_vectors(World& world, const vecfuncT& mo) {
        vecfuncT ao = project_ao_basis(world);
        Tensor<double> Saomo = matrix_inner(world, ao, mo);
        Tensor<double> Saoao = matrix_inner(world, ao, ao, true);

        ao.clear(); world.gop.fence(); // free memory

        Tensor<double> C;
        gesv(Saoao, Saomo, &C);
        C = transpose(C); 

        // psi_i(x) approx = sum (mu) C(i,mu)*chi_mu(x)

        long nmo = mo.size();
        for (long i=0; i<nmo; i++) {
            printf("  MO%4ld:", i);
            aobasis.print_anal(molecule, C(i,_));
        }
    }

    void localize_boys(World& world, vecfuncT& mo) {
        // Should modify this to broadcast the result from 0 to all nodes
        // just to be sure all have bitwise the same data

        const bool doprint=false;
        START_TIMER;
#define DIP(i,j,k,l) (dip(0,i,j)*dip(0,k,l)+dip(1,i,j)*dip(1,k,l)+dip(2,i,j)*dip(2,k,l))
        // Form the dipole moment matrices

        truncate(world,mo); /// <<<<<<<<<<<<<<<<<<<<<<<  FOR DEBUG OF LOCALIZATION ONLY

        long nmo = mo.size();
        Tensor<double> dip(3, nmo, nmo);

        for (int axis=0; axis<3; axis++) {
            functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
            dip(axis,_,_) = matrix_inner(world, mo, mul(world, fdip, mo), true);
        }

        Tensor<double> U(nmo,nmo);
        for (long i=0; i<nmo; i++) U(i,i) = 1.0;

        const double thresh = 1e-9; // Final convergence test
        double tol = 0.1;           // Current test
        long ndone=0;               // Total no. of rotations performed

        bool converged = false;
        for (long iter=0; iter<100; iter++) {
            // Compute the objective function to track convergence
            double sum = 0.0;
            for (long i=0; i<nmo; i++) {
                sum += DIP(i,i,i,i);
            }

            tol *= 0.5;
            if (tol < thresh) tol = thresh;
            long ndone_iter=0; // No. of rotations done this iteration
            double screen = tol*sum;

            if (doprint) print("iteration", iter, sum, ndone);
            for (long i=0; i<nmo; i++) {
                for (long j=0; j<i; j++) {
                    double g = 0.0;
                    double h = 0.0;
                    double sij = 0.0; // Estimate of the overlap of |i> and |j>
                    g += DIP(i,j,j,j) - DIP(i,j,i,i);
                    h += 4.0*DIP(i,j,i,j) + 2.0*DIP(i,i,j,j) - DIP(i,i,i,i) - DIP(j,j,j,j);
                    sij += DIP(i,j,i,j);

                    double thetamax = 0.5;  // circa 30 degrees
                    if (h == 0.0) {
                        if (doprint) print("     exactly zero h", i, j, h);
                        h = -1.0;
                        thetamax *= 0.5;
                    }
                    else if (h > 0.0) {
                        if (doprint) print("     forcing negative h", i, j, h);
                        h = -1.0;
                        thetamax *= 0.5;
                    }

                    double theta = -g/h;

                    if (fabs(theta) > thetamax) {
                        if (doprint) print("  restricting", i, j);
                        if (g < 0) theta = -thetamax;
                        else       theta =  thetamax*0.8; // Breaking symmetry is good
                    }

                    bool randomized=false;
                    if (iter == 0 && sij > 0.01 && fabs(theta) < 0.01) { //
                        randomized=true;
                        if (doprint) print("   randomizing", i, j, h, g, sij, theta);
                        theta += (RandomNumber<double>() - 0.5);
                    }
                 
                    if (fabs(theta) > screen || randomized) {
                        ndone_iter++;
                        if (doprint) print("    ", i, j, g, h, theta);
                     
                        double c = cos(theta);
                        double s = sin(theta);
                        for (int axis=0; axis<3; axis++) {
                            drot(nmo, &dip(axis,i,0), &dip(axis,j,0), s, c, 1);
                            // This next line is the expensive one due to the non-unit stride.
                            // It can be eliminated at the expense of adding more unit
                            // stride work when computing the integrals above (inside DIP).
                            drot(nmo, &dip(axis,0,i), &dip(axis,0,j), s, c, nmo);
                        }
                        drot(nmo, &U(i,0), &U(j,0), s, c, 1);
                    }
                }
            }
            ndone += ndone_iter;
            if (ndone_iter==0 && tol==thresh) {
                if (doprint) print("Converged!", ndone);
                converged = true;
                break;
            }
        }
        U = transpose(U);

        if (!converged && world.rank() == 0) {
            print("warning: boys localization did not fully converge");
        }

        if (doprint) {
            vecfuncT monew = transform(world, mo, U);
            truncate(world,monew);
            
            Tensor<double> Snew = matrix_inner(world, monew, monew);
            print("Snew");
            print(Snew);
            
            print("Original vector analysis");
            analyze_vectors(world, mo);
            print("Localized vector analysis");
            analyze_vectors(world, monew);
            
            long nold=0, nnew=0;
            for (int i=0; i<nmo; i++) {
                int nn = monew[i].size();
                int no = mo[i].size();
                nnew += nn;
                nold += no;
                print(i, no, nn);
            }
            print("totals", nold, nnew);
        }

#undef DIP
        END_TIMER("Boys localize");
    }

    Tensor<double> kinetic_energy_matrix(World& world, const vecfuncT& v) {
        reconstruct(world, v);
        int n = v.size();
        Tensor<double> r(n,n);
        for (int axis=0; axis<3; axis++) {
            vecfuncT dv = diff(world,v,axis);
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
        // We use the density to help with load balance and also
        // for computing the initial coulomb potential.
        START_TIMER;
        functionT rho = factoryT(world).functor(functorT(new GuessDensity(molecule, aobasis)));
        END_TIMER("guess density");
        if (world.rank() == 0) print("\ndistribution of rho");
        rho.print_info();
        if (world.rank() == 0) print("\ndistribution of vnuc");
        vnuc.print_info();
        double nel = rho.trace();
        if (world.rank() == 0) print("guess dens trace", nel);

//         if (world.size() > 1) {
//             LoadBalImpl<3> lb(rho, lbcost<double,3>);
//             lb.load_balance();
//             FunctionDefaults<3>::set_pmap(lb.load_balance());
//             world.gop.fence();
//             rho = copy(rho, FunctionDefaults<3>::get_pmap(), false);
//             vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//         }

        functionT vlocal;
        if (param.nalpha+param.nbeta > 1) {
            START_TIMER;
            vlocal = vnuc + apply(*coulop, rho); //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
            END_TIMER("guess Coulomb potn");

            // Shove the closed-shell LDA potential on there also
            bool save = param.spin_restricted;
            param.spin_restricted = true;
            rho.scale(0.5);
            vlocal = vlocal + make_lda_potential(world, rho, rho, functionT(), functionT());
            vlocal.truncate();
            param.spin_restricted = save;
        }
        else {
            vlocal = vnuc;
        }
        rho.clear();
        vlocal.reconstruct();

        START_TIMER;
        vecfuncT ao   = project_ao_basis(world);
        END_TIMER("project ao basis");
//         for (unsigned int i=0; i<ao.size(); i++) {
//             if (world.rank() == 0) print("\ndistribution of ao",i);
//             ao[i].print_info();
//         }

//         if (world.size() > 1) {
//             LoadBalImpl<3> lb(ao[0], lbcost<double,3>);
//             for (unsigned int i=1; i<ao.size(); i++) {
//                 lb.add_tree(ao[i],lbcost<double,3>);
//             }
//             lb.load_balance();
//             FunctionDefaults<3>::set_pmap(lb.load_balance());
//             world.gop.fence();
//             vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//             vlocal = copy(vlocal, FunctionDefaults<3>::get_pmap(), false);
//             for (unsigned int i=0; i<amo.size(); i++) {
//                 ao[i] = copy(ao[i], FunctionDefaults<3>::get_pmap(), false);
//             }
//         }

        START_TIMER;
        Tensor<double> overlap = matrix_inner(world,ao,ao,true);
        END_TIMER("make overlap");
        START_TIMER;
        Tensor<double> kinetic = kinetic_energy_matrix(world, ao);
        END_TIMER("make KE matrix");

        reconstruct(world, ao);
        START_TIMER;
        //vecfuncT vpsi = mul(world, vlocal, ao);
        vecfuncT vpsi = mul_sparse(world, vlocal, ao, vtol);
        world.gop.fence();
        END_TIMER("make V*psi");
        world.gop.fence();
        START_TIMER;
        compress(world,vpsi);
        END_TIMER("Compressing Vpsi");
        START_TIMER;
        truncate(world, vpsi);
        END_TIMER("Truncating vpsi");
        START_TIMER;
        compress(world, ao);
        END_TIMER("Compressing AO");

        START_TIMER;
        Tensor<double> potential = matrix_inner(world, vpsi, ao, true); 
        world.gop.fence();
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
            print("initial eigenvalues");
             print(e);
         }

        compress(world,ao);
        world.gop.fence();
        START_TIMER;
        amo = transform(world, ao, c(_,Slice(0,param.nmo_alpha-1)));
        world.gop.fence();
        END_TIMER("transform initial MO");
        truncate(world, amo);
        normalize(world, amo);

        aeps = e(Slice(0,param.nmo_alpha-1));

        aocc = Tensor<double>(param.nmo_alpha);
        for (int i=0; i<param.nalpha; i++) aocc[i] = 1.0;

        print("AOCC", aocc);

        if (param.nbeta && !param.spin_restricted) {
            bmo = transform(world, ao, c(_,Slice(0,param.nmo_beta-1)));
            truncate(world, bmo);
            normalize(world, bmo);
            beps = e(Slice(0,param.nmo_beta-1));
            bocc = Tensor<double>(param.nmo_beta);
            for (int i=0; i<param.nbeta; i++) bocc[i] = 1.0;
        }
        
        localize_boys(world, amo);
    }

    functionT make_density(World& world, const Tensor<double>& occ, const vecfuncT& v) {
        vecfuncT vsq = square(world, v);
        compress(world,vsq);
        functionT rho = factoryT(world).thresh(vtol);
        rho.compress();
        for (unsigned int i=0; i<vsq.size(); i++) {
            if (occ[i]) rho.gaxpy(1.0,vsq[i],occ[i],false);
        }
        world.gop.fence();
        vsq.clear();
        world.gop.fence();

        double dtrace = rho.trace();
        if (world.rank() == 0) print("err in trace of density", dtrace-occ.sum());
        return rho;
    }

    vector<poperatorT> make_bsh_operators(World& world, const Tensor<double>& evals) {
        int nmo = evals.dim[0];
        vector<poperatorT> ops(nmo);
        int k = FunctionDefaults<3>::get_k();
        double tol = FunctionDefaults<3>::get_thresh();
        for (int i=0; i<nmo; i++) {
            double eps = evals(i);
            if (eps > 0) {
                if (world.rank() == 0) {
                    print("bsh: warning: positive eigenvalue", i, eps);
                }
                eps = -0.1;
            }
            ops[i] = poperatorT(BSHOperatorPtr<double,3>(world, sqrt(-2.0*eps), k, param.lo, tol));
        }
        return ops;
    }

    vecfuncT apply_hf_exchange(World& world, 
                                        const Tensor<double>& occ,
                                        const vecfuncT& psi,
                                        const vecfuncT& f) {
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

        vecfuncT Kf = zero_functions<double,3>(world, nf);

        compress(world,Kf);
        reconstruct(world,psi);
        for (int i=0; i<nocc; i++) {
            if (occ[i] > 0.0) {
                //vecfuncT psif = mul(world, psi[i], f);
                vecfuncT psif = mul_sparse(world, psi[i], f, vtol);
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
        
    static double munge(double r) {
        if (r < 1e-12) r = 1e-12;
        return r;
    }

    static void ldaop(const Key<3>& key, Tensor<double>& t) {
        UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
    }

    static void ldaeop(const Key<3>& key, Tensor<double>& t) {
        UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q1; double q2; double dq;x_rks_s__(&r, &q1, &dq);c_rks_vwn5__(&r, &q2, &dq); *_p0 = q1+q2);
    }


    functionT
    make_lda_potential(World& world, 
                       const functionT& arho,
                       const functionT& brho,
                       const functionT& adelrhosq,
                       const functionT& bdelrhosq)
    {
        MADNESS_ASSERT(param.spin_restricted);
        functionT vlda = copy(arho);
        vlda.reconstruct();
        vlda.unaryop(&ldaop);
        return vlda;
    }

    double
    make_lda_energy(World& world, 
                    const functionT& arho,
                    const functionT& brho,
                    const functionT& adelrhosq,
                    const functionT& bdelrhosq)
    {
        MADNESS_ASSERT(param.spin_restricted);
        functionT vlda = copy(arho);
        vlda.reconstruct();
        vlda.unaryop(&ldaeop);
        return vlda.trace();
    }

    vecfuncT 
    apply_potential(World& world, 
                    const Tensor<double>& occ,
                    const vecfuncT& amo, 
                    const functionT& arho,
                    const functionT& brho,
                    const functionT& adelrhosq,
                    const functionT& bdelrhosq,
                    const functionT& vlocal,
                    double& exc) 
    {

        functionT vloc = vlocal;
        if (param.lda) {
            exc = make_lda_energy(world, arho, brho, adelrhosq, bdelrhosq);
            vloc = vloc + make_lda_potential(world, arho, brho, adelrhosq, bdelrhosq);
        }

        vecfuncT Vpsi = mul_sparse(world, vloc, amo, vtol);

        if (!param.lda) {
            vecfuncT Kamo = apply_hf_exchange(world, occ, amo, amo);
            Tensor<double> excv = inner(world, Kamo, amo);
            exc = 0.0;
            for (unsigned long i=0; i<amo.size(); i++) {
                exc -= 0.5*excv[i]*occ[i];
            }
            gaxpy(world, 1.0, Vpsi, -1.0, Kamo);
            Kamo.clear();
        }

        truncate(world,Vpsi);
        world.gop.fence(); // ensure memory is cleared
        return Vpsi;
    }


    /// Updates the orbitals and eigenvalues of one spin
    void update(World& world, 
                Tensor<double>& occ, 
                Tensor<double>& eps,
                vecfuncT& psi,
                vecfuncT& Vpsi) 
    {
        int nmo = psi.size();
        vector<double> fac(nmo,-2.0);
        scale(world, Vpsi, fac);
        
        vector<poperatorT> ops = make_bsh_operators(world, eps);
        set_thresh(world, Vpsi, FunctionDefaults<3>::get_thresh());  // <<<<< Since cannot set in apply
        
        START_TIMER;
        vecfuncT new_psi = apply(world, ops, Vpsi);
        END_TIMER("Apply BSH");

        ops.clear();            // free memory
        
        vecfuncT r = sub(world, psi, new_psi); // residuals
        vector<double> rnorm = norm2(world, r);
        vector<double> new_norm = norm2(world, new_psi);
        Tensor<double> Vpr = inner(world, Vpsi, r); // Numerator of delta_epsilon
        Tensor<double> deps(nmo);
        for (int i=0; i<nmo; i++) deps(i) = 0.5*Vpr(i)/(new_norm[i]*new_norm[i]);

        Vpsi.clear(); 
        normalize(world, new_psi);

        Tensor<double> new_eps(nmo);

        for (int i=0; i<nmo; i++) {
            double step = (rnorm[i] < param.maxrotn) ? 1.0 : param.maxrotn;
            if (step!=1.0 && world.rank()==0) {
                print("  restricting step for orbital ", i, step);
            }
            psi[i].gaxpy(1.0-step, new_psi[i], step, false);
            double dd = deps[i]*step;
            if (abs(dd) > abs(0.1*eps[i])) dd = (dd/abs(dd))*0.1*abs(eps[i]);
            new_eps[i] = eps[i] + dd;
            if (new_eps[i] > -0.05) new_eps = -0.05;
        }
        world.gop.fence();


        if (world.rank() == 0) {
            print("residual norms");
            print(rnorm);
        }

        eps = new_eps;

        world.gop.fence();
        new_psi.clear(); // free memory

        truncate(world, psi);

//         START_TIMER;
//         // Orthog the new orbitals using sqrt(overlap).
//         // Try instead an energy weighted orthogonalization
//         Tensor<double> c = energy_weighted_orthog(matrix_inner(world, psi, psi, true), eps);
//         psi = transform(world, psi, c);
//         truncate(world, psi);
//         normalize(world, psi);
//         END_TIMER("Eweight orthog");

        return;
    }

    /// Diagonalizes the fock matrix and also returns energy contribution from this spin

    /// Also computes energy contributions for this spin (consistent with the input 
    /// orbitals not the output)
    void diag_fock_matrix(World& world, 
                          vecfuncT& psi, 
                          vecfuncT& Vpsi,
                          Tensor<double>& occ,
                          Tensor<double>& evals,
                          const functionT& arho, 
                          const functionT& brho, 
                          const functionT& adelrhosq, 
                          const functionT& bdelrhosq, 
                          double& ekinetic)
    {
        // This is unsatisfactory for large molecules, but for now will have to suffice.
        Tensor<double> overlap = matrix_inner(world, psi, psi, true);
        Tensor<double> ke = kinetic_energy_matrix(world, psi);
        Tensor<double> pe = matrix_inner(world, Vpsi, psi, true);

        int nocc = occ.size;
        ekinetic = 0.0;
        for (int i=0; i<nocc; i++) {
            ekinetic += occ[i]*ke(i,i);
        }

        Tensor<double> fock = ke + pe;
        pe = Tensor<double>();

        fock.gaxpy(0.5,transpose(fock),0.5);

        Tensor<double> c;
        sygv(fock, overlap, 1, &c, &evals);

        print("NEW EVALS", evals);

        Vpsi = transform(world, Vpsi, c);
        psi = transform(world, psi, c);

        truncate(world, psi);
        normalize(world, psi);
    }

    void loadbal(World& world) {
//         if (world.size() == 1) return;
//         LoadBalImpl<3> lb(vnuc, lbcost<double,3>);
//         for (unsigned int i=0; i<amo.size(); i++) {
//             lb.add_tree(amo[i],lbcost<double,3>);
//         }
//         if (param.nbeta && !param.spin_restricted) {
//             for (unsigned int i=0; i<bmo.size(); i++) {
//                 lb.add_tree(bmo[i],lbcost<double,3>);
//             }
//         }
//         lb.load_balance();
//         FunctionDefaults<3>::set_pmap(lb.load_balance());
//         world.gop.fence();
//         vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//         for (unsigned int i=0; i<amo.size(); i++) {
//             amo[i] = copy(amo[i], FunctionDefaults<3>::get_pmap(), false);
//         }
//         if (param.nbeta && !param.spin_restricted) {
//             for (unsigned int i=0; i<bmo.size(); i++) {
//                 bmo[i] = copy(bmo[i], FunctionDefaults<3>::get_pmap(), false);
//             }
//         }
//         world.gop.fence();
    }

    void solve(World& world) {
        functionT arho_old, brho_old;
        functionT adelrhosq, bdelrhosq; // placeholders for GGAs

        for (int iter=0; iter<param.maxiter; iter++) {
            if (world.rank()==0) print("\nIteration", iter,"\n");
            START_TIMER;
            loadbal(world);
            if (iter > 0) {
                arho_old = copy(arho_old, FunctionDefaults<3>::get_pmap(), false);
               if (param.nbeta && !param.spin_restricted) {
                   brho_old = copy(brho_old, FunctionDefaults<3>::get_pmap(), false);
               }
            }
            END_TIMER("Load balancing");

            START_TIMER;
            functionT arho = make_density(world, aocc, amo);
            functionT brho;
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

            double enuclear = inner(rho, vnuc);

            START_TIMER;
            functionT vcoul = apply(*coulop, rho);
            END_TIMER("Coulomb");

            double ecoulomb = 0.5*inner(rho, vcoul);
            rho.clear(false);

            functionT vlocal = vcoul + vnuc;
            vcoul.clear(false);
            vlocal.truncate(); 

            START_TIMER;
            double exca=0.0, excb=0.0;
            vecfuncT Vpsia = apply_potential(world, aocc, amo, arho, brho, adelrhosq, bdelrhosq, vlocal, exca);
            vecfuncT Vpsib;
            if (param.spin_restricted) {
                if (!param.lda) excb = exca;
            }
            else if (param.nbeta) {
                Vpsib = apply_potential(world, bocc, bmo, brho, arho, bdelrhosq, adelrhosq, vlocal, excb);
            }
            END_TIMER("Apply potential");

            double ekina=0.0, ekinb=0.0;;
            START_TIMER;
            diag_fock_matrix(world, amo, Vpsia, aocc, aeps, arho, brho, adelrhosq, bdelrhosq, ekina);
            if (param.spin_restricted) {
                ekinb = ekina;
            }
            else if (param.nbeta) {
                diag_fock_matrix(world, bmo, Vpsib, bocc, beps, brho, arho, bdelrhosq, adelrhosq, ekinb);
            }
            END_TIMER("Diag and transform");

            if (world.rank() == 0) {
                double enrep = molecule.nuclear_repulsion_energy();
                double ekinetic = ekina + ekinb;
                double exc = exca + excb;
                double etot = ekinetic + enuclear + ecoulomb + exc + enrep;
                printf("\n              kinetic %16.8f\n", ekinetic);
                printf("   nuclear attraction %16.8f\n", enuclear);
                printf("              coulomb %16.8f\n", ecoulomb);
                printf(" exchange-correlation %16.8f\n", exc);
                printf("    nuclear-repulsion %16.8f\n", enrep);
                printf("                total %16.8f\n\n", etot);
            }

            if (iter > 0) {
                double dconv = max(FunctionDefaults<3>::get_thresh(), param.dconv);
                if (da<dconv && db<dconv) {
                    if (world.rank()==0) {
                        print("\nConverged!\n");
                        print(" ");
                        print("alpha eigenvalues");
                        print(aeps);
                        if (param.nbeta && !param.spin_restricted) {
                            print("beta eigenvalues");
                            print(beps);
                        }
                    }

                    return;
                }
            }

            update(world, aocc, aeps, amo, Vpsia); 
            if (param.nbeta && !param.spin_restricted) update(world, bocc, beps, bmo, Vpsib);

        }
    }

};

//#include <sched.h>

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
//     cpu_set_t mask;
//     CPU_ZERO(&mask);
//     CPU_SET(MPI::COMM_WORLD.Get_rank(), &mask);
//     if( sched_setaffinity( 0, sizeof(mask), &mask ) == -1 ) {
//         printf("WARNING: Could not set CPU Affinity, continuing...\n");
//     }

    World world(MPI::COMM_WORLD);

    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));
    
        
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

