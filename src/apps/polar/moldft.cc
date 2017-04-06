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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <type_traits>
#include <madness/mra/mra.h>

#include <polar/molecule.h>
#include <polar/molecularbasis.h>
#include <polar/corepotential.h>
#include <polar/xcfunctional.h>
#include <polar/potentialmanager.h>

#include <madness/tensor/elem.h>
#include <madness/mra/lbdeux.h>
#include <madness/mra/qmprop.h>
#include <madness/misc/misc.h>
#include <madness/misc/ran.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/distributed_matrix.h>
#include <madness/world/worldmem.h>

using namespace madness;

//#include <jacob/abinitdftsolventsolver.h>
//#include <examples/molecularmask.h>

template<int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<double_complex>& t) const {
        //vzExp(t.size, t.ptr(), t.ptr());
        UNARY_OPTIMIZED_ITERATOR(double_complex, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

// Returns exp(-I*t*V)
Function<double_complex,3> make_exp(double t, const Function<double,3>& v) {
    v.reconstruct();
    Function<double_complex,3> expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<3>());
    //expV.truncate(); expV.reconstruct();
    return expV;
}

extern void drot(long n, double* restrict a, double* restrict b, double s, double c, long inc);

void drot3(long n, double* restrict a, double* restrict b, double s, double c, long inc) {
    if (inc == 1) {
        n*=3;
        for (long i=0; i<n; i+=3) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
    else {
        inc*=3;
        n*=inc;
        for (long i=0; i<n; i+=inc) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
}

// class NuclearDensityFunctor : public FunctionFunctorInterface<double,3> {
//   Molecule molecule;
//   std::vector<coord_3d> specialpts;
// public:
//   NuclearDensityFunctor(const Molecule& molecule) : molecule(molecule) {}

//   double operator()(const Vector<double,3>& r) const {
//     return molecule.mol_nuclear_charge_density(r[0], r[1], r[2]);
//   }

//   std::vector<coord_3d> special_points() const{
//     return molecule.get_all_coords_vec();
//   }

//   Level special_level() {
//     return 15;
//   }
// };


typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>,3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;

extern tensorT distributed_localize_PM(World & world, 
                                const vecfuncT & mo, 
                                const vecfuncT & ao, 
                                const std::vector<int> & set, 
                                const std::vector<int> & at_to_bf,
                                const std::vector<int> & at_nbf, 
                                const double thresh = 1e-9, 
                                const double thetamax = 0.5, 
                                const bool randomize = true, 
                                       const bool doprint = false);

static double ttt, sss;
void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}


inline double mask1(double x) {
    /* Iterated first beta function to switch smoothly
       from 0->1 in [0,1].  n iterations produce 2*n-1
       zero derivatives at the end points. Order of polyn
       is 3^n.

       Currently use one iteration so that first deriv.
       is zero at interior boundary and is exactly representable
       by low order multiwavelet without refinement */

    x = (x*x*(3.-2.*x));
    return x;
}

double mask3(const coordT& ruser) {
    coordT rsim;
    user_to_sim(ruser, rsim);
    double x= rsim[0], y=rsim[1], z=rsim[2];
    double lo = 0.0625, hi = 1.0-lo, result = 1.0;
    double rlo = 1.0/lo;

    if (x<lo)
        result *= mask1(x*rlo);
    else if (x>hi)
        result *= mask1((1.0-x)*rlo);
    if (y<lo)
        result *= mask1(y*rlo);
    else if (y>hi)
        result *= mask1((1.0-y)*rlo);
    if (z<lo)
        result *= mask1(z*rlo);
    else if (z>hi)
        result *= mask1((1.0-z)*rlo);

    return result;
}

// template<int NDIM>
// class DensityIsosurfaceCharacteristic {
//     const double rho0, twobeta;

//     double f(double rho) const {
//         double r = std::max(rho,1e-12)/rho0;
//         double r2b = std::pow(r,twobeta);
//         return r2b/(1.0+r2b);
//     }

// public:
//     DensityIsosurfaceCharacteristic(double rho0, double beta)
//         : rho0(rho0), twobeta(2.0*beta)
//     {}

//     void operator()(const Key<NDIM>& key, Tensor<double>& t) const {
//         UNARY_OPTIMIZED_ITERATOR(double, t, *_p0 = f(*_p0););
//     }
//     template <typename Archive>
//     void serialize(Archive& ar) {}
// };

// template<int NDIM>
// class DensityIsosurfaceCharacteristicDerivative {
//     const double rho0, twobeta;

//     double f(double rho) const {
//         double r = std::max(rho,1e-12)/rho0;
//         double r2b = std::pow(r,twobeta);
//         return twobeta*r2b/(rho0*r*(1.0+r2b)*(1.0+r2b));
//     }

// public:
//     DensityIsosurfaceCharacteristicDerivative(double rho0, double beta)
//         : rho0(rho0), twobeta(2.0*beta)
//     {}

//     void operator()(const Key<NDIM>& key, Tensor<double>& t) const {
//         UNARY_OPTIMIZED_ITERATOR(double, t, *_p0 = f(*_p0););
//     }
//     template <typename Archive>
//     void serialize(Archive& ar) {}
// };


/// simple projector class for 1- and 2-particle projectors
template<typename T, std::size_t NDIM>
class Projector {

    int particle_;
    std::vector<Function<T,NDIM> > p_;

public:

    Projector() : p_(std::vector<Function<T,NDIM> >()) {}

    /// simple constructor with only one orbital to project out
    Projector(const Function<T,NDIM>& p, const int particle=0)
            : particle_(particle), p_(std::vector<Function<T,NDIM> >(1,p)) {
        MADNESS_ASSERT(particle_==0 or particle_==1);
        MADNESS_ASSERT(p_.size()>0);
    }

    /// constructor with a set of orbitals to project out
    Projector(const std::vector<Function<T,NDIM> >& p, const int particle=0) : particle_(particle), p_(p) {
        MADNESS_ASSERT(particle_==0 or particle_==1);
        MADNESS_ASSERT(p_.size()>0);
    }

    int& particle() {return particle_;}
    const int& particle() const {return particle_;}

    /// get a const reference to the orbitals
    const std::vector<Function<T,NDIM> >& p() const {return p_;}

    /// project f on p: |result> =  | p><p | f>
    template<std::size_t FDIM>
    typename std::enable_if<NDIM==FDIM, Function<T,FDIM> >::type
    operator()(const Function<T,FDIM>& f) const {

        const double ovlp=inner(f,p_[0]);
        Function<T,NDIM> sum=ovlp*p_[0];

        for (unsigned int i=1; i<p_.size(); ++i) {
            const double ovlp2=inner(f,p_[i]);
            sum=(sum+ovlp2*p_[i]).truncate().reduce_rank();
        }
        return sum;
    }

//prod    /// project p out of f: |result(1,2)> = sum_p | p(1)><p(1) | f(1,2)>
//prod    template<std::size_t FDIM>
//prod    typename std::enable_if<2*NDIM==FDIM, Function<T,FDIM> >::type
//prod    operator()(const Function<T,FDIM>& f) const {
//prod        real_function_6d sum=real_factory_6d(p_.begin()->world());
//prod        for (unsigned int i=0; i<p_.size(); ++i) {
//prod            const real_function_3d pf2=f.project_out(p_[i],particle_);
//prod            real_function_6d tmp;
//prod            MADNESS_EXCEPTION("Projector class: the hartree product is inaccurate -- don't use it",1);
//prod            if (particle_==0) tmp=hartree_product(p_[i],pf2);
//prod            else tmp=hartree_product(pf2,p_[i]);
//prod            sum=(sum+tmp);
//prod        }
//prod        sum.truncate();
//prod        return sum;
//prod    }
};


class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis)
        : molecule(molecule), aobasis(aobasis) {}

    double operator()(const coordT& x) const {
        return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
    }

    std::vector<coordT> special_points() const {return molecule.get_all_coords_vec();}
};


class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
    const AtomicBasisFunction aofunc;

public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc)
        : aofunc(aofunc)
    {}

    double operator()(const coordT& x) const {
        return aofunc(x[0], x[1], x[2]);
    }

    std::vector<coordT> special_points() const {
        return std::vector<coordT>(1,aofunc.get_coords_vec());
    }
};

class MolecularDerivativeFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int atom;
    const int axis;

public:
    MolecularDerivativeFunctor(const Molecule& molecule, int atom, int axis)
        : molecule(molecule), atom(atom), axis(axis)
    {}

    double operator()(const coordT& x) const {
        return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0], x[1], x[2]);
    }

    std::vector<coordT> special_points() const {
        return std::vector<coordT>(1,molecule.get_atom(atom).get_coords());
    }
};

class CorePotentialDerivativeFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int atom;
    const int axis;
    std::vector<coordT> specialpt;
public:
    CorePotentialDerivativeFunctor(const Molecule& molecule, int atom, int axis)
        : molecule(molecule), atom(atom), axis(axis) {}

    double operator()(const coordT& r) const {
        return molecule.core_potential_derivative(atom, axis, r[0], r[1], r[2]);
    }
};

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coordT& x) const {
        return x[axis];
    }
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

/// A MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class MomentFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int i, j, k;
public:
    MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    MomentFunctor(const std::vector<int>& x) : i(x[0]), j(x[1]), k(x[2]) {}
    double operator()(const coordT& r) const {
        double xi=1.0, yj=1.0, zk=1.0;
        for (int p=0; p<i; ++p) xi *= r[0];
        for (int p=0; p<j; ++p) yj *= r[1];
        for (int p=0; p<k; ++p) zk *= r[2];
        return xi*yj*zk;
    }
};

/// A generic functor to compute external potential for TDDFT
template<typename T>
class VextCosFunctor {
  double _omega;
  Function<T,3> _f;
public:
    VextCosFunctor(World& world,
//        const std::shared_ptr<FunctionFunctorInterface<T,3> >& functor,
        const FunctionFunctorInterface<T,3>* functor,
        double omega) : _omega(omega)
    {
//      _f = factoryT(world).functor(functor);
      _f = factoryT(world).functor(functorT(new DipoleFunctor(2)));
    }
    Function<T,3> operator()(const double t) const {
        return std::cos(_omega * t) * _f;
    }
};


/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
    tensorT Q = inner(s,s);
    Q.gaxpy(0.2,s,-2.0/3.0);
    for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.0;
    return Q.scale(15.0/8.0);
}

/// Computes matrix square root (not used any more?)
tensorT sqrt(const tensorT& s, double tol=1e-8) {
    int n=s.dim(0), m=s.dim(1);
    MADNESS_ASSERT(n==m);
    tensorT c, e;
    //s.gaxpy(0.5,transpose(s),0.5); // Ensure exact symmetry
    syev(s, c, e);
    for (int i=0; i<n; ++i) {
        if (e(i) < -tol) {
            MADNESS_EXCEPTION("Matrix square root: negative eigenvalue",i);
        }
        else if (e(i) < tol) { // Ugh ..
            print("Matrix square root: Warning: small eigenvalue ", i, e(i));
            e(i) = tol;
        }
        e(i) = 1.0/sqrt(e(i));
    }
    for (int j=0; j<n; ++j) {
        for (int i=0; i<n; ++i) {
            c(j,i) *= e(i);
        }
    }
    return c;
}

template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() < 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};


struct CalculationParameters {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!                                                                   !!!
    // !!! If you add more data don't forget to add them to serialize method !!!
    // !!!                                                                   !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // First list input parameters
    double charge;              ///< Total molecular charge
    double smear;               ///< Smearing parameter
    double dconv;               ///< Density convergence
    double L;                   ///< User coordinates box size
    double maxrotn;             ///< Step restriction used in autoshift algorithm
    int nvalpha;                ///< Number of alpha virtuals to compute
    int nvbeta;                 ///< Number of beta virtuals to compute
    int nopen;                  ///< Number of unpaired electrons = napha-nbeta
    int maxiter;                ///< Maximum number of iterations
    int nio;                    ///< No. of io servers to use
    bool spin_restricted;       ///< True if spin restricted
    int plotlo,plothi;          ///< Range of MOs to print (for both spins if polarized)
    bool plotdens;              ///< If true print the density at convergence
    bool plotcoul;              ///< If true plot the total coulomb potential at convergence
    bool localize;              ///< If true solve for localized orbitals
    bool localize_pm;           ///< If true use PM for localization
    bool restart;               ///< If true restart from orbitals on disk
    bool save;                  ///< If true save orbitals to disk
    unsigned int maxsub;        ///< Size of iterative subspace ... set to 0 or 1 to disable
    int npt_plot;               ///< No. of points to use in each dim for plots
    tensorT plot_cell;          ///< lo hi in each dimension for plotting (default is all space)
    std::string aobasis;        ///< AO basis used for initial guess (6-31g or sto-3g)
    std::string core_type;      ///< core potential type ("" or "mcp")
    bool derivatives;           ///< If true calculate derivatives
    bool dipole;                ///< If true calculatio dipole moment
    bool conv_only_dens;        ///< If true remove bsh_residual from convergence criteria   how ugly name is...
    // Next list inferred parameters
    int nalpha;                 ///< Number of alpha spin electrons
    int nbeta;                  ///< Number of beta  spin electrons
    int nmo_alpha;              ///< Number of alpha spin molecular orbitals
    int nmo_beta;               ///< Number of beta  spin molecular orbitals
    double lo;                  ///< Smallest length scale we need to resolve
    std::string xc_data;         ///< XC input line
    std::vector<double> protocol_data;  ///< Calculation protocol
    bool gopt;                  ///< geometry optimizer
    double gtol;                ///< geometry tolerance
    bool gtest;                 ///< geometry tolerance
    double gval;                ///< value precision
    double gprec;               ///< gradient precision
    int  gmaxiter;               ///< optimization maxiter
    std::string algopt;         ///< algorithm used for optimization
    bool tdksprop;               ///< time-dependent Kohn-Sham equation propagate
    //bool absolvent;             ///< If true calculate solvation effects
    //double epsilon_2;           ///< dielectric constant of solvent
    //double Gamma;               ///< surface tension of solvent
    //double beta;                ///switching parameter controles boundary conditions of solvent cavity
    //double rho_0;               /// threshold density--determines size of molecular cavity
    //double sigma;               ///switching parameter controles boundary conditions of solvent cavity-(SVPE)
    bool polar;                 ///< If true calculate polarizability on
    int polar_ds;               ///< Static or Dynamic
    bool response;              ///< If true calculate response on
    int response_axis;          ///< The axis for calculation response function
    double response_freq;       ///< Frequency for calculation response function
    bool nonrotate;             ///< If true do not molcule orient
    double rconv;               ///< Response convergence
    double epsf;                ///< eps for finite field polarizability
    
    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & dconv & L & maxrotn & nvalpha & nvbeta & nopen & maxiter & nio & spin_restricted;
        ar & plotlo & plothi & plotdens & plotcoul & localize & localize_pm & restart & save & maxsub & npt_plot & plot_cell & aobasis;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
        ar & core_type & derivatives & conv_only_dens & dipole;
        ar & xc_data & protocol_data;
        ar & gopt & gtol & gtest & gval & gprec & gmaxiter & algopt & tdksprop;
        //ar & absolvent & epsilon_2 & Gamma & Gamma & beta & rho_0 & sigma;
        ar & polar & polar_ds;
        ar & response & response_axis & response_freq;
        ar & nonrotate;
        ar & rconv;
        ar & epsf;
    }

    CalculationParameters()
        : charge(0.0)
        , smear(0.0)
        , dconv(1e-4)
        , L(0.0)
        , maxrotn(0.25)
        , nvalpha(0)
        , nvbeta(0)
        , nopen(0)
        , maxiter(20)
        , nio(1)
        , spin_restricted(true)
        , plotlo(0)
        , plothi(-1)
        , plotdens(false)
        , plotcoul(false)
        , localize(true)
        , localize_pm(true)
        , restart(false)
        , save(true)
        , maxsub(8)
        , npt_plot(101)
        , aobasis("6-31g")
        , core_type("")
        , derivatives(false)
        , dipole(false)
        , conv_only_dens(false)
        , nalpha(0)
        , nbeta(0)
        , nmo_alpha(0)
        , nmo_beta(0)
        , lo(1e-10)
        , xc_data("lda")
        , protocol_data(madness::vector_factory(1e-4, 1e-6))
        , gopt(false)
        , gtol(1e-3)
        , gtest(false)
        , gval(1e-5)
        , gprec(1e-4)
        , gmaxiter(20)
        , algopt("BFGS")
        , tdksprop(false)
          //, absolvent(false)
          //, epsilon_2(78.304)
          //, Gamma(0.0719)
          //, beta(1.3)
          //, rho_0(0.00048)
          //, sigma(0.3)
        , polar(false)
        , polar_ds(0)
        , response(false)
        , response_axis(0)
        , response_freq(0.0)
        , nonrotate(false)
        , rconv(1e-6)
        , epsf(0.0)
    {}


    void read_file(const std::string& filename) {
        std::ifstream f(filename.c_str());
        position_stream(f, "dft");
        std::string s;
        xc_data = "lda";
        protocol_data = madness::vector_factory(1e-4, 1e-6);

        while (f >> s) {
            if (s == "end") {
                break;
            }
            else if (s == "charge") {
                f >> charge;
            }
            else if (s == "smear") {
                f >> smear;
            }
            else if (s == "dconv") {
                f >> dconv;
            }
            else if (s == "L") {
                f >> L;
            }
            else if (s == "maxrotn") {
                f >> maxrotn;
            }
            else if (s == "nvalpha") {
                f >> nvalpha;
            }
            else if (s == "nvbeta") {
                f >> nvbeta;
            }
            else if (s == "nopen") {
                f >> nopen;
            }
            else if (s == "unrestricted") {
                spin_restricted = false;
            }
            else if (s == "restricted") {
                spin_restricted = true;
            }
            else if (s == "maxiter") {
                f >> maxiter;
            }
            else if (s == "nio") {
                f >> nio;
            }
            else if (s == "xc") {
                char buf[1024];
                f.getline(buf,sizeof(buf));
                xc_data = buf;
            }
            else if (s == "protocol") {
                std::string buf;
                std::getline(f,buf);
                protocol_data = std::vector<double>();
                double d;
                std::stringstream s(buf);
                while (s >> d) protocol_data.push_back(d);
            }
            else if (s == "plotmos") {
                f >> plotlo >> plothi;
            }
            else if (s == "plotdens") {
                plotdens = true;
            }
            //            else if (s == "absolvent") {
            //                absolvent = true;
            //            }
            //            else if (s == "dielec") {
            //                f >> epsilon_2;
            //            }
            //            else if (s == "Gamma") {
            //                f >> Gamma;
            //            }
            //            else if (s == "beta") {
            //                f >> beta;
            //            }
            //            else if (s == "rho_0") {
            //                f >> rho_0;
            //            }
            //            else if (s == "sigma") {
            //                f >> sigma;
            //            }
            else if (s == "plotcoul") {
                plotcoul = true;
            }
            else if (s == "plotnpt") {
                f >> npt_plot;
            }
            else if (s == "plotcell") {
                plot_cell = tensorT(3L,2L);
                f >> plot_cell(0,0) >> plot_cell(0,1) >> plot_cell(1,0) >> plot_cell(1,1) >> plot_cell(2,0) >> plot_cell(2,1);
            }
            else if (s == "aobasis") {
                f >> aobasis;
                if (aobasis!="sto-3g" && aobasis!="sto-6g" && aobasis!="6-31g") {
                    std::cout << "moldft: unrecognized aobasis (sto-3g or sto-6g or 6-31g only): " << aobasis << std::endl;
                    MADNESS_EXCEPTION("input_error", 0);
                }
            }
            else if (s == "canon") {
                localize = false;
            }
            else if (s == "local") {
                localize = true;
            }
            else if (s == "pm") {
                localize_pm = true;
            }
            else if (s == "boys") {
                localize_pm = false;
            }
            else if (s == "restart") {
                restart = true;
            }
            else if (s == "save") {
                f >> save;
            }
            else if (s == "maxsub") {
                f >> maxsub;
                if (maxsub <= 0) maxsub = 1;
                if (maxsub > 20) maxsub = 20;
            }
            else if (s == "core_type") {
                f >> core_type;
            }
            else if (s == "derivatives") {
                derivatives = true;
            }
            else if (s == "dipole") {
                dipole = true;
            }
            else if (s == "convonlydens") {
                conv_only_dens = true;
            }
            else if (s == "gopt") {
               gopt = true;
            }
            else if (s == "gtol") {
                f >> gtol;
            }
            else if (s == "gtest") {
               gtest = true;
            }
            else if (s == "gval") {
                f >> gval;
            }
            else if (s == "gprec") {
                f >> gprec;
            }
            else if (s == "gmaxiter") {
                f >> gmaxiter;
            }
            else if (s == "algopt") {
                char buf[1024];
                f.getline(buf,sizeof(buf));
                algopt = buf;
            }
            else if (s == "tdksprop") {
              tdksprop = true;
            }
            else if (s == "polar") {
              polar = true;
              std::string ds;
              f >> ds;
              if(ds == "static")
                  polar_ds = 0;
              else if (ds == "dynamic")
                  polar_ds = 1;
            }
            else if (s == "response") {
              response = true;
              std::string axis;
              f >> axis;
              if(axis == "x")
                  response_axis = 0;
              else if (axis == "y") 
                  response_axis = 1;
              else if (axis == "z") 
                  response_axis = 2;
            }
            else if (s == "response_freq") {
              double freq;
              f >> freq;
                response_freq = freq;
            }
            else if (s == "nonrotate") {
              nonrotate = true; 
            }
            else if (s == "rconv") {
                f >> rconv;
            }
            else if (s == "eps_finite") {
                f >> epsf;
            }
            else {
                std::cout << "moldft: unrecognized input keyword " << s << std::endl;
                MADNESS_EXCEPTION("input error",0);
            }
            if (nopen != 0) spin_restricted = false;
        }
    }

    void set_molecular_info(const Molecule& molecule, const AtomicBasisSet& aobasis, unsigned int n_core) {
        double z = molecule.total_nuclear_charge();
        int nelec = int(z - charge - n_core*2);
        if (fabs(nelec+charge+n_core*2-z) > 1e-6) {
            error("non-integer number of electrons?", nelec+charge+n_core*2-z);
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
        nmo_alpha = std::min(nbf,nmo_alpha);
        nmo_beta = std::min(nbf,nmo_beta);
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
        //time_t t = time((time_t *) 0);
        //char *tmp = ctime(&t);
        //tmp[strlen(tmp)-1] = 0; // lose the trailing newline

        //madness::print(" date of calculation ", tmp);
        madness::print("             restart ", restart);
        madness::print(" number of processes ", world.size());
        madness::print("   no. of io servers ", nio);
        madness::print("     simulation cube ", -L, L);
        madness::print("        total charge ", charge);
        madness::print("            smearing ", smear);
        madness::print(" number of electrons ", nalpha, nbeta);
        madness::print("  number of orbitals ", nmo_alpha, nmo_beta);
        madness::print("     spin restricted ", spin_restricted);
        madness::print("       xc functional ", xc_data);
#ifdef MADNESS_HAS_LIBXC
        madness::print("         xc libraray ", "libxc");
#else
        madness::print("         xc libraray ", "default (lda only)");
#endif
        if (core_type != "")
            madness::print("           core type ", core_type);
        madness::print(" initial guess basis ", aobasis);
        madness::print(" max krylov subspace ", maxsub);
        madness::print("    compute protocol ", protocol_data);
        madness::print(" density convergence ", dconv);
        madness::print(" response convergence ", rconv);
        madness::print("    maximum rotation ", maxrotn);
        if (conv_only_dens)
            madness::print(" Convergence criterion is only density delta.");
        else
            madness::print(" Convergence criteria are density delta & BSH residual.");
        madness::print("        plot density ", plotdens);
        madness::print("        plot coulomb ", plotcoul);
        madness::print("        plot orbital ", plotlo, plothi);
        madness::print("        plot npoints ", npt_plot);
        if (plot_cell.size() > 0)
            madness::print("        plot  volume ", plot_cell(0,0), plot_cell(0,1),
                           plot_cell(1,0), plot_cell(1,1), plot_cell(2,0), plot_cell(2,1));
        else
            madness::print("        plot  volume ", "default");

        std::string loctype = "pm";
        if (!localize_pm) loctype = "boys";
        if (localize)
            madness::print("  localized orbitals ", loctype);
        else
            madness::print("  canonical orbitals ");
        if (derivatives)
            madness::print("    calc derivatives ");
        if (dipole)
            madness::print("         calc dipole ");
        //        if(absolvent){
        //            madness::print("       isodensity solvation ", absolvent);
        //            madness::print("       surface tension      ", Gamma);
        //            madness::print("       switching param(beta)", beta);
        //            madness::print("       dielectric constant  ", epsilon_2);
        //            madness::print("       threshold density    ", rho_0);
        //        }
    }
//};
    void gprint(World& world) const {
        madness::print(" Optimizer parameters:           ");
        madness::print("   Maximum iterations (gmaxiter) ", gmaxiter);
        madness::print("                Tolerance (gtol) ", gtol);
        madness::print("           Gradient value (gval) ", gval);
        madness::print("      Gradient precision (gprec) ", gprec);
        madness::print(" Optimization algorithm (algopt) ", algopt);
        madness::print(" Gradient numerical test (gtest) ", gtest);
    }
};

struct Calculation {
    std::shared_ptr<PotentialManager> potentialmanager;
    Molecule molecule;
    CalculationParameters param;
    XCfunctional xc;
    AtomicBasisSet aobasis;
    functionT vacuo_rho;
    functionT rhoT;
    functionT rho_elec;
    functionT rhon;
    functionT mol_mask;
    functionT Uabinit;
    functionT mask;
    vecfuncT amo, bmo;
    vecfuncT amo_p, bmo_p;
    vecfuncT amo_m, bmo_m;
    std::vector<int> aset, bset;
    std::vector<int> aset_mp, bset_mp;
    vecfuncT ao;
    std::vector<int> at_to_bf, at_nbf;
    tensorT aocc, bocc;
    tensorT aocc_mp, bocc_mp;
    tensorT aeps, beps;
    tensorT aeps_mp, beps_mp;
    poperatorT coulop;
    std::vector< std::shared_ptr<real_derivative_3d> > gradop;
    double vtol;
    double current_energy;
    double esol;//etot;
    double vacuo_energy;
    static const int vnucextra = 12; // load balance parameter for nuclear pot.

    Calculation(World & world, const char *filename)
    {
        if(world.rank() == 0) {
            molecule.read_file(filename);
            param.read_file(filename);
            unsigned int n_core = 0;
            if (param.core_type != "") {
                molecule.read_core_file(param.core_type);
                param.aobasis = molecule.guess_file();
                n_core = molecule.n_core_orb_all();
            }
            // if nonrotate : not rotate the axis
            if (param.nonrotate){
               molecule.center();
            }
            else {
                 molecule.orient();
            } 
            aobasis.read_file(param.aobasis);
            param.set_molecular_info(molecule, aobasis, n_core);
        }
        world.gop.broadcast_serializable(molecule, 0);
        world.gop.broadcast_serializable(param, 0);
        world.gop.broadcast_serializable(aobasis, 0);

        xc.initialize(param.xc_data, !param.spin_restricted, world);
        //xc.plot();

        FunctionDefaults<3>::set_cubic_cell(-param.L, param.L);
        set_protocol(world, 1e-4);

        potentialmanager = std::shared_ptr<PotentialManager>(new PotentialManager(molecule, param.core_type));
    }

    void set_protocol(World & world, double thresh)
    {
        int k;
        // Allow for imprecise conversion of threshold
        if(thresh >= 0.9e-2)
            k = 4;
        else if(thresh >= 0.9e-4)
            k = 6;
        else if(thresh >= 0.9e-6)
            k = 8;
        else if(thresh >= 0.9e-8)
            k = 10;
        else
            k = 12;

        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_initial_level(2);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(false);
        GaussianConvolution1DCache<double>::map.clear();
        double safety = 0.1;
        vtol = FunctionDefaults<3>::get_thresh() * safety;
        coulop = poperatorT(CoulombOperatorPtr(world, param.lo, thresh));
        gradop = gradient_operator<double,3>(world);
        mask = functionT(factoryT(world).f(mask3).initial_level(4).norefine());
        if(world.rank() == 0){
            print("\nSolving with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "   dconv", std::max(thresh, param.dconv), "\n");
        }
    }

    void save_mos(World& world) {
        archive::ParallelOutputArchive ar(world, "restartdata", param.nio);
        ar & param.spin_restricted;
        ar & (unsigned int)(amo.size());
        ar & aeps & aocc & aset;
        for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
        if (!param.spin_restricted) {
            ar & (unsigned int)(bmo.size());
            ar & beps & bocc & bset;
            for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];
        }
    }

    void load_mos(World& world) {
        const double trantol = vtol / std::min(30.0, double(param.nalpha));
        const double thresh = FunctionDefaults<3>::get_thresh();
        const int k = FunctionDefaults<3>::get_k();
        unsigned int nmo = 0;
        bool spinrest = false;
        amo.clear(); bmo.clear();

        archive::ParallelInputArchive ar(world, "restartdata");

        /*
          File format:

          bool spinrestricted --> if true only alpha orbitals are present

          unsigned int nmo_alpha;
          Tensor<double> aeps;
          Tensor<double> aocc;
          vector<int> aset;
          for i from 0 to nalpha-1:
          .   Function<double,3> amo[i]

          repeat for beta if !spinrestricted

        */

        // LOTS OF LOGIC MISSING HERE TO CHANGE OCCUPATION NO., SET,
        // EPS, SWAP, ... sigh

        ar & spinrest;

        ar & nmo;
        MADNESS_ASSERT(nmo >= unsigned(param.nmo_alpha));
        ar & aeps & aocc & aset;
        amo.resize(nmo);
        for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
        unsigned int n_core = molecule.n_core_orb_all();
        if (nmo > unsigned(param.nmo_alpha)) {
            aset = vector<int>(aset.begin()+n_core, aset.begin()+n_core+param.nmo_alpha);
            amo = vecfuncT(amo.begin()+n_core, amo.begin()+n_core+param.nmo_alpha);
            aeps = copy(aeps(Slice(n_core, n_core+param.nmo_alpha-1)));
            aocc = copy(aocc(Slice(n_core, n_core+param.nmo_alpha-1)));
        }

        if (amo[0].k() != k) {
            reconstruct(world,amo);
            for(unsigned int i = 0;i < amo.size();++i) amo[i] = madness::project(amo[i], k, thresh, false);
            world.gop.fence();
        }
        normalize(world, amo);
        amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)), trantol, true);
        truncate(world, amo);
        normalize(world, amo);

        if (!param.spin_restricted) {

            if (spinrest) { // Only alpha spin orbitals were on disk
                MADNESS_ASSERT(param.nmo_alpha >= param.nmo_beta);
                bmo.resize(param.nmo_beta);
                bset.resize(param.nmo_beta);
                beps = copy(aeps(Slice(0,param.nmo_beta-1)));
                bocc = copy(aocc(Slice(0,param.nmo_beta-1)));
                for (int i=0; i<param.nmo_beta; ++i) bmo[i] = copy(amo[i]);
            }
            else {
                ar & nmo;
                ar & beps & bocc & bset;

                bmo.resize(nmo);
                for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];

                if (nmo > unsigned(param.nmo_beta)) {
                    bset = vector<int>(bset.begin()+n_core, bset.begin()+n_core+param.nmo_beta);
                    bmo = vecfuncT(bmo.begin()+n_core, bmo.begin()+n_core+param.nmo_beta);
                    beps = copy(beps(Slice(n_core, n_core+param.nmo_beta-1)));
                    bocc = copy(bocc(Slice(n_core, n_core+param.nmo_beta-1)));
                }

                if (bmo[0].k() != k) {
                    reconstruct(world,bmo);
                    for(unsigned int i = 0;i < bmo.size();++i) bmo[i] = madness::project(bmo[i], k, thresh, false);
                    world.gop.fence();
                }

                normalize(world, bmo);
                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo, bmo)), trantol, true);
                truncate(world, bmo);
                normalize(world, bmo);

            }
        }
    }

    void do_plots(World& world) {
        START_TIMER(world);

        std::vector<long> npt(3,param.npt_plot);

        if (param.plot_cell.size() == 0)
            param.plot_cell = copy(FunctionDefaults<3>::get_cell());

        if (param.plotdens || param.plotcoul) {
            functionT rho;
            rho = make_density(world, aocc, amo);

            if (param.spin_restricted) {
                rho.scale(2.0);
            }
            else {
                functionT rhob = make_density(world, bocc, bmo);
                functionT rho_spin = rho - rhob;
                rho += rhob;
                plotdx(rho_spin, "spin_density.dx", param.plot_cell, npt, true);

            }
            plotdx(rho, "total_density.dx", param.plot_cell, npt, true);
            if (param.plotcoul) {
                real_function_3d vnuc = potentialmanager->vnuclear();
                functionT vlocl = vnuc + apply(*coulop, rho);
                vlocl.truncate();
                vlocl.reconstruct();
                plotdx(vlocl, "coulomb.dx", param.plot_cell, npt, true);
            }
        }

        for (int i=param.plotlo; i<=param.plothi; ++i) {
            char fname[256];
            if (i < param.nalpha) {
                sprintf(fname, "amo-%5.5d.dx", i);
                plotdx(amo[i], fname, param.plot_cell, npt, true);
            }
            if (!param.spin_restricted && i < param.nbeta) {
                sprintf(fname, "bmo-%5.5d.dx", i);
                plotdx(bmo[i], fname, param.plot_cell, npt, true);
            }
        }
        END_TIMER(world, "plotting");
    }

    void project(World & world)
    {
        reconstruct(world, amo);
        for(unsigned int i = 0;i < amo.size();++i){
            amo[i] = madness::project(amo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
        }
        world.gop.fence();
        truncate(world, amo);
        normalize(world, amo);
        if(param.nbeta && !param.spin_restricted){
            reconstruct(world, bmo);
            for(unsigned int i = 0;i < bmo.size();++i){
                bmo[i] = madness::project(bmo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
            }
            world.gop.fence();
            truncate(world, bmo);
            normalize(world, bmo);
        }
    }

    void make_nuclear_potential(World & world)
    {
        START_TIMER(world);
        potentialmanager->make_nuclear_potential(world);
        END_TIMER(world, "Project vnuclear");
    }

    void project_ao_basis(World & world)
    {
        // Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
        aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

        START_TIMER(world);
        ao = vecfuncT(aobasis.nbf(molecule));
        for(int i = 0;i < aobasis.nbf(molecule);++i){
            functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule, i)));
            ao[i] = factoryT(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(1);
        }
        world.gop.fence();
        truncate(world, ao);
        normalize(world, ao);
        END_TIMER(world, "project ao basis");
	print_meminfo(world.rank(), "project ao basis");
    }

    double PM_q(const tensorT & S, const double * restrict Ci, const double * restrict Cj, int lo, int nbf)
    {
        double qij = 0.0;
        if (nbf == 1) { // H atom in STO-3G ... often lots of these!
            qij = Ci[lo]*S(0,0)*Cj[lo];
        }
        else {
            for(int mu = 0;mu < nbf;++mu){
                double Smuj = 0.0;
                for(int nu = 0;nu < nbf;++nu){
                    Smuj += S(mu, nu) * Cj[nu + lo];
                }
                qij += Ci[mu + lo] * Smuj;
            }
        }

        return qij;
    }


    void localize_PM_ij(const int seti, const int setj, 
                        const double tol, const double thetamax,
                        const int natom, const int nao,  const int nmo,
                        const std::vector<tensorT>& Svec, 
                        const std::vector<int>& at_to_bf, const std::vector<int>& at_nbf, 
                        long& ndone_iter, double& maxtheta, 
                        double * restrict Qi, double * restrict Qj,  
                        double * restrict Ci, double * restrict Cj, 
                        double * restrict Ui, double * restrict Uj)
    {
        if(seti == setj){
            // Q could be recomputed for each ij, but by saving it we reduce computation by a factor
            // of 2 when far from convergence, and as we approach convergence we save more since
            // most work is associated with computing ovij.

            double ovij = 0.0;
            for(long a = 0;a < natom;++a)
                ovij += Qi[a] * Qj[a];
            
            print("ovij", ovij);
            if(fabs(ovij) > tol * tol){
                double aij = 0.0;
                double bij = 0.0;
                for(long a = 0;a < natom;++a){
                    double qiia = Qi[a];
                    double qija = PM_q(Svec[a], Ci, Cj, at_to_bf[a], at_nbf[a]);
                    double qjja = Qj[a];
                    double d = qiia - qjja;
                    aij += qija * qija - 0.25 * d * d;
                    bij += qija * d;
                }
                double theta = 0.25 * acos(-aij / sqrt(aij * aij + bij * bij));

                if(bij > 0.0)
                    theta = -theta;
                
                print("theta", theta);
                if(theta > thetamax)
                    theta = thetamax;
                else
                    if(theta < -thetamax)
                        theta = -thetamax;
                
                maxtheta = std::max(fabs(theta), maxtheta);
                if(fabs(theta) >= tol){
                    ++ndone_iter;
                    double c = cos(theta);
                    double s = sin(theta);
                    print(c,s);
                    drot(nao, Ci, Cj, s, c, 1);
                    drot(nmo, Ui, Uj, s, c, 1);
                    for(long a = 0;a < natom;++a){
                        Qi[a] = PM_q(Svec[a], Ci, Ci, at_to_bf[a], at_nbf[a]);
                        Qj[a] = PM_q(Svec[a], Cj, Cj, at_to_bf[a], at_nbf[a]);
                    }
                }
            }
        }
    }
    


    void localize_PM_task_kernel(tensorT & Q, std::vector<tensorT> & Svec, tensorT & C,
                                 const bool & doprint, const std::vector<int> & set,
                                 const double thetamax, tensorT & U, const double thresh)
    {
        long nmo = C.dim(0);
        long nao = C.dim(1);
        long natom = molecule.natom();

        for(long i = 0;i < nmo;++i){
            for(long a = 0;a < natom;++a){
                Q(i, a) = PM_q(Svec[a], &C(i,0), &C(i,0), at_to_bf[a], at_nbf[a]);
            }
        }

        print("Q\n", Q);

        double tol = 0.1;
        long ndone = 0;
        for(long iter = 0;iter < 100;++iter){

            // Diagnostics at beginning of iteration
            double sum = 0.0;
            for(long i = 0;i < nmo;++i){
                for(long a = 0;a < natom;++a){
                    double qiia = Q(i, a);
                    sum += qiia * qiia;
                }
            }
            if(doprint)
                printf("iteration %ld sum=%.4f ndone=%ld tol=%.2e\n", iter, sum, ndone, tol);
            // End diagnostics at beginning of iteration

            // Initialize variables for convergence test
            long ndone_iter = 0;
            double maxtheta = 0.0;
            for(long i = 0;i < nmo;++i){
                for(long j = 0;j < i;++j){

                    localize_PM_ij(set[i], set[j], 
                                   tol, thetamax, 
                                   natom, nao, nmo, 
                                   Svec, 
                                   at_to_bf, at_nbf, 
                                   ndone_iter, maxtheta, 
                                   &Q(i,0), &Q(j,0),
                                   &C(i,0), &C(j,0),
                                   &U(i,0), &U(j,0));

                }
            }

            ndone += ndone_iter;
            if(ndone_iter == 0 && tol == thresh){
                if(doprint)
                    print("PM localization converged in", ndone,"steps");

                break;
            }
            tol = std::max(0.1 * std::min(maxtheta, tol), thresh);
        }
    }

    tensorT localize_PM(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true, const bool doprint = false)
    {
        START_TIMER(world);
        tensorT UT = distributed_localize_PM(world, mo, ao, set, at_to_bf, at_nbf, thresh, thetamax, randomize, doprint);
        END_TIMER(world, "Pipek-Mezy distributed ");
        //print(UT);

        return UT;


        START_TIMER(world);
        long nmo = mo.size();
        long natom = molecule.natom();

        tensorT S = matrix_inner(world, ao, ao, true);
        std::vector<tensorT> Svec(natom);
        for(long a = 0;a < natom;++a){
            Slice as(at_to_bf[a], at_to_bf[a] + at_nbf[a] - 1);
            Svec[a] = copy(S(as, as));
        }
        S = tensorT();
        tensorT C = matrix_inner(world, mo, ao);
        tensorT U(nmo, nmo);
        tensorT Q(nmo, natom);
        if(world.rank() == 0){
            for(long i = 0;i < nmo;++i)
                U(i, i) = 1.0;

            localize_PM_task_kernel(Q, Svec, C, doprint, set, thetamax, U, thresh);
            U = transpose(U);

            // Fix orbital orders
	    bool switched = true;
	    while (switched) {
	      switched = false;
	      for (int i=0; i<nmo; i++) {
		for (int j=i+1; j<nmo; j++) {
		  if (set[i] == set[j]) {
		    double sold = U(i,i)*U(i,i) + U(j,j)*U(j,j);
		    double snew = U(i,j)*U(i,j) + U(j,i)*U(j,i);
		    if (snew > sold) {
		      tensorT tmp = copy(U(_,i));
		      U(_,i) = U(_,j);
		      U(_,j) = tmp;
		      switched = true;
		    }
		  }
		}
	      }
	    }

            // Fix phases.
            for (long i=0; i<nmo; ++i) {
                if (U(i,i) < 0.0) U(_,i).scale(-1.0);
            }
        }
        world.gop.broadcast(U.ptr(), U.size(), 0);
        END_TIMER(world, "Pipek-Mezy localize");
	print_meminfo(world.rank(), "Pipek-Mezy localize");
        return U;
    }

    void analyze_vectors(World & world, const vecfuncT & mo, const tensorT & occ = tensorT(), const tensorT & energy = tensorT(), const std::vector<int> & set = std::vector<int>())
    {
        tensorT Saomo = matrix_inner(world, ao, mo);
        tensorT Saoao = matrix_inner(world, ao, ao, true);
        int nmo = mo.size();
        tensorT rsq, dip(3, nmo);
        {
            functionT frsq = factoryT(world).f(rsquared).initial_level(4);
            rsq = inner(world, mo, mul_sparse(world, frsq, mo, vtol));
            for(int axis = 0;axis < 3;++axis){
                functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
                dip(axis, _) = inner(world, mo, mul_sparse(world, fdip, mo, vtol));
                for(int i = 0;i < nmo;++i)
                    rsq(i) -= dip(axis, i) * dip(axis, i);

            }
        }
            tensorT C;

            START_TIMER(world);
            gesvp(world, Saoao, Saomo, C);
            END_TIMER(world, " compute eigen gesv analize vectors");
        if(world.rank() == 0){
            C = transpose(C);
            long nmo = mo.size();
            for(long i = 0;i < nmo;++i){
                printf("  MO%4ld : ", i);
                if(set.size())
                    printf("set=%d : ", set[i]);

                if(occ.size())
                    printf("occ=%.2f : ", occ(i));

                if(energy.size())
                    printf("energy=%13.8f : ", energy(i));

                printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i), dip(1, i), dip(2, i), sqrt(rsq(i)));
                aobasis.print_anal(molecule, C(i, _));
            }
        }

    }

    inline double DIP(const tensorT & dip, int i, int j, int k, int l)
    {
        return dip(i, j, 0) * dip(k, l, 0) + dip(i, j, 1) * dip(k, l, 1) + dip(i, j, 2) * dip(k, l, 2);
    }

    tensorT localize_boys(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true)
    {
        START_TIMER(world);
        const bool doprint = false;
        long nmo = mo.size();
        tensorT dip(nmo, nmo, 3);
        for(int axis = 0;axis < 3;++axis){
            functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
            dip(_, _, axis) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
        }
        tensorT U(nmo, nmo);
        if(world.rank() == 0){
            for(long i = 0;i < nmo;++i)
                U(i, i) = 1.0;

            double tol = thetamax;
            long ndone = 0;
            bool converged = false;
            for(long iter = 0;iter < 300;++iter){
                double sum = 0.0;
                for(long i = 0;i < nmo;++i){
                    sum += DIP(dip, i, i, i, i);
                }
                long ndone_iter = 0;
                double maxtheta = 0.0;
                if(doprint)
                    printf("iteration %ld sum=%.4f ndone=%ld tol=%.2e\n", iter, sum, ndone, tol);

                for(long i = 0;i < nmo;++i){
                    for(long j = 0;j < i;++j){
                        if (set[i] == set[j]) {
                            double g = DIP(dip, i, j, j, j) - DIP(dip, i, j, i, i);
                            double h = 4.0 * DIP(dip, i, j, i, j) + 2.0 * DIP(dip, i, i, j, j) - DIP(dip, i, i, i, i) - DIP(dip, j, j, j, j);
                            double sij = DIP(dip, i, j, i, j);
                            bool doit = false;
                            if(h >= 0.0){
                                doit = true;
                                if(doprint)
                                    print("             forcing negative h", i, j, h);

                                h = -1.0;
                            }
                            double theta = -g / h;
                            maxtheta = std::max<double>(std::abs(theta), maxtheta);
                            if(fabs(theta) > thetamax){
                                doit = true;
                                if(doprint)
                                    print("             restricting", i, j);

                                if(g < 0)
                                    theta = -thetamax;

                                else
                                    theta = thetamax * 0.8;

                            }
                            bool randomized = false;
                            if(randomize && iter == 0 && sij > 0.01 && fabs(theta) < 0.01){
                                randomized = true;
                                if(doprint)
                                    print("             randomizing", i, j);

                                theta += (RandomValue<double>() - 0.5);
                            }
                            if(fabs(theta) >= tol || randomized || doit){
                                ++ndone_iter;
                                if(doprint)
                                    print("     rotating", i, j, theta);

                                double c = cos(theta);
                                double s = sin(theta);
                                drot3(nmo, &dip(i, 0, 0), &dip(j, 0, 0), s, c, 1);
                                drot3(nmo, &dip(0, i, 0), &dip(0, j, 0), s, c, nmo);
                                drot(nmo, &U(i, 0), &U(j, 0), s, c, 1);
                            }
                        }
                    }
                }

                ndone += ndone_iter;
                if(ndone_iter == 0 && tol == thresh){
                    if(doprint)
                        print("Boys localization converged in", ndone,"steps");

                    converged = true;
                    break;
                }
                tol = std::max(0.1 * maxtheta, thresh);
            }

            if(!converged){
                print("warning: boys localization did not fully converge: ", ndone);
            }
            U = transpose(U);

	    bool switched = true;
	    while (switched) {
	      switched = false;
	      for (int i=0; i<nmo; i++) {
		for (int j=i+1; j<nmo; j++) {
		  if (set[i] == set[j]) {
		    double sold = U(i,i)*U(i,i) + U(j,j)*U(j,j);
		    double snew = U(i,j)*U(i,j) + U(j,i)*U(j,i);
		    if (snew > sold) {
		      tensorT tmp = copy(U(_,i));
		      U(_,i) = U(_,j);
		      U(_,j) = tmp;
		      switched = true;
		    }
		  }
		}
	      }
	    }

        // Fix phases.
        for (long i=0; i<nmo; ++i) {
            if (U(i,i) < 0.0) U(_,i).scale(-1.0);
        }

        }

        world.gop.broadcast(U.ptr(), U.size(), 0);
        END_TIMER(world, "Boys localize");
        return U;
    }

    tensorT kinetic_energy_matrix(World & world, const vecfuncT & v)
    {
        reconstruct(world, v);
        int n = v.size();
        tensorT r(n, n);
        for(int axis = 0;axis < 3;++axis){
            vecfuncT dv = apply(world, *(gradop[axis]), v);
            r += matrix_inner(world, dv, dv, true);
            dv.clear();
        }
        return r.scale(0.5);
    }


    vecfuncT core_projection(World & world, const vecfuncT & psi, const bool include_Bc = true)
    {
        int npsi = psi.size();
        if (npsi == 0) return psi;
        int natom = molecule.natom();
        vecfuncT proj = zero_functions<double,3>(world, npsi);
        tensorT overlap_sum(static_cast<long>(npsi));

        for (int i=0; i<natom; ++i) {
            Atom at = molecule.get_atom(i);
            unsigned int atn = at.atomic_number;
            unsigned int nshell = molecule.n_core_orb(atn);
            if (nshell == 0) continue;
            for (unsigned int c=0; c<nshell; ++c) {
                unsigned int l = molecule.get_core_l(atn, c);
                int max_m = (l+1)*(l+2)/2;
                nshell -= max_m - 1;
                for (int m=0; m<max_m; ++m) {
                    functionT core = factoryT(world).functor(functorT(new CoreOrbitalFunctor(molecule, i, c, m)));
                    tensorT overlap = inner(world, core, psi);
                    overlap_sum += overlap;
                    for (int j=0; j<npsi; ++j) {
                        if (include_Bc) overlap[j] *= molecule.get_core_bc(atn, c);
                        proj[j] += core.scale(overlap[j]);
                    }
                }
            }
            world.gop.fence();
        }
        if (world.rank() == 0) print("sum_k <core_k|psi_i>:", overlap_sum);
        return proj;
    }

    double core_projector_derivative(World & world, const vecfuncT & mo, const tensorT & occ, int atom, int axis)
    {
        vecfuncT cores, dcores;
        std::vector<double> bc;
        unsigned int atn = molecule.get_atom(atom).atomic_number;
        unsigned int ncore = molecule.n_core_orb(atn);

        // projecting core & d/dx core
        for (unsigned int c=0; c<ncore; ++c) {
            unsigned int l = molecule.get_core_l(atn, c);
            int max_m = (l+1)*(l+2)/2;
            for (int m=0; m<max_m; ++m) {
                functorT func = functorT(new CoreOrbitalFunctor(molecule, atom, c, m));
                cores.push_back(functionT(factoryT(world).functor(func).truncate_on_project()));
                func = functorT(new CoreOrbitalDerivativeFunctor(molecule, atom, axis, c, m));
                dcores.push_back(functionT(factoryT(world).functor(func).truncate_on_project()));
                bc.push_back(molecule.get_core_bc(atn, c));
            }
        }

        // calc \sum_i occ_i <psi_i|(\sum_c Bc d/dx |core><core|)|psi_i>
        double r = 0.0;
        for (unsigned int c=0; c<cores.size(); ++c) {
            double rcore= 0.0;
            tensorT rcores = inner(world, cores[c], mo);
            tensorT rdcores = inner(world, dcores[c], mo);
            for (unsigned int i=0; i<mo.size(); ++i) {
                rcore += rdcores[i] * rcores[i] * occ[i];
            }
            r += 2.0 * bc[c] * rcore;
        }

        return r;
    }

    void initial_guess(World & world)
    {
        START_TIMER(world);
        if (param.restart) {
            load_mos(world);
        }
        else {
            // Use the initial density and potential to generate a better process map
            functionT rho = factoryT(world).functor(functorT(new MolecularGuessDensityFunctor(molecule, aobasis))).truncate_on_project();
            END_TIMER(world, "guess density");
            double nel = rho.trace();
            if(world.rank() == 0)
                print("guess dens trace", nel);

            if(world.size() > 1) {
                START_TIMER(world);
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0), false);
                lb.add_tree(rho, lbcost<double,3>(1.0, 8.0), true);

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
                END_TIMER(world, "guess loadbal");
            }

            // Diag approximate fock matrix to get initial mos
            functionT vlocal;
            if(param.nalpha + param.nbeta > 1){
                START_TIMER(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                vlocal = vnuc + apply(*coulop, rho);
                END_TIMER(world, "guess Coulomb potn");
                bool save = param.spin_restricted;
                param.spin_restricted = true;
                vlocal = vlocal + make_lda_potential(world, rho);
                vlocal.truncate();
                param.spin_restricted = save;
            } else {
                real_function_3d vnuc = potentialmanager->vnuclear();
                vlocal = vnuc;
            }
            rho.clear();
            vlocal.reconstruct();
            if(world.size() > 1){
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0), false);
                for(unsigned int i = 0;i < ao.size();++i){
                    lb.add_tree(ao[i], lbcost<double,3>(1.0, 8.0), false);
                }

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
            }

            tensorT overlap = matrix_inner(world, ao, ao, true);
            START_TIMER(world);
            tensorT kinetic = kinetic_energy_matrix(world, ao);
            END_TIMER(world, "guess Kinet potn");
            reconstruct(world, ao);
            vlocal.reconstruct();
            vecfuncT vpsi = mul_sparse(world, vlocal, ao, vtol);
            compress(world, vpsi);
            truncate(world, vpsi);
            compress(world, ao);
            tensorT potential = matrix_inner(world, vpsi, ao, true);
            vpsi.clear();
            tensorT fock = kinetic + potential;
            fock = 0.5 * (fock + transpose(fock));
            tensorT c, e;

            START_TIMER(world);
            sygvp(world, fock, overlap, 1, c, e);
            END_TIMER(world, "guess eigen sol");
	    print_meminfo(world.rank(), "guess eigen sol");

	    // NAR 7/5/2013
            // commented out because it generated a lot of output
            // if(world.rank() == 0 && 0){
            //   print("initial eigenvalues");
            //   print(e);
            //   print("\n\nWSTHORNTON: initial eigenvectors");
            //   print(c);
            // }

            compress(world, ao);

            unsigned int ncore = 0;
            if (param.core_type != "") {
                ncore = molecule.n_core_orb_all();
            }
            amo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha - 1)), 0.0, true);
            truncate(world, amo);
            normalize(world, amo);
            aeps = e(Slice(ncore, ncore + param.nmo_alpha - 1));

            aocc = tensorT(param.nmo_alpha);
            for(int i = 0;i < param.nalpha;++i)
                aocc[i] = 1.0;

            aset = std::vector<int>(param.nmo_alpha,0);
            if(world.rank() == 0)
                std::cout << "alpha set " << 0 << " " << 0 << "-";

            for(int i = 1;i < param.nmo_alpha;++i) {
                aset[i] = aset[i - 1];
                //vamastd::cout << "aeps -" << i << "- " << aeps[i] << std::endl;
                if(aeps[i] - aeps[i - 1] > 1.5 || aocc[i] != 1.0){
                    ++(aset[i]);
                    if(world.rank() == 0){
                        std::cout << i - 1 << std::endl;
                        std::cout << "alpha set " << aset[i] << " " << i << "-";
                    }
                }
            }
            if(world.rank() == 0)
                std::cout << param.nmo_alpha - 1 << std::endl;

            if(param.nbeta && !param.spin_restricted){
                bmo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta - 1)), 0.0, true);
                truncate(world, bmo);
                normalize(world, bmo);
                beps = e(Slice(ncore, ncore + param.nmo_beta - 1));
                bocc = tensorT(param.nmo_beta);
                for(int i = 0;i < param.nbeta;++i)
                    bocc[i] = 1.0;

                bset = std::vector<int>(param.nmo_beta,0);
                if(world.rank() == 0)
                    std::cout << " beta set " << 0 << " " << 0 << "-";

                for(int i = 1;i < param.nmo_beta;++i) {
                    bset[i] = bset[i - 1];
                    if(beps[i] - beps[i - 1] > 1.5 || bocc[i] != 1.0){
                        ++(bset[i]);
                        if(world.rank() == 0){
                            std::cout << i - 1 << std::endl;
                            std::cout << " beta set " << bset[i] << " " << i << "-";
                        }
                    }
                }
                if(world.rank() == 0)
                    std::cout << param.nmo_beta - 1 << std::endl;
            }
        }
    }

    void initial_guess_mp(World & world)
    {
        START_TIMER(world);
        if (param.restart) {
            load_mos(world);
        }
        else {
            // Use the initial density and potential to generate a better process map
            functionT rho = factoryT(world).functor(functorT(new MolecularGuessDensityFunctor(molecule, aobasis))).truncate_on_project();
            END_TIMER(world, "guess density");
            double nel = rho.trace();
            if(world.rank() == 0)
                print("guess dens trace", nel);

            if(world.size() > 1) {
                START_TIMER(world);
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0), false);
                lb.add_tree(rho, lbcost<double,3>(1.0, 8.0), true);

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
                END_TIMER(world, "guess loadbal");
            }

            // Diag approximate fock matrix to get initial mos
            functionT vlocal;
            if(param.nalpha + param.nbeta > 1){
                START_TIMER(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                vlocal = vnuc + apply(*coulop, rho);
                END_TIMER(world, "guess Coulomb potn");
                bool save = param.spin_restricted;
                param.spin_restricted = true;
                vlocal = vlocal + make_lda_potential(world, rho);
                vlocal.truncate();
                param.spin_restricted = save;
            } else {
                real_function_3d vnuc = potentialmanager->vnuclear();
                vlocal = vnuc;
            }
            rho.clear();
            vlocal.reconstruct();
            if(world.size() > 1){
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc = potentialmanager->vnuclear();
                lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0), false);
                for(unsigned int i = 0;i < ao.size();++i){
                    lb.add_tree(ao[i], lbcost<double,3>(1.0, 8.0), false);
                }

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
            }

            tensorT overlap = matrix_inner(world, ao, ao, true);
            START_TIMER(world);
            tensorT kinetic = kinetic_energy_matrix(world, ao);
            END_TIMER(world, "guess Kinet potn");
            reconstruct(world, ao);
            vlocal.reconstruct();
            vecfuncT vpsi = mul_sparse(world, vlocal, ao, vtol);
            compress(world, vpsi);
            truncate(world, vpsi);
            compress(world, ao);
            tensorT potential = matrix_inner(world, vpsi, ao, true);
            vpsi.clear();
            tensorT fock = kinetic + potential;
            fock = 0.5 * (fock + transpose(fock));
            tensorT c, e;

            START_TIMER(world);
            sygvp(world, fock, overlap, 1, c, e);
            END_TIMER(world, "guess eigen sol");
	    print_meminfo(world.rank(), "guess eigen sol");

	    // NAR 7/5/2013
            // commented out because it generated a lot of output
            // if(world.rank() == 0 && 0){
            //   print("initial eigenvalues");
            //   print(e);
            //   print("\n\nWSTHORNTON: initial eigenvectors");
            //   print(c);
            // }

            compress(world, ao);

            unsigned int ncore = 0;
            if (param.core_type != "") {
                ncore = molecule.n_core_orb_all();
            }
            amo_p = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha - 1)), 0.0, true);
            amo_m = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha - 1)), 0.0, true);
            truncate(world, amo_p);
            truncate(world, amo_m);
            normalize(world, amo_p);
            normalize(world, amo_m);
            aeps_mp = e(Slice(ncore, ncore + param.nmo_alpha - 1));

            aocc_mp = tensorT(param.nmo_alpha);
            for(int i = 0;i < param.nalpha;++i)
                aocc_mp[i] = 1.0;

            aset_mp = std::vector<int>(param.nmo_alpha,0);
            if(world.rank() == 0)
                std::cout << "alpha set mp " << 0 << " " << 0 << "-";

            for(int i = 1;i < param.nmo_alpha;++i) {
                aset_mp[i] = aset_mp[i - 1];
                //vamastd::cout << "aeps -" << i << "- " << aeps[i] << std::endl;
                if(aeps_mp[i] - aeps_mp[i - 1] > 1.5 || aocc_mp[i] != 1.0){
                    ++(aset_mp[i]);
                    if(world.rank() == 0){
                        std::cout << i - 1 << std::endl;
                        std::cout << "alpha set mp " << aset_mp[i] << " " << i << "-";
                    }
                }
            }
            if(world.rank() == 0)
                std::cout << param.nmo_alpha - 1 << std::endl;

            if(param.nbeta && !param.spin_restricted){
                bmo_p = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta - 1)), 0.0, true);
                bmo_m = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta - 1)), 0.0, true);
                truncate(world, bmo_p);
                truncate(world, bmo_m);
                normalize(world, bmo_p);
                normalize(world, bmo_m);
                beps_mp = e(Slice(ncore, ncore + param.nmo_beta - 1));
                bocc_mp = tensorT(param.nmo_beta);
                for(int i = 0;i < param.nbeta;++i)
                    bocc_mp[i] = 1.0;

                bset_mp = std::vector<int>(param.nmo_beta,0);
                if(world.rank() == 0)
                    std::cout << " beta set mp " << 0 << " " << 0 << "-";

                for(int i = 1;i < param.nmo_beta;++i) {
                    bset_mp[i] = bset_mp[i - 1];
                    if(beps_mp[i] - beps_mp[i - 1] > 1.5 || bocc_mp[i] != 1.0){
                        ++(bset_mp[i]);
                        if(world.rank() == 0){
                            std::cout << i - 1 << std::endl;
                            std::cout << " beta set mp " << bset_mp[i] << " " << i << "-";
                        }
                    }
                }
                if(world.rank() == 0)
                    std::cout << param.nmo_beta - 1 << std::endl;
            }
        }
    }

    void initial_load_bal(World & world)
    {
        LoadBalanceDeux<3> lb(world);
        real_function_3d vnuc = potentialmanager->vnuclear();
        lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0));

        FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));
    }

    functionT make_density(World & world, const tensorT & occ, const vecfuncT & v)
    {
        vecfuncT vsq = square(world, v);
        compress(world, vsq);
        functionT rho = factoryT(world);
        rho.compress();
        for(unsigned int i = 0;i < vsq.size();++i){
            if(occ[i])
                rho.gaxpy(1.0, vsq[i], occ[i], false);

        }
        world.gop.fence();
        vsq.clear();
        return rho;
    }

    functionT make_density(World & world, const tensorT & occ, const cvecfuncT & v)
    {
      reconstruct(world, v); // For max parallelism
      std::vector<functionT> vsq(v.size());
      for (unsigned int i=0; i < v.size(); i++) {
          vsq[i] = abssq(v[i], false);
      }
      world.gop.fence();

      compress(world, vsq); // since will be using gaxpy for accumulation
      functionT rho = factoryT(world);
      rho.compress();

      for(unsigned int i = 0; i < vsq.size();++i) {
          if(occ[i])
              rho.gaxpy(1.0, vsq[i], occ[i], false);

      }
      world.gop.fence();
      vsq.clear();
      rho.truncate();

      return rho;
    }

    std::vector<poperatorT> make_bsh_operators(World & world, const tensorT & evals)
    {
        int nmo = evals.dim(0);
        std::vector<poperatorT> ops(nmo);
        double tol = FunctionDefaults<3>::get_thresh();
        for(int i = 0;i < nmo;++i){
            double eps = evals(i);
            if(eps > 0){
                if(world.rank() == 0){
                    print("bsh: warning: positive eigenvalue", i, eps);
                }
                eps = -0.1;
            }

            ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  param.lo, tol));
        }

        return ops;
    }

    vecfuncT apply_hf_exchange(World & world, const tensorT & occ, const vecfuncT & psi, const vecfuncT & f)
    {
        const bool same = (&psi == &f);
        int nocc = psi.size();
        int nf = f.size();
        double tol = FunctionDefaults<3>::get_thresh(); /// Important this is consistent with Coulomb
        vecfuncT Kf = zero_functions_compressed<double,3>(world, nf);
        reconstruct(world, psi);
        norm_tree(world, psi);
        if (!same) {
            reconstruct(world, f);
            norm_tree(world, f);
        }

//         // Smaller memory algorithm ... possible 2x saving using i-j sym
//         for(int i=0; i<nocc; ++i){
//             if(occ[i] > 0.0){
//                 vecfuncT psif = mul_sparse(world, psi[i], f, tol); /// was vtol
//                 truncate(world, psif);
//                 psif = apply(world, *coulop, psif);
//                 truncate(world, psif);
//                 psif = mul_sparse(world, psi[i], psif, tol); /// was vtol
//                 gaxpy(world, 1.0, Kf, occ[i], psif);
//             }
//         }

        // Larger memory algorithm ... use i-j sym if psi==f
        vecfuncT psif;
        for (int i=0; i<nocc; ++i) {
            int jtop = nf;
            if (same) jtop = i+1;
            for (int j=0; j<jtop; ++j) {
                psif.push_back(mul_sparse(psi[i], f[j], tol, false));
            }
        }

        world.gop.fence();
        truncate(world, psif);
        psif = apply(world, *coulop, psif);
        truncate(world, psif, tol);
        reconstruct(world, psif);
        norm_tree(world, psif);
        vecfuncT psipsif = zero_functions<double,3>(world, nf*nocc);
        int ij = 0;
        for (int i=0; i<nocc; ++i) {
            int jtop = nf;
            if (same) jtop = i+1;
            for (int j=0; j<jtop; ++j,++ij) {
                psipsif[i*nf+j] = mul_sparse(psif[ij],psi[i],false);
                if (same && i!=j) {
                    psipsif[j*nf+i] = mul_sparse(psif[ij],psi[j],false);
                }
            }
        }
        world.gop.fence();
        psif.clear();
        world.gop.fence();
        compress(world, psipsif);
        for (int i=0; i<nocc; ++i) {
            for (int j=0; j<nf; ++j) {
                Kf[j].gaxpy(1.0,psipsif[i*nf+j],occ[i],false);
            }
        }
        world.gop.fence();
        psipsif.clear();
        world.gop.fence();

        truncate(world, Kf, tol);
        return Kf;
    }

    // Used only for initial guess that is always spin-restricted LDA
    functionT make_lda_potential(World & world, const functionT & arho)
    {
        functionT vlda = copy(arho);
        vlda.reconstruct();
        vlda.unaryop(xc_lda_potential());
        return vlda;
    }


    functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin, int what)
    {
        return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
    }

    double make_dft_energy(World & world, const vecfuncT& vf, int ispin)
    {
        functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc, ispin), vf);
        return vlda.trace();
    }

    vecfuncT apply_potential(World & world, const tensorT & occ, const vecfuncT & amo,
                             const vecfuncT& vf, const vecfuncT& delrho, const functionT & vlocal, double & exc, int ispin)
    {
        functionT vloc = vlocal;
        exc = 0.0;

        //print("DFT", xc.is_dft(), "LDA", xc.is_lda(), "GGA", xc.is_gga(), "POLAR", xc.is_spin_polarized());
        if (xc.is_dft() && !(xc.hf_exchange_coefficient()==1.0)) {
            START_TIMER(world);
#ifdef MADNESS_HAS_LIBXC
            exc = make_dft_energy(world, vf, ispin);
#else
            if (ispin == 0) exc = make_dft_energy(world, vf, ispin);
#endif
            vloc = vloc + make_dft_potential(world, vf, ispin, 0);
            //print("VLOC1", vloc.trace(), vloc.norm2());

#ifdef MADNESS_HAS_LIBXC
            if (xc.is_gga() ) {
                if (world.rank() == 0) print(" WARNING GGA XC functionals must be used with caution in this version \n");
                real_function_3d vsig = make_dft_potential(world, vf, ispin, 1);
                //print("VSIG", vsig.trace(), vsig.norm2());
                real_function_3d vr(world);
                for (int axis=0; axis<3; axis++) {
                     vr += (*gradop[axis])(vsig);
                 //, int flagprint("VR", vr.trace(), vr.norm2());
                }
                vloc = vloc - vr;
            }
#endif
            END_TIMER(world, "DFT potential");
        }


        START_TIMER(world);
        vecfuncT Vpsi = mul_sparse(world, vloc, amo, vtol);
        END_TIMER(world, "V*psi");
	print_meminfo(world.rank(), "V*psi");
        if(xc.hf_exchange_coefficient()){
            START_TIMER(world);
            vecfuncT Kamo = apply_hf_exchange(world, occ, amo, amo);
            tensorT excv = inner(world, Kamo, amo);
            double exchf = 0.0;
            for(unsigned long i = 0;i < amo.size();++i){
                exchf -= 0.5 * excv[i] * occ[i];
            }
            if (!xc.is_spin_polarized()) exchf *= 2.0;
            gaxpy(world, 1.0, Vpsi, -xc.hf_exchange_coefficient(), Kamo);
            Kamo.clear();
            END_TIMER(world, "HF exchange");
            exc = exchf* xc.hf_exchange_coefficient() + exc;
        }
        potentialmanager->apply_nonlocal_potential(world, amo, Vpsi);

        if (param.core_type.substr(0,3) == "mcp") {
            START_TIMER(world);
            gaxpy(world, 1.0, Vpsi, 1.0, core_projection(world, amo));
            END_TIMER(world, "MCP Core Projector");
        }

        START_TIMER(world);
        truncate(world, Vpsi);
        END_TIMER(world, "Truncate Vpsi");
	print_meminfo(world.rank(), "Truncate Vpsi");
        world.gop.fence();
        return Vpsi;
    }

    vecfuncT apply_potential_response(World & world, const tensorT & occ, const vecfuncT & dmo,
                             const vecfuncT& vf, const functionT & vlocal, double & exc, int ispin)
    {
        functionT vloc = vlocal;
        exc = 0.0;

        //print("DFT", xc.is_dft(), "LDA", xc.is_lda(), "GGA", xc.is_gga(), "POLAR", xc.is_spin_polarized());
        if (xc.is_dft() && !(xc.hf_exchange_coefficient()==1.0)) {
            START_TIMER(world);
#ifdef MADNESS_HAS_LIBXC
            exc = make_dft_energy(world, vf, ispin);
#else
            if (ispin == 0) exc = make_dft_energy(world, vf, ispin);
#endif
            vloc = vloc + make_dft_potential(world, vf, ispin, 0);
            //print("VLOC1", vloc.trace(), vloc.norm2());

#ifdef MADNESS_HAS_LIBXC
            if (xc.is_gga() ) {
                if (world.rank() == 0) print(" WARNING GGA XC functionals must be used with caution in this version \n");
                real_function_3d vsig = make_dft_potential(world, vf, ispin, 1);
                //print("VSIG", vsig.trace(), vsig.norm2());
                real_function_3d vr(world);
                for (int axis=0; axis<3; axis++) {
                     vr += (*gradop[axis])(vsig);
                 //, int flagprint("VR", vr.trace(), vr.norm2());
                }
                vloc = vloc - vr;
            }
#endif
            END_TIMER(world, "DFT potential");
        }


        START_TIMER(world);
        vecfuncT Vdmo = mul_sparse(world, vloc, dmo, vtol);
        END_TIMER(world, "V*dmo");
	print_meminfo(world.rank(), "V*dmo");
        if(xc.hf_exchange_coefficient()){
            START_TIMER(world);
            vecfuncT Kdmo;
            if(ispin == 0)
                Kdmo= apply_hf_exchange(world, occ, amo, dmo);
            if(ispin == 1)
                Kdmo= apply_hf_exchange(world, occ, bmo, dmo);
            //tensorT excv = inner(world, Kdmo, dmo);
            //double exchf = 0.0;
            //for(unsigned long i = 0;i < dmo.size();++i){
            //    exchf -= 0.5 * excv[i] * occ[i];
            //}
            //if (!xc.is_spin_polarized()) exchf *= 2.0;
            gaxpy(world, 1.0, Vdmo, -xc.hf_exchange_coefficient(), Kdmo);
            Kdmo.clear();
            END_TIMER(world, "HF exchange");
            //exc = exchf* xc.hf_exchange_coefficient() + exc;
        }
        potentialmanager->apply_nonlocal_potential(world, dmo, Vdmo);

        if (param.core_type.substr(0,3) == "mcp") {
            START_TIMER(world);
            gaxpy(world, 1.0, Vdmo, 1.0, core_projection(world, dmo));
            END_TIMER(world, "MCP Core Projector");
        }

        START_TIMER(world);
        truncate(world, Vdmo);
        END_TIMER(world, "Truncate Vdmo");
	print_meminfo(world.rank(), "Truncate Vdmo");
        world.gop.fence();
        return Vdmo;
    }

    vecfuncT apply_potential_mp(World & world, const tensorT & occ, const vecfuncT & amo_mp,
                             const vecfuncT& vf, const vecfuncT& delrho, const functionT & vlocal,
                             double & exc, int ispin, int flag, int & axis, double & ep)
    {
        functionT vloc = vlocal;
        exc = 0.0;

        //print("DFT", xc.is_dft(), "LDA", xc.is_lda(), "GGA", xc.is_gga(), "POLAR", xc.is_spin_polarized());
        if (xc.is_dft() && !(xc.hf_exchange_coefficient()==1.0)) {
            START_TIMER(world);
#ifdef MADNESS_HAS_LIBXC
            exc = make_dft_energy(world, vf, ispin);
#else
            if (ispin == 0) exc = make_dft_energy(world, vf, ispin);
#endif
            vloc = vloc + make_dft_potential(world, vf, ispin, 0);
            //print("VLOC1", vloc.trace(), vloc.norm2());

#ifdef MADNESS_HAS_LIBXC
            if (xc.is_gga() ) {
                if (world.rank() == 0) print(" WARNING GGA XC functionals must be used with caution in this version \n");
                real_function_3d vsig = make_dft_potential(world, vf, ispin, 1);
                //print("VSIG", vsig.trace(), vsig.norm2());
                real_function_3d vr(world);
                for (int axis=0; axis<3; axis++) {
                     vr += (*gradop[axis])(vsig);
                 //print("VR", vr.trace(), vr.norm2());
                }
                vloc = vloc - vr;
            }
#endif
            END_TIMER(world, "DFT potential");
        }

        //vecfuncT perturbation = zero_functions<double,3>(world, 3) ;

                print("Vloc1", vloc.trace(), vloc.norm2());
        std::vector<int> f(3, 0);
        f[axis] = true;
        functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
        dipolefunc.scale(ep);
                print("Vdip", dipolefunc.trace(), dipolefunc.norm2());
        //perturbation = dipolefunc;
        //if (flag == 0) vloc = vloc + perturbation[i];
        if (flag == 0) vloc = vloc + dipolefunc;
        //else if(flag == 1) vloc = vloc - perturbation[i];
        else if(flag == 1) vloc = vloc - dipolefunc;
        //vloc = vloc - dipolefunc;
                print("Vloc2", vloc.trace(), vloc.norm2());

        START_TIMER(world);
        vecfuncT Vpsi = mul_sparse(world, vloc, amo_mp, vtol);
        END_TIMER(world, "V*psi");
	print_meminfo(world.rank(), "V*psi");
        if(xc.hf_exchange_coefficient()){
            START_TIMER(world);
            vecfuncT Kamo = apply_hf_exchange(world, occ, amo_mp, amo_mp);
            tensorT excv = inner(world, Kamo, amo_mp);
            double exchf = 0.0;
            for(unsigned long i = 0;i < amo_mp.size();++i){
                exchf -= 0.5 * excv[i] * occ[i];
            }
            if (!xc.is_spin_polarized()) exchf *= 2.0;
            gaxpy(world, 1.0, Vpsi, -xc.hf_exchange_coefficient(), Kamo);
            Kamo.clear();
            END_TIMER(world, "HF exchange");
            exc = exchf* xc.hf_exchange_coefficient() + exc;
        }
        potentialmanager->apply_nonlocal_potential(world, amo_mp, Vpsi);

        if (param.core_type.substr(0,3) == "mcp") {
            START_TIMER(world);
            gaxpy(world, 1.0, Vpsi, 1.0, core_projection(world, amo_mp));
            END_TIMER(world, "MCP Core Projector");
        }

        START_TIMER(world);
        truncate(world, Vpsi);
        END_TIMER(world, "Truncate Vpsi");
	print_meminfo(world.rank(), "Truncate Vpsi");
        world.gop.fence();
        return Vpsi;
    }

    tensorT derivatives(World & world)
    {
        START_TIMER(world);

        functionT rho = make_density(world, aocc, amo);
        functionT brho = rho;
        if (!param.spin_restricted) brho = make_density(world, bocc, bmo);
        rho.gaxpy(1.0, brho, 1.0);

        vecfuncT dv(molecule.natom() * 3);
        vecfuncT du = zero_functions<double,3>(world, molecule.natom() * 3);
        tensorT rc(molecule.natom() * 3);
        for(int atom = 0;atom < molecule.natom();++atom){
            for(int axis = 0;axis < 3;++axis){
                functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
                dv[atom * 3 + axis] = functionT(factoryT(world).functor(func).nofence().truncate_on_project());
                if (param.core_type != "" && molecule.is_potential_defined_atom(atom)) {
                    // core potential contribution
                    func = functorT(new CorePotentialDerivativeFunctor(molecule, atom, axis));
                    du[atom * 3 + axis] = functionT(factoryT(world).functor(func).truncate_on_project());

                    // core projector contribution
                    rc[atom * 3 + axis] = potentialmanager->core_projector_derivative(world, amo, aocc, atom, axis);
                    if (!param.spin_restricted) {
                        if (param.nbeta) rc[atom * 3 + axis] += potentialmanager->core_projector_derivative(world, bmo, bocc, atom, axis);
                    }
                    else {
                        rc[atom * 3 + axis] *= 2 * 2;
                            // because of 2 electrons in each valence orbital bra+ket
                    }
                }
            }
        }

        world.gop.fence();
        tensorT r = inner(world, rho, dv);
        world.gop.fence();
        tensorT ru = inner(world, rho, du);
        dv.clear();
        du.clear();
        world.gop.fence();
        tensorT ra(r.size());
        for(int atom = 0;atom < molecule.natom();++atom){
            for(int axis = 0;axis < 3;++axis){
                ra[atom * 3 + axis] = molecule.nuclear_repulsion_derivative(atom, axis);
            }
        }
        //if (world.rank() == 0) print("derivatives:\n", r, ru, rc, ra);
        r +=  ra + ru + rc;
        END_TIMER(world,"derivatives");

        if (world.rank() == 0) {
            print("\n Derivatives (a.u.)\n -----------\n");
            print("  atom        x            y            z          dE/dx        dE/dy        dE/dz");
            print(" ------ ------------ ------------ ------------ ------------ ------------ ------------");
            for (int i=0; i<molecule.natom(); ++i) {
                const Atom& atom = molecule.get_atom(i);
                printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                       i, atom.x, atom.y, atom.z,
                       r[i*3+0], r[i*3+1], r[i*3+2]);
            }
        }
        return r;
    }

    tensorT dipole(World & world)
    {
        START_TIMER(world);
        tensorT mu(3);
        for (unsigned int axis=0; axis<3; ++axis) {
            std::vector<int> x(3, 0);
            x[axis] = true;
            functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(x)));
            functionT rho = make_density(world, aocc, amo);
            if (!param.spin_restricted) {
                if (param.nbeta) rho += make_density(world, bocc, bmo);
            }
            else {
                rho.scale(2.0);
            }
            mu[axis] = -dipolefunc.inner(rho);
            mu[axis] += molecule.nuclear_dipole(axis);
        }

        if (world.rank() == 0) {
            print("\n Dipole Moment (a.u.)\n -----------\n");
            print("     x: ", mu[0]);
            print("     y: ", mu[1]);
            print("     z: ", mu[2]);
            print(" Total Dipole Moment: ", mu.normf());
        }
        END_TIMER(world, "dipole");

        return mu;
    }

    void vector_stats(const std::vector<double> & v, double & rms, double & maxabsval)
    {
        rms = 0.0;
        maxabsval = v[0];
        for(unsigned int i = 0;i < v.size();++i){
            rms += v[i] * v[i];
            maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
        }
        rms = sqrt(rms / v.size());
    }

    vecfuncT compute_residual(World & world, tensorT & occ, tensorT & fock, const vecfuncT & psi, vecfuncT & Vpsi, double & err)
    {
        double trantol = vtol / std::min(30.0, double(psi.size()));
        int nmo = psi.size();

        tensorT eps(nmo);
        for(int i = 0;i < nmo;++i){
            eps(i) = std::min(-0.05, fock(i, i));
            fock(i, i) -= eps(i);
        }
        vecfuncT fpsi = transform(world, psi, fock, trantol, true);

        for(int i = 0;i < nmo;++i){ // Undo the damage
            fock(i, i) += eps(i);
        }

        gaxpy(world, 1.0, Vpsi, -1.0, fpsi);
        fpsi.clear();
        std::vector<double> fac(nmo, -2.0);
        scale(world, Vpsi, fac);
        std::vector<poperatorT> ops = make_bsh_operators(world, eps);
        set_thresh(world, Vpsi, FunctionDefaults<3>::get_thresh());
        if(world.rank() == 0)
            std::cout << "entering apply\n";

        START_TIMER(world);
        vecfuncT new_psi = apply(world, ops, Vpsi);
        END_TIMER(world, "Apply BSH");
        ops.clear();
        Vpsi.clear();
        world.gop.fence();

        // Thought it was a bad idea to truncate *before* computing the residual
        // but simple tests suggest otherwise ... no more iterations and
        // reduced iteration time from truncating.
        START_TIMER(world);
        truncate(world, new_psi);
        END_TIMER(world, "Truncate new psi");

        vecfuncT r = sub(world, psi, new_psi);
        std::vector<double> rnorm = norm2s(world, r);
        if (world.rank() == 0) print("residuals", rnorm);
        double rms, maxval;
        vector_stats(rnorm, rms, maxval);
        err = maxval;
        if(world.rank() == 0)
            print("BSH residual: rms", rms, "   max", maxval);

        return r;
    }

    tensorT make_fock_matrix(World & world, const vecfuncT & psi, const vecfuncT & Vpsi, const tensorT & occ, double & ekinetic)
    {
        START_TIMER(world);
        tensorT pe = matrix_inner(world, Vpsi, psi, true);
        END_TIMER(world, "PE matrix");
        START_TIMER(world);
        tensorT ke = kinetic_energy_matrix(world, psi);
        END_TIMER(world, "KE matrix");
	START_TIMER(world);
        int nocc = occ.size();
        ekinetic = 0.0;
        for(int i = 0;i < nocc;++i){
            ekinetic += occ[i] * ke(i, i);
        }
        ke += pe;
        pe = tensorT();
        ke.gaxpy(0.5, transpose(ke), 0.5);
	END_TIMER(world, "Make fock matrix rest");
        return ke;
    }

    /// Compute the two-electron integrals over the provided set of orbitals

    /// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
    /// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
    Tensor<double> twoint(World& world, const vecfuncT& psi) {
        double tol = FunctionDefaults<3>::get_thresh(); /// Important this is consistent with Coulomb
        reconstruct(world, psi);
        norm_tree(world, psi);

        // Efficient version would use mul_sparse vector interface
        vecfuncT pairs;
        for (unsigned int i=0; i<psi.size(); ++i) {
            for (unsigned int j=0; j<=i; ++j) {
                pairs.push_back(mul_sparse(psi[i], psi[j], tol, false));
            }
        }

        world.gop.fence();
        truncate(world, pairs);
        vecfuncT Vpairs = apply(world, *coulop, pairs);

        return matrix_inner(world, pairs, Vpairs, true);
    }

    tensorT matrix_exponential(const tensorT& A) {
        const double tol = 1e-13;
        MADNESS_ASSERT(A.dim((0) == A.dim(1)));

        // Scale A by a power of 2 until it is "small"
        double anorm = A.normf();
        int n = 0;
        double scale = 1.0;
        while (anorm*scale > 0.1) {
            ++n;
            scale *= 0.5;
        }
        tensorT B = scale*A;    // B = A*2^-n

        // Compute exp(B) using Taylor series
        tensorT expB = tensorT(2, B.dims());
        for (int i=0; i<expB.dim(0); ++i) expB(i,i) = 1.0;

        int k = 1;
        tensorT term = B;
        while (term.normf() > tol) {
            expB += term;
            term = inner(term,B);
            ++k;
            term.scale(1.0/k);
        }

        // Repeatedly square to recover exp(A)
        while (n--) {
            expB = inner(expB,expB);
        }

        return expB;
    }

    tensorT diag_fock_matrix(World & world, tensorT& fock, vecfuncT & psi, vecfuncT & Vpsi, tensorT & evals, const tensorT & occ, double thresh)
    {
        long nmo = psi.size();
        tensorT overlap = matrix_inner(world, psi, psi, true);

        START_TIMER(world);
        tensorT U;
        sygvp(world, fock, overlap, 1, U, evals);
        END_TIMER(world, "Diagonalization Fock-mat w sygv");

        START_TIMER(world);
        // Within blocks with the same occupation number attempt to
        // keep orbitals in the same order (to avoid confusing the
        // non-linear solver).
	// !!!!!!!!!!!!!!!!! NEED TO RESTRICT TO OCCUPIED STATES?
	bool switched = true;
	while (switched) {
	  switched = false;
	  for (int i=0; i<nmo; i++) {
	    for (int j=i+1; j<nmo; j++) {
	      if (occ(i) == occ(j)) {
		double sold = U(i,i)*U(i,i) + U(j,j)*U(j,j);
		double snew = U(i,j)*U(i,j) + U(j,i)*U(j,i);
		if (snew > sold) {
		  tensorT tmp = copy(U(_,i));
		  U(_,i) = U(_,j);
		  U(_,j) = tmp;
		  std::swap(evals[i],evals[j]);
		  switched = true;
		}
	      }
	    }
	  }
	}
        // Fix phases.
        for (long i=0; i<nmo; ++i) {
            if (U(i,i) < 0.0) U(_,i).scale(-1.0);
        }

        // Rotations between effectively degenerate states confound
        // the non-linear equation solver ... undo these rotations
        long ilo = 0; // first element of cluster
        while (ilo < nmo-1) {
            long ihi = ilo;
            while (fabs(evals[ilo]-evals[ihi+1]) < thresh*10.0*std::max(fabs(evals[ilo]),1.0)) {
                ++ihi;
                if (ihi == nmo-1) break;
            }
            long nclus = ihi - ilo + 1;
            if (nclus > 1) {
                //print("   found cluster", ilo, ihi);
                tensorT q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
                //print(q);
                // Special code just for nclus=2
                // double c = 0.5*(q(0,0) + q(1,1));
                // double s = 0.5*(q(0,1) - q(1,0));
                // double r = sqrt(c*c + s*s);
                // c /= r;
                // s /= r;
                // q(0,0) = q(1,1) = c;
                // q(0,1) = -s;
                // q(1,0) = s;

                // Iteratively construct unitary rotation by
                // exponentiating the antisymmetric part of the matrix
                // ... is quadratically convergent so just do 3
                // iterations
                tensorT rot = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(q,rot);
                tensorT rot2 = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(q,rot2);
                tensorT rot3 = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(rot,inner(rot2,rot3));
                U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
            }
            ilo = ihi+1;
        }

	// if (world.rank() == 0) {
	//   print("Fock");
	//   print(fock);
	//   print("Evec");
	//   print(U);;
	//   print("Eval");
	//   print(evals);
	// }

        world.gop.broadcast(U.ptr(), U.size(), 0);
        world.gop.broadcast(evals.ptr(), evals.size(), 0);

        fock = 0;
        for (unsigned int i=0; i<psi.size(); ++i) fock(i,i) = evals(i);

        Vpsi = transform(world, Vpsi, U, vtol / std::min(30.0, double(psi.size())), false);
        psi = transform(world, psi, U, FunctionDefaults<3>::get_thresh() / std::min(30.0, double(psi.size())), true);
        truncate(world, Vpsi, vtol, false);
        truncate(world, psi);
        normalize(world, psi);

        END_TIMER(world, "Diagonalization rest");
        return U;
    }

    void loadbal(World & world, functionT & arho, functionT & brho, functionT & arho_old, functionT & brho_old, subspaceT & subspace)
    {
        if(world.size() == 1)
            return;

        LoadBalanceDeux<3> lb(world);
        real_function_3d vnuc = potentialmanager->vnuclear();
        lb.add_tree(vnuc, lbcost<double,3>(vnucextra*1.0, vnucextra*8.0), false);
        lb.add_tree(arho, lbcost<double,3>(1.0, 8.0), false);
        for(unsigned int i = 0;i < amo.size();++i){
            lb.add_tree(amo[i], lbcost<double,3>(1.0, 8.0), false);
        }
        if(param.nbeta && !param.spin_restricted){
            lb.add_tree(brho, lbcost<double,3>(1.0, 8.0), false);
            for(unsigned int i = 0;i < bmo.size();++i){
                lb.add_tree(bmo[i], lbcost<double,3>(1.0, 8.0), false);
            }
        }
	world.gop.fence();

        FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0)); // 6.0 needs retuning after vnucextra
    }

    void rotate_subspace(World& world, const tensorT& U, subspaceT& subspace, int lo, int nfunc, double trantol) {
        for (unsigned int iter=0; iter<subspace.size(); ++iter) {
            vecfuncT& v = subspace[iter].first;
            vecfuncT& r = subspace[iter].second;
            transform(world, vecfuncT(&v[lo],&v[lo+nfunc]), U, trantol, false);
            transform(world, vecfuncT(&r[lo],&r[lo+nfunc]), U, trantol, true);
        }
    }

    void update_response_subspace(World & world,
                         vecfuncT & ax, vecfuncT & ay,
                         vecfuncT & bx, vecfuncT & by,
                         vecfuncT & rax, vecfuncT & ray,
                         vecfuncT & rbx, vecfuncT & rby,
                         subspaceT & subspace, tensorT & Q, double & update_residual)
    {
        vecfuncT vm = ax;
        vm.insert(vm.end(), ay.begin(), ay.end());

        vecfuncT rm = rax;
        rm.insert(rm.end(), ray.begin(), ray.end());

        if(param.nbeta != 0 && !param.spin_restricted){
            vm.insert(vm.end(), bx.begin(), bx.end());
            vm.insert(vm.end(), by.begin(), by.end());
            rm.insert(rm.end(), rbx.begin(), rbx.end());
            rm.insert(rm.end(), rby.begin(), rby.end());
        }

        compress(world, vm, false);
        compress(world, rm, false);
        world.gop.fence();
        subspace.push_back(pairvecfuncT(vm, rm));
        int m = subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for(int s = 0;s < m;++s){
            const vecfuncT & vs = subspace[s].first;
            const vecfuncT & rs = subspace[s].second;
            for(unsigned int i = 0;i < vm.size();++i){
                ms[s] += vm[i].inner_local(rs[i]);
                sm[s] += vs[i].inner_local(rm[i]);
            }
        }

        world.gop.sum(ms.ptr(), m);
        world.gop.sum(sm.ptr(), m);
        tensorT newQ(m, m);
        if(m > 1)
            newQ(Slice(0, -2), Slice(0, -2)) = Q;

        newQ(m - 1, _) = ms;
        newQ(_, m - 1) = sm;
        Q = newQ;
        //if (world.rank() == 0) { print("kain Q"); print(Q); }
        tensorT c;
        if(world.rank() == 0){
            double rcond = 1e-12;
            while(1){
                c = KAIN(Q, rcond);
                //if (world.rank() == 0) print("kain c:", c);
                if(std::abs(c[m - 1]) < 3.0){
                    break;
                } else if(rcond < 0.01){
                    print("Increasing subspace singular value threshold ", c[m - 1], rcond);
                    rcond *= 100;
                } else {
                    print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m - 1] = 1.0;
                    break;
                }
            }
        }

        world.gop.broadcast_serializable(c, 0);
        if(world.rank() == 0){
            print("Response Subspace solution", c);
        }
        START_TIMER(world);
        vecfuncT ax_new = zero_functions_compressed<double,3>(world, ax.size());
        vecfuncT ay_new = zero_functions_compressed<double,3>(world, ay.size());
        vecfuncT bx_new = zero_functions_compressed<double,3>(world, bx.size());
        vecfuncT by_new = zero_functions_compressed<double,3>(world, by.size());

        for(unsigned int m = 0;m < subspace.size();++m){
            const vecfuncT & vm = subspace[m].first;
            const vecfuncT & rm = subspace[m].second;
            const vecfuncT vmax(vm.begin(), vm.begin() + ax.size());
            const vecfuncT rmax(rm.begin(), rm.begin() + rax.size());
            const vecfuncT vmay(vm.begin() + ax.size(), vm.begin() + ax.size() + ay.size());
            const vecfuncT rmay(rm.begin() + rax.size(), rm.begin() + rax.size() + ray.size());
            gaxpy(world, 1.0, ax_new, c(m), vmax, false);
            gaxpy(world, 1.0, ax_new, -c(m), rmax, false);
            gaxpy(world, 1.0, ay_new, c(m), vmay, false);
            gaxpy(world, 1.0, ay_new, -c(m), rmay, false);
            //if(param.nbeta != 0 && !param.spin_restricted){
                const vecfuncT vmbx(vm.end() - by.size() - bx.size(), vm.end() - by.size());
                const vecfuncT rmbx(rm.end() - rby.size() - rbx.size(), rm.end() - rby.size());
                const vecfuncT vmby(vm.end() - by.size(), vm.end());
                const vecfuncT rmby(rm.end() - rby.size(), rm.end());
                gaxpy(world, 1.0, bx_new, c(m), vmbx, false);
                gaxpy(world, 1.0, bx_new, -c(m), rmbx, false);
                gaxpy(world, 1.0, by_new, c(m), vmby, false);
                gaxpy(world, 1.0, by_new, -c(m), rmby, false);
            //}
        }
        world.gop.fence();
        END_TIMER(world, "Subspace transform");
        if(param.maxsub <= 1){
            subspace.clear();
        } else if(subspace.size() == param.maxsub){
            subspace.erase(subspace.begin());
            Q = Q(Slice(1, -1), Slice(1, -1));
        }

        std::vector<double> axnorm = norm2s(world, sub(world, ax, ax_new));
        std::vector<double> aynorm = norm2s(world, sub(world, ay, ay_new));
        std::vector<double> bxnorm = norm2s(world, sub(world, bx, bx_new));
        std::vector<double> bynorm = norm2s(world, sub(world, by, by_new));
        int nres = 0;
        for(unsigned int i = 0;i < ax.size();++i){
            if(axnorm[i] > param.maxrotn){
                double s = param.maxrotn / axnorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for alpha orbitals:");

                    printf(" %d", i);
                }
                ax_new[i].gaxpy(s, ax[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");
        
        nres = 0;
        for(unsigned int i = 0;i < ay.size();++i){
            if(aynorm[i] > param.maxrotn){
                double s = param.maxrotn / aynorm[i];
                nres++;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for alpha orbitals:");

                    printf(" %d", i);
                }
                ay_new[i].gaxpy(s, ay[i], 1.0 - s, false);
            }
        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        //if(param.nbeta != 0 && !param.spin_restricted){
            nres = 0;
            for(unsigned int i = 0;i < bx.size();++i){
                if(bxnorm[i] > param.maxrotn){
                    double s = param.maxrotn / bxnorm[i];
                    nres++;
                    if(world.rank() == 0){
                        if(nres == 1)
                            printf("  restricting step for  beta orbitals:");

                        printf(" %d", i);
                    }
                    bx_new[i].gaxpy(s, bx[i], 1.0 - s, false);
                }
            }
            if(nres > 0 && world.rank() == 0)
                printf("\n");

            nres = 0;
            for(unsigned int i = 0;i < by.size();++i){
                if(bynorm[i] > param.maxrotn){
                    double s = param.maxrotn / bynorm[i];
                    ++nres;
                    if(world.rank() == 0){
                        if(nres == 1)
                            printf("  restricting step for  beta orbitals:");

                        printf(" %d", i);
                    }
                    by_new[i].gaxpy(s, by[i], 1.0 - s, false);
                }

            }
            if(nres > 0 && world.rank() == 0)
                printf("\n");
        //}

        world.gop.fence();
        double rms, maxval_x, maxval_y, maxval_b;
        vector_stats(axnorm, rms, maxval_x);
        vector_stats(aynorm, rms, maxval_y);

        update_residual = std::max(maxval_x, maxval_y);

        if(bxnorm.size()){
            vector_stats(bxnorm, rms, maxval_x);
            vector_stats(bynorm, rms, maxval_y);
            
            maxval_b = std::max(maxval_x, maxval_y);
            update_residual = std::max(update_residual, maxval_b);
        }
        //START_TIMER(world);
        //double trantol = vtol / std::min(30.0, double(ax.size()));
        //normalize(world, ax_new);
        //normalize(world, ay_new);
        //ax_new = transform(world, ax_new, Q3(matrix_inner(world, ax_new, ax_new)), trantol, true);
        //ay_new = transform(world, ay_new, Q3(matrix_inner(world, ay_new, ay_new)), trantol, true);
        truncate(world, ax_new);
        truncate(world, ay_new);
        //normalize(world, ax_new);
        //normalize(world, ay_new);
        if(param.nbeta != 0  && !param.spin_restricted){
            //normalize(world, bx_new);
            //normalize(world, by_new);
            //bx_new = transform(world, bx_new, Q3(matrix_inner(world, bx_new, bx_new)), trantol, true);
            //by_new = transform(world, by_new, Q3(matrix_inner(world, by_new, by_new)), trantol, true);
            truncate(world, bx_new);
            truncate(world, by_new);
            //normalize(world, bx_new);
            //normalize(world, by_new);
        }
        //END_TIMER(world, "Orthonormalize");
        ax = ax_new;
        ay = ay_new;
        bx = bx_new;
        by = by_new;
    }

    void update_subspace(World & world,
                         vecfuncT & Vpsia, vecfuncT & Vpsib,
                         tensorT & focka, tensorT & fockb,
                         subspaceT & subspace, tensorT & Q,
                         double & bsh_residual, double & update_residual)
    {
        double aerr = 0.0, berr = 0.0;
        vecfuncT vm = amo;

        // Orbitals with occ!=1.0 exactly must be solved for as eigenfunctions
        // so zero out off diagonal lagrange multipliers
        for (int i=0; i<param.nmo_alpha; i++) {
            if (aocc[i] != 1.0) {
                double tmp = focka(i,i);
                focka(i,_) = 0.0;
                focka(_,i) = 0.0;
                focka(i,i) = tmp;
            }
        }

        vecfuncT rm = compute_residual(world, aocc, focka, amo, Vpsia, aerr);
        if(param.nbeta != 0 && !param.spin_restricted){
            for (int i=0; i<param.nmo_beta; i++) {
                if (bocc[i] != 1.0) {
                    double tmp = fockb(i,i);
                    fockb(i,_) = 0.0;
                    fockb(_,i) = 0.0;
                    fockb(i,i) = tmp;
                }
            }

            vecfuncT br = compute_residual(world, bocc, fockb, bmo, Vpsib, berr);
            vm.insert(vm.end(), bmo.begin(), bmo.end());
            rm.insert(rm.end(), br.begin(), br.end());
        }
        bsh_residual = std::max(aerr, berr);
        world.gop.broadcast(bsh_residual, 0);
        compress(world, vm, false);
        compress(world, rm, false);
        world.gop.fence();
        subspace.push_back(pairvecfuncT(vm, rm));
        int m = subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for(int s = 0;s < m;++s){
            const vecfuncT & vs = subspace[s].first;
            const vecfuncT & rs = subspace[s].second;
            for(unsigned int i = 0;i < vm.size();++i){
                ms[s] += vm[i].inner_local(rs[i]);
                sm[s] += vs[i].inner_local(rm[i]);
            }
        }

        world.gop.sum(ms.ptr(), m);
        world.gop.sum(sm.ptr(), m);
        tensorT newQ(m, m);
        if(m > 1)
            newQ(Slice(0, -2), Slice(0, -2)) = Q;

        newQ(m - 1, _) = ms;
        newQ(_, m - 1) = sm;
        Q = newQ;
        //if (world.rank() == 0) { print("kain Q"); print(Q); }
        tensorT c;
        if(world.rank() == 0){
            double rcond = 1e-12;
            while(1){
                c = KAIN(Q, rcond);
                //if (world.rank() == 0) print("kain c:", c);
                if(std::abs(c[m - 1]) < 3.0){
                    break;
                } else  if(rcond < 0.01){
                    print("Increasing subspace singular value threshold ", c[m - 1], rcond);
                    rcond *= 100;
                } else {
                    print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m - 1] = 1.0;
                    break;
                }
            }
        }

        world.gop.broadcast_serializable(c, 0);
        if(world.rank() == 0){
            print("Subspace solution", c);
        }
        START_TIMER(world);
        vecfuncT amo_new = zero_functions_compressed<double,3>(world, amo.size());
        vecfuncT bmo_new = zero_functions_compressed<double,3>(world, bmo.size());
        for(unsigned int m = 0;m < subspace.size();++m){
            const vecfuncT & vm = subspace[m].first;
            const vecfuncT & rm = subspace[m].second;
            const vecfuncT vma(vm.begin(), vm.begin() + amo.size());
            const vecfuncT rma(rm.begin(), rm.begin() + amo.size());
            const vecfuncT vmb(vm.end() - bmo.size(), vm.end());
            const vecfuncT rmb(rm.end() - bmo.size(), rm.end());
            gaxpy(world, 1.0, amo_new, c(m), vma, false);
            gaxpy(world, 1.0, amo_new, -c(m), rma, false);
            gaxpy(world, 1.0, bmo_new, c(m), vmb, false);
            gaxpy(world, 1.0, bmo_new, -c(m), rmb, false);
        }
        world.gop.fence();
        END_TIMER(world, "Subspace transform");
        if(param.maxsub <= 1){
            subspace.clear();
        } else if(subspace.size() == param.maxsub){
            subspace.erase(subspace.begin());
            Q = Q(Slice(1, -1), Slice(1, -1));
        }

        std::vector<double> anorm = norm2s(world, sub(world, amo, amo_new));
        std::vector<double> bnorm = norm2s(world, sub(world, bmo, bmo_new));
        int nres = 0;
        for(unsigned int i = 0;i < amo.size();++i){
            if(anorm[i] > param.maxrotn){
                double s = param.maxrotn / anorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for alpha orbitals:");

                    printf(" %d", i);
                }
                amo_new[i].gaxpy(s, amo[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        nres = 0;
        for(unsigned int i = 0;i < bmo.size();++i){
            if(bnorm[i] > param.maxrotn){
                double s = param.maxrotn / bnorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for  beta orbitals:");

                    printf(" %d", i);
                }
                bmo_new[i].gaxpy(s, bmo[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        world.gop.fence();
        double rms, maxval;
        vector_stats(anorm, rms, maxval);
        if(world.rank() == 0)
            print("Norm of vector changes alpha: rms", rms, "   max", maxval);

        update_residual = maxval;
        if(bnorm.size()){
            vector_stats(bnorm, rms, maxval);
            if(world.rank() == 0)
                print("Norm of vector changes  beta: rms", rms, "   max", maxval);

            update_residual = std::max(update_residual, maxval);
        }
        START_TIMER(world);
        double trantol = vtol / std::min(30.0, double(amo.size()));
        normalize(world, amo_new);
        amo_new = transform(world, amo_new, Q3(matrix_inner(world, amo_new, amo_new)), trantol, true);
        truncate(world, amo_new);
        normalize(world, amo_new);
        if(param.nbeta != 0  && !param.spin_restricted){
            normalize(world, bmo_new);
            bmo_new = transform(world, bmo_new, Q3(matrix_inner(world, bmo_new, bmo_new)), trantol, true);
            truncate(world, bmo_new);
            normalize(world, bmo_new);
        }
        END_TIMER(world, "Orthonormalize");
        amo = amo_new;
        bmo = bmo_new;
    }

    void update_subspace_mp(World & world,
                         vecfuncT & Vpsia, vecfuncT & Vpsib,
                         tensorT & focka, tensorT & fockb,
                         subspaceT & subspace, tensorT & Q,
                         double & bsh_residual, double & update_residual,
                         int & mflag,
                         vecfuncT & amo_mp, vecfuncT & bmo_mp)
    {
        double aerr = 0.0, berr = 0.0;
        #if 0
        
        if(mflag == 1) {
            amo_mp = copy(world, amo_m);
            bmo_mp = copy(world, bmo_m);
        }
        else {
            amo_mp = copy(world, amo_p);
            bmo_mp = copy(world, bmo_p);
        }
        #endif
        vecfuncT vm = amo_mp;

        // Orbitals with occ!=1.0 exactly must be solved for as eigenfunctions
        // so zero out off diagonal lagrange multipliers
        for (int i=0; i<param.nmo_alpha; i++) {
            if (aocc_mp[i] != 1.0) {
                double tmp = focka(i,i);
                focka(i,_) = 0.0;
                focka(_,i) = 0.0;
                focka(i,i) = tmp;
            }
        }

        vecfuncT rm = compute_residual(world, aocc_mp, focka, amo_mp, Vpsia, aerr);
        if(param.nbeta != 0 && !param.spin_restricted){
            for (int i=0; i<param.nmo_beta; i++) {
                if (bocc_mp[i] != 1.0) {
                    double tmp = fockb(i,i);
                    fockb(i,_) = 0.0;
                    fockb(_,i) = 0.0;
                    fockb(i,i) = tmp;
                }
            }

            vecfuncT br = compute_residual(world, bocc, fockb, bmo_mp, Vpsib, berr);
            vm.insert(vm.end(), bmo_mp.begin(), bmo_mp.end());
            rm.insert(rm.end(), br.begin(), br.end());
        }
        bsh_residual = std::max(aerr, berr);
        world.gop.broadcast(bsh_residual, 0);
        compress(world, vm, false);
        compress(world, rm, false);
        world.gop.fence();
        subspace.push_back(pairvecfuncT(vm, rm));
        int m = subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for(int s = 0;s < m;++s){
            const vecfuncT & vs = subspace[s].first;
            const vecfuncT & rs = subspace[s].second;
            for(unsigned int i = 0;i < vm.size();++i){
                ms[s] += vm[i].inner_local(rs[i]);
                sm[s] += vs[i].inner_local(rm[i]);
            }
        }

        world.gop.sum(ms.ptr(), m);
        world.gop.sum(sm.ptr(), m);
        tensorT newQ(m, m);
        if(m > 1)
            newQ(Slice(0, -2), Slice(0, -2)) = Q;

        newQ(m - 1, _) = ms;
        newQ(_, m - 1) = sm;
        Q = newQ;
        //if (world.rank() == 0) { print("kain Q"); print(Q); }
        tensorT c;
        if(world.rank() == 0){
            double rcond = 1e-12;
            while(1){
                c = KAIN(Q, rcond);
                //if (world.rank() == 0) print("kain c:", c);
                if(std::abs(c[m - 1]) < 3.0){
                    break;
                } else  if(rcond < 0.01){
                    print("Increasing subspace singular value threshold ", c[m - 1], rcond);
                    rcond *= 100;
                } else {
                    print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m - 1] = 1.0;
                    break;
                }
            }
        }

        world.gop.broadcast_serializable(c, 0);
        if(world.rank() == 0){
            print("Subspace solution", c);
        }
        START_TIMER(world);
        vecfuncT amo_new = zero_functions_compressed<double,3>(world, amo_mp.size());
        vecfuncT bmo_new = zero_functions_compressed<double,3>(world, bmo_mp.size());
        for(unsigned int m = 0;m < subspace.size();++m){
            const vecfuncT & vm = subspace[m].first;
            const vecfuncT & rm = subspace[m].second;
            const vecfuncT vma(vm.begin(), vm.begin() + amo_mp.size());
            const vecfuncT rma(rm.begin(), rm.begin() + amo_mp.size());
            const vecfuncT vmb(vm.end() - bmo_mp.size(), vm.end());
            const vecfuncT rmb(rm.end() - bmo_mp.size(), rm.end());
            gaxpy(world, 1.0, amo_new, c(m), vma, false);
            gaxpy(world, 1.0, amo_new, -c(m), rma, false);
            gaxpy(world, 1.0, bmo_new, c(m), vmb, false);
            gaxpy(world, 1.0, bmo_new, -c(m), rmb, false);
        }
        world.gop.fence();
        END_TIMER(world, "Subspace transform");
        if(param.maxsub <= 1){
            subspace.clear();
        } else if(subspace.size() == param.maxsub){
            subspace.erase(subspace.begin());
            Q = Q(Slice(1, -1), Slice(1, -1));
        }

        std::vector<double> anorm = norm2s(world, sub(world, amo_mp, amo_new));
        std::vector<double> bnorm = norm2s(world, sub(world, bmo_mp, bmo_new));
        int nres = 0;
        for(unsigned int i = 0;i < amo_mp.size();++i){
            if(anorm[i] > param.maxrotn){
                double s = param.maxrotn / anorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for alpha orbitals:");

                    printf(" %d", i);
                }
                amo_new[i].gaxpy(s, amo_mp[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        nres = 0;
        for(unsigned int i = 0;i < bmo_mp.size();++i){
            if(bnorm[i] > param.maxrotn){
                double s = param.maxrotn / bnorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for  beta orbitals:");

                    printf(" %d", i);
                }
                bmo_new[i].gaxpy(s, bmo_mp[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        world.gop.fence();
        double rms, maxval;
        vector_stats(anorm, rms, maxval);
        if(world.rank() == 0)
            print("Norm of vector changes alpha: rms", rms, "   max", maxval);

        update_residual = maxval;
        if(bnorm.size()){
            vector_stats(bnorm, rms, maxval);
            if(world.rank() == 0)
                print("Norm of vector changes  beta: rms", rms, "   max", maxval);

            update_residual = std::max(update_residual, maxval);
        }
        START_TIMER(world);
        double trantol = vtol / std::min(30.0, double(amo_mp.size()));
        normalize(world, amo_new);
        amo_new = transform(world, amo_new, Q3(matrix_inner(world, amo_new, amo_new)), trantol, true);
        truncate(world, amo_new);
        normalize(world, amo_new);
        if(param.nbeta != 0  && !param.spin_restricted){
            normalize(world, bmo_new);
            bmo_new = transform(world, bmo_new, Q3(matrix_inner(world, bmo_new, bmo_new)), trantol, true);
            truncate(world, bmo_new);
            normalize(world, bmo_new);
        }
        END_TIMER(world, "Orthonormalize");
        amo_mp = amo_new;
        bmo_mp = bmo_new;
    }

//vama1//    template <typename Func>
//vama1//    void propagate(World& world, const Func& Vext, int step0)
//vama1    void propagate(World& world, double omega, int step0)
//vama1    {
//vama1      // Load molecular orbitals
//vama1      set_protocol(world,1e-4);
//vama1      make_nuclear_potential(world);
//vama1      initial_load_bal(world);
//vama1      load_mos(world);
//vama1
//vama1      int nstep = 1000;
//vama1      double time_step = 0.05;
//vama1
//vama1      double strength = 0.1;
//vama1
//vama1      // temporary way of doing this for now
//vama1//      VextCosFunctor<double> Vext(world,new DipoleFunctor(2),omega);
//vama1      functionT fdipx = factoryT(world).functor(functorT(new DipoleFunctor(0))).initial_level(4);
//vama1      functionT fdipy = factoryT(world).functor(functorT(new DipoleFunctor(1))).initial_level(4);
//vama1      functionT fdipz = factoryT(world).functor(functorT(new DipoleFunctor(2))).initial_level(4);
//vama1
//vama1      world.gop.broadcast(time_step);
//vama1      world.gop.broadcast(nstep);
//vama1
//vama1      // Need complex orbitals :(
//vama1      double thresh = 1e-4;
//vama1      cvecfuncT camo = zero_functions<double_complex,3>(world, param.nalpha);
//vama1      cvecfuncT cbmo = zero_functions<double_complex,3>(world, param.nbeta);
//vama1      for (int iorb = 0; iorb < param.nalpha; iorb++)
//vama1      {
//vama1        camo[iorb] = std::exp(double_complex(0.0,2*constants::pi*strength))*amo[iorb];
//vama1        camo[iorb].truncate(thresh);
//vama1      }
//vama1      if (!param.spin_restricted && param.nbeta) {
//vama1        for (int iorb = 0; iorb < param.nbeta; iorb++)
//vama1        {
//vama1          cbmo[iorb] = std::exp(double_complex(0.0,2*constants::pi*strength))*bmo[iorb];
//vama1          cbmo[iorb].truncate(thresh);
//vama1        }
//vama1      }
//vama1
//vama1      // Create free particle propagator
//vama1      // Have no idea what to set "c" to
//vama1      double c = 20.0;
//vama1      printf("Creating G\n");
//vama1      Convolution1D<double_complex>* G = qm_1d_free_particle_propagator(FunctionDefaults<3>::get_k(), c, 0.5*time_step, 2.0*param.L);
//vama1      printf("Done creating G\n");
//vama1
//vama1      // Start iteration over time
//vama1      for (int step = 0; step < nstep; step++)
//vama1      {
//vama1//        if (world.rank() == 0) printf("Iterating step %d:\n\n", step);
//vama1        double t = time_step*step;
//vama1//        iterate_trotter(world, G, Vext, camo, cbmo, t, time_step);
//vama1        iterate_trotter(world, G, camo, cbmo, t, time_step, thresh);
//vama1        functionT arho = make_density(world,aocc,camo);
//vama1        functionT brho = (!param.spin_restricted && param.nbeta) ?
//vama1            make_density(world,aocc,camo) : copy(arho);
//vama1        functionT rho = arho + brho;
//vama1        double xval = inner(fdipx,rho);
//vama1        double yval = inner(fdipy,rho);
//vama1        double zval = inner(fdipz,rho);
//vama1        if (world.rank() == 0) printf("%15.7f%15.7f%15.7f%15.7f\n", t, xval, yval, zval);
//vama1      }
//vama1
//vama1
//vama1    }

    complex_functionT APPLY(const complex_operatorT* q1d, const complex_functionT& psi) {
        complex_functionT r = psi;  // Shallow copy violates constness !!!!!!!!!!!!!!!!!
        coordT lo, hi;
        lo[2] = -10;
        hi[2] = +10;

        r.reconstruct();
        r.broaden();
        r.broaden();
        r.broaden();
        r.broaden();
        r = apply_1d_realspace_push(*q1d, r, 2); r.sum_down();
        r = apply_1d_realspace_push(*q1d, r, 1); r.sum_down();
        r = apply_1d_realspace_push(*q1d, r, 0); r.sum_down();

        return r;
    }

    void iterate_trotter(World& world,
                         Convolution1D<double_complex>* G,
                         cvecfuncT& camo,
                         cvecfuncT& cbmo,
                         double t,
                         double time_step,
                         double thresh)
    {

      // first kinetic energy apply
      cvecfuncT camo2 = zero_functions<double_complex,3>(world, param.nalpha);
      cvecfuncT cbmo2 = zero_functions<double_complex,3>(world, param.nbeta);
      for (int iorb = 0; iorb < param.nalpha; iorb++)
      {
//        if (world.rank()) printf("Apply free-particle Green's function to alpha orbital %d\n", iorb);
        camo2[iorb] = APPLY(G,camo[iorb]);
        camo2[iorb].truncate(thresh);
      }
      if(!param.spin_restricted && param.nbeta)
      {
        for (int iorb = 0; iorb < param.nbeta; iorb++)
        {
          cbmo2[iorb] = APPLY(G,cbmo[iorb]);
          cbmo2[iorb].truncate(thresh);
        }
      }
      // Construct new density
//      START_TIMER(world);
      functionT arho = make_density(world, aocc, amo), brho;

      if (param.nbeta) {
          if (param.spin_restricted) {
              brho = arho;
          }
          else {
              brho = make_density(world, bocc, bmo);
          }
      }
      else {
          brho = functionT(world); // zero
      }
      functionT rho = arho + brho;
//      END_TIMER(world, "Make densities");

      // Do RPA only for now
      real_function_3d vnuc = potentialmanager->vnuclear();
      functionT vlocal = vnuc;
//      START_TIMER(world);
      functionT vcoul = apply(*coulop, rho);
//      END_TIMER(world, "Coulomb");
//      vlocal += vcoul + Vext(t+0.5*time_step);
//      vlocal += vcoul + std::cos(0.1*(t+0.5*time_step))*fdip;

      // exponentiate potential
//      if (world.rank()) printf("Apply Kohn-Sham potential to orbitals\n");
      complex_functionT expV = make_exp(time_step, vlocal);
      cvecfuncT camo3 = mul_sparse(world,expV,camo2,vtol,false);
      world.gop.fence();


      // second kinetic energy apply
      for (int iorb = 0; iorb < param.nalpha; iorb++)
      {
//        if (world.rank() == 0) printf("Apply free-particle Green's function to alpha orbital %d\n", iorb);
        camo3[iorb].truncate(thresh);
        camo[iorb] = APPLY(G,camo3[iorb]);
        camo[iorb].truncate();
      }
      if (!param.spin_restricted && param.nbeta)
      {
        cvecfuncT cbmo3 = mul_sparse(world,expV,cbmo2,vtol,false);

        // second kinetic energy apply
        for (int iorb = 0; iorb < param.nbeta; iorb++)
        {
          cbmo[iorb] = APPLY(G,cbmo3[iorb]);
          cbmo[iorb].truncate();
        }
      }
    }

    void solve(World & world)
    {
        functionT arho_old, brho_old;
        const double dconv = std::max(FunctionDefaults<3>::get_thresh(), param.dconv);
        const double trantol = vtol / std::min(30.0, double(amo.size()));
        const double tolloc = 1e-3;
        double update_residual = 0.0, bsh_residual = 0.0;
        subspaceT subspace;
        tensorT Q;
        bool do_this_iter = true;
        // Shrink subspace until stop localizing/canonicalizing
        int maxsub_save = param.maxsub;
        param.maxsub = 2;

        for(int iter = 0;iter < param.maxiter;++iter){
            if(world.rank() == 0)
                printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());

            if (iter > 0 && update_residual < 0.1) {
                //do_this_iter = false;
                param.maxsub = maxsub_save;
            }

            if(param.localize && do_this_iter) {
                tensorT U;
                if (param.localize_pm) {
                    U = localize_PM(world, amo, aset, tolloc, 0.25, iter == 0, true);
                }
                else {
                    U = localize_boys(world, amo, aset, tolloc, 0.25, iter==0);
                }
                amo = transform(world, amo, U, trantol, true);
                truncate(world, amo);
                normalize(world, amo);
                rotate_subspace(world, U, subspace, 0, amo.size(), trantol);
                if(!param.spin_restricted && param.nbeta != 0 ){
                    if (param.localize_pm) {
                        U = localize_PM(world, bmo, bset, tolloc, 0.25, iter == 0, true);
                    }
                    else {
                        U = localize_boys(world, bmo, bset, tolloc, 0.25, iter==0);
                    }
                    bmo = transform(world, bmo, U, trantol, true);
                    truncate(world, bmo);
                    normalize(world, bmo);
                    rotate_subspace(world, U, subspace, amo.size(), bmo.size(), trantol);
                }
            }

            START_TIMER(world);
            functionT arho = make_density(world, aocc, amo), brho;

            if (param.nbeta) {
                if (param.spin_restricted) {
                    brho = arho;
                }
                else {
                    brho = make_density(world, bocc, bmo);
                }
            }
            else {
                brho = functionT(world); // zero
            }
            END_TIMER(world, "Make densities");
	    print_meminfo(world.rank(), "Make densities");

            if(iter < 2 || (iter % 10) == 0){
                START_TIMER(world);
                loadbal(world, arho, brho, arho_old, brho_old, subspace);
                END_TIMER(world, "Load balancing");
		print_meminfo(world.rank(), "Load balancing");
            }
            double da = 0.0, db = 0.0;
            if(iter > 0){
                da = (arho - arho_old).norm2();
                db = (brho - brho_old).norm2();
                if(world.rank() == 0)
                    print("delta rho", da, db, "residuals", bsh_residual, update_residual);

            }

            arho_old = arho;
            brho_old = brho;
            functionT rho = arho + brho;
            rho.truncate();
            real_function_3d vnuc = potentialmanager->vnuclear();
            double enuclear = inner(rho, vnuc);


            START_TIMER(world);
            functionT vcoul = apply(*coulop, rho);
            functionT vlocal;
            END_TIMER(world, "Coulomb");
	    print_meminfo(world.rank(), "Coulomb");

            double ecoulomb = 0.5 * inner(rho, vcoul);
            rho.clear(false);
            vlocal = vcoul + vnuc ;
            
            vcoul.clear(false);
            vlocal.truncate();
            double exca = 0.0, excb = 0.0;

            vecfuncT vf, delrho;
            if (xc.is_dft()) {
                arho.reconstruct();
                if (param.nbeta != 0 && xc.is_spin_polarized()) brho.reconstruct();
                // brho.reconstruct();

                vf.push_back(arho);

                if (xc.is_spin_polarized()) vf.push_back(brho);

                if (xc.is_gga()) {

                    for (int axis=0; axis<3; ++axis) delrho.push_back((*gradop[axis])(arho,false)); // delrho
                    if (xc.is_spin_polarized())
                        for (int axis=0; axis<3; ++axis) delrho.push_back((*gradop[axis])(brho,false));


                    world.gop.fence(); // NECESSARY

                    vf.push_back(delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2]);     // sigma_aa

                    if (xc.is_spin_polarized())
                        vf.push_back(delrho[0]*delrho[3]+delrho[1]*delrho[4]+delrho[2]*delrho[5]); // sigma_ab
                    if (xc.is_spin_polarized())
                        vf.push_back(delrho[3]*delrho[3]+delrho[4]*delrho[4]+delrho[5]*delrho[5]); // sigma_bb

                    for (int axis=0; axis<3; ++axis) vf.push_back(delrho[axis]);        // dda_x

                    if (xc.is_spin_polarized())
                        for (int axis=0; axis<3; ++axis) vf.push_back(delrho[axis + 3]); // ddb_x
                    world.gop.fence(); // NECESSARY
                }
                if (vf.size()) {
                    reconstruct(world, vf);
//                    arho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
                    refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
                }
            }

            vecfuncT Vpsia = apply_potential(world, aocc, amo, vf, delrho, vlocal, exca, 0);

            vecfuncT Vpsib;
            if(!param.spin_restricted && param.nbeta) {
                Vpsib = apply_potential(world, bocc, bmo, vf, delrho, vlocal, excb, 1);
            }

            double ekina = 0.0, ekinb = 0.0;
            tensorT focka = make_fock_matrix(world, amo, Vpsia, aocc, ekina);
            tensorT fockb = focka;

            if (!param.spin_restricted && param.nbeta != 0)
                fockb = make_fock_matrix(world, bmo, Vpsib, bocc, ekinb);
            else if (param.nbeta != 0) {
                ekinb = ekina;
            }

            if (!param.localize && do_this_iter) {
                tensorT U = diag_fock_matrix(world, focka, amo, Vpsia, aeps, aocc, FunctionDefaults<3>::get_thresh());
                rotate_subspace(world, U, subspace, 0, amo.size(), trantol);
                if (!param.spin_restricted && param.nbeta != 0) {
                    U = diag_fock_matrix(world, fockb, bmo, Vpsib, beps, bocc, FunctionDefaults<3>::get_thresh());
                    rotate_subspace(world, U, subspace, amo.size(), bmo.size(), trantol);
                }
            }
            
            double enrep = molecule.nuclear_repulsion_energy();
            double ekinetic = ekina + ekinb;
            double exc = exca + excb;
            //double etot = ekinetic + enuclear + ecoulomb + exc + enrep + edipole;
            double etot = ekinetic + enuclear + ecoulomb + exc + enrep;
            current_energy = etot;
            esol = etot;

            if(world.rank() == 0){
                printf("\n              kinetic %16.8f\n", ekinetic);
                printf("   nuclear attraction %16.8f\n", enuclear);
                printf("              coulomb %16.8f\n", ecoulomb);
                printf(" exchange-correlation %16.8f\n", exc);
                printf("    nuclear-repulsion %16.8f\n", enrep);
                printf("                total %16.8f\n\n", etot);
            }

            if(iter > 0){
                //print("##convergence criteria: density delta=", da < dconv * molecule.natom() && db < dconv * molecule.natom(), ", bsh_residual=", (param.conv_only_dens || bsh_residual < 5.0*dconv));
                if(da < dconv * molecule.natom() && db < dconv * molecule.natom() && (param.conv_only_dens || bsh_residual < 5.0*dconv)){
                    if(world.rank() == 0) {
                        print("\nConverged!\n");
                    }

                    // Diagonalize to get the eigenvalues and if desired the final eigenvectors
                    tensorT U;
                    tensorT overlap = matrix_inner(world, amo, amo, true);

                    START_TIMER(world);
                    sygvp(world, focka, overlap, 1, U, aeps);
                    END_TIMER(world, "focka eigen sol");

                    if (!param.localize) {
                        amo = transform(world, amo, U, trantol, true);
                        truncate(world, amo);
                        normalize(world, amo);
                    }
                    if(param.nbeta != 0 && !param.spin_restricted){
                        overlap = matrix_inner(world, bmo, bmo, true);

                        START_TIMER(world);
                        sygvp(world, fockb, overlap, 1, U, beps);
                        END_TIMER(world, "fockb eigen sol");

                        if (!param.localize) {
                            bmo = transform(world, bmo, U, trantol, true);
                            truncate(world, bmo);
                            normalize(world, bmo);
                        }
                    }

                    if(world.rank() == 0) {
                        print(" ");
                        print("alpha eigenvalues");
                        print(aeps);
                        if(param.nbeta==0.0 && !param.spin_restricted){
                            print("beta eigenvalues");
                            print(beps);
                        }
                    }

                    if (param.localize) {
                        // Restore the diagonal elements for the analysis
                        for (unsigned int i=0; i<amo.size(); ++i) aeps[i] = focka(i,i);
                        for (unsigned int i=0; i<bmo.size(); ++i) beps[i] = fockb(i,i);
                    }

                    break;
                }

            }

            update_subspace(world, Vpsia, Vpsib, focka, fockb, subspace, Q, bsh_residual, update_residual);
        }

        if (world.rank() == 0) {
            if (param.localize) print("Orbitals are localized - energies are diagonal Fock matrix elements\n");
            else print("Orbitals are eigenvectors - energies are eigenvalues\n");
            print("Analysis of alpha MO vectors");
        }

        analyze_vectors(world, amo, aocc, aeps);
        if (param.nbeta != 0 && !param.spin_restricted) {
            if (world.rank() == 0)
                print("Analysis of beta MO vectors");

            analyze_vectors(world, bmo, bocc, beps);
        }


    }// end solve function

//dsol    void solve_finite_mo(World & world) {
//dsol
//dsol        double ep = param.epsf;
//dsol
//dsol        if(param.response_freq != 0.0) 
//dsol            error("Frequency should be ZERO.");
//dsol        if(world.rank() == 0) { 
//dsol            print("**********************");
//dsol            print("   start solve plus   ");
//dsol            print("**********************");
//dsol        }
//dsol        // amo_p, bmo_p
//dsol        solve_mp(world, amo_p, bmo_p, 0, ep);
//dsol
//dsol        if(world.rank() == 0) { 
//dsol            print("**********************");
//dsol            print("   start solve minus  ");
//dsol            print("**********************");
//dsol        }
//dsol        // amo_m, bmo_m
//dsol        solve_mp(world, amo_m, bmo_m, 1, ep);
//dsol    }
//dsol
//dsol    void solve_spolar(World & world)
//dsol    {
//dsol        double ep = param.epsf;
//dsol        subspaceT subspace;
//dsol        tensorT Q;
//dsol        //bool do_this_iter = true;
//dsol        // Shrink subspace until stop localizing/canonicalizing
//dsol
//dsol        
//dsol        // X:axis=0, Y:axis=1, Z:axis=2
//dsol        double omega = 0.0;
//dsol        int axis = param.response_axis;
//dsol
//dsol        response_frequency(world, axis);
//dsol        if(world.rank() == 0)  
//dsol            print(" Frequency for response function = ", omega);
//dsol        
//dsol        const double rconv = std::max(FunctionDefaults<3>::get_thresh(), param.rconv);
//dsol        vecfuncT ax = zero_functions<double, 3>(world, param.nalpha);
//dsol        vecfuncT ay = zero_functions<double, 3>(world, param.nalpha);
//dsol        vecfuncT bx = zero_functions<double, 3>(world, param.nbeta);
//dsol        vecfuncT by = zero_functions<double, 3>(world, param.nbeta);
//dsol
//dsol        // new response function
//dsol        vecfuncT ax_old = zero_functions<double,3>(world, param.nalpha);
//dsol        vecfuncT ay_old = zero_functions<double,3>(world, param.nalpha);
//dsol        vecfuncT bx_old = zero_functions<double,3>(world, param.nbeta);
//dsol        vecfuncT by_old = zero_functions<double,3>(world, param.nbeta);
//dsol        
//dsol        if(world.rank() == 0) { 
//dsol            print("\n\n\n");
//dsol            print(" ------------------------------------------------------------------------------");
//dsol            print(" |                MADNESS FINITE FIELD SOLVATION MODULE                  |");
//dsol            print(" ------------------------------------------------------------------------------");
//dsol            print(" \n\n");
//dsol        }
//dsol
//dsol        for(int iter = 0;iter < param.maxiter;++iter){
//dsol            if(world.rank() == 0)
//dsol                printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());
//dsol            double residual = 0.0;
//dsol
//dsol            for(int p=0; p<param.nalpha; ++p) {
//dsol                ax[p] = amo_p[p] - amo_m[p];
//dsol                ax[p].scale(1/(2*ep));
//dsol            }
//dsol            if(!param.spin_restricted && param.nbeta != 0) {
//dsol                for(int p=0; p<param.nbeta; ++p) {
//dsol                    bx[p] = bmo_p[p] - bmo_m[p];
//dsol                    bx[p].scale(1/(2*ep));
//dsol                }
//dsol            }
//dsol
//dsol            if(iter > 0) {
//dsol#if 0
//dsol                vecfuncT rax = zero_functions<double,3>(world, param.nalpha); //residual alpha x
//dsol                vecfuncT ray = zero_functions<double,3>(world, param.nalpha); //residual alpha y
//dsol                vecfuncT rbx = zero_functions<double,3>(world, param.nbeta);  //residual beta x
//dsol                vecfuncT rby = zero_functions<double,3>(world, param.nbeta);  //residual beta y
//dsol
//dsol                rax = sub(world, ax, ax_old);
//dsol                ray = sub(world, ay, ay_old);
//dsol
//dsol                if(!param.spin_restricted && param.nbeta != 0) {
//dsol                    // bxerr = (bx_new - bx).norm2() 
//dsol                    rbx = sub(world, bx, bx_old);
//dsol                    rby = sub(world, by, by_old);
//dsol                }
//dsol                double update_residual = 0.0;
//dsol                    update_response_subspace(world, ax, ay, bx, by, rax, ray, rbx, rby, subspace, Q, update_residual); 
//dsol#endif
//dsol                orthogonalize_response(ax, 0);
//dsol                orthogonalize_response(ay, 0);
//dsol                orthogonalize_response(bx, 1);
//dsol                orthogonalize_response(by, 1);
//dsol
//dsol                if(!param.spin_restricted && param.nbeta != 0) residual = norm2(world, sub(world, ax, ax_old)) + norm2(world, sub(world, bx, bx_old));
//dsol                else if (param.nbeta != 0) residual = 2 * norm2(world, sub(world, ax, ax_old));
//dsol                else residual = norm2(world, sub(world, ax, ax_old)); 
//dsol
//dsol                print("\nOLD response function (plus, alpha_spin) = ", norm2(world, ax_old));
//dsol                print("NEW response function (plus, alpha_spin) = ", norm2(world, ax));
//dsol                print("\nresiduals_response (final) = ", residual);
//dsol                print("rconv *(param.nalpha + param.nbeta)*2", rconv *(param.nalpha + param.nbeta)*2);
//dsol
//dsol                // 収束判定
//dsol                if( residual < (rconv *(param.nalpha + param.nbeta)*2))
//dsol                {
//dsol                    print("\n\n\n");
//dsol                    print(" ------------------------------------------------------------------------------");
//dsol                    print(" |                  MADNESS CALCULATION POLARIZABILITY                        |");
//dsol                    print(" ------------------------------------------------------------------------------");
//dsol                    print(" \n\n");
//dsol
//dsol                    break; 
//dsol                }
//dsol            }
//dsol            ax_old = ax;
//dsol            ay_old = ay;
//dsol            bx_old = bx;
//dsol            by_old = by;
//dsol
//dsol        } // end iteration
//dsol
//dsol        print("\nConverged response function!!\n");
//dsol
//dsol        double polar_alpha_p, polar_alpha_m, polar_beta_p, polar_beta_m;
//dsol        double total_polar;
//dsol        polar_alpha_p = polar_alpha_m = polar_beta_p = polar_beta_m = 0.0;
//dsol
//dsol        std::vector<int> f(3, 0);
//dsol        f[axis] = true;
//dsol        functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
//dsol
//dsol        vecfuncT damo = zero_functions<double,3>(world, param.nalpha); // dipolefunc * amo
//dsol        vecfuncT dbmo = zero_functions<double,3>(world, param.nbeta);  // dipolefunc * bmo
//dsol
//dsol        for(int i = 0; i<param.nalpha; ++i) {
//dsol            damo[i] = dipolefunc * amo[i];
//dsol
//dsol            polar_alpha_p += -2 * inner(ax[i], damo[i]);
//dsol            polar_alpha_m += -2 * inner(ay[i], damo[i]);
//dsol        }
//dsol        if(!param.spin_restricted && param.nbeta != 0) {
//dsol            for(int i = 0; i<param.nbeta; ++i) {
//dsol                dbmo[i] = dipolefunc * bmo[i];
//dsol                
//dsol                polar_beta_p += -2 * inner(bx[i], dbmo[i]);
//dsol                polar_beta_m += -2 * inner(by[i], dbmo[i]);
//dsol            }
//dsol        }
//dsol
//dsol        if(!param.spin_restricted && param.nbeta != 0) {
//dsol            total_polar = polar_alpha_p + polar_beta_p;
//dsol        }
//dsol        else if(param.nbeta != 0) total_polar = 2 * polar_alpha_p;
//dsol        else total_polar = polar_alpha_p;
//dsol
//dsol        if(!param.spin_restricted && param.nbeta != 0) {
//dsol            if(axis == 0) {
//dsol                print("Static Polarizability alpha [X][X]", polar_alpha_p);
//dsol                print("Static Polarizability beta  [X][X]", polar_beta_p);
//dsol                print("Static Polarizability TOTAL [X][X]", total_polar);
//dsol            }
//dsol            if(axis == 1) {
//dsol                print("Static Polarizability alpha [Y][Y]", polar_alpha_p);
//dsol                print("Static Polarizability beta  [Y][Y]", polar_beta_p);
//dsol                print("Static Polarizability TOTAL [Y][Y]", total_polar);
//dsol            }
//dsol            if(axis == 2) {
//dsol                print("Static Polarizability alpha [Z][Z]", polar_alpha_p);
//dsol                print("Static Polarizability beta  [Z][Z]", polar_beta_p);
//dsol                print("Static Polarizability TOTAL [Z][Z]", total_polar);
//dsol            }
//dsol        }
//dsol        else {
//dsol            if(axis == 0) {
//dsol                print("Static Polarizability alpha [X][X]", polar_alpha_p);
//dsol                print("Static Polarizability TOTAL [X][X]", total_polar);
//dsol            }
//dsol            if(axis == 1) {
//dsol                print("Static Polarizability alpha [Y][Y]", polar_alpha_p);
//dsol                print("Static Polarizability TOTAL [Y][Y]", total_polar);
//dsol            }
//dsol            if(axis == 2) { 
//dsol                print("Static Polarizability alpha [Z][Z]", polar_alpha_p);
//dsol                print("Static Polarizability TOTAL [Z][Z]", total_polar);
//dsol            }
//dsol        }
//dsol
//dsol    } // end solve_polar function

//dsolmp    void solve_mp(World & world, vecfuncT & amo_mp, vecfuncT & bmo_mp, int mflag, double & ep) 
//dsolmp    {
//dsolmp
//dsolmp        functionT arho_old, brho_old;
//dsolmp        const double dconv = std::max(FunctionDefaults<3>::get_thresh(), param.dconv);
//dsolmp        const double trantol = vtol / std::min(30.0, double(amo_mp.size()));
//dsolmp        const double tolloc = 1e-3;
//dsolmp        double update_residual = 0.0, bsh_residual = 0.0;
//dsolmp        subspaceT subspace;
//dsolmp        tensorT Q;
//dsolmp        bool do_this_iter = true;
//dsolmp        int maxsub_save = param.maxsub;
//dsolmp        param.maxsub = 2;
//dsolmp
//dsolmp        for(int iter = 0;iter < param.maxiter;++iter){
//dsolmp            if(world.rank() == 0)
//dsolmp                printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());
//dsolmp
//dsolmp            if (iter > 0 && update_residual < 0.1) {
//dsolmp                //do_this_iter = false;
//dsolmp                param.maxsub = maxsub_save;
//dsolmp            }
//dsolmp            if(param.localize && do_this_iter) {
//dsolmp                tensorT U;
//dsolmp                if (param.localize_pm) {
//dsolmp                    U = localize_PM(world, amo_mp, aset_mp, tolloc, 0.25, iter == 0, true);
//dsolmp                }
//dsolmp                else {
//dsolmp                    U = localize_boys(world, amo_mp, aset_mp, tolloc, 0.25, iter==0);
//dsolmp                }
//dsolmp
//dsolmp                amo_mp = transform(world, amo_mp, U, trantol, true);
//dsolmp                truncate(world, amo_mp);
//dsolmp                normalize(world, amo_mp);
//dsolmp                rotate_subspace(world, U, subspace, 0, amo_mp.size(), trantol);
//dsolmp                if(!param.spin_restricted && param.nbeta != 0 ){
//dsolmp                    if (param.localize_pm) {
//dsolmp                        U = localize_PM(world, bmo_mp, bset_mp, tolloc, 0.25, iter == 0, true);
//dsolmp                    }
//dsolmp                    else {
//dsolmp                        U = localize_boys(world, bmo_mp, bset_mp, tolloc, 0.25, iter==0);
//dsolmp                    }
//dsolmp                    bmo_mp = transform(world, bmo_mp, U, trantol, true);
//dsolmp                    truncate(world, bmo_mp);
//dsolmp                    normalize(world, bmo_mp);
//dsolmp                    rotate_subspace(world, U, subspace, amo_mp.size(), bmo_mp.size(), trantol);
//dsolmp                }
//dsolmp            }
//dsolmp
//dsolmp            START_TIMER(world);
//dsolmp            functionT arho, brho;
//dsolmp            if(iter == 0) arho = make_density(world, aocc, amo);
//dsolmp            else arho = make_density(world, aocc_mp, amo_mp);
//dsolmp
//dsolmp            if (param.nbeta) {
//dsolmp                if (param.spin_restricted) {
//dsolmp                    brho = arho;
//dsolmp                }
//dsolmp                else {
//dsolmp                    if(iter == 0) brho = make_density(world, bocc, bmo);
//dsolmp                    else brho = make_density(world, bocc_mp, bmo_mp);
//dsolmp                }
//dsolmp            }
//dsolmp            else {
//dsolmp                brho = functionT(world); // zero
//dsolmp            }
//dsolmp            END_TIMER(world, "Make densities");
//dsolmp            print_meminfo(world.rank(), "Make densities");
//dsolmp
//dsolmp            if(iter < 2 || (iter % 10) == 0){
//dsolmp                START_TIMER(world);
//dsolmp                loadbal(world, arho, brho, arho_old, brho_old, subspace);
//dsolmp                END_TIMER(world, "Load balancing");
//dsolmp                print_meminfo(world.rank(), "Load balancing");
//dsolmp            }
//dsolmp            double da = 0.0, db = 0.0;
//dsolmp            if(iter > 0){
//dsolmp                da = (arho - arho_old).norm2();
//dsolmp                db = (brho - brho_old).norm2();
//dsolmp                if(world.rank() == 0)
//dsolmp                    print("delta rho", da, db, "residuals", bsh_residual, update_residual);
//dsolmp
//dsolmp            }
//dsolmp
//dsolmp            arho_old = arho;
//dsolmp            brho_old = brho;
//dsolmp            functionT rho = arho + brho;
//dsolmp            rho.truncate();
//dsolmp            print("rho = ", rho.norm2());
//dsolmp            real_function_3d vnuc = potentialmanager->vnuclear();
//dsolmp            double enuclear = inner(rho, vnuc);
//dsolmp
//dsolmp            START_TIMER(world);
//dsolmp            functionT vcoul = apply(*coulop, rho);
//dsolmp            functionT vlocal;
//dsolmp            END_TIMER(world, "Coulomb");
//dsolmp            print_meminfo(world.rank(), "Coulomb");
//dsolmp
//dsolmp            double ecoulomb = 0.5 * inner(rho, vcoul);
//dsolmp            //rho.clear(false);
//dsolmp            vlocal = vcoul + vnuc ;
//dsolmp
//dsolmp            print("Vloc1", vlocal.trace(), vlocal.norm2());
//dsolmp
//dsolmp            std::vector<int> f(3, 0);
//dsolmp            f[param.response_axis] = true; // 2 -> axis = z
//dsolmp            functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
//dsolmp            dipolefunc.scale(ep);
//dsolmp            dipolefunc.truncate();
//dsolmp            
//dsolmp            if (mflag == 1) vlocal = vlocal - dipolefunc;
//dsolmp            else vlocal = vlocal + dipolefunc;
//dsolmp            print("Vloc2", vlocal.trace(), vlocal.norm2());
//dsolmp            
//dsolmp            double edipole = inner(rho, dipolefunc);
//dsolmp            print("Energy of finite field = ", edipole);
//dsolmp
//dsolmp            rho.clear(false);
//dsolmp
//dsolmp            vcoul.clear(false);
//dsolmp            vlocal.truncate();
//dsolmp            double exca = 0.0, excb = 0.0;
//dsolmp
//dsolmp            vecfuncT vf, delrho;
//dsolmp            if (xc.is_dft()) {
//dsolmp                arho.reconstruct();
//dsolmp                if (param.nbeta != 0 && xc.is_spin_polarized()) brho.reconstruct();
//dsolmp                // brho.reconstruct();
//dsolmp
//dsolmp                vf.push_back(arho);
//dsolmp
//dsolmp                if (xc.is_spin_polarized()) vf.push_back(brho);
//dsolmp
//dsolmp                if (xc.is_gga()) {
//dsolmp
//dsolmp                    for (int axis=0; axis<3; ++axis) delrho.push_back((*gradop[axis])(arho,false)); // delrho
//dsolmp                    if (xc.is_spin_polarized())
//dsolmp                        for (int axis=0; axis<3; ++axis) delrho.push_back((*gradop[axis])(brho,false));
//dsolmp
//dsolmp
//dsolmp                    world.gop.fence(); // NECESSARY
//dsolmp
//dsolmp                    vf.push_back(delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2]);     // sigma_aa
//dsolmp
//dsolmp                    if (xc.is_spin_polarized())
//dsolmp                        vf.push_back(delrho[0]*delrho[3]+delrho[1]*delrho[4]+delrho[2]*delrho[5]); // sigma_ab
//dsolmp                    if (xc.is_spin_polarized())
//dsolmp                        vf.push_back(delrho[3]*delrho[3]+delrho[4]*delrho[4]+delrho[5]*delrho[5]); // sigma_bb
//dsolmp
//dsolmp                    for (int axis=0; axis<3; ++axis) vf.push_back(delrho[axis]);        // dda_x
//dsolmp
//dsolmp                    if (xc.is_spin_polarized())
//dsolmp                        for (int axis=0; axis<3; ++axis) vf.push_back(delrho[axis + 3]); // ddb_x
//dsolmp                    world.gop.fence(); // NECESSARY
//dsolmp                }
//dsolmp                if (vf.size()) {
//dsolmp                    reconstruct(world, vf);
//dsolmp                    arho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
//dsolmp                }
//dsolmp            }
//dsolmp
//dsolmp            vecfuncT Vpsia = apply_potential(world, aocc_mp, amo_mp, vf, delrho, vlocal, exca, 0);
//dsolmp
//dsolmp            vecfuncT Vpsib;
//dsolmp            if(!param.spin_restricted && param.nbeta) {
//dsolmp                Vpsib = apply_potential(world, bocc_mp, bmo_mp, vf, delrho, vlocal, excb, 1);
//dsolmp            }
//dsolmp
//dsolmp            double ekina = 0.0, ekinb = 0.0;
//dsolmp            tensorT focka = make_fock_matrix(world, amo_mp, Vpsia, aocc_mp, ekina);
//dsolmp            tensorT fockb = focka;
//dsolmp
//dsolmp            if (!param.spin_restricted && param.nbeta != 0)
//dsolmp                fockb = make_fock_matrix(world, bmo_mp, Vpsib, bocc_mp, ekinb);
//dsolmp            else if (param.nbeta != 0) {
//dsolmp                ekinb = ekina;
//dsolmp            }
//dsolmp
//dsolmp            if (!param.localize && do_this_iter) {
//dsolmp                tensorT U = diag_fock_matrix(world, focka, amo_mp, Vpsia, aeps_mp, aocc_mp, FunctionDefaults<3>::get_thresh());
//dsolmp                rotate_subspace(world, U, subspace, 0, amo_mp.size(), trantol);
//dsolmp                if (!param.spin_restricted && param.nbeta != 0) {
//dsolmp                    U = diag_fock_matrix(world, fockb, bmo_mp, Vpsib, beps_mp, bocc_mp, FunctionDefaults<3>::get_thresh());
//dsolmp                    rotate_subspace(world, U, subspace, amo_mp.size(), bmo_mp.size(), trantol);
//dsolmp                }
//dsolmp            }
//dsolmp
//dsolmp            double enrep = molecule.nuclear_repulsion_energy();
//dsolmp            double ekinetic = ekina + ekinb;
//dsolmp            double exc = exca + excb;
//dsolmp            double etot = ekinetic + enuclear + ecoulomb + exc + enrep + edipole;
//dsolmp            current_energy = etot;
//dsolmp            esol = etot;
//dsolmp
//dsolmp            if(world.rank() == 0){
//dsolmp                printf("\n              kinetic %16.8f\n", ekinetic);
//dsolmp                printf("   nuclear attraction %16.8f\n", enuclear);
//dsolmp                printf("              coulomb %16.8f\n", ecoulomb);
//dsolmp                printf(" exchange-correlation %16.8f\n", exc);
//dsolmp                printf("    nuclear-repulsion %16.8f\n", enrep);
//dsolmp                printf("                total %16.8f\n\n", etot);
//dsolmp            }
//dsolmp
//dsolmp            if(iter > 0){
//dsolmp                //print("##convergence criteria: density delta=", da < dconv * molecule.natom() && db < dconv * molecule.natom(), ", bsh_residual=", (param.conv_only_dens || bsh_residual < 5.0*dconv));
//dsolmp                if(da < dconv * molecule.natom() && db < dconv * molecule.natom() && (param.conv_only_dens || bsh_residual < 5.0*dconv)){
//dsolmp                    if(world.rank() == 0) {
//dsolmp                        print("\nConverged!\n");
//dsolmp                    }
//dsolmp
//dsolmp                    // Diagonalize to get the eigenvalues and if desired the final eigenvectors
//dsolmp                    tensorT U;
//dsolmp                    tensorT overlap = matrix_inner(world, amo_mp, amo_mp, true);
//dsolmp
//dsolmp                    START_TIMER(world);
//dsolmp                    sygvp(world, focka, overlap, 1, U, aeps_mp);
//dsolmp                    END_TIMER(world, "focka eigen sol");
//dsolmp
//dsolmp                    if (!param.localize) {
//dsolmp                        amo_mp = transform(world, amo_mp, U, trantol, true);
//dsolmp                        truncate(world, amo_mp);
//dsolmp                        normalize(world, amo_mp);
//dsolmp                    }
//dsolmp                    if(param.nbeta != 0 && !param.spin_restricted){
//dsolmp                        overlap = matrix_inner(world, bmo_mp, bmo_mp, true);
//dsolmp
//dsolmp                        START_TIMER(world);
//dsolmp                        sygvp(world, fockb, overlap, 1, U, beps_mp);
//dsolmp                        END_TIMER(world, "fockb eigen sol");
//dsolmp
//dsolmp                        if (!param.localize) {
//dsolmp                            bmo_mp = transform(world, bmo_mp, U, trantol, true);
//dsolmp                            truncate(world, bmo_mp);
//dsolmp                            normalize(world, bmo_mp);
//dsolmp                        }
//dsolmp                    }
//dsolmp
//dsolmp                    if(world.rank() == 0) {
//dsolmp                        print(" ");
//dsolmp                        print("alpha mp eigenvalues");
//dsolmp                        print(aeps_mp);
//dsolmp                        if(param.nbeta==0.0 && !param.spin_restricted){
//dsolmp                            print("beta mp eigenvalues");
//dsolmp                            print(beps_mp);
//dsolmp                        }
//dsolmp                    }
//dsolmp
//dsolmp                    if (param.localize) {
//dsolmp                        // Restore the diagonal elements for the analysis
//dsolmp                        for (unsigned int i=0; i<amo_mp.size(); ++i) aeps_mp[i] = focka(i,i);
//dsolmp                        for (unsigned int i=0; i<bmo_mp.size(); ++i) beps_mp[i] = fockb(i,i);
//dsolmp                    }
//dsolmp
//dsolmp                    break;
//dsolmp                }
//dsolmp            }
//dsolmp            update_subspace_mp(world, Vpsia, Vpsib, focka, fockb, subspace, Q, bsh_residual, update_residual, mflag, amo_mp, bmo_mp);
//dsolmp        }
//dsolmp        if (world.rank() == 0) {
//dsolmp            if (param.localize) print("Orbitals are localized - energies are diagonal Fock matrix elements\n");
//dsolmp            else print("Orbitals are eigenvectors - energies are eigenvalues\n");
//dsolmp            print("Analysis of alpha MO vectors");
//dsolmp        }
//dsolmp
//dsolmp        analyze_vectors(world, amo_mp, aocc_mp, aeps_mp);
//dsolmp        if (param.nbeta != 0 && !param.spin_restricted) {
//dsolmp            if (world.rank() == 0)
//dsolmp                print("Analysis of beta MO vectors");
//dsolmp
//dsolmp            analyze_vectors(world, bmo_mp, bocc_mp, beps_mp);
//dsolmp        }
//dsolmp
//dsolmp    }// end solve mp function

    void response_frequency(World & world, int & axis)
    {
        if (world.rank() == 0) { 
            print("\n");
            if(axis == 0) 
                print(" AXIS of frequency = x");

            else if(axis == 1) 
                print(" AXIS of frequency = y");

            else if(axis == 2) 
                print(" AXIS of frequency = z");
        }
    }

    vecfuncT calc_dipole_mo(World & world,  vecfuncT & mo, int & axis, int  & size)
    {
        START_TIMER(world);

        vecfuncT dipolemo = zero_functions<double,3>(world, size);
        
        std::vector<int> f(3, 0);
        f[axis] = 1;
        //print("f = ", f[0]," ",  f[1], " ", f[2]);
        functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));

        // dipolefunc * mo[iter]
        for(int p=0; p<size; ++p)
            dipolemo[p] = dipolefunc * mo[p];

        END_TIMER(world, "Make perturbation");
        print_meminfo(world.rank(), "Make perturbation");

        return dipolemo;
    }

    void calc_freq(World & world, double & omega, tensorT & ak, tensorT & bk, int sign)
    {
        if (world.rank() == 0) { 
            print(" eps_alpha");
            print(aeps);
            if(!param.spin_restricted && param.nbeta != 0) print(" eps_beta  = ", beps);

            print(" frequency = ", omega);
            print(" number of alpha orbital = ", param.nalpha);
            print(" number of beta orbital = ", param.nbeta);
        }

        for(int i=0; i<param.nalpha; ++i){
            ak[i] = sqrt(-2.0 * (aeps[i] + sign * omega));
            if (world.rank() == 0)  
                print(" kxy(alpha) [", i, "] : sqrt(-2 * (eps +/- omega)) = ", ak[i]);
        }
        if(!param.spin_restricted && param.nbeta != 0) {
            for(int i=0; i<param.nbeta; ++i){
                bk[i] = sqrt(-2.0 * (beps[i] + sign * omega));
                    if (world.rank() == 0)  
                        print(" kxy(beta) [", i, "]: sqrt(-2 * (eps +/- omega)) = ", bk[i]);
    }
    }
    }

    void make_BSHOperatorPtr(World & world, tensorT & ak, tensorT & bk,
            std::vector<poperatorT> & aop, std::vector<poperatorT> & bop)
    {
        START_TIMER(world);
        double tol = FunctionDefaults < 3 > ::get_thresh();

        for(int i=0; i<param.nalpha; ++i) {
            // thresh tol : 1e-6
            aop[i] = poperatorT(BSHOperatorPtr3D(world, ak[i], param.lo , tol));
        }
        if(!param.spin_restricted && param.nbeta != 0) {
            for(int i=0; i<param.nbeta; ++i) {
                bop[i] = poperatorT(BSHOperatorPtr3D(world, bk[i], param.lo, tol));
            }
        } 

        END_TIMER(world, "Make BSHOp");
        print_meminfo(world.rank(), "Make BSHOp");
    }

    vecfuncT initial_guess_response(World & world, vecfuncT & dipolemo,
            functionT & vlocal, vecfuncT & vf, int  spin, int & axis) {
        
        START_TIMER(world);

        vecfuncT rhs, Vmo, djkmo;
        double exca = 0.0, excb = 0.0;
        
        if (world.rank() == 0)  
            print("\nmake initial guess for response function");
        
      functionT arho = make_density(world, aocc, amo);
      functionT brho;

      if (param.nbeta) {
          if (!param.spin_restricted) {
              brho = make_density(world, bocc, bmo);
          }
          else {
              brho = arho;
          }
      }
      else {
          brho = functionT(world); // zero
      }
      functionT rho = arho + brho;

        if(spin == 0){
          // initial guess >>  x,y_MOs = dipoleamo
            Vmo = apply_potential_response(world, aocc, dipolemo, vf, vlocal, exca, 0 );
            djkmo = calc_djkmo(world, dipolemo, dipolemo, rho, 0);
            rhs = calc_rhs(world, amo, Vmo, dipolemo, djkmo);
        }
        else {
            Vmo = apply_potential_response(world, bocc, dipolemo, vf, vlocal, excb, 1);
            djkmo = calc_djkmo(world, dipolemo, dipolemo, rho, 1);
            rhs = calc_rhs(world, bmo, Vmo, dipolemo, djkmo);
        }

        //rhs = add(world, Vmo, add(world, dipolemo, djkmo));

        END_TIMER(world, "Initial for response");
        print_meminfo(world.rank(), "Initial for response");
        
        //truncate(world, rhs);
        return rhs;
    }

    functionT make_density_ground(World & world, functionT & arho, functionT & brho)
    {
        functionT rho = factoryT(world);

        START_TIMER(world);            

        arho = make_density(world, aocc, amo);
        if (!param.spin_restricted) {
            brho = make_density(world, bocc, bmo);
        }
        else {
            brho = arho;
        }

        rho = arho + brho;
        rho.truncate();

        END_TIMER(world, "Make densities");
        print_meminfo(world.rank(), "Make densities");

        return rho;
    }

    functionT make_derivative_density(World & world, vecfuncT & x, vecfuncT & y,
            int  flag, int & size)
    {
        functionT drho = factoryT(world);

        if(flag == 0){
            for(int i=0; i<size; ++i) {
                drho += (amo[i] * x[i]) + (amo[i] * y[i]);
            }
            //drho.scale(2.0);
        }
        else{
            for(int i=0; i<size; ++i) {
                drho += (bmo[i] * x[i]) + (bmo[i] * y[i]);
            }
        }
        return drho;
    }

    functionT calc_derivative_Jmo(World & world, int & p, 
            functionT & drho, 
            vecfuncT & mo, int  spin, int & size)
    {
        functionT dJ = factoryT(world);
        functionT dJmo = factoryT(world);
        functionT j = factoryT(world);
        
        //START_TIMER(world);

        dJ = apply(*coulop,   drho  );
//vam2        for(int i=0; i<size; ++i) {
//vam2            j = apply(*coulop,  ( dmo1[i] * mo[i] ) + ( dmo2[i] * mo[i] ) );
//vam2            dJ = dJ + j;
//vam2            j.clear(false);
//vam2        }
        dJmo = dJ * mo[p];
        
        //END_TIMER(world, "Calc derivative J");
        return dJmo;
    }

    functionT calc_exchange_function(World & world, int & p,
            vecfuncT & dmo1, vecfuncT & dmo2,
            vecfuncT & mo, int & size)
    {
        functionT dKmo = factoryT(world);
        functionT k1 = factoryT(world);
        functionT k2 = factoryT(world);

        //START_TIMER(world);
        for(int i=0; i<size; ++i) {
            k1 = apply(*coulop, ( mo[i] * mo[p] )) * dmo1[i];
            k2 = apply(*coulop, ( mo[p] * dmo2[i] )) * mo[i];
            dKmo = dKmo + (k1 + k2);
            k1.clear(false);
            k2.clear(false);
        }
        //END_TIMER(world, "Calc derivative xc");
        return dKmo;
    }
    
    vecfuncT calc_djkmo(World & world, vecfuncT & dmo1, 
            vecfuncT & dmo2,  functionT drho, int  spin)
    {
        START_TIMER(world);

        int size = 0;
        vecfuncT mo; 

        if(spin == 0) {
           size = param.nalpha;
           mo = amo;
        }
        else {
           size = param.nbeta;
           mo = bmo;
        }

        vecfuncT djkmo = zero_functions<double,3>(world, size);
        
        for(int p=0; p<size; ++p) {
            functionT dJmo = calc_derivative_Jmo(world, p, drho, mo,  spin, size);        
            functionT dKmo = calc_exchange_function(world, p, dmo1, dmo2, mo, size);

            djkmo[p] = dJmo - dKmo;

            dJmo.clear(false);
            dKmo.clear(false);
        }

        END_TIMER(world, "Calc derivative V");
        print_meminfo(world.rank(), "Calc derivative V");
        return djkmo;
    }

//prod    vecfuncT calc_rhs(World & world, vecfuncT & Vdmo, vecfuncT & dipolemo, vecfuncT & djkmo) 
//prod    {
//prod        vecfuncT rhs;
//prod
//prod        START_TIMER(world);
//prod        //dmo_rhs = Vdmo + dipolemo + djkmo;
//prod        rhs = add(world, Vdmo, add(world, dipolemo, djkmo));
//prod        END_TIMER(world, "Sum rhs response");
//prod        print_meminfo(world.rank(), "Sum rhs response");
//prod        truncate(world, rhs);
//prod
//prod        return rhs;
//prod    }

    vecfuncT calc_rhs(World & world, vecfuncT & mo , vecfuncT & Vdmo, vecfuncT & dipolemo, vecfuncT & djkmo )
    {
        vecfuncT rhs = zero_functions<double,3>(world, Vdmo.size());

        START_TIMER(world);
        //dmo_rhs = Vdmo + dipolemo + djkmo;

        // the projector on the unperturbed density
        Projector<double,3> rho0(mo);

       vecfuncT gp = add(world, dipolemo, djkmo);
//        vecfuncT gp1= mul_sparse(world, rho0, gp, vtol);
//prod        vecfuncT Gampsi;
        for (int i=0; i<Vdmo.size(); ++i) {
            functionT gp1 =  gp[i];
            gp1 = gp1 - rho0(gp1);
//prod            functionT a = rho0*gp;
//prod            Gampsi.push_back(gp);
            gp1 = Vdmo[i] + gp1 ;
            rhs.push_back(gp1);
        }
//prod
//prod        //rhs = add(world, Vdmo, add(world, dipolemo, djkmo));
//prod        rhs = add(world, Vdmo, Gampsi);
//prod        rhs = add(world, Vdmo, add(world, dipolemo, djkmo));

        END_TIMER(world, "Sum rhs response");
        print_meminfo(world.rank(), "Sum rhs response");
        truncate(world, rhs);

        return rhs;
    }


    void calc_response_function(World & world, vecfuncT & dmo, 
            std::vector<poperatorT> & op, vecfuncT & rhs)
    {
        // new response function
        // BSHOperatorPrt3D : op
        dmo = apply(world, op, rhs);
        scale(world, dmo, -2.0);
        truncate(world, dmo);
    }

    // orthogonalization
    void orthogonalize_response(vecfuncT & dmo, int flag)
    {
        if(flag == 0){ 
            for(int i=0; i<param.nalpha; ++i){
                for (int j=0; j<param.nalpha; ++j){
                    // new_x = new_x - < psi | new_x > * psi
                    dmo[i] = dmo[i] - dmo[i].inner(amo[j])*amo[j];
                }
            }
        }
        else{
            if(!param.spin_restricted && param.nbeta != 0) {
                for(int i=0; i<param.nbeta; ++i){
                    for (int j=0; j<param.nbeta; ++j){
                        // new_x = new_x - < psi | new_x > * psi
                        dmo[i] = dmo[i] - dmo[i].inner(bmo[j])*bmo[j];
                    }
                }
            }
        }
    }


//vama ugly ! alpha_ij(w) = - sum(m occ) [<psi_m(0)|r_i|psi_mj(1)(w)> + <psi_mj(1)(-w)|r_i|psi_m(0)>]

    void dpolar(World & world, tensorT & polar, functionT & drho, int & axis)
    {
        functionT dipolefunc;

        for(int i=0; i<3; ++i) {
            std::vector<int> f(3, 0);
            f[i] = 1;
            dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(f)));
            polar(axis, i) = -2 * dipolefunc.inner(drho);

            dipolefunc.clear(false);
        }
    }

    void calc_dpolar(World & world, double & Dpolar_average, 
            vecfuncT & ax, vecfuncT & ay, 
            vecfuncT & bx, vecfuncT & by, 
            int & axis, tensorT & Dpolar_total,
            tensorT & Dpolar_alpha, tensorT & Dpolar_beta)
    {
        int fflag = 0;

        START_TIMER(world);
        // derivative density matrix
        functionT drhoa = make_derivative_density(world, ax, ay, 0, param.nalpha);
        functionT drhob;
        if(!param.spin_restricted) 
            drhob = make_derivative_density(world, bx, by ,1, param.nbeta);
        else 
            drhob = drhoa;

        functionT drho = drhoa + drhob;
        
        dpolar(world, Dpolar_alpha, drhoa, axis);

        dpolar(world, Dpolar_beta, drhob, axis);

        dpolar(world, Dpolar_total, drho, axis);
        
        for(int i=0; i<3; ++i) {
            Dpolar_total(axis, i) = 0.5 * Dpolar_total(axis, i);

            if (world.rank() == 0) { 
                print("Dynamic Polarizability alpha [", axis, "][", i, "]", Dpolar_alpha(axis, i));
                if(param.nbeta != 0) 
                    print("Dynamic Polarizability beta  [", axis, "][", i, "]", Dpolar_beta(axis, i));
                print("Dynamic Polarizability TOTAL [", axis, "][", i, "]", Dpolar_total(axis, i), "\n");
            }
            
            // final flag
            if((i == 2) && (axis == 2)) fflag = 1;
        }

        drhoa.clear(false);
        drhob.clear(false);

        //diagonalize
        tensorT V, epolar, eapolar, ebpolar;
        if(fflag != 0) {
            syev(Dpolar_alpha, V, eapolar);
            syev(Dpolar_total, V, epolar);
                if(param.nbeta != 0) 
                    syev(Dpolar_beta, V, ebpolar);
            for(unsigned int i=0; i<3; ++i) {
                Dpolar_average = Dpolar_average + epolar[i];
            }
        }

        if (world.rank() == 0) { 
            if(fflag != 0) {
                for(int j=0; j<3; ++j) {
                    for(int k=0; k<3; ++k) {
                        if(Dpolar_alpha(j, k) < 0.001 && Dpolar_alpha(j,k) < -0.001)
                            Dpolar_alpha(j, k) = 0.0;
                        if(Dpolar_beta(j, k) < 0.001 && Dpolar_beta(j, k) < -0.001)
                            Dpolar_beta(j, k) = 0.0;
                        if(Dpolar_total(j, k) < 0.001 && Dpolar_total(j, k) < -0.001)
                            Dpolar_total(j, k) = 0.0;
                    }
                }
                    print("\n");
                print("************************************************************************");
                print("\tEigenvalues of Dynamic Polarizability ");
                print("************************************************************************");
                print("\n\t", epolar[0], "\t", epolar[1], "\t", epolar[2]);
                print("\n************************************************************************");
                
                print("\n\n");
                print("************************************************************************");
                print("\tDynamic Polarizability alpha ( Frequency = ", param.response_freq, ")");
                print("************************************************************************");
                print("\n\t", Dpolar_alpha(0,0), "\t", Dpolar_alpha(0,1), "\t", Dpolar_alpha(0,2));
                print("\t", Dpolar_alpha(1,0), "\t", Dpolar_alpha(1,1), "\t", Dpolar_alpha(1,2));
                print("\t", Dpolar_alpha(2,0), "\t", Dpolar_alpha(2,1), "\t", Dpolar_alpha(2,2));
                print("\n************************************************************************");

                
                if(param.nbeta != 0) {
                    print("\n");
                    print("************************************************************************");
                    print("\tDynamic Polarizability beta ( Frequency = ", param.response_freq, ")");
                    print("************************************************************************");
                    print("\n\t", Dpolar_beta(0,0), "\t", Dpolar_beta(0,1), "\t", Dpolar_beta(0,2));
                    print("\t", Dpolar_beta(1,0), "\t", Dpolar_beta(1,1), "\t", Dpolar_beta(1,2));
                    print("\t", Dpolar_beta(2,0), "\t", Dpolar_beta(2,1), "\t", Dpolar_beta(2,2));
                    print("\n************************************************************************");
               }
                print("\n");
                print("\tTotal Dynamic Polarizability ( Frequency = ", param.response_freq, ")");
                print("************************************************************************");
                print("\n\t", Dpolar_total(0,0), "\t", Dpolar_total(0,1), "\t", Dpolar_total(0,2));
                print("\t", Dpolar_total(1,0), "\t", Dpolar_total(1,1), "\t", Dpolar_total(1,2));
                print("\t", Dpolar_total(2,0), "\t", Dpolar_total(2,1), "\t", Dpolar_total(2,2));
                print("\n************************************************************************");
                
                print("\nThe average of polarizability = ", Dpolar_average/3.0, "\n");
            }
        }
        END_TIMER(world, "Calc D polar");
        print_meminfo(world.rank(), "Calc D polar");

    // end of solving polarizability
    }

    double residual_response(World & world, vecfuncT & x, vecfuncT & y,
            vecfuncT & x_old, vecfuncT & y_old,
            vecfuncT & rx, vecfuncT & ry)
    {
        double residual = 0.0;
        
        START_TIMER(world);
        rx = sub(world, x_old, x);
        ry = sub(world, y_old, y);
        std::vector<double> rnormx = norm2s(world, rx);
        std::vector<double> rnormy = norm2s(world, ry);

        double rms, maxval_x, maxval_y;
        vector_stats(rnormx, rms, maxval_x);
        vector_stats(rnormy, rms, maxval_y);
        residual = std::max(maxval_x, maxval_y);
        
        END_TIMER(world, "Residual X,Y");
        print_meminfo(world.rank(), "Residual X,Y");

        return residual;
    }

    // solve response function
    void solve_response(World & world)
    {
        if(world.rank() == 0) {
            print("\n\n\n");
            print(" ------------------------------------------------------------------------------");
            print(" |                MADNESS RESPONSE                                            |");
            print(" ------------------------------------------------------------------------------");
            print(" \n\n");
        }

        const double rconv = std::max(FunctionDefaults<3>::get_thresh(), param.rconv);

        subspaceT subspace;
        tensorT Q;

        double update_residual = 0.0;
        int maxsub_save = param.maxsub;

                    if(world.rank() == 0) {
                        print(" ");
                        print("alpha eigenvalues 2");
                        print(aeps);
                    }
        
        START_TIMER(world);
        // Green's function
        tensorT akx(param.nalpha);
        tensorT aky(param.nalpha);
        tensorT bkx(param.nbeta);
        tensorT bky(param.nbeta);
        
        // X:axis=0, Y:axis=1, Z:axis=2
        double omega = param.response_freq;
        
        START_TIMER(world);
        //calculate frequency term
        calc_freq(world, omega, akx, bkx, 1);

        if(omega != 0.0)
            calc_freq(world, omega, aky, bky, -1);
        END_TIMER(world, "Make frequency term");
        print_meminfo(world.rank(), "Make frequency term");

        // make density matrix
        functionT arho ; // = factoryT(world);
        functionT brho ; // = factoryT(world);
        functionT rho = make_density_ground(world, arho, brho);

        vecfuncT vf, delrho;
        if (xc.is_dft()) {
            START_TIMER(world);
            arho.reconstruct();
            if (param.nbeta != 0 && xc.is_spin_polarized())
                brho.reconstruct();
            // brho.reconstruct();

            vf.push_back(arho);

            if (xc.is_spin_polarized())
                vf.push_back(brho);

            if (xc.is_gga()) {

                for (int axis = 0; axis < 3; ++axis)
                    delrho.push_back((*gradop[axis])(arho, false)); // delrho
                    if (xc.is_spin_polarized() && param.nbeta != 0)
                    for (int axis = 0; axis < 3; ++axis)
                        delrho.push_back((*gradop[axis])(brho, false));

                world.gop.fence(); // NECESSARY

                vf.push_back(
                             delrho[0] * delrho[0] + delrho[1] * delrho[1]
                             + delrho[2] * delrho[2]);     // sigma_aa

                if (xc.is_spin_polarized() && param.nbeta != 0)
                    vf.push_back(
                                 delrho[0] * delrho[3] + delrho[1] * delrho[4]
                                 + delrho[2] * delrho[5]); // sigma_ab
                if (xc.is_spin_polarized() && param.nbeta != 0)
                    vf.push_back(
                                 delrho[3] * delrho[3] + delrho[4] * delrho[4]
                                 + delrho[5] * delrho[5]); // sigma_bb

                world.gop.fence(); // NECESSARY
            }
            if (vf.size()) {
                reconstruct(world, vf);
//                arho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
                refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
            }

            // this is a nasty hack, just adding something so that make_libxc_args receives 5 arguments
            // has to be here or refine_to_common_level(vf) above hangs, but we really need a better solution for when nbeta=0
            if (xc.is_spin_polarized() && param.nbeta == 0 && xc.is_gga()){
                    vf.push_back(brho);
                    vf.push_back(brho);}
            END_TIMER(world, "DFT setup");
        } 


        print_meminfo(world.rank(), "Make delrho");

        // vlocal = vnuc + 2*J  

        functionT vnuc;
        vnuc = potentialmanager->vnuclear();
        START_TIMER(world);
        functionT vcoul = apply(*coulop, rho);
        END_TIMER(world, "Coulomb");
        functionT vlocal = vcoul + vnuc;
        
        vcoul.clear(false);
        vlocal.truncate();

        // BSHOperatorPtr
        std::vector<poperatorT> aopx(param.nalpha); 
        std::vector<poperatorT> bopx(param.nbeta); 
        std::vector<poperatorT> aopy(param.nalpha); 
        std::vector<poperatorT> bopy(param.nbeta); 

        make_BSHOperatorPtr(world, akx, bkx, aopx, bopx);
        if(omega != 0.0)
            make_BSHOperatorPtr(world, aky, bky, aopy, bopy);
       
        tensorT Dpolar_total(3, 3), Dpolar_alpha(3, 3), Dpolar_beta(3, 3);
        double Dpolar_average = 0.0;

        for(int axis = 0; axis<3; ++axis) {
        
            //int axis = param.response_axis;
            response_frequency(world, axis);
            if (world.rank() == 0) 
                print(" Frequency for response function = ", omega);

            // perturbation
            vecfuncT dipoleamo;// = zero_functions<double,3>(world, param.nalpha);
            vecfuncT dipolebmo;// = zero_functions<double,3>(world, param.nbeta);

            // make response function x, y
            vecfuncT ax;// = zero_functions<double,3>(world, param.nalpha);
            vecfuncT ay;// = zero_functions<double,3>(world, param.nalpha);
            vecfuncT bx;// = zero_functions<double,3>(world, param.nbeta);
            vecfuncT by;// = zero_functions<double,3>(world, param.nbeta);

            // old response function
            vecfuncT ax_old;// = zero_functions<double,3>(world, param.nalpha);
            vecfuncT ay_old; // = zero_functions<double,3>(world, param.nalpha);
            vecfuncT bx_old; // = zero_functions<double,3>(world, param.nbeta);
            vecfuncT by_old; // = zero_functions<double,3>(world, param.nbeta);

            vecfuncT axrhs; // = zero_functions<double,3>(world, param.nalpha);
            vecfuncT ayrhs; //= zero_functions<double,3>(world, param.nalpha);
            vecfuncT bxrhs; // = zero_functions<double,3>(world, param.nbeta);
            vecfuncT byrhs; // = zero_functions<double,3>(world, param.nbeta);

            dipoleamo = calc_dipole_mo(world, amo, axis, param.nalpha);
            if(!param.spin_restricted && param.nbeta != 0) {
                dipolebmo = calc_dipole_mo(world, bmo, axis, param.nbeta);
            }
            else {
                dipolebmo = dipoleamo;
            }
//goo            dipoleamo = calc_dipole_mo(world, amo, axis, param.nalpha);
//goo            if(!param.spin_restricted && param.nbeta != 0) {
//goo                dipolebmo = calc_dipole_mo(world, bmo, axis, param.nbeta);
//goo            }
//goo            else {
//goo                dipolebmo = dipoleamo;
//goo            }

            for(int iter = 0; iter < param.maxiter; ++iter) {
                if(world.rank() == 0)
                    printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());
                
                double residual = 0.0;
                double exca = 0.0, excb = 0.0;
            
                if (iter > 0 && update_residual < 0.1) {
                    //do_this_iter = false;
                    param.maxsub = maxsub_save;
                }
                
                if(iter == 0) {
                    // iter = 0 initial_guess

                    axrhs = initial_guess_response(world, dipoleamo, vlocal, vf, 0, axis);
                    if(!param.spin_restricted && param.nbeta != 0) { 
                        bxrhs = initial_guess_response(world, dipolebmo, vlocal, vf, 1, axis);
                    } 
 
                    if(omega != 0.0) {
                        ayrhs = initial_guess_response(world, dipoleamo, vlocal, vf, 0, axis);
                        if(!param.spin_restricted && param.nbeta != 0) 
                            byrhs = initial_guess_response(world, dipolebmo, vlocal, vf, 1, axis);
                    } 
                }

                else{
                    // make potential * wave function
                    vecfuncT aVx;// = zero_functions<double,3>(world, param.nalpha);
                    vecfuncT bVx; // = zero_functions<double,3>(world, param.nbeta); 

                    vecfuncT aVy;// = zero_functions<double,3>(world, param.nalpha);
                    vecfuncT bVy;// = zero_functions<double,3>(world, param.nbeta);

                    // make (dJ-dK)*2*mo
                    vecfuncT djkamox; // = zero_functions<double,3>(world, param.nalpha);
                    vecfuncT djkamoy; // = zero_functions<double,3>(world, param.nalpha);

                    vecfuncT djkbmox; // = zero_functions<double,3>(world, param.nbeta);
                    vecfuncT djkbmoy; // = zero_functions<double,3>(world, param.nbeta);


                    // calculate (dJ-dK)*2*mo
                    functionT drhoa = make_derivative_density( world, ax_old, ay_old, 0, param.nalpha );
                    functionT drhob;
                    if(!param.spin_restricted && param.nbeta != 0) {
                       drhob = make_derivative_density( world, bx_old, by_old, 1, param.nbeta );
                    } else {
                       drhob = drhoa;
                    } 
                    functionT drho = drhoa + drhob; 

                    aVx = apply_potential_response(world, aocc, ax_old, vf, vlocal, exca, 0);
                    djkamox = calc_djkmo(world, ax_old, ay_old, drho, 0);
                    // axrhs = -2.0 * (aVx + dipoleamo + duamo)
                    axrhs = calc_rhs(world, amo,  aVx, dipoleamo, djkamox);

                    if(!param.spin_restricted && param.nbeta != 0) {
                        bVx = apply_potential_response(world, bocc, bx_old, vf, vlocal, excb, 1);
                        djkbmox = calc_djkmo(world, bx_old, by_old, drho, 1);
                        bxrhs = calc_rhs(world, bmo, bVx, dipolebmo, djkbmox);
                    }
//goo                    else {
//goo                       bxrhs = axrhs;
//goo                    }

                    if(omega != 0.0) {
                        aVy = apply_potential_response(world, aocc, ay_old, vf, vlocal, exca, 0);
                        djkamoy = calc_djkmo(world, ay_old, ax_old,  drho, 0);
                        // bxrhs = -2.0 * (bVx + dipolebmo + dubmo)
                        ayrhs = calc_rhs(world, amo, aVy, dipoleamo, djkamoy);

                        if(!param.spin_restricted && param.nbeta != 0) {
                            bVy = apply_potential_response(world, bocc, by_old, vf, vlocal, excb, 1);
                            djkbmoy = calc_djkmo(world, by_old, bx_old, drho, 1);
                            byrhs = calc_rhs(world, bmo, bVy, dipolebmo, djkbmoy);
                        }
//goo                        else {
//goo                           byrhs = ayrhs;
//goo                         }
                    }
                    aVx.clear();
                    bVx.clear(); 
                    aVy.clear();
                    bVy.clear();
                    djkamox.clear();
                    djkamoy.clear();
                    djkbmox.clear();
                    djkbmoy.clear();

                }

                START_TIMER(world);
                // ax_new = G * axrhs;
                calc_response_function(world, ax, aopx, axrhs);
                orthogonalize_response(ax, 0);
                truncate(world, ax);
                axrhs.clear();

                if(!param.spin_restricted && param.nbeta != 0) {
                    // bx_new = G * bxrhs;
                    calc_response_function(world, bx, bopx, bxrhs);
                    orthogonalize_response(bx, 1);
                    truncate(world, bx);
                    bxrhs.clear();
                }
                else {
                    bx = ax;
                }

                if(omega != 0.0){
                    calc_response_function(world, ay, aopy, ayrhs);
                    orthogonalize_response(ay, 0);
                    truncate(world, ay);
                    ayrhs.clear();

                    if(!param.spin_restricted && param.nbeta != 0) {
                        calc_response_function(world, by, bopy, byrhs);
                        orthogonalize_response(by, 1);
                        truncate(world, by);
                        byrhs.clear();
                    }
                    else {
                       by = ay;
                    }
                }     
                else {
                       ay = ax;
                       by = bx;
                }
                END_TIMER(world, "Make response func");
                print_meminfo(world.rank(), "Make response func");

                if(iter > 0) {
                   // START_TIMER(world);
                    residual = 0.0;
                    
                    vecfuncT rax = zero_functions<double,3>(world, param.nalpha); //residual alpha x
                    vecfuncT ray = zero_functions<double,3>(world, param.nalpha); //residual alpha y
                    vecfuncT rbx = zero_functions<double,3>(world, param.nbeta);  //residual beta x
                    vecfuncT rby = zero_functions<double,3>(world, param.nbeta);  //residual beta y

                    double aresidual =  residual_response(world, ax, ay, ax_old, ay_old, rax, ray);
                    double bresidual = 0.0;
                    world.gop.fence(); 
                    if(!param.spin_restricted && param.nbeta != 0) { 
                        bresidual = aresidual + residual_response(world, bx, by, bx_old, by_old, rbx, rby);
                        residual = std::max(aresidual, bresidual);
                        world.gop.fence(); 
                    }
                    else { 
                      residual = aresidual;
                    } 

                    if (world.rank() == 0)  
                        print("\nresiduals_response (first) = ", residual);
#if 0
                    rax = sub(world, ax_old, ax);
                    ray = sub(world, ay_old, ay);
                    if(!param.spin_restricted && param.nbeta != 0) {
                        // bxerr = (bx_new - bx).norm2() 
                        rbx = sub(world, bx_old, bx);
                        rby = sub(world, by_old, by);
                    }
#endif
                    residual = 0.0;
#if 0
                    residual = norm2(world, sub(world, ax_old, ax)) + norm2(world, sub(world, ay_old, ay));
                    if(!param.spin_restricted && param.nbeta != 0) 
                        residual = residual + norm2(world, sub(world, bx_old, bx)) + norm2(world, sub(world, by_old, by));
                    else if(param.nbeta != 0) residual = 2.0 * residual; 
#endif
                    double nx,ny;
                    ////////UPDATE
                    nx=norm2(world, ax);
                    if (world.rank() == 0) 
                        print("CURRENT_X_norm2() = ", nx);
                    update_response_subspace(world, ax, ay, bx, by, rax, ray, rbx, rby, subspace, Q, update_residual); 

                    nx=norm2(world, ax);
                    ny=norm2(world, ay);
                    if (world.rank() == 0) { 
                        print("new X (alpha) norm2() = ", nx);
                        print("new Y (alpha) norm2() = ", ny);
                    }

                    aresidual = residual_response(world, ax, ay, ax_old, ay_old, rax, ray);
                    bresidual = 0.0;
                    
                    if(!param.spin_restricted && param.nbeta != 0) { 
                        bresidual = residual_response(world, bx, by, bx_old, by_old, rbx, rby);
                        residual = std::max(aresidual, bresidual);
                    }
                    else residual = aresidual;

                    double thresh = rconv *(param.nalpha + param.nbeta)*2; 
                    if (world.rank() == 0) { 
                        print("\nresiduals_response (final) = ", residual);
                        print("rconv *(param.nalpha + param.nbeta)*2", thresh);
                    }

                  //  END_TIMER(world, "Update response func");
                    print_meminfo(world.rank(), "Update response func");
                    
                    // 収束判定
                    if( residual < (rconv *(param.nalpha + param.nbeta)*2))
                    {
                        if (world.rank() == 0) { 
                            print("\nConverged response function!!\n");
                            print("\n\n\n");
                            print(" ------------------------------------------------------------------------------");
                            print(" |                  MADNESS CALCULATION POLARIZABILITY                        |");
                            print(" ------------------------------------------------------------------------------");
                            print(" \n\n");
                        } 
                        break; 
                    }
                }

                ax_old = ax;
                bx_old = bx;

                ay_old = ay;
                by_old = by;
                    
            } //end iteration
            END_TIMER(world, "Make response func");
            print_meminfo(world.rank(), "Make response func");

            calc_dpolar(world, Dpolar_average, ax, ay, bx, by, axis, Dpolar_total, Dpolar_alpha, Dpolar_beta);
            
            ax.clear();
            ay.clear();
            bx.clear();
            by.clear();
            ax_old.clear();
            ay_old.clear();
            bx_old.clear();
            by_old.clear();

            dipoleamo.clear();
            dipolebmo.clear();

            if(axis == 2) break;
        }
    } // end solve response function
};
// Computes molecular energy as a function of the geometry
// This is cludgy ... need better factorization of functionality
// between calculation, main program and this ... or just merge it all.
class MolecularEnergy : public OptimizationTargetInterface {
    World& world;
    Calculation& calc;
    mutable double coords_sum;     // sum of square of coords at last solved geometry
    mutable double E; //< Current energy

    public:
    MolecularEnergy(World& world, Calculation& calc)
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
                    if(calc.param.polar == true) {
                        if(calc.param.polar_ds == 0)
                            calc.initial_guess_mp(world);
                    }
                    //calc.param.restart = true;
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
            //if(proto==0){
            // if(calc.param.absolvent)//||calc.param.svpe)
            //     calc.solve_gas_phase(world);// computes vacuo density and energy
            //     //}//  calc.save_mos(world); //debugging

            calc.solve(world);
                    if(world.rank() == 0) {
                        print(" ");
                        print("alpha eigenvalues 1");
                        print(calc.aeps);
                    }
            if(calc.param.polar == true) {
                if(calc.param.polar_ds == 0) {
//dsolmp                    calc.solve_finite_mo(world);
//dsolmp                    calc.solve_spolar(world);
                }
                else if (calc.param.polar_ds == 1) {
                    calc.solve_response(world);
                }
            }
            if (calc.param.save)
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
            Calculation calc(world, inpname);

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
            //vama
            calc.set_protocol(world,calc.param.protocol_data[0]);


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
//vama1            else if (calc.param.tdksprop) {
//vama1                print("\n\n Propagation of Kohn-Sham equation                      ");
//vama1                print(" ----------------------------------------------------------\n");
//vama1                //          calc.propagate(world,VextCosFunctor<double>(world,new DipoleFunctor(2),0.1),0);
//vama1                calc.propagate(world,0.1,0);
//vama1            }
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
