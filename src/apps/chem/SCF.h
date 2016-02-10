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


#ifndef MADNESS_CHEM_SCF_H__INCLUDED
#define MADNESS_CHEM_SCF_H__INCLUDED

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/mra/mra.h>

#include <chem/molecule.h>
#include <chem/molecularbasis.h>
#include <chem/corepotential.h>
#include <chem/xcfunctional.h>
#include <chem/potentialmanager.h>
#include <chem/gth_pseudopotential.h> 

#include <madness/tensor/solvers.h>
#include <madness/tensor/distributed_matrix.h>


namespace madness {

typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef DistributedMatrix<double> distmatT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>,3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;

    extern distmatT distributed_localize_PM(World & world,
                                const vecfuncT & mo,
                                const vecfuncT & ao,
                                const std::vector<int> & set,
                                const std::vector<int> & at_to_bf,
                                const std::vector<int> & at_nbf,
                                const double thresh = 1e-9,
                                const double thetamax = 0.5,
                                const bool randomize = true,
                                       const bool doprint = false);

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

static double mask3(const coordT& ruser) {
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


class AtomicAttractionFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int iatom;

public:
    AtomicAttractionFunctor(const Molecule& molecule, int iatom)
        : molecule(molecule), iatom(iatom) {}

    double operator()(const coordT& x) const {
        const Atom& atom=molecule.get_atom(iatom);
        const coordT coord={atom.x,atom.y,atom.z};
        double r = (x-coord).normf();
        return -atom.q * smoothed_potential(r*molecule.get_rcut()[iatom])
                *molecule.get_rcut()[iatom];
    }

    std::vector<coordT> special_points() const {
        return std::vector<coordT>(1,molecule.get_atom(iatom).get_coords());
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

class MolecularSecondDerivativeFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int atom;
    const int iaxis, jaxis;

public:
    MolecularSecondDerivativeFunctor(const Molecule& molecule, int atom,
            int iaxis, int jaxis)
        : molecule(molecule), atom(atom),iaxis(iaxis), jaxis(jaxis)
    {}

    double operator()(const coordT& x) const {
        return molecule.nuclear_attraction_potential_second_derivative(atom,
                iaxis, jaxis, x[0], x[1], x[2]);
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
    double econv;               ///< Energy convergence
    double dconv;               ///< Density convergence
    int k;						///< polynomial order
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
    bool no_compute;            ///< If true use orbitals on disk, set value to computed
    bool no_orient;				///< If true the molecule coordinates will not be reoriented
    bool save;                  ///< If true save orbitals to disk
    unsigned int maxsub;        ///< Size of iterative subspace ... set to 0 or 1 to disable
    double orbitalshift;		///< scf orbital shift: shift the occ orbitals to lower energies
    int npt_plot;               ///< No. of points to use in each dim for plots
    tensorT plot_cell;          ///< lo hi in each dimension for plotting (default is all space)
    std::string aobasis;        ///< AO basis used for initial guess (6-31g or sto-3g)
    std::string core_type;      ///< core potential type ("" or "mcp")
    bool derivatives;           ///< If true calculate derivatives
    bool dipole;                ///< If true calculate dipole moment
    bool conv_only_dens;        ///< If true remove bsh_residual from convergence criteria   how ugly name is...
    bool psp_calc;              ///< pseudopotential calculation for all atoms
    bool print_dipole_matels;   ///< If true output dipole matrix elements
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
    int  gmaxiter;              ///< optimization maxiter
    bool ginitial_hessian;      ///< compute inital hessian for optimization
    std::string algopt;         ///< algorithm used for optimization
    bool hessian;               ///< compute the hessian matrix
    bool read_cphf;             ///< read the orbital response for nuclear displacements from file
    bool tdksprop;               ///< time-dependent Kohn-Sham equation propagate
    std::string nuclear_corrfac;	///< nuclear correlation factor
    bool pure_ae;                 ///< pure all electron calculation with no pseudo-atoms
    int nv_factor;              ///< factor to multiply number of virtual orbitals with when automatically decreasing nvirt

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & k & L & maxrotn & nvalpha & nvbeta
           & nopen & maxiter & nio & spin_restricted;
        ar & plotlo & plothi & plotdens & plotcoul & localize & localize_pm
           & restart & save & no_compute &no_orient & maxsub & orbitalshift & npt_plot & plot_cell & aobasis;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
        ar & core_type & derivatives & conv_only_dens & dipole;
        ar & xc_data & protocol_data;
        ar & gopt & gtol & gtest & gval & gprec & gmaxiter & ginitial_hessian & algopt & tdksprop
            & nuclear_corrfac & psp_calc & print_dipole_matels & pure_ae & hessian & read_cphf;
    }

    CalculationParameters()
        : charge(0.0)
        , smear(0.0)
        , econv(1e-5)
        , dconv(1e-4)
    	, k(-1)
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
    	, no_compute(false)
    	, no_orient(false)
        , save(true)
        , maxsub(5)
    	, orbitalshift(0.0)
        , npt_plot(101)
        , aobasis("6-31g")
        , core_type("")
        , derivatives(false)
        , dipole(false)
        , conv_only_dens(false)
        , psp_calc(false)
        , print_dipole_matels(false)
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
        , ginitial_hessian(false)
        , algopt("BFGS")
        , hessian(false)
        , read_cphf(false)
        , tdksprop(false)
        , nuclear_corrfac("none")
        , pure_ae(true)
        , nv_factor(1)
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
            else if (s == "econv") {
                f >> econv;
            }
            else if (s == "dconv") {
                f >> dconv;
            }
            else if (s == "k") {
                f >> k;
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
                //can't redirect true/false as with other variables so create temporary variable
                std::string tmp_save;
                f >> tmp_save;
                if (tmp_save=="true"){
                    save=true;
                }
                else if (tmp_save=="false"){
                    save=false;
                }
                else{
                    std::cout << "moldft: unrecognized value for save (true or false only): " << tmp_save << std::endl;
                    MADNESS_EXCEPTION("input_error", 0);
                }
            }
            else if (s == "no_compute") {
                no_compute = true;
            }
            else if (s == "no_orient") {
            	no_orient = true;
            }
            else if (s == "maxsub") {
                f >> maxsub;
                if (maxsub <= 0) maxsub = 1;
                if (maxsub > 20) maxsub = 20;
            }
            else if (s == "orbitalshift") {
                f >> orbitalshift;
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
            else if (s == "ginitial_hessian") {
                ginitial_hessian = true;
            }
            else if (s == "algopt") {
                f >> algopt;

//                char buf[1024];
//                f.getline(buf,sizeof(buf));
//                algopt = buf;
            }
            else if (s == "hessian") {
               hessian = true;
            }
            else if (s == "read_cphf") {
                read_cphf = true;
            }
            else if (s == "tdksprop") {
              tdksprop = true;
            }
            else if (s == "nuclear_corrfac") {
                std::string str;
                std::getline(f,str);
            	nuclear_corrfac=str;
            }
            else if (s == "psp_calc") {
              psp_calc = true;
              pure_ae = false;
            }
            else if (s == "print_dipole_matels") {
              print_dipole_matels = true;
            }
            else if (s == "nv_factor") {
              f >> nv_factor;
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
        madness::print("         xc library  ", "libxc");
#else
        madness::print("         xc library  ", "default (lda only)");
#endif
        if (core_type != "")
            madness::print("           core type ", core_type);
        madness::print(" initial guess basis ", aobasis);
        madness::print(" max krylov subspace ", maxsub);
        madness::print("    compute protocol ", protocol_data);
        madness::print("  energy convergence ", econv);
        madness::print(" density convergence ", dconv);
        madness::print("    maximum rotation ", maxrotn);
        madness::print("    polynomial order ", k);
        madness::print("       truncate mode ", FunctionDefaults<3>::get_truncate_mode());
        madness::print("  maximum iterations ", maxiter);
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
        if (psp_calc)
            madness::print(" psp or all electron ", "pseudopotential");
        else if (pure_ae)
            madness::print(" psp or all electron ", "all electron");
        else
            madness::print(" psp or all electron ", "mixed psp/AE");
    }

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

class SCF {
public:
    std::shared_ptr<PotentialManager> potentialmanager;
    std::shared_ptr<GTHPseudopotential <double> > gthpseudopotential;
    Molecule molecule;
    CalculationParameters param;
    XCfunctional xc;
    AtomicBasisSet aobasis;
    functionT mask;

    /// alpha and beta molecular orbitals
    vecfuncT amo, bmo;

    /// sets of orbitals grouped by their orbital energies (for localization?)
    /// only orbitals within the same set will be mixed to localize
    std::vector<int> aset, bset;

    /// MRA projection of the minimal basis set
    vecfuncT ao;

    std::vector<int> at_to_bf, at_nbf;

    /// occupation numbers for alpha and beta orbitals
    tensorT aocc, bocc;

    /// orbital energies for alpha and beta orbitals
    tensorT aeps, beps;
    poperatorT coulop;
    std::vector< std::shared_ptr<real_derivative_3d> > gradop;
    double vtol;
    double current_energy;
    //double esol;//etot;
    //double vacuo_energy;
    static const int vnucextra = 12; // load balance parameter for nuclear pot.
    static const int loadbalparts = 2.0; // was 6.0

    SCF(World & world, const char *filename);

    template<std::size_t NDIM>
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

        // k defaults to make sense with thresh, override by providing k in input file
        if (param.k == -1) {
        	FunctionDefaults<NDIM>::set_k(k);
//        	param.k=k;
        } else {
        	FunctionDefaults<NDIM>::set_k(param.k);
        }
        // don't forget to adapt the molecular smoothing parameter!!
//        molecule.set_eprec(std::min(thresh,molecule.get_eprec()));
        FunctionDefaults<NDIM>::set_thresh(thresh);
        FunctionDefaults<NDIM>::set_refine(true);
        FunctionDefaults<NDIM>::set_initial_level(2);
//        FunctionDefaults<NDIM>::set_truncate_mode(1);
        FunctionDefaults<NDIM>::set_autorefine(false);
        FunctionDefaults<NDIM>::set_apply_randomize(false);
        FunctionDefaults<NDIM>::set_project_randomize(false);
        FunctionDefaults<NDIM>::set_cubic_cell(-param.L, param.L);
        GaussianConvolution1DCache<double>::map.clear();
        double safety = 0.1;
        vtol = FunctionDefaults<NDIM>::get_thresh() * safety;
        coulop = poperatorT(CoulombOperatorPtr(world, param.lo, thresh));
        gradop = gradient_operator<double,3>(world);
        mask = functionT(factoryT(world).f(mask3).initial_level(4).norefine());
        if(world.rank() == 0){
            print("\nSolving NDIM=",NDIM," with thresh", thresh, "    k",
            	FunctionDefaults<NDIM>::get_k(), "   conv", std::max(thresh, param.dconv), "\n");
        }
    }

    /// getter for the molecular orbitals, alpha spin
    const vecfuncT& get_amo() const {return amo;}

    /// getter for the molecular orbitals, beta spin
    const vecfuncT& get_bmo() const {return bmo;}

    /// getter for the occupation numbers, alpha spin
    const tensorT& get_aocc() const {return aocc;}

    /// getter for the occupation numbers, alpha spin
    const tensorT& get_bocc() const {return bocc;}

    bool is_spin_restricted() const {return param.spin_restricted;}

    void save_mos(World& world);

    void load_mos(World& world);

    void do_plots(World& world);

    void project(World & world);

    void make_nuclear_potential(World & world);

    void project_ao_basis(World & world);

    /// group orbitals into sets of similar orbital energies for localization

    /// @param[in]	eps	orbital energies
    /// @param[in]	occ	occupation numbers
    /// @param[in]	nmo number of MOs for the given spin
    /// @return		vector of length nmo with the set index for each MO
    std::vector<int> group_orbital_sets(World& world, const tensorT& eps,
    		const tensorT& occ, const int nmo) const;

    /// compute the unitary localization matrix according to Pipek-Mezey

    /// @param[in]	world	the world
    /// @param[in]	mo		the MOs
    /// @param[in]	set		only orbitals within the same set will be mixed
    /// @param[in]	thresh	the localization threshold
    /// @param[in]	thetamax	??
    /// @param[in]	randomize	??
    distmatT localize_PM(World & world, const vecfuncT & mo, const std::vector<int> & set,
    		const double thresh = 1e-9, const double thetamax = 0.5,
    		const bool randomize = true, const bool doprint = false) const;


    void analyze_vectors(World & world, const vecfuncT & mo, const tensorT & occ = tensorT(),
    		const tensorT & energy = tensorT(), const std::vector<int> & set = std::vector<int>());

    inline double DIP(const tensorT & dip, int i, int j, int k, int l) {
        return dip(i, j, 0) * dip(k, l, 0) + dip(i, j, 1) * dip(k, l, 1) + dip(i, j, 2) * dip(k, l, 2);
    }

    // tensorT localize_boys(World & world, const vecfuncT & mo, const std::vector<int> & set,
    // 		const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true);


    distmatT kinetic_energy_matrix(World & world, const vecfuncT & v) const;
    distmatT kinetic_energy_matrix(World & world, const vecfuncT & vbra, const vecfuncT & vket) const;

    vecfuncT core_projection(World & world, const vecfuncT & psi, const bool include_Bc = true);

    double core_projector_derivative(World & world, const vecfuncT & mo,
    			const tensorT & occ, int atom, int axis);

    void initial_guess(World & world);

    void initial_load_bal(World & world);

    functionT make_density(World & world, const tensorT & occ, const vecfuncT & v) const;

    functionT make_density(World & world, const tensorT & occ, const cvecfuncT & v);

    real_function_3d make_sigma(const real_function_3d& rho1,
            const real_function_3d& rho2) const;


    std::vector<poperatorT> make_bsh_operators(World & world, const tensorT & evals) const;
    std::vector<poperatorT> make_gradbsh_operators(World & world,
            const tensorT & evals, const int axis) const;

    /// apply the HF exchange on a set of orbitals

    /// @param[in]  world   the world
    /// @param[in]  occ     occupation numbers
    /// @param[in]  psi     the orbitals in the exchange operator
    /// @param[in]  f       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT apply_hf_exchange(World & world, const tensorT & occ,
    		const vecfuncT & psi, const vecfuncT & f) const ;

    // Used only for initial guess that is always spin-restricted LDA
    functionT make_lda_potential(World & world, const functionT & arho);


//    functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin, int what)
//    {
//        return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
//    }

    double make_dft_energy(World & world, const vecfuncT& vf, int ispin)
    {
        functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc), vf);
        return vlda.trace();
    }

	vecfuncT apply_potential(World & world, const tensorT & occ,
			const vecfuncT & amo,
			const functionT & vlocal, double & exc, double & enl, int ispin);

    tensorT derivatives(World & world, const functionT& rho) const;


    /// compute the three components of \f[\nabla f\f]
    std::vector<real_function_3d> nabla(const real_function_3d& f,
            const bool do_refine=false) const {
        World& world=f.world();

        f.reconstruct();
        if (do_refine) f.refine();      // refine to make result more precise

        std::vector< std::shared_ptr< Derivative<double,3> > > grad=
                gradient_operator<double,3>(world);

        std::vector<real_function_3d> result(3);
        for (int i=0; i<3; ++i) {
            result[i]=apply(*(grad[i]),f,false);
        }
        world.gop.fence();
        return result;

    }

    /// compute the total dipole moment of the molecule

    /// @param[in]  rho the total (alpha + beta) density
    /// @return     the x,y,z components of the el. + nucl. dipole moment
    tensorT dipole(World & world, const functionT& rho) const;

    void dipole_matrix_elements(World & world, const vecfuncT & mo, const tensorT & occ = tensorT(),
    		const tensorT & energy = tensorT(), int spin=0);

    void vector_stats(const std::vector<double> & v, double & rms,
    		double & maxabsval) const;

	vecfuncT compute_residual(World & world, tensorT & occ, tensorT & fock,
			const vecfuncT & psi, vecfuncT & Vpsi, double & err);

	tensorT make_fock_matrix(World & world, const vecfuncT & psi,
			const vecfuncT & Vpsi, const tensorT & occ,
			double & ekinetic) const;


    /// make the Coulomb potential given the total density
    functionT make_coulomb_potential(const functionT& rho) const {
        return apply(*coulop, rho);
    }

    /// Compute the two-electron integrals over the provided set of orbitals

    /// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
    /// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
    Tensor<double> twoint(World& world, const vecfuncT& psi) const;

    tensorT matrix_exponential(const tensorT& A) const ;

    /// compute the unitary transformation that diagonalizes the fock matrix

    /// @param[in]	world	the world
    /// @param[in]	overlap	the overlap matrix of the orbitals
    /// @param[in,out]	fock	the fock matrix; diagonal upon exit
    /// @param[out]	evals	the orbital energies
    /// @param[in]	occ	the occupation numbers
    /// @param[in]	thresh_degenerate	threshold for orbitals being degenerate
    /// @return		the unitary matrix U: U^T F U = evals
    tensorT get_fock_transformation(World& world, const tensorT& overlap,
    		tensorT& fock, tensorT& evals, const tensorT& occ,
    		const double thresh_degenerate) const;


    /// diagonalize the fock matrix, taking care of degenerate states

    /// Vpsi is passed in to make sure orbitals and Vpsi are in phase
    /// @param[in]	world	the world
    /// @param[in,out]	fock	the fock matrix (diagonal upon exit)
    /// @param[in,out]	psi		the orbitals
    /// @param[in,out]	Vpsi	the orbital times the potential
    /// @param[out]	evals	the orbital energies
    /// @param[in]	occ		occupation numbers
    /// @param[in]	thresh	threshold for rotation and truncation
    /// @return		the unitary matrix U: U^T F U = evals
    tensorT diag_fock_matrix(World& world, tensorT& fock,
    		vecfuncT& psi, vecfuncT& Vpsi, tensorT& evals,
    		const tensorT& occ, const double thresh) const;


    void loadbal(World & world, functionT & arho, functionT & brho, functionT & arho_old,
    		functionT & brho_old, subspaceT & subspace);


    void rotate_subspace(World& world, const tensorT& U, subspaceT& subspace,
    		int lo, int nfunc, double trantol) const;

    void rotate_subspace(World& world, const distmatT& U, subspaceT& subspace,
    		int lo, int nfunc, double trantol) const;

    void update_subspace(World & world,
                         vecfuncT & Vpsia, vecfuncT & Vpsib,
                         tensorT & focka, tensorT & fockb,
                         subspaceT & subspace, tensorT & Q,
                         double & bsh_residual, double & update_residual);

    /// perform step restriction following the KAIN solver

    /// undo the rotation from the KAIN solver if the rotation exceeds the
    /// maxrotn parameter
    /// @param[in]		world	the world
    /// @param[in]		mo		vector of orbitals from previous iteration
    /// @param[in,out]	mo_new	vector of orbitals from the KAIN solver
    /// @param[in]		spin	"alpha" or "beta" for user information
    /// @return			max residual
    double do_step_restriction(World& world, const vecfuncT& mo,
    		vecfuncT& mo_new, std::string spin) const;


    /// orthonormalize the vectors

    /// @param[in]		world	the world
    /// @param[in,out]	amo_new	the vectors to be orthonormalized
    void orthonormalize(World& world, vecfuncT& amo_new) const;

    void orthonormalize(World& world, vecfuncT& amo_new, int nocc) const;



    void propagate(World& world, double omega, int step0);


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


    void iterate_trotter(World& world, Convolution1D<double_complex>* G,
			cvecfuncT& camo, cvecfuncT& cbmo, double t, double time_step,
			double thresh);


    // For given protocol, solve the DFT/HF/response equations
    void solve(World & world);

};


// Computes molecular energy as a function of the geometry
// This is cludgy ... need better factorization of functionality
// between calculation, main program and this ... or just merge it all.
class MolecularEnergy : public OptimizationTargetInterface {
    World& world;
    SCF& calc;
    mutable double coords_sum;     // sum of square of coords at last solved geometry

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

	  // read converged wave function from disk if there is one
	  if (calc.param.no_compute) {
		calc.load_mos(world);
		calc.make_nuclear_potential(world);
		calc.project_ao_basis(world);
		return calc.current_energy;
	  }

	  int nvalpha = calc.param.nmo_alpha - calc.param.nalpha;
	  int nvbeta = calc.param.nmo_beta - calc.param.nbeta;
	  int nvalpha_start, nv_old;

        // The below is missing convergence test logic, etc.

        // Make the nuclear potential, initial orbitals, etc.
        for (unsigned int proto=0; proto<calc.param.protocol_data.size(); proto++) {

            //repeat with gradually decreasing nvirt, only for first protocol
            if (proto == 0 && nvalpha > 0){
                nvalpha_start = nvalpha * calc.param.nv_factor;}
            else{
                nvalpha_start = nvalpha;}

            nv_old = nvalpha_start;

            for (int nv=nvalpha_start;nv>=nvalpha;nv-=nvalpha){

                if (nv > 0 && world.rank() == 0) std::cout << "Running with " << nv << " virtual states" << std::endl;
            
                calc.param.nmo_alpha = calc.param.nalpha + nv;
                // check whether this is sensible for spin restricted case
                if (calc.param.nbeta && !calc.param.spin_restricted){
                    if (nvbeta == nvalpha){
                        calc.param.nmo_beta = calc.param.nbeta + nv;}
                    else{
                        calc.param.nmo_beta = calc.param.nbeta + nv + nvbeta - nvalpha;}
                }

                calc.set_protocol<3>(world,calc.param.protocol_data[proto]);
                calc.make_nuclear_potential(world);

                calc.project_ao_basis(world);

                if (proto == 0 && nv == nvalpha_start) {
                    if (calc.param.restart) {
                        calc.load_mos(world);
                    }
                    else {
                        calc.initial_guess(world);
                        //calc.param.restart = true;
                    }
                }
                else {
                   if (nv != nv_old){
                       calc.amo.resize(calc.param.nmo_alpha);
                       calc.bmo.resize(calc.param.nmo_beta);

                       calc.aocc = tensorT(calc.param.nmo_alpha);
                       for (int i = 0; i < calc.param.nalpha; ++i)
                           calc.aocc[i] = 1.0;

                       calc.bocc = tensorT(calc.param.nmo_beta);
                       for (int i = 0; i < calc.param.nbeta; ++i)
                           calc.bocc[i] = 1.0;

                       // might need to resize aset, bset, but for the moment this doesn't seem to be necessary

                   }
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

                if (calc.param.save)
                  calc.save_mos(world);

                nv_old=nv;
                // exit loop over decreasing nvirt if nvirt=0
                if (nv==0) break;

            }

        }
        return calc.current_energy;
    }

    madness::Tensor<double> gradient(const Tensor<double>& x) {
        value(x); // Ensures DFT equations are solved at this geometry


        functionT rho = calc.make_density(world, calc.aocc, calc.amo);
        functionT brho = rho;
        if (!calc.param.spin_restricted)
            brho = calc.make_density(world, calc.bocc, calc.bmo);
        rho.gaxpy(1.0, brho, 1.0);

        return calc.derivatives(world,rho);
    }
};
}

#endif /* SCF_H_ */

