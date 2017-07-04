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

/// \file CalculationParameters
/// \brief solution parameters for SCF calculations



#ifndef MADNESS_CHEM_CALCULATIONPARAMETERS_H__INCLUDED
#define MADNESS_CHEM_CALCULATIONPARAMETERS_H__INCLUDED

#include <chem/molecule.h>
#include <chem/molecularbasis.h>

namespace madness {


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
    int k;                      ///< polynomial order
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
    bool restartao;             ///< If true restart from orbitals projected into AO basis (STO3G) on disk
    bool no_compute;            ///< If true use orbitals on disk, set value to computed
    bool no_orient;             ///< If true the molecule coordinates will not be reoriented
    bool save;                  ///< If true save orbitals to disk
    unsigned int maxsub;        ///< Size of iterative subspace ... set to 0 or 1 to disable
    double orbitalshift;        ///< scf orbital shift: shift the occ orbitals to lower energies
    int npt_plot;               ///< No. of points to use in each dim for plots
    Tensor<double> plot_cell;   ///< lo hi in each dimension for plotting (default is all space)
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
    bool read_cphf;             ///< read the converged orbital response for nuclear displacements from file
    bool restart_cphf;          ///< read the guess orbital response for nuclear displacements from file
    bool purify_hessian;        ///< symmetrize the hessian matrix based on atomic charges
    bool tdksprop;               ///< time-dependent Kohn-Sham equation propagate
    std::string nuclear_corrfac;    ///< nuclear correlation factor
    bool pure_ae;                 ///< pure all electron calculation with no pseudo-atoms
    int nv_factor;              ///< factor to multiply number of virtual orbitals with when automatically decreasing nvirt
    int vnucextra; // load balance parameter for nuclear pot.
    int loadbalparts = 2; // was 6
    std::string pcm_data;            ///< do a PCM (solvent) calculation
    std::string ac_data;             ///< do a calculation with asymptotic correction (see ACParameters class in chem/AC.h for details)

    // Next list for response code from a4v4
    bool response;                    ///< response function calculation
    double response_freq;             ///< Frequency for calculation response function
    std::vector<bool> response_axis;  ///< Calculation protocol
    bool nonrotate;                   ///< If true do not molcule orient
    double rconv;                     ///< Response convergence
    double efield;                    ///< eps for finite field
    double efield_axis;               ///< eps for finite field axis


    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & k & L & maxrotn & nvalpha & nvbeta
        & nopen & maxiter & nio & spin_restricted;
        ar & plotlo & plothi & plotdens & plotcoul & localize & localize_pm
        & restart & restartao & save & no_compute &no_orient & maxsub & orbitalshift & npt_plot & plot_cell & aobasis;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
        ar & core_type & derivatives & conv_only_dens & dipole;
        ar & xc_data & protocol_data;
        ar & gopt & gtol & gtest & gval & gprec & gmaxiter & ginitial_hessian & algopt & tdksprop
        & nuclear_corrfac & psp_calc & print_dipole_matels & pure_ae & hessian & read_cphf & restart_cphf
        & purify_hessian & vnucextra & loadbalparts & pcm_data & ac_data;
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
    , restartao(false)
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
    , restart_cphf(false)
    , purify_hessian(false)
    , tdksprop(false)
    , nuclear_corrfac("none")
    , pure_ae(true)
    , nv_factor(1)
    , vnucextra(12)
    , loadbalparts(2)
    , pcm_data("none")
    , ac_data("none")
    , response(false)
    , response_freq(0.0)
    , response_axis(madness::vector_factory(true, true, true))
    , nonrotate(false)
    , rconv(1e-6)
    , efield(0.0)
    , efield_axis(0)
    {}

    // initializes CalculationParameters using the contents of file \c filename
    void read_file(const std::string& filename) {
        std::ifstream f(filename.c_str());
        read(f);
    }

    // initializes CalculationParameters using the contents of stream \c f
    void read(std::istream& f) {
        position_stream(f, "dft");
        std::string s;
        xc_data = "lda";
        protocol_data = madness::vector_factory(1e-4, 1e-6);

        while (f >> s) {
            if (s == "end") {
                break;
            }
            else if (s == "loadbal") {
                f >> vnucextra >> loadbalparts;
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
                plot_cell = Tensor<double>(3L,2L);
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
            else if (s == "restartao") {
                restartao = true;
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
            else if (s == "restart_cphf") {
                restart_cphf = true;
            }
            else if (s == "purify_hessian") {
                purify_hessian = true;
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
            else if (s=="pcm") {
                char buf[1024];
                f.getline(buf,sizeof(buf));
                pcm_data = buf;
            }
            else if (s=="ac") {
                char buf[1024];
                f.getline(buf,sizeof(buf));
                ac_data = buf;
            }
            else if (s == "print_dipole_matels") {
                print_dipole_matels = true;
            }
            else if (s == "nv_factor") {
                f >> nv_factor;
            }
            else if (s == "response") {
                response = true;
            }
            else if (s == "response_freq") {
                double freq;
                f >> freq;
                response_freq = freq;
            }
            else if (s == "response_axis") {
                std::string buf;
                std::getline(f,buf);
                response_axis = std::vector<bool>();
                bool d;
                std::stringstream s(buf);
                while (s >> d) response_axis.push_back(d);
            }
            else if (s == "nonrotate") {
                nonrotate = true;
            }
            else if (s == "rconv") {
                f >> rconv;
            }
            else if (s == "efield") {
                f >> efield;
            }
            else if (s == "efield_axis") {
                std::string axis;
                f >> axis;
                if(axis == "x")
                    efield_axis = 0;
                else if (axis == "y")
                    efield_axis = 1;
                else if (axis == "z")
                    efield_axis = 2;
                else if (axis == "none")
                    efield_axis = -1;
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
        madness::print("    restart from AOs ", restartao);
        madness::print(" number of processes ", world.size());
        madness::print("   no. of io servers ", nio);
        madness::print("   vnuc load bal fac ", vnucextra);
        madness::print("      load bal parts ", loadbalparts);
        madness::print("     simulation cube ", -L, L);
        madness::print("        total charge ", charge);
        madness::print("            smearing ", smear);
        madness::print(" number of electrons ", nalpha, nbeta);
        madness::print("  number of orbitals ", nmo_alpha, nmo_beta);
        madness::print("     spin restricted ", spin_restricted);
        madness::print("       xc functional ", xc_data);
#ifdef MADNESS_HAS_LIBXC
        madness::print("          xc library ", "libxc");
#else
        madness::print("          xc library ", "default (lda only)");
#endif
#ifdef MADNESS_HAS_PCM
            madness::print("          PCM solver ", pcm_data);
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

} // namespace madness

#endif /* MADNESS_CHEM_CALCULATIONPARAMETERS_H__INCLUDED */
