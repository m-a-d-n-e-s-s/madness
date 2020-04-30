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
#include <chem/QCCalculationParametersBase.h>


namespace madness {

#if 1
struct CalculationParameters : public QCCalculationParametersBase {

	CalculationParameters(const CalculationParameters& other) : QCCalculationParametersBase(other) {
	}

	/// ctor reading out the input file
	CalculationParameters() {

		initialize<double>("charge",0.0,"total molecular charge");
		initialize<std::string> ("xc","hf","XC input line");
		initialize<double>("smear",0.0,"smearing parameter");
		initialize<double>("econv",1.e-5,"energy convergence");
		initialize<double>("dconv",1.e-4,"density convergence");
		initialize<bool>  ("converge_each_energy",false,"converge all fock operator components");
		initialize<int>   ("k",-1,"polynomial order");
		initialize<double>("l",20,"user coordinates box size");
		initialize<std::string>("deriv","abgv","derivative method",{"abgv","bspline","ble"});
		initialize<std::string>("dft_deriv","abgv","derivative method for gga potentials",{"abgv","bspline","ble"});
		initialize<double>("maxrotn",0.25,"step restriction used in autoshift algorithm");
		initialize<int>   ("nvalpha",0,"number of alpha virtuals to compute");
		initialize<int>   ("nvbeta",0,"number of beta virtuals to compute");
		initialize<int>   ("nopen",0,"number of unpaired electrons = nalpha-nbeta");
		initialize<int>   ("maxiter",25,"maximum number of iterations");
		initialize<int>   ("nio",1,"no. of io servers to use");
		initialize<bool>  ("spin_restricted",true,"true if spin restricted");
		initialize<int>   ("plotlo",0,"range of MOs to print (for both spins if polarized");
		initialize<int>   ("plothi",-1,"range of MOs to print (for both spins if polarized");
		initialize<bool>  ("plotdens",false,"If true print the density at convergence");
		initialize<bool>  ("plotcoul",false,"If true plot the total coulomb potential at convergence");
		initialize<std::string> ("localize","new","localization method",{"pm","boys","new","canon"});
//		initialize<bool localize_pm;           ///< If true use PM for localization
//		initialize<bool localize_boys;         ///< If true use boys for localization
//		initialize<bool localize_new;          ///< If true use new for localization
		initialize<std::string> ("pointgroup","c1","use point (sub) group symmetry if not localized",{"c1","c2","ci","cs","c2v","c2h","d2","d2h"});
		initialize<bool>  ("restart",false,"if true restart from orbitals on disk");
		initialize<bool>  ("restartao",false,"if true restart from orbitals projected into AO basis (STO3G) on disk");
		initialize<bool>  ("no_compute",false,"if true use orbitals on disk, set value to computed");
		initialize<bool>  ("no_orient",false,"if true the molecule coordinates will not be reoriented");
		initialize<bool>  ("save",true,"if true save orbitals to disk");
		initialize<int>   ("maxsub",5,"size of iterative subspace ... set to 0 or 1 to disable");
		initialize<double> ("orbitalshift",0.0,"scf orbital shift: shift the occ orbitals to lower energies");
		initialize<int>    ("npt_plot",101,"no. of points to use in each dim for plots");
//		initialize<Tensor<double> > ("plot_cell",Tensor<double>(),"lo hi in each dimension for plotting (default is all space)");
		initialize<std::vector<double> > ("plot_cell",std::vector<double>(),"lo hi in each dimension for plotting (default is all space)");
		initialize<std::string> ("aobasis","6-31g","AO basis used for initial guess (6-31g or sto-3g)");
		initialize<std::string> ("core_type","none","core potential type",{"none","mpc"});
		initialize<bool> ("derivatives",false,"if true calculate nuclear derivatives");
		initialize<bool> ("dipole",false,"if true calculate dipole moment");
		initialize<bool> ("conv_only_dens",false,"if true remove bsh_residual from convergence criteria");
		initialize<bool> ("psp_calc",false,"pseudopotential calculation for all atoms");
		initialize<bool> ("print_dipole_matels",false,"if true output dipole matrix elements");
		initialize<std::string> ("pcm_data","none","do a PCM (solvent) calculation");
		initialize<std::string> ("ac_data","none","do a calculation with asymptotic correction (see ACParameters class in chem/AC.h for details)");
		initialize<bool> ("pure_ae",true,"pure all electron calculation with no pseudo-atoms");
		initialize<int>  ("print_level",3,"0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");

		// Next list inferred parameters
		initialize<int> ("nalpha",-1,"number of alpha spin electrons");
		initialize<int> ("nbeta",-1,"number of beta  spin electrons");
		initialize<int> ("nmo_alpha",-1,"number of alpha spin molecular orbitals");
		initialize<int> ("nmo_beta",-1,"number of beta spin molecular orbitals");
		initialize<double> ("lo",1.e10,"smallest length scale we need to resolve");
		initialize<std::vector<double> > ("protocol",{1.e-4,1.e-6},"calculation protocol");

		// geometry optimization parameters
		// @TODO: need to be moved to molecular optimizer class
		initialize<bool> ("gopt",false,"geometry optimizer");
		initialize<double> ("gtol",1.e-4,"geometry tolerance");
		initialize<bool> ("gtest",false,"geometry tolerance");
		initialize<double> ("gval",1.e-5,"value precision");
		initialize<double> ("gprec",1.e-4,"gradient precision");
		initialize<int> ("gmaxiter",20,"optimization maxiter");
		initialize<bool> ("ginitial_hessian",false,"compute inital hessian for optimization");
		initialize<std::string> ("algopt","bfgs","algorithm used for optimization",{"bfgs","cg"});
		initialize<bool> ("tdksprop",false,"time-dependent Kohn-Sham equation propagate");
		initialize<int> ("nv_factor",1,"factor to multiply number of virtual orbitals with when automatically decreasing nvirt");
		initialize<int> ("vnucextra",2,"load balance parameter for nuclear pot");
		initialize<int> ("loadbalparts",2,"??");

		// Next list for response code from a4v4
		initialize<bool> ("response",false,"response function calculation");
		initialize<double> ("response_freq",0.0,"frequency for calculation response function");
		initialize<std::vector<bool> > ("response_axis",{true,true,true},"response axis");
		initialize<bool> ("nonrotate",false,"if true do not molcule orient (redundant with no_orient");
		initialize<double> ("rconv",1.e-6,"Response convergence");
		initialize<double> ("efield",0.0,"eps for finite field");
		initialize<int> ("efield_axis",0,"finite field axis",{0l,1,2});
//		initialize<std::map<std::string,std::string> generalkeyval;  ///< general new key/value pair

	}

	public:
	using QCCalculationParametersBase::read;


	double econv() const {return get<double>("econv");}
	double dconv() const {return get<double>("dconv");}
	bool converge_each_energy() {return get<bool>("converge_each_energy");}

	int nopen() const {return get<int>("nopen");}
	int nalpha() const {return get<int>("nalpha");}
	int nbeta() const {return get<int>("nbeta");}

	int nvalpha() const {return get<int>("nvalpha");}
	int nvbeta() const {return get<int>("nvbeta");}
	int nv_factor() const {return get<int>("nv_factor");}

	int nmo_alpha() const {return get<int>("nmo_alpha");}
	int nmo_beta() const {return get<int>("nmo_beta");}

	bool have_beta() const {return (nbeta()>0) and (not spin_restricted());}

	bool spin_restricted() const {return get<bool>("spin_restricted");}
	bool no_compute() const {return get<bool>("no_compute");}

	double lo() const {return get<double>("lo");}
	double L() const {return get<double>("l");}
	int k() const {return get<int>("k");}

	std::string localize_method() const {return get<std::string>("localize");}
	bool do_localize() const {return (localize_method()!="canon");}
	bool localize_pm() const {return (localize_method()=="pm");}

	std::string pointgroup() const {return get<std::string>("pointgroup");}
	bool do_symmetry() const {return (pointgroup()!="c1");}
	bool no_orient() const {return get<bool>("no_orient");}
	double charge() const {return get<double>("charge");}
	int print_level() const {return get<int>("print_level");}

	int maxiter() const {return get<int>("maxiter");}
	double orbitalshift() const {return get<double>("orbitalshift");}

	std::string deriv() const {return get<std::string>("deriv");}
	std::string dft_deriv() const {return get<std::string>("dft_deriv");}
	std::string pcm_data() const {return get<std::string>("pcm_data");}
	std::string ac_data() const {return get<std::string>("ac_data");}
	std::string xc() const {return get<std::string>("xc");}

	std::string aobasis() const {return get<std::string>("aobasis");}
	std::string core_type() const {return get<std::string>("core_type");}
	bool psp_calc() const {return get<bool>("psp_calc");}
	bool pure_ae() const {return get<bool>("pure_ae");}

	std::vector<double> protocol() const {return get<std::vector<double> >("protocol");}
	bool save() const {return get<bool>("save");}
	bool restart() const {return get<bool>("restart");}
	bool restartao() const {return get<bool>("restartao");}
	bool restart_cphf() const {return get<bool>("restart_cphf");}

	int maxsub() const {return get<int>("maxsub");}
	double maxrotn() const {return get<double>("maxrotn");}

	int vnucextra() const {return get<int>("vnucextra");}
	int loadbalparts() const {return get<int>("loadbalparts");}

	double response_freq() const {return get<double>("response_freq");}
	std::vector<bool> response_axis() const {return get<std::vector<bool> >("response_axis");}

	bool derivatives() const {return get<bool>("derivatives");}
	bool response() const {return get<bool>("response");}
	bool tdksprop() const {return get<bool>("tdksprop");}
	bool dipole() const {return get<bool>("dipole");}

	bool gopt() const {return get<bool>("gopt");}
	std::string algopt() const {return get<std::string>("algopt");}
	int gmaxiter() const {return get<int>("gmaxiter");}
	double gtol() const {return get<double>("gtol");}
	double gval() const {return get<double>("gval");}
	double gprec() const {return get<double>("gprec");}
	bool ginitial_hessian() const {return get<bool>("ginitial_hessian");}

	Tensor<double> plot_cell() const {
		std::vector<double> vcell=get<std::vector<double> >("plot_cell");
		if (vcell.size()==0) return Tensor<double>();
		Tensor<double> cell(3,2);
		cell(0,0)=vcell[0];
		cell(0,1)=vcell[1];
		cell(1,0)=vcell[2];
		cell(1,1)=vcell[3];
		cell(2,0)=vcell[4];
		cell(2,1)=vcell[5];
		return cell;
	}


	void set_derived_values(const Molecule& molecule, const AtomicBasisSet& aobasis) {

        for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
            if (molecule.get_pseudo_atom(iatom)){
                set_derived_value("pure_ae",false);
                continue;
            }
        }
        set_derived_value("aobasis",molecule.guess_file());;
        const int n_core = molecule.n_core_orb_all();


        std::vector<double> proto=get<std::vector<double> >("protocol");
	// No ... The accuracy of computation is INDEPENDENT of the convergence requirement
	// --- actually need more precision than convergence threshold in order to have
	// variational principle working and for robust convergence
        //proto.back()=get<double>("econv");
	set_derived_value("protocol",proto);
	// No ... the energy is variational!  Don't need more accuracy in dconv --- in fact the opposite is true.
	// set_derived_value("dconv",sqrt(get<double>("econv"))*0.1);

        double z = molecule.total_nuclear_charge();
        const double charge=get<double>("charge");
        int nelec = int(z - charge - n_core*2);
        if (fabs(nelec+charge+n_core*2-z) > 1e-6) {
            error("non-integer number of electrons?", nelec+charge+n_core*2-z);
        }

        set_derived_value("nalpha",(nelec + nopen())/2);
        set_derived_value("nbeta",(nelec - nopen())/2);

        if (nalpha() < 0) error("negative number of alpha electrons?", nalpha());
        if (nbeta() < 0) error("negative number of beta electrons?", nbeta());
        if ((nalpha()+nbeta()) != nelec) error("nalpha+nbeta != nelec", nalpha()+nbeta());
        if (nalpha() != nbeta()) set_derived_value("spin_restricted",false);

        set_derived_value("nmo_alpha",nalpha() + nvalpha());
        set_derived_value("nmo_beta",nbeta() + nvbeta());

        // Ensure we have enough basis functions to guess the requested
        // number of states ... a minimal basis for a closed-shell atom
        // might not have any functions for virtuals.
        int nbf = aobasis.nbf(molecule);
        if ((nmo_alpha()>nbf) or (nmo_beta()>nbf)) error("too few basis functions?", nbf);
//        nmo_alpha = std::min(nbf,nmo_alpha);
//        nmo_beta = std::min(nbf,nmo_beta);
//        if (nalpha>nbf || nbeta>nbf) error("too few basis functions?", nbf);
//        nvalpha = nmo_alpha - nalpha;
//        nvbeta = nmo_beta - nbeta;

        // Unless overridden by the user use a cell big enough to
        // have exp(-sqrt(2*I)*r) decay to 1e-6 with I=1ev=0.037Eh
        // --> need 50 a.u. either side of the molecule
        set_derived_value("l",molecule.bounding_cube() + 50.0);

        set_derived_value("lo",molecule.smallest_length_scale());

        // set highest possible point group for symmetry
        if (do_localize()) set_derived_value("pointgroup",std::string("c1"));
        else set_derived_value("pointgroup",molecule.pointgroup_);

        // above two lines will not override user input, so check input is sane
        if (do_localize() and do_symmetry()) {
        	error("\n\nsymmetry and localization cannot be used at the same time\n"
        			"switch from local to canonical orbitals (keyword canon)\n\n");
        }
	}

};


#else
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
    bool localize_boys;         ///< If true use boys for localization
    bool localize_new;          ///< If true use new for localization
    std::string symmetry;		///< use point group symmetry for all orbitals: default/full/schoenflies
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
    std::map<std::string,std::string> generalkeyval;  ///< general new key/value pair

    // Different derivatives can be used
    std::string deriv;       ///< Which method of derivative should be used for the KE matrix
    std::string dft_deriv;    ///< Which method of derivative should be used for dft functional

    static bool stringtobool(std::string str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str=="true" or str=="1" or str=="yes") return true;
        if (str=="false" or str=="0" or str=="no") return false;
        madness::print("unknown boolean ",str);
        return 0;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & k & L & maxrotn & nvalpha & nvbeta
        & nopen & maxiter & nio & spin_restricted;
        ar & plotlo & plothi & plotdens & plotcoul & localize & localize_pm & localize_boys & localize_new & symmetry
        & restart & restartao & save & no_compute &no_orient & maxsub & orbitalshift & npt_plot & plot_cell & aobasis;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
        ar & core_type & derivatives & conv_only_dens & dipole;
        ar & xc_data & protocol_data;
        ar & gopt & gtol & gtest & gval & gprec & gmaxiter & ginitial_hessian & algopt & tdksprop
        & nuclear_corrfac & psp_calc & print_dipole_matels & pure_ae & hessian & read_cphf & restart_cphf
        & purify_hessian & vnucextra & loadbalparts & pcm_data & ac_data & deriv & dft_deriv;
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
    , localize_pm(false)
    , localize_boys(false)
    , localize_new(true)
    , symmetry("default")
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
    , deriv("abgv")
    , dft_deriv("abgv")
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
                    std::cout << "moldft: unrecognized aobasis (sto-3g or sto-6g or 6-31g only): proceed with care " << aobasis << std::endl;
                    //MADNESS_EXCEPTION("input_error", 0);
                }
            }
            else if (s == "canon") {
                localize = false;
            }
            else if (s == "local") {
                localize = true;
            }
            else if (s == "pm") {
                localize = localize_pm = true;
                localize_boys = localize_new = false;
            }
            else if (s == "boys") {
                localize = localize_boys = true;
                localize_pm = localize_new = false;
            }
            else if (s == "newloc") {
                localize = localize_new = true;
                localize_boys = localize_pm = false;
            }
            else if (s == "symmetry") {
            	symmetry="full";
            	std::string buf,buf1;
                std::getline(f,buf);
            	std::stringstream ff(buf);
            	ff >> buf1;
                if (buf1.size()>0) symmetry=buf1;
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
            else if (s == "keyval") {
                std::string key, val;
                f >> key;
                f >> val;
                generalkeyval.insert(std::make_pair(key,val));
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
            else if (s == "deriv") {
               f >> deriv;
               if (deriv!="abgv" && deriv!="bspline" && deriv!="ble") {
                  throw "deriv must be \"abgv\", \"bspline\", or \"ble\"";
               }
            }
            else if (s == "dft_deriv") {
               f >> dft_deriv;
               if (dft_deriv!="abgv" && dft_deriv!="bspline" && dft_deriv!="ble") {
                  throw "dft_deriv must be \"abgv\", \"bspline\", or \"ble\"";
               }
            }
            else {
                std::cout << "moldft: unrecognized input keyword " << s << std::endl;
                MADNESS_EXCEPTION("input error",0);
            }
        }
        if (nopen != 0) spin_restricted = false;
        if (nvalpha || nvbeta) localize = false; // must use canonical orbitals if computing virtuals
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

        // set molecular and computational point groups
        // use highest point group unless specified by user

        // complain if symmetry has been set to anything other than c1
        if ((symmetry!="default" and symmetry!="c1") and localize) {
        	error("\n\nsymmetry and localization cannot be used at the same time\n"
        			"switch from local to canonical orbitals (keyword canon)\n\n");
        }

        // no symmetry keyword specified
    	if (symmetry=="default") {
    		if (localize) symmetry="c1";
    		else symmetry=molecule.pointgroup_;

    	// symmetry keyword specified without pointgroup
    	} else if (symmetry=="full") {
    		symmetry=molecule.pointgroup_;
    	}

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
        madness::print("  KE derivative type ", deriv);
        madness::print(" DFT derivative type ", dft_deriv);
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
        if (localize_pm) loctype = "pm";
        else if (localize_new) loctype = "new";
        else if (localize_boys) loctype = "boys";
        else loctype = "none";
        if (localize)
            madness::print("  localized orbitals ", loctype);
        else
            madness::print("  canonical orbitals ");
        madness::print("   comp. point group ", symmetry);
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
#endif

} // namespace madness

#endif /* MADNESS_CHEM_CALCULATIONPARAMETERS_H__INCLUDED */
