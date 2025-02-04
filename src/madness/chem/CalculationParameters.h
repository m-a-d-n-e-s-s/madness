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

#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/mra/QCCalculationParametersBase.h>
#include<madness/mra/commandlineparser.h>


namespace madness {

struct CalculationParameters : public QCCalculationParametersBase {

	CalculationParameters(const CalculationParameters& other) = default;

	CalculationParameters(World& world, const commandlineparser& parser) : CalculationParameters() {
		read_input_and_commandline_options(world, parser, "dft");
        // convenience option -- needs to be moved to the MolecularOptimizer class
        if (parser.key_exists("optimize")) set_user_defined_value("gopt",true);
    }

	/// ctor reading out the input file
	CalculationParameters() {

        initialize<std::string>("prefix","mad","prefixes your output/restart/json/plot/etc files");
		initialize<double>("charge",0.0,"total molecular charge");
		initialize<std::string> ("xc","hf","XC input line");
		initialize<std::string> ("hfexalg","multiworld","hf exchange algorithm: choose from multiworld (default), multiworld_row, smallmem, largemem");
		initialize<double>("smear",0.0,"smearing parameter");
		initialize<double>("econv",1.e-5,"energy convergence");
		initialize<double>("dconv",1.e-4,"density convergence");
		initialize<std::vector<std::string> >("convergence_criteria",{"bsh_residual","total_energy"},"possible values are: bsh_residual, total_energy, each_energy, density");
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
		initialize<std::string> ("pointgroup","c1","use point (sub) group symmetry if not localized",{"c1","c2","ci","cs","c2v","c2h","d2","d2h"});
		initialize<bool>  ("restart",false,"if true restart from orbitals on disk");
		initialize<bool>  ("restartao",false,"if true restart from orbitals projected into AO basis (STO3G) on disk");
		initialize<bool>  ("no_compute",false,"if true use orbitals on disk, set value to computed");
		initialize<bool>  ("save",true,"if true save orbitals to disk");
		initialize<int>   ("maxsub",10,"size of iterative subspace ... set to 0 or 1 to disable");
		initialize<double> ("orbitalshift",0.0,"scf orbital shift: shift the occ orbitals to lower energies");
		initialize<int>    ("npt_plot",101,"no. of points to use in each dim for plots");
//		initialize<Tensor<double> > ("plot_cell",Tensor<double>(),"lo hi in each dimension for plotting (default is all space)");
		initialize<std::vector<double> > ("plot_cell",std::vector<double>(),"lo hi in each dimension for plotting (default is all space)");
		initialize<std::string> ("aobasis","6-31g","AO basis used for initial guess (6-31gss, 6-31g, 3-21g, sto-6g, sto-3g)");
		initialize<bool> ("derivatives",false,"if true calculate nuclear derivatives");
		initialize<bool> ("dipole",false,"if true calculate dipole moment");
		initialize<bool> ("conv_only_dens",false,"if true remove bsh_residual from convergence criteria (deprecated)");
		initialize<bool> ("psp_calc",false,"pseudopotential calculation for all atoms");
		initialize<std::string> ("pcm_data","none","do a PCM (solvent) calculation");
		initialize<std::string> ("ac_data","none","do a calculation with asymptotic correction (see ACParameters class in chem/AC.h for details)");
		initialize<bool> ("pure_ae",true,"pure all electron calculation with no pseudo-atoms");
		initialize<int>  ("print_level",3,"0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");
		initialize<std::string>  ("molecular_structure","inputfile","where to read the molecule from: inputfile or name from the library");

		// Next list inferred parameters
		initialize<int> ("nalpha",-1,"number of alpha spin electrons");
		initialize<int> ("nbeta",-1,"number of beta  spin electrons");
		initialize<int> ("nmo_alpha",-1,"number of alpha spin molecular orbitals");
		initialize<int> ("nmo_beta",-1,"number of beta spin molecular orbitals");
		initialize<double> ("lo",1.e-10,"smallest length scale we need to resolve");
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
		initialize<int> ("nv_factor",1,"factor to multiply number of virtual orbitals with when automatically decreasing nvirt");
		initialize<int> ("vnucextra",2,"load balance parameter for nuclear pot");
		initialize<int> ("loadbalparts",2,"??");

          //Keyword to use nwchem output for initial guess
          initialize<std::string> ("nwfile","none","Base name of nwchem output files (.out and .movecs extensions) to read from");

	}

	public:
	using QCCalculationParametersBase::read_input_and_commandline_options;

    std::string prefix() const {return get<std::string>("prefix");}

	double econv() const {return get<double>("econv");}
	double dconv() const {return get<double>("dconv");}

	bool converge_density() const {
		std::vector<std::string> criteria=get<std::vector<std::string> >("convergence_criteria");
		return std::find(criteria.begin(),criteria.end(),"density")!=criteria.end();
	}
	bool converge_bsh_residual() const {
		std::vector<std::string> criteria=get<std::vector<std::string> >("convergence_criteria");
		return std::find(criteria.begin(),criteria.end(),"bsh_residual")!=criteria.end();
	}
	bool converge_total_energy() const {
		std::vector<std::string> criteria=get<std::vector<std::string> >("convergence_criteria");
		return std::find(criteria.begin(),criteria.end(),"total_energy")!=criteria.end();
	}
	bool converge_each_energy() const {
		std::vector<std::string> criteria=get<std::vector<std::string> >("convergence_criteria");
		return std::find(criteria.begin(),criteria.end(),"each_energy")!=criteria.end();
	}

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
	double charge() const {return get<double>("charge");}
	int print_level() const {return get<int>("print_level");}

	int maxiter() const {return get<int>("maxiter");}
	double orbitalshift() const {return get<double>("orbitalshift");}

	std::string deriv() const {return get<std::string>("deriv");}
	std::string dft_deriv() const {return get<std::string>("dft_deriv");}
	std::string pcm_data() const {return get<std::string>("pcm_data");}
	std::string ac_data() const {return get<std::string>("ac_data");}
	std::string xc() const {return get<std::string>("xc");}
        std::string hfexalg() const {return get<std::string>("hfexalg");}

	std::string aobasis() const {return get<std::string>("aobasis");}

	std::vector<double> protocol() const {return get<std::vector<double> >("protocol");}
	bool save() const {return get<bool>("save");}
	bool restart() const {return get<bool>("restart");}
	bool restartao() const {return get<bool>("restartao");}
	bool restart_cphf() const {return get<bool>("restart_cphf");}

	int maxsub() const {return get<int>("maxsub");}
	double maxrotn() const {return get<double>("maxrotn");}

	int vnucextra() const {return get<int>("vnucextra");}
	int loadbalparts() const {return get<int>("loadbalparts");}


	bool derivatives() const {return get<bool>("derivatives");}
	bool dipole() const {return get<bool>("dipole");}

	bool gopt() const {return get<bool>("gopt");}
	std::string algopt() const {return get<std::string>("algopt");}
	int gmaxiter() const {return get<int>("gmaxiter");}
	double gtol() const {return get<double>("gtol");}
	double gval() const {return get<double>("gval");}
	double gprec() const {return get<double>("gprec");}
	bool ginitial_hessian() const {return get<bool>("ginitial_hessian");}

     std::string nwfile() const {return get<std::string>("nwfile");}

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


	void set_derived_values(const Molecule& molecule, const AtomicBasisSet& aobasis, const commandlineparser& parser) {
        std::string inputfile=parser.value("input");
        std::string prefix=commandlineparser::remove_extension(commandlineparser::base_name(inputfile));
        if (prefix!="input") set_derived_value("prefix",prefix);

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
        else set_derived_value("pointgroup",molecule.get_pointgroup());

        // above two lines will not override user input, so check input is sane
        if (do_localize() and do_symmetry()) {
        	error("\n\nsymmetry and localization cannot be used at the same time\n"
        			"switch from local to canonical orbitals (keyword canon)\n\n");
        }

        //NWChem interface doesn't support geometry optimization
        if (get<bool>("gopt") && nwfile() != "none") error("NWchem initialization only supports single point energy calculations.");

        //NWChem only supports Boys localization (or canonical)
        if (nwfile() != "none") {
             set_derived_value("localize",std::string("boys"));
             //Error if user requested something other than Boys
             if(localize_method() != "boys" and localize_method() != "canon") error("NWchem initialization only supports Boys localization");
        }
	}
};


} // namespace madness

#endif /* MADNESS_CHEM_CALCULATIONPARAMETERS_H__INCLUDED */
