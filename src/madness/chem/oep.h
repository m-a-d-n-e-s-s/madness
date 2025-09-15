/*
 * oep.h
 *
 *  Created on: Nov 6, 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_OEP_H_
#define SRC_APPS_CHEM_OEP_H_



#include<madness/chem/nemo.h>
#include<madness/chem/SCFOperators.h>
#include<madness/chem/projector.h>
#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/Results.h>

namespace madness {

/// Class to compute terms of the potential
struct divide_add_interpolate {

	double thresh_high=1.e-5;
	double thresh_low=1.e-7;
	double eps_regularize=1.e-8;
	double log_high, log_low;
	bool square_denominator=false;

	divide_add_interpolate(double hi, double lo, double eps_regularize) :
		thresh_high(hi), thresh_low(lo), eps_regularize(eps_regularize),
			log_high(log10(thresh_high)), log_low(log10(thresh_low)) {}

    std::size_t get_result_size() const {
    	return 1;
    }

	/// @param[in]	t	vector containing: oep, refdensity, longrange
    /// @return		(num1/denom1 - num2/denom2) * mask + (1-mask)*longrange
	std::vector<Tensor<double> > operator()(const Key<3> & key,
			const std::vector<Tensor<double> >& t) const {

		const Tensor<double>& refdens=t[0];
		const Tensor<double>& num1=t[1];
		const Tensor<double>& denom1=t[2];
		const Tensor<double>& num2=t[3];
		const Tensor<double>& denom2=t[4];
		const Tensor<double>& longrange=t[5];

    	// set 1 if above munge interval (dens > thresh_high) or 0 if below munge interval (dens < thresh_low)
    	// interpolate inbetween with logarithmic form of (r - lo)/(hi - lo) because this yields near linear behavior
    	Tensor<double> U = copy(refdens); // copy refdens in oder to have the same dimension in the tensor
    	U.fill(1.0); // start with a fuction that is 1.0 everywhere
        ITERATOR(
            U,
			double r = refdens(IND);
        	double result=num1(IND)/(denom1(IND)+eps_regularize)
        			- num2(IND)/(denom2(IND)+eps_regularize);
        	if (square_denominator) result=num1(IND)/(denom1(IND)*denom1(IND)+eps_regularize)
        			- num2(IND)/(denom2(IND)*denom2(IND)+eps_regularize);
            if (r > thresh_high) {
            	U(IND) = result;
            } else if (r < thresh_low) {
            	U(IND) = longrange(IND);
            } else {
            	// mask = 1 if refdens>hi, 0 if refdens < lo, interpolate in between
            	double mask=(log10(refdens(IND)) - log_low)/(log_high - log_low);
            	U(IND)=mask*result + (1.0-mask)*longrange(IND);
            }
        );
        return std::vector<Tensor<double> > (1,U);
	}
};

/// TODO:
///  - change algorithm to handle localized orbitals
///  - do not recompute the OEP potentials that don't depend on the orbital energies/Fock matrix elements
///  - do not recompute the HF contributions to the OEP potentials
///  - think about the long-range part of the Slater potential (or medium-range)
class OEP_Parameters : public QCCalculationParametersBase {
public:
	static constexpr char const *tag = "oep";

    OEP_Parameters() {
        initialize<std::vector<std::string> >("model",{"dcep"},"comment on this: oaep ocep dcep mrks");
        initialize<unsigned int>("maxiter",150,"maximum number of iterations in OEP algorithm");
        initialize<bool>("restart",false,"restart from previous OEP calculation");
        initialize<bool>("no_compute",false,"read from previous OEP calculation, no computation");
        initialize<double>("levelshift",0.0,"shift occupied orbital energies in the BSH operator");
//		initialize<double>("conv_threshold",1.e-5,"comment on this");
        initialize<double>("density_threshold_high",1.e-6,"comment on this");
        initialize<double>("density_threshold_low",1.e-8,"comment on this");
        initialize<double>("density_threshold_inv",1.e-9,"comment on this");
        initialize<std::vector<double> >("kain_param",{1.0e-8, 3.0},"comment on this");

//		std::vector<bool> oep_model = {false, false, false, false};
        initialize<unsigned int>("saving_amount",0,"choose level 0, 1, 2 or 3 for saving functions");
        initialize<unsigned int>("save_iter_orbs",0,"if > 0 save all orbitals every ... iterations (needs a lot of storage!");
        initialize<unsigned int>("save_iter_density",0,"if > 0 save KS density every ... iterations");
        initialize<unsigned int>("save_iter_IKS",0,"if > 0 save IKS every ... iterations");
        initialize<unsigned int>("save_iter_kin_tot_KS",0,"if > 0 save kin_tot_KS every ... iterations");
        initialize<unsigned int>("save_iter_kin_P_KS",0,"if > 0 save kin_P_KS every ... iterations");
        initialize<unsigned int>("save_iter_corrections",0,"if > 0 save OEP correction(s) every ... iterations");
        initialize<unsigned int>("save_iter_effective_potential",0,"if > 0 save effective potential every ... iterations");
    }

	OEP_Parameters(World& world, const commandlineparser& parser) : OEP_Parameters() {
        read_input_and_commandline_options(world, parser, "oep");
	}

    OEP_Parameters(const OEP_Parameters& other) = default;


	std::string get_tag() const override {
		return std::string("oep");
	}


    void set_derived_values(const CalculationParameters& cparam) {
    	set_derived_value("density_threshold_high",10.0*cparam.econv());
    	set_derived_value("density_threshold_low",0.01*get<double>("density_threshold_high"));
		if (dens_thresh_hi()<dens_thresh_lo()) {
			MADNESS_EXCEPTION("confused thresholds for the long-range transition",1);
		}
    	set_derived_value("density_threshold_inv",0.1*get<double>("density_threshold_low"));

    	if (levelshift()>0.0) {
    		print("levelshift > 0.0 in oep parameters\n\n");
    		MADNESS_EXCEPTION("levelshift > 0.0 in oep parameters",1);
    	}
	}

	// convenience functions
	std::vector<std::string> model() const {return get<std::vector<std::string> >("model");}

	bool restart() const {return get<bool>("restart");}
    bool no_compute() const {return get<bool>("no_compute");}
	unsigned int maxiter() const {return get<unsigned int>("maxiter");}
	double levelshift() const {return get<double>("levelshift");}
//	double conv_thresh() const {return get<double>("conv_threshold");}
	double dens_thresh_hi() const {return get<double>("density_threshold_high");}
	double dens_thresh_lo() const {return get<double>("density_threshold_low");}
	double dens_thresh_inv() const{return get<double>("density_threshold_inv");}
	unsigned int saving_amount() const {return get<unsigned int>("saving_amount");}
	unsigned int save_iter_corrections        () const {return get<unsigned int>("save_iter_corrections");}

	std::vector<double> kain_param() const {return get<std::vector<double> >("kain_param");}

};


class OEP : public Nemo {


protected:
    /// parameters for this OEP calculation
	OEP_Parameters oep_param;

private:
	/// the wave function reference that determines the local potential
	std::shared_ptr<const Nemo> reference;

	/// the final local potential
	real_function_3d Vfinal;

	/// collection of results
	std::vector<OEPResults> results;

public:

	OEP(World& world, const OEP_Parameters& oepparam, const std::shared_ptr<const Nemo>& reference)
			: Nemo(world, reference->get_calc_param(), reference->get_nemo_param(), reference->molecule()),		// OEP is a nemo calculation, ctor will set param from reference
			  oep_param(oepparam),
			  reference(reference) {

		// add tight convergence criteria
		std::vector<std::string> convergence_crit=get_calc_param().get<std::vector<std::string> >("convergence_criteria");
		if (std::find(convergence_crit.begin(),convergence_crit.end(),"each_energy")==convergence_crit.end()) {
			convergence_crit.push_back("each_energy");
		}
		calc->param.set_derived_value("convergence_criteria",convergence_crit);

		oep_param.set_derived_values(get_calc_param());

		// check that reference is well converged
		MADNESS_CHECK_THROW(reference->get_calc_param().converge_each_energy(),"need tightly converged reference for OEP calculation");

	}


    OEP(World& world, const commandlineparser& parser)
            : Nemo(world, parser),
              oep_param(world, parser) {

        // add tight convergence criteria
        std::vector<std::string> convergence_crit=get_calc_param().get<std::vector<std::string> >("convergence_criteria");
        if (std::find(convergence_crit.begin(),convergence_crit.end(),"each_energy")==convergence_crit.end()) {
            convergence_crit.push_back("each_energy");
        }
        calc->param.set_derived_value("convergence_criteria",convergence_crit);

        // set reference
		auto ref=std::make_shared<Nemo>(world,parser);
        ref->get_calc()->param.set_derived_value("convergence_criteria",convergence_crit);
        set_reference(ref);

        oep_param.set_derived_values(get_calc_param());
    }

    std::string name() const override {return "oep";}

    static void help() {
        print_header2("help page for oep");
        print("The oep code computes local exchange potentials based on a Hartree-Fock calculation from nemo");
        print("oep --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("oep --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        OEP_Parameters param;
        print("default parameters for the oep program are");
        param.print("oep", "end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    void set_reference(const std::shared_ptr<Nemo> reference1) {
	    reference=reference1;
	}

    std::shared_ptr<const Nemo> get_reference() const {
        return reference;
    }

    void print_parameters(std::vector<std::string> what) const {
        for (auto w : what) {
            if (w=="oep") oep_param.print("oep");
            else if (w=="reference") reference->get_calc_param().print("dft");
            else if (w=="oep_calc") get_calc_param().print("oep_calc");
            else {MADNESS_EXCEPTION(std::string("unknown parameter set to print "+w).c_str(),1);}
        }
    }

	real_function_3d get_final_potential() const {
        return Vfinal;
    }

    double value() override {
        return value(calc->molecule.get_all_coords());
    }

    /// update the json file with calculation input and output
    void output_calc_info_schema(const double& energy) const;


    virtual double value(const Tensor<double>& x) override {
    	MADNESS_CHECK_THROW(reference->check_converged(x),"OEP reference is not converged at this geometry");
        set_protocol(get_calc_param().econv());
        calc->copy_data(world,*(reference->get_calc()));
        double energy=0.0;
        Tensor<double> fock;
        bool load_mos=(oep_param.restart() or oep_param.no_compute());
        if (load_mos) load_restartdata(fock);

        if (not oep_param.no_compute())  energy=solve(reference->get_calc()->get_amo());
        output_calc_info_schema(energy);
        return energy;
	};

    /// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
	/// for other functionals, slater potential must be modified
	/// HF orbitals and eigenvalues are used as the guess here
	/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
	/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps
	/// converged if norm, total energy difference and orbital energy differences (if not OAEP) are converged
    double solve(const vecfuncT& HF_nemo);

	/// results are computed in compute_and_print_final_energies
    nlohmann::json analyze() const override {
    	// turn Results vector into a json object
    	nlohmann::json results_json;
    	for (const auto& r : results) {
			results_json.push_back(r.to_json());
		}
    	return  results_json;
    };

    OEPResults iterate(const std::string model, const vecfuncT& HF_nemo, const tensorT& HF_eigvals,
    		vecfuncT& KS_nemo, tensorT& KS_Fock, real_function_3d& Voep,
			const real_function_3d Vs) const;

    std::shared_ptr<Fock<double,3>> make_fock_operator() const override;

//    MolecularOrbitals<double,3> to_MO() const {
//    	std::vector<std::string> str_irreps;
//    	vecfuncT aaa=symmetry_projector(calc->amo,R_square,str_irreps);
//    	return MolecularOrbitals<double,3>(aaa,this->get_calc()->aeps,str_irreps,this->get_calc()->aocc,this->get_calc()->aset);
//    }

    void save_restartdata(const Tensor<double>& fock) const;

    void load_restartdata(Tensor<double>& fock);

    std::tuple<Tensor<double>, vecfuncT> recompute_HF(const vecfuncT& HF_nemo) const;

    /// The following function tests all essential parts of the OEP program qualitatively and some also quantitatively
    bool selftest() override;

    bool need_ocep_correction(const std::string& model) const {
    	return (model=="ocep") or (model=="dcep") or (model=="mrks");
    }

    bool need_dcep_correction(const std::string& model) const {
    	return (model=="dcep");
    }

    bool need_mrks_correction(const std::string& model) const {
    	return (model=="mrks");
    }

    /// print orbital energies in reverse order with optional shift
    void print_orbens(const tensorT orbens) const {
		for (long i = orbens.size() - 1; i >= 0; i--) {
			printf(" e%2.2lu = %12.8f Eh\n", i, orbens(i));
		}
    }

    OEPResults compute_and_print_final_energies(const std::string model, const real_function_3d& Voep,
    		const vecfuncT& KS_nemo, const tensorT& KS_Fock,
			const vecfuncT& HF_nemo, const tensorT& HF_Fock) const;

    /// compute density from orbitals with ragularization (Bischoff, 2014_1, equation (19))
    real_function_3d compute_density(const vecfuncT& nemo) const {
    	real_function_3d density = 2.0*R_square*dot(world, nemo, nemo); // 2 because closed shell
    	return density;
    }

    /// compute Delta rho as an indicator for the result's quality
    double compute_delta_rho(const real_function_3d rho_HF, const real_function_3d rho_KS) const {
    	// from Ospadov_2017, equation (26)
    	real_function_3d rho_diff = abs(rho_KS - rho_HF);
    	double Drho = rho_diff.trace();
    	return Drho;
    }

     /// compute Slater potential (Kohut, 2014, equation (15))
    real_function_3d compute_slater_potential(const vecfuncT& nemo) const {

        Exchange<double,3> K(world,reference->get_calc()->param.lo());
         K.set_bra_and_ket(R_square * nemo, nemo);
        const vecfuncT Knemo = K(nemo);
        // 2.0*R_square in numerator and density (rho) cancel out upon division
        real_function_3d numerator = -1.0*dot(world, nemo, Knemo);
        real_function_3d rho = dot(world, nemo, nemo);

        // long-range asymptotic behavior for Slater potential is \int 1/|r-r'| * |phi_HOMO|^2 dr'
        // in order to compute this lra, use Coulomb potential with only HOMO density (= |phi_HOMO|^2)
//        Coulomb J(world);
//        J.reset_poisson_operator_ptr(param.lo(),param.econv());
//        real_function_3d lra = -1.0*J.compute_potential(R_square*square(nemo[homo_ind]));
        Coulomb<double,3> J(world,this);
        real_function_3d lra=-1.0/(get_calc_param().nalpha()+get_calc_param().nbeta())*J.compute_potential(this);
//        print("compute long-range part of the Slater potential from the full molecular density");
        if (oep_param.saving_amount() >= 3) save(lra, "lra_slater");

        real_function_3d zero=real_factory_3d(world);
        real_function_3d one=real_factory_3d(world).functor([](const coord_3d& r){return 1.0;});
    	std::vector<real_function_3d> args={R_square*rho,numerator,rho,zero,one,lra};
        refine_to_common_level(world,args);

        divide_add_interpolate op(oep_param.dens_thresh_hi(), oep_param.dens_thresh_lo(),
        		oep_param.dens_thresh_inv());
        real_function_3d VSlater=multi_to_multi_op_values(op,args)[0];
        return VSlater;
    }


    /// compute the total kinetic energy density of equation (6) from Kohut

    /// return without the NCF and factor 2 for closed shell !
    real_function_3d compute_total_kinetic_density(const vecfuncT& nemo) const {

	    // compute tau = 1/2 * sum |\nabla phi_i|^2
	    // = 1/2 * sum {(\nabla R)^2 * nemo_i^2 + 2 * R * nemo_i * (\nabla R) * (\nabla nemo_i)) + R^2 * (\nabla nemo_i)^2}
	    // = 1/2 * R^2 * sum {U1dot * nemo_i^2 + 2 * nemo_i * U1 * (\nabla nemo_i)) + (\nabla nemo_i)^2}

    	// get \nabla R and (\nabla R)^2 via and U1 = -1/R \nabla R and U1dot = (1/R \nabla R)^2
    	const vecfuncT U1 = this->ncf->U1vec();
	    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(ncf.get());
	    const real_function_3d U1dot = real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

	    real_function_3d u1u1term=U1dot*dot(world,nemo,nemo);
	    real_function_3d u1dnterm=real_factory_3d(world).compressed();
	    real_function_3d dndnterm=real_factory_3d(world).compressed();

	    for (int idim=0; idim<3; ++idim) {
	    	real_derivative_3d D(world,idim);
	    	if(get_calc_param().dft_deriv() == "bspline") D.set_bspline1();
	    	vecfuncT nemo_copy=copy(world,nemo);
	    	refine(world,nemo_copy);
	    	std::vector<real_function_3d> dnemo=apply(world,D,nemo_copy);
	    	u1dnterm+=U1[idim]*dot(world,nemo_copy,dnemo);
	    	refine(world,dnemo);
	    	dndnterm+=dot(world,dnemo,dnemo);
	    }

	    real_function_3d tau=0.5*(u1u1term - 2.0*u1dnterm + dndnterm);

		return tau;

    }

    /// compute the Pauli kinetic energy density divided by the density tau_P/rho with equation (16) from Ospadov, 2017
    real_function_3d compute_Pauli_kinetic_density(const vecfuncT& nemo) const {

    	// compute the numerator tau_P

	    // get \nabla nemo
	    std::vector<vecfuncT> grad_nemo(nemo.size());
	    for (size_t i = 0; i < nemo.size(); i++) {
	    	vecfuncT nemo_copy=copy(world,nemo);
	    	refine(world,nemo_copy);

	    	if(get_calc_param().dft_deriv() == "bspline") grad_nemo[i] = grad_bspline_one(nemo_copy[i]);  // gradient using b-spline
	    	else grad_nemo[i] = grad(nemo_copy[i]);  // default gradient using abgv
	    }

	    vecfuncT grad_nemo_term;
		for (size_t i = 0; i < nemo.size(); i++) {
			for (size_t j = i + 1; j < nemo.size(); j++) {
				vecfuncT tmp = nemo[i]*grad_nemo[j] - nemo[j]*grad_nemo[i];
				grad_nemo_term.push_back(dot(world, tmp, tmp));
			}
		}
		real_function_3d numerator = sum(world, grad_nemo_term); // numerator = tau_P * 2 * rho / R^4

		return numerator;

    }

    /// compute correction of the given model

    /// shift of the diagonal elements of the fock matrix results in a global shift
    /// of the potential
    /// \f[
    /// \frac{\bar \epsilon}{\rho} = \frac{1}{\rho}\sum_{ij}\phi_i(F_ij+\delta_ij s)\phi_j
    ///        = \frac{1}{\rho} ( \sum_{ij}\phi_i F_ij \phi_j + s\sum_i\phi_i\phi_i )
    ///        = s + \frac{1}{\rho} \sum_{ij}\phi_i F_ij \phi_j
    /// \f]
    real_function_3d compute_ocep_correction(const real_function_3d& ocep_numerator_HF,
    		const vecfuncT& nemoHF, const vecfuncT& nemoKS,
			const tensorT& fockHF, const tensorT& fockKS) const {

    	if (fockKS.normf()<1.e-10) {
    		real_function_3d zero=real_factory_3d(world);
    		return zero;
    	}

    	// compute HOMO energies and the long-range asymptotics
    	auto [eval, evec] = syev(fockKS);
       	double homoKS = -eval.max();

       	auto [eval1, evec1] = syev(fockHF);
       	double homoHF = -eval1.max();

        double longrange=homoHF-homoKS;

        tensorT fock1=copy(fockKS);
        for (int i=0; i<fock1.dim(0); ++i) fock1(i,i)-=longrange;

		// 2.0*R_square in numerator and density (rho) cancel out upon division
//    	real_function_3d numeratorHF=-1.0*compute_energy_weighted_density_local(nemoHF,fockHF);
//    	real_function_3d numeratorHF=-1.0*compute_energy_weighted_density(nemoHF,eigvalsHF);
    	real_function_3d numeratorKS=-1.0*compute_energy_weighted_density_local(nemoKS,fock1);
//    	real_function_3d numeratorKS=-1.0*compute_energy_weighted_density(nemoKS,eval);

		// 2.0*R_square in numerator and density (rho) cancel out upon division
        real_function_3d densityKS = dot(world, nemoKS, nemoKS);
        real_function_3d densityHF = dot(world, nemoHF, nemoHF);

    	real_function_3d lra_func=real_factory_3d(world);
    	std::vector<real_function_3d> args={densityKS,ocep_numerator_HF,densityHF,
    			numeratorKS,densityKS,lra_func};
        refine_to_common_level(world,args);

        divide_add_interpolate op(oep_param.dens_thresh_hi(), oep_param.dens_thresh_lo(),
        		oep_param.dens_thresh_inv());
        real_function_3d correction=multi_to_multi_op_values(op,args)[0];


    	return correction;

    }

    /// return without the NCF and factor 2 for closed shell !

    /// follow Eq. (4) of Ryabinkin2014
    real_function_3d compute_energy_weighted_density_local(const vecfuncT& nemo,
    		const tensorT& fock) const {
    	return dot(world,nemo,transform(world,nemo,fock));
    }

    /// compute correction of the given model
    real_function_3d compute_dcep_correction(const real_function_3d& dcep_numerator_HF,
    		const vecfuncT& nemoHF, const vecfuncT& nemoKS) const {

		// 2.0*R_square in numerator and density (rho) cancel out upon division
//    	real_function_3d numeratorHF=compute_total_kinetic_density(nemoHF);
    	real_function_3d numeratorKS=compute_total_kinetic_density(nemoKS);

		// 2.0*R_square in numerator and density (rho) cancel out upon division
        real_function_3d densityKS = dot(world, nemoKS, nemoKS);
        real_function_3d densityHF = dot(world, nemoHF, nemoHF);

    	real_function_3d lra_func=real_factory_3d(world).functor([](const coord_3d& r)
    			{return 0.0;});

    	std::vector<real_function_3d> args={densityKS,dcep_numerator_HF,densityHF,
    			numeratorKS,densityKS,lra_func};
        refine_to_common_level(world,args);

        divide_add_interpolate op(oep_param.dens_thresh_hi(), oep_param.dens_thresh_lo(),
        		oep_param.dens_thresh_inv());
        real_function_3d correction=multi_to_multi_op_values(op,args)[0];

    	return correction;
    }

    /// compute correction of the given model
    real_function_3d compute_mrks_correction(const real_function_3d& mrks_numerator_HF,
    		const vecfuncT& nemoHF, const vecfuncT& nemoKS) const {

		// 2.0*R_square in numerator and density (rho) cancel out upon division
//    	real_function_3d numeratorHF=compute_Pauli_kinetic_density(nemoHF);
    	real_function_3d numeratorKS=compute_Pauli_kinetic_density(nemoKS);

		// 2.0*R_square in numerator and density (rho) cancel out upon division
        real_function_3d densityKS = dot(world, nemoKS, nemoKS);
        real_function_3d densityHF = dot(world, nemoHF, nemoHF);

        // longrange correction is zero
    	real_function_3d lra_func=real_factory_3d(world).functor([](const coord_3d& r) {return 0.0;});

    	real_function_3d denssq=square(densityKS);
    	std::vector<real_function_3d> args={densityKS,mrks_numerator_HF,densityHF,
    			numeratorKS,densityKS,lra_func};
        refine_to_common_level(world,args);

        divide_add_interpolate op(oep_param.dens_thresh_hi(), oep_param.dens_thresh_lo(),
        		oep_param.dens_thresh_inv());
        op.square_denominator=true;
        real_function_3d correction=0.5*multi_to_multi_op_values(op,args)[0];

    	return correction;
    }

    /// compute all potentials from given nemos except kinetic energy
    void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& Jnemo, vecfuncT& Unemo) const {

    	// compute Coulomb part
    	compute_coulomb_potential(nemo, Jnemo);

    	// compute nuclear potential part
    	Nuclear<double,3> Unuc(world, this->ncf);
    	Unemo = Unuc(nemo);

    }

    /// compute Coulomb potential
    void compute_coulomb_potential(const vecfuncT& nemo, vecfuncT& Jnemo) const {

    	Coulomb<double,3> J(world, this);
    	real_function_3d density=this->compute_density(nemo);
    	J.potential()=J.compute_potential(density);
    	Jnemo = J(nemo);
    	truncate(world, Jnemo);

    }

    /// compute exchange potential (needed for Econv)
    void compute_exchange_potential(const vecfuncT& nemo, vecfuncT& Knemo) const {

    	Exchange<double,3> K(world,this->get_calc()->param.lo());
        K.set_bra_and_ket(R_square * nemo, nemo);
    	Knemo = K(nemo);
    	truncate(world, Knemo);

    }

    /// compute Evir using Levy-Perdew virial relation (Kohut_2014, (43) or Ospadov_2017, (25))
    double compute_exchange_energy_vir(const vecfuncT& nemo, const real_function_3d Vx) const {

    	// make vector of functions r = (x, y, z)
    	auto monomial_x = [] (const coord_3d& r) {return r[0];};
    	auto monomial_y = [] (const coord_3d& r) {return r[1];};
    	auto monomial_z = [] (const coord_3d& r) {return r[2];};
    	vecfuncT r(3);
    	r[0]=real_factory_3d(world).functor(monomial_x);
    	r[1]=real_factory_3d(world).functor(monomial_y);
    	r[2]=real_factory_3d(world).functor(monomial_z);

    	real_function_3d rhonemo = 2.0*dot(world, nemo, nemo); // 2 because closed shell

    	real_function_3d rhoterm = 3*rhonemo + dot(world, r, -2.0*ncf->U1vec()*rhonemo + grad(rhonemo,true));
    	double Ex = inner(Vx, R_square*rhoterm);
    	return Ex;

    }

    /// compute exchange energy using the expectation value of the exchange operator
    double compute_exchange_energy_conv(const vecfuncT phi, const vecfuncT Kphi) const {

    	double Ex = -1.0*inner(world, phi, Kphi).sum();
    	return Ex;

    }

    /// compute all energy contributions except exchange and sum up for total energy
    /// the exchange energy must be computed priorly with the compute_exchange_energy_... methods
    std::vector<double> compute_energy(const vecfuncT& nemo, const double E_X) const {

    	// compute kinetic energy
    	double E_kin = compute_kinetic_energy(nemo);

    	real_function_3d density=compute_density(nemo);

    	// compute external potential (nuclear attraction)
    	real_function_3d Vext = calc->potentialmanager->vnuclear();
    	Coulomb<double,3> J(world,this);
    	real_function_3d Jpotential=J.compute_potential(density);

    	// compute remaining energies: nuclear attraction, Coulomb, nuclear repulsion
    	// computed as expectation values (see Szabo, Ostlund (3.81))
    	const double E_ext = inner(Vext,density); // 2 because closed shell
    	const double E_J = 0.5*inner(density,Jpotential);
    	const double E_nuc = calc->molecule.nuclear_repulsion_energy();
    	double energy = E_kin + E_ext + E_J + E_X + E_nuc;

    	if (world.rank() == 0) {
    		printf("\n                      kinetic energy %15.8f Eh\n", E_kin);
    		printf("  electron-nuclear attraction energy %15.8f Eh\n", E_ext);
    		printf("                      Coulomb energy %15.8f Eh\n", E_J);
    		printf("                     exchange energy %15.8f Eh\n", E_X);
    		printf("    nuclear-nuclear repulsion energy %15.8f Eh\n", E_nuc);
    		printf("                        total energy %15.8f Eh\n\n", energy);
    	}
    	return {energy,E_kin,E_ext,E_J,E_X};

    }

    /// compute diagonal elements of Fock matrix
    Tensor<double> compute_fock_diagonal_elements(const Tensor<double>& KS_eigvals,
    		const vecfuncT& phi, const vecfuncT& Kphi, const real_function_3d& Vx) const {
    	return KS_eigvals - inner(world, phi, Kphi) - inner(world, phi, Vx*phi);
    }

    /// cumpute E^(0) = \sum_i \epsilon_i^KS
    double compute_E_zeroth(const tensorT eigvals) const {
    	double E_0 = 0.0;
    	for (long i = 0; i < eigvals.size(); i++) {
    		E_0 += eigvals(i);
    	}
    	E_0 *= 2.0; // closed shell: every orbital energy must be counted twice
    	return E_0;
    }

    /// compute E^(1) = 1/2*\sum_ij <ij||ij> - \sum_i <i|J + Vx|i> = \sum_i <i|- 0.5*J - 0.5*K - Vx|i>
    double compute_E_first(const vecfuncT phi, const vecfuncT Jphi, const vecfuncT Kphi, const real_function_3d Vx) const {

    	//compute expectation values:
    	const double E_J = inner(world, phi, Jphi).sum();
    	const double E_K = inner(world, phi, Kphi).sum();
    	const double E_Vx = inner(world, phi, Vx*phi).sum();
    	printf("  E_J   =  %15.8f Eh\n", E_J);
    	printf("  E_K   =  %15.8f Eh\n", E_K);
    	printf("  E_Vx  =  %15.8f Eh\n", E_Vx);
    	const double E_nuc = calc->molecule.nuclear_repulsion_energy();

    	double E_1 = -1.0*(E_J + E_K + 2.0*E_Vx) + E_nuc;
    	return E_1;

    }

};



} /* namespace madness */

#endif /* SRC_APPS_CHEM_OEP_H_ */
