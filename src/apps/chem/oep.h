/*
 * oep.h
 *
 *  Created on: Nov 6, 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_OEP_H_
#define SRC_APPS_CHEM_OEP_H_



#include <chem/nemo.h>
#include <chem/cheminfo.h>
#include <chem/SCFOperators.h>
#include <chem/projector.h>

using namespace madness;
namespace madness {

struct dens_inv{

	double threshold;

	// default value for threshold is 1.0e-8
	dens_inv(const double thresh = 1.0e-8) {
		threshold = thresh;
	}

    /// @param[out] U   result
    /// @param[in]  t   numerator
    /// @param[in]  inv density to be inverted >0.0
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& inv) const {
        ITERATOR(
            U, double d = t(IND);
            double p = std::max(inv(IND), threshold);
            U(IND) = d/p;
        );
   }
    template <typename Archive>
    void serialize(Archive& ar) {}

};


//struct binary_munge{
//
//	double longrangevalue;
//	double threshold;
//
//	// default value for threshold is 1.0e-8
//	binary_munge(const double thresh = 1.0e-8, const double lrv = 0.0) {
//		longrangevalue = lrv;
//		threshold = thresh;
//	}
//
//    /// @param[out] U   result
//    /// @param[in]  f   function to be munged
//    /// @param[in]  r   reference density
//    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& f,
//            const Tensor<double>& refdens) const {
//        ITERATOR(
//            U, double r = refdens(IND);
//            double ff = f(IND);
//            U(IND) = (r > threshold) ? ff : longrangevalue;
//        );
//    }
//
//    template <typename Archive>
//    void serialize(Archive& ar) {}
//
//};


struct interpolate_munge_refdens{

	double thresh_high, thresh_low, log_high, log_low;
	typedef double resultT;

	// same default values as for dens_thresh_hi and dens_thresh_lo
	interpolate_munge_refdens(const double hi = 1.0e-4, const double lo = 1.0e-7) {
		thresh_high = hi;
		thresh_low = lo;
		log_high = log10(thresh_high);
		log_low = log10(thresh_low);
	}

    /// @param[out] U   result
    /// @param[in]  r   reference density
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& refdens) const {

    	// set 1 if above munge interval (dens > thresh_high) or 0 if below munge interval (dens < thresh_low)
    	// interpolate inbetween with logarithmic form of (r - lo)/(hi - lo) because this yields near linear behavior
    	Tensor<double> U = copy(refdens); // copy refdens in oder to have the same dimension in the tensor
    	U.fill(1.0); // start with a fuction that is 1.0 everywhere
        ITERATOR(
            U, double r = refdens(IND);
            if (r > thresh_high) {
            	U(IND) = 1.0;
            } else if (r < thresh_low) {
            	U(IND) = 0.0;
            } else {
            	U(IND) = (log10(refdens(IND)) - log_low)/(log_high - log_low);
            }
        );
        return U;

    }

    template <typename Archive>
    void serialize(Archive& ar) {}

};

class OEP_Parameters : public QCCalculationParametersBase {
public:
	OEP_Parameters(World& world, std::string inputfile) {

		initialize<std::string>("model","mrks","comment on this",{"oaep","ocep","dcep","mrks"});
		initialize<unsigned int>("maxiter",250,"maximum number of iterations in OEP algorithm");
		initialize<double>("conv_threshold",1.e-5,"comment on this");
		initialize<double>("density_threshold_high",1.e-4,"comment on this");
		initialize<double>("density_threshold_low",1.e-7,"comment on this");
		initialize<double>("density_threshold_inv",1.e-8,"comment on this");
		initialize<bool>("set_thresh_inv",false,"comment on this");
		initialize<std::vector<double> >("kain_param",{1.0e-8, 3.0},"comment on this");
		initialize<std::vector<double> >("damp_coeff",{1.0},"set coefficients for the new and a number of old potentials for damping");

//		std::vector<bool> oep_model = {false, false, false, false};
		initialize<unsigned int>("saving_amount",1,"choose level 0, 1, 2 or 3 for saving functions");
		initialize<unsigned int>("save_iter_orbs",0,"if > 0 save all orbitals every ... iterations (needs a lot of storage!");
		initialize<unsigned int>("save_iter_density",0,"if > 0 save KS density every ... iterations");
		initialize<unsigned int>("save_iter_IKS",0,"if > 0 save IKS every ... iterations");
		initialize<unsigned int>("save_iter_kin_tot_KS",0,"if > 0 save kin_tot_KS every ... iterations");
		initialize<unsigned int>("save_iter_kin_P_KS",0,"if > 0 save kin_P_KS every ... iterations");
		initialize<unsigned int>("save_iter_ocep_correction",0,"if > 0 save OCEP correction every ... iterations");
		initialize<unsigned int>("save_iter_dcep_correction",0,"if > 0 save DCEP correction every ... iterations");
		initialize<unsigned int>("save_iter_mrks_correction",0,"if > 0 save mRKS correction every ... iterations");
		initialize<unsigned int>("save_iter_total_correction",0,"if > 0 save total correction (OCEP + DCEP) every ... iterations");
		initialize<unsigned int>("save_iter_effective_potential",0,"if > 0 save effective potential every ... iterations");

		read(world,inputfile,"oep");

	}

	void set_derived_values(const Nemo::NemoCalculationParameters& nparam) {
    	set_derived_value("density_threshold_inv",0.1*get<double>("density_threshold_low"));
    	set_derived_value("conv_threshold",nparam.econv());

	}

	// convenience functions
	std::string model() const {return get<std::string>("model");}
	bool is_oaep() const {return (get<std::string>("model")=="oaep");}
	bool is_ocep() const {return (get<std::string>("model")=="ocep");}
	bool is_dcep() const {return (get<std::string>("model")=="dcep");}
	bool is_mrks() const {return (get<std::string>("model")=="mrks");}

	long damp_num() const {return get<std::vector<double> >("damp_coeff").size();}
	bool do_damping() const {return damp_num() > 1;}
	unsigned int maxiter() const {return get<unsigned int>("maxiter");}
	double conv_thresh() const {return get<double>("conv_threshold");}
	double dens_thresh_hi() const {return get<double>("density_threshold_high");}
	double dens_thresh_lo() const {return get<double>("density_threshold_low");}
	double dens_thresh_inv() const{return get<double>("density_threshold_inv");}
	unsigned int saving_amount() const {return get<unsigned int>("saving_amount");}

	unsigned int save_iter_orbs               () const {return get<unsigned int>("save_iter_orbs");}
	unsigned int save_iter_density            () const {return get<unsigned int>("save_iter_density");}
	unsigned int save_iter_IKS                () const {return get<unsigned int>("save_iter_iks");}
	unsigned int save_iter_kin_tot_KS         () const {return get<unsigned int>("save_iter_kin_tot_ks");}
	unsigned int save_iter_kin_P_KS           () const {return get<unsigned int>("save_iter_kin_p_ks");}
	unsigned int save_iter_ocep_correction    () const {return get<unsigned int>("save_iter_ocep_correction");}
	unsigned int save_iter_dcep_correction    () const {return get<unsigned int>("save_iter_dcep_correction");}
	unsigned int save_iter_mrks_correction    () const {return get<unsigned int>("save_iter_mrks_correction");}
	unsigned int save_iter_total_correction   () const {return get<unsigned int>("save_iter_total_correction");}
	unsigned int save_iter_effective_potential() const {return get<unsigned int>("save_iter_effective_potential");}

	std::vector<double> kain_param() const {return get<std::vector<double> >("kain_param");}
	std::vector<double> damp_coeff() const {return get<std::vector<double> >("damp_coeff");}

};


class OEP : public Nemo {

	typedef std::shared_ptr<real_convolution_3d> poperatorT;

private:
//
//	double dens_thresh_hi = 1.0e-4;                   // default 1.0e-4
//	double dens_thresh_lo = 1.0e-7;                   // default 1.0e-7
//	double dens_thresh_inv = 1.0e-8;                  // default 0.1*dens_thresh_lo
//	bool set_thresh_inv = false;                      // to check if default or custom
//	double conv_thresh = calc->param.econv();           // default convergence threshold is same as econv
//	std::vector<double> kain_param = {1.0e-8, 3.0};   // default KAIN settings for rcondtol and cabsmax (see nonlinsol.h)
//	unsigned int damp_num = 0;                        // default 0 (no damping)
//	std::vector<double> damp_coeff;
//	std::string model;
//	std::vector<bool> oep_model = {false, false, false, false};
//	unsigned int saving_amount = 1;                   // choose level 0, 1, 2 or 3 for saving functions
//	unsigned int save_iter_orbs = 0;                  // if > 0 save all orbitals every ... iterations (needs a lot of storage!)
//	unsigned int save_iter_density = 0;               // if > 0 save KS density every ... iterations
//	unsigned int save_iter_IKS = 0;                   // if > 0 save IKS every ... iterations
//	unsigned int save_iter_kin_tot_KS = 0;            // if > 0 save kin_tot_KS every ... iterations
//	unsigned int save_iter_kin_P_KS = 0;              // if > 0 save kin_P_KS every ... iterations
//	unsigned int save_iter_ocep_correction = 0;       // if > 0 save OCEP correction every ... iterations
//	unsigned int save_iter_dcep_correction = 0;       // if > 0 save DCEP correction every ... iterations
//	unsigned int save_iter_mrks_correction = 0;       // if > 0 save mRKS correction every ... iterations
//	unsigned int save_iter_total_correction = 0;      // if > 0 save total correction (OCEP + DCEP) every ... iterations
//	unsigned int save_iter_effective_potential = 0;   // if > 0 save effective potential every ... iterations
//
//	void set_model_oaep() {oep_model[0] = true;}
//	bool is_oaep() const {return oep_model[0];}
//	void set_model_ocep() {oep_model[1] = true;}
//	bool is_ocep() const {return oep_model[1];}
//	void set_model_dcep() {oep_model[2] = true;}
//	bool is_dcep() const {return oep_model[2];}
//	void set_model_mrks() {oep_model[3] = true;}
//	bool is_mrks() const {return oep_model[3];}

	/// returns true if all members of a vector of booleans are true, otherwise false
	bool IsAlltrue(std::vector<bool> vec) {
		for (int i = 0; i < vec.size(); i++) {
			if (!vec[i]) return false;
		}
		return true;
	}

	OEP_Parameters oep_param;

public:

	OEP(World& world, const std::shared_ptr<SCF> calc, std::string inputfile)
		: Nemo(world, calc, inputfile), oep_param(world, inputfile) {
		oep_param.set_derived_values(param);

		if (param.localize_method()!="canon") {
			MADNESS_EXCEPTION("use localized orbitals for OEP calculations",1);
		}

		oep_param.print("oep","end");

	}

//	void read_oep_param(std::istream& in) {
//        position_stream(in, "oep");
//        std::string str;
//
//        while (in >> str) {
//            if (str == "end") {
//                break;
//            }
//            else if (str == "model") {
//            	in >> model; // default model is mRKS
//            }
//            else if (str == "density_threshold_high") {
//            	in >> dens_thresh_hi;
//            }
//            else if (str == "density_threshold_low") {
//            	in >> dens_thresh_lo;
//            }
//            else if (str == "density_threshold_inverting") {
//            	in >> dens_thresh_inv;
//            	set_thresh_inv = true;
//            }
//            else if (str == "conv_threshold") {
//            	in >> conv_thresh;
//            }
//            else if (str == "kain_parameters") {
//            	in >> kain_param[0];
//            	in >> kain_param[1];
//            }
//            else if (str == "damping") {
//            	in >> damp_num;
//            	for (unsigned int i = 0; i < damp_num + 1; i++) {
//            		double coeff;
//            		in >> coeff;
//            		damp_coeff.push_back(coeff);
//            	}
//            }
//            else if (str == "saving_amount") {
//            	in >> saving_amount;
//            }
//            else if (str == "save_orbitals") {
//            	in >> save_iter_orbs;
//            }
//            else if (str == "save_density") {
//            	in >> save_iter_density;
//            }
//            else if (str == "save_IKS") {
//            	in >> save_iter_IKS;
//            }
//            else if (str == "save_kin_tot_KS") {
//            	in >> save_iter_kin_tot_KS;
//            }
//            else if (str == "save_kin_P_KS") {
//            	in >> save_iter_kin_P_KS;
//            }
//            else if (str == "save_OCEP_correction") {
//            	in >> save_iter_ocep_correction;
//            }
//            else if (str == "save_DCEP_correction") {
//            	in >> save_iter_dcep_correction;
//            }
//            else if (str == "save_mRKS_correction") {
//            	in >> save_iter_mrks_correction;
//            }
//            else if (str == "save_total_correction") {
//            	in >> save_iter_total_correction;
//            }
//            else if (str == "save_effective_potential") {
//            	in >> save_iter_effective_potential;
//            }
//            else {
//                print("oep: unrecognized input keyword:", str);
//                MADNESS_EXCEPTION("input error",0);
//            }
//        }
//
//        // set variables from input and print notes in output
//
//    	if (model == "oaep" or model == "OAEP") {
//    		set_model_oaep();
//    		model = "OAEP";
//    	} else if (model == "ocep" or model == "OCEP") {
//    		set_model_ocep();
//    		model = "OCEP";
//    	} else if (model == "dcep" or model == "DCEP") {
//    		set_model_dcep();
//    		model = "DCEP";
//    	} else if (model == "mrks" or model == "mRKS" or model == "MRKS" or model.empty()) { // default model is mRKS
//    		set_model_mrks();
//    		model = "mRKS";
//    	} else {
//            print("oep: no available OEP model selected, please choose one of the following: oaep/ocep/dcep/mrks");
//            MADNESS_EXCEPTION("input error",0);
//    	}
//
//    	if (!set_thresh_inv) dens_thresh_inv = 0.1*dens_thresh_lo; // take default dens_thresh_inv
//
//    	print("using", model, "model as approximation to OEP");
//    	print("using upper density threshold =", dens_thresh_hi);
//    	print("using lower density threshold =", dens_thresh_lo);
//    	print("using density threshold for inverting =", dens_thresh_inv);
//    	print("using convergence threshold for optimized potential =", conv_thresh);
//    	print("using KAIN parameters rcondtol =", kain_param[0], "and cabsmax =", kain_param[1]);
//    	if (damp_num == 0) {
//    		damp_coeff.push_back(1.0);
//    		print("using no damping");
//    	}
//    	else {
//    		print("using damping with", damp_num, "old potential(s) and the following coefficients:");
//    		print("         new potential =", damp_coeff[0]);
//    		for (unsigned int i = 1; i < damp_num + 1; i++) {
//    			print("  previous potential", i, "=", damp_coeff[i]);
//    		}
//    	}
//    	if (is_dcep() or is_mrks()) {
//        	if (calc->param.dft_deriv() == "bspline") print("unsing b-spline gradient operator for DCEP correction");
//        	else print("using default abgv gradient operator for", model, "correction");
//    	}
//    	print("");
//
//        // check some common mistakes in input file
//
//        if (dens_thresh_hi <= dens_thresh_lo) {
//            print("oep: density_threshold_high must always be larger than density_threshold_low!");
//            MADNESS_EXCEPTION("input error",0);
//        }
//
//        double all_coeffs = 0.0;
//        for (unsigned int i = 0; i < damp_num + 1; i++) {
//        	all_coeffs += damp_coeff[i];
//        }
//        if (all_coeffs != 1.0) {
//            print("oep: sum of damping coefficients does not equal 1.0, please check the input file!");
//            MADNESS_EXCEPTION("input error",0);
//        }
//	}

    /// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
	/// for other functionals, slater potential must be modified
	/// HF orbitals and eigenvalues are used as the guess here
	/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
	/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps
	/// converged if norm, total energy difference and orbital energy differences (if not OAEP) are converged
    void solve_oep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals);

    /// The following function tests all essential parts of the OEP program qualitatively and some also quantitatively
    void test_oep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals);

    /// get index of HOMO from a given set of orbital energies
    long homo_ind(const tensorT orbens) const {
    	long index;
    	double en_homo = orbens.max(&index);
    	return index;
    }

    /// get difference of HF and KS HOMO energies as HOMO_KS - HOMO_HF
    double homo_diff(const tensorT ev1, const tensorT ev2) const {
    	return ev1(homo_ind(ev1)) - ev2(homo_ind(ev2));
    }

    /// print orbital energies in reverse order with optional shift
    void print_orbens(const tensorT orbens, const double shift = 0.0) const {
		for (long i = orbens.size() - 1; i >= 0; i--) {
			printf(" e%2.2lu = %12.8f Eh\n", i, orbens(i) + shift);
		}
    }

    /// compute density from orbitals with ragularization (Bischoff, 2014_1, equation (19))
    real_function_3d compute_density(const vecfuncT& nemo) const {
    	real_function_3d density = 2.0*R_square*dot(world, nemo, nemo); // 2 because closed shell
    	return density;
    }

    /// get function that is 1 for density > dens_thresh_hi and 0 for density < dens_thresh_lo
    /// interpolate logarithmic density inbetween
    real_function_3d compute_weighting_function(const vecfuncT& nemo) const {
    	real_function_3d density = compute_density(nemo);
    	real_function_3d weighting = unary_op(density, interpolate_munge_refdens(
    			oep_param.dens_thresh_hi(), oep_param.dens_thresh_lo()));
    	return weighting;
    }

    /// compute \Delta rho as an indicator for the result's quality
    double compute_delta_rho(const real_function_3d rho_HF, const real_function_3d rho_KS) const {
    	// from Ospadov_2017, equation (26)
    	real_function_3d rho_diff = abs(rho_KS - rho_HF);
    	double Drho = rho_diff.trace();
    	return Drho;
    }

    /// compute Slater potential (Kohut, 2014, equation (15))
    real_function_3d compute_slater_potential(const vecfuncT& nemo, const long homo_ind) const {

        Exchange K(world, this, 0);
        vecfuncT Knemo = K(nemo);
        // 2.0*R_square in numerator and density (rho) cancel out upon division
        real_function_3d numerator = dot(world, nemo, Knemo);
        real_function_3d rho = dot(world, nemo, nemo);
//        save(numerator, "Slaterpotential_numerator");

        // dividing by rho: the minimum value for rho is dens_thresh_lo
        real_function_3d Vs = -1.0*binary_op(numerator, rho, dens_inv(oep_param.dens_thresh_inv())); // in Kohut_2014, -1.0 is included in K
        if (oep_param.saving_amount() >= 3) save(Vs, "Slaterpotential_nolra");

        // long-range asymptotic behavior for Slater potential is \int 1/|r-r'| * |phi_HOMO|^2 dr'
        // in order to compute this lra, use Coulomb potential with only HOMO density (= |phi_HOMO|^2)
        Coulomb J(world, this);
        real_function_3d lra = -1.0*J.compute_potential(R_square*square(nemo[homo_ind]));
        if (oep_param.saving_amount() >= 3) save(lra, "lra_slater");

//        // these are not yet forgotten tests about the correct long-range behavior
//        real_function_3d lra = (-0.5/nemo.size())*J.compute_potential(this);

//        for (int i = 0; i < nemo.size(); i++) {
//        	real_function_3d product = R_square*square(nemo[i]);
//        	real_function_3d potential = -1.0*J.compute_potential(R_square*square(nemo[i]));
//        	save(R*nemo[i], "phi"+stringify(i));
//        	save(product, "phi"+stringify(i)+"phi"+stringify(i));
//        	save(potential, "int_phi"+stringify(i)+"phi"+stringify(i));
//        }

        // interpolate in interval between explicit calculation and long range asymptotics
        real_function_3d weight = compute_weighting_function(nemo);
//        save(weight, "weight_slater");
        Vs = Vs*weight + lra*(1.0 - weight);

        return Vs;

    }

    /// compute average ionization energy I like Kohut, 2014, equations (21) and (25)
    real_function_3d compute_average_I(const vecfuncT& nemo, const tensorT eigvals) const {

    	// transform 1d tensor eigvals to vector epsilon
		std::vector<double> epsilon(eigvals.size());
		for (int i = 0; i < eigvals.size(); i++) epsilon[i] = eigvals(i);

		vecfuncT nemo_square = square(world, nemo); // |nemo|^2
		scale(world, nemo_square, epsilon); // epsilon*|nemo|^2
		// 2.0*R_square in numerator and density (rho) cancel out upon division
		real_function_3d numerator = sum(world, nemo_square);
        real_function_3d rho = dot(world, nemo, nemo);

        // like Kohut, 2014, equations (21) and (25)
        real_function_3d I = -1.0*binary_op(numerator, rho, dens_inv(oep_param.dens_thresh_inv()));

        // munge I for long-range asymptotic behavior which is -epsilon_HOMO
       	print("computing I: index of HOMO is", homo_ind(eigvals));
       	double lra = -1.0*eigvals(homo_ind(eigvals));
       	real_function_3d weight = compute_weighting_function(nemo);
       	I = I*weight + lra*(1.0 - weight);

        return I;

    }

    /// compute the total kinetic energy density devided by the density tau/rho with equation (6) from Kohut
    real_function_3d compute_total_kinetic_density(const vecfuncT& nemo, const tensorT eigvals) const {

    	// compute the denominator rho (density)
    	real_function_3d rho = compute_density(nemo);

    	// compute the numerator tau

	    // get \nabla R and (\nabla R)^2 via and U1 = -1/R \nabla R and U1dot = (1/R \nabla R)^2
    	const vecfuncT U1 = this->nuclear_correlation->U1vec();
	    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
	    const real_function_3d U1dot = real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

	    // get \nabla nemo
	    std::vector<vecfuncT> grad_nemo(nemo.size());
	    for (long i = 0; i < nemo.size(); i++) {
	    	if(calc->param.dft_deriv() == "bspline") grad_nemo[i] = grad_bspline_one(nemo[i]);  // gradient using b-spline
	    	else grad_nemo[i] = grad(nemo[i]);  // default gradient using abgv
	    }

	    // compute tau = 1/2 * sum |phi_i|^2
	    // = 1/2 * sum {(\nabla R)^2 * nemo_i^2 + 2 * R * nemo_i * (\nabla R) * (\nabla nemo_i)) + R^2 * (\nabla nemo_i)^2}
	    // = 1/2 * R^2 * sum {U1dot * nemo_i^2 + 2 * nemo_i * U1 * (\nabla nemo_i)) + (\nabla nemo_i)^2}
	    vecfuncT grad_nemo_squared(nemo.size());
		for (long i = 0; i < nemo.size(); i++) {
			grad_nemo_squared[i] = U1dot*square(nemo[i])
								   - 2.0*nemo[i]*dot(world, U1, grad_nemo[i])
								   + dot(world, grad_nemo[i], grad_nemo[i]);
		}
		real_function_3d tau = R_square*sum(world, grad_nemo_squared); // 1/2 * sum |Nabla nemo|^2 (* 2 because closed shell)

		// calculate quotient = tau/rho
		real_function_3d quotient = binary_op(tau, rho, dens_inv(oep_param.dens_thresh_inv()));

        // munge quotient for long-range asymptotic behavior which is -epsilon_HOMO
       	print("computing tau/rho: index of HOMO is", homo_ind(eigvals));
       	double lra = -1.0*eigvals(homo_ind(eigvals));
       	real_function_3d weight = compute_weighting_function(nemo);
       	quotient = quotient*weight + lra*(1.0 - weight);

    	return quotient;

    }

    /// compute the Pauli kinetic energy density devided by the density tau_P/rho with equation (16) from Ospadov, 2017
    real_function_3d compute_Pauli_kinetic_density(const vecfuncT& nemo, const tensorT eigvals) const {

    	// compute the denominator rho (density)
    	real_function_3d rho_square = square(dot(world, nemo, nemo)); // density squared without 2.0*R_square

    	// compute the numerator tau_P

	    // get \nabla nemo
	    std::vector<vecfuncT> grad_nemo(nemo.size());
	    for (long i = 0; i < nemo.size(); i++) {
	    	if(calc->param.dft_deriv() == "bspline") grad_nemo[i] = grad_bspline_one(nemo[i]);  // gradient using b-spline
	    	else grad_nemo[i] = grad(nemo[i]);  // default gradient using abgv
	    }

	    vecfuncT grad_nemo_term;
		for (long i = 0; i < nemo.size(); i++) {
			for (long j = i + 1; j < nemo.size(); j++) {
				vecfuncT tmp = nemo[i]*grad_nemo[j] - nemo[j]*grad_nemo[i];
				grad_nemo_term.push_back(dot(world, tmp, tmp));
			}
		}
		real_function_3d numerator = sum(world, grad_nemo_term); // numerator = tau_P * 2 * rho / R^4

		// calculate quotient = tau_P/rho
		real_function_3d quotient = 0.5*binary_op(numerator, rho_square, dens_inv(oep_param.dens_thresh_inv()));

        // munge quotient for long-range asymptotic behavior which is 0
       	print("computing tau_P/rho: index of HOMO is", homo_ind(eigvals));
       	real_function_3d weight = compute_weighting_function(nemo);
       	quotient = quotient*weight; // + lra*(1.0 - weight) but lra is 0 here

    	return quotient;

    }

    /// compute correction of the given model
    real_function_3d compute_oep_correction(const std::string model, const real_function_3d corrHF, const vecfuncT& nemoKS,
    		const tensorT eigvalsKS) const {

    	// Kohn-Sham correction of given model
    	real_function_3d corrKS;
    	if (model == "ocep") corrKS = compute_average_I(nemoKS, eigvalsKS);
    	if (model == "dcep") corrKS = compute_total_kinetic_density(nemoKS, eigvalsKS);
    	if (model == "mrks") corrKS = compute_Pauli_kinetic_density(nemoKS, eigvalsKS);

    	// calculate correction corrHF - corrKS like Kohut, 2014, equation (26) or (33)
    	real_function_3d correction = corrHF - corrKS;

    	return correction;

    }

    /// compute all potentials from given nemos except kinetic energy
    void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& Jnemo, vecfuncT& Unemo,
    		const real_function_3d V, vecfuncT& Vnemo) const {

    	// compute Coulomb part
    	compute_coulomb_potential(nemo, Jnemo);

    	// compute nuclear potential part
    	Nuclear Unuc(world, this->nuclear_correlation);
    	Unemo = Unuc(nemo);

    	// compute approximate OEP exchange potential part
    	Vnemo = V*nemo;

    }

    /// compute Coulomb potential
    void compute_coulomb_potential(const vecfuncT& nemo, vecfuncT& Jnemo) const {

    	Coulomb J = Coulomb(world, this);
    	Jnemo = J(nemo);
    	truncate(world, Jnemo);

    }

    /// compute exchange potential (needed for Econv)
    void compute_exchange_potential(const vecfuncT& nemo, vecfuncT& Knemo) const {

    	Exchange K = Exchange(world, this, 0);
    	Knemo = K(nemo);
    	truncate(world, Knemo);

    }

    /// compute kinetic energy needed for total energy
    double compute_kinetic_energy(const vecfuncT phi) const {

    	// it's ok to use phi here, no regularization necessary for this eigenvalue
    	double E_kin = 0.0;
    	for (int axis = 0; axis < 3; axis++) {
    		real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    		const vecfuncT dphi = apply(world, D, phi);
    		E_kin += 0.5 * (inner(world, dphi, dphi)).sum();
    		// -1/2 sum <Psi|Nabla^2|Psi> = 1/2 sum <NablaPsi|NablaPsi>   (integration by parts)
    	}
    	E_kin *= 2.0; // 2 because closed shell
    	return E_kin;

    }

    /// compute Evir using Levy-Perdew virial relation (Kohut_2014, (43) or Ospadov_2017, (25))
    double compute_exchange_energy_vir(const vecfuncT phi, const real_function_3d Vx) const {

    	// make vector of functions r = (x, y, z)
    	auto monomial_x = [] (const coord_3d& r) {return r[0];};
    	auto monomial_y = [] (const coord_3d& r) {return r[1];};
    	auto monomial_z = [] (const coord_3d& r) {return r[2];};
    	vecfuncT r(3);
    	r[0]=real_factory_3d(world).functor(monomial_x);
    	r[1]=real_factory_3d(world).functor(monomial_y);
    	r[2]=real_factory_3d(world).functor(monomial_z);

    	// for density note that phi is R*nemo, so no R_square is needed
    	real_function_3d rho = 2.0*dot(world, phi, phi); // 2 because closed shell
    	real_function_3d rhoterm = 3*rho + dot(world, r, grad(rho));
    	double Ex = inner(Vx, rhoterm);
    	return Ex;

    }

    /// compute exchange energy using the expectation value of the exchange operator
    double compute_exchange_energy_conv(const vecfuncT phi, const vecfuncT Kphi) const {

    	double Ex = -1.0*inner(world, phi, Kphi).sum();
    	return Ex;

    }

    /// compute all energy contributions except exchange and sum up for total energy
    /// the exchange energy must be computed priorly with the compute_exchange_energy_... methods
    double compute_energy(const vecfuncT phi, const vecfuncT Jphi, const double E_X) const {

    	// compute kinetic energy
    	double E_kin = compute_kinetic_energy(phi);

    	// compute external potential (nuclear attraction)
    	real_function_3d Vext = calc->potentialmanager->vnuclear();
    	const vecfuncT Vextphi = Vext*phi;

    	// compute remaining energies: nuclear attraction, Coulomb, nuclear repulsion
    	// computed as expectation values (see Szabo, Ostlund (3.81))
    	const double E_ext = 2.0*inner(world, phi, Vextphi).sum(); // 2 because closed shell
    	const double E_J = inner(world, phi, Jphi).sum();
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

    	return energy;

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
