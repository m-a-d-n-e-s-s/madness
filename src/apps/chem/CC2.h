/*
 * CC2.h
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#ifndef CC2_H_
#define CC2_H_

#include <chem/projector.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/CCOperators.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include <chem/CCOperators.h>

#include <examples/nonlinsol.h>
namespace madness {




class CC2 {
public:
	CC2(World &world_,const std::string &inputFileName, const Nemo &nemo_):
		world(world_),
		parameters(inputFileName, nemo_.get_calc() -> param.lo),
		nemo(nemo_),
		active_mo(nemo_.get_calc()->amo),
		//correlationfactor(world,parameters.corrfac_gamma,parameters.thresh_f12,nemo.get_calc()->molecule),
		correlationfactor(world,1.0,1.e-7,nemo.get_calc()->molecule),
		CCOPS(world,nemo,correlationfactor,parameters)
{
		output_section("CC2 Class has been initialized with the following parameters");
		// set the threshholds
		// Set Protocoll
		output("Set Protocol 3D");
		nemo_.get_calc() -> set_protocol<3>(world,parameters.thresh_3D);
		output("Set Protocol 6D");
		nemo_.get_calc() -> set_protocol<6>(world,parameters.thresh_6D);

		FunctionDefaults<3>::set_thresh(parameters.thresh_3D);
		FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
		// Make shure that k is the same in 3d and 6d functions
		FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
		// by default SCF sets the truncate_mode to 1
		FunctionDefaults<3>::set_truncate_mode(3);
        FunctionDefaults<6>::set_truncate_mode(3);
        parameters.information(world);
        parameters.sanity_check(world);
        //output_section("Testing Section in Constructor");
        //CCOPS.test_fill_tree();
}
	void plot(const real_function_3d &f, const std::string &msg = "unspecified function")const{
		plot_plane(world,f,msg);
		output("Plotted " + msg);
	}
	/// Check energy convergence: Creates the difference between two vectors and compares against given thresh in parameters
	bool check_energy_convergence(const std::vector<double> &current, const std::vector<double> &updated)const{
		if(current.size()!=updated.size())MADNESS_EXCEPTION("error in energy convergence check: different sizes in vectors",1);
		bool conv = true;
		std::vector<double> diff(current.size(),0.0);
		for(size_t i=0;i<current.size();i++){
			double diffi = updated[i] - current[i];
			diff[i] = diffi;
			if(diffi > parameters.econv) conv=false;
		}
		if(world.rank()==0){
			std::cout << "\n\n";
			std::cout << "Pair Correlation Energies: New, Old, Diff\n";
			for(size_t i=0;i<current.size();i++) std::cout << updated[i] << ", " << current[i] << ", " << diff[i] << std::endl;
			std::cout << "\n\n";
		}
		return conv;
	}
	/// make consistency tests
	bool test()const;
	/// The World
	World &world;
	/// Structure holds all the parameters used in the CC2 calculation
	const CC_Parameters parameters;
	/// The SCF Calculation
	const Nemo &nemo;
	/// Active MO
	const vecfuncT &active_mo;
	/// The electronic Correlation Factor
	CorrelationFactor correlationfactor;
	/// The CC Operator Class
	CC_Operators CCOPS;

	/// solve the CC2 ground state equations, returns the correlation energy
	double solve()const;
	bool solve_CCS()const;
	/// solve the MP2 equations (uncoupled -> Canonical Orbitals)
	double solve_uncoupled_mp2(Pairs<CC_Pair> &u)const;
	/// solve the coupled CC2 equations
	double solve_cc2(Pairs<CC_Pair> &u, CC_Singles &tau)const;
	bool iterate_cc2_singles(const Pairs<CC_Pair> &doubles, CC_Singles &singles)const;
	bool iterate_cc2_doubles( Pairs<CC_Pair> &doubles, const CC_Singles &singles)const;
	/// Compute the pair correlation energy of an electron pair function at mp2/CCD level (no singles contributions)
	double compute_mp2_pair_energy(CC_Pair &u)const;
	CC_Singles initialize_cc2_singles(const Pairs<CC_Pair> &doubles)const;
	/// Initialize an electron pair
	/// Calculate the constant Term, and the <ij|gQf|ij> etc terms
	void initialize_electron_pair(CC_Pair &u)const;
	/// Calculate the current CC2 correlation energy
	double compute_correlation_energy(const vecfuncT &singles, const Pairs<real_function_6d> &doubles)const;
	/// update the pair energies of cc2
	std::vector<double> update_cc2_pair_energies(const Pairs<CC_Pair> &doubles, const CC_Singles &singles)const;
	/// Iterates the CC2 singles equations
	void iterate_singles(vecfuncT &singles, const Pairs<real_function_6d> &doubles)const;
	/// Iterates the CC2 doubles equations
	void iterate_doubles(const vecfuncT &singles, Pairs<real_function_6d> &doubles)const;

	/// Create formated output, std output with world rank 0
	void output(const std::string &msg)const{
		if(world.rank()==0) std::cout << msg << "\n";
	}
	/// Create formated output, New programm section
	void output_section(const std::string&msg)const{
		if(world.rank()==0){
			std::cout << std::setw(100) << std::setfill('#') << std::endl;
			std::cout << "\n" << msg << "\n";
			std::cout << std::setw(100) << std::setfill('#') << "\n" << std::endl;
		}
	}
	/// Create formated output, New programm subsection
	void output_subsection(const std::string&msg)const{
		if(world.rank()==0){
			std::cout << std::setw(50) << std::setfill('*') << std::endl;
			std::cout << "\n" << msg << "\n";
			std::cout << std::setw(50) << std::setfill('*') << "\n" << std::endl;
		}
	}

};

} /* namespace madness */

#endif /* CC2_H_ */
