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
	/// solve the MP2 equations (uncoupled -> Canonical Orbitals)
	double solve_uncoupled_mp2(Pairs<CC_Pair> &u)const;
	/// Compute the pair correlation energy of an electron pair function at mp2/CCD level (no singles contributions)
	double compute_mp2_pair_energy(CC_Pair &u)const;
	/// Initialize an electron pair
	/// Calculate the constant Term, and the <ij|gQf|ij> etc terms
	void initialize_electron_pair(CC_Pair &u)const;
	/// Calculate the current CC2 correlation energy
	double compute_correlation_energy(const vecfuncT &singles, const Pairs<real_function_6d> &doubles)const;
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
			std::cout << std::setw(100) << std::setfill('*') << std::endl;
			std::cout << "\n" << msg << "\n";
			std::cout << std::setw(100) << std::setfill('*') << "\n" << std::endl;
		}
	}
	/// Create formated output, New programm subsection
	void output_subsection(const std::string&msg)const{
		if(world.rank()==0){
			std::cout << std::setw(20) << std::setfill(' ') << std::setw(60) << std::setfill('*') << std::setw(20) << std::setfill(' ') << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "\n" << msg << "\n";
			std::cout << std::setw(20) << std::setfill(' ') << std::setw(60) << std::setfill('*') << std::setw(20) << std::setfill(' ') << std::endl;
		}
	}

};

} /* namespace madness */

#endif /* CC2_H_ */
