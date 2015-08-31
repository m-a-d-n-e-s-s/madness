/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#include "CC2.h"

namespace madness {

/// solve the CC2 ground state equations, returns the correlation energy
double CC2::solve()const{
	// Initialize the Pair functions (uij, i>=j)
	output_section("Initialize Electron Pairs: First guess will be the constant Term of MP2");
	CC_Timer timer_init(world,"Initialization of all pairs");
	Pairs<CC_Pair> pairs;
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=i;j<active_mo.size();j++){
			CC_Pair u(i,j);
			initialize_electron_pair(u);
			pairs.insert(i,j,u);
		}
	}
	timer_init.info();

	output_section("Solve the uncoupled MP2 equations");
	double time_mp2_start = wall_time();
	double correlation_energy_mp2 = solve_uncoupled_mp2(pairs);
	output("Solving of MP2 ended at " + stringify(wall_time()) + "\n it took " + stringify(wall_time()-time_mp2_start) + " seconds");

	output_section("Nothing more Implemented right now");
	return correlation_energy_mp2;
}

double CC2::solve_uncoupled_mp2(Pairs<CC_Pair> &pairs)const{
	// Loop over all Pair functions uij (i=<j)
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=i;j<active_mo.size();j++){
			output_subsection("Solving uncoupled MP2 equations for pair |u" + stringify(i) + stringify(j) + ">");
			CC_Timer timer_bsh_setup(world,"Setup the BSH-Operator for the pair function");
			real_convolution_6d G = BSHOperator<6>(world, sqrt(-2 * CCOPS.get_epsilon(i,j)),
							parameters.lo, parameters.thresh_bsh_6D);
			timer_bsh_setup.info();

			// Beginn the iterations
			for(size_t iter=0;iter<30;iter++){
			CC_Timer timer_iteration(world,"Iteration "+ stringify(iter));
			double current_energy=(pairs(i,j).e_singlet + pairs(i,j).e_triplet);
			double current_error=99.9;
			// Compute the non constant part of the MP2 equations which is the regularized 6D Fock Residue
			//and apply the G Operator G[(2J-K(R)+Un)|uij>]
			{
				CC_Timer timer_mp2_residue(world,"Calculate MP2 Residue (2J-K(R)+Un)|uij>");
				real_function_6d mp2_residue = CCOPS.get_MP2_potential_residue(pairs(i,j)).truncate();
				timer_mp2_residue.info();

				CC_Timer timer_apply_bsh(world,"Apply BSH Operator to MP2 Residue");
				real_function_6d Gresidue = G(mp2_residue);
				timer_apply_bsh.info();

				// Add the constant part and the residue and make the new u function
				// |u_new> = -2.0*(constant_part + Gresidue)
				CC_Timer timer_addition(world,"Add the constant_term and the MP2 Residue");
				real_function_6d unew = (pairs(i,j).constant_term + Gresidue).truncate();
				unew.scale(-2.0);
				timer_addition.info();
				// Get the error
				CC_Timer timer_make_bsh_residue(world,"Make the BSH-Residue");
				real_function_6d bsh_residue = (unew - pairs(i,j).function);
				timer_make_bsh_residue.info();
				current_error = bsh_residue.norm2();
				// update the pair function
				pairs(i,j).function = unew;
			}
			// evaluate the current mp2 energy
			double new_energy = compute_mp2_pair_energy(pairs(i,j));
			double delta = new_energy - current_energy;
			output("End of Iteration " + stringify(iter)+ "at time: " + stringify(wall_time()));
			output("Norm of BSH Residue: " + stringify(current_error));
			output("MP2 Energy: New, Old, Difference : " + stringify(new_energy) + ", " + stringify(current_energy) + ", " + stringify(delta));
			timer_iteration.info();
			current_energy = new_energy;
			}
		}
	}
}

// Unnecessary function
double CC2::compute_mp2_pair_energy(CC_Pair &u)const{
	return CCOPS.compute_mp2_pair_energy(u);
}

/// Initialize an electron pair
/// Calculate the constant Term, and the <ij|gQf|ij> etc terms
void CC2::initialize_electron_pair(CC_Pair &u)const{
	output("Initialize Electron Pair |u" + stringify(u.i) + stringify(u.j)+">");
	CC_Timer timer_integrals(world,"Make constant energy Integrals");
	double tmp = CCOPS.make_ij_gQf_ij(u.i,u.j,u);
	timer_integrals.info();
	output("\n<"+stringify(u.i)+stringify(u.j)+"|gQf|(2.0| "+stringify(u.i)+stringify(u.j) + "> - |" +stringify(u.i)+stringify(u.j)+">) = " + stringify(tmp) + "\n");

	double epsij = CCOPS.get_epsilon(u.i,u.j);
	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * epsij), parameters.lo, parameters.thresh_bsh_6D);

	CC_Timer timer_const(world,"Calculation of constant MP2 part");
	real_function_6d mp2_constant_part = CCOPS.get_MP2_potential_constant_part(u).truncate();
	timer_const.info();

	CC_Timer timer_Gconst(world,"Apply BSH to constant MP2 part");
	u.constant_term = G(-2.0*mp2_constant_part).truncate();
	if(world.rank()==0) std::cout << "||G(-2.0*QConstPart)|| = " << u.constant_term.norm2();
	timer_Gconst.info();

	// Make the first guess for the mp2 pair function
	u.function = u.constant_term;

	// Calculate the pair energy
	double test_energy = compute_mp2_pair_energy(u);
	output("Initialized Electron Pair: |u" + stringify(u.i) + stringify(u.j) + "> with pair energy: " + stringify(test_energy) + "\n");
	u.info();
}

/// Calculate the current CC2 correlation energy
double CC2::compute_correlation_energy(const vecfuncT &singles, const Pairs<real_function_6d> &doubles)const{

}
/// Iterates the CC2 singles equations
void CC2::iterate_singles(vecfuncT &singles, const Pairs<real_function_6d> &doubles)const{

}
/// Iterates the CC2 doubles equations
void CC2::iterate_doubles(const vecfuncT &singles, Pairs<real_function_6d> &doubles)const{

}

} /* namespace madness */
