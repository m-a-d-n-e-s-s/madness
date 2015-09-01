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
	std::vector<double> pair_energies;
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=i;j<active_mo.size();j++){

			output_subsection("Solving uncoupled MP2 equations for pair |u" + stringify(i) + stringify(j) + ">");

			NonlinearSolverND<6> solver(parameters.kain_subspace);
			solver.do_print = (world.rank() == 0);

			output_subsection("Setup the BSH Operator");
			CC_Timer timer_bsh_setup(world,"Setup the BSH-Operator for the pair function");
			real_convolution_6d G = BSHOperator<6>(world, sqrt(-2 * CCOPS.get_epsilon(i,j)),
							parameters.lo, parameters.thresh_bsh_6D);
			G.destructive_ = true;
			output("Constructed Green Operator is destructive ? : " + stringify(G.destructive()));
			timer_bsh_setup.info();

			// Beginn the iterations
			for(size_t iter=0;iter<30;iter++){
			CC_Timer timer_iteration(world,"\n\nIteration "+ stringify(iter)+"\n\n");
			double current_energy=(pairs(i,j).e_singlet + pairs(i,j).e_triplet);
			double current_error=99.9;
			// Compute the non constant part of the MP2 equations which is the regularized 6D Fock Residue
			//and apply the G Operator G[(2J-K(R)+Un)|uij>]
			{
				output_subsection("Calculate MP2 Residue");
				CC_Timer timer_mp2_residue(world,"\n\nCalculate MP2 Residue (2J-K(R)+Un)|uij>\n\n");
				real_function_6d mp2_residue = CCOPS.get_MP2_potential_residue(pairs(i,j)).truncate();
				mp2_residue.print_size("Vpsi");
				timer_mp2_residue.info();

				output_subsection("Apply the Green's Operator");
				CC_Timer timer_apply_bsh(world,"\n\nApply BSH Operator to MP2 Residue\n\n");
				mp2_residue.scale(-2.0);
				real_function_6d Gresidue = G(mp2_residue);
				Gresidue.print_size("G(J+U-K)|u>");
				timer_apply_bsh.info();

				// Add the constant part and the residue and make the new u function
				// |u_new> = (constant_part + Gresidue)
				output_subsection("Add the Constant Term");
				CC_Timer timer_addition(world,"\n\nAdd the constant_term and the MP2 Residue\n\n");
				real_function_6d unew = (pairs(i,j).constant_term + Gresidue).truncate();
				unew.print_size("unew");
				unew = CCOPS.Q12(unew);
				unew.print_size("Q12(unew)");
				timer_addition.info();
				// Get the error
				CC_Timer timer_make_bsh_residue(world,"\n\nMake the BSH-Residue\n\n");
				real_function_6d bsh_residue = (pairs(i,j).function-unew);
				timer_make_bsh_residue.info();
				current_error = bsh_residue.norm2();
				// update the pair function
				if(parameters.kain) pairs(i,j).function = CCOPS.Q12(solver.update(unew, bsh_residue));
				else pairs(i,j).function = unew;
			}
			// evaluate the current mp2 energy
			double new_energy = compute_mp2_pair_energy(pairs(i,j));
			double delta = new_energy - current_energy;
			output("End of Iteration " + stringify(iter)+ "at time: " + stringify(wall_time()));
			output("Norm of BSH Residue: " + stringify(current_error));
			output("MP2 Energy: New, Old, Difference : " + stringify(new_energy) + ", " + stringify(current_energy) + ", " + stringify(delta));

			output_subsection("End of Iteration " + stringify(iter));
			if(world.rank()==0){
				std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) <<    "current correlation energy:"  << new_energy << std::endl;
				std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) <<    "previous correlation energy:"  << current_energy << std::endl;
				std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) <<    "correlation energy difference:"  << delta << std::endl;
				std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) <<    "current wavefunction error:"  << current_error << std::endl;
				std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) <<    "wavefunction norm:"  << pairs(i,j).function.norm2() << std::endl;
			}

			timer_iteration.info();
			current_energy = new_energy;
			if(current_error < parameters.dconv_6D){
				output("Wavefunction convergence fullfilled");
				if(delta < parameters.econv){
					output("Energy connverged");
					pair_energies.push_back(current_energy);
					break;
				}
			}
			}
		}
	}
	output_section("All Pair Energies Converged");
	output("Converged Pair Energies are:");
	if(world.rank()==0){
		for(auto x:pair_energies) std::cout << std::setprecision(parameters.output_prec) << std::fixed << x << std::endl;
	}
	double correlation_energy=0.0;
	for(auto x:pair_energies) correlation_energy += x;
	output("Correlation Energy is: " + stringify(correlation_energy));
	return correlation_energy;
}

// Unnecessary function
double CC2::compute_mp2_pair_energy(CC_Pair &u)const{
	return CCOPS.compute_mp2_pair_energy(u);
}

/// Initialize an electron pair
/// Calculate the constant Term, and the <ij|gQf|ij> etc terms
void CC2::initialize_electron_pair(CC_Pair &u)const{
	output_subsection("Initialize Electron Pair |u" + stringify(u.i) + stringify(u.j)+">");
	output("\n\nCheck for saved pairs...");
	if(u.load_pair(world)){
		output("...Found saved pair\n\n");
		u.info();
		return;
	}else{
	output("...No saved pair found... recalculate\n\n");
	CC_Timer timer_integrals(world,"\n\nMake constant energy Integrals\n\n");
	double tmp = CCOPS.make_ij_gQf_ij(u.i,u.j,u);
	timer_integrals.info();
	output("\n<"+stringify(u.i)+stringify(u.j)+"|gQf|(2.0| "+stringify(u.i)+stringify(u.j) + "> - |" +stringify(u.i)+stringify(u.j)+">) = " + stringify(tmp) + "\n");

	double epsij = CCOPS.get_epsilon(u.i,u.j);
	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * epsij), parameters.lo, parameters.thresh_bsh_6D);

	output_subsection("Calculation of constant MP2 potential");
	CC_Timer timer_const(world,"\n\nCalculation of constant MP2 part\n\n");
	real_function_6d mp2_constant_part = CCOPS.get_MP2_potential_constant_part(u).truncate();
	timer_const.info();

	output_subsection("Apply the Green's Operator");
	CC_Timer timer_Gconst(world,"Apply BSH to constant MP2 part");
	real_function_6d GVPhi = G(-2.0*mp2_constant_part).truncate();
	GVPhi.print_size("G(-2.0*Q(ConstantTerm))");
	u.constant_term = CCOPS.Q12(GVPhi);
	u.constant_term.print_size("Q(G(-2.0*Q(ConstantTerm)))");
	timer_Gconst.info();

	// Make the first guess for the mp2 pair function
	u.function = copy(u.constant_term);

	// Calculate the pair energy
	double test_energy = compute_mp2_pair_energy(u);
	output("Initialized Electron Pair: |u" + stringify(u.i) + stringify(u.j) + "> with pair energy: " + stringify(test_energy) + "\n");
	u.info();
	u.store_pair(world);
	}
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
