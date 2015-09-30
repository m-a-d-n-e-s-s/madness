/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#include "CC2.h"

namespace madness {

bool CC2::test()const{
	bool sanity = true;
	output_section("CC2 CONSISTENCY CHECK");

	output_subsection("Testing Integal Routines:");

	bool integrals_sane = true;
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=i;j<active_mo.size();j++){
			// make empty electron pair
			CC_Pair u(i,j);

			// test the integral routine
			double ij_gQf_ij = CCOPS.make_ijgQfxy(u.i,u.j,CCOPS.mo_ket(i),CCOPS.mo_ket(j));
			double ij_gQf_ji = CCOPS.make_ijgQfxy(u.i,u.j,CCOPS.mo_ket(j),CCOPS.mo_ket(i));
			double ji_gQf_ij = CCOPS.make_ijgQfxy(u.j,u.i,CCOPS.mo_ket(i),CCOPS.mo_ket(j));
			double ji_gQf_ji = CCOPS.make_ijgQfxy(u.j,u.i,CCOPS.mo_ket(j),CCOPS.mo_ket(i));
			output("------------------------------------");
			if(world.rank()==0){
				std::cout << "<" << i << j << "|gQf|" << i << j <<  "> = " << ij_gQf_ij << std::endl;
				std::cout << "<" << j << i << "|gQf|" << j << i <<  "> = " << ji_gQf_ji << std::endl;
				std::cout << "<" << i << j << "|gQf|" << j << i <<  "> = " << ij_gQf_ji << std::endl;
				std::cout << "<" << j << i << "|gQf|" << i << j <<  "> = " << ji_gQf_ij << std::endl;
			}
			double diff1 = fabs(ij_gQf_ij - ji_gQf_ji);
			double diff2 = fabs(ij_gQf_ji - ji_gQf_ij);
			if(diff1 > FunctionDefaults<3>::get_thresh()){
				output("Error in Integrals ij = " +stringify(i) + stringify(j));
				integrals_sane = false;
			}
			if(diff2 > FunctionDefaults<3>::get_thresh()){
				output("Error in Exchange Integrals ij = " +stringify(i) + stringify(j));
				integrals_sane = false;
			}

		}
	}
	sanity = integrals_sane;
	if(integrals_sane) output("Integrals are sane");


	return sanity;
}

/// solve the CC2 ground state equations, returns the correlation energy
double CC2::solve()const{
	// Check if HF is converged
	//if(parameters.debug) solve_CCS();
	//if(parameters.debug) test();
	// Initialize the Pair functions (uij, i>=j)
	if(parameters.restart) output_section("Initialize Electron Pairs: Loading Stored Pairs");
	else output_section("Initialize Electron Pairs: First guess will be the constant Term of MP2");
	CC_Timer timer_init(world,"Initialization of all pairs");
	Pairs<CC_Pair> pairs;
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=i;j<active_mo.size();j++){
			CC_Pair u(i,j);
			if (parameters.restart == true){
				if(u.load_pair(world)){
					output("...Found saved pair\n\n");
					u.info();
					double mp2_energy = CCOPS.compute_mp2_pair_energy(u);
					output("Current MP2 Energy of the Pair is: " + stringify(mp2_energy));
				}else MADNESS_EXCEPTION(("No Restartdata found for pair " + stringify(i)+stringify(j)).c_str(),1);
			}else initialize_electron_pair(u);
			pairs.insert(i,j,u);
		}
	}
	timer_init.info();

	double correlation_energy_mp2=0.0;
	double correlation_energy_cc2=0.0;
	if(not parameters.restart){
		output_section("Solve the uncoupled MP2 equations");
		CC_Timer timer_mp2(world,"Solve MP2 equations");
		correlation_energy_mp2 = solve_uncoupled_mp2(pairs);
		timer_mp2.info();
		output("Solving of MP2 ended at " + stringify(wall_time()) + "s (wall), " +  stringify(cpu_time()) + "s (cpu)");
		timer_mp2.info();
	}
	{
		output_section("Solve the CC2 equations");
		CC_Timer timer_cc2(world,"Solve CC2 equations");
		// Make empty CC_Singles Vector (will be initialized in the solve routine)
		CC_Singles singles;
		pairs(0,0).function.print_size("\t Guess Doubles");
		correlation_energy_cc2 = solve_cc2(pairs,singles);
		timer_cc2.info();
		output("Solving of CC2 ended at " + stringify(wall_time()) + "s (wall), " +  stringify(cpu_time()) + "s (cpu)");
	}

	output_section("Solve CC2 ended");
	if(world.rank()==0) std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << correlation_energy_mp2 << "\nCC2 Correlation Energy is: " << correlation_energy_cc2 << std::endl;
	output_section("Nothing more Implemented right now");
	return correlation_energy_cc2;
}

// Solve the CCS equations for the ground state (debug potential and check HF convergence)
bool CC2::solve_CCS()const{
	output_section("SOLVE CCS");
	// since the symmetry should not change use the projected aos from moldft as guess
	CC_Singles singles;
	real_function_3d guessi = real_factory_3d(world);
	for(size_t i=0;i<nemo.get_calc()->ao.size();i++) guessi += nemo.get_calc()->ao[i];
	CCOPS.Q(guessi);
	for(size_t i=0;i<active_mo.size();i++){
		CC_Single tmp(i,guessi);
		singles.push_back(tmp);
	}

	for(size_t iter=0;iter<30;iter++){
		output_subsection("Iterate CCS: Iteration "+stringify(iter+1));
		CCOPS.update_intermediates(singles);
		CC_Timer timer_potential(world,"CCS Potential");
		vecfuncT potential = CCOPS.get_CCS_potential(singles);
		timer_potential.info();

		output_subsection("Apply the Green's Operator");
		CC_Timer timer_G(world,"Apply the Green's Operator");
		vecfuncT G_potential = zero_functions<double,3>(world,potential.size());
		scale(world,potential,-2.0);
		for(size_t i=0;i<potential.size();i++){
			double epsi = CCOPS.get_orbital_energies()[i];
			real_convolution_3d G = BSHOperator<3>(world, sqrt(-2.0 * epsi), parameters.lo, parameters.thresh_bsh_3D);
			real_function_3d tmp = (G(potential[i])).truncate();
			CCOPS.Q(tmp);
			G_potential[i] = tmp;
		}
		timer_G.info();

		std::vector<double> errors;
		bool converged = true;
		std::vector<double> omega;
		for(size_t i=0;i<potential.size();i++){
			if(world.rank()==0) std::cout << "|| |tau" + stringify(i)+">|| =" << G_potential[i].norm2() << std::endl;
			real_function_3d residue = singles[i].function() - G_potential[i];
			double error = residue.norm2();
			errors.push_back(error);
			if(world.rank()==0) std::cout << "|| residue" + stringify(i)+">|| =" << error << std::endl;
			CCOPS.Q(G_potential[i]);
			singles[i].update(G_potential[i]);
			singles[i].iterations ++;
			if(fabs(error) > parameters.dconv_3D) converged = false;
			omega.push_back(CCOPS.compute_ccs_correlation_energy(singles[i],singles[i]));
		}
		// print out the norms
		output("Current CCS Correlation energies (Diagonal Part)");
		if(world.rank()==0) std::cout << omega << std::endl;
		output("CCS Norms");
		for(size_t i=0;i<singles.size();i++){
			double norm = singles[i].function().norm2();
			if(world.rank()==0) std::cout << norm << std::endl;
		}
		if(converged) break;

	}


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
				CC_Timer timer_iteration(world,"Iteration "+ stringify(iter));
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
					if(fabs(delta) < parameters.econv){
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

double CC2::solve_cc2(Pairs<CC_Pair> &doubles, CC_Singles &singles)const{
	output_section("Initialize CC2 Singles from the MP2 Doubles");
	CC_Timer init_singles(world,"Initialize CC2 Singles");
	doubles(0,0).function.print_size("doubles 1");
	singles = initialize_cc2_singles(doubles);
	init_singles.info();

	// Calculate pair energies with the initialized singles
	CC_Timer timer_cc2_energies(world,"Update CC2 Pair Energies");
	doubles(0,0).function.print_size("doubles 2");
	update_cc2_pair_energies(doubles,singles);
	doubles(0,0).function.print_size("doubles 3");
	timer_cc2_energies.info();

	output_section("Beginn the CC2 Iterations");
	bool singles_converged = false;
	bool doubles_converged = false;
	for(size_t iter=0;iter<parameters.iter_max_6D;iter++){
		CC_Timer timer_iter_all(world,"Iteration " + stringify(iter));
		CCOPS.update_intermediates(singles);
		output_subsection("Iteration "+stringify(iter));
		CC_Timer timer_iter_singles(world,"Iteration " + stringify(iter) + " Singles");
		singles_converged=iterate_cc2_singles(doubles,singles);
		timer_iter_singles.info();
		CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Doubles");
		doubles_converged = iterate_cc2_doubles(doubles,singles);
		timer_iter_doubles.info();
		if(singles_converged and doubles_converged){
			output("Singles and Doubles Converged (!)");
			timer_iter_all.info();
			break;
		}
		timer_iter_all.info();

	}

	return 0.0;
}

std::vector<double> CC2::update_cc2_pair_energies(const Pairs<CC_Pair> &doubles, const CC_Singles &singles)const{
	std::vector<double> omegas;
	for(size_t i=0;i<active_mo.size();i++){
		for(size_t j=0;j<active_mo.size();j++){
			double tmp = CCOPS.compute_cc2_pair_energy(doubles(i,j), singles[i].function(), singles[j].function());
			omegas.push_back(tmp);
		}
	}
	if(world.rank()==0){
		std::cout << "Updated CC2 pair energies:\n";
		for(auto x:omegas) std::cout << std::scientific << std::setprecision(parameters.output_prec) << x << std::endl;
	}
	return omegas;
}

bool CC2::iterate_cc2_singles(const Pairs<CC_Pair> &doubles, CC_Singles &singles)const{
	output_subsection("Iterate CC2 Singles");
	CC_Timer timer_potential(world,"CC2 Singles Potential");
	doubles(0,0).function.print_size("doubles 5");
	vecfuncT potential = CCOPS.get_CC2_singles_potential(singles,doubles);
	timer_potential.info();

	output_subsection("Apply the Green's Operator");
	CC_Timer timer_G(world,"Apply the Green's Operator");
	vecfuncT G_potential = zero_functions<double,3>(world,potential.size());
	scale(world,potential,-2.0);
	for(size_t i=0;i<potential.size();i++){
		double epsi = CCOPS.get_orbital_energies()[i];
		real_convolution_3d G = BSHOperator<3>(world, sqrt(-2.0 * epsi), parameters.lo, parameters.thresh_bsh_3D);
		real_function_3d tmp = (G(potential[i])).truncate();
		CCOPS.Q(tmp);
		G_potential[i] = tmp;
	}
	timer_G.info();

	std::vector<double> errors;
	bool converged = true;
	for(size_t i=0;i<potential.size();i++){
		if(world.rank()==0) std::cout << "|| |tau" + stringify(i)+">|| =" << G_potential[i].norm2() << std::endl;
		real_function_3d residue = singles[i].function() - G_potential[i];
		double error = residue.norm2();
		errors.push_back(error);
		if(world.rank()==0) std::cout << "|| residue" + stringify(i)+">|| =" << error << std::endl;
		CCOPS.Q(G_potential[i]);
		singles[i].update(G_potential[i]);
		singles[i].iterations ++;
		if(fabs(error) > parameters.dconv_3D) converged = false;
	}
	return converged;
}

bool CC2::iterate_cc2_doubles(Pairs<CC_Pair> &doubles, const CC_Singles &singles)const{
	return false;
}

CC_Singles CC2::initialize_cc2_singles(const Pairs<CC_Pair> &doubles)const{

	output_subsection("Calculate the singles guess potential: S2b+X + S2c+X");
	CC_Timer timer_guess_potential(world,"Calculate the singles guess potential");
	vecfuncT guess_potential = CCOPS.get_CC2_singles_initial_potential(doubles);
	timer_guess_potential.info();

	output_subsection("Apply the Green's Operator");
	CC_Timer timer_G(world,"Apply the Green's Operator");
	vecfuncT G_guess_potential = zero_functions<double,3>(world,guess_potential.size());
	scale(world,guess_potential,-2.0);
	for(size_t i=0;i<guess_potential.size();i++){
		double epsi = CCOPS.get_orbital_energies()[i];
		real_convolution_3d G = BSHOperator<3>(world, sqrt(-2.0 * epsi), parameters.lo, parameters.thresh_bsh_3D);
		real_function_3d tmp = (G(guess_potential[i])).truncate();
		CCOPS.Q(tmp);
		G_guess_potential[i] = tmp;
	}
	timer_G.info();

	output_section("Initialized CC2 Singles");
	CC_Singles singles;
	for(size_t i=0;i<guess_potential.size();i++){
		if(world.rank()==0) std::cout << "|| |tau" + stringify(i)+">|| =" << G_guess_potential[i].norm2() << std::endl;
		CC_Single taui(i,G_guess_potential[i]);
		singles.push_back(taui);
	}

	return singles;
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

		output("...No saved pair found... recalculate\n\n");
		CC_Timer timer_integrals(world,"Make constant energy Integrals");
		u.ij_gQf_ij = CCOPS.make_ijgQfxy(u.i,u.j,active_mo[u.i],active_mo[u.j]);
		u.ji_gQf_ij = CCOPS.make_ijgQfxy(u.i,u.j,active_mo[u.j],active_mo[u.i]);
		timer_integrals.info();

		double epsij = CCOPS.get_epsilon(u.i,u.j);
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * epsij), parameters.lo, parameters.thresh_bsh_6D);

		output_subsection("Calculation of constant MP2 potential");
		CC_Timer timer_const(world,"Calculation of constant MP2 part");
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

/// Calculate the current CC2 correlation energy
double CC2::compute_correlation_energy(const vecfuncT &singles, const Pairs<real_function_6d> &doubles)const{
	MADNESS_EXCEPTION("CC2::compute_correlation_energy Not Yet Implemented",1);
return 0.0;
}
/// Iterates the CC2 singles equations
void CC2::iterate_singles(vecfuncT &singles, const Pairs<real_function_6d> &doubles)const{
	MADNESS_EXCEPTION("CC2::iterate_singles Not Yet Implemented",1);
}
/// Iterates the CC2 doubles equations
void CC2::iterate_doubles(const vecfuncT &singles, Pairs<real_function_6d> &doubles)const{
	MADNESS_EXCEPTION("CC2::iterate_doubles Not Yet Implemented",1);
}

} /* namespace madness */
