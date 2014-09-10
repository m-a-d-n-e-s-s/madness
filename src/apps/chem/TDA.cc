/*
 * TDA.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: kottmanj
 */

#include <chem/TDA.h>

using namespace madness;

/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};

void TDA::solve(xfunctionsT &xfunctions) {

	if (plot_ == true)
		plot_vecfunction(active_mo_, "active_mo_");

	// check for saved xfunctions
	if (read_) {
		std::cout << std::setw(40) << "Found read keyword..." << " : checking for saved xfunctions" << std::endl;
		xfunctionsT read_xfunctions;
		for (size_t i = 0; i < excitations_; i++) {
			xfunction tmp(world);
			vecfuncT x(active_mo_.size());
			tmp.x = x;
			tmp.delta.push_back(10.0);
			tmp.error.push_back(10.0);
			std::cout << std::setw(40) << "Reading xfunction " << i;
			if (load_vecfunction(world, tmp.x, stringify(i))){
			std::cout << " : found" << std::endl;
			tmp.omega = guess_omega_;
			converged_xfunctions_.push_back(tmp);
			}
			else std::cout << " : not found" << std::endl;
		}
		// Iterate one time to check level of convergence
		print("\n----ITERATE READ XFUNCTIONS------\n\n");
		normalize(converged_xfunctions_);
		orthonormalize_fock(converged_xfunctions_,false);
		normalize(converged_xfunctions_);
		for(size_t i=0;i<converged_xfunctions_.size();i++){
			converged_xfunctions_[i].omega = converged_xfunctions_[i].expectation_value.back();
			iterate_one(converged_xfunctions_[i],false,false);
			if(converged_xfunctions_[i].error.back() > dconv_){
				print("Warning: Read xfunctions may not converged");
			}
			print_xfunction(converged_xfunctions_[i]);
		}
		print("\n----END ITERATE READ XFUNCTIONS-----\n\n");

	}


		// Create the excitation_function vector and initialize
		std::cout << "---Start Initialize Guess Functions---" << "\n\n " << std::endl;
		TDA_TIMER init(world,"\nfinished to initialize guess excitations ...");
		if(xfunctions.empty())initialize(xfunctions);
		iterate_guess(xfunctions);
		std::cout << std::setw(100) << "---End Initialize Guess Functions---" << " " << std::endl;
		init.info();

		// make shure everything is fine
		if(kain_) {
			kain_solvers.increase_subspace(world,xfunctions);
			kain_solvers.sanity_check(xfunctions);
		}

		// Iterate till convergence is reached
		iterate(xfunctions);

		bool converged = check_convergence(xfunctions);
		if(converged==true) {
			if(world.rank()==0)std::cout << std::setw(100) << stringify(excitations_)+" excitations CONVERGED!" << std::endl;
		}
		if(converged==false) {
			if(world.rank()==0)std::cout <<std::setw(100) << stringify(excitations_)+"excitations NOT CONVERGED!" << std::endl;
		}

		// Save to disk and plot
		for(size_t i=0;i<converged_xfunctions_.size();i++) {
			save_vecfunction(world,converged_xfunctions_[i].x,stringify(i));
			plot_vecfunction(converged_xfunctions_[i].x,"final_guess_excitation_"+stringify(i),true);
		}
		// print final result
		print_status(xfunctions);
		print_performance(xfunctions,"pre-");

		// Analyze
		//analyze(xfunctions);


}

void TDA::solve_sequential(xfunctionsT xfunctions) {

	print("\n\n\n\n-------------------------------------------------------");
	print("BEGINNING THE FINAL ITERATIONS TO AN ACCURACY OF ", hard_dconv_);
	print("-------------------------------------------------------\n\n\n\n");

	// The given xfunctions should be sorted by energy in ascending order:
	for (size_t iroot = 0; iroot < xfunctions.size(); iroot++) {
		print("\n\n-----xfunction ", iroot, " ------\n\n");
		// create the kain solver for the current xfunction
		solverT solver(allocator(world, xfunctions[iroot].x.size()));
		solver.set_maxsub(3);
		kain_ = true;

		// begin the final iterations
		for (int iter = 0; iter < 100; iter++) {
			TDA_TIMER iteration_timer(world, "");
			normalize(xfunctions[iroot]);
			// on the fly or not makes no sense here, but since the same input file is used for both solve functions this has to be here
			if (not on_the_fly_)
				on_the_fly_ = true;
			// first false: expectation value is calculated and saved, second false: no guess calculation
			iterate_one(xfunctions[iroot], false, false);
			normalize(xfunctions[iroot]);

			// update with kain
			xfunction tmp = solver.update(xfunction(world, xfunctions[iroot].x),
					xfunction(world, xfunctions[iroot].current_residuals));
			xfunctions[iroot].x = tmp.x;
			// need to pack another vector to proejct out the converged functions (function only takes vectors)
			xfunctionsT courier;
			courier.push_back(xfunctions[iroot]);
			project_out_converged_xfunctions(courier);
			if (fabs(xfunctions[iroot].delta.back())
					< xfunctions[iroot].error.back() * 1.e-2
					or fabs(xfunctions[iroot].delta.back()) < 8.e-4) {
				xfunctions[iroot].omega += xfunctions[iroot].delta.back();
			} else
				xfunctions[iroot].omega =
						xfunctions[iroot].expectation_value.back();

			// print update

			//print("Iteration ", iter ," on xfunction " ,iroot, " starts at ",wall_time());
			print_xfunction(xfunctions[iroot]);
			std::cout << "time: "; iteration_timer.info(); std::cout << std::endl;

			// check convergence
			if (xfunctions[iroot].error.back() < hard_dconv_) {
				std::cout << "\n ------xfunction " << iroot << " converged!!! -----\n" << std::endl;
				xfunctions[iroot].converged =true;
				converged_xfunctions_.push_back(xfunctions[iroot]);
				break;
			}
			if (iter > iter_max_) {
				std::cout << "\n ------xfunction " << iroot << " did not converge ------\n" << std::endl;
				xfunctions[iroot].converged =false;
				converged_xfunctions_.push_back(xfunctions[iroot]);
				break;
			}

		}

	}
	// Save to disk and plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		save_vecfunction(world,converged_xfunctions_[i].x,stringify(i));
		plot_vecfunction(converged_xfunctions_[i].x,"final_excitation_"+stringify(i),true);
	}

	print("Final result :");
	print_status(xfunctions);
	print_performance(converged_xfunctions_,"final-");
	analyze(converged_xfunctions_);

}

void TDA::print_status(const xfunctionsT & xfunctions) const {
	std::cout << "\n" <<std::setw(5) << " #" << std::setw(20) << "omega" << std::setw(20) << "delta" << std::setw(20)
	<< "error"<<std::setw(20)
	<<"expv" << std::setw(7) <<"iter"<< std::setw(7)<< "conv" << std::endl;
	for(size_t i=0;i<converged_xfunctions_.size();i++) print_xfunction(converged_xfunctions_[i]);
	std::cout <<"--------"<< std::endl;
	for(size_t i=0;i<xfunctions.size();i++)print_xfunction(xfunctions[i]);
	std::cout << "\n" << std::endl;
}

void TDA::print_xfunction(const xfunction &x) const {
	std::cout << std::setw(5) << x.number;
	std::cout << std::scientific << std::setprecision(10) << std::setw(20) << x.omega << std::setw(20)<< x.delta.back()
	<< std::setw(20)<< x.error.back()<< std::setw(20) << x.expectation_value.back();
	std::cout << std::fixed <<std::setw(7)<< x.iterations << "   " << std::setw(7)<<x.converged << std::endl;
}

void TDA::initialize(xfunctionsT & xfunctions) {
	if (guess_ == "physical") {
		guess_physical(xfunctions);
	} else if (guess_ == "...")
		print("dummy");
	else {
		if (world.rank() == 0)
			print("unknown keyword for guess: ", guess_);
		MADNESS_EXCEPTION("Reached end of initialize", 1);
	}

}

void TDA::guess_physical(xfunctionsT & xfunctions) {

	// when dft is used add diffuse functions to the mos to increase delocalisation
	vecfuncT diffuse_mos;
	if (dft_) {
		for (size_t i = 0; i < active_mo_.size(); i++) {
			real_function_3d tmp = copy(active_mo_[i]);
			diffuse_mos.push_back(tmp);
		}
		add_diffuse_functions(diffuse_mos);
	}

	vecfuncT xoperators = make_excitation_operators();

	real_function_3d all_orbitals = real_factory_3d(world);
	if (guess_mode_ == "all_orbitals") {
		for (size_t i = 0; i < mos_.size(); i++) {
			all_orbitals += mos_[i];
		}
		all_orbitals.truncate(truncate_thresh_);
	}

	for (size_t j = 0; j < xoperators.size(); j++) {

		vecfuncT x;
		for (size_t i = 0; i < active_mo_.size(); i++) {
			real_function_3d tmp;
			if (guess_mode_ == "all_orbitals") {
				tmp = all_orbitals * xoperators[j];
			} else
				tmp = active_mo_[i] * xoperators[j];

			//if(dft_)tmp += diffuse_mos[i]*xoperators[j];
			double norm = tmp.norm2();
			tmp.scale(1.0 / norm);
			x.push_back(copy(tmp));
		}
		xfunction xtmp(world, x);
		xtmp.number = xfunctions.size();
		xtmp.omega = guess_omega_;
		xtmp.kain = false;
		normalize(xtmp);
		plot_vecfunction(xtmp.x,
				"guess_excitation_" + stringify(xtmp.number) + "_");
		xfunctions.push_back(xtmp);
	}
	normalize(xfunctions);

	// Initalize energies
	double factor = 1.0;
	for (size_t i = 0; i < xfunctions.size(); i++) {
		xfunctions[i].omega = guess_omega_ * factor;
		//factor += 0.05;
	}

	for (size_t i = 0; i < xfunctions.size(); i++)
		plot_vecfunction(xfunctions[i].x,
				"guess_excitation_after_fock" + stringify(i) + "_");

}

void TDA::add_diffuse_functions(vecfuncT &mos) {

	// mos must be reconstructed for the eval function
	reconstruct(world, mos);

	// get the box size L and the molecule bounding cube bc
	double L = calc_.param.L;

	// Get all coordinates and pack them into the right format
	Tensor<double> Coord = calc_.molecule.get_all_coords();
	std::vector<coord_3d> coord;
	for (size_t i = 0; i < calc_.molecule.natom(); i++) {
		coord_3d tmp;
		for (int j = 0; j < 3; j++)
			tmp[j] = Coord(i, j);
		coord.push_back(tmp);
	}

	// Add diffuse 1s functions with the right sign to each mo
	for (size_t k = 0; k < mos.size(); k++) {

		// Check the signs of the orbital on the nuclei
		std::vector<int> signs;
		double tol = 10.0 * FunctionDefaults<3>::get_thresh();
		for (size_t p = 0; p < coord.size(); p++) {
			double value = mos[k].eval(coord[p]);
			if (fabs(value) < tol)
				signs.push_back(0.0);
			else if (value < 0.0)
				signs.push_back(-1.0);
			else if (value > 0.0)
				signs.push_back(1.0);
		}
		std::cout << "signs\n" << signs << std::endl;

		// make the diffuse functions with the right signs
		real_function_3d diffuse_tmp;
		if (dft_) {
			std::shared_ptr<FunctionFunctorInterface<double, 3> > rydberg_functor(
					new diffuse_functions(L, coord, signs, mos.size(), 1));
			diffuse_tmp = real_factory_3d(world).functor(rydberg_functor);
		} else {
			std::shared_ptr<FunctionFunctorInterface<double, 3> > rydberg_functor(
					new diffuse_functions(L, coord, signs, mos.size(), 0));
			diffuse_tmp = real_factory_3d(world).functor(rydberg_functor);
		}

		// add them to the corresponding mo
		mos[k] += diffuse_tmp;

	}

}

vecfuncT TDA::make_excitation_operators() const {

	std::string key;
	if (guess_exop_ == "symmetry") {
		key = calc_.molecule.pointgroup_;
		std::cout << std::setw(40) << "symmetry_guess" << " : point group is " << key;
		if (key == "C1") {
			key = "quadrupole+";
			std::cout<< "and not yet implemented ...";}
		if (key == "Cs") {
			key = "quadrupole+";
			std::cout<< "and not yet implemented ...";}
		if (key == "Ci") {
			key = "quadrupole+";
			std::cout<< "and not yet implemented ...";}
		if (key == "D2") {
			key = "quadrupole+";
			std::cout<< "and not yet implemented ...";}
		std::cout << std::endl;
	}
	else key = guess_exop_;

	exoperators exops;
	vecfuncT xoperators = exops.get_exops(world, key);

	return xoperators;
}

void TDA::iterate_guess(xfunctionsT &xfunctions) {
	std::cout << "\n" << std::endl;
	std::cout << "---Start Guess Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,true);
	// set back iteration counter
	for(size_t i=0;i<xfunctions.size();i++)xfunctions[i].iterations = 0;
	std::cout <<std::setw(100)<< "---End Guess Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate(xfunctionsT &xfunctions) {
	std::cout << "\n" << std::endl;
	std::cout << "---Start Main Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,false);
	std::cout <<std::setw(100)<< "---End Main Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate_all(xfunctionsT &xfunctions, bool guess) {
	size_t max;
	if (guess == true) {max = guess_iter_;}
	if (guess == false) {max = 1000;}	//iter_max_; // using iter_max now for single roots to kick out}

	// start iteration cycle
	for (size_t current_iteration = 0; current_iteration < max;
			current_iteration++) {

		// get information about the used memory
		double xfunction_size = memwatch(xfunctions,debug_);
		double conv_xfunction_size = memwatch(converged_xfunctions_,debug_);
		if(debug_)print("mos ",get_size(world,mos_)," (GB)\n");
		std::cout << std::scientific << std::setprecision(0);
		std::cout << std::setw(40) << " memory information..." << " : ";
		std::cout << "xfunctions " << xfunction_size << " (GB), ";
		std::cout <<" converged xfunctions "<< conv_xfunction_size <<" (GB)" << std::endl;


		int iter = current_iteration;
		if (guess == false)
			iter += guess_iter_;
		std::cout<< std::setw(40)<<"Starting iteration"+stringify(iter)+" at time   "
		<<" : "<<wall_time() << std::endl;
		// Check memory management
		for (size_t i = 0; i < xfunctions.size(); i++) {
			if (on_the_fly_ and not xfunctions[i].Vx.empty())
				print(
						"MEMORY PROBLEM: On the fly calculation of the applied potential is ",
						on_the_fly_, " but Vx of xfunction ", i,
						" is not empty, size is ", xfunctions[i].Vx.size());
			if (not kain_ and not xfunctions[i].current_residuals.empty())
				print("MEMORY PROBLEM: Kain is ", kain_,
						" but Residual-vector of xfunction-struct ", i,
						" is not empty, size is ",
						xfunctions[i].current_residuals.size());
		}
		bool ptfock_orthonormalization = true;
		if (only_fock_)
			ptfock_orthonormalization = true;
		if (only_GS_)
			ptfock_orthonormalization = false;

		TDA_TIMER Iterationtime(world,
				"\nIteration " + stringify(current_iteration) + " time: ");

		TDA_TIMER poo(world, "project out occ. space...");
		for (size_t p = 0; p < xfunctions.size(); p++) {
			project_out_occupied_space(xfunctions[p].x);
		}
		poo.info();

		normalize(xfunctions);

		// Truncate xfunctions
		TDA_TIMER truncate1(world, "truncate |x> and |Vx>...");
		truncate_xfunctions(xfunctions);
		truncate1.info();

		if (not on_the_fly_) {
			// Update the perturbed vectorpotentials
			TDA_TIMER update_potentials(world, "update all potentials...");
			for (size_t i = 0; i < xfunctions.size(); i++) {
				if (xfunctions[i].converged == false)
					xfunctions[i].Vx = apply_perturbed_potential(xfunctions[i]);
			}
			update_potentials.info();

			// Truncate xfunctions
			TDA_TIMER truncate2(world, "truncate |x> and |Vx>...");
			truncate_xfunctions(xfunctions);
			truncate2.info();
		}

		if (ptfock_orthonormalization) {
			TDA_TIMER orthonormalize_fock_matrix(world, "orthogonalize_fock : ");

			ptfock_orthonormalization = orthonormalize_fock(xfunctions, guess);
			orthonormalize_fock_matrix.info();

			// Update the energies (not during the guess calculations)
			if (guess == false) update_energies(xfunctions);
			else std::cout << std::setw(40) << "update energies..." << " : " << "not during guess iterations" << std::endl;

			// Truncate xfunctions
			TDA_TIMER truncate3(world, "truncate |x> and |Vx>...");
			truncate_xfunctions(xfunctions);
			truncate3.info();
		}
		// if there is no fock update demanded: just update the energy
		else if (guess == false) update_energies(xfunctions);

		TDA_TIMER apply_bsh(world, "apply BSH operators...");
		// apply BSH operator
		std::vector<size_t> frozen_excitations;
		for (size_t i = 0; i < xfunctions.size(); i++) {iterate_one(xfunctions[i], ptfock_orthonormalization, guess);}
		apply_bsh.info();

		if (not only_fock_)
			orthonormalize_GS(xfunctions);
		else {
			normalize(xfunctions);
		}

		// Update with the KAIN solver (not during guess iterations)
		if (kain_ and not guess) {
			TDA_TIMER kain(world, "update with kain solver...");
			print("\n");
			kain_solvers.update(xfunctions);
			kain.info();
		}

		//check convergence (but not for guess iterations)
		if (guess == false) {
			bool converged = check_convergence(xfunctions);
			project_out_converged_xfunctions(xfunctions);
			//if(kain_) kain_solvers.reduce_subspace(xfunctions);
			std::cout<<"converged : "<<converged<<std::endl;
			if (converged == true)
				break;
			if (converged_xfunctions_.size() >= excitations_)
				break;

		} else
			project_out_converged_xfunctions(xfunctions);

		Iterationtime.info();
		print_status(xfunctions);
		for (size_t p = 0; p < xfunctions.size(); p++) {
			plot_vecfunction(xfunctions[p].x,
					"current_iteration_" + stringify(p) + "_", plot_);
		}
	}
}

void TDA::iterate_one(xfunction & xfunction, bool ptfock, bool guess) {
	xfunction.iterations += 1;

	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("DFT is off but potential shift is not zero", 1);

	vecfuncT Vpsi;
	if (on_the_fly_)
		Vpsi = apply_perturbed_potential(xfunction);
	else if (not on_the_fly_)
		Vpsi = xfunction.Vx;

	double omega = xfunction.omega;

	// If orthonormalize_fock was not used previously the expectation value must be calculated:
	if (not ptfock){
		xfunction.expectation_value.push_back(
				expectation_value(xfunction, Vpsi));
	}
	scale(world, Vpsi, -2.0);
	truncate(world, Vpsi, truncate_thresh_); // no fence

	std::vector<poperatorT> bsh(active_mo_.size());
	for (size_t p = 0; p < active_mo_.size(); p++) {
		double eps = get_calc().aeps[p + nfreeze_] + omega + shift_; // BSH parametrization
		if (eps > 0) {
			if (world.rank() == 0)
				print("bsh: warning: positive eigenvalue", p + nfreeze_, eps);
			eps = get_calc().aeps[p + nfreeze_] + guess_omega_;
		}
		bsh[p] = poperatorT(
				BSHOperatorPtr3D(world, sqrt(-2.0 * eps), lo, bsh_eps_));
	}

	world.gop.fence();

	vecfuncT GVpsi = apply(world, bsh, Vpsi);

	// Residual as error estimation
	project_out_occupied_space(GVpsi);
	vecfuncT residual = sub(world, xfunction.x, GVpsi);
	Tensor<double> inner_tmp = inner(world, residual, residual);
	double error = sqrt(inner_tmp.sum());
	xfunction.error.push_back(error);
	if (error < 8.e-2)
		xfunction.kain = true;

	// Calculate 2nd order update:
	// Inner product of Vpsi and the residual (Vspi is scaled to -2.0 --> multiply later with 0.5)
	Tensor<double> tmp = inner(world, Vpsi, residual);

	// Norm of GVpsi (Psi_tilde)
	double tmp2 = 0.0;
	for (size_t i = 0; i < GVpsi.size(); ++i) {
		double n = GVpsi[i].norm2();
		tmp2 += n * n;
	}

	// Factor 0.5 removes the factor 2 from the scaling before
	xfunction.delta.push_back(0.5 * tmp.sum() / tmp2);

	// Update x-function (only if kain_ is false or during the guess iterations where kain is not used)
	if (not kain_ or guess) {
		xfunction.x = GVpsi;
		project_out_occupied_space(xfunction.x);
	}

	// Use tmp2 to normalize the new x-functions GVPsi:
	scale(world, xfunction.x, 1.0 / tmp2);

	// when kain update should be used we will need the residual vectors later
	if (kain_ and not guess) {
		//print("do not rescale residual wiht -1.0");
		//scale(world,residual,-1.0); //???????
		truncate(world, residual,truncate_thresh_);
		xfunction.current_residuals = residual;
	}
	// when kain is not used the current_residual vector of the xfunction struct should be empty all the time (check this here to avoid memory overflow)
	if (not kain_) {
		if (xfunction.current_residuals.size() != 0)
			MADNESS_EXCEPTION(
					"KAIN is turned off, but the current_residuals vector of the xfunction structures is not empty",
					1);
	}

}

void TDA::update_energies(xfunctionsT &xfunctions) {
	std::cout << std::setw(40) << "update energies..." << " : ";
	for(size_t k=0;k<xfunctions.size();k++) {

		//failsafe: make shure the delta and expectation values vectors are not empty to avoid segmentation faults
		if(not xfunctions[k].delta.empty() and not xfunctions[k].expectation_value.empty() and not xfunctions[k].error.empty()) {
			if(xfunctions[k].expectation_value.back() < highest_excitation_) {
				if(fabs(xfunctions[k].delta.back()) < xfunctions[k].error.back()*1.e-2 or fabs(xfunctions[k].delta.back())< 8.e-4) {
					xfunctions[k].omega +=xfunctions[k].delta.back();
					std::cout << k << "(2nd), ";
				} else {
					xfunctions[k].omega = xfunctions[k].expectation_value.back();
					std::cout << k << "(exp), ";
				}
			} else {
				// set the last converged value of the same type of guess as the default
				double new_omega = highest_excitation_*0.9;// default
				xfunctions[k].omega = new_omega;
				std::cout << k << "(setback), ";
			}

		}

	}std::cout<<std::endl;

}

void TDA::normalize(xfunctionsT & xfunctions) {
	for (size_t i = 0; i < xfunctions.size(); i++)
		normalize(xfunctions[i]);
}
void TDA::normalize(xfunction & xfunction) {
	const Tensor<double> self_overlaps = inner(world, xfunction.x, xfunction.x);
	const double squared_norm = self_overlaps.sum();
	scale(world, xfunction.x, 1.0 / sqrt(squared_norm));
}

void TDA::project_out_converged_xfunctions(xfunctionsT & xfunctions) {
	TDA_TIMER project(world, "project out converged xfunctions...");

	for (size_t p = 0; p < xfunctions.size(); p++) {
		for (size_t k = 0; k < converged_xfunctions_.size(); k++) {
			Tensor<double> overlap = inner(world, xfunctions[p].x,
					converged_xfunctions_[k].x);
			double c = overlap.sum();
			for (size_t i = 0; i < xfunctions[p].x.size(); i++) {
				xfunctions[p].x[i] -= c * converged_xfunctions_[k].x[i];
			}
		}
	}
	project.info(low_calculation_);
}

void TDA::orthonormalize_GS(xfunctionsT &xfunctions) {
	TDA_TIMER gs(world, "Gram-Schmidt-Orthonormalization...");
	// First normalize
	normalize(xfunctions);

	// Project converged roots out of unconverged
	project_out_converged_xfunctions(xfunctions);

	//orthonormalize the rest with gram schmidt procedure
	for (size_t r = 0; r < xfunctions.size(); r++) {
		if (xfunctions[r].converged == false) {
			for (size_t p = 0; p < r; p++) {
				if (xfunctions[p].converged == false) {
					Tensor<double> overlap = inner(world, xfunctions[r].x,
							xfunctions[p].x);
					const double c = overlap.sum();
					for (size_t i = 0; i < xfunctions[r].x.size(); i++)
						xfunctions[r].x[i] -= c * xfunctions[p].x[i];
				}
			}
		}
	}
	// and again to be shure
	normalize(xfunctions);
	gs.info();

}
//
bool TDA::orthonormalize_fock(xfunctionsT &xfunctions, bool guess) {

	Tensor<double> overlap(xfunctions.size(), xfunctions.size());
	for (size_t p = 0; p < xfunctions.size(); p++) {
		for (size_t k = 0; k < xfunctions.size(); k++) {
			Tensor<double> overlap_vec = inner(world, xfunctions[p].x,
					xfunctions[k].x);
			overlap(p, k) = overlap_vec.sum();
		}
	}

	std::cout<< "\n overlap matrix\n" << overlap << std::endl;
	// if the overlap matrix is already the unit matrix then no orthogonalization is needed
	double overlap_offdiag = measure_offdiagonality(overlap, xfunctions.size());
	if (fabs(overlap_offdiag) < FunctionDefaults<3>::get_thresh()) {
		std::cout << " already orthogonal: perturbed fock matrix will not be calculated \n" <<std ::endl;
		return false;
	}

	Tensor<double> F(3L, xfunctions.size());
	if (dft_) {
		F = make_perturbed_fock_matrix(xfunctions);
	} else
		F = make_perturbed_fock_matrix(xfunctions);

	// Diagonalize the perturbed Fock matrix
	Tensor<double> U, evals;

	madness::Tensor<double> dummy(xfunctions.size());
	U = calc_.get_fock_transformation(world, overlap, F, evals, dummy,
			1.5 * econv_);
	//}
	print("\n\n");
	std::cout<<std::setw(40)<< "Transformation-Matrix-U \n" << U << std::endl;
	print("\n\n");

	//Prevent printout when expectation value is calculated
	if (xfunctions.size() > 1) {
		// Make an estimation how "strong" the xfunctions are linear combinated
		double offdiagonal = measure_offdiagonality(U, xfunctions.size());
		std::cout << std::setw(40) << "offdiagonal transformation part..." << " : " << fabs(offdiagonal) << std::endl;
	}

	// Transform the xfunctions, and if demanded the kain subspace and the potentials
	// X-functions
	std::vector<vecfuncT> old_x;
	for (size_t i = 0; i < xfunctions.size(); i++) {
		old_x.push_back(xfunctions[i].x);
	}
	std::vector<vecfuncT> new_x = transform_vecfunctions(old_x, U);
	for (size_t i = 0; i < xfunctions.size(); i++) {
		xfunctions[i].x = new_x[i];
		xfunctions[i].expectation_value.push_back(evals(i));
	}
	std::cout << std::setw(40) << "Transforming..." << " : xfunctions... ";

	// potentials
	if (not on_the_fly_) {
		std::vector<vecfuncT> old_Vx;
		for (size_t i = 0; i < xfunctions.size(); i++) {
			old_Vx.push_back(xfunctions[i].Vx);
		}
		std::vector<vecfuncT> new_Vx = transform_vecfunctions(old_Vx, U);
		for (size_t i = 0; i < xfunctions.size(); i++) {
			xfunctions[i].Vx = new_Vx[i];
		}
		std::cout << "potentials... ";
	}

	// Kain subspace
	if (kain_ and not guess) {
		kain_solvers.transform_subspace(world, U);
		std::cout << "Kain subspace...";
	}
	std::cout << std::endl;

	return true;

}

double TDA::measure_offdiagonality(const madness::Tensor<double> &U,
		const size_t size) const {
	if(size ==0) return 0.0;
	std::vector<double> offdiag_elements;
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			if(i!=j)offdiag_elements.push_back(fabs(U(i, j)));
		}
	}
	std::sort(offdiag_elements.begin(),offdiag_elements.end());

	return offdiag_elements.back();
}

std::vector<vecfuncT> TDA::transform_vecfunctions(
		const std::vector<vecfuncT> &xfunctions,
		const madness::Tensor<double> U) const {

	std::vector<vecfuncT> new_xfunctions(xfunctions.size());
	for (std::size_t i = 0; i < xfunctions.size(); i++) {
		new_xfunctions[i] = zero_functions<double, 3>(world,
				xfunctions[i].size());
		compress(world, new_xfunctions[i]);
		compress(world, xfunctions[i]);
	}

	for (size_t i = 0; i < xfunctions.size(); i++) {
		for (size_t j = 0; j < xfunctions.size(); ++j) {
			gaxpy(world, 1.0, new_xfunctions[i], U(j, i), xfunctions[j]);
		}
	}

	// Return the transformed vector of vecfunctions
	return new_xfunctions;

}

void TDA::project_out_occupied_space(vecfuncT &x) {
	if (x.size() != active_mo_.size())
		MADNESS_EXCEPTION(
				"ERROR IN PROJECTOR: Size of xfunctions and active orbitals is not equal",
				1);
	for (size_t p = 0; p < active_mo_.size(); ++p)
		x[p] -= rho0(x[p]);
}

double TDA::perturbed_fock_matrix_element(const vecfuncT &xr,
		const vecfuncT &Vxp, const vecfuncT &xp) const {
	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("No DFT calculation but shift is not zero", 1);
	if (xr.size() != Vxp.size())
		MADNESS_EXCEPTION(
				"ERROR IN EXPECTATION VALUE: Function vectors with different sizes",
				1);
	Tensor<double> tmp = inner(world, xr, Vxp);
	double value = tmp.sum();

	// Epsilon Part
	double weighted_sum = 0.0;
	Tensor<double> overlaps = inner(world, xr, xp);
	for (size_t i = 0; i < xr.size(); i++) {
		weighted_sum += (get_calc().aeps[nfreeze_ + i] + shift_) * overlaps[i];
	}

	return value - weighted_sum;
}

double TDA::expectation_value(const xfunction &x, const vecfuncT &Vx) {
	Tensor<double> pot = inner(world, x.x, Vx);
	double expv = pot.sum();
	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	for (int axis = 0; axis < 3; axis++) {
		const vecfuncT dx = apply(world, *(gradop[axis]), x.x);
		Tensor<double> kin = inner(world, dx, dx);
		expv += 0.5 * kin.sum();
	}
	for (size_t i = 0; i < x.x.size(); i++) {
		Tensor<double> overlap = inner(world, x.x, x.x);
		expv -= (get_calc().aeps(nfreeze_ + i) + shift_) * overlap(i);
	}
	return expv;
}

Tensor<double> TDA::make_perturbed_fock_matrix(
		const xfunctionsT &xfunctions) const {
	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("No DFT calculation but shift is not zero", 1);
	std::cout << std::setw(40) << "perturbed fock matrix dimension..." << " : " << xfunctions.size() << "x" << xfunctions.size() << std::endl;

	//Tensor<double> F(xfunctions.size(),xfunctions.size());
	Tensor<double> F(xfunctions.size(), xfunctions.size());

	//The potential part
	// if bool on_the_fly_ is true the potential is calculated, if not the potential has been calculated before and is stored in the
	// xfunction structure via : xfunction.Vx[i]
	// the later is too costly in memory for our clusters right know
	if (on_the_fly_) {
		for (std::size_t p = 0; p < xfunctions.size(); p++) {
			vecfuncT Vxp = apply_perturbed_potential(xfunctions[p]);
			for (std::size_t k = 0; k < xfunctions.size(); k++) {
				Tensor<double> fpk_i = inner(world, xfunctions[k].x, Vxp);
				F(p, k) = fpk_i.sum();
			}
		}

	} else {
		for (std::size_t p = 0; p < xfunctions.size(); p++) {
			for (std::size_t r = 0; r < xfunctions.size(); r++) {
				Tensor<double> fpr_i = inner(world, xfunctions[p].x,
						xfunctions[r].Vx);
				F(p, r) = fpr_i.sum();
			}
		}
	}
	//The kinetic part -1/2<xip|nabla^2|xir> = +1/2 <nabla xip||nabla xir>
	for (std::size_t iroot = 0; iroot < xfunctions.size(); ++iroot) {
		const vecfuncT& xp = xfunctions[iroot].x;
		reconstruct(world, xp);
	}

	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);

	for (int axis = 0; axis < 3; ++axis) {

		std::vector<vecfuncT> dxp;
		for (std::size_t iroot = 0; iroot < xfunctions.size(); ++iroot) {

			const vecfuncT& xp = xfunctions[iroot].x;
			const vecfuncT d = apply(world, *(gradop[axis]), xp);
			dxp.push_back(d);
		}
		for (std::size_t iroot = 0; iroot < xfunctions.size(); ++iroot) {
			for (std::size_t jroot = 0; jroot < xfunctions.size(); ++jroot) {
				Tensor<double> xpi_Txqi = inner(world, dxp[iroot], dxp[jroot]);
				F(iroot, jroot) += 0.5 * xpi_Txqi.sum();
			}
		}
	}

	// The epsilon part
	for (std::size_t p = 0; p < xfunctions.size(); p++) {
		for (std::size_t r = 0; r < xfunctions.size(); r++) {
			Tensor<double> eij = inner(world, xfunctions[p].x, xfunctions[r].x);
			for (size_t ii = 0; ii < xfunctions[p].x.size(); ++ii) {
				F(p, r) -= (get_calc().aeps[ii + nfreeze_] + shift_) * eij[ii];
			}
		}
	}

	return F;

}

vecfuncT TDA::apply_perturbed_potential(const xfunction & xfunction) const {

	vecfuncT Gamma;
	TDA_TIMER gammatimer(world, "make gamma...");
	if (not dft_)
		Gamma = apply_gamma(xfunction);
	if (dft_)
		Gamma = apply_gamma_dft(xfunction);
	gammatimer.info(debug_);

	TDA_TIMER vxctimer(world, "apply the unperturbed potential...");
	vecfuncT V0 = get_V0(xfunction.x);
	vxctimer.info(debug_);

	vecfuncT Vpsi = add(world, V0, Gamma);
	return Vpsi;

}

vecfuncT TDA::apply_gamma(const xfunction &xfunction) const {

	TDA_TIMER hartree(world, "apply perturbed hartree potential...");
	vecfuncT gamma = apply_hartree_potential(xfunction.x);
	hartree.info(debug_);

	TDA_TIMER exchange(world, "apply hf-exchange potential kernel...");\
	for (std::size_t p = 0; p < xfunction.x.size(); p++) {

		vecfuncT x_Ppi = mul(world, xfunction.x, exchange_intermediate_[p]);
		for (std::size_t i = 0; i < xfunction.x.size(); i++)
			gamma[p] -= x_Ppi[i];

		// Project out occupied space (can not use the function here, because of const ... and dont want to give it up)
		for (size_t p = 0; p < active_mo_.size(); ++p)
			gamma[p] -= rho0(gamma[p]);
	}
	exchange.info(debug_);

	return gamma;
}

vecfuncT TDA::apply_gamma_dft(const xfunction &xfunction) const {

	TDA_TIMER hartree(world, "apply perturbed hartree potential...");
	vecfuncT gamma = apply_hartree_potential(xfunction.x);
	hartree.info(debug_);

	TDA_TIMER rhoprime(world, "make perturbed density...");

	// Make the perturbed density for closed shell molecules
	real_function_3d perturbed_density = real_factory_3d(world);
	for (size_t i = 0; i < active_mo_.size(); i++) {
		perturbed_density += 2.0 * xfunction.x[i] * active_mo_[i];
	}
	rhoprime.info(debug_);
	//
	//	TDA_TIMER applyit(world,"apply vxc...");
	//
	// Get the perturbed xc potential from the dft class
	real_function_3d vxc = xclib_interface_.convolution_with_kernel(
			perturbed_density);

	for (size_t i = 0; i < gamma.size(); i++) {
		gamma[i] += vxc * active_mo_[i];
	}

	// Alternative way (more expensive, but avoid the unprecise kernel)
	// for small test molecules this seems to bring no improvement
	// when using this:
	// 1.comment out the line before
	// 2.return add(world,gamma,gamma2)
	// 3. dont forget to project out occupied space also from gamma2 (below here)
	// THIS DOES NOT WORK FOR GGA
	//	vecfuncT gamma2=xclib_interface_.apply_kernel(xfunction.x);
	//	for (int p=0; p<active_mo_.size(); ++p) gamma2[p] -= rho0(gamma2[p]);

	// project out occupied space
	for (size_t p = 0; p < active_mo_.size(); ++p)
		gamma[p] -= rho0(gamma[p]);

	//applyit.info(debug_);
	plot_vecfunction(gamma, "gamma", plot_);
	//return add(world,gamma,gamma2);
	truncate(world, gamma,truncate_thresh_);
	return gamma;
}

vecfuncT TDA::apply_hartree_potential(const vecfuncT &x) const {

	// Make the perturbed density (Factor 2 is for closed shell)
	real_function_3d perturbed_density = real_factory_3d(world);
	for (size_t i = 0; i < x.size(); i++) {
		perturbed_density += 2.0 * x[i] * active_mo_[i];
	}

	real_convolution_3d J = CoulombOperator(world, lo, bsh_eps_);
	real_function_3d Jrhoprime = J(perturbed_density);

	return mul(world, Jrhoprime, active_mo_);

}

vecfuncT TDA::get_V0(const vecfuncT& x) const {

	real_function_3d rho = density_;

	// the local potential V^0 of Eq. (4)
	real_function_3d coulomb;
	real_function_3d vlocal = get_calc().potentialmanager->vnuclear()
			+ get_coulomb_potential();

	// make the potential for V0*xp
	vecfuncT Vx = mul(world, vlocal, x);

	// and the exchange potential is K xp
	vecfuncT Kx;
	if (not dft_)
		Kx = get_calc().apply_hf_exchange(world, get_calc().aocc, mos_, x);
	if (dft_) {
		real_function_3d vxc = xclib_interface_.get_unperturbed_vxc();
		// Shift the unperturbed potential down
		vxc = vxc + shift_;
		Kx = mul(world, vxc, x);
	}

	// sum up: V0 xp = V_loc xp - K xp // this is 2J - K (factor 2 is included in get_coulomb_potential())
	vecfuncT V0;
	if (not dft_)
		V0 = sub(world, Vx, Kx);
	if (dft_)
		V0 = add(world, Vx, Kx);

	return V0;

}

real_function_3d TDA::get_coulomb_potential() const {
	MADNESS_ASSERT(get_calc().param.spin_restricted);
	if (coulomb_.is_initialized())
		return copy(coulomb_);
	functionT rho = get_calc().make_density(world, get_calc().aocc, mos_).scale(
			2.0); // here is the factor 2
	coulomb_ = get_calc().make_coulomb_potential(rho);
	return copy(coulomb_);
}

std::vector<vecfuncT> TDA::make_exchange_intermediate() const {
	const vecfuncT & active_mo = active_mo_;
	const vecfuncT & amo = active_mo_;
	// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr<
			real_convolution_3d>(CoulombOperatorPtr(world, lo, bsh_eps_));

	std::vector<vecfuncT> intermediate(amo.size());

	for (std::size_t p = 0; p < active_mo.size(); ++p) {
		intermediate[p] = apply(world, (*poisson),
				mul(world, active_mo[p], amo));
	}
	return intermediate;
}

void TDA::plot_vecfunction(const vecfuncT &x, std::string msg = "name_",
		bool plot) const {
	if (plot == true) {
		TDA_TIMER plot(world, "plotting...");
		for (size_t i = 0; i < x.size(); i++) {
			plot_plane(world, x[i], msg + stringify(i));
		}
		plot.info();
	}
}

bool TDA::check_convergence(xfunctionsT &xfunctions) {
	// assign convergence to every excitation
	std::cout<< std::setw(40) << "converged excitations" << " : ";
	for (size_t i = 0; i < xfunctions.size(); i++) {
		if(xfunctions[i].iterations == 0) return false;

		bool push_to_converged = false;
		bool replace = false;
		if(xfunctions[i].error.back() < dconv_
				and fabs(xfunctions[i].delta.back()) < econv_) {
			xfunctions[i].converged = true;
			push_to_converged = true;
		}
		else if(xfunctions[i].iterations>iter_max_ and fabs(xfunctions[i].delta.back()) < econv_) {
			push_to_converged = true;
			xfunctions[i].converged = false;
		}
		else if(xfunctions[i].iterations>iter_max_) replace = true;

		if (push_to_converged) {
			std::cout << i;
			//xfunctions[i].converged = true;
			xfunctions[i].current_residuals.clear();
			xfunctions[i].Vx.clear();
			//if(kain_) kain_solvers.reduce_subspace(i);
			converged_xfunctions_.push_back(xfunctions[i]);
			//will be substituted with guess function now//xfunctions.erase(xfunctions.begin()+i);

			// replace the converged xfunctions with the guess function it once was
			//(there should be no sorting on the way, e.g in the fock orthonormalization)
			xfunctionsT new_xfunctions;
			initialize(new_xfunctions);
			xfunctions[i]=new_xfunctions[i];
			// be shure
			xfunctions[i].x = copy(world,new_xfunctions[i].x);
			xfunctions[i].omega = guess_omega_;
			xfunctions[i].error.push_back(100);
			xfunctions[i].delta.push_back(100);
			xfunctions[i].expectation_value.push_back(guess_omega_);
			new_xfunctions.clear();
			if(kain_) {
				// all subspaces have to be erased because of the fock transformation
				for(size_t j=0;j<xfunctions.size();j++) kain_solvers.erase_subspace(j);
			}
		} else {
			// sort out xfunctions that iterated to long and showed no progress
			if(replace) {
				if(xfunctions[i].error.back() > 5.e-2) {
					xfunctionsT new_xfunctions;
					initialize(new_xfunctions);
					xfunctions[i]=new_xfunctions[i];
					// be shure
					xfunctions[i].x = copy(world,new_xfunctions[i].x);
					xfunctions[i].omega = guess_omega_;
					xfunctions[i].error.push_back(100);
					xfunctions[i].delta.push_back(100);
					xfunctions[i].expectation_value.push_back(guess_omega_);
					new_xfunctions.clear();
					if(kain_) {
						// all subspaces have to be erased because of the fock transformation
						for(size_t j=0;j<xfunctions.size();j++) kain_solvers.erase_subspace(j);
					}
				}
			}
		}

	}std::cout<<std::endl;

	// Check overall convergence
	bool converged = true;
	for (size_t i = 0; i < xfunctions.size(); i++) {
		if (xfunctions[i].converged == false) {
			converged = false;
			return false;
		}
	}
	return converged;
}

void TDA::print_performance(const xfunctionsT &xfunctions,const std::string prename) const {
	for (size_t i = 0; i < xfunctions.size(); i++) {
		std::fstream expv;
		expv.open("excitation_" + stringify(i) + "_expectation_value",
				std::ios::out);
		for (size_t k = 0; k < xfunctions[i].expectation_value.size(); k++) {
			expv << xfunctions[i].expectation_value[k] << endl;}
		expv.close();
		std::fstream err;
		err.open("excitation_" + stringify(i) + "_error", std::ios::out);
		for (size_t l = 0; l < xfunctions[i].error.size(); l++) {
			err << xfunctions[i].error[l] << endl;}
		err.close();
	}
	// Print results in TeX format
	std::fstream results;
	results.open(prename+"results.tex", std::ios::out);
	results << "% bsh_eps: " << bsh_eps_ << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << "\n";
	results << "\\begin{tabular}{lll}"<< "\n";
	results << "\\toprule" << " \\\\ \n";
	results << "\\multicolumn{1}{c}{omega}" << "&" << "\\multicolumn{1}{c}{error}"<< " \\\\ " << " \n";
	results << "\\midrule" << " \\\\ \n";
	for (size_t i = 0; i < xfunctions.size(); i++) {
		results << "\\num{" << std::setprecision(10) << xfunctions[i].omega
		<< "} & " << std::scientific << std::setprecision(1) << xfunctions[i].error.back() << " \\\\" << "\n";
	}
	results << "\\bottomrule" << " \\\\ \n";
	results << "\\end{tabular} \n";
	results.close();
}

void TDA::truncate_xfunctions(xfunctionsT &xfunctions) {
	for (size_t k = 0; k < xfunctions.size(); k++) {
		for (size_t i = 0; i < xfunctions[k].x.size(); i++) {
			xfunctions[k].x[i].truncate(truncate_thresh_);
			// The potential or residual vectors can be empty, not failsafe for the x-vector because it should never be empty
			if (not xfunctions[k].Vx.empty())
				xfunctions[k].Vx[i].truncate(truncate_thresh_);
			if (not xfunctions[k].current_residuals.empty())
				xfunctions[k].current_residuals[i].truncate(truncate_thresh_);
		}
	}
}

bool TDA::load_vecfunction(World& world, vecfuncT &v,
		const std::string filename_end) const {
	std::string filename = "vecfunc_" + filename_end;
	bool exists = archive::ParallelInputArchive::exists(world,
			filename.c_str());
	if (not exists)
		return false;

	archive::ParallelInputArchive ar(world, filename.c_str(), 1);
	for (std::size_t i = 0; i < v.size(); ++i) {
		ar & v[i];
		v[i].set_thresh(FunctionDefaults<3>::get_thresh());
	}
	return true;
}

void TDA::save_vecfunction(World& world, const vecfuncT &v,
		const std::string filename_end) const {
	std::string filename = "vecfunc_" + filename_end;
	archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
	for (std::size_t i = 0; i < v.size(); ++i)
		ar & v[i];
}

double TDA::oscillator_strength_length(const xfunction& xfunction) const {
	Tensor<double> mu_if(3);
	for (int idim=0; idim<3; idim++) {
		real_function_3d ri = real_factory_3d(world).functor2(xyz(idim));
		vecfuncT amo_times_x=mul(world,ri,active_mo_);
		Tensor<double> a=inner(world,amo_times_x,xfunction.x);
		mu_if(idim)=a.sum();
	}
	const double f= 2.0/3.0 * xfunction.omega * mu_if.sumsq() * 2.0;
	return f;
}

double TDA::oscillator_strength_velocity(const xfunction& root) const {
	Tensor<double> p_if(3);
	// compute the derivatives of the MOs in all 3 directions
	for (int idim=0; idim<3; idim++) {
		real_derivative_3d D = free_space_derivative<double,3>(world, idim);
		vecfuncT Damo=apply(world,D,active_mo_);
		Tensor<double> a=inner(world,Damo,root.x);
		p_if(idim)=a.sum();
	}
	const double f= 2.0/(3.0 * root.omega) * p_if.sumsq() * 2.0;
	return f;
}

void TDA::analyze(const xfunctionsT& roots) const {

	const size_t noct=active_mo_.size();

	std::vector<xfunction>::const_iterator it;
	int iroot=0;
	for (it=roots.begin(); it!=roots.end(); ++it, ++iroot) {
		std::vector<double> norms=norm2s(world,it->x);

		// compute the oscillator strengths
		double osl=this->oscillator_strength_length(*it);
		double osv=this->oscillator_strength_velocity(*it);

		std::cout << std::scientific << std::setprecision(10) << std::setw(20);
		if (world.rank()==0) {
			std::cout.width(10); std::cout.precision(8);
			print("excitation energy for root ",iroot,": ",it->omega);
			print("oscillator strength (length)    ", osl);
			print("oscillator strength (velocity)  ", osv);

			// print out the most important amplitudes
			print("\n  dominant contributions ");

			for (std::size_t p=0; p<noct; ++p) {
				const double amplitude=norms[p]*norms[p];
				if (amplitude > 0.1) {
					std::cout << "  norm(x_"<<p<<") **2  ";
					std::cout.width(10); std::cout.precision(6);
					std::cout << amplitude << std::endl;
				}
			}
		}
	}
}

