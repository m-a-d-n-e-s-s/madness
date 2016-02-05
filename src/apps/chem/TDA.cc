/*
 * TDA.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: kottmanj
 */

#include "TDA.h"

//#include <chem/nemo.h>
//#include <madness/mra/mra.h>
#include <cmath>
#include "../../madness/mra/derivative.h"
#include "../../madness/mra/funcdefaults.h"
#include "../../madness/mra/funcplot.h"
#include "../../madness/mra/function_interface.h"
#include "../../madness/mra/functypedefs.h"
#include "../../madness/mra/vmra1.h"
#include "../../madness/tensor/slice.h"
#include "../../madness/tensor/srconf.h"
#include "../../madness/tensor/tensor.h"
#include "../../madness/tensor/tensor_lapack.h"
#include "../../madness/world/madness_exception.h"
#include "../../madness/world/print.h"
#include "../../madness/world/timers.h"
#include "../../madness/world/world.h"
#include "../../madness/world/worldgop.h"
#include "SCFOperators.h"

using namespace madness;

/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};

void TDA::solve_guess(xfunctionsT &xfunctions) {
	output_section("SOLVE GUESS START");
	if(world.rank()==0) std::cout << "Solving with thresh " << FunctionDefaults<3>::get_thresh() << std::endl;
	if(world.rank()==0) std::cout << "dconv=" << guess_dconv_ << "\neconv=" << guess_econv_ << std::endl;
	converged_xfunctions_.clear();
	if(guess_ =="koala"){
		guess_koala(xfunctions);
		return;
	}

	plot_vecfunction(active_mo_, "active_mo_");

	// Create the excitation_function vector and initialize
	if(world.rank()==0)std::cout << "\n\n---Start Initialize Guess Functions---" << "\n\n " << std::endl;
	TDA_TIMER init(world,"\nfinished to initialize guess excitations ...");
	// if xfunctions were read in before then xfunctions.empty() will be false
	if(xfunctions.empty())initialize(xfunctions);
	for(size_t i=0;i<10;i++){
		output("Guess Iteration " + stringify(i));
		converged_xfunctions_.clear();
		iterate_guess(xfunctions);
		std::sort(converged_xfunctions_.begin(),converged_xfunctions_.end());
		size_t convergedc = 0;
		for(auto x:converged_xfunctions_){
			if(x.converged) convergedc++;
		}
		if(convergedc>=excitations_) break;
		else xfunctions = converged_xfunctions_;
	}
	output("\nEnd Initialize Guess Functions\n");
	init.info();

	// now sort the pre-converged xfunctions
	std::sort(converged_xfunctions_.begin(),converged_xfunctions_.end());
	if(converged_xfunctions_.size()>excitations_)converged_xfunctions_.erase(converged_xfunctions_.begin()+excitations_,converged_xfunctions_.end());

}

void TDA::solve(xfunctionsT &xfunctions) {
	output_section("SOLVE START");
	if(world.rank()==0) std::cout << "Solving with thresh " << FunctionDefaults<3>::get_thresh() << std::endl;
	if(world.rank()==0) std::cout << "dconv=" << dconv_ << "\neconv=" << econv_ << std::endl;
	converged_xfunctions_.clear();
	// Iterate till convergence is reached
	for(size_t iter=0;iter<300;iter++){
		output("\nIteration " + stringify(iter));
		converged_xfunctions_.clear();
		iterate(xfunctions);
		size_t convergedc=0;
		for(auto x:converged_xfunctions_){
			if(x.converged) convergedc++;
		}
		if(convergedc>=excitations_) break;
		else{
			xfunctions = converged_xfunctions_;
		}

	}

	// print final result
	print_status(xfunctions);
	print_performance(xfunctions,"pre-");

}

void TDA::solve_sequential(xfunctionsT &xfunctions) {
        converged_xfunctions_.clear();
	// manual now
	bool kain=true;

	if(world.rank()==0) std::cout << "\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE_SEQUENTIAL START " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n"
			"The first " << excitations_ << " of the following xfunctions will be solved sequentially "<< std::endl;

	if(xfunctions.empty()) initialize(xfunctions);

	std::sort(xfunctions.begin(),xfunctions.end());
	xfunctions.erase(xfunctions.begin()+excitations_,xfunctions.end());
	print_status(xfunctions);

	print("\n\n\n\n-------------------------------------------------------");
	print("BEGINNING THE FINAL ITERATIONS TO AN ACCURACY OF ... no kain right now", hard_dconv_);
	print("-------------------------------------------------------\n\n\n\n");

	// Failsafe for comming fock orthonormalization
	for(size_t i=0;i<xfunctions.size();i++)xfunctions[i].smooth_potential.clear();

	orthonormalize_fock(xfunctions);

	size_t max = xfunctions.size();

	print_status(xfunctions);

	// The given xfunctions should be sorted by energy in ascending order:
	for (size_t iroot = 0; iroot < max; iroot++) {
		print("\n\n-----xfunction ", iroot, " ------\n\n");
		xfunction current_root = xfunctions[iroot];
		current_root.converged = false;

		// make kain solver
		sequential_kain_solver kain_solver(TDA_allocator(world,active_mo_.size()),true);
		kain_solver.set_maxsub(kain_subspace_);

		for(size_t iter =0;iter<iter_max_+1;iter++){
			TDA_TIMER seq_iter(world,"\nSequential iteration "+ stringify(iter) + " time:");
			normalize(current_root);
			project_out_occupied_space(current_root.x);
			xfunctionsT carrier(1,current_root);
			project_out_converged_xfunctions(carrier);
			current_root = carrier.front();
			normalize(current_root);
			current_root.smooth_potential = apply_smooth_potential(current_root);
			current_root.expectation_value.push_back(expectation_value(current_root,current_root.smooth_potential));
			truncate(world,current_root.x);
			truncate(world,current_root.smooth_potential);
			update_energy(current_root);

			if(not kain){
				current_root.x = iterate_one(current_root);
			}else{
				vecfuncT updated_x = iterate_one(current_root);
				// The kain structure needs the residual like this: old - new
				vecfuncT residual = sub(world,current_root.x,updated_x);
				xfunction new_x = kain_solver.update(xfunction(world,current_root.x),xfunction(world,residual));
				current_root.x = new_x.x;
			}
			seq_iter.info();
			print_xfunction(current_root);


			if(fabs(current_root.error.back())<hard_dconv_ and fabs(current_root.delta.back())<hard_econv_){
				current_root.converged =true;
				converged_xfunctions_.push_back(current_root);
				if(world.rank()==0) std::cout << "\nexcitation " << iroot << " converged!!!\n" << std::endl;
				break;
			}else if(iter==iter_max_){
				current_root.converged=false;
				converged_xfunctions_.push_back(current_root);
				if(world.rank()==0) std::cout << "\nexcitation " << iroot << " converged not ...\n" << std::endl;
				break;
			}

		}
	}


	plot_vecfunction(active_mo_, "active_mo_");
	// plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		plot_vecfunction(converged_xfunctions_[i].x,"final_excitation_"+stringify(i),true);
	}

	// calculate densities and plot
	real_function_3d rho = CCOPS_.make_density(active_mo_,active_mo_);
	plot_plane(world,rho,"rho");

	for(size_t p=0;p<converged_xfunctions_.size();p++){
		real_function_3d rhoprime = CCOPS_.make_density(active_mo_,converged_xfunctions_[p].x);
		plot_plane(world,rhoprime,"rho_pert_"+stringify(p));
	}

	print("Final result :");
	print_status(converged_xfunctions_);
	analyze(converged_xfunctions_);
	print_performance(converged_xfunctions_,"final-");
	print_status(xfunctions);

}

void TDA::print_status(const xfunctionsT & xfunctions) const {
	if(world.rank()==0){
		if(compute_virtuals_) std::cout << "\nVirtual orbitals: Energy is omega + epsilon\n";
		std::cout << "\n" <<std::setw(5) << " " << std::setw(20) << "omega" << std::setw(20) << "delta" << std::setw(20)
		<< "error"<<std::setw(20)
		<<"expv" << std::setw(7) <<"iter"<< std::setw(7)<< "conv" << std::endl;
		for(size_t i=0;i<converged_xfunctions_.size();i++) print_xfunction(converged_xfunctions_[i]);
		std::cout <<"--------"<< std::endl;
		for(size_t i=0;i<xfunctions.size();i++)print_xfunction(xfunctions[i]);
		std::cout << "\n" << std::endl;
	}
}

void TDA::print_xfunction(const xfunction &x) const {
	if(world.rank()==0){
		std::cout << std::setw(5) << x.number;
		std::cout << std::scientific << std::setprecision(10) << std::setw(20) << x.omega << std::setw(20)<< x.delta.back()
																																	<< std::setw(20)<< x.error.back()<< std::setw(20) << x.expectation_value.back();
		std::cout << std::fixed <<std::setw(7)<< x.iterations << "   " << std::setw(7)<<x.converged << std::endl;
	}
}

void TDA::initialize(xfunctionsT & xfunctions){
	if(world.rank()==0) std::cout << "\nDemanded guess is " << guess_ <<"\n"<< std::endl;
	std::vector<std::string> exop_strings;
	if(guess_ == "custom"){
		exop_strings = custom_exops_;
	}
	else{
		polynomial_guess guess_structure;
		exop_strings =  guess_structure.make_predefined_guess_strings(guess_);
	}
	if(exop_strings.empty()) MADNESS_EXCEPTION("ERROR in intializing xfunctions, no guess strings were created ... unknown keyword ?",1);
	for(size_t i=0;i<exop_strings.size();i++){
		xfunction tmp(world);
		tmp.number=i;
		tmp.guess_excitation_operator = exop_strings[i];
		tmp.omega = guess_omega_;
		xfunctions.push_back(tmp);
	}
	for(size_t i=0;i<xfunctions.size();i++){
		xfunctions[i].x = make_guess_vector(xfunctions[i].guess_excitation_operator);
	}
	make_big_fock_guess(xfunctions);
	normalize(xfunctions);
	for(size_t i=0;i<xfunctions.size();i++){
		project_out_occupied_space(xfunctions[i].x);
		normalize(xfunctions[i]);
		plot_vecfunction(xfunctions[i].x,"guess_excitation_" + stringify(i) + "_",true);
	}

}

void TDA::make_big_fock_guess(xfunctionsT &xfunctions)const{
	if(world.rank()==0) std::cout << "\n\nMaking big perturbed fock guess...\n\n" << std::endl;

	TDA_TIMER big_ortho(world,"Big Perturbed Fock Matrix: ");
	normalize(xfunctions);
	truncate_xfunctions(xfunctions);

	// Diagonalize the big fock guess
	{
		Tensor<double> overlap(xfunctions.size(), xfunctions.size());
		for (size_t p = 0; p < xfunctions.size(); p++) {
			for (size_t k = 0; k < xfunctions.size(); k++) {
				overlap(p, k) = CCOPS_.make_inner_product(xfunctions[p].x,xfunctions[k].x);
			}
		}

		Tensor<double> F = make_perturbed_fock_matrix(xfunctions);
		if(parameters.debug) std::cout<< "guess overlap matrix is\n" << overlap << "\n";
		if(parameters.debug) std::cout<<"guess pert. Fock Matrix is\n" << F << "\n";
		Tensor<double> U, evals, dummy(xfunctions.size());
		U = get_nemo().get_calc() -> get_fock_transformation(world, overlap, F, evals, dummy,
				1.5 * econv_);
		if(parameters.debug) std::cout<<"Transformation Matrix is\n" << U << "\n";
		std::vector<vecfuncT> old_x;
		std::vector<std::string> old_exop_strings;
		for (size_t i = 0; i < xfunctions.size(); i++) {
			old_x.push_back(xfunctions[i].x);
			old_exop_strings.push_back(xfunctions[i].guess_excitation_operator);
		}
		std::vector<vecfuncT> new_x = transform_vecfunctions(old_x, U);
		std::vector<std::string> new_exop_strings(new_x.size());

		for (size_t i = 0; i < xfunctions.size(); i++) {
			size_t counter =0;
			for (size_t j = 0; j < xfunctions.size(); ++j) {
				if(fabs(U(i,j)) > FunctionDefaults<3>::get_thresh()*10.0){
					if(counter != 0) new_exop_strings[i] += " , ";
					new_exop_strings[i] +=old_exop_strings[j] + " c " + stringify(U(i,j));
					counter++;
				}
			}
		}
		for (size_t i = 0; i < xfunctions.size(); i++) {
			xfunctions[i].x = new_x[i];
			xfunctions[i].expectation_value.push_back(evals(i));
			xfunctions[i].guess_excitation_operator = new_exop_strings[i];
			xfunctions[i].number = i;
		}
	}
	std::sort(xfunctions.begin(),xfunctions.end());
	// safe memory
	if(xfunctions.size()>guess_excitations_) xfunctions.erase(xfunctions.begin()+guess_excitations_,xfunctions.end());
	big_ortho.info();
	if(world.rank()==0) std::cout << "\nthe following guess functions have been created:\n " << std::endl;
	print_status(xfunctions);
	if(world.rank()==0 and parameters.debug){
		std::cout << "\nCorresponding excitation operators are:\n" << std::endl;
		for(size_t i=0;i<xfunctions.size();i++){
			xfunctions[i].number =i ;
			std::cout << std::setprecision(2) << xfunctions[i].guess_excitation_operator << std::endl;
		}
	}

}

vecfuncT TDA::make_guess_vector(const std::string &input_string)const{
	if(world.rank()==0) std::cout << "\nMaking guess excitation with excitation operator \n " << input_string << std::endl;
	std::shared_ptr<FunctionFunctorInterface<double, 3> > exop_functor(
			new polynomial_exop_functor(input_string));
	real_function_3d exop = real_factory_3d(world).functor(exop_functor);
	vecfuncT guess = mul(world,exop,active_mos_for_guess_calculation_);
	project_out_occupied_space(guess);
	return guess;
}


void TDA::guess_koala(xfunctionsT &roots)const{

	// for convenience
	const std::size_t nmo=orbital_energies_.size();

	// first we need to determine the rotation of the external orbitals
	// to the MRA orbitals, so that we can subsequently rotate the guess
	Tensor<double> guess_phases_;
	vecfuncT koala_mo;

	// read koala's orbitals from disk
	for (std::size_t i=parameters.freeze; i<nmo; ++i) {
		real_function_3d x_i=real_factory_3d(world).empty();
		const std::string valuefile="grid.koala.orbital"+stringify(i);
		x_i.get_impl()->read_grid2<3>(valuefile,functorT());
		koala_mo.push_back(x_i);
	}
	if(koala_mo.size()!=active_mo_.size()) MADNESS_EXCEPTION("ERROR in Koala guess: not the same number of Koala mos and active_mos of MRA",1);
	// this is the transformation matrix for the rotation
	guess_phases_=matrix_inner(world,koala_mo,mos_);
	guess_phases_=guess_phases_(_,Slice(parameters.freeze,nmo-1));

	// compute the inverse of the overlap matrix
	Tensor<double> S=(guess_phases_+transpose(guess_phases_)).scale(0.5);
	Tensor<double> U, eval;
	syev(S,U,eval);
	Tensor<double> Sinv=copy(U);
	for (int i=0; i<U.dim(0); ++i) {
		for (int j=0; j<U.dim(1); ++j) {
			Sinv(i,j)/=eval(j);
		}
	}
	Sinv=madness::inner(Sinv,U,-1,-1);

	vecfuncT transformed_koala_mos = transform(world,koala_mo,Sinv);
	Tensor<double> check = matrix_inner(world, transformed_koala_mos,active_mo_);
	if(world.rank()==0 and active_mo_.size()<6) std::cout << "Overlap between transformed Koala MOs and MRA MOs \n" << check << std::endl;
	double check_2 = measure_offdiagonality(check,active_mo_.size());
	double size_tmp = (double) active_mo_.size();
	if(check_2 > size_tmp*FunctionDefaults<3>::get_thresh()){
		if(world.rank()==0) std::cout << "WARNING: Koala and MRA MOs are maybe rotated" << check << std::endl;
	}


	for(size_t iroot=0;iroot<guess_excitations_;iroot++){


		xfunction root(world);

		// read the actual external guess from file
		for (std::size_t i=parameters.freeze; i<nmo; ++i) {

			// this is the file where the guess is on disk
			const std::string valuefile="grid.koala.orbital"+stringify(i)
	    																							+".excitation"+stringify(iroot);
			real_function_3d x_i=real_factory_3d(world).empty();
			x_i.get_impl()->read_grid2<3>(valuefile,functorT());
			root.x.push_back(x_i);
		}

		// now rotate the active orbitals from the guess to conform with
		// the MRA orbitals
		// Sinv=Sinv(Slice(parameters.freeze,nmo-1),Slice(parameters.freeze,nmo-1));
		root.x=transform(world,root.x,Sinv);

		root.omega = guess_omega_;
		root.number = iroot;

		for(size_t i=0;i<roots.size();i++){
			for(size_t j=0;j<roots[i].x.size();j++){
				plot_plane(world,roots[i].x[j],"Koala_guess_root_"+stringify(i)+"_"+stringify(j));
			}
		}
		roots.push_back(root);
	}

	// Print information about the read roots
	if(world.rank()==0) std::cout << "Read " << roots.size() << " excitation vectors from the koala calculation" << std::endl;
	if(roots.empty()) MADNESS_EXCEPTION("ERROR: Koala guess demanded, but no KOALA results found ...",1);
}


void TDA::iterate_guess(xfunctionsT &xfunctions) {
	if(world.rank()==0)std::cout << "\n" << std::endl;
	if(world.rank()==0)std::cout << "---Start Guess Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,true);
	// set back iteration counter
	//for(size_t i=0;i<xfunctions.size();i++)xfunctions[i].iterations = 0;
	if(world.rank()==0)std::cout <<std::setw(100)<< "---End Guess Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate(xfunctionsT &xfunctions) {
	if(world.rank()==0)std::cout << "\n" << std::endl;
	if(world.rank()==0)std::cout << "---Start Main Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,false);
	if(world.rank()==0)std::cout <<std::setw(100)<< "---End Main Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate_all(xfunctionsT &all_xfunctions, bool guess) {

	// Bool checks if the perturbed fock matrix has been calculated
	//(if not the expencation value has to be calculated in the iterate_one routine)
	bool pert_fock_applied = false;

	// Restrict the number of parallel iterating guess functions
	for(size_t i=0;i<all_xfunctions.size();i++){
		all_xfunctions[i].current_residuals.clear();
		all_xfunctions[i].smooth_potential.clear();
		all_xfunctions[i].converged = false;
	}

	// make big fock diagonalization
	orthonormalize_fock(all_xfunctions);
	if(iterating_excitations_ > all_xfunctions.size()) iterating_excitations_=all_xfunctions.size();
	std::sort(all_xfunctions.begin(),all_xfunctions.end());
	print_status(all_xfunctions);
	xfunctionsT xfunctions(all_xfunctions.begin(),all_xfunctions.begin()+iterating_excitations_);
	xfunctionsT remaining_xfunctions(all_xfunctions.begin()+iterating_excitations_,all_xfunctions.end());

	size_t guess_iter_counter =1;
	std::vector<size_t> iteration_counters(xfunctions.size(),0);
	size_t maxiter = 100;
	if(guess) maxiter = 10;
	for(size_t i=0;i<maxiter;i++){
		TDA_TIMER iteration_time(world,"\nEnd of iteration " + stringify(i) +": ");

		{
			TDA_TIMER update_potentials(world,"Update smooth potentials: ");
			for(size_t i=0;i<xfunctions.size();i++) xfunctions[i].smooth_potential=apply_smooth_potential(xfunctions[i]);
			update_potentials.info();
		}{
			memory_information(xfunctions);
		}{
			TDA_TIMER normalization(world,"Normalization: ");
			normalize(xfunctions);
			normalization.info();
		}{
			TDA_TIMER orthonormalization(world,"Orthonormalization: ");
			pert_fock_applied= orthonormalize_fock(xfunctions);
			orthonormalization.info();
		}{
			TDA_TIMER truncation(world,"Truncate: ");
			truncate_xfunctions(xfunctions);
			truncation.info();
		}{
			memory_information(xfunctions);
		}{
			TDA_TIMER apply_bsh(world,"Apply Green's operator: ");
			if(world.rank()==0) std::cout << std::setw(40) << "update energy..." << " : ";
			for(size_t i=0;i<xfunctions.size();i++){
				// if the fock matrix has not been calculated we need to caluclate the expectation value
				if(not pert_fock_applied) xfunctions[i].expectation_value.push_back(expectation_value(xfunctions[i], xfunctions[i].smooth_potential));
				if(world.rank()==0) std::cout << " " <<  update_energy(xfunctions[i]) << " ";
				xfunctions[i].x=iterate_one(xfunctions[i]);
				xfunctions[i].smooth_potential.clear();
			}
			if(world.rank()==0)std::cout << std::endl;
			apply_bsh.info();
		}{
			TDA_TIMER normalization(world,"Normalization: ");
			normalize(xfunctions);
			normalization.info();
		}{
			TDA_TIMER truncation(world,"Truncate: ");
			truncate_xfunctions(xfunctions);
			truncation.info();
		}{
			check_convergence(xfunctions,guess);
			if(not guess){
				for(size_t ix=0;ix<xfunctions.size();ix++){

					if(xfunctions[ix].converged or iteration_counters[ix]>=solve_iter_){
						converged_xfunctions_.push_back(xfunctions[ix]);
						if(remaining_xfunctions.empty()){
							xfunctions.erase(xfunctions.begin()+ix);
							iteration_counters.erase(iteration_counters.begin()+ix);
						}
						else{
							xfunctions[ix] = remaining_xfunctions.front();
							remaining_xfunctions.erase(remaining_xfunctions.begin());
							iteration_counters[ix]=0;
						}
					}else iteration_counters[ix]++;
				}
			}else{
				if(guess_iter_counter==guess_iter_){
					for(size_t j=0;j<xfunctions.size();j++)converged_xfunctions_.push_back(xfunctions[j]);
					xfunctions.clear();
					if(remaining_xfunctions.empty()) break;
					else{
						for(size_t j=0;j<iterating_excitations_;j++){
							if(remaining_xfunctions.empty()) break;
							xfunctions.push_back(remaining_xfunctions.front());
							remaining_xfunctions.erase(remaining_xfunctions.begin());
						}
					}
					guess_iter_counter = 0;
				}
			}
			TDA_TIMER project_time(world,"Project out converged excitations and occupied space: ");
			project_out_converged_xfunctions(xfunctions);
			for(size_t i=0;i<xfunctions.size();i++) project_out_occupied_space(xfunctions[i].x);
			project_time.info();
			print_status(xfunctions);

			size_t break_condition = excitations_;
			if(guess) break_condition = guess_excitations_;
			if(converged_xfunctions_.size()>=break_condition){
				// push all current iterating roots which are not freshly re-initialized to the converged roots
				for(size_t i=0;i<xfunctions.size();i++){
					if(xfunctions[i].iterations>2) converged_xfunctions_.push_back(xfunctions[i]);
				}
				break;
			}
		}

		iteration_time.info();
		if(world.rank()==0) std::cout << "\n\n" << std::endl;
		guess_iter_counter++;
	}

}

vecfuncT TDA::iterate_one(xfunction & xfunction)const {
	xfunction.iterations += 1;


	// Add the nuclear potential
	TDA_TIMER VNUC(world,"Adding nuclear potential to xfunction ");
	vecfuncT Vpsi = add(world,apply_nuclear_potential(xfunction),xfunction.smooth_potential);
	plot_plane(world,Vpsi.back(),"Vpsi_homo");
	VNUC.info();

	double omega = xfunction.omega;
	if(Vpsi.empty()) MADNESS_EXCEPTION("ERROR in iterate_one function: Applied potential of xfunction is empty",1);

	if(not use_omega_for_bsh_){
	truncate(world, Vpsi); // no fence
	vecfuncT omegapsi = xfunction.x;
	scale(world,omegapsi,-omega);
	Vpsi = add(world,Vpsi,omegapsi);
	}
	truncate(world, Vpsi); // no fence
	scale(world, Vpsi, -2.0);



	std::vector<poperatorT> bsh(active_mo_.size());
	for (size_t p = 0; p < active_mo_.size(); p++) {
		double eps = active_eps(p);// + omega;
		if(use_omega_for_bsh_) eps += omega;
		if (eps > 0) {
			if (world.rank() == 0)
				print("bsh: warning: positive eigenvalue", p + parameters.freeze, eps);
			eps = active_eps(p) + guess_omega_;
		}
		bsh[p] = poperatorT(
				BSHOperatorPtr3D(world, sqrt(-2.0 * eps), parameters.lo, parameters.thresh_bsh_3D));
	}

	world.gop.fence();

	vecfuncT GVpsi = apply(world, bsh, Vpsi);

	// Residual as error estimation
	project_out_occupied_space(GVpsi);
	vecfuncT residual = sub(world, xfunction.x, GVpsi);

	double resinner = CCOPS_.make_inner_product(residual, residual);
	double error = sqrt(resinner);
	xfunction.error.push_back(error);

	// Calculate 2nd order update:
	// Inner product of Vpsi and the residual (Vspi is scaled to -2.0 --> multiply later with 0.5)
	double tmp = CCOPS_.make_inner_product(residual,Vpsi);

	// squared norm of GVpsi (Psi_tilde)
	double tmp2 = CCOPS_.make_inner_product(GVpsi,GVpsi);

	// Factor 0.5 removes the factor 2 from the scaling before
	xfunction.delta.push_back(0.5 * tmp / tmp2);

	// return Updated x-function
	return GVpsi;

}


std::string TDA::update_energy(xfunction &xfunction)const {
	double thresh = FunctionDefaults<3>::get_thresh();
	if(use_omega_for_bsh_){
	//failsafe: make shure the delta and expectation values vectors are not empty to avoid segmentation faults
	if(not xfunction.delta.empty() and not xfunction.expectation_value.empty() and not xfunction.error.empty()) {
		if(xfunction.expectation_value.back() < highest_excitation_) {
			if(fabs(xfunction.delta.back()) < thresh*20.0) {
				xfunction.omega +=xfunction.delta.back();
				return "(2nd)";
			} else {
				xfunction.omega = xfunction.expectation_value.back();
				return "(exp)";
			}
		} else {
			// get the highest converged energy, of no xfunction converged already use the guess_omega_ energy
			double setback_energy = guess_omega_;
			// set the last converged value of the same type of guess as the default
			//double new_omega = highest_excitation_*0.8;// default
			xfunction.omega = setback_energy;
			return "(-)";
		}

	}else return "(no energies)";

	}else{
		xfunction.omega = xfunction.expectation_value.back();
		return "(omega not used for BSH, no second order update possible -> use expectation value)";
	}
}

void TDA::normalize(xfunctionsT & xfunctions)const {
	for (size_t i = 0; i < xfunctions.size(); i++)
		normalize(xfunctions[i]);
}
void TDA::normalize(xfunction & xfunction)const {
	double self_overlap = CCOPS_.make_inner_product(xfunction.x, xfunction.x);
	const double norm = sqrt(self_overlap);
	scale(world, xfunction.x, 1.0 / norm);

}

void TDA::project_out_converged_xfunctions(xfunctionsT & xfunctions)const {
	for (size_t p = 0; p < xfunctions.size(); p++) {
		for (size_t k = 0; k < converged_xfunctions_.size(); k++) {
			double overlap = CCOPS_.make_inner_product(xfunctions[p].x,converged_xfunctions_[k].x);
			for (size_t i = 0; i < xfunctions[p].x.size(); i++) {
				xfunctions[p].x[i] -= overlap * converged_xfunctions_[k].x[i];
			}
		}
	}
}

bool TDA::orthonormalize_fock(xfunctionsT &xfunctions)const {
	normalize(xfunctions);

	if(xfunctions.size()<2){
		if(world.rank()==0) std::cout << "no Fock diagonalization just " << xfunctions.size() << " xfunctions" << std::endl;
		return false;
	}

	Tensor<double> overlap(xfunctions.size(), xfunctions.size());
	for (size_t p = 0; p < xfunctions.size(); p++) {
		for (size_t k = 0; k < xfunctions.size(); k++) {
			overlap(p, k)  = CCOPS_.make_inner_product(xfunctions[p].x,xfunctions[k].x);
		}
	}

	if(world.rank()==0) std::cout << "Overlap Matrix\n" << overlap <<"\n";

	Tensor<double> F(3L, xfunctions.size());
	if (dft_) {
		F = make_perturbed_fock_matrix(xfunctions);
	} else
		F = make_perturbed_fock_matrix(xfunctions);

	// Diagonalize the perturbed Fock matrix
	Tensor<double> U, evals;

	madness::Tensor<double> dummy(xfunctions.size());
	U = get_nemo().get_calc() -> get_fock_transformation(world, overlap, F, evals, dummy,
			1.5 * econv_);
	//}
	print("\n\n");
	if(world.rank()==0 and xfunctions.size()<6) std::cout << std::setprecision(3) << std::setw(40)<< "Transformation-Matrix-U \n" << U << std::endl;
	print("\n\n");

	//Prevent printout when expectation value is calculated
	if (xfunctions.size() > 1) {
		// Make an estimation how "strong" the xfunctions are linear combinated
		double offdiagonal = measure_offdiagonality(U, xfunctions.size());
		if(world.rank()==0) std::cout << std::setw(40) << "offdiagonal transformation part..." << " : " << fabs(offdiagonal) << std::endl;
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
	if(world.rank()==0) std::cout << std::setw(40) << "Transforming..." << " : xfunctions... ";

	// potentials
	{
		std::vector<vecfuncT> old_Vx;
		for (size_t i = 0; i < xfunctions.size(); i++) {
			old_Vx.push_back(xfunctions[i].smooth_potential);
		}
		std::vector<vecfuncT> new_Vx = transform_vecfunctions(old_Vx, U);
		for (size_t i = 0; i < xfunctions.size(); i++) {
			xfunctions[i].smooth_potential = new_Vx[i];
		}
		if(world.rank()==0)std::cout << "potentials... ";
	}

	if(world.rank()==0) std::cout << std::endl;

	return true;

}

double TDA::measure_offdiagonality(const madness::Tensor<double> &U,
		const size_t size) const {
	if(size ==0) return 0.0;
	if(size ==1) return 0.0;
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
		new_xfunctions[i] = zero_functions_compressed<double, 3>(world,
				xfunctions[i].size());
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

void TDA::project_out_occupied_space(vecfuncT &x)const {
	CCOPS_.Q(x);
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
		weighted_sum += (orbital_energies_(parameters.freeze + i) + shift_) * overlaps[i];
	}

	return value - weighted_sum;
}

double TDA::expectation_value(const xfunction &x, const vecfuncT &smooth_potential)const {
	// The part from the smooth potential
	double expv = CCOPS_.make_inner_product(x.x, smooth_potential);

	// The Nuclear part
	real_function_3d vnuc = get_nemo().get_calc() -> potentialmanager->vnuclear();
	vecfuncT vnuci = mul(world,vnuc,x.x);
	Tensor<double> vnuc_pot = inner(world, x.x, vnuci);
	double exp_vnuc = CCOPS_.make_inner_product(x.x,apply_nuclear_potential(x));
	expv += exp_vnuc;

	// The kinetic part
	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	for (int axis = 0; axis < 3; axis++) {
		const vecfuncT dx = apply(world, *(gradop[axis]), x.x);
		vecfuncT bra_dx;
		bra_dx = apply(world, *(gradop[axis]), mul(world,get_nemo().nuclear_correlation -> square(),x.x));
		Tensor<double> kin;
		kin= inner(world, bra_dx, dx);
		expv += 0.5 * kin.sum();
	}

	// The epsilon part to get the excitation energy
	for (size_t i = 0; i < x.x.size(); i++) {
		double overlap = CCOPS_.make_inner_product(x.x[i],x.x[i]);
		expv -= active_eps(i)* overlap;
	}
	return expv;
}

Tensor<double> TDA::make_perturbed_fock_matrix(
		const xfunctionsT &xfunctions) const {
	TDA_TIMER MAKE_FOCK_MATRIX_TIMER2(world,"Make new perturbed fock matrix");
	output("Make Perturbed Fock Matrix new version");
	TDA_TIMER VPART(world,"Potential Part");
	Tensor<double> new_F(xfunctions.size(),xfunctions.size());
	for(size_t q=0;q<xfunctions.size();q++){
		vecfuncT Vxq;
		if(xfunctions[q].smooth_potential.empty()){
			output("smooth potential needs to be recalculated");
			Vxq = add(world,apply_smooth_potential(xfunctions[q]),apply_nuclear_potential(xfunctions[q]));
		}else Vxq = add(world,xfunctions[q].smooth_potential,apply_nuclear_potential(xfunctions[q]));
		for(size_t p=0;p<xfunctions.size();p++){
			new_F(p,q) = CCOPS_.make_inner_product(xfunctions[p].x,Vxq);
		}}
	VPART.info();
	TDA_TIMER TPART(world,"Kinetic Part");
	// pack the xfunctions in a vector
	std::vector<vecfuncT> XVEC;
//	for(size_t q=0;q<xfunctions.size();q++){
//		for(size_t p=0;p<xfunctions.size();p++){
//			//new_F(p,q) += CCOPS_.get_matrix_element_kinetic_energy(xfunctions[p].x,xfunctions[q].x);
//		}}
	for(size_t i=0;i<xfunctions.size();i++) XVEC.push_back(xfunctions[i].x);
	new_F += CCOPS_.get_matrix_kinetic(XVEC);
	TPART.info();
	TDA_TIMER EPART(world,"Orbital Energy Part");
	for(size_t q=0;q<xfunctions.size();q++){
		for(size_t p=0;p<xfunctions.size();p++){
			// substract the weighted orbital energies
			for(size_t i=0;i<xfunctions[p].x.size();i++){
				new_F(p,q) -= active_eps(i)*CCOPS_.make_inner_product(xfunctions[p].x[i],xfunctions[q].x[i]);
			}
		}
	}
	EPART.info();
	MAKE_FOCK_MATRIX_TIMER2.info();

	return new_F;
}



// The smooth potential is the potential without the nuclear potential
// The nuclear potential has to be calculated sepparately
vecfuncT TDA::apply_smooth_potential(const xfunction&xfunction)const{
	//std::cout << "\nPotential for virtuals: J-K , Vnuc is added later\n";
	vecfuncT smooth_potential;
	if(compute_virtuals_){
		if(dft_) smooth_potential = CCOPS_.KS_residue_closed_shell(xfunction.x);
		else smooth_potential=CCOPS_.fock_residue_closed_shell(xfunction.x);
	}else{
		// for now
		if(dft_) smooth_potential = CCOPS_.get_TDA_potential_singlet(xfunction.x);
		else smooth_potential = CCOPS_.get_CIS_potential_singlet(xfunction.x);
	}
	truncate(world,smooth_potential);
	plot_plane(world,smooth_potential.back(),"smooth_potential_homo");
	return smooth_potential;
}

vecfuncT TDA::apply_nuclear_potential(const xfunction &xfunction) const{
	Nuclear U(world,&get_nemo());
	vecfuncT Ux = U(xfunction.x);
	truncate(world,Ux);
	return Ux;
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

void TDA::check_convergence(xfunctionsT &xfunctions, const bool guess)const{
	// get the convergence criteria
	double econv = econv_;
	double dconv = dconv_;
	if(guess){
		econv = guess_econv_;
		dconv = guess_dconv_;
	}
	// check convergence for every root
	for(size_t i=0;i<xfunctions.size();i++){
		if(fabs(xfunctions[i].error.back()) < dconv and fabs(xfunctions[i].delta.back()) < econv){
			xfunctions[i].converged = true;
		}else xfunctions[i].converged = false;
	}
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
	results << "% bsh_eps: " <<  parameters.thresh_bsh_3D << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << " \t \t \n";
	results << "\\begin{tabular}{llll}"<< "\n";
	results << "\\toprule" << "  \n";
	results << "\\multicolumn{1}{c}{$\\omega$}" << "&" << "\\multicolumn{1}{c}{error}" << "& \\multicolumn{1}{c}{$\\Delta$}" <<"\\multicolumn{1}{c}{iter}"" \\\\ " << " \n";
	results << "\\midrule" << " \n";
	for (size_t i = 0; i < xfunctions.size(); i++) {
		results << "\\num{" << std::setprecision(10) << xfunctions[i].omega
				<< "} & " << std::scientific << std::setprecision(0)
		<< xfunctions[i].error.back()
		<< " & "  <<xfunctions[i].delta.back()
		<< " & (" <<xfunctions[i].iterations << ")" << " \\\\" << "\n";
	}
	results << "\\bottomrule" << "  \n";
	results << "\\end{tabular} \n";
	results.close();

	// Print results with oscillator strength in Tex format
	std::fstream results2;
	results2.open(prename+"results_full.tex", std::ios::out);
	results2 << "% bsh_eps: " <<  parameters.thresh_bsh_3D << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << "\n";
	results2 << "\\begin{tabular}{ll}"<< "\n";
	results2 << "\\toprule" << " \n";
	for (size_t i = 0; i < xfunctions.size(); i++) {
		results2 << "Excitation " << i+1 << " & " << "\\multicolumn{1}{c}{MRA}" << " \\\\ \n ";
		results2 << "\\midrule" << "  \n";
		results2 << "Excitation energy "<< " & " << "\\num{" << std::setprecision(10) << xfunctions[i].omega << "}" << " \\\\ \n ";
		results2 << "f\\textsubscript{osc}(length) "<< " & " << "\\num{" << std::setprecision(10) << xfunctions[i].f_length << "}" << " \\\\ \n ";
		results2 << "f\\textsubscript{osc}(velocity) "<< " & " << "\\num{" << std::setprecision(10) << xfunctions[i].f_velocity << "}" << " \\\\ \n ";
		results2 << "\\midrule" << " \\\\ \n";
	}
	results2 << "\\bottomrule" << " \n";
	results2 << "\\end{tabular} \n";
	results2 << "\\multicolumn{2}{l}{ bsh_eps: " <<  parameters.thresh_bsh_3D << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << "} \\\\ \n" ;
	results2.close();
}

void TDA::truncate_xfunctions(xfunctionsT &xfunctions)const {
	for (size_t k = 0; k < xfunctions.size(); k++) {
		for (size_t i = 0; i < xfunctions[k].x.size(); i++) {
			xfunctions[k].x[i].truncate();
			// The potential or residual vectors can be empty, not failsafe for the x-vector because it should never be empty
			if (not xfunctions[k].smooth_potential.empty())
				xfunctions[k].smooth_potential[i].truncate();
			if (not xfunctions[k].current_residuals.empty())
				xfunctions[k].current_residuals[i].truncate();
		}
	}
}

double TDA::oscillator_strength_length(const xfunction& xfunction) const {
	Tensor<double> mu_if(3);
	for (int idim=0; idim<3; idim++) {
		real_function_3d ri = real_factory_3d(world).functor(xyz(idim));
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

void TDA::analyze(xfunctionsT& roots) const {

	std::cout << "\n\n!!!!!WARNING: Analyze not correct if nuclear correlation factors are used!!!!!!\n\n" << std::endl;

	const size_t noct=active_mo_.size();

	std::vector<xfunction>::iterator it;
	int iroot=0;
	for (it=roots.begin(); it!=roots.end(); ++it, ++iroot) {
		std::vector<double> norms=norm2s(world,it->x);

		// compute the oscillator strengths
		double osl=this->oscillator_strength_length(*it);
		double osv=this->oscillator_strength_velocity(*it);

		it->f_length = osl;
		it->f_velocity = osv;

		if(world.rank()==0) std::cout << std::scientific << std::setprecision(10) << std::setw(20);
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


	// Get overlap with chosen components of guess operators
}

vecfuncT TDA::project_to_ao_basis(const vecfuncT & active_mo, const vecfuncT &ao_basis)const{
	if(world.rank()==0) std::cout << "Projecting to AO basis ..." << std::endl;
	Tensor<double> Saa = matrix_inner(world,ao_basis,ao_basis,true);
	// Make AOMO overlap Matrix
	Tensor<double> Sam = matrix_inner(world,ao_basis,active_mo,false);
	// Get the coefficients : solve Sam = C*Sao
	Tensor<double> C;
	gesv(Saa, Sam, C);
	C = transpose(C);

	// make projected mos
	vecfuncT projected_mos;
	for(size_t i=0;i<active_mo.size();i++){
		//if(world.rank()==0)std::cout << "Projecting MO " << i << " to AO Basis ...\n Coefficients: ";
		real_function_3d tmp = real_factory_3d(world).compressed();
		for(size_t j=0;j<ao_basis.size();j++){
			tmp += C(i,j)*ao_basis[j];
			if(world.rank()==0 and fabs(C(i,j))>FunctionDefaults<3>::get_thresh()*100.0 and active_mo.size()<6)std::cout << C(i,j) <<  " ";
		}
		if(world.rank()==0)std::cout << std::endl;
		tmp.truncate();
		projected_mos.push_back(tmp);
		plot_plane(world,tmp,"projected_mo_"+stringify(i));
	}
	return projected_mos;
}

void TDA::memory_information(const xfunctionsT &xfunctions)const{
	double xmem = 0.0, potmem=0.0;
	for(size_t i=0;i<xfunctions.size();i++){
		xmem += get_size(world,xfunctions[i].x);
		potmem += get_size(world,xfunctions[i].smooth_potential);
	}
	if(world.rank()==0) std::cout << std::scientific << std::setprecision(2);
	if(world.rank()==0) std::cout << std::setw(40) << "xfunctions size..." << " : " << xmem << " GB" << std::endl;
	if(world.rank()==0) std::cout << std::setw(40) << "smoothed potentials..." << " : " << potmem << " GB" << std::endl;
}


