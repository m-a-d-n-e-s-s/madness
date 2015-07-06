/*
 * TDA.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: kottmanj
 */

#include "TDA.h"

//#include <chem/nemo.h>
//#include <madness/mra/mra.h>
#include "/usr/include/math.h"
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
	if(world.rank()==0) std::cout << "\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE_GUESS START " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n\n\n" << std::endl;

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
		if(world.rank()==0) std::cout << "\n\n\n" << "Guess Iteration Cycle " << i << "\n\n\n"<< std::endl;
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
	if(world.rank()==0)std::cout << std::setw(100) << "---End Initialize Guess Functions---" << " " << std::endl;
	init.info();

	// plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		plot_vecfunction(converged_xfunctions_[i].x,"final_guess_excitation_"+stringify(i),true);
	}

	// now sort the pre-converged xfunctions
	std::sort(converged_xfunctions_.begin(),converged_xfunctions_.end());
	if(converged_xfunctions_.size()>excitations_)converged_xfunctions_.erase(converged_xfunctions_.begin()+excitations_,converged_xfunctions_.end());


	if(world.rank()==0) std::cout << "\n\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE_GUESS ENDED " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n\n\n\n" << std::endl;

}

void TDA::solve(xfunctionsT &xfunctions) {
	if(world.rank()==0) std::cout << "\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE START " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n" << std::endl;

	if(world.rank()==0) std::cout << "\n------------------------------------------\n"
			<< "The following pre-converged guess_xfunctions will be used from now on: " << "\n------------------------------------------\n"  << std::endl;
	print_status(xfunctions);

	// Iterate till convergence is reached
	for(size_t iter=0;iter<300;iter++){
		if(world.rank()==0) std::cout << "\n\n\n" << "Solve Iteration Cycle " << iter << "\n\n\n"<< std::endl;
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


	// plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		plot_vecfunction(converged_xfunctions_[i].x,"final_par_excitation_"+stringify(i),true);
	}
	// print final result
	print_status(xfunctions);
	print_performance(xfunctions,"pre-");

	// Analyze
	//analyze(xfunctions);

	if(world.rank()==0) std::cout << "\n\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE ENDED " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n\n\n\n" << std::endl;

}

void TDA::solve_sequential(xfunctionsT &xfunctions) {

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
	real_function_3d rho = real_factory_3d(world);
	for(size_t i=0;i<active_mo_.size();i++){
		rho += active_mo_[i]*active_mo_[i];
	}
	plot_plane(world,rho,"rho");

	for(size_t p=0;p<converged_xfunctions_.size();p++){
		real_function_3d rhoprime = real_factory_3d(world);
		for(size_t i=0;i<active_mo_.size();i++){
			rhoprime += active_mo_[i]*converged_xfunctions_[p].x[i];
		}
		plot_plane(world,rhoprime,"rhoprime_"+stringify(p));
	}

	print("Final result :");
	print_status(converged_xfunctions_);
	analyze(converged_xfunctions_);
	print_performance(converged_xfunctions_,"final-");
	print_status(xfunctions);

}

void TDA::print_status(const xfunctionsT & xfunctions) const {
	if(world.rank()==0){
		std::cout << "\n" <<std::setw(5) << " #" << std::setw(20) << "omega" << std::setw(20) << "delta" << std::setw(20)
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
		if(use_nemo_){
			plot_vecfunction(mul(world,get_nemo().nuclear_correlation -> function(),xfunctions[i].x),"R_times_guess_function_"+stringify(i)+"_",true);
			plot_plane(world,get_nemo().nuclear_correlation -> function(),"nuclear_correlation_factor");
			plot_plane(world,get_calc().amo[0],"NEMO1");
			plot_plane(world,get_calc().amo[0]*get_nemo().nuclear_correlation -> function(),"RNEMO1");
			std::cout << "Overlap between reconstructed guess and reconstructed occupied space\n" << inner(world,mul(world,get_nemo().nuclear_correlation -> square(),xfunctions[i].x),get_calc().amo) << "\n";
		}

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
				Tensor<double> overlap_vec;
				if(use_nemo_){
					overlap_vec = inner(world, mul(world,get_nemo().nuclear_correlation ->square(),xfunctions[p].x),xfunctions[k].x);
				}else overlap_vec = inner(world, xfunctions[p].x,xfunctions[k].x);
				overlap(p, k) = overlap_vec.sum();
			}
		}

		Tensor<double> F = make_perturbed_fock_matrix(xfunctions);
		if(debug_) std::cout<< "guess overlap matrix is\n" << overlap << "\n";
		if(debug_) std::cout<<"guess pert. Fock Matrix is\n" << F << "\n";
		Tensor<double> U, evals, dummy(xfunctions.size());
		U = calc_ -> get_fock_transformation(world, overlap, F, evals, dummy,
				1.5 * econv_);
		if(debug_) std::cout<<"Transformation Matrix is\n" << U << "\n";
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
	if(world.rank()==0 and debug_){
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
	const std::size_t nmo=get_calc().amo.size();

	// first we need to determine the rotation of the external orbitals
	// to the MRA orbitals, so that we can subsequently rotate the guess
	Tensor<double> guess_phases_;
	vecfuncT koala_mo;

	// read koala's orbitals from disk
	for (std::size_t i=nfreeze_; i<nmo; ++i) {
		real_function_3d x_i=real_factory_3d(world).empty();
		const std::string valuefile="grid.koala.orbital"+stringify(i);
		x_i.get_impl()->read_grid2<3>(valuefile,functorT());
		koala_mo.push_back(x_i);
	}
	if(koala_mo.size()!=active_mo_.size()) MADNESS_EXCEPTION("ERROR in Koala guess: not the same number of Koala mos and active_mos of MRA",1);
	// this is the transformation matrix for the rotation
	guess_phases_=matrix_inner(world,koala_mo,get_calc().amo);
	guess_phases_=guess_phases_(_,Slice(nfreeze_,nmo-1));

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
	Sinv=inner(Sinv,U,-1,-1);

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
		for (std::size_t i=nfreeze_; i<nmo; ++i) {

			// this is the file where the guess is on disk
			const std::string valuefile="grid.koala.orbital"+stringify(i)
	    																			+".excitation"+stringify(iroot);
			real_function_3d x_i=real_factory_3d(world).empty();
			x_i.get_impl()->read_grid2<3>(valuefile,functorT());
			root.x.push_back(x_i);
		}

		// now rotate the active orbitals from the guess to conform with
		// the MRA orbitals
		// Sinv=Sinv(Slice(nfreeze_,nmo-1),Slice(nfreeze_,nmo-1));
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

	// Bool checks if the perturbed fock matrix has been calculated (if not the expencation value has to be calculated in the iterate_one routine)
	bool pert_fock_applied = false;

	// Restrict the number of parallel iterating guess functions
	for(size_t i=0;i<all_xfunctions.size();i++){
		all_xfunctions[i].current_residuals.clear();
		all_xfunctions[i].smooth_potential.clear();
		all_xfunctions[i].converged = false;
	}

	// make big fock diagonalization
	//make_big_fock_guess(all_xfunctions);
	orthonormalize_fock(all_xfunctions);

	std::sort(all_xfunctions.begin(),all_xfunctions.end());
	print_status(all_xfunctions);
	xfunctionsT xfunctions(all_xfunctions.begin(),all_xfunctions.begin()+iterating_excitations_);
	xfunctionsT remaining_xfunctions(all_xfunctions.begin()+iterating_excitations_,all_xfunctions.end());



	//	if(world.rank()==0) std::cout << "\nremaining guess functions are\n " << std::endl;
	//	print_status(remaining_xfunctions);

	size_t guess_iter_counter =1;
	std::vector<size_t> iteration_counters(xfunctions.size(),0);
	for(size_t i=0;i<1000;i++){
		TDA_TIMER iteration_time(world,"\nEnd of iteration " + stringify(i) +": ");

		{
			TDA_TIMER update_potentials(world,"Update potentials: ");
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
			for(size_t i=0;i<xfunctions.size();i++) project_out_occupied_space(xfunctions[i].x);
			project_out_converged_xfunctions(xfunctions);
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

	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("DFT is off but potential shift is not zero", 1);

	vecfuncT Vpsi;
	if(not use_nemo_){
	real_function_3d vn = get_calc().potentialmanager -> vnuclear();
	vn.truncate();
	Vpsi=add(world,xfunction.smooth_potential,mul(world,vn,xfunction.x));
	}else if(use_nemo_){
		Vpsi = xfunction.smooth_potential;
	}
	double omega = xfunction.omega;
	if(Vpsi.empty()) MADNESS_EXCEPTION("ERROR in iterate_one function: Applied potential of xfunction is empty",1);

	truncate(world, Vpsi); // no fence
	scale(world, Vpsi, -2.0);

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

	vecfuncT residual_bra = residual;
	if(use_nemo_) residual_bra = mul(world,get_nemo().nuclear_correlation -> square(),residual);
	Tensor<double> inner_tmp = inner(world, residual_bra, residual);
	if(debug_ and world.rank()==0) std::cout << "\n residual-self-overlaps \n" << inner_tmp << std::endl;
	double error = sqrt(inner_tmp.sum());
	xfunction.error.push_back(error);

	// Calculate 2nd order update:
	// Inner product of Vpsi and the residual (Vspi is scaled to -2.0 --> multiply later with 0.5)
	Tensor<double> tmp = inner(world, Vpsi, residual_bra);

	// Norm of GVpsi (Psi_tilde)
	double tmp2 = 0.0;
	for (size_t i = 0; i < GVpsi.size(); ++i) {
		double n = GVpsi[i].norm2();
		tmp2 += n * n;
	}

	// Factor 0.5 removes the factor 2 from the scaling before
	xfunction.delta.push_back(0.5 * tmp.sum() / tmp2);

	// return Updated x-function
	return GVpsi;

}


std::string TDA::update_energy(xfunction &xfunction)const {
	double thresh = FunctionDefaults<3>::get_thresh();
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
}

void TDA::normalize(xfunctionsT & xfunctions)const {
	for (size_t i = 0; i < xfunctions.size(); i++)
		normalize(xfunctions[i]);
}
void TDA::normalize(xfunction & xfunction)const {
	Tensor<double> self_overlaps;
	if(use_nemo_){
		vecfuncT RX = mul(world,get_nemo().nuclear_correlation -> function(),xfunction.x);
		self_overlaps = inner(world,RX,RX);
	}else self_overlaps = inner(world, xfunction.x, xfunction.x);
	const double squared_norm = self_overlaps.sum();
	const double norm = sqrt(squared_norm);
	scale(world, xfunction.x, 1.0 / norm);

	/// DEBUG
	if(use_nemo_){
	//std::cout << " Norm of regularized x-vector is " << sqrt(inner(world,xfunction.x,xfunction.x).sum()) <<"\n";
	//std::cout << " Norm of reconstructed x-vector is" << sqrt(inner(world,mul(world,get_nemo().nuclear_correlation -> square(),xfunction.x),xfunction.x).sum()) << "\n";
	}
}

void TDA::project_out_converged_xfunctions(xfunctionsT & xfunctions)const {
	for (size_t p = 0; p < xfunctions.size(); p++) {
		compress(world,xfunctions[p].x);
		for (size_t k = 0; k < converged_xfunctions_.size(); k++) {
			compress(world,converged_xfunctions_[k].x);
			Tensor<double> overlap;
			if(use_nemo_) overlap = inner(world, mul(world,get_nemo().nuclear_correlation -> square(),xfunctions[p].x),
					converged_xfunctions_[k].x);
			else overlap = inner(world, xfunctions[p].x,
					converged_xfunctions_[k].x);
			double c = overlap.sum();
			for (size_t i = 0; i < xfunctions[p].x.size(); i++) {
				xfunctions[p].x[i] -= c * converged_xfunctions_[k].x[i];
			}
		}
	}
}

void TDA::orthonormalize_GS(xfunctionsT &xfunctions)const {
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
bool TDA::orthonormalize_fock(xfunctionsT &xfunctions)const {
	normalize(xfunctions);

	if(xfunctions.size()<2){
		if(world.rank()==0) std::cout << "no Fock diagonalization just " << xfunctions.size() << " xfunctions" << std::endl;
		return false;
	}

	Tensor<double> overlap(xfunctions.size(), xfunctions.size());
	for (size_t p = 0; p < xfunctions.size(); p++) {
		for (size_t k = 0; k < xfunctions.size(); k++) {
			Tensor<double> overlap_vec;
			if(use_nemo_){
				overlap_vec = inner(world,mul(world,get_nemo().nuclear_correlation->square(),xfunctions[k].x),xfunctions[p].x);
			}
			else overlap_vec = inner(world, xfunctions[p].x,xfunctions[k].x);
			overlap(p, k) = overlap_vec.sum();
		}
	}

	if(world.rank()==0 and xfunctions.size()<6) std::cout << std::setprecision(3)<< "\n overlap matrix\n" << overlap << std::endl;
	// if the overlap matrix is already the unit matrix then no orthogonalization is needed
//	double overlap_offdiag = measure_offdiagonality(overlap, xfunctions.size());
//	if (fabs(overlap_offdiag) < FunctionDefaults<3>::get_thresh()) {
//		if(world.rank()==0) std::cout << " already orthogonal: perturbed fock matrix will not be calculated \n" <<std ::endl;
//		return false;
//	}

	Tensor<double> F(3L, xfunctions.size());
	if (dft_) {
		F = make_perturbed_fock_matrix(xfunctions);
	} else
		F = make_perturbed_fock_matrix(xfunctions);

	/// DEBUG
	std::cout << "\n\nPerturbed Fock Matrix F is \n" << F << "\n\n";

	// Diagonalize the perturbed Fock matrix
	Tensor<double> U, evals;

	madness::Tensor<double> dummy(xfunctions.size());
	U = calc_ -> get_fock_transformation(world, overlap, F, evals, dummy,
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
	if (x.size() != active_mo_.size())
		MADNESS_EXCEPTION(
				"ERROR IN PROJECTOR: Size of xfunctions and active orbitals is not equal",
				1);
	for (size_t p = 0; p < active_mo_.size(); ++p)
		if(use_nemo_){real_function_3d R2X = x[p]*get_nemo().nuclear_correlation -> square();
		x[p] -= rho0(R2X);
		}else{
			x[p] -= rho0(x[p]);
		}
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

double TDA::expectation_value(const xfunction &x, const vecfuncT &smooth_potential)const {
	Tensor<double> smooth_pot;
	if(use_nemo_) smooth_pot= inner(world, mul(world,get_nemo().nuclear_correlation -> square(),x.x), smooth_potential);
	else smooth_pot= inner(world,x.x, smooth_potential);
	double exp_smooth_v = smooth_pot.sum();

	double expv = exp_smooth_v;

	// The Vnuc part
	if(not use_nemo_){
	real_function_3d vnuc = get_calc().potentialmanager->vnuclear();
	vecfuncT vnuci = mul(world,vnuc,x.x);
	Tensor<double> vnuc_pot = inner(world, x.x, vnuci);
	double exp_vnuc = vnuc_pot.sum();
	expv += exp_vnuc;
	}

	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	for (int axis = 0; axis < 3; axis++) {
		const vecfuncT dx = apply(world, *(gradop[axis]), x.x);
		vecfuncT bra_dx;
		if(use_nemo_) bra_dx = apply(world, *(gradop[axis]), mul(world,get_nemo().nuclear_correlation -> square(),x.x));
		Tensor<double> kin;
		if(use_nemo_) kin= inner(world, bra_dx, dx);
		else kin= inner(world, dx, dx);
		expv += 0.5 * kin.sum();
	}
	for (size_t i = 0; i < x.x.size(); i++) {
		Tensor<double> overlap;
		if(use_nemo_) overlap = inner(world, x.x, x.x);
		else overlap = inner(world, x.x, mul(world,get_nemo().nuclear_correlation -> square(),x.x));
		expv -= (get_calc().aeps(nfreeze_ + i) + shift_) * overlap(i);
	}
	return expv;
}

Tensor<double> TDA::make_perturbed_fock_matrix(
		const xfunctionsT &xfunctions) const {
	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("No DFT calculation but shift is not zero", 1);
	if(world.rank()==0) std::cout << std::setw(40) << "perturbed fock matrix dimension..." << " : " << xfunctions.size() << "x" << xfunctions.size() << std::endl;

	//Tensor<double> F(xfunctions.size(),xfunctions.size());
	Tensor<double> F(xfunctions.size(), xfunctions.size());

/// NEMO: TO USE NEMO MAKE SHURE THE SMOOTH POTENTIAL IS THE TRANSFORMED NEMO POTENTIAL
	//The smooth potential part

	// The bra_element for nemos has to be multiplied with R^2 because of the changed metric

	for (std::size_t p = 0; p < xfunctions.size(); p++) {
		vecfuncT Vxp;
		if(xfunctions[p].smooth_potential.empty()) Vxp = apply_smooth_potential(xfunctions[p]);
		else Vxp = xfunctions[p].smooth_potential;
		for (std::size_t k = 0; k < xfunctions.size(); k++) {
			vecfuncT bra_x = xfunctions[k].x;
			if(use_nemo_){bra_x= mul(world,get_nemo().nuclear_correlation->square(),xfunctions[k].x);}
			Tensor<double> fpk_i = inner(world, bra_x, Vxp);
			F(p, k) = fpk_i.sum();
		}
	}

	if(debug_) std::cout << "Fock Matrix: smooth potential part:\n" << F << "\n";

/// NEMO: TURN THIS OF IF NEMO IS USED
	// The nuclear potential part
	if(not use_nemo_){
	for(size_t i=0; i < xfunctions.size();i++){
		real_function_3d vn = get_calc().potentialmanager->vnuclear();
		vn.truncate();
		vecfuncT vnuci = mul(world,vn,xfunctions[i].x);
		truncate(world,vnuci);
		for(size_t j=0; j < xfunctions.size();j++){
			Tensor<double> fij = inner(world, xfunctions[j].x, vnuci);
			F(i,j) += fij.sum();
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

		std::vector<vecfuncT> dxp, bra_dxp;
		for (std::size_t iroot = 0; iroot < xfunctions.size(); ++iroot) {

			const vecfuncT& xp = xfunctions[iroot].x;
			vecfuncT d = apply(world, *(gradop[axis]), xp);
			vecfuncT bra_d;
			if(use_nemo_){
				bra_d = apply(world, *(gradop[axis]), mul(world,get_nemo().nuclear_correlation -> square(),xp));
			}
			truncate(world,d);
			dxp.push_back(d);
			if(use_nemo_) bra_dxp.push_back(bra_d);
		}
		for (std::size_t iroot = 0; iroot < xfunctions.size(); ++iroot) {
			for (std::size_t jroot = 0; jroot < xfunctions.size(); ++jroot) {
				Tensor<double> xpi_Txqi;
				if(use_nemo_)xpi_Txqi = inner(world, bra_dxp[iroot], dxp[jroot]);
				else xpi_Txqi = inner(world, dxp[iroot], dxp[jroot]);
				F(iroot, jroot) += 0.5 * xpi_Txqi.sum();
			}
		}
	}

	// The epsilon part
	for (std::size_t p = 0; p < xfunctions.size(); p++) {
		for (std::size_t r = 0; r < xfunctions.size(); r++) {
			vecfuncT bra_x = xfunctions[p].x;
			if(use_nemo_){bra_x= mul(world,get_nemo().nuclear_correlation->square(),xfunctions[p].x);}
			Tensor<double> eij = inner(world, bra_x, xfunctions[r].x);
			for (size_t ii = 0; ii < xfunctions[p].x.size(); ++ii) {
				F(p, r) -= (get_calc().aeps[ii + nfreeze_] + shift_) * eij[ii];
			}
		}
	}

	return F;

}

Tensor<double> TDA::make_perturbed_fock_matrix_for_guess_functions(const xfunctionsT &xguess)const{
MADNESS_EXCEPTION("perturbed fock matrix for guess functions is not used",1);
//	Tensor<double> F(xguess.size(),xguess.size());
//
//	TDA_TIMER V0T(world,"unperturbed part: ");
//	// Unperturbed potential
//	for(size_t i=0;i<xguess.size();i++){
//		vecfuncT v0i = get_V0(xguess[i].x);
//		for(size_t j=0;j<xguess.size();j++){
//			Tensor<double> tmp = inner(world,xguess[j].x,v0i);
//			F(i,j) = tmp.sum();
//		}
//	}
//	V0T.info();
//
//	TDA_TIMER HT(world,"perturbed hartree part: ");
//	//	// perturbed hartree
//	for(size_t i=0;i<xguess.size();i++){
//		vecfuncT Ji = apply_hartree_potential(xguess[i].x);
//		for(size_t j=0;j<xguess.size();j++){
//			Tensor<double> tmp = inner(world,xguess[j].x,Ji);
//			F(i,j) += tmp.sum();
//		}
//	}
//	HT.info();
//
//	TDA_TIMER KT(world,"perturbed exchange part: ");
//	// perturbed exchange
//	vecfuncT intermediate= zero_functions<double,3>(world,active_mo_.size());
//	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr<
//			real_convolution_3d>(CoulombOperatorPtr(world, lo, bsh_eps_));
//	for(size_t i=0;i<intermediate.size();i++){
//		real_function_3d tmp =real_factory_3d(world);
//		for(size_t j=0;j<active_mo_.size();j++){
//			tmp += active_mo_[j] * (*poisson)( active_mo_[j]*active_mo_[i]);
//		}
//		intermediate[i] = tmp;
//	}
//	for(size_t i=0;i<xguess.size();i++){
//		std::shared_ptr<FunctionFunctorInterface<double, 3> >
//		polynom_functor(new polynomial_exop_functor(xguess[i].guess_excitation_operator));
//		for(size_t j=0;j<xguess.size();j++){
//			std::shared_ptr<FunctionFunctorInterface<double, 3> >
//			polynom_functor(new polynomial_exop_functor(xguess[j].guess_excitation_operator));
//			real_function_3d exop_j = real_factory_3d(world).functor(polynom_functor);
//			vecfuncT tmp = mul(world,exop_j,intermediate);
//			project_out_occupied_space(tmp);
//			Tensor<double> tensor_tmp = inner(world,xguess[i].x,tmp);
//			F(i,j) += tensor_tmp.sum();
//		}
//	}
//	KT.info();
//
//	//The kinetic part -1/2<xip|nabla^2|xir> = +1/2 <nabla xip||nabla xir>
//	for (std::size_t iroot = 0; iroot < xguess.size(); ++iroot) {
//		const vecfuncT& xp = xguess[iroot].x;
//		reconstruct(world, xp);
//	}
//
//	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
//	gradop = gradient_operator<double, 3>(world);
//
//	for (int axis = 0; axis < 3; ++axis) {
//
//		std::vector<vecfuncT> dxp;
//		for (std::size_t iroot = 0; iroot <xguess.size(); ++iroot) {
//
//			const vecfuncT& xp =xguess[iroot].x;
//			const vecfuncT d = apply(world, *(gradop[axis]), xp);
//			dxp.push_back(d);
//		}
//		for (std::size_t iroot = 0; iroot < xguess.size(); ++iroot) {
//			for (std::size_t jroot = 0; jroot < xguess.size(); ++jroot) {
//				Tensor<double> xpi_Txqi = inner(world, dxp[iroot], dxp[jroot]);
//				F(iroot, jroot) += 0.5 * xpi_Txqi.sum();
//			}
//		}
//	}
//
//	// The epsilon part
//	for (std::size_t p = 0; p < xguess.size(); p++) {
//		for (std::size_t r = 0; r < xguess.size(); r++) {
//			Tensor<double> eij = inner(world, xguess[p].x, xguess[r].x);
//			for (size_t ii = 0; ii < xguess[p].x.size(); ++ii) {
//				F(p, r) -= (get_calc().aeps[ii + nfreeze_] + shift_) * eij[ii];
//			}
//		}
//	}
//
//	return F;
}

vecfuncT TDA::apply_smooth_potential(const xfunction&xfunction) const{
/// NEMO: MOVE nemo_ BOOL TO THE .h FILE
	if(use_nemo_){
		real_function_3d vlocal = get_coulomb_potential();
		vecfuncT JX = mul(world,vlocal,xfunction.x);
		truncate(world,JX);
		Nuclear Uop(world,&get_nemo());
		vecfuncT UX = Uop(xfunction.x);
		truncate(world,UX);
		Exchange KOp(world,&get_nemo(),0);
		vecfuncT KX =KOp(xfunction.x);
		truncate(world,KX);
		vecfuncT TMP1 = sub(world,JX,KX);
		vecfuncT smooth_V0 = add(world,TMP1,UX);
		truncate(world,smooth_V0);
		vecfuncT gamma = apply_gamma(xfunction);
		truncate(world,gamma);
		vecfuncT smooth_V = add(world,smooth_V0,gamma);
		return smooth_V;
	}
	if(not use_nemo_){
	real_function_3d vlocal = get_coulomb_potential();
	vecfuncT J = mul(world,vlocal,xfunction.x);
	truncate(world,J);
//	vecfuncT K = get_calc().apply_hf_exchange(world, get_calc().aocc, mos_, xfunction.x);
	Exchange KOp(world,&*calc_,0);
	vecfuncT K=KOp(xfunction.x);
	truncate(world,K);
	vecfuncT smooth_V0 = sub(world,J,K);
	vecfuncT gamma = apply_gamma(xfunction);
	truncate(world,gamma);
	vecfuncT smooth_V = add(world,smooth_V0,gamma);
	return smooth_V;
	}
	MADNESS_EXCEPTION("ERROR in apply_smooth_potential... we should not reach the bottom",1);
}

vecfuncT TDA::apply_perturbed_potential(const xfunction & xfunction) const {

	// Take out the singularity in the nuclear potential with a nucleii correlation factor
	bool use_nuclear_correlation_factor_=true;
	if(use_nuclear_correlation_factor_){

	}
	else{
	// The smooth potential is the unperturbed potential without the nuclear potential plus the perturbed potential
	vecfuncT smooth_V  = apply_smooth_potential(xfunction);
	real_function_3d vnuc = get_calc().potentialmanager->vnuclear();
	vecfuncT vnucx = mul(world,vnuc,xfunction.x);
	vecfuncT Vpsi = add(world,vnucx,smooth_V);

	if(world.rank()==0) std::cout << std::scientific << std::setprecision(2) << "\n---Memory information for the potentials---" << std::endl;
	if(world.rank()==0) std::cout << std::setw(40) << "Applied Vnuc in GB: " << get_size(world,vnucx)  <<   std::endl;
	if(world.rank()==0) std::cout << std::setw(40) << "Smooth potential in GB: " << get_size(world,smooth_V) << std::endl;
	if(world.rank()==0) std::cout << std::setw(40) << "Overall potential in GB: " << get_size(world,Vpsi) << std::endl;
	if(world.rank()==0) std::cout << "\n\n"<<  std::setw(40) <<  std::endl;
	return Vpsi;
	}
}

vecfuncT TDA::apply_gamma(const xfunction &xfunction) const {

	TDA_TIMER hartree(world, "apply perturbed hartree potential...");
	vecfuncT gamma ;
	if(not triplet_) gamma = apply_hartree_potential(xfunction.x);
	else gamma = zero_functions<double,3>(world,xfunction.x.size());
	hartree.info(debug_);
	std::string the_message = "apply hf-exchange potential kernel, localization is " + stringify(localize_exchange_intermediate_) + " ...";
	TDA_TIMER exchange(world, the_message);\
	// transform xfunction to lmo basis
	if(localize_exchange_intermediate_){
		// not used and not consistent with nemo
		MADNESS_EXCEPTION("Localization was demanded ... not consistent with nemo ... and does not work anyway right now",1);
//		std::vector<int> set=calc_.group_orbital_sets(world,calc_.aeps,calc_.aocc,active_mo_.size());
//		distmatT dmo2lmo=calc_.localize_PM(world,active_mo_,set);
//		tensorT mo2lmo(active_mo_.size(),active_mo_.size());
//		dmo2lmo.copy_to_replicated(mo2lmo);
//		vecfuncT lmo_x=madness::transform(world,xfunction.x,mo2lmo,true);
//		vecfuncT gamma_ex_lmo=zero_functions<double,3>(world,lmo_x.size());
//		for (std::size_t p = 0; p < xfunction.x.size(); p++) {
//
//			vecfuncT x_Ppi=zero_functions<double,3>(world,lmo_x.size());
//			norm_tree(world,lmo_x);
//			norm_tree(world,exchange_intermediate_[p]);
//			for(size_t i=0 ; i<x_Ppi.size();i++){
//				x_Ppi[i]=mul_sparse(lmo_x[i],exchange_intermediate_[p][i],FunctionDefaults<3>::get_thresh());
//			}
//			compress(world, x_Ppi);
//			for (std::size_t i = 0; i < lmo_x.size(); i++) {
//				gamma_ex_lmo[p] -= x_Ppi[i];
//			}
//
//			// Project out occupied space
//			gamma_ex_lmo[p] -= rho0(gamma[p]);
//		}
//		// Transform gamma to canonical basis
//		vecfuncT gamma_ex_canon=madness::transform(world,gamma_ex_lmo,transpose(mo2lmo),true);
//		gamma = add(world,gamma,gamma_ex_canon);
	}
	else{
		for (std::size_t p = 0; p < xfunction.x.size(); p++) {

			const vecfuncT x_Ppi = mul(world, xfunction.x, exchange_intermediate_[p]);
			compress(world, x_Ppi);
			for (std::size_t i = 0; i < xfunction.x.size(); i++) {
				gamma[p] -= x_Ppi[i];
			}

			// Project out occupied space
			if(use_nemo_){
				real_function_3d R2gamma = get_nemo().nuclear_correlation -> square()*gamma[p];
				gamma[p] -= rho0(R2gamma);
			}else gamma[p] -= rho0(gamma[p]);
		}
	}

	exchange.info(debug_);

	return gamma;
}

vecfuncT TDA::apply_gamma_dft(const xfunction &xfunction) const {
MADNESS_EXCEPTION("NO TDDFT AVAILABLE RIGHT NOW",1);
//	TDA_TIMER hartree(world, "apply perturbed hartree potential...");
//	vecfuncT gamma = apply_hartree_potential(xfunction.x);
//	hartree.info(debug_);
//
//	TDA_TIMER rhoprime(world, "make perturbed density...");
//
//	// Make the perturbed density for closed shell molecules
//	real_function_3d perturbed_density = real_factory_3d(world);
//	for (size_t i = 0; i < active_mo_.size(); i++) {
//		perturbed_density += 2.0 * xfunction.x[i] * active_mo_[i];
//	}
//	rhoprime.info(debug_);
//	//
//	//	TDA_TIMER applyit(world,"apply vxc...");
//	//
//	// Get the perturbed xc potential from the dft class
//	//		real_function_3d vxc = xclib_interface_.convolution_with_kernel(
//	//				perturbed_density);
//
//	//		for (size_t i = 0; i < gamma.size(); i++) {
//	//			gamma[i] += vxc * active_mo_[i];
//	//		}
//
//	// Alternative way (more expensive, but avoid the unprecise kernel)
//	// for small test molecules this seems to bring no improvement
//	// when using this:
//	// 1.comment out the line before
//	// 2.return add(world,gamma,gamma2)
//	// 3. dont forget to project out occupied space also from gamma2 (below here)
//	// THIS DOES NOT WORK FOR GGA
//	//	vecfuncT gamma2=xclib_interface_.apply_kernel(xfunction.x);
//	//	for (int p=0; p<active_mo_.size(); ++p) gamma2[p] -= rho0(gamma2[p]);
//
//	//	vecfuncT gamma2 = mul(world,perturbed_density,lda_intermediate_);
//
//	for(size_t p=0;p<gamma.size();p++){
//		gamma[p] += xclib_interface_.get_fxc()*perturbed_density*active_mo_[p];
//	}
//	plot_vecfunction(gamma,"complete_gamma_");
//
//	// project out occupied space
//	for (size_t p = 0; p < active_mo_.size(); ++p)
//		gamma[p] -= rho0(gamma[p]);
//
//	//applyit.info(debug_);
//	//plot_vecfunction(gamma, "gamma", plot_);
//	//return add(world,gamma,gamma2);
//	truncate(world, gamma);
//	return gamma;
}

vecfuncT TDA::apply_hartree_potential(const vecfuncT &x) const {
	// Make the perturbed density (Factor 2 is for closed shell)
	real_function_3d perturbed_density = real_factory_3d(world);
	for (size_t i = 0; i < x.size(); i++) {
		perturbed_density += 2.0 * x[i] * active_mo_[i];
	}
	// if nemo is used we need the R2 metric so we just have to multiply the perturbed density with R2 before the coulomb operator is applied
	if(use_nemo_){
		perturbed_density= perturbed_density*get_nemo().nuclear_correlation->square();
	}
	real_convolution_3d J = CoulombOperator(world, lo, bsh_eps_);
	real_function_3d Jrhoprime = J(perturbed_density);

	return mul(world, Jrhoprime, active_mo_);

}

vecfuncT TDA::get_V0(const vecfuncT& x) const {
	// This function is not consistent with nemo right now
	MADNESS_EXCEPTION("get_V0 was called ... this should not happen anymore",1);
	// the local potential V^0 of Eq. (4)
	real_function_3d vlocal = get_calc().potentialmanager->vnuclear()
	        + get_coulomb_potential();

	// make the potential for V0*xp
	vecfuncT Vx = mul(world, vlocal, x);

	// and the exchange potential is K xp
	vecfuncT Kx;
	if (not dft_) {
//		Kx = get_calc().apply_hf_exchange(world, get_calc().aocc, mos_, x);
	    Exchange K(world,&*calc_,0);
	    Kx=K(x);
	}
	if (dft_) {
		MADNESS_EXCEPTION("NO TDDFT AVAILABLE RIGHT NOW",1);
		//real_function_3d vxc = xclib_interface_.get_unperturbed_vxc();
		// Shift the unperturbed potential down
		//vxc = vxc + shift_;
		//Kx = mul(world, vxc, x);
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
	if(use_nemo_) rho = rho*get_nemo().nuclear_correlation->square();
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
		// Integrant
		vecfuncT Integrant = mul(world, active_mo[p], amo);
		if(use_nemo_) Integrant = mul(world,get_nemo().nuclear_correlation->square(),Integrant);
		truncate(world,Integrant);
		intermediate[p] = apply(world, (*poisson),Integrant);
	}
	double immem=0.0;
	for(size_t i=0;i<intermediate.size();i++){
		immem += get_size(world,intermediate[i]);
	}
	if(world.rank()==0) std::cout << "\n\n----Memory information in GB for the exchange intermediate---- "<<  std::endl;
	if(world.rank()==0) std::cout << std::scientific << std::setprecision(2) << immem <<"\n\n"<< std::endl;
	return intermediate;
}

std::vector<vecfuncT> TDA::make_localized_exchange_intermediate() const {
	vecfuncT lmo=madness::transform(world,active_mo_,mo2lmo_,true);

	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr<
			real_convolution_3d>(CoulombOperatorPtr(world, lo, bsh_eps_));

	std::vector<vecfuncT> intermediate(lmo.size());

	for (std::size_t p = 0; p < lmo.size(); ++p) {
		intermediate[p] = apply(world, (*poisson),
				mul(world, lmo[p], lmo));
	}
	if(world.rank()==0) std::cout << std::setw(40) << "\n\n----Memory information in GB for the exchange intermediate---- "<<  std::endl;
	if(world.rank()==0) std::cout << std::scientific << std::setprecision(2) << get_size(world,intermediate[0]) << " x " << intermediate.size() << std::setw(40) <<  std::endl;
	return intermediate;
}

vecfuncT TDA::make_lda_intermediate()const{
	MADNESS_EXCEPTION("NO TDDFT AVAILABLE RIGHT NOW",1);
//	vecfuncT mo;
//	for(size_t i=0;i<active_mo_.size();i++) mo.push_back(copy(active_mo_[i]));
//	vecfuncT result = xclib_interface_.get_lda_intermediate(mo);
//	plot_vecfunction(result,"lda_intermediate_");
//	mo.clear();
//	return result;
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
	results << "% bsh_eps: " << bsh_eps_ << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << " \t \t \n";
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
	results2 << "% bsh_eps: " << bsh_eps_ << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << "\n";
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
	results2 << "\\multicolumn{2}{l}{ bsh_eps: " << bsh_eps_ << " thresh: " << FunctionDefaults<3>::get_thresh() << " time: " << wall_time() << "} \\\\ \n" ;
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

void TDA::analyze(xfunctionsT& roots) const {

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
