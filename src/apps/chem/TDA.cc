/*
 * TDA.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: kottmanj
 */

#include "TDA.h"

#include "/usr/include/math.h"
#include "../../madness/mra/derivative.h"
#include "../../madness/mra/funcdefaults.h"
#include "../../madness/mra/funcplot.h"
#include "../../madness/mra/function_interface.h"
#include "../../madness/mra/functypedefs.h"
#include "../../madness/mra/vmra1.h"
#include "../../madness/tensor/distributed_matrix.h"
#include "../../madness/tensor/srconf.h"
#include "../../madness/tensor/tensor.h"
#include "../../madness/world/archive.h"
#include "../../madness/world/parar.h"
#include "../../madness/world/print.h"
#include "../../madness/world/shared_ptr_bits.h"
#include "../../madness/world/worldexc.h"
#include "../../madness/world/worldfwd.h"
#include "../../madness/world/worldgop.h"
#include "../../madness/world/worldtime.h"
#include "molecule.h"
#include "nemo.h"
#include "potentialmanager.h"

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
	plot_vecfunction(active_mo_, "active_mo_");
	// check for saved xfunctions
	read_xfunctions(xfunctions);
	// read keyword : only read xfunctions and analyze them (this happens in solve_sequential function)
	if (read_ or only_sequential_){
		std::cout << "\n\n ----- found read keyword ... skipping iterations \n\n" << std::endl;
		return;
	}

	// Create the excitation_function vector and initialize
	std::cout << "\n\n---Start Initialize Guess Functions---" << "\n\n " << std::endl;
	TDA_TIMER init(world,"\nfinished to initialize guess excitations ...");
	// if xfunctions were read in before then xfunctions.empty() will be false
	if(xfunctions.empty())initialize(xfunctions);
	iterate_guess(xfunctions);
	std::cout << std::setw(100) << "---End Initialize Guess Functions---" << " " << std::endl;
	init.info();

	// plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		plot_vecfunction(converged_xfunctions_[i].x,"final_guess_excitation_"+stringify(i),true);
	}

	// now sort the pre-converged xfunctions
	std::sort(converged_xfunctions_.begin(),converged_xfunctions_.end());

	if(world.rank()==0) std::cout << "\n\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE_GUESS ENDED " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n\n\n\n" << std::endl;

}

void TDA::solve(xfunctionsT &xfunctions) {
	// read keyword : only read xfunctions and analyze them (this happens in solve_sequential function)
	if (read_ or only_sequential_){
		std::cout << "\n\n ----- found read keyword ... skipping iterations \n\n" << std::endl;
		return;
	}
	if(world.rank()==0) std::cout << "\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE START " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n" << std::endl;

	// use only the demanded number of xfunctions
	std::sort(xfunctions.begin(),xfunctions.end());
	if(xfunctions.size()>guess_excitations_) xfunctions.erase(xfunctions.begin()+guess_excitations_,xfunctions.end());

	if(world.rank()==0) std::cout << "\n------------------------------------------\n"
			<< "The following pre-converged guess_xfunctions will be used from now on: " << "\n------------------------------------------\n"  << std::endl;
	print_status(xfunctions);

	// Iterate till convergence is reached
	iterate(xfunctions);

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
	if(world.rank()==0) std::cout << "\n\n\n\n------------------------------------------------------------------------------------------------------------------------\n"
			<< "SOLVE_SEQUENTIAL START " << "\n------------------------------------------------------------------------------------------------------------------------\n\n\n"
			"The first " << excitations_ << " of the following xfunctions will be solved sequentially "<< std::endl;
	std::sort(xfunctions.begin(),xfunctions.end());
	xfunctions.erase(xfunctions.begin()+excitations_,xfunctions.end());
	print_status(xfunctions);

	// on the fly or not makes no sense here, but since the same input file is used for both solve functions this has to be here
	if (not on_the_fly_)
		on_the_fly_ = true;

	if(not read_){

		print("\n\n\n\n-------------------------------------------------------");
		print("BEGINNING THE FINAL ITERATIONS TO AN ACCURACY OF ", hard_dconv_);
		print("-------------------------------------------------------\n\n\n\n");

		if(xfunctions.size()!=excitations_) {
			print("Wrong size in xfunctions!!!!");
			//MADNESS_EXCEPTION("Wrong size in xfunction vector",1);
		}
		// The fock matrix routine needs a dummy kain solver
		kain_solver_helper_struct dummy_kain(world,1,active_mo_.size(),xfunctions.size(),false);
		orthonormalize_fock(xfunctions,true,dummy_kain);

		size_t max = xfunctions.size();

		// failsafe
		if(excitations_ > xfunctions.size()){
			print("\nDemanded ", excitations_, " excitations, but only ",
					xfunctions.size(), " pre-converged\n");
			max = xfunctions.size();
		}

		// The given xfunctions should be sorted by energy in ascending order:
		for (size_t iroot = 0; iroot < max; iroot++) {
			print("\n\n-----xfunction ", iroot, " ------\n\n");
			// create the kain solver for the current xfunction
			sequential_kain_solver solver(TDA_allocator(world, xfunctions[iroot].x.size()));
			solver.set_maxsub(kain_subspace_);
			kain_ = true;

			// begin the final iterations
			for (int iter = 0; iter < 100; iter++) {
				TDA_TIMER iteration_timer(world, "");
				normalize(xfunctions[iroot]);
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
				if (xfunctions[iroot].error.back() < hard_dconv_ and fabs(xfunctions[iroot].delta.back()) < hard_econv_) {
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
	}else{
		for(size_t i=0;i<xfunctions.size();i++) converged_xfunctions_.push_back(xfunctions[i]);
		std::cout << "\n\n ----- found read keyword ... skipping iterations \n\n" << std::endl;
	}

	plot_vecfunction(active_mo_, "active_mo_");
	// plot
	for(size_t i=0;i<converged_xfunctions_.size();i++) {
		plot_vecfunction(converged_xfunctions_[i].x,"final_excitation_"+stringify(i),true);
	}
	// save converged xfunctions to disc
	if(not read_) save_xfunctions(xfunctions);

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
	print_status(xfunctions);
	analyze(converged_xfunctions_);
	print_performance(converged_xfunctions_,"final-");
	print_status(xfunctions);

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

void TDA::initialize(xfunctionsT & xfunctions)const{
	plot_vecfunction(active_mo_,"active_mo_",true);
	if (guess_ == "physical") {
		guess_physical(xfunctions);
		if (not guess_omegas_.empty() and guess_omegas_.size()<=xfunctions.size()){
			for(size_t i=0;i<guess_omegas_.size();i++){
				xfunctions[i].omega = guess_omegas_[i];
			}
		}
	}
	else if(guess_ == "custom") {
		guess_custom(xfunctions);
	}
	else if(guess_ == "atomic_excitation"){
		guess_atomic_excitation(xfunctions);
	}
	else if(guess_ == "custom_2"){
		guess_custom_2(xfunctions);
	}
    else if(guess_ == "local"){
        guess_local(xfunctions);
    }
	else {
		if (world.rank() == 0)
			print("unknown keyword for guess: ", guess_);
		MADNESS_EXCEPTION("Reached end of initialize", 1);
	}
	// failsafe
	if(guess_xfunctions_.size()<guess_excitations_ and not replace_guess_functions_) MADNESS_EXCEPTION("Too many guess_excitation demanded ... replace_guess_functions_ is false",1);

}

void TDA::guess_custom(xfunctionsT & xfunctions)const{

	std::shared_ptr<FunctionFunctorInterface<double, 3> > smooth_functor(
			new guess_smoothing(guess_box_));

	real_function_3d smoothing_function = real_factory_3d(world).functor(smooth_functor);
	plot_plane(world,smoothing_function,"smoothing_function");

	exoperators exops(world,excitation_point_);

	//	// project the active_mos on ao functions
	//	vecfuncT projected_mos = zero_functions<double,3>(world,active_mo_.size());
	//	for(size_t i=0;i<active_mo_.size();i++){
	//		projected_mos[i] = exops.project_function_on_aos(world,active_mo_[i]);
	//	}
	//	plot_vecfunction(projected_mos,"projected_mos_");

	std::vector<vecfuncT> xguess_inner = exops.make_custom_guess(world,custom_exops_,active_mo_, smoothing_function);
	//std::vector<vecfuncT> xguess_outer = exops.make_custom_guess(world,custom_exops_,projected_mos, anti_smoothing_function);
	for(size_t i=0;i<xguess_inner.size();i++){
		xfunction tmp(world);
		tmp.omega = guess_omegas_[i];
		tmp.x = xguess_inner[i];
		plot_vecfunction(tmp.x,"custom_guessfunction_"+stringify(i));
		//plot_vecfunction(xguess_inner[i],"custom_guessfunction_inner"+stringify(i));
		//plot_vecfunction(xguess_outer[i],"custom_guessfunction_outer"+stringify(i));
		xfunctions.push_back(tmp);
	}
	//	//TESTING DEBUG
	//	double c=inner(xfunctions[0].x[0],projected_mos[0]);
	//	xfunctions[0].x[0] -= c*projected_mos[0];
	//	plot_vecfunction(xfunctions[0].x,"custom_guessfunction_woocc");



}

void TDA::guess_physical(xfunctionsT & xfunctions)const{
	std::shared_ptr<FunctionFunctorInterface<double, 3> > smooth_functor(
			new guess_smoothing(guess_box_));

	real_function_3d smoothing_function = real_factory_3d(world).functor(smooth_functor);
	smoothing_function.print_size("smoothing_function size");
	plot_plane(world,smoothing_function,"smoothing_function");

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
				tmp = all_orbitals * xoperators[j]*smoothing_function;
			} else
				tmp = active_mo_[i] * xoperators[j]*smoothing_function;

			double norm = tmp.norm2();
			tmp.scale(1.0 / norm);
			x.push_back(copy(tmp));
		}
		xfunction xtmp(world, x);
		xtmp.number = xfunctions.size();
		xtmp.omega = guess_omega_;
		xtmp.kain = false;
		normalize(xtmp);
		project_out_occupied_space(xtmp.x);
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

}

void TDA::guess_custom_2(xfunctionsT & xfunctions)const{
	if(world.rank()==0) std::cout << "Making custom_2 guess: " << custom_exops_.size() << " custom excitations given" << std::endl;
	guess guess_structure(world,calc_.param.L,guess_mode_,active_mo_,calc_.ao);
	for(size_t i=0;i<custom_exops_.size();i++){
		xfunction tmp(world,guess_omega_);
		vecfuncT xtmp = guess_structure.make_custom_guess(world,custom_exops_[i]);
		project_out_occupied_space(xtmp);
		tmp.x=xtmp;
		tmp.number=i;
		xfunctions.push_back(tmp);
		plot_vecfunction(tmp.x,"guess_excitation_" + stringify(tmp.number) + "_");
	}
	normalize(xfunctions);
	plot_vecfunction(active_mo_,"active_mo_",true);
}

void TDA::guess_local(xfunctionsT & xfunctions)const{

    if(world.rank()==0) print("Making local guess");

    // localize the active orbitals

    std::vector<int> set=calc_.group_orbital_sets(world,calc_.aeps,calc_.aocc,active_mo_.size());
    distmatT dmo2lmo=calc_.localize_PM(world,active_mo_,set);
    tensorT mo2lmo(active_mo_.size(),active_mo_.size());
    dmo2lmo.copy_to_replicated(mo2lmo);
    vecfuncT lmo=madness::transform(world,active_mo_,mo2lmo,true);

    guess guess_structure(world,calc_.param.L,guess_mode_,lmo,calc_.ao);

    for(size_t i=0;i<custom_exops_.size();i++){
        xfunction tmp(world,guess_omega_);
        vecfuncT xtmp = guess_structure.make_custom_guess(world,custom_exops_[i]);
        project_out_occupied_space(xtmp);

        tmp.x=transform(world,xtmp,transpose(mo2lmo),true);
        tmp.number=i;
        xfunctions.push_back(tmp);
        plot_vecfunction(tmp.x,"guess_excitation_" + stringify(tmp.number) + "_");
    }
    normalize(xfunctions);
    plot_vecfunction(active_mo_,"active_mo_",true);
}

void TDA::guess_atomic_excitation(xfunctionsT & xfunctions)const{
	guess guess_structure(world,calc_.param.L,guess_mode_,active_mo_,calc_.ao);
	for(size_t i=0;i<custom_exops_.size();i++){
		xfunction tmp(world,guess_omega_);
		vecfuncT xtmp = guess_structure.make_atomic_guess(world,custom_exops_[i],active_mo_.size());
		project_out_occupied_space(xtmp);
		tmp.x=xtmp;
		tmp.number=i;
		plot_vecfunction(tmp.x,"guess_excitation_"+stringify(i)+"_",true);
		normalize(tmp);
		xfunctions.push_back(tmp);
	}
	normalize(xfunctions);
}

void TDA::guess_koala(World &world, xfunctionsT &roots)const{

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
	if(world.rank()==0) std::cout << "Overlap between transformed Koala MOs and MRA MOs \n" << check << std::endl;
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


vecfuncT TDA::make_excitation_operators()const {
	if(guess_exop_ == "custom"){

		exoperators exops(world,excitation_point_);
		vecfuncT xoperators = exops.get_custom_exops(world,custom_exops_);
		return xoperators;

	}else{
		std::string key;
		if (guess_exop_ == "symmetry") {
			key = calc_.molecule.pointgroup_;
			std::cout << std::setw(40) << "symmetry_guess" << " : point group is " << key;
			if (key == "C1") {
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

		exoperators exops(world,excitation_point_);
		vecfuncT xoperators = exops.get_exops(world, key);

		return xoperators;
	}
}

void TDA::iterate_guess(xfunctionsT &xfunctions) {
	std::cout << "\n" << std::endl;
	std::cout << "---Start Guess Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,true);
	// set back iteration counter
	//for(size_t i=0;i<xfunctions.size();i++)xfunctions[i].iterations = 0;
	std::cout <<std::setw(100)<< "---End Guess Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate(xfunctionsT &xfunctions) {
	std::cout << "\n" << std::endl;
	std::cout << "---Start Main Iterations---" << "\n\n " << std::endl;
	iterate_all(xfunctions,false);
	std::cout <<std::setw(100)<< "---End Main Iterations---" << "\n\n " << std::endl;
}

void TDA::iterate_all(xfunctionsT &xfunctions, bool guess) {
	size_t max = iter_max_ * excitations_ + guess_iter_;

	// Initialize the kain solver (not during guess iteraions)
	bool kain = kain_;
	if(guess) kain = false;
	kain_solver_helper_struct kain_solver(world,kain_subspace_,active_mo_.size(),xfunctions.size(),kain_);

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
		std::cout<< std::setw(40)<<"Starting iteration"+stringify(iter)+" at time   "
				<<" : "<<std::setprecision(6)<<wall_time() << std::endl;
		//  memory management
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
				xfunctions[i].Vx = apply_perturbed_potential(xfunctions[i]);
				//xfunctions[i].Vx.back().print_size("size of homo pert pot "+stringify(i));
			}
			update_potentials.info();

			// Truncate xfunctions
			TDA_TIMER truncate2(world, "truncate |x> and |Vx>...");
			truncate_xfunctions(xfunctions);
			truncate2.info();
		}

		if (ptfock_orthonormalization) {
			TDA_TIMER orthonormalize_fock_matrix(world, "orthogonalize_fock : ");

			ptfock_orthonormalization = orthonormalize_fock(xfunctions, guess, kain_solver);
			orthonormalize_fock_matrix.info();

			// Update the energies
			update_energies(xfunctions);

			// Truncate xfunctions
			TDA_TIMER truncate3(world, "truncate |x> and |Vx>...");
			truncate_xfunctions(xfunctions);
			truncate3.info();
		}
		// if there is no fock update demanded: just update the energy
		update_energies(xfunctions);

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
			kain_solver.update(xfunctions);
			kain.info();
		}
		//check convergence
		check_convergence(xfunctions,guess);
		size_t counter = 0;
		std::cout << std::setw(40) << "converged excitations ..." << " : ";
		for(size_t i=0;i<xfunctions.size();i++){
			std::cout << " " << i << " (" << xfunctions[i].converged << ") ";
			if(xfunctions[i].converged){
				if(guess){
					converged_xfunctions_.push_back(xfunctions[i]);
					if(guess_xfunctions_.empty()) initialize(guess_xfunctions_);
					xfunctions[i]=guess_xfunctions_[i];
					i=-1; // set back i for the case that more than one xfunction converged
				}else if(kain_){
					counter++;
				}else{
					converged_xfunctions_.push_back(xfunctions[i]);
					xfunctions.erase(xfunctions.begin()+i);
				}
			}
		}std::cout << std::endl;
		if(not kain_) project_out_converged_xfunctions(xfunctions);

		// convergence criterium
		if(guess or not kain_){
			if (converged_xfunctions_.size() >= guess_excitations_){
				// push back all the xfunctions to converged (for the case that lower energy solutions are there)
				for(size_t i=0;i<xfunctions.size();i++){
					if(xfunctions[i].iterations > 1) converged_xfunctions_.push_back(xfunctions[i]);
				}
				break;
			}
		}
		if(kain_){
			if (counter >= excitations_){
				converged_xfunctions_.clear();
				for(size_t i=0;i<xfunctions.size();i++){
					if(xfunctions[i].converged) converged_xfunctions_.push_back(xfunctions[i]);
				}
				break;
			}
		}

		Iterationtime.info();
		print_status(xfunctions);
		for (size_t p = 0; p < xfunctions.size(); p++) {
			plot_vecfunction(xfunctions[p].x,
					"current_iteration_" + stringify(p) + "_", plot_);
		}
	}
}

void TDA::iterate_one(xfunction & xfunction, bool ptfock, bool guess)const {
	xfunction.iterations += 1;

	if (not dft_ and shift_ != 0.0)
		MADNESS_EXCEPTION("DFT is off but potential shift is not zero", 1);

	vecfuncT Vpsi;
	if (on_the_fly_)
		Vpsi = apply_perturbed_potential(xfunction);
	else if (not on_the_fly_){
		if(xfunction.Vx.empty()) Vpsi = apply_perturbed_potential(xfunction);
		else Vpsi = xfunction.Vx;
	}

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
	if(debug_ and world.rank()==0) std::cout << "\n residual-self-overlaps \n" << inner_tmp << std::endl;
	double error = sqrt(inner_tmp.sum());
	xfunction.error.push_back(error);
	if (error < kain_conv_thresh_)
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
	scale(world, xfunction.x, 1.0 / sqrt(tmp2));

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


void TDA::update_energies(xfunctionsT &xfunctions)const {
	double thresh = FunctionDefaults<3>::get_thresh();
	std::cout << std::setw(40) << "update energies..." << " : ";
	for(size_t k=0;k<xfunctions.size();k++) {

		//failsafe: make shure the delta and expectation values vectors are not empty to avoid segmentation faults
		if(not xfunctions[k].delta.empty() and not xfunctions[k].expectation_value.empty() and not xfunctions[k].error.empty()) {
			if(xfunctions[k].expectation_value.back() < highest_excitation_) {
				if(fabs(xfunctions[k].delta.back()) < thresh*2.0) {
					xfunctions[k].omega +=xfunctions[k].delta.back();
					std::cout << k << "(2nd), ";
				} else {
					xfunctions[k].omega = xfunctions[k].expectation_value.back();
					std::cout << k << "(exp), ";
				}
			} else {
				// get the highest converged energy, of no xfunction converged already use the guess_omega_ energy
				double setback_energy = guess_omega_;
				// set the last converged value of the same type of guess as the default
				//double new_omega = highest_excitation_*0.8;// default
				xfunctions[k].omega = setback_energy;
				std::cout << k << "(-) , ";
			}

		}

	}std::cout<<std::endl;

}

void TDA::normalize(xfunctionsT & xfunctions)const {
	for (size_t i = 0; i < xfunctions.size(); i++)
		normalize(xfunctions[i]);
}
void TDA::normalize(xfunction & xfunction)const {
	const Tensor<double> self_overlaps = inner(world, xfunction.x, xfunction.x);
	const double squared_norm = self_overlaps.sum();
	const double norm = sqrt(squared_norm);
	//std::cout << "Norm if xfunction " << xfunction.number << " is " << norm << " squared norm is " << squared_norm << std::endl;
	scale(world, xfunction.x, 1.0 / norm);
}

void TDA::project_out_converged_xfunctions(xfunctionsT & xfunctions)const {
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
	project.info();
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
bool TDA::orthonormalize_fock(xfunctionsT &xfunctions, bool guess, kain_solver_helper_struct & kain_solver)const {
	normalize(xfunctions);
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

	// rotate Kain subspace
	if (kain_ and not guess) {
		kain_solver.transform_subspace(world, U);
		std::cout << "Kain subspace is transformed...";
	}
	std::cout << std::endl;

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

double TDA::expectation_value(const xfunction &x, const vecfuncT &Vx)const {
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
	vecfuncT gamma ;
	if(not triplet_) gamma = apply_hartree_potential(xfunction.x);
	else gamma = zero_functions<double,3>(world,xfunction.x.size());
	hartree.info(debug_);
	std::string the_message = "apply hf-exchange potential kernel, localization is " + stringify(localize_exchange_intermediate_) + " ...";
	TDA_TIMER exchange(world, the_message);\
	// transform xfunction to lmo basis
	if(localize_exchange_intermediate_){
		std::vector<int> set=calc_.group_orbital_sets(world,calc_.aeps,calc_.aocc,active_mo_.size());
		distmatT dmo2lmo=calc_.localize_PM(world,active_mo_,set);
		tensorT mo2lmo(active_mo_.size(),active_mo_.size());
		dmo2lmo.copy_to_replicated(mo2lmo);
		vecfuncT lmo_x=madness::transform(world,xfunction.x,mo2lmo,true);
		vecfuncT gamma_ex_lmo=zero_functions<double,3>(world,lmo_x.size());
		for (std::size_t p = 0; p < xfunction.x.size(); p++) {

			vecfuncT x_Ppi=zero_functions<double,3>(world,lmo_x.size());
			norm_tree(world,lmo_x);
			norm_tree(world,exchange_intermediate_[p]);
			for(size_t i=0 ; i<x_Ppi.size();i++){
				x_Ppi[i]=mul_sparse(lmo_x[i],exchange_intermediate_[p][i],FunctionDefaults<3>::get_thresh());
			}
			compress(world, x_Ppi);
			for (std::size_t i = 0; i < lmo_x.size(); i++) {
				gamma_ex_lmo[p] -= x_Ppi[i];
			}

			// Project out occupied space
			gamma_ex_lmo[p] -= rho0(gamma[p]);
		}
		// Transform gamma to canonical basis
		vecfuncT gamma_ex_canon=madness::transform(world,gamma_ex_lmo,transpose(mo2lmo),true);
		gamma = add(world,gamma,gamma_ex_canon);
	}
	else{
		for (std::size_t p = 0; p < xfunction.x.size(); p++) {

			const vecfuncT x_Ppi = mul(world, xfunction.x, exchange_intermediate_[p]);
			compress(world, x_Ppi);
			for (std::size_t i = 0; i < xfunction.x.size(); i++) {
				gamma[p] -= x_Ppi[i];
			}

			// Project out occupied space
			gamma[p] -= rho0(gamma[p]);
		}
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
	//		real_function_3d vxc = xclib_interface_.convolution_with_kernel(
	//				perturbed_density);

	//		for (size_t i = 0; i < gamma.size(); i++) {
	//			gamma[i] += vxc * active_mo_[i];
	//		}

	// Alternative way (more expensive, but avoid the unprecise kernel)
	// for small test molecules this seems to bring no improvement
	// when using this:
	// 1.comment out the line before
	// 2.return add(world,gamma,gamma2)
	// 3. dont forget to project out occupied space also from gamma2 (below here)
	// THIS DOES NOT WORK FOR GGA
	//	vecfuncT gamma2=xclib_interface_.apply_kernel(xfunction.x);
	//	for (int p=0; p<active_mo_.size(); ++p) gamma2[p] -= rho0(gamma2[p]);

	//	vecfuncT gamma2 = mul(world,perturbed_density,lda_intermediate_);

	for(size_t p=0;p<gamma.size();p++){
		gamma[p] += xclib_interface_.get_fxc()*perturbed_density*active_mo_[p];
	}
	plot_vecfunction(gamma,"complete_gamma_");

	// project out occupied space
	for (size_t p = 0; p < active_mo_.size(); ++p)
		gamma[p] -= rho0(gamma[p]);

	//applyit.info(debug_);
	//plot_vecfunction(gamma, "gamma", plot_);
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

std::vector<vecfuncT> TDA::make_localized_exchange_intermediate() const {
	vecfuncT lmo=madness::transform(world,active_mo_,mo2lmo_,true);

	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr<
			real_convolution_3d>(CoulombOperatorPtr(world, lo, bsh_eps_));

	std::vector<vecfuncT> intermediate(lmo.size());

	for (std::size_t p = 0; p < lmo.size(); ++p) {
		intermediate[p] = apply(world, (*poisson),
				mul(world, lmo[p], lmo));
	}
	return intermediate;
}

vecfuncT TDA::make_lda_intermediate()const{
	vecfuncT mo;
	for(size_t i=0;i<active_mo_.size();i++) mo.push_back(copy(active_mo_[i]));
	vecfuncT result = xclib_interface_.get_lda_intermediate(mo);
	plot_vecfunction(result,"lda_intermediate_");
	mo.clear();
	return result;
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
			xfunctions[k].x[i].truncate(truncate_thresh_);
			// The potential or residual vectors can be empty, not failsafe for the x-vector because it should never be empty
			if (not xfunctions[k].Vx.empty())
				xfunctions[k].Vx[i].truncate(truncate_thresh_);
			if (not xfunctions[k].current_residuals.empty())
				xfunctions[k].current_residuals[i].truncate(truncate_thresh_);
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
	// get the overlap with the guess exops
	guess guess_structure(world,calc_.param.L,guess_mode_,active_mo_,calc_.ao);
	double tol = 10.0 * FunctionDefaults<3>::get_thresh();
	for(size_t i=0;i<roots.size();i++){
		guess_structure.project_to_guess_basis(world,roots[i].x,tol);
	}

	//	// get the overlaps with exops
	//	exoperators exops(world,excitation_point_);
	//	for(size_t i=0;i<roots.size();i++){
	//		functorT smooth_functor(
	//				new guess_smoothing(guess_box_));
	//		real_function_3d smoothing_function = real_factory_3d(world).functor(smooth_functor);
	//
	//		// TESTING DEBUG
	//		// project the active_mos on ao functions
	//		//		vecfuncT projected_mos = zero_functions<double,3>(world,active_mo_.size());
	//		//		for(size_t i=0;i<active_mo_.size();i++){
	//		//			projected_mos[i] = exops.project_function_on_aos(world,active_mo_[i]);
	//		//		}
	//
	//		std::vector<double> overlap_tmp = exops.get_overlaps_with_guess(world,roots[i].x,active_mo_,smoothing_function);
	//		std::vector<std::string> key = exops.key_;
	//		if(world.rank()==0){
	//			//			std::cout <<"\n\n----excitation "<< i << "----"<< std::endl;
	//			//			std::cout <<"\n dipole contributions"<< std::endl;
	//			//			for(size_t k=0;k<3;k++) std::cout << std::fixed << std::setprecision(2) << key[k]<<" " <<overlap_tmp[k]<<" ";
	//			//			std::cout <<"\n quadrupole contributions"<< std::endl;
	//			//			for(size_t k=3;k<9;k++) std::cout << std::fixed << std::setprecision(2) << key[k]<<" "<< overlap_tmp[k]<<" ";
	//			//			std::cout <<"\n cubic contributions"<< std::endl;
	//			//			for(size_t k=9;k<19;k++) std::cout << std::fixed << std::setprecision(2) << key[k]<<" "<< overlap_tmp[k] << " ";
	//
	//			std::cout <<"\n\n all significant contributions " << std::endl;
	//			for(size_t k=0;k<overlap_tmp.size();k++){
	//				if(fabs(overlap_tmp[k]) > 1.e-4) std::cout << std::fixed << std::setprecision(4) << key[k]<<" "<< overlap_tmp[k]<<" ";
	//			}
	//			std::cout << std::endl;
	//		}
	//		// Get overlap with for each mo
	//		//		for(size_t j=0;j<roots[i].x.size();j++){
	//		//			double xtmp_norm = roots[i].x[j].norm2();
	//		//			std::cout << "Contributions for excitation " << i << " and MO " << j << " ... norm is " << xtmp_norm << std::endl;
	//		//			vecfuncT xtmp; xtmp.push_back(roots[i].x[j]);
	//		//			vecfuncT motmp; motmp.push_back(active_mo_[j]);
	//		//
	//		//			std::vector<double> mo_overlap_tmp = exops.get_overlaps_with_guess(world,xtmp,motmp,smoothing_function);
	//		//			for(size_t k=0;k<mo_overlap_tmp.size();k++){
	//		//				if(world.rank()==0)if(fabs(mo_overlap_tmp[k]) > 1.e-4) std::cout << std::fixed << std::setprecision(4) << key[k]<<" "<< mo_overlap_tmp[k]<<" ";
	//		//			}
	//		//		}
	//	}
}

void TDA::save_xfunctions(const xfunctionsT &xfunctions)const{
	const std::string name_x = "xfunctions_current";
	const std::string name_xconv = "xfunctions_converged";
	// save current xfunctions
	for(size_t i=0;i<xfunctions.size();i++){
		std::string filename = name_x+stringify(i);
		archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
		ar & xfunctions[i].omega;
		ar & xfunctions[i].expectation_value;
		ar & xfunctions[i].delta;
		ar & xfunctions[i].number;
		ar & xfunctions[i].error;
		ar & xfunctions[i].converged;
		ar & xfunctions[i].iterations;
		for(size_t j=0;j<xfunctions[i].x.size();j++){ar & xfunctions[i].x[j];}
	}
	// save converged xfunctions
	if(not converged_xfunctions_.empty()){
		for(size_t i=0;i<converged_xfunctions_.size();i++){
			std::string filename = name_xconv+stringify(i);
			archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
			ar & converged_xfunctions_[i].omega;
			ar & converged_xfunctions_[i].expectation_value;
			ar & converged_xfunctions_[i].delta;
			ar & converged_xfunctions_[i].number;
			ar & converged_xfunctions_[i].error;
			ar & converged_xfunctions_[i].converged;
			ar & converged_xfunctions_[i].iterations;
			for(size_t j=0;j<converged_xfunctions_[i].x.size();j++){ar & converged_xfunctions_[i].x[j];}
		}
	}
}

bool TDA::read_xfunctions(xfunctionsT &xfunctions){
	size_t noct = active_mo_.size();
	const std::string name_x = "xfunctions_current";
	const std::string name_xconv = "xfunctions_converged";
	// check for converged_xfunctions
	std::cout << "check for saved converged functions " << std::endl;
	for(size_t i=0;i<guess_excitations_;i++){
		std::string filename = name_xconv+stringify(i);
		std::cout << "converged xfunction " << i;
		if(archive::ParallelInputArchive::exists(world,filename.c_str())){
			xfunction dummy(world,active_mo_);
			archive::ParallelInputArchive ar(world, filename.c_str(), 1);
			std::cout << "...found"<<std::endl;
			ar & dummy.omega;
			ar & dummy.expectation_value;
			ar & dummy.delta;
			ar & dummy.number;
			ar & dummy.error;
			ar & dummy.converged;
			ar & dummy.iterations;
			for(size_t j=0;j<noct;j++){ar & dummy.x[j];}
			converged_xfunctions_.push_back(dummy);
		}else std::cout << "...not found" << std::endl;
	}

	// check for unconverged xfunctions
	for(size_t i=0;i<guess_excitations_;i++){
		std::string filename = name_x+stringify(i);
		std::cout << "xfunction " << i;
		if(archive::ParallelInputArchive::exists(world,filename.c_str())){
			xfunction dummy(world,active_mo_);
			archive::ParallelInputArchive ar(world, filename.c_str(), 1);
			std::cout << "...found"<<std::endl;
			ar & dummy.omega;
			ar & dummy.expectation_value;
			ar & dummy.delta;
			ar & dummy.number;
			ar & dummy.error;
			ar & dummy.converged;
			ar & dummy.iterations;
			for(size_t j=0;j<noct;j++){ar & dummy.x[j];}
			xfunctions.push_back(dummy);
		}else std::cout << "...not found" << std::endl;
	}
	if(not xfunctions.empty()) return true;
	else return false;
}
