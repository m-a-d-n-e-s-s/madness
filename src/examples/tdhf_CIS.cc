/*
 * tdhfCIS.cc
 *
 *  Created on: May 5, 2014
 *      Author: kottmanj
 */

#include "tdhf_CIS.h"

using namespace madness;


static double ttt, sss;
void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss;
    if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

/// Print information of root vector
void CIS::print_roots(const std::vector<root> &roots,const int iter) const{

	if (world.rank()==0) {
		print(" root   excitation energy   energy correction     error		converged");
		for(size_t i=0;i<roots.size();i++){
			// $ for better grep
			std::cout << " $ " << i << " $ " ;
			print_root(roots[i]);
			std::cout << "  iter:" <<iter << std::endl;
		}
		std::cout << std::endl;
	}

}

void CIS::print_roots(const std::vector<root> &roots) const{

	if (world.rank()==0) {
		print(" root   excitation energy   energy correction     error		converged");
		for(size_t i=0;i<roots.size();i++){
			// $ for better grep
			std::cout << " " << i << " "; ;
			print_root(roots[i]);
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

}

void CIS::print_root(const root &root) const{
	if (world.rank()==0) {
		std::cout << std::scientific << std::setprecision(10);
		std::cout << std::setw(20) << root.omega	<< std::setw(20)<< root.delta
				<< std::setw(20) << root.err << "   " << root.converged;
	}
}

void CIS::sort_roots(std::vector<root> & roots){
	if(world.rank()==0) print("Sorting roots ...");
	std::sort(roots.begin(),roots.end(),compare_roots);
}


/// solve the CIS equations for n roots
void CIS::solve() {
	// Set convergence criteria for the actual CIS calculation
	// tight: Final iterations with second order and Gram-Schmidt
	// loose: Pre iteration with Fock Matrix

	prec tight;
	tight.thresh=thresh_;
	tight.econv=econv_;
	tight.dconv=dconv_;

	prec loose;
	loose.thresh=guess_thresh_;
	loose.econv=guess_econv_;
	loose.dconv=guess_dconv_;

	prec MO_guess;
	MO_guess.thresh=guess_thresh_;
	MO_guess.econv=guess_econv_;
	MO_guess.dconv=guess_dconv_;


	//plot the MOS
	vecfuncT MOs = get_calc().amo;
	for(size_t i=0;i<MOs.size();i++)plot_plane(world,MOs[i],"MO_"+stringify(i));

	set_prot(world,loose);

	// The roots vector
	std::vector<root> roots;
	initialize_roots(world,roots);

	if(read_and_save_koala_ == true) return;

	// Pre iteration with loose prec, Add noise if demanded
	if(guess_ == "physical"){
		if(world.rank()==0)print("Starting iteration of all ",guess_roots_," Guess roots with Fock update and loose convergence");
		if(noise_ == true){
			if(world.rank()==0)print("Noise will be added all ", guess_iter_," iterations");
			for(int cycle=0;cycle<40;cycle++){
				add_noise(world,roots);
				solve_internal_par("fock",roots,guess_iter_);
			}
		}
		else solve_internal_par("fock",roots,guess_iter_);
	}


	// Iterate with lose convergence criteria
	int cycles_ = 10; // make this member variable later
	for(int cycles=0; cycles<cycles_;cycles++){
	solve_internal_par("second order",roots,iter_max_);
	for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}
	solve_internal_par("second order",roots,1);
	if(check_convergence(roots))break;
	}



	//change to tight convergence criteria and second order update (use kain)
	set_prot(world,tight);
	if(world.rank()==0)print("Changing convergence criteria to tight");// and the update mode to second order - KAIN will be used form now on");
	// Iterate all roots of interest till convergence (converged roots are not iterated anymore)
	if(world.rank()==0) print("Iterate only nroot_ roots",nroot_);
	// erase all roots that are not of interest
	int old_size = roots.size();
	for(int i=nroot_;i<old_size;i++){roots.pop_back();}

	// Convergence criteria changed: set all roots back to non-converged (else the first iteration will do nothing)
	for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}

	// Iterate with tight convergence criteria
	if(world.rank()==0) print("using KAIN now");
	solve_internal_par("KAIN",roots,iter_max_*cycles_);
	for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}
	solve_internal_par("second order",roots,1);



	analyze(roots);

	// plot final roots
	for(size_t iroot=0;iroot<roots.size();iroot++){
		for(size_t i=0;i<roots[iroot].x.size();i++){
			plot_plane(world,roots[iroot].x[i],"Final_root"+stringify(iroot)+"_"+stringify(i));
		}
	}

}


bool CIS::solve_internal_par(const std::string mode,std::vector<root> &roots,const int iter_max) {


	// for convenience
	const vecfuncT& amo=get_calc().amo;
	const std::size_t nmo=amo.size();
	const int noct=nmo-nfreeze_;	// # of active orbitals

	// this is common for all roots in all iterations
	exchange_intermediate_=make_exchange_intermediate(active_mo(),active_mo());


	// one KAIN solver for each root
	typedef  XNonlinearSolver<root,double,allocator> solverT;
	solverT onesolver(allocator(world,noct));
	onesolver.set_maxsub(3);
	//!!!!!!!!!!!!!!! for now this works because the second order update uses KAIN and fock update doesnt
	// if you want to update guess roots (more than nroot_) with KAIN this will fail
	std::vector<solverT> solver(nroot_,onesolver);

	bool converged=iterate_all_CIS_roots(world,solver,roots,mode,iter_max);
	if (converged) {
		//analyze(roots);
		if (world.rank()==0) print(" CIS iterations converged ");
		return true;
	} else {
		if (world.rank()==0) print(" CIS iterations not converged ");
		return false;
	}

	// We should never be here
	print("END OF SOLVE_INTERNAL_PAR REACHED .... this shouldnt happen");
	MADNESS_EXCEPTION("Reached end of solve_internal_par",1);
	return false;
}

//*/
/// return the roots of the response equation
std::vector<root>& CIS::roots() {return roots_;}

/// are we supposed to print the grid for an external guess
bool CIS::print_grid() const {return print_grid_;}

// Excitation functions for the guess (x,y,z,r,x^2,xy,xz ...)

std::vector<real_function_3d> CIS::excitation_functions(const std::string exf) const {

	real_function_3d fx = real_factory_3d(world).f(x);
	fx.truncate();
	real_function_3d fy = real_factory_3d(world).f(y);
	fy.truncate();
	real_function_3d fz = real_factory_3d(world).f(z);
	fz.truncate();
	real_function_3d fr = real_factory_3d(world).f(rfunction);
	fz.truncate();
	real_function_3d fr2 = real_factory_3d(world).f(rfunction2);
	fz.truncate();


	std::vector<real_function_3d> exfunctions;

	if(exf == "dipole"){

		//2p
		exfunctions.push_back(fx);
		exfunctions.push_back(fy);
		exfunctions.push_back(fz);
		//2s
		exfunctions.push_back(fr2);
	}


	//	//3s (1+r)+r^2
	//	exfunctions.push_back(fr2+fr*fr);
	//	//3p as 2s*x
	//	exfunctions.push_back(fr2*fx);
	//	exfunctions.push_back(fr2*fy);
	//	exfunctions.push_back(fr2*fz);

	if(exf == "quadrupole"){
		//3d
		exfunctions.push_back(fx*fx);
		exfunctions.push_back(fx*fy);
		exfunctions.push_back(fx*fz);
		exfunctions.push_back(fy*fy);
		exfunctions.push_back(fy*fz);
		exfunctions.push_back(fz*fz);
	}
	//	//4d as 2s*x^2 (bzw xy,yz ...)
	//	exfunctions.push_back(fr2*fx*fx);
	//	exfunctions.push_back(fr2*fx*fy);
	//	exfunctions.push_back(fr2*fx*fz);
	//	exfunctions.push_back(fr2*fy*fz);
	//	exfunctions.push_back(fr2*fy*fy);
	//	exfunctions.push_back(fr2*fz*fz);

	if(exf == "benzene"){
		real_function_3d tmp = fx*fx + fy*fy;
		real_function_3d tmp2 = fx*fx-fy*fy;
		exfunctions.push_back(fx);
		exfunctions.push_back(fy);
		exfunctions.push_back(fz);
		exfunctions.push_back(tmp);
		exfunctions.push_back(fz *fz);
		exfunctions.push_back(fx*fz);
		exfunctions.push_back(fy*fz);
		exfunctions.push_back(tmp2);
		exfunctions.push_back(fx*fy);

	}

	return exfunctions;
}

void CIS::initialize_roots(World &world,std::vector<root> &roots){
	if(world.rank()==0) {
		printf("\n");print("Initialize roots...");print("guess is",guess_);print("guess omega is ",guess_omega_);
	}

	// Simplest guess: Start with the MOs iterate a few times and create the next guess root
	if(guess_ == "MO")guess_MO(world,roots);

	// Read saved roots from previous calculations
	else if(guess_ == "read")guess_read(world,roots);

	// Create a guess and iterate it a few times
	else if(guess_ =="physical")guess_physical(world,roots);

	// get the guess from koala
	else if (guess_=="koala")guess_koala(world,roots);

	else print("Reached the end of initialize_roots function ... this should not happen");
}

void CIS::guess_MO(World &world,std::vector<root> &roots){
	// for convenience
	const std::size_t nmo=get_calc().amo.size();
	const size_t noct=nmo-nfreeze_;
	root all_orbitals_root(world);
	all_orbitals_root.amplitudes_=std::vector<double>(nmo,1.0);


	// Create the all_orbitals guess

	root orbitals_root(world);
	real_function_3d all_orbitals=real_factory_3d(world);
	for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
		all_orbitals+=get_calc().ao[ivir];
		orbitals_root.x.push_back(get_calc().ao[ivir]);
	}


	for(size_t i=nfreeze_;i<nmo;i++){

		real_function_3d tmp = all_orbitals;
		double norm = tmp.norm2();
		tmp.scale(1.0/norm);
		all_orbitals_root.x.push_back(copy(tmp));
		if(plot_==true) plot_plane(world,tmp,"MO_guess_"+stringify(i));
	}

//	orbitals_root.omega = guess_omega_;
//	all_orbitals_root.omega = guess_omega_;


	// Preoptimize the all_orbital guess one by one
	for(int i=0;i<guess_roots_;i++){
		// Read roots if roots are there
		root tmp(world); bool check = false;
		tmp.x.resize(noct);
		check = load_root(world,i,tmp);
		if(check==true){
			if(world.rank()==0) print("Root ",i,"found and loaded");
			if(omega_[i]<100.0){
				if(world.rank()==0)print("Custom value for root ",i,"found (",omega_[i],")");
				tmp.omega = omega_[i];
			}
			else{
				if(world.rank()==0)print("Using default value for read root");
				//tmp.omega= -0.9*get_calc().aeps(nmo-1);
			}
			roots.push_back(tmp);
		}
		if(check==false){
			if(world.rank()==0) print("No saved root found, use all_orbital guess for next root ...");
			//if(guess_mode_=="all_orbitals")roots.push_back(all_orbitals_root);
			//if(guess_mode_=="orbitals")roots.push_back(orbitals_root);
			//if(guess_mode_=="noise")roots.push_back();
			roots.push_back(all_orbitals_root);
		}


		// these are the active orbitals in a vector (shallow-copied)
		//vecfuncT act=this->active_mo();
		//const vecfuncT& amo=get_calc().amo;
		// the projector on the unperturbed density
		//Projector<double,3> rho0(amo);
		//make_perturbed_fock_matrix(roots,act,rho0);



		if(guess_damp_ == true){
			if(check==false){
				std::vector<double> starting_energies;
				for(size_t iroot=0;iroot<roots.size();iroot++){ starting_energies.push_back(roots[iroot].omega);}
				if(world.rank()==0)print("Pre iteration without energy update:");
				for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}
				for(int instances=0;instances<guess_damp_iter_;instances++){

					// Fock Matrix update for energy guess
					solve_internal_par("fock",roots,1);

					if(world.rank()==0)print("Damped Iteration ",instances);
					solve_internal_par("second_order",roots,1);

					if(instances<2){
						for(size_t iroot=0;iroot<roots.size();iroot++){roots[iroot].omega=starting_energies[iroot];}
					}
					if(instances<4){
						for(size_t iroot=0;iroot<roots.size();iroot++){
							roots[iroot].omega=starting_energies[iroot]+0.25*roots[iroot].delta;
						}
					}
					if(instances<6){
						for(size_t iroot=0;iroot<roots.size();iroot++)
						{roots[iroot].omega=starting_energies[iroot]+0.50*roots[iroot].delta;
						}
					}
					if(instances<8){
						for(size_t iroot=0;iroot<roots.size();iroot++){
							roots[iroot].omega=starting_energies[iroot]+0.75*roots[iroot].delta;
						}
					}

				}



				// set back convergence
				for(size_t iroot=0;iroot<roots.size();iroot++){roots[iroot].converged = false;}
				solve_internal_par("second order",roots,guess_iter_);
				for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}

			}
		}
		// Iterate the guess MOs
		solve_internal_par("second order",roots,guess_iter_);



	}
}

void CIS::guess_read(World &world,std::vector<root> &roots){
	// for convenience
	const std::size_t nmo=get_calc().amo.size();
	const size_t noct=nmo-nfreeze_;

	for(int iroot=0;iroot<guess_roots_;iroot++){
		if(world.rank()==0) print("Trying to load root ",iroot," out of ",guess_roots_," roots");
		root tmp(world); bool check=false;
		// Correct size of X-Vector is necessary for the load_root function
		tmp.x.resize(noct);
		//Check if root is there
		check = load_root(world,iroot,tmp);
		if(check == true) {
			roots.push_back(tmp);
			if(world.rank()==0) print("Root ",iroot," found and loaded");
		}
		if(check == false){
			if(world.rank()==0) print("Root ",iroot," not found, use guess MO if you do not have all the data for the roots to read");
		}

	}
	if(world.rank()==0) print("Read roots are:");
	print_roots(roots);

	// Sort the roots
	sort_roots(roots);

}

void CIS::guess_koala(World &world,std::vector<root> &roots){

	// for convenience
	const std::size_t nmo=get_calc().amo.size();

	for(size_t iroot=0;iroot<guess_roots_;iroot++){

		vecfuncT koala_amo;
		root root(world);

		// first we need to determine the rotation of the external orbitals
		// to the MRA orbitals, so that we can subsequently rotate the guess
		if (not guess_phases_.has_data()) {

			// read koala's orbitals from disk
			for (std::size_t i=nfreeze_; i<nmo; ++i) {
				real_function_3d x_i=real_factory_3d(world).empty();
				const std::string valuefile="grid.koala.orbital"+stringify(i);
				x_i.get_impl()->read_grid2<3>(valuefile,functorT());
				koala_amo.push_back(x_i);
			}
			// this is the transformation matrix for the rotation
			guess_phases_=matrix_inner(world,koala_amo,get_calc().amo);
			guess_phases_=guess_phases_(_,Slice(nfreeze_,nmo-1));



		}

		// read the actual external guess from file
		for (std::size_t i=nfreeze_; i<nmo; ++i) {

			// this is the file where the guess is on disk
			const std::string valuefile="grid.koala.orbital"+stringify(i)
	    																					+".excitation"+stringify(iroot);
			real_function_3d x_i=real_factory_3d(world).empty();
			x_i.get_impl()->read_grid2<3>(valuefile,functorT());
			root.x.push_back(x_i);
		}

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

		// now rotate the active orbitals from the guess to conform with
		// the MRA orbitals
		//	    	Sinv=Sinv(Slice(nfreeze_,nmo-1),Slice(nfreeze_,nmo-1));
		root.x=transform(world,root.x,Sinv);
		save_root(world,iroot,root);

		if( omega_[iroot]<100.0) root.omega=omega_[iroot];
		else root.omega=get_calc().aeps(nmo-1);

		for(size_t i=0;i<roots.size();i++){
			for(size_t j=0;j<roots[i].x.size();j++){
				plot_plane(world,roots[i].x[j],"Koala_guess_root_"+stringify(i)+"_"+stringify(j));
			}
		}

		roots.push_back(root);
	}

	if(guess_damp_ == true){
		std::vector<double> starting_energies;
		for(size_t iroot=0;iroot<roots.size();iroot++){ starting_energies.push_back(roots[iroot].omega);}
		if(world.rank()==0)print("Pre iteration without energy update:");
		for(int instances=0;instances<guess_damp_iter_;instances++){
			solve_internal_par("fock",roots,1);
			for(size_t iroot=0;iroot<roots.size();iroot++){roots[iroot].omega=starting_energies[iroot];}

		}
	}
}

void CIS::guess_physical(World &world,std::vector<root> &roots){
	if(world.rank()==0) print("Guess is physical...");

	if(world.rank()==0)print("Starting Solve with lose convergence and Fock-Matrix update without solver");
	if(world.rank()==0)print("Guess is physical: x,y,z,r");
	if(world.rank()==0)print("20 pre iterations: x,y,z,r");

	// The guess roots
	std::vector<root> guess_roots=guess_big("dipole");

	int iterations = 2;

	// Create first 4 Guess fuctions (MO*x,y,z) and iterate them 5 times

	std::vector<root> guess = guess_roots;
	for(int i=0;i<4;i++){roots.push_back(guess[i]);}
	orthonormalize(world,roots);
	if(guess_damp_ == true){
		std::vector<double> starting_energies;
		for(size_t iroot=0;iroot<roots.size();iroot++){ starting_energies.push_back(roots[iroot].omega);}
		if(world.rank()==0)print("Pre iteration without energy update:");
		for(int instances=0;instances<guess_damp_iter_;instances++){
			solve_internal_par(0,roots,1);
			for(size_t iroot=0;iroot<roots.size();iroot++){roots[iroot].omega=starting_energies[iroot];}
		}
	}
	solve_internal_par("second order",roots,iterations+1);

	if(noise_ == true){
		if(world.rank()==0)print("Adding noise to the guess function ...");
		add_noise(world,roots);
	}
	// iterate the first 4 guess functions 3 times
	solve_internal_par("second order",roots,iterations+1);

	// Cont. if more than 4 roots are demanded for the guess
	if(guess_roots_>4){

		// Create the next 3+1 guess functions and pre iterate 5 times

		if(guess_exf_ == "dipole"){
			for(int i=0;i<4;i++){roots.push_back(guess[i]);}
		}
		else{
			if(world.rank()==0) print("Create quadrupole guess");
			std::vector<root> quadrupole_guess = guess_big("quadrupole");
			for(size_t iroot=0;iroot<quadrupole_guess.size();iroot++){
				roots.push_back(quadrupole_guess[iroot]);
			}
		}
		orthonormalize(world,roots);
		solve_internal_par("second order",roots,iterations);
	}

	// Save the guess roots if demanded
	if(guess_save_ == true){
		for(size_t iroot=0;iroot<roots.size();iroot++){
			save_root(world,iroot,roots[iroot],"Guess_");
		}
	}
}


void CIS::add_noise(World &world,std::vector<root> &roots)const{
	if(world.rank()==0)print("Adding noise to the X-functions, noise_box_ is ",noise_box_);

	int number=noise_gaussnumber_;
	double comp_factor=noise_comp_;
	real_function_3d noise_function=real_factory_3d(world);

	vecfuncT gauss;

	// Generate number random gauss functions and compress

	for(int i=0;i<number;i++){
		double size; double width=noise_width_;
		if(i<number/3.0) size=noise_box_; if(i<2.0*number/3.0) size=noise_box_*1.5; else size=noise_box_*2.0;
		functorT random_noise=functorT(new noise(size,width));
		real_function_3d tmp =real_factory_3d(world).functor(random_noise);
		double norm = tmp.norm2(); tmp.scale(comp_factor/norm);
		gauss.push_back(tmp);
		random_noise.reset();
	}

	// Sum the gauss functions
	for(size_t i=0;i<gauss.size();i++){
		if(plot_==true)plot_plane(world,gauss[i],"gauss"+stringify(i));
		noise_function+=gauss[i];
	}

	double noisenorm = noise_function.norm2();
	noise_function.scale(comp_factor/noisenorm);

	// Control output
	if(world.rank()==0)print("The norm of the added noise is: ", noisenorm);
	if(world.rank()==0)print("The norm of the rescaled noise function is: ", noise_function.norm2());

	// Add the noise to the x=functions of the roots which has not converged yet
	for(size_t iroot=0;iroot<roots.size();iroot++){
		for(size_t i=0;i<roots[iroot].x.size();i++){
			if(roots[iroot].converged ==false){
				roots[iroot].x[i] += noise_function;
				roots[iroot].x[i].truncate();
			}
		}
	}

	if(plot_==true)plot_plane(world,noise_function,"noise");

	orthonormalize(world,roots);

}



// Creates an (xyz and r)*MO guess
std::vector<root> CIS::guess_big(const std::string exf){
	if(world.rank()==0){
		printf("\n \n");print("Create guess functions...");
		print("guess_mode is: ",guess_mode_);
		printf("\n\n");
	}
	// default empty root vector
	std::vector<root> roots;

	// for convenience
	const std::size_t nmo=get_calc().amo.size();	// all orbitals

	// Sort the eigenvalues from the HF calculation and pick the smallest
	std::vector<double> orbital_energies;
	if(world.rank()==0)printf("\n----------------------------------------\n");
	if(world.rank()==0)print("Molecular Orbital energies will be sorted");
	if(world.rank()==0)print("They are:");
	for (std::size_t i=0; i<get_calc().amo.size(); ++i){
		double tmp = get_calc().aeps(i);
		orbital_energies.push_back(tmp);
		if(world.rank()==0) print(tmp);
	}

	// Now sort
	std::sort(orbital_energies.begin(),orbital_energies.end());
	double guess_energy=-0.9*orbital_energies.back();
	// Be sure
	if(world.rank()==0)printf("\n----------------------------------------\n");
	if(world.rank()==0)print("HOMO has energy: ",orbital_energies.back());
	if(world.rank()==0)printf("So the guess energy will be:%.3f \n\n",guess_energy);
	if(guess_energy<0){
		if(world.rank()==0) print("WARNING: Guess_Energy is negative");
	}

	// construct the projector on the occupied space
	const vecfuncT& amo=get_calc().amo;
	Projector<double,3> rho0(amo);

	real_function_3d all_orbitals=real_factory_3d(world);
	for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
		all_orbitals+=get_calc().ao[ivir];
	}

	std::vector<real_function_3d> exfunctions = excitation_functions(exf);

	// for all roots
	for(int iroot=0;iroot<4;iroot++){
		root root(world);
		root.amplitudes_=std::vector<double>(nmo,1.0);
		// For all MOs
		for(size_t i=nfreeze_;i<amo.size();i++){

			real_function_3d tmp = real_factory_3d(world);

			if(guess_mode_ == "mo"){
				tmp = get_calc().ao[iroot]*exfunctions[iroot];
			}

			if(guess_mode_ == "homo"){
				tmp = get_calc().ao[nmo-1]*exfunctions[iroot];
			}

			if(guess_mode_ == "all_orbitals"){
				tmp = all_orbitals*exfunctions[iroot];
			}
			// construct the projector on the occupied space
			const vecfuncT& amo=get_calc().amo;
			Projector<double,3> rho0(amo);
			// Project out the occupied space
			tmp -= rho0(tmp);
			double norm = tmp.norm2();
			tmp.scale(1.0/norm);
			root.x.push_back(copy(tmp));
			if(plot_==true) plot_plane(world,amo[i],"MO"+stringify(i));
			if(plot_==true) plot_plane(world,tmp,"Guess_root_"+stringify(iroot)+"_"+stringify(i));
		}

		// guess an excitation energy: 0.9* HOMO or use input from file (0.9 was too high)
		root.omega=guess_energy;
		roots.push_back(root);
	}

	orthonormalize(world,roots);

	return roots;
}


/// Just used to read roots
/// note that the orbitals from an external guess might be rotated wrt the
/// MRA orbitals, which would spoil the guess. Therefore we rotate the
/// guess with the overlap matrix so that the guess conforms with the MRA
/// phases. The overlap to the external orbitals (if applicable) is computed
/// only once
/// @param[in]	iroot	guess for root iroot
root CIS::guess_amplitudes(const int iroot) {

	// for convenience
	const std::size_t nmo=get_calc().amo.size();	// all orbitals
	const size_t noct=nmo-nfreeze_;						// active orbitals

	// default empty root
	root root(world);
	root.amplitudes_=std::vector<double>(nmo,1.0);

	// guess an excitation energy: 0.9* HOMO or use input from file
	root.omega=-0.9*get_calc().aeps(nmo-1);
	if (omega_[iroot]!=100.0) root.omega=omega_[iroot];

	// check if there's a root on disk, if so return those amplitudes
	root.x.resize(noct);
	if (load_root(world,iroot,root)) {
		for(size_t i=0;i<root.x.size();i++){
			root.x[i].set_thresh(thresh_);
		}
		return root;
	}
	root.x.clear();

	// get the guess from koala
	if (guess_=="koala") {

		vecfuncT koala_amo;

		// first we need to determine the rotation of the external orbitals
		// to the MRA orbitals, so that we can subsequently rotate the guess
		if (not guess_phases_.has_data()) {

			// read koala's orbitals from disk
			for (std::size_t i=nfreeze_; i<nmo; ++i) {
				real_function_3d x_i=real_factory_3d(world).empty();
				const std::string valuefile="grid.koala.orbital"+stringify(i);
				x_i.get_impl()->read_grid2<3>(valuefile,functorT());
				koala_amo.push_back(x_i);
			}
			// this is the transformation matrix for the rotation
			guess_phases_=matrix_inner(world,koala_amo,get_calc().amo);
			guess_phases_=guess_phases_(_,Slice(nfreeze_,nmo-1));

		}

		// read the actual external guess from file
		for (std::size_t i=nfreeze_; i<nmo; ++i) {

			// this is the file where the guess is on disk
			const std::string valuefile="grid.koala.orbital"+stringify(i)
    	    																				+".excitation"+stringify(iroot);
			real_function_3d x_i=real_factory_3d(world).empty();
			x_i.get_impl()->read_grid2<3>(valuefile,functorT());
			root.x.push_back(x_i);
		}

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

		// now rotate the active orbitals from the guess to conform with
		// the MRA orbitals
		//	    	Sinv=Sinv(Slice(nfreeze_,nmo-1),Slice(nfreeze_,nmo-1));
		root.x=transform(world,root.x,Sinv);
		save_root(world,iroot,root);

		// assign the excitation energies (need to be given in the input as omega0 omega1 ... etc)
		if(omega_.size() >= guess_roots_){
			root.omega = omega_[iroot];
		}
		else{root.omega = get_calc().aeps(nmo-1);}


	} else if (guess_=="all_orbitals") {
		/// New Try

		if(world.rank()==0){
			print("Try new guess 2 with all orbitals!!!....");}
		// construct the projector on the occupied space
		const vecfuncT& amo=get_calc().amo;
		Projector<double,3> rho0(amo);

		real_function_3d all_orbitals=real_factory_3d(world);
		for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
			all_orbitals+=get_calc().ao[ivir];
		}
		real_function_3d fx = real_factory_3d(world).f(x);
		fx.truncate();
		real_function_3d fy = real_factory_3d(world).f(y);
		fy.truncate();
		real_function_3d fz = real_factory_3d(world).f(z);
		fz.truncate();
		real_function_3d fr = real_factory_3d(world).f(rfunction);
		fz.truncate();
		real_function_3d fr2 = real_factory_3d(world).f(rfunction2);
		fz.truncate();


		std::vector<real_function_3d> exfunctions;
		//2p
		exfunctions.push_back(fx);
		exfunctions.push_back(fy);
		exfunctions.push_back(fz);
		//2s
		exfunctions.push_back(fr2);
		//3d
		exfunctions.push_back(fx*fx);
		exfunctions.push_back(fx*fy);
		exfunctions.push_back(fx*fz);
		exfunctions.push_back(fy*fy);
		exfunctions.push_back(fy*fz);
		exfunctions.push_back(fz*fz);
		//3s (1+r)+r^2
		exfunctions.push_back(fr2+fr*fr);
		//3p as 2s*x
		exfunctions.push_back(fr2*fx);
		exfunctions.push_back(fr2*fy);
		exfunctions.push_back(fr2*fz);
		//3d as 2s*x^2 (bzw xy,yz ...)
		exfunctions.push_back(fr2*fx*fx);
		exfunctions.push_back(fr2*fx*fy);
		exfunctions.push_back(fr2*fx*fz);
		exfunctions.push_back(fr2*fy*fz);
		exfunctions.push_back(fr2*fy*fy);
		exfunctions.push_back(fr2*fz*fz);


		// For all Molecules
		for(size_t i=nfreeze_;i<amo.size();i++){

			real_function_3d tmp = all_orbitals*exfunctions[iroot];
			if(plot_==true) plot_plane(world,tmp,"Guess_root_"+stringify(iroot)+"_"+stringify(i)+"before_P");

			// construct the projector on the occupied space
			const vecfuncT& amo=get_calc().amo;
			Projector<double,3> rho0(amo);
			// Project out the occupied space
			tmp -= rho0(tmp);
			double norm = tmp.norm2();
			tmp.scale(1.0/norm);
			root.x.push_back(copy(tmp));
			if(plot_==true) plot_plane(world,amo[i],"MO"+stringify(i));
			if(plot_==true) plot_plane(world,tmp,"Guess_root_"+stringify(iroot)+"_"+stringify(i));
		}



	} else {
		MADNESS_EXCEPTION("unknown source to guess CIS amplitudes",1);
	}
	MADNESS_ASSERT(root.x.size()==noct);

	// normalize this root
	normalize(world,root);

	return root;
}

/// @param[in]	world	the world
/// @param[in]	solver	the KAIN solver (unused right now)
/// @param[inout]	roots	on entry: guess for the roots
///                         on successful exit: converged root
/// @param[guess] for guess optimization with shifted HF
/// @return	convergence reached or not

template<typename solverT>
bool CIS::iterate_all_CIS_roots(World& world, std::vector<solverT>& solver,
		std::vector<root>& roots,const std::string mode,const int iter_max) const {

	// check orthogonality of the roots
	orthonormalize(world,roots);
	Tensor<double> ovlp=overlap(roots,roots);
	for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
	if (ovlp.normf()/ovlp.size()>econv_) {
		print(ovlp);
		MADNESS_EXCEPTION("non-orthogonal roots",1);
	}



	// start iterations
	//Iterate theif(guess==true) iter_max=guess_iter_hf_;
	for (int iter=0; iter<iter_max; ++iter) {

		// Add some noise (fock==true: only in guess iteration))
		//if(noise_==true) {if(fock_==true) add_noise(world,roots);}

		std::vector<double> error(roots.size());
		// print progress on the computation
		if (world.rank()==0) {
			printf("\n\n-----------------------------------------------\n\n");
			print("starting iteration ", iter, " at time ", wall_time());
			print(" root   excitation energy   energy correction         error");
		}

		// apply the BSH operator on each root
		START_TIMER(world);
		for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {
			std::cout << std::setw(4) << iroot << " ";
			if(roots[iroot].converged != true){
				roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],mode);

			}
			if(plot_==true) {
				for(size_t i=0;i<roots[iroot].x.size();i++)
					plot_plane(world,roots[iroot].x[i],"root_"+stringify(iroot)+"_0_iter_"+stringify(iter));
			}

		}
		orthonormalize(world,roots);
		END_TIMER(world,"BSH step");


		// Orthonormalization using Gram-Schmidt (fock_==false) or the perturbed fock-matrix (fock_==true)
		// Orthonormalization with the fock matrix updates the energy
		// if fock==false the energy has been updated during the second order step (if fock==true this has been suppresed)
		if(mode == "fock"){
			orthonormalize_fock(world,roots);
		}
		else{
			if(mode != "KAIN")std::sort(roots.begin(),roots.end(),compare_roots);
			orthonormalize(world,roots);
		}

		///save the roots
		for (std::size_t i=0; i<roots.size(); ++i) save_root(world,i,roots[i]);

		// check energy and density convergence (max residual norm)
		// global convergence is checked
		// convergence of single roots is checked
		if (world.rank()==0) printf("checking convergence ... \n   with econv:%f , dconv:%f , thresh:%f \n",econv_,dconv_,thresh_);
		bool converged=check_convergence(roots);

		//Print information

		print("");print("------Iteration-",iter,"---------------------------------");
		print_roots(roots,iter);print("");
		if (converged==true){print("converged!");return true;}
		if (world.rank()==0) print("not converged yet ...");
	}
	print("NOT CONVERGED");
	return false;
}

bool CIS::check_convergence(std::vector<root> &roots)const{
	bool converged=true;
			for (std::size_t i=0; i<roots.size(); ++i) {
				roots[i].converged = true;
				if (std::fabs(roots[i].delta)>econv_){
					converged=false;
					roots[i].converged = false;

				}
				if (roots[i].err>dconv_){
					converged = false;
					roots[i].converged = false;

				}

			}
	return converged;
}

/// follow Eq (4) of
/// T. Yanai, R. J. Harrison, and N. Handy,
/// Multiresolution quantum chemistry in multiwavelet bases: time-dependent
/// density functional theory with asymptotically corrected potentials in
/// local density and generalized gradient approximations,
/// Mol. Phys., vol. 103, no. 2, pp. 413-424, 2005.
///
/// The convergence criterion is that the excitation amplitudes don't change
/// @param[in]		world	the world
/// @param[in]		solver	the KAIN solver (not used right now..)
/// @param[inout]	root	the current root that we solve
/// @return			the residual error
template<typename solverT>
double CIS::iterate_one_CIS_root(World& world, solverT& solver, root& thisroot,const std::string mode) const {

	const vecfuncT& amo=get_calc().amo;
	const int nmo=amo.size();		// # of orbitals in the HF calculation
	const int noct=nmo-nfreeze_;	// # of active orbitals

	vecfuncT& x=thisroot.x;
	double& omega=thisroot.omega;

	// these are the active orbitals in a vector (shallow-copied)
	vecfuncT active_mo=this->active_mo();

	// the projector on the unperturbed density
	Projector<double,3> rho0(amo);

	// apply the unperturbed potential
	const vecfuncT V0=apply_fock_potential(x);

	// apply the gamma potential on the active MOs; Eq. (4-6)
	const vecfuncT Gamma=apply_gamma(x,active_mo,rho0);

	// add the local potential and the electron interaction term
	vecfuncT Vphi=add(world,V0,Gamma);
	scale(world,Vphi,-2.0);
	truncate(world,Vphi);


	// the bound-state helmholtz function for omega < orbital energy

	std::vector<poperatorT> bsh(noct);
	for(int p = 0; p<noct; ++p){
		double eps = get_calc().aeps[p+nfreeze_] + omega;	// orbital energy
		if(eps > 0){
			if(world.rank() == 0)
				print("bsh: warning: positive eigenvalue", p+nfreeze_, eps);
			eps = -0.03;
		}
		bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), lo, bsh_eps_));
	}
	const vecfuncT GVphi=apply(world,bsh,Vphi);

	// compute the residual aka the difference between to old and the new solution vector
	const vecfuncT residual=sub(world,x,GVphi);

	// update the excitation energy omega

	Tensor<double> t1=inner(world,Vphi,residual);
	double t2=0.0;
	for (int i=0; i<noct; ++i) {
		double n=GVphi[i].norm2();
		t2+=n*n;
	}

	// remove factor 2 from Vphi coming from the BSH application
	double delta=0.5*t1.sum()/t2;
	thisroot.delta = delta;

	double error=sqrt(inner(root(world,residual),root(world,residual)));

	// some update on the progress for the user
	if (world.rank()==0) {
		std::cout << std::scientific << std::setprecision(10);
		std::cout << std::setw(20) << omega	<< std::setw(20)<< delta
				<< std::setw(19) << error << std::endl;

	}


	// Update Energy and X-Function
	// Use the solver only when using the second order update
	// For the Fock Matrix procedure do not use the solver because the X-Vectors will be sorted
	// update the x vector: orthogonalize against the occupied space
	if (mode == "fock") {
		x=GVphi;
	}else {
		omega+=delta;
		if(mode=="KAIN"){root ff=solver.update(root(world,x),root(world,residual));x=ff.x;}
		else {x=GVphi;}
	}

	for (int p=0; p<noct; ++p) x[p] -= rho0(GVphi[p]);

	// BYPASS THE KAIN SOLVER




	return error;
}




vecfuncT CIS::apply_gamma(const vecfuncT& x, const vecfuncT& act,
		const Projector<double,3>& rho0) const {

	// for convenience
	const std::size_t noct=act.size();

	// now construct the two-electron contribution Gamma, cf Eqs. (7,8)
	real_function_3d rhoprime=real_factory_3d(world);

	vecfuncT Gamma;
	if(triplet_ == false){
		// a poisson solver
		std::shared_ptr<real_convolution_3d> poisson
		=std::shared_ptr<real_convolution_3d>
		(CoulombOperatorPtr(world,lo,bsh_eps_));

		// the Coulomb part Eq. (7) and Eq. (3)
		for (size_t i=0; i<noct; ++i) rhoprime+=x[i]*act[i];
		real_function_3d pp=2.0*(*poisson)(rhoprime); // here is the factor 2

		// the Coulomb part Eq. (7) and Eq. (3)
		Gamma=mul(world,pp,act);
	}

	// set Gamma to zero
	if(triplet_==true){
		//create emtpy function
		real_function_3d empty = real_factory_3d(world);
		for (std::size_t p=0; p<noct; ++p) Gamma.push_back(empty);
	}

	// the exchange part Eq. (8)
	for (std::size_t p=0; p<noct; ++p) {

		// this is x_i * \int 1/r12 \phi_i \phi_p
		vecfuncT x_Ppi=mul(world,x,exchange_intermediate_[p]);
		for (std::size_t i=0; i<noct; ++i) Gamma[p]-=x_Ppi[i];

		// project out the zeroth-order density Eq. (4)
		Gamma[p]-=rho0(Gamma[p]);
	}
	return Gamma;
}

/// make the zeroth-order potential J-K+V_nuc

/// note the use of all (frozen & active) orbitals in the computation
/// @param[in]	x	the response vector
/// @return		(J-K+V_nuc) x
vecfuncT CIS::apply_fock_potential(const vecfuncT& x) const {

	// the local potential V^0 of Eq. (4)
	real_function_3d coulomb;
	real_function_3d vlocal = get_calc().potentialmanager->vnuclear() +
			get_coulomb_potential();

	// make the potential for V0*xp
	vecfuncT Vx=mul(world,vlocal,x);

	// and the exchange potential is K xp
	vecfuncT Kx=get_calc().apply_hf_exchange(world,get_calc().aocc,
			get_calc().amo,x);

	// sum up: V0 xp = V_loc xp - K xp // this is 2J - K (factor 2 is included in get_coulomb_potential())
	vecfuncT V0=sub(world,Vx,Kx);

	return V0;

}

/// return the Coulomb potential
real_function_3d CIS::get_coulomb_potential() const {
	MADNESS_ASSERT(get_calc().param.spin_restricted);
	if (coulomb_.is_initialized()) return copy(coulomb_);
	functionT rho = get_calc().make_density(world, get_calc().aocc,
			get_calc().amo).scale(2.0); // here is the factor 2
	coulomb_=get_calc().make_coulomb_potential(rho);
	return copy(coulomb_);
}

/// the intermediate is the same for all roots:
/// \f[
///   Int[i,p] = \int \frac{1}{r_{12}} \phi_i(1) * \phi_p(1)
/// \f]
/// both i and p are active MOs
/// @param[in]	active_mo	active orbitals in the CIS computation
/// @param[in]	amo			all MOs of the HF calculation
/// @return		a vector of vectors of functions: [noct][nocc]
std::vector<vecfuncT> CIS::make_exchange_intermediate(const vecfuncT& active_mo,
		const vecfuncT& amo) const {
	// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson
	=std::shared_ptr<real_convolution_3d>
	(CoulombOperatorPtr(world,lo,bsh_eps_));

	std::vector<vecfuncT> intermediate(amo.size());

	for (std::size_t p=0; p<active_mo.size(); ++p) {
		intermediate[p]=apply(world,(*poisson),mul(world,active_mo[p],amo));
	}
	return intermediate;
}

/// compute the perturbed fock matrix

/// the matrix is given by
/// \f[
///   F_{pq} = < x^p_i | F | x^q_i > + < x^p_i | \Gamma^q | \phi_i >
/// \f]
/// where an amplitude is given by its components x_i, such that and
/// summation over the occupied orbitals i is implied
/// \f[
///   < \vec x^p | \vec x^q > = \sum_i x^p_i x^q_i
/// \f]
/// and similar for the Fock matrix elements.
/// Note this is NOT the response matrix A
Tensor<double> CIS::make_perturbed_fock_matrix(const std::vector<root>& roots,
		const vecfuncT& act, const Projector<double,3>& rho0) const {

	const std::size_t nroot=roots.size();
	Tensor<double> Fock_pt(nroot,nroot);

	START_TIMER(world);
	// apply the unperturbed Fock operator and the Gamma potential on all
	// components of all roots
	for (std::size_t iroot=0; iroot<nroot; ++iroot) {

		const vecfuncT& xp=roots[iroot].x;
		// apply the unperturbed potential and the Gamma potential
		// on the response amplitude x^p
		const vecfuncT V0=apply_fock_potential(xp);
		const vecfuncT Gamma=apply_gamma(xp,act,rho0);

		const vecfuncT V=add(world,V0,Gamma);

		// compute the matrix element V_pq
		// <p | V | q> = \sum_i <x^p_i | V | x^q_i>
		for (std::size_t jroot=0; jroot<nroot; ++jroot) {
			const vecfuncT& xq=roots[jroot].x;
			Tensor<double> xpi_Vxqi=inner(world,xq,V);
			Fock_pt(iroot,jroot)=xpi_Vxqi.sum();
		}
	}
	END_TIMER(world,"Fock_pt, potential");
	START_TIMER(world);

	// add kinetic part
	for (std::size_t iroot=0; iroot<nroot; ++iroot) {
		const vecfuncT& xp=roots[iroot].x;
		reconstruct(world, xp);
	}

	std::vector< std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double,3>(world);

	for(int axis = 0;axis < 3;++axis) {

		std::vector<vecfuncT> dxp;
		for (std::size_t iroot=0; iroot<nroot; ++iroot) {

			const vecfuncT& xp=roots[iroot].x;
			const vecfuncT d = apply(world, *(gradop[axis]), xp);
			dxp.push_back(d);
		}
		for (std::size_t iroot=0; iroot<nroot; ++iroot) {
			for (std::size_t jroot=0; jroot<nroot; ++jroot) {
				Tensor<double> xpi_Txqi=inner(world,dxp[iroot],dxp[jroot]);
				Fock_pt(iroot,jroot)+=0.5*xpi_Txqi.sum();
			}
		}
	}
	END_TIMER(world,"Fock_pt, kinetic");
	return Fock_pt;
}

/// load a converged root from disk

/// @param[in]	world 	the world
/// @param[in]	iroot	the i-th root
/// @param[inout]	x	the x-vector for the i-th root
/// @param[out]	omega	the excitation energy
/// @return	successfully loaded a root or not
/// @param[in] filename_end end of filename (for guess roots)
bool CIS::load_root(World& world, const int i, root& root)  const {
	return load_root(world,i,root,"");
}
bool CIS::load_root(World& world, const int i, root& root,const std::string filename_end)  const {
	std::string filename="root_"+stringify(i)+filename_end;
	bool exists=archive::ParallelInputArchive::exists(world,filename.c_str());
	if (not exists) return false;

	archive::ParallelInputArchive ar(world, filename.c_str(), 1);
	ar & root.omega;
	for (std::size_t i=0; i<root.x.size(); ++i) ar & root.x[i];
	return true;
}

/// save a converged root to disk

/// @param[in]	world 	the world
/// @param[in]	iroot	the i-th root
/// @param[inout]	x	the x-vector for the i-th root
/// @param[in]	omega	the excitation energy
/// @param[in] filename for guess roots to be saved
void CIS::save_root(World& world, const int i, const root& root) const {
	save_root(world,i,root,"");
}
void CIS::save_root(World& world, const int i, const root& root, const std::string filename_end) const {
	std::string filename="root_"+stringify(i)+filename_end;
	archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
	ar & root.omega;
	for (std::size_t i=0; i<root.x.size(); ++i) ar & root.x[i];
}

/// normalize the excitation amplitudes

/// normalize the set of excitation amplitudes, such that the sum of square
/// of all amplitudes equals 1.
/// @param[in]		world the world
/// @param[inout]	x	the excitation vector
void CIS::normalize(World& world, root& x) const {
	const double n2=inner(x,x);
	scale(world,x.x,1.0/sqrt(n2));
}



/// orthonormalize all roots using the perturbed fock matrix (Energy update included)
void CIS::orthonormalize_fock(World &world,std::vector<root> &roots)const{

	// construct the projector on the occupied space
	const vecfuncT& amo=get_calc().amo;
	vecfuncT act=this->active_mo();
	const std::size_t noct=act.size();
	Projector<double,3> rho0(amo);

	// compute the Fock matrix elements
	Tensor<double> Fock_pt=make_perturbed_fock_matrix(roots,act,rho0);

	// add the term with the orbital energies:
	//  F_pq -= \sum_i\epsilon_i <x^p_i | x^q_i>
	for (std::size_t i=0; i<roots.size(); ++i) {
		for (std::size_t j=0; j<roots.size(); ++j) {
			Tensor<double> eij=inner(world,roots[i].x,roots[j].x);
			for (size_t ii=0; ii<noct; ++ii) {
				Fock_pt(i,j)-=get_calc().aeps[ii+nfreeze_]*eij[ii];
			}
		}
	}

	// diagonalize the roots in their subspace
	Tensor<double> U, evals;
	syev(Fock_pt,U,evals);

	// Print out the details:
	printf("\n Perturbed Fock Matrix eigenvectors and values \n\n");
	print("eval(F)",evals);
	print("evec(F)"); std::cout << U; std::cout<<std::endl;


	// diagonalize the amplitudes in the space of the perturbed Fock
	// matrix
	std::vector<vecfuncT> vc(roots.size());
	for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {
		vc[iroot]=zero_functions<double,3>(world,noct);
		compress(world, vc[iroot]);
		compress(world,roots[iroot].x);
	}


	for (size_t i=0; i<roots.size(); ++i) {
		for (size_t j=0; j<roots.size(); ++j) {
			gaxpy(world,1.0,vc[i],U(j,i),roots[j].x);
		}
	}
	if (world.rank()==0){printf("\n\n******Update Process from perturbed Fock Matrix**************\n");}
	for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {
		roots[iroot].delta =evals[iroot] - roots[iroot].omega;
		roots[iroot].x=vc[iroot];
		normalize(world,roots[iroot]);
		// Update the Energy
		roots[iroot].omega=evals[iroot];
	}
	if (world.rank()==0)print_roots(roots);
	if (world.rank()==0){printf("*****************************************************************\n\n\n");}
}


/// orthonormalize all roots using Gram-Schmidt
void CIS::orthonormalize(World& world, std::vector<root>& roots) const {

	// first normalize
	for (std::size_t r=0; r<roots.size(); ++r) {
		normalize(world,roots[r]);
	}

	// project out the converged roots
	root converged_roots(world);
	for(std::size_t i=0;i<roots.size();i++){
		if(roots[i].converged==false){
			for(std::size_t j=0;j<roots.size();j++){
				vecfuncT& ucx =roots[i].x;
				if(roots[j].converged==true){
					vecfuncT& cx=roots[j].x;
					const double overlap=inner(world,cx,ucx).sum();
					compress(world,cx,false);
					compress(world,ucx,true);

					for (unsigned int p=0; p<ucx.size(); ++p) {
						ucx[p].gaxpy(1.0, cx[p], -overlap, false);
					}
					world.gop.fence();

				}
			}
		}
	}


	// orthogonalize the unconverged roots wrt each other
	// Do not operate on converged roots
	for (std::size_t r=0; r<roots.size(); ++r) {
		if(roots[r].converged==false){
		vecfuncT& x=roots[r].x;
		for (std::size_t rr=0; rr<r; ++rr) {
			if(roots[rr].converged==false){

			const vecfuncT& lower=roots[rr].x;

			const double ovlp=inner(world,lower,x).sum();
			compress(world,lower,false);
			compress(world,x,true);

			for (unsigned int p=0; p<x.size(); ++p) {
				x[p].gaxpy(1.0, lower[p], -ovlp, false);
			}
			world.gop.fence();
		}
		normalize(world,roots[r]);
	}}
}}

/// compute the overlap between 2 sets of roots
Tensor<double> CIS::overlap(const std::vector<root>& r1,
		const std::vector<root>& r2) const {

	Tensor<double> ovlp(r1.size(),r2.size());
	for (std::size_t i=0; i<r1.size(); ++i) {
		const vecfuncT& x1=r1[i].x;
		for (std::size_t j=0; j<r2.size(); ++j) {
			const vecfuncT& x2=r2[j].x;
			ovlp(i,j)=inner(world,x1,x2).sum();
		}
	}
	return ovlp;
}

/// compute the oscillator strength in the length representation

/// the oscillator strength is given by
/// \f[
/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]	root	a converged root
double CIS::oscillator_strength_length(const root& root) const {
	Tensor<double> mu_if(3);
	for (int idim=0; idim<3; idim++) {
		real_function_3d ri = real_factory_3d(world).functor2(xyz(idim));
		vecfuncT amo_times_x=mul(world,ri,active_mo());
		Tensor<double> a=inner(world,amo_times_x,root.x);
		mu_if(idim)=a.sum();
	}
	const double f= 2.0/3.0 * root.omega * mu_if.sumsq() * 2.0;
	return f;
}

/// compute the oscillator strength in the velocity representation

/// the oscillator strength is given by
/// \f[
/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]	root	a converged root
double CIS::oscillator_strength_velocity(const root& root) const {
	Tensor<double> p_if(3);
	// compute the derivatives of the MOs in all 3 directions
	for (int idim=0; idim<3; idim++) {
		real_derivative_3d D = free_space_derivative<double,3>(world, idim);
		vecfuncT Damo=apply(world,D,active_mo());
		Tensor<double> a=inner(world,Damo,root.x);
		p_if(idim)=a.sum();
	}
	const double f= 2.0/(3.0 * root.omega) * p_if.sumsq() * 2.0;
	return f;
}

/// analyze the root: oscillator strength and contributions from occ
void CIS::analyze(const std::vector<root>& roots) const {

	const size_t noct=active_mo().size();

	std::vector<root>::const_iterator it;
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

/// return the active MOs only, note the shallow copy
const vecfuncT CIS::active_mo() const {
	const vecfuncT& amo=get_calc().amo;
	vecfuncT actmo;
	for (int i=nfreeze_; i<amo.size(); ++i) actmo.push_back(amo[i]);
	return actmo;
}

