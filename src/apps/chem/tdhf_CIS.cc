/*
 * tdhfCIS.cc
 *
 *  Created on: May 5, 2014
 *      Author: kottmanj
 */
/*!
  \file examples/tdhf_CIS.cc
  \brief definitions for CIS.h

  \class CIS
  \brief CIS class provides all necessary function to do a CIS calculation (currently only HF exchange)


  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

  ... moved to github (source tree is the same)

*/
#include <chem/tdhf_CIS.h>

using namespace madness;

// Timer
static double ttt, sss;
static void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

static void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct allocator {
    World& world;
    const int n;

    /// @param[in]	world	the world
    /// @param[in]	nn		the number of functions in a given vector
    allocator(World& world, const int nn) : world(world), n(nn) {}

    /// allocate a vector of n empty functions
    root operator()() {
        return root(world,zero_functions<double,3>(world,n));
    }
};

/// for convenience
// This is used by many functions in tdhf_CIS.cc even if eclipse says it is unused
static double inner(const root& a, const root& b) {
	if (a.x.size()==0) return 0.0;
	return madness::inner(a.x[0].world(),a.x,b.x).sum();
}

/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};



/// Print information of root vector
void CIS::print_roots(const std::vector<root> &roots,const int iter) const{

	// Create a copy so that it can be sorted
	std::vector<root> copyroot;
	for(size_t i=0;i<roots.size();i++){
		const root &a= roots[i];
		root tmp(a.world,a.omega,a.expv,a.delta,a.err,a.converged,a.iter,a.number);
		copyroot.push_back(tmp);
	}

	sort_roots(copyroot,"energy");

	if (world.rank()==0) {
		print(" root init     excitation energy     energy correction      error   all-iter   conv     expv");
		for(size_t i=0;i<copyroot.size();i++){
			// $ for better grep
			std::cout << " $" << i ;
			print_root(copyroot[i]);
			std::cout << "   current-iter:" <<iter << std::endl;
		}
		std::cout << std::endl;
	}

}

void CIS::print_roots(const std::vector<root> &roots) const{

	if (world.rank()==0) {
		print(" root   excitation energy   energy correction     error		converged");
		for(size_t i=0;i<roots.size();i++){
			// $ for better grep
			//std::cout << " " << i << " "; ;
			print_root(roots[i]);
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

}

void CIS::print_root(const root &root) const{
	if (world.rank()==0) {
		std::cout << std::scientific << std::setprecision(10);
		std::cout << std::setw(5)<< "&" << root.number << std::setw(20) << root.omega	<< std::setw(20)<< root.delta
				<< std::setw(20) << root.err;
		printf(" %i  %s   ",root.iter,root.converged ? "true" : "false");
		printf("   %.4f   ",root.expv);
	}
}

void CIS::sort_roots(std::vector<root> & roots,std::string criterium)const{
	// identified a possible bug (iter and number will not be sorted ... so maybe x will also not be sorted)
	// do nothing for now
	if(world.rank()==0) print("Sorting roots ...");
	if(criterium=="energy")std::sort(roots.begin(),roots.end());
	if(criterium=="error")std::sort(roots.begin(),roots.end(),compare_roots_error);
}

/// Read and analyze roots
void CIS::Analyze(){

	/// read_roots
	std::vector<root> roots;
	for(int i=0;i<nroot_;i++){
		root tmp(world); tmp.number=i;
		load_root(world,i,tmp);
		roots.push_back(tmp);
	}

	// Print out the norm of the roots
	for(size_t iroot=0;iroot<roots.size();iroot++){
		// Print out the energy of the roots
		print("Root",iroot," omega is ",roots[iroot].omega);
		print("roots.size() ",roots.size()," roots[",iroot,"].size() ",roots[iroot].x.size());
		std::vector<double> amplitudes;
		for(size_t i=0;i<roots[iroot].x.size();i++){
			double tmp = roots[iroot].x[i].norm2();
			amplitudes.push_back(tmp);
		}
		print("Amplitudes Root",iroot);
		std::cout << amplitudes;
		std::cout << std::endl;
	}



}

/// solve the CIS equations for n roots
void CIS::solve() {
	// Set convergence criteria for the actual CIS calculation
	// tight: Final iterations with second order and Gram-Schmidt
	// loose: Pre iteration with Fock Matrix


	//plot the MOS
	vecfuncT MOs = get_calc().amo;
	for(size_t i=0;i<MOs.size();i++)plot_plane(world,MOs[i],"MO_"+stringify(i));

	// The roots vector
	std::vector<root> roots;
	//set_prot(world,loose); // THIS MADE HUGE PROBLEMS (apparently .... ) dont know why
	initialize_roots(world,roots);


	if(read_and_save_koala_ == true) return;

	// Use only number of demanded roots
	for(size_t i=nroot_;i<roots.size();i++) roots.pop_back();

	// Iterate till convergence
	for(int cycle=0;cycle<100;cycle++){
		for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}
		solve_internal_par("expectation_value",roots,guess_iter_fock_);
		solve_internal_par("fock",roots,1);
		for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}
		bool converged=check_convergence(roots);
		if(converged==true) break;
	}
	solve_internal_par("expectation_value",roots,iter_max_);

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

	bool converged=false;
	converged = iterate_all_CIS_roots(world,solver,roots,mode,iter_max);

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

	if(exf == "quadrupole"){
		//3d
		exfunctions.push_back(fx*fx);
		exfunctions.push_back(fx*fy);
		exfunctions.push_back(fx*fz);
		exfunctions.push_back(fy*fy);
		exfunctions.push_back(fy*fz);
		exfunctions.push_back(fz*fz);
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

	/// Not recommended
	else if (guess_=="noise")guess_noise(world,roots);

	else if (guess_=="active_space" or guess_=="aspace") guess_aspace(world,roots);
	else if (guess_=="forced") guess_forced(world,roots);

	else if (guess_=="benzene") guess_benzene_custom(world,roots);

	else print("Reached the end of initialize_roots function ... this should not happen");
}

void CIS::guess_MO(World &world,std::vector<root> &roots){
	// for convenience
	const std::size_t nmo=get_calc().amo.size();
	const size_t noct=nmo-nfreeze_;
	root all_orbitals_root(world);
	all_orbitals_root.amplitudes_=std::vector<double>(nmo,1.0);

	// Test smthg
	for(size_t i=0;i<get_calc().ao.size();i++){
		plot_plane(world,get_calc().ao[i],"AO"+stringify(i));
	}
	for(size_t i=0;i<get_calc().amo.size();i++){
		plot_plane(world,get_calc().amo[i],"AMO"+stringify(i));
	}



	// Create the all_orbitals guess

	root orbitals_root(world);
	real_function_3d all_orbitals=real_factory_3d(world);
	if(guess_mode_ =="mo"){
		for (std::size_t ivir=0; ivir<get_calc().amo.size(); ++ivir) {
			all_orbitals+=get_calc().amo[ivir];
			orbitals_root.x.push_back(get_calc().amo[ivir]);
		}
	}
	if(guess_mode_=="ao"){
		for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
			all_orbitals+=get_calc().ao[ivir];
			orbitals_root.x.push_back(get_calc().ao[ivir]);
		}
	}
	//double size=1.0; double width=2.0;
	//functorT random_noise=functorT(new noise(size,width));
	//real_function_3d tmp =real_factory_3d(world).functor(random_noise);

	if(noise_==true){
		// Add a wide gauss function to the all_orbital guess
		double width = get_calc().molecule.bounding_cube();
		typedef std::shared_ptr<FunctionFunctorInterface<double,3> > functorT;
		functorT wide_gauss=functorT(new gauss_function(width));
		real_function_3d gauss = real_factory_3d(world).functor(wide_gauss);
		//gauss.scale(0.01);
		all_orbitals+=gauss;
	}

	for(size_t i=nfreeze_;i<nmo;i++){

		real_function_3d tmp = all_orbitals;
		double norm = tmp.norm2();
		tmp.scale(1.0/norm);
		all_orbitals_root.x.push_back(copy(tmp));
		if(plot_==true) plot_plane(world,tmp,"MO_guess_"+stringify(i));
	}

	// Calculate the Energy of the guess root using the perturbed fock matrix:

	// these are the active orbitals in a vector (shallow-copied)
	vecfuncT act=this->active_mo();
	const vecfuncT& amo=get_calc().amo;
	// the projector on the unperturbed density
	Projector<double,3> rho0(amo);
	// this is common for all roots in all iterations
	exchange_intermediate_=make_exchange_intermediate(active_mo(),active_mo());
	std::vector<root> groots;
	groots.push_back(all_orbitals_root);
	orthonormalize(world,groots);
	orthonormalize_fock(world,groots);
	all_orbitals_root = groots[0];

	//	orbitals_root.omega = guess_omega_;
	//	all_orbitals_root.omega = guess_omega_;


	// Preoptimize the all_orbital guess one by one
	int number_counter=0;
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
			tmp.number=number_counter;
			number_counter++;
			roots.push_back(tmp);
		}
		if(check==false){
			if(world.rank()==0) print("No saved root found, use all_orbital guess for next root ...");
			all_orbitals_root.number=number_counter;
			number_counter++;
			roots.push_back(all_orbitals_root);
		}

		// Iterate the guess MOs
		for(size_t iroot=0;iroot<roots.size();iroot++){roots[iroot].converged = false;}

			// Try to pull the roots to a lower level

			if(guess_pull_ ==true){
				solve_internal_par("pull",roots,guess_iter_);
				orthonormalize_fock(world,roots);
			}
			if(guess_pull_ ==false){
				solve_internal_par("fock",roots,1);
				solve_internal_par("expectation_value",roots,guess_iter_);
			}


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
			tmp.number=iroot;
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
	sort_roots(roots,"energy");

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
		root.number=iroot;
		roots.push_back(root);
	}

	// Print information about the read roots
	printf("\n\n\n");print("Roots from Koala are:");printf("\n");
	print_roots(roots);


	if(guess_mode_ == "fock"){
		solve_internal_par("fock",roots,guess_iter_);
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

	// The guess roots
	std::vector<root> guess_roots=guess_big("dipole");

	// Create first 4 Guess fuctions (MO*x,y,z) etc

	int number_counter=0;
	for(int instances=0;instances<100;instances++){
		std::vector<root> guess = guess_roots;
		int choice=3;
		if(guess_exf_ == "r") choice =4;
		for(int i=0;i<choice;i++){
			root tmp=guess[i];
			tmp.number=number_counter; tmp.iter=0;
			number_counter +=1;
			roots.push_back(tmp);
		}
		orthonormalize(world,roots);
		for(size_t i=0;i<roots.size();i++){roots[i].converged=false;}

		if(guess_pull_ ==true){
			// Try to pull the roots to a lower level
			for(int i=0;i<guess_iter_;i++){
			solve_internal_par("pull",roots,guess_iter_fock_);
			orthonormalize_fock(world,roots);
			}
		}
		if(guess_pull_ == false){
			solve_internal_par("expectation_value",roots,guess_iter_);
			solve_internal_par("fock",roots,1);
		}


		if(noise_ == true){
			if(world.rank()==0)print("Adding noise to the guess function ...");
			add_noise(world,roots);
		}
		if(roots.size()>=guess_roots_) break;
	}

	// Save the guess roots if demanded
	if(guess_save_ == true){
		for(size_t iroot=0;iroot<roots.size();iroot++){
			save_root(world,iroot,roots[iroot],"Guess_");
		}
	}
}

/// Not used guess functions (for tests)
void CIS::guess_noise(World &world,std::vector<root> & roots){


	// Make a diffuse gauss function and then add noise

	// This will not work for atoms (bounding cube is 0)
	double width=get_calc().molecule.bounding_cube();
	width=width*1.5;
	if(width==0) width=5.0;

	typedef std::shared_ptr<FunctionFunctorInterface<double,3> > functorT;
	functorT wide_gauss=functorT(new gauss_function(width));
	real_function_3d gauss = real_factory_3d(world).functor(wide_gauss);

	// Print it to test
	plot_plane(world,gauss,"wide_gauss_function_width_"+stringify(width));

	for(int iroot=0;iroot<guess_roots_;iroot++){
		root tmp(world);
		for(size_t i=nfreeze_;i<get_calc().amo.size();i++){
			real_function_3d mo_gauss = get_calc().amo[i] * gauss;
			tmp.x.push_back(mo_gauss);
			plot_plane(world,mo_gauss,"MO"+stringify(i)+"_gauss");
		}
		roots.push_back(tmp);

		//add_noise(world,roots);

		for(size_t i=0;i<roots[0].x.size();i++){plot_plane(world,roots[0].x[i],"Guess_root0_"+stringify(i));}

		// iterate one by one
		solve_internal_par("expectation_value",roots,guess_iter_);
	}

}

void CIS::guess_aspace(World &world,std::vector<root> & roots){
	//shortcuts
	const vecfuncT &mos = get_calc().amo; size_t nmo=mos.size();
	int iter=guess_iter_;

	// Print out information
	if(world.rank()==0){
		printf("\n\n********Active Space Guess********\n");
		print(" ", active_mo_, " active homos");
		print(" ", "Guess mode is: ",guess_mode_);
		if(active_mo_ ==nmo-nfreeze_) print("Active Space was not defined ... std may is too big !!!! ");
		printf("**********************************\n\n");

	}

	// Prepare the excitation function(s)
	std::vector<real_function_3d> exfunctions;
	if(guess_mode_ == "xyz" or guess_mode_ == "xyzr" or guess_mode_ == "r"){
		iter = 3;
		// Get the excitation functions: x,y,z,r, + quadrupoles (will not be used now)
		exfunctions=excitation_functions("dipole");
		if(guess_mode_ == "xyz") exfunctions.pop_back(); // delete r function
		if(guess_mode_ == "r"){exfunctions[3]=exfunctions[0]; exfunctions.pop_back(); exfunctions.pop_back(); iter=guess_iter_;}
	}
	if(guess_mode_ == "r2"){
		real_function_3d r2 = real_factory_3d(world).f(monopole);
		exfunctions.push_back(r2);
	}
	else if(guess_mode_ == "gauss"){
		iter = guess_iter_;
		// This will not work for atoms (bounding cube is 0)
		double width=get_calc().molecule.bounding_cube();
		width=width*1.5;
		if(width==0) width=5.0;

		typedef std::shared_ptr<FunctionFunctorInterface<double,3> > functorT;
		functorT wide_gauss=functorT(new gauss_function(width));
		real_function_3d gauss = real_factory_3d(world).functor(wide_gauss);
		exfunctions.push_back(gauss);
	}


	int number_counter =0;
	for(size_t instances=0;instances<exfunctions.size();instances++){
		for(int i=0;i<active_mo_;i++){

			print("Build a guess root");

			root tmp(world);
			tmp.number = number_counter; number_counter++;
			tmp.omega = guess_omega_;

			// Choose the right exfunction
			real_function_3d exfunction = real_factory_3d(world);

			exfunction = exfunctions[instances];

			// Use the active homos and multiply with the chosen exfunctions
			real_function_3d exmo = mos[nmo-1-i]*exfunction;

			// Normalize
			double norm =exmo.norm2();
			exmo.scale(1.0/norm);
			exmo.truncate();

			// Project out the occupied space
			Projector<double,3> rho0(mos);
			exmo -= rho0(exmo);

			// Normalize
			norm =exmo.norm2();
			exmo.scale(1.0/norm);
			exmo.truncate();

			// Put the "x*homo" function in the root x-vector
			for(size_t j=0;j<nmo-nfreeze_;j++){tmp.x.push_back(copy(exmo));}



			// Push back the tmp root and iterate before adding a new one
			roots.push_back(tmp);
			orthonormalize(world,roots);

			// debug
			//std::vector<root> guess_roots=guess_big("dipole");
			//roots[0].x=guess_roots[0].x;
			print("size of roots is ",roots.size());
			std::vector<double> normen;
			for(size_t k=0;k<roots[0].x.size();k++){
				real_function_3d asd = roots[0].x[k];
				double normus =asd.norm2();
				normen.push_back(normus);
			}
			print("Norm is "); std::cout << normen << std::endl;

			orthonormalize(world,roots);
			solve_internal_par("damped",roots,iter);
		}
		solve_internal_par("damped",roots,guess_iter_);
	}




}

void CIS::guess_forced(World &world,std::vector<root> & roots){

	// All Atomic Orbitals for all Symmetries
	real_function_3d allao = real_factory_3d(world);
	for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
		allao+=get_calc().ao[ivir];
	}

	// Normalize
	double norm =allao.norm2();
	allao.scale(1.0/norm);
	allao.truncate();

	// Project out the occupied space
	Projector<double,3> rho0(allao);
	allao -= rho0(allao);

	// Normalize
	norm =allao.norm2();
	allao.scale(1.0/norm);
	allao.truncate();

	// Create the guess root
	root tmp(world);
	tmp.iter=0; tmp.number=0; tmp.omega=guess_omega_;
	for(int k=0;k<get_calc().amo.size()-nfreeze_;k++){tmp.x.push_back(copy(allao));}

	for(size_t iroot=0;iroot<guess_roots_;iroot++){

		roots.push_back(tmp);
		orthonormalize(world,roots);


		// Iterate one time and set back the energy
		for(int i=0;i<guess_iter_;i++){
			double forced_energy;
			(iroot==0) ? forced_energy=guess_omega_ : forced_energy = roots[iroot-1].omega*1.1;
			solve_internal_par("second order",roots,1);
			if(roots[iroot].err < guess_dconv_ ) break;
			else roots[iroot].omega=forced_energy ;
		}
	}

}

void CIS::guess_benzene_custom(World &world,std::vector<root> & roots){

	// Try to create an b1u guess
	// use the cubic function: x(x^2-3y^2)
	// source: http://www.pci.tu-bs.de/aggericke/PC2/Punktgruppen/D6h.htm
	// Care with the orientation of moldft

	print("Entering Benzene Custom guess function for a b1u excitation");
	print("Be shure that you only try to optimize one to seven roots and that moldft did not reorient the molecule");

	vecfuncT exfunctions;

	real_function_3d f1 = real_factory_3d(world).f(b1u);
	exfunctions.push_back(f1);
	real_function_3d f2 = real_factory_3d(world).f(b2u);
	exfunctions.push_back(f2);

	real_function_3d f3 = real_factory_3d(world).f(e1g1); exfunctions.push_back(f3);
	real_function_3d f4 = real_factory_3d(world).f(e1g2); exfunctions.push_back(f4);
	real_function_3d f5 = real_factory_3d(world).f(a2u);  exfunctions.push_back(f5);
	real_function_3d f6 = real_factory_3d(world).f(e2u1);  exfunctions.push_back(f6);
	real_function_3d f7 = real_factory_3d(world).f(e2u2);  exfunctions.push_back(f7);

	// construct the projector on the occupied space
	const vecfuncT& amo=get_calc().amo;
	Projector<double,3> rho0(amo);

	// The molecular orbitals
	vecfuncT orbitals;
	size_t nmo = amo.size();
	print("nfreeze is ",nfreeze_," and amo.size ",nmo );
	for(size_t i=nfreeze_;i<nmo;i++){orbitals.push_back(amo[i]);}


	for(int j=0;j<guess_roots_;j++){
		root guess(world);
		guess.amplitudes_=std::vector<double>(nmo,1.0);
	for(size_t i=nfreeze_;i<amo.size();i++){
		real_function_3d tmp = exfunctions[j]*amo[i];

		// Project out the occupied space
		tmp -= rho0(tmp);

		double tmpnorm = tmp.norm2();
		tmp.scale(1.0/tmpnorm);
		tmp.truncate();
		guess.x.push_back(copy(tmp));
	}

	guess.number = j;
	guess.iter =0;
	guess.omega=0.23;
	normalize(world,guess);
	roots.push_back(guess);
	orthonormalize(world,roots);
	solve_internal_par("expectation_value",roots,1);
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
	for (std::size_t ivir=0; ivir<get_calc().amo.size(); ++ivir) {
		all_orbitals+=get_calc().amo[ivir];
	}

	std::vector<real_function_3d> exfunctions = excitation_functions(exf);

	// for all roots
	for(int iroot=0;iroot<4;iroot++){
		root root(world);
		root.amplitudes_=std::vector<double>(nmo,1.0);
		// For all MOs
		for(size_t i=nfreeze_;i<amo.size();i++){

			//create empty function
			real_function_3d tmp = real_factory_3d(world);

			if(guess_mode_ == "mo"){
				tmp = get_calc().amo[i]*exfunctions[iroot];
			}

			if(guess_mode_ =="ao"){
				real_function_3d aos=real_factory_3d(world);
				for (std::size_t ivir=0; ivir<get_calc().amo.size(); ++ivir) {
					aos+=get_calc().amo[ivir];

					tmp = aos*exfunctions[iroot];
				}
			}

			if(guess_mode_ == "homo"){
				tmp = get_calc().amo[nmo-1]*exfunctions[iroot];
			}

			if(guess_mode_ == "active_space"){
				for(int i=0;i<active_mo_;i++){
					tmp += get_calc().amo[nmo-1-i]*exfunctions[iroot];
				}
			}

			if(guess_mode_ == "all_orbitals"){
				tmp = all_orbitals*exfunctions[iroot];
			}

			if(guess_mode_ == "ALL"){
				real_function_3d aos=real_factory_3d(world);
				for (std::size_t ivir=0; ivir<get_calc().amo.size(); ++ivir) {
					aos+=get_calc().amo[ivir];
				}
				real_function_3d all = all_orbitals + aos;
				all.truncate();
				tmp = all * exfunctions[iroot];
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
	bool flip =false;
	for (int iter=0; iter<iter_max; ++iter) {

		// Add noise every second iteration

		if(noise_==true){
			if(flip==true){add_noise(world,roots);print("MAKE SOME NOISE"); flip=false;}
			else flip=true;
		}


		std::vector<double> error(roots.size());
		// print progress on the computation

		if (world.rank()==0) {
			printf("\n\n-----------------------------------------------\n\n");
			print("starting iteration ", iter," at time ", wall_time(), "(Update mode: ",mode,")");
		}

		if(mode == "pull"){

			//iterate
			if(guess_=="physical"){
			// Pull the first 3 roots down to guess_omega_, pull the others down to the values of the 3 before
				// dont do anything for root6 and higher

				for(size_t iroot=0;iroot<3;iroot++){

					roots[iroot].omega = guess_omega_;
					if(roots[iroot].err > guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"pull");
					if(roots[iroot].err < guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"expectation_value");
				}
				if(roots.size()>3){
				for(size_t iroot=3;iroot<6;iroot++){

					roots[iroot].omega = guess_omega_*1.1;
					if(roots[iroot].err > guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"pull");
					if(roots[iroot].err < guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"expectation_value");
				}
				}
				if(roots.size()>6){
				for(size_t iroot=6;iroot<roots.size();iroot++){

					roots[iroot].omega = guess_omega_*1.2;
					if(roots[iroot].err > guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"pull");
					if(roots[iroot].err < guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"expectation_value");
				}
				}


			}
			if(guess_=="MO"){
				sort_roots(roots,"energy");
				for(size_t iroot=0;iroot<roots.size();iroot++){
					if(iroot==0) roots[iroot].omega = guess_omega_;
					if(iroot!=0) roots[iroot].omega = 1.05*roots[iroot-1].omega;
					if(roots[iroot].err > guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"pull");
					if(roots[iroot].err < guess_dconv_)roots[iroot].err=iterate_one_CIS_root(world,solver[iroot],roots[iroot],"expectation_value");
				}
			}
			orthonormalize(world,roots);
		}


		if(mode!="pull"){

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
		}


		// Orthonormalization using Gram-Schmidt (fock_==false) or the perturbed fock-matrix (fock_==true)
		// Orthonormalization with the fock matrix updates the energy
		// if fock==false the energy has been updated during the second order step (if fock==true this has been suppresed)
		if(mode == "fock"){
			orthonormalize_fock(world,roots);
		}
		else{
			if(mode != "KAIN")sort_roots(roots,"error");
			orthonormalize(world,roots);
		}

		///save the roots
		for (std::size_t i=0; i<roots.size(); ++i) save_root(world,i,roots[i]);

		if (world.rank()==0) printf("checking convergence ... \n   with econv:%f , dconv:%f , thresh:%f \n",econv_,dconv_,thresh_);
		bool converged=check_convergence(roots);

		//Print information
		print("");print("------Iteration-",iter,"---------------------------------");
		print_roots(roots,iter);print("");
		if (converged==true){if(world.rank()==0)print("converged!");return true;}
		if (world.rank()==0) print("not converged yet ...");
	}
	if(world.rank()==0)print("NOT CONVERGED");
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
/// ÒMultiresolution quantum chemistry in multiwavelet bases: time-dependent
/// density functional theory with asymptotically corrected potentials in
/// local density and generalized gradient approximations,Ó
/// Mol. Phys., vol. 103, no. 2, pp. 413Ð424, 2005.
///
/// The convergence criterion is that the excitation amplitudes don't change
/// @param[in]		world	the world
/// @param[in]		solver	the KAIN solver (not used right now..)
/// @param[inout]	root	the current root that we solve
/// @return			the residual error
template<typename solverT>
double CIS::iterate_one_CIS_root(World& world, solverT& solver, root& thisroot,const std::string mode) const {
	if(world.rank()==0){ printf("\n---------------------\n");print("Iterate Root ",thisroot.number," mode: ",mode);printf("\n\n");}

	// Make the settings
	bool second_order = false;
	if(mode == "expectation_value"){
		// Switch to second order if the accuracy of the energy is beyond the scope of the expectation value
		(fabs(thisroot.delta)<guess_econv_) ? second_order =true : second_order=false;

	}
	if(mode == "damped"){
		(thisroot.err<guess_dconv_) ? second_order =true : second_order =false ;
	}
	if(mode == "second order"){
		second_order = true;
	}
	if(mode == "pull"){
		second_order = true;
	}


	// Update iteration number in root
	thisroot.iter+=1;

	const vecfuncT& amo=get_calc().amo;
	vecfuncT& x=thisroot.x;

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

	// Update the energy with expectation value in the beginning
	double exp_value=expectation_value(world,thisroot,Vphi);
	thisroot.expv=exp_value;
	double exp_delta = exp_value - thisroot.omega;
	if(world.rank()==0)print("delta: ",thisroot.delta);
	if(second_order==false and mode !="pull"){
		thisroot.omega=exp_value;
		thisroot.delta=exp_delta;
		if(fabs(exp_delta) < guess_econv_){second_order =true;}
		else{if(world.rank()==0)print("Expectation value is used");}
	}

	const int nmo=amo.size();		// # of orbitals in the HF calculation
	const int noct=nmo-nfreeze_;	// # of active orbitals

	double &omega=thisroot.omega;

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

	double error=sqrt(inner(root(world,residual),root(world,residual)));

	// some update on the progress for the user
	if (world.rank()==0) {
		print("BSH Step:");
		std::cout << thisroot.number;
		std::cout << std::scientific << std::setprecision(10);
		std::cout << std::setw(20) << omega	<< std::setw(20)<< delta
				<< std::setw(19) << error << std::endl;

	}


	// Update Energy and X-Function
	// Use the solver only when using the second order update
	// For the Fock Matrix procedure do not use the solver because the X-Vectors will be sorted
	// update the x vector: orthogonalize against the occupied space

	if (mode == "fock")x=GVphi;
	if(mode=="KAIN"){
		thisroot.delta = delta;
		omega+=delta;
		root ff=solver.update(root(world,x),root(world,residual));x=ff.x;
	}
	else if( second_order==true) {
		x=GVphi;
		print("Use second order update");
		thisroot.delta = delta;
		omega+=thisroot.delta;
	}
	else x=GVphi;
	if(mode == "pull") thisroot.omega = thisroot.expv;

	for (int p=0; p<noct; ++p) x[p] -= rho0(GVphi[p]);

	// BYPASS THE KAIN SOLVER




	return error;
	if(world.rank()==0){ printf("\n---------------------\n");}
}

vecfuncT CIS::gamma_update(World &world, const root root)const{

	// construct the projector on the occupied space
	const vecfuncT& amo=get_calc().amo;
	vecfuncT act=this->active_mo();
	Projector<double,3> rho0(amo);

	const vecfuncT& xp=root.x;
	// apply the unperturbed potential and the Gamma potential
	// on the response amplitude x^p
	const vecfuncT V0=apply_fock_potential(xp);
	const vecfuncT Gamma=apply_gamma(xp,act,rho0);

	vecfuncT V_iroot=add(world,V0,Gamma);

	return V_iroot;

}

vecfuncT CIS::apply_gamma(const vecfuncT& x, const vecfuncT& act,
		const Projector<double,3>& rho0) const {
	print("apply_gamma");
	START_TIMER(world);
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
	if(exchange_ == "hf"){
		for (std::size_t p=0; p<noct; ++p) {

			// this is x_i * \int 1/r12 \phi_i \phi_p
			vecfuncT x_Ppi=mul(world,x,exchange_intermediate_[p]);
			for (std::size_t i=0; i<noct; ++i) Gamma[p]-=x_Ppi[i];

			// project out the zeroth-order density Eq. (4)
			Gamma[p]-=rho0(Gamma[p]);
		}
	}
	if(exchange_ == "dirac"){
		print("Use dirac exchange");
		double c = -3.0/4.0*pow(1.0/3.0,3.0/constants::pi);
		// Integral over rhoprime
		double dirac=0;
		for (size_t i=0; i<noct; ++i) dirac+=inner(x[i],act[i]);
		dirac *=c;
		for (std::size_t p=0; p<noct; ++p) {Gamma[p]+=dirac*act[p];}

		// project out the zeroth-order density Eq. (4)
		for (std::size_t p=0; p<noct; ++p) {
			Gamma[p]-=rho0(Gamma[p]);
		}
	}
	END_TIMER(world,"Apply Gamma");
	return Gamma;
}

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

/// Compute the expectation value from the perturbed fock operator
double CIS::expectation_value(World &world,const root &thisroot,const vecfuncT &Vx)const{

	std::vector<root> roots;
	roots.push_back(thisroot);
	std::vector<vecfuncT> V;
	V.push_back(Vx);

	orthonormalize_fock(world,roots,V);

	return roots[0].omega;

	/*
	START_TIMER(world);
	// <p | V | p> = \sum_i <x^p_i | V | x^q_i>
	const vecfuncT& x=thisroot.x;
	Tensor<double> tmp=inner(world,x,Vx);
	double omega=tmp.sum();
	END_TIMER(world,"Expectation value, potential");

	START_TIMER(world);
	// Kinetic part
	reconstruct(world, x);

	// The gradient operator
	std::vector< std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double,3>(world);

	// loop over all axes
	for(int axis = 0;axis < 3;++axis) {
		const vecfuncT dx = apply(world, *(gradop[axis]), x);
		Tensor<double> x_i_T_x_i=inner(world,dx,dx);
		omega+=0.5*x_i_T_x_i.sum();
		}
	END_TIMER(world,"Expectation value, kinetic");

	// add the term with the orbital energies:
	//  F_pq -= \sum_i\epsilon_i <x^p_i | x^q_i>
	vecfuncT act=this->active_mo();
	const std::size_t noct=act.size();

			Tensor<double> eij=inner(world,x,x);
			for (size_t ii=0; ii<noct; ++ii) {
				omega-=get_calc().aeps[ii+nfreeze_]*eij[ii];
			}



	return omega;//*/

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
Tensor<double> CIS::make_perturbed_fock_matrix(const std::vector<root>& roots, const std::vector<vecfuncT> &V) const{
	//const vecfuncT& act, const Projector<double,3>& rho0) const {

	const std::size_t nroot=roots.size();
	Tensor<double> Fock_pt(nroot,nroot);

	START_TIMER(world);
	// apply the unperturbed Fock operator and the Gamma potential on all
	// components of all roots
	for (std::size_t iroot=0; iroot<nroot; ++iroot) {

		//const vecfuncT& xp=roots[iroot].x;
		// apply the unperturbed potential and the Gamma potential
		// on the response amplitude x^p
		//const vecfuncT V0=apply_fock_potential(xp);
		//const vecfuncT Gamma=apply_gamma(xp,act,rho0);

		//const vecfuncT V=add(world,V0,Gamma);

		// compute the matrix element V_pq
		// <p | V | q> = \sum_i <x^p_i | V | x^q_i>
		for (std::size_t jroot=0; jroot<nroot; ++jroot) {
			const vecfuncT& xq=roots[jroot].x;
			Tensor<double> xpi_Vxqi=inner(world,xq,V[iroot]);
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

// To calculate the guess energy (Gamma must be created)
void CIS::orthonormalize_fock(World &world,std::vector<root> &roots)const{

	// construct the projector on the occupied space
	const vecfuncT& amo=get_calc().amo;
	vecfuncT act=this->active_mo();
	Projector<double,3> rho0(amo);

	std::vector<vecfuncT> Vphi;

	for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {

		const vecfuncT& xp=roots[iroot].x;
		// apply the unperturbed potential and the Gamma potential
		// on the response amplitude x^p
		const vecfuncT V0=apply_fock_potential(xp);
		const vecfuncT Gamma=apply_gamma(xp,act,rho0);

		const vecfuncT V=add(world,V0,Gamma);

		Vphi.push_back(V);
	}

	orthonormalize_fock(world,roots,Vphi);
}

/// orthonormalize all roots using the perturbed fock matrix (Energy update included)
void CIS::orthonormalize_fock(World &world,std::vector<root> &roots, const std::vector<vecfuncT> &Vphi)const{



	// construct the projector on the occupied space
	//const vecfuncT& amo=get_calc().amo;
	vecfuncT act=this->active_mo();
	const std::size_t noct=act.size();
	//Projector<double,3> rho0(amo);

	// compute the Fock matrix elements
	Tensor<double> Fock_pt=make_perturbed_fock_matrix(roots,Vphi);

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
		vc[iroot]=zero_functions_compressed<double,3>(world,noct);
		compress(world,roots[iroot].x);
	}


	for (size_t i=0; i<roots.size(); ++i) {
		for (size_t j=0; j<roots.size(); ++j) {
			gaxpy(world,1.0,vc[i],U(j,i),roots[j].x);
		}
	}
	// Suppress the update for the expectation value
	if(roots.size()>1){
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
	if(roots.size()==1){
		roots[0].omega=evals[0];
		if(world.rank()==0)print("Expectation value: ",roots[0].omega);
	}

}


/// orthonormalize all roots using Gram-Schmidt
void CIS::orthonormalize(World& world, std::vector<root>& roots) const {
	if(world.rank()==0)print("orthonormalize GS");

	// TEST THIS CAREFULLY
	sort_roots(roots,"error");
	if(world.rank()==0){
		std::cout<< "Control output: Sorting criterium is error, roots are: ";
		for(size_t i=0;i<roots.size();i++) std::cout<< roots[i].number<<" " ;
	}


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
	}

	// TEST THIS CAREFULLY ... sort back
	sort_roots(roots,"energy");

}

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


