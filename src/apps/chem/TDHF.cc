/*
 * TDHF.cc
 *
 *  Created on: Aug 11, 2016
 *      Author: kottmanj
 */

#include "TDHF.h"

namespace madness {

// KAIN allocator for vectorfunctions
struct TDHF_allocator{
	World& world;
	const int noct;

	/// @param[in]	world	the world
	/// @param[in]	nnoct	the number of functions in a given vector
	/// @todo validate doxygen on `nnoct`
	TDHF_allocator(World& world, const int nnoct) : world(world), noct(nnoct) {}

	vecfuncT operator()(){
		return zero_functions<double,3>(world,noct);
	}
	TDHF_allocator operator=(const TDHF_allocator &other){
		TDHF_allocator tmp(world,other.noct);
		return tmp;
	}

};





// conveniece to interface functions
TDHF::TDHF(World &world, const Nemo & nemo_, const Parameters& param):
								    																		  world(world),
																											  nemo(nemo_),
																											  parameters(param),
																											  g12(world,OT_G12,parameters.get_ccc_parameters()),
																											  mo_ket_(make_mo_ket(nemo_)),
																											  mo_bra_(make_mo_bra(nemo_)),
																											  Q(world,mo_bra_.get_vecfunction(),mo_ket_.get_vecfunction()),
																											  msg(world) {
	msg.section("TDHF initialized without the usual initialization routine (no intermediates, no parameters read)");
}

TDHF::TDHF(World &world, const Nemo & nemo_, const std::string& input):
						    																		  world(world),
																									  nemo(nemo_),
																									  parameters(nemo_.get_calc(),input),
																									  g12(world,OT_G12,parameters.get_ccc_parameters()),
																									  mo_ket_(make_mo_ket(nemo_)),
																									  mo_bra_(make_mo_bra(nemo_)),
																									  Q(world,mo_bra_.get_vecfunction(),mo_ket_.get_vecfunction()),
																									  msg(world) {
	msg.section("Initialize TDHF Class");

	msg.debug = parameters.debug;

	msg.subsection("General Information about settings from SCF object:\n");
	msg << " is_dft() = " << nemo.get_calc()->xc.is_dft() << "\n";
	msg << " hf_coeff = " << nemo.get_calc()->xc.hf_exchange_coefficient() << "\n";
	msg << " do_pcm() = " << nemo.do_pcm() << "\n";
	msg << " do_ac()  = " << nemo.do_ac() << "\n";

	parameters.print(world);
	const double old_thresh = FunctionDefaults<3>::get_thresh();
	if(old_thresh>parameters.thresh*0.1 and old_thresh>1.e-5){
		msg.warning("Threshold of Reference might be too loose |  Response thresh="+std::to_string(parameters.thresh)+ " and Reference thresh="+std::to_string(old_thresh) + ". Be careful, reference should be tight");
	}
	FunctionDefaults<3>::set_thresh(parameters.thresh);
	msg << "MRA Threshold is set to: " << FunctionDefaults<3>::get_thresh() << " with k=" << FunctionDefaults<3>::get_k() << "\n";

	if (not parameters.no_compute) {

		if(nemo.get_calc()->xc.hf_exchange_coefficient()!=0.0){
			msg.subsection("Computing Exchange Intermediate");
			CCTimer timer(world,"Computing ExIm");
			g12.update_elements(mo_bra_,mo_ket_);
			timer.info();
		}else msg.output("No Exchange Intermediate Computed\n");

		msg.output("Orbital Energies of Reference");
		const Tensor<double> eps=nemo.get_calc()->aeps;
		msg << eps << "\n";

	}
	if(nemo.get_calc()->param.localize){
		Fock F(world, &nemo);
		F_occ = F(get_active_mo_bra(),get_active_mo_ket());
		for(size_t i=0;i<get_active_mo_ket().size();++i){
			msg << std::scientific << std::setprecision(10);
			msg << "F(" << i << "," << i << ")=" << F_occ(i,i) << "\n";
			if(std::fabs(get_orbital_energy(i+parameters.freeze)-F_occ(i,i))>1.e-5){
				msg << "eps(" << i << ")=" << get_orbital_energy(i) << " | diff=" << get_orbital_energy(i+parameters.freeze)-F_occ(i,i) << "\n";
			}
		}
	}else{
		F_occ = Tensor<double>(get_active_mo_bra().size(),get_active_mo_ket().size());
		F_occ*=0.0;
		for(size_t i=0;i<get_active_mo_ket().size();++i){
			F_occ(i,i)=get_orbital_energy(i+parameters.freeze);
		}
	}
}

/// plot planes and cubes
void TDHF::plot(const vecfuncT& vf, const std::string& name)const{
	if(parameters.plot){
		CCTimer timer(world,"plot planes and cubes");
		madness::plot(vf,name,nemo.get_calc()->molecule.cubefile_header());
		timer.print();
	}
}

/// sort the xfunctions according to their excitation energy and name the excitation energies accordingly
std::vector<CC_vecfunction> TDHF::sort_xfunctions(std::vector<CC_vecfunction> x)const{
	std::sort(x.begin(),x.end());
	for(size_t i=0;i<x.size();++i){
		x[i].excitation=i;
	}
	return x;
}

/// print information
void TDHF::print_xfunctions(const std::vector<CC_vecfunction> & f, const bool& fullinfo)const{
	for(const auto& x:f){
		const double mem=get_size(world,x.get_vecfunction());

			msg << "ex. vector "
					<< x.excitation << " | "
					<< x.omega << " (ex. energy) | "
					<< x.current_error << " (error) | "
					<< x.delta << " (Edelta) | "
					<< mem << " (Gbyte)" << "\n";

		if(fullinfo){
			print_size(world,x.get_vecfunction(),"ex. vector "+std::to_string(x.excitation));
		}
	}
}

void TDHF::initialize(std::vector<CC_vecfunction> &start)const{



	msg.subsection("Calculate Guess");
	std::vector<CC_vecfunction> guess;
	bool use_old_guess=(parameters.generalkeyval.find("use_old_guess")!=parameters.generalkeyval.end());
	if(use_old_guess) use_old_guess=(parameters.generalkeyval.find("use_old_guess")->second=="true" or parameters.generalkeyval.find("use_old_guess")->second=="1");
	if(use_old_guess){
		guess=make_old_guess(get_active_mo_ket());
	}
	else{
		guess= make_guess_from_initial_diagonalization();
	}
	// combine guess and start vectors
	for(const auto& tmp:start) guess.push_back(tmp);
	std::vector<vecfuncT> empty;
	//orthonormalize(guess,empty);


	// failsafe (works in most cases)
	if(guess.size()<parameters.guess_excitations){
		std::string message=("WARNING: You demanded: " + std::to_string(parameters.guess_excitations)
		+ " Guess vectors, but your demanded guess has only "
		+ std::to_string(guess.size()) + "vectors. So we will not iterate the first vectors and then do the same guess again ... this might be unstable").c_str();
		msg.output(message);
		if(parameters.guess_maxiter==0){
			msg.output("In this case you demanded guess_maxiter=0 (which is also the default), so this can not work!");
			MADNESS_EXCEPTION("Faulty combinations of parameters given",1);
		}
		iterate_cis_guess_vectors(guess);
		initialize(guess);
	}

	//sort guess (according to excitation energies)
	std::sort(guess.begin(),guess.end());
	//truncate the guess
	std::vector<CC_vecfunction> guess_vectors;
	for(size_t i=0;i<parameters.guess_excitations;i++) guess_vectors.push_back(guess[i]);
	// this is the return value
	start = guess_vectors;
}

/// @param[in/out] CC_vecfunction
/// on input the guess functions (if empty or not enough the a guess will be generated)
/// on output the solution
std::vector<CC_vecfunction> TDHF::solve_cis()const{
	std::vector<CC_vecfunction> ccs;
	// look for restart options
	for(size_t k=0;k<parameters.restart.size();k++){
		CC_vecfunction tmp;
		const bool found= initialize_singles(tmp,RESPONSE,parameters.restart[k]);
		if(found) ccs.push_back(tmp);
	}

	for(size_t macrocycle=0;macrocycle<1;++macrocycle){
		//msg.section("CIS Macroiteration " + std::to_string(macrocycle));
		ccs=solve_cis(ccs);
		if(converged_roots.size()>=parameters.excitations) break;
	}

	return converged_roots;
}

std::vector<CC_vecfunction> TDHF::solve_cis(std::vector<CC_vecfunction> &start)const{
	msg.section("SOLVING CIS EQUATIONS");

	mo_ket_.plot("MOS_");

	CCTimer time(world,"TDHF/CIS");
	// decide if a guess calculation is needed
	bool need_guess =false;
	if(start.size()<parameters.guess_excitations) need_guess=true;
	std::vector<CC_vecfunction> guess_vectors;
	if(need_guess){
		initialize(start);
		guess_vectors = start;
	}else guess_vectors=start;

	msg.output("====Guess-Vectors=====");
	print_xfunctions(guess_vectors);

	std::vector<CC_vecfunction> final_vectors;
	msg.subsection("Iterate Guess Vectors");
	{
		// do guess iterations
		iterate_cis_guess_vectors(guess_vectors);
		// sort according to excitation energies
		std::sort(guess_vectors.begin(),guess_vectors.end());
		// save
		for(size_t i=0;i<guess_vectors.size();i++) guess_vectors[i].save_functions(std::to_string(i));
		// prepare final iterations
		for(size_t i=0;i<guess_vectors.size();i++){
			if(i<parameters.iterating_excitations) final_vectors.push_back(guess_vectors[i]);
			else guess_roots.push_back(guess_vectors[i]);
		}
		// sort guess_roots backwards in order to feed them into the cycle easily with pop_back
		std::sort(guess_roots.rbegin(),guess_roots.rend());

	}

	msg.subsection("Iterate Final Vectors");

	msg.output("Vectors in Iteration Cycle");
	print_xfunctions(final_vectors);
	msg.output("Remaining Guess Vectors");
	print_xfunctions(guess_roots);

	iterate_cis_final_vectors(final_vectors);
	msg.section("CIS CALCULATIONS ENDED");
	std::sort(converged_roots.begin(),converged_roots.end());
	for(size_t i=0;i<converged_roots.size();++i) converged_roots[i].excitation=i;
	std::sort(final_vectors.begin(),final_vectors.end());
	for(size_t i=0;i<final_vectors.size();++i) final_vectors[i].excitation=converged_roots.size()+i;
	// information
	msg.output("\n\nCONVERGED ROOTS\n");
	print_xfunctions(converged_roots);
	msg.output("\n\nRest\n");
	print_xfunctions(final_vectors);
	time.info();
	// plot final vectors
	for(size_t i=0;i<final_vectors.size();i++) final_vectors[i].plot(std::to_string(i)+"_converged_cis");
	// save final vectors
	for(size_t i=0;i<final_vectors.size();i++) final_vectors[i].save_functions(std::to_string(i));
	// assign solution vectors
	//start = final_vectors;\

	//give back as much as demanded
	std::vector<CC_vecfunction> result=converged_roots;
	for(const auto& x:final_vectors){
		result.push_back(x);
		if(result.size()==parameters.excitations) break;
	}
	result=sort_xfunctions(result);
	msg << "writing final functions to disc...\n";
	for(size_t i=0;i<result.size();i++) result[i].save_functions(std::to_string(i));

	// print out all warnings
	msg.print_warnings();
	return result;
}

void TDHF::solve_tdhf(std::vector<CC_vecfunction> &x)const{
	msg.section("SOLVING TDHF EQUATIONS");
	MADNESS_EXCEPTION("TDHF NOT IMPLEMENTED",1);
}

bool TDHF::iterate_cis_guess_vectors(std::vector<CC_vecfunction> &x)const{
	std::vector<CC_vecfunction> dummy;
	return iterate_vectors(x,dummy,false,parameters.guess_dconv,parameters.guess_econv,parameters.guess_maxiter, false);
}
bool TDHF::iterate_cis_final_vectors(std::vector<CC_vecfunction> &x)const{
	std::vector<CC_vecfunction> dummy;
	return iterate_vectors(x,dummy,false,parameters.dconv,parameters.econv,parameters.maxiter, parameters.kain_subspace>0);
}
bool TDHF::iterate_vectors(std::vector<CC_vecfunction> &x,const std::vector<CC_vecfunction> &y,
		const bool iterate_y,const double dconv, const double econv, const double iter_max, const bool kain)const{
	if(iterate_y) MADNESS_ASSERT(x.size()==y.size());

	msg.subsection("Iterating excitation vectors with the following parameters");
	msg		<< "dconv    = " << dconv << "\n"
			<< "econv    = " << econv << "\n"
			<< "iter_max = " << iter_max << "\n"
			<< "kain     = " << kain << "\n";

	// set up the kain solvers ... if needed or not
	std::vector<std::shared_ptr< XNonlinearSolver<vecfuncT,double,TDHF_allocator> > > solvers(x.size());
	// initialize solvers
	if(kain){
		for(size_t i=0;i<x.size();i++){
			solvers[i]=std::make_shared<XNonlinearSolver<vecfuncT,double,TDHF_allocator> >(TDHF_allocator(world,x[i].size()),true);
			solvers[i]->set_maxsub(parameters.kain_subspace);
		}
	}

	bool converged = true;

	// get potentials (if demanded)
	std::vector<vecfuncT> V;

	// for TDHF the potential is always stored
	if(parameters.store_potential and y.empty()) V = make_potentials(x);
	else if(not y.empty()) V = make_tdhf_potentials(x,y);


	for(size_t iter=0;iter<iter_max;iter++){
		CCTimer timer(world,"iteration "+std::to_string(iter));
		const bool this_is_a_guess_iteration=not (kain); // better readable code

		// if this is the first iteration, the potentials are not calculated yet
		if(iter==0) orthonormalize(x,V);

		// if we solve for the y functions we need to switch the sign of omega
		if(iterate_y){
			for(auto& tmp:x) tmp.omega *=-1.0;
		}
		// apply Greens Functions
		for(size_t i=0;i<x.size();i++) if(iterate_y and x[i].omega>=0.0) x[i].omega = -1.0*y[i].omega;
		std::vector<vecfuncT> residuals=apply_G(x,V);
		// check convergence

		converged = true;
		// error (vector norm)
		std::vector<double> errors;
		// largest individual error
		std::vector<std::pair<int,double> > largest_errors;
		for(size_t i=0;i<x.size();i++){
			double av_error = norm2(world,residuals[i]);
			std::vector<double> ind_err = norm2s(world,residuals[i]);
			auto it=max_element(ind_err.begin(),ind_err.end());
			const double largest_error=(*it);
			errors.push_back(av_error);
			largest_errors.push_back(std::make_pair(std::distance(ind_err.begin(),it),*it));
			// convergece criteria are the individual functions (orhterwise we have a dependence on the number of functions)
			if((*it)>dconv) converged=false;
			if(fabs(x[i].delta)>econv) converged=false;
			x[i].current_error=(largest_error);
		}
		// update, store old omegas in the deltas vector
		for(size_t i=0;i<x.size();i++){
			if(std::fabs(x[i].current_error)>dconv or std::fabs(x[i].delta)>econv or this_is_a_guess_iteration){
				vecfuncT new_x0;
				if(kain){
					new_x0 = solvers[i]->update(x[i].get_vecfunction(),residuals[i],10.0*parameters.thresh,5.0);
				}else{
					new_x0 = sub(world,x[i].get_vecfunction(),residuals[i]);
				}

				// Project out the converged roots
				if(converged_roots.size()>0){
					vecfuncT bra_new_x0=make_bra(new_x0);
					CCTimer timeP(world,"project out converged roots");
					for(const auto it:converged_roots){
						const double overlap = inner(world,it.get_vecfunction(),bra_new_x0).sum();
						new_x0 -= overlap*it.get_vecfunction();
					}
					timeP.print();
				}

				vecfuncT Q_new_x0 = Q(new_x0);
				truncate(world,Q_new_x0);
				x[i].set_functions(Q_new_x0,x[i].type,parameters.freeze);
			} else msg.output("Root " +std::to_string(i) + " converged");
		}

		// remove converged roots
		if(this_is_a_guess_iteration){
			// not in the guess iteration
		}else{
			std::vector<CC_vecfunction> unconverged_roots;
			std::vector<std::shared_ptr<XNonlinearSolver<vecfuncT,double,TDHF_allocator> > > corresponding_solvers;
			for(size_t i=0;i<x.size();++i){
				if(x[i].current_error<dconv and std::fabs(x[i].delta)<econv){
					// put converged roots into converged roots and replace it by one of the remaining guess functions
					converged_roots.push_back(x[i]);
					if(guess_roots.empty()) msg.output("Ran out of guess functions"); // no replacement
					else{
						// fill the unconverged roots with new guess roots
						unconverged_roots.push_back(guess_roots.back());
						// allocate a new solver for the new guess function
						auto sol=std::make_shared<XNonlinearSolver<vecfuncT,double,TDHF_allocator> >(TDHF_allocator(world,x[i].size()),true);
						sol->set_maxsub(parameters.kain_subspace);
						corresponding_solvers.push_back(sol);
						// delete the used guess root
						guess_roots.pop_back();
					}
				}else{
					unconverged_roots.push_back(x[i]);
					corresponding_solvers.push_back(solvers[i]);
				}
			}
			x=unconverged_roots;
			solvers=corresponding_solvers;
			sort_xfunctions(converged_roots);
		}

		// for TDHF the potential is always stored
		if(parameters.store_potential and y.empty()) V = make_potentials(x);
		else if(not y.empty()) V = make_tdhf_potentials(x,y);
		// orthonormalize
		orthonormalize(x,V);
		// save functions (every 5 iterations)
		if(iter%5==0){
			if(not iterate_y)for(size_t i=0;i<x.size();i++) x[i].save_functions(std::to_string(i));
			else for(size_t i=0;i<y.size();i++) y[i].save_functions(std::to_string(i));
		}
		// make plots if demanded
		if(parameters.plot){
			if(iterate_y) for(size_t i=0;i<y.size();i++) y[i].plot(std::to_string(i)+"_y_iter_"+std::to_string(iter)+"_");
			else for(size_t i=0;i<x.size();i++) x[i].plot(std::to_string(i)+"_iter_"+std::to_string(iter)+"_");
		}
		// information
		msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
		msg << "Iteration " << iter <<": omega, largest error, delta, size "<< "\n";

		msg.output("===========CONVERGED=ROOTS=====================");
		print_xfunctions(converged_roots);
		msg.output("=========NON-CONVERGED=ROOTS====================");
		print_xfunctions(x);

		msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');

		timer.stop().print();
		if(converged) break;
		if(converged_roots.size()==parameters.excitations) break;
		if(x.empty()) break;
	}
	return converged;
}



std::vector<vecfuncT>  TDHF::apply_G(std::vector<CC_vecfunction> &x,std::vector<vecfuncT> &V)const{

	std::string msg1 = "Applying Greens Function to vectors";
	if(V.empty()) msg1+=", with recalculated Potentials";
	CCTimer time(world,msg1);
	std::vector<vecfuncT> result;
	for(size_t i=0;i<x.size();i++){

		vecfuncT Vi;
		if(V.empty()) Vi=get_tda_potential(x[i]);
		else Vi=V[i];
		double omega = x[i].omega;
		if(x[i].type==RESPONSE) MADNESS_ASSERT(omega>0.0);
		else MADNESS_ASSERT(omega<0.0);
		if(x[i].type==UNDEFINED and V.empty()) msg.warning("Empty V but x is y state from TDHF");

		CCTimer time_N(world,"add nuclear potential");
		// the potentials still need the nuclear potential
		const Nuclear V(world,&nemo);
		vecfuncT VNi = V(x[i].get_vecfunction());
		Vi += VNi;
		time_N.info(parameters.debug);

		// scale potential
		scale(world,Vi,-2.0);

		// make bsh operators as pointers in order to apply them in parallel
		std::vector<std::shared_ptr<SeparatedConvolution<double,3> > > bsh(x[i].size());
		for (size_t p = 0; p < bsh.size(); p++) {
			double eps = get_orbital_energy(p+parameters.freeze) + omega;
			// if eps is above zero we have an unbound state (or are early in the iteration) however this needs a shift of the potential
			// we shift to -0.05 (same as moldft does, no particular reason)
			if(eps>0.0){
				msg.output("potential shift needed for V" + std::to_string(p+parameters.freeze));
				double shift = eps+0.05;
				eps = eps - shift;
				Vi[p] -= (-2.0*shift*x[i].get_vecfunction()[p]);
			}
			if(eps>0.0){
				msg.warning("eps is " + std::to_string(eps) + "... should not happen ... setting to zero");
				eps = -1.0*parameters.thresh;
			}
			MADNESS_ASSERT(not(eps>0.0));
			bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), parameters.lo, parameters.thresh_op));
		}
		world.gop.fence();

		vecfuncT GV = Q(apply(world, bsh, Vi));

		vecfuncT residual = sub(world,x[i].get_vecfunction(),GV);
		result.push_back(residual);

		// Calculate Second Order Energy Update
		const vecfuncT bra_GV=make_bra(GV);
		{
			// Inner product of Vpsi and the residual (Vi is scaled to -2.0 --> multiply later with 0.5)
			double tmp = inner(world,make_bra(residual),Vi).sum();
			// squared norm of GVpsi (Psi_tilde)
			double tmp2 = inner(world,make_bra(GV),GV).sum();

			// Factor 0.5 removes the factor 2 from the scaling before
			const double sou= (0.5 * tmp / tmp2);
			msg << "FYI: second order update would be: " << sou << " norm after QG is " << tmp2 << "\n";
		}
		// clear potential
		Vi.clear();
	}
	// clear the potentials
	V.clear();
	time.info();
	return result;
}

std::vector<vecfuncT> TDHF::make_potentials(const std::vector<CC_vecfunction> &x)const{
	CCTimer time(world,"Make Potentials");
	std::vector<vecfuncT> V;
	for(auto& xi:x){
		if(parameters.debug) msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
		const vecfuncT pot = get_tda_potential(xi);
		V.push_back(pot);
		if(parameters.debug) msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
	}
	time.info();
	MADNESS_ASSERT(V.size()==x.size());
	return V;
}

vecfuncT TDHF::get_tda_potential(const CC_vecfunction &x)const{
	// XC information
	const std::string xc_data = nemo.get_calc()->param.xc_data;
	// HF exchange Coefficient
	double hf_coeff = nemo.get_calc()->xc.hf_exchange_coefficient();

	// Use the PCMSolver
	bool pcm=nemo.do_pcm();
	if(parameters.debug){
		msg<< "TDA Potential is " << xc_data << ", hf_coeff=" << hf_coeff << ", pcm is=" << pcm <<  "\n";
	}
	if(hf_coeff<0.0) msg.warning("hf_exchange_coefficient is negative");

	// Occupation numbers
	const Tensor<double> occ=nemo.get_calc()->get_aocc();
	// Closed shell full density of the nemo orbitals (without any nuclear cusps)
	const real_function_3d nemo_density=2.0*nemo.make_density(occ,mo_ket_.get_vecfunction());
	// Real Alpha density (with nuclear cusps)
	const real_function_3d alpha_density=0.5*nemo.R_square*nemo_density;





	// Apply Ground State Potential to x-states
	vecfuncT Vpsi1;
	{
		// construct unperturbed operators
		const Coulomb J(world,&nemo);
		// const Nuclear V(world,&nemo); // not included in the TDA potential anymore


		std::string xc_data=nemo.get_calc()->param.xc_data;
		xc_data = xc_data.erase(0,xc_data.find_first_not_of(" "));
		xc_data = xc_data.erase(xc_data.find_last_not_of(" ")+1);

		// Nuclear Potential applied to x
		//CCTimer timeN(world,"Nx");
		//const vecfuncT Nx=V(x.get_vecfunction());
		//timeN.info(parameters.debug);
		// Applied Hartree Potential (J|x>) -> factor two is absorbed into the density for the J Operator
		CCTimer timeJ(world,"Jx");
		const vecfuncT Jx=J(x.get_vecfunction());
		timeJ.info(parameters.debug);

		if(nemo.get_calc()->xc.is_dft()){
			// XC Potential
			const XCOperator xc(world,xc_data, not nemo.get_calc()->param.spin_restricted,alpha_density,alpha_density);

			// Applied XC Potential
			CCTimer timeXCx(world,"XCx");
			real_function_3d xc_pot = xc.make_xc_potential();

			// compute the asymptotic correction of exchange-correlation potential
			if(nemo.do_ac()) {
				double charge = double(nemo.molecule().total_nuclear_charge());
				real_function_3d scaledJ = -1.0/charge*J.potential()*(1.0-hf_coeff);
				xc_pot = nemo.get_ac().apply(xc_pot, scaledJ);
			}

			const vecfuncT XCx=mul(world, xc_pot, x.get_vecfunction());
			// Ground State Potential applied to x, without exchange
			Vpsi1 = Jx+XCx; // Nx removed
		}else Vpsi1=Jx;
		// add exchange if demanded
		if(hf_coeff!=0.0){
			CCTimer timeKx(world,"Kx");
			Exchange K=Exchange(world,&nemo,0).small_memory(false);
			K.set_parameters(mo_bra_.get_vecfunction(),mo_ket_.get_vecfunction(),occ,parameters.lo,parameters.thresh_op);
			vecfuncT Kx =K(x.get_vecfunction());
			scale(world,Kx,hf_coeff);
			Vpsi1 = sub(world, Vpsi1, Kx);
			timeKx.info(parameters.debug);
		}
		// compute the solvent (PCM) contribution to the potential
		if (pcm) {
			CCTimer timepcm(world,"pcm:gs");
			const real_function_3d vpcm = nemo.get_pcm().compute_pcm_potential(J.potential(),false);
			if(parameters.plot or parameters.debug) plot_plane(world,vpcm,"vpcm_gs");
			const vecfuncT pcm_x=vpcm*x.get_vecfunction();
			timepcm.info(parameters.debug);
			Vpsi1 = add(world,Vpsi1,pcm_x);
		}
	}

	// Apply the Perturbed Potential to the Active Ground State Orbitals
	vecfuncT Vpsi2;
	{
		// active mo
		const vecfuncT active_mo = get_active_mo_ket();
		const vecfuncT active_bra = get_active_mo_bra();
		// construct perturbed operators
		CCTimer timeJ(world,"pXC");
		Coulomb Jp(world);
		real_function_3d density_pert=2.0*nemo.make_density(occ,active_bra,x.get_vecfunction());
		Jp.potential()=Jp.compute_potential(density_pert);

		vecfuncT XCp=zero_functions<double,3>(world,get_active_mo_ket().size());
		if(nemo.get_calc()->xc.is_dft()){
			// XC Potential
			const XCOperator xc(world,xc_data, not nemo.get_calc()->param.spin_restricted,alpha_density,alpha_density);
			// reconstruct the full perturbed density: do not truncate!
			real_function_3d gamma=xc.apply_xc_kernel(density_pert);
			vecfuncT XCp=mul(world,gamma,active_mo);
			truncate(world,XCp);
		}

		if(parameters.triplet){
			if(norm2(world,XCp)!=0.0) MADNESS_EXCEPTION("Triplets only for CIS",1);
			Vpsi2 = XCp;
		}
		else Vpsi2 = Jp(active_mo)+XCp;
		timeJ.info(parameters.debug);
		// Exchange Part
		if(hf_coeff>0.0){
			CCTimer timeK(world,"pK");
			vecfuncT Kp;
			// summation over all active indices
			for(const auto itmp:x.functions){
				const size_t i=itmp.first;
				real_function_3d Ki=real_factory_3d(world);
				for(const auto ktmp:x.functions){
					const size_t k=ktmp.first;
					Ki+=(g12(mo_bra_(k),mo_ket_(i))*x(k).function).truncate();
				}
				Kp.push_back(Ki);
			}
			scale(world,Kp,hf_coeff);
			Vpsi2 = sub(world,Vpsi2,Kp);
			timeK.info(parameters.debug);
			truncate(world,Vpsi2);
		}

		// compute the solvent (PCM) contribution to the kernel
		if (pcm) {
			CCTimer timepcm(world,"pcm:ex");
			const real_function_3d vpcm = nemo.get_pcm().compute_pcm_potential(Jp.potential(),true);
			if(parameters.plot or parameters.debug) plot_plane(world,vpcm,"vpcm_ex");
			const vecfuncT pcm_orbitals=vpcm*active_mo;
			timepcm.info(parameters.debug);
			Vpsi2 = add(world,Vpsi2,pcm_orbitals);
		}

		truncate(world,Vpsi2);
	}
	// whole tda potential
	vecfuncT Vpsi = Vpsi1 + Q(Vpsi2);
	// if the ground state is localized add the coupling terms
	// canonical: -ei|xi> (part of greens function)
	// local:  -fik|xk> (fii|xi> part of greens functions, rest needs to be added)

	if(nemo.get_calc()->param.localize){
		const vecfuncT vx=x.get_vecfunction();
		vecfuncT fock_coupling=madness::transform(world,vx,F_occ);
		// subtract the diagonal terms
		for(size_t i=0;i<fock_coupling.size();++i){
			fock_coupling[i]=(fock_coupling[i]-F_occ(i,i)*vx[i]);
		}
		Vpsi-=fock_coupling;
	}


	truncate(world,Vpsi);



	// debug output
	if(parameters.debug or parameters.plot){
		plot_plane(world,Vpsi,"Vpsi");
		plot_plane(world,Vpsi1,"Vpsi1");
		plot_plane(world,Vpsi2,"Vpsi2");
	}

	return Vpsi;

}

std::vector<vecfuncT> TDHF::make_tdhf_potentials(std::vector<CC_vecfunction> &x,const std::vector<CC_vecfunction> &y)const{
	MADNESS_EXCEPTION("NOT IMPLEMENTED",1);
}



void TDHF::orthonormalize(std::vector<CC_vecfunction> &x,std::vector<vecfuncT> &V)const{
	if(x.empty()) return;
	CCTimer time(world,"Orthonormalization");

	// make the overlap matrix
	Tensor<double> S = make_overlap_matrix(x);
	if(parameters.debug) msg << "The Overlap Matrix\n " << S << "\n";

	//make the Hamilton matrix for the vectorfunctions
	Tensor<double> F = make_perturbed_fock_matrix(x,V);

	// Diagonalize the F Matrix
	Tensor<double> U, evals;
	Tensor<double> dummy(x.size());
	U = nemo.get_calc() -> get_fock_transformation(world, S, F, evals, dummy, 2.0*parameters.thresh);

	if(parameters.debug){
		msg << "Perturbed Fock-Matrix Eigenvalues\n";
		for(int x=0;x<evals.size();++x) msg << evals(x) << "\n";
	}

	// Transform the states
	x = transform(x,U);

	// Transform the potentials (if V is empty nothing will happen)
	V = transform(V,U);

	// assign new energies and  get energy differences
	for(size_t i=0;i<x.size();i++){
		const double old_omega = x[i].omega;
		const double new_omega = evals(i);
		const double delta = new_omega-old_omega;
		x[i].omega=new_omega;
		x[i].delta=delta;
	}

	time.info();
}

std::vector<CC_vecfunction> TDHF::transform(const std::vector<CC_vecfunction> &x,const madness::Tensor<double> U) const {
	std::vector<CC_vecfunction> transformed;
	for(size_t k=0;k<x.size();k++){
		vecfuncT new_x = zero_functions_compressed<double,3>(world,x[k].size());
		compress(world,x[k].get_vecfunction());
		for(size_t l=0;l<x.size();l++){
			gaxpy(world,1.0,new_x,U(l,k),x[l].get_vecfunction()); // gaxpy(alpha,a,beta,b) -> a[i]=alpha*a[i] + beta*b[i], since there is no += for vectorfunctions implemented
		}
		CC_vecfunction tmp(x[k]);
		tmp.set_functions(new_x,tmp.type,parameters.freeze);
		transformed.push_back(tmp);
	}
	MADNESS_ASSERT(transformed.size()==x.size());
	return transformed;

}

Tensor<double> TDHF::make_overlap_matrix(const std::vector<CC_vecfunction> &x)const{
	CCTimer time(world,"Make Overlap Matrix");
	Tensor<double> S(x.size(),x.size());
	for(size_t k=0;k<x.size();k++){
		const vecfuncT kbra = make_bra(x[k]);
		for(size_t l=0;l<x.size();l++){
			S(l,k) = inner(world,kbra,x[l].get_vecfunction()).sum();
		}
	}
	time.info();
	if(parameters.debug)msg << std::fixed << std::setprecision(5) << "\nOverlap Matrix\n" << S << "\n";
	return S;
}

Tensor<double> TDHF::make_perturbed_fock_matrix(const std::vector<CC_vecfunction> &x, const std::vector<vecfuncT> &V)const{
	// Make formated timings
	CCTimer timeF(world,"Matrix: F");
	CCTimer timeT(world,"Matrix: T+Vn");
	CCTimer timeV(world,"Matrix: V");
	CCTimer timeR(world,"Matrix: e");

	// bra elements of x
	std::vector<vecfuncT> xbra;

	{
		CCTimer time_bra(world,"Make bra elements");
		for(size_t k=0;k<x.size();k++){
			const vecfuncT xbrak = make_bra(x[k]);
			xbra.push_back(xbrak);
		}
		MADNESS_ASSERT(xbra.size()==x.size());
		time_bra.info(parameters.debug);
	}

	timeF.start();
	Tensor<double> F(x.size(),x.size());
	{
		Tensor<double> T(x.size(),x.size());
		{
			timeT.start();
			// gradient operator
			std::vector < std::shared_ptr<real_derivative_3d> > D =gradient_operator<double, 3>(world);

			real_function_3d Vnuc = nemo.get_calc()->potentialmanager->vnuclear();
			if(not Vnuc.is_initialized()){
				msg.output("Compute Nuclear Potential");
				nemo.get_calc()->potentialmanager->make_nuclear_potential(world);
				Vnuc = nemo.get_calc()->potentialmanager->vnuclear();
			}

			const real_function_3d R = nemo.nuclear_correlation -> function();
			std::vector<vecfuncT> Rx(x.size(),zero_functions<double,3>(world,x.front().size()));
			CCTimer timeR(world,"make Rx");
			for(size_t k=0;k<x.size();k++){
				Rx[k] = mul(world,R,x[k].get_vecfunction(),false);
			}
			world.gop.fence();
			timeR.info(parameters.debug);
			std::vector<vecfuncT> dx(x.size(),zero_functions<double,3>(world,x.front().size()));
			std::vector<vecfuncT> dy(x.size(),zero_functions<double,3>(world,x.front().size()));
			std::vector<vecfuncT> dz(x.size(),zero_functions<double,3>(world,x.front().size()));
			CCTimer timeD(world,"make Grad(Rx)");
			for(size_t k=0;k<x.size();k++){
				dx[k] = apply(world,*(D[0]),Rx[k],false);
				dy[k] = apply(world,*(D[1]),Rx[k],false);
				dz[k] = apply(world,*(D[2]),Rx[k],false);
			}
			world.gop.fence();
			timeD.info(parameters.debug);

			CCTimer time_mat(world,"T+V Mat");
			for(size_t k=0;k<x.size();k++){
				const vecfuncT Vxk = mul(world,Vnuc,x[k].get_vecfunction());
				for(size_t l=0;l<x.size();l++){
					T(l,k) = inner(world,xbra[l],Vxk).sum();
					T(l,k) += 0.5*inner(world,dx[l],dx[k]).sum();
					T(l,k) += 0.5*inner(world,dy[l],dy[k]).sum();
					T(l,k) += 0.5*inner(world,dz[l],dz[k]).sum();
				}

			}
			time_mat.info(parameters.debug);
			timeT.stop();
		}
		Tensor<double> MV(x.size(),x.size());
		{
			timeV.start();
			bool recompute_V = V.empty();
			for(size_t k=0;k<x.size();k++){
				vecfuncT Vk;
				//if(recompute_V) Vk = CCOPS.get_CIS_potential(x[k]);
				if(recompute_V){
					msg.output("Recompute V");
					Vk = get_tda_potential(x[k]);
				}
				else Vk = V[k];
				for(size_t l=0;l<x.size();l++){
					MV(l,k) = inner(world,xbra[l],Vk).sum();
				}
			}
			timeV.stop();
		}

		// now set the fock matrix together: F(l,k) = T(l,k) + MV(l,k) - eps(l,k)
		// with eps(l,k) = \sum_i eps_i <xl_i|xk_i>
		// first get active eps
		timeR.start();
		std::vector<double> eps;
		for(size_t i=0;i<mo_ket_.size();i++) eps.push_back(get_orbital_energy(i));
		std::vector<double> active_eps;
		for(size_t i=parameters.freeze;i<eps.size();i++) active_eps.push_back(eps[i]);
		for(size_t k=0;k<x.size();k++){
			for(size_t l=0;l<x.size();l++){
				Tensor<double> xlk = inner(world,xbra[l],x[k].get_vecfunction());
				MADNESS_ASSERT(xlk.size()==active_eps.size());
				double eps_part =0.0;
				for(size_t i=0;i<active_eps.size();i++) eps_part += xlk(i)*active_eps[i];
				F(l,k) = T(l,k) + MV(l,k) - eps_part;
			}
		}
		timeR.stop();
		if(parameters.debug)msg<< std::fixed << std::setprecision(5) << "\n(T+V) Matrix\n" << T << "\n";
		if(parameters.debug)msg<< std::fixed << std::setprecision(5) << "\nPotential Matrix\n" << MV << "\n";
		if(parameters.debug)msg<< std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
	}
	timeF.stop();
	if(parameters.debug)msg << std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
	//formated timings output
	timeT.print();
	timeV.print();
	timeR.print();
	timeF.print();

	// symmetryze
	F = 0.5*(F + transpose<double>(F));
	if(parameters.debug) msg << std::fixed << std::setprecision(5) << "\nSymmetrized Perturbed Fock Matrix\n" << F << "\n";

	return F;
}

/// Makes the (old) guess functions by exciting active orbitals with excitation operators
std::vector<CC_vecfunction> TDHF::make_old_guess(const vecfuncT& f)const{
	CCTimer time(world,"Making Guess Functions: " + parameters.guess_virtuals);
	std::vector<std::string> exop_strings;
	if(parameters.guess_virtuals=="custom"){
		exop_strings = parameters.exops;
		msg << "Custom Excitation Operators Demanded:\n";
      	msg << exop_strings << "\n";
	}
	else exop_strings = make_predefined_exop_strings(parameters.guess_virtuals);

	// make the excitation operators
	vecfuncT exops;
	for(const auto& exs:exop_strings){
		std::shared_ptr<FunctionFunctorInterface<double, 3> > exop_functor(new polynomial_functor(exs));
		real_function_3d exop = real_factory_3d(world).functor(exop_functor);
		// do damp
		if(parameters.damping_width > 0.0){
			std::shared_ptr<FunctionFunctorInterface<double, 3> > damp_functor(new gauss_functor(parameters.damping_width));
			real_function_3d damp = real_factory_3d(world).functor(damp_functor);
			plot_plane(world,damp,"damping_function");
			exop = (exop*damp).truncate();
		}
		exops.push_back(exop);
	}

	// Excite the last N unfrozen MOs

	size_t N = std::min(parameters.guess_occ_to_virt,int(f.size()));
	// if N was not assigned we use all orbitals
	if(N==0){
		N=f.size();
	}

	// making the guess
	std::vector<CC_vecfunction> guess;
	for(size_t i=0;i<exops.size();i++){
		const vecfuncT& vm = f;
		reconstruct(world,vm);
		reconstruct(world,exops);
		MADNESS_ASSERT(not(N>vm.size()));
		vecfuncT tmp= zero_functions<double,3>(world,vm.size());
		// exciting the first N orbitals (from the homo to the homo-N)
		for(size_t k=0;k<N;k++){
			real_function_3d xmo = (exops[i]*vm[vm.size()-1-k]).truncate();
			tmp[tmp.size()-1-k]=xmo;
			plot_plane(world,xmo,std::to_string(i)+"_cis_guess_"+"_"+std::to_string(vm.size()-k-1+parameters.freeze));
		}
		{
			const double norm = sqrt(inner(world,make_bra(tmp),tmp).sum());
			scale(world,tmp,1.0/norm);
		}
		tmp = Q(tmp);
		{
			const double norm = sqrt(inner(world,make_bra(tmp),tmp).sum());
			scale(world,tmp,1.0/norm);
		}
		CC_vecfunction guess_tmp(tmp,RESPONSE,parameters.freeze);
		guess.push_back(guess_tmp);
	}

	time.info();
	return guess;
}

vecfuncT TDHF::make_virtuals() const {
	CCTimer time(world, "make virtuals");
	// create virtuals
	vecfuncT virtuals;
	if (parameters.guess_virtuals == "external") {
		madness::load_function(world, virtuals, "mybasis");
		//virtuals=Q(virtuals);
		for (auto& x : virtuals) {
			const double norm = sqrt(inner(make_bra(x), x));
			x.scale(1.0 / norm);
		}
	} else if (parameters.guess_virtuals == "scf") {
		// use the ao basis set from the scf calculations as virtuals (like projected aos)
		virtuals = (nemo.get_calc()->ao);
		for (auto& x : virtuals) {
			const double norm = sqrt(inner(make_bra(x), x));
			x.scale(1.0 / norm);
		}
	} else{
		// create the seeds
		vecfuncT xmo;
		for(size_t i=0;i<parameters.guess_occ_to_virt;++i) xmo.push_back(get_active_mo_ket()[get_active_mo_ket().size()-1-i]);

		bool use_trigo=true;
		if(parameters.generalkeyval.find("polynomial_exops")!=parameters.generalkeyval.end()) use_trigo = (std::stoi(parameters.generalkeyval.find("polynomial_exops")->second)==0);
		virtuals = apply_excitation_operators(xmo,use_trigo);

	}

	if(parameters.guess_cm>0.0){
		// add center of mass diffuse functions
		const double factor=parameters.guess_cm;
		const double width = (-1.0*factor/(get_orbital_energy(mo_ket_.size()-1)));
		msg.subsection("adding center of mass functions with exponent homo/c and c=" + std::to_string(factor));
		msg.output("width="+std::to_string(width));

		Tensor<double> cm = nemo.get_calc()->molecule.center_of_mass();
		msg << "center of mass is " << cm << "\n";
		polynomial_functor px("x 1.0",width,cm);
		polynomial_functor py("y 1.0",width,cm);
		polynomial_functor pz("z 1.0",width,cm);
		gauss_functor s(width,cm);
		real_function_3d vpx=real_factory_3d(world).functor(px);
		real_function_3d vpy=real_factory_3d(world).functor(py);
		real_function_3d vpz=real_factory_3d(world).functor(pz);
		real_function_3d vs=real_factory_3d(world).functor(s);
		virtuals.push_back(vpx);
		virtuals.push_back(vpy);
		virtuals.push_back(vpz);
		virtuals.push_back(vs);
	}
	if (world.rank() == 0)
		msg << virtuals.size() << " virtuals\n";

	plot(virtuals,"virtuals");
	CCTimer timerQ(world,"virt=Qvirt");
	virtuals=Q(virtuals);
	timerQ.print();
	plot(virtuals,"Qvirtuals");

	if (parameters.debug) {
		for (const auto& x : virtuals)
			x.print_size("virtual");
	}
	world.gop.fence();
	time.print();
	return virtuals;
}

vecfuncT TDHF::apply_excitation_operators(const vecfuncT& seed, const bool& use_trigo) const {
	//const int nvirt = seed.size() * ((2 * order * 2 * order * 2 * order) - 1);
	//	msg.subsection("creating a set of " + std::to_string(nvirt) + " virtuals by multiplying functions with plane waves");
	// compute the centers of the seed functions
	CCTimer time_centers(world,"compute centers");
	std::vector<coord_3d> centers = compute_centroids(seed);
	time_centers.print();


	// prepare the list of excitation operators and copied seeds
	CCTimer time_init_exop(world,"initialize excitation operators");
	std::vector<std::pair<vecfuncT, std::string> > exlist;
	{
		std::vector<std::string> exop_strings=parameters.exops;
		if(parameters.guess_virtuals!="custom") exop_strings=(make_predefined_exop_strings(parameters.guess_virtuals));
		for(const auto ex: exop_strings){
			vecfuncT cseed=copy(world,seed,false);
			exlist.push_back(std::make_pair(cseed,ex));
		}
	}
	world.gop.fence();
	time_init_exop.print();
	msg << "will create " << exlist.size()*seed.size() << " virtuals, from " << seed.size() << " seeds and " << exlist.size() << " excitation operators"   <<" \n";

	// create the virtuals by unary operations: multiply excitation operators with seeds
	CCTimer time_create_virtuals(world,"create virtuals");
	vecfuncT virtuals;
	for(auto it:exlist){
		if(use_trigo) virtuals=append(virtuals,apply_trigonometric_exop(it.first,it.second,centers,false));
		else virtuals=append(virtuals,apply_polynomial_exop(it.first,it.second,centers,false));
	}
	world.gop.fence();
	time_create_virtuals.print();

	return virtuals;

}

/// make the initial guess by explicitly diagonalizing a CIS matrix with virtuals from the make_virtuals routine
vector<CC_vecfunction> TDHF::make_guess_from_initial_diagonalization() const {
	CCTimer time(world, "make_guess_from_initial_diagonalization");
	//convenience
	const int nact = get_active_mo_ket().size();
	// create virtuals
	vecfuncT virtuals = make_virtuals();
	// canonicalize virtuals
	virtuals = canonicalize(virtuals);
	// compute the CIS matrix
	Tensor<double> MCIS = make_cis_matrix(virtuals);
	// initialize the guess functions
	if (world.rank() == 0 && MCIS.dim(0) < parameters.guess_excitations) {
		msg.warning(std::to_string(parameters.guess_excitations) + " guess vectors where demanded, but with the given options only " + std::to_string(MCIS.dim(0)) + " can be created\n");
	}
	std::vector<CC_vecfunction> xfunctions;
	for (int x = 0; x < MCIS.dim(0); ++x) {
		if (x >= parameters.guess_excitations)
			break;

		CC_vecfunction init(zero_functions<double, 3>(world, nact), RESPONSE, parameters.freeze);
		xfunctions.push_back(init);
	}
	{
		Tensor<double> U, evals;
		CCTimer time_diag(world, "cis-matrix diagonalization");
		syev(MCIS, U, evals);
		time_diag.print();

			msg << "Initial Diagonalization of CIS Matrix:\n Lowest Eigenvalues are \n";
			for (int x = 0; x < std::max(std::min(int(evals.size()), 10), parameters.guess_excitations); ++x)
				msg << evals(x) << "\n";

		CCTimer time_assemble(world, "assemble guess vectors");
		// make x functions from amplitudes and virtuals
		// todo: probably not optimal for highly parallel applications
		const int nvirt = virtuals.size();
		auto get_com_idx = [nvirt](int i, int a) {
			return i * nvirt + a;
		};
		auto get_vir_idx = [nvirt](int I) {
			return I % nvirt;
		};
		auto get_occ_idx = [nvirt](int I) {
			return I / nvirt;
		};
		for (size_t I = 0; I < MCIS.dim(0); ++I) {
			if (I >= parameters.guess_excitations)
				break;

			const int a = get_vir_idx(I);
			const int i = get_occ_idx(I);
			if (evals(I) < 0.0 && world.rank() == 0)
				msg.warning("NEGATIVE EIGENVALUE IN INITIAL DIAGONALIZATION: CHECK YOUR REFERENCE!\n");

			if (evals(I) < 1.e-5) {

					msg << "skipping root " << evals(I) << " \n";

				continue;
			}
			for (size_t J = 0; J < MCIS.dim(1); ++J) {
				const int b = get_vir_idx(J);
				const int j = get_occ_idx(J);
				const double xjb = U(J, I);
				xfunctions[I].get_vecfunction()[j] += xjb * virtuals[b];
				xfunctions[I].omega = evals(I);
				xfunctions[I].excitation = I;
			}
		}
		time_assemble.print();
		CCTimer time_truncate(world, "truncate guess");
		for (auto& x : xfunctions) {
			vecfuncT tmp = x.get_vecfunction();
			truncate(world, tmp, parameters.thresh);
			// should be truncated by shallow copy anyways ... but just to be sure
			x.set_functions(tmp, x.type, parameters.freeze);
		}
		time_truncate.print();
	}
	if (parameters.debug) {
		Tensor<double> S = make_overlap_matrix(xfunctions);

			msg << "Overlap matrix of guess:\n" << S << "\n";
	}
	print_xfunctions(xfunctions, true);

		msg << "created " << xfunctions.size() << " guess vectors\n";

	time.print();
	return xfunctions;
}
/// canonicalize a set of orbitals (here the virtuals for the guess)
vecfuncT TDHF::canonicalize(const vecfuncT& v)const{
	CCTimer time(world,"canonicalize");
	Fock F(world, &nemo);
	const vecfuncT vbra=make_bra(v);
	Tensor<double> Fmat = F(vbra,v);
	Tensor<double> S = matrix_inner(world, vbra, v);
	Tensor<double> occ(v.size());
	occ=1.0;
	Tensor<double> evals;
	if(parameters.debug) msg << "Canonicalize: Fock Matrix\n" << Fmat(Slice(0,std::min(10,int(v.size()))-1),Slice(0,std::min(10,int(v.size()))-1));
	if(parameters.debug) msg << "Canonicalize: Overlap Matrix\n" << S(Slice(0,std::min(10,int(v.size()))-1),Slice(0,std::min(10,int(v.size()))-1));
	Tensor<double> U = nemo.get_calc()->get_fock_transformation(world, S, Fmat, evals, occ, std::min(parameters.thresh,1.e-4));
	vecfuncT result = madness::transform(world, v, U);
	time.print();
	return result;
}
/// compute the CIS matrix for a given set of virtuals
Tensor<double> TDHF::make_cis_matrix(const vecfuncT virtuals)const{

	// make bra elements
	const vecfuncT virtuals_bra = make_bra(virtuals);
	// make Fock Matrix of virtuals for diagonal elements
	Fock F(world, &nemo);
	Tensor<double> Fmat = F(virtuals_bra, virtuals);


	if (parameters.debug) {
		const int dim = std::min(10,int(virtuals.size()));
			msg << "Debug Part of Virtual Fock Matrix\n" << Fmat(Slice(0,dim-1),Slice(0,dim-1)) << "\n";

		Tensor<double> S = matrix_inner(world, virtuals_bra, virtuals);
			msg << "Debug Overlap of virtuals\n" << S(Slice(0,dim-1),Slice(0,dim-1)) << "\n";
	}


	CCTimer time_cis(world, "make CIS matrix");

	// the cis matrix is indexed by ij and ab
	// we will use the combined indixes from ia and jb named I and J
	// in order to not be confused we use the following helper functions
	const int nocc = get_active_mo_ket().size();
	// determines for which orbitals (couting from the HOMO downwards) the off-diagonal elements will be computed
	// this simplifies the guess
	int active_guess_orbitals = parameters.guess_active_orbitals;
	const int nvirt = virtuals.size();
	auto get_com_idx = [nvirt](int i, int a) { return i*nvirt+a; };
	auto get_vir_idx = [nvirt](int I) {return I%nvirt;};
	auto get_occ_idx = [nvirt](int I) {return I/nvirt;};
	auto delta = [](int x, int y) {if (x==y) return 1; else return 0;};


	const int dim=(virtuals.size()*nocc);
	msg << "CIS-Matrix for guess calculation will be of size " << dim << "x" << dim << "\n";
	// the number of the matrix where elements which are not determined by orbital energies and the fock matrix are computed (controlled over active_guess_orbitals parameter)
	const int dim2=(virtuals.size()*active_guess_orbitals);
	if(dim2<dim) msg << "Effective size through neglect of some orbitals will be: " << dim2 << "x" << dim2 << "\n";
	const int start_ij = nocc-active_guess_orbitals;
	Tensor<double> MCIS(dim,dim);

	// make CIS matrix
	// first do the "diagonal" entries
	if(nemo.get_calc()->param.localize){
		Tensor<double> Focc = F(get_active_mo_bra(),get_active_mo_ket());
		for(int I=0;I<dim;++I){
			const int a=get_vir_idx(I);
			const int i=get_occ_idx(I);
			for(int J=0;J<dim;++J){
				const int b=get_vir_idx(J);
				const int j=get_occ_idx(J);
				MCIS(I,J) = Fmat(a,b)*delta(i,j)-Focc(i,j)*delta(a,b);
			}
		}
	}else{
		for(int I=0;I<dim;++I){
			const int a=get_vir_idx(I);
			const int i=get_occ_idx(I);
			MCIS(I,I) = Fmat(a,a)-get_orbital_energy(i+parameters.freeze);
		}
	}

	if(not parameters.guess_diag){
		int I = -1; // combined index from i and a, start is -1 so that initial value is 0 (not so important anymore since I dont use ++I)
		for (int i = start_ij; i < get_active_mo_ket().size(); ++i) {
			const real_function_3d brai = get_active_mo_bra()[i];
			const vecfuncT igv = g12(brai * virtuals);
			for (int a = 0; a < virtuals.size(); ++a) {
				I=get_com_idx(i,a);
				int J =-1;
				for (int j = start_ij; j < get_active_mo_ket().size(); ++j) {
					const real_function_3d braj =get_active_mo_bra()[j];
					for (int b = 0; b < virtuals.size(); ++b) {
						J=get_com_idx(j,b);
						if(J<=I){
							const real_function_3d igj = g12(mo_bra_(i+parameters.freeze),mo_ket_(j+parameters.freeze)); // use exchange intermediate
							const double rIJ = 2.0 * inner(braj * virtuals[b], igv[a]) - inner(virtuals_bra[a] * virtuals[b],igj);
							MCIS(J,I) += rIJ;
							MCIS(I,J) += rIJ;
						}
					}
				}
			}
		}
	}
		int sdim=std::min(int(MCIS.dim(0)),10);
		msg << "Part of the CIS Matrix:\n" << MCIS(Slice(dim-sdim,-1),Slice(dim-sdim,-1)) << "\n";
		if(parameters.debug) msg << "Debug: Full CIS Matrix:\n" << MCIS<< "\n";


	// test if symmetric
	if (parameters.debug) {
		const double symm_norm = (MCIS - transpose(MCIS)).normf();
			msg << "Hermiticity of CIS Matrix:\n" << "||MCIS-transpose(MCIS)||=" << symm_norm << "\n";

		if (symm_norm > 1.e-4) {
			int sliced_dim = 8;
			if (8 > MCIS.dim(0))
				sliced_dim = MCIS.dim(0);

				msg << "first " << sliced_dim << "x" << sliced_dim << " block of MCIS Matrix\n" << MCIS(_, Slice(sliced_dim - 1, sliced_dim - 1));
		}
	}
	time_cis.info();

	return MCIS;
}

bool TDHF::initialize_singles(CC_vecfunction &singles,const FuncType type,const int ex) const {
	MADNESS_ASSERT(singles.size()==0);
	bool restarted = false;

	std::vector<CCFunction> vs;
	for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
		CCFunction single_i;
		single_i.type=type;
		single_i.i = i;
		std::string name;
		if(ex<0) name = single_i.name();
		else name = std::to_string(ex)+"_"+single_i.name();
		real_function_3d tmpi = real_factory_3d(world);
		print("trying to load function ",name);
		const bool found=load_function<double,3>(tmpi,name);
		if(found) restarted = true;
		else msg.output("Initialized " + single_i.name()+" of type " + assign_name(type) +" as zero-function");
		single_i.function = copy(tmpi);
		vs.push_back(single_i);
	}

	singles = CC_vecfunction(vs,type);
	if(type==RESPONSE) singles.excitation=ex;

	return restarted;
}





/// compute the oscillator strength in the length representation

/// the oscillator strength is given by
/// \f[
/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]  root    a converged root
double TDHF::oscillator_strength_length(const CC_vecfunction& x) const {
	Tensor<double> mu_if(3);
	for (int idim=0; idim<3; idim++) {
		real_function_3d ri = real_factory_3d(world).functor(xyz(idim));
		vecfuncT amo_times_x=ri*get_active_mo_bra();
		Tensor<double> a=inner(world,amo_times_x,x.get_vecfunction());
		mu_if(idim)=a.sum();
	}
	const double f= 2.0/3.0 * x.omega * mu_if.sumsq() * 2.0;
	return f;
}

/// compute the oscillator strength in the velocity representation

/// the oscillator strength is given by
/// \f[
/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]  root    a converged root
double TDHF::oscillator_strength_velocity(const CC_vecfunction& x) const {
	Tensor<double> p_if(3);
	// compute the derivatives of the MOs in all 3 directions
	const vecfuncT Rroot=nemo.R*x.get_vecfunction();
	const vecfuncT Rnemo=nemo.R*get_active_mo_ket();

	for (int idim=0; idim<3; idim++) {
		real_derivative_3d D = free_space_derivative<double,3>(world, idim);
		vecfuncT Damo=apply(world,D,Rnemo);
		Tensor<double> a=inner(world,Damo,Rroot);
		p_if(idim)=a.sum();
	}
	const double f= 2.0/(3.0 * x.omega) * p_if.sumsq() * 2.0;
	return f;
}


/// analyze the root: oscillator strength and contributions from occ
void TDHF::analyze(const std::vector<CC_vecfunction> &x) const {

	const size_t noct=get_active_mo_ket().size();

	for (const CC_vecfunction& root : x) {

		const vecfuncT Rroot=nemo.R*root.get_vecfunction(); // reintroduce the nuclear correlation factor
		std::vector<double> norms=norm2s(world,Rroot);

		// compute the oscillator strengths and dominant contributions
		double osl=this->oscillator_strength_length(root);
		double osv=this->oscillator_strength_velocity(root);

		msg << std::scientific << std::setprecision(10) << std::setw(20);

			msg<< "excitation energy for root "
					<< std::fixed << std::setprecision(1) << root.excitation <<": "
					<< std::fixed << std::setprecision(10) << root.omega << " Eh         "
					<< root.omega*constants::hartree_electron_volt_relationship << " eV\n";
			msg << std::scientific;
			if(world.rank()==0)print("  oscillator strength (length)    ", osl);
			if(world.rank()==0)print("  oscillator strength (velocity)  ", osv);
			// print out the most important amplitudes
			if(world.rank()==0)print("  dominant contributions ");

		for (std::size_t p=0; p<noct; ++p) {
			const double amplitude=norms[p]*norms[p];
			if ((amplitude > 0.1)) {
				msg << "    norm(x_"<<p+parameters.freeze<<") **2  ";
				std::cout.width(10); std::cout.precision(6);
				msg << amplitude << "\n";
			}
		}
		if (world.rank()==0) print(" ");
	}

	// compute the transition densities
	const vecfuncT bra_oct=get_active_mo_bra();
	for (std::size_t i=0; i<x.size(); ++i) {
		const vecfuncT root=x[i].get_vecfunction();
		const real_function_3d td=dot(world,root,bra_oct);
		const double trace=td.trace();
		if (world.rank()==0) print("trace over transition density",i,trace);
		save(td,"transition_density_"+std::to_string(i));
	}
}

/// todo: read_from_file compatible with dist. memory computation
TDHF::Parameters::Parameters(const std::shared_ptr<SCF>& scf,const std::string& input) :
													lo(scf->param.lo) {
	read_from_file(input, "response");
	complete_with_defaults(scf);
}
TDHF::Parameters::Parameters(const std::shared_ptr<SCF>& scf) :
													lo(scf->param.lo) {
	complete_with_defaults(scf);
}

/// auto assigns all parameters which where not explicitly given and which depend on other parameters of the reference calculation
void TDHF::Parameters::complete_with_defaults(const std::shared_ptr<SCF>& scf) {
	if (econv == -1.0)
		econv = thresh;

	if (dconv == -1.0)
		dconv = thresh*10.0;

	if (guess_econv == -1.0)
		guess_econv = std::max(1.e-1, 10 * econv);

	if (guess_dconv == -1.0)
		guess_dconv = std::max(1.0, 10 * dconv);

	if (thresh_op == -1.0)
		thresh_op = std::min(1.e-4, thresh);

	if (iterating_excitations == -1)
		iterating_excitations = std::min(excitations,4);

	if (guess_excitations == -1)
		guess_excitations = std::min(excitations + iterating_excitations,2*excitations);

	if(guess_occ_to_virt<0) guess_occ_to_virt=(scf->amo.size()-freeze);
	guess_occ_to_virt = std::min(size_t(guess_occ_to_virt),(scf->amo.size()-freeze));
	if(guess_active_orbitals<0) guess_active_orbitals=(scf->amo.size()-freeze);
	guess_active_orbitals = std::min(size_t(guess_active_orbitals),(scf->amo.size()-freeze));

}
/// todo: read_from_file compatible with dist. memory computation
void TDHF::Parameters::read_from_file(const std::string input, const std::string& key) {
	{
		// get the parameters from the input file
		std::ifstream f(input.c_str());
		position_stream(f, key);
		std::string s;
		while (f >> s) {
			//std::cout << "input tag is: " << s << std::endl;
			std::transform(s.begin(), s.end(), s.begin(), ::tolower);
			//std::cout << "transformed input tag is: " << s << std::endl;
			if (s == "end")
				break;
			else if (s == "thresh") f >> thresh;
			else if (s == "freeze") f >> freeze;
			else if (s == "no_compute")
				f >> no_compute;
			else if (s == "debug")
				f >> debug;
			else if (s == "plot")
				f >> plot;
			else if (s == "restart") {
				size_t tmp;
				f >> tmp;
				restart.push_back(tmp);
			} else if (s == "exop" || s == "exop") {
				std::string tmp;
				char buf[1024];
				f.getline(buf, sizeof(buf));
				tmp = buf;
				exops.push_back(tmp);
			} else if (s == "calculation") {
				f >> calculation;
			} else if (s == "guess_occ_to_virt")
				f >> guess_occ_to_virt;
			else if(s =="guess_cm")
				f >> guess_cm;
			else if (s == "guess_active_orbitals")
				f >> guess_active_orbitals;
			else if (s == "guess_diag")
				f >> guess_diag;
			else if (s == "guess_excitations")
				f >> guess_excitations;
			else if (s == "excitations")
				f >> excitations;
			else if (s == "iterating_excitations")
				f >> iterating_excitations;
			else if (s == "guess_virtuals")
				f >> guess_virtuals;
			else if (s == "dconv_guess" or s=="guess_dconv")
				f >> guess_dconv;
			else if (s == "dconv")
				f >> dconv;
			else if (s == "econv_guess" or s=="guess_econv")
				f >> guess_econv;
			else if (s == "econv")
				f >> econv;
			else if (s == "store_potential")
				f >> store_potential;
			else if (s == "iter_max" or s =="maxiter")
				f >> maxiter;
			else if (s == "iter_guess" or s=="guess_maxiter" or s=="maxiter_guess")
				f >> guess_maxiter;
			else if (s == "damping_width")
				f >> damping_width;
			else if (s == "triplet")
				f >> triplet;
			else if (s == "kain_subspace" or s=="kain")
				f>>kain_subspace;
			else if (s == "keyval") {
				std::string key, val;
				f >> key;
				f >> val;
				generalkeyval.insert(std::make_pair(key, val));
			} else {
				std::cout << "UNKNOWN KEYWORD: " << s << "\n";
				continue;
			}
		}
	}
}

void TDHF::Parameters::print(World& world) const {
	if (world.rank() == 0) {
		std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
		std::cout << std::setfill(' ');
		std::cout << "TDHF PARAMETERS:\n";
		std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
		std::cout << std::setfill(' ');
		std::cout << std::scientific << std::setprecision(2);
		std::cout << "thresh               :" << thresh << std::endl;
		std::cout << "calculation          :" << calculation << std::endl;
		std::cout << "freeze               :" << freeze << std::endl;
		std::cout << "excitations          :" << excitations << std::endl;
		std::cout << "guess_excitations    :" << guess_excitations << std::endl;
		std::cout << "iterating_excitations:" << iterating_excitations << std::endl;
		std::cout << "dconv                :" << dconv << std::endl;
		std::cout << "econv                :" << econv << std::endl;
		std::cout << "store_potential      :" << store_potential << std::endl;
		std::cout << "maxiter              :" << maxiter << std::endl;
		std::cout << "triplet              :" << triplet << std::endl;
		if (kain_subspace > 0)
			std::cout << "kain_subspace        :" << kain_subspace << std::endl;
		else
			std::cout << "no kain is used\n";
		std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
		std::cout << "Parameters for the guess:\n";
		std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
		std::cout << "guess_virtuals       :" << guess_virtuals << std::endl;
		std::cout << "guess_diag           :" << guess_diag << std::endl;
		if(guess_maxiter==0){
			std::cout << "guess_maxiter    :"  << "no iteration of guess" << std::endl;
		}else{
			std::cout << "guess_dconv          :" << guess_dconv << std::endl;
			std::cout << "guess_econv          :" << guess_econv << std::endl;
			std::cout << "guess_maxiter        :" << guess_maxiter << std::endl;
		}
		if(guess_diag==false) std::cout << "guess_active_orbitals:" << guess_active_orbitals << std::endl;
		// if this is the case then virtuals are generated
		if(guess_virtuals!="scf" and guess_virtuals!="extern"){
			std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
			std::cout << "Parameters for the generation of virtuals\nby multiplying excitation operators to the occupied orbitals\n";
			std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
			std::cout << "guess_occ_to_virt    :" << guess_occ_to_virt << std::endl;
			std::cout << "damping_width        :" << damping_width << std::endl;
			std::cout << "guess_cm             :" << guess_cm << std::endl;
			if(guess_virtuals!="custom") 		std::cout << "predefined excitation operators are " << guess_virtuals << std::endl;
			else std::cout << "custom excitation operators are \n" << exops << "\n";
		}

		std::cout << std::setfill('-') << std::setw(50) << std::setfill('-') << "\n";
		std::cout << std::setfill(' ');
		if (generalkeyval.size() > 0) {
			std::cout << "Also found the following general key-value pairs:\n";
			for (const auto x : generalkeyval)
				std::cout << x.first << " : " << x.second;
			std::cout << "\n\n";
		}
	}
}


TDHF::~TDHF() {
	// TODO Auto-generated destructor stub
}

} /* namespace madness */
