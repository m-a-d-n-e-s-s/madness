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

/// GaussFunctor so let the exciation operators go to zero at the boundaries
/// 1S symmetry: totally symmetric
struct gauss_functor : public FunctionFunctorInterface<double,3> {
public:
	gauss_functor();
	gauss_functor(const double& width): width_(width){
		MADNESS_ASSERT(not(width<0.0));
	}
	gauss_functor(const double& width, const coord_3d c): width_(width), center(c){
		MADNESS_ASSERT(not(width<0.0));
	}
	gauss_functor(const double& width, const Tensor<double> c): width_(width), center(tensor_to_coord(c)){
		MADNESS_ASSERT(not(width<0.0));
	}
	const double width_;
	const coord_3d center=coord_3d();
	coord_3d tensor_to_coord(const Tensor<double>& t)const{
		coord_3d result;
		MADNESS_ASSERT(t.size()>=3);
		for(size_t i=0;i<3;++i) result[i]=t[i];
		return result;
	}
	double operator()(const coord_3d &rr)const{
		coord_3d r;
		r[0]=rr[0]-center[0];
		r[1]=rr[1]-center[1];
		r[2]=rr[2]-center[2];
		if(width_<=0.0) return 1.0;
		else{
			const double prefactor = 0.06349363593424097/(width_*width_*width_);
			const double exponent=0.5/(width_*width_);
			return prefactor*exp(-exponent*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
		}
		return 1.0;
	}


};

/// Project a general 3D polynomial to the MRA Grid
struct polynomial_functor : public FunctionFunctorInterface<double,3> {
public :
	polynomial_functor(const std::string input, const double& damp_width=0.0, const coord_3d& c=coord_3d()) : input_string_(input), data_(read_string(input)), dampf(damp_width), center(c) {}
	polynomial_functor(const std::string input, const double& damp_width, const Tensor<double>& c) : input_string_(input), data_(read_string(input)), dampf(damp_width), center(dampf.tensor_to_coord(c)) {}

	double operator()(const coord_3d &rr)const{
		coord_3d r;
		r[0]=rr[0]-center[0];
		r[1]=rr[1]-center[1];
		r[2]=rr[2]-center[2];
		return dampf(r)*make_polynom(r);
	}
	double make_polynom(const coord_3d& r)const{
		double result =0.0;
		for(size_t i=0;i<data_.size();i++){
			if(data_[i].size()!=4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, empty data_ entry",1);
			result += ( data_[i][3]*pow(r[0],data_[i][0])*pow(r[1],data_[i][1])*pow(r[2],data_[i][2]) );
		}
		return result;
	}
private:
	const std::string input_string_;
	/// The data for the construction of the polynomial chain
	/// every entry of data_ is vector containing the threee exponents and the coefficient of a monomial dx^ay^bz^c , data_[i] = (a,b,c,d)
	const std::vector<std::vector<double>> data_;
	gauss_functor dampf;
	coord_3d center=coord_3d();
public:
	std::vector<std::vector<double> > read_string(const std::string string)const{
		std::stringstream line(string);
		std::string name;
		size_t counter = 0;
		std::vector<double> current_data = vector_factory(0.0,0.0,0.0,1.0);
		std::vector<std::vector<double> > read_data;
		while(line>>name){
			if(name=="c") line>>current_data[3];
			else if(name=="x") line>>current_data[0];
			else if(name=="y") line>>current_data[1];
			else if(name=="z") line>>current_data[2];
			else if(name==","){
				counter++; read_data.push_back(current_data); current_data = vector_factory(0.0,0.0,0.0,1.0);
			}
		}
		// dont forget the last read polynomial
		read_data.push_back(current_data);
		return read_data;
	}
	void test(){
		std::cout << "Test polynomial functor " << "\n input string is " << input_string_ << std::endl;
		std::cout << "\n read data is \n" << data_ << std::endl;
	}
	std::vector<std::vector<double> > give_data(){return data_;}
};



/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};


coord_3d compute_centroid(const real_function_3d& f){
	coord_3d result(0.0);
	real_function_3d density = f*f;
	const double integral= density.trace();
	for(size_t x=0;x<3;++x){
		const auto mf = xyz(x);
		real_function_3d m=real_factory_3d(f.world()).functor(mf);
		result[x]=(m*density).trace()/integral;
	}
	return result;
}
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
	if(world.rank()==0) std::cout << " is_dft() = " << nemo.get_calc()->xc.is_dft() << "\n";
	if(world.rank()==0) std::cout << " hf_coeff = " << nemo.get_calc()->xc.hf_exchange_coefficient() << "\n";
	if(world.rank()==0) std::cout << " do_pcm() = " << nemo.do_pcm() << "\n";
	if(world.rank()==0) std::cout << " do_ac()  = " << nemo.do_ac() << "\n";

	parameters.print(world);

	FunctionDefaults<3>::set_thresh(parameters.thresh);
	if(world.rank()==0) std::cout << "MRA Threshold is set to: " << FunctionDefaults<3>::get_thresh() << " with k=" << FunctionDefaults<3>::get_k() << "\n";

	if (not parameters.no_compute) {

		if(nemo.get_calc()->xc.hf_exchange_coefficient()!=0.0){
			msg.subsection("Computing Exchange Intermediate");
			CCTimer timer(world,"Computing ExIm");
			g12.update_elements(mo_bra_,mo_ket_);
			timer.info();
		}else msg.output("No Exchange Intermediate Computed\n");

		msg.output("Orbital Energies of Reference");
		const Tensor<double> eps=nemo.get_calc()->aeps;
		if(world.rank()==0) std::cout << eps << "\n";

	}
	if(nemo.get_calc()->param.localize){
		Fock F(world, &nemo);
		F_occ = F(get_active_mo_bra(),get_active_mo_ket());
		for(size_t i=0;i<get_active_mo_ket().size();++i){
			std::cout << std::scientific << std::setprecision(10);
			if(world.rank()==0) std::cout << "F(" << i << "," << i << ")=" << F_occ(i,i) << "\n";
			if(std::fabs(get_orbital_energy(i+parameters.freeze)-F_occ(i,i))>1.e-5){
				if(world.rank()==0) std::cout << "eps(" << i << ")=" << get_orbital_energy(i) << " | diff=" << get_orbital_energy(i+parameters.freeze)-F_occ(i,i) << "\n";
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
	orthonormalize(guess,empty);


	// failsafe
	if(guess.size()<parameters.guess_excitations){
		std::string msg=("WARNING: You demanded: " + std::to_string(parameters.guess_excitations)
		+ " Guess vectors, but your demanded guess has only "
		+ std::to_string(guess.size()) + "vectors. So we will not iterate the first vectors and then do the same guess again ... this might be unstable").c_str();
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

	std::vector<CC_vecfunction> final_vectors;
	msg.subsection("Iterate Guess Vectors");
	{
		// do guess iterations
		iterate_cis_guess_vectors(guess_vectors);
		// sort according to excitation energies
		std::sort(guess_vectors.begin(),guess_vectors.end());
		// save
		for(size_t i=0;i<guess_vectors.size();i++) guess_vectors[i].save_functions(std::to_string(i));
		// truncate again
		for(size_t i=0;i<parameters.excitations;i++) final_vectors.push_back(guess_vectors[i]);
	}

	msg.subsection("Iterate Final Vectors");
	iterate_cis_final_vectors(final_vectors);
	msg.section("CIS CALCULATIONS ENDED");
	std::sort(final_vectors.begin(),final_vectors.end());
	// information
	if(world.rank()==0) std::cout << std::setfill('-') << std::setw(25) << "\n" << std::setfill(' ');
	if(world.rank()==0) std::cout << "Results of CIS Calculation: Excitation, Excitation Energy, WF-Error \n";
	for(size_t i=0;i<final_vectors.size();i++){
		if(world.rank()==0) std::cout << "Excitation " << std::fixed << std::setprecision(1) << i <<": "
				<< std::fixed << std::setprecision(10) << final_vectors[i].omega
				<< std::scientific << std::setprecision(2)  << ", " << final_vectors[i].current_error << "\n";
	}
	if(world.rank()==0) std::cout << std::setfill('-') << std::setw(25) << "\n" << std::setfill(' ');
	time.info();
	// plot final vectors
	for(size_t i=0;i<final_vectors.size();i++) final_vectors[i].plot(std::to_string(i)+"_converged_cis");
	// save final vectors
	for(size_t i=0;i<final_vectors.size();i++) final_vectors[i].save_functions(std::to_string(i));
	// assign solution vectors
	//start = final_vectors;

	return final_vectors;
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
	if(world.rank()==0) std::cout
			<< "dconv    = " << dconv << "\n"
			<< "econv    = " << econv << "\n"
			<< "iter_max = " << iter_max << "\n"
			<< "kain     = " << kain << "\n";

	// set up the kain solvers ... if needed or not
	std::vector<XNonlinearSolver<vecfuncT,double,TDHF_allocator> > solvers;
	// initialize solvers
	if(kain){
		for(size_t i=0;i<x.size();i++){
			XNonlinearSolver<vecfuncT,double,TDHF_allocator>  solver(TDHF_allocator(world,x[i].size()),true);
			solver.set_maxsub(parameters.kain_subspace);
			solvers.push_back(solver);
		}
	}

	bool converged = true;

	// get potentials (if demanded)
	std::vector<vecfuncT> V;

	// for TDHF the potential is always stored
	if(parameters.store_potential and y.empty()) V = make_potentials(x);
	else if(not y.empty()) V = make_tdhf_potentials(x,y);
	// orthonormalize
	orthonormalize(x,V);

	for(size_t iter=0;iter<iter_max;iter++){
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
			errors.push_back(av_error);
			largest_errors.push_back(std::make_pair(std::distance(ind_err.begin(),it),*it));
			// convergece criteria are the individual functions (orhterwise we have a dependence on the number of functions)
			if((*it)>dconv) converged=false;
		}
		// update, store old omegas in the deltas vector
		std::vector<double> deltas;
		for(size_t i=0;i<x.size();i++){
			vecfuncT new_x0;
			if(kain){
				new_x0 = solvers[i].update(x[i].get_vecfunction(),residuals[i],10.0*parameters.thresh,5.0);
			}else{
				new_x0 = sub(world,x[i].get_vecfunction(),residuals[i]);
			}
			vecfuncT Q_new_x0 = Q(new_x0);
			truncate(world,Q_new_x0);
			CC_vecfunction new_x1(Q_new_x0,x[i].type,parameters.freeze);

			new_x1.omega = x[i].omega;
			x[i]=new_x1;
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
		if(world.rank()==0) std::cout << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
		if(world.rank()==0) std::cout << "Iteration " << iter <<": omega, largest error (number) , vector error, delta "<< "\n";
		bool econv=true;
		for(size_t i=0;i<errors.size();i++){
			if(world.rank()==0) std::cout << "Excitation " << std::fixed << std::setprecision(1) << i <<": "
					<< std::fixed << std::setprecision(10)
			<< x[i].omega
			<< std::scientific << std::setprecision(2)  << ", "
			<< largest_errors[i].second << "(" << largest_errors[i].first  <<"), "
			<<errors[i] << ", "
			<< x[i].delta
			<< "\n";
			x[i].current_error=errors[i];

			// check energy convergence
			if(fabs(x[i].delta)>econv) converged=false;
		}

		if(world.rank()==0) std::cout << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');

		if(converged) break;
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
		time_N.info();

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

		vecfuncT GV = apply(world, bsh, Vi);

		vecfuncT residual = sub(world,x[i].get_vecfunction(),GV);
		result.push_back(residual);

		// Calculate Second Order Energy Update
		{
			// Inner product of Vpsi and the residual (Vi is scaled to -2.0 --> multiply later with 0.5)
			double tmp = inner(world,make_bra(residual),Vi).sum();
			// squared norm of GVpsi (Psi_tilde)
			double tmp2 = inner(world,make_bra(GV),GV).sum();

			// Factor 0.5 removes the factor 2 from the scaling before
			x[i].delta = (0.5 * tmp / tmp2);
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
		if(world.rank()==0 and parameters.debug) std::cout << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
		const vecfuncT pot = get_tda_potential(xi);
		V.push_back(pot);
		if(world.rank()==0 and parameters.debug) std::cout << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
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
		if(world.rank()==0) std::cout << "TDA Potential is " << xc_data << ", hf_coeff=" << hf_coeff << ", pcm is=" << pcm <<  "\n";
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
	CCTimer time(world,"Orthonormalization");

	// make the overlap matrix
	Tensor<double> S = make_overlap_matrix(x);
	if(parameters.debug) std::cout << "The Overlap Matrix\n " << S << "\n";

	//make the Hamilton matrix for the vectorfunctions
	Tensor<double> F = make_perturbed_fock_matrix(x,V);

	// Diagonalize the F Matrix
	Tensor<double> U, evals;
	Tensor<double> dummy(x.size());
	U = nemo.get_calc() -> get_fock_transformation(world, S, F, evals, dummy, 2.0*parameters.thresh);

	if(parameters.debug){
		std::cout << "Perturbed Fock-Matrix Eigenvalues\n";
		for(int x=0;x<evals.size();++x) if(world.rank()==0) std::cout << evals(x) << "\n";
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
		CC_vecfunction tmp(new_x,x[k].type,parameters.freeze);
		tmp.omega=x[k].omega;
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
	if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\nOverlap Matrix\n" << S << "\n";
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
		if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\n(T+V) Matrix\n" << T << "\n";
		if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\nPotential Matrix\n" << MV << "\n";
		if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
	}
	timeF.stop();
	if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
	//formated timings output
	timeT.print();
	timeV.print();
	timeR.print();
	timeF.print();

	// symmetryze
	F = 0.5*(F + transpose<double>(F));
	if(parameters.debug and world.rank()==0) std::cout << std::fixed << std::setprecision(5) << "\nSymmetrized Perturbed Fock Matrix\n" << F << "\n";

	return F;
}

/// Makes the guess for the y-states from empty y states
std::vector<CC_vecfunction> TDHF::make_y_guess(const std::vector<CC_vecfunction> & x, std::vector<CC_vecfunction> & y)const{
	std::vector<vecfuncT> V;
	std::vector<CC_vecfunction> result;
	MADNESS_EXCEPTION("Not Implemented",1);
	return result;
}

/// Makes the (old) guess functions by exciting active orbitals with excitation operators
std::vector<CC_vecfunction> TDHF::make_old_guess(const vecfuncT& f)const{
	CCTimer time(world,"Making Guess Functions: " + parameters.guess_virtuals);
	std::vector<std::string> exop_strings;
	if(parameters.guess_virtuals=="custom"){
		exop_strings = parameters.exops;
		if(world.rank()==0){
			std::cout << "Custom Excitation Operators Demanded:\n";
			std::cout << exop_strings << "\n";
		}
	}
	else exop_strings = make_predefined_guess_strings(parameters.guess_virtuals);

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
	} else if(parameters.guess_virtuals == "ppw"){
		// product plane wave (similar as in Baker, Burke, White 2018
		int order = 1;
		if(parameters.generalkeyval.find("ppw_order")!=parameters.generalkeyval.end()) order = std::stoi(parameters.generalkeyval.find("ppw_order")->second);
		if(order!=1) msg.output("Found ppw_order keyval\nTake Care: a lot of virtuals will be created!");
		vecfuncT xmo;
		for(size_t i=0;i<parameters.guess_occ_to_virt;++i) xmo.push_back(get_active_mo_ket()[get_active_mo_ket().size()-1-i]);
		virtuals = make_ppw_virtuals(xmo,order);

	} else	// use the old guess format and convert to normalized virtuals (result is essentially the same as ppw but with polynomials instead of plane waves;
	{
		std::vector<CC_vecfunction> cguess = make_old_guess(get_active_mo_ket());
		for (const auto& xg : cguess) {
			for (const auto& f : xg.get_vecfunction()) {
				if (f.norm2() > 0.1) {
					const double norm = sqrt(inner(make_bra(f), f));
					virtuals.push_back(1.0 / norm * f);
				}
			}
		}
	}

	if(parameters.guess_cm>0.0){
		// add center of mass diffuse functions
		const double factor=parameters.guess_cm;
		const double width = (-1.0*factor/(get_orbital_energy(mo_ket_.size()-1)));
		msg.subsection("adding center of mass functions with exponent homo/c and c=" + std::to_string(factor));
		msg.output("width="+std::to_string(width));

		Tensor<double> cm = nemo.get_calc()->molecule.center_of_mass();
		if(world.rank()==0) std::cout << "center of mass is " << cm << "\n";
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
		std::cout << "created " << virtuals.size() << " virtuals\n";

	virtuals = Q(virtuals);
	if (parameters.plot) {
		CCTimer timep(world,"plot");
		for(size_t i=0;i<virtuals.size();++i){
			plot_plane(world, virtuals[i], "virt_"+std::to_string(i));
			save(virtuals[i],"virt_"+std::to_string(i));
		}
		timep.print();
	}
	if (parameters.debug) {
		for (const auto& x : virtuals)
			x.print_size("virtual");
	}
	world.gop.fence();
	time.print();
	return virtuals;
}

madness::vecfuncT TDHF::make_ppw_virtuals(const vecfuncT& seed, const int& order) const {
	//const int nvirt = seed.size() * ((2 * order * 2 * order * 2 * order) - 1);
	//	msg.subsection("creating a set of " + std::to_string(nvirt) + " virtuals by multiplying functions with plane waves");

	// compute the centers of the seed functions
	CCTimer time_centers(world,"compute centers");
	std::vector<coord_3d> centers;
	for(const auto& s:seed){
		const coord_3d c=compute_centroid(s);
		centers.push_back(c);
	}
	time_centers.print();
	if(world.rank()==0) std::cout << "Centers of the seed functions are:\n" << centers << "\n";

	CCTimer timec(world, "initialize exops");
	vecfuncT virtuals;// = zero_functions<double, 3>(world, nvirt);
	// create the exops (plane waves) and initialize the virtuals with the right seed functions
	std::vector<ExopUnaryOpStructure> exops;

	// make sin and cos guess instead of x and x2
	for (size_t i=0;i<seed.size();++i) {
		const coord_3d& c=centers[i];
		std::vector<double> n1 = { 1.0,0.0,0.0 };
		std::vector<double> n2 = { 0.0,1.0,0.0 };
		std::vector<double> n3 = { 0.0,0.0,1.0 };

		std::vector<double> n4 = { 0.0,1.0,1.0 };
		std::vector<double> n5 = { 1.0,0.0,1.0 };
		std::vector<double> n6 = { 1.0,1.0,0.0 };;

		std::vector<double> n7 = { 1.0,1.0,1.0 };

		std::vector<bool> cosinus = {true,true,true };
		std::vector<bool> sinusx = {false,true,true };
		std::vector<bool> sinusy = {true,false,true };
		std::vector<bool> sinusz = {true,true,false };
		std::vector<bool> sinusxz = {false,true,false };
		std::vector<bool> sinusyz = {true,false,false };
		std::vector<bool> sinusxy = {false,false,true };

		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n1, sinusx, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n2, sinusy, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n3, sinusz, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n4, sinusyz, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n5, sinusxz, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n6, sinusxy, c))));
		exops.push_back(ExopUnaryOpStructure(std::shared_ptr<PlaneWaveFunctor>(new PlaneWaveFunctor(n7, cosinus, c))));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
		virtuals.push_back(copy(seed[i],false));
	}



	if(world.rank()==0) std::cout << " will create " << virtuals.size() << " virtuals\n";
	timec.print();

	MADNESS_ASSERT(exops.size()==virtuals.size());

	world.gop.fence();

	timec.print();
	CCTimer timeex(world, "excite");
	const size_t nvirt=virtuals.size();
	for(size_t v=0;v<nvirt;++v){
		virtuals[v].unaryop(exops[v],false);
	}
	world.gop.fence();
	timeex.print();
	CCTimer timep(world, "apply Q");
	virtuals = Q(virtuals);
	for(size_t v=0;v<nvirt;++v){plot_plane(world,virtuals[v],"Qvirt_"+std::to_string(v));}
	timep.print();
	return virtuals;
}

vector<CC_vecfunction> TDHF::make_guess_from_initial_diagonalization() const {
	CCTimer time(world, "make_guess_from_initial_diagonalization");
	//convenience
	const int nact = get_active_mo_ket().size();

	// create virtuals
	vecfuncT virtuals = make_virtuals();
	// canonicalize virtuals
	virtuals = canonicalize(virtuals);

	// compute the CIS matrix
	Tensor<double> MCIS=make_cis_matrix(virtuals);

	// initialize the guess functions
	if(world.rank()==0 and MCIS.dim(0)<parameters.guess_excitations){
		msg.warning("asd"+std::to_string(1));
		msg.warning(std::to_string(parameters.guess_excitations)
		+ " guess vectors where demanded, but with the given options only "
		+std::to_string(MCIS.dim(0))+" can be created\n");
	}
	std::vector<CC_vecfunction> xfunctions;
	for (int x = 0; x < MCIS.dim(0); ++x) {
		if(x>=parameters.guess_excitations) break;
		CC_vecfunction init(zero_functions<double, 3>(world, nact), RESPONSE, parameters.freeze);
		xfunctions.push_back(init);
	}

	{
		Tensor<double> U, evals;
		CCTimer time_diag(world,"cis-matrix diagonalization");
		syev(MCIS, U, evals);
		time_diag.print();
		if (world.rank() == 0) {
			std::cout << "Initial Diagonalization of CIS Matrix:\n Lowest Eigenvalues are \n";
			for (int x = 0; x < std::max(std::min(int(evals.size()),10),parameters.guess_excitations); ++x) std::cout << evals(x) << "\n";
		}
		CCTimer time_assemble(world,"assemble guess vectors");
		// make x functions from amplitudes and virtuals
		// todo: probably not optimal for highly parallel applications
		const int nvirt = virtuals.size();
		auto get_com_idx = [nvirt](int i, int a) { return i*nvirt+a; };
		auto get_vir_idx = [nvirt](int I) {return I%nvirt;};
		auto get_occ_idx = [nvirt](int I) {return I/nvirt;};
		for(size_t I=0;I<MCIS.dim(0);++I){
			if(I>=parameters.guess_excitations) break;
			const int a=get_vir_idx(I);
			const int i=get_occ_idx(I);
			if(evals(I) < 0.0 and world.rank()==0) msg.warning("NEGATIVE EIGENVALUE IN INITIAL DIAGONALIZATION: CHECK YOUR REFERENCE!\n");
			if (evals(I) < 1.e-5) {
				if (world.rank() == 0)
					std::cout << "skipping root " << evals(I) << " \n";
				continue;
			}
			for(size_t J=0;J<MCIS.dim(1);++J){
				const int b=get_vir_idx(J);
				const int j=get_occ_idx(J);
				const double xjb = U(J, I);
				xfunctions[I].get_vecfunction()[j] += xjb * virtuals[b];
				xfunctions[I].omega = evals(I);
				xfunctions[I].excitation = I;
			}
		}
		time_assemble.print();
		CCTimer time_truncate(world,"truncate guess");
		for(auto& x:xfunctions){
			vecfuncT tmp = x.get_vecfunction();
			truncate(world,tmp,parameters.thresh);
			// should be truncated by shallow copy anyways ... but just to be sure
			x.set_functions(tmp,x.type,parameters.freeze);
		}
		time_truncate.print();
	}
	if(parameters.debug){
		Tensor<double> S = make_overlap_matrix(xfunctions);
		if (world.rank() == 0) std::cout << "Overlap matrix of guess:\n" << S << "\n";
	}

	print_xfunctions(xfunctions,true);

	if(world.rank()==0) std::cout << "created " << xfunctions.size() << " guess vectors\n";

	time.print();
	return xfunctions;
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

/// Makes an automated excitation operator string for the excitation operators needed to create virtuals from the reference orbitals
std::vector<std::string> TDHF::make_auto_polynom_guess(const size_t order)const{
	std::vector<std::string> exop_strings;
	for(size_t i=0; i<order+1; i++){
		for(size_t j=0; j<order+1 ; j++){
			for(size_t k=0;k<order+1 ; k++){
				if(i+j+k > order) MADNESS_ASSERT(i+j+k >order) ; // do nothing
				else if(i==0 and j==0 and k==0) MADNESS_ASSERT(i==0) ; // do nothing
				else{
					if(i==0 and j!=0 and k!=0) exop_strings.push_back(" y " + madness::stringify(j) + " z " + madness::stringify(k) );
					else if(j==0 and i!=0 and k!=0) exop_strings.push_back("x " + madness::stringify(i) + " z " + madness::stringify(k) );
					else if(k==0 and i!=0 and j!=0) exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j));
					else if(i==0 and j==0) exop_strings.push_back(" z " + madness::stringify(k) );
					else if(i==0 and k==0) exop_strings.push_back(" y " + madness::stringify(j));
					else if(j==0 and k==0) exop_strings.push_back("x " + madness::stringify(i));
					else exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j) + " z " + madness::stringify(k) );
				}
			}
		}
	}
	return exop_strings;
}
/// Makes an excitation operator string based on predefined keywords
std::vector<std::string> TDHF::make_predefined_guess_strings(const std::string what)const{
	std::vector<std::string> exop_strings;
	if(what == "dipole"){
		exop_strings.resize(3);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
	}else if(what == "x"){
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0";
	}else if(what == "y"){
		exop_strings.resize(1);
		exop_strings[0] = "y 1.0";
	}else if(what == "z"){
		exop_strings.resize(1);
		exop_strings[0] = "z 1.0";
	}else if(what == "r2"){
		exop_strings.resize(1);
		exop_strings[0] = "x 2.0 , y 2.0 , z 2.0";
	}else if(what == "quadrupole"){
		exop_strings.resize(9);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0";
		exop_strings[4] = "y 2.0";
		exop_strings[5] = "z 2.0";
		exop_strings[6] = "x 1.0 y 1.0";
		exop_strings[7] = "x 1.0 z 1.0";
		exop_strings[8] = "y 1.0 z 1.0";
	}else if(what == "dipole+"){
		exop_strings.resize(4);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
	}else if(what == "dipole+diffuse"){
		exop_strings.resize(4);
		exop_strings[0] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[1] = "x 2.0 y 1.0 , y 3.0 , y 1.0 z 2.0";
		exop_strings[2] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[3] = "x 4.0 , y 4.0 , z 4.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	}else if(what == "dipole+diffuse_big"){
		exop_strings.resize(8);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
		exop_strings[4] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[5] = "x 2.0 y 1.0 , y 3.0 ,y 1.0 z 2.0";
		exop_strings[6] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[7] = "x 4.0 , 4 2.0 , 4 2.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	}else if(what == "c2v"){
		exop_strings.resize(4);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0 , x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[1] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0 , x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[2] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0 , x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[3] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0, x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0 ";
	}else if(what == "water_first"){
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0 y 1.0, x 1.0 y 1.0 z 1.0";
	}else if(what == "c2v_big"){
		exop_strings.resize(8);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0";
		exop_strings[1] = "x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[2] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0";
		exop_strings[3] = "x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[4] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0";
		exop_strings[5] = "x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[6] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0";
		exop_strings[7] = "x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0";
	}else if(what == "big_fock"){exop_strings = make_auto_polynom_guess(6);
	}else if(what == "small_fock"){exop_strings = make_auto_polynom_guess(4);
	}else if(what == "big_fock_2"){exop_strings = make_auto_polynom_guess(2);
	}else if(what == "big_fock_3"){exop_strings = make_auto_polynom_guess(3);
	}else if(what == "big_fock_4"){exop_strings = make_auto_polynom_guess(4);
	}else if(what == "big_fock_5"){exop_strings = make_auto_polynom_guess(5);
	}else if(what == "big_fock_6"){exop_strings = make_auto_polynom_guess(6);
	}else if(what == "big_fock_7"){exop_strings = make_auto_polynom_guess(7);
	}else if(what == "big_fock_8"){exop_strings = make_auto_polynom_guess(8);
	}else if(what == "big_fock_9"){exop_strings = make_auto_polynom_guess(9);
	}else std::cout << "Keyword " << what << " is not known" << std::endl;
	return exop_strings;
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

		std::cout << std::scientific << std::setprecision(10) << std::setw(20);
		if (world.rank()==0) {
			std::cout << "excitation energy for root "
					<< std::fixed << std::setprecision(1) << root.excitation <<": "
					<< std::fixed << std::setprecision(10) << root.omega << " Eh         "
					<< root.omega*constants::hartree_electron_volt_relationship << " eV\n";
			std::cout << std::scientific;
			print("  oscillator strength (length)    ", osl);
			print("  oscillator strength (velocity)  ", osv);
			// print out the most important amplitudes
			print("  dominant contributions ");
		}
		for (std::size_t p=0; p<noct; ++p) {
			const double amplitude=norms[p]*norms[p];
			if (world.rank()==0 and (amplitude > 0.1)) {
				std::cout << "    norm(x_"<<p+parameters.freeze<<") **2  ";
				std::cout.width(10); std::cout.precision(6);
				std::cout << amplitude << std::endl;
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
	if (guess_econv == -1.0)
		guess_econv = std::max(1.e-1, 10 * econv);

	if (guess_dconv == -1.0)
		guess_dconv = std::max(1.0, 10 * dconv);

	if (thresh_op == -1.0)
		thresh_op = std::min(1.e-4, thresh);

	if (iterating_excitations == -1)
		iterating_excitations = std::min(excitations, 10);

	if (guess_excitations == -1)
		guess_excitations = excitations + 2;

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
		std::cout << "dconv_guess          :" << guess_dconv << std::endl;
		std::cout << "dconv                :" << dconv << std::endl;
		std::cout << "econv_guess          :" << guess_econv << std::endl;
		std::cout << "econv                :" << econv << std::endl;
		std::cout << "store_potential      :" << store_potential << std::endl;
		std::cout << "iter_max             :" << maxiter << std::endl;
		std::cout << "iter_guess           :" << guess_maxiter << std::endl;
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
