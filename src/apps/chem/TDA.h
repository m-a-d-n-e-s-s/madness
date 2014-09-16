/*
 * TDA.h
 *
 *  Created on: Jul 14, 2014
 *      Author: kottmanj
 */

#ifndef TDA_H_
#define TDA_H_

//#include<examples/dft_solver.h>
#include <chem/projector.h>
//#include<examples/nonlinsol.h> not used anymore
#include<chem/SCF.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include <chem/TDA_XC.h>
#include <chem/TDA_exops.h>

// Kain solver
#include<examples/nonlinsol.h>

// std::sort
#include <algorithm>

namespace madness {

typedef std::vector<Function<double,3> > vecfuncT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

/// Strucutre for TIMER
struct TDA_TIMER{
	/// TDA_TIMER contructor
	/// @param[in] world the world
	/// @param[in] msg	a string that contains the desired printout when info function is called
	TDA_TIMER(World &world,std::string msg) : world(world),start_wall(wall_time()),start_cpu(cpu_time()),operation(msg){}
	World & world;
	const double start_wall;
	const double start_cpu;
	std::string operation;
	double end_wall;
	double end_cpu;
	void update_time(){
		end_wall = wall_time()-start_wall;
		end_cpu = cpu_time()-start_cpu;
	}
public:
	/// print out information about the passed time since the TDA_TIMER object was created
	void info(bool debug = true){
		if(debug==true){
			update_time();
			if(world.rank()==0) std::cout<< std::setw(40) << operation << " : " << std::scientific << std::setprecision(1)
			<< end_wall << " (wall) "<< end_cpu << " (cpu)" << std::endl;
		}
	}

};

/// The Root structure is needed by the TDA class
struct xfunction{
	/// default constructor
	/// @param[in] world	the world is needed as a reference
	xfunction(World &world) :world(world),omega(0.00001),converged(false),number(100),iterations(0),kain(false),f_length(999),f_velocity(999) {}
	/// constructs a xfunctions object and initializes the x-vecfunction (response orbitals)
	/// @param[in] world	the world is needed
	/// @param[in] x1	vectorfunction of response orbitals
	xfunction(World& world, const vecfuncT& x1) : world(world), x(x1),omega(0.00001),converged(false),number(100),iterations(0),kain(true),f_length(999),f_velocity(999) {}
	/// the copy contructor
	xfunction(const xfunction &other) : world(other.world),x(other.x),Vx(other.Vx),omega(other.omega),expectation_value(other.expectation_value),error(other.error),
			delta(other.delta),converged(other.converged),number(other.number),iterations(other.iterations),kain(other.kain),f_length(other.f_length),f_velocity(other.f_velocity){}

	World & world;
	/// the response orbitals
	vecfuncT x;
	/// the applied potentials (to save memory this will mostly be empty)
	vecfuncT Vx;
	/// the currrent excitation energy used to parametrize the BSH operator
	double omega;
	/// the expectation values (as vector so the conergence can be plotted)
	std::vector<double> expectation_value;
	/// the errors after each bsh step
	std::vector<double> error;
	/// the errors in omega after each bsh step
	std::vector<double> delta;
	/// true if the xfunctions has converged
	bool converged;
	/// the number of the xfunction
	size_t number;
	/// number of iterations already taken
	size_t iterations;
	/// true if the kain update should be used, false if a full step update should be forced
	bool kain;
	/// the residuals of the last bsh step, is needed if the kain solver should be used
	vecfuncT current_residuals;

	/// Oscillator strenght in length and velocity gauge
	/// will be calcualted after convergece, default is 999
	double f_length;
	double f_velocity;

	/// assignment operator (needed by kain)
	xfunction& operator=(const xfunction &other){
		x=other.x;
		Vx=other.Vx;
		omega = other.omega;
		expectation_value = other.expectation_value;
		error=other.error;
		delta=other.delta;
		converged=other.converged;
		number=other.number;
		iterations=other.iterations;
		kain=other.kain;
		f_length = other.f_length;
		f_velocity = other.f_velocity;

		return *this;
	}

	// Operators needed by the KAIN solver
	xfunction operator-=(const xfunction& b) {
		x=sub(world,x,b.x);
		return *this;
	}

	xfunction operator-(const xfunction &b)const {
		return xfunction(world,sub(world,x,b.x));
	}

	xfunction operator+=(const xfunction& b) { // Operator+= necessary
		x=add(world,x,b.x);
		return *this;
	}

	xfunction operator*(double a) const { // Scale by a constant necessary

		PROFILE_BLOCK(Vscale);
		xfunction result(*this);
		for (unsigned int i=0; i<x.size(); ++i) {
			result.x[i]=mul(a,x[i],false);
		}
		world.gop.fence();

		//		scale(world, x, a);
		return result;
	}

	// finally an operator that the result can be sorted after the energy
	// sorting of xfunctions should not happen during the iterations
	// therefore a warning is installed
	bool operator<=(const xfunction &b)const{return omega<=b.omega;std::cout << "WARNING XFUNCTIONS ARE SORTED" << std::endl;}
	bool operator< (const xfunction &b)const{return omega<b.omega; std::cout << "WARNING XFUNCTIONS ARE SORTED" << std::endl;}



};

struct allocator{
	World& world;
	const int noct;

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(World& world, const int nnoct) : world(world), noct(nnoct) {}

	xfunction operator()(){
		return xfunction(world,zero_functions<double,3>(world,noct));
	}
	allocator operator=(const allocator &other){
		allocator tmp(world,other.noct);
		return tmp;
	}
};

/// An inner product for the xfunction class also needed by the KAIN solver
static double inner(const xfunction &a, const xfunction &b) {
	if (a.x.size()!=b.x.size()) MADNESS_EXCEPTION("ERROR :Inner product of two xfunction structures: Different sizes in x-vectors",1);
	if (a.x.size()==0) return 0.0;
	return madness::inner(a.x[0].world(),a.x,b.x).sum();
}

// TYPEDEFS
typedef std::vector<xfunction> xfunctionsT;
typedef XNonlinearSolver<xfunction,double,allocator> solverT;

/// The structure needed if the kain solver shall be used
struct kain_solver_helper_struct{
	kain_solver_helper_struct(){}
private:
	/// number of occupied orbitals (non frozen)
	size_t noct;
	/// is kain used ?, this bool is needed because the struct has to be initialized if kain is used or not
	bool is_used;
	/// a vector of kain solvers (for each xfunction one solver) -> used when all_at_once is false
	std::vector<solverT> solver;
public:

	/// Check if everything is fine within the kain solver (use this when adding or deleting xfunctions)
	/// @param[in] xfunctions, the xfunctions of the current iteration
	void sanity_check(const xfunctionsT &xfunctions)const{
		if(xfunctions.size()!=solver.size()){
			std::cout << "KAIN SOLVER SANITY CHECK FAILED: Unequal sizes " << xfunctions.size() << " xfunctions and " << solver.size() << " solvers" << std::endl;
			MADNESS_EXCEPTION("Error in kain solver helper struct",1);
		}
		if(xfunctions[0].x.size()!=noct){
			MADNESS_EXCEPTION("Error in kain solver helper struct, unequal sizes in occupied orbitals (check freeze_ keyword)",1);
		}
	}

	/// Initialize the kain solver helper structure
	/// @param[in] world 	the world
	/// @param[in] excitations	the number of excitations planned to calculate
	/// @param[in] nnoct 	the number of occupied (non frozen) orbitals = the number of response orbitals
	/// @param[in] kain		true if kain should be used, false if not (when false it is just the default initialization)
	/// @param[in] allatonce	should all xfunctions be solved ad once or not (currently: not)
	/// @param[in] guess_iter_	number of guess iterations (not needed anymore)
	void initialize(World &world,const size_t excitations,const size_t nnoct,const bool kain){
		noct=nnoct;
		is_used = kain;
		if(kain){
			for(size_t i=0;i<excitations; i++){
				solverT onesolver(allocator(world,noct));
				onesolver.set_maxsub(3);
				solver.push_back(onesolver);
			}
		}
	}
	/// Update the current response orbitals x of the xfunctions structures
	/// @param[in] xfunctions	a vector if xfunctions structures of the current iteration (each xfunction should contain the x and curren_residuals)
	void update(xfunctionsT &xfunctions){
		if(not is_used) MADNESS_EXCEPTION("Kain update was requested, but kain should not be used",1);
		if(is_used){
			for(size_t k=0;k<xfunctions.size();k++){
				// we need ot push this back if we update or not, because of the subspace transformation
				xfunction tmp = solver[k].update(xfunction(xfunctions[k].world,xfunctions[k].x),xfunction(xfunctions[k].world,xfunctions[k].current_residuals));
				if(xfunctions[k].kain){
					xfunctions[k].x = tmp.x;
				}else{
					std::cout << "(convergence too low) forcing full step for  " << k << std::endl;
					xfunctions[k].x = sub(xfunctions[k].world,xfunctions[k].x,xfunctions[k].current_residuals);
				}
			}
		}
	}
	/// Transform the Kain subspace(s)
	/// @param[in] world 	the world
	/// @param[in] U		The transformation matrix
	void transform_subspace(World &world,const madness::Tensor<double> U){
		if(solver[0].get_ulist().empty()){ std::cout<<"kain: empty ulist, no transformation in "; return;}
		if(solver[0].get_rlist().empty()){ std::cout<<"kain: empty rlist, no transformation in "; return;}
		for(int k=0;k<solver[0].get_ulist().size();k++){
			std::vector<vecfuncT> all_x;
			std::vector<vecfuncT> all_r;
			for(size_t i=0;i<solver.size();i++){
				all_x.push_back(solver[i].get_ulist()[k].x);
				all_r.push_back(solver[i].get_rlist()[k].x);
			}
			std::vector<vecfuncT> new_x = transform_vecfunctions(world,all_x,U);
			std::vector<vecfuncT> new_r = transform_vecfunctions(world,all_r,U);
			for(size_t i=0;i<solver.size();i++){
				solver[i].get_ulist()[k].x = new_x[i];
				solver[i].get_rlist()[k].x = new_r[i];
			}
		}
	}
	/// reduce the subspace: delete whole solvers
	/// Use this if you delete xfunctions
	/// @param[in] i	the number of the xfunction that was deleted
	void reduce_subspace(const size_t i){
		solver.erase(solver.begin()+i);
	}
	/// Only erase the subspace, not delete the solver
	/// use this if you replace xfunctions
	/// @param[in] i 	the number of the solver/xfunction which subspace should be erased
	void erase_subspace(const size_t i){
		if(solver.size()<i) MADNESS_EXCEPTION("Tried to delete subspace of nonexisting solver",1);
		solver[i].get_ulist().clear();
		solver[i].get_rlist().clear();
	}
	/// increase the subspace: adding new solvers
	/// use this if you add new xfunctions
	/// @param[in] world 	the world
	/// @param[in] xfunctions the current xfunctions (which have been extended)
	void increase_subspace(World &world,const xfunctionsT &xfunctions){
		// if the subspace is increased (means more solvers are added) the subspaces of the existing solvers must be erased
		// else there will be problems with the next subspace transformation
		// the simplest solution is to add cimpletely new solvers for all xfunctions

		// clean up: erase all old solvers
		solver.clear();
		// add new solvers
		solverT onesolver(allocator(world,noct));
		onesolver.set_maxsub(3);
		for(size_t j=0;j<xfunctions.size();j++){
			solver.push_back(onesolver);
		}
		sanity_check(xfunctions);
	}
	/// Helper function for the transform_subspace function of the structure
	/// @param[in] world the world
	/// @param[in] xfunctions a vector of vectorfunctions
	/// @param[in] U the transformation matrix
	/// @param[out] U*xfunctions : the transformed vector of vectorfunctions
	std::vector<vecfuncT> transform_vecfunctions(World &world,const std::vector<vecfuncT> &xfunctions,const madness::Tensor<double> U)const{

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
};

/// Functor that adds diffuse 1s functions on given coordinates with given signs (phases)
struct diffuse_functions : public FunctionFunctorInterface<double,3> {
public:
	/// constructor
	/// @param[in] L the box size
	/// @param[in] bc the molecular bounding cube
	/// @param[in] coord the coordinates (in most cases of nuclei) where the diffuse functions shall be centered
	/// @param[in] signs the signs the diffuse functions should have on the corresponding coordinates
	/// @param[in] natoms the number of atoms (if the coordinates are of nuclei) or just the size of the coord vector
	/// @param[in] mode the level of diffuseness (0 and 1 possible)
	diffuse_functions(const double L,const std::vector<coord_3d> coord,const std::vector<int> signs, const size_t natoms, const size_t mode):
		L(L), coord(coord), signs(signs),natoms(natoms),mode(mode)
{if(mode>1) MADNESS_EXCEPTION("mode in diffuse_function struct is not 0,1 or 2",1);}

	/// Make a rydberg guess function where the near range area is 0.0 and the long range is a diffuse 2s or 2p function
	double operator()(const coord_3d &r)const{
		std::vector<double> diffuse_1s_functions;
		double exponent =1.0;
		if(mode==0) exponent = 20.0;
		if(mode==1) exponent = 15.0;
		for(size_t i=0;i<natoms;i++){
			if (signs[i]!=0){
				// make diffuse 1s functions localized at the atoms
				double x = r[0]-coord[i][0];
				double y = r[1]-coord[i][0];
				double z = r[2]-coord[i][0];
				double rad = sqrt(x*x+y*y+z*z);
				double diffuse_tmp = signs[i]*exp(-exponent/L*rad);
				diffuse_1s_functions.push_back(diffuse_tmp);
			}else diffuse_1s_functions.push_back(0.0);
		}
		double result=0;
		for(size_t i=0;i<diffuse_1s_functions.size();i++) result+=diffuse_1s_functions[i];
		return result;
	}

private:
	/// Box size
	const double L;
	/// The coordinates of the atoms
	const std::vector<coord_3d> coord;
	/// The signs of the respoinse function at the atom coordinates (0 is node)
	const std::vector<int> signs;
	/// number of atoms
	const size_t natoms;
	/// the mode of diffuseness (from 0 to 2)
	const size_t mode;
	/// The diffuse 2s and 2p functions
	double diffuse_2s (const coord_3d &r)const{return r_function(r)*exp(-10.0/L*r_function(r));}
	double diffuse_2px(const coord_3d &r)const{return          r[0]*exp(-10.0/L*r_function(r));}
	double diffuse_2py(const coord_3d &r)const{return          r[1]*exp(-10.0/L*r_function(r));}
	double diffuse_2pz(const coord_3d &r)const{return          r[2]*exp(-10.0/L*r_function(r));}
	double diffuse_1s (const coord_3d &r)const{return               exp(-10.0/L*r_function(r));}
	/// Helper function to evaluate the radius
	double r_function(const coord_3d &r)const{return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);}

};



/// The TDA class: computes TDA and CIS calculations
class TDA {
public:
	/// the constructor
	/// @param[in] world	the world
	/// @param[in] calc 	the SCF calcualtion
	/// @param[in] mos		the occupied molecular orbitals from the scf calculation
	/// @param[in] input	name of the input file
	/// @param[in] lowt		will be used later to ditinguish between low and high threshold computations (not yet implemented)
	TDA(World &world,const SCF &calc,const vecfuncT &mos,const std::string input):
		world(world),
		dft_(false),
		calc_(calc),
		mos_(mos),
		print_grid_(false),
		guess_("physical"),
		guess_iter_(15),
		guess_mode_("physical"),
		guess_exop_("quadrupole"),
		guess_excitations_(6),
		excitations_(4),
		bsh_eps_(1.e-5),
		iter_max_(100),
		econv_(1.e-4),
		dconv_(1.e-3),
		hard_dconv_(1.e-3),
		nfreeze_(0),
		plot_(false),
		debug_(false),
		only_fock_(false),
		only_GS_(false),
		on_the_fly_(true),
		read_(false),
		xclib_interface_(world,calc),
		ipot_(0.0),
		kain_(false),
		shift_(0.0),
		safety_(1.0)
{
		setup(mos,input);
}
	/// reads the input file and calculates needed functions
	void setup(const vecfuncT &mos,const std::string input){


		// so that the thresh can be changed from the outside
		mos_ = mos;

		size_t noct = calc_.aeps.size();
		// The highest possible excitation (-homo_energy)
		double highest_excitation_default = -calc_.aeps(noct-1);
		highest_excitation_ = highest_excitation_default;
		ipot_ = -calc_.aeps(noct-1)*2.0;


		// The guessed lowest excitation (if no guess_omega_ is in the input)
		double guess_omega_default = -0.1*calc_.aeps[noct-1];
		guess_omega_ = guess_omega_default;

		std::ifstream f(input.c_str());
		position_stream(f, "TDA");
		std::string s, tag;
		while (std::getline(f,s)) {
			std::istringstream ss(s);
			ss >> tag;
			if (tag == "end") break;
			else if (tag == "dft") dft_=true;
			else if (tag == "excitations") ss >> excitations_;
			else if (tag == "guess") ss >> guess_;
			else if (tag == "hard_dconv") ss >> hard_dconv_;
			else if (tag == "guess_iter") ss >> guess_iter_;
			else if (tag == "guess_omega") ss >> guess_omega_;
			else if (tag == "guess_mode") ss >> guess_mode_;
			else if (tag == "guess_exop") ss >> guess_exop_;
			else if (tag == "guess_excitations") ss >> guess_excitations_;
			else if (tag == "bsh_eps") ss >> bsh_eps_;
			else if (tag == "iter_max") ss >> iter_max_;
			else if (tag == "econv") ss >> econv_;
			else if (tag == "dconv") ss >> dconv_;
			else if (tag == "freeze") ss >> nfreeze_;
			else if (tag == "print_grid") print_grid_=true;
			else if (tag == "plot") plot_=true;
			else if (tag == "debug") debug_=true;
			else if (tag == "only_fock") only_fock_=true;
			else if (tag == "only_GS") only_GS_=true;
			else if (tag == "highest_excitation") ss >> highest_excitation_;
			else if (tag == "no_otf") on_the_fly_=false;
			else if (tag == "read") read_ = true;
			else if (tag == "ipot") ss >> ipot_;
			else if (tag == "kain") kain_=true;
			else if (tag == "exop1") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop2") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop3") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop4") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop5") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop6") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop7") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop8") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop9") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "exop10") {std::string tmp; ss >> tmp; custom_exops_.push_back(tmp);}
			else if (tag == "truncate_safety") ss>>safety_;
			else continue;
		}

		// make potential shift = -ipot - homo
		if(dft_) shift_= -ipot_ - get_calc().aeps[noct-1];
		highest_excitation_=highest_excitation_-shift_;


		if (world.rank() == 0) {
			std::cout<< std::setw(60) <<"\n\n\n\n ======= TDA info =======\n\n\n" << std::endl;
			if (nfreeze_==0) std::cout<< std::setw(40) <<"# frozen orbitals : "<<"none" << std::endl;
			if (nfreeze_>0) std::cout<< std::setw(40) <<"# frozen orbitals : " <<  "0 to " << nfreeze_-1 << std::endl;
			std::cout<< std::setw(40) <<"active orbitals : " << nfreeze_ << " to " << calc_.param.nalpha-1 << std::endl;
			std::cout<< std::setw(40) << "guess from : " << guess_ << std::endl;
			std::cout<< std::setw(40) << "Gram-Schmidt is used : " << !only_fock_ << std::endl;
			std::cout<< std::setw(40) << "threshold 3D : " << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout<< std::setw(40) << "energy convergence : " << econv_ << std::endl;
			std::cout<< std::setw(40) << "max residual (dconv) : " << dconv_ << std::endl;
			std::cout<< std::setw(40) << "number of excitations : " << excitations_ << std::endl;
			std::cout<< std::setw(40) << "number of guess excitations : " << guess_excitations_ << std::endl;
			std::cout<< std::setw(40) << "guessed lowest extitation energy : " << guess_omega_ << std::endl;
			std::cout<< std::setw(40) << "guessed excitation operators : " << guess_exop_ << std::endl;
			std::cout<< std::setw(40) << "highest possible excitation : " << highest_excitation_default << std::endl;
			std::cout<< std::setw(40) << "used highest possible excitation : " << highest_excitation_ << std::endl;
			std::cout<< std::setw(40) << "guessed ionization potential is : " << ipot_ << std::endl;
			std::cout<< std::setw(40) << "potential shift is : " << shift_ << std::endl;
			std::cout<< std::setw(40) << "guessed lowest excitation : " << guess_omega_default << std::endl;
			std::cout<< std::setw(40) << "chosen lowest excitation : " << guess_omega_ << std::endl;
			std::cout<< std::setw(40) << "orthonormalization : ";
			if(only_fock_) std::cout << "only perturbed fock matrix"<< std::endl;
			else if(only_GS_) std::cout << "only Gram-Schmidt"<< std::endl;
			else std::cout << "use both"<< std::endl;
			std::cout<< std::setw(40) << "potential calculation : " << "on_the_fly is " << on_the_fly_ << std::endl;
			std::cout<< std::setw(40) << "use KAIN : " << kain_ << std::endl;
		}



		lo=get_calc().param.lo;
		bsh_eps_ = FunctionDefaults<3>::get_thresh()*0.1;
		for(size_t i=nfreeze_;i<mos_.size();i++){active_mo_.push_back(mos_[i]);}

		// Initialize the projector on the occupied space
		Projector<double,3> projector(mos_);
		rho0 = projector;

		// make the unperturbed density (closed shell)
		real_function_3d active_density = real_factory_3d(world);
		for(size_t i=0;i<active_mo_.size();i++){active_density += 2.0*active_mo_[i]*active_mo_[i];}
		active_density_ = active_density;
		real_function_3d density = real_factory_3d(world);
		for(size_t i=0;i<mos_.size();i++){density+=2.0*mos_[i]*mos_[i];}
		density_=density;
		std::cout <<std::setw(40) <<"Norm of unperturbed density is : " << density_.norm2() << std::endl;

		// Initialize the exchange intermediate
		if(not dft_) {
			exchange_intermediate_ = make_exchange_intermediate();
			std::cout << std::setw(40) << "CIS is used" << " : LIBXC Interface is not initialized" << std::endl;
		}

		// Initialize the KAIN solvers
		kain_solvers.initialize(world,excitations_,noct,kain_);


		// Prevent misstakes:
		if(shift_>0){MADNESS_EXCEPTION("Potential shift is positive",1);}
		if(not dft_ and shift_ !=0.0){MADNESS_EXCEPTION("Non zero potential shift in TDHF calculation",1);}

		if(only_fock_ and only_GS_){
			print("\nWARNING: only_fock and only_GS demanded ...use both");
			only_fock_ = false;
			only_GS_ = false;
		}

		// make the truncate thresh
		truncate_thresh_ = FunctionDefaults<3>::get_thresh() * safety_;
		std::cout << "Truncate threshold is set to " << truncate_thresh_ << std::endl;

		// Truncate the current mos
		truncate(world,mos_,truncate_thresh_);
		std::cout << "truncate molecular orbitals to " << truncate_thresh_ << std::endl;

		std::cout << "setup of TDA class ended\n" << std::endl;
	}

	/// try to gain a little bit information about the used memory

	double memwatch(const xfunctionsT &xfunctions,const bool printout)const{
		// sanity_check
		if(xfunctions.empty())return 0.0;
		double allx=0.0; double allVx=0.0; double allr=0.0;
		if(printout)std::cout << "\n\n#" << "  " << "     x " << "     Vx " << "     r " << std::endl;
		if(printout)print("-------------------------------");
		for(size_t i=0;i<xfunctions.size();i++){
			if(on_the_fly_ and not xfunctions[i].Vx.empty()) MADNESS_EXCEPTION("on the fly calculation used but Vx not empty",1);
			if(not kain_ and not xfunctions[i].current_residuals.empty()) MADNESS_EXCEPTION("no kain is used but current residuals are not empty",1);

			// mem information
			double x_size = get_size(world,xfunctions[i].x);
			double Vx_size= get_size(world,xfunctions[i].Vx);
			double r_size=get_size(world,xfunctions[i].current_residuals);
			allx+=x_size; allVx+=Vx_size; allr=r_size;
			if(printout)std::cout << i << "  " << x_size <<" "<< Vx_size <<" "<< r_size << " (GB)" <<  std::endl;
		}
		if(printout)print("-------------------------------");
		if(printout)std::cout << "all" << "  " << allx <<" "<< allVx <<" "<< allr << " (GB)\n\n" <<  std::endl;
		return allx+allVx+allr;
	}

	//virtual ~TDA();

	/// Solves the CIS or TDA equations
	void solve(xfunctionsT &xfunctions);

	/// Solves the CIS or TDA equations sequentially for a set of preconverged xfunctions
	void solve_sequential(xfunctionsT xfunctions);

	/// Returns the MolDFT calulation
	const SCF &get_calc() const {return calc_;}

	// Print out grid (e.g for Koala or other external program)
	const bool print_grid_TDA() const {return print_grid_;}

	// returns a shallow copy the converged xfunctions
	xfunctionsT get_converged_xfunctions(){return converged_xfunctions_;}

private:

	/// The World
	World & world;

	/// DFT or HF Calculation
	/// for TDA calculations currently only LDA works
	bool dft_;

	/// The SCF calculation of MolDFT
	const SCF &calc_;

	/// The molecular orbitals of moldft
	/// extra member variable that the thresh can be changed without changing calc_
	vecfuncT mos_;

	/// Print grid option
	bool print_grid_;

	/// Options for the guess calculation

	/// guess == physical is the only implementation left
	/// new guess functions can be implemented and called in the intialize function
	std::string guess_;
	/// guess iterations are the first iterations where the energy is kept fixed at the guess_omega energy
	size_t guess_iter_;
	double guess_omega_;

	/// mode is either mo or all_orbitals (decides on which of the two functions the excitation operators act)
	/// mo is the default, all_orbitals mode can increase the freedom (if there are convergence problems) of the guess functions
	std::string guess_mode_;

	/// Excitation operator for the guess functions (bsp "dipole" or "quadrupole" which will be dipole + quadrupole operators)
	std::string guess_exop_;
	/// how many excitations should pre_converge (recommended: 1-2 more than demanded in the end)
	size_t guess_excitations_;
	std::vector<std::string> custom_exops_;

	/// Number of excitations to be caluclated
	size_t excitations_;

	/// Thresholds and convergence cirteria
	double bsh_eps_;

	/// maximal iterations per guess_function
	size_t iter_max_;
	/// energy convergence level for the guess functions in the solve routine
	double econv_;
	/// maximal residual for the guess_functions in the solve routine
	double dconv_;
	// Convergence criteria (residual norm) for the high thresh sequential iterations in the end
	double hard_dconv_;

	/// Frozen Orbitals
	size_t nfreeze_;

	/// Many Plots
	bool plot_;

	/// More output
	bool debug_;

	/// use only the fock orthonormalization procedure (default)
	bool only_fock_;
	/// use only Gram-Schmidt orthonormalization (not recommended)
	bool only_GS_;

	/// The highest possible excitation to calculate (higher values will result in positive eigenvalues for the BSH operator)
	double highest_excitation_;

	double lo;

	/// Vector of active molecular orbitals
	vecfuncT active_mo_;

	// the projector on the unperturbed density
	Projector<double,3> rho0;

	/// Active density (closed shell)
	real_function_3d active_density_;

	/// Complete density (includes frozen orbitals), closed shell
	real_function_3d density_;

	/// The potential is calculated when needed and then deleted (saves memory but the potential has to be calculated more often)
	bool on_the_fly_;

	/// only read and analyze functions
	bool read_;

	/// Iterate the read xfunctions sequentially
	bool sequential_;

	/// The interface to XCLIB library
	TDA_DFT xclib_interface_;

	/// Ionization potential for the potential shift used in TDDFT calculations to get bound states for the first excitations (default is -2.0*homo)
	double ipot_;

	/// Use the Kain solver to update the functions
	bool kain_;
	kain_solver_helper_struct kain_solvers;

	/// The potential shift for the unperturbed DFT potential when using TDDFT (shift = -ipot_ -homo)
	double shift_;

	/// The truncate threshold (default is the default threshold)
	double truncate_thresh_;

	/// Truncate threshold as factor ot the detault thresh (safety)
	double safety_;

	/// The unperturbed dft potential;
	real_function_3d unperturbed_vxc_;

	/// the intermediate is the same for all roots:
	/// \[
	///   int[p,i] = \int 1/r12 \phi_i(1) * \phi_p(1)
	/// \]
	/// with p \in noct, i \in nocc
	std::vector<vecfuncT> exchange_intermediate_;

	/// the coulomb potential
	mutable real_function_3d coulomb_;

	/// The converged xfunctions
	std::vector<xfunction> converged_xfunctions_;

	/// Print the current xfunctions in a formated way
	/// @param[in] xfunctions, a vector of xfunction structures
	void print_status(const xfunctionsT & xfunctions)const;

	/// just a helper function for print_status and others
	/// @param[in] xfunction, a single xfunction structure
	/// the function will print out the information of the xfunction structure (energy, iterations, convergence ...) in a formated way
	void print_xfunction(const xfunction &x)const;

	/// Takes an empty vector of excitation functions and passes it to one of the guess functions
	/// @param[in] xfunctions, empty vector of xfunctions (no necessarily empty)
	void initialize(xfunctionsT & xfunctions);

	/// Creates physical guess functions (x,y,z excitations - depending on the input file, see make_excitation_operators function)
	void guess_physical(xfunctionsT & xfunctions);

	/// Add diffuse 1s functions to the molecular orbitals (for LDA calculations)
	void add_diffuse_functions(vecfuncT &mos);

	/// Create excitation operators (e.g x,y,z for dipole excitations bzw symmetry operators)
	/// @param[out] gives back a vectorfunction of excitation operators (specified in the input file)
	/// e.g the keyword: "dipole+" in the TDA section of the input file will give back a vecfunction containing (x,y,z,r)
	vecfuncT make_excitation_operators()const;

	/// iterates the guess_functions
	// Calls iterate_all with the right settings
	void iterate_guess(xfunctionsT &xfunctions);

	/// iterations after guess initialisation
	// Calls iterate_all with the right settings
	void iterate(xfunctionsT &xfunctions);

	/// iterates the xfunctions
	// guess: guess inititialisation or final iterations
	/// @param[in] xfunctions, all excitation structures
	/// @param[in] guess, for the first iterations (no energy update, no kain update)
	void iterate_all(xfunctionsT &xfunctions,bool guess);

	/// Update process for one xfunction
	/// @param[in] xfunction, a single xfunction structure (contains the response orbitals as vecfunc x)
	/// @param[in] ptfock, this should be true if orthonormalize_fock was used before (if false, the function will calculate the expectation value of the xfunction)
	/// @param[in] guess, true if this is one of the very first iterations (no energy update, no kain update)
	void iterate_one(xfunction & xfunction,bool ptfock,bool guess);

	/// basicaly the same than iterate_one
	// instead of x = G(-2VPsi) and G parametrized with epsilon+omega
	// here: x= G(-2VPsi + 2omega x) and G aprametrized only with epsilon (always bound)
	void iterate_one_unbound(xfunction & xfunction,bool ptfock,bool guess);

	void iterate_one_yanai(xfunction & xfunction,bool guess);

	/// Update energies (decide if second order or expectation value should be used)
	void update_energies(xfunctionsT &xfunctions);

	/// Normalize one or all excitation functions
	void normalize(xfunctionsT &xfunctions);
	void normalize(xfunction &xfunction);

	/// Project out the converged xfunctions
	void project_out_converged_xfunctions(xfunctionsT & xfunctions);

	/// Orthonormalize the exfunctions with Gram-Schmidt
	void orthonormalize_GS(xfunctionsT &xfunctions);

	/// Orthonormalize the xfunction with the perturbed Fock Matrix
	// 1. Call make_perturbed_fock_matrix(xfunctions)
	// 2. Diagonalize
	// 3. Update Energy and xfunctions
	/// @param[in] the xfunctions
	/// @param[in] guess : is it a guess iterations (the first iterations where the energy is fixed or not
	/// @param[out] true if the fock update was carried out, false if not
	bool orthonormalize_fock(xfunctionsT &xfunctions,const bool guess);

	/// a little helper routine to measure the degree of offdiagonality in a 2d tensor
	double measure_offdiagonality(const madness::Tensor<double> &U,const size_t size)const;

	std::vector<vecfuncT> transform_vecfunctions(const std::vector<vecfuncT> &xfunctions,const madness::Tensor<double> U)const;

	/// Projects out the occupied space
	// 1. Make projector
	// 2. apply
	void project_out_occupied_space(vecfuncT &x);

	/// Calculate offdiagonal elements of the perturbed fock matrix
	double perturbed_fock_matrix_element(const vecfuncT &xr, const vecfuncT &Vxp,const vecfuncT &xp)const;

	/// Calculate the expectation value and update xfunction.expectation_value
	// Can also be used to calculate diagonal elements of the fock matrix
	double expectation_value(const xfunction &x,const vecfuncT &Vx);

	/// Make the perturbed Fock Matrix
	// 1. Calculate potential (V0 + Gamma) --> Call perturbed_potential(exfunctions)
	// 2. Calculate Matrix Elements
	Tensor<double> make_perturbed_fock_matrix(const xfunctionsT &xfunctions)const;

	/// Calculate the perturbed Potential (V0 + Gamma)
	// 1. Call get_V0
	// 2. Call apply_gamma or apply_gamma_dft
	// 3. return V0 + Gamma
	vecfuncT apply_perturbed_potential(const xfunction & xfunction)const;

	/// Calculate the perturbed vector potential for CIS or TDA calculations
	/// @param[in] xfunction, a single xfunction structure which contains the response orbital vector x
	/// @param[out] the applied perturbed potential (applied on the unperturbed MOs) given back as vector function
	vecfuncT apply_gamma(const xfunction &xfunction)const;
	vecfuncT apply_gamma_dft(const xfunction &xfunction)const;

	/// The perturbed Hartree potential is the same for TDA and CIS
	/// @param[in] x vectorfunction of response orbitals
	/// @param[out] the applied perturbed hartree potential
	/// the function will evaluate the perturbed density and then calculate the hartree potential
	vecfuncT apply_hartree_potential(const vecfuncT &x)const;

	/// Create the exchange intermediate
	// This has to be done just one time because only the unperturbed orbitals are needed
	std::vector<vecfuncT> make_exchange_intermediate()const;

	/// Return the unperturbed fock potential as vector potential : V0*x_p
	vecfuncT get_V0(const vecfuncT &x)const;

	/// Return the coulomb potential of the moldft calculation
	real_function_3d get_coulomb_potential() const;

	/// Return the unperturbed exchange-correlation functional (for dft calculations)
	real_function_3d get_vxc_potential()const;

	/// Plot vectorfunction (for convenience)
	/// @param[in] x the vectorfunction to plot
	/// @param[in] msg the name for the vecfunction plot
	/// @param[in] plot, if false nothing is done
	void plot_vecfunction(const vecfuncT &x,std::string msg, bool plot = true)const;

	/// Check convergence
	/// checks if the xfunctions have converged
	/// if so then the converged flag will be set so zero
	/// also the converged xfunction is pushed into the converged_xfunctions_ vector and is replaced by a new guess function
	/// @param[in] xfunctions a vector of xfunction structures
	bool check_convergence(xfunctionsT &xfunctions);

	/// Print performance: Values of expectation values and errors of each iteration into a file
	/// @param[in] xfunctions a vector of xfunction structures
	/// @param[in] string the saved file will be prename+results.tex
	void print_performance(const xfunctionsT &xfunctions,const std::string prename)const;

	/// Truncate the xfunctions structure:
	/// @param[in] xfunctions a vector of xfunction structures
	/// Truncates the vecfunctions: x, Vx and the current_residual (if not empty)
	void truncate_xfunctions(xfunctionsT &xfunctions);

	/// load a converged root from disk

	/// compute the oscillator strength in the length representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_length(const xfunction& xfunction) const;

	/// compute the oscillator strength in the velocity representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_velocity(const xfunction& root) const;

	/// analyze the root: oscillator strength and contributions from occ
	void analyze(xfunctionsT& roots) const;

	void save_xfunctions(const xfunctionsT &xfunctions)const;
	bool read_xfunctions(xfunctionsT &xfunctions);
};

} /* namespace madness */

#endif /* TDA_H_ */
