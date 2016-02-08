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
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/CISOperators.h>
#include <chem/CCOperators.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
//#include <chem/TDA_XC.h>
//#include <madness/world/print.h>

//#include <chem/TDA_exops.h>
#include <chem/TDA_guess.h>

// Kain solver
#include <examples/nonlinsol.h>

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
	xfunction(World &world) :world(world),omega(0.00001),converged(false),number(100),iterations(0),kain(false),
			f_length(999),f_velocity(999) {error.push_back(999);delta.push_back(999);expectation_value.push_back(999);guess_excitation_operator="not initialized";}
	xfunction(World &world, const double in_omega) :world(world),omega(in_omega),converged(false),number(100),iterations(0),kain(false),
			f_length(999),f_velocity(999) {error.push_back(999);delta.push_back(999);expectation_value.push_back(999);guess_excitation_operator="not initialized";}
	/// constructs a xfunctions object and initializes the x-vecfunction (response orbitals)
	/// @param[in] world	the world is needed
	/// @param[in] x1	vectorfunction of response orbitals
	xfunction(World& world, const vecfuncT& x1) : world(world), x(x1),omega(0.00001),converged(false),number(100),iterations(0),kain(true),
			f_length(999),f_velocity(999) {error.push_back(999);delta.push_back(999);expectation_value.push_back(999);guess_excitation_operator="not initialized";}
	/// the copy contructor
	xfunction(const xfunction &other) : world(other.world),x(other.x),smooth_potential(other.smooth_potential),omega(other.omega),expectation_value(other.expectation_value),error(other.error),
			delta(other.delta),converged(other.converged),number(other.number),iterations(other.iterations),kain(other.kain),
			f_length(other.f_length),f_velocity(other.f_velocity),guess_excitation_operator(other.guess_excitation_operator){}

	World & world;
	/// the response orbitals
	vecfuncT x;
	/// the applied potentials (to save memory the nuclear potential is missing)
	vecfuncT smooth_potential;
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

	/// A string which determines the guess excitation operator
	std::string guess_excitation_operator;

	/// assignment operator (needed by kain)
	xfunction& operator=(const xfunction &other){
		x=other.x;
		smooth_potential=other.smooth_potential;
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
		guess_excitation_operator = other.guess_excitation_operator;

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
	bool operator<=(const xfunction &b)const{return expectation_value.back()<=b.expectation_value.back() ;std::cout << "WARNING XFUNCTIONS ARE SORTED" << std::endl;}
	bool operator< (const xfunction &b)const{return expectation_value.back()<b.expectation_value.back(); std::cout << "WARNING XFUNCTIONS ARE SORTED" << std::endl;}



};

///// This is a structure to perform operations on a vector of xfunctions
//struct vector_of_xfunctions{
//public:
//	vector_of_xfunctions(const size_t active_element_size): active_element_size(active_element_size){}
//	vector_of_xfunctions(const vector_of_xfunctions &other) : active_element_size(other.active_element_size),active_elements(other.active_elements),remaining_elements(other.remaining_elements),converged_elements(other.converged_elements){}
//	void sort(){
//		std::sort(active_elements.begin(),active_elements.end());
//		std::sort(remaining_elements.begin(),remaining_elements.end());
//		std::sort(converged_elements.begin(),converged_elements.end());
//	}
//	std::vector<xfunction> get_converged_xfunctions(){return converged_elements;}
//	std::vector<xfunction> get_active_xfunctions(){return active_elements;}
//	std::vector<xfunction> get_all_xfunctions(){
//		std::vector<xfunction> all_xfunctions = converged_elements;
//		for(auto x:active_elements) all_xfunctions.push_back(x);
//		for(auto x:remaining_elements) all_xfunctions.push_back(x);
//		return all_xfunctions;
//	}
//	void push_to_converged(){
//		std::vector<xfunction> not_converged;
//		for(auto x:active_elements){
//			if(x.converged)converged_elements.push_back(x);
//			else not_converged.push_back(x);
//		}
//		active_elements = not_converged;
//	}
//	void fill_up(const size_t i){
//		if(active_elements.size()>=i) return;
//		else if(remaining_elements.empty()) return;
//		else {
//			for(size_t k=0;k<i;k++){
//				active_elements.push_back(remaining_elements.front());
//				remaining_elements.erase(remaining_elements.begin());
//				if(remaining_elements.empty())break;
//				if(active_elements.size()==i)break;
//			}
//		}
//	}
//	void clear(){
//		active_elements.clear();
//		converged_elements.clear();
//		remaining_elements.clear();
//	}
//	void set(const std::vector<xfunction> &x){
//		clear();
//		for(size_t i=0;i<active_element_size;i++) active_elements.push_back(x[i]);
//		for(size_t i=active_element_size;i<x.size();i++) remaining_elements.push_back(x[i]);
//	}
//	void reset(){
//		std::vector<xfunction> all_xfunctions = get_all_xfunctions();
//		std::sort(all_xfunctions.begin(),all_xfunctions.end());
//		set(all_xfunctions);
//
//	}
//	void print_status(){
//		std::cout << "\n" <<std::setw(5) << " #" << std::setw(20) << "omega" << std::setw(20) << "delta" << std::setw(20)
//		<< "error"<<std::setw(20)
//		<<"expv" << std::setw(7) <<"iter"<< std::setw(7)<< "conv" << std::endl;
//		std::cout << "_._._._(pre) converged xfunctions" << std::endl;
//		for(auto x:converged_elements) print_xfunction(x);
//		std::cout << "_._._._active xfunctions"<< std::endl;
//		for(auto x:active_elements) print_xfunction(x);
//		std::cout << "_._._._remaining xfunctions" << std::endl;
//		for(auto x:remaining_elements) print_xfunction(x);
//	}
//	void print_xfunction(const xfunction &x){
//		std::cout << std::setw(5) << x.number;
//		std::cout << std::scientific << std::setprecision(10) << std::setw(20) << x.omega << std::setw(20)<< x.delta.back()
//																											<< std::setw(20)<< x.error.back()<< std::setw(20) << x.expectation_value.back();
//		std::cout << std::fixed <<std::setw(7)<< x.iterations << "   " << std::setw(7)<<x.converged << std::endl;
//	}
//	xfunction operator()(const size_t i){return active_elements[i];}
//
//private:
//	const size_t active_element_size;
//	std::vector<xfunction> active_elements;
//	std::vector<xfunction> remaining_elements;
//	std::vector<xfunction> converged_elements;
//};
//
//
//
/// Kain allocator for single roots
struct TDA_allocator{
	World& world;
	const int noct;

	/// @param[in]	world	the world
	/// @param[in]	nnoct	the number of functions in a given vector
	/// @todo validate doxygen on `nnoct`
	TDA_allocator(World& world, const int nnoct) : world(world), noct(nnoct) {}

	xfunction operator()(){
		return xfunction(world,zero_functions<double,3>(world,noct));
	}
	TDA_allocator operator=(const TDA_allocator &other){
		TDA_allocator tmp(world,other.noct);
		return tmp;
	}
};
//
//
//
/// An inner product for the xfunction class also needed by the KAIN solver
static double inner(const xfunction &a, const xfunction &b) {
	if (a.x.size()!=b.x.size()) MADNESS_EXCEPTION("ERROR :Inner product of two xfunction structures: Different sizes in x-vectors",1);
	if (a.x.size()==0) return 0.0;
	return madness::inner(a.x[0].world(),a.x,b.x).sum();
}

// TYPEDEFS
typedef std::vector<xfunction> xfunctionsT;
typedef XNonlinearSolver<xfunction,double,TDA_allocator> sequential_kain_solver;

/// The structure needed if the kain solver shall be used

/// Functor that smoothes guess functions with the error functions (no fluctuations at the box borders)
//struct guess_smoothing : public FunctionFunctorInterface<double,3> {
//private:
//	/// The size of the smoothing box (rectangular function, borders must be at dyadic points)
//	const double box_size_;
//public:
//	guess_smoothing(const double box_size) : box_size_(box_size) {}
//	// Smoothing function
//	//	double operator()(const coord_3d &r)const{
//	//		return 0.5*(erf(-(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])-box_size_))+1.0);
//	//	}
//
//	double operator()(const coord_3d &r)const{
//		if(fabs(r[0])>box_size_) return 0.0;
//		else if(fabs(r[1])>box_size_) return 0.0;
//		else if(fabs(r[2])>box_size_) return 0.0;
//		else return 1.0;
//	}
//};

/// Structure that makes the excitation operators in polynomial form from strings
struct polynomial_exop_functor : public FunctionFunctorInterface<double,3> {
public :
	polynomial_exop_functor(const std::string input) : input_string_(input), data_(read_string(input)) {}

	double operator()(const coord_3d &r)const{
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

/// Structure that multiplicates two polynomial strings (needed for fock matrix of guess functions)
struct multiply_polynomials_functor : public FunctionFunctorInterface<double,3> {
public:
	multiply_polynomials_functor(const std::string polynomial1, const std::string polynomial2){
		polynomial_exop_functor dummy(polynomial1);
		std::vector<std::vector<double>> data1 = dummy.read_string(polynomial1);
		std::vector<std::vector<double>> data2 = dummy.read_string(polynomial2);
		for(size_t i=0;i<data1.size();i++){
			for(size_t j=0;j<data2.size();j++){
				std::vector<double> new_data(4);
				new_data[0] = data1[i][0]+data2[j][0];
				new_data[1] = data1[i][1]+data2[j][1];
				new_data[2] = data1[i][2]+data2[j][2];
				new_data[3] = data1[i][3]*data2[j][3];
				data_.push_back(new_data);
			}
		}
		std::cout << " \ndata1\n " << data1 << " \ndata2\n " << data2 << std::endl;
		test();

	}
	double operator()(const coord_3d &r)const{
		double result =0.0;
		for(size_t i=0;i<data_.size();i++){
			if(data_[i].size()!=4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, empty data_ entry",1);
			result += ( data_[i][3]*pow(r[0],data_[i][0])*pow(r[1],data_[i][1])*pow(r[2],data_[i][2]) );
		}
		return result;
	}
	void test(){
		std::cout << "\n Multiplied functor made polynomial\n " << data_ << std::endl;
	}
private :
	std::vector<std::vector<double>> data_;
};


/// The TDA class: computes TDA and CIS calculations
class TDA {
public:
	/// the constructor
	/// @param[in] world	the world
	/// @param[in] mos		the occupied molecular orbitals from the scf calculation
	/// @param[in] input	name of the input file
	TDA(World &world,const Nemo &nemo,const vecfuncT &mos,const std::string input):
		orbital_energies_(nemo.get_calc()->aeps),
		parameters(input, 1.e-7),
		world(world),
		dft_(false),
		nemo_(nemo),
		mos_(mos),
		CCOPS_(CIS_Operators(world,nemo,mos)),
		print_grid_(false),
		guess_("dipole+"),
		solve_iter_(5),
		guess_iter_(3),
		guess_mode_("projected"),
		replace_guess_functions_(true),
		guess_excitations_(0),
		excitations_(8),
		iterating_excitations_(0),
		iter_max_(20),
		econv_(parameters.thresh_3D),
		guess_econv_(10.0*parameters.dconv_3D), // econv_6D is always looser than econv 3D
		dconv_(10.0*parameters.dconv_3D),
		guess_dconv_(100.0*parameters.dconv_3D), // same as with econv
		hard_dconv_(parameters.dconv_3D),
		hard_econv_(0.1*parameters.thresh_3D),
		plot_(false),
		only_fock_(false),
		only_GS_(false),
		ipot_(0.0),
		kain_(false),
		kain_subspace_(3),
		shift_(0.0),
		triplet_(false),
		use_omega_for_bsh_(true),
		compute_virtuals_(false)
{
		setup(mos,input);
}
	/// reads the input file and calculates needed functions
	void setup(const vecfuncT &mos,const std::string input){

		// so that the thresh can be changed from the outside
		mos_ = mos;

		size_t noct = orbital_energies_.size();
		// The highest possible excitation (-homo_energy)
		double highest_excitation_default = -orbital_energies_(noct-1);
		highest_excitation_ = highest_excitation_default;
		ipot_ = -orbital_energies_(noct-1)*1.2;


		// The guessed lowest excitation (if no guess_omega_ is in the input)
		double guess_omega_default = -0.99*orbital_energies_[noct-1];
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
			else if (tag == "iterating_excitations") ss >> iterating_excitations_;
			else if (tag == "guess") ss >> guess_;
			else if (tag == "hard_dconv") ss >> hard_dconv_;
			else if (tag == "hard_econv") ss >> hard_econv_;
			else if (tag == "solve_iter") ss >> solve_iter_;
			else if (tag == "guess_iter") ss >> guess_iter_;
			else if (tag == "guess_omega") ss >> guess_omega_;
			else if (tag == "guess_mode") ss >> guess_mode_;
			else if (tag == "replace_guess_functions") ss >> replace_guess_functions_;
			else if (tag == "guess_excitations") ss >> guess_excitations_;
			else if (tag == "iter_max") ss >> iter_max_;
			else if (tag == "econv") ss >> econv_;
			else if (tag == "guess_econv") ss >> guess_econv_;
			else if (tag == "dconv") ss >> dconv_;
			else if (tag == "guess_dconv") ss >> guess_dconv_;
			else if (tag == "print_grid") print_grid_=true;
			else if (tag == "plot") plot_=true;
			else if (tag == "debug") parameters.debug=true;
			else if (tag == "only_fock") only_fock_=true;
			else if (tag == "only_GS") only_GS_=true;
			else if (tag == "highest_excitation") ss >> highest_excitation_;
			else if (tag == "ipot") ss >> ipot_;
			else if (tag == "kain") kain_=true;
			else if (tag == "kain_subspace") ss>> kain_subspace_;
			else if (tag == "exop") {std::string tmp;char buf[1024];ss.getline(buf,sizeof(buf));tmp=buf; custom_exops_.push_back(tmp);}
			else if (tag == "triplet") triplet_=true;
			else if (tag == "compute_virtuals") compute_virtuals_ = true;
			else if (tag == "dft") dft_ = true;
			else continue;
		}
		if(compute_virtuals_){
			if(parameters.freeze != mos.size()-1){
				if(world.rank()==0) std::cout << "Virtual orbital calculation demanded: Freeze Key is set to number_of_mos -1 which is " << mos.size()-1 << std::endl;
				MADNESS_EXCEPTION("Freeze key is set wrong",1);
			}
		}



		// this will be the case if guess_excitations are not assigned
		if(guess_excitations_ == 0) guess_excitations_ = excitations_;
		if(guess_excitations_ < excitations_){
			if(world.rank()==0) std::cout << "WARNING: More converged excitations than guess excitations demanded ... correcting that " << std::endl;
			guess_excitations_ = excitations_;
		}

		if(iterating_excitations_==0) iterating_excitations_ = guess_excitations_;

		// make potential shift = -ipot - homo
		if(dft_ and ipot_>0.0) shift_= -ipot_ - orbital_energies_[noct-1];
		highest_excitation_=highest_excitation_-shift_;

		if(guess_ =="koala"){
			if(replace_guess_functions_){
			if(world.rank()==0) std::cout << "For the koala guess the guess functions will not be replaced after convergece \n"
					<< std::endl;
			replace_guess_functions_ = false;
			}
		}

		if (world.rank() == 0) {
			std::cout << std::scientific << std::endl;
			std::cout<< std::setw(60) <<"\n\n\n\n ======= TDA info =======\n\n\n" << std::endl;
			if (parameters.freeze==0) std::cout<< std::setw(40) <<"frozen orbitals : "<<"none" << std::endl;
			if (parameters.freeze>0) std::cout<< std::setw(40) <<"frozen orbitals : " <<  "0 to " << parameters.freeze-1 << std::endl;
			std::cout<< std::setw(40) <<"active orbitals : " << parameters.freeze << " to " << mos_.size()-1 << std::endl;
			std::cout<< std::setw(40) << "guess from : " << guess_ << std::endl;
			std::cout<< std::setw(40) << "Gram-Schmidt is used : " << !only_fock_ << std::endl;
			std::cout<< std::setw(40) << "threshold 3D : " << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout<< std::setw(40) << "energy convergence : " << econv_ << std::endl;
			std::cout<< std::setw(40) << "max residual (dconv) : " << dconv_ << std::endl;
			std::cout<< std::setw(40) << "number of final excitations : " << excitations_ << std::endl;
			std::cout<< std::setw(40) << "number of guess excitations : " << guess_excitations_ << std::endl;
			std::cout<< std::setw(40) << "number of parallel iterating excitations : " << iterating_excitations_ << std::endl;
			std::cout<< std::setw(40) << "guess_iter : " << guess_iter_<< std::endl;
			std::cout<< std::setw(40) << "solve_iter : " << solve_iter_ << std::endl;
			std::cout<< std::setw(40) << "guessed lowest extitation energy : " << guess_omega_ << std::endl;
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
			//std::cout<< std::setw(40) << "potential calculation : " << "on_the_fly is " << on_the_fly_ << std::endl;
			std::cout<< std::setw(40) << "use KAIN : " << kain_ << std::endl;
			std::cout<< std::setw(40) << "triplet is " << triplet_ << std::endl;
			std::cout<< std::setw(40) << "dft is " << dft_ << std::endl;
		}

		// Make the active_mos_ vector
		for(size_t i=parameters.freeze;i<mos_.size();i++){active_mo_.push_back(mos_[i]);}

		// project the mos if demanded (default is true)
		if(guess_mode_ != "numerical"){
			active_mos_for_guess_calculation_ = project_to_ao_basis(active_mo_,get_nemo().get_calc() -> ao);
		}else active_mos_for_guess_calculation_ = active_mo_;

		/// Make transformation matrix from cannical to localized MOs
		std::vector<int> set=get_nemo().get_calc() -> group_orbital_sets(world,get_nemo().get_calc() -> aeps,get_nemo().get_calc() -> aocc,active_mo_.size());
		distmatT dmo2lmo=get_nemo().get_calc() -> localize_PM(world,active_mo_,set);
		tensorT mo2lmo(active_mo_.size(),active_mo_.size());
		dmo2lmo.copy_to_replicated(mo2lmo);
		mo2lmo_ = mo2lmo;

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

		// Prevent misstakes:
		if(shift_>0){MADNESS_EXCEPTION("Potential shift is positive",1);}
		if(not dft_ and shift_ !=0.0){MADNESS_EXCEPTION("Non zero potential shift in TDHF calculation",1);}

		if(only_fock_ and only_GS_){
			print("\nWARNING: only_fock and only_GS demanded ...use both");
			only_fock_ = false;
			only_GS_ = false;
		}

		// Truncate the current mos
		truncate(world,mos_);
		if(world.rank()==0)std::cout << "setup of TDA class ended\n" << std::endl;

		if(compute_virtuals_){
			if(world.rank()==0){
				std::cout << "\nCOMPUTE VIRTUAL ORBITALS\n";
				if(active_mo_.size()!=1) std::cout << "\nWARNING: Active MOs are larger than one, for virtuals only one entry in the excitation vector is needed "
						"-> save time and use the freeze keyword to freeze the rest\n";
			}
		}
		std::cout << "TESTING SECTION\n";
		CCOPS_.test_tda(dft_,nemo_);
		std::cout << "shift is " << shift_ << std::endl;
		std::cout << "ipot is " << ipot_ << std::endl;
	}

	//virtual ~TDA();

	/// Creates and solves guess_xfunctions till pre_convergence is reached
	void solve_guess(xfunctionsT &xfunctions);

	/// Solves the CIS or TDA equations
	void solve(xfunctionsT &xfunctions);

	/// Solves the CIS or TDA equations sequentially for a set of preconverged xfunctions
	void solve_sequential(xfunctionsT &xfunctions);

	/// Returns the MolDFT calulation
//	const SCF get_calc() const {return *calc_;}

	/// Return the Nemo based scf calculations
	const Nemo & get_nemo() const {return nemo_;}

	// Print out grid (e.g for Koala or other external program)
	bool print_grid_TDA() const {return print_grid_;}

	// returns a shallow copy the converged xfunctions
	xfunctionsT get_converged_xfunctions(){return converged_xfunctions_;}

	/// The orbital energies of the SCF calculation
	Tensor<double> const orbital_energies_;

	// returns the orbital energies of the SCF calculation
	Tensor<double>get_orbital_energies()const{
		return orbital_energies_;
	}
	// returns demanded orbital energy of the scf calculation with the right shift (the first unfrozen orbital is number 0)
	double active_eps(const size_t & i)const{
		return orbital_energies_(i+parameters.freeze);
	}
private:

	/// global parameters
	CC_Parameters parameters;

	/// The World
	World & world;

	/// MO to LMO transformation matrix
	Tensor<double> mo2lmo_;

	/// DFT or HF Calculation
	/// for TDA calculations currently only LDA works
	bool dft_;

	/// The SCF calculation of MolDFT
//	const std::shared_ptr<SCF> &calc_;



	/// The SCF calculation using NEMOs (MOs without Nuclear Cusp): MO = R*NEMO, R= Nuclear correlation factor
	//bool use_nemo_;
	const Nemo & nemo_;

	/// The molecular orbitals of moldft
	/// extra member variable that the thresh can be changed without changing calc_
	vecfuncT mos_;

	/// The Operators for the Potential, structure contains also an exchange intermediate
	CIS_Operators CCOPS_;

	/// The molecular orbitals that are used to calcualte the guess excitation vectors
	/// Theese are either the projected numerical mos (projected to minimal AO basis) or just a reference to the numerical mos from moldft
	vecfuncT active_mos_for_guess_calculation_;

	/// Print grid option
	bool print_grid_;

	/// Options for the guess calculation

	/// guess == physical is the only implementation left
	/// new guess functions can be implemented and called in the intialize function
	std::string guess_;
	size_t solve_iter_;
	size_t guess_iter_;
	double guess_omega_;

	/// if guess_mode_ is "numerical" the MOs from moldft will not be projected to the ao basis to form the guess functions
	std::string guess_mode_;

	/// Determine if guess functions should be replaced after pre convergence
	bool replace_guess_functions_;

	/// Excitation operators given in string form
	/// bsp for the excitationoperatr: 1.0*x^2z^3 - 2.0y the string c 1.0 x 2.0 z 3.0 , c -2.0 y 1.0 is needed
	std::vector<std::string> custom_exops_;

	/// Excitation operators for the guess used by big_fock guess
	std::vector<std::string> guess_exops_;

	/// Number of guess excitations to be calculated
	size_t guess_excitations_;
	/// Number of excitations to be caluclated
	size_t excitations_;
	/// Number of parallel iterating excitations
	size_t iterating_excitations_;

	/// Thresholds and convergence cirteria
	//double bsh_eps_;

	/// maximal iterations per guess_function
	size_t iter_max_;

	/// energy convergence level for the guess functions in the solve routine
	double econv_;
	double guess_econv_;
	/// maximal residual for the guess_functions in the solve routine
	double dconv_;
	double guess_dconv_;
	// Convergence criteria (residual norm) for the high thresh sequential iterations in the end
	double hard_dconv_;
	double hard_econv_;
	/// Frozen Orbitals
	//size_t nfreeze_;

	/// Many Plots
	bool plot_;

	/// More output
	//bool parameters.debug;

	/// use only the fock orthonormalization procedure (default)
	bool only_fock_;
	/// use only Gram-Schmidt orthonormalization (not recommended)
	bool only_GS_;

	/// The highest possible excitation to calculate (higher values will result in positive eigenvalues for the BSH operator)
	double highest_excitation_;

	//double lo;

	/// Vector of active molecular orbitals
	vecfuncT active_mo_;

	// the projector on the unperturbed density
	Projector<double,3> rho0;

	/// Active density (closed shell)
	real_function_3d active_density_;

	/// Complete density (includes frozen orbitals), closed shell
	real_function_3d density_;

	/// The potential is calculated when needed and then deleted (saves memory but the potential has to be calculated more often)
	//bool on_the_fly_;

	/// The interface to XCLIB library
	//TDA_DFT xclib_interface_;

	/// Ionization potential for the potential shift used in TDDFT calculations to get bound states for the first excitations (default is -2.0*homo)
	double ipot_;

	/// Kain solver used or not
	bool kain_;

	/// Kain subspace size for the sequential iterations
	size_t kain_subspace_;

	/// The potential shift for the unperturbed DFT potential when using TDDFT (shift = -ipot_ -homo)
	double shift_;

	/// the coulomb potential
	mutable real_function_3d coulomb_;

	/// The converged xfunctions
	std::vector<xfunction> converged_xfunctions_;

	/// Calculate triplets
	bool triplet_;

	/// Use the excitation energy in the BSH operator (if not it is added to the potential)
	bool use_omega_for_bsh_;

	/// Compute virtual orbitals, the freeze_ key should then be set to (all_mos)-1
	bool compute_virtuals_;

	/// The thresholds for the guess, solve and sequential calculation
	//double guess_thresh_;
	//double solve_thresh_;
	//double solve_sequential_thresh_;

	/// Print the current xfunctions in a formated way
	/// @param[in] xfunctions a vector of xfunction structures
	void print_status(const xfunctionsT & xfunctions)const;

	/// just a helper function for print_status and others
	/// @param[in] x a single xfunction structure
	/// the function will print out the information of the xfunction structure (energy, iterations, convergence ...) in a formated way
	void print_xfunction(const xfunction &x)const;

	/// Takes an empty vector of excitation functions and passes it to one of the guess functions
	/// @param[in] xfunctions empty vector of xfunctions (no necessarily empty)
	void initialize(xfunctionsT & xfunctions);

	void make_big_fock_guess(xfunctionsT &xfunctions)const;

	/// Creates physical guess functions (x,y,z excitations - depending on the input file, see make_excitation_operators function)
	void guess_physical(xfunctionsT & xfunctions)const;

	/// guess_ao_excitation
	void guess_custom_2(xfunctionsT &xfunctions)const;

    /// guess: localize MOs and excite with a dipole
    void guess_local(xfunctionsT &xfunctions)const;

	/// Make a huge guess: Excite on every non hydrogen atom
	void guess_atomic_excitation(xfunctionsT & xfunctions)const;

	vecfuncT make_guess_vector(const std::string &input)const;

	void guess_koala(xfunctionsT &roots)const;

	/// Create excitation operators (e.g x,y,z for dipole excitations bzw symmetry operators)
	/// @return gives back a vectorfunction of excitation operators (specified in the input file)
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
	/// @param[in] xfunctions all excitation structures
	/// @param[in] guess for the first iterations (no energy update, no kain update)
	void iterate_all(xfunctionsT &xfunctions,bool guess);

	/// Applies the greens operator and calcualtes the updated xfunction for one xfunction
	/// @param[in] xfunction a single xfunction structure (contains the response orbitals as vecfunc x)
	/// @return The updated xfunction
	vecfuncT iterate_one(xfunction & xfunction)const;

	/// Update energies (decide if second order or expectation value should be used)
	/// @return the update method (2nd order, expectation value, setback)
	std::string update_energy(xfunction &xfunction)const;

	/// Normalize one or all excitation functions
	void normalize(xfunctionsT &xfunctions)const;
	void normalize(xfunction &xfunction)const;

	/// Project out the converged xfunctions
	void project_out_converged_xfunctions(xfunctionsT & xfunctions)const;

	/// Orthonormalize the xfunction with the perturbed Fock Matrix
	// 1. Call make_perturbed_fock_matrix(xfunctions)
	// 2. Diagonalize
	// 3. Update Energy and xfunctions
	/// @param[in] xfunctions the xfunctions
	/// @return true is fock matrix was calculated (if not that means no energy was calculated and that the expectation value needs to be calculated in the iterate_one procedure)
	bool orthonormalize_fock(xfunctionsT &xfunctions)const;

	/// a little helper routine to measure the degree of offdiagonality in a 2d tensor
	double measure_offdiagonality(const madness::Tensor<double> &U,const size_t size)const;

	std::vector<vecfuncT> transform_vecfunctions(const std::vector<vecfuncT> &xfunctions,const madness::Tensor<double> U)const;

	/// Projects out the occupied space
	// 1. Make projector
	// 2. apply
	void project_out_occupied_space(vecfuncT &x)const;

	/// Calculate offdiagonal elements of the perturbed fock matrix
	double perturbed_fock_matrix_element(const vecfuncT &xr, const vecfuncT &Vxp,const vecfuncT &xp)const;

	/// Calculate the expectation value and update xfunction.expectation_value
	// Can also be used to calculate diagonal elements of the fock matrix
	double expectation_value(const xfunction &x,const vecfuncT &Vx)const;

	/// Make the perturbed Fock Matrix
	// 1. Calculate potential (V0 + Gamma) --> Call perturbed_potential(exfunctions)
	// 2. Calculate Matrix Elements
	Tensor<double> make_perturbed_fock_matrix(const xfunctionsT &xfunctions)const;

	/// The CIS or TDA Potential without the nuclear potential is applied to one xfunction
	/// @param[in] xfunction one
	vecfuncT apply_smooth_potential(const xfunction&xfunction)const;

	vecfuncT apply_nuclear_potential(const xfunction &xfunction) const;

	/// Plot vectorfunction (for convenience)
	/// @param[in] x the vectorfunction to plot
	/// @param[in] msg the name for the vecfunction plot
	/// @param[in] plot if false nothing is done
	void plot_vecfunction(const vecfuncT &x,std::string msg, bool plot = true)const;

	/// Check convergence
	/// checks if the xfunctions have converged
	/// @param[in] xfunctions a vector of xfunction structures
	/// @param[in] guess : decide if guess criteria for convergece should be used
	void check_convergence(xfunctionsT &xfunctions, const bool guess)const ;

	/// Print performance: Values of expectation values and errors of each iteration into a file
	/// @param[in] xfunctions a vector of xfunction structures
	/// @param[in] prename the saved file will be prename+results.tex
	void print_performance(const xfunctionsT &xfunctions,const std::string prename)const;

	/// Truncate the xfunctions structure:
	/// @param[in] xfunctions a vector of xfunction structures
	/// Truncates the vecfunctions: x, Vx and the current_residual (if not empty)
	void truncate_xfunctions(xfunctionsT &xfunctions)const;

	/// load a converged root from disk

	/// compute the oscillator strength in the length resentation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	xfunction	a converged root
	double oscillator_strength_length(const xfunction& xfunction) const;

	/// compute the oscillator strength in the velocity representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_velocity(const xfunction& root) const;

	void memory_information(const xfunctionsT &xfunctions)const;
public:
	/// get the threshholds
	//double get_guess_thresh(){return guess_thresh_;}
	//double get_solve_thresh(){return solve_thresh_;}
	//double get_solve_sequential_thresh(){return solve_sequential_thresh_;}
	/// analyze the root: oscillator strength and contributions from occ
	void analyze(xfunctionsT& roots) const;
	/// Project a vecfuncT to the ao basis (used to create projected MOs for the guess calculation)
	vecfuncT project_to_ao_basis(const vecfuncT & mos, const vecfuncT& ao_basis)const;

	void output_section(const std::string &msg)const{
		if(world.rank()==0){
			std::cout << "\n\n------------------------------------------------\n";
			std::cout << msg;
			std::cout << "\n------------------------------------------------\n\n";
		}
	}

	void output(const std::string &msg)const{
		if(world.rank()==0){
			std::cout << msg << std::endl;
		}
	}

	void sanitycheck(const xfunctionsT &xfunctions)const{
		for(auto xf:xfunctions){
			for(auto x:xf.x){
				if(x.thresh()!=FunctionDefaults<3>::get_thresh()) MADNESS_EXCEPTION("ERROR: WRONG THRESH IN XFUNCTIONS DETECTED",1);
			}
		}
	}

};

} /* namespace madness */

#endif /* TDA_H_ */
