/*
 * CCOperators.h
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */

////// TODO LIST
// Make two_electron integrals calculation parallel
// S5c diagram: Test alternative with perturbed density (not 3 for loops anymore) -> see time difference (but keep in mind that integral intermediate has to be calculated anyway for the exchange part
// Or do S5c and S5c_X together
#ifndef CCOPERATORS_H_
#define CCOPERATORS_H_

// Operators for coupled cluster and CIS
#include <chem/SCFOperators.h>
#include <chem/electronic_correlation_factor.h>

#include <algorithm> // tolower function for strings
//#include <string>

// to debug
#include<chem/mp2.h>

namespace madness {

// Timer Structure
struct CC_Timer{
	/// TDA_TIMER contructor
	/// @param[in] world the world
	/// @param[in] msg	a string that contains the desired printout when info function is called
	CC_Timer(World &world,std::string msg) : world(world),start_wall(wall_time()),start_cpu(cpu_time()),operation(msg) {}
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
			if(world.rank()==0) std::cout<< std::setw(40) << std::setfill(' ') << "Timer: " << operation << " : " << std::scientific << std::setprecision(1)
			<< end_wall << "s (wall) "<< end_cpu << "s (cpu)" << std::endl;
		}
	}

};

struct CC_Parameters{
	// default constructor
	CC_Parameters():
	lo(FunctionDefaults<3>::get_thresh()),
	thresh_3D(FunctionDefaults<3>::get_thresh()),
	thresh_6D(FunctionDefaults<6>::get_thresh()),
	thresh_6D_tight(thresh_6D*0.1),
	thresh_bsh_3D(FunctionDefaults<3>::get_thresh()),
	thresh_bsh_6D(FunctionDefaults<6>::get_thresh()),
	thresh_poisson_3D(FunctionDefaults<3>::get_thresh()),
	thresh_poisson_6D(FunctionDefaults<6>::get_thresh()),
	thresh_f12(FunctionDefaults<6>::get_thresh()),
	thresh_Ue(FunctionDefaults<6>::get_thresh()),
	econv(1.e-4),
	dconv_3D(1.e-2),
	dconv_6D(1.e-2),
	nfreeze(0),
	iter_max_3D(30),
	iter_max_6D(30),
	restart(false),
	corrfac_gamma(2.0),
	output_prec(8),
	debug(false)
	{}

	// read parameters from input
	/// ctor reading out the input file
	CC_Parameters(const std::string& input,const double &low) :
		lo(1.e-6),
		thresh_3D(FunctionDefaults<3>::get_thresh()),
		thresh_6D(FunctionDefaults<6>::get_thresh()),
		thresh_6D_tight(thresh_6D*0.1),
		thresh_bsh_3D(FunctionDefaults<3>::get_thresh()*0.1),
		thresh_bsh_6D(FunctionDefaults<6>::get_thresh()*0.1),
		thresh_poisson_3D(FunctionDefaults<3>::get_thresh()*0.1),
		thresh_poisson_6D(FunctionDefaults<6>::get_thresh()*0.1),
		thresh_f12(FunctionDefaults<6>::get_thresh()*0.1),
		econv(1.e-4),
		dconv_3D(1.e-2),
		dconv_6D(1.e-2),
		nfreeze(0),
		iter_max_3D(30),
		iter_max_6D(30),
		restart(false),
		corrfac_gamma(2.0),
		output_prec(8),
		debug(false),
		kain(false),
		kain_subspace(3)
	{
		// get the parameters from the input file
        std::ifstream f(input.c_str());
        position_stream(f, "cc2");
        std::string s;

        // minimum operator thresh
        double minopthresh = 1.e-4;

        while (f >> s) {
        	std::transform(s.begin(),s.end(),s.begin(), ::tolower);
            if (s == "end") break;
            else if (s == "debug") debug=true;
            else if (s == "lo") f >> lo;
            else if (s == "econv"){
            	f >> econv;
            }

            else if (s == "dconv"){
            	double tmp = 0.0;
            	f >> tmp;
            	dconv_3D = tmp; dconv_6D=tmp;
            }
            else if (s == "dconv_3d"){
            	f >> dconv_3D;
            }
            else if (s == "dconv_6d"){
            	f >> dconv_6D;
            }
            else if (s == "thresh"){
            	double tmp = 0.0;
            	f >> tmp;
            	double opthresh = tmp*0.1;
            	double opthresh_3D = tmp*0.01;
            	if(opthresh > minopthresh) opthresh = minopthresh;
            	thresh_3D 		  = 0.01*tmp;
            	thresh_6D 		  = tmp;
            	thresh_poisson_3D = opthresh_3D;
            	thresh_poisson_6D = opthresh;
            	thresh_bsh_3D     = opthresh_3D;
            	thresh_bsh_6D     = opthresh;
            	thresh_f12        = opthresh;
            	thresh_Ue		  = opthresh;
            }
            else if (s == "thresh_operators" or s == "thresh_operator"){
            	double tmp =0.0;
            	f >> tmp;
            	if(tmp>minopthresh) tmp = minopthresh;
            	thresh_poisson_3D = tmp;
            	thresh_poisson_6D = tmp;
            	thresh_bsh_3D     = tmp;
            	thresh_bsh_6D     = tmp;
            	thresh_f12        = tmp;
            	thresh_Ue		  = tmp;
            }
            else if (s == "thresh_operators_3d" or s == "thresh_operator_3d"){
            	double tmp =0.0;
            	f >> tmp;
            	if(tmp>minopthresh) tmp = minopthresh;
            	thresh_poisson_3D = tmp;
            	thresh_bsh_3D     = tmp;
            }
            else if (s == "thresh_operators_6d" or s == "thresh_operator_3d"){
            	double tmp =0.0;
            	f >> tmp;
            	if(tmp>minopthresh) tmp = minopthresh;
            	thresh_poisson_6D = tmp;
            	thresh_bsh_6D     = tmp;
            	thresh_f12        = tmp;
            	thresh_Ue		  = tmp;
            }
            else if (s == "thresh_3d") f >> thresh_3D;
            else if (s == "thresh_6d") f >> thresh_6D;
            else if (s == "thresh_bsh_3d") f >> thresh_bsh_3D;
            else if (s == "thresh_bsh_6d") f >> thresh_bsh_6D;
            else if (s == "thresh_poisson_3d") f >> thresh_poisson_3D;
            else if (s == "thresh_poisson_6d") f >> thresh_poisson_6D;
            else if (s == "thresh_f12") f >> thresh_f12;
            else if (s == "thresh_Ue") f >> thresh_Ue;
            else if (s == "freeze" or s=="nfreeze") f >> nfreeze;
            else if (s == "iter_max"){
            	f >> iter_max_3D;
            	iter_max_6D = iter_max_3D;
            }
            else if (s == "iter_max_3d") f >> iter_max_3D;
            else if (s == "iter_max_6d") f >> iter_max_6D;
            else if (s == "restart") restart=true;
            else if (s == "corrfac_gamma" or "gamma") f>>corrfac_gamma;
            else if (s == "kain") kain=true;
            else if (s == "kain_subspace") f>>kain_subspace;
            else continue;
        }

        thresh_6D_tight = thresh_6D*0.1;

        if(not kain) kain_subspace = 0;

        // set the thresholds
        FunctionDefaults<3>::set_thresh(thresh_3D);
        FunctionDefaults<6>::set_thresh(thresh_6D);
        if(econv < 1.e-1) output_prec = 2;
        if(econv < 1.e-2) output_prec = 3;
        if(econv < 1.e-3) output_prec = 4;
        if(econv < 1.e-4) output_prec = 5;
        if(econv < 1.e-5) output_prec = 6;
        if(econv < 1.e-6) output_prec = 7;
	}

	double lo;
	// function thresh 3D
	double thresh_3D;
	// function thresh 6D
	double thresh_6D;
	// tight thresh for 6D functions (for addition)
	double thresh_6D_tight;
	// BSH thresh
	double thresh_bsh_3D;
	double thresh_bsh_6D;
	// Poisson thresh
	double thresh_poisson_3D;
	double thresh_poisson_6D;
	// f12 thresh
	double thresh_f12;
	// Ue thresh
	double thresh_Ue;
	// Convergence for Correlation Energy
	double econv;
	// Convergence for CC-singles
	double dconv_3D;
	// Convergence for CC-Doubles
	double dconv_6D;
	// Number of frozen Orbitals
	size_t nfreeze;
	// iterations
	size_t iter_max_3D;
	size_t iter_max_6D;
	// restart
	bool restart;
	// Exponent for the correlation factor
	double corrfac_gamma;
	// for formated output
	size_t output_prec;
	// debug mode
	bool debug;
	// use kain or not
	bool kain;
	size_t kain_subspace;

	// print out the parameters
	void information(World &world)const{
		if(world.rank()==0){
			//std::cout << "CC2 Parameters:\n";
			std::cout << std::setw(20) << std::setfill(' ') << "restart :"           << restart << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "lo :"                << lo << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "k (3D) :"                << FunctionDefaults<3>::get_k() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "k (6D) :"                << FunctionDefaults<6>::get_k() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_3D demanded :"         << thresh_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_3D set :"         << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_6D demanded :"         << thresh_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_6D set :"         << FunctionDefaults<6>::get_thresh() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_bsh_3D :"     << thresh_bsh_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_bsh_6D :"     << thresh_bsh_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_poisson_3D :" << thresh_poisson_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_poisson_6D :" << thresh_poisson_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_f12 :"        << thresh_f12 << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_Ue :"        << thresh_Ue << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "econv :"             << econv << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "dconv_3D :"          << dconv_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "dconv_6D :"          << dconv_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "nfreeze :"           << nfreeze << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "iter_max_3D :"           << iter_max_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "iter_max_6D :"           << iter_max_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "truncation mode 3D :" << FunctionDefaults<3>::get_truncate_mode()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "truncation mode 6D :" << FunctionDefaults<6>::get_truncate_mode()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "tensor type: " << FunctionDefaults<6>::get_tensor_type()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "facReduce:" << GenTensor<double>::fac_reduce()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "max. displacement:" << Displacements<6>::bmax_default()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "apply randomize:" << FunctionDefaults<6>::get_apply_randomize()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "Cell min width (3D, 6D) :" << FunctionDefaults<6>::get_cell_min_width() << ", " << FunctionDefaults<3>::get_cell_min_width()  <<std::endl;
			//std::cout << std::setw(20) << std::setfill(' ') << "Cell widths (3D) :" << FunctionDefaults<3>::get_cell_width()  <<std::endl;
			//std::cout << std::setw(20) << std::setfill(' ') << "Cell widths (6D) :" << FunctionDefaults<6>::get_cell_width()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "Autorefine (3D, 6D) :" << FunctionDefaults<6>::get_autorefine() << ", " << FunctionDefaults<3>::get_autorefine()  <<std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "debug mode is: " << debug  <<std::endl;
			if(kain) std::cout << std::setw(20) << std::setfill(' ') << "Kain is used with subspace " << kain_subspace <<std::endl;
			else std::cout << std::setw(20) << std::setfill(' ') << "Kain is not used" << debug  <<std::endl;

		}
	}

	void sanity_check(World &world)const{
		size_t warnings = 0;
		if(FunctionDefaults<3>::get_thresh() > 0.01*FunctionDefaults<6>::get_thresh()) warning(world,"3D Thresh is too low, should be 0.01*6D_thresh");
		if(FunctionDefaults<3>::get_thresh() > 0.1*FunctionDefaults<6>::get_thresh()) error(world,"3D Thresh is way too low, should be 0.01*6D_thresh");
		if(FunctionDefaults<3>::get_cell_min_width() != FunctionDefaults<6>::get_cell_min_width()) error(world,"3D and 6D Cell sizes differ");
		if(FunctionDefaults<3>::get_k() != FunctionDefaults<6>::get_k()) error(world, "k-values of 3D and 6D differ ");
		if(FunctionDefaults<3>::get_truncate_mode()!=3) warnings+=warning(world,"3D Truncate mode is not 3");
		if(FunctionDefaults<6>::get_truncate_mode()!=3) warnings+=warning(world,"6D Truncate mode is not 3");
		if(dconv_3D < FunctionDefaults<3>::get_thresh()) error(world,"Demanded higher convergence than threshold for 3D");
		if(dconv_6D < FunctionDefaults<6>::get_thresh()) error(world,"Demanded higher convergence than threshold for 6D");
		if(thresh_3D != FunctionDefaults<3>::get_thresh()) error(world,"3D thresh set unequal 3D thresh demanded");
		if(thresh_6D != FunctionDefaults<6>::get_thresh()) error(world,"6D thresh set unequal 6D thresh demanded");
		if(econv < FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 3D");
		if(econv < FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 6D");
		if(econv < 0.1*FunctionDefaults<3>::get_thresh()) warning(world,"Demanded higher energy convergence than threshold for 3D (more than factor 10 difference)");
		if(econv < 0.1*FunctionDefaults<6>::get_thresh()) warning(world,"Demanded higher energy convergence than threshold for 6D (more than factor 10 difference)");
		// Check if the 6D thresholds are not too high
		if(thresh_6D < 1.e-3) warnings+=warning(world,"thresh_6D is smaller than 1.e-3");
		if(thresh_Ue < 1.e-4) warnings+=warning(world,"thresh_Ue is smaller than 1.e-4");
		if(thresh_3D > 0.01*thresh_6D) warnings+=warning(world,"Demanded 6D thresh is to precise compared with the 3D thresh");
		if(thresh_3D > 0.1*thresh_6D) error(world,"Demanded 6D thresh is to precise compared with the 3D thresh");
		if(kain and kain_subspace !=0) warnings+=warning(world,"Demanded Kain solver but the size of the iterative subspace is set to zero");
		if(warnings >0){
			if(world.rank()==0) std::cout << warnings <<"Warnings in parameters sanity check!\n\n";
		}else{
			if(world.rank()==0) std::cout << "Sanity check for parameters passed\n\n" << std::endl;
		}
	}

	void error(World& world,const std::string &msg)const{
		if(world.rank()==0) std::cout << "\n\n\n\n\n!!!!!!!!!\n\nERROR IN CC_PARAMETERS:\n    ERROR MESSAGE IS: " << msg << "\n\n\n!!!!!!!!" << std::endl;
		MADNESS_EXCEPTION("ERROR IN CC_PARAMETERS",1);
	}
	size_t warning(World& world,const std::string &msg)const{
		if(world.rank()==0) std::cout << "WARNING IN CC_PARAMETERS!: " << msg << std::endl;
		return 1;
	}
};

// TAKEN FROM MP2.h
/// POD holding all electron pairs with easy access
template<typename T>
struct Pairs {

	typedef std::map<std::pair<int, int>, T> pairmapT;
	pairmapT allpairs;

	/// getter
	const T& operator()(int i, int j) const {
		return allpairs.find(std::make_pair(i, j))->second;
	}

	/// getter
	T& operator()(int i, int j) {
		return allpairs[std::make_pair(i, j)];
	}

	/// setter
	void insert(int i, int j, T pair) {
		std::pair<int, int> key = std::make_pair(i, j);
		allpairs.insert(std::make_pair(key, pair));
	}
};

/// enhanced POD for the pair functions
class CC_Pair: public archive::ParallelSerializableObject {

public:
	/// default ctor; initialize energies with a large number
	CC_Pair() :
			i(-1), j(-1), e_singlet(uninitialized()), e_triplet(
					uninitialized()), ij_gQf_ij(uninitialized()), ji_gQf_ij(
					uninitialized()), iteration(0), converged(false) {
	}

	/// ctor; initialize energies with a large number
	CC_Pair(const int i, const int j) :
			i(i), j(j), e_singlet(uninitialized()), e_triplet(uninitialized()), ij_gQf_ij(
					uninitialized()), ji_gQf_ij(uninitialized()), iteration(0), converged(
					false) {
	}

	/// print the pair's energy
	void print_energy() const {
		if (function.world().rank() == 0) {
			printf("final correlation energy %2d %2d %12.8f %12.8f\n", i, j,
					e_singlet, e_triplet);
		}
	}

	// print information
	void info()const{
		if(function.world().rank()==0){
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " Current Information about Electron Pair |u" << i << j << ">"  << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " corelation energy: " << e_singlet + e_triplet << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " e_singlet: " << e_singlet << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " e_triplet: " << e_triplet << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ij_gQf_ij: " << ij_gQf_ij << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ji_gQf_ij: " << ji_gQf_ij << std::endl;
		}
	}

	static double uninitialized() {
		return 1.e10;
	}

	int i, j;                       ///< orbitals i and j
	real_function_6d function; ///< pair function for a specific pair w/o correlation factor part
	real_function_6d constant_term;	///< the first order contribution to the MP1 wave function

	double e_singlet;				///< the energy of the singlet pair ij
	double e_triplet;				///< the energy of the triplet pair ij

	double ij_gQf_ij;         	  	///< <ij | g12 Q12 f12 | ij>
	double ji_gQf_ij;          		///< <ji | g12 Q12 f12 | ij>

	int iteration;					///< current iteration for restart
	bool converged;					///< is the pair function converged

	/// serialize this CC_Pair

	/// store the function only if it has been initialized
	/// load the function only if there is one
	/// don't serialize recomputable intermediates r12phi, Uphi, KffKphi
	template<typename Archive> void serialize(Archive& ar) {
		bool fexist = function.is_initialized();
		bool cexist = constant_term.is_initialized();
		ar & ij_gQf_ij & ji_gQf_ij & e_singlet & e_triplet & converged
				& iteration & fexist & cexist;
		if (fexist)
			ar & function;
		if (cexist)
			ar & constant_term;
	}

	bool load_pair(World& world) {
		std::string name = "pair_" + stringify(i) + stringify(j);
		bool exists = archive::ParallelInputArchive::exists(world,
				name.c_str());
		if (exists) {
			if (world.rank() == 0)
				printf("loading matrix elements %s", name.c_str());
			archive::ParallelInputArchive ar(world, name.c_str(), 1);
			ar & *this;
			if (world.rank() == 0)
				printf(" %s\n", (converged) ? " converged" : " not converged");
			function.set_thresh(FunctionDefaults<6>::get_thresh());
			constant_term.set_thresh(FunctionDefaults<6>::get_thresh());
		} else {
			if (world.rank() == 0)
				print("could not find pair ", i, j, " on disk");
		}
		return exists;
	}

	void store_pair(World& world) {
		std::string name = "pair_" + stringify(i) + stringify(j);
		if (world.rank() == 0)
			printf("storing matrix elements %s\n", name.c_str());
		archive::ParallelOutputArchive ar(world, name.c_str(), 1);
		ar & *this;
	}
};

class ElectronSingles{
public:
	// Constructors
	ElectronSingles() : converged(false),iterations(0){}
	ElectronSingles(const size_t ii): i(ii), converged(false),iterations(0) {}
	ElectronSingles(const size_t i, const real_function_3d &f): i(i), converged(false),iterations(0), function_(f){}
	ElectronSingles(const ElectronSingles &other): i(other.i), converged(other.converged),iterations(other.iterations), function_(other.function()){}

	real_function_3d function()const{return function_;}

	size_t i;
	bool converged;
	mutable size_t iterations;
private:
	real_function_3d function_;
	static size_t uninitialized(){return 999;}


};

static double unitfunction(const coord_3d &r) {
	return 1.0;
}

// forward declaration
//class SCF;
//class Nemo;
//class NuclearCorrelationFactor;
//class XCfunctional;
//class Nuclear;
typedef std::vector<Function<double, 3> > vecfuncT;

/// Structure that holds the CC intermediates and is able to refresh them
struct CC_Intermediates {
public:
	CC_Intermediates(World&world, const vecfuncT &bra, const vecfuncT &ket,
			const Nemo&nemo, const CC_Parameters &param) :
			world(world), parameters(param), mo_bra_(bra), mo_ket_(ket), poisson(
					std::shared_ptr < real_convolution_3d
							> (CoulombOperatorPtr(world,
									parameters.lo,
									parameters.thresh_poisson_3D))), density_(
					make_density(bra, ket)), exchange_intermediate_(
					make_exchange_intermediate(bra, ket)), hartree_potential_(
							make_hartree_potential(density_)), integrals_hf_(
					make_two_electron_integrals_hf()) {

	}

	/// Get the intermediates
	real_function_3d get_density() {
		return density_;
	}
	real_function_3d get_perturbed_density() const {
		return perturbed_density_;
	}
	real_function_3d get_hartree_potential() const {
		return hartree_potential_;
	}
	real_function_3d get_J() {
		return hartree_potential_;
	}
	real_function_3d get_perturbed_hartree_potential() const {
		return perturbed_hartree_potential_;
	}
	real_function_3d get_pJ() {
		return perturbed_hartree_potential_;
	}
	std::vector<vecfuncT> get_exchange_intermediate() const {
		return exchange_intermediate_;
	}
	std::vector<vecfuncT> get_EX() {
		return exchange_intermediate_;
	}
	std::vector<vecfuncT> get_perturbed_exchange_intermediate() const {
		return perturbed_exchange_intermediate_;
	}
	std::vector<vecfuncT> get_pEX() const {
		return perturbed_exchange_intermediate_;
	}
	Tensor<double> get_intergrals_hf() const {
		return integrals_hf_;
	}
	Tensor<double> get_integrals_mixed_t1() const {
		return integrals_mixed_t1_;
	}
	Tensor<double> get_integrals_t1() const {
		return integrals_t1_;
	}
	/// refresh the intermediates that depend on the \tau functions
	void update(const vecfuncT &tau) {
		if (world.rank() == 0)
			std::cout << "Update Intermediates:\n";
		perturbed_density_ = make_density(mo_bra_, tau);
		perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
		perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,
				tau);
		integrals_mixed_t1_ = make_two_electron_integrals_mixed_t1(tau);
		integrals_t1_ = make_two_electron_integrals_t1(tau);
	}

	/// make a density from two input functions
	/// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
	/// @param[in] vecfuncT_bra
	/// @param[in] vecfuncT_ket
	/// @param[out] \sum_i bra_i * ket_i
	real_function_3d make_density(const vecfuncT &bra,
			const vecfuncT &ket) const {
		if (bra.size() != ket.size())
			error(
					"error in make density: unequal sizes ("
							+ stringify(bra.size()) + " and "
							+ stringify(ket.size()) + ")");
		if (bra.empty())
			error("error in make_density: bra_element is empty");
		// make the density
		real_function_3d density = real_factory_3d(world);
		for (size_t i = 0; i < bra.size(); i++)
			density += bra[i] * ket[i];
		density.truncate();
		return density;
	}
	/// Poisson operator
	std::shared_ptr<real_convolution_3d> get_poisson()const{return poisson;}

private:
	World &world;
	const CC_Parameters &parameters;
	const vecfuncT &mo_bra_;
	const vecfuncT &mo_ket_;
	const std::shared_ptr<real_convolution_3d> poisson;
	/// const intermediates
	const real_function_3d density_;
	/// Exchange intermediate: EX(i,j) = <i|g|j>
	const std::vector<vecfuncT> exchange_intermediate_;
	/// Hartree_Potential  = J = \sum_k <k|g|k> = Poisson(density)
	const real_function_3d hartree_potential_;
	/// intermediates that need to be recalculated before every iteration
	/// Perturbed Density = \sum_k |k><\tau_k|
	real_function_3d perturbed_density_;
	/// Perturbed Hartree Poptential PJ = \sum_k <k|g|\tau_k> = Poisson(perturbed_density)
	real_function_3d perturbed_hartree_potential_;
	/// Perturbed Exchange Intermediate: PEX(i,j) = <i|g|\tau_j>
	std::vector<vecfuncT> perturbed_exchange_intermediate_;

	/// Two electron integrals
	/// The Integrals which consist of the hartree-fock ground state
	/// <ij|g|kl> = <ji|g|lk>
	const Tensor<double> integrals_hf_;
	/// The Integrals which consist of ground state and t1 amplitudes
	/// <ij|g|k\tau_l> = <ji|g|\tau_lk>
	Tensor<double> integrals_mixed_t1_;
	/// The Integrals from the t1 functions and the hf orbitals
	/// <ij|g|\tau_k\tau_l> = <ji|g|\tau_l\tau_k>
	Tensor<double> integrals_t1_;

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_INTERMEDIATES!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
				1);
	}

public:
	/// Make the exchange intermediate: EX[j][i] <bra[i](r2)|1/r12|ket[j](r2)>
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT &bra,
			const vecfuncT &ket) const {
		if (bra.size() != ket.size() or bra.empty())
			error(
					"in make_exchange_intermediate, bra and ket empty or unequal sizes:\n bra_size: "
							+ stringify(bra.size()) + ", ket_size: "
							+ stringify(ket.size()));
		std::vector<vecfuncT> EX;
		EX.resize(bra.size());
		for (size_t i = 0; i < bra.size(); i++) {
			EX[i].resize(ket.size());
			for (size_t j = 0; j < ket.size(); j++) {
				EX[i][j] = (*poisson)(bra[j] * ket[i]);
			}
			truncate(world, EX[i]);
		}
		return EX;
	}
	/// Calculates the hartree potential Poisson(density)
	/// @param[in] density: a 3d function on which the poisson operator is applied (can be the occupied density and the perturbed density)
	/// @param[out] poisson(density) = \int 1/r12 density(r2) dr2
	real_function_3d make_hartree_potential(
			const real_function_3d &density) const {
		real_function_3d hartree = (*poisson)(density);
		hartree.truncate();
		return hartree;
	}

	/// Calculates two electron integrals
	/// <ij|g|kl>
	Tensor<double> make_two_electron_integrals_hf() const {
		Tensor<double> result(mo_bra_.size(), mo_bra_.size(), mo_ket_.size(),
				mo_ket_.size());
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			for (size_t j = 0; j < mo_bra_.size(); j++) {
				for (size_t k = 0; k < mo_ket_.size(); k++) {
					for (size_t l = 0; l < mo_ket_.size(); l++) {
						result(i, j, k, l) = (mo_bra_[i] * mo_ket_[k]).inner(
								exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}
	/// <ij|g|k\tau_l>
	Tensor<double> make_two_electron_integrals_mixed_t1(
			const vecfuncT &tau) const {
		Tensor<double> result(mo_bra_.size(), mo_bra_.size(), mo_ket_.size(),
				tau.size());
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			for (size_t j = 0; j < mo_bra_.size(); j++) {
				for (size_t k = 0; k < mo_ket_.size(); k++) {
					for (size_t l = 0; l < tau.size(); l++) {
						result(i, j, k, l) = (mo_bra_[i] * mo_ket_[k]).inner(
								perturbed_exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}
	// <ij|g|\tau_k \tau_l>
	Tensor<double> make_two_electron_integrals_t1(const vecfuncT &tau) const {
		Tensor<double> result(mo_bra_.size(), mo_bra_.size(), tau.size(),
				tau.size());
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			for (size_t j = 0; j < mo_bra_.size(); j++) {
				for (size_t k = 0; k < tau.size(); k++) {
					for (size_t l = 0; l < tau.size(); l++) {
						result(i, j, k, l) = (mo_bra_[i] * tau[k]).inner(
								perturbed_exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}

	bool use_timer_;
	mutable double ttt, sss;
	void START_TIMER() const {
		if (use_timer_)
			world.gop.fence();
		ttt = wall_time();
		sss = cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if (use_timer_)
			END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if (use_timer_) {
			ttt = wall_time() - ttt;
			sss = cpu_time() - sss;
			if (world.rank() == 0)
				printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}
};

/// Coupled Cluster Operators (all closed shell)
class CC_Operators {
public:
	/// Constructor
	CC_Operators(World& world, const Nemo &nemo,
			const CorrelationFactor &correlationfactor, const CC_Parameters &param) : Q12(world),
			world(world), nemo(nemo), corrfac(correlationfactor),parameters(param), mo_bra_(
					make_mo_bra(nemo)), mo_ket_(nemo.get_calc()->amo), intermediates_(
					world, mo_bra_, mo_ket_, nemo, param), use_timer_(true) {

		// initialize the Q12 projector
		Q12.set_spaces(mo_bra_,mo_ket_,mo_bra_,mo_ket_);
	}

	StrongOrthogonalityProjector<double,3> Q12;

	vecfuncT get_CIS_potential(const vecfuncT &tau) {
		START_TIMER();
		intermediates_.update(tau);
		END_TIMER("update intermediates");
		vecfuncT result = add(world, S3c(tau), S3c_X(tau));
		Q(result);
		return add(world, result, fock_residue_closed_shell(tau));
	}

	vecfuncT get_CC2_singles_potential(const vecfuncT &singles, const Pairs<real_function_6d> &doubles)const{
		START_TIMER();
		vecfuncT result = fock_residue_closed_shell(singles);
		END_TIMER("Singles Potential: Fock Residue");
		START_TIMER();
		result =add(world, S3c(singles),result);
		END_TIMER("Singles Potential: S3c");
		START_TIMER();
		result =add(world, S3c_X(singles),result);
		END_TIMER("Singles Potential: S3cX");
		START_TIMER();
		result =add(world, S5b(singles),result);
		END_TIMER("Singles Potential: S5b");
		START_TIMER();
		result =add(world, S5b_X(singles),result);
		END_TIMER("Singles Potential: S5bX");
		START_TIMER();
		result =add(world, S5c(singles),result);
		END_TIMER("Singles Potential: S5c");
		START_TIMER();
		result =add(world, S5c_X(singles),result);
		END_TIMER("Singles Potential: S5cX");
		START_TIMER();
		result =add(world, S6(singles),result);
		END_TIMER("Singles Potential: S6");
		START_TIMER();
		result =add(world, S2b(doubles),result);
		END_TIMER("Singles Potential: S2b+X");
		START_TIMER();
		result =add(world, S2c(doubles),result);
		END_TIMER("Singles Potential: S2c+X");
		START_TIMER();
		result =add(world, S4a(doubles,singles),result);
		END_TIMER("Singles Potential: S4a+X");
		START_TIMER();
		START_TIMER();
		result =add(world, S4b(doubles,singles),result);
		END_TIMER("Singles Potential: S4b+X");
		START_TIMER();
		result =add(world, S4c(doubles,singles),result);
		END_TIMER("Singles Potential: S4c+X");
		Q(result);
		truncate(world,result);
		return result;
	}

	real_function_6d get_CC2_doubles_potential(const vecfuncT &singles, const CC_Pair &u)const{
	MADNESS_EXCEPTION("CC2 doubles potential not yet implemented",1);
	return real_factory_6d(world);
	}

	real_function_6d get_MP2_potential_constant_part(CC_Pair &u)const{
		CC_Timer timer_U(world,"Ue(R)|ij>");
		real_function_6d UePart = apply_transformed_Ue(mo_ket_[u.i],mo_ket_[u.j],u.i,u.j,u);
		UePart.print_size("Ue|"+stringify(u.i)+stringify(u.j)+">");
		timer_U.info();

		CC_Timer timer_KffK(world,"Kf|ij>");
		real_function_6d KffKPart = apply_exchange_commutator(mo_ket_[u.i],mo_ket_[u.j],"occupied",u.i,u.j);
		KffKPart.print_size("[K,f]|"+stringify(u.i)+stringify(u.j)+">");
		timer_KffK.info();

		real_function_6d unprojected_result = (UePart - KffKPart).truncate();
		unprojected_result.print_size("Ue - [K,f]|"+stringify(u.i)+stringify(u.j)+">");
		CC_Timer timer_Q(world,"Apply Q12");
		real_function_6d result = Q12(unprojected_result);
		result.print_size("Q12(Ue - [K,f]|"+stringify(u.i)+stringify(u.j)+">)");
		timer_Q.info();
		return result;

	}

	/// returns the non constant part of the MP2 potential which is
	/// (2J-K+Un)|uij>
	real_function_6d get_MP2_potential_residue(const CC_Pair &u)const{
		START_TIMER();
		real_function_6d result = fock_residue_6d(u);
		END_TIMER("(2J-K(R)+Un)|uij>");
		START_TIMER();

		if(parameters.debug){
		// DEBUG
		real_function_6d debug_result=multiply_with_0th_order_Hamiltonian(u.function,u.i,u.j);
		real_function_6d diff = debug_result - result;
		std::cout << "\n\nMP2/CC2 DEBUG: DIFFERENCE BETWEEN MP2 and CC2 IMPLEMENTATION OF MP2 RESIDUE is: "
				<< diff.norm2() << "\n\n" << std::endl;
		// DEBUG END
		END_TIMER("MP2 Potential Debug");
		}
		return result;
	}

	real_function_6d multiply_with_0th_order_Hamiltonian(
			const real_function_6d& f, const int i, const int j) const {

		real_function_6d vphi;

		START_TIMER();
			// the purely local part: Coulomb and U2
			real_function_3d v_local = 2.0*intermediates_.get_hartree_potential()
										+ nemo.nuclear_correlation->U2();

			v_local.print_size("vlocal");
			f.print_size("u");

			// screen the construction of Vphi: do only what is needed to
			// get an accurate result of the BSH operator
			const double eps = get_epsilon(i, j);
			real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), parameters.lo,
					parameters.thresh_bsh_6D);
			op_mod.modified() = true;
			vphi = CompositeFactory<double, 6, 3>(world).ket(copy(f)).V_for_particle1(
					copy(v_local)).V_for_particle2(copy(v_local));
			vphi.fill_tree(op_mod);

			vphi.print_size("vphi: local parts");

			// the part with the derivative operators: U1
			for (int axis = 0; axis < 6; ++axis) {
				real_derivative_6d D = free_space_derivative<double, 6>(world,
						axis);
				const real_function_6d Drhs = D(f).truncate();

				// note integer arithmetic
				if (world.rank() == 0)
					print("axis, axis^%3, axis/3+1", axis, axis % 3, axis / 3 + 1);
				const real_function_3d U1_axis =
						nemo.nuclear_correlation->U1(axis % 3);
				//                    real_function_6d x=multiply(copy(Drhs),copy(U1_axis),axis/3+1).truncate();

				double tight_thresh = std::min(FunctionDefaults<6>::get_thresh(), 1.e-4);
				real_function_6d x;
				if (axis / 3 + 1 == 1) {
					x =CompositeFactory<double, 6, 3>(world).ket(Drhs)
												.V_for_particle1(copy(U1_axis))
												.thresh(tight_thresh);

				} else if (axis / 3 + 1 == 2) {
					x =CompositeFactory<double, 6, 3>(world).ket(Drhs)
												.V_for_particle2(copy(U1_axis))
												.thresh(tight_thresh);
				}
				x.fill_tree(op_mod);
				x.set_thresh(FunctionDefaults<6>::get_thresh());
				vphi += x;
				vphi.truncate().reduce_rank();

			}
			vphi.print_size("(U_nuc + J) |ket>:  made V tree");

		END_TIMER("apply (U + J) |ket>");

		// and the exchange
		START_TIMER();
		vphi = (vphi - K(f, i == j)).truncate().reduce_rank();
		vphi.print_size("(U_nuc + J - K) |ket>:  made V tree");
		END_TIMER("apply K |ket>");

		return vphi;
	}

	// right now this is all copied from mp2.cc
	double compute_mp2_pair_energy(CC_Pair &pair)const{

	    START_TIMER();
		// this will be the bra space
	real_function_6d eri = TwoElectronFactory(world).dcut(parameters.lo);
		real_function_6d ij_g =
				CompositeFactory<double, 6, 3>(world).particle1(
						copy(mo_bra_[pair.i])).particle2(
						copy(mo_bra_[pair.j])).g12(eri);
		real_function_6d ji_g =
				CompositeFactory<double, 6, 3>(world).particle1(
						copy(mo_bra_[pair.i])).particle2(
						copy(mo_bra_[pair.j])).g12(eri);

		// compute < ij | g12 | psi >
		const double ij_g_uij = inner(pair.function, ij_g);
		if (world.rank() == 0)
			printf("<ij | g12       | psi^1>  %12.8f\n", ij_g_uij);

		// compute < ji | g12 | psi > if (i/=j)
		const double ji_g_uij = (pair.i == pair.j) ? 0 : inner(pair.function, ji_g);
		if (world.rank() == 0)
			printf("<ji | g12       | psi^1>  %12.8f\n", ji_g_uij);

		// the singlet and triplet triplet pair energies
		if (pair.i == pair.j) {
			pair.e_singlet = ij_g_uij + pair.ij_gQf_ij;
			pair.e_triplet = 0.0;
		} else {
			pair.e_singlet = (ij_g_uij + pair.ij_gQf_ij)
					+ (ji_g_uij + pair.ji_gQf_ij);
			pair.e_triplet = 3.0
					* ((ij_g_uij - ji_g_uij) + (pair.ij_gQf_ij - pair.ji_gQf_ij));
		}

		// print the pair energies
		if (world.rank() == 0) {
			printf("current energy %2d %2d %12.8f %12.8f\n", pair.i, pair.j,
					pair.e_singlet, pair.e_triplet);
		}

		END_TIMER("compute MP2 energy");
		// return the total energy of this pair
		return pair.e_singlet + pair.e_triplet;
	}



	/// Projectors to project out the occupied space
	// 3D on vector of functions
	void Q(vecfuncT &f) const {
		for (size_t i = 0; i < f.size(); i++)
			Q(f[i]);
	}
	// 3D on single function
	void Q(real_function_3d &f) const {
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			f -= mo_bra_[i].inner(f) * mo_ket_[i];
		}
	}

	/// CCSD/CC2 singles potential parts

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_i = \sum_k <k|r12|k> |tau_i>
	// K_i = \sum_k <k|r12|tau_i> |k>
	vecfuncT fock_residue_closed_shell(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT J = mul(world, intermediates_.get_hartree_potential(), tau);
		truncate(world, J);
		scale(world, J, 2.0);
		END_TIMER("J");
		START_TIMER();
		vecfuncT K;
		for (size_t i = 0; i < tau.size(); i++) {
			real_function_3d tmp = real_factory_3d(world);
			vecfuncT vectmp = mul(world,
					intermediates_.get_perturbed_exchange_intermediate()[i],
					mo_ket_);
			for (size_t j = 0; j < tau.size(); j++)
				tmp += vectmp[j];
			tmp.truncate();
			K.push_back(tmp);
		}
		truncate(world, K);
		scale(world, K, -1);
		END_TIMER("K");
		return add(world, J, K);
	}

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3c(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT result = mul(world,
				intermediates_.get_perturbed_hartree_potential(), mo_ket_);
		Q(result);
		truncate(world, result);
		scale(world, result, 2.0);
		END_TIMER("S3c");
		return result;
	}
	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	//	vecfuncT S3c(const vecfuncT &tau)const{
	//		START_TIMER();
	//		vecfuncT result = mul(world,(*intermediates_.poisson)(intermediates_.make_density(mo_bra_,tau)),mo_ket_);
	//		Q(result);
	//		truncate(world,result);
	//		scale(world,result,2.0);
	//		END_TIMER("S3C_C");
	//		return result;
	//	}

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = -Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3c_X(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT result;
		for (size_t i = 0; i < tau.size(); i++) {
			real_function_3d tmp = real_factory_3d(world);
			vecfuncT vectmp = mul(world,
					intermediates_.get_exchange_intermediate()[i], tau);
			for (size_t j = 0; j < tau.size(); j++)
				tmp += vectmp[j];
			tmp.truncate();
			result.push_back(tmp);
		}
		Q(result);
		truncate(world, result);
		scale(world, result, -1.0);
		END_TIMER("S3c_X");
		return result;
	}

	/// The S5b term
	//[i]    [Q]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// 2\sum_k <k|g|\tau_k> |\tau_i>
	// No Q is applied yet !
	vecfuncT S5b(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT result = mul(world,
				intermediates_.get_perturbed_hartree_potential(), mo_ket_);
		truncate(world, result);
		scale(world, result, 2.0);
		return result;
		END_TIMER("S5b");
	}

	/// The S5b Exchange Term
	//[i]         [Q]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_k <k|g|\tau_i> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5b_X(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT tmp;
		vecfuncT result = zero_functions_compressed<double, 3>(world,
				(tau.size()));
		for (size_t i = 0; i < tau.size(); i++) {
			tmp = mul(world,
					intermediates_.get_perturbed_exchange_intermediate()[i],
					tau);
			for (size_t k = 0; k < tau.size(); k++) {
				result[i] += tmp[k];
			}
		}
		truncate(world, result);
		scale(world, result, -1);
		END_TIMER("S5b_X");
		return result;
	}

	/// The S5c term
	//[Q]    [i]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// -2\sum_kl <kl|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	// May use alteriative algorithm with perturbed density intermediate
	vecfuncT S5c(const vecfuncT&tau) const {
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double, 3>(world,
				tau.size());
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			for (size_t k = 0; k < mo_bra_.size(); k++) {
				for (size_t l = 0; l < mo_bra_.size(); l++) {
					result[i] += intermediates_.get_integrals_mixed_t1()(k, l,
							i, l) * tau[k];
				}
			}
		}
		truncate(world, result);
		scale(world, result, -2.0);
		END_TIMER("S5c");
		return result;
	}

	/// The S5c_X echange term
	//[Q]         [i]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_kl <lk|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5c_X(const vecfuncT&tau) const {
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double, 3>(world,
				tau.size());
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			for (size_t k = 0; k < mo_bra_.size(); k++) {
				for (size_t l = 0; l < mo_bra_.size(); l++) {
					result[i] += intermediates_.get_integrals_mixed_t1()(l, k,
							i, l) * tau[k];
				}
			}
		}
		truncate(world, result);
		scale(world, result, -1.0);
		END_TIMER("S5c_X");
		return result;
	}

	/// The S6+X Term
	// \    /\    /...
	//  \  /  \  /   /\
	//  _\/_  _\/_  _\/_
	// -Q \sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>
	// Q is not applied yet!
	vecfuncT S6(const vecfuncT &tau) const {
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double, 3>(world,
				tau.size());
		for (size_t i = 0; i < tau.size(); i++) {
			for (size_t k = 0; k < mo_bra_.size(); k++) {
				for (size_t l = 0; l < mo_bra_.size(); l++) {
					result[i] += (-2
							* intermediates_.get_integrals_t1()(k, l, k, i)
							- intermediates_.get_integrals_t1()(k, l, i, k))
							* tau[l];
				}
			}
		}
		truncate(world, result);
		END_TIMER("S6+X");
		return result;
	}

	/// CC2 singles diagrams with 6d functions as input
	/// Use GFInterface in function_interface.h as kernel (f*g) and do not reconstruct \tau = f12u(1,2) if possible
	/// Since the correlation factor of CC2 has Slater form like in MP2: g12f12 = g12(1-exp(-mu*r12)/r12) = g12 - exp(-mu*r12)/r12 = Coulomb_Operator - BSH_Operator(mu)

	/// S2b + X Term
	// [i]   [Q]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// Current procedure:
	/// use g12 = \int \delta(1-3) g32 d3
	/// <k(2)|g12|u(1,2)> = \int d2[ g12x(1,2 ] with x(1,2) = k(2)u(1,2)
	/// = int d2 [ int d3[ \delta(1-3) g32 ] x(1,2) ]
	/// = \int d3[\delta(1-3) \int d2 [ g32 x(1,2 ] ]
	/// = \int d3[\delta(1-3) h(1,3)] with h(1,3) = \int d2 g23 x(1,2)
	vecfuncT S2b(const Pairs<real_function_6d> u) const {
		//double prefactor = 1.0 / (2.0 * corrfac.gamma());
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(mo_ket_.size());
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			real_function_3d resulti = real_factory_3d(world);
			for (size_t k = 0; k < mo_ket_.size(); k++) {
				// calculate x(1,2) from u(1,2) and k(2), --> F.A.B uses multiply(copy(f),copy(bra) ...) deep copy of functions (dont know why)
				real_function_6d xik = multiply(u(i, k), mo_bra_[k], 2);
				real_function_6d xki = multiply(u(k, i), mo_bra_[k], 2);
				// calculate the convolution with fg = 1/(2gamma)*(Coulomb - 4pi*BSH(gamma))
				real_function_6d hik;
				real_function_6d hki;
				{
//					real_function_6d CoulombTerm = ((*poisson)(xik)).truncate();
//					real_function_6d fBSHTerm = 4.0*constants::pi*((*fBSH)(xik)).truncate();
//					hik = prefactor * (CoulombTerm + fBSHTerm);
					hik = apply_gf(xik,2);
				}
				{
//					real_function_6d CoulombTerm = ((*poisson)(xki)).truncate();
//					real_function_6d fBSHTerm = 4.0*constants::pi*((*fBSH)(xki)).truncate();
//					hki = prefactor * (CoulombTerm + fBSHTerm);
					hki = apply_gf(xki,2);
				}
				// Make the projection to 3D with the unit function
				real_function_3d resultik = hik.project_out(unity, 3);
				real_function_3d resultki = hki.project_out(unity, 3);
				resultik.truncate();
				resulti += (2.0 * resultik - resultki);
			}
			result[i] = resulti;
		}
		Q(result);
		truncate(world, result);
		return result;
	}

	/// S2c + X Term
	// [Q]   [i]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// = \sum <k(3)l(4)|g34 f31| u_{lk}(1,3) i(4)> = \sum <k(3)|f13| X_{lk,li}(1,3) > with X_{lk,li}(1,3) = u_lk(1,3) * (<l(4)|g34>|i(4)>_4)(3)
	/// = \sum \int d5 \delta(5-1) \int d3 f53 k(3)*X_{lk,li}(5,3)
	vecfuncT S2c(const Pairs<real_function_6d> u) const {
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(mo_ket_.size());
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			real_function_3d resulti = real_factory_3d(world);
			for (size_t k = 0; k < mo_ket_.size(); i++) {
				for (size_t l = 0; l < mo_ket_.size(); l++) {
					// make X_lkli(5,3) = ulk(5,3) * (<l(4)|g34>|i(4)>_4)(3) and X
					real_function_6d xlkli = multiply(u(l, k),
							intermediates_.get_exchange_intermediate()[i][l],
							2);
					real_function_6d xlkki = multiply(u(l, k),
							intermediates_.get_exchange_intermediate()[i][k],
							2);
					// make Atmp =  k(3) * X_lkli(5,3) and X
					real_function_6d Aklkli = multiply(xlkli, mo_bra_[k], 2);
					real_function_6d Allkki = multiply(xlkki, mo_bra_[l], 2);
					// make f12 convolution
					real_function_6d tmpklkli = (*f12op)(Aklkli);
					real_function_6d tmpllkki = (*f12op)(Allkki);
					// Project to 3d
					resulti -= (2.0 * tmpklkli.project_out(unity, 0)
							- tmpllkki.project_out(unity, 0));
				}
			}
			resulti.truncate();
			result[i] = resulti;
		}
		Q(result);
		return result;
	}

	/// The S4a + X diagram
	//[Q]       [i]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k>
	vecfuncT S4a(const Pairs<real_function_6d> u, const vecfuncT & tau) const {
		vecfuncT result(mo_ket_.size());
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			real_function_3d resulti = real_factory_3d(world);
			for (size_t k = 0; k < mo_ket_.size(); k++) {
				for (size_t l = 0; l < mo_ket_.size(); l++) {
					// Coulomb Part of f12g12 = g12 - BSH
					{
						real_function_6d eri = TwoElectronFactory(world).dcut(
								FunctionDefaults<3>::get_thresh());
						real_function_6d kl_g = CompositeFactory<double, 6, 3>(
								world).particle1(copy(mo_bra_[k])).particle2(
								copy(mo_bra_[l])).g12(eri);
						resulti -= (2.0 * inner(u(i, l), kl_g) * tau[k]
								- inner(u(i, k), kl_g) * tau[l]);
						resulti.truncate();
					}
					// BSH part of f12g12
					{
						real_function_6d bsh_kernel =
								TwoElectronFactory(world).BSH().dcut(
										FunctionDefaults<3>::get_thresh());
						real_function_6d kl_bsh =
								CompositeFactory<double, 6, 3>(world).particle1(
										copy(mo_bra_[k])).particle2(
										copy(mo_bra_[l])).g12(bsh_kernel);
						resulti -= (2.0 * inner(u(i, l), kl_bsh) * tau[k]
								- inner(u(i, k), kl_bsh) * tau[l]);
						resulti.truncate();
					}
				}
			}
			result[i] = resulti;
		}
		Q(result);
		return result;
	}

	/// The S4b
	//[i]       [Q]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum_{kl} (2<k(3)l(4)|g34f14|\tau_{i}(3)u_{kl}(1,4)>  // exchange part - <k(4)l(3)|g34f14|\tau_i(3)u_{lk}(1,4)>)
	// 1. make exchange intermedaite X(4) = <k(3)|g34|\tau_i(3)>_3 *  l(4)			Exchange part : Xx(4) = <l(3)|g34|\tau_i(3)>(4) * k(4)
	// 2. make 6d intermediate Y(1,4) = X(4)* u_{kl}(1,4)							Exchange part : Yx(1,4) = X(4)*u_{lk}(1,4)
	// 3. make f14 integration via delta function trick: result(1) = \int f14 Y(1,4) d4 = \int delta(5-1) (\int f54 Y(1,4) d4)d5
	// 3.1 do the convolution Z(1,5) = \int f54 Y(1,4) d4							Exchange part: Zx(1,5) = int f54 Yx(1,4)d4
	// 3.2 project out the unit function: result(1) = <I(5)|Z(1,5)>_5				Exchange part: resultx(1) = <I(5)|Zx(1,5>_5
	vecfuncT S4b(const Pairs<real_function_6d> u, const vecfuncT & tau) const {
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(mo_ket_.size());
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			real_function_3d resulti = real_factory_3d(world);
			for (size_t k = 0; i < mo_ket_.size(); k++) {
				for (size_t l = 0; l < mo_ket_.size(); l++) {
					// Make 3d Intermediates
					real_function_3d X =
							intermediates_.get_perturbed_exchange_intermediate()[i][k]
									* mo_bra_[l];
					real_function_3d Xx =
							intermediates_.get_perturbed_exchange_intermediate()[i][l]
									* mo_bra_[k];
					// Make 6d Intermediates
					real_function_6d Y = (multiply(u(k, l), X, 2)).truncate();
					real_function_6d Yx = (multiply(u(l, k), Xx, 2)).truncate();
					// Do the convolution with f
					real_function_6d Z = (*f12op)(Y);
					real_function_6d Zx = (*f12op)(Yx);
					// Project out the second particle
					resulti -= 2.0 * Z.project_out(unity, 2);
					resulti += Zx.project_out(unity, 2);
				}
			}
			result[i] = resulti;
		}
		Q(result);
		return result;
	}

	/// The S4c + X + X + X + X Diagrams
	//            [i]   [Q]
	//   .......   \    /
	//  /\     /\   \  /
	// _\/_   _\/____\/_
	/// Q\sum_{kl}[ 4*<k(3)l(4)|g34 f14| \tau_k(3) u_{il}(1,4)> - 2* <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)>
	/// - 2* <k(3)l(4)|g34 f14| \tau_k(3) U_{li}(1,4)> + <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)>  ]
	// First and third Terms are solved like this:
	// 1. X(4) = \sum_k (<k(3)|g34|\tau_k(3)>_3(4)) * l(4) = perturbed_hartree_potential(4) * l(4)
	// 2. Y(1,4) = X(4) u_{il}(1,4)			Exchange Part: Yx(4,1) = X(4) u_{li}(4,1)
	// 3.1 Z(1,5) = \int f54 Y(1,4) d4		Exchange Part: Zx(5,1) = \int f54 Yx(4,1) d4
	// 3.2 result(1) = -4 <I(5)|Z(1,5)>_5 -2 <I(5)|Zx(1,5)>_5
	// Second and fourth terms can not use the perturbed hartree potential
	vecfuncT S4c(const Pairs<real_function_6d> u, const vecfuncT & tau) const {
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(tau.size());
		for (size_t i = 0; i < result.size(); i++) {
			real_function_3d resulti = real_factory_3d(world);
			// Term 1 and 3
			for (size_t l = 0; l < mo_ket_.size(); l++) {
				real_function_3d X =
						intermediates_.get_perturbed_hartree_potential()
								* mo_bra_[l];
				real_function_6d Y = (multiply(u(i, l), X, 2)).truncate();
				real_function_6d Yx = (multiply(u(l, i), X, 1)).truncate();
				real_function_6d Z = (*f12op)(Y);
				real_function_6d Zx = (*f12op)(Yx);
				resulti += (4.0 * Z.project_out(unity, 2)
						- 2.0 * Zx.project_out(unity, 1));
			}
			// Term 2 and 4
			for (size_t k = 0; k < mo_ket_.size(); k++) {
				for (size_t l = 0; l < mo_ket_.size(); l++) {
					real_function_3d X =
							intermediates_.get_perturbed_exchange_intermediate()[k][l]
									* mo_bra_[k];
					real_function_6d Y = (multiply(u(i, l), X, 2)).truncate();
					real_function_6d Yx = (multiply(u(i, l), X, 1)).truncate();
					real_function_6d Z = (*f12op)(Y);
					real_function_6d Zx = (*f12op)(Yx);
					resulti += (-2.0 * Z.project_out(unity, 2)
							+ Zx.project_out(unity, 1));
				}
			}
			result[i] = resulti;
		}
		Q(result);
		truncate(world, result);
		return result;
	}

	// CC2 Doubles Potential

	// MP2 Terms are
	// fock_residue_6d = (2J - Kn + Un) |u_{ij}> , KR = R12^{-1}*K*R12 (nuclear tranformed K)
	// Uen|ij> = R12{-1}*U_e*R12 |ij>

	/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
	real_function_6d fock_residue_6d(const CC_Pair &u) const {
		const double eps = get_epsilon(u.i, u.j);
		// make the coulomb and local Un part with the composite factory
		real_function_3d local_part = (2.0
				* intermediates_.get_hartree_potential()
				+ nemo.nuclear_correlation->U2());
		local_part.print_size("vlocal");
		u.function.print_size("u");

		// Contruct the BSH operator in order to screen

		real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps),
				parameters.lo, parameters.thresh_bsh_6D);
		// apparently the modified_NS form is necessary for the screening procedure
		op_mod.modified() = true;
		// Make the CompositeFactory
		real_function_6d vphi =
				CompositeFactory<double, 6, 3>(world).ket(copy(u.function)).V_for_particle1(
						copy(local_part)).V_for_particle2(copy(local_part));
		// Screening procedure
		vphi.fill_tree(op_mod);

		vphi.print_size("vlocal|u>");

		// the part with the derivative operators: U1
		for (int axis = 0; axis < 6; ++axis) {
			real_derivative_6d D = free_space_derivative<double, 6>(world,
					axis);
			// Partial derivative of the pari function
			const real_function_6d Du = D(u.function).truncate();

			// % operator gives division rest (modulo operator)
			if (world.rank() == 0)
				print("axis, axis^%3, axis/3+1", axis, axis % 3, axis / 3 + 1);
			const real_function_3d U1_axis = nemo.nuclear_correlation->U1(
					axis % 3);

			double tight_thresh = parameters.thresh_Ue;
			real_function_6d x;
			if (axis / 3 + 1 == 1) {
				x =
						CompositeFactory<double, 6, 3>(world).ket(Du).V_for_particle1(
								copy(U1_axis)).thresh(tight_thresh);

			} else if (axis / 3 + 1 == 2) {
				x =
						CompositeFactory<double, 6, 3>(world).ket(Du).V_for_particle2(
								copy(U1_axis)).thresh(tight_thresh);
			}
			x.fill_tree(op_mod);
			x.set_thresh(FunctionDefaults<6>::get_thresh());
			vphi += x;
			vphi.truncate().reduce_rank();
		}

		vphi.print_size("(Un + J1 + J2)|u>");

		// Exchange Part
		vphi = (vphi - K(u.function, u.i == u.j)).truncate().reduce_rank();
		vphi.print_size("(Un + J1 + J2 - K1 - K2)|U>");
		vphi.truncate();
		vphi.print_size("truncated: (Un + J1 + J2 - K1 - K2)|U>");
		return vphi;

	}

	/// Echange Operator on 3D function
	/// !!!!Prefactor (-1) is not included
	real_function_3d K(const real_function_3d &f,const size_t &i, const bool hc=false)const{
		real_function_3d result = real_factory_3d(world);
		if(hc==true) MADNESS_EXCEPTION("ERROR in K, hc=true not implemented",1);

		for (std::size_t k = 0; k < mo_ket_.size(); ++k) {
			result += mo_ket_[k]*intermediates_.get_exchange_intermediate()[i][k];
		}

		// Sanity Check (expensive when not helium)
		if(mo_ket_.size()<3){
			real_function_3d result2 = real_factory_3d(world);
			// multiply rhs with R2orbitals (the bra space)
			vecfuncT R2rhs = mul(world, f, mo_bra_);
			for (std::size_t k = 0; k < mo_ket_.size(); ++k) {
				result2 += mo_ket_[k] * (*poisson)(R2rhs[k]);
			}
			double sanity = (result-result2).norm2();
			if(sanity < FunctionDefaults<3>::get_thresh()){
				if(world.rank()==0) std::cout << "Sanity Check of K passed\n";
			}else{
				if(world.rank()==0) std::cout << "Sanity Check of K NOT passed\n";
			}
		}

		return result;
	}

	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// if i==j in uij then the symmetry will be exploited
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d K(const real_function_6d &u,
			const bool symmetric = false) const {
		/// DEBUG
		if(world.rank()==0)std::cout << "Entering K" << std::endl;
		/// DEBUG END

		/// TEST IF THIS WILL WORK FOR THE += Operator
		real_function_6d result = real_factory_6d(world).compressed();
		// K(1) Part
		result += apply_K(u, 1);
		// K(2) Part
		if (symmetric)
			result += swap_particles(result);
		else
			result += apply_K(u, 2);

		return (result.truncate());
	}

	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
	/// 1. X(3,2) = bra_k(3)*u(3,2)
	/// 2. Y(1,2) = \int X(3,2) g13 d3
	/// 3. result = Y(1,2)*ket_k(1)
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d apply_K(const real_function_6d &u,
			const size_t &particle) const {
		/// DEBUG
		if(world.rank()==0)std::cout << "Entering apply_K" << std::endl;
		/// DEBUG END
		MADNESS_ASSERT(particle == 1 or particle == 2);
		poisson->particle() = particle;
		/// WARNING: CHECK IF THIS WORKS -> bc of the += operator later
		real_function_6d result = real_factory_6d(world).compressed();
		for (size_t k = 0; k < mo_ket_.size(); k++) {
			real_function_6d X = (multiply(copy(u), copy(mo_bra_[k]), particle)).truncate();
			real_function_6d Y = (*poisson)(X);
			result += multiply(copy(Y), copy(mo_ket_[k]), particle).truncate();
		}
		return result;
	}

	/// Apply Ue on a tensor product of two 3d functions: Ue(1,2) |x(1)y(2)> (will be either |ij> or |\tau_i\tau_j> or mixed forms)
	/// The Transformed electronic regularization potential (Kutzelnigg) is R_{12}^{-1} U_e R_{12} with R_{12} = R_1*R_2
	/// It is represented as: R_{12}^{-1} U_e R_{12} = U_e + R^-1[Ue,R]
	/// where R^-1[Ue,R] = R^-1 [[T,f],R] (see: Regularizing the molecular potential in electronic structure calculations. II. Many-body
	/// methods, F.A.Bischoff)
	/// The double commutator can be evaluated as follows:  R^-1[[T,f],R] = -Ue_{local}(1,2)*(Un_{local}(1) - Un_{local}(2))
	/// @param[in] x the 3D function for particle 1
	/// @param[in] y the 3D function for particle 2
	/// @param[in] i the first index of the current pair function (needed to construct the BSH operator for screening)
	/// @param[in] j the second index of the current pair function
	/// @param[out]  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
	real_function_6d apply_transformed_Ue(const real_function_3d x,
			const real_function_3d y, const size_t &i, const size_t &j, CC_Pair &u) const {
		real_function_6d Uxy = real_factory_6d(world);
		// Apply the untransformed U Potential
		const double eps = get_epsilon(i, j);
		Uxy = corrfac.apply_U(x, y, eps);

		// Get the 6D BSH operator in modified-NS form for screening
		real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps),
				parameters.lo,
				parameters.thresh_Ue);
		op_mod.modified() = true;

		// make shure the thresh is high enough
		double tight_thresh = parameters.thresh_6D_tight;

		// Apply the double commutator R^{-1}[[T,f,R]
		for (size_t axis = 0; axis < 3; axis++) {
			// Make the local parts of the Nuclear and electronic U potentials
			const real_function_3d Un_local = nemo.nuclear_correlation->U1(
					axis);
			const real_function_3d Un_local_x = (Un_local * x).truncate();
			const real_function_3d Un_local_y = (Un_local * y).truncate();
			const real_function_6d Ue_local = corrfac.U1(axis);
			// Now add the Un_local_x part to the first particle of the Ue_local potential
			real_function_6d UeUnx = CompositeFactory<double, 6, 3>(world).g12(
					Ue_local).particle1(Un_local_x).particle2(copy(y)).thresh(
					tight_thresh);
			// Fill the Tree were it will be necessary
			UeUnx.fill_tree(op_mod);
			// Set back the thresh
			UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());

			UeUnx.print_size("UeUnx");

			// Now add the Un_local_y part to the second particle of the Ue_local potential
			real_function_6d UeUny = CompositeFactory<double, 6, 3>(world).g12(
					Ue_local).particle1(copy(x)).particle2(Un_local_y).thresh(
					tight_thresh);
			// Fill the Tree were it will be necessary
			UeUny.fill_tree(op_mod);
			// Set back the thresh
			UeUny.set_thresh(FunctionDefaults<6>::get_thresh());

			UeUny.print_size("UeUny");

			// Construct the double commutator part and add it to the Ue part
			real_function_6d diff = (UeUnx - UeUny).scale(-1.0);
			diff.truncate();
			Uxy = (Uxy+diff).truncate();
		}

		// sanity check: <xy|R2 [T,g12] |xy> = <xy |R2 U |xy> - <xy|R2 g12 | xy> = 0
		real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
				copy(mo_bra_[i])).particle2(copy(mo_bra_[j]));

		const double a = inner(Uxy, tmp);
		const real_function_3d xx = (x * x*nemo.nuclear_correlation -> square());
		const real_function_3d yy = (y * y*nemo.nuclear_correlation -> square());
		const real_function_3d gxx = (*poisson)(xx);
		const double aa = inner(yy, gxx);
		const double error = std::fabs(a - aa);
		if (world.rank() == 0) {
			printf("< phi0 | U_R   | phi0 >  %12.8f\n", a);
			printf("< phi0 | 1/r12 | phi0 >  %12.8f\n", aa);
			if (error > FunctionDefaults<6>::get_thresh())
				print("WARNING : Kutzelnigg's potential inaccurate (box size, thresh ?)");
			//if (error > FunctionDefaults<6>::get_thresh() * 10.0)
			//	MADNESS_EXCEPTION("Kutzelnigg's potential plain wrong (box size, thresh ?)", 1);
		}
		Uxy.print_size("Uphi0");

		return Uxy;
	}

	/// Apply the Exchange Commutator [K,f]|xy>
	real_function_6d apply_exchange_commutator(const real_function_3d &x, const real_function_3d &y,const std::string &type, const size_t &i, const size_t &j)const{
		MADNESS_ASSERT(
				type == "occupied" or type == "mixed" or type == "virtual");


		// make first part of commutator
		real_function_6d Kfxy = apply_Kf(x,y,type,i,j).truncate();

		// for sanity check:
		double expv_first_part = 0.0;
		double expv_second_part = 0.0;
		if(type=="occupied"){
			real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[i])).particle2(copy(mo_bra_[j]));
			expv_first_part = inner(Kfxy,tmp);
		}


		// make the second part of the commutator
		real_function_6d fKxy = apply_fK(x,y,type,i,j).truncate();

		// fot the sanity check
		if(type=="occupied"){
			real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[i])).particle2(copy(mo_bra_[j]));
			expv_second_part = inner(fKxy,tmp);
		}

		if(type=="occupied"){
		if(world.rank()==0){
			std::cout << "Apply [K,f]|x,y> sanity check:";
			std::cout <<  "\n<ij|Kf|ij> =" << expv_first_part;
			std::cout <<  "\n<ij|fK|ij> =" << expv_second_part;
		}
		}

		real_function_6d result = (Kfxy - fKxy);

		if(type=="occupied"){
		// sanity check: The Expectation value of the Kommutator must vanish (symmetry)
		// <0|[A,B]|0> = <0|AB|0> - <0|BA|0> since A=f and B=K -> both are hermitian
			real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[i])).particle2(copy(mo_bra_[j]));
		const double a = inner(result, tmp);
		if (world.rank() == 0) {
			printf("< nemo0 | R^2 R-1 [K,f] R | nemo0 >  %12.8f\n", a);
			if (std::fabs(a) > FunctionDefaults<6>::get_thresh())
				print("WARNING : exchange commutator inaccurate");
			if (std::fabs(a) > FunctionDefaults<6>::get_thresh() * 10.0)
				MADNESS_EXCEPTION("exchange commutator plain wrong", 1);
		}
		}

		return result;
	}

	/// Apply the Exchange operator on a tensor product multiplied with f12
	/// !!! Prefactor of (-1) is not inclued in K here !!!!
	real_function_6d apply_Kf(const real_function_3d &x,
			const real_function_3d &y, const std::string &type, const size_t &i, const size_t &j) const {
		MADNESS_ASSERT(
				type == "occupied" or type == "mixed" or type == "virtual");

		bool symmetric = false;
		if(type == "occupied" and i==j) symmetric = true;
		if(type == "virtual" and i==j) symmetric = true;
		if(type == "mixed") symmetric = false;

		START_TIMER();
		// First make the 6D function f12|x,y>
		real_function_6d f12xy = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(x)).particle2(copy(y));
		f12xy.fill_tree().truncate().reduce_rank();
		END_TIMER("Constructed f12|xy>");
		// Apply the Exchange Operator
		real_function_6d result = K(f12xy, symmetric);
		return result.truncate();
	}

	/// Apply fK on a tensor product of two 3D functions
	/// fK|xy> = fK_1|xy> + fK_2|xy>
	/// @param[in] x, the first 3D function in |xy>
	/// @param[in] y, the second 3D function in |xy>
	/// @param[in] type, specifies if |xy> = |ij> (occupied), |xy> = |\tau_i,j> (mixed) or |xy> = |\tau_i\tau_j> (virtual)
	/// @param[in] i, the number of the function: bsp if occupied then x_i = |i>, if virtual then x_i = \tau_i etc
	/// @param[in] j , index of the second function
	real_function_6d apply_fK(const real_function_3d &x,
			const real_function_3d &y, const std::string &type, const size_t &i,
			const size_t &j) const {
		MADNESS_ASSERT(type == "occupied" or type == "mixed" or type == "virtual");

		const real_function_3d& phi_i = x;
		const real_function_3d& phi_j = y;

		const real_function_3d Kphi_i = K(phi_i,i,false);
		const real_function_3d Kphi_j = K(phi_j,j,false);

		real_function_6d fKphi0a = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(phi_i)).particle2(copy(Kphi_j));
		fKphi0a.fill_tree().truncate();
		real_function_6d fKphi0b = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(Kphi_i)).particle2(copy(phi_j));
		fKphi0b.fill_tree().truncate();

		real_function_6d fKphi0 = (fKphi0a + fKphi0b).truncate();
		return fKphi0;

	}

	// gives back \epsilon_{ij} = \epsilon_i + \epsilon_j
	double get_epsilon(const size_t &i, const size_t &j) const {
		double eps = (nemo.get_calc()->aeps(i) + nemo.get_calc()->aeps(j));
		if(world.rank()==0)std::cout << "Epsilon:" << eps << std::endl;
		return (nemo.get_calc()->aeps(i) + nemo.get_calc()->aeps(j));
	}

	/// swap particles 1 and 2

	/// param[in]	f	a function of 2 particles f(1,2)
	/// return	the input function with particles swapped g(1,2) = f(2,1)
	real_function_6d swap_particles(const real_function_6d& f) const {
		CC_Timer timer_swap(world,"swap particles");
		// this could be done more efficiently for SVD, but it works decently
		std::vector<long> map(6);
		map[0] = 3;
		map[1] = 4;
		map[2] = 5;	// 2 -> 1
		map[3] = 0;
		map[4] = 1;
		map[5] = 2;	// 1 -> 2
		timer_swap.info();
		return mapdim(f, map);
	}

	// Calculate the CC2 energy equation which is
	// \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
	// with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
	double get_CC2_correlation_energy() const {
		MADNESS_EXCEPTION("get_cc2_correlation_energy not implemented yet",1);
		return 0.0;
	}
	double get_CC2_pair_energy(const CC_Pair &u,
			const real_function_3d &taui, const real_function_3d &tauj) const {
		double omega = 0.0;
		const size_t i = u.i;
		const size_t j = u.j;
		// Contribution from u itself, we will calculate <uij|g|ij> instead of <ij|g|uij> and then just make the inner product (see also mp2.cc)
		{
			real_function_6d coulomb = TwoElectronFactory(world).dcut(
					FunctionDefaults<6>::get_thresh());
			real_function_6d g_ij =
					CompositeFactory<double, 6, 3>(world).particle1(
							copy(mo_bra_[i])).particle2(copy(mo_bra_[j])).g12(
							coulomb);
			real_function_6d g_ji =
					CompositeFactory<double, 6, 3>(world).particle1(
							copy(mo_bra_[j])).particle2(copy(mo_bra_[i])).g12(
							coulomb);
			const double uij_g_ij = inner(u.function, g_ij);
			const double uij_g_ji = inner(u.function, g_ji); // =uji_g_ij
			omega += 2.0 * uij_g_ij - uij_g_ji;
		}
		// Contribution from the mixed f12(|\tau_i,j>+|i,\tau_j>) part
		{

		}
		// Contribution from the f12|ij> part, this should be calculated in the beginning
		{
			omega += (2.0*u.ij_gQf_ij - u.ji_gQf_ij );
		}
		// Contribution from the f12|\tau_i\tau_j> part
		{

		}
		// Singles Contribution
		{
			// will be added in the correlation energy function (perturbed density can be used with the summation)
		}
		return omega;
	}

	/// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
	// the bra elements are always the R2orbitals
	// the ket elements can be \tau_i , or orbitals dependet n the type given
	double make_ij_gQf_ij(const size_t &i, const size_t &j,CC_Pair &u)const{
		double result=0.0;
		double exchange_result = 0.0;
		real_function_3d ii = (mo_bra_[i]*mo_ket_[i]).truncate();
		real_function_3d jj = (mo_bra_[j]*mo_ket_[j]).truncate();
		real_function_3d ji = (mo_bra_[j]*mo_ket_[i]).truncate();
		real_function_3d ij = (mo_bra_[i]*mo_ket_[j]).truncate();

		// the pure part
		real_function_3d tmp = apply_gf(ii);
		real_function_3d tmp_ex = apply_gf(ij);
		result += tmp.inner(jj);
		exchange_result += tmp_ex.inner(ji);

		// the O1 part: \sum_k<ij|g|k><k|f|ij> and exchange part <ij|g|k><k|f|ji>
		for(size_t k=0;k<mo_ket_.size();k++){
			real_function_3d igk = intermediates_.get_exchange_intermediate()[k][i];
			real_function_3d kfi = (*f12op)((mo_bra_[k]*mo_ket_[i]).truncate());
			real_function_3d kfj = (*f12op)((mo_bra_[k]*mo_ket_[j]).truncate());
			real_function_3d igkkfi = (igk*kfi).truncate();
			real_function_3d igkkfj = (igk*kfj).truncate();
			result -= igkkfi.inner(jj);
			exchange_result -= igkkfj.inner(ji);
		}
		// the O2 part: (can be speed up by using kfi and kfj from the O1 part)
		for(size_t k=0;k<mo_ket_.size();k++){
			real_function_3d jgk = intermediates_.get_exchange_intermediate()[k][i];
			real_function_3d kfj = (*f12op)((mo_bra_[k]*mo_ket_[j]).truncate());
			real_function_3d kfi = (*f12op)((mo_bra_[k]*mo_ket_[i]).truncate());
			real_function_3d jgkkfj = (jgk*kfj).truncate();
			real_function_3d jgkkfi = (jgk*kfi).truncate();
			result -= jgkkfj.inner(ii);
			exchange_result -= jgkkfi.inner(ij);
		}
		// the O12 part (may include O1 and O2 parts in the first loop and avoid recalulating of kfj and kfi
		// and also the mo_bra_[k]*i etc products
		for(size_t k=0;k<mo_ket_.size();k++){
			for(size_t l=0;l<mo_ket_.size();l++){
				// make <ij|g|kl>
				real_function_3d jgk = intermediates_.get_exchange_intermediate()[k][j];
				double ijgkl = jgk.inner(mo_bra_[i]*mo_ket_[l]);
				// make <kl|f|ij>
				real_function_3d kfi = (*f12op)((mo_bra_[k]*mo_ket_[i]).truncate());
				double klfij = kfi.inner(mo_bra_[l]*mo_ket_[j]);
				// make <kl|f|ji> for exchange part
				real_function_3d kfj = (*f12op)((mo_bra_[k]*mo_ket_[j]).truncate());
				double klfji = kfj.inner(mo_bra_[l]*mo_ket_[i]);
				// Get result
				result += (ijgkl*klfij);
				exchange_result += (ijgkl*klfji);
			}
		}
		if(world.rank()==0) std::cout << "\nCalculated: 2.0 <" << stringify(i) << stringify(j) << "|gQf|" << stringify(i) << stringify(j) << "> = " << result ;
		if(world.rank()==0) std::cout << "\nCalculated: - <" << stringify(i) << stringify(j) << "|gQf|" << stringify(j) << stringify(i) << "> = " << exchange_result<<"\n" ;
		u.ij_gQf_ij = result;
		u.ji_gQf_ij = exchange_result;
		return result+exchange_result;
	}

	/// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
	/// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
	real_function_3d apply_gf(const real_function_3d &f)const{
		double bsh_prefactor = 4.0 * constants::pi;
		double prefactor = 1.0/(2.0*corrfac.gamma());
		return prefactor*((*poisson)(f) - bsh_prefactor*(*fBSH)(f)).truncate();
	}
	real_function_6d apply_gf(const real_function_6d &f,const size_t &particle)const{
		poisson->particle()=particle;
		fBSH->particle()=particle;
		double bsh_prefactor = 4.0 * constants::pi;
		double prefactor = 1.0/(2.0*corrfac.gamma());
		return prefactor*((*poisson)(f) - bsh_prefactor*(*fBSH)(f)).truncate();
	}

	real_function_6d test_fill_tree()const{
		return real_factory_6d(world);
	}

private:
	/// The World
	World &world;
	/// Nemo
	const Nemo &nemo;
	/// Thresh for the bsh operator
	double bsh_eps = std::min(FunctionDefaults<6>::get_thresh(), 1.e-4);
	/// Electronic correlation factor
	CorrelationFactor corrfac;
	/// All necessary parameters
	const CC_Parameters &parameters;
	/// The ket and the bra element of the occupied space
	/// if a  nuclear correlation factor is used the bra elements are the MOs multiplied by the squared nuclear correlation factor (done in the constructor)
	const vecfuncT mo_bra_;
	const vecfuncT mo_ket_;
	/// Helper function to initialize the const mo_bra and ket elements
	vecfuncT make_mo_bra(const Nemo &nemo) const {
		START_TIMER();
		return mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		END_TIMER("Initialized molecular orbital bra_elements");
	}
	/// The poisson operator (Coulomb Operator)
	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr
			< real_convolution_3d
			> (CoulombOperatorPtr(world,parameters.lo,
					parameters.thresh_poisson_3D));
	/// The BSH Operator for the f12g12 convolution which is with f12= 1/(2gamma)[1-exp(-gamma*r12)], f12g12 = 1/(2gamma) [CoulombOp - BSHOp(gamma)]
	std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr
			< real_convolution_3d
			> (BSHOperatorPtr3D(world, corrfac.gamma(),
					parameters.lo,
					parameters.thresh_f12));
	/// The f12 convolution operator
	std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr
			< real_convolution_3d
			> (SlaterF12OperatorPtr(world, corrfac.gamma(),
					parameters.lo,
					parameters.thresh_f12));
	/// Intermediates (some need to be refreshed after every iteration)
	CC_Intermediates intermediates_;

	/// Take the time of operations
	bool use_timer_;
	mutable double ttt, sss;
	void START_TIMER() const {
		if (use_timer_)
			world.gop.fence();
		ttt = wall_time();
		sss = cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if (use_timer_)
			END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if (use_timer_) {
			ttt = wall_time() - ttt;
			sss = cpu_time() - sss;
			if (world.rank() == 0)
				printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}
};

} /* namespace madness */

#endif /* CCOPERATORS_H_ */
