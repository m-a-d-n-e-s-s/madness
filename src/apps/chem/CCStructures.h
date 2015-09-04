/*
 * CCStructures.h
 *
 *  Created on: Sep 3, 2015
 *      Author: kottmanj
 */


/// File holds all helper structures necessary for the CC_Operator and CC2 class
#ifndef CCSTRUCTURES_H_
#define CCSTRUCTURES_H_

#include <chem/SCFOperators.h>
#include <chem/electronic_correlation_factor.h>
#include <algorithm> // tolower function for strings

namespace madness{



typedef std::vector<Function<double, 3> > vecfuncT;

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
        	//std::cout << "input tag is: " << s << std::endl;
        	std::transform(s.begin(),s.end(),s.begin(), ::tolower);
        	//std::cout << "transformed input tag is: " << s << std::endl;
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
            else if ((s == "corrfac_gamma") or (s== "gamma")) f>>corrfac_gamma;
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
			std::cout << std::setw(20) << std::setfill(' ') << "Kain is: " << kain << std::endl;
			if(kain) std::cout << std::setw(20) << std::setfill(' ') << "Kain subspace: " << kain_subspace << std::endl;



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

class CC_Single{
public:
	// Constructors
	CC_Single() : converged(false),iterations(0){}
	CC_Single(const size_t ii): i(ii), converged(false),iterations(0) {}
	CC_Single(const size_t i, const real_function_3d &f): i(i), converged(false),iterations(0), function_(f){}
	CC_Single(const CC_Single &other): i(other.i), converged(other.converged),iterations(other.iterations), function_(other.function()){}

	real_function_3d function()const{return function_;}

	size_t i;
	mutable bool converged;
	mutable size_t iterations;
private:
	real_function_3d function_;
	static size_t uninitialized(){return 999;}


};

typedef  std::vector<CC_Single> CC_Singles;

}

#endif /* CCSTRUCTURES_H_ */
