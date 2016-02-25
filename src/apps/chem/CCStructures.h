/*
 * CCStructures.h
 *
 *  Created on: Sep 3, 2015
 *      Author: kottmanj
 */


/// File holds all helper structures necessary for the CC_Operator and CC2 class
#ifndef CCSTRUCTURES_H_
#define CCSTRUCTURES_H_

//#include <chem/SCFOperators.h>
#include <chem/electronic_correlation_factor.h>
#include <algorithm> // tolower function for strings
#include <examples/nonlinsol.h>


namespace madness{

enum optype {g12_,f12_};
enum calctype {MP2_, CC2_, CCS_response, CC2_response_, CISpD_, experimental_};
enum functype {HOLE,PARTICLE,MIXED,RESPONSE,UNDEFINED};
enum pairtype {GROUND_STATE,EXCITED_STATE};
enum potentialtype_s {pot_F3D_, pot_S2b_u_, pot_S2c_u_, pot_S4a_u_, pot_S4b_u_, pot_S4c_u_,pot_S2b_r_, pot_S2c_r_, pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_ccs_};
enum potentialtype_d {pot_F6D_, pot_cc2_coulomb_,pot_cc2_residue_};
// The pair function is:  \tau = u + Qf(|titj>), FULL means that \tau is calculated in 6D form, DECOMPOSED means that u is used in 6D and the rest is tried to solve in 3D whenever possible
enum pair_function_form{DECOMPOSED, FULL};
static std::string assign_name(const optype &input){
  switch(input){
    case g12_ : return "g12";
    case f12_ : return "f12";
  }
  MADNESS_EXCEPTION("bad string assign name",1);
}
static calctype assign_calctype(const std::string name){
  if(name=="mp2") return MP2_;
  else if(name=="cc2") return CC2_;
  else if(name=="cc2_response") return CC2_response_;
  else if(name=="cispd") return CISpD_;
  else if(name=="cis" or name=="ccs" or name=="ccs_response") return CCS_response;
  else if(name=="experimental") return experimental_;
  else{
    std::string msg= "CALCULATION OF TYPE: " + name + " IS NOT KNOWN!!!!";
    MADNESS_EXCEPTION(msg.c_str(),1);
  }
}
static std::string assign_name(const calctype &inp){
	switch(inp){
	case CC2_ : return "CC2";
	case MP2_ : return "MP2";
	case CC2_response_ : return "CC2-Response";
	case CISpD_ : return "CIS(D)";
	case CCS_response: return "CCS/CIS";
	}
	return "unknown";
}
static std::string assign_name(const potentialtype_s &inp){
	switch(inp){
	case pot_F3D_ : return "Fock-Residue-3D";
	case pot_S2b_u_ : return "S2b_u_part";
	case pot_S2c_u_ : return "S2c_u_part";
	case pot_S4a_u_ : return "S4a_u_part";
	case pot_S4b_u_ : return "S4b_u_part";
	case pot_S4c_u_ : return "S4c_u_part";
	case pot_S2b_r_ : return "S2b_r_part";
	case pot_S2c_r_ : return "S2c_r_part";
	case pot_S4a_r_ : return "S4a_r_part";
	case pot_S4b_r_ : return "S4b_r_part";
	case pot_S4c_r_ : return "S4c_r_part";
	case pot_ccs_ : return "ccs-potential";
	}
	return "undefined";
}

static std::string assign_name(const potentialtype_d &inp){
	switch(inp){
	case pot_F6D_ : return "Fock-Residue-6D";
	case pot_cc2_coulomb_ : return "CC2-Coulomb";
	case pot_cc2_residue_ : return "CC2-Residue";
	}
	return "undefined";
}

static std::string assign_name(const functype &inp){
	switch(inp){
	case HOLE : return "Hole";
	case PARTICLE : return "Particle";
	case MIXED : return "Mixed";
	case RESPONSE: return "Response";
	case UNDEFINED : return "Undefined";
	}
	return "???";
}


typedef std::vector<Function<double, 3> > vecfuncT;

// Timer Structure
struct CC_Timer{
	/// TDA_TIMER contructor
	/// @param[in] world the world
	/// @param[in] msg	a string that contains the desired printout when info function is called
	CC_Timer(World &world,std::string msg) : world(world),start_wall(wall_time()),start_cpu(cpu_time()),operation(msg),end_wall(0.0), end_cpu(0.0) {}
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
			if(world.rank()==0) std::cout<< std::setw(17) << std::setfill(' ') << std::setw(60) << "Timer:"+operation+" : "<< std::setfill(' ') << std::scientific << std::setprecision(1)
			<< end_wall << "s (wall) "<< end_cpu << "s (cpu)" << std::endl;
		}
	}

	std::pair<double,double> current_time(bool printout = false){
		update_time();
		info(printout);
		return std::make_pair(end_wall,end_cpu);
	}

	void print(const std::pair<double,double> &times)const{
		if(world.rank()==0) std::cout<< std::setw(20) << std::setfill(' ')  << "Timer: " << std::setw(60)<< operation+" : "<< std::setfill(' ') << std::scientific << std::setprecision(1)
		<< times.first << "s (wall) "<< times.second << "s (cpu)" << std::endl;
	}


};

struct CC_Parameters{
	// default constructor
//	CC_Parameters():
//	{return CC_Parameters("default constructor")}

	const double uninitialized = 123.456;

	// read parameters from input
	/// ctor reading out the input file
	CC_Parameters(const std::string& input,const double &low) :
		calculation(CISpD_),
		lo(uninitialized),
		thresh_3D(uninitialized),
		tight_thresh_3D(uninitialized),
		thresh_6D(uninitialized),
		tight_thresh_6D(uninitialized),
		thresh_bsh_3D(uninitialized),
		thresh_bsh_6D(uninitialized),
		thresh_poisson(uninitialized),
		thresh_f12(uninitialized),
		thresh_Ue(uninitialized),
		econv(uninitialized),
		dconv_3D(uninitialized),
		dconv_6D(uninitialized),
		iter_max_3D(10),
		iter_max_6D(10),
		restart(false),
		no_compute(false),
		corrfac_gamma(1.0),
		output_prec(8),
		debug(false),
		kain(false),
		freeze(0),
		test(false)
	{
		// get the parameters from the input file
		std::ifstream f(input.c_str());
		position_stream(f, "cc2");
		std::string s;

		// general operators thresh
		double thresh_operators=uninitialized;
		double thresh_operators_3D=uninitialized;
		double thresh_operators_6D=uninitialized;

		while (f >> s) {
			//std::cout << "input tag is: " << s << std::endl;
			std::transform(s.begin(),s.end(),s.begin(), ::tolower);
			//std::cout << "transformed input tag is: " << s << std::endl;
			if (s == "end") break;
			else if (s == "calculation"){
			  std::string tmp;
			  f>>tmp;
			  calculation = assign_calctype(tmp);
			}
			else if (s == "lo") f >> lo;
			else if (s == "thresh") f >> thresh_6D;
			else if (s == "thresh_3d") f >> thresh_3D;
			else if (s == "tight_thresh_3d") f >> tight_thresh_3D;
			else if (s == "thresh_6d") f >> thresh_6D;
			else if (s == "tight_thresh_6d") f >> tight_thresh_6D;
			else if (s == "debug") debug = true;
			else if (s == "econv")f >> econv;
			else if (s == "dconv") f >> dconv_6D;
			else if (s == "dconv_3d")f >> dconv_3D;
			else if (s == "dconv_6d")f >> dconv_6D;
			else if (s == "thresh_operators" or s == "thresh_operator") f>> thresh_operators;
			else if (s == "thresh_operators_3d" or s == "thresh_operator_3d") f >> thresh_operators_3D;
			else if (s == "thresh_operators_6d" or s == "thresh_operator_6d") f >> thresh_operators_6D;
			else if (s == "thresh_bsh_3d") f >> thresh_bsh_3D;
			else if (s == "thresh_bsh_6d") f >> thresh_bsh_6D;
			else if (s == "thresh_poisson") f >> thresh_poisson;
			else if (s == "thresh_f12") f >> thresh_f12;
			else if (s == "thresh_ue") f >> thresh_Ue;
			else if (s == "freeze") f >> freeze;
			else if (s == "iter_max_3d") f >> iter_max_3D;
			else if (s == "iter_max_6d") f >> iter_max_6D;
			else if (s == "restart") restart=true;
			else if (s == "no_compute") no_compute=true;
			else if (s == "kain") kain=true;
			else if (s == "kain_subspace") f>>kain_subspace;
			else if (s == "freeze") f>>freeze;
			else if (s == "test") test =true;
			else if (s == "corrfac") f>>corrfac_gamma;
			else{
			  std::cout << "Unknown Keyword: " << s << "\n";
			  continue;
			}
		}

		// set defaults
		if(not kain) kain_subspace = 0;

		// set all parameters that were not explicitly given
		if(lo==uninitialized) lo = 1.e-7;
		if(thresh_6D==uninitialized) thresh_6D = 1.e-3;
		if(tight_thresh_6D==uninitialized) tight_thresh_6D = thresh_6D*0.1;
		if(thresh_3D==uninitialized) thresh_3D = thresh_6D*0.01;
		if(tight_thresh_3D==uninitialized) tight_thresh_3D = thresh_3D*0.1;
		if(thresh_operators==uninitialized) thresh_operators = 1.e-6;
		if(thresh_operators_3D==uninitialized) thresh_operators_3D = thresh_operators;
		if(thresh_operators_6D==uninitialized) thresh_operators_6D = thresh_operators;
		if(thresh_bsh_3D==uninitialized) thresh_bsh_3D = thresh_operators_3D;
		if(thresh_bsh_6D==uninitialized) thresh_bsh_6D = thresh_operators_6D;
		if(thresh_poisson==uninitialized) thresh_poisson = thresh_operators_3D;
		if(thresh_f12==uninitialized) thresh_f12 = thresh_operators_3D;
		if(thresh_Ue==uninitialized) thresh_Ue = tight_thresh_6D;
		if(dconv_6D==uninitialized) dconv_6D = thresh_6D;
		if(dconv_3D==uninitialized) dconv_3D = dconv_6D;
		if(econv ==uninitialized) econv = 0.1*dconv_6D;
		if(iter_max_6D==uninitialized) iter_max_6D = 10;
		if(iter_max_3D==uninitialized) iter_max_3D = iter_max_6D;

		// set the thresholds
		FunctionDefaults<3>::set_thresh(thresh_3D);
		FunctionDefaults<6>::set_thresh(thresh_6D);
		if(econv < 1.e-1) output_prec = 2;
		if(econv < 1.e-2) output_prec = 3;
		if(econv < 1.e-3) output_prec = 4;
		if(econv < 1.e-4) output_prec = 5;
		if(econv < 1.e-5) output_prec = 6;
		if(econv < 1.e-6) output_prec = 7;
		std::cout.precision(output_prec);
	}


	// the demanded calculation: possibilities are MP2_, CC2_, CIS_, CCS_ (same as CIS), CISpD_
	calctype calculation;
	double lo;
	// function thresh 3D
	double thresh_3D;
	double tight_thresh_3D;
	// function thresh 6D
	double thresh_6D;
	double tight_thresh_6D;
	// BSH thresh
	double thresh_bsh_3D;
	double thresh_bsh_6D;
	// Poisson thresh
	double thresh_poisson;
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
	// iterations
	size_t iter_max_3D;
	size_t iter_max_6D;
	// restart
	bool restart;
	bool no_compute;
	// Exponent for the correlation factor
	double corrfac_gamma;
	// for formated output
	size_t output_prec;
	// debug mode
	bool debug;
	// use kain
	bool kain;
	size_t kain_subspace;
	// freeze MOs
	size_t freeze;
	// Gamma of the correlation factor
	double gamma()const{
		if(corrfac_gamma<0) MADNESS_EXCEPTION("ERROR in CC_PARAMETERS: CORRFAC_GAMMA WAS NOT INITIALIZED",1);
		return corrfac_gamma;
	}
	bool test;

	// print out the parameters
	void information(World &world)const{
		if(world.rank()==0){
//			std::cout << "Defaults for 6D and 3D Functions:\n";
//			FunctionDefaults<3>::print();
//			FunctionDefaults<6>::print();
			std::cout << "THE DEMANDED CALCULATION IS " << assign_name(calculation) << std::endl;
			if(no_compute) warning(world,"no computation demanded");
			std::cout << "\n\nThe" <<  assign_name(calculation) << " Parameters are:\n";
			std::cout << std::setw(20) << std::setfill(' ') << "Corrfac. Gamma :"           << corrfac_gamma << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "freeze :"           << freeze << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "restart :"           << restart << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "lo :"                << lo << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "k (3D) :"                << FunctionDefaults<3>::get_k() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "k (6D) :"                << FunctionDefaults<6>::get_k() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_3D demanded :"         << thresh_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_3D set :"         << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_6D demanded :"         << thresh_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_6D set :"         << FunctionDefaults<6>::get_thresh() << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "tight_thresh_6D :"     << tight_thresh_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "tight_thresh_3D :"     << tight_thresh_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_bsh_3D :"     << thresh_bsh_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_bsh_6D :"     << thresh_bsh_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_poisson :" << thresh_poisson << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_f12 :"        << thresh_f12 << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "thresh_Ue :"        << thresh_Ue << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "econv :"             << econv << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "dconv_3D :"          << dconv_3D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "dconv_6D :"          << dconv_6D << std::endl;
			std::cout << std::setw(20) << std::setfill(' ') << "freeze :"           << freeze << std::endl;
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
			if(kain) std::cout << std::setw(20) << std::setfill(' ') << "Kain subspace: " << kain_subspace << std::endl;
			if(test) std::cout << "\n\n\t\t\t!Test Mode is on!\n\n" << std::endl;
		}
	}

	void sanity_check(World &world)const{
		size_t warnings = 0;
		if(FunctionDefaults<3>::get_thresh() > 0.01*FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"3D Thresh is too low, should be 0.01*6D_thresh");
		if(FunctionDefaults<3>::get_thresh() > 0.1*FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"3D Thresh is way too low, should be 0.01*6D_thresh");
		if(FunctionDefaults<3>::get_cell_min_width() != FunctionDefaults<6>::get_cell_min_width()) warnings+=warning(world,"3D and 6D Cell sizes differ");
		if(FunctionDefaults<3>::get_k() != FunctionDefaults<6>::get_k()) warnings+=warning(world, "k-values of 3D and 6D differ ");
		if(FunctionDefaults<3>::get_truncate_mode()!=3) warnings+=warning(world,"3D Truncate mode is not 3");
		if(FunctionDefaults<6>::get_truncate_mode()!=3) warnings+=warning(world,"6D Truncate mode is not 3");
		if(dconv_3D < FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher convergence than threshold for 3D");
		if(dconv_6D < FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher convergence than threshold for 6D");
		if(thresh_3D != FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"3D thresh set unequal 3D thresh demanded");
		if(thresh_6D != FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"6D thresh set unequal 6D thresh demanded");
		if(econv < FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 3D");
		if(econv < FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 6D");
		if(econv < 0.1*FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 3D (more than factor 10 difference)");
		if(econv < 0.1*FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 6D (more than factor 10 difference)");
		// Check if the 6D thresholds are not too high
		if(thresh_6D < 1.e-3) warnings+=warning(world,"thresh_6D is smaller than 1.e-3");
		if(thresh_6D < tight_thresh_6D) warnings+=warning(world,"tight_thresh_6D is larger than thresh_6D");
		if(thresh_6D < tight_thresh_3D) warnings+=warning(world,"tight_thresh_3D is larger than thresh_3D");
		if(thresh_6D < 1.e-3) warnings+=warning(world,"thresh_6D is smaller than 1.e-3");
		if(thresh_Ue < 1.e-4) warnings+=warning(world,"thresh_Ue is smaller than 1.e-4");
		if(thresh_Ue > 1.e-4) warnings+=warning(world,"thresh_Ue is larger than 1.e-4");
		if(thresh_3D > 0.01*thresh_6D) warnings+=warning(world,"Demanded 6D thresh is to precise compared with the 3D thresh");
		if(thresh_3D > 0.1*thresh_6D) warnings+=warning(world,"Demanded 6D thresh is to precise compared with the 3D thresh");
		if(kain and kain_subspace ==0) warnings+=warning(world,"Demanded Kain solver but the size of the iterative subspace is set to zero");
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

/// enhanced POD for the pair functions
class CC_Pair: public archive::ParallelSerializableObject {

public:

	/// default ctor; initialize energies with a large number
	CC_Pair(const pairtype type) : type(type),
		i(-1), j(-1), e_singlet(uninitialized()), e_triplet(
				uninitialized()), ij_gQf_ij(uninitialized()), ji_gQf_ij(
						uninitialized()), converged(false), current_error(uninitialized()), current_energy_difference(uninitialized()), epsilon(uninitialized()){
	}

	/// ctor; initialize energies with a large number
	CC_Pair(const int i, const int j,const pairtype type) : type(type),
		i(i), j(j), e_singlet(uninitialized()), e_triplet(uninitialized()), ij_gQf_ij(
				uninitialized()), ji_gQf_ij(uninitialized()), converged(
						false), current_error(uninitialized()),epsilon(uninitialized()){
	}
	/// ctor; initialize energies with a large number
	CC_Pair(const real_function_6d &f,const int i, const int j,const pairtype type) : type(type),
		i(i), j(j),function(f), e_singlet(uninitialized()), e_triplet(uninitialized()), ij_gQf_ij(
				uninitialized()), ji_gQf_ij(uninitialized()), converged(
						false), current_error(uninitialized()),epsilon(uninitialized()) {
	}

	/// print the pair's energy
	void print_energy() const {
		if (function.world().rank() == 0) {
			printf("final correlation energy %2ld %2ld %12.8f %12.8f\n", i, j,
					e_singlet, e_triplet);
		}
	}

	// print information
	void info()const{
		if(function.world().rank()==0){
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " Current Information about Electron Pair " << name() << std::endl;
			//std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ij_gQf_ij: " << ij_gQf_ij << std::endl;
			//std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ji_gQf_ij: " << ji_gQf_ij << std::endl;
			if(function.impl_initialized()) std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ||u||    : "      <<std::setprecision(4)<<std::scientific<< function.norm2() << std::endl;
			if(constant_term.impl_initialized()) std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " ||const||: " <<std::setprecision(4)<<std::scientific<< constant_term.norm2() << std::endl;
			if(current_error != uninitialized()) std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " |error|  : " << current_error << std::endl;
			if(current_energy_difference != uninitialized()) std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " |deltaE|  : " << current_energy_difference << std::endl;
			if(current_energy != uninitialized()) std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << "  omega   : " <<std::setprecision(10)<<std::fixed<< current_energy << std::endl;
			//if(epsilon == uninitialized()) std::cout << "WARNING: BSH-epsilon is not initialized" << std::endl;
		}
	}

	std::string name()const{
	    std::string name = "???";
	    if(type==GROUND_STATE) name ="u";
	    if(type==EXCITED_STATE) name = "x";
		return name+stringify(i)+stringify(j);
	}

	static double uninitialized() {
		return 1.e10;
	}

	const pairtype type;

	const size_t i, j;                       ///< orbitals i and j
	real_function_6d function; ///< pair function for a specific pair w/o correlation factor part
	real_function_6d constant_term;	///< the first order contribution to the MP1 wave function

	double e_singlet;				///< the energy of the singlet pair ij
	double e_triplet;				///< the energy of the triplet pair ij

	double ij_gQf_ij;         	  	///< <ij | g12 Q12 f12 | ij>
	double ji_gQf_ij;          		///< <ji | g12 Q12 f12 | ij>

	bool converged;					///< is the pair function converged

	double current_error;			///< error of the last iteration: ||function_old - function||_L2
	double current_energy_difference;/// difference of current_energy and energy of the last iteration
	double current_energy = uninitialized(); /// < the correlation energy of the last iteration
	double epsilon;					///< the summed up orbital energies corresponding to the pair function indices: epsilon_i + epsilon_j


	/// serialize this CC_Pair

	/// store the function only if it has been initialized
	/// load the function only if there is one
	/// don't serialize recomputable intermediates r12phi, Uphi, KffKphi
	template<typename Archive> void serialize(Archive& ar) {
		bool fexist = function.is_initialized();
		bool cexist = constant_term.is_initialized();
		ar & ij_gQf_ij & ji_gQf_ij & e_singlet & e_triplet & converged
		& fexist & cexist;
		if (fexist)
			ar & function;
		if (cexist)
			ar & constant_term;
	}

	bool load_pair(World& world) {
		std::string name_ = name();
		bool exists = archive::ParallelInputArchive::exists(world,
				name_.c_str());
		if (exists) {
			if (world.rank() == 0)
				printf("loading pair %s", name_.c_str());
			archive::ParallelInputArchive ar(world, name_.c_str(), 1);
			ar & *this;
			if (world.rank() == 0)
			function.set_thresh(FunctionDefaults<6>::get_thresh());
			constant_term.set_thresh(FunctionDefaults<6>::get_thresh());
		} else {
			if (world.rank() == 0) std::cout << "pair " << name_ << " not found " << std::endl;
		}
		return exists;
	}

	void store_pair(World& world, const std::string &msg = "") {
		std::string name_ = msg+name();
		if (world.rank() == 0)
			printf("storing CC_Pair %s\n", name_.c_str());
		archive::ParallelOutputArchive ar(world, name_.c_str(), 1);
		ar & *this;
	}

};


// TAKEN FROM MP2.h
/// POD holding all electron pairs with easy access
template<typename T>
struct Pairs {


	typedef std::map<std::pair<int, int>, T> pairmapT;
	pairmapT allpairs;


	/// getter
	const T & operator()(int i,int j)const{
		return allpairs.at(std::make_pair(i, j));
	}

	/// getter
	// at instead of [] operator bc [] inserts new element if nothing is found while at throws out of range error
	T& operator()(int i, int j) {
		return allpairs.at(std::make_pair(i, j));
	}

	/// setter
	void insert(int i, int j, T pair) {
		std::pair<int, int> key = std::make_pair(i, j);
		allpairs.insert(std::make_pair(key, pair));
	}
};

typedef Pairs<real_function_3d> intermediateT;
static double size_of(const intermediateT &im){
  double size=0.0;
  for(const auto & tmp:im.allpairs){
	size += get_size<double,3>(tmp.second);
  }
  return size;
}


// structure for a CC Function 3D which holds an index and a type
struct CC_function{
	CC_function(): current_error(99),i(99), type(UNDEFINED){};
	CC_function(const real_function_3d &f): current_error(99),function(f), i(99),type(UNDEFINED){};
	CC_function(const real_function_3d &f,const size_t &ii): current_error(99), function(f), i(ii), type(UNDEFINED){};
	CC_function(const real_function_3d &f,const size_t &ii, const functype &type_): current_error(99),function(f), i(ii), type(type_){};
	CC_function(const CC_function &other): current_error(other.current_error),function(other.function), i(other.i), type(other.type){};
	double current_error;
	real_function_3d function;
	real_function_3d get()const{return function;}
	real_function_3d f()const{return function;}
	void set(const real_function_3d &other){function=other;}
	size_t i;
	functype type;
	void info(World &world,const std::string &msg = " ")const{
		if(world.rank()==0){
			std::cout <<"Information about 3D function: " << name() << " " << msg << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " |f|    : " << function.norm2() << std::endl;
			std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " |error|: " << current_error << std::endl;
		}
	}
	std::string name()const{
		if(type==HOLE){return "phi"+stringify(i);}
		else if(type==PARTICLE){return "tau"+stringify(i);}
		else if(type==MIXED){return "t"+stringify(i);}
		else if(type==RESPONSE){return "x"+stringify(i);}
		else{return "function"+stringify(i);}
	}
	double inner(const CC_function &f)const{
		return inner(f.function);
	}
	double inner(const real_function_3d &f)const{
		return function.inner(f);
	}

	CC_function operator*(const CC_function &f)const{
		real_function_3d product = function*f.function;
		return CC_function(product,999,UNDEFINED);
	}
	CC_function operator+(const CC_function &f)const{
		real_function_3d sum = function+f.function;
		return CC_function(sum,i,combine_types(f));
	}

	functype combine_types(const CC_function &f)const{
		if(type == UNDEFINED or f.type == UNDEFINED) return UNDEFINED;
		if(i==f.i){
			if(type == f.type) return type;
			else return MIXED;
		}
		else return UNDEFINED;
	}


};


// structure for CC Vectorfunction
struct CC_vecfunction{

	CC_vecfunction(): type(UNDEFINED){}
	CC_vecfunction(const functype type_): type(type_){}
	CC_vecfunction(const vecfuncT &v,const functype &type): type(type){
		for(size_t i=0;i<v.size();i++){
			CC_function tmp(v[i],i,type);
			functions.insert(std::make_pair(i,tmp));
		}
	}
	CC_vecfunction(const vecfuncT &v,const functype &type,const size_t &freeze): type(type){
		for(size_t i=0;i<v.size();i++){
			CC_function tmp(v[i],freeze+i,type);
			functions.insert(std::make_pair(freeze+i,tmp));
		}
	}
	CC_vecfunction(const std::vector<CC_function> &v,const functype type_): type(type_){
		for(auto x:v){
			functions.insert(std::make_pair(x.i,x));
		}
	}
	CC_vecfunction(const CC_vecfunction &other) : functions(other.functions),type(other.type) {}

	typedef std::map<std::size_t, CC_function> CC_functionmap;
	CC_functionmap functions;

	functype type;
	std::string name()const{
	  if (type==PARTICLE) return "singles_gs";
	  else if(type==HOLE) return "mos_gs";
	  else if(type==MIXED) return "t_gs";
	  else if(type==RESPONSE) return "singles_response";
	  else return "UNKNOWN";
	}

	/// getter
	const CC_function& operator()(const CC_function &i) const {
		return functions.find(i.i)->second;
	}

	/// getter
	const CC_function& operator()(const size_t &i) const {
		return functions.find(i)->second;
	}

	/// getter
	CC_function& operator()(const CC_function &i) {
		return functions[i.i];
	}

	/// getter
	CC_function& operator()(const size_t &i) {
		return functions[i];
	}

	/// setter
	void insert(const size_t &i, const CC_function &f) {
		functions.insert(std::make_pair(i, f));
	}

	vecfuncT get_vecfunction()const{
		vecfuncT tmp;
		for(auto x:functions) tmp.push_back(x.second.function);
		return tmp;
	}

	size_t size()const{
		return functions.size();
	}

	void print_size(const std::string &msg="!?not assigned!?")const{
		if(functions.size()==0){
			std::cout << "CC_vecfunction " << msg << " is empty\n";
		}else{
			std::string msg2;
			if(msg=="!?not assigned!?") msg2 = "";
			else msg2 = "_("+msg+")";
			for(auto x:functions){
				x.second.function.print_size(x.second.name()+msg2);
			}
		}
	}

};

// data structure which contains information about performances of a functions
struct CC_data{
	CC_data(): name("UNDEFINED"), time(std::make_pair(999.999,999.999)), result_size(999.999), result_norm(999.999){}
	CC_data(const std::string &name_):name(name_), time(std::make_pair(999.999,999.999)), result_size(999.999), result_norm(999.999){}
	CC_data(const potentialtype_s &name_):name(assign_name(name_)), time(std::make_pair(999.999,999.999)), result_size(999.999), result_norm(999.999){}
	CC_data(const CC_data &other) : name(other.name), time(other.time), result_size(other.result_size), result_norm(other.result_norm), warnings(other.warnings){}
	const std::string name;
	std::pair<double,double> time; // overall time
	double result_size;
	double result_norm;
	std::vector<std::string> warnings;

	void info(World & world)const{
		if(world.rank()==0) info();
	}
	void info(const bool &x)const{
		if(x) info();
		else return;
	}
	void info()const{
		std::cout << std::setw(25) <<name << std::setfill(' ') << ", ||f||=" << result_norm << ", (" << result_size << ") GB, " << time.first << "s (Wall), " << time.second << "s (CPU)\n";
		if(not warnings.empty()){
			std::cout << "!!!Problems were detected in " << name <<"!!!\n";
			std::cout << warnings << std::endl;
		}
	}
};

// structure which holds all CC_data structures sorted by name of the function and iteration
struct CC_performance{

	CC_performance():current_iteration(0){}

	typedef std::map<std::pair<std::string, std::size_t>, CC_data> datamapT;
	datamapT data;

	/// getter
	const CC_data& operator()(const std::string &name, const size_t &iter) const {
		return data.find(std::make_pair(name, iter))->second;
	}

	/// getter
	const CC_data& operator()(const std::string &name) const {
		return data.find(std::make_pair(name, current_iteration))->second;
	}


	/// getter
	CC_data& operator()(const std::string &name, const std::size_t &iter) {
		return data[std::make_pair(name, iter)];
	}

	/// getter
	CC_data& operator()(const std::string &name) {
		return data[std::make_pair(name, current_iteration)];
	}

	/// setter
	void insert(const std::string &name, const CC_data &new_data) {
		std::pair<std::string, std::size_t> key = std::make_pair(name, current_iteration);
		data.insert(std::make_pair(key, new_data));
	}

	mutable std::size_t current_iteration;

	void info()const{
		std::cout << "CC2 Performance information: Iteration \n";
		for(auto x:data) x.second.info();
	}

	void info(const std::size_t &iter)const{
		std::cout << "CC2 Performance information: Iteration" << iter << "\n";
		for(auto x:data){
			if(x.first.second == iter) x.second.info();
		}
	}

	void info_last_iter()const {
		if(current_iteration !=0)info(current_iteration -1);
		else info(0);
	}

	std::pair<double,double> get_average_time(const std::string &name)const{
		double overall_time_cpu = 0.0;
		double overall_time_wall = 0.0;
		size_t iterations = 0;
		for(auto x:data){
			if(x.first.first == name){
				overall_time_wall += x.second.time.first;
				overall_time_cpu += x.second.time.second;
				iterations++;
			}
		}
		if(iterations==0) return std::make_pair(0.0,0.0);
		double iter = (double) iterations;
		return std::make_pair(overall_time_wall/iter,overall_time_cpu/iter);
	}
};

}//namespace madness

#endif /* CCSTRUCTURES_H_ */
