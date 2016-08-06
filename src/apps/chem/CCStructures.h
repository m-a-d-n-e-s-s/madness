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
#include "electronic_correlation_factor.h"
#include <algorithm> // tolower function for strings
#include <examples/nonlinsol.h>


namespace madness{

  enum functype_6d {pure_, decomposed_, op_decomposed_};
  enum optype {g12_,f12_};
  enum calctype {MP2_, CC2_, LRCCS_, LRCC2_, CISpD_,ADC2_ , experimental_};
  enum functype {HOLE,PARTICLE,MIXED,RESPONSE,UNDEFINED};
  enum pairtype {GROUND_STATE,EXCITED_STATE};
  enum potentialtype_s {pot_F3D_, pot_s2b_, pot_s2c_, pot_s4a_, pot_s4b_, pot_s4c_, pot_S2b_u_, pot_S2c_u_, pot_S4a_u_, pot_S4b_u_, pot_S4c_u_,pot_S2b_r_, pot_S2c_r_, pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_ccs_,pot_cis_,pot_singles_};
  enum potentialtype_d {pot_F6D_, pot_cc2_coulomb_,pot_cc2_residue_};
  // The pair function is:  \tau = u + Qf(|titj>), FULL means that \tau is calculated in 6D form, DECOMPOSED means that u is used in 6D and the rest is tried to solve in 3D whenever possible
  enum pair_function_form{DECOMPOSED, FULL};

  static std::string assign_name(const pairtype &input){
    switch(input){
      case GROUND_STATE : return "Ground State";
      case EXCITED_STATE: return "Excited State";
    }
    MADNESS_EXCEPTION("assign_name:pairtype, should not end up here",1);
    return "unknown pairtype";
  }

  static std::string assign_name(const optype &input){
    switch(input){
      case g12_ : return "g12";
      case f12_ : return "f12";
    }
    MADNESS_EXCEPTION("assign_name:optype, should not end up here",1);
    return "unknown operatortype";
  }

  static calctype assign_calctype(const std::string name){
    if(name=="mp2") return MP2_;
    else if(name=="cc2") return CC2_;
    else if(name=="lrcc2" or name=="cc2_response") return LRCC2_;
    else if(name=="cispd") return CISpD_;
    else if(name=="cis" or name=="ccs" or name=="ccs_response" or name=="lrccs") return LRCCS_;
    else if(name=="experimental") return experimental_;
    else if(name=="adc2" or name=="adc(2)") return ADC2_;
    else{
      std::string msg= "CALCULATION OF TYPE: " + name + " IS NOT KNOWN!!!!";
      MADNESS_EXCEPTION(msg.c_str(),1);
    }
  }
  static std::string assign_name(const calctype &inp){
    switch(inp){
      case CC2_ : return "CC2";
      case MP2_ : return "MP2";
      case LRCC2_ : return "LRCC2";
      case CISpD_ : return "CISpD";
      case LRCCS_: return "LRCCS";
      case ADC2_: return "ADC2";
      case experimental_: return "experimental";
    }
    return "unknown";
  }
  static std::string assign_name(const potentialtype_s &inp){
    switch(inp){
      case pot_F3D_ : return "F3D";
      case pot_s2b_ : return "s2b";
      case pot_s2c_ : return "s2c";
      case pot_s4a_ : return "s4a";
      case pot_s4b_ : return "s4b";
      case pot_s4c_ : return "s4c";
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
      case pot_ccs_ : return "ccs";
      case pot_cis_ : return "cis-potential";
      case pot_singles_: return "singles potential";
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
    CC_Timer(World &world,std::string msg) : world(world),start_wall(wall_time()),start_cpu(cpu_time()),operation(msg),end_wall(0.0), end_cpu(0.0), time_wall(-1.0), time_cpu(-1.0) {}
    World & world;
    double start_wall;
    double start_cpu;
    std::string operation;
    double end_wall;
    double end_cpu;
    double time_wall;
    double time_cpu;
    void update_time(){
      time_wall = wall_time()-start_wall;
      time_cpu = cpu_time()-start_cpu;
    }
  public:
    /// print out information about the passed time since the TDA_TIMER object was created
    void info(const bool debug = true, const double norm=12345.6789){
      if(debug==true){
	update_time();
	std::string s_norm = "";
	if(norm!=12345.6789) s_norm=", ||result||="+std::to_string(norm);
	if(world.rank()==0){
	  std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
	  << "Timer: " << time_wall << " (Wall), " << time_cpu << " (CPU)" << s_norm << ", (" +operation+")" <<  "\n";
	}
      }
    }

    void start(){
      start_wall = wall_time();
      start_cpu  = cpu_time();
    }
    void stop(){
      end_wall = wall_time();
      end_cpu = cpu_time();
      time_wall = end_wall - start_wall;
      time_cpu = end_cpu - start_cpu;
    }


    double get_wall_time_diff()const{return end_wall;}
    double get_cpu_time_diff()const{return end_cpu;}

    std::pair<double,double> current_time(bool printout = false){
      if(time_wall<0.0 or time_cpu<0.0) stop();
      return std::make_pair(time_wall,time_cpu);
    }

    double current_wall(){return current_time().first;}
    double current_cpu(){return current_time().second;}

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

    // copy constructor
    CC_Parameters(const CC_Parameters& other) :
      calculation(other.calculation),
      lo(other.lo),
      dmin(other.dmin),
      thresh_3D(other.thresh_3D),
      tight_thresh_3D(other.tight_thresh_3D),
      thresh_6D(other.thresh_6D),
      tight_thresh_6D(other.tight_thresh_6D),
      thresh_bsh_3D(other.thresh_bsh_3D),
      thresh_bsh_6D(other.thresh_bsh_6D),
      thresh_poisson(other.thresh_poisson),
      thresh_f12(other.thresh_f12),
      thresh_Ue(other.thresh_Ue),
      econv(other.econv),
      econv_pairs(other.econv_pairs),
      dconv_3D(other.dconv_3D),
      dconv_6D(other.dconv_6D),
      iter_max_3D(other.iter_max_3D),
      iter_max_6D(other.iter_max_6D),
      restart(other.restart),
      no_compute(other.no_compute),
      no_compute_gs(other.no_compute_gs),
      no_compute_response(other.no_compute_response),
      no_compute_mp2(other.no_compute_response),
      no_compute_cc2(other.no_compute_mp2),
      no_compute_cispd(other.no_compute_cispd),
      no_compute_lrcc2(other.no_compute_lrcc2),
      corrfac_gamma(other.corrfac_gamma),
      output_prec(other.output_prec),
      debug(other.debug),
      kain(other.kain),
      freeze(other.freeze),
      test(other.test),
      decompose_Q(other.decompose_Q),
      QtAnsatz(other.QtAnsatz),
      excitations_(other.excitations_),
      tda_guess_mode("uninitialized"),
      tda_excitations(0),
      tda_guess_excitations(0),
      tda_iterating_excitations(0),
      tda_guess("uninitialized"),
      tda_energy_guess_factor(uninitialized),
      tda_dconv_guess(uninitialized),
      tda_dconv(uninitialized),
      tda_dconv_hard(uninitialized),
      tda_econv_guess(uninitialized),
      tda_econv(uninitialized),
      tda_econv_hard(uninitialized)
    {}

    // read parameters from input
    /// ctor reading out the input file
    CC_Parameters(const std::string& input,const double &low) :
      calculation(LRCC2_),
      lo(uninitialized),
      dmin(1.0),
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
      econv_pairs(uninitialized),
      dconv_3D(uninitialized),
      dconv_6D(uninitialized),
      iter_max_3D(10),
      iter_max_6D(10),
      restart(false),
      no_compute(false),
      no_compute_gs(false),
      no_compute_response(false),
      no_compute_mp2(false),
      no_compute_cc2(false),
      no_compute_cispd(false),
      no_compute_lrcc2(false),
      corrfac_gamma(1.0),
      output_prec(8),
      debug(false),
      kain(false),
      freeze(0),
      test(false),
      decompose_Q(false),
      QtAnsatz(false),
      excitations_(0),
      tda_guess_mode("uninitialized"),
      tda_excitations(0),
      tda_guess_excitations(0),
      tda_iterating_excitations(0),
      tda_guess("uninitialized"),
      tda_energy_guess_factor(uninitialized),
      tda_dconv_guess(uninitialized),
      tda_dconv(uninitialized),
      tda_dconv_hard(uninitialized),
      tda_econv_guess(uninitialized),
      tda_econv(uninitialized),
      tda_econv_hard(uninitialized)
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
	else if (s == "dmin") f>>dmin;
	else if (s == "thresh") f >> thresh_6D;
	else if (s == "thresh_3d") f >> thresh_3D;
	else if (s == "tight_thresh_3d") f >> tight_thresh_3D;
	else if (s == "thresh_6d") f >> thresh_6D;
	else if (s == "tight_thresh_6d") f >> tight_thresh_6D;
	else if (s == "debug") debug = true;
	else if (s == "econv")f >> econv;
	else if (s == "econv_pairs")f >> econv_pairs;
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
	else if (s == "kain") kain=true;
	else if (s == "kain_subspace") f>>kain_subspace;
	else if (s == "freeze") f>>freeze;
	else if (s == "test") test =true;
	else if (s == "corrfac" or s=="corrfac_gamma" or s=="gamma") f>>corrfac_gamma;
	else if (s == "decompose_q") decompose_Q=true;
	else if (s == "restart")restart = true;
	else if (s == "no_compute"){
	  no_compute = true;
	  no_compute_gs=true;
	  no_compute_mp2=true;
	  no_compute_cispd=true;
	  no_compute_response=true;
	}
	else if (s == "no_compute_gs"){
	  no_compute_gs=true;
	  no_compute_mp2=true;
	  no_compute_cc2=true;
	}
	else if (s == "no_compute_response"){
	  no_compute_response=true;
	  no_compute_cispd=true;
	  no_compute_lrcc2=true;
	}
	else if (s == "no_compute_cc2") no_compute_cc2=true;
	else if (s == "no_compute_cispd") no_compute_cispd=true;
	else if (s == "no_compute_lrcc2") no_compute_lrcc2=true;
	else if (s == "no_compute_mp2") no_compute_mp2=true;
	else if (s == "excitation"){
	  size_t tmp;
	  f>>tmp;
	  excitations_.push_back(tmp);
	}
	else if ( s == "qtansatz") QtAnsatz=true;
	else if ( s == "tda_guess_mode") f>>tda_guess_mode;
	else if ( s == "tda_guess_excitations") f>>tda_guess_excitations;
	else if ( s == "tda_excitations") f>>tda_excitations;
	else if ( s == "tda_iterating_excitations") f>>tda_iterating_excitations;
	else if ( s == "tda_guess") f >> tda_guess;
	else if ( s == "tda_energy_guess_factor") f >> tda_energy_guess_factor;
	else if ( s == "tda_dconv_guess") f >> tda_dconv_guess;
	else if ( s == "tda_dconv") f >> tda_dconv;
	else if ( s == "tda_dconv_hard")f >> tda_dconv_hard;
	else if ( s == "tda_econv_guess") f >> tda_econv_guess;
	else if ( s == "tda_econv") f >> tda_econv;
	else if ( s == "tda_econv_hard") f >> tda_econv_hard;
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
      if(econv_pairs ==uninitialized) econv_pairs = econv;
      if(iter_max_6D==uninitialized) iter_max_6D = 10;
      if(iter_max_3D==uninitialized) iter_max_3D = iter_max_6D;

      // set the thresholds
      FunctionDefaults<3>::set_thresh(thresh_3D);
      FunctionDefaults<6>::set_thresh(thresh_6D);
      if(thresh_3D < 1.1e-1) output_prec = 3;
      if(thresh_3D < 1.1e-2) output_prec = 4;
      if(thresh_3D < 1.1e-3) output_prec = 5;
      if(thresh_3D < 1.1e-4) output_prec = 6;
      if(thresh_3D < 1.1e-5) output_prec = 7;
      if(thresh_3D < 1.1e-6) output_prec = 8;
      std::cout.precision(output_prec);

      // set the default TDA parameters
      if(tda_guess=="uninitialized") tda_guess = "big_fock_3";
      if(tda_guess_mode=="uninitialized") tda_guess_mode = "projected";
      if(tda_excitations==0) tda_excitations = 4;
      if(tda_guess_excitations==0) tda_guess_excitations = tda_excitations;
      if(tda_iterating_excitations==0) tda_iterating_excitations=tda_guess_excitations;
      if(tda_energy_guess_factor==uninitialized) tda_energy_guess_factor=0.99;
      if(tda_dconv_guess==uninitialized) tda_dconv_guess = 1.0;
      if(tda_dconv==uninitialized) tda_dconv = 1.0;
      if(tda_dconv_hard==uninitialized) tda_dconv_hard = 10.0*thresh_3D;
      if(tda_econv_guess==uninitialized) tda_econv_guess = 1.e-1;
      if(tda_econv==uninitialized) tda_econv =1.e-2;
      if(tda_econv_hard==uninitialized) tda_econv_hard = thresh_3D;

      if(no_compute==true and restart ==false) restart = true;
    }


    // the demanded calculation: possibilities are MP2_, CC2_, CIS_, CCS_ (same as CIS), CISpD_
    calctype calculation;
    double lo;
    // the finest length to be resolved by 6D operators which needs special refinement
    // this will define the depth of the special level (default is 1.0 bohr)
    double dmin;
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
    // Convergence for Correlation Energy (overall and pairs)
    double econv;
    double econv_pairs;
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
    bool no_compute_gs;
    bool no_compute_response;
    bool no_compute_mp2;
    bool no_compute_cc2;
    bool no_compute_cispd;
    bool no_compute_lrcc2;
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
    // choose if Q for the constant part of MP2 and related calculations should be decomposed: GQV or GV - GO12V
    bool decompose_Q;
    // if true the ansatz for the CC2 ground state pairs is |tau_ij> = |u_ij> + Qtf12|titj>, with Qt = Q - |tau><phi|
    // if false the ansatz is the same with normal Q projector
    // the response ansatz is the corresponding response of the gs ansatz
    bool QtAnsatz;

    /// a vector containing the excitations which shall be optizmized later (with CIS(D) or CC2)
    std::vector<size_t> excitations_;

    // Parameters for the TDA Algorithm

    /// Guess mode for TDA:
    /// "numerical" use the std numerical occupied orbitals to create the guess for the excited states
    /// "projected" use the projected occupied orbitals (projected to guess gauss basis) and avoid noise for high guess polynomials
    std::string tda_guess_mode;

    /// The number of excitation vectors for which the alorithm will solve
    std::size_t tda_excitations;
    /// The number of guess_excitation vectors for the first iterations
    std::size_t tda_guess_excitations;
    /// The number of excitation vectors which will be iterated parallel
    std::size_t tda_iterating_excitations;

    /// the guess which will be applied
    /// see the file guess.h
    /// can be "dipole", "dipole+", "quadrupole", "qualdrupole+" , " big_fock_3", "big_fock_4"
    std::string tda_guess;

    /// the guess factor for the first energy guess which is: omega = - factor*HOMO
    /// the factor has to be between ]0,1[
    double tda_energy_guess_factor;

    /// convergence for the excitation vectors in the guess, solve and solve_sequential mode
    double tda_dconv_guess;
    double tda_dconv;
    double tda_dconv_hard;
    /// convergence for the excitation energy in the guess, solve and solve_sequential mode
    double tda_econv_guess;
    double tda_econv;
    double tda_econv_hard;


    // print out the parameters
    // the TDA parameters are printed out in the TDA section of the program, so no need here
    void information(World &world)const{
      if(world.rank()==0){
	//			std::cout << "Defaults for 6D and 3D Functions:\n";
	//			FunctionDefaults<3>::print();
	FunctionDefaults<6>::print();
	std::cout << "THE DEMANDED CALCULATION IS " << assign_name(calculation) << std::endl;
	if(no_compute) warning(world,"no computation demanded");
	std::cout << "The Ansatz for the Pair functions |tau_ij> is: ";
	if(QtAnsatz) std::cout << "(Qt)f12|titj> and response: (Qt)f12(|tixj> + |xitj>) - (OxQt + QtOx)f12|titj>";
	else std::cout << "Qf12|titj> and response: Qf12(|xitj> + |tixj>)";
	std::cout << "\n\nThe" <<  assign_name(calculation) << " Parameters are:\n";
	std::cout << std::setw(20) << std::setfill(' ') << "Corrfac. Gamma :"           << corrfac_gamma << std::endl;
	std::cout << std::setw(20) << std::setfill(' ') << "freeze :"           << freeze << std::endl;
	std::cout << std::setw(20) << std::setfill(' ') << "restart :"           << restart << std::endl;
	std::cout << std::setw(20) << std::setfill(' ') << "lo :"                << lo << std::endl;
	std::cout << std::setw(20) << std::setfill(' ') << "dmin :"                << dmin << std::endl;
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
	std::cout << std::setw(20) << std::setfill(' ') << "econv_pairs :"             << econv_pairs << std::endl;
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
	if(restart) std::cout << "restart is on";
	if(no_compute) std::cout << "no_compute for all";
	if(no_compute_gs) std::cout << "no_compute for ground-state";
	if(no_compute_response) std::cout << "no_compute for excited-state";
	std::cout << "Excitations to optimize are:\n";
	if(excitations_.empty()) std::cout << "All" << "\n";
	else std::cout << excitations_ <<"\n";
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
      if(restart == false and no_compute==true){
	warnings+=warning(world,"no_compute flag detected but no restart flag");
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
    CC_Pair(const pairtype type) : type(type),i(-1), j(-1){}

    /// ctor; initialize energies with a large number
    CC_Pair(const int i, const int j,const pairtype type) : type(type),i(i), j(j){}
    /// ctor; initialize energies with a large number
    CC_Pair(const real_function_6d &f,const int i, const int j,const pairtype type) : type(type),i(i), j(j),function(f) {}


    // print information
    void info()const{
      if(function.world().rank()==0){
	std::cout <<std::setw(10) << std::setfill(' ')<<std::setw(50) << " Current Information about Electron Pair " << name() << std::endl;
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
      if(type==EXCITED_STATE) name = "chi";
      return name+stringify(i)+stringify(j);
    }

    static double uninitialized() {
      return 1.e10;
    }

    const pairtype type;

    const size_t i, j;                       ///< orbitals i and j
    real_function_6d function; ///< pair function for a specific pair w/o correlation factor part
    real_function_6d constant_term;	///< the first order contribution to the MP1 wave function

    double constant_energy= uninitialized(); /// the energy from <ij|gQf|ij> in MP2 or corresponsing terms with singles in CC2 (singles are kept constant during doubles interations)

    double current_error= uninitialized();;			///< error of the last iteration: ||function_old - function||_L2
    double current_energy_difference= uninitialized();;/// difference of current_energy and energy of the last iteration
    double current_energy = uninitialized(); /// < the correlation energy of the last iteration


    /// serialize this CC_Pair

    /// store the function only if it has been initialized
    /// load the function only if there is one
    /// don't serialize recomputable intermediates r12phi, Uphi, KffKphi
    template<typename Archive> void serialize(Archive& ar) {
      bool fexist = function.is_initialized();
      bool cexist = constant_term.is_initialized();
      ar & fexist & cexist;
      if (fexist)
	ar & function;
      if (cexist)
	ar & constant_term;
    }

    bool load_pair(World& world, const std::string &msg = "") {
      const std::string name_ = msg+name();
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

    void store_pair(World& world, const std::string &msg ="")const{
      const std::string name_ = msg+name();
      CC_Timer time(world,"Storing pair " + name_);
      if (world.rank() == 0)
	printf("storing CC_Pair %s\n", name_.c_str());
      archive::ParallelOutputArchive ar(world, name_.c_str(), 1);
      ar & *this;
      time.info();
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
    /// can NOT replace elements (for this construct new pair map and swap the content)
    void insert(int i, int j, const T& pair) {
      std::pair<int, int> key = std::make_pair(i, j);
      allpairs.insert(std::make_pair(key, pair));
    }

    /// swap the contant of the pairmap
    void swap(Pairs<T>& other){
      allpairs.swap(other.allpairs);
    }

    bool empty()const{
      if(allpairs.size()==0) return true;
      else return false;
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

    CC_vecfunction(): type(UNDEFINED),omega(0.0){}
    CC_vecfunction(const functype type_): type(type_),omega(0.0){}
    CC_vecfunction(const vecfuncT &v): type(UNDEFINED),omega(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],i,type);
	functions.insert(std::make_pair(i,tmp));
      }
    }
    CC_vecfunction(const std::vector<CC_function> &v): type(UNDEFINED),omega(0.0){
      for(size_t i=0;i<v.size();i++){
	functions.insert(std::make_pair(v[i].i,v[i]));
      }
    }
    CC_vecfunction(const vecfuncT &v,const functype &type): type(type),omega(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],i,type);
	functions.insert(std::make_pair(i,tmp));
      }
    }
    CC_vecfunction(const vecfuncT &v,const functype &type,const size_t &freeze): type(type),omega(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],freeze+i,type);
	functions.insert(std::make_pair(freeze+i,tmp));
      }
    }
    CC_vecfunction(const std::vector<CC_function> &v,const functype type_): type(type_),omega(0.0){
      for(auto x:v){
	functions.insert(std::make_pair(x.i,x));
      }
    }
    CC_vecfunction(const CC_vecfunction &other) : functions(other.functions),type(other.type), omega(other.omega) {}

    typedef std::map<std::size_t, CC_function> CC_functionmap;
    CC_functionmap functions;

    /// returns a deep copy (void shallow copy errors)
    CC_vecfunction copy()const{
      std::vector<CC_function> vn;
      for(auto x:functions){
	const CC_function fn(madness::copy(x.second.function),x.second.i,x.second.type);
	vn.push_back(fn);
      }
      return CC_vecfunction(vn,type);
    }


    functype type;
    double omega; // excitation energy
    std::string name()const{
      if (type==PARTICLE) return "tau";
      else if(type==HOLE) return "phi";
      else if(type==MIXED) return "t";
      else if(type==RESPONSE) return "x";
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
      std::cout << std::setw(25) <<name << std::setfill(' ') << ", ||f||=" <<std::fixed << std::setprecision(6) << result_norm
	  <<std::scientific << std::setprecision(2) << ", (" << result_size << ") GB, " << time.first << "s (Wall), " << time.second << "s (CPU)\n";
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
