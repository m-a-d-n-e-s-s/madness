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

namespace madness{

//  enum functype_6d {pure_, decomposed_, op_decomposed_};
  enum optype {g12_,f12_};
  enum calctype {MP2_, CC2_, LRCCS_, LRCC2_, CISpD_,ADC2_, TDHF_ , experimental_};
  enum functype {HOLE,PARTICLE,MIXED,RESPONSE,UNDEFINED};
//  enum pairtype {GROUND_STATE,EXCITED_STATE};
//  enum potentialtype_s {pot_F3D_,pot_s3a_,pot_s3b_, pot_s3c_,pot_s5a_,pot_s5b_,pot_s5c_,pot_s6_, pot_s2b_, pot_s2c_, pot_s4a_, pot_s4b_, pot_s4c_, pot_S2b_u_, pot_S2c_u_, pot_S4a_u_, pot_S4b_u_, pot_S4c_u_,pot_S2b_r_, pot_S2c_r_, pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_ccs_,pot_cis_,pot_singles_};
//  enum potentialtype_d {pot_F6D_, pot_cc2_coulomb_,pot_cc2_residue_};
//  // The pair function is:  \tau = u + Qf(|titj>), FULL means that \tau is calculated in 6D form, DECOMPOSED means that u is used in 6D and the rest is tried to solve in 3D whenever possible
//  enum pair_function_form{DECOMPOSED, FULL};

  std::string
  assign_name(const optype& input);

  calctype
  assign_calctype(const std::string name);
  std::string
  assign_name(const calctype& inp);

  std::string
  assign_name(const functype& inp);

  // Little structure for formated and controllable output
  struct messenger{
    messenger(World &world) : world(world), output_prec(10), scientific(true), debug(false){}
    World & world;
    size_t output_prec;
    bool scientific;
    bool debug;
    void debug_output(const std::string &msg)const{
      if(debug) output(msg);
    }
    void output(const std::string &msg)const{
      if(scientific)std::cout << std::scientific;
      else std::cout << std::fixed;
      std::cout << std::setprecision(output_prec);
      if(world.rank()==0) std::cout << msg << std::endl;
    }
    void output_section(const std::string &msg)const{
      std::cout <<"\n"<<std::setw(msg.size()+10) << std::setfill('*') << "\n";
      std::setfill(' ');
      output(msg);
      std::cout << std::setw(msg.size()+10) << std::setfill('*') << "\n\n";
      std::setfill(' ');
    }
    void output_subsection(const std::string &msg)const{
      std::cout <<"\n"<< std::setw(msg.size()+5) << std::setfill('-') <<"\n";
      std::setfill(' ');
      output(msg);
      std::cout << std::setw(msg.size()+5) << std::setfill('-') <<"\n";
      std::setfill(' ');
    }
    void warning(const std::string&msg)const{
      std::string tmp = "!!!!!WARNING:" + msg + "!!!!!!";
      output(tmp);
      warnings.push_back(msg);
    }
    mutable std::vector<std::string> warnings;
  };

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

    // print timings and information about the calculated function/vecfunction
    void info(const real_function_3d &f){
      vecfuncT dummy(1,f);
      info(dummy);
    }
    void info(const vecfuncT& f){
      const double size = get_size(world,f);
      const double norm = norm2(world,f);
      update_time();
	if(world.rank()==0){
	  std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
	  << "Timer: " << time_wall << " (Wall), " << time_cpu << " (CPU), " << "||f||=" << norm << ", " << size << " (GB), "  << "(" +operation+")" <<  "\n";
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

    void print(){
      print(current_time());
    }

    void print()const{
      print(std::make_pair(time_wall,time_cpu));
    }

    void print(const std::pair<double,double> &times)const{
      if(world.rank()==0){
      	  std::cout << std::setfill(' ') << std::scientific << std::setprecision(2)
      	  << "Timer: " << times.first << " (Wall), " << times.second << " (CPU)" << ", (" +operation+")" <<  "\n";
      	}
    }


  };

  struct CC_Parameters{
    // default constructor
    //	CC_Parameters():
    //	{return CC_Parameters("default constructor")}

    const double uninitialized = 123.456;

    // copy constructor
    CC_Parameters(const CC_Parameters& other);

    // read parameters from input
    /// ctor reading out the input file
    CC_Parameters(const std::string& input,const double& low);


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
    size_t iter_max;
    size_t iter_max_3D;
    size_t iter_max_6D;
    // restart
    bool restart;
    bool no_compute;
    std::pair<std::size_t,std::size_t> only_pair;
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
    // make additional plots
    bool plot;
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

    /// The number of orbitals which are not zero in the guess (default is one, so the guess is just homo)
    std::size_t tda_guess_orbitals;

    /// Guess mode for TDA:
    /// "numerical" use the std numerical occupied orbitals to create the guess for the excited states
    /// "projected" use the projected occupied orbitals (projected to guess gauss basis) and avoid noise for high guess polynomials
    std::string tda_guess_mode;

    /// The number of excitation vectors for which the alorithm will solve
    size_t tda_excitations;
    /// The number of guess_excitation vectors for the first iterations
    size_t tda_guess_excitations;
    /// The number of excitation vectors which will be iterated parallel
    size_t tda_iterating_excitations;

    /// the guess which will be applied
    /// see the file guess.h
    /// can be "dipole", "dipole+", "quadrupole", "qualdrupole+" , " big_fock_3", "big_fock_4"
    std::string tda_guess;

    /// the guess factor for the first energy guess which is: omega = - factor*HOMO
    /// the factor has to be between ]0,1[
    double tda_energy_guess_factor;

    /// convergence for the excitation vectors
    double tda_dconv_guess;
    double tda_dconv;
    double tda_dconv_hard;
    /// convergence for the excitation energy
    double tda_econv_guess;
    double tda_econv;
    double tda_econv_hard;

    /// store the potential for orthogonalizations or recalculate it (doubles the time but saves memory)
    bool tda_store_potential;

    /// maximum number of iterations in the final iterations
    size_t tda_iter_max;
    /// maximum number of guess iterations (mostly more than the final ones and always without KAIN)
    size_t tda_iter_guess;
    /// specify if for the guess only the homo orbitals (or Homo untill homo-N controlled over guess_orbitals parameter) are used
    bool tda_homo_guess;
    /// Vector of strings which contains the polynomial excitation operators
    /// For this to be used the tda_guess key has to be "custom"
    /// The strings are given in a format like: "c c1 x x1 y y1 z z1, c c2 x x2 y y2 z z2, ..." which will be interpreted as: c1*x^x1*y^y1*z^z1 + c2*x^x2*y^y2*z^z2 + ....
    std::vector<std::string> tda_exops;

    // print out the parameters
    // the TDA parameters are printed out in the TDA section of the program, so no need here
    void
    information(World& world) const;

    void
    print_tda_parameters(World& world) const;
    void
    sanity_check(World& world) const;

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

  /// Intermediates for the CC_Convolution_Operator Structure
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

    /// scalar multiplication
    CC_function operator*(const double &fac)const{
      real_function_3d fnew = fac*function;
      return CC_function(fnew,i,type);
    }

    /// plotting
    void plot(const std::string &msg="")const{
      plot_plane(function.world(),function,msg+name());
    }

  };

  // structure for CC Vectorfunction
  struct CC_vecfunction{

    CC_vecfunction(): type(UNDEFINED),omega(0.0),excitation(-1), current_error(99.9), delta(0.0){}
    CC_vecfunction(const functype type_): type(type_),omega(0.0),excitation(-1), current_error(99.9),delta(0.0){}
    CC_vecfunction(const vecfuncT &v): type(UNDEFINED),omega(0.0),excitation(-1), current_error(99.9),delta(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],i,type);
	functions.insert(std::make_pair(i,tmp));
      }
    }
    CC_vecfunction(const std::vector<CC_function> &v): type(UNDEFINED),omega(0.0),excitation(-1), current_error(99.9),delta(0.0){
      for(size_t i=0;i<v.size();i++){
	functions.insert(std::make_pair(v[i].i,v[i]));
      }
    }
    CC_vecfunction(const vecfuncT &v,const functype &type): type(type),omega(0.0),excitation(-1), current_error(99.9),delta(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],i,type);
	functions.insert(std::make_pair(i,tmp));
      }
    }
    CC_vecfunction(const vecfuncT &v,const functype &type,const size_t &freeze): type(type),omega(0.0),excitation(-1), current_error(99.9),delta(0.0){
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],freeze+i,type);
	functions.insert(std::make_pair(freeze+i,tmp));
      }
    }
    CC_vecfunction(const std::vector<CC_function> &v,const functype type_): type(type_),omega(0.0),excitation(-1),current_error(99.9),delta(0.0){
      for(auto x:v){
	functions.insert(std::make_pair(x.i,x));
      }
    }
    CC_vecfunction(const CC_vecfunction &other) : functions(other.functions),type(other.type), omega(other.omega),excitation(other.excitation),current_error(other.current_error),delta(other.delta) {}

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
    double omega; /// excitation energy
    int excitation; /// the excitation number
    double current_error;
    double delta; // Last difference in Energy

    std::string name()const{
      if (type==PARTICLE) return "tau";
      else if(type==HOLE) return "phi";
      else if(type==MIXED) return "t";
      else if(type==RESPONSE){
	if(excitation<0) MADNESS_EXCEPTION("EXCITATION VECTOR HAS NO NUMBER ASSIGNED!",1);
	return std::to_string(excitation)+"_"+"x";
      }
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

    /// setter
    void set_functions(const vecfuncT & v, const functype& type, const size_t& freeze){
      functions.clear();
      for(size_t i=0;i<v.size();i++){
	CC_function tmp(v[i],freeze+i,type);
	functions.insert(std::make_pair(freeze+i,tmp));
      }
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

    // scalar multiplication
    CC_vecfunction operator*(const double &fac)const{
      vecfuncT vnew = fac*get_vecfunction();
      const size_t freeze = functions.cbegin()->first;
      return CC_vecfunction(vnew,type,freeze);
    }

    // scaling (inplace)
    void scale(const double& factor){
      for(auto& ktmp:functions){
	ktmp.second.function.scale(factor);
      }
    }

    // plotting
    void plot(const std::string &msg="")const{
      for(auto& ktmp:functions){
	ktmp.second.plot(msg);
      }
    }

    // saving the functions
    void save(const std::string& msg="")const{
    std::string pre_name = "";
    if(msg!="") pre_name = msg+"_";
    for(const auto&tmp:functions) madness::save<double,3>(tmp.second.function,pre_name+tmp.second.name());
    }

    // operator needed for sort operation (sorted by omega values)
    bool operator<=(const CC_vecfunction &b)const{return omega<=b.omega;}
    bool operator< (const CC_vecfunction &b)const{return omega<b.omega;}

  };


  /// Helper Structure that carries out operations on CC_functions
  /// The structure can hold intermediates for g12 and f12 of type : <mo_bra_k|op|type> with type=HOLE,PARTICLE or RESPONSE
  /// some 6D operations are also included
  /// The structure does not know if nuclear correlation facors are used, so the corresponding bra states have to be prepared beforehand
  struct CC_convolution_operator{
    /// @param[in] world
    /// @param[in] optype: the operatortype (can be g12_ or f12_)
    /// @param[in] param: the parameters of the current CC-Calculation (including function and operator thresholds and the exponent for f12)
    CC_convolution_operator(World &world,const optype type, const CC_Parameters param):world(world),operator_type(type),op(init_op(type,param)){}

    /// @param[in] f: a 3D function
    /// @param[out] the convolution op(f), no intermediates are used
    real_function_3d operator()(const real_function_3d &f)const {return ((*op)(f)).truncate();}

    // @param[in] f: a vector of 3D functions
    // @param[out] the convolution of op with each function, no intermeditates are used
    vecfuncT operator()(const vecfuncT &f)const{
      return apply<double,double,3>(world,(*op),f);
    }

    // @param[in] bra: a 3D CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] ket: a 3D CC_function,
    // @param[in] use_im: default is true, if false then no intermediates are used
    // @param[out] the convolution <bra|op|ket> = op(bra*ket), if intermediates were calculated before the operator uses them
    real_function_3d operator()(const CC_function &bra, const CC_function &ket, const bool use_im=true)const;

    // @param[in] u: a 6D-function
    // @param[out] the convolution \int g(r,r') u(r,r') dr' (if particle==2) and g(r,r') u(r',r) dr' (if particle==1)
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    real_function_6d operator()(const real_function_6d &u, const size_t particle)const;

    // @param[in] bra: a 3D-CC_function, if nuclear-correlation factors are used they have to be applied before
    // @param[in] u: a 6D-function
    // @param[in] particle: specifies on which particle of u the operator will act (particle ==1 or particle==2)
    // @param[out] the convolution <bra|g12|u>_particle
    real_function_3d operator()(const CC_function &bra,const real_function_6d &u, const size_t particle)const;

    /// @param[in] bra: a vector of CC_functions, the type has to be HOLE
    /// @param[in] ket: a vector of CC_functions, the type can be HOLE,PARTICLE,RESPONSE
    /// updates intermediates of the type <bra|op|ket>
    void update_elements(const CC_vecfunction &bra, const CC_vecfunction &ket);

    /// @param[out] prints the name of the operator (convenience) which is g12 or f12 or maybe other things like gf in the future
    std::string name()const{return assign_name(operator_type);}

    /// @param[in] the type of which intermediates will be deleted
    /// e.g if(type==HOLE) then all intermediates of type <mo_bra_k|op|HOLE> will be deleted
    void clear_intermediates(const functype &type);
    /// name speaks for itself
    void clear_all_intermediates(){
      clear_intermediates(HOLE);
      clear_intermediates(PARTICLE);
      clear_intermediates(RESPONSE);
    }
    /// prints out information (operatorname, number of stored intermediates ...)
    size_t info()const;

    /// sanity check .. doens not do so much
    void sanity()const{print_intermediate(HOLE);}

    /// @param[in] type: the type of intermediates which will be printed, can be HOLE,PARTICLE or RESPONSE
    void print_intermediate(const functype type)const{
      if(type==HOLE)for(const auto& tmp:imH.allpairs)tmp.second.print_size("<H"+std::to_string(int(tmp.first.first))+"|"+assign_name(operator_type)+"|H"+std::to_string(int(tmp.first.second))+"> intermediate");
      else if(type==PARTICLE)for(const auto& tmp:imP.allpairs)tmp.second.print_size("<H"+std::to_string(int(tmp.first.first))+"|"+assign_name(operator_type)+"|P"+std::to_string(int(tmp.first.second))+"> intermediate");
      else if(type==RESPONSE)for(const auto& tmp:imR.allpairs)tmp.second.print_size("<H"+std::to_string(int(tmp.first.first))+"|"+assign_name(operator_type)+"|R"+std::to_string(int(tmp.first.second))+"> intermediate");
    }

  private:
    /// the world
    World &world;
    /// the operatortype, currently this can be g12_ or f12_
    const optype operator_type;
    /// @param[in] optype: can be f12_ or g12_ depending on which operator shall be intitialzied
    /// @param[in] parameters: parameters (thresholds etc)
    /// initializes the operators
    SeparatedConvolution<double,3>* init_op(const optype &type,const CC_Parameters &parameters)const;
    const std::shared_ptr<real_convolution_3d> op;
    intermediateT imH;
    intermediateT imP;
    intermediateT imR;
    /// @param[in] msg: output message
    /// the function will throw an MADNESS_EXCEPTION
    void error(const std::string &msg)const{
      if(world.rank()==0) std::cout <<"\n\n!!!!ERROR in CC_convolution_operator: " << msg <<"!!!!!\n\n"<< std::endl;
      MADNESS_EXCEPTION(msg.c_str(),1);
    }
  };



}//namespace madness

#endif /* CCSTRUCTURES_H_ */
