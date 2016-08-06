/*
 * CCOperators.h
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */

#ifndef CCOPERATORS_H_
#define CCOPERATORS_H_

// Operators for coupled cluster and CIS

#include <chem/CCStructures.h>
#include <chem/projector.h>
#include <chem/nemo.h>
#include <chem/SCFOperators.h>
//#include <string>o

// to debug
//#include<chem/mp2.h>

namespace madness {

  template<size_t NDIM>
  static double rsquare(const Vector<double,NDIM> &r){
    double result = 0.0;
    for(size_t i=0;i<r.size();i++){
      result += r[i]*r[i];
    }
    return result;
  }
  template<size_t NDIM>
  static double gauss_ND(const Vector<double,NDIM> &r){
    const double r2=rsquare<NDIM>(r);
    const double double_N = (double) NDIM;
    const double c = pow(1.0/(sqrt(2.0*M_PI)),double_N);
    return c*exp(-0.5*r2);
  }
  template<size_t NDIM>
  static double unitfunction(const Vector<double,NDIM>&r){
    return 1.0;
  }

  // functors for gauss function
  static double f_gauss(const coord_3d &r) {
    return exp(-((r[0]) * (r[0]) + (r[1]) * (r[1]) + (r[2]) * (r[2])));
  }
  static double f_r2(const coord_3d &r) {
    return (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  }
  static double f_r(const coord_3d &r){
    return sqrt(f_r2(r));
  }
  static double f_laplace_gauss(const coord_3d&r) {
    return -6.0 * f_gauss(r) + 4.0 * f_r2(r) * f_gauss(r);
  }

  typedef std::vector<Function<double, 3> > vecfuncT;

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


  /// Helper structure for the coupling potential of CC Singles and Doubles
  /// because of the regularization of the CC-Wavefunction (for CC2: |tauij> = |uij> + Qt12*f12*|titj>)
  /// we have 6D-functions in std format |u> : type==pure_
  /// we have 6D-functions in sepparated format: type==decomposed_ (e.g O1*f12*|titj> = |xy> with x=|k> and y=<k|f12|ti>*|tj>)
  /// we have 6D-function like f12|xy> which are not needed to be represented on the 6D MRA-Grid, type==op_decomposed_
  struct CC_function_6d{

  public:
    CC_function_6d(World&world,const real_function_6d &ket):world(world), type(pure_), a(),b(), op(0),u(ket) {}
    CC_function_6d(World&world,const vecfuncT &f1,const vecfuncT &f2):world(world), type(decomposed_), a(f1),b(f2), op(0),u() {}
    CC_function_6d(World&world,const std::pair<vecfuncT,vecfuncT> &f):world(world), type(decomposed_), a(f.first),b(f.second), op(0),u() {}
    CC_function_6d(World&world,const CC_convolution_operator *op_,const CC_function &f1, const CC_function &f2):world(world), type(op_decomposed_), a(),b(), op(op_),x(f1),y(f2),u() {}

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D)
    real_function_3d project_out(const CC_function &f,const size_t particle)const;

    // result is: <x|op12|f>_particle
    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    real_function_3d dirac_convolution(const CC_function &x, const CC_convolution_operator &op, const size_t particle)const;

    /// @param[out] particles are interchanged, if the function was u(1,2) the result is u(2,1)
    CC_function_6d swap_particles()const;

    /// @param[in] the Greens operator
    /// @param[out] the Greens operator is applied to the function: G(u)
    real_function_6d apply_G(const real_convolution_6d &G)const;


  private:
    /// the 3 types of 6D-function that occur in the CC potential which coupled doubles to singles
    enum functype_6d {pure_, decomposed_, op_decomposed_};
    World &world;
    /// the type of the given 6D-function
    const functype_6d type;
    /// if type==decomposed this is the first particle
    const vecfuncT a;
    /// if type==decomposed this is the second particle
    const vecfuncT b;
    /// if type==op_decomposed_ this is the symmetric 6D-operator (g12 or f12) in u=op12|xy>
    const CC_convolution_operator* op;
    /// if type==op_decomposed_ this is the first particle in u=op12|xy>
    const CC_function x;
    /// if type==op_decomposed_ this is the second particle in u=op12|xy>
    const CC_function y;
    /// if type=pure_ this is just the MRA 6D-function
    const real_function_6d u;

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=|ab> so <f|u>_particle = <f|a>*|b> if particle==1
    real_function_3d project_out_decomposed(const real_function_3d &f,const size_t particle)const;

    /// @param[in] f: a 3D-CC_function
    /// @param[in] particle: the particle on which the operation acts
    /// @param[out] <f|u>_particle (projection from 6D to 3D) for the case that u=op|xy> so <f|u>_particle = <f|op|x>*|y> if particle==1
    real_function_3d project_out_op_decomposed(const CC_function &f, const size_t particle)const;

    /// @param[in] x: a 3D-CC_function
    /// @param[in] op: a CC_convoltion_operator which is currently either f12 or g12
    /// @param[in] particle: the particle on which the operation acts (can be 1 or 2)
    /// @param[out] the operator is applied and afterwards a convolution with the delta function makes a 3D-function: <x|op|u>_particle
    /// in this case u=|ab> and the result is <x|op|u>_1 = <x|op|a>*|b> for particle==1
    real_function_3d dirac_convolution_decomposed(const CC_function &x, const CC_convolution_operator &op, const size_t particle)const;

    /// small helper function that gives back (a,b) or (b,a) depending on the value of particle
    const std::pair<vecfuncT,vecfuncT> assign_particles(const size_t particle)const;

    /// swap particle function if type==pure_
    CC_function_6d swap_particles_pure() const;
    /// swap particle function if type==decomposed_
    CC_function_6d swap_particles_decomposed()const;
    /// swap particle function if type==op_decomposed_ (all ops are assumed to be symmetric)
    CC_function_6d swap_particles_op_decomposed()const;
  };

  /// Coupled Cluster Operators (all closed shell)
  class CC_Operators {
  public:
    /// Constructor
    CC_Operators(World& world, const Nemo &nemo,
		 //const CorrelationFactor &correlationfactor,
		 const CC_Parameters &param) :
		   world(world),
		   nemo(nemo),
		   corrfac(world,param.gamma(),1.e-7,nemo.get_calc()->molecule),
		   parameters(param),
		   mo_bra_(make_mo_bra(nemo)),
		   mo_ket_(make_mo_ket(nemo)),
		   orbital_energies(init_orbital_energies(nemo)),
		   g12(world,g12_,param),
		   f12(world,f12_,param),
		   //intermediates_(world, mo_bra_, mo_ket_, nemo, param),
		   projector_Q12(world),
		   projector_Q(world,mo_bra_.get_vecfunction(),mo_ket_.get_vecfunction()){
      // make operators

      // make the active mo vector (ket nemos, bra is not needed for that)
      MADNESS_ASSERT(mo_ket_.size() == mo_bra_.size());
      // initialize the Q12 projector
      projector_Q12.set_spaces(mo_bra_.get_vecfunction(), mo_ket_.get_vecfunction(),
			       mo_bra_.get_vecfunction(), mo_ket_.get_vecfunction());
      // make exchange intermediate for ground state
      g12.update_elements(mo_bra_,mo_ket_);
      g12.sanity();
      f12.update_elements(mo_bra_,mo_ket_);
      f12.sanity();
      performance_S.current_iteration = 0;
      performance_D.current_iteration = 0;

    }



    // collect all the data for every function and every iteration
    mutable CC_performance performance_S;
    mutable CC_performance performance_D;
    // collect all the warnings that are put out over a calculation
    mutable std::vector<std::string> warnings;

    // interface to mul_sparse which sets the same tolerance for all operations and makes shure the functions are reconstructed;
    real_function_3d msparse(const real_function_3d &l,const real_function_3d &r)const{
      l.reconstruct();
      r.reconstruct();
      return mul_sparse(l,r,parameters.thresh_3D);
    }


    void save_functions(const CC_vecfunction &f)const{
      CC_Timer time(world,"Saving CC_vecfunction");
      for(const auto&tmp:f.functions) save_function<double,3>(tmp.second.function,tmp.second.name());
      time.info();
    }

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const Function<T, NDIM>& f,const std::string name) const {
      if(world.rank() == 0) print("saving function: ",name);
      f.print_size(name);
      archive::ParallelOutputArchive ar(world,name.c_str(),1);
      ar & f;
    }

    template<typename T, size_t NDIM>
    bool load_function(Function<T, NDIM>& f, const std::string name) const {
      bool exists = archive::ParallelInputArchive::exists(world,name.c_str());
      if(exists){
	if (world.rank() == 0) print("loading function", name);
	archive::ParallelInputArchive ar(world, name.c_str());
	ar & f;
	f.print_size(name);
	return true;
      }else return false;
    }

    void plot(const real_function_3d &f, const std::string &msg) const {
      CC_Timer plot_time(world, "plotting " + msg);
      plot_plane(world, f, msg);
      plot_time.info();
    }

    void error(const std::string &msg) const {
      std::cout << "\n\n\nERROR IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
      MADNESS_EXCEPTION(
	  "\n\n!!!!ERROR IN CC_OPERATORS!!!!\n\n\n\n\n\n\n\n\n\n\n\n", 1);
    }
    void warning(const std::string &msg, CC_data &data) const {
      warning(msg);
      data.warnings.push_back(msg);
    }
    void warning(const std::string &msg) const {
      if(world.rank()==0){
	std::cout << "\n\n" << std::endl;
	std::cout << "|| || || || || || || || || || || || || " << std::endl;
	std::cout << "\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/" << std::endl;
	std::cout << "=> WARNING IN CC_OPERATORS" << std::endl;
	std::cout << "=> " <<msg << std::endl;
	std::cout << "/\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ " << std::endl;
	std::cout << "|| || || || || || || || || || || || || " << std::endl;
	std::cout << "\n\n" << std::endl;
      }
      warnings.push_back(msg);
      std::ofstream file;
      file.open ("warnings");
      file << msg << "\n";
      file.close();
    }

    void output_section(const std::string&msg) const {
      if (world.rank() == 0) {
	std::cout << "\n\n--------------\n";
	std::cout << msg << std::endl;
	std::cout << "\n";
      }
    }
    void output(const std::string &msg) const {
      if (world.rank() == 0) {
	std::cout << msg << std::endl;
      }
    }

    void update_intermediates(const CC_vecfunction &singles) {
      CC_Timer update(world, "Update Intermediates");
      g12.update_elements(mo_bra_,singles);
      f12.update_elements(mo_bra_,singles);
      update.info();
    }

    void update_response_intermediates(const CC_vecfunction &singles) {
      CC_Timer update(world, "Update Intermediates");
      g12.update_elements(mo_bra_,singles);
      f12.update_elements(mo_bra_,singles);
      update.info();
    }

    const vecfuncT get_active_mo_ket()const{
      vecfuncT result;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++) result.push_back(mo_ket_(i).function);
      return result;
    }

    const vecfuncT get_active_mo_bra()const{
      vecfuncT result;
      for(size_t i=parameters.freeze;i<mo_bra_.size();i++) result.push_back(mo_bra_(i).function);
      return result;
    }

    const CC_function mo_ket(const size_t &i) const {
      return mo_ket_(i);
    }
    const CC_vecfunction mo_ket() const {
      return mo_ket_;
    }
    const CC_function mo_bra(const size_t &i) const {
      return mo_bra_(i);
    }
    const CC_vecfunction mo_bra() const {
      return mo_bra_;
    }

    /// makes the t intermediate which is defined as: \f$ |t_i> = |\tau_i> + |i> \f$
    /// for response functions the t intermediate will be just the hole states (CIS(D) -> GS singles are zero)
    CC_function make_t_intermediate(const CC_function &tau) const {
      if(tau.type==PARTICLE){
	CC_function t(mo_ket_(tau.i).function + tau.function, tau.i, MIXED);
	return t;
      }else if(tau.type==RESPONSE){
	return mo_ket_(tau.i);
      }else if(tau.type==HOLE){
	return mo_ket_(tau.i);
      }else{
	error("Wrong type for t_intermediate: " + assign_name(tau.type));
	return CC_function();
      }
    }
    CC_vecfunction make_t_intermediate(const CC_vecfunction &tau) const {
      if(tau.type==PARTICLE){
	CC_vecfunction result(MIXED);
	for (auto x : tau.functions) {
	  CC_function tmpi = make_t_intermediate(x.second);
	  result.insert(tmpi.i, tmpi);
	}
	return result;
      } else if(tau.type==RESPONSE or tau.type==HOLE){
	CC_vecfunction result(HOLE);
	for (auto x : tau.functions) {
	  result.insert(x.second.i, mo_ket_(x.second.i));
	}
	return result;
      }else{
	error("Wrong type for t_intermediate: " + assign_name(tau.type));
	return CC_vecfunction();
      }
    }
    // returns frozen orbitals + t intermedaite
    CC_vecfunction make_t_intermediate_full(const CC_vecfunction &tau)const{
      if(tau.type==HOLE) return mo_ket_;
      CC_vecfunction result(MIXED);
      for(size_t i=0;i<mo_ket_.size();i++){
	if(i<parameters.freeze){
	  result.insert(i,mo_ket_(i));
	}
	else{
	  result.insert(i,make_t_intermediate(tau(i)));
	}
      }
      return result;
    }

    double make_norm(const CC_function &f)const{return make_norm(f.function);}
    // makes vectorfunction norm = sqrt( \sum_i <fi|fi> )
    double make_norm(const CC_vecfunction &f)const{
      double norm2 = 0.0;
      for(const auto& ktmp:f.functions){
	double tmp = make_norm(ktmp.second);
	norm2+= tmp*tmp;
      }
      return sqrt(norm2);
    }
    double make_norm(const real_function_3d &f)const{
      const real_function_3d bra = f*nemo.nuclear_correlation -> square();
      const double norm2 = bra.inner(f);
      return sqrt(norm2);
    }

    void test_singles_potential();


    vecfuncT get_LRCCS_potential(const CC_vecfunction &x){
      if(x.type!=RESPONSE) error("get_LRCCS_potential: Wrong type of input singles");
      Pairs<CC_Pair> empty_doubles;
      CC_vecfunction empty_singles(PARTICLE);
      const vecfuncT fock_residue = response_potential_singles(empty_singles,empty_doubles,x,empty_doubles,pot_F3D_);
      vecfuncT potential =          response_potential_singles(empty_singles,empty_doubles,x,empty_doubles,pot_cis_);
      potential = apply_Q(potential,"CCS-Response-Singles-Potential");
      truncate(world,potential);
      current_singles_potential_response = copy(world,potential);
      return add(world,fock_residue,potential);
    }

    vecfuncT get_CC2_singles_response_potential(const CC_vecfunction &gs_singles,
						const Pairs<CC_Pair> &gs_doubles,const CC_vecfunction &response_singles,
						const Pairs<CC_Pair> &response_doubles) {
      if(gs_singles.type != PARTICLE) error("cc2_singles_response_potential: gs_singles have wrong type");
      if(response_singles.type != RESPONSE) error("cc2_singles_response_potential: response_singles have wrong type");
      if(gs_doubles(parameters.freeze,parameters.freeze).type != GROUND_STATE) error("cc2_singles_response_potential: gs_doubles have wrong type");
      if(response_doubles(parameters.freeze,parameters.freeze).type != EXCITED_STATE) error("cc2_singles_response_potential: response_doubles have wrong type");

      const vecfuncT fock_residue = response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_F3D_);
      vecfuncT potential =          response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_ccs_);
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S2b_u_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S2c_u_));

      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4a_u_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4b_u_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4c_u_));

      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S2b_r_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S2c_r_));

      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4a_r_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4b_r_));
      potential=add(world,potential,response_potential_singles(gs_singles,gs_doubles,response_singles,response_doubles,pot_S4c_r_));

      potential = apply_Q(potential,"CC2-Response-Singles-Potential");
      truncate(world,potential);

      current_singles_potential_response = copy(world,potential);
      vecfuncT result = add(world,fock_residue,potential);
      performance_S.current_iteration++;
      return result;
    }

    vecfuncT get_CC2_singles_potential(const CC_vecfunction &singles,
				       const Pairs<CC_Pair> &doubles) {

      const double norm = make_norm(singles);

      vecfuncT fock_residue =    potential_singles(doubles, singles,pot_F3D_);
      vecfuncT result =          potential_singles(doubles, singles, pot_ccs_);
      result = add(world, result,potential_singles(doubles, singles, pot_S2b_u_));
      result = add(world, result,potential_singles(doubles, singles, pot_S2c_u_));

      if(norm > parameters.thresh_3D){
	result = add(world, result,potential_singles(doubles, singles, pot_S4a_u_));
	result = add(world, result,potential_singles(doubles, singles, pot_S4b_u_));
	result = add(world, result,potential_singles(doubles, singles, pot_S4c_u_));

	result = add(world, result,potential_singles(doubles, singles, pot_S2b_r_));
	result = add(world, result,potential_singles(doubles, singles, pot_S2c_r_));
	result = add(world, result,potential_singles(doubles, singles, pot_S4a_r_));
	result = add(world, result,potential_singles(doubles, singles, pot_S4b_r_));
	result = add(world, result,potential_singles(doubles, singles, pot_S4c_r_));
      }else output("Norm of Singles Vector is below threshold ||singles||=" + std::to_string(norm));

      // the fock residue does not get projected, but all the rest
      result=apply_Q(result,"CC2-Singles-Potential");
      truncate(world, result);
      // need to store this for application of Fock oerator on singles ( F|taui> = -Singles_potential[i] + \epsilon_i|taui>)
      current_singles_potential_gs = copy(world,result);
      result = add(world, result, fock_residue);
      performance_S.current_iteration++;
      return result;
    }


    /// returns: $\f - G(Q*kgj*moi,xk) - G(xk,kgi*moj)   $\f
    real_function_6d make_cispd_coulomb_parts(const CC_function & moi, const CC_function &moj, const CC_vecfunction &x)const{
      if(FunctionDefaults<6>::get_thresh()>parameters.tight_thresh_6D) warning("Thresh has not been tightened for coulomb parts of CIS(D)");
      const double omega = x.omega;
      if(world.rank()==0) std::cout << "current omega for BSH is " << omega << "\n";
      const double bsh_eps = get_epsilon(moi.i,moj.i)+omega;
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 *bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive() = true;

      real_function_6d result = real_factory_6d(world);
      for(const auto& ktmp:x.functions){
	CC_function mokbra = mo_bra_(ktmp.first);
	CC_function xk = ktmp.second;

	const real_function_3d p1 = apply_Q(g12(mokbra,moj)*moi.function);
	const real_function_3d p2 = apply_Q(g12(mokbra,moi)*moj.function);

	const std::string p1_name = "Qg(" + mokbra.name() +"*"+ moj.name() + ") *"+ moi.name();
	const std::string p2_name = "Qg(" + mokbra.name() +"*"+ moj.name() + ") *"+ moi.name();
	p1.print_size("p1="+p1_name);
	p2.print_size("p2="+p2_name);
	const real_function_3d xk1 = copy(xk.function); // G is destructive
	const real_function_3d xk2 = copy(xk.function);
	output("now doing G|p1,"+xk.name()+">");
	result -= -2.0*(G(p1,xk2));
	output("now doing G|"+xk.name()+",p2>");
	result -= -2.0*(G(xk1,p2));

      }

      return result;
    }

    real_function_6d make_cc2_coulomb_parts(const CC_function &taui, const CC_function &tauj, const CC_vecfunction &singles, const double omega=0.0) const;



    real_function_6d make_nonorthogonal_mp2_coulomb_parts(const size_t i,const size_t j)const{
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(i,j)),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;
      const CC_function ti = mo_ket_(i);
      const CC_function tj = mo_ket_(j);
      // first do the O1 and O2 parts which are
      // Otau1(g|titj) = |tauk><k|(1)g|titj> = kgti(2)|tauktj>
      // same for Otau2 = kgtj(1)|titauk>
      real_function_6d G_O1tau_part=real_factory_6d(world);
      real_function_6d G_O2tau_part=real_factory_6d(world);
      for(const auto& ktmp : mo_ket_.functions){
	const size_t k=ktmp.first;
	const CC_function &tauk=ktmp.second;

	real_function_3d kgti_tj=g12(mo_bra_(k),ti) * tj.function;
	real_function_3d kgtj_ti=g12(mo_bra_(k),tj) * ti.function;

	real_function_3d l1=real_factory_3d(world);
	real_function_3d l2=real_factory_3d(world);
	for(const auto& ltmp : mo_ket_.functions){
	  const CC_function& mo_bra_l=mo_bra_(ltmp.first);
	  const CC_function& taul=ltmp.second;
	  l1+=mo_bra_l.function.inner(kgtj_ti) * taul.function;
	  l2+=mo_bra_l.function.inner(kgti_tj) * taul.function;
	}

	real_function_3d part1=-1.0 * kgtj_ti + 0.5 * l1;
	real_function_3d part2=-1.0 * kgti_tj + 0.5 * l2;

	G_O1tau_part+=-2.0 * G(copy(tauk.function),part2);
	G_O2tau_part+=-2.0 * G(part1,copy(tauk.function));
      }

      return G_O1tau_part + G_O2tau_part;
    }
    // computes: G(f(F-eij)|titj> + Ue|titj> - [K,f]|titj>) and uses G-operator screening
    real_function_6d make_cc2_residue_sepparated(const CC_function &taui, const CC_function &tauj, const double omega=0.0)const;



    // return <kl|A|xy>
    // the given operator a should act as: A(x,y,z) = <x|A|yz>_1
    Tensor<double> make_matrix_elements(const CC_function &x, const CC_function &y,real_function_3d (*A)(const CC_function&,const CC_function&, const CC_function&))const{
      Tensor<double> result(mo_bra_.size(),mo_bra_.size());
      for(const auto & ktmp:mo_bra_.functions){
	const CC_function & k=ktmp.second;
	const real_function_3d kfxy = A(k,x,y);
	for(const auto & ltmp:mo_bra_.functions){
	  const CC_function & l=ltmp.second;
	  result(k.i,l.i) = l.inner(kfxy);
	}
      }
      return result;
    }

    /// returns the non constant part of the MP2 potential which is
    /// \f$ (2J-K+Un)|uij> \f$
    real_function_6d get_MP2_potential_residue(const CC_Pair &pair) const {
      CC_Timer timer(world, "(2J-K(R)+Un)|uij>");
      CC_data data("mp2_residue");
      real_function_6d result = fock_residue_6d(pair);

      // make sanity test
      // (T-eps_ij)|uij> = -result
      // do <ij|T|uij> and <ij|result>
      // and <uij|T|uij> = - <uij|result> + eps_ij <uij|uij>
      if(parameters.debug){
	CC_Timer time(world,"MP2-Potential-Sanity-Check");
	output("\n MP2-Potential Sanity Check");
	//	const CC_function moi = mo_bra_(u.i);
	//	const CC_function moj = mo_bra_(u.j);
	std::vector<real_function_6d> grad_u;
	for(size_t axis=0;axis<6;axis++){
	  real_derivative_6d D = free_space_derivative<double,6>(world, axis);
	  const real_function_6d Du = D(pair.function);
	  Du.print_size("d_"+std::to_string(axis)+"("+pair.name()+")");
	  grad_u.push_back(Du);
	}
	const double ur =pair.function.inner(result);
	const double eps_uij = get_epsilon(pair.i,pair.j)*pair.function.inner(pair.function);
	double uTu =0.0;
	for(const auto& k:grad_u){
	  uTu += k.inner(k);
	}
	uTu = 0.5*uTu;
	const double diff = uTu + ur - eps_uij;
	if(world.rank()==0){
	  std::cout << "Current error is " << pair.current_error << "\n";
	  std::cout << "<"<<pair.name()<<"|T|"<<pair.name() << "=" << uTu << "\n";
	  std::cout << "<"<<pair.name()<<"|V_MP2>=" << ur << "\n";
	  std::cout << "<"<<pair.name()<<"|" << pair.name() << "*epsij=" << eps_uij << "\n";
	  std::cout << "0 = " << diff << "\n";
	}
	if(fabs(diff)>fabs(pair.current_error)) warning("MP2 Potential Inaccurate");
	else output("MP2 Potential seems to be sane");
	time.info();
      }

      data.result_size = get_size(result);
      data.result_norm = result.norm2();
      data.time = timer.current_time();
      performance_D.insert(data.name, data);
      timer.info();
      return result;
    }


    double compute_cis_expectation_value(const CC_vecfunction &x, const vecfuncT &V = vecfuncT()){

      const vecfuncT xket = x.get_vecfunction();
      const vecfuncT xbra = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());
      const double norm = sqrt(inner(world,xbra,xket).sum());
      Kinetic<double,3> T(world);
      double kinetic = 0.0;
      for(size_t k=0;k<xket.size();k++) kinetic += T(xbra[k],xket[k]);

      vecfuncT pot;
      if(V.empty()) pot = get_LRCCS_potential(x);
      else pot = V;

      double eps = 0.0;
      for(size_t k=0;k<xket.size();k++) eps -= get_orbital_energies()[k+parameters.freeze]*xbra[k].inner(xket[k]);

      double potential = inner(world,xbra,V).sum();

      const double result = 1.0/norm*(potential + kinetic + eps);
      if(world.rank()==0){
	std::cout << "CCS Expectation Value:\n--------\n";
	std::cout << "Kinetic-Energy  =" << std::fixed << std::setprecision(8) <<  kinetic << "\n";
	std::cout << "Potential-Energy=" << std::fixed << std::setprecision(8) <<  potential << "\n";
	std::cout << "ei*<xi|xi>      =" << std::fixed << std::setprecision(8) <<  eps << "\n";
	std::cout << "||x||           =" << std::fixed << std::setprecision(8) <<  norm << "\n";
	std::cout << "Expectationvalue=" << std::fixed << std::setprecision(8)<< result << "\n--------\n";
      }

      return result;

    }

    double compute_cispd_energy_correction(const CC_vecfunction &x, const Pairs<CC_Pair> &u, const Pairs<CC_Pair> &chi )const{
//      // make empty singles for the ground state
//      const CC_vecfunction tau(zero_functions<double,3>(world,x.size()),PARTICLE,parameters.freeze);
//
//      // now calculate the singles-response potential without the CIS potential (so just the S2 and S4 terms since tau is zero anyways)
//      vecfuncT V = response_potential_singles(tau,u,x,chi,pot_S2b_u_);
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S2c_u_));
//
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4a_u_));
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4b_u_));
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4c_u_));
//
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S2b_r_));
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S2c_r_));
//
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4a_r_));
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4b_r_));
//      V = add(world,V,response_potential_singles(tau,u,x,chi,pot_S4c_r_));
//
//      // apply Q-projector
//      V = Q(V,mo_ket_);
//      truncate(world,V);
//      const double Vnorm = norm2(world,V);
//
//      // make inner product with cis_states
//      vecfuncT xbra = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());
//      Tensor<double> xV = inner(world,xbra,V);
//      const double result = xV.sum();
      if(world.rank()==0){
	std::cout << "Current S2b_reg_part_gs\n";
      }
      print_size(world,current_s2b_reg_part_gs,"current_s2b_reg_part_gs");
      current_s2b_reg_part_gs.clear();
      double result = 0.0;
      result += compute_cispd_energy_V(x,u,chi,pot_S2b_u_);
      result += compute_cispd_energy_V(x,u,chi,pot_S2c_u_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4a_u_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4b_u_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4c_u_);
      result += compute_cispd_energy_V(x,u,chi,pot_S2b_r_);
      result += compute_cispd_energy_V(x,u,chi,pot_S2c_r_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4a_r_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4b_r_);
      result += compute_cispd_energy_V(x,u,chi,pot_S4c_r_);

      if(world.rank()==0){
	std::cout << "CIS(D) Energy Correction to CCS/CIS Excitation Energy: " << x.omega << "\n";
	//std::cout << "||V||     =" << Vnorm  << "\n";
	std::cout << "CIS       =" << std::fixed << std::setprecision(parameters.output_prec+3)<< x.omega << "\n";
	std::cout << "CIS(D)    =" << std::fixed << std::setprecision(parameters.output_prec+3)<< x.omega+result << "\n";
	std::cout << "correction=" << std::fixed << std::setprecision(parameters.output_prec+3)<< result << "\n";
      }
      print_size(world,current_s2b_reg_part_gs,"current_s2b_reg_part_gs");
      return result;
    }
    double compute_cispd_energy_V(const CC_vecfunction &x, const Pairs<CC_Pair> u, const Pairs<CC_Pair> &chi, const potentialtype_s &name)const{
      const CC_vecfunction tau(zero_functions<double,3>(world,x.size()),PARTICLE,parameters.freeze);
      const vecfuncT V = response_potential_singles(tau,u,x,chi,name);
      const vecfuncT QV = Q(V,mo_ket_);
      const double Vnorm = norm2(world,V);
      const double QVnorm = norm2(world,QV);
      vecfuncT xbra = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());
      const double result = inner(world,xbra,QV).sum();

      if(world.rank()==0){
	std::cout << "CIS(D) Contribution from " + assign_name(name) + ":\n";
	std::cout << "||V|| =" <<  Vnorm << "\n";
	std::cout << "||QV||=" << QVnorm << "\n";
	std::cout << " delta-omega = " << result << "\n";
      }
      return result;
    }



    // right now this s4cuis all copied from mp2.cc
    double compute_mp2_pair_energy(CC_Pair &pair) const;

    /// returns <mo_bra|f> |g>
    vecfuncT P(const vecfuncT &f, const CC_vecfunction &g)const{
      vecfuncT active_mo_bra;
      for(const auto& gtmp:g.functions) active_mo_bra.push_back(mo_bra_(gtmp.first).function);
      Projector<double,3> P(active_mo_bra,g.get_vecfunction());
      return P(f);
    }
    /// returns (1-factor*P)f
    vecfuncT Q(const vecfuncT &f, const CC_vecfunction &g, const double factor=1.0)const{
      if(world.rank()==0) std::cout << "Applying Q-Projector with " << g.name() << " functions and factor " << factor << "\n";
      vecfuncT pf=P(f,g);
      if(factor!=1.0) scale(world,pf,factor);
      vecfuncT result = sub(world,f,pf);
      return result;
    }
    real_function_3d Q(const real_function_3d &f, const CC_vecfunction &g)const{
      vecfuncT tmp; tmp.push_back(f); return Q(tmp,g).front();
    }
    real_function_3d P(const real_function_3d &f, const CC_vecfunction &g)const{
      vecfuncT tmp; tmp.push_back(f); return P(tmp,g).front();
    }

    /// CC2 singles potential parts of the ground state
    /// only for CC2 ground state right now
    /// Genereal function which evaluates a CC_singles potential
    vecfuncT potential_singles(const Pairs<CC_Pair>& u,
			       const CC_vecfunction & singles,
			       const potentialtype_s &name) const {
      //output_section("Now doing Singles Potential " + assign_name(name));
      if (singles.functions.size() != mo_ket_.size() - parameters.freeze)
	warning("Somethings wrong: Size of singles unequal to size of orbitals minus freeze parameter");
      const std::string full_name = assign_name(CC2_)+":"+assign_name(name);
      CC_Timer timer(world, full_name);
      CC_data data(full_name);

      vecfuncT result;

      switch (name) {
	case pot_F3D_:
	  result = fock_residue_closed_shell(singles);
	  break;
	case pot_ccs_:{
	  result = ccs_potential_new(singles);
	  std::cout <<"Make new ccs potential\n";
	  vecfuncT old = ccs_potential(singles);
	  std::cout <<"Make old ccs potential\n";
	  vecfuncT diff = sub(world,result,old);
	  double norm2 = inner(world,diff,diff).sum();
	  std::cout <<"Difference between new and old ccs potential: " << sqrt(norm2) << std::endl;
	  if(sqrt(norm2)>parameters.thresh_3D) warning("Warning: New CCS Potential not good");
	}
	break;
	case pot_cis_: error("No Ground State Singles Potential for CIS");
	break;
	case pot_S2b_u_:{
	  vecfuncT s2b = current_s2b_u_part_gs;
	  if(s2b.empty()){
	    s2b = S2b_u_part(make_pairs(u),singles);
	    vecfuncT tmp = copy(world,s2b);
	    truncate(world,tmp);
	    current_s2b_u_part_gs = tmp;
	  }else output("Found u-part of S2b-potential");
	  result = s2b;
	}
	break;
	case pot_S2c_u_:{
	  vecfuncT s2c = current_s2c_u_part_gs;
	  if(s2c.empty()){
	    s2c=S2c_u_part(make_pairs(u), singles);
	    vecfuncT tmp = copy(world,s2c);
	    truncate(world,s2c);
	    current_s2c_u_part_gs = tmp;
	  }else output("Found u-part of S2c-potential");
	  result = s2c;
	}
	break;
	case pot_S4a_u_:{
	  const vecfuncT s2b = current_s2b_u_part_gs;
	  result = S4a_from_S2b(s2b,singles);
	  //result = S4a_u_part(u, singles);
	}
	break;
	case pot_S4b_u_:
	  result = S4b_u_part(make_pairs(u), singles);
	  break;
	case pot_S4c_u_:
	  result = S4c_u_part(make_pairs(u), singles);
	  break;
	case pot_S2b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  CC_vecfunction tfull = make_t_intermediate_full(singles);
	  //result = S2b_reg_part(t,t,tfull);
	  vecfuncT f12_part = S2b_gf_part(t,t);
	  vecfuncT tmp1 = S2b_u_part(make_regularization_tails_O1(tfull,tfull,t,t,0.5),singles);
	  vecfuncT tmp2 = S2b_u_part(make_regularization_tails_O2(tfull,tfull,t,t,0.5),singles);
	  vecfuncT projected_part = add(world,tmp1,tmp2);
	  result=sub(world,f12_part,projected_part);
	  // save for s4a_reg_ potential
	  vecfuncT tmp = copy(world,result);
	  truncate(world,tmp);
	  current_s2b_reg_part_gs = tmp;
	  break;
	}
	case pot_S2c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  CC_vecfunction tfull = make_t_intermediate_full(singles);
	  //	  result = S2c_reg_part(t,t,tfull);
	  vecfuncT f12_part = S2c_u_part(make_regularization_tails_f12(t,t),singles);
	  vecfuncT tmp1 = S2c_u_part(make_regularization_tails_O1(tfull,tfull,t,t,0.5),singles);
	  vecfuncT tmp2 = S2c_u_part(make_regularization_tails_O2(tfull,tfull,t,t,0.5),singles);
	  vecfuncT projected_part = add(world,tmp1,tmp2);
	  result=sub(world,f12_part,projected_part);
	  break;
	}
	case pot_S4a_r_:
	{
	  //	  CC_vecfunction t = make_t_intermediate(singles);
	  //	  CC_vecfunction tfull = make_t_intermediate_full(singles);
	  //	  result = S4a_reg_part(t,t,singles,tfull);
	  const vecfuncT s2b = current_s2b_reg_part_gs;
	  result = S4a_from_S2b(s2b,singles);
	  break;
	}
	case pot_S4b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  CC_vecfunction tfull = make_t_intermediate_full(singles);
	  //result = S4b_reg_part(t,t,singles,tfull);
	  vecfuncT tmp1 = S4b_u_part(make_regularization_tails_O1(tfull,tfull,t,t,0.5),singles);
	  vecfuncT tmp2 = S4b_u_part(make_regularization_tails_O2(tfull,tfull,t,t,0.5),singles);
	  vecfuncT projected_part=add(world,tmp1,tmp2);
	  vecfuncT f12_part = S4b_u_part(make_regularization_tails_f12(t,t),singles);
	  result = sub(world,f12_part,projected_part);
	  break;
	}
	case pot_S4c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  CC_vecfunction tfull = make_t_intermediate_full(singles);
	  //result = S4c_reg_part(t,t,singles,tfull);
	  vecfuncT tmp1 = S4c_u_part(make_regularization_tails_O1(tfull,tfull,t,t,0.5),singles);
	  vecfuncT tmp2 = S4c_u_part(make_regularization_tails_O2(tfull,tfull,t,t,0.5),singles);
	  vecfuncT projected_part=add(world,tmp1,tmp2);
	  vecfuncT f12_part = S4c_u_part(make_regularization_tails_f12(t,t),singles);
	  result = sub(world,f12_part,projected_part);
	  break;
	}
      }

      data.result_size = get_size(world, result);
      data.result_norm = norm2(world, result);
      data.time = timer.current_time();
      data.info();
      performance_S.insert(data.name, data);
      return result;
    }

    // same thing as potential_singles for the response potential of CC2
    vecfuncT response_potential_singles(
	const CC_vecfunction & tau,
	const Pairs<CC_Pair>& u,
	const CC_vecfunction & x,
	const Pairs<CC_Pair>& chi,
	const potentialtype_s &name) const {

      const std::string full_name = assign_name(LRCC2_)+":"+assign_name(name);
      CC_Timer timer(world, full_name);
      CC_data data(full_name);
      if(x.type!=RESPONSE) error("Error in response_potential_singles " + assign_name(name)+" response_singles have wrong type");
      if(tau.type!=PARTICLE) error("Error in response_potential_singles " + assign_name(name)+" gs_singles have wrong type");
      vecfuncT result;

      switch (name) {
	case pot_F3D_:
	  result = fock_residue_closed_shell(x);
	  break;
	case pot_ccs_:
	  result = LRCCS_potential(tau,x);
	  break;
	case pot_cis_:
	  result = cis_response_potential(x);
	  break;
	case pot_S2b_u_:{
	  vecfuncT s2b = current_s2b_u_part_response;
	  if(s2b.empty()){
	    s2b = S2b_u_part(make_pairs(chi),tau);
	    vecfuncT tmp = copy(world,s2b);
	    truncate(world,tmp);
	    current_s2b_u_part_response = tmp;
	  }else output("Found u-part of S2b-potential");
	  result = s2b;
	}
	break;
	case pot_S2c_u_:{
	  vecfuncT s2c = current_s2c_u_part_gs;
	  if(s2c.empty()){
	    s2c=S2c_u_part(make_pairs(chi),tau);
	    vecfuncT tmp = copy(world,s2c);
	    truncate(world,s2c);
	    current_s2c_u_part_gs = tmp;
	  }else output("Found u-part of S2c-potential");
	  result = s2c;
	}break;
	case pot_S4a_u_:{
	  vecfuncT s2b_response = current_s2b_u_part_response;
	  if(s2b_response.empty()) s2b_response = S2b_u_part(make_pairs(chi),tau);
	  const vecfuncT part1 = S4a_from_S2b(s2b_response,tau);
	  vecfuncT s2b_gs = current_s2b_u_part_gs;
	  if(s2b_gs.empty()) s2b_gs = S2b_u_part(make_pairs(u),x);
	  const vecfuncT part2 = S4a_from_S2b(s2b_gs,x);
	  print_size(world,part2,assign_name(name)+" ux-part"); // for CIS(D) this is the nonzero part
	  result = add(world,part1,part2);
	}break;
	case pot_S4b_u_:{
	  const vecfuncT part1 = S4b_u_part(make_pairs(u), x);
	  print_size(world,part1,assign_name(name)+" ux-part"); // for CIS(D) this is the nonzero part
	  const vecfuncT part2 = S4b_u_part(make_pairs(chi), tau);
	  result = add(world,part1,part2);
	}break;
	case pot_S4c_u_:{
	  const vecfuncT part1 = S4c_u_part(make_pairs(u), x);
	  print_size(world,part1,assign_name(name)+" ux-part"); // for CIS(D) this is the nonzero part
	  const vecfuncT part2 = S4c_u_part(make_pairs(chi), tau);
	  result = add(world,part1,part2);
	}break;
	case pot_S2b_r_:{
	  CC_vecfunction t = make_t_intermediate(tau);
	  CC_vecfunction tfull = make_t_intermediate_full(tau);

	  const Pairs<CC_function_6d> reg_O1_xt = make_regularization_tails_O1(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> reg_O2_xt = make_regularization_tails_O2(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> reg_O1_tx = make_regularization_tails_O1(tfull,tfull,t,x,0.5);
	  const Pairs<CC_function_6d> reg_O2_tx = make_regularization_tails_O2(tfull,tfull,t,x,0.5);

	  const Pairs<CC_function_6d> reg_O1_tt = make_regularization_tails_O1(x,tfull,t,t,1.0);
	  const Pairs<CC_function_6d> reg_O2_tt = make_regularization_tails_O2(tfull,x,t,t,1.0);

	  vecfuncT f12_part_xt = S2b_gf_part(x,t);
	  vecfuncT f12_part_tx = S2b_gf_part(t,x);
	  vecfuncT f12_part = add(world,f12_part_xt,f12_part_tx);

	  vecfuncT function_response_projected = S2b_u_part(reg_O1_xt,tau);
	  function_response_projected=add(world,function_response_projected,S2b_u_part(reg_O1_tx,tau));
	  function_response_projected=add(world,function_response_projected,S2b_u_part(reg_O2_xt,tau));
	  function_response_projected=add(world,function_response_projected,S2b_u_part(reg_O2_tx,tau));

	  vecfuncT function_response = sub(world,f12_part,function_response_projected);

	  vecfuncT projector_response = add(world,S2b_u_part(reg_O1_tt,tau),S2b_u_part(reg_O2_tt,tau));

	  result = sub(world,function_response,projector_response);
	  // save for s4a_reg_ potential
	  vecfuncT tmp = copy(world,result);
	  truncate(world,tmp);
	  current_s2b_reg_part_response = tmp;
	}break;
	case pot_S2c_r_:{
	  CC_vecfunction t = make_t_intermediate(tau);
	  CC_vecfunction tfull = make_t_intermediate_full(tau);

	  // NO GS REG PAIRS FOR S2C -> ONLY RESPONSE FUNCTION (SINCE NO SINGLES IN THIS POTENTIAL)

	  // response reg_pairs: Qtf|xitj+tixj> - (OxQt+QtOx)f|titj>
	  const Pairs<CC_function_6d> re_reg_f12_xt=make_regularization_tails_f12(x,t);
	  const Pairs<CC_function_6d> re_reg_f12_tx=make_regularization_tails_f12(t,x);
	  const Pairs<CC_function_6d> re_reg_O1_xt=make_regularization_tails_O1(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O1_tx=make_regularization_tails_O1(tfull,tfull,t,x,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_xt=make_regularization_tails_O2(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_tx=make_regularization_tails_O2(tfull,tfull,t,x,0.5);

	  const Pairs<CC_function_6d> re_reg_O1_tt=make_regularization_tails_O1(x,tfull,t,t,1.0);
	  const Pairs<CC_function_6d> re_reg_O2_tt=make_regularization_tails_O2(tfull,x,t,t,1.0);

	  // chi with tau function response
	  result = 		    S2c_u_part(re_reg_f12_xt,tau);
	  result = add(world,result,S2c_u_part(re_reg_f12_tx,tau));
	  result = sub(world,result,S2c_u_part(re_reg_O1_tx,tau));
	  result = sub(world,result,S2c_u_part(re_reg_O1_xt,tau));
	  result = sub(world,result,S2c_u_part(re_reg_O2_tx,tau));
	  result = sub(world,result,S2c_u_part(re_reg_O2_xt,tau));
	  // chi with tau projector response
	  result = sub(world,result,S2c_u_part(re_reg_O1_tt,tau));
	  result = sub(world,result,S2c_u_part(re_reg_O2_tt,tau));
	}break;
	case pot_S4a_r_:{
	  vecfuncT s2b_response = current_s2b_reg_part_response;
	  if(s2b_response.empty()){
	    output("recalculating s2b_reg_part_response");
	    s2b_response = response_potential_singles(tau,u,x,chi,pot_S2b_r_);
	  }
	  const vecfuncT part1 = S4a_from_S2b(s2b_response,tau);
	  vecfuncT s2b_gs = current_s2b_reg_part_gs;
	  print_size(world,current_s2b_reg_part_gs,"current_s2b_reg_part_gs");
	  if(s2b_gs.empty()){
	    output("recalculating s2b_reg_part_gs");
	    s2b_gs = potential_singles(u,tau,pot_S2b_r_);
	  }
	  print_size(world,current_s2b_reg_part_gs,"current_s2b_reg_part_gs");
	  const vecfuncT part2 = S4a_from_S2b(s2b_gs,x);
	  print_size(world,part2,"S4a ux part");
	  result = add(world,part1,part2);
	}break;
	case pot_S4b_r_:{
	  CC_vecfunction t = make_t_intermediate(tau);
	  CC_vecfunction tfull = make_t_intermediate_full(tau);

	  // ground state reg_pairs: Qtf|titj>
	  const Pairs<CC_function_6d> gs_reg_f12=make_regularization_tails_f12(t,t);
	  const Pairs<CC_function_6d> gs_reg_O1=make_regularization_tails_O1(tfull,tfull,t,t,0.5);
	  const Pairs<CC_function_6d> gs_reg_O2=make_regularization_tails_O2(tfull,tfull,t,t,0.5);
	  // response reg_pairs: Qtf|xitj+tixj> - (OxQt+QtOx)f|titj>
	  const Pairs<CC_function_6d> re_reg_f12_xt=make_regularization_tails_f12(x,t);
	  const Pairs<CC_function_6d> re_reg_f12_tx=make_regularization_tails_f12(t,x);
	  const Pairs<CC_function_6d> re_reg_O1_xt=make_regularization_tails_O1(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O1_tx=make_regularization_tails_O1(tfull,tfull,t,x,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_xt=make_regularization_tails_O2(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_tx=make_regularization_tails_O2(tfull,tfull,t,x,0.5);

	  const Pairs<CC_function_6d> re_reg_O1_tt=make_regularization_tails_O1(x,tfull,t,t,1.0);
	  const Pairs<CC_function_6d> re_reg_O2_tt=make_regularization_tails_O2(tfull,x,t,t,1.0);

	  // u with x
	  result = S4b_u_part(gs_reg_f12,x);
	  result = sub(world,result,S4b_u_part(gs_reg_O1,x));
	  result = sub(world,result,S4b_u_part(gs_reg_O2,x));
	  print_size(world,result,assign_name(name)+" ux-part"); // for CIS(D) this is the nonzero part

	  // chi with tau function response
	  result = add(world,result,S4b_u_part(re_reg_f12_xt,tau));
	  result = add(world,result,S4b_u_part(re_reg_f12_tx,tau));
	  result = sub(world,result,S4b_u_part(re_reg_O1_tx,tau));
	  result = sub(world,result,S4b_u_part(re_reg_O1_xt,tau));
	  result = sub(world,result,S4b_u_part(re_reg_O1_tx,tau));
	  result = sub(world,result,S4b_u_part(re_reg_O2_xt,tau));
	  // chi with tau projector response
	  result = sub(world,result,S4b_u_part(re_reg_O1_tt,tau));
	  result = sub(world,result,S4b_u_part(re_reg_O2_tt,tau));
	}break;
	case pot_S4c_r_:{
	  CC_vecfunction t = make_t_intermediate(tau);
	  CC_vecfunction tfull = make_t_intermediate_full(tau);

	  // ground state reg_pairs: Qtf|titj>
	  const Pairs<CC_function_6d> gs_reg_f12=make_regularization_tails_f12(t,t);
	  const Pairs<CC_function_6d> gs_reg_O1=make_regularization_tails_O1(tfull,tfull,t,t,0.5);
	  const Pairs<CC_function_6d> gs_reg_O2=make_regularization_tails_O2(tfull,tfull,t,t,0.5);
	  // response reg_pairs: Qtf|xitj+tixj> - (OxQt+QtOx)f|titj>
	  const Pairs<CC_function_6d> re_reg_f12_xt=make_regularization_tails_f12(x,t);
	  const Pairs<CC_function_6d> re_reg_f12_tx=make_regularization_tails_f12(t,x);
	  const Pairs<CC_function_6d> re_reg_O1_xt=make_regularization_tails_O1(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O1_tx=make_regularization_tails_O1(tfull,tfull,t,x,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_xt=make_regularization_tails_O2(tfull,tfull,x,t,0.5);
	  const Pairs<CC_function_6d> re_reg_O2_tx=make_regularization_tails_O2(tfull,tfull,t,x,0.5);

	  const Pairs<CC_function_6d> re_reg_O1_tt=make_regularization_tails_O1(x,tfull,t,t,1.0);
	  const Pairs<CC_function_6d> re_reg_O2_tt=make_regularization_tails_O2(tfull,x,t,t,1.0);

	  // u with x
	  result = S4c_u_part(gs_reg_f12,x);
	  result = sub(world,result,S4c_u_part(gs_reg_O1,x));
	  result = sub(world,result,S4c_u_part(gs_reg_O2,x));
	  print_size(world,result,assign_name(name)+" ux-part"); // for CIS(D) this is the nonzero part

	  // chi with tau function response
	  result = add(world,result,S4c_u_part(re_reg_f12_xt,tau));
	  result = add(world,result,S4c_u_part(re_reg_f12_tx,tau));
	  result = sub(world,result,S4c_u_part(re_reg_O1_tx,tau));
	  result = sub(world,result,S4c_u_part(re_reg_O1_xt,tau));
	  result = sub(world,result,S4c_u_part(re_reg_O1_tx,tau));
	  result = sub(world,result,S4c_u_part(re_reg_O2_xt,tau));
	  // chi with tau projector response
	  result = sub(world,result,S4c_u_part(re_reg_O1_tt,tau));
	  result = sub(world,result,S4c_u_part(re_reg_O2_tt,tau));
	}break;
      }
      data.result_size = get_size(world, result);
      data.result_norm = norm2(world, result);
      data.time = timer.current_time();
      data.info();
      performance_S.insert(data.name, data);
      return result;

    }


    vecfuncT fock_residue_closed_shell(const CC_vecfunction &tau) const;

    /// brilloin terms of ccs
    vecfuncT S1(const CC_vecfunction &tau) const;
    vecfuncT S5a(const CC_vecfunction &tau) const;

    /// The CCS Potential without Brillouin terms and Fock residue
    vecfuncT ccs_potential(const CC_vecfunction &tau) const;
    vecfuncT ccs_potential_new(const CC_vecfunction &tau)const;
    vecfuncT ccs_unprojected(const CC_vecfunction &ti, const CC_vecfunction &tk)const;
    vecfuncT LRCCS_potential(const CC_vecfunction &singles, const CC_vecfunction &response)const;

    vecfuncT cis_response_potential(const CC_vecfunction &x)const{
      vecfuncT result;
      for(const auto& itmp:x.functions){
	const size_t i=itmp.first;
	real_function_3d resulti = real_factory_3d(world);
	for(const auto& ktmp:x.functions){
	  const size_t k=ktmp.first;
	  resulti+= 2.0*g12(mo_bra_(k),x(k))*mo_ket(i).function - g12(mo_bra_(k),mo_ket_(i))*x(k).function;
	}
	result.push_back(resulti);
      }
      return result;
    }


    // result: \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 )
    // singles are not needed explicitly but to determine if it is response or ground state
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S2b_u_part(const Pairs<CC_function_6d> &doubles,
			const CC_vecfunction &singles) const;

    // the <k|gf|xy> part of S2b where the f|xy> term comes from the unprojected regularization tail
    vecfuncT S2b_gf_part(const CC_vecfunction &x, const CC_vecfunction &y)const{
      vecfuncT result;
      for(const auto &itmp:x.functions){
	const CC_function &ti = x(itmp.first);
	real_function_3d resulti = real_factory_3d(world);

	for(const auto &ktmp:y.functions){
	  const CC_function &tk = y(ktmp.first);
	  const CC_function &k = mo_bra_(ktmp.first);
	  const real_function_3d part1 = apply_gf(k.function*tk.function)*ti.function;
	  const real_function_3d partx = apply_gf(k.function*ti.function)*tk.function;
	  resulti += 2.0*part1 - partx;
	}
	result.push_back(resulti);
      }
      return result;
    }
    // result: -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1)
    // singles are not needed explicitly but to determine if it is response or ground state
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1) \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S2c_u_part(const Pairs<CC_function_6d> &doubles,
			const CC_vecfunction &singles) const;

    vecfuncT
    S2c_u_part_old(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
      vecfuncT result;
      if(singles.type==PARTICLE) result=copy(world,current_s2c_u_part_gs);
      else if(singles.type==RESPONSE) result=copy(world,current_s2c_u_part_response);
      else warning("singles of type " + assign_name(singles.type) +" in S2c_u_part");

      if(not result.empty()){
	output("S2c_u_part already calculated");
      }else{
	for(const auto& itmp : singles.functions){
	  const size_t i=itmp.first;
	  real_function_3d resulti=real_factory_3d(world);
	  for(const auto& ktmp : singles.functions){
	    const size_t k=ktmp.first;
	    const real_function_3d kgi=g12(mo_bra_(k),mo_ket_(i));
	    for(const auto& ltmp : singles.functions){
	      const size_t l=ltmp.first;
	      const real_function_6d ukl=get_pair_function(doubles,k,l);
	      const real_function_3d l_kgi=mo_bra_(l).function * kgi;
	      resulti+=-2.0 * ukl.project_out(l_kgi,1);     // 1 means second particle
	      resulti+=ukl.project_out(l_kgi,0);
	    }
	  }
	  result.push_back(resulti);
	}
	if(singles.type==PARTICLE) current_s2c_u_part_gs=copy(world,result);
	else if(singles.type==RESPONSE) current_s2c_u_part_response = copy(world,result);
      }
      return result;
    }

    /// The Part of the CC2 singles potential which depends on singles and doubles (S4a, S4b, S4c)

    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -|singles_l><l|S2b_i> = -|singles_l>(2.0*<lk|g|uik>-<kl|g|uik> \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4a_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    /// the S4a potential can be calcualted from the S2b potential
    /// result is \f$ s4a_i = - <l|s2b_i>*|tau_l> \f$
    vecfuncT S4a_from_S2b(const vecfuncT &s2b, const CC_vecfunction&singles)const{
      if(s2b.empty())warning("S2b-potential is empy --> S4a will be zero");
      vecfuncT result;
      for(size_t i=0;i<s2b.size();i++){
	real_function_3d resulti = real_factory_3d(world);
	const Tensor<double> ls2bi = inner(world,s2b[i],get_active_mo_bra());
	for(const auto& ltmp:singles.functions){
	  resulti -= ls2bi[ltmp.first-parameters.freeze]*singles(ltmp.first).function;
	}
	result.push_back(resulti);
      }
      return result;
    }

    // result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4b_u_part(const Pairs<CC_function_6d> &doubles,
			const CC_vecfunction &singles) const;

    vecfuncT
    S4b_u_part_old(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
      vecfuncT result;
      const vecfuncT active_mo_bra = get_active_mo_bra();
      for(const auto& itmp : singles.functions){
	const size_t i=itmp.first;
	real_function_3d resulti=real_factory_3d(world);
	for(const auto& ktmp : singles.functions){
	  const size_t k=ktmp.first;
	  const real_function_3d kgi = g12(mo_bra_(k),singles(i));
	  vecfuncT l_kgi = mul_sparse(world,kgi,active_mo_bra,parameters.thresh_3D);
	  truncate(world,l_kgi);
	  for(const auto& ltmp : singles.functions){
	    const size_t l=ltmp.first;
	    const real_function_6d ukl=get_pair_function(doubles,k,l);
	    //const real_function_3d l_kgi=mo_bra_(l).function * kgi;
	    resulti+=-2.0 * ukl.project_out(l_kgi[l-parameters.freeze],1);     // 1 means second particle
	    resulti+=ukl.project_out(l_kgi[l-parameters.freeze],0);
	  }
	}
	resulti.print_size("s4b_"+std::to_string(i));
	result.push_back(resulti);
      }
      return result;
    }
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ ( 4<l|kgtauk|uil>_2 - 2<l|kgtauk|uil>_1 - 2<k|lgtauk|uil>_2 + <k|lgtauk|uil>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4c_u_part(const Pairs<CC_function_6d> &doubles,
			const CC_vecfunction &singles) const;

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[out] \f$ 2<k|gQf|reg1_i,reg2_k>_2 - <k|gQf|reg1_i,reg2_k>_1 \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S2b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2, const CC_vecfunction &projector_t) const;


    /// for response equations the pair function looks like: tauij = uij + Qtf|xitj> + Qtf|tixj> - (OxQt + QtOx)f|titj>
    vecfuncT S2b_response_reg_part(const CC_vecfunction x, const CC_vecfunction &t, const CC_vecfunction& projector_t) const{
      CC_Timer time1(world,"Qtf|xt>");
      const vecfuncT Qtfxt = S2b_reg_part(x,t,projector_t);
      time1.info();
      CC_Timer time2(world,"Qtf|tx>");
      const vecfuncT Qtftx = S2b_reg_part(t,x,projector_t);
      time2.info();
      const vecfuncT function_response = add(world,Qtfxt,Qtftx);

      vecfuncT OxQt_part;
      vecfuncT QtOx_part;
      CC_Timer time3(world,"(OxQt-QtOx)f|tt>");
      for(const auto& itmp:t.functions){
	const size_t i=itmp.first;
	const CC_function ti = t(i);

	real_function_3d OxQt_i = real_factory_3d(world);
	real_function_3d QtOx_i = real_factory_3d(world);

	for(const auto& ktmp:t.functions){
	  const size_t k=ktmp.first;
	  const CC_function tk=t(k);

	  const std::pair<vecfuncT,vecfuncT> Ox1_ftitk = make_O1t_op_xy(x,f12,ti,tk);
	  const std::pair<vecfuncT,vecfuncT> OxQt_ftitk = std::make_pair(Ox1_ftitk.first,Q(Ox1_ftitk.second,projector_t));

	  const std::pair<vecfuncT,vecfuncT> Ox2_ftitk = make_O2t_op_xy(x,f12,ti,tk);
	  const std::pair<vecfuncT,vecfuncT>  QtOx_ftitk = std::make_pair(Q(Ox2_ftitk.first,projector_t),Ox2_ftitk.second);

	  for(std::size_t m=0;m<QtOx_ftitk.first.size();m++){
	    OxQt_i += 2.0*(g12(mo_bra_(k).function*OxQt_ftitk.second[m])*OxQt_ftitk.first[m]) - (g12(mo_bra_(k).function*OxQt_ftitk.first[m])*OxQt_ftitk.second[m]);
	    QtOx_i += 2.0*(g12(mo_bra_(k).function*QtOx_ftitk.second[m])*QtOx_ftitk.first[m]) - (g12(mo_bra_(k).function*QtOx_ftitk.first[m])*QtOx_ftitk.second[m]);
	  }


	}
	OxQt_part.push_back(OxQt_i);
	QtOx_part.push_back(QtOx_i);
      }

      const vecfuncT projector_response = add(world,OxQt_part,QtOx_part);

      const vecfuncT result = sub(world,function_response,projector_response);
      time3.info();
      return result;

    }

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[out] \f$ 2<l*kgi|Qf|reg1_k,reg2_l>_2 - <l*kgi|gQf|reg1_k,reg2_l>_1, kgi = <k|g|i> \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S2c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2, const CC_vecfunction &projector_t) const;

    /// for response equations the pair function looks like: tauij = uij + Qtf|xitj> + Qtf|tixj> - (OxQt + QtOx)f|titj>
    vecfuncT S2c_response_reg_part(const CC_vecfunction x, const CC_vecfunction &t, const CC_vecfunction& projector_t) const{
      CC_Timer time1(world,"S2c-Qt|xt>");
      const vecfuncT Qtxt = S2c_reg_part(x,t,projector_t);
      time1.info();
      CC_Timer time2(world,"S2c_Qt|tx>");
      const vecfuncT Qttx = S2c_reg_part(t,x,projector_t);
      time2.info();
      const vecfuncT function_response = add(world,Qtxt,Qttx);



      CC_Timer time3(world,"S2c-(OxQt+QtOx)|tt>");

      vecfuncT OxQt_part;
      vecfuncT QtOx_part;
      for(const auto& itmp:t.functions){
	const size_t i=itmp.first;
	const CC_function ti = t(i);

	real_function_3d OxQt_i = real_factory_3d(world);
	real_function_3d QtOx_i = real_factory_3d(world);

	for(const auto& ktmp:t.functions){
	  const size_t k=ktmp.first;
	  const CC_function tk=t(k);

	}
	OxQt_part.push_back(OxQt_i);
	QtOx_part.push_back(QtOx_i);
      }


      const vecfuncT projector_response = add(world,OxQt_part,QtOx_part);
      const vecfuncT result = sub(world,function_response,projector_response);
      time3.info();
      return result;


    }

    /// result: -\sum_{kl}( 2 <l|kgtaui|Qftktl> - <l|kgtaui|Qftltk>
    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ -|singles_l>(2.0*<lk|gQf|reg1_i,reg2_k>-<kl|gQf|reg1_i,reg2_k> \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4a_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles, const CC_vecfunction &projector) const;

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ -( <l*kgtaui|Qf|reg1_k,reg2_l>_2 - <l*kgtaui|Qf|reg1_k,reg2_l>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles, const CC_vecfunction & projector) const;

    /// result: 4<l|kgtauk|Qftitl> - 2<l|kgtauk|Qftlti> - 2<k|lgtauk|Qftitl> + <k|lgtauk|Qftlti>
    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ ( 4<l*kgtauk|Qf|reg1_i,reg2_l>_2 - 2<l*kgtauk|Qf|reg1_i,reg2_l>_1 - 2<k*lgtauk|Qf|reg1_i,reg2_l>_2 + <k*lgtauk|Qf|reg1_i,reg2_l>_1 ) \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles, const CC_vecfunction &projector) const;

    // tests if functions are suitable for a projector (same size as all active orbitals)
    void projector_functions_consistency(const CC_vecfunction &f)const{
      if(f.size()!=mo_ket_.size()){
	warning(f.name()+" functions not suitable for projector because the size differs from that of all MOs, maybe you gave a frozen intermediate");
      }
    }

    /// CC2 singles diagrams with 6d functions as input
    /// Use GFInterface in function_interface.h as kernel (f*g) and do not reconstruct \tau = f12u(1,2) if possible
    /// Since the correlation factor of CC2 has Slater form like in MP2: g12f12 = g12(1-exp(-mu*r12)/r12) = g12 - exp(-mu*r12)/r12 = Coulomb_Operator - BSH_Operator(mu)


    // CC2 Doubles Potential

    /// smalll helper function to track the time for the projetor
    void apply_Q12(real_function_6d &f,
		   const std::string &msg = "6D-function") const {
      CC_Timer Q12_time(world, "Applying Q12 to " + msg);
      f = projector_Q12(f);
      Q12_time.info();
    }
    vecfuncT apply_Q(const vecfuncT &f, const std::string msg="vectorfunction")const{
      CC_Timer time(world,"Q("+msg+")");
      vecfuncT result= projector_Q(f);
      time.info();
      return result;
    }
    real_function_3d apply_Q(const real_function_3d &f, const std::string msg="function")const{
      CC_Timer time(world,"Q("+msg+")");
      real_function_3d result= projector_Q(f);
      time.info();
      return result;
    }

    real_function_6d apply_regularization_potential(const CC_function &a, const CC_function &b, const double omega)const;

    // uses ansatz: tau = u + Qtf12(|xitj>+|tixj>) - (OxQt + QtOx)f12|titj>
    real_function_6d make_constant_part_cc2_response(const std::pair<size_t,size_t> &u,const CC_vecfunction &x,const CC_vecfunction&t ){
      const size_t i=u.first;
      const size_t j=u.second;
      const CC_function xi = x(i);
      const CC_function xj = x(j);
      const CC_function ti = t(i);
      const CC_function tj = t(j);
      const double omega = x.omega;
      const double bsh_eps = get_epsilon(i,j)+omega;

      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 *bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;

      output_section("Computing Commutator-Response of Constant-Part: (OVxQt + QtOVx - OxOVt - OVtOx)f|"+ti.name()+tj.name()+">");

      output_section("Computing Function-Response of Constant-Part: (Qt*Vreg + [F,Qt])|"+xi.name()+tj.name()+"+"+ti.name()+xj.name()+">");
      // QtVreg|titj> + [F,Qt]f|tifj> function response: |titj> --> >|xitj>+|tixj>
      const real_function_6d function_response_xt = make_constant_part_cc2_gs(xi,tj,t,omega);
      const real_function_6d function_response_tx = make_constant_part_cc2_gs(ti,xj,t,omega);
      const real_function_6d function_response = function_response_xt + function_response_tx;

      // Projector  and commutator response
      // -(OxQt+QtOx)Vregf|titj> + (OVxQt - OxOVt + QtOVx - OVtOx)f|titj>

      output_section("Computing Projector-Response of Constant-Part: -(QtOx+OxQt)*Vreg|"+ti.name()+tj.name()+">"+ " and projector-ket-states: " + t.name() +"\n");
      // this operation can be made into a complete 3D operation because of the projectors Ox
      const real_function_6d Vreg_tt = apply_regularization_potential(ti,tj,omega);
      // two intermediates
      const std::pair<vecfuncT,vecfuncT> Ox1_Vreg_tt = make_O1t_u(x,CC_function_6d(world,Vreg_tt));
      const std::pair<vecfuncT,vecfuncT> Ox2_Vreg_tt = make_O2t_u(x,CC_function_6d(world,Vreg_tt));
      // those are the needed functions
      const std::pair<vecfuncT,vecfuncT> OxQt_Vreg_tt = std::make_pair(Ox1_Vreg_tt.first,Q(Ox1_Vreg_tt.second,t));
      const std::pair<vecfuncT,vecfuncT> QtOx_Vreg_tt = std::make_pair(Q(Ox2_Vreg_tt.first,t),Ox2_Vreg_tt.second);

      // now apply the Greens Operator
      output("Raising 6D-thresh to " + std::to_string(parameters.tight_thresh_6D) + " in order to avoid accumulation of errors from 6D low-rank addition");
      FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);
      const real_function_6d projector_response_xt = apply_G_decomposed(G,OxQt_Vreg_tt);
      const real_function_6d projector_response_tx = apply_G_decomposed(G,QtOx_Vreg_tt);
      const real_function_6d projector_response = -1.0*projector_response_xt - projector_response_tx;
      projector_response_xt.print_size("OxQt*Vreg|"+ti.name()+tj.name()+">");
      projector_response_tx.print_size("QtOx*Vreg|"+ti.name()+tj.name()+">");
      projector_response.print_size("-(QtOx+OxQt)*Vreg|"+ti.name()+tj.name()+">");

      // now the Response from the Commutator;
      const vecfuncT Vx = current_singles_potential_response;
      vecfuncT Vt = current_singles_potential_gs;

      // for (CIS(D) there is no singles potential for the ground state (it is zero) so the parts including Vt will not be calculated
      const bool cispd = (t.type==HOLE);
      if(Vx.empty())warning("Constant-Part of CC2-Response: Vx is empty");
      if(not cispd and Vt.empty())warning("Constant-Part of CC2-Response: Vt is empty");
      if(cispd) output("t==HOLE: CIS(D)-Calculation!");

      // OVxQt parts
      const std::pair<vecfuncT,vecfuncT> OVx1_ftitj = make_O1t_op_xy(Vx,f12,ti,tj);
      const std::pair<vecfuncT,vecfuncT> OVx2_ftitj = make_O2t_op_xy(Vx,f12,ti,tj);
      const std::pair<vecfuncT,vecfuncT> OVxQt_ftitj = std::make_pair(OVx1_ftitj.first,Q(OVx1_ftitj.second,t));
      const std::pair<vecfuncT,vecfuncT> QtOVx_ftitj = std::make_pair(Q(OVx2_ftitj.first,t),OVx2_ftitj.second);

      // OVtOx parts (zero for CIS(D))
      std::pair<vecfuncT,vecfuncT> OxOVt_ftitj;
      std::pair<vecfuncT,vecfuncT> OVtOx_ftitj;
      if(not cispd){
      const std::pair<vecfuncT,vecfuncT> OVt1_ftitj = make_O1t_op_xy(Vt,f12,ti,tj);
      const std::pair<vecfuncT,vecfuncT> OVt2_ftitj = make_O2t_op_xy(Vt,f12,ti,tj);
      OxOVt_ftitj = std::make_pair(P(OVt1_ftitj.first,x),OVt1_ftitj.second);
      OVtOx_ftitj = std::make_pair(OVt2_ftitj.first,P(OVt2_ftitj.second,x));
      }else output("Skipping Vt parts bc. they are zero in CIS(D)");

      // apply the greens operator
      real_function_6d commutator_response=real_factory_6d(world);
      commutator_response += apply_G_decomposed(G,OVxQt_ftitj);
      commutator_response += apply_G_decomposed(G,QtOVx_ftitj);
      if(not cispd) commutator_response -= apply_G_decomposed(G,OxOVt_ftitj); // for CIS(D) the given pair is empty, the return is then just real_factory_6d(world) so it is initialized to zero
      if(not cispd) commutator_response -= apply_G_decomposed(G,OVtOx_ftitj); // same here

      // add all together and apply Q12
      real_function_6d constant_term = function_response+projector_response+commutator_response;
      constant_term.print_size("ConstTermResponse");
      function_response.print_size("function-response");
      commutator_response.print_size("commutator-response");
      apply_Q12(constant_term,"ConstTermResponse");
      constant_term.print_size("Q12ConstTermResponse");

      output("Lowering 6D-thresh back to " + std::to_string(parameters.thresh_6D));
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      constant_term.set_thresh(parameters.thresh_6D);

      return constant_term;
    }

    void test_on_mp2_constant_part(const CC_function &reg1, const CC_function &reg2){
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(reg1.i,reg2.i)),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(reg1.i,reg2.i)),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.destructive()=true;
      Gscreen.modified()=true;

      //fK|r1r2> part
      real_function_6d fK = real_factory_6d(world);
      CC_Timer time_fK(world,"fK-part");
      {
	CC_Timer time_K(world,"fK: Apply K");
      const real_function_3d Kreg1 = K(reg1);
      const real_function_3d Kreg2 = K(reg2);
      time_K.info();

      CC_Timer time_fill(world,"fK: fill-functions");
      real_function_6d V1 = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(reg1.function)).particle2(copy(Kreg2));
      real_function_6d V2 = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Kreg2)).particle2(copy(reg2.function));
      V1.fill_tree(Gscreen).truncate().reduce_rank();
      V2.fill_tree(Gscreen).truncate().reduce_rank();
      V1.print_size("fK: V1");
      V2.print_size("fK: V2");
      time_fill.info();

      CC_Timer time_G(world,"fK: Apply G");
      real_function_6d GV1 = -2.0*G(V1);
      real_function_6d GV2 = -2.0*G(V2);
      time_G.info();

      fK = GV1 + GV2;
      fK.print_size("fK: result");
      }time_fK.info();

      // Kf|r1r2> part
      real_function_6d Kf = real_factory_6d(world);
      CC_Timer time_Kf(world,"Kf-part");
      {
	CC_Timer time_fill(world,"Kf: fill-functions");
	real_function_6d freg1reg2 = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(reg1.function)).particle2(copy(reg2.function));
	freg1reg2.fill_tree(Gscreen).truncate().reduce_rank();
	time_fill.info();

	CC_Timer time_K(world,"Kf: Apply K");
	real_function_6d V = K(freg1reg2).truncate().reduce_rank();
	V.print_size("Kf: V");
	time_K.info();

	CC_Timer time_G(world,"Kf: Apply G");
	Kf = -2.0*G(V);
	time_G.info();

	Kf.print_size("Kf: result");
      }time_Kf.info();

      // Ue part
      CC_Timer time_Ue(world,"Ue-part");
      real_function_6d Ue = real_factory_6d(world);
      {
	CC_Timer time_U(world,"Ue: apply U");
	real_function_6d V = corrfac.apply_U(reg1.function,reg2.function,get_epsilon(reg1.i,reg2.i));
	V.print_size("Ue");
	time_U.info();

	CC_Timer time_G(world,"Ue: apply G");
	Ue = -2.0*G(V);
	Ue.print_size("Ue: result");
	time_G.info();
      }time_Ue.info();




    }



    // new version which uses: tau = u + Qt*f12*|titj>
    real_function_6d make_constant_part_cc2_gs(const CC_function &reg1, const CC_function &reg2, const CC_vecfunction &t, const double omega=0.0){
      // output convenience
      const std::string name = "|" + reg1.name()+reg2.name()+">";
      output("Calculating Constant-Part for CC2-Ground-State with functions: " + name + " and projector-ket-states: " + t.name() +"\n");
      if(parameters.test) test_on_mp2_constant_part(reg1,reg2);
      // consistency check
      if(omega==0.0 and (reg1.type==RESPONSE or reg2.type==RESPONSE)) warning("make_constant_part_cc2_gs: omega is zero but a response function was given");
      if(omega!=0.0 and (reg1.type!=RESPONSE and reg2.type!=RESPONSE)) warning("make_constant_part_cc2_gs: omega is not zero, but only ground state functions were given");
      if(t.size()!=mo_ket_.size())warning("make_constant_part_CC2_gs: The t-functions for the projector have not the same size as the unfrozen MOs");

      const bool symmetric=(reg1.i==reg2.i and reg1.type==reg2.type);

      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(reg1.i,reg2.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;

      CC_Timer time_V(world,"Apply Vreg to " + name);
      real_function_6d Vreg = apply_regularization_potential(reg1,reg2,0.0);
      Vreg.truncate().reduce_rank();
      Vreg.print_size("(Vreg).truncate()");
      time_V.info();

      // unprojected part: G(Vreg)
      CC_Timer time_GV(world,"GVreg"+name);
      const real_function_6d cVreg = copy(Vreg); // needed bc G is destructive
      const real_function_6d GV = -2.0*G(cVreg);
      GV.print_size("GVreg"+name);
      time_GV.info();

      output("Raising 6D-thresh to " + std::to_string(parameters.tight_thresh_6D) + " to avoid accumulation of addition errors");
      FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);

      CC_Timer timeG2(world,"GOt12Vreg"+name);
      // projected part is : -O1t(1-0.5O2t) - O2t(1-0.5O1t)
      const std::pair<vecfuncT,vecfuncT> O1t_Vreg = make_O1t_u(t,CC_function_6d(world,Vreg));
      const std::pair<vecfuncT,vecfuncT> O1_part = std::make_pair(O1t_Vreg.first,Q(O1t_Vreg.second,t,0.5));
      const real_function_6d GOV_O1 = apply_G_decomposed(G,O1_part);


      real_function_6d GOV_O2;
      if(symmetric) GOV_O2 = swap_particles(GOV_O1);
      else{
	const std::pair<vecfuncT,vecfuncT> O2t_Vreg = make_O2t_u(t,CC_function_6d(world,Vreg));
	const std::pair<vecfuncT,vecfuncT> O2_part = std::make_pair(Q(O2t_Vreg.first,t,0.5),O2t_Vreg.second);
	GOV_O2 = apply_G_decomposed(G,O2_part);
      }

      // negativ sign because this is projected out: Qt12 = (1-Ot1)(x)(1-Ot2)
      const real_function_6d GOV = -1.0*(GOV_O1 + GOV_O2);

      GV.print_size("GV"+name);
      GOV.print_size("GO12V"+name);
      timeG2.info();

      // fock commutator part [F,Qt12]f12|reg1reg2> = |Vk> (x) Qt <k|f|ti>*|tj> + Qt <k|f|tj>*|ti> (x) |Vk> , Vk = current singles potential: -Vk = (F-ek)|tauk>
      // the negative sign from the singles equation and the negative sign from the projector cancel each other

      real_function_6d GFQQF=real_factory_6d(world);
      if(t.type==HOLE){
	output("t-functions for projector are of typ: HOLE, the commutator [F,Qt] is therefore zero");
      }else{
	// little sanity check
	if(current_singles_potential_gs.empty())warning("Constant_part of cc2 needs the singles potential but it is empty");
	if((reg1.type==RESPONSE or reg2.type==RESPONSE) and current_singles_potential_response.empty()) warning("Constant_part of cc2-response needs the response-singles potential but it is empty");


	real_function_6d GFQQF_part1 = real_factory_6d(world);
	real_function_6d GFQQF_part2 = real_factory_6d(world);
	// raise thresh for many additions
	GFQQF_part1.set_thresh(parameters.tight_thresh_6D);
	GFQQF_part2.set_thresh(parameters.tight_thresh_6D);
	CC_Timer time_FQQF(world,"apply [F,Qt12]f12"+name);
	const CC_vecfunction OV(current_singles_potential_gs,UNDEFINED,parameters.freeze);
	const std::pair<vecfuncT,vecfuncT> O1Vfxy = make_O1t_op_xy(OV,f12,reg1,reg2);
	const std::pair<vecfuncT,vecfuncT> O1VQ2tfxy = std::make_pair(O1Vfxy.first,Q(O1Vfxy.second,t));
	GFQQF_part1 = apply_G_decomposed(G,O1VQ2tfxy);

	if(symmetric)  GFQQF_part2 = swap_particles(GFQQF_part1);
	else{
	  const std::pair<vecfuncT,vecfuncT> O2Vfxy = make_O2t_op_xy(OV,f12,reg1,reg2);
	  const std::pair<vecfuncT,vecfuncT> Q1tO2Vfxy = std::make_pair(Q(O2Vfxy.first,t),O2Vfxy.second); // here the particle swap is made
	  GFQQF_part2 = apply_G_decomposed(G,Q1tO2Vfxy);
	}
	GFQQF =  GFQQF_part1 +  GFQQF_part2;
	time_FQQF.info();
      }

      // sign of GOV has been applied before
      real_function_6d result = GV + GOV + GFQQF;
      result.print_size("(GV-GOV+GFQQF)");
      result.truncate().reduce_rank();
      GV.print_size("GV"+name);
      GOV.print_size("GOV"+name);
      GFQQF.print_size("GFQQF"+name);
      result.print_size("(GV-GOV+GFQQF).truncate()"+name);
      apply_Q12(result);
      result.print_size("Q12(GV-GOV+GFQQF)"+name);

      output("Lowering 6D-thresh back to " + std::to_string(parameters.thresh_6D));
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);

      return result;
    }

    /// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\tau i>+|i>
    // @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
    // @param[in] \tau_j
    /// \todo Parameter descriptions.
    /// @return Equation: \f$ Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj> \f$  with \f$ |ti> = |\tau i>+|i> \f$
    /// Right now Calculated in the decomposed form: \f$ |titj> = |i,j> + |\tau i,\tau j> + |i,\tau j> + |\tau i,j> \f$
    /// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
    real_function_6d make_regularization_residue(const CC_function &taui,
						 const CC_function &tauj, const calctype &type, const double omega=0.0) const;

    real_function_6d make_response_regularization_residue(const CC_function &ti, const CC_function &xi, const CC_function &tj, const CC_function &xj,const double &omega) const {
      output("Calculating GV_response_reg(|"+ti.name()+xj.name()+"> + |"+xi.name()+tj.name()+">)");
      // consistency check
      if(xi.type!=RESPONSE)error("response_Vreg: xi has wrong type");
      if(xj.type!=RESPONSE)error("response_Vreg: xj has wrong type");
      if(ti.type!=MIXED)error("response_Vreg: ti has wrong type");
      if(tj.type!=MIXED)error("response_Vreg: tj has wrong type");
      if(xi.i!=ti.i)error("response_Vreg: ti.i!=xi.i");
      if(xj.i!=tj.i)error("response_Vreg: tj.i!=xj.i");

      const bool symmetric(xi.i==xj.i and ti.i==tj.i and xi.i==tj.i);

      real_function_6d Vreg;
      const real_function_6d Vreg1 = apply_regularization_potential(ti,xj,omega);
      if(not symmetric) Vreg = Vreg1 + apply_regularization_potential(xi,tj,omega);
      else{
	output("Exploiting Symmetry for Diagonal Pairs");
	Vreg = Vreg1 + swap_particles(Vreg1);
      }

      Vreg.scale(-2.0);
      apply_Q12(Vreg,"Vreg");
      Vreg.truncate().reduce_rank();
      Vreg.print_size("-2.0*Q12Vreg");

      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(ti.i,tj.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;
      real_function_6d GV=G(Vreg);
      apply_Q12(GV,"CC2-Residue:G(V)");
      return GV;
    }



    real_function_6d make_nonorthogonal_regularization_residue(size_t i, size_t j)const{
      output("Calculating nonorthogonal GV_reg");
      const bool symmetric(i==j);

      real_function_6d Vreg;
      Vreg = apply_regularization_potential(mo_ket_(i),mo_ket_(i),0.0);
      Vreg.scale(-2.0);
      Vreg.truncate().reduce_rank();
      Vreg.print_size("-2.0*Vreg");

      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(i,j)),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;
      real_function_6d GV=G(Vreg);
      //apply_Q12(GV,"MP2-Residue:G(V)");
      return GV;
    }
    // apply the kinetic energy operator to a decomposed 6D function
    /// @param[in] y a 3d function x (will be particle 1 in the decomposed 6d function)
    /// @param[in] x a 3d function y (will be particle 2 in the decomposed 6d function)
    /// @return a 6d function: G(f12*T*|xy>)
    real_function_6d make_GQfT_xy(const real_function_3d &x,
				  const real_function_3d &y, const size_t &i, const size_t &j) const;

    // MP2 Terms are
    // fock_residue_6d = (2J - Kn + Un) |u_{ij}> , KR = R12^{-1}*K*R12 (nuclear tranformed K)
    // Uen|ij> = R12{-1}*U_e*R12 |ij>

    /// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
    real_function_6d fock_residue_6d(const CC_Pair &u) const;


    /// Exchange Operator on 3D function
    /// !!!!Prefactor (-1) is not included
    real_function_3d K(const CC_function &f) const;
    real_function_3d K(const real_function_3d &f) const;

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// if i==j in uij then the symmetry will be exploited
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d K(const real_function_6d &u,
		       const bool symmetric = false, const double thresh = FunctionDefaults<6>::get_thresh()) const;

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
    /// 1. X(3,2) = bra_k(3)*u(3,2)
    /// 2. Y(1,2) = \int X(3,2) g13 d3
    /// 3. result = Y(1,2)*ket_k(1)
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d apply_K(const real_function_6d &u,
			     const size_t &particle, const double thresh = FunctionDefaults<6>::get_thresh()) const;

    /// returns \sum_k (<k|g|f> *|k>).truncate()
    real_function_3d apply_K(const CC_function &f) const;

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
    /// @return  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
    real_function_6d apply_transformed_Ue(const CC_function &x,
					  const CC_function &y, const double &omega=0.0) const;

    /// Apply the Exchange Commutator [K,f]|xy>
    real_function_6d apply_exchange_commutator(const CC_function &x,
					       const CC_function &y, const double thresh = FunctionDefaults<6>::get_thresh()) const;

    /// Apply the Exchange operator on a tensor product multiplied with f12
    /// !!! Prefactor of (-1) is not inclued in K here !!!!
    real_function_6d apply_Kf(const CC_function &x, const CC_function &y, const double thresh = FunctionDefaults<6>::get_thresh()) const;

    /// Apply fK on a tensor product of two 3D functions
    /// fK|xy> = fK_1|xy> + fK_2|xy>
    /// @param[in] x the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    /// @param[in] y the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    real_function_6d apply_fK(const CC_function &x, const CC_function &y, const double thresh = FunctionDefaults<6>::get_thresh()) const;

    real_function_3d apply_laplacian(const real_function_3d &x) const;

    real_function_3d apply_F(const CC_function &x) const;
    real_function_3d apply_reduced_F(const CC_function &x)const;

    /// little helper function to pack a vector of CC_3D_functions (just structures which hold the function the index and the type)
    std::vector<CC_function> make_CC_3D_function(const vecfuncT &f,
						 const functype &type) {
      std::vector<CC_function> result(f.size());
      for (size_t i = 0; i < f.size(); i++) {
	CC_function tmp(f[i], i, type);
	result[i] = tmp;
      }
      return result;
    }

    // gives back \epsilon_{ij} = \epsilon_i + \epsilon_j
    double get_epsilon(const size_t &i, const size_t &j) const {
      return (orbital_energies[i] + orbital_energies[j]);
    }
    // gives back the orbital energies
    std::vector<double> get_orbital_energies() const {
      return orbital_energies;
    }
    /// swap particles 1 and 2

    /// param[in] all CC_Pairs
    /// param[in] the i index
    /// param[in] the j index
    /// param[out] a 6d function correspoding to electron pair ij
    /// if i>j the pair will be created via: fij(1,2) = fji(2,1)
    real_function_6d get_pair_function(const Pairs<CC_Pair> &pairs,
				       const size_t i, const size_t j) const {
      if (i > j) {
	const real_function_6d & function = pairs(j, i).function;
	const real_function_6d & swapped_function = swap_particles(
	    function);
	return swapped_function;

      } else {
	return pairs(i, j).function;
      }
    }

    CC_function_6d get_pair_function(const Pairs<CC_function_6d> &pairs,
				     const size_t i, const size_t j) const {
      if (i > j) {
	const CC_function_6d & function = pairs(j, i);
	const CC_function_6d & swapped_function = swap_particles(function);
	return swapped_function;
      } else {
	return pairs(i, j);
      }
    }

    /// param[in]	f	a function of 2 particles f(1,2)
    /// return	the input function with particles swapped g(1,2) = f(2,1)
    real_function_6d swap_particles(const real_function_6d& f) const;

    CC_function_6d swap_particles(const CC_function_6d& f) const{
      return f.swap_particles();
    }

    // Calculate the CC2 energy equation which is
    // \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
    // with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
    double get_CC2_correlation_energy() const;
    double compute_cc2_pair_energy(const CC_Pair &u, const CC_function &taui,
				   const CC_function &tauj) const;
    double compute_pair_energy(CC_Pair &u, const CC_vecfunction &singles)const{
      if(singles.type!=PARTICLE)error("compute_pair_energy with singles of type " + assign_name(singles.type) + " is not possible");
      CC_vecfunction t = make_t_intermediate_full(singles);

      // convenience
      const CC_function ti = t(u.i);
      const CC_function tj = t(u.j);
      const CC_function taui = singles(u.i);
      const CC_function tauj = singles(u.j);
      const CC_function mobrai = mo_bra_(u.i);
      const CC_function mobraj = mo_bra_(u.j);

      // 2<ij|g|taui,tauj> - <ji|g|taui,tauj> = i.inner(g12(j*tauj)*taui) - i.inner(g12(j*taui)*tauj)
      const double pure_singles_part = 2.0*mobrai.inner(g12(mobraj,tauj)*taui.function) - mobrai.inner(g12(mobraj,taui)*tauj.function);

      const double pure_u_part = 2.0*make_ijgu(mo_ket_(u.i),mo_ket_(u.j),u) - make_ijgu(mo_ket_(u.j),mo_ket_(u.i),u);

      const double regularization_part = 2.0*make_ijgQtfxy(mobrai,mobraj,t,ti,tj) - make_ijgQtfxy(mobraj,mobrai,t,ti,tj) ;

      const double result = pure_singles_part + pure_u_part + regularization_part;
      if(world.rank()==0){
	std::cout << "Correlation Energy of Pair " << u.name() << "\n";
	std::cout << std::fixed << std::setprecision(parameters.output_prec);
	std::cout << "<" << u.i << u.j << "|g|" << u.name() << ">_as         =" << pure_u_part << "\n";
	std::cout << "<" << u.i << u.j << "|g|" << taui.name() << tauj.name() << ">_as    =" << pure_singles_part << "\n";
	std::cout << "<" << u.i << u.j << "|gQtf|" << ti.name() << tj.name() << ">_as     =" << regularization_part << "\n";
	std::cout << "-------------------------\n";
	std::cout << "overall = " << result << "\n\n";
      }
      u.current_energy = result;
      u.constant_energy = pure_singles_part + regularization_part;
      return result;
    }
    /// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
    // the bra elements are always the R2orbitals
    // the ket elements can be \tau_i , or orbitals dependet n the type given
    double make_ij_gQf_ij(const size_t &i, const size_t &j, CC_Pair &u) const;
    double make_ijgQfxy(const size_t &i, const size_t &j, const CC_function &x,
			const CC_function &y) const;
    double make_ijgQfxy(const  CC_function &i, const  CC_function &j, const CC_function &x,
			const CC_function &y) const;
    // most general functon
    double make_ijgQtfxy(const CC_function &i, const CC_function &j, const CC_vecfunction &t, const CC_function &x, const CC_function &y)const{
      // the t vector must contain frozen orbitals
      if(t.size()!=mo_ket_.size()) error("Error in make_ijgQtfxy t-vector misses frozen orbitals -> important for projector");
      // part1: No projector, <ij|gf|xy>
      const real_function_3d ix = i.function*x.function;
      const real_function_3d jy = j.function*y.function;
      const real_function_3d gfix = apply_gf(ix);
      const double part1 = jy.inner(gfix);

      // part2: O1-part, <ij|gO1f|xy> = <j|igm*mfx|y>
      // part3: O2-part, <ij|gO3f|xy> = <i|jgm*mfy|x>
      // make part 4 from part2 ims
      // part4: O12-part,<ij|gO12f|xy>= ijgmn*mnfyx = <j|igm|n>*<n|mfx|y>
      double part2 = 0.0;
      double part3 = 0.0;
      double part4 = 0.0;
      for(const auto& mtmp:t.functions){
	const CC_function& mket =t(mtmp.first);
	const CC_function& mbra =mo_bra_(mtmp.first);
	const real_function_3d igm = g12(i,mket);
	const real_function_3d mfx = f12(mbra,x);
	part2 += jy.inner(igm*mfx);
	part3 += ix.inner(g12(j,mket)*f12(mbra,y));
	for(const auto& ntmp:t.functions){
	  const CC_function &nket = t(ntmp.first);
	  const CC_function &nbra = mo_bra_(ntmp.first);
	  const real_function_3d jn = j.function*nket.function;
	  const real_function_3d ny = nbra.function*y.function;
	  part4 += jn.inner(igm)*ny.inner(mfx);
	}
      }
      return part1-part2-part3+part4;
    }
    double make_ijgfxy(const size_t &i, const size_t &j,
		       const real_function_3d &x, const real_function_3d &y) const;
    /// Make two electron integral (expensive without intermediates) use just for debugging
    // double make_ijgxy(const size_t &i, const size_t &j,
    //const CC_function &x, const CC_function &y) const;
    double make_integral(const size_t &i, const size_t &j, const CC_function &x,
			 const CC_function&y, const optype type) const;
    /// Make two electron integral with the pair function
    double make_ijgu(const size_t &i, const size_t &j, const CC_Pair &u) const;
    double make_ijgu(const CC_function &phi_i, const CC_function &phi_j, const CC_Pair &u) const;
    double make_ijgu(const size_t &i, const size_t &j,
		     const real_function_6d &u) const;
    /// Make two electron integral with BSH operator
    double make_ijGu(const size_t &i, const size_t &j, const CC_Pair &u) const;
    /// apply the operator \f$ gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma) \f$
    /// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
    real_function_3d apply_gf(const real_function_3d &f) const;


    /// @param[in] x Function which is convoluted with (it is assumed that x is already multiplied with R2)
    /// @param[in] y function over which is not integrated
    /// @param[in] z function which is correlated with y over Q12f12
    /// @return <x(2)|Q12f12|y(1)z(2)>_2
    /// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
    /// 1. <x|f12|z>*|y>
    /// 2. - \sum_m mxfyz |m>
    /// 3. - \sum_m <x|m>*mfz*|y>
    /// 4. +\sum_{mn} <x|m> nmfyz |n>
    /// Description: Similar to convolute x_gQf_yz just without coulomb operator
    real_function_3d convolute_x_Qf_yz(const CC_function &x,
				       const CC_function &y, const CC_function &z) const;
    real_function_3d convolute_x_Qtf_yz(const CC_function &x,const CC_vecfunction &t,
					const CC_function &y, const CC_function &z) const{
      if(t.size()!=mo_ket_.size()) error("convolute_x_Qtf_yz, t for Qt has not the size of all orbitals -> maybe used frozen intermediate ?");
      real_function_3d xfz=f12(msparse(x.function,z.function));
      xfz.truncate();
      xfz.reconstruct();
      const real_function_3d xfz_y=msparse(xfz,y.function).truncate();
      const real_function_3d part1=msparse(xfz,y.function);

      real_function_3d part2=real_factory_3d(world);
      real_function_3d part3tmp=real_factory_3d(world);
      real_function_3d part4=real_factory_3d(world);
      for(const auto& mtmp : t.functions){
	const CC_function& mom=mtmp.second;
	const double mxfyz=mo_bra_(mom).function.inner(xfz_y);
	part2-=mxfyz * mom.function;

	const double xm=x.function.inner(mom.function);

	const real_function_3d mfz=f12(mo_bra_(mom),z);
	const real_function_3d mfz_y=msparse(mfz,y.function);

	part3tmp-=xm * mfz;

	for(const auto& ntmp : t.functions){
	  const CC_function& mon=ntmp.second;
	  const double nmfyz=mo_bra_(mon).function.inner(mfz_y);
	  part4+=xm * nmfyz * mon.function;
	}

      }
      const real_function_3d part3=msparse(part3tmp,y.function);
      real_function_3d result=part1 + part2 + part3 + part4;
      result.truncate();
      return result;

    }
    /// Doubles potentials
    /// G_D4b = G(Q12\sum_k <k|g|j>(1) |i\tau_k> + <k|g|i>(2) |\tau_k,j>)
    /// use that Q12 = Q1*Q2
    /// need to apply G to every single term in order to avoid entangelment


    real_function_6d make_xy(const CC_function &x, const CC_function &y, const double thresh=FunctionDefaults<6>::get_thresh() ) const;

    real_function_6d make_f_xy(const CC_function &x,
			       const CC_function &y,const double thresh = FunctionDefaults<6>::get_thresh()) const;

    real_function_6d make_f_xy_screened(const CC_function &x,
					const CC_function &y, const real_convolution_6d &screenG) const;


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
    const CC_vecfunction mo_bra_;
    const CC_vecfunction mo_ket_;
    /// The orbital energies
    const std::vector<double> orbital_energies;
    std::vector<double> init_orbital_energies(const Nemo &nemo) const {
      std::vector<double> eps;
      if (world.rank() == 0)
	std::cout << "SCF Orbital Energies are:\n";
      for (size_t i = 0; i < mo_ket_.size(); i++) {
	eps.push_back(nemo.get_calc()->aeps(i));
	if (world.rank() == 0)
	  std::cout << nemo.get_calc()->aeps(i);
      }
      if (world.rank() == 0)
	std::cout << "\n" << std::endl;
      return eps;
    }
    /// Helper function to initialize the const mo_bra and ket elements
    CC_vecfunction make_mo_bra(const Nemo &nemo) const {
      vecfuncT tmp = mul(world, nemo.nuclear_correlation->square(),
			 nemo.get_calc()->amo);
      set_thresh(world, tmp, parameters.thresh_3D);
      truncate(world,tmp);
      reconstruct(world,tmp);
      CC_vecfunction mo_bra(tmp, HOLE);
      return mo_bra;
    }

    CC_vecfunction make_mo_ket(const Nemo&nemo) const {
      vecfuncT tmp = nemo.get_calc()->amo;
      set_thresh(world, tmp, parameters.thresh_3D);
      truncate(world,tmp);
      reconstruct(world,tmp);
      CC_vecfunction mo_ket(tmp, HOLE);
      return mo_ket;
    }
    /// The poisson operator (Coulomb Operator)
    //    std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr
    //	< real_convolution_3d
    //	> (CoulombOperatorPtr(world, parameters.lo,
    //			      parameters.thresh_poisson));
    /// The BSH Operator for the f12g12 convolution which is with f12= 1/(2gamma)[1-exp(-gamma*r12)], f12g12 = 1/(2gamma) [CoulombOp - BSHOp(gamma)]
    std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr
	< real_convolution_3d
	> (BSHOperatorPtr3D(world, corrfac.gamma(), parameters.lo,
			    parameters.thresh_poisson));
    /// The f12 convolution operator
    //    std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr
    //	< real_convolution_3d
    //	> (SlaterF12OperatorPtr(world, corrfac.gamma(), parameters.lo,
    //				parameters.thresh_poisson));
    /// Intermediates (some need to be refreshed after every iteration)
    //CC_Intermediates intermediates_;
    /// The current singles potential (Q\sum singles_diagrams) , needed for application of the fock opeerator on a singles function
    vecfuncT current_singles_potential_gs;
    vecfuncT current_singles_potential_response;
    /// The 6D part of S2c and S2b (only the parts which depends on the u-function, not the regularization tail (since the singles change)
    mutable vecfuncT current_s2b_u_part_gs;
    mutable vecfuncT current_s2b_u_part_response;
    mutable vecfuncT current_s2c_u_part_gs;
    mutable vecfuncT current_s2c_u_part_response;
    mutable vecfuncT current_s2b_reg_part_gs;
    mutable vecfuncT current_s2b_reg_part_response;

    CC_convolution_operator g12;
    CC_convolution_operator f12;
    StrongOrthogonalityProjector<double, 3> projector_Q12;
    QProjector<double,3> projector_Q;

  public:
    void check_stored_singles_potentials() {

    }
    void remove_stored_singles_potentials() {
      output("Removing stored singles potentials\n");
      current_s2b_u_part_gs.clear();
      current_s2c_u_part_gs.clear();
      current_singles_potential_gs.clear();
    }
    void remove_stored_response_singles_potentials() {
      output("Removing stored singles potentials\n");
      current_s2b_u_part_response.clear();
      current_s2c_u_part_response.clear();
      current_singles_potential_response.clear();
    }



    void screening(const real_function_3d &x, const real_function_3d &y) const {
      double normx = x.norm2();
      double normy = y.norm2();
      double norm_xy = normx * normy;
      if (world.rank() == 0)
	std::cout << "Screening |xy> 6D function, norm is: " << norm_xy
	<< std::endl;
      //return norm_xy;
    }

    double guess_thresh(const CC_function &x, const CC_function &y) const {
      double norm = x.function.norm2() * y.function.norm2();
      double thresh = parameters.thresh_6D;
      if (norm > parameters.thresh_6D)
	thresh = parameters.thresh_6D;
      else thresh = parameters.tight_thresh_6D;
      return thresh;
    }

    // Debug function, content changes from time to time
    void test_potentials(const int k, const double thresh =
	FunctionDefaults<3>::get_thresh(), const bool refine = true) const {
      output_section("Little Debug and Testing Session");
      // testing laplace operator with different thresholds
      const size_t old_k = FunctionDefaults<3>::get_k();
      {
	CC_Timer time(world, "time");
	FunctionDefaults<3>::set_thresh(thresh);
	FunctionDefaults<3>::set_k(k);
	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	std::string name = "_" + stringify(k) + "_" + stringify(thresh);
	if (world.rank() == 0)
	  std::cout
	  << "Testing Laplace operator with threshold  "
	  + stringify(FunctionDefaults<3>::get_thresh())
	  + " and k="
	  + stringify(FunctionDefaults<3>::get_k())
	  + " and refinement=" + stringify(refine) + "\n";

	real_function_3d gauss = real_factory_3d(world).f(f_gauss);
	real_function_3d laplace_gauss_analytical =
	    real_factory_3d(world).f(f_laplace_gauss);
	real_function_3d laplace_gauss_analytical_old_k = project(
	    laplace_gauss_analytical, old_k);
	plot_plane(world, gauss, "gauss" + name);
	plot_plane(world, laplace_gauss_analytical,
		   "laplace_gauss_analytical_old_k" + name);

	real_function_3d laplace_gauss_numerical = real_factory_3d(world);
	for (size_t i = 0; i < 3; i++) {
	  real_function_3d tmp = (*gradop[i])(gauss);
	  real_function_3d tmp2 = (*gradop[i])(tmp);
	  laplace_gauss_numerical += tmp2;
	}

	real_function_3d laplace_gauss_diff = laplace_gauss_analytical
	    - laplace_gauss_numerical;
	plot_plane(world, laplace_gauss_diff, "laplace_gauss_diff" + name);
	laplace_gauss_diff.print_size(
	    "||laplace on gauss num and ana      ||");

	real_function_3d projected_numerical = project(
	    laplace_gauss_numerical, old_k);
	real_function_3d laplace_gauss_diff2 =
	    laplace_gauss_analytical_old_k - projected_numerical;
	plot_plane(world, laplace_gauss_diff2,
		   "laplace_gauss_diff_old_k" + name);
	laplace_gauss_diff.print_size(
	    "||laplace on gauss num and ana old k||");

	FunctionDefaults<3>::set_thresh(parameters.thresh_3D);
	FunctionDefaults<3>::set_k(old_k);
	world.gop.fence();
	time.info();
      }
    }

    real_function_3d smooth_function(const real_function_3d &f,
				     const size_t mode) const {
      size_t k = f.get_impl()->get_k();
      real_function_3d fproj = project(f, k - 1);
      real_function_3d freproj = project(fproj, k);
      real_function_3d smoothed2 = 0.5 * (f + freproj);
      // double diff = (freproj - f).norm2();
      // double diff2 = (smoothed2 - f).norm2();
      // if(world.rank()==0) std::cout << "||f - f_smoothed|| =" << diff << std::endl;
      // if(world.rank()==0) std::cout << "||f - f_smoothed2||=" << diff2 << std::endl;
      if (mode == 1)
	return freproj;
      else if (mode == 2)
	return smoothed2;
      else {
	std::cout << "Unknown smoothing mode, returning unsmoothed function"
	    << std::endl;
	return f;
      }
    }



    template<size_t NDIM>
    void test_greens_operators(const double thresh, const size_t k, const double eps, const double bsh_thresh = 1.e-6)const{
      std::cout << "\n\nTesting " << NDIM << "-dimensional Greens Operator with thresh=" << thresh << " BSH-thresh=" << bsh_thresh << " and epsilon=" << eps<< std::endl;
      FunctionDefaults<NDIM>::set_k(k);
      FunctionDefaults<NDIM>::set_thresh(thresh);
      Function<double,NDIM> f = FunctionFactory<double,NDIM>(world).f(gauss_ND<NDIM>);
      f.print_size("TestGaussFunction");
      //Function<double,NDIM> one = FunctionFactory<double,NDIM>(world).f(unitfunction<NDIM>);
      SeparatedConvolution<double,NDIM> G = BSHOperator<NDIM>(world, sqrt(-2.0 * eps),1.e-6, 1.e-5);
      Function<double,NDIM> Lf = general_apply_laplacian<NDIM>(f);
      Lf.print_size("Laplace(f)");
      // Helmholtz eq: (Delta + 2eps)f = ...
      Lf += 2.0*eps*f;
      Lf.truncate();
      Lf.print_size("(Laplace +2eps)f");
      Function<double,NDIM> GLf= G(Lf);

      //const double integral_f = one.inner(f);
      //const double integral_GLf = one.inner(GLf);

      //std::cout << "integral(f)  =" << integral_f   << std::endl;
      //std::cout << "integral(GLf)=" << integral_GLf << std::endl;
      std::cout << "<f|f>     = " << f.inner(f)     << std::endl;
      std::cout << "<GLf|f>   = " << GLf.inner(f)   << std::endl;
      std::cout << "<f|GLf>   = " << f.inner(GLf)   << std::endl;
      std::cout << "<GLf|GLf> = " << GLf.inner(GLf) << std::endl;
      std::cout << "\n\n";


    }



    template<size_t NDIM>
    Function<double,NDIM> general_apply_laplacian(const Function<double,NDIM> &f)const{
      Function<double,NDIM> laplace_f = FunctionFactory<double,NDIM>(world);
      for (int axis = 0; axis < NDIM; ++axis) {
	Derivative<double,NDIM> D = free_space_derivative<double, NDIM>(world,axis);
	const Function<double,NDIM> Df = D(f);
	const Function<double,NDIM> D2f= D(Df).truncate();
	laplace_f += D2f;
      }
      return laplace_f;
    }

    double size_of(const CC_vecfunction &f)const{
      double size = 0.0;
      for(const auto &tmp:f.functions){
	size += get_size<double,3>(tmp.second.function);
      }
      return size;
    }

    double size_of(const Pairs<CC_Pair> &pairs)const{
      double size =0.0;
      for(const auto & pair:pairs.allpairs){
	size += get_size<double,6>(pair.second.function);
      }
      return size;
    }



    void print_memory_information(const CC_vecfunction &singles, const Pairs<CC_Pair> &doubles)const{
      const double moket_size = size_of(mo_ket_);
      const double mobra_size = size_of(mo_bra_);
      const double singles_size = size_of(singles);
      const double doubles_size = size_of(doubles);
      const double intermediates_size = g12.info() + f12.info();
      const double all = singles_size + doubles_size + intermediates_size;
      if(world.rank()==0){
	std::cout << "MO-ket       :" << moket_size      <<" (GByte)"<< std::endl;
	std::cout << "MO-bra       :" << mobra_size      <<" (GByte)"<< std::endl;
	std::cout << "Singles      :" << singles_size    <<" (GByte)"<< std::endl;
	std::cout << "Singles      :" << singles_size    <<" (GByte)"<< std::endl;
	std::cout << "Doubles      :" << doubles_size    <<" (GByte)"<< std::endl;
	std::cout << "Intermediates:" << intermediates_size    <<" (GByte)"<< std::endl;
	std::cout << "--------      " << "----------------------\n";
	std::cout << "all           " << all              <<" (GByte)"<< std::endl;
	std::cout << "\n----------------------------------\n";
      }
    }


    void plot(const CC_vecfunction &vf)const{
      CC_Timer time(world,"Plotting " +vf.name());
      for(const auto& tmp:vf.functions){
	plot(tmp.second);
      }
      time.info();
    }
    void plot(const CC_function &f)const{
      plot_plane(world,f.function,f.name());
    }
    void plot(const Pairs<CC_Pair> &pairs)const{
      for(const auto& tmp:pairs.allpairs){
	plot(tmp.second);
      }
    }
    void plot(const CC_Pair &f)const{
      plot_plane(world,f.function,f.name());
    }

    // omega is an excitation energy
    // this function checks if omega should be zero (ground state calculation)
    void consistency_check(const CC_function &a, const CC_function&b, const double omega)const{
      if((a.type!=RESPONSE and b.type!=RESPONSE) and omega !=0.0){
	error("Inconsistency detected: " + a.name() + b.name() + " and omega nonzero " +std::to_string(omega) );
      }
      if((a.type==RESPONSE or b.type==RESPONSE) and omega ==0.0){
	error("Inconsistency detected: " + a.name() + b.name() + " and omega zero " +std::to_string(omega) );
      }
    }

    // make: result = Q12u - O1f|xy> - O2f|xy> + O12f|xy>
    real_function_6d do_f12_projection(const real_function_6d &u,const CC_function &x, const CC_function &y)const{
      CC_Timer time(world,"f12-projection");
      output("Now Doing f12_projection");
      const real_function_6d Qu = projector_Q12(u);

      real_function_6d  Ofxy= real_factory_6d(world);
      Ofxy.set_thresh(parameters.tight_thresh_6D);
      for(const auto& mtmp:mo_ket_.functions){
	const size_t m = mtmp.first;
	const real_function_3d mfx = f12(mo_bra_(m),x);
	const real_function_3d my = mo_bra_(m).function*y.function;
	const real_function_3d mfx_y = mfx*y.function;
	const real_function_3d mfy = f12(mo_bra_(m),y);
	const real_function_3d mfy_x = mfy*x.function;

	real_function_3d im1 = real_factory_3d(world);
	real_function_3d im2 = real_factory_3d(world);
	for(const auto& ntmp:mo_ket_.functions){
	  const size_t n = ntmp.first;
	  const real_function_3d ny = mo_bra_(n).function*y.function;
	  const real_function_3d nx = mo_bra_(n).function*x.function;
	  const double mnfxy = ny.inner(mfx);
	  const double nmfxy = nx.inner(mfy);
	  im2 += mnfxy*mo_ket_(n).function;
	  im1 += nmfxy*mo_ket_(n).function;
	}
	im1.scale(-0.5);
	im2.scale(-0.5);
	const real_function_3d particle2 = mfx_y + im2;
	const real_function_3d particle1 = mfy_x + im1;
	const real_function_6d O1part = make_xy(mo_ket_(m),CC_function(particle2,99,UNDEFINED),parameters.tight_thresh_6D);
	const real_function_6d O2part = make_xy(CC_function(particle1,99,UNDEFINED),mo_ket_(m),parameters.tight_thresh_6D);
	Ofxy += O1part + O2part;
      }
      const real_function_6d result = Qu - Ofxy;
      if(world.rank()==0){
	std::cout << "\n\n u=Qu-O12f12|"+x.name()+y.name()+"> results:\n";
	std::cout << "||u||     =" << u.norm2() << "\n";
	std::cout << "||Qu||    ="<< Qu.norm2() << "\n";
	std::cout << "||Ofxy||  ="<< Ofxy.norm2() << "\n";
	std::cout << "||QfU||   ="<< result.norm2() << "\n";
	std::cout << "||diff||  ="<< (u-result).norm2() << "\n\n";
      }
      time.info();
      return result;
    }


    real_function_6d make_G_P_g_xy(const CC_vecfunction &p, const CC_function &x, const CC_function &y, const double omega =0.0)const{
      return make_G_P_op_xy(p,g12,x,y,omega);
    }
    real_function_6d make_G_P_op_xy(const CC_vecfunction &p, const CC_convolution_operator &op, const CC_function &x, const CC_function &y, const double omega =0.0)const{

      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(x.i,y.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;

      real_function_6d result = real_factory_6d(world);
      result.set_thresh(parameters.tight_thresh_6D);
      for(const auto &ptmp:p.functions){
	const size_t k=ptmp.first;
	const CC_function& pk = ptmp.second;

	real_function_3d particle1 = pk.function;
	real_function_3d particle2 = op(mo_bra_(k),x)*y.function;
	particle1 = apply_Q(particle1);
	particle2 = apply_Q(particle2);
	result += -2.0*G(particle1,particle2);
      }
      result.set_thresh(parameters.thresh_6D);
      return result;
    }

    real_function_6d make_G_P1P2_g_xy(const CC_vecfunction &p1, const CC_vecfunction &p2,const CC_function &x, const CC_function &y, const double omega =0.0)const{
      return make_G_P1P2_op_xy(p1,p2,g12,x,y,omega);
    }

    real_function_6d make_G_P1P2_op_xy(const CC_vecfunction &p1, const CC_vecfunction &p2, const CC_convolution_operator &op, const CC_function &x, const CC_function &y, const double omega =0.0)const{
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * get_epsilon(x.i,y.i)+omega),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;

      Tensor<double> kl_op_xy = make_matrix_mn_Op_xy(mo_bra_,mo_bra_,op,x,y);

      real_function_6d part1 = real_factory_6d(world);
      for(const auto& p1tmp:p1.functions){
	real_function_3d particle1 = p1tmp.second.function;
	real_function_3d particle2 = real_factory_3d(world);
	for(const auto& p2tmp:p2.functions){
	  particle2 += 0.5*kl_op_xy(p1tmp.first,p2tmp.first)*p2tmp.second.function;
	}
	part1 += -2.0*G(particle1,particle2);
      }

      real_function_6d part2 = real_factory_6d(world);
      for(const auto& p2tmp:p1.functions){
	real_function_3d particle2 = p2tmp.second.function;
	real_function_3d particle1 = real_factory_3d(world);
	for(const auto& p1tmp:p2.functions){
	  particle1 += 0.5*kl_op_xy(p1tmp.first,p2tmp.first)*p1tmp.second.function;
	}
	part2 += -2.0*G(particle1,particle2);
      }

      return part1 + part2;
    }

    // GO1Q2op(xy) = GO1op(xy) - GO1O2op(xy)
    real_function_6d make_G_O1Q2_op_xy(const CC_vecfunction &O1, const CC_vecfunction&Q2, const CC_convolution_operator &op, const CC_function &x, const CC_function &y, const double omega=0.0)const{
      const real_function_6d GO1opxy = make_G_P_op_xy(O1,op,x,y,omega);
      const real_function_6d GO1O2opxy = make_G_P1P2_op_xy(O1,Q2,op,x,y,omega);
      return GO1opxy - GO1O2opxy;
    }

    real_function_6d make_O12_op_xy(const CC_convolution_operator &op, const CC_function &x, const CC_function &y)const{
      Tensor<double> mn_op_xy = make_matrix_mn_Op_xy(mo_bra_,mo_bra_,op,x,y);
      // make intermediate for particle2 of O1(1-0.5O2) part
      vecfuncT O1_particle1 = mo_ket_.get_vecfunction();
      vecfuncT O1_particle2;
      vecfuncT O2_particle2 = mo_ket_.get_vecfunction();
      vecfuncT O2_particle1;
      for(size_t m=0;m<mo_ket_.size();m++){
	real_function_3d O1p2_m = op(mo_bra_(m),x)*y.function;
	real_function_3d O2p1_m = op(mo_bra_(m),y)*x.function;
	for(size_t n=0;n<mo_ket_.size();n++){
	  O1p2_m -= 0.5*mn_op_xy(m,n)*mo_ket_(n).function;
	  O2p1_m -= 0.5*mn_op_xy(n,m)*mo_ket_(n).function;
	}
	O1_particle2.push_back(O1p2_m);
	O2_particle1.push_back(O2p1_m);
      }
      if(O1_particle1.size()!=O1_particle2.size())error("inconsistent sizes in O1-part of O12"+op.name()+x.name()+y.name());
      if(O2_particle1.size()!=O2_particle2.size())error("inconsistent sizes in O2-part of O12"+op.name()+x.name()+y.name());
      if(O1_particle1.size()!=O2_particle1.size())error("inconsitent sizes in "+op.name()+x.name()+y.name());
      real_function_6d result = real_factory_6d(world);
      result.set_thresh(parameters.tight_thresh_6D);
      for(size_t i=0;i<O1_particle1.size();i++){
	result += make_xy(O1_particle1[i],O1_particle2[i],parameters.tight_thresh_6D);
	result += make_xy(O2_particle1[i],O2_particle2[i],parameters.tight_thresh_6D);
      }
      return result;
    }

    Tensor<double> make_matrix_mn_Op_xy(const CC_vecfunction &m, const CC_vecfunction &n, const CC_convolution_operator &op ,const CC_function &x, const CC_function &y)const{
      bool mhole = (m.type==HOLE);
      bool nhole = (n.type==HOLE);
      bool same_size =(m.size()==n.size());
      if(not(mhole and nhole and same_size)) error("Matrix for m,n not HOLE states or different sizes not possible");

      Tensor<double> result(m.size(),n.size());
      for(size_t i=0;i<m.size();i++){
	const real_function_3d mfx = op(mo_bra_(i),x);
	const real_function_3d mfx_y = mfx*y.function;
	for(size_t j=0;j<n.size();j++){
	  const double mnfxy = mo_bra_(j).inner(mfx_y);
	  result(i,j)=mnfxy;
	}
      }
      return result;
    }

    // Ot(1)f|xy> = |tk(1)>(x) <k|f|x>*|y>
    // t has always the size of all mos (since the frozen t-functions are just the mos)
    std::pair<vecfuncT,vecfuncT> make_O1t_op_xy(const CC_vecfunction &t, const CC_convolution_operator &op, const CC_function &x, const CC_function &y)const{
      if(t.size()==0) output("Ot1|"+x.name()+y.name()+"> : FUNCTIONS FOR PROJECTOR ARE EMPTY");
      vecfuncT particle1 = copy(world,t.get_vecfunction()); // avoid shallow copy errors (in most cases a destructive G-operator will be applied to this)
      vecfuncT particle2;
      for(const auto & ktmp:t.functions){
	const size_t k=ktmp.first;
	const real_function_3d p2tmp = op(mo_bra_(k),x)*y.function;
	particle2.push_back(p2tmp);
      }

      return std::make_pair(particle1,particle2);
    }
    std::pair<vecfuncT,vecfuncT> make_O2t_op_xy(const CC_vecfunction &t, const CC_convolution_operator &op, const CC_function &x, const CC_function &y)const{
      if(t.size()==0) output("Ot2|"+x.name()+y.name()+"> : FUNCTIONS FOR PROJECTOR ARE EMPTY");
      vecfuncT particle1;
      vecfuncT particle2= copy(world,t.get_vecfunction()); // avoid shallow copy errors (in most cases a destructive G-operator will be applied to this)
      for(const auto & ktmp:t.functions){
	const size_t k=ktmp.first;
	const real_function_3d p1tmp = op(mo_bra_(k),y)*x.function;
	particle1.push_back(p1tmp);
      }

      return std::make_pair(particle1,particle2);
    }

    std::pair<vecfuncT,vecfuncT> make_O1t_u(const CC_vecfunction &t, const CC_function_6d &u)const{
      if(t.size()==0) output("Ot1|V6D> : FUNCTIONS FOR PROJECTOR ARE EMPTY");
      vecfuncT particle1 = copy(world,t.get_vecfunction()); // avoid shallow copy errors (in most cases a destructive G-operator will be applied to this)
      vecfuncT particle2;
      for(const auto &ktmp:t.functions){
	const size_t k=ktmp.first;
	const real_function_3d p2tmp = u.project_out(mo_bra_(k),1);
	particle2.push_back(p2tmp);
      }
      return std::make_pair(particle1,particle2);
    }
    std::pair<vecfuncT,vecfuncT> make_O2t_u(const CC_vecfunction &t, const CC_function_6d &u)const{
      if(t.size()==0) output("Ot2|V6D> : FUNCTIONS FOR PROJECTOR ARE EMPTY");
      vecfuncT particle2 = copy(world,t.get_vecfunction()); // avoid shallow copy errors (in most cases a destructive G-operator will be applied to this)
      vecfuncT particle1;
      for(const auto &ktmp:t.functions){
	const size_t k=ktmp.first;
	const real_function_3d p1tmp = u.project_out(mo_bra_(k),2);
	particle1.push_back(p1tmp);
      }
      return std::make_pair(particle1,particle2);
    }

    real_function_6d apply_G_decomposed(const real_convolution_6d &G,const std::pair<vecfuncT,vecfuncT> &f)const{
      CC_Timer time(world,"Apply_G_decomposed");
      if(f.first.size()!=f.second.size()) error("Apply_G_decomposed: Different sizes in given vectorfunctions");
      // tighten thresh from the outside so that we dont loose const
      if(FunctionDefaults<6>::get_thresh()!=parameters.tight_thresh_6D) warning("Thresh for apply_G_decomposed was not tightened");
      real_function_6d result = real_factory_6d(world);
      for(size_t i=0;i<f.first.size();i++){
	result += G(f.first[i],f.second[i]);
      }
      result.scale(-2.0);
      time.info();
      return result;
    }

    const Pairs<CC_function_6d> make_pairs(const Pairs<CC_Pair> &u)const{
      Pairs<CC_function_6d> result;
      for(const auto & tmp:u.allpairs){
	CC_function_6d tmpa(world,tmp.second.function);
	result.insert(tmp.second.i,tmp.second.j,tmpa);
      }
      return result;
    }

    // gives back pairs like this: f12|reg1,reg2> (not as full 6d function)
    Pairs<CC_function_6d> make_regularization_tails_f12(const CC_vecfunction &reg1, const CC_vecfunction &reg2)const{
      Pairs<CC_function_6d> result;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
	for(size_t j=i;j<mo_ket_.size();j++){
	  CC_function_6d tmp(world,&f12,reg1(i),reg2(j));
	  result.insert(i,j,tmp);
	}
      }
      return result;
    }

    // gives back pairs like this:  O1(1-factor*O2)f12|reg1,reg2>
    Pairs<CC_function_6d> make_regularization_tails_O1(const CC_vecfunction &projector1,const CC_vecfunction &projector2, const CC_vecfunction &reg1, const CC_vecfunction &reg2, const double factor)const{
      // factor should be 0.5 or 1
      if(factor!=0.5 and factor!=1.0) warning("Factor of " + std::to_string(factor) + " in make_regularization_tails_O1");
      Pairs<CC_function_6d> result;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
	for(size_t j=i;j<mo_ket_.size();j++){
	  const std::pair<vecfuncT,vecfuncT> O1fxy = make_O1t_op_xy(projector1,f12,reg1(i),reg2(j));
	  const vecfuncT p1 = O1fxy.first;
	  const vecfuncT p2 = Q(O1fxy.second,projector2,factor);
	  CC_function_6d tmp(world,p1,p2);
	  result.insert(i,j,tmp);
	}
      }
      return result;
    }

    // gives back pairs like this:  (1-factor*O1)O2f12|reg1,reg2>
    Pairs<CC_function_6d> make_regularization_tails_O2(const CC_vecfunction &projector1,const CC_vecfunction &projector2, const CC_vecfunction &reg1, const CC_vecfunction &reg2, const double factor)const{
      // factor should be 0.5 or 1
      if(factor!=0.5 and factor!=1.0) warning("Factor of " + std::to_string(factor) + " in make_regularization_tails_O1");
      Pairs<CC_function_6d> result;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
	for(size_t j=i;j<mo_ket_.size();j++){
	  const std::pair<vecfuncT,vecfuncT> O2fxy = make_O2t_op_xy(projector2,f12,reg1(i),reg2(j));
	  const vecfuncT p1 = Q(O2fxy.first,projector1,factor);
	  const vecfuncT p2 = O2fxy.second;
	  CC_function_6d tmp(world,p1,p2);
	  result.insert(i,j,tmp);
	}
      }
      return result;
    }


    bool test_f12_projections()const{
      const real_function_6d fij = make_f_xy(mo_ket_(parameters.freeze),mo_ket_(parameters.freeze));
      real_function_6d Qfij_1 = copy(fij);
      apply_Q12(Qfij_1);
      const real_function_6d zero = real_factory_6d(world);
      const real_function_6d Ofij_1 = do_f12_projection(zero,mo_ket_(parameters.freeze),mo_ket_(parameters.freeze)); // gets 0 - Ofij back
      const real_function_6d Ofij_2 = make_O12_op_xy(f12,mo_ket_(parameters.freeze),mo_ket_(parameters.freeze));
      const real_function_6d diff_Ofij = Ofij_1 + Ofij_2;
      const real_function_6d Qfij_2 = fij + Ofij_1; // minus sign included in function before
      const real_function_6d diff = Qfij_1 - Qfij_2;
      Qfij_1.print_size("Qfij_1");
      Qfij_2.print_size("Qfij_2");
      diff.print_size("difference");
      if(world.rank()==0){
	std::cout << "\n\nEnd of f12-projection Test:\n";
	std::cout << "||Qfij_1||=" << Qfij_1.norm2() << "\n";
	std::cout << "||Qfij_2||=" << Qfij_2.norm2() << "\n";
	std::cout << "||differ||=" << diff.norm2() << "\n\n";
	std::cout << "difference between Ofij from different functions: " << diff_Ofij.norm2() <<"\n\n"<< std::endl;
      }
      if(diff.norm2()<parameters.thresh_6D) return true;
      else return false;
    }

    bool test_inverse_correlation_factor()const{
      if(corrfac.gamma()!=0.5){
	output("Gamma of correlationfactor is not 1/2!");
	return false;
      }

      CorrelationFactor2 corrfac2(world);
      real_function_3d x = mo_ket_(parameters.freeze).function;
      real_function_3d y = copy(x);
      const real_function_6d xy = make_xy(x,y,parameters.thresh_6D);
      CC_Timer time1(world,"make fxy");
      real_function_6d fxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(y));
      fxy.fill_tree().truncate().reduce_rank();
      time1.info();
      CC_Timer time2(world,"make f2xy");
      real_function_6d f2xy = CompositeFactory<double,6,3>(world).g12(corrfac2.function()).particle1(copy(x)).particle2(copy(y));
      f2xy.fill_tree().truncate().reduce_rank();
      time2.info();
      real_function_6d f2xy_2 = 0.5*(xy+fxy);
      real_function_6d diff_f2 = f2xy_2 - f2xy;
      xy.print_size("|xy>");
      fxy.print_size("f|xy>");
      f2xy.print_size("f2|xy>");
      f2xy_2.print_size("0.5(1+f)|xy>");
      diff_f2.print_size("diff_f2");

      CC_Timer time3(world,"make (1/f2)*f2|xy>");
      real_function_6d xy_2 = CompositeFactory<double,6,3>(world).g12(corrfac2.inverse()).ket(f2xy);
      xy_2.fill_tree().truncate().reduce_rank();
      time3.info();

      real_function_6d diff_xy = xy - xy_2;
      xy.print_size("|xy>");
      xy_2.print_size("(1/f2)*f2|xy>");
      diff_xy.print_size("diff_|xy>");

      if(world.rank()==0){
	std::cout << "Test corrfactor2 and inverse of it:\n";
	std::cout << "diff_f2 =" << diff_f2.norm2() << "\n";
	std::cout << "diff_inv=" << diff_xy.norm2() << "\n";
      }

      if(diff_f2.norm2()<parameters.thresh_6D and diff_xy.norm2()<parameters.thresh_6D) return true;
      else return false;

    }

    bool test_energy_functions()const{

      const double testa1 = make_ijgQfxy(mo_bra_(0),mo_bra_(0),mo_ket_(0),mo_ket_(0));
      const double testa2 = make_ijgQfxy(0,0,mo_ket_(0),mo_ket_(0));
      const double testa3 = make_ijgQtfxy(mo_bra_(0),mo_bra_(0),mo_ket_,mo_ket_(0),mo_ket_(0));
      const double diff12 = testa1 - testa2;
      const double diff13 = testa1 - testa3;
      if(world.rank()==0){
	std::cout << "Test of energy functions:\n";
	std::cout << std::fixed << std::setprecision(parameters.output_prec);
	std::cout << "<00|gQf|00>=" << testa1 << "\n";
	std::cout << "<00|gQf|00>=" << testa2 << "\n";
	std::cout << "<00|gQf|00>=" << testa3 << "\n";
      }
      if(fabs(diff12)>parameters.thresh_3D or fabs(diff13)>parameters.thresh_3D){
	warning("Energy test failed!");
	return false;
      } else output("Energy Test passed\n\n");

      output("Testing convolute_x_Qtf_yz functions\n");

      const real_function_3d testb1 = convolute_x_Qf_yz(mo_ket_(0),mo_ket_(0),mo_ket_(0));
      const real_function_3d testb2 = convolute_x_Qtf_yz(mo_ket_(0),mo_ket_,mo_ket_(0),mo_ket_(0));
      const double diff=(testb1-testb2).norm2();
      if(world.rank()==0) std::cout << "||old-new||=" << diff << "\n";
      if(diff>parameters.thresh_3D){
	warning("convolute_x_Qtf_yz test failed");
	return false;
      }else output("conolute_x_Qtf_yz test passed");


      return true;
    }

    bool test_projectors()const{



      bool result=true;

      output("Testing general project_out function");
      {
	real_function_6d f = make_xy(mo_ket_(0),mo_ket_(0));
	vecfuncT tmp; tmp.push_back(mo_ket_(0).function);
	CC_function_6d f_decomp(world,tmp,tmp);
	CC_function_6d f2(world,f);

	real_function_3d test1 = f.project_out(mo_bra_(0).function,0);
	real_function_3d test2 = f_decomp.project_out(mo_bra_(0).function,1);
	real_function_3d test3 = f2.project_out(mo_bra_(0).function,1);

	double diff12 = (test1-test2).norm2();
	double diff13 = (test1-test3).norm2();

	if(world.rank()==0){
	  std::cout << "||00_6d-00_decomposed||=" <<diff12 <<"\n";
	  std::cout << "||00_6d-00_6d||        =" <<diff13 <<"\n";
	}

	if(diff12>parameters.thresh_6D){
	  result =false;
	  warning("decomp_ proejct out failed");
	}else output("decomp_ project out passed");

	real_function_6d f00 = make_f_xy(mo_ket_(0),mo_ket_(0));
	CC_function_6d f00_decomp(world,&f12,mo_ket_(0).function,mo_ket_(0).function);
	CC_function_6d f002(world,f00);

	real_function_3d testa1 = f00.project_out(mo_bra_(0).function,0);
	real_function_3d testa2 = f00_decomp.project_out(mo_bra_(0).function,1);
	real_function_3d testa3 = f002.project_out(mo_bra_(0).function,1);

	double diffa12 = (testa1-testa2).norm2();
	double diffa13 = (testa1-testa3).norm2();

	if(world.rank()==0){
	  std::cout << "||f00_6d-f00_decomposed||=" <<diffa12 <<"\n";
	  std::cout << "||f00_6d-f00_6d||        =" <<diffa13 <<"\n";
	}

	if(diffa12>parameters.thresh_6D){
	  result =false;
	  warning("op_decomp_ proejct out failed");
	}else output("op_decomp_ project out passed");

      }



      output("Testing S2c with new 6d_function");


      Pairs<CC_function_6d> fmomo;
      Pairs<CC_function_6d> O1fmomo;
      Pairs<CC_function_6d> O2fmomo;
      Pairs<CC_function_6d> O1O2fmomo;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
	for(size_t j=i;j<mo_ket_.size();j++){
	  CC_function_6d tmp1(world,&f12,mo_ket_(i).function,mo_ket_(j).function);
	  fmomo.insert(i,j,tmp1);

	  std::pair<vecfuncT,vecfuncT> O1f00 = make_O1t_op_xy(mo_ket_,f12,mo_ket_(i),mo_ket_(j));
	  std::pair<vecfuncT,vecfuncT> O2f00 = make_O2t_op_xy(mo_ket_,f12,mo_ket_(i),mo_ket_(j));
	  std::pair<vecfuncT,vecfuncT> O1O2f00 = std::make_pair(P(O2f00.first,mo_ket_),O2f00.second);

	  CC_function_6d O1f(world,O1f00.first,O1f00.second);
	  O1fmomo.insert(i,j,O1f);
	  CC_function_6d O2f(world,O2f00.first,O2f00.second);
	  O2fmomo.insert(i,j,O2f);
	  CC_function_6d O1O2f(world,O1O2f00.first,O1O2f00.second);
	  O1O2fmomo.insert(i,j,O1O2f);


	}
      }

      vecfuncT old_reg = S2c_reg_part(CC_vecfunction(get_active_mo_ket(),HOLE,parameters.freeze),CC_vecfunction(get_active_mo_ket(),HOLE,parameters.freeze),mo_ket_);
      vecfuncT new_reg_1 = S2c_u_part(fmomo,mo_ket_);
      vecfuncT new_reg_2 = S2c_u_part(O1fmomo,mo_ket_);
      vecfuncT new_reg_3 = S2c_u_part(O2fmomo,mo_ket_);
      vecfuncT new_reg_4 = S2c_u_part(O1O2fmomo,mo_ket_);
      vecfuncT new_reg = sub(world,add(world,new_reg_1,new_reg_4),add(world,new_reg_2,new_reg_3));

      const double testx1 = norm2(world,old_reg);
      const double testx2 = norm2(world,new_reg);
      const double diffx = norm2(world,sub(world,old_reg,new_reg));
      if(world.rank()==0) std::cout << "||Vold||=" << testx1 << "\n";
      if(world.rank()==0) std::cout << "||Vnew||=" << testx2 << "\n";
      if(world.rank()==0) std::cout << "||diff||=" << diffx << "\n";
      if(diffx>parameters.thresh_3D){
	result = false;
	warning("Test for S2c_new failed!");
      }else output("Test for S2c passed");

      real_function_6d tmp_fxy=make_f_xy(mo_ket_(0),mo_ket_(0));
      apply_Q12(tmp_fxy);
      Pairs<CC_Pair> fxy_pairs;
      CC_Pair fxy_pair(tmp_fxy,0,0,GROUND_STATE);
      fxy_pairs.insert(0,0,fxy_pair);
      vecfuncT old_reg_u = S2c_u_part_old(fxy_pairs,mo_ket_);
      if(world.rank()==0) std::cout << "||Vold with Q12fxy||=" << norm2(world,old_reg_u) << "\n";




      output("Testing G_O1Q2_f12_xy function");
      const CC_function x=mo_ket_(mo_ket_.size()-1);
      const CC_function y=mo_ket_(mo_ket_.size()-1);
      const CC_vecfunction Oprojector = mo_ket_;
      const CC_vecfunction Qprojector = mo_ket_;

      const std::pair<vecfuncT,vecfuncT> O1fxy_decomposed = make_O1t_op_xy(Oprojector,f12,x,y);
      real_function_6d O1fxy_a = real_function_6d(world);
      for(size_t i=0;i<O1fxy_decomposed.first.size();i++){
	O1fxy_a += make_xy(O1fxy_decomposed.first[i],O1fxy_decomposed.second[i]);
      }

      real_function_6d O1Q2fxy_a = real_function_6d(world);
      vecfuncT Q2part = Q(O1fxy_decomposed.second,Qprojector);
      for(size_t i=0;i<O1fxy_decomposed.first.size();i++){
	O1Q2fxy_a += make_xy(O1fxy_decomposed.first[i],Q2part[i]);
      }

      real_function_6d O1O2fxy_a = real_function_6d(world);
      vecfuncT O2part = P(O1fxy_decomposed.second,Qprojector);
      for(size_t i=0;i<O1fxy_decomposed.first.size();i++){
	O1O2fxy_a += make_xy(O1fxy_decomposed.first[i],O2part[i]);
      }

      const real_function_6d fxy = make_f_xy(x,y);
      real_function_6d O1fxy = real_factory_6d(world);
      for(const auto & ktmp:Oprojector.functions){
	const real_function_3d tmp = fxy.project_out(mo_bra_(ktmp.first).function,0);
	O1fxy += make_xy(ktmp.second.function,tmp);
      }
      real_function_6d O1O2fxy = real_factory_6d(world);
      for(const auto & ktmp:Qprojector.functions){
	const real_function_3d tmp = O1fxy.project_out(mo_bra_(ktmp.first).function,1);
	O1O2fxy += make_xy(tmp,ktmp.second.function);
      }
      const real_function_6d O1Q2fxy = O1fxy - O1O2fxy;

      const double diff1 = (O1fxy-O1fxy_a).norm2();
      const double diff2 = (O1Q2fxy-O1Q2fxy_a).norm2();
      const double diff3 = (O1O2fxy - O1O2fxy_a).norm2();
      if(world.rank()==0){
	std::cout << "||O1fxy - O1fxy||=" << diff1 << "\n";
	std::cout << "||O1O2fxy - O1O2fxy||=" << diff3 << "\n";
	std::cout << "||O1Q2fxy - O1Q2fxy||=" << diff2 << "\n";
      }
      if(diff1>parameters.thresh_6D or diff3>parameters.thresh_6D){
	warning("Test failed");
	result=false;
      }output("Test passed");

      return result;
    }






  };



} /* namespace madness */

#endif /* CCOPERATORS_H_ */
