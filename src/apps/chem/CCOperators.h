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

  /// Small Helper Structure which applies either the Coulomb or f12 operator and can store intermediates <H|g|X> where the bra element is always a HOLE state and X can be HOLE, PARTICLE or RESPONSE
  /// The strucutre is initialized with the CC_Parameters to make sure all operations are carried out with the same parametrization in the operators
  struct CC_convolution_operator{
    CC_convolution_operator(World &world,const optype type, const CC_Parameters param):world(world),operator_type(type),op(init_op(type,param)){}

    real_function_3d operator()(const real_function_3d &f)const {return ((*op)(f)).truncate();}

    real_function_3d operator()(const CC_function &bra, const CC_function &ket, const bool use_im=true)const{
      real_function_3d result;
      if(not use_im){
	if(world.rank()==0) std::cout <<"Recalculating <" << bra.name()<<"|"<< assign_name(operator_type) <<"|"<<ket.name() <<">\n";
	result = ((*op)(bra.function*ket.function)).truncate();
      }
      else if(bra.type==HOLE and ket.type==HOLE and not imH.allpairs.empty()) result = imH(bra.i,ket.i);
      else if(bra.type==HOLE and ket.type==RESPONSE and not imR.allpairs.empty()) result = imR(bra.i,ket.i);
      else if(bra.type==HOLE and ket.type==PARTICLE and not imP.allpairs.empty()) result = imP(bra.i,ket.i);
      else if(bra.type==HOLE and ket.type==MIXED and (not imP.allpairs.empty() and not imH.allpairs.empty())) result = (imH(bra.i,ket.i)+imP(bra.i,ket.i));
      else{
	if(world.rank()==0) std::cout <<"No Intermediate found for <" << bra.name()<<"|"<<assign_name(operator_type) <<"|"<<ket.name() <<"> ... recalculate \n";
	result = ((*op)(bra.function*ket.function)).truncate();
      }
      return result;
    }

    real_function_6d operator()(const real_function_6d &u, const size_t particle)const{
      MADNESS_ASSERT(particle==1 or particle==2);
      MADNESS_ASSERT(operator_type == g12_);
      op->particle()=particle;
      return (*op)(u);
    }

    real_function_3d operator()(const CC_function &bra,const real_function_6d &u, const size_t particle)const{
      MADNESS_ASSERT(particle==1 or particle==2);
      MADNESS_ASSERT(operator_type == g12_);
      const real_function_6d tmp = multiply(copy(u),copy(bra.function),particle);
      op->particle()=particle;
      const real_function_6d g_tmp = (*op)(tmp);
      const real_function_3d result = g_tmp.dirac_convolution<3>();
      return result;
    }


    void update_elements(const CC_vecfunction &bra, const CC_vecfunction &ket){
      const  std::string operation_name = "<"+assign_name(bra.type)+"|"+name()+"|"+assign_name(ket.type)+">";
      if(world.rank()==0) std::cout << "updating operator elements: " << operation_name << " (" << bra.size() <<"x"<<ket.size() <<")"<< std::endl;
      if(bra.type != HOLE) error("Can not create intermediate of type "+operation_name+" , bra-element has to be of type HOLE");
      intermediateT xim;
      for(auto tmpk : bra.functions){
	const CC_function & k=tmpk.second;
	for(auto tmpl : ket.functions){
	  const CC_function& l=tmpl.second;
	  real_function_3d kl=(bra(k).function * l.function);
	  real_function_3d result=((*op)(kl)).truncate();
	  result.reconstruct(); // for sparse multiplication
	  xim.insert(k.i,l.i,result);
	}
      }
      if(ket.type==HOLE) imH=xim;
      else if(ket.type==PARTICLE) imP=xim;
      else if(ket.type==RESPONSE) imR=xim;
      else error("Can not create intermediate of type <"+assign_name(bra.type)+"|op|"+assign_name(ket.type)+">");
    }
    std::string name()const{return assign_name(operator_type);}

    void clear_intermediates(const functype &type){
      if(world.rank()==0) std::cout <<"Deleting all <HOLE|" << name() <<"|" << assign_name(type) << "> intermediates \n";
      switch(type){
	case HOLE : {imH.allpairs.clear(); break;}
	case PARTICLE:{imP.allpairs.clear(); break;}
	case RESPONSE:{imR.allpairs.clear(); break;}
	default: error("intermediates for " + assign_name(type) + " are not defined");
      }
    }
    void clear_all_intermediates(){
      clear_intermediates(HOLE);
      clear_intermediates(PARTICLE);
      clear_intermediates(RESPONSE);
    }
    size_t info()const{
      const size_t size_imH = size_of(imH);
      const size_t size_imP = size_of(imP);
      const size_t size_imR = size_of(imR);
      if(world.rank()==0){
	std::cout <<"Size of " << name() <<" intermediates:\n";
	std::cout <<std::setw(5)<<"("<<imH.allpairs.size() << ") x <H|"+name()+"H>=" << std::scientific << std::setprecision(1) << size_imH << " (Gbyte)\n";
	std::cout <<std::setw(5)<<"("<<imP.allpairs.size() << ") x <H|"+name()+"P>=" << std::scientific << std::setprecision(1) << size_imH << " (Gbyte)\n";
	std::cout <<std::setw(5)<<"("<<imR.allpairs.size() << ") x <H|"+name()+"R>=" << std::scientific << std::setprecision(1) << size_imH << " (Gbyte)\n";
      }
      return size_imH+size_imP + size_imR;
    }
    void sanity()const{
      print_intermediate(HOLE);
    }
    void print_intermediate(const functype type)const{
      if(type==HOLE)for(const auto& tmp:imH.allpairs)tmp.second.print_size("<H"+std::to_string(tmp.first.first)+"|"+assign_name(operator_type)+"|H"+std::to_string(tmp.first.second)+"> intermediate");
      else if(type==PARTICLE)for(const auto& tmp:imP.allpairs)tmp.second.print_size("<H"+std::to_string(tmp.first.first)+"|"+assign_name(operator_type)+"|P"+std::to_string(tmp.first.second)+"> intermediate");
      else if(type==RESPONSE)for(const auto& tmp:imR.allpairs)tmp.second.print_size("<H"+std::to_string(tmp.first.first)+"|"+assign_name(operator_type)+"|R"+std::to_string(tmp.first.second)+"> intermediate");
    }
  private:
    World &world;
    const optype operator_type;
    SeparatedConvolution<double,3>* init_op(const optype &type,const CC_Parameters &parameters)const{
      switch(type){
	case g12_ : {
	  if(world.rank()==0) std::cout << "Creating " << assign_name(type) <<" Operator with thresh=" << parameters.thresh_poisson <<" and lo=" << parameters.lo << std::endl;
	  return CoulombOperatorPtr(world, parameters.lo,parameters.thresh_poisson);
	}
	case f12_ : {
	  if(world.rank()==0) std::cout << "Creating " << assign_name(type) <<" Operator with thresh=" << parameters.thresh_poisson <<" and lo=" << parameters.lo << " and Gamma=" << parameters.gamma() << std::endl;
	  return SlaterF12OperatorPtr(world, parameters.gamma(),parameters.lo, parameters.thresh_poisson);
	}
	default : {
	  error("Unknown operatorype " + assign_name(type));
	  MADNESS_EXCEPTION("error",1);
	}
      }

    }
    const std::shared_ptr<real_convolution_3d> op;
    intermediateT imH;
    intermediateT imP;
    intermediateT imR;
    void error(const std::string &msg)const{
      if(world.rank()==0) std::cout <<"\n\n!!!!ERROR in CC_convolution_operator: " << msg <<"!!!!!\n\n"<< std::endl;
      MADNESS_EXCEPTION(msg.c_str(),1);
    }

  };


  //  /// Structure that holds the CC intermediates and is able to refresh them
  //  struct CC_Intermediates {
  //  public:
  //    CC_Intermediates(World&world, const CC_vecfunction &bra,
  //		     const CC_vecfunction &ket, const Nemo&nemo,
  //		     const CC_Parameters &param) :
  //		       world(world),
  //		       parameters(param),
  //		       mo_bra_(bra),
  //		       mo_ket_(ket),
  //		       poisson(std::shared_ptr < real_convolution_3d> (CoulombOperatorPtr(world, parameters.lo,parameters.thresh_poisson))),
  //		       f12op(std::shared_ptr < real_convolution_3d> (SlaterF12OperatorPtr(world, parameters.gamma(),parameters.lo, parameters.thresh_poisson))),
  //		       density_(make_density(bra, ket)),
  //		       exchange_intermediate_(make_exchange_intermediate(bra, ket)),
  //		       f12_exchange_intermediate_(make_f12_exchange_intermediate(bra, ket)),
  //		       hartree_potential_(make_hartree_potential(density_)) { }
  //
  //    // return the whole intermediate (to estimate size)
  //    const intermediateT& get_EX()const{return exchange_intermediate_;}
  //    const intermediateT& get_pEX()const{return perturbed_exchange_intermediate_;}
  //    const intermediateT& get_fEX()const{return f12_exchange_intermediate_;}
  //    const intermediateT& get_pfEX()const{return perturbed_f12_exchange_intermediate_;}
  //    const intermediateT& get_rEX()const{return response_exchange_intermediate_;}
  //    const intermediateT& get_rfEX()const{return response_f12_exchange_intermediate_;}
  //
  //    /// Get the intermediates
  //    const real_function_3d get_density() const {
  //      return density_;
  //    }
  //    const real_function_3d get_perturbed_density() const {
  //      return perturbed_density_;
  //    }
  //    const real_function_3d get_hartree_potential() const {
  //      return hartree_potential_;
  //    }
  //    const real_function_3d get_perturbed_hartree_potential() const {
  //      return perturbed_hartree_potential_;
  //    }
  //    /// returns <k|g|l>
  //    const real_function_3d get_EX(const size_t &k, const size_t &l) const {
  //      return exchange_intermediate_(k, l);
  //    }
  //    const real_function_3d get_EX(const CC_function &k,
  //				  const CC_function &l) const {
  //      return exchange_intermediate_(k.i, l.i);
  //    }
  //    const real_function_3d get_fEX(const size_t &k, const size_t &l) const {
  //      return f12_exchange_intermediate_(k, l);
  //    }
  //    const real_function_3d get_fEX(const CC_function &k,
  //				   const CC_function &l) const {
  //      return f12_exchange_intermediate_(k.i, l.i);
  //    }
  //
  //    /// returns <k|g|\tau_l>
  //    const real_function_3d get_pEX(const CC_function &k,
  //				   const CC_function &l) const {
  //      return perturbed_exchange_intermediate_(k.i, l.i);
  //    }
  //    const real_function_3d get_pEX(const size_t &k, const size_t &l) const {
  //      return perturbed_exchange_intermediate_(k, l);
  //    }
  //    /// returns <k|f|\tau_l>
  //    const real_function_3d get_pfEX(const CC_function &k,
  //				    const CC_function &l) const {
  //      return perturbed_f12_exchange_intermediate_(k.i, l.i);
  //    }
  //    const real_function_3d get_pfEX(const size_t &k, const size_t &l) const {
  //      return perturbed_f12_exchange_intermediate_(k, l);
  //    }
  //    const real_function_3d get_rEX(const CC_function &k,
  //				   const CC_function &l) const {
  //      return response_exchange_intermediate_(k.i, l.i);
  //    }
  //    const real_function_3d get_rEX(const size_t &k, const size_t &l) const {
  //      return response_exchange_intermediate_(k, l);
  //    }
  //    const real_function_3d get_rfEX(const CC_function &k,
  //				    const CC_function &l) const {
  //      return response_f12_exchange_intermediate_(k.i, l.i);
  //    }
  //    const real_function_3d get_rfEX(const size_t &k, const size_t &l) const {
  //      return response_f12_exchange_intermediate_(k, l);
  //    }
  //
  //    /// refresh the intermediates that depend on the \tau functions
  //    void update(const CC_vecfunction &tau) {
  //      if (world.rank() == 0)
  //	std::cout << "Update Intermediates:\n";
  //      perturbed_density_ = make_density(mo_bra_, tau);
  //      perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
  //      perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,tau);
  //      perturbed_f12_exchange_intermediate_ = make_f12_exchange_intermediate(mo_bra_, tau);
  //    }
  //
  //    void update_response(const CC_vecfunction &x){
  //      std::cout << "Update Response Intermediates:\n";
  //      response_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,x);
  //      response_f12_exchange_intermediate_ = make_f12_exchange_intermediate(mo_bra_,x);
  //    }
  //
  //    /// make a density from two input functions
  //    /// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
  //    /// @param[in] vecfuncT_bra
  //    /// @param[in] vecfuncT_ket
  //    /// @param[out] \sum_i bra_i * ket_i
  //    //real_function_3d make_density(const vecfuncT &bra,const vecfuncT &ket) const;
  //    real_function_3d make_density(const CC_vecfunction &bra,
  //				  const CC_vecfunction &ket) const;
  //    /// Poisson operator
  //    std::shared_ptr<real_convolution_3d> get_poisson() const {
  //      return poisson;
  //    }
  //
  //  private:
  //    World &world;
  //    const CC_Parameters &parameters;
  //    const CC_vecfunction &mo_bra_;
  //    const CC_vecfunction &mo_ket_;
  //    const std::shared_ptr<real_convolution_3d> poisson;
  //    const std::shared_ptr<real_convolution_3d> f12op;
  //    /// const intermediates
  //    const real_function_3d density_;
  //    /// Exchange intermediate: \f$EX(i,j) = <i|g|j>\f$
  //    intermediateT exchange_intermediate_;
  //    /// The f12 exchange intermediate \f$fEX(i,j) = <i|f12|j>\f$
  //    intermediateT f12_exchange_intermediate_;
  //    /// Hartree_Potential  \f$ = J = \sum_k <k|g|k> = \f$ Poisson(density)
  //    const real_function_3d hartree_potential_;
  //    /// intermediates that need to be recalculated before every iteration
  //    /// Perturbed Density \f$= \sum_k |k><\tau_k| \f$
  //    real_function_3d perturbed_density_;
  //    /// Perturbed Hartree Potential PJ \f$ = \sum_k <k|g|\tau_k> = \f$ Poisson(perturbed_density)
  //    real_function_3d perturbed_hartree_potential_;
  //    /// Perturbed Exchange Intermediate: \f$ PEX(i,j) = <i|g|\tau_j> \f$
  //    intermediateT perturbed_exchange_intermediate_;
  //    /// Perturbed f12-exchange-intermediate: \f$ pfEX(i,j) = <i|f12|tau_j> \f$
  //    intermediateT perturbed_f12_exchange_intermediate_;
  //    /// response exchange-intermediate: \f$ rEX(i,j) = <i|g12|x_j> \f$
  //    intermediateT response_exchange_intermediate_;
  //    /// response f12-exchange-intermediate: \f$ rEX(i,j) = <i|f12|x_j> \f$
  //    intermediateT response_f12_exchange_intermediate_;
  //
  //    void error(const std::string &msg) const {
  //      std::cout << "\n\n\nERROR IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
  //      MADNESS_EXCEPTION(
  //	  "\n\n!!!!ERROR IN CC_INTERMEDIATES!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
  //	  1);
  //    }
  //    void warning(const std::string &msg) const {
  //      std::cout << "\n\n\nWARNING IN CC_INTERMEDIATES:\n" << msg
  //	  << "\n\n\n!!!";
  //    }
  //  public:
  //    /// Make the exchange intermediate: EX[j][i] \f$ <bra[i](r2)|1/r12|ket[j](r2)> \f$
  //    intermediateT make_exchange_intermediate(const CC_vecfunction &bra,
  //					     const CC_vecfunction &ket) const;
  //    intermediateT make_f12_exchange_intermediate(const CC_vecfunction &bra,
  //						 const CC_vecfunction &ket) const;
  //    /// Calculates the hartree potential Poisson(density)
  //    /// @param[in] density A 3d function on which the poisson operator is applied (can be the occupied density and the perturbed density)
  //    /// @return poisson(density) \f$ = \int 1/r12 density(r2) dr2 \f$
  //    real_function_3d make_hartree_potential(
  //	const real_function_3d &density) const {
  //      real_function_3d hartree = (*poisson)(density);
  //      hartree.truncate();
  //      return hartree;
  //    }
  //  };

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
      //g12.print_intermediate(PARTICLE);
      //f12.print_intermediate(PARTICLE);
      //      intermediates_.update(singles);
      update.info();
    }

    void update_response_intermediates(const CC_vecfunction &singles) {
      CC_Timer update(world, "Update Intermediates");
      g12.update_elements(mo_bra_,singles);
      f12.update_elements(mo_bra_,singles);
      //      intermediates_.update(singles);
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


    vecfuncT get_CCS_response_potential(const CC_vecfunction &x){
      if(x.type!=RESPONSE) error("get_CCS_response_potential: Wrong type of input singles");
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
      const double bsh_eps = get_epsilon(moi.i,moj.i)+omega;
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 *bsh_eps),parameters.lo,parameters.thresh_bsh_6D);

      real_function_6d result = real_factory_6d(world);
      for(const auto& ktmp:x.functions){
	CC_function mokbra = mo_bra_(ktmp.first);
	CC_function xk = ktmp.second;

	const real_function_3d p1 = apply_Q(g12(mokbra,moj)*moi.function);
	const real_function_3d p2 = apply_Q(g12(mokbra,moi)*moj.function);

	result -= -2.0*(G(p1,xk.function));
	result -= -2.0*(G(xk.function,p2));

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

    /// reconstructs the full pair function from the regularized pair functions
    /// used to compute norms of the doubles to compare with LCAO codes
    /// used to debug the singles potential
    /// @param[in] u the regularized function
    /// @return Equation: \f$ \tau = u + Q12f12(|ij> + |taui,j> + |i,tauj> + |taui,tauj>) = u + Q12f12|titj> \f$ with \f$ ti = taui + i \f$
    real_function_6d make_full_pair_function(const CC_Pair &u,
					     const CC_function &taui, const CC_function &tauj) const {
      const size_t i = u.i;
      const size_t j = u.j;
      MADNESS_ASSERT(i == taui.i);
      MADNESS_ASSERT(j == tauj.i);
      real_function_3d ti = mo_ket_(i).function + taui.function;
      real_function_3d tj = mo_ket_(j).function + tauj.function;
      real_function_6d Q12f12titj = make_f_xy(ti, tj);
      apply_Q12(Q12f12titj);
      real_function_6d result = u.function + Q12f12titj;
      return result;
    }
    CC_Pair make_full_pair(const CC_Pair &u, const CC_function &taui,
			   const CC_function &tauj) const {
      real_function_6d full_pair_function = make_full_pair_function(u, taui,
								    tauj);
      CC_Pair result(full_pair_function, u.i, u.j,u.type);
      return result;
    }
    Pairs<CC_Pair> make_full_pairs(const Pairs<CC_Pair> &pairs,
				   const CC_vecfunction &singles) const {
      Pairs<CC_Pair> result;
      for (auto utmp : pairs.allpairs) {
	CC_Pair u = utmp.second;
	CC_Pair ufull = make_full_pair(u, singles(u.i), singles(u.j));
	result.insert(ufull.i, ufull.j, ufull);
      }
      return result;
    }

    Pairs<CC_Pair> make_reg_residues(const CC_vecfunction &singles) const {
      CC_Timer time(world,"Making Regularization-Tails of Pair-Functions");
      Pairs<CC_Pair> result;
      CC_vecfunction t = make_t_intermediate(singles);
      for (size_t i=parameters.freeze;i<mo_ket_.size();i++) {
	for(size_t j=i;j<mo_ket_.size();j++){
	  real_function_6d Qftitj = make_f_xy(t(i),t(j));
	  apply_Q12(Qftitj);
	  Qftitj.print_size("Q12|t"+stringify(i)+"t"+stringify(j)+">");
	  CC_Pair pair_tmp(Qftitj,i,j,GROUND_STATE);
	  result.insert(pair_tmp.i, pair_tmp.j, pair_tmp);
	}
      }
      time.info();
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
      if(V.empty()) pot = get_CCS_response_potential(x);
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

      const vecfuncT xbra = mul(world,nemo.nuclear_correlation->square(),x.get_vecfunction());

      // S2 Part;
      const vecfuncT amo_tmp = get_active_mo_ket();
      const CC_vecfunction amo(amo_tmp,HOLE,parameters.freeze);
      const double s2bu = inner(world,xbra,S2b_u_part(chi,x)).sum();
      const double s2cu = inner(world,xbra,S2c_u_part(chi,x)).sum();
      const double s2br = inner(world,xbra,add(world,S2b_reg_part(x,amo),S2b_reg_part(amo,x))).sum();
      const double s2cr = inner(world,xbra,add(world,S2c_reg_part(x,amo),S2c_reg_part(amo,x))).sum();

      // S4 Part;
      const double s4au = inner(world,xbra,S4a_u_part(u,x)).sum();
      const double s4ar = inner(world,xbra,S4a_reg_part(amo,amo,x)).sum();
      const double s4bu = inner(world,xbra,S4b_u_part(u,x)).sum();
      const double s4br = inner(world,xbra,S4b_reg_part(amo,amo,x)).sum();
      const double s4cu = inner(world,xbra,S4c_u_part(u,x)).sum();
      const double s4cr = inner(world,xbra,S4c_reg_part(amo,amo,x)).sum();

      double result = 0.0;
      {
	if(world.rank()==0)std::cout << " CIS(D) Energy Correction:\n";
	if(world.rank()==0)std::cout <<"s2bu-part="<< s2bu << "\n";result+=s2bu;
	if(world.rank()==0)std::cout <<"s2cu-part="<< s2cu << "\n";result+=s2cu;
	if(world.rank()==0)std::cout <<"s2br-part="<< s2br << "\n";result+=s2br;
	if(world.rank()==0)std::cout <<"s2cr-part="<< s2cr << "\n";result+=s2cr;
	if(world.rank()==0)std::cout <<"s4au-part="<< s4au << "\n";result+=s4au;
	if(world.rank()==0)std::cout <<"s4ar-part="<< s4ar << "\n";result+=s4ar;
	if(world.rank()==0)std::cout <<"s4bu-part="<< s4bu << "\n";result+=s4bu;
	if(world.rank()==0)std::cout <<"s4br-part="<< s4br << "\n";result+=s4br;
	if(world.rank()==0)std::cout <<"s4cu-part="<< s4cu << "\n";result+=s4cu;
	if(world.rank()==0)std::cout <<"s4cr-part="<< s4cr << "\n";result+=s4cr;
	if(world.rank()==0)std::cout <<"All together = " << result << "\n";
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

    /// CCSD/CC2 singles potential parts
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
	case pot_S2b_u_:
	  result = S2b_u_part(u, singles);
	  break;
	case pot_S2c_u_:
	  result = S2c_u_part(u, singles);
	  break;
	case pot_S4a_u_:
	  result = S4a_u_part(u, singles);
	  break;
	case pot_S4b_u_:
	  result = S4b_u_part(u, singles);
	  break;
	case pot_S4c_u_:
	  result = S4c_u_part(u, singles);
	  break;
	case pot_S2b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  result = S2b_reg_part(t,t);
	  break;
	}
	case pot_S2c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  result = S2c_reg_part(t,t);
	  break;
	}
	case pot_S4a_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  result = S4a_reg_part(t,t,singles);
	  break;
	}
	case pot_S4b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  result = S4b_reg_part(t,t,singles);
	  break;
	}
	case pot_S4c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(singles);
	  result = S4c_reg_part(t,t,singles);
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

    vecfuncT response_potential_singles(
	const CC_vecfunction & gs_singles,
	const Pairs<CC_Pair>& gs_u,
	const CC_vecfunction & response_singles,
	const Pairs<CC_Pair>& response_u,
	const potentialtype_s &name) const {

      const std::string full_name = assign_name(CC2_response_)+":"+assign_name(name);
      CC_Timer timer(world, full_name);
      CC_data data(full_name);
      if(response_singles.type!=RESPONSE) error("Error in response_potential_singles " + assign_name(name)+" response_singles have wrong type");
      if(gs_singles.type!=PARTICLE) error("Error in response_potential_singles " + assign_name(name)+" gs_singles have wrong type");
      vecfuncT result;

      switch (name) {
	case pot_F3D_:
	  result = fock_residue_closed_shell(response_singles);
	  break;
	case pot_ccs_:
	  result = ccs_response_potential(gs_singles,response_singles);
	  break;
	case pot_cis_:
	  result = cis_response_potential(response_singles);
	  break;
	case pot_S2b_u_:
	  result = S2b_u_part(response_u, response_singles);
	  break;
	case pot_S2c_u_:
	  result = S2c_u_part(response_u, response_singles);
	  break;
	case pot_S4a_u_:
	{
	  const vecfuncT part1  = S4a_u_part(response_u, gs_singles);
	  const vecfuncT part2  = S4a_u_part(gs_u, response_singles);
	  result = add(world,part1,part2);
	  break;
	}
	case pot_S4b_u_:
	{
	  const vecfuncT part1  = S4b_u_part(response_u, gs_singles);
	  const vecfuncT part2  = S4b_u_part(gs_u, response_singles);
	  result = add(world,part1,part2);
	  break;
	}
	case pot_S4c_u_:
	{
	  const vecfuncT part1  = S4c_u_part(response_u, gs_singles);
	  const vecfuncT part2  = S4c_u_part(gs_u, response_singles);
	  result = add(world,part1,part2);
	  break;
	}
	case pot_S2b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(gs_singles);
	  const vecfuncT tx_part = S2b_reg_part(t,response_singles);
	  const vecfuncT xt_part = S2b_reg_part(response_singles,t);
	  result = add(world,tx_part,xt_part);
	  break;
	}
	case pot_S2c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(gs_singles);
	  const vecfuncT tx_part = S2c_reg_part(t,response_singles);
	  const vecfuncT xt_part = S2c_reg_part(response_singles,t);
	  result = add(world,tx_part,xt_part);
	  break;
	}
	case pot_S4a_r_:
	{
	  CC_vecfunction t = make_t_intermediate(gs_singles);
	  // part with gs_doubles and response_singles
	  const vecfuncT part1 = S4a_reg_part(t,t,response_singles);
	  // part with response_doubles and gs_singles
	  const vecfuncT xt_part = S4a_reg_part(response_singles,t,gs_singles);
	  const vecfuncT tx_part = S4a_reg_part(t,response_singles,gs_singles);
	  const vecfuncT part2 = add(world,xt_part,tx_part);
	  result = add(world,part1,part2);
	  break;
	}
	case pot_S4b_r_:
	{
	  CC_vecfunction t = make_t_intermediate(gs_singles);
	  // part with gs_doubles and response_singles
	  const vecfuncT part1 = S4b_reg_part(t,t,response_singles);
	  // part with response_doubles and gs_singles
	  const vecfuncT xt_part = S4b_reg_part(response_singles,t,gs_singles);
	  const vecfuncT tx_part = S4b_reg_part(t,response_singles,gs_singles);
	  const vecfuncT part2 = add(world,xt_part,tx_part);
	  result = add(world,part1,part2);
	  break;
	}
	case pot_S4c_r_:
	{
	  CC_vecfunction t = make_t_intermediate(gs_singles);
	  // part with gs_doubles and response_singles
	  const vecfuncT part1 = S4c_reg_part(t,t,response_singles);
	  // part with response_doubles and gs_singles
	  const vecfuncT xt_part = S4c_reg_part(response_singles,t,gs_singles);
	  const vecfuncT tx_part = S4c_reg_part(t,response_singles,gs_singles);
	  const vecfuncT part2 = add(world,xt_part,tx_part);
	  result = add(world,part1,part2);
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


    vecfuncT fock_residue_closed_shell(const CC_vecfunction &tau) const;

    /// brilloin terms of ccs
    vecfuncT S1(const CC_vecfunction &tau) const;
    vecfuncT S5a(const CC_vecfunction &tau) const;

    /// The CCS Potential without Brillouin terms and Fock residue
    vecfuncT ccs_potential(const CC_vecfunction &tau) const;
    vecfuncT ccs_potential_new(const CC_vecfunction &tau)const;
    vecfuncT ccs_unprojected(const CC_vecfunction &ti, const CC_vecfunction &tk)const;
    vecfuncT ccs_response_potential(const CC_vecfunction &singles, const CC_vecfunction &response)const;

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
    vecfuncT S2b_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    // result: -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1)
    // singles are not needed explicitly but to determine if it is response or ground state
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1) \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S2c_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    /// The Part of the CC2 singles potential which depends on singles and doubles (S4a, S4b, S4c)

    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -|singles_l><l|S2b_i> = -|singles_l>(2.0*<lk|g|uik>-<kl|g|uik> \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4a_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    // result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ -( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4b_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[out] \f$ ( 4<l|kgtauk|uil>_2 - 2<l|kgtauk|uil>_1 - 2<k|lgtauk|uil>_2 + <k|lgtauk|uil>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    vecfuncT S4c_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const;

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[out] \f$ 2<k|gQf|reg1_i,reg2_k>_2 - <k|gQf|reg1_i,reg2_k>_1 \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S2b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2) const;

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[out] \f$ 2<l*kgi|Qf|reg1_k,reg2_l>_2 - <l*kgi|gQf|reg1_k,reg2_l>_1, kgi = <k|g|i> \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S2c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2) const;

    /// result: -\sum_{kl}( 2 <l|kgtaui|Qftktl> - <l|kgtaui|Qftltk>
    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ -|singles_l>(2.0*<lk|gQf|reg1_i,reg2_k>-<kl|gQf|reg1_i,reg2_k> \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4a_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles) const;

    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ -( <l*kgtaui|Qf|reg1_k,reg2_l>_2 - <l*kgtaui|Qf|reg1_k,reg2_l>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4b_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles) const;

    /// result: 4<l|kgtauk|Qftitl> - 2<l|kgtauk|Qftlti> - 2<k|lgtauk|Qftitl> + <k|lgtauk|Qftlti>
    /// Regularization-Part of electron Pair: Qf|reg1,reg2>
    /// @param[in] reg1: 3D function for particle 1 in regularization-part of electron-pair
    /// @param[in] reg2: 3D function for particle 2 in regularization-part of electron-pair
    /// @param[in] singles: the singles functions (particle or response) which interact with the Regularization-Part
    /// @param[out] \f$ ( 4<l*kgtauk|Qf|reg1_i,reg2_l>_2 - 2<l*kgtauk|Qf|reg1_i,reg2_l>_1 - 2<k*lgtauk|Qf|reg1_i,reg2_l>_2 + <k*lgtauk|Qf|reg1_i,reg2_l>_1 ) \f$
    /// Q-Projector is not applied (to the result), sign is correct
    vecfuncT S4c_reg_part(const CC_vecfunction &reg1, const CC_vecfunction &reg2,const CC_vecfunction &singles) const;

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

    /// param[in]	f	a function of 2 particles f(1,2)
    /// return	the input function with particles swapped g(1,2) = f(2,1)
    real_function_6d swap_particles(const real_function_6d& f) const;

    // Calculate the CC2 energy equation which is
    // \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
    // with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
    double get_CC2_correlation_energy() const;
    double compute_cispd_energy(const Pairs<CC_Pair> &u, const Pairs<CC_Pair> mp2_doubles, const CC_vecfunction x);
    double compute_cispd_energy_constant_part(const Pairs<CC_Pair> &u, const CC_vecfunction x)const;
    double compute_cc2_pair_energy(const CC_Pair &u, const CC_function &taui,
				   const CC_function &tauj) const;
    /// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
    // the bra elements are always the R2orbitals
    // the ket elements can be \tau_i , or orbitals dependet n the type given
    double make_ij_gQf_ij(const size_t &i, const size_t &j, CC_Pair &u) const;
    double make_ijgQfxy(const size_t &i, const size_t &j, const CC_function &x,
			const CC_function &y) const;
    double make_ijgQfxy(const  CC_function &i, const  CC_function &j, const CC_function &x,
			const CC_function &y) const;
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




  };



} /* namespace madness */

#endif /* CCOPERATORS_H_ */
