/*
 * CCStructures.cc
 *
 *  Created on: 6 Jan 2017
 *      Author: kottmanj
 */

#include"CCStructures.h"

namespace madness{

    real_function_3d CC_convolution_operator::operator()(const CC_function &bra, const CC_function &ket, const bool use_im)const{
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
	if(world.rank()==0 and ket.type==HOLE) std::cout <<"No Intermediate found for <" << bra.name()<<"|"<<assign_name(operator_type) <<"|"<<ket.name() <<"> ... recalculate \n";
	result = ((*op)(bra.function*ket.function)).truncate();
      }
      return result;
    }

    real_function_6d CC_convolution_operator::operator()(const real_function_6d &u, const size_t particle)const{
      MADNESS_ASSERT(particle==1 or particle==2);
      MADNESS_ASSERT(operator_type == g12_);
      op->particle()=particle;
      return (*op)(u);
    }

    real_function_3d CC_convolution_operator::operator()(const CC_function &bra,const real_function_6d &u, const size_t particle)const{
      MADNESS_ASSERT(particle==1 or particle==2);
      MADNESS_ASSERT(operator_type == g12_);
      const real_function_6d tmp = multiply(copy(u),copy(bra.function),particle);
      op->particle()=particle;
      const real_function_6d g_tmp = (*op)(tmp);
      const real_function_3d result = g_tmp.dirac_convolution<3>();
      return result;
    }


    void CC_convolution_operator::update_elements(const CC_vecfunction &bra, const CC_vecfunction &ket){
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


    void CC_convolution_operator::clear_intermediates(const functype &type){
      if(world.rank()==0) std::cout <<"Deleting all <HOLE|" << name() <<"|" << assign_name(type) << "> intermediates \n";
      switch(type){
	case HOLE : {imH.allpairs.clear(); break;}
	case PARTICLE:{imP.allpairs.clear(); break;}
	case RESPONSE:{imR.allpairs.clear(); break;}
	default: error("intermediates for " + assign_name(type) + " are not defined");
      }
    }

    size_t CC_convolution_operator::info()const{
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

    SeparatedConvolution<double,3>* CC_convolution_operator::init_op(const optype &type,const CC_Parameters &parameters)const{
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

  //  enum pairtype {GROUND_STATE,EXCITED_STATE};
  //  enum potentialtype_s {pot_F3D_,pot_s3a_,pot_s3b_, pot_s3c_,pot_s5a_,pot_s5b_,pot_s5c_,pot_s6_, pot_s2b_, pot_s2c_, pot_s4a_, pot_s4b_, pot_s4c_, pot_S2b_u_, pot_S2c_u_, pot_S4a_u_, pot_S4b_u_, pot_S4c_u_,pot_S2b_r_, pot_S2c_r_, pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_ccs_,pot_cis_,pot_singles_};
  //  enum potentialtype_d {pot_F6D_, pot_cc2_coulomb_,pot_cc2_residue_};
  //  // The pair function is:  \tau = u + Qf(|titj>), FULL means that \tau is calculated in 6D form, DECOMPOSED means that u is used in 6D and the rest is tried to solve in 3D whenever possible
  //  enum pair_function_form{DECOMPOSED, FULL};
  std::string
  assign_name(const optype& input) {
    switch(input){
      case g12_:
	return "g12";
      case f12_:
	return "f12";
    }
    MADNESS_EXCEPTION("assign_name:optype, should not end up here",1);
    return "unknown operatortype";
  }

  calctype
  assign_calctype(const std::string name) {
    if(name == "mp2") return MP2_;
    else if(name == "cc2") return CC2_;
    else if(name == "lrcc2" or name == "cc2_response") return LRCC2_;
    else if(name == "cispd") return CISpD_;
    else if(name == "cis" or name == "ccs" or name == "ccs_response" or name == "lrccs") return LRCCS_;
    else if(name == "experimental") return experimental_;
    else if(name == "adc2" or name == "adc(2)") return ADC2_;
    else if(name == "tdhf") return TDHF_;
    else{
      std::string msg="CALCULATION OF TYPE: " + name + " IS NOT KNOWN!!!!";
      MADNESS_EXCEPTION(msg.c_str(),1);
    }
  }

  std::string
  assign_name(const calctype& inp) {
    switch(inp){
      case CC2_:
	return "CC2";
      case MP2_:
	return "MP2";
      case LRCC2_:
	return "LRCC2";
      case CISpD_:
	return "CISpD";
      case LRCCS_:
	return "LRCCS";
      case ADC2_:
	return "ADC2";
      case TDHF_:
	return "TDHF";
      case experimental_:
	return "experimental";
    }
    return "unknown";
  }

  std::string
  assign_name(const functype& inp) {
    switch(inp){
      case HOLE:
	return "Hole";
      case PARTICLE:
	return "Particle";
      case MIXED:
	return "Mixed";
      case RESPONSE:
	return "Response";
      case UNDEFINED:
	return "Undefined";
    }
    return "???";
  }

  void
  CC_Parameters::information(World& world) const {
    if(world.rank() == 0){
      //			std::cout << "Defaults for 6D and 3D Functions:\n";
      //			FunctionDefaults<3>::print();
      //                      FunctionDefaults<6>::print();
      std::cout << "Demanded Calculation is " << assign_name(calculation) << std::endl;
      if(calculation != LRCCS_ && calculation != TDHF_){
	std::cout << "The Ansatz for the Pair functions |tau_ij> is: ";
	if(QtAnsatz) std::cout << "(Qt)f12|titj> and response: (Qt)f12(|tixj> + |xitj>) - (OxQt + QtOx)f12|titj>";
	else std::cout << "Qf12|titj> and response: Qf12(|xitj> + |tixj>)";

	std::cout << "Gamma of correlation factor is " << corrfac_gamma << std::endl;
      }
      if(test) std::cout << "\n\n\t\t\t!Test Mode is on!\n\n" << std::endl;

      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      std::cout << "\nMRA-Parameters:\n";
      std::cout << "lo                         :" << lo << std::endl;
      std::cout << "dmin                       :" << dmin << std::endl;
      std::cout << "k (3D)                     :" << FunctionDefaults<3>::get_k() << std::endl;
      std::cout << "k (6D)                     :" << FunctionDefaults<6>::get_k() << std::endl;
      std::cout << "thresh_3D demanded         :" << thresh_3D << std::endl;
      std::cout << "thresh_3D set              :" << FunctionDefaults<3>::get_thresh() << std::endl;
      std::cout << "thresh_6D demanded         :" << thresh_6D << std::endl;
      std::cout << "thresh_6D set              :" << FunctionDefaults<6>::get_thresh() << std::endl;
      std::cout << "tight_thresh_6D            :" << tight_thresh_6D << std::endl;
      std::cout << "tight_thresh_3D            :" << tight_thresh_3D << std::endl;
      std::cout << "thresh_bsh_3D              :" << thresh_bsh_3D << std::endl;
      std::cout << "thresh_bsh_6D              :" << thresh_bsh_6D << std::endl;
      std::cout << "thresh_poisson             :" << thresh_poisson << std::endl;
      std::cout << "thresh_f12                 :" << thresh_f12 << std::endl;
      std::cout << "thresh_Ue                  :" << thresh_Ue << std::endl;
      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      std::cout << "\nAdvanced-MRA-Parameters:\n";
      std::cout << "truncation mode 3D         :" << FunctionDefaults<3>::get_truncate_mode() << std::endl;
      std::cout << "truncation mode 6D         :" << FunctionDefaults<6>::get_truncate_mode() << std::endl;
      std::cout << "tensor type                :" << FunctionDefaults<6>::get_tensor_type() << std::endl;
      std::cout << "facReduce                  :" << GenTensor<double>::fac_reduce() << std::endl;
      std::cout << "max. displacement          :" << Displacements<6>::bmax_default() << std::endl;
      std::cout << "apply randomize            :" << FunctionDefaults<6>::get_apply_randomize() << std::endl;
      std::cout << "Cell min width (3D, 6D)    :" << FunctionDefaults<6>::get_cell_min_width() << ", " << FunctionDefaults<3>::get_cell_min_width() << std::endl;
      std::cout << "Autorefine (3D, 6D)        :" << FunctionDefaults<6>::get_autorefine() << ", " << FunctionDefaults<3>::get_autorefine() << std::endl;
      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      std::cout << "\nCC-Parameters:\n";
      std::cout << "freeze      :" << freeze << std::endl;
      std::cout << "restart     :" << restart << std::endl;
      std::cout << "econv       :" << econv << std::endl;
      std::cout << "econv_pairs :" << econv_pairs << std::endl;
      std::cout << "dconv_3D    :" << dconv_3D << std::endl;
      std::cout << "dconv_6D    :" << dconv_6D << std::endl;
      std::cout << "freeze      :" << freeze << std::endl;
      std::cout << "iter_max    :" << iter_max << std::endl;
      std::cout << "iter_max_3D :" << iter_max_3D << std::endl;
      std::cout << "iter_max_6D :" << iter_max_6D << std::endl;
      std::cout << "debug mode  :" << debug << std::endl;
      std::cout << "Kain        :";
      if(kain) std::cout << " on , Subspace size " << kain_subspace << "\n";
      else std::cout << "off\n";

      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      //std::cout << std::setw(20) << std::setfill(' ') << "Cell widths (3D) :" << FunctionDefaults<3>::get_cell_width()  <<std::endl;
      //std::cout << std::setw(20) << std::setfill(' ') << "Cell widths (6D) :" << FunctionDefaults<6>::get_cell_width()  <<std::endl;
      if(restart) std::cout << "restart is on";

      if(no_compute) std::cout << "no_compute for all";

      if(no_compute_gs) std::cout << "no_compute for ground-state";

      if(no_compute_response) std::cout << "no_compute for excited-state";

      std::cout << "Excitations to optimize are:\n";
      if(excitations_.empty()) std::cout << "All" << "\n";
      else std::cout << excitations_ << "\n\n\n";
    }
  }

  void
  CC_Parameters::print_tda_parameters(World& world) const {
    if(world.rank() == 0){
      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      std::cout << "TDA PARAMETERS:\n";
      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
      std::cout << std::scientific << std::setprecision(2);
      std::cout << "tda_guess_orbitals       :" << tda_guess_orbitals << std::endl;
      //std::cout << "tda_guess_mode           :" << tda_guess_mode              << std::endl;
      std::cout << "tda_excitations          :" << tda_excitations << std::endl;
      std::cout << "tda_guess_excitations    :" << tda_guess_excitations << std::endl;
      //std::cout << "tda_iterating_excitations:" << tda_iterating_excitations   << std::endl;
      std::cout << "tda_guess                :" << tda_guess << std::endl;
      std::cout << "tda_energy_guess_factor  :" << tda_energy_guess_factor << std::endl;
      std::cout << "tda_dconv_guess          :" << tda_dconv_guess << std::endl;
      std::cout << "tda_dconv                :" << tda_dconv << std::endl;
      //std::cout << "tda_dconv_hard           :" << tda_dconv_hard              << std::endl;
      std::cout << "tda_econv_guess          :" << tda_econv_guess << std::endl;
      std::cout << "tda_econv                :" << tda_econv << std::endl;
      //std::cout << "tda_econv_hard           :" << tda_econv_hard              << std::endl;
      std::cout << "tda_store_potential      :" << tda_store_potential << std::endl;
      std::cout << "tda_iter_max             :" << tda_iter_max << std::endl;
      std::cout << "tda_iter_guess           :" << tda_iter_guess << std::endl;
      std::cout << "tda_homo_guess           :" << tda_homo_guess << std::endl;
      std::cout << std::setfill('-') << std::setw(35) << std::setfill('-') << "\n";
      std::cout << std::setfill(' ');
    }
  }

  void
  CC_Parameters::sanity_check(World& world) const {
    size_t warnings=0;
    if(FunctionDefaults<3>::get_thresh() > 0.01 * FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"3D Thresh is too low, should be 0.01*6D_thresh");

    if(FunctionDefaults<3>::get_thresh() > 0.1 * FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"3D Thresh is way too low, should be 0.01*6D_thresh");

    if(FunctionDefaults<3>::get_cell_min_width() != FunctionDefaults<6>::get_cell_min_width()) warnings+=warning(world,"3D and 6D Cell sizes differ");

    if(FunctionDefaults<3>::get_k() != FunctionDefaults<6>::get_k()) warnings+=warning(world,"k-values of 3D and 6D differ ");

    if(FunctionDefaults<3>::get_truncate_mode() != 3) warnings+=warning(world,"3D Truncate mode is not 3");

    if(FunctionDefaults<6>::get_truncate_mode() != 3) warnings+=warning(world,"6D Truncate mode is not 3");

    if(dconv_3D < FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher convergence than threshold for 3D");

    if(dconv_6D < FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher convergence than threshold for 6D");

    if(thresh_3D != FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"3D thresh set unequal 3D thresh demanded");

    if(thresh_6D != FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"6D thresh set unequal 6D thresh demanded");

    if(econv < FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 3D");

    if(econv < FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 6D");

    if(econv < 0.1 * FunctionDefaults<3>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 3D (more than factor 10 difference)");

    if(econv < 0.1 * FunctionDefaults<6>::get_thresh()) warnings+=warning(world,"Demanded higher energy convergence than threshold for 6D (more than factor 10 difference)");

    // Check if the 6D thresholds are not too high
    if(thresh_6D < 1.e-3) warnings+=warning(world,"thresh_6D is smaller than 1.e-3");

    if(thresh_6D < tight_thresh_6D) warnings+=warning(world,"tight_thresh_6D is larger than thresh_6D");

    if(thresh_6D < tight_thresh_3D) warnings+=warning(world,"tight_thresh_3D is larger than thresh_3D");

    if(thresh_6D < 1.e-3) warnings+=warning(world,"thresh_6D is smaller than 1.e-3");

    if(thresh_Ue < 1.e-4) warnings+=warning(world,"thresh_Ue is smaller than 1.e-4");

    if(thresh_Ue > 1.e-4) warnings+=warning(world,"thresh_Ue is larger than 1.e-4");

    if(thresh_3D > 0.01 * thresh_6D) warnings+=warning(world,"Demanded 6D thresh is to precise compared with the 3D thresh");

    if(thresh_3D > 0.1 * thresh_6D) warnings+=warning(world,"Demanded 6D thresh is to precise compared with the 3D thresh");

    if(kain && kain_subspace == 0) warnings+=warning(world,"Demanded Kain solver but the size of the iterative subspace is set to zero");

    if(warnings > 0){
      if(world.rank() == 0) std::cout << warnings << "Warnings in parameters sanity check!\n\n";
    }else{
      if(world.rank() == 0) std::cout << "Sanity check for parameters passed\n\n" << std::endl;
    }
    if(restart == false && no_compute == true){
      warnings+=warning(world,"no_compute flag detected but no restart flag");
    }
  }


  CC_Parameters::CC_Parameters(const std::string& input,const double& low)
      : calculation(LRCC2_), lo(uninitialized), dmin(1.0), thresh_3D(uninitialized), tight_thresh_3D(uninitialized), thresh_6D(uninitialized), tight_thresh_6D(uninitialized), thresh_bsh_3D(
	  uninitialized), thresh_bsh_6D(uninitialized), thresh_poisson(uninitialized), thresh_f12(uninitialized), thresh_Ue(uninitialized), econv(uninitialized), econv_pairs(uninitialized), dconv_3D(
	  uninitialized), dconv_6D(uninitialized), iter_max(4), iter_max_3D(5), iter_max_6D(5), restart(false), no_compute(false), no_compute_gs(false), no_compute_response(false), no_compute_mp2(
	  false), no_compute_cc2(false), no_compute_cispd(false), no_compute_lrcc2(false), corrfac_gamma(1.0), output_prec(8), debug(false), plot(false), kain(false), freeze(0), test(false), decompose_Q(
	  true), QtAnsatz(false), excitations_(0), tda_guess_orbitals(0), tda_guess_mode("uninitialized"), tda_excitations(0), tda_guess_excitations(0), tda_iterating_excitations(0), tda_guess(
	  "uninitialized"), tda_energy_guess_factor(uninitialized), tda_dconv_guess(uninitialized), tda_dconv(uninitialized),     // tda_dconv_hard(uninitialized),
      tda_econv_guess(uninitialized), tda_econv(uninitialized),     //tda_econv_hard(uninitialized),
      tda_store_potential(true), tda_iter_max(25), tda_iter_guess(10), tda_homo_guess(false) {
    // get the parameters from the input file
    std::ifstream f(input.c_str());
    position_stream(f,"cc2");
    std::string s;
    // general operators thresh
    double thresh_operators=uninitialized;
    double thresh_operators_3D=uninitialized;
    double thresh_operators_6D=uninitialized;
    while(f >> s){
      //std::cout << "input tag is: " << s << std::endl;
      std::transform(s.begin(),s.end(),s.begin(),::tolower);
      //std::cout << "transformed input tag is: " << s << std::endl;
      if(s == "end") break;
      else if(s == "calculation"){
	std::string tmp;
	f >> tmp;
	calculation=assign_calctype(tmp);
      }else if(s == "lo") f >> lo;
      else if(s == "dmin") f >> dmin;
      else if(s == "thresh") f >> thresh_6D;
      else if(s == "thresh_3d") f >> thresh_3D;
      else if(s == "tight_thresh_3d") f >> tight_thresh_3D;
      else if(s == "thresh_6d") f >> thresh_6D;
      else if(s == "tight_thresh_6d") f >> tight_thresh_6D;
      else if(s == "debug") debug=true;
      else if(s == "plot") plot=true;
      else if(s == "econv") f >> econv;
      else if(s == "econv_pairs") f >> econv_pairs;
      else if(s == "dconv") f >> dconv_6D;
      else if(s == "dconv_3d") f >> dconv_3D;
      else if(s == "dconv_6d") f >> dconv_6D;
      else if(s == "thresh_operators" || s == "thresh_operator") f >> thresh_operators;
      else if(s == "thresh_operators_3d" || s == "thresh_operator_3d") f >> thresh_operators_3D;
      else if(s == "thresh_operators_6d" || s == "thresh_operator_6d") f >> thresh_operators_6D;
      else if(s == "thresh_bsh_3d") f >> thresh_bsh_3D;
      else if(s == "thresh_bsh_6d") f >> thresh_bsh_6D;
      else if(s == "thresh_poisson") f >> thresh_poisson;
      else if(s == "thresh_f12") f >> thresh_f12;
      else if(s == "thresh_ue") f >> thresh_Ue;
      else if(s == "freeze") f >> freeze;
      else if(s == "iter_max") f >> iter_max;
      else if(s == "iter_max_3d") f >> iter_max_3D;
      else if(s == "iter_max_6d") f >> iter_max_6D;
      else if(s == "kain") kain=true;
      else if(s == "kain_subspace") f >> kain_subspace;
      else if(s == "freeze") f >> freeze;
      else if(s == "test") test=true;
      else if(s == "corrfac" || s == "corrfac_gamma" || s == "gamma") f >> corrfac_gamma;
      else if(s == "decompose_q") decompose_Q=true;
      else if(s == "restart") restart=true;
      else if(s == "no_compute"){
	no_compute=true;
	no_compute_gs=true;
	no_compute_mp2=true;
	no_compute_cispd=true;
	no_compute_response=true;
      }else if(s == "no_compute_gs"){
	no_compute_gs=true;
	no_compute_mp2=true;
	no_compute_cc2=true;
      }else if(s == "no_compute_response"){
	no_compute_response=true;
	no_compute_cispd=true;
	no_compute_lrcc2=true;
      }else if(s == "no_compute_cc2") no_compute_cc2=true;
      else if(s == "no_compute_cispd") no_compute_cispd=true;
      else if(s == "no_compute_lrcc2") no_compute_lrcc2=true;
      else if(s == "no_compute_mp2") no_compute_mp2=true;
      else if(s == "only_pair"){
	size_t tmp1, tmp2;
	f >> tmp1;
	f >> tmp2;
	std::cout << "found only pair in the world: " << tmp1 << ", " << tmp2 << "\n";
	only_pair=std::make_pair(tmp1,tmp2);
	MADNESS_ASSERT(!(tmp1 > tmp2));
      }else if(s == "excitation"){
	size_t tmp;
	f >> tmp;
	excitations_.push_back(tmp);
      }else if(s == "qtansatz") QtAnsatz=true;
      else if(s == "full_residue") decompose_Q=false;
      else if(s == "tda_guess_orbitals") f >> tda_guess_orbitals;
      else if(s == "tda_guess_mode") f >> tda_guess_mode;
      else if(s == "tda_guess_excitations") f >> tda_guess_excitations;
      else if(s == "tda_excitations") f >> tda_excitations;
      else if(s == "tda_iterating_excitations") f >> tda_iterating_excitations;
      else if(s == "tda_guess") f >> tda_guess;
      else if(s == "tda_energy_guess_factor") f >> tda_energy_guess_factor;
      else if(s == "tda_dconv_guess") f >> tda_dconv_guess;
      else if(s == "tda_dconv") f >> tda_dconv;
      else
      //else if ( s == "tda_dconv_hard")f >> tda_dconv_hard;
      if(s == "tda_econv_guess") f >> tda_econv_guess;
      else if(s == "tda_econv") f >> tda_econv;
      else
      //else if ( s == "tda_econv_hard") f >> tda_econv_hard;
      if(s == "tda_store_potential") f >> tda_store_potential;
      else if(s == "tda_iter_max") f >> tda_iter_max;
      else if(s == "tda_iter_guess") f >> tda_iter_guess;
      else if(s == "tda_homo_guess") tda_homo_guess=true;
      else if(s == "tda_exop" || s == "exop"){
	std::string tmp;
	char buf[1024];
	f.getline(buf,sizeof(buf));
	tmp=buf;
	tda_exops.push_back(tmp);
      }else{
	std::cout << "Unknown Keyword: " << s << "\n";
	continue;
      }
    }
    // set defaults
    if(!kain) kain_subspace=0;

    // set all parameters that were not explicitly given
    if(lo == uninitialized) lo=1.e-7;

    if(thresh_6D == uninitialized) thresh_6D=1.e-3;

    if(tight_thresh_6D == uninitialized) tight_thresh_6D=thresh_6D * 0.1;

    if(thresh_3D == uninitialized) thresh_3D=thresh_6D * 0.01;

    if(tight_thresh_3D == uninitialized) tight_thresh_3D=thresh_3D * 0.1;

    if(thresh_operators == uninitialized) thresh_operators=1.e-6;

    if(thresh_operators_3D == uninitialized) thresh_operators_3D=thresh_operators;

    if(thresh_operators_6D == uninitialized) thresh_operators_6D=thresh_operators;

    if(thresh_bsh_3D == uninitialized) thresh_bsh_3D=thresh_operators_3D;

    if(thresh_bsh_6D == uninitialized) thresh_bsh_6D=thresh_operators_6D;

    if(thresh_poisson == uninitialized) thresh_poisson=thresh_operators_3D;

    if(thresh_f12 == uninitialized) thresh_f12=thresh_operators_3D;

    if(thresh_Ue == uninitialized) thresh_Ue=tight_thresh_6D;

    if(dconv_6D == uninitialized) dconv_6D=thresh_6D;

    if(dconv_3D == uninitialized) dconv_3D=dconv_6D;

    if(econv == uninitialized) econv=0.1 * dconv_6D;

    if(econv_pairs == uninitialized) econv_pairs=econv;

    if(iter_max_6D == uninitialized) iter_max_6D=10;

    if(iter_max_3D == uninitialized) iter_max_3D=iter_max_6D;

    // set the thresholds
    FunctionDefaults<3>::set_thresh(thresh_3D);
    FunctionDefaults<6>::set_thresh(thresh_6D);
    if(thresh_3D < 1.1e-1) output_prec=3;

    if(thresh_3D < 1.1e-2) output_prec=4;

    if(thresh_3D < 1.1e-3) output_prec=5;

    if(thresh_3D < 1.1e-4) output_prec=6;

    if(thresh_3D < 1.1e-5) output_prec=7;

    if(thresh_3D < 1.1e-6) output_prec=8;

    std::cout.precision(output_prec);
    // set the default TDA parameters
    if(tda_guess == "uninitialized") tda_guess="big_fock_3";

    //if(tda_guess_mode=="uninitialized") tda_guess_mode = "projected";
    if(tda_excitations == 0) tda_excitations=1;

    if(tda_guess_excitations == 0) tda_guess_excitations=tda_excitations;

    if(tda_iterating_excitations == 0) tda_iterating_excitations=tda_guess_excitations;

    if(tda_energy_guess_factor == uninitialized) tda_energy_guess_factor=0.99;

    if(tda_dconv_guess == uninitialized) tda_dconv_guess=1.0;

    if(tda_dconv == uninitialized) tda_dconv=thresh_3D * 10.0;

    if(tda_econv_guess == uninitialized) tda_econv_guess=1.e-1;

    if(tda_econv == uninitialized) tda_econv=thresh_3D;

    if(tda_dconv_hard == uninitialized) tda_econv_hard=tda_econv;

    ;
    if(tda_econv_hard == uninitialized) tda_dconv_hard=tda_dconv;

    ;
    if(no_compute == true && restart == false) restart=true;
  }


  CC_Parameters::CC_Parameters(const CC_Parameters& other)
      : calculation(other.calculation), lo(other.lo), dmin(other.dmin), thresh_3D(other.thresh_3D), tight_thresh_3D(other.tight_thresh_3D), thresh_6D(other.thresh_6D), tight_thresh_6D(
	  other.tight_thresh_6D), thresh_bsh_3D(other.thresh_bsh_3D), thresh_bsh_6D(other.thresh_bsh_6D), thresh_poisson(other.thresh_poisson), thresh_f12(other.thresh_f12), thresh_Ue(
	  other.thresh_Ue), econv(other.econv), econv_pairs(other.econv_pairs), dconv_3D(other.dconv_3D), dconv_6D(other.dconv_6D), iter_max(other.iter_max), iter_max_3D(other.iter_max_3D), iter_max_6D(
	  other.iter_max_6D), restart(other.restart), no_compute(other.no_compute), no_compute_gs(other.no_compute_gs), no_compute_response(other.no_compute_response), no_compute_mp2(
	  other.no_compute_response), no_compute_cc2(other.no_compute_mp2), no_compute_cispd(other.no_compute_cispd), no_compute_lrcc2(other.no_compute_lrcc2), corrfac_gamma(other.corrfac_gamma), output_prec(
	  other.output_prec), debug(other.debug), plot(other.plot), kain(other.kain), kain_subspace(other.kain_subspace), freeze(other.freeze), test(other.test), decompose_Q(other.decompose_Q), QtAnsatz(
	  other.QtAnsatz), excitations_(other.excitations_), tda_guess_orbitals(0), tda_guess_mode(other.tda_guess_mode), tda_excitations(other.tda_excitations), tda_guess_excitations(
	  other.tda_guess_excitations), tda_iterating_excitations(other.tda_iterating_excitations), tda_guess(other.tda_guess), tda_energy_guess_factor(other.tda_energy_guess_factor), tda_dconv_guess(
	  other.tda_dconv_guess), tda_dconv(other.tda_dconv), tda_dconv_hard(other.tda_dconv_hard), tda_econv_guess(other.tda_econv_guess), tda_econv(other.tda_econv), tda_econv_hard(
	  other.tda_econv_hard), tda_store_potential(other.tda_store_potential), tda_iter_max(other.tda_iter_max), tda_iter_guess(other.tda_iter_guess), tda_homo_guess(other.tda_homo_guess), tda_exops(
	  other.tda_exops) {
  }

}; // end of namespace madness
