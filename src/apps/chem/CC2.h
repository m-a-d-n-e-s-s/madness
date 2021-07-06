/*
 * CC2.h
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#ifndef CC2_H_
#define CC2_H_

#include <chem/projector.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/CCPotentials.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include <chem/TDHF.h>
#include <madness/mra/nonlinsol.h>

namespace madness {

  class CC2{
  public:

	  CC2(World &world_,const CCParameters& param,const std::shared_ptr<Nemo> nemo_)
	 	: world(world_),
	 	parameters(param),
	 	nemo(nemo_),
	 	CCOPS(world,*nemo,parameters),
	 	output(CCOPS.output)
	   {
		  parameters.sanity_check(world);
	   }

    CC2(World &world_,const std::string &inputFileName,const std::shared_ptr<Nemo> nemo_)
	: world(world_),
	parameters(inputFileName,nemo_->get_calc()->param.lo()),
	nemo(nemo_),
	CCOPS(world,*nemo,parameters),
	output(CCOPS.output)
  {

      output.section("CC2 Class has been initialized with the following parameters");
      // set the threshholds
      // Set Protocoll
      output("Set Protocol 3D");
      nemo_->get_calc()->set_protocol < 3 > (world, parameters.thresh_3D);
      output("Set Protocol 6D");
      nemo_->get_calc()->set_protocol < 6 > (world, parameters.thresh_6D);

      FunctionDefaults<3>::set_thresh(parameters.thresh_3D);
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      // Make shure that k is the same in 3d and 6d functions
      FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
      // by default SCF sets the truncate_mode to 1
      FunctionDefaults<3>::set_truncate_mode(3);
      FunctionDefaults<6>::set_truncate_mode(3);
      // set initial level to 3 (important to avoid undersampling issues)
      // the "real initial level" is then 2 since the initial level gets substracted by one if refine is true (which is the default)
      FunctionDefaults<6>::set_initial_level(2);
      FunctionDefaults<3>::set_initial_level(2);
      FunctionDefaults<6>::set_special_level(FunctionDefaults<6>::set_length_scale(parameters.dmin,FunctionDefaults<6>::get_k()));
      FunctionDefaults<3>::set_special_level(FunctionDefaults<3>::set_length_scale(parameters.dmin,FunctionDefaults<3>::get_k()));
      parameters.information(world);
      parameters.sanity_check(world);

      std::string nuc="???";
      if(nemo->ncf->type() == NuclearCorrelationFactor::None) nuc="None";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::GaussSlater) nuc="GaussSlater";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::GradientalGaussSlater) nuc="GradientalGaussSlater";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::LinearSlater) nuc="LinearSlater";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::Polynomial) nuc="Polynomial";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::Slater) nuc="Slater";
      else if(nemo->ncf->type() == NuclearCorrelationFactor::Two) nuc="Two";
      if(world.rank() == 0) std::cout << "Nuclear Correlation Factor is " << nuc << std::endl;

    }

    void
    plot(const real_function_3d &f,const std::string &msg="unspecified function") const {
      plot_plane(world,f,msg);
      output("Plotted " + msg);
    }
    /// Check energy convergence: Creates the difference between two vectors and compares against given thresh in parameters
    bool
    check_energy_convergence(const std::vector<double> &current,const std::vector<double> &updated) const {
      if(current.size() != updated.size()) MADNESS_EXCEPTION("error in energy convergence check: different sizes in vectors",1);
      bool conv=true;
      std::vector<double> diff(current.size(),0.0);
      for(size_t i=0; i < current.size(); i++){
	double diffi=updated[i] - current[i];
	diff[i]=diffi;
	if(diffi > parameters.econv) conv=false;
      }
      if(world.rank() == 0){
	std::cout << "\n\n";
	std::cout << "Pair Correlation Energies: New, Old, Diff\n";
	for(size_t i=0; i < current.size(); i++)
	  std::cout << updated[i] << ", " << current[i] << ", " << diff[i] << std::endl;
	std::cout << "\n\n";
      }
      return conv;
    }
    /// make consistency tests
    bool
    test() const;
    /// The World
    World &world;
    /// The electronic Correlation Factor, has to be initialized before parameters so that parameters has the right gamma value
    //CorrelationFactor correlationfactor;
    /// Structure holds all the parameters used in the CC2 calculation
    const CCParameters parameters;
    /// The SCF Calculation
    const std::shared_ptr<Nemo> nemo;
    /// The CC Operator Class
    CCPotentials CCOPS;
    /// Formated Output (same as used in CC2Potentials structure)
    CCMessenger& output;

    /// solve the CC2 ground state equations, returns the correlation energy
    void
    solve();
    std::vector<CC_vecfunction>
    solve_ccs();
    /// solve the MP2 equations (uncoupled -> Canonical Orbitals)
    double
    solve_mp2(Pairs<CCPair> &doubles);
    double
    solve_cc2(CC_vecfunction &tau,Pairs<CCPair> &u);
    double
    solve_cispd(Pairs<CCPair> &doubles,const Pairs<CCPair> &mp2_pairs,const CC_vecfunction & cis_singles);

    /// convencience function to iterate the CC2 ground state singles,
    /// makes the right call on the iterate_singles functions
    bool
    iterate_cc2_singles(CC_vecfunction &singles,Pairs<CCPair> &doubles) {
      CCOPS.clear_potentials(singles);
      Pairs<CCPair> empty;
      return iterate_singles(singles,CC_vecfunction(RESPONSE),doubles,empty,CT_CC2,parameters.iter_max_3D);
    }
    bool
    iterate_adc2_singles(Pairs<CCPair> &mp2, CC_vecfunction &singles,Pairs<CCPair> &x){
      MADNESS_ASSERT(singles.type==RESPONSE);
      CCOPS.clear_potentials(singles);
      return iterate_singles(singles, CC_vecfunction(UNDEFINED),mp2,x,CT_ADC2,parameters.iter_max_3D);
    }
    bool
    iterate_lrcc2_singles(CC_vecfunction &cc2_s, Pairs<CCPair> & cc2_d, CC_vecfunction &lrcc2_s, Pairs<CCPair> lrcc2_d){
      MADNESS_ASSERT(cc2_s.type==PARTICLE);
      MADNESS_ASSERT(lrcc2_s.type==RESPONSE);
      CCOPS.clear_potentials(lrcc2_s);
      return iterate_singles(lrcc2_s,cc2_s,cc2_d,lrcc2_d,CT_LRCC2,parameters.iter_max_3D);
    }
    /// convencience function to iterate the CCS Response singles,
    /// makes the right call on the iterate_singles functions
    bool
    iterate_ccs_singles(CC_vecfunction &x){
      Pairs<CCPair> empty;
      CCOPS.clear_potentials(x);
      return iterate_singles(x,CC_vecfunction(PARTICLE),empty,empty,CT_LRCCS,1);
    }
    bool
    /// Iterates the singles equations for CCS, CC2, LRCC2
    /// The corresponding regulairzation tails of the doubles are updated in every iteration (therefore doubles are not marked as const)
    /// @param[in] : singles, the singles that are iterated
    /// @param[in] : singles2, additional set of singles for potential (in LRCC2 this are the Ground State singles)
    /// @param[in] : gs_doubles, Ground State doubles (Needed for CC2 and LRCC2)
    /// @param[in] : ex_doubles, Excited State doubles (Needed for LRCC2)
    /// @param[in] : ctype: the calculation type: CCS, CC2, CC2_response_
    /// @param[in] : maxiter: maxmial number of iterations
    /// @param[out]: true if the overall change of the singles is below 10*donv_6D
    iterate_singles(CC_vecfunction &singles,const CC_vecfunction singles2,Pairs<CCPair> &gs_doubles, Pairs<CCPair> &ex_doubles,const CalcType ctype,const std::size_t maxiter) {
      output.subsection("Iterate " + assign_name(ctype) + "-Singles");
      CCTimer time_all(world,"Overall Iteration of " + assign_name(ctype) + "-Singles");
      bool converged=true;

      CC_vecfunction old_singles(singles);
      for(auto& tmp : singles.functions)
	old_singles(tmp.first).function=copy(tmp.second.function);

      // KAIN solver
      typedef allocator<double, 3> allocT;
      typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
      allocT alloc(world,singles.size());
      solverT solver(allocT(world,singles.size()));
      solver.do_print=(world.rank() == 0);

      for(size_t iter=0; iter < maxiter; iter++){
	output.subsection("Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
	CCTimer time(world,"Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
	double omega=0.0;
	if(ctype == CT_LRCC2) omega=singles.omega;
	else if(ctype == CT_LRCCS) omega=singles.omega;
	else if(ctype == CT_ADC2) omega=singles.omega;

	// consistency check
	switch(ctype){
	  case CT_CC2:
	    if(singles.type != PARTICLE) output.warning("iterate_singles: CC2 demanded but singles are not of type PARTICLE");
	    break;
	  case CT_MP2:
	    MADNESS_EXCEPTION("Demanded Singles Calculation for MP2 ????",1);
	    break;
	  case CT_LRCC2:
	    if(singles.type != RESPONSE or singles2.type != PARTICLE) output.warning("iterate_singles: CC2_response_ singles have wrong types");
	    break;
	  case CT_LRCCS:
	    if(singles.type != RESPONSE) output.warning("iterate_singles: CCS_response_ singles have wrong types");
	    break;
	  case CT_CISPD:
	    MADNESS_EXCEPTION("Demanded Singles Calculation for CIS(D)",1);
	    break;
	  case CT_ADC2:
	    MADNESS_ASSERT(singles.type==RESPONSE);
	    break;
	  case CT_TEST:
	    MADNESS_EXCEPTION("Iterate Singles not implemented for Experimental calculation",1);
	    break;
	  default:
	    MADNESS_EXCEPTION(("Unknown calculation type in iterate singles: " + assign_name(ctype)).c_str(),1);
	}

	// get potentials
	CCTimer time_V(world,assign_name(ctype) + "-Singles Potential");
	vector_real_function_3d V;
	if(ctype == CT_CC2) V=CCOPS.get_CC2_singles_potential_gs(singles,gs_doubles);
	else if(ctype == CT_LRCC2) V=CCOPS.get_CC2_singles_potential_ex(singles2,gs_doubles,singles,ex_doubles);
	else if(ctype == CT_LRCCS) V=CCOPS.get_CCS_potential_ex(singles);
	else if(ctype == CT_ADC2) V=CCOPS.get_ADC2_singles_potential(gs_doubles,singles,ex_doubles);
	else  MADNESS_EXCEPTION("iterate singles: unknown type",1);
	time_V.info(true,norm2(world,V));

	if(ctype == CT_LRCCS or ctype==CT_LRCC2 or ctype==CT_ADC2){
//	  const double expv=CCOPS.compute_cis_expectation_value(singles,V);
//	  if(world.rank() == 0) std::cout << "Current CCS/CIS or CC2 Expectation Value " << expv << "\n";
//	  if(world.rank() == 0) std::cout << "using expectation-value for bsh-operator\n";
//	  singles.omega = expv;
	  omega = singles.omega; // computed with the potential
	}

	scale(world,V,-2.0);
	truncate(world,V);

	// make bsh operators
	CCTimer time_makebsh(world,"Make G-Operators");
	std::vector < std::shared_ptr<SeparatedConvolution<double, 3> > > G(singles.size());
	for(size_t i=0; i < G.size(); i++){
	  const double bsh_eps=CCOPS.get_orbital_energies()[i + parameters.freeze] + omega;
	  G[i]=std::shared_ptr < SeparatedConvolution<double, 3> > (BSHOperatorPtr3D(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_3D));
	}
	world.gop.fence();
	time_makebsh.info();

	// apply bsh operators
	CCTimer time_applyG(world,"Apply G-Operators");
	vector_real_function_3d GV=apply<SeparatedConvolution<double, 3>, double, 3>(world,G,V);
	world.gop.fence();
	time_applyG.info();

	// apply Q-Projector to result
	GV=CCOPS.apply_Qt(GV,CCOPS.mo_ket());

	// Normalize Singles if it is excited state
	if(ctype==CT_LRCCS or ctype==CT_LRCC2 or ctype==CT_ADC2){
	  output("Normalizing new singles");
	  const vector_real_function_3d x = GV;
	  const vector_real_function_3d xbra = mul(world,nemo->ncf->square(),GV);
	  const double norm = sqrt(inner(world,xbra,x).sum());
	  if(world.rank()==0) std::cout << " Norm was " <<std::fixed<< std::setprecision(parameters.output_prec) << norm << "\n";
	  scale(world,GV,1.0/norm);
	} else output("Singles not normalized");

	// residual
	const vector_real_function_3d residual=sub(world,singles.get_vecfunction(),GV);

	// information
	const Tensor<double> R2xinnerx=inner(world,mul(world,nemo->ncf->square(),singles.get_vecfunction()),singles.get_vecfunction());
	const Tensor<double> R2GVinnerGV=inner(world,mul(world,nemo->ncf->square(),GV),GV);
	const Tensor<double> R2rinnerr=inner(world,mul(world,nemo->ncf->square(),residual),residual);
	const double R2vector_error=sqrt(R2rinnerr.sum());

	// print information
	if(world.rank() == 0) std::cout << "\n\n-----Results of current interation:-----\n";
	if(world.rank() == 0) std::cout << "\nName: ||" << singles.name() << "||, ||GV" << singles.name() << ", ||residual||" << "\n";
	if(world.rank() == 0)
	std::cout << singles.name() << ": " << std::scientific << std::setprecision(parameters.output_prec)
				    << sqrt(R2xinnerx.sum()) <<", " << sqrt(R2GVinnerGV.sum())<<", "<< sqrt(R2rinnerr.sum())
				    << "\n----------------------------------------\n";
	for(size_t i=0; i < GV.size(); i++){
	  if(world.rank() == 0) std::cout << singles(i + parameters.freeze).name() << ": " << std::scientific << std::setprecision(parameters.output_prec)
	  << sqrt(R2xinnerx(i))<<", " << sqrt(R2GVinnerGV(i))<<", "  << sqrt(R2rinnerr(i)) << "\n";
	}
	if(world.rank() == 0) std::cout << "\n----------------------------------------\n\n";

	// make second order update (only for response)
	if(ctype == CT_LRCC2 or ctype == CT_LRCCS){
	  output("\nMake 2nd order energy update:");
	  // include nuclear factors
	  {
	    vector_real_function_3d bra_res=mul(world,nemo->ncf->square(),residual);
	    vector_real_function_3d bra_GV=mul(world,nemo->ncf->square(),GV);
	    double Rtmp=inner(world,bra_res,V).sum();
	    double Rtmp2=inner(world,bra_GV,GV).sum();
	    const double Rdelta=(0.5 * Rtmp / Rtmp2);
	    double old_omega=omega;
	    output("Delta-Update is not used");
	    if(world.rank() == 0)
	      std::cout << "omega, old_omega, delta" << std::fixed << std::setprecision(parameters.output_prec + 2) << omega << ", " << old_omega << ", " << Rdelta << "\n\n";
	  }

	}

	// update singles
	singles.omega=omega;
	vector_real_function_3d new_singles=GV;
	if(parameters.kain) new_singles=solver.update(singles.get_vecfunction(),residual);
	print_size(world,new_singles,"new_singles");
	truncate(world,new_singles);
	print_size(world,new_singles,"new_singles");
	for(size_t i=0; i < GV.size(); i++){
	  singles(i + parameters.freeze).function=copy(new_singles[i]);
	}

	// update intermediates
	CCOPS.update_intermediates(singles);

	// update reg_residues of doubles
	//if(ctype==CC2_) update_reg_residues_gs(singles,gs_doubles);
	//else if(ctype==LRCC2_) update_reg_residues_ex(singles2,singles,ex_doubles);

	converged=(R2vector_error < parameters.dconv_3D);

	time.info();
	if(converged) break;
	if(ctype==CT_LRCCS) break; // for CCS just one iteration to check convergence
      }
      time_all.info();

      // Assign the overall changes
      bool no_change=true;
      if(world.rank() == 0) std::cout << "Change in Singles functions after all the CC2-Single-Microiterations" << std::endl;
      for(auto& tmp : singles.functions){
	const double change = (tmp.second.function - old_singles(tmp.first).function).norm2();
	tmp.second.current_error= change;
	if(change>parameters.dconv_3D) no_change=false;
	if(world.rank() == 0) std::cout << "Change of " << tmp.second.name() << "=" << tmp.second.current_error << std::endl;
      }
      // update reg_residues of doubles
      if(ctype==CT_CC2) update_reg_residues_gs(singles,gs_doubles);
      else if(ctype==CT_LRCC2) update_reg_residues_ex(singles2,singles,ex_doubles);

      //CCOPS.plot(singles);
      singles.save_functions();
      if(no_change) output("Change of Singles was below  = " +std::to_string(parameters.dconv_3D)+"!");
      return no_change;
    }


    bool initialize_singles(CC_vecfunction &singles,const FuncType type,const int ex= -1) const;
    bool initialize_pairs(Pairs<CCPair> &pairs,const CCState ftype, const CalcType ctype, const CC_vecfunction &tau, const CC_vecfunction &x, const size_t extitation=0)const;
    void update_reg_residues_gs(const CC_vecfunction& singles,Pairs<CCPair> &doubles)const;
    void update_reg_residues_ex(const CC_vecfunction& singles, const CC_vecfunction& response,Pairs<CCPair> &doubles)const;

    /// Iterates a pair of the CC2 doubles equations
    bool
    iterate_pair(CCPair & pair,const CC_vecfunction &singles=CC_vecfunction(UNDEFINED)) const;
    bool
    iterate_adc2_pairs(Pairs<CCPair> &cispd,const CC_vecfunction & ccs);
    bool
    iterate_lrcc2_pairs(const CC_vecfunction & cc2_s,const Pairs<CCPair>& cc2_d,const CC_vecfunction lrcc2_s,Pairs<CCPair> &lrcc2_d);

    bool update_constant_part_mp2(CCPair &pair){
      MADNESS_ASSERT(pair.ctype==CT_MP2);
      MADNESS_ASSERT(pair.type==GROUND_STATE);
      if(pair.constant_part.is_initialized()) return false;

      // make screening Operator
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * CCOPS.get_epsilon(pair.i,pair.j)),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.modified()=true;

      const CCFunction& moi = CCOPS.mo_ket(pair.i);
      const CCFunction& moj = CCOPS.mo_ket(pair.j);

      pair.constant_part = CCOPS.make_constant_part_mp2(moi,moj,&Gscreen);
      save(pair.constant_part,pair.name()+"_const");
      return true;
    }

    bool update_constant_part_cc2_gs(const CC_vecfunction &tau,CCPair &pair){
      MADNESS_ASSERT(pair.ctype==CT_CC2);
      MADNESS_ASSERT(pair.type==GROUND_STATE);
      // make screening Operator
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * pair.bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.modified()=true;

      if(parameters.QtAnsatz)pair.constant_part = CCOPS.make_constant_part_cc2_Qt_gs(pair,tau,&Gscreen);
      else pair.constant_part = CCOPS.make_constant_part_cc2_gs(pair,tau,&Gscreen);
      save(pair.constant_part,pair.name()+"_const");
      return true;
    }

    bool update_constant_part_cispd(const CC_vecfunction &ccs, CCPair &pair){
      MADNESS_ASSERT(pair.ctype==CT_CISPD);
      MADNESS_ASSERT(pair.type==EXCITED_STATE);
      MADNESS_ASSERT(pair.bsh_eps==CCOPS.get_epsilon(pair.i,pair.j)+ccs.omega);
      if(pair.constant_part.is_initialized()) return false; // the CIS(D) constant part does not change because there is no singles iteration (like MP2)
      // make screening Operator
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * pair.bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.modified()=true;

      if(parameters.QtAnsatz) pair.constant_part = CCOPS.make_constant_part_cispd_Qt(pair,ccs,&Gscreen);
      else  pair.constant_part = CCOPS.make_constant_part_cispd(pair,ccs,&Gscreen);
      save(pair.constant_part,pair.name()+"_const");
      return true;

    }

    bool update_constant_part_adc2(const CC_vecfunction &ccs, CCPair &pair){
      std::cout << assign_name(pair.ctype);
      MADNESS_ASSERT(pair.ctype==CT_ADC2);
      MADNESS_ASSERT(pair.type==EXCITED_STATE);
      MADNESS_ASSERT(pair.bsh_eps==CCOPS.get_epsilon(pair.i,pair.j)+ccs.omega);
      // make screening Operator
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * pair.bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.modified()=true;

      if(parameters.QtAnsatz) pair.constant_part = CCOPS.make_constant_part_cispd_Qt(pair,ccs,&Gscreen);
      else  pair.constant_part = CCOPS.make_constant_part_cispd(pair,ccs,&Gscreen);
      save(pair.constant_part,pair.name()+"_const");
      return true;

    }

    bool update_constant_part_lrcc2(CCPair &pair, const CC_vecfunction &tau,const CC_vecfunction &x){
      MADNESS_ASSERT(pair.ctype==CT_LRCC2);
      MADNESS_ASSERT(tau.type==PARTICLE);
      MADNESS_ASSERT(x.type==RESPONSE);

      // make screening Operator
      real_convolution_6d Gscreen=BSHOperator<6>(world,sqrt(-2.0 * pair.bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
      Gscreen.modified()=true;

      if(parameters.QtAnsatz)pair.constant_part = CCOPS.make_constant_part_cc2_Qt_ex(pair,tau,x,&Gscreen);
      else pair.constant_part = CCOPS.make_constant_part_cc2_ex(pair,tau,x,&Gscreen);
      save(pair.constant_part,pair.name()+"_const");
      return true;
    }

  };



} /* namespace madness */

#endif /* CC2_H_ */
