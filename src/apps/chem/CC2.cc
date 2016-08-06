/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#include "CC2.h"
namespace madness {

  /// solve the CC2 ground state equations, returns the correlation energy
  void
  CC2::solve() {
    if(parameters.calculation!=LRCCS_) MADNESS_EXCEPTION("Only Linear-Response for CCS possible right now, type keyword: calculation lrccs in the cc2 section",1);
    solve_ccs();
  }

  // Solve the CCS equations for the ground state (debug potential and check HF convergence)
  std::vector<std::pair<CC_vecfunction,double> > CC2::solve_ccs() {
    output_section("SOLVE CCS");
    // this is that the output of CIS is not screwed
    std::cout << std::setw(5) << std::setfill(' ') << std::endl;
    TDA CCS_Solver(world,nemo,nemo.get_calc()->amo,"input");
    {
      xfunctionsT guess;
      CCS_Solver.solve_guess(guess);
      guess.clear();
    }
        {
          xfunctionsT solve = CCS_Solver.get_converged_xfunctions();
          CCS_Solver.solve(solve);
          solve.clear();
        }
        {
          xfunctionsT final = CCS_Solver.get_converged_xfunctions();
          CCS_Solver.solve_sequential(final);
          final.clear();
        }
    output_section("SOLVE CCS finished");
    std::vector<std::pair<CC_vecfunction,double>> ccs_vectors;
    {
      xfunctionsT ccs_vectors_tmp = CCS_Solver.get_converged_xfunctions();
      for(size_t i=0;i<ccs_vectors_tmp.size();i++){
	const vecfuncT xtmp = ccs_vectors_tmp[i].x;
	const CC_vecfunction x(copy(world,xtmp),RESPONSE,parameters.freeze);
	const double omega = ccs_vectors_tmp[i].omega;
	ccs_vectors.push_back(std::make_pair(x,omega));
      }
    }
    return ccs_vectors;
  }

  double
  CC2::solve_mp2(Pairs<CC_Pair> &doubles) {
    output_section("Solve MP2");
    if(not parameters.no_compute){
      bool mp2_converged=true;
      for(auto& tmp_pair : doubles.allpairs){
	update_constant_part_mp2(tmp_pair.second);
	bool pair_converged=iterate_pair(tmp_pair.second);
	if(not pair_converged) mp2_converged=false;
	else tmp_pair.second.store_pair(world,"mp2_converged_");
      }
      print_results(doubles,initialize_cc2_singles());
    }else output("keyword: no_compute found");
    return get_correlation_energy(doubles);
  }

  double CC2::solve_mp2_nonorthogonal(Pairs<CC_Pair> &doubles) {
    output_section("Test f12-projection");
    CCOPS.test_inverse_correlation_factor();
    CCOPS.test_f12_projections();
    output_section("Solve MP2");
    if(not parameters.no_compute){
      bool mp2_converged=true;
      for(auto& tmp_pair : doubles.allpairs){
	const double ij_gf_ij = CCOPS.make_ijgfxy(tmp_pair.second.i,tmp_pair.second.j,CCOPS.mo_ket(tmp_pair.second.i).function,CCOPS.mo_ket(tmp_pair.second.j).function);
	const double ji_gf_ij = CCOPS.make_ijgfxy(tmp_pair.second.j,tmp_pair.second.i,CCOPS.mo_ket(tmp_pair.second.i).function,CCOPS.mo_ket(tmp_pair.second.j).function);
	tmp_pair.second.constant_energy = 2.0*ij_gf_ij - ji_gf_ij;
	if(world.rank()==0){
	  std::cout << "Pair " << tmp_pair.second.name() << "\n";
	  std::cout << "ij_gf_ij=" << ij_gf_ij << "\n";
	  std::cout << "ji_gf_ij=" << ji_gf_ij << "\n";
	}
	bool pair_converged=iterate_nonorthogonal_pair(tmp_pair.second);
	if(not pair_converged) mp2_converged=false;
      }
      print_results(doubles,initialize_cc2_singles());
    }else output("keyword: no_compute found");
    return get_correlation_energy(doubles);
  }

  double
  CC2::solve_cispd(Pairs<CC_Pair> &doubles,const Pairs<CC_Pair> &mp2_pairs, const CC_vecfunction & cis_singles, const double cis_omega) {
    output_section("Solve CIS(D) for CIS-Exctation energy " + std::to_string(cis_omega));

    CCOPS.update_response_intermediates(cis_singles);
    bool cispd_converged=true;
    output("\n\n Making empty CC2-gs-singles for CIS(D) calculation");
    CC_vecfunction empty_singles(zero_functions<double,3>(world,cis_singles.size()),PARTICLE,parameters.freeze);
    CCOPS.update_intermediates(empty_singles);

    if(not parameters.no_compute){
    for(auto& tmp_pair : doubles.allpairs){
      tmp_pair.second.current_energy = cis_omega;
      update_constant_part_cc2_response(tmp_pair.second,CC_vecfunction(PARTICLE),cis_singles);
      bool pair_converged=iterate_pair(tmp_pair.second);
      if(not pair_converged) cispd_converged=false;
    }
    }

    output_section("Calculating the CIS(D) Excitation Energies");
    const double result = CCOPS.compute_cispd_energy_correction(cis_singles,mp2_pairs,doubles);
    if(world.rank()==0){
      std::cout <<"Excitation Energy CIS   ="<<std::fixed << std::setprecision(parameters.output_prec) <<  cis_omega << "\n";
      std::cout <<"CIS(D) Correction       ="<<std::fixed << std::setprecision(parameters.output_prec) << result << "\n";
      std::cout <<"Excitation Energy CIS(D)="<<std::fixed << std::setprecision(parameters.output_prec) << result+ cis_omega << "\n";
    }
    return result+ cis_omega;
  }

  double CC2::solve_cc2_response(const CC_vecfunction &tau,const Pairs<CC_Pair> &u,CC_vecfunction x,Pairs<CC_Pair> &chi){
    output_section("Beginning the CC2-Response Iterations");
    CCOPS.plot(tau);
    double omega = x.omega;
    CC_Timer time(world,"CC2-Response");
    CCOPS.update_response_intermediates(x);
    bool singles_converged = true;
    bool doubles_converged = true;
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      output_subsection("Mactoiteration " + std::to_string(iter) + " with excitation energy " + std::to_string(omega));
      CC_Timer timer_iter_all(world,"Macroiteration " + stringify(iter));

      CCOPS.remove_stored_response_singles_potentials();
      std::vector<Pairs<CC_Pair> > doubles;
      doubles.push_back(u);
      doubles.push_back(chi);
      singles_converged=iterate_singles(x,tau,doubles,LRCC2_);
      const double new_omega = x.omega;
      const double diff = new_omega - omega;
      bool energy_converged = (fabs(diff)<parameters.econv);

      CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Response of Doubles");
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
      for(auto& tmp_pair : chi.allpairs){
	update_constant_part_cc2_response(tmp_pair.second,tau,x);
	bool pair_converged=iterate_pair(tmp_pair.second);
	if(not pair_converged) doubles_converged=false;
      }

      output("Performance Overview of Iteration " + stringify(iter));
      if(world.rank() == 0) CCOPS.performance_S.info_last_iter();
      if(world.rank() == 0) CCOPS.performance_D.info_last_iter();

      CCOPS.remove_stored_response_singles_potentials();


      if(world.rank() == 0){
	output("End of Macroiteration " + stringify(iter));
	std::cout << "singles converged: " << singles_converged << std::endl;
	std::cout << "doubles converged: " << doubles_converged << std::endl;
	std::cout << "Excitation energy converged: " << energy_converged << std::endl;
      }



      if(singles_converged and doubles_converged and energy_converged){
	output("Singles and Doubles Converged (!)");
	output("Testing if singles do not change anymore");
	vecfuncT old_x=x.get_vecfunction();
	std::vector<Pairs<CC_Pair> > doubles;
	doubles.push_back(u);
	doubles.push_back(chi);
	iterate_singles(x,tau,doubles,LRCC2_);
	vecfuncT new_x=x.get_vecfunction();
	vecfuncT difference=sub(world,old_x,new_x);
	bool full_convergence=true;

	for(auto x : difference){
	  if(x.norm2() > parameters.dconv_6D*0.1) full_convergence=false;
	}
	if(full_convergence){
	  output("Change in Response-Singles is under: 0.1*dconv_6D="+std::to_string(parameters.dconv_6D*0.1));
	  output_section("Response of CC2 CONVERGED!!!");
	  timer_iter_all.info();
	  if(world.rank()==0) std::cout << "Current omega is " << x.omega << std::endl;
	  timer_iter_all.info();
	  break;
	}else output("Overall convergence not yet reached ... starting cycle again");
      }

      timer_iter_all.info();
    }

    time.info();
    return omega;
  }

  double
  CC2::solve_cc2(Pairs<CC_Pair> &doubles,CC_vecfunction &singles) {
    if(singles.size()==0) CCOPS.error("Forgot to initialize Singles");
    CCOPS.update_intermediates(singles);
    output_section("Beginn the CC2 Iterations");
    bool singles_converged=true;
    bool doubles_converged=true;

    // iterate singles for the first time (initialization)
    if(not parameters.restart)output_subsection("Initialize Singles with MP2 Doubles");
    else output_section("Initialize Singles from Restart");
    singles_converged=iterate_cc2_singles(doubles,singles);

    std::vector<double> current_energies=update_cc2_pair_energies(doubles,singles);

    double cc2_correlation_energy = get_correlation_energy(doubles);

    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      CC_Timer timer_iter_all(world,"Macroiteration " + stringify(iter));
      output_subsection("Macroiteration " + stringify(iter));
      CCOPS.print_memory_information(singles,doubles);

      CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Doubles");
      for(auto& tmp_pair : doubles.allpairs){
	update_constant_part_cc2_gs(tmp_pair.second,singles);
	bool pair_converged=iterate_pair(tmp_pair.second);
	if(not pair_converged) doubles_converged=false;
      }

      std::vector<double> updated_energies=update_cc2_pair_energies(doubles,singles);
      bool pair_energies_converged=check_energy_convergence(current_energies,updated_energies);
      current_energies=updated_energies;

      output("Pair Correlation Energies of Macro-Iteration " + stringify(iter));
      if(world.rank() == 0) std::cout << std::setprecision(parameters.output_prec) << current_energies << std::endl;
      timer_iter_doubles.info();

      output("Performance Overview of Iteration " + stringify(iter));
      if(world.rank() == 0) CCOPS.performance_S.info_last_iter();
      if(world.rank() == 0) CCOPS.performance_D.info_last_iter();

      output("\nPair Energies");
      if(world.rank() == 0) std::cout << std::setprecision(parameters.output_prec) << current_energies << std::endl;

      output("End performance Overview\n");

      const double new_correlation_energy= get_correlation_energy(doubles);
      const double delta_correlation = new_correlation_energy -cc2_correlation_energy;
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(parameters.output_prec) << "Overall Correlation_Energy\n new, old , diff\n";
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(parameters.output_prec) << new_correlation_energy << ", " << cc2_correlation_energy <<", " << delta_correlation <<"\n";
      bool energy_converged = true;
      if(fabs(delta_correlation>parameters.econv)) energy_converged = false;
      if(world.rank() == 0){
	output("End of Macroiteration " + stringify(iter));
	std::cout << "singles converged: " << singles_converged << std::endl;
	std::cout << "doubles converged: " << doubles_converged << std::endl;
	std::cout << "overall energy  converged: " << energy_converged << std::endl;
	std::cout << "all pair-energies converged:" << pair_energies_converged << std::endl;
      }
      cc2_correlation_energy = new_correlation_energy;

      print_results(doubles,singles);

      if(doubles_converged and energy_converged and pair_energies_converged){
	output("Doubles Converged (!)");
	output("Testing if singles do not change anymore");
	bool no_change_in_singles=iterate_cc2_singles(doubles,singles);
	if(no_change_in_singles){
	  output("");
	  output_section("CC2 CONVERGED!!!");
	  timer_iter_all.info();
	  print_results(doubles,singles);
	  break;
	}else output("Overall convergence not yet reached ... starting cycle again");
      }
      timer_iter_all.info();

    }
    CCOPS.plot(singles);
    return get_correlation_energy(doubles);
  }

  std::vector<double>
  CC2::update_cc2_pair_energies(Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    std::vector<double> omegas;
    double correlation_energy = 0.0;
    for(size_t i=parameters.freeze; i < mo.size(); i++){
      for(size_t j=i; j < mo.size(); j++){
	double tmp= CCOPS.compute_pair_energy(doubles(i,j),singles);//CCOPS.compute_cc2_pair_energy(doubles(i,j),singles(i),singles(j)); //CHANGED
	doubles(i,j).current_energy = tmp;
	omegas.push_back(tmp);
	correlation_energy += tmp;
	if(i!=j) correlation_energy +=tmp;
      }
    }
    if(world.rank() == 0){
      std::cout << "Updated CC2 pair energies:\n";
      for(auto x : omegas) std::cout << std::fixed << std::setprecision(parameters.output_prec) << x << std::endl;
      std::cout <<"Current Correlation Energy = " << std::fixed << std::setprecision(10) << correlation_energy << std::endl;
    }
    return omegas;
  }

  bool
  CC2::iterate_cc2_singles(const Pairs<CC_Pair> &doubles,CC_vecfunction &singles) {
    const std::vector<Pairs<CC_Pair>> u(1,doubles);
    CC_vecfunction empty(UNDEFINED);
    CCOPS.remove_stored_singles_potentials();
    return iterate_singles(singles,empty,u,CC2_);
  }
  
bool CC2::iterate_pair(CC_Pair &pair,const double omega) const {
    output_section("Iterate Pair " + pair.name());
    if(omega!=0 and pair.type!=EXCITED_STATE) CCOPS.warning("Iterate_pair: Omega is not zero but pair is part of the excited-state");

    const real_function_6d constant_part = pair.constant_term;

    output_subsection("Converge pair " + pair.name() + " on constant singles potential");

    CC_Timer make_BSH_time(world,"Make Detructive Greens Operator");
    double bsh_eps = CCOPS.get_epsilon(pair.i,pair.j)+omega;
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    make_BSH_time.info();

    NonlinearSolverND<6> solver(parameters.kain_subspace);
    solver.do_print=(world.rank() == 0);

    bool converged=false;
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      output_subsection("MP2-Microiteration");
      CC_Timer timer_mp2(world,"MP2-Microiteration of pair " + pair.name());

      CC_Timer timer_mp2_potential(world,"MP2-Potential of pair " + pair.name());
      real_function_6d mp2_potential=-2.0 * CCOPS.get_MP2_potential_residue(pair);
      mp2_potential.truncate().reduce_rank();
      timer_mp2_potential.info();

      CC_Timer timer_G(world,"Apply Greens Operator on MP2-Potential of pair " + pair.name());
      const real_function_6d GVmp2=G(mp2_potential);
      timer_G.info();

      CC_Timer timer_addup(world,"Add constant parts and update pair " + pair.name());
      real_function_6d unew=GVmp2 + constant_part;
      CCOPS.apply_Q12(unew,"new-pair-function");
      const real_function_6d residue=pair.function - unew;
      const double error=residue.norm2();
      if(parameters.kain){
	output("Update with KAIN");
	real_function_6d kain_update=copy(solver.update(unew,residue));
	CCOPS.apply_Q12(kain_update,"Kain-Update-Function");
	kain_update.truncate();
	kain_update.print_size("Kain-Update-Function");
	pair.function=copy(kain_update);
      }else{
	output("Update without KAIN");
	pair.function=unew;
      }

      const double old_energy=pair.current_energy;
      if(pair.type==GROUND_STATE)pair.current_energy=CCOPS.compute_mp2_pair_energy(pair);
      const double delta=pair.current_energy - old_energy;
      pair.current_error=error;
      pair.current_energy_difference=delta;

      timer_addup.info();

      output("\n\n Iteration " + stringify(iter) + " ended");
      pair.info();
      pair.store_pair(world);
      timer_mp2.info();
      if(fabs(error) < parameters.dconv_6D){
	output(pair.name() + " converged!");
	if(fabs(delta) < parameters.econv_pairs){
	  converged=true;
	  break;
	}else output("Energy not yet converged");
      }else output("Convergence for pair " + pair.name() + " not reached yet");
    }

    // store current functions for restart later
    pair.store_pair(world);

    return converged;
  }

  bool CC2::iterate_nonorthogonal_pair(CC_Pair &pair) {
    output("Iterate Nonorthogonal Pair " + pair.name());

    if(not parameters.restart) pair.function = real_factory_6d(world);
    output_subsection("Get Nonorthogonal Regularization Potential of Pair " + pair.name());
    CC_Timer timer_cc2_regular(world,"Get Nonorthogonal Regularization Potential of Pair " + pair.name());
    const real_function_6d regular_part=CCOPS.make_nonorthogonal_regularization_residue(pair.i,pair.j);
    regular_part.print_size("Regularization part of pair " + pair.name());
    timer_cc2_regular.info();



    real_function_6d coulomb_part=real_factory_6d(world);
    CC_Timer timer_cc2_coulomb(world,"Get Screened Coulomb Potentials");
    output("Increasing 6D thresh for screened Coulomb parts");
    FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);
    coulomb_part=CCOPS.make_nonorthogonal_mp2_coulomb_parts(pair.i,pair.j);
    coulomb_part.print_size("Screened Coulomb part for pair " + pair.name());
    FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
    timer_cc2_coulomb.info();


    real_function_6d constant_part=(regular_part + coulomb_part);
    constant_part = CCOPS.do_f12_projection(constant_part,CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j));
    constant_part = CCOPS.do_f12_projection(constant_part,CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j));
    // if(ctype == CISpD_)  pair.constant_term = copy(constant_part);
    constant_part.print_size("constant part of pair " + pair.name());
    pair.constant_term=copy(constant_part);

    output_subsection("Converge pair " + pair.name() + " on constant singles potential");

    CC_Timer make_BSH_time(world,"Make Detructive Greens Operator");
    double bsh_eps = CCOPS.get_epsilon(pair.i,pair.j);
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    make_BSH_time.info();

    NonlinearSolverND<6> solver(parameters.kain_subspace);
    solver.do_print=(world.rank() == 0);

    bool converged=false;
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      output_subsection("MP2-Iteration");
      CC_Timer timer_mp2(world,"MP2-Iteration of pair " + pair.name());

      CC_Timer timer_mp2_potential(world,"MP2-Potential of pair " + pair.name());
      real_function_6d mp2_potential=-2.0 * CCOPS.get_MP2_potential_residue(pair);
      mp2_potential.print_size("mp2-potential");
      mp2_potential.truncate().reduce_rank();
      timer_mp2_potential.info();

      CC_Timer timer_G(world,"Apply Greens Operator on MP2-Potential of pair " + pair.name());
      const real_function_6d GVmp2=G(mp2_potential);
      timer_G.info();

      CC_Timer timer_addup(world,"Add constant parts and update pair " + pair.name());
      real_function_6d unew=GVmp2 + constant_part;
      unew = CCOPS.do_f12_projection(unew,CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j));
      unew = CCOPS.do_f12_projection(unew,CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j));
      const real_function_6d residue=pair.function - unew;
      const double error=residue.norm2();
      if(parameters.kain){
	output("Update with KAIN");
	real_function_6d kain_update=copy(solver.update(unew,residue));
	kain_update.truncate();
	kain_update.print_size("Kain-Update-Function");
	pair.function=copy(kain_update);
      }else{
	output("Update without KAIN");
	pair.function=unew;
      }
      const double old_energy=pair.current_energy;
      pair.current_energy=CCOPS.compute_mp2_pair_energy(pair);
      const double delta=pair.current_energy - old_energy;
      pair.current_error=error;
      pair.current_energy_difference=delta;
      timer_addup.info();

      output("\n\n Iteration " + stringify(iter) + " ended");
      pair.info();

      timer_mp2.info();
      if(fabs(error) < parameters.dconv_6D){
	output(pair.name() + " converged!");
	if(fabs(delta) < parameters.econv){
	  converged=true;
	  break;
	}else output("Energy not yet converged");
      }else output("Convergence for pair " + pair.name() + " not reached yet");
    }

    // store current functions for restart later
    if(true){
      std::string msg="nono_mp2_";
      if(converged) msg+="converged_pair_";
      else msg+="not_converged_pair_";
      pair.store_pair(world,msg);
    }else{
      pair.store_pair(world,"current_pair_");
    }

    return converged;
  }



  CC_vecfunction
  CC2::initialize_cc2_singles() const {

    std::vector<CC_function> singles;
    for(size_t i=parameters.freeze;i<mo.size();i++){
      CC_function single_i;
      single_i.type=PARTICLE;
      single_i.i = i;
      real_function_3d tmpi = real_factory_3d(world);
      if(parameters.restart){
	if(CCOPS.load_function<double,3>(tmpi,single_i.name())) output("found " + single_i.name());
	else CCOPS.warning("Restart demanded but single " + single_i.name()+" not found! .... initialize as zero-function");
      }
      else output("Initialized ground state single " + single_i.name()+" as zero-function");
      single_i.function = copy(tmpi);
      singles.push_back(single_i);
    }

    return CC_vecfunction(singles,PARTICLE);
  }


  // Unnecessary function
  double
  CC2::compute_mp2_pair_energy(CC_Pair &u) const {
    return CCOPS.compute_mp2_pair_energy(u);
  }

  Pairs<CC_Pair> CC2::initialize_pairs(const pairtype type,const double omega)const{
    if(type==GROUND_STATE and omega!=0.0) CCOPS.warning("Initialize Pairs: Ground state demanded but omega is not zero");
    if(type==EXCITED_STATE and omega==0.0) CCOPS.warning("Initialize Pairs: Excited state demanded but omega is zero");
    Pairs<CC_Pair> pairs;
    for(size_t i=parameters.freeze;i<mo.size();i++){
      for(size_t j=i;j<mo.size();j++){
	CC_Pair u(i,j,type);
	u.function=real_factory_6d(world);
	u.current_energy =0.0;
	u.current_energy_difference = CC_Pair::uninitialized();
	if(parameters.restart){
	  output("Looking for Restartdata for pair " + u.name());
	  if(u.load_pair(world)){
	    u.function.print_size("loaded pair "+ u.name());
	    u.info();
	  }
	}
	if(u.type==GROUND_STATE){
	  const double ij_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.i],u.i,HOLE),CC_function(mo[u.j],u.j,HOLE));
	  const double ji_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.j],u.j,HOLE),CC_function(mo[u.i],u.i,HOLE));
	  u.constant_energy = 2.0*ij_gQf_ij - ji_gQf_ij;
	  u.current_energy = CCOPS.compute_mp2_pair_energy(u);
	}else if(u.type==EXCITED_STATE){
	  u.current_energy = omega;
	  u.constant_energy = omega;
	}

	pairs.insert(i,j,u);
      }
    }
    return pairs;
  }

  double
  CC2::get_correlation_energy(const Pairs<CC_Pair> &doubles) const {
    double result=0.0;
    output("Getting Correlation Energies:");
    for(const auto& tmp_pair : doubles.allpairs){
      const CC_Pair& pair=tmp_pair.second;
      const double omega=pair.current_energy;
      const bool symmetric=(pair.i == pair.j);
      result+=omega;
      if(not symmetric) result+=omega;     // off diagonal pairs count twice because only one of them was computed
      std::cout << std::fixed << std::setprecision(parameters.output_prec);
      if(world.rank() == 0) std::cout << "pair-correlation-energy of pair " << pair.name() << " = " << omega << std::endl;
    }
    std::cout << std::fixed << std::setprecision(parameters.output_prec);
    if(world.rank() == 0) std::cout << "Overall correlation energy = " << result << std::endl;
    return result;
  }

} /* namespace madness */
