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
    calctype type = parameters.calculation;
    if(type==CCS_response_)solve_ccs();
    if(type==MP2_){
      Pairs<CC_Pair> pairs = initialize_pairs(GROUND_STATE);
      const double mp2_correlation_energy = solve_mp2(pairs);
      output_section("MP2 Ended");
      print_results(pairs,initialize_cc2_singles());
      if(world.rank()==0) std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << mp2_correlation_energy << "\n";
    }
    if(type==experimental_){
      Pairs<CC_Pair> pairs = initialize_pairs(GROUND_STATE);
      const double mp2_correlation_energy = solve_mp2_nonorthogonal(pairs);
      output_section("Nonorthogonal-MP2 Ended");
      print_results(pairs,initialize_cc2_singles());
      if(world.rank()==0) std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << mp2_correlation_energy << "\n";
    }
    if(type==CC2_){
      CC_Timer time_cc2(world,"CC2-Overall-Time");
      Pairs<CC_Pair> pairs = initialize_pairs(GROUND_STATE);
      double mp2_correlation_energy = 0.0;
      if(not parameters.no_compute) mp2_correlation_energy = solve_mp2(pairs);
      CC_vecfunction singles=initialize_cc2_singles();
      print_results(pairs,singles);
      CCOPS.update_intermediates(singles);
      const double cc2_correlation_energy = solve_cc2(pairs,singles);
      output_section("CC2 Ended");
      print_results(pairs,singles);
      if(world.rank()==0){
	std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << mp2_correlation_energy << "\n";
	std::cout << "CC2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << cc2_correlation_energy << "\n";
      }
      time_cc2.info();
    }
    if(type==CISpD_){
      Pairs<CC_Pair> mp2_pairs = initialize_pairs(GROUND_STATE);
      output("Test ijgu");
      const CC_Pair test = mp2_pairs.allpairs.find(std::make_pair(0,0))->second;
      CCOPS.make_ijgu(CCOPS.mo_ket(0),CCOPS.mo_ket(0),test);
      output("Test finished");
      const double mp2_correlation_energy = solve_mp2(mp2_pairs);
      output_section("MP2 Ended");
      if(world.rank()==0) std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << mp2_correlation_energy << "\n";

      std::vector<std::pair<CC_vecfunction,double> > cis_results = solve_ccs();
      std::vector<std::pair<double,double>> result;
      size_t i=0;
      for(const auto& cistmp:cis_results){
	const double cis_omega = cistmp.second;
	const CC_vecfunction cis_singles = cistmp.first;
	cis_singles.print_size();
	Pairs<CC_Pair> cispd_pairs = initialize_pairs(EXCITED_STATE,cis_omega);
	const double cispd_omega = solve_cispd(cispd_pairs,mp2_pairs,cis_singles,cis_omega);

	if(world.rank()==0) std::cout << "CIS(D) for Excitation " << i << "\n"
	    << std::fixed << std::setprecision(parameters.output_prec)
	<< cis_omega <<" (CIS), " << cispd_omega << " (CIS(D))\n";
	result.push_back(std::make_pair(cis_omega,cispd_omega));
	i++;
      }

      output_section("CIS(D) Ended:");
      i=0;
      for(const auto & res:result){
	if(world.rank()==0) std::cout << "CIS(D) for Excitation " << i << "\n"
	    << std::fixed << std::setprecision(parameters.output_prec)
	<< res.first <<" (CIS), " << res.second << " (CIS(D))\n";
	i++;
      }

    }
    if(type==CC2_response_){
      CC_Timer time_cc2(world,"CC2-Overall-Time");
      CC_Timer time_cc2_gs(world,"CC2-Ground-State-Time");
      Pairs<CC_Pair> pairs = initialize_pairs(GROUND_STATE);
      double mp2_correlation_energy = 0.0;
      if(not parameters.no_compute) mp2_correlation_energy = solve_mp2(pairs);
      CC_vecfunction singles=initialize_cc2_singles();
      print_results(pairs,singles);
      CCOPS.update_intermediates(singles);
      double cc2_correlation_energy =0.0;
      if(not parameters.no_compute) cc2_correlation_energy = solve_cc2(pairs,singles);
      else cc2_correlation_energy = get_correlation_energy(pairs);
      output_section("CC2 Ground-State Calculation Ended");
      print_results(pairs,singles);
      if(world.rank()==0){
	std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << mp2_correlation_energy << "\n";
	std::cout << "CC2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << cc2_correlation_energy << "\n";
      }
      time_cc2_gs.info();

      std::vector<std::pair<CC_vecfunction,double> > cis_results = solve_ccs();
      size_t i=0;
      std::vector<std::pair<CC_vecfunction,double> > cc2_results;
      for(const auto& cistmp:cis_results){
	const double omega = cistmp.second;
	output_section("Solving Excitation " + std::to_string(i) + " with CIS excitation energy " + std::to_string(omega));
	CC_vecfunction x(cistmp.first);
	x.omega = omega;
	Pairs<CC_Pair> chi = initialize_pairs(EXCITED_STATE,omega);
	CCOPS.update_response_intermediates(x);

	// reiterate ccs/cis (make shure it converged and make shure the singles potential is stored (needed for fock application)
	std::vector<Pairs<CC_Pair>> empty_pair_vector;
	CC_vecfunction empty_vector;
	if(not iterate_singles(x,empty_vector,empty_pair_vector, CCS_response_)) CCOPS.warning("CCS/CIS not converged!");

	// make CIS(D) as first guess for doubles
	const double omega_cispd = solve_cispd(chi,pairs,x,x.omega);
	singles.omega = omega_cispd;

	// reiterate cc2 singles (make shure they converged and that the singles potential is stored
	if(iterate_cc2_singles(pairs,singles)) output("CC2 Singles Converged and Potential is stored");
	else CCOPS.warning("CC2 Singles are not fully converged --> doubles maybe also not");

	const double cc2_omega = solve_cc2_response(singles,pairs,x,chi);
	cc2_results.push_back(std::make_pair(x,cc2_omega));
	if(world.rank()==0) std::cout << " Excitation " << i << "\n" << "Excitation Energy (CIS) " << omega << "\n" << "Excitation Energy (CC2)" << cc2_omega << "\n";
	i++;
      }
      output_section("CC2 Response Ended");
      if(world.rank()==0) std::cout << "\n\n" << "Results of CC2 Response: CIS | CC2 Excitation Energies:\n";
      for(size_t i=0;i<cc2_results.size();i++){
	if(world.rank()==0) std::cout << "Excitation " << i << ": " << cis_results[i].second << " | " << cc2_results[i].second << "\n";
      }



      solve_ccs();

      time_cc2.info();
    }
  }

  // Solve the CCS equations for the ground state (debug potential and check HF convergence)
  std::vector<std::pair<CC_vecfunction,double> > CC2::solve_ccs() {
    output_section("SOLVE CCS");
    // this is that the output of CIS is not screwed
    std:cout << std::setw(5) << std::setfill(' ') << std::endl;
    TDA CCS_Solver(world,nemo,nemo.get_calc()->amo,"input");
    {
      xfunctionsT guess;
      CCS_Solver.solve_guess(guess);
      guess.clear();
    }
    //    {
    //      xfunctionsT solve = CCS_Solver.get_converged_xfunctions();
    //      CCS_Solver.solve(solve);
    //      solve.clear();
    //    }
    //    {
    //      xfunctionsT final = CCS_Solver.get_converged_xfunctions();
    //      CCS_Solver.solve_sequential(final);
    //      final.clear();
    //    }
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
	bool pair_converged=iterate_pair(tmp_pair.second,CC_vecfunction());
	if(not pair_converged) mp2_converged=false;
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
	tmp_pair.second.ij_gQf_ij = CCOPS.make_ijgfxy(tmp_pair.second.i,tmp_pair.second.j,CCOPS.mo_ket(tmp_pair.second.i).function,CCOPS.mo_ket(tmp_pair.second.j).function);
	tmp_pair.second.ji_gQf_ij = CCOPS.make_ijgfxy(tmp_pair.second.j,tmp_pair.second.i,CCOPS.mo_ket(tmp_pair.second.i).function,CCOPS.mo_ket(tmp_pair.second.j).function);
	if(world.rank()==0){
	  std::cout << "Pair " << tmp_pair.second.name() << "\n";
	  std::cout << "ij_gf_ij=" << tmp_pair.second.ij_gQf_ij << "\n";
	  std::cout << "ji_gf_ij=" << tmp_pair.second.ji_gQf_ij << "\n";
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
    for(auto& tmp_pair : doubles.allpairs){
      tmp_pair.second.current_energy = cis_omega;
      bool pair_converged=iterate_pair(tmp_pair.second,CC_vecfunction(),cis_singles,CISpD_);
      if(not pair_converged) cispd_converged=false;
    }

    output_section("Calculating the CIS(D) Excitation Energies");
    const double result = CCOPS.compute_cispd_energy_correction(cis_singles,mp2_pairs,doubles);
    if(world.rank()==0){
      std::cout <<"Excitation Energy CIS   =" <<  cis_omega << "\n";
      std::cout <<"CIS(D) Correction       =" << result << "\n";
      std::cout <<"Excitation Energy CIS(D)=" << result+ cis_omega << "\n";
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
      singles_converged=iterate_singles(x,tau,doubles,CC2_response_);
      const double new_omega = x.omega;
      const double diff = new_omega - omega;
      bool energy_converged = (fabs(diff)<parameters.econv);

      CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Response of Doubles");
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
      for(auto& tmp_pair : chi.allpairs){
	bool pair_converged=iterate_pair(tmp_pair.second,tau,x,CC2_response_);
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
	iterate_singles(x,tau,doubles,CC2_response_);
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
    std::vector<double> current_energies;
    for(const auto &tmp:doubles.allpairs) current_energies.push_back(tmp.second.current_energy);
    double cc2_correlation_energy = get_correlation_energy(doubles);
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      CC_Timer timer_iter_all(world,"Macroiteration " + stringify(iter));
      output_subsection("Macroiteration " + stringify(iter));
      CCOPS.print_memory_information(singles,doubles);

      // Iterate singles
      CC_vecfunction old_singles(singles);
      for(auto& tmp : singles.functions) old_singles(tmp.first).function=copy(tmp.second.function);
      singles_converged=iterate_cc2_singles(doubles,singles);

      CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Doubles");
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
      for(auto& tmp_pair : doubles.allpairs){
	bool pair_converged=iterate_pair(tmp_pair.second,singles);
	if(not pair_converged) doubles_converged=false;
      }
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
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

      CCOPS.remove_stored_singles_potentials();
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

      if(singles_converged and doubles_converged and energy_converged and pair_energies_converged){
	output("Singles and Doubles Converged (!)");
	output("Testing if singles do not change anymore");
	vecfuncT old_singles=singles.get_vecfunction();
	iterate_cc2_singles(doubles,singles);
	vecfuncT new_singles=singles.get_vecfunction();
	vecfuncT difference=sub(world,old_singles,new_singles);
	bool full_convergence=true;

	for(auto x : difference){
	  if(x.norm2() > parameters.dconv_6D*0.1) full_convergence=false;
	}
	if(full_convergence){
	  output("Change in Singles is under: 0.1*dconv_6D="+std::to_string(parameters.dconv_6D*0.1));
	  output_section("CC2 CONVERGED!!!");
	  timer_iter_all.info();
	  //print_results(doubles,singles);
	  break;
	}else output("Overall convergence not yet reached ... starting cycle again");
      }
      timer_iter_all.info();

    }
    CCOPS.plot(singles);
    return get_correlation_energy(doubles);
  }

  std::vector<double>
  CC2::update_cc2_pair_energies(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    std::vector<double> omegas;
    double correlation_energy = 0.0;
    for(size_t i=parameters.freeze; i < mo.size(); i++){
      for(size_t j=i; j < mo.size(); j++){
	double tmp=CCOPS.compute_cc2_pair_energy(doubles(i,j),singles(i),singles(j));
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
    return iterate_singles(singles,empty,u,CC2_);

    //    if(singles.functions.size() != active_mo.size())
    //      MADNESS_EXCEPTION(("Wrong size of singles at beginning of iterations " + stringify(singles.functions.size())).c_str(),1);
    //    output_subsection("Iterate CC2 Singles");
    //    CC_Timer timer_potential(world,"CC2 Singles Potential");
    //    vecfuncT potential=CCOPS.get_CC2_singles_potential(singles,doubles);
    //    timer_potential.info();
    //
    //    output_subsection("Apply the Green's Operator");
    //    CC_Timer timer_G(world,"Apply the Green's Operator");
    //    vecfuncT G_potential=zero_functions<double, 3>(world,potential.size());
    //    scale(world,potential,-2.0);
    //    for(size_t i=0; i < potential.size(); i++){
    //      double epsi=CCOPS.get_orbital_energies()[i + parameters.freeze];
    //      output("Make Greens Operator for single " + stringify(i + parameters.freeze));
    //      real_convolution_3d G=BSHOperator<3>(world,sqrt(-2.0 * epsi),parameters.lo,parameters.thresh_bsh_3D);
    //      real_function_3d tmp=(G(potential[i])).truncate();
    //      G_potential[i]=tmp;
    //    }
    //    G_potential=CCOPS.apply_Q(G_potential,"G_potential");
    //    timer_G.info();
    //
    //    std::vector<double> errors;
    //    bool converged=true;
    //    for(size_t i=0; i < potential.size(); i++){
    //      MADNESS_ASSERT(singles(i + parameters.freeze).i == i + parameters.freeze);
    //      if(world.rank() == 0) std::cout << "|| |tau" + stringify(i + parameters.freeze) + ">|| =" << G_potential[i].norm2() << std::endl;
    //      real_function_3d residue=singles(i + parameters.freeze).function - G_potential[i];
    //      double error=residue.norm2();
    //      errors.push_back(error);
    //      if(world.rank() == 0) std::cout << "|| residue" + stringify(i + parameters.freeze) + ">|| =" << error << std::endl;
    //      CC_function new_single(G_potential[i],singles(i + parameters.freeze).i,PARTICLE);
    //      new_single.current_error=error;
    //      singles(i + parameters.freeze)=new_single;
    //      if(fabs(error) > parameters.dconv_3D) converged=false;
    //    }
    //    if(singles.functions.size() != active_mo.size())
    //      MADNESS_EXCEPTION(("Wrong size of singles at the end of the iteration " + stringify(singles.functions.size())).c_str(),1);
    //    if(converged) output("singles converged");
    //    else output("No convergence in singles");
    //    return converged;
  }

  // the constant part of CC2 will be calculated and stored in the pair, after that the pair is iterated with the MP2 alg. till it converges
  // Note: if the singles are initialized to 0 then this is just MP2
  bool CC2::iterate_pair(CC_Pair &pair,const CC_vecfunction &singles) const {
    if(singles.size()==0) return iterate_pair(pair,singles,CC_vecfunction(PARTICLE),MP2_);
    else if(singles.type==PARTICLE) return iterate_pair(pair,singles,CC_vecfunction(RESPONSE),CC2_);
    else if(singles.type==RESPONSE) return iterate_pair(pair,CC_vecfunction(PARTICLE),singles,CISpD_);
    else error("Calles iterate Pair with singles of undefined type");
    MADNESS_EXCEPTION("cannot get here but need to make compiler happy",1);
  }
  bool CC2::iterate_pair(CC_Pair &pair,const CC_vecfunction &singles, const CC_vecfunction &response_singles,const calctype ctype) const {
    output("Iterate " + assign_name(ctype) + " Pair " + pair.name());

    double omega = 0.0;
    if(ctype==CC2_response_ or CISpD_){
      omega = response_singles.omega;
      pair.current_energy = omega;
    }

    // check if the constant part has to be recalculated
    bool recalc_const=true;
    if(ctype == MP2_ and parameters.restart) recalc_const=false;
    if(ctype == CC2_){
      if(singles(pair.i).current_error < parameters.dconv_3D and singles(pair.j).current_error < parameters.dconv_3D) recalc_const=false;
    }
    if(ctype == CC2_response_){
      if(response_singles(pair.i).current_error < parameters.dconv_3D and response_singles(pair.j).current_error < parameters.dconv_3D) recalc_const=false;
    }
    if(ctype == CISpD_) recalc_const=false;


    if(not pair.constant_term.impl_initialized()){
      output("semi-constant-term was not stored -> recalculation");
      recalc_const=true;
    }


    real_function_6d regular_part;
    if(recalc_const){
      output_subsection("Get" + assign_name(ctype)+" Regularization Potential of Pair " + pair.name());
      CC_Timer timer_cc2_regular(world,"Get Regularization Potential of Pair " + pair.name());
      if(ctype == CC2_)        regular_part=CCOPS.make_regularization_residue(singles(pair.i),singles(pair.j),ctype,0.0);           //   CCOPS.make_cc2_residue_sepparated(singles(pair.i),singles(pair.j));
      else if(ctype == MP2_)   regular_part=CCOPS.make_regularization_residue(CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j),ctype,0.0);//   CCOPS.make_cc2_residue_sepparated(CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j));
      else if(ctype == CISpD_) regular_part=CCOPS.make_regularization_residue(response_singles(pair.i),response_singles(pair.j),ctype,omega);
      else if(ctype == CC2_response_){
	regular_part=CCOPS.make_response_regularization_residue(CCOPS.make_t_intermediate(singles(pair.i)),response_singles(pair.i),CCOPS.make_t_intermediate(singles(pair.j)),response_singles(pair.j),omega);
      }
      regular_part.print_size("Regularization part of pair " + pair.name());
      timer_cc2_regular.info();
      pair.constant_term=copy(regular_part);
    }else{
      if(ctype==CC2_){
	output("Recalculation of Regularized-Potential of Pair " + pair.name() + " is not neccesary since the involved singles have not changed");
	output("Last change in " + singles(pair.i).name() + "=" + stringify(singles(pair.i).current_error));
	output("Last change in " + singles(pair.j).name() + "=" + stringify(singles(pair.j).current_error));
      }else if(ctype==CC2_response_){
	output("Recalculation of Regularized-Potential of Pair " + pair.name() + " is not neccesary since the involved singles have not changed");
	output("Last change in " + response_singles(pair.i).name() + "=" + stringify(response_singles(pair.i).current_error));
	output("Last change in " + response_singles(pair.j).name() + "=" + stringify(response_singles(pair.j).current_error));
      }
      else output("Loading constant_term from saved pair");
      regular_part=pair.constant_term;
    }

    real_function_6d coulomb_part=real_factory_6d(world);
    if(ctype == CC2_ or ctype==CISpD_ or ctype==CC2_response_){
      CC_Timer timer_cc2_coulomb(world,"Get Screened Coulomb Potentials of singles");
      output("Increasing 6D thresh for screened Coulomb parts of singles");
      FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);
      if(ctype==CC2_) coulomb_part=CCOPS.make_cc2_coulomb_parts(CCOPS.make_t_intermediate(singles(pair.i)),CCOPS.make_t_intermediate(singles(pair.j)),singles);
      else if(ctype==CISpD_) coulomb_part=CCOPS.make_cispd_coulomb_parts(CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j),response_singles);
      else if(ctype==CC2_response_){
	bool symmetric = (pair.i==pair.j);
	const CC_function ti = CCOPS.make_t_intermediate(singles(pair.i));
	const CC_function tj = CCOPS.make_t_intermediate(singles(pair.j));
	const CC_function xi = response_singles(pair.i);
	const CC_function xj = response_singles(pair.j);

	FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);
	real_function_6d part1;
	real_function_6d part2;

	const real_function_6d xt_part = -1.0*CCOPS.make_G_P_g_xy(singles,xi,tj,omega) - CCOPS.make_G_P_g_xy(singles,tj,xi,omega) + CCOPS.make_G_P1P2_g_xy(singles,singles,xi,tj,omega);
	real_function_6d tx_part;
	if(symmetric) tx_part = CCOPS.swap_particles(xt_part);
	else tx_part = -1.0*CCOPS.make_G_P_g_xy(singles,ti,xj,omega) - CCOPS.make_G_P_g_xy(singles,xj,ti,omega) + CCOPS.make_G_P1P2_g_xy(singles,singles,ti,xj,omega);
	part1 = xt_part + tx_part;

	if(symmetric){
	  const real_function_6d tmp1 = CCOPS.make_G_P_g_xy(response_singles,ti,tj,omega);
	  const real_function_6d tmp2 = CCOPS.make_G_P1P2_g_xy(response_singles,singles,ti,tj,omega);
	  part2 = -1.0*tmp1 - CCOPS.swap_particles(tmp1) + tmp2 + CCOPS.swap_particles(tmp2);
	}else{
	  part2 = -1.0*CCOPS.make_G_P_g_xy(response_singles,ti,tj,omega)
		        		  -  CCOPS.make_G_P_g_xy(response_singles,tj,ti,omega)
					  +  CCOPS.make_G_P1P2_g_xy(response_singles,singles,ti,tj,omega)
					  +  CCOPS.make_G_P1P2_g_xy(singles,response_singles,ti,tj,omega);
	}


	coulomb_part = part1 + part2;
	FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      }
      coulomb_part.print_size("Coulomb part of Singles for pair " + pair.name());
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      timer_cc2_coulomb.info();
    }

    real_function_6d constant_part=(regular_part + coulomb_part);
    CCOPS.apply_Q12(constant_part,"constant_part+screened_coulomb_part");
    // if(ctype == CISpD_)  pair.constant_term = copy(constant_part);
    constant_part.print_size("semi-constant-cc2-part of pair " + pair.name());

    output_subsection("Converge pair " + pair.name() + " on constant singles potential");

    CC_Timer make_BSH_time(world,"Make Detructive Greens Operator");
    double bsh_eps = CCOPS.get_epsilon(pair.i,pair.j);
    if(ctype==CISpD_) bsh_eps = CCOPS.get_epsilon(pair.i,pair.j) + omega;
    if(ctype==CC2_response_) bsh_eps = CCOPS.get_epsilon(pair.i,pair.j) + omega;
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    make_BSH_time.info();

    NonlinearSolverND<6> solver(parameters.kain_subspace);
    solver.do_print=(world.rank() == 0);

    bool converged=false;
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      if(ctype==MP2_)output_subsection("MP2-Microiteration");
      if(ctype==CC2_)output_subsection("MP2-Microiteration with Frozen CC2-Singles");
      if(ctype==CC2_response_)output_subsection("CC2-Response: MP2-Microiteration with Frozen CC2-Singles");
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
      if(ctype == CC2_) pair.current_energy=CCOPS.compute_cc2_pair_energy(pair,singles(pair.i),singles(pair.j));
      else if(ctype == MP2_) pair.current_energy=CCOPS.compute_mp2_pair_energy(pair);
      const double delta=pair.current_energy - old_energy;
      pair.current_error=error;
      pair.current_energy_difference=delta;
      timer_addup.info();

      output("\n\n Iteration " + stringify(iter) + " ended");
      pair.info();

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
    if(ctype == MP2_){
      std::string msg="mp2_";
      if(converged) msg+="converged_pair_";
      else msg+="not_converged_pair_";
      pair.store_pair(world,msg);
    }else{
      pair.store_pair(world,"current_pair_");
    }

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
	  u.ij_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.i],u.i,HOLE),CC_function(mo[u.j],u.j,HOLE));
	  u.ji_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.j],u.j,HOLE),CC_function(mo[u.i],u.i,HOLE));
	  u.epsilon = CCOPS.get_epsilon(u.i,u.j);
	  u.current_energy = CCOPS.compute_mp2_pair_energy(u);
	}else if(u.type==EXCITED_STATE){
	  u.epsilon = CCOPS.get_epsilon(u.i,u.j)+omega;
	  u.current_energy = omega;
	}

	pairs.insert(i,j,u);
      }
    }
    return pairs;
  }

  /// Initialize an electron pair
  void
  CC2::initialize_electron_pair(CC_Pair &u) const {
    output_subsection("Initialize Electron Pair |u" + stringify(u.i) + stringify(u.j) + ">");
    CC_Timer timer_integrals(world,"Make constant energy Integrals");
    u.ij_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.i],u.i,HOLE),CC_function(mo[u.j],u.j,HOLE));     // the u.i u.j have the right numbers (noo freeze problems should occur)
    u.ji_gQf_ij=CCOPS.make_ijgQfxy(u.i,u.j,CC_function(mo[u.j],u.j,HOLE),CC_function(mo[u.i],u.i,HOLE));
    timer_integrals.info();

    double epsij=CCOPS.get_epsilon(u.i,u.j);
    u.epsilon=epsij;

    u.function=real_factory_6d(world);
    u.constant_term=real_factory_6d(world);
    u.current_energy=0.0;
    u.current_error=0.0;

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
