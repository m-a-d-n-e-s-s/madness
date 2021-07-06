/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#include <chem/CC2.h>
#include <chem/commandlineparser.h>

namespace madness {

  /// solve the CC2 ground state equations, returns the correlation energy
  void
  CC2::solve() {
    if(parameters.test){
      CCOPS.test();
    }
    const CalcType ctype = parameters.calculation;

    if(ctype==CT_TDHF){
//      TDHF tdhf(world,nemo);
        commandlineparser parser;
        TDHF tdhf(world,parser);
        tdhf.set_reference(nemo);
        MADNESS_EXCEPTION("check if the previous 4 lines are OK",1);
      std::vector<CC_vecfunction> ccs;
      for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) ccs.push_back(tmp);
      }
      tdhf.solve_tdhf(ccs);
    }
    else if(ctype==CT_MP2){
      Pairs<CCPair> pairs;
      initialize_pairs(pairs,GROUND_STATE,CT_MP2,CC_vecfunction(PARTICLE),CC_vecfunction(RESPONSE),0);
      const double mp2_correlation_energy = solve_mp2(pairs);
      output.section(assign_name(ctype)+" Calculation Ended !");
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << " MP2 Correlation Energy =" << mp2_correlation_energy << "\n";
    }else if(ctype==CT_CC2){
      double mp2_correlation_energy=0.0;
      double cc2_correlation_energy=0.0;

      Pairs<CCPair> doubles;
      CC_vecfunction singles;

      // check if singles or/and doubles to restart are there
      initialize_singles(singles,PARTICLE);
      const bool load_doubles = initialize_pairs(doubles,GROUND_STATE,CT_CC2,singles,CC_vecfunction(RESPONSE),0);

      // nothing to restart -> make MP2
      if(not load_doubles){
	output("\n--- No Restartdata found: Solving MP2 as first guess");
	Pairs<CCPair> mp2_pairs;
	initialize_pairs(mp2_pairs,GROUND_STATE,CT_MP2,CC_vecfunction(PARTICLE),CC_vecfunction(RESPONSE),0);
	mp2_correlation_energy = solve_mp2(mp2_pairs);
	// use mp2 as cc2 guess
	for(auto& tmp:mp2_pairs.allpairs){
	  const size_t i=tmp.second.i;
	  const size_t j=tmp.second.j;
	  doubles(i,j).update_u(tmp.second.function());
	}

      }

      cc2_correlation_energy = solve_cc2(singles,doubles);

      output.section(assign_name(ctype)+" Calculation Ended !");
      if(world.rank()==0 and mp2_correlation_energy!=0.0) std::cout << std::fixed << std::setprecision(10) << " MP2 Correlation Energy =" << mp2_correlation_energy << "\n";
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << " CC2 Correlation Energy =" << cc2_correlation_energy << "\n";

    }else if(ctype==CT_LRCCS){
//      TDHF tdhf(world,nemo);
        commandlineparser parser;
        TDHF tdhf(world,parser);
        tdhf.set_reference(nemo);
        MADNESS_EXCEPTION("check if the previous 4 lines are OK",1);
      std::vector<CC_vecfunction> ccs;
      for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) ccs.push_back(tmp);
      }
      tdhf.solve_cis(ccs);

    }else if(ctype==CT_CISPD){
      CCTimer time(world,"whole CIS(D) Calculation");
      CCTimer time_mp2(world,"MP2 Calculation");
      Pairs<CCPair> mp2;
      initialize_pairs(mp2,GROUND_STATE,CT_MP2,CC_vecfunction(PARTICLE),CC_vecfunction(RESPONSE),0);
      const double mp2_correlation_energy = solve_mp2(mp2);
      output.section("MP2 Ground State Calculation Ended !");
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << " MP2 Correlation Energy =" << mp2_correlation_energy << "\n";
      time_mp2.info();

      std::vector<CC_vecfunction> vccs;
      // try to load cis functions
      std::vector<std::pair<double,double> > cispd_results;
      for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) vccs.push_back(tmp);
      }

      CCTimer time_ccs(world,"Time CCS");
      if(vccs.empty()) vccs = solve_ccs();
      time_ccs.info();

      CCTimer time_cispd(world,"Time CIS(D) Response");

      for(size_t k=0;k<parameters.excitations_.size();k++){

	CC_vecfunction& ccs = vccs[k];
	const size_t excitation = parameters.excitations_[k];
	CCTimer time_ex(world,"CIS(D) for Excitation " + std::to_string(int(excitation)));

	// check the convergence of the cis function (also needed to store the ccs potential) and to recalulate the excitation energy
	iterate_ccs_singles(ccs);

	Pairs<CCPair> cispd;
	initialize_pairs(cispd,EXCITED_STATE,CT_CISPD,CC_vecfunction(PARTICLE),ccs,excitation);

	const double ccs_omega = ccs.omega;
	const double cispd_omega =solve_cispd(cispd,mp2,ccs);

	cispd_results.push_back(std::make_pair(ccs_omega,cispd_omega));
	time_ex.info();
      }

      output.section("CIS(D) Calculation Ended");
      if(world.rank()==0){
	std::cout << std::fixed << std::setprecision(10) << "\n"
	    << "MP2 Correlation Energy: " << mp2_correlation_energy <<"\n\n";
      }
      for(size_t i=0;i<cispd_results.size();i++){
	if(world.rank()==0){
	  std::cout << std::fixed << std::setprecision(10) << "\n"
	      << "--------------------------------\n"
	      << "Excitation " << parameters.excitations_[i] << "\n"
	      << "CIS   =" << cispd_results[i].first << "\n"
	      << "CIS(D)=" << cispd_results[i].second << "\n"
	      << "Delta =" << cispd_results[i].second - cispd_results[i].first << "\n"
	      << "--------------------------------\n";
	}
      }
      time_cispd.info();
      time.info();

    } else if(ctype==CT_ADC2){
      // we will never need the GS singles, but we use the CC2 potential functions so we initialize all gs singles potentials to zero
      CCOPS.update_intermediates(CC_vecfunction(zero_functions<double,3>(world,CCOPS.get_active_mo_ket().size()),PARTICLE,parameters.freeze));
      output.section("ADC(2) Calculation");
      CCTimer time(world,"Whole ADC(2) Calculation");

      CCTimer time_gs(world,"ADC(2): MP2 Ground State");
      output.section("ADC(2): Calculating MP2 Ground State");

      Pairs<CCPair> mp2;
      initialize_pairs(mp2,GROUND_STATE,CT_MP2,CC_vecfunction(PARTICLE),CC_vecfunction(RESPONSE),0);
      const double mp2_correlation_energy = solve_mp2(mp2);
      output.section("MP2 Ground State Calculation Ended !");
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << " MP2 Correlation Energy =" << mp2_correlation_energy << "\n";
      time_gs.info();

      CCTimer time_cis(world,"ADC(2): CIS Guess");
      output.section("ADC(2): Calculating CIS Guess for Excited State");
      std::vector<CC_vecfunction> vccs;
      // try to load cis functions
      std::vector<std::pair<double,double> > cispd_results;
      for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) vccs.push_back(tmp);
      }
      if(vccs.empty()) vccs = solve_ccs(); // only solve if no CIS vectors where given
      time_cis.info();

      CCTimer time_ex(world,"ADC(2) Calculation");
      output.section("ADC(2): Calculating ADC(2) Correction to CIS");
      std::vector<std::vector<double> > adc2_results;
      for(size_t k=0;k<parameters.excitations_.size();k++){

	CC_vecfunction& ccs = vccs[k];
	const size_t excitation = parameters.excitations_[k];
	CCTimer time_ex(world,"ADC(2) for Excitation " + std::to_string(int(excitation)));

	// check the convergence of the cis function (also needed to store the ccs potential) and to recalulate the excitation energy
	CC_vecfunction dummy = ccs.copy();
	iterate_ccs_singles(dummy);
	ccs.omega = dummy.omega; // will be overwritten soon
	output("Changes not stored!");

	Pairs<CCPair> xpairs;
	const bool restart =initialize_pairs(xpairs,EXCITED_STATE,CT_ADC2,CC_vecfunction(PARTICLE),ccs,excitation);

	// if no restart: Calculate CIS(D) as first guess
	const double ccs_omega = ccs.omega;
	double cispd_omega = 0.0;
	if(not restart){
	  output.section("No Restart-Pairs found: Calculating CIS(D) as first Guess");
	  Pairs<CCPair> cispd;
	  initialize_pairs(cispd,EXCITED_STATE,CT_CISPD,CC_vecfunction(PARTICLE),ccs,excitation);
	  cispd_omega =solve_cispd(cispd,mp2,ccs);
	  for(auto& tmp:cispd.allpairs){
	    const size_t i=tmp.first.first;
	    const size_t j=tmp.first.second;
	    xpairs(i,j).update_u(cispd(i,j).function());
	  }
	}

	iterate_adc2_singles(mp2,ccs,xpairs);
	for(size_t iter=0;iter<10;iter++){
	  bool dconv = iterate_adc2_pairs(xpairs,ccs);
	  bool sconv = iterate_adc2_singles(mp2,ccs,xpairs);
	  if(sconv and dconv){
	    output("ADC(2) Converged");
	    break;
	  } else output("Not yet converged");
	}

	output.section("ADC(2) For Excitation " + std::to_string(int(excitation)) + " ended");
	const double adc2_omega = ccs.omega;
	std::vector<double> resulti;
	resulti.push_back(ccs_omega);
	resulti.push_back(cispd_omega);
	resulti.push_back(adc2_omega);
	adc2_results.push_back(resulti);
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< std::setfill(' ') << std::setw(12) << "CIS" << std::setw(12) << "CIS(D)" << std::setw(12) << "ADC(2)" << "\n"
	<< ccs_omega << ", " << cispd_omega << ", " << adc2_omega << "\n";


	time_ex.info();
      }

      output.section("ADC(2) Ended!");
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
      << std::setfill(' ') << std::setw(12) << "CIS" << std::setw(12) << "CIS(D)" << std::setw(12) << "ADC(2)" << "\n";
      for(size_t i=0;i<adc2_results.size();i++){
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< adc2_results[i][0] << ", " << adc2_results[i][1] << ", " << adc2_results[i][2] << "\n";
      }


      time_ex.info();

      time.info();
    } else if(ctype==CT_LRCC2){
      CCTimer time(world,"Whole LRCC2 Calculation");

      std::vector<std::pair<std::string,double> > results;
      std::vector<std::pair<std::string,std::pair<double,double> > > timings;

      double mp2_correlation_energy = 0.0;
      double cc2_correlation_energy = 0.0;

      CC_vecfunction cc2_s;
      initialize_singles(cc2_s,PARTICLE);
      Pairs<CCPair> cc2_d;
      bool found_cc2d = initialize_pairs(cc2_d,GROUND_STATE,CT_CC2,cc2_s,CC_vecfunction(RESPONSE));

      if(not found_cc2d){
	Pairs<CCPair> mp2;
	initialize_pairs(mp2,GROUND_STATE,CT_MP2,CC_vecfunction(PARTICLE),CC_vecfunction(RESPONSE));
	CCTimer time_mp2(world,"MP2 Calculation");
	mp2_correlation_energy = solve_mp2(mp2);
	results.push_back(std::make_pair("MP2 correlation energy",mp2_correlation_energy));
	timings.push_back(std::make_pair("MP2",time_mp2.current_time(true)));
	// use mp2 as cc2 guess
	for(auto tmp:mp2.allpairs){
	  const size_t i=tmp.second.i;
	  const size_t j=tmp.second.j;
	  cc2_d(i,j).update_u(tmp.second.function());
	}
      }

      CCTimer time_cc2(world,"CC2 Calculation");
      cc2_correlation_energy = solve_cc2(cc2_s,cc2_d);
      results.push_back(std::make_pair("CC2 correlation energy",cc2_correlation_energy));
      timings.push_back(std::make_pair("CC2",time_cc2.current_time(true)));

      output.section("Ground State Calculation of LRCC2 ended");
      for(const auto& res:results){
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< res.first << "=" << res.second << "\n";
      }
      for(const auto& time:timings){
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< time.first << ": " << time.second.first << " (Wall), " << time.second.second << " (CPU)" << "\n";
      }

      std::vector<CC_vecfunction> vccs;
      // try to load cis functions
      std::vector<std::pair<double,double> > cispd_results;
      for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) vccs.push_back(tmp);
      }


      if(vccs.empty()){
	CCTimer time_ccs(world,"CCS Calculation");
	vccs = solve_ccs();
	timings.push_back(std::make_pair("CCS",time_ccs.current_time(true)));
      }

      std::vector<std::pair<std::string,std::pair<double,double> > >results_ex;
      for(size_t xxx=0;xxx<vccs.size();xxx++){
	const size_t excitation = parameters.excitations_[xxx];
	CCTimer time_ex(world,"LRCC2 Calculation for Excitation " + std::to_string(int(excitation)));
	CC_vecfunction lrcc2_s = vccs[xxx];
	// needed to assign an omega
	const vector_real_function_3d backup = copy(world,lrcc2_s.get_vecfunction());
	CC_vecfunction test(backup,RESPONSE,parameters.freeze);
	test.excitation = lrcc2_s.excitation;
	iterate_ccs_singles(test);
	lrcc2_s.omega=test.omega;
	output("CCS Iteration: Changes are not applied (just omega)!");


	Pairs<CCPair> lrcc2_d;
	bool found_lrcc2d = initialize_pairs(lrcc2_d,EXCITED_STATE,CT_LRCC2,cc2_s,lrcc2_s,excitation);

	if(found_lrcc2d) iterate_lrcc2_singles(cc2_s,cc2_d,lrcc2_s,lrcc2_d);
	else iterate_ccs_singles(lrcc2_s);
	const double omega_cis = lrcc2_s.omega;

	for(size_t iter=0;iter<parameters.iter_max;iter++){
	  output.section("Macroiteration " + std::to_string(int(iter)) + " of LRCC2");
	  bool dconv = iterate_lrcc2_pairs(cc2_s,cc2_d,lrcc2_s,lrcc2_d);
	  bool sconv = iterate_lrcc2_singles(cc2_s,cc2_d,lrcc2_s,lrcc2_d);
	  if(dconv and sconv) break;
	}
	const double omega_cc2 = lrcc2_s.omega;
	const std::string msg = "Excitation " + std::to_string(int(excitation));
	results_ex.push_back(std::make_pair(msg,std::make_pair(omega_cis,omega_cc2)));
	timings.push_back(std::make_pair(msg,time_ex.current_time(true)));

      }

      timings.push_back(std::make_pair("Whole LRCC2",time.current_time(true)));
      output.section("LRCC2 Finished");
      output("Ground State Results:");
      for(const auto& res:results){
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< res.first << "=" << res.second << "\n";
      }
      output("Response Results:");
      for(const auto& res:results_ex){
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10)
	<< res.first << ": "<< res.second.first << " (CIS)*, " << res.second.second << " (CC2)\n";
      }
      if(world.rank()==0) std::cout << "*only if CIS vectors where given in the beginning (not for CC2 restart)\n";
      output("\nTimings");
      for(const auto& time:timings){
	if(world.rank()==0) std::cout << std::scientific << std::setprecision(2)
	<< std::setfill(' ') <<std::setw(15) << time.first
	<< ": " << time.second.first << " (Wall), " << time.second.second << " (CPU)" << "\n";
      }




    } else MADNESS_EXCEPTION(("Unknown Calculation Type: " + assign_name(ctype)).c_str(),1);

  }

  // Solve the CCS equations for the ground state (debug potential and check HF convergence)
  std::vector<CC_vecfunction> CC2::solve_ccs() {
    output.section("SOLVE CCS");
//    TDHF tdhf(world,nemo);
      commandlineparser parser;
      TDHF tdhf(world,parser);
      tdhf.set_reference(nemo);
      MADNESS_EXCEPTION("check if the previous 4 lines are OK",1);
    std::vector<CC_vecfunction> excitations;
    for(size_t k=0;k<parameters.excitations_.size();k++){
	CC_vecfunction tmp;
	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
	if(found) excitations.push_back(tmp);
    }
    tdhf.solve_cis(excitations);

    // return only those functions which are demanded
    std::vector<CC_vecfunction> result;
    for(const auto& x:parameters.excitations_){
      if(excitations.size()-1<x) MADNESS_EXCEPTION("Not Enough CIS Vectors to solve for the demanded CC2 vector",1);
      excitations[x].excitation = x;
      result.push_back(excitations[x]);
    }
    return result;
  }

  double
  CC2::solve_mp2(Pairs<CCPair> &doubles) {
    output.section("Solve MP2");
    double omega =0.0;
    Pairs<double> pair_energies;
    for(auto& tmp_pair : doubles.allpairs){
      MADNESS_ASSERT(tmp_pair.second.type==GROUND_STATE);
      MADNESS_ASSERT(tmp_pair.second.ctype==CT_MP2);

      if(parameters.no_compute_mp2) output("Found no_compute_mp2 keyword");
      else {
	update_constant_part_mp2(tmp_pair.second);
	iterate_pair(tmp_pair.second);
      }
      save(tmp_pair.second.function(),tmp_pair.second.name());
      const double pair_energy = CCOPS.compute_pair_correlation_energy(tmp_pair.second);
      pair_energies.insert(tmp_pair.second.i,tmp_pair.second.j,pair_energy);
      omega += pair_energy;

    }
    if(world.rank()==0) std::cout << "\nMP2 Pair Correlation Energies:\n";
    for(auto& a : pair_energies.allpairs){
      if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << "omega_" << a.first.first << a.first.second << "=" << a.second << "\n";
    }
    if(world.rank()==0) std::cout << "sum     =" << omega << "\n";
    return omega;
  }

  double
  CC2::solve_cispd(Pairs<CCPair> &cispd,const Pairs<CCPair> &mp2,const CC_vecfunction & ccs) {
    output.section("Solve CIS(D) for CIS Excitation energy " + std::to_string(double(ccs.omega)));
    MADNESS_ASSERT(ccs.type==RESPONSE);
    CCOPS.update_intermediates(ccs);

    for(auto &pairs:cispd.allpairs){
      CCPair& pair = pairs.second;
      pair.bsh_eps = CCOPS.get_epsilon(pair.i,pair.j) + ccs.omega;
      if(size_t(parameters.only_pair.first) == pair.i and size_t(parameters.only_pair.second) == pair.j){
	output("Found only_pair exception");
	update_constant_part_cispd(ccs,pair);
	iterate_pair(pair,ccs);
      }
      else if(parameters.no_compute_cispd) output("Found no_compute_cispd key");
      else{
	update_constant_part_cispd(ccs,pair);
	iterate_pair(pair,ccs);
      }
      // test consitency of the two approaches
      if(parameters.debug and parameters.thresh_6D > 1.e-4) CCOPS.test_pair_consistency(pair.functions[0],pair.i,pair.j,ccs);
    }

    const double diff = CCOPS.compute_cispd_energy(ccs,mp2,cispd);
    CC_vecfunction empty(zero_functions<double,3>(world,ccs.size()),PARTICLE,parameters.freeze);
    const double omega_cc2 = CCOPS.compute_cc2_excitation_energy(empty,ccs,mp2,cispd);
    output.section("CIS(D) Calculation for CIS Excitation " + std::to_string(double(ccs.omega)) + " ended");
    if(world.rank()==0){
      std::cout << std::fixed << std::setprecision(10)
      << "CIS   =" << ccs.omega << "\n"
      << "CIS(D)=" << ccs.omega+diff << "\n"
      << "Diff  =" << diff
      << "\nomega_cc2 =" << omega_cc2 << "\n\n\n";
    }

    return ccs.omega+diff;
  }

  bool
  CC2::iterate_adc2_pairs(Pairs<CCPair> &cispd,const CC_vecfunction & ccs) {
    output.section("Solve ADC(2) for Excitation energy " + std::to_string(double(ccs.omega)));
    MADNESS_ASSERT(ccs.type==RESPONSE);
    CCOPS.update_intermediates(ccs);

    bool conv = true;
    for(auto &pairs:cispd.allpairs){
      CCPair& pair = pairs.second;
      pair.bsh_eps = CCOPS.get_epsilon(pair.i,pair.j) + ccs.omega;
      update_constant_part_adc2(ccs,pair);
      conv=iterate_pair(pair,ccs);
    }

    return conv;
  }
  bool
  CC2::iterate_lrcc2_pairs(const CC_vecfunction & cc2_s,const Pairs<CCPair>& cc2_d,const CC_vecfunction lrcc2_s,Pairs<CCPair> &lrcc2_d) {
    output.section("Solve LRCC2 for Excitation energy " + std::to_string(double(lrcc2_s.omega)));
    MADNESS_ASSERT(lrcc2_s.type==RESPONSE);
    CCOPS.update_intermediates(lrcc2_s);

    bool conv = true;
    for(auto &tmp:lrcc2_d.allpairs){
      CCPair& pair = tmp.second;
      const size_t i=pair.i;
      const size_t j=pair.j;
      // check if singles have significantly changed
      if(lrcc2_s(i).current_error<0.1*parameters.thresh_6D and lrcc2_s(j).current_error<0.1*parameters.thresh_6D) output("Skipping Pair Iteration, No significant Change in Singles");
      else{
	pair.bsh_eps = CCOPS.get_epsilon(pair.i,pair.j) + lrcc2_s.omega;
	update_constant_part_lrcc2(pair,cc2_s,lrcc2_s);
	conv=iterate_pair(pair,lrcc2_s);
      }
    }

    return conv;
  }


  double
  CC2::solve_cc2(CC_vecfunction &singles,Pairs<CCPair> &doubles) {

    output.section("Solving CC2 Ground State");

    MADNESS_ASSERT(singles.type==PARTICLE);
    CCOPS.update_intermediates(singles);
    output.section("Solve CC2 Ground State");
    CCTimer time(world,"CC2 Ground State");

    double omega = CCOPS.compute_cc2_correlation_energy(singles,doubles);
    if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";

    if(not parameters.no_compute_cc2){
      // first singles iteration
      output.section("Initialize Singles to the Doubles");
      iterate_cc2_singles(singles,doubles);
      // update correlation energy
      omega = CCOPS.compute_cc2_correlation_energy(singles,doubles);

      for(size_t iter=0;iter<parameters.iter_max;iter++){
	CCTimer time_miter(world,"Macroiteration " + std::to_string(int(iter)) + " of CC2");
	output.section("Macroiteration " + std::to_string(int(iter)) + " of CC2");

	// iterate doubles
	bool doubles_converged = true;
	for(auto& pairs: doubles.allpairs){
	  CCPair& pair = pairs.second;
	  update_constant_part_cc2_gs(singles,pair);
	  bool pair_converged=iterate_pair(pair,singles);
	  save(pair.function(),pair.name());
	  if(not pair_converged) doubles_converged=false;
	}

	// new omega
	omega = CCOPS.compute_cc2_correlation_energy(singles,doubles);

	// check if singles converged
	const bool singles_converged = iterate_cc2_singles(singles,doubles);

	// check if energy converged
	const double omega_new = CCOPS.compute_cc2_correlation_energy(singles,doubles);
	const double delta = omega_new - omega;
	const bool omega_converged(delta<parameters.econv);
	omega = omega_new;
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";
	if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << "Difference                  = " << delta << "\n";

	if(doubles_converged and singles_converged and omega_converged) break;

	time_miter.info();
      }
      omega= CCOPS.compute_cc2_correlation_energy(singles,doubles);
      output.section("CC2 Iterations Eneded");
    }else{
      output.section("Found no_compute_cc2 Key: Reiterating Singles to check convergence");
      // need the singles potential for the constant part of LRCC2 so we recompute it (also good to check if it is converged)
      bool sconv = iterate_cc2_singles(singles,doubles);
      if(not sconv ) output.warning("Singles not Converged");
    }

    if(world.rank()==0) std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";
    time.info();
    return omega;

  }


  bool CC2::iterate_pair(CCPair &pair,const CC_vecfunction &singles) const {
    output.section("Iterate Pair " + pair.name());
    if(pair.ctype==CT_CC2) MADNESS_ASSERT(singles.type==PARTICLE);
    if(pair.ctype==CT_CISPD) MADNESS_ASSERT(singles.type==RESPONSE);
    if(pair.ctype==CT_MP2) MADNESS_ASSERT(singles.get_vecfunction().empty());
    if(pair.ctype==CT_ADC2)MADNESS_ASSERT(singles.type==RESPONSE);

    real_function_6d constant_part = pair.constant_part;
    constant_part.truncate().reduce_rank();
    pair.function().truncate().reduce_rank();

    output.subsection("Converge pair " + pair.name() + " on constant singles potential");

    double bsh_eps = pair.bsh_eps; //CCOPS.get_epsilon(pair.i,pair.j)+omega;
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;

    NonlinearSolverND<6> solver(parameters.kain_subspace);
    solver.do_print=(world.rank() == 0);

    bool converged=false;

    double omega =0.0;
    if(pair.type==GROUND_STATE) omega= CCOPS.compute_pair_correlation_energy(pair,singles);
    if(pair.type==EXCITED_STATE) omega= CCOPS.compute_excited_pair_energy(pair,singles);

    if(world.rank()==0) std::cout << "Correlation Energy of Pair " << pair.name() <<" =" << std::fixed << std::setprecision(10) << omega << "\n";

    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      output.subsection(assign_name(pair.ctype)+"-Microiteration");
      CCTimer timer_mp2(world,"MP2-Microiteration of pair " + pair.name());



      CCTimer timer_mp2_potential(world,"MP2-Potential of pair " + pair.name());
      real_function_6d mp2_potential=-2.0 * CCOPS.fock_residue_6d(pair);
      if(parameters.debug) mp2_potential.print_size(assign_name(pair.ctype)+" Potential");
      mp2_potential.truncate().reduce_rank();
      timer_mp2_potential.info(true,mp2_potential.norm2());

      CCTimer timer_G(world,"Apply Greens Operator on MP2-Potential of pair " + pair.name());
      const real_function_6d GVmp2=G(mp2_potential);
      timer_G.info(true,GVmp2.norm2());

      CCTimer timer_addup(world,"Add constant parts and update pair " + pair.name());
      real_function_6d unew=GVmp2 + constant_part;
      unew.print_size("unew");
      unew = CCOPS.apply_Q12t(unew,CCOPS.mo_ket());
      unew.print_size("Q12unew");
      //unew.truncate().reduce_rank(); // already done in Q12 application at the end
      if(parameters.debug)unew.print_size("truncated-unew");
      const real_function_6d residue=pair.function() - unew;
      const double error=residue.norm2();
      if(parameters.kain){
	output("Update with KAIN");
	real_function_6d kain_update=copy(solver.update(pair.function(),residue));
	kain_update=CCOPS.apply_Q12t(kain_update,CCOPS.mo_ket());
	kain_update.truncate().reduce_rank();
	kain_update.print_size("Kain-Update-Function");
	pair.update_u(copy(kain_update));
      }else{
	output("Update without KAIN");
	pair.update_u(unew);
      }

      timer_addup.info(true,pair.function().norm2());

      double omega_new = 0.0;
      double delta =0.0;
      if(pair.type==GROUND_STATE) omega_new = CCOPS.compute_pair_correlation_energy(pair,singles);
      else if(pair.type==EXCITED_STATE) omega_new = CCOPS.compute_excited_pair_energy(pair,singles);
      delta = omega - omega_new;

      const double current_norm = pair.function().norm2();

      omega = omega_new;
      if(world.rank()==0){
	  std::cout<< std::fixed
	  << std::setw(50) << std::setfill('#')
	  << "\n" << "Iteration " << iter <<" of pair " << pair.name()
	  << std::setprecision(4) << "||u|| = " << current_norm
	  << "\n"<< std::setprecision(10) << "error = " << error <<"\nomega = " << omega << "\ndelta = " << delta << "\n"
	  << std::setw(50) << std::setfill('#') << "\n";
      }



      output("\n--Iteration " + stringify(iter) + " ended--");
      save(pair.function(),pair.name());
      timer_mp2.info();
      if(fabs(error) < parameters.dconv_6D){
	output(pair.name() + " converged!");
	if(fabs(delta) < parameters.econv_pairs){
	  converged=true;
	  break;
	}else output("Energy not yet converged");
      }else output("Convergence for pair " + pair.name() + " not reached yet");
    }

    return converged;
  }


  bool
  CC2::initialize_singles(CC_vecfunction &singles,const FuncType type,const int ex) const {
    MADNESS_ASSERT(singles.size()==0);
    bool restarted = false;

    std::vector<CCFunction> vs;
    for(size_t i=parameters.freeze;i<CCOPS.mo_ket().size();i++){
      CCFunction single_i;
      single_i.type=type;
      single_i.i = i;
      std::string name;
      if(ex<0) name = single_i.name();
      else name = std::to_string(ex)+"_"+single_i.name();
      real_function_3d tmpi = real_factory_3d(world);
      const bool found=CCOPS.load_function<double,3>(tmpi,name);
      if(found) restarted = true;
      else output("Initialized " + single_i.name()+" of type " + assign_name(type) +" as zero-function");
      single_i.function = copy(tmpi);
      vs.push_back(single_i);
    }

    singles = CC_vecfunction(vs,type);
    if(type==RESPONSE) singles.excitation=ex;

    return restarted;
  }

  bool
  CC2::initialize_pairs(Pairs<CCPair> &pairs,const CCState ftype, const CalcType ctype, const CC_vecfunction &tau, const CC_vecfunction &x, const size_t excitation)const{
    MADNESS_ASSERT(tau.type==PARTICLE);
    MADNESS_ASSERT(x.type==RESPONSE);
    MADNESS_ASSERT(pairs.empty());
    output("Initialize " + assign_name(ctype) + " Pairs for " + assign_name(ftype));

    bool restarted = false;

    for(size_t i=parameters.freeze;i<CCOPS.mo_ket().size();i++){
      for(size_t j=i;j<CCOPS.mo_ket().size();j++){

	std::string name=CCPair(i,j,ftype,ctype).name();
	if(ftype==GROUND_STATE){
	  real_function_6d utmp = real_factory_6d(world);
	  const bool found =CCOPS.load_function(utmp,name);
	  if(found) restarted =true; // if a single pair was found then the calculation is not from scratch
	  real_function_6d const_part;
	  CCOPS.load_function(const_part,name+"_const");
	  CCPair tmp = CCOPS.make_pair_gs(utmp,tau,i,j);
	  tmp.constant_part = const_part;
	  pairs.insert(i,j,tmp);

	  //const double omega = CCOPS.compute_pair_correlation_energy(tmp);
	  //if(world.rank()==0) std::cout << "initialized pair " << tmp.name() << " with correlation energy=" << std::fixed << std::setprecision(10) << omega << "\n";

	}else if(ftype==EXCITED_STATE){
	  name = std::to_string(int(excitation))+"_"+name;
	  real_function_6d utmp = real_factory_6d(world);
	  const bool found =CCOPS.load_function(utmp,name);
	  if(found) restarted =true;
	  real_function_6d const_part;
	  CCOPS.load_function(const_part,name+"_const");
	  CCPair tmp = CCOPS.make_pair_ex(utmp,tau,x,i,j,ctype);
	  tmp.excitation=excitation;
	  tmp.constant_part = const_part;
	  pairs.insert(i,j,tmp);
	}else error("Unknown pairtype");
      }
    }
    return restarted;
  }

  void CC2::update_reg_residues_gs(const CC_vecfunction& singles,Pairs<CCPair> &doubles)const{
    CCTimer time(world,"Updated Regularization Residues of the Ground State");
    MADNESS_ASSERT(singles.type==PARTICLE);
    Pairs<CCPair> updated_pairs;
    //    output("Correlation energy with old pairs");
    //    CCOPS.compute_cc2_correlation_energy(singles,doubles);
    for(auto& tmp:doubles.allpairs){
      MADNESS_ASSERT(tmp.second.type==GROUND_STATE);
      CCPair& pair = tmp.second;
      const size_t i=pair.i;
      const size_t j=pair.j;
      const CCPair updated_pair = CCOPS.make_pair_gs(pair.function(),singles,i,j);
      updated_pairs.insert(i,j,updated_pair);
    }
    //    output("Correlation energy with updated pairs");
    //    CCOPS.compute_cc2_correlation_energy(singles,updated_pairs);
    doubles.swap(updated_pairs);
    //    output("Correlation energy with swapped pairs");
    //    CCOPS.compute_cc2_correlation_energy(singles,updated_pairs);
    time.info();
  }
  void CC2::update_reg_residues_ex(const CC_vecfunction& singles,const CC_vecfunction& response,Pairs<CCPair> &doubles)const{
    CCTimer time(world,"Updated Regularization Residues of the Excited State");
    MADNESS_ASSERT(singles.type==PARTICLE);
    MADNESS_ASSERT(response.type==RESPONSE);
    Pairs<CCPair> updated_pairs;
    for(auto& tmp:doubles.allpairs){
      MADNESS_ASSERT(tmp.second.type==EXCITED_STATE);
      CCPair& pair = tmp.second;
      const size_t i=pair.i;
      const size_t j=pair.j;
      CCPair updated_pair = CCOPS.make_pair_ex(pair.function(),singles,response,i,j,pair.ctype);
      updated_pairs.insert(i,j,updated_pair);
    }
    doubles.swap(updated_pairs);
    time.info();
  }


} /* namespace madness */
