/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#include "CC2.h"
namespace madness {

/// solve the CC2 ground state equations, returns the correlation energy
  double
  CC2::solve() {
    // Check if HF is converged
    //if(parameters.test) decompose_constant_part();
    if(parameters.ccs) solve_CCS();
    if(parameters.test) CCOPS.test_singles_potential();
    // Initialize the Pair functions (uij, i>=j)
    if(parameters.restart) output_section("Initialize Electron Pairs: Loading Stored Pairs");
    else output_section("Initialize Electron Pairs: First guess will be the constant Term of MP2");
    CC_Timer timer_init(world,"Initialization of all pairs");
    Pairs<CC_Pair> pairs;
    double mp2_energy=0.0;
    for(size_t i=parameters.freeze; i < mo.size(); i++){
      for(size_t j=i; j < mo.size(); j++){
	CC_Pair u(i,j);
	if(parameters.restart == true){
	  if(u.load_pair(world)){
	    u.function.print_size("loaded pair u" + stringify(i) + stringify(j));
	    output("...Found saved pair\n\n");
	    u.current_energy = CC_Pair::uninitialized();
	    u.current_error = CC_Pair::uninitialized();
	    u.current_energy_difference = CC_Pair::uninitialized();
	    mp2_energy=CCOPS.compute_mp2_pair_energy(u);
	    u.current_energy = mp2_energy;
	    output("Current MP2 Energy of the Pair is: " + stringify(mp2_energy));
	    u.info();
	  }else
	  MADNESS_EXCEPTION(("No Restartdata found for pair " + stringify(i) + stringify(j)).c_str(),1);
	}else initialize_electron_pair(u);
	pairs.insert(i,j,u);
      }
    }
    timer_init.info();

    output_section("Solve the CC2 equations");

    CC_vecfunction singles=initialize_cc2_singles();
    CCOPS.update_intermediates(singles);

    double correlation_energy_mp2=0.0;
    double correlation_energy_cc2=0.0;
    if(not parameters.restart) correlation_energy_mp2=solve_mp2(pairs);
    else if(parameters.mp2) correlation_energy_mp2=solve_mp2(pairs);

    if(not parameters.mp2_only) correlation_energy_cc2=solve_cc2(pairs,singles);
    output("Solving of CC2 ended at " + stringify(wall_time()) + "s (wall), " + stringify(cpu_time()) + "s (cpu)");

    output_section("Solve CC2 ended");
    if(world.rank() == 0)
      std::cout << "MP2 Correlation Energy is: " << std::fixed << std::setprecision(parameters.output_prec) << correlation_energy_mp2 << "\nCC2 Correlation Energy is: " << correlation_energy_cc2
	  << std::endl;
    output_section("Nothing more Implemented right now");
    return correlation_energy_cc2;
  }

// Solve the CCS equations for the ground state (debug potential and check HF convergence)
/// \todo Matt added "return false" at the end; otherwise, the end was reached without a return statement. please verify!
  bool
  CC2::solve_CCS() {
    output_section("SOLVE CCS");
    // since the symmetry should not change use the projected aos from moldft as guess
    // calculate HF energy = sum_i 2.0*\epsilon_i + \sum_ij 2.0*<ij|g|ij> - <ij|g|ji>
    double HF_energy=nemo.get_calc()->molecule.nuclear_repulsion_energy();
    if(world.rank() == 0) std::cout << "Nuclear repulsion is " << HF_energy << std::endl;
    double Corr_energy=0.0;
    for(size_t i=0; i < active_mo.size(); i++){
      HF_energy+=2.0 * CCOPS.get_orbital_energies()[i];
      for(size_t j=0; j < active_mo.size(); j++){
	HF_energy-=2.0 * CCOPS.make_integral(i,j,CC_function(active_mo[i],i,HOLE),CC_function(active_mo[j],j,HOLE));
	HF_energy+=CCOPS.make_integral(j,i,CC_function(active_mo[i],i,HOLE),CC_function(active_mo[j],j,HOLE));
      }
    }
    if(world.rank() == 0) std::cout << "HARTREE-FOCK ENERGY IS: " << HF_energy << std::endl;
    real_function_3d guessi=real_factory_3d(world);
    //for(size_t i=0;i<nemo.get_calc()->ao.size();i++) guessi += nemo.get_calc()->ao[i];
    //CCOPS.Q(guessi);
    vecfuncT guess(active_mo.size(),guessi);
    CC_vecfunction singles(guess,PARTICLE,parameters.freeze);

    std::vector<double> omega;
    for(size_t iter=0; iter < 30; iter++){
      output_subsection("Iterate CCS: Iteration " + stringify(iter + 1));
      CCOPS.update_intermediates(singles);
      CC_Timer timer_potential(world,"CCS Potential");
      vecfuncT potential=CCOPS.get_CCS_potential(singles);
      timer_potential.info();

      output_subsection("Apply the Green's Operator");
      CC_Timer timer_G(world,"Apply the Green's Operator");
      vecfuncT G_potential=zero_functions<double, 3>(world,potential.size());
      scale(world,potential,-2.0);
      for(size_t i=0; i < potential.size(); i++){
	double epsi=CCOPS.get_orbital_energies()[i];
	real_convolution_3d G=BSHOperator<3>(world,sqrt(-2.0 * epsi),parameters.lo,parameters.thresh_bsh_3D);
	real_function_3d tmp=(G(potential[i])).truncate();
	CCOPS.Q(tmp);
	G_potential[i]=tmp;
      }
      timer_G.info();

      std::vector<double> errors;
      bool converged=true;
      for(size_t i=0; i < potential.size(); i++){
	if(world.rank() == 0) std::cout << "|| |tau" + stringify(i) + ">|| =" << G_potential[i].norm2() << std::endl;
	real_function_3d residue=singles(i).function - G_potential[i];
	double error=residue.norm2();
	errors.push_back(error);
	if(world.rank() == 0) std::cout << "|| residue" + stringify(i) + ">|| =" << error << std::endl;
	CCOPS.Q(G_potential[i]);
	real_function_3d tau_before=singles(i).function;
	CC_function new_single_i(G_potential[i],singles(i).i,PARTICLE);
	singles(i)=new_single_i;
	real_function_3d tau_after=singles(i).function;
	double difference=(tau_before - tau_after).norm2();
	double difference2=(tau_after - G_potential[i]).norm2();
	if(world.rank() == 0) std::cout << "tau replacement debug (difference should be NOT zero): difference=" << difference << std::endl;
	if(world.rank() == 0) std::cout << "tau replacement debug (difference should be zero): difference=" << difference2 << std::endl;
	if(fabs(error) > parameters.dconv_3D) converged=false;
	omega.push_back(CCOPS.compute_ccs_correlation_energy(singles(i),singles(i)));
      }
      // print out the norms
      output("Performance Overview of Iteration " + stringify(iter));
      if(world.rank() == 0) CCOPS.performance_S.info_last_iter();
      if(world.rank() == 0) CCOPS.performance_D.info_last_iter();
      output("\nNorm of Singles\n");
      for(auto x : singles.functions)
	x.second.function.print_size("|tau_" + stringify(x.first) + ">");
      output("End performance Overview\n");

      output("Current CCS Correlation energies (Diagonal Part)");
      if(world.rank() == 0) std::cout << omega << std::endl;
      output("CCS Norms");
      for(size_t i=0; i < singles.size(); i++){
	double norm=singles(i).function.norm2();
	if(world.rank() == 0) std::cout << norm << std::endl;
      }
      output("Current <ti|F|ti> values");
      const CC_vecfunction t=CCOPS.make_t_intermediate(singles);
      const real_function_3d R2=nemo.nuclear_correlation->square();
      for(const auto itmp : t.functions){
	const real_function_3d brai=itmp.second.function * R2;
	const real_function_3d Fketi=CCOPS.apply_F(CC_function(itmp.second.function,itmp.first,UNDEFINED));
	const double omegai=brai.inner(Fketi);
	std::cout << "<" << itmp.second.name() << "|F|" << itmp.second.name() << ">=" << omegai << std::endl;
      }
      if(converged) break;

    }
    Corr_energy=omega.back();
    if(world.rank() == 0) std::cout << "Correlation Energy Convergence:\n" << omega << std::endl;
    if(world.rank() == 0)
      std::cout << "CCS finished\n" << "Hartree_Fock energy is " << HF_energy << "\nCorrelation Energy is " << Corr_energy << "\nTotel Energy is " << HF_energy + Corr_energy << std::endl;

    return false;
  }

  double
  CC2::solve_MP2_alternative(Pairs<CC_Pair> &pairs) const {
    for(auto pairtmp : pairs.allpairs){
      CC_Pair & pair=pairtmp.second;
      const size_t i=pair.i;
      const size_t j=pair.j;
      pair.epsilon=CCOPS.get_epsilon(i,j);
      real_function_6d u_final=copy(pair.constant_term);
      double omega=CCOPS.compute_mp2_pair_energy(pair);
      std::string gvstring="GVc";

      CC_Timer timer_makeG(world,"Construct destructive Greens Operator");
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * pair.epsilon),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive()=true;
      timer_makeG.info();

      for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
	if(0.1 * pair.function.norm2() < FunctionDefaults<6>::get_thresh()){
	  output("Doing adaptive thresh increase");
	  FunctionDefaults<6>::set_thresh(0.1 * pair.function.norm2());
	  output("6D Thresh is now " + stringify(FunctionDefaults<6>::get_thresh()));
	}
	output("Now doing " + gvstring);
	CC_Timer timer_potential(world,"MP2:Potential");
	real_function_6d potential=CCOPS.get_MP2_potential_residue(pair);
	potential.scale(-2.0);
	timer_potential.info();

	CC_Timer timer_applyG(world,"Applying Greens Operator");
	real_function_6d GVpsi=G(potential);
	timer_applyG.info();
	CCOPS.apply_Q12(GVpsi);
	u_final+=GVpsi;
	double delta=(2.0 * CCOPS.make_ijgu(pair.i,pair.j,GVpsi) - CCOPS.make_ijgu(pair.j,pair.i,GVpsi));
	omega+=delta;
	pair.function=GVpsi;
	pair.current_error=GVpsi.norm2();

	CC_Timer timer_omega(world,"Compute current pair correlation energy");
	{
	  CC_Pair wrapper(u_final,i,j);
	  wrapper.ij_gQf_ij=pair.ij_gQf_ij;
	  wrapper.ji_gQf_ij=pair.ji_gQf_ij;
	  wrapper.constant_term=pair.constant_term;
	  pair.current_energy=CCOPS.compute_mp2_pair_energy(wrapper);
	}
	timer_omega.info();

	output("\n\n");
	pair.info();
	output("alternative omega = " + stringify(omega));
	output("delta = " + stringify(delta));
	u_final.print_size("full u-function");
	output("\n\nIteration " + stringify(iter) + " ended\n\n");

	if(fabs(pair.current_error) < parameters.dconv_6D){
	  output("Pair " + pair.name() + " converged!");
	  if(fabs(delta) < parameters.econv){
	    pair.function=u_final;
	    break;
	  }else{
	    output("Energy did not converge yet");
	  }
	}
	gvstring="GV" + gvstring;
      }

    }
    double mp2_correlation_energy=0.0;
    for(auto pairtmp : pairs.allpairs){
      mp2_correlation_energy+=pairtmp.second.current_energy;
    }
    return mp2_correlation_energy;
  }

  double
  CC2::solve_uncoupled_mp2(Pairs<CC_Pair> &pairs) const {
    // Loop over all Pair functions uij (i=<j)
    std::vector<double> pair_energies;
    for(auto& utmp : pairs.allpairs){
      CC_Pair& pair=utmp.second;
      output_subsection("Solving uncoupled MP2 equations for pair " + pair.name());

      output_subsection("Setup the BSH Operator");
      CC_Timer timer_bsh_setup(world,"Setup the BSH-Operator for the pair function");
      real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * CCOPS.get_epsilon(pair.i,pair.j)),parameters.lo,parameters.thresh_bsh_6D);
      G.destructive_=false;
      output("Constructed Green Operator is destructive ? : " + stringify(G.destructive()));
      output("eps in Green's Operator is " + stringify(CCOPS.get_epsilon(pair.i,pair.j)));
      timer_bsh_setup.info();

      NonlinearSolverND<6> solver(parameters.kain_subspace);
      solver.do_print=(world.rank() == 0);

      double current_energy=CCOPS.compute_mp2_pair_energy(pair);

      // Beginn the iterations
      for(size_t iter=0; iter < 30; iter++){
	CC_Timer timer_iteration(world,"Iteration " + stringify(iter));

	double current_error=99.9;
	// Compute the non constant part of the MP2 equations which is the regularized 6D Fock Residue
	//and apply the G Operator G[(2J-K(R)+Un)|uij>]

	output_subsection("Calculate MP2 Residue");
	CC_Timer timer_mp2_residue(world,"\n\nCalculate MP2 Residue (2J-K(R)+Un)|uij>\n\n");
	real_function_6d mp2_residue=CCOPS.get_MP2_potential_residue(pair).truncate().reduce_rank();
	mp2_residue.print_size("Vpsi");
	timer_mp2_residue.info();

	output_subsection("Apply the Green's Operator");
	CC_Timer timer_apply_bsh(world,"\n\nApply BSH Operator to MP2 Residue\n\n");
	mp2_residue.scale(-2.0);
	real_function_6d Gresidue=G(mp2_residue);
	Gresidue.print_size("G(J+U-K)|u>");
	timer_apply_bsh.info();

	// Add the constant part and the residue and make the new u function
	// |u_new> = (constant_part + Gresidue)
	output_subsection("Add the Constant Term");
	CC_Timer timer_addition(world,"\n\nAdd the constant_term and the MP2 Residue\n\n");
	real_function_6d unew=(pair.constant_term + Gresidue);
	unew.print_size("unew");
	CCOPS.apply_Q12(unew);
	unew.print_size("Q12(unew)");
	timer_addition.info();
	// Get the error
	CC_Timer timer_make_bsh_residue(world,"\n\nMake the BSH-Residue\n\n");
	real_function_6d bsh_residue=(pair.function - unew);
	timer_make_bsh_residue.info();
	current_error=bsh_residue.norm2();
	// update the pair function
	real_function_6d updated_function=unew;
	if(parameters.kain) solver.update(unew,bsh_residue);
	pair.function=updated_function;
	pair.function.truncate();
	pair.store_pair(world);
	// evaluate the current mp2 energy
	double new_energy=compute_mp2_pair_energy(pair);
	double delta=new_energy - current_energy;
	current_energy=new_energy;
	output("End of Iteration " + stringify(iter) + "at time: " + stringify(wall_time()));
	output("Norm of BSH Residue: " + stringify(current_error));
	output("MP2 Energy: New, Old, Difference : " + stringify(new_energy) + ", " + stringify(current_energy) + ", " + stringify(delta));

	output_subsection("End of Iteration " + stringify(iter));
	if(world.rank() == 0){
	  std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << "current correlation energy:" << std::fixed << std::setprecision(parameters.output_prec) << new_energy << std::endl;
	  std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << "previous correlation energy:" << current_energy << std::endl;
	  std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << "correlation energy difference:" << delta << std::endl;
	  std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << "current wavefunction error:" << current_error << std::endl;
	  std::cout << std::setw(10) << std::setfill(' ') << std::setw(50) << "wavefunction norm:" << pair.function.norm2() << std::endl;
	}

	timer_iteration.info();
	current_energy=new_energy;
	if(current_error < parameters.dconv_6D){
	  output("Wavefunction convergence fullfilled");
	  if(fabs(delta) < parameters.econv){
	    output("Energy converged");
	    pair_energies.push_back(current_energy);
	    pair.store_pair(world,"converged_mp2_");
	    break;
	  }
	}
      }

    }
    output_section("All Pair Energies Converged");
    output("Converged Pair Energies are:");
    if(world.rank() == 0){
      for(auto x : pair_energies)
	std::cout << std::setprecision(parameters.output_prec) << std::fixed << x << std::endl;
    }
    double correlation_energy=0.0;
    for(auto x : pair_energies)
      correlation_energy+=x;
    output("Correlation Energy is: " + stringify(correlation_energy));
    return correlation_energy;
  }

  double
  CC2::solve_mp2(Pairs<CC_Pair> &doubles) {
    output_section("Solve MP2");
    bool mp2_converged=true;
    for(auto& tmp_pair : doubles.allpairs){
      bool pair_converged=iterate_pair(tmp_pair.second);
      if(not pair_converged) mp2_converged=false;
    }
    print_results(doubles,initialize_cc2_singles());
    return get_correlation_energy(doubles);
  }

  double
  CC2::solve_cc2(Pairs<CC_Pair> &doubles,CC_vecfunction &singles) {
    if(singles.size() == 0){
      output_section("Initialize CC2 with Zero functions");
      CC_Timer init_singles(world,"Initialize CC2 Singles");
      singles=initialize_cc2_singles();
      if(singles.size() != active_mo.size())
      MADNESS_EXCEPTION(("Singles have wrong size after initialization: " + stringify(singles.size())).c_str(),1);
      CCOPS.update_intermediates(singles);     // remeber that singles are zero right now
    }

    output_section("Beginn the CC2 Iterations");
    bool singles_converged=true;
    bool doubles_converged=true;
    std::vector<double> current_energies=update_cc2_pair_energies(doubles,singles);
    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      CC_Timer timer_iter_all(world,"Macroiteration " + stringify(iter));
      CCOPS.update_intermediates(singles);
      output_subsection("Macroiteration " + stringify(iter));
      CCOPS.print_memory_information(singles,doubles);

      // Iterate singles
      CCOPS.check_stored_singles_potentials();
      CC_vecfunction old_singles(singles);
      for(auto& tmp : singles.functions)
	old_singles(tmp.first).function=copy(tmp.second.function);
      for(size_t mis=0; mis < parameters.iter_max_3D; mis++){
	CC_Timer timer_iter_singles(world,"Iteration " + stringify(iter) + " Singles" + " Microiteration #" + stringify(mis));
	singles_converged=iterate_cc2_singles(doubles,singles);
	timer_iter_singles.info();
	if(singles_converged == true) break;
      }
      // assign errors to singles (the difference between the old_singles and the new_ones and not between the last iteration and the new_ones (iter_cc2_singles makes more than one iteration)
      if(world.rank() == 0) std::cout << "Change in Singles functions after all the CC2-Single-Microiterations" << std::endl;
      for(auto& tmp : singles.functions){
	tmp.second.current_error=(tmp.second.function - old_singles(tmp.first).function).norm2();
	if(world.rank() == 0) std::cout << "Change of " << tmp.second.name() << "=" << tmp.second.current_error << std::endl;
      }

      CCOPS.update_intermediates(singles);
      CC_Timer timer_iter_doubles(world,"Iteration " + stringify(iter) + " Doubles");
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
      for(auto& tmp_pair : doubles.allpairs){
	bool pair_converged=iterate_pair(tmp_pair.second,singles);
	if(not pair_converged) doubles_converged=false;
      }
      //doubles_converged = iterate_cc2_doubles(doubles,singles);
      std::vector<double> updated_energies=update_cc2_pair_energies(doubles,singles);
      bool energy_converged=check_energy_convergence(current_energies,updated_energies);
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

      if(world.rank() == 0){
	output("End of Macroiteration " + stringify(iter));
	std::cout << "singles converged: " << singles_converged << std::endl;
	std::cout << "doubles converged: " << doubles_converged << std::endl;
	std::cout << "energy  converged: " << energy_converged << std::endl;
      }

      if(singles_converged and doubles_converged and energy_converged){
	output("Singles and Doubles Converged (!)");
	output("Testing if singles do not change anymore");
	vecfuncT old_singles=singles.get_vecfunction();
	iterate_cc2_singles(doubles,singles);
	vecfuncT new_singles=singles.get_vecfunction();
	vecfuncT difference=sub(world,old_singles,new_singles);
	bool full_convergence=true;

	for(auto x : difference){
	  if(x.norm2() > parameters.dconv_6D) full_convergence=false;
	}
	if(full_convergence){
	  output_section("CC2 CONVERGED!!!");
	  timer_iter_all.info();
	  print_results(doubles,singles);
	  break;
	}else output("Overall convergence not yet reached ... starting cycle again");
      }
      timer_iter_all.info();

    }

    return get_correlation_energy(doubles);
  }

  std::vector<double>
  CC2::update_cc2_pair_energies(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    std::vector<double> omegas;
    for(size_t i=parameters.freeze; i < mo.size(); i++){
      for(size_t j=i; j < mo.size(); j++){
	double tmp=CCOPS.compute_cc2_pair_energy(doubles(i,j),singles(i),singles(j));
	omegas.push_back(tmp);
      }
    }
    if(world.rank() == 0){
      std::cout << "Updated CC2 pair energies:\n";
      for(auto x : omegas)
	std::cout << std::scientific << std::setprecision(parameters.output_prec) << x << std::endl;
    }
    return omegas;
  }

  bool
  CC2::iterate_cc2_singles(const Pairs<CC_Pair> &doubles,CC_vecfunction &singles) {
    if(singles.functions.size() != active_mo.size())
    MADNESS_EXCEPTION(("Wrong size of singles at beginning of iterations " + stringify(singles.functions.size())).c_str(),1);
    output_subsection("Iterate CC2 Singles");
    CC_Timer timer_potential(world,"CC2 Singles Potential");
    vecfuncT potential=CCOPS.get_CC2_singles_potential(singles,doubles);
    timer_potential.info();

    output_subsection("Apply the Green's Operator");
    CC_Timer timer_G(world,"Apply the Green's Operator");
    vecfuncT G_potential=zero_functions<double, 3>(world,potential.size());
    scale(world,potential,-2.0);
    for(size_t i=0; i < potential.size(); i++){
      double epsi=CCOPS.get_orbital_energies()[i + parameters.freeze];
      output("Make Greens Operator for single " + stringify(i + parameters.freeze));
      real_convolution_3d G=BSHOperator<3>(world,sqrt(-2.0 * epsi),parameters.lo,parameters.thresh_bsh_3D);
      real_function_3d tmp=(G(potential[i])).truncate();
      CCOPS.Q(tmp);
      G_potential[i]=tmp;
    }
    timer_G.info();

    std::vector<double> errors;
    bool converged=true;
    for(size_t i=0; i < potential.size(); i++){
      MADNESS_ASSERT(singles(i + parameters.freeze).i == i + parameters.freeze);
      if(world.rank() == 0) std::cout << "|| |tau" + stringify(i + parameters.freeze) + ">|| =" << G_potential[i].norm2() << std::endl;
      CCOPS.Q(G_potential[i]);
      real_function_3d residue=singles(i + parameters.freeze).function - G_potential[i];
      double error=residue.norm2();
      errors.push_back(error);
      if(world.rank() == 0) std::cout << "|| residue" + stringify(i + parameters.freeze) + ">|| =" << error << std::endl;
      CC_function new_single(G_potential[i],singles(i + parameters.freeze).i,PARTICLE);
      new_single.current_error=error;
      singles(i + parameters.freeze)=new_single;
      if(fabs(error) > parameters.dconv_3D) converged=false;
    }
    if(singles.functions.size() != active_mo.size())
    MADNESS_EXCEPTION(("Wrong size of singles at the end of the iteration " + stringify(singles.functions.size())).c_str(),1);
    if(converged) output("singles converged");
    else output("No convergence in singles");
    return converged;
  }

// the constant part of CC2 will be calculated and stored in the pair, after that the pair is iterated with the MP2 alg. till it converges
// Note: if the singles are initialized to 0 then this is just MP2
  bool
  CC2::iterate_pair(CC_Pair &pair,const CC_vecfunction &singles) const {
    bool converged=false;
    calctype type=CC2_;
    if(singles.size() == 0){
      output("No Singles: Current Iterations are MP2");
      type=MP2_;
    }

// check if the constant part has to be recalculated
    bool recalc_const=true;
    if(type == MP2_ and parameters.restart) recalc_const=false;
    if(type == CC2_){
      if(singles(pair.i).current_error < parameters.dconv_3D and singles(pair.j).current_error < parameters.dconv_3D) recalc_const=false;
    }
    if(not pair.constant_term.impl_initialized()){
      output("semi-constant-term was not stored -> recalculation");
      recalc_const=true;
    }

    output_section("Iterate " + pair.name());

    real_function_6d regular_part;
    if(recalc_const){
      output_subsection("Get Regularized-Potential of Pair " + pair.name());
      CC_Timer timer_cc2_regular(world,"Get Semi-Constant CC2 Part of Pair " + pair.name());
      if(type == CC2_) regular_part=CCOPS.make_cc2_residue_sepparated(singles(pair.i),singles(pair.j));
      else if(type == MP2_) regular_part=CCOPS.make_cc2_residue_sepparated(CC_function(real_factory_3d(world),pair.i,HOLE),CC_function(real_factory_3d(world),pair.j,HOLE));
      regular_part.print_size("Regularization part of pair " + pair.name());
      timer_cc2_regular.info();
      pair.constant_term=copy(regular_part);
    }else{
      output("Recalculation of Regularized-Potential of Pair " + pair.name() + " is not neccesary since the involved singles have not changed");
      output("Last change in " + singles(pair.i).name() + "=" + stringify(singles(pair.i).current_error));
      output("Last change in " + singles(pair.j).name() + "=" + stringify(singles(pair.j).current_error));
      regular_part=pair.constant_term;
    }

    real_function_6d coulomb_part=real_factory_6d(world);
    if(type == CC2_){
      CC_Timer timer_cc2_coulomb(world,"Get Screened Coulomb Potentials of CC2 singles");
      const double thresh=CCOPS.guess_thresh(singles(pair.i),singles(pair.j));
      if(thresh < FunctionDefaults<6>::get_thresh()){
	output("Increasing 6D thresh for screened Coulomb parts of singles");
	FunctionDefaults<6>::set_thresh(parameters.tight_thresh_6D);
      }
      coulomb_part=CCOPS.make_cc2_coulomb_parts(singles(pair.i),singles(pair.j),singles);
      coulomb_part.print_size("Coulomb part of Singles for pair " + pair.name());
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      timer_cc2_coulomb.info();
    }

    const real_function_6d constant_part=(regular_part + coulomb_part);
    constant_part.print_size("semi-constant-cc2-part of pair " + pair.name());

    output_subsection("Converge pair " + pair.name() + " on constant singles potential");

    CC_Timer make_BSH_time(world,"Make Detructive Greens Operator");
    real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * CCOPS.get_epsilon(pair.i,pair.j)),parameters.lo,parameters.thresh_bsh_6D);
    G.destructive()=true;
    make_BSH_time.info();

    NonlinearSolverND<6> solver(parameters.kain_subspace);
    solver.do_print=(world.rank() == 0);

    for(size_t iter=0; iter < parameters.iter_max_6D; iter++){
      output_subsection("MP2-Microiteration with Frozen CC2-Singles");
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
      if(type == CC2_) pair.current_energy=CCOPS.compute_cc2_pair_energy(pair,singles(pair.i),singles(pair.j));
      else if(type == MP2_) pair.current_energy=CCOPS.compute_mp2_pair_energy(pair);
      const double delta=pair.current_energy - old_energy;
      pair.current_error=error;
      pair.current_energy_difference=delta;
      pair.iteration++;
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
    if(type == MP2_){
      std::string msg="mp2_";
      if(converged) msg+="converged_pair_";
      else msg+="not_converged_pair_";
      pair.store_pair(world,msg);
    }else{
      pair.store_pair(world,"current_pair_");
    }

    return converged;
  }

  bool
  CC2::iterate_cc2_doubles(Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
    output_subsection("Iterate CC2 Doubles");

    bool converged=true;
    for(size_t i=parameters.freeze; i < mo.size(); i++){
      for(size_t j=i; j < mo.size(); j++){
	CC_Timer whole_potential(world,"whole doubles potential");

	CC_Timer make_BSH_time(world,"Make BSH Operator");
	real_convolution_6d G=BSHOperator<6>(world,sqrt(-2.0 * CCOPS.get_epsilon(i,j)),parameters.lo,parameters.thresh_bsh_6D);
	make_BSH_time.info();

	real_function_6d doubles_potential=CCOPS.get_CC2_doubles_potential(doubles(i,j),singles);

	CC_Timer G_time(world,"Apply Greens Operator to doubles potential");
	doubles_potential.scale(-2.0);
	const real_function_6d G_doubles_potential=G(doubles_potential);
	G_doubles_potential.print_size("-2.0*G(doubles_potential)");
	G_time.info();

	// add the constant term
	output("Do NOT Add constant MP2 Term, since it has been recalculated");
	//real_function_6d u_new = doubles(i,j).constant_term + G_doubles_potential;
	real_function_6d u_new=G_doubles_potential;

	u_new.print_size("u_new");
	CCOPS.apply_Q12(u_new,"u_new");
	u_new.print_size("Q12(u_New)");
	real_function_6d BSH_residue=doubles(i,j).function - u_new;
	const double error=BSH_residue.norm2();
	const double omega=CCOPS.compute_cc2_pair_energy(doubles(i,j),singles(i),singles(j));
	CCOPS.performance_D.current_iteration++;
	std::cout << "Iteration of pair |u" << i << j << "> completed, current error is:" << error << ", current correlation energy is " << std::setprecision(parameters.output_prec) << omega
	    << std::endl;
	if(error > parameters.dconv_6D) converged=false;
	if(world.rank() == 0){
	  if(converged) std::cout << "\n\n\n\t\t\tPair" << i << j << " converged!\n\n\n" << std::endl;
	}
	real_function_6d updated_function=u_new;
	doubles(i,j).function=updated_function;
	doubles(i,j).current_error=error;
	doubles(i,j).current_energy=omega;
	BSH_residue.print_size("Residue");
	doubles(i,j).function.truncate();
	doubles(i,j).store_pair(world);
	whole_potential.info();
      }
    }

    return converged;
  }

  CC_vecfunction
  CC2::initialize_cc2_singles() const {
    output("Initialize CC2-Singles as Zero-functions");
    vecfuncT G_guess_potential=zero_functions<double, 3>(world,mo.size() - parameters.freeze);
    output_section("Initialized CC2 Singles");
    std::vector<CC_function> tmp;
    for(size_t i=0; i < G_guess_potential.size(); i++){
      if(world.rank() == 0) std::cout << "|| |tau" + stringify(i) + ">|| =" << G_guess_potential[i].norm2() << std::endl;
      CC_function taui(G_guess_potential[i],i + parameters.freeze,PARTICLE);
      tmp.push_back(taui);
    }
    return CC_vecfunction(tmp);
  }
// Unnecessary function
  double
  CC2::compute_mp2_pair_energy(CC_Pair &u) const {
    return CCOPS.compute_mp2_pair_energy(u);
  }

/// Initialize an electron pair
  void
  CC2::initialize_electron_pair(CC_Pair &u) const {
    output_subsection("Initialize Electron Pair |u" + stringify(u.i) + stringify(u.j) + ">");
    output("\n\nCheck for saved pairs...");

    output("...No saved pair found... recalculate\n\n");
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
    u.iteration=0;

//	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * epsij), parameters.lo, parameters.thresh_bsh_6D);
//
//	output_subsection("Calculation of constant MP2 potential");
//	CC_Timer timer_const(world,"Calculation of constant MP2 part");
//	real_function_6d mp2_constant_part = CCOPS.get_MP2_potential_constant_part(u).truncate();
//	mp2_constant_part.print_size("mp2_constant_part");
//	timer_const.info();
//
//	output_subsection("Apply the Green's Operator");
//	CC_Timer timer_Gconst(world,"Apply BSH to constant MP2 part");
//	real_function_6d GVPhi = G(-2.0*mp2_constant_part).truncate();
//	GVPhi.print_size("G(-2.0*Q(ConstantTerm))");
//	real_function_6d unprojected_GVPhi = copy(GVPhi);
//	CCOPS.apply_Q12(GVPhi);
//	u.constant_term = copy(GVPhi);
//	u.constant_term.print_size("Q(G(-2.0*Q(ConstantTerm)))");
//	timer_Gconst.info();
//
//	// Make the first guess for the mp2 pair function
//	u.function = copy(GVPhi);
//
//	// Calculate the pair energy
//	double test_energy = compute_mp2_pair_energy(u);
//	output("Initialized Electron Pair: |u" + stringify(u.i) + stringify(u.j) + "> with pair energy: " + stringify(test_energy) + "\n");
//	u.info();
//	u.store_pair(world,"const_");

  }

  static double
  gauss_2s(const coord_3d &r) {

  }
  void
  CC2::decompose_constant_part() {
// works only for one orbital molecules
    output_section("Try to Decompose Constant Part");
    if(mo.size() != 1) output("Will not work, too many orbitals");
    CC_Pair u(0,0);
    const real_function_6d const_part=CCOPS.get_MP2_potential_constant_part(u);
    const double const_norm=const_part.norm2();

// make greens operator
    const double epsi=CCOPS.get_orbital_energies()[0];
    real_convolution_3d Gi=BSHOperator<3>(world,sqrt(-2.0 * epsi),parameters.lo,parameters.thresh_bsh_3D);
    const double epsj=CCOPS.get_orbital_energies()[0];
    real_convolution_3d Gj=BSHOperator<3>(world,sqrt(-2.0 * epsj),parameters.lo,parameters.thresh_bsh_3D);

// make guess basis
    vecfuncT basis1;
    vecfuncT basis2;
    real_function_3d guess1=real_factory_3d(world);
    vecfuncT aos=nemo.get_calc()->ao;
    std::cout << "AO_Basis size is " << aos.size() << std::endl;
    for(const auto ao : aos){
      guess1+=ao;
    }
    CCOPS.Q(guess1);
    guess1.print_size("Guess1");
    plot_plane(world,guess1,"Guess1");
    basis1.push_back(guess1);
    basis2.push_back(guess1);

// test reduced fock matrix
    CC_vecfunction moket(mo,HOLE);
    CCOPS.update_intermediates(moket);
    Tensor<double> test_redF=CCOPS.make_reduced_fock_matrix(moket,epsi);
    std::cout << " Reduced Fock Matrix is: " << test_redF << std::endl;

// for more than one guess the orthonormalization is missing
    for(size_t iter=0; iter < 100; iter++){
      std::cout << "\n\n------\nIter" << iter << std::endl;
      const CC_vecfunction ccbas1(basis1,UNDEFINED);
      const CC_vecfunction ccbas2(basis2,UNDEFINED);
      CCOPS.update_intermediates(ccbas1);
      Tensor<double> fock_matrix_basis1=CCOPS.make_reduced_fock_matrix(ccbas1,epsi);
      const vecfuncT fock_residue_basis1=CCOPS.fock_residue_closed_shell(ccbas1);
      CCOPS.update_intermediates(ccbas2);
      Tensor<double> fock_matrix_basis2=CCOPS.make_reduced_fock_matrix(ccbas2,epsj);
      const vecfuncT fock_residue_basis2=CCOPS.fock_residue_closed_shell(ccbas2);
      vecfuncT basis1_const_part;
      for(const auto &b : basis1)
	basis1_const_part.push_back(const_part.project_out((b * nemo.nuclear_correlation->square()),0));
      vecfuncT basis2_const_part;
      for(const auto &b : basis2)
	basis2_const_part.push_back(const_part.project_out((b * nemo.nuclear_correlation->square()),0));

      for(size_t i=0; i < basis1.size(); i++){
	real_function_3d pot1i=fock_residue_basis1[i] + basis2_const_part[i];
	for(size_t j=0; j < basis2.size(); j++){
	  pot1i+=fock_matrix_basis2(i,j) * basis1[j];
	}
	pot1i.scale(-2.0);
	real_function_3d new_basis1_i=Gi(pot1i);
	CCOPS.Q(new_basis1_i);
	const double diff=(basis1[i] - new_basis1_i).norm2();
	basis1[i]=new_basis1_i;
	std::cout << "error basis1_" << i << "=" << diff << std::endl;
      }

      for(size_t i=0; i < basis2.size(); i++){
	real_function_3d pot2i=fock_residue_basis2[i] + basis1_const_part[i];
	for(size_t j=0; j < basis1.size(); j++){
	  pot2i+=fock_matrix_basis1(i,j) * basis2[j];
	}
	pot2i.scale(2.0);
	real_function_3d new_basis2_i=Gi(pot2i);
	CCOPS.Q(new_basis2_i);
	const double diff=(basis2[i] - new_basis2_i).norm2();
	basis1[i]=new_basis2_i;
	std::cout << "error basis2_" << i << "=" << diff << std::endl;
      }

      double omega_u=0.0;
      for(size_t i=0; i < basis1.size(); i++){
	omega_u+=CCOPS.make_ijgxy(0,0,basis1[i],basis2[i]);
      }
      std::cout << "omega=" << omega_u << std::endl;

    }

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
