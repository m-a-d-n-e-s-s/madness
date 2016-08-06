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
#include <chem/CCOperators.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include <chem/TDA.h>
#include <examples/nonlinsol.h>

namespace madness {

  class CC2{
  public:

    CC2(World &world_,const std::string &inputFileName,const Nemo &nemo_)
	: world(world_),
	//correlationfactor(world,1.0,1.e-7,nemo_.get_calc()->molecule),
	parameters(inputFileName,nemo_.get_calc()->param.lo), nemo(nemo_), mo(nemo_.get_calc()->amo), active_mo(make_active_mo()), CCOPS(world,nemo,parameters) {
      output_section("CC2 Class has been initialized with the following parameters");
      // set the threshholds
      // Set Protocoll
      output("Set Protocol 3D");
      nemo_.get_calc()->set_protocol < 3 > (world, parameters.thresh_3D);
      output("Set Protocol 6D");
      nemo_.get_calc()->set_protocol < 6 > (world, parameters.thresh_6D);

      FunctionDefaults<3>::set_thresh(parameters.thresh_3D);
      FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
      // Make shure that k is the same in 3d and 6d functions
      FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
      // by default SCF sets the truncate_mode to 1
      FunctionDefaults<3>::set_truncate_mode(3);
      FunctionDefaults<6>::set_truncate_mode(3);
      parameters.information(world);
      parameters.sanity_check(world);
      // sanity checks
      if(active_mo.size() + parameters.freeze != CCOPS.mo_ket().size()) CCOPS.error("active_mo + freeze of CC2 class and mo_ket_ of CC_Operators have not the same size");
      if(active_mo.size() + parameters.freeze != CCOPS.mo_bra().size()) CCOPS.error("active_mo + freeze of CC2 class and mo_bra_ of CC_Operators have not the same size");
      output("Active molecular orbitals have been created...");
      if(world.rank() == 0) std::cout << mo.size() << " MOs\n " << active_mo.size() << " Active MOs\n" << parameters.freeze << "frozen MOs\n";

      std::string nuc="???";
      if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::None) nuc="None";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::GaussSlater) nuc="GaussSlater";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::GradientalGaussSlater) nuc="GradientalGaussSlater";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::LinearSlater) nuc="LinearSlater";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::Polynomial) nuc="Polynomial";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::Slater) nuc="Slater";
      else if(nemo.nuclear_correlation->type() == NuclearCorrelationFactor::Two) nuc="Two";
      if(world.rank() == 0) std::cout << "Nuclear Correlation Factor is " << nuc << std::endl;

      //output_section("Testing Section in Constructor");
      //CCOPS.test_fill_tree();
    }
    vecfuncT
    make_active_mo() {
      if(mo.empty()) MADNESS_EXCEPTION("Tried to init. active MOs, but MO vector is empty",1);
      if(parameters.freeze != 0){
	output("Make Active MOs from " + stringify(parameters.freeze) + " to " + stringify(mo.size()));
	vecfuncT tmp;
	for(size_t i=parameters.freeze; i < mo.size(); i++){
	  tmp.push_back(mo[i]);
	}
	return tmp;
      }else{
	output("No freezing demanded, active_mo = mo");
	return mo;
      }
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
    const CC_Parameters parameters;
    /// The SCF Calculation
    const Nemo &nemo;
    /// Molecular orbitals (all of them, NEMOS!!! )
    const vecfuncT mo;
    /// Active MO
    const vecfuncT active_mo;
    /// The CC Operator Class
    CC_Operators CCOPS;

    /// solve the CC2 ground state equations, returns the correlation energy
    void
    solve();
    std::vector<std::pair<CC_vecfunction, double> >
    solve_ccs();
    /// solve the MP2 equations (uncoupled -> Canonical Orbitals)
    double
    solve_mp2(Pairs<CC_Pair> &doubles);
    double
    solve_mp2_nonorthogonal(Pairs<CC_Pair> &doubles);
    double
    solve_cc2(Pairs<CC_Pair> &u,CC_vecfunction &tau);
    double
    solve_cc2_response(const CC_vecfunction &tau,const Pairs<CC_Pair> &u,CC_vecfunction x,Pairs<CC_Pair> &chi);
    double
    solve_cispd();
    double
    solve_cispd(Pairs<CC_Pair> &doubles,const Pairs<CC_Pair> &mp2_pairs,const CC_vecfunction & cis_singles,const double cis_omega);

    // doubles[0] = gs_doubles, doubles[1] = response_doubles
    bool
    iterate_singles(CC_vecfunction &singles,const CC_vecfunction singles2,const std::vector<Pairs<CC_Pair>> &doubles,const calctype ctype) {
      output_subsection("Iterate " + assign_name(ctype) + "-Singles");
      CC_Timer time_all(world,"Overall Iteration of " + assign_name(ctype) + "-Singles");
      bool converged=true;

      CC_vecfunction old_singles(singles);
      for(auto& tmp : singles.functions)
	old_singles(tmp.first).function=copy(tmp.second.function);

      // KAIN solver
      typedef allocator<double, 3> allocT;
      typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
      allocT alloc(world,singles.size());
      solverT solver(allocT(world,singles.size()));
      solver.do_print=(world.rank() == 0);

      for(size_t iter=0; iter < parameters.iter_max_3D; iter++){
	output_subsection("Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
	CC_Timer time(world,"Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
	double omega=0.0;
	if(ctype == LRCC2_) omega=singles.omega;
	else if(ctype == LRCCS_) omega=singles.omega;

	// consistency check
	switch(ctype){
	  case CC2_:
	    if(singles.type != PARTICLE) CCOPS.warning("iterate_singles: CC2 demanded but singles are not of type PARTICLE");
	    break;
	  case MP2_:
	    CCOPS.error("Demanded Singles Calculation for MP2 ????");
	    break;
	  case LRCC2_:
	    if(singles.type != RESPONSE or singles2.type != PARTICLE) CCOPS.warning("iterate_singles: LRCC2_ singles have wrong types");
	    break;
	  case LRCCS_:
	    if(singles.type != RESPONSE) CCOPS.warning("iterate_singles: LRCCS_ singles have wrong types");
	    break;
	  case CISpD_:
	    CCOPS.error("Demanded Singles Calculation for CIS(D)");
	    break;
	  case experimental_:
	    CCOPS.error("Iterate Singles not implemented for Experimental calculation");
	    break;
	  default:
	    CCOPS.error("Unknown calculation type in iterate singles: " + assign_name(ctype));
	}

	// get potentials
	CC_Timer time_V(world,assign_name(ctype) + "-Singles Potential");
	vecfuncT V;
	if(ctype == CC2_) V=CCOPS.get_CC2_singles_potential(singles,doubles.front());
	else if(ctype == LRCC2_) V=CCOPS.get_CC2_singles_response_potential(singles2,doubles.front(),singles,doubles.back());
	else if(ctype == LRCCS_) V=CCOPS.get_LRCCS_potential(singles);
	else CCOPS.error("iterate singles: unknown type");
	time_V.info();

	if(ctype == LRCCS_){
	  const double expv=CCOPS.compute_cis_expectation_value(singles,V);
	  if(world.rank() == 0) std::cout << "Current CCS/CIS Expectation Value " << expv << "\n";
	  if(world.rank() == 0) std::cout << "using expectation-value for bsh-operator\n";
	  singles.omega = expv;
	  omega = expv;
	}

	scale(world,V,-2.0);
	truncate(world,V);

	// make bsh operators
	CC_Timer time_makebsh(world,"Make G-Operators");
	std::vector < std::shared_ptr<SeparatedConvolution<double, 3> > > G(singles.size());
	for(size_t i=0; i < G.size(); i++){
	  const double bsh_eps=CCOPS.get_orbital_energies()[i + parameters.freeze] + omega;
	  G[i]=std::shared_ptr < SeparatedConvolution<double, 3> > (BSHOperatorPtr3D(world,sqrt(-2.0 * bsh_eps),parameters.lo,parameters.thresh_bsh_3D));
	}
	world.gop.fence();
	time_makebsh.info();

	// apply bsh operators
	CC_Timer time_applyG(world,"Apply G-Operators");
	vecfuncT GV=apply<SeparatedConvolution<double, 3>, double, 3>(world,G,V);
	world.gop.fence();
	time_applyG.info();

	// apply Q-Projector to result
	GV=CCOPS.apply_Q(GV);

	if(ctype==LRCCS_){
	  output("Normalizing new singles");
	  const vecfuncT x = GV;
	  const vecfuncT xbra = mul(world,nemo.nuclear_correlation->square(),GV);
	  const double norm = sqrt(inner(world,xbra,x).sum());
	  if(world.rank()==0) std::cout << " Norm was " <<std::fixed<< std::setprecision(parameters.output_prec) << norm << "\n";
	  scale(world,GV,1.0/norm);
	}

	// residual
	const vecfuncT residual=sub(world,singles.get_vecfunction(),GV);

	// information with and without nuclear correlation factor
	const Tensor<double> xinnerx=inner(world,singles.get_vecfunction(),singles.get_vecfunction());
	const Tensor<double> R2xinnerx=inner(world,mul(world,nemo.nuclear_correlation->square(),singles.get_vecfunction()),singles.get_vecfunction());
	const Tensor<double> GVinnerGV=inner(world,GV,GV);
	const Tensor<double> R2GVinnerGV=inner(world,mul(world,nemo.nuclear_correlation->square(),GV),GV);
	const Tensor<double> rinnerr=inner(world,residual,residual);
	const Tensor<double> R2rinnerr=inner(world,mul(world,nemo.nuclear_correlation->square(),residual),residual);
	const double R2vector_error=sqrt(R2rinnerr.sum());

	// print information
	if(world.rank() == 0) std::cout << "\n\n-----Results of current interation:-----\nresult with nuclear correlation factor (result without)\n";
	if(world.rank() == 0) std::cout << "\nName: ||" << singles.name() << "||, ||GV" << singles.name() << ", ||residual||" << "\n";
	if(world.rank() == 0)
	  std::cout << singles.name() << ": " << std::scientific << std::setprecision(parameters.output_prec) << sqrt(R2xinnerx.sum()) << " (" << sqrt(xinnerx.sum()) << "), "
	      << sqrt(R2GVinnerGV.sum()) << " (" << sqrt(GVinnerGV.sum()) << "), " << sqrt(R2rinnerr.sum()) << " (" << sqrt(rinnerr.sum()) << "), \n----------------------------------------\n";
	for(size_t i=0; i < GV.size(); i++){
	  if(world.rank() == 0)
	    std::cout << singles(i + parameters.freeze).name() << ": " << std::scientific << std::setprecision(parameters.output_prec) << sqrt(R2xinnerx(i)) << " (" << sqrt(xinnerx(i)) << "), "
		<< sqrt(R2GVinnerGV(i)) << " (" << sqrt(GVinnerGV(i)) << "), " << sqrt(R2rinnerr(i)) << " (" << sqrt(rinnerr(i)) << "), ";
	}
	if(world.rank() == 0) std::cout << "\n----------------------------------------\n\n";

	// make second order update (only for response)
	if(ctype == LRCC2_ or ctype == LRCCS_){
	  output("\nMake 2nd order energy update:");
	  double tmp=inner(world,residual,V).sum();
	  double tmp2=inner(world,GV,GV).sum();
	  const double delta=(0.5 * tmp / tmp2);
	  // include nuclear factors
	  {
	    vecfuncT bra_res=mul(world,nemo.nuclear_correlation->square(),residual);
	    vecfuncT bra_GV=mul(world,nemo.nuclear_correlation->square(),GV);
	    double Rtmp=inner(world,bra_res,V).sum();
	    double Rtmp2=inner(world,bra_GV,GV).sum();
	    const double Rdelta=(0.5 * Rtmp / Rtmp2);
	    double old_omega=omega;
	    if(fabs(delta) < 0.1) omega+=Rdelta;
	    if(world.rank() == 0)
	      std::cout << "omega, old_omega, Rdelta, (delta)" << std::fixed << std::setprecision(parameters.output_prec + 2) << omega << ", " << old_omega << ", " << Rdelta << ", (" << delta
		  << ")\n\n";
	  }

	}

	// update singles
	singles.omega=omega;
	vecfuncT new_singles=GV;
	if(parameters.kain) new_singles=solver.update(singles.get_vecfunction(),residual).x;
	for(size_t i=0; i < GV.size(); i++){
	  singles(i + parameters.freeze).function=copy(new_singles[i]);
	}

	// update intermediates
	if(singles.type == RESPONSE) CCOPS.update_response_intermediates(singles);
	else if(singles.type == PARTICLE) CCOPS.update_intermediates(singles);

	converged=(R2vector_error < parameters.dconv_3D);

	time.info();
	if(converged) break;
      }
      time_all.info();

      // Assign the overall changes
      bool no_change=true;
      if(world.rank() == 0) std::cout << "Change in Singles functions after all the CC2-Single-Microiterations" << std::endl;
      for(auto& tmp : singles.functions){
	const double change = (tmp.second.function - old_singles(tmp.first).function).norm2();
	tmp.second.current_error= change;
	if(change>parameters.dconv_6D*0.1) no_change=false;
	if(world.rank() == 0) std::cout << "Change of " << tmp.second.name() << "=" << tmp.second.current_error << std::endl;
      }

      CCOPS.plot(singles);
      CCOPS.save_functions(singles);
      if(no_change) output("Change of Singles was below (0.1*thresh_6D) = " +std::to_string(parameters.dconv_6D*0.1)+"!");
      return no_change;
    }

    bool
    iterate_cc2_singles(const Pairs<CC_Pair> &doubles,CC_vecfunction &singles);

    bool
    iterate_cc2_doubles(Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const;
    /// Compute the pair correlation energy of an electron pair function at mp2/CCD level (no singles contributions)
    double
    compute_mp2_pair_energy(CC_Pair &u) const;
    CC_vecfunction
    initialize_cc2_singles() const;
    Pairs<CC_Pair>
    initialize_pairs(const pairtype type,const double omega=0.0) const;
    /// Initialize an electron pair
    void
    initialize_electron_pair(CC_Pair &u) const;
    /// Calculate the current CC2 correlation energy
    double
    get_correlation_energy(const Pairs<CC_Pair> &doubles) const;
    /// update the pair energies of cc2
    std::vector<double>
    update_cc2_pair_energies(Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const;
    /// Iterates a pair of the CC2 doubles equations
    bool
    iterate_pair(CC_Pair & pair,const double omega=0.0) const;
    bool
    iterate_nonorthogonal_pair(CC_Pair &pair);

    bool update_constant_part_mp2(CC_Pair &pair){
      if(pair.constant_term.is_initialized()) return false;
      else pair.constant_term = CCOPS.make_regularization_residue(CCOPS.mo_ket(pair.i),CCOPS.mo_ket(pair.j),MP2_,0.0);
      return true;
    }

    bool update_constant_part_cc2_gs(CC_Pair &pair, const CC_vecfunction &tau){
      if(pair.type!=GROUND_STATE)CCOPS.error("asked for constant part of ground state, but given pair is not a ground state pair");
      const CC_function taui = tau(pair.i);
      const CC_function tauj = tau(pair.j);
      bool recalc=(taui.current_error > parameters.dconv_3D or tauj.current_error>parameters.dconv_3D);
      if(not pair.constant_term.is_initialized()) recalc=true;

      if(recalc){
	output_section("(Re)-Calculating the (Semi-) Constant-Terms of the CC2-Ground-State");
	const CC_vecfunction t = CCOPS.make_t_intermediate_full(tau);
	pair.constant_term = CCOPS.make_constant_part_cc2_gs(t(pair.i),t(pair.j),t,0.0);
      }else output("Singles did not change significantly: Constant part is not recalculated");
      return recalc;
    }

    bool update_constant_part_cc2_response(CC_Pair &pair, const CC_vecfunction &tau,const CC_vecfunction &x){
      if(pair.type!=EXCITED_STATE)CCOPS.error("asked for constant part of the response, but given pair is not a response pair");
      if(x.type!=RESPONSE) error("update_constant_part_response: x!=RESPONSE");
      if(tau.type!=PARTICLE) error("update_constant_part_response: x!=PARTICLE");

      bool recalc=(x(pair.i).current_error > parameters.dconv_3D or x(pair.j).current_error>parameters.dconv_3D);
      if(not pair.constant_term.is_initialized()) recalc=true;

      if(recalc){
	output_section("(Re)-Calculating the (Semi-) Constant-Terms of the CC2-Response");
	pair.constant_term = CCOPS.make_constant_part_cc2_response(std::make_pair(pair.i,pair.j),x,CCOPS.mo_ket());
      }else output("Singles did not change significantly: Constant part is not recalculated");
      return recalc;
    }
    /// Create formated output, std output with world rank 0
    void
    output(const std::string &msg) const {
      if(world.rank() == 0) std::cout << msg << "\n";
    }
    /// Create formated output, New programm section
    void
    output_section(const std::string&msg) const {
      if(world.rank() == 0){
	std::cout << std::setw(100) << std::setfill('#') << std::endl;
	std::cout << "\n" << msg << "\n";
	std::cout << std::setw(100) << std::setfill('#') << "\n" << std::endl;
      }
    }
    /// Create formated output, New programm subsection
    void
    output_subsection(const std::string&msg) const {
      if(world.rank() == 0){
	std::cout << std::setw(50) << std::setfill('*') << std::endl;
	std::cout << "\n" << msg << "\n";
	std::cout << std::setw(50) << std::setfill('*') << "\n" << std::endl;
      }
    }
    void
    decompose_constant_part();

    void
    print_results(const Pairs<CC_Pair> &doubles,const CC_vecfunction &singles) const {
      const double Ecorr=get_correlation_energy(doubles);
      output("\n Results:\n");
      const size_t prec=std::max<double>(6,parameters.output_prec);
      std::cout << std::setw(5) << std::setfill(' ') << "Pair" << "|" << std::setw(prec + 1) << std::setfill(' ') << "omega" << "|" << std::setw(prec + 1) << std::setfill(' ') << "omega*2" << "|"
	  << std::setw(7) << std::setfill(' ') << "error" << "|" << std::setw(7) << std::setfill(' ') << "deltaE" << "|" << std::setw(7) << std::setfill(' ') << "||uij||" << "|" << std::setw(7)
	  << std::setfill(' ') << "||ti||" << "|" << std::setw(7) << std::setfill(' ') << "||tj||" << "\n";
      for(const auto utmp : doubles.allpairs){
	const CC_Pair & u=utmp.second;
	double omega=u.current_energy;
	if(u.i != u.j) omega=2.0 * u.current_energy;
	double delta=u.current_energy_difference;
	double delta_sign = delta/fabs(delta);
	std::string sign = "+"; if(delta_sign<0.0) sign="-";
	const std::string sdelta = sign+std::to_string(fabs(delta));
	if(world.rank() == 0){
	  std::cout << std::fixed << std::setprecision(prec) << std::setw(5) << std::setfill(' ') << u.name() << "|" << std::setw(prec + 1) << std::setfill(' ') << u.current_energy << "|"
	      << std::setw(prec + 1) << std::setfill(' ') << omega << "|" << std::scientific << std::setprecision(3) << std::setw(7) << std::setfill(' ') << u.current_error << "|" << std::setw(7)
	      << std::setfill(' ') << sdelta << "|" << std::setw(7) << std::setfill(' ') << u.function.norm2() << "|" << std::setw(7) << std::setfill(' ')
	      << singles(u.i).function.norm2() << "|" << std::setw(7) << std::setfill(' ') << singles(u.j).function.norm2() << "\n";
	}
      }
      if(world.rank() == 0) std::cout << "\n ---> overall correlation energy: " << std::fixed << std::setprecision(parameters.output_prec) << Ecorr << std::endl;
    }
  };

} /* namespace madness */

#endif /* CC2_H_ */
