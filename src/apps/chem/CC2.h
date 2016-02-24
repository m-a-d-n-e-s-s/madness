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




  class CC2 {
  public:



    CC2(World &world_,const std::string &inputFileName, const Nemo &nemo_):
      world(world_),
      //correlationfactor(world,1.0,1.e-7,nemo_.get_calc()->molecule),
      parameters(inputFileName, nemo_.get_calc() -> param.lo),
      nemo(nemo_),
      mo(nemo_.get_calc()->amo),
      active_mo(make_active_mo()),
      CCOPS(world,nemo,parameters)
  {
      output_section("CC2 Class has been initialized with the following parameters");
      // set the threshholds
      // Set Protocoll
      output("Set Protocol 3D");
      nemo_.get_calc() -> set_protocol<3>(world,parameters.thresh_3D);
      output("Set Protocol 6D");
      nemo_.get_calc() -> set_protocol<6>(world,parameters.thresh_6D);

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
      if (active_mo.size()+parameters.freeze != CCOPS.mo_ket().size()) CCOPS.error("active_mo + freeze of CC2 class and mo_ket_ of CC_Operators have not the same size");
      if (active_mo.size()+parameters.freeze != CCOPS.mo_bra().size()) CCOPS.error("active_mo + freeze of CC2 class and mo_bra_ of CC_Operators have not the same size");
      output("Active molecular orbitals have been created...");
      if(world.rank()==0) std::cout << mo.size() << " MOs\n " << active_mo.size() << " Active MOs\n" << parameters.freeze << "frozen MOs\n";


      std::string nuc = "???";
      if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::None) nuc="None";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::GaussSlater) nuc="GaussSlater";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::GradientalGaussSlater) nuc="GradientalGaussSlater";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::LinearSlater) nuc="LinearSlater";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::Polynomial) nuc="Polynomial";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::Slater) nuc="Slater";
      else if(nemo.nuclear_correlation -> type()==NuclearCorrelationFactor::Two) nuc="Two";
      if(world.rank()==0) std::cout << "Nuclear Correlation Factor is " << nuc << std::endl;

      //output_section("Testing Section in Constructor");
      //CCOPS.test_fill_tree();
  }
    vecfuncT make_active_mo(){
      if(mo.empty()) MADNESS_EXCEPTION("Tried to init. active MOs, but MO vector is empty",1);
      if(parameters.freeze != 0){
	output("Make Active MOs from " + stringify(parameters.freeze) + " to " + stringify(mo.size()));
	vecfuncT tmp;
	for(size_t i=parameters.freeze; i<mo.size();i++){
	  tmp.push_back(mo[i]);
	}
	return tmp;
      }else{
	output("No freezing demanded, active_mo = mo");
	return mo;
      }
    }
    void plot(const real_function_3d &f, const std::string &msg = "unspecified function")const{
      plot_plane(world,f,msg);
      output("Plotted " + msg);
    }
    /// Check energy convergence: Creates the difference between two vectors and compares against given thresh in parameters
    bool check_energy_convergence(const std::vector<double> &current, const std::vector<double> &updated)const{
      if(current.size()!=updated.size())MADNESS_EXCEPTION("error in energy convergence check: different sizes in vectors",1);
      bool conv = true;
      std::vector<double> diff(current.size(),0.0);
      for(size_t i=0;i<current.size();i++){
	double diffi = updated[i] - current[i];
	diff[i] = diffi;
	if(diffi > parameters.econv) conv=false;
      }
      if(world.rank()==0){
	std::cout << "\n\n";
	std::cout << "Pair Correlation Energies: New, Old, Diff\n";
	for(size_t i=0;i<current.size();i++) std::cout << updated[i] << ", " << current[i] << ", " << diff[i] << std::endl;
	std::cout << "\n\n";
      }
      return conv;
    }
    /// make consistency tests
    bool test()const;
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
    void solve();
    std::vector<std::pair<CC_vecfunction,double> > solve_ccs();
    /// solve the MP2 equations (uncoupled -> Canonical Orbitals)
    double solve_mp2(Pairs<CC_Pair> &doubles);
    double solve_mp2_nonorthogonal(Pairs<CC_Pair> &doubles);
    double solve_cc2(Pairs<CC_Pair> &u, CC_vecfunction &tau);
    double solve_cc2_response(const CC_vecfunction &tau, const Pairs<CC_Pair> &u,CC_vecfunction x, Pairs<CC_Pair> &chi);
    double solve_cispd();
    double solve_cispd(Pairs<CC_Pair> &doubles,const Pairs<CC_Pair> &mp2_pairs, const CC_vecfunction & cis_singles, const double cis_omega);
    bool iterate_cc2_singles(const Pairs<CC_Pair> &doubles, CC_vecfunction &singles);
    bool iterate_cc2_singles_response(const CC_vecfunction &tau, const Pairs<CC_Pair> &u,CC_vecfunction x, const Pairs<CC_Pair> &chi) {
      output_subsection("Iterate Response of CC2 Singles");
      const double current_omega = x.omega;
      output("Current Omega is " + std::to_string(current_omega));
      CC_Timer timer_potential(world,"Response of CC2 Singles Potential");
      vecfuncT potential=CCOPS.get_CC2_singles_response_potential(tau,u,x,chi);
      timer_potential.info();

      output_subsection("Apply the Green's Operator");
      CC_Timer timer_G(world,"Apply the Green's Operator");
      vecfuncT G_potential=zero_functions<double, 3>(world,potential.size());
      scale(world,potential,-2.0);
      for(size_t i=0; i < potential.size(); i++){
        double epsi=CCOPS.get_orbital_energies()[i + parameters.freeze]+current_omega;
        output("Make Greens Operator for single " + stringify(i + parameters.freeze));
        real_convolution_3d G=BSHOperator<3>(world,sqrt(-2.0 * epsi),parameters.lo,parameters.thresh_bsh_3D);
        real_function_3d tmp=(G(potential[i])).truncate();
        G_potential[i]=tmp;
      }
      G_potential=CCOPS.apply_Q(G_potential,"G_potential");
      timer_G.info();

      const vecfuncT residue = sub(world,x.get_vecfunction(),G_potential);
      const vecfuncT res_bra = mul(world,nemo.nuclear_correlation->square(),residue);
      const vecfuncT GV_bra = mul(world,nemo.nuclear_correlation->square(),G_potential);
      const double tmp1 = inner(world,res_bra,G_potential).sum();
      const double tmp2 = inner(world,GV_bra,G_potential).sum();

      const double delta = 0.5*tmp1/tmp2;
      if(world.rank()==0) std::cout << " Delta is " << delta << "\n";
      const double new_omega = current_omega + delta;
      x.omega = new_omega;

      std::vector<double> errors;
      bool converged=true;
      for(size_t i=0; i < potential.size(); i++){
        MADNESS_ASSERT(x(i + parameters.freeze).i == i + parameters.freeze);
        if(world.rank() == 0) std::cout << "|| |tau" + stringify(i + parameters.freeze) + ">|| =" << G_potential[i].norm2() << std::endl;
        real_function_3d residue=x(i + parameters.freeze).function - G_potential[i];
        double error=residue.norm2();
        errors.push_back(error);
        if(world.rank() == 0) std::cout << "|| residue" + stringify(i + parameters.freeze) + ">|| =" << error << std::endl;
        CC_function new_x(G_potential[i],x(i + parameters.freeze).i,RESPONSE);
        new_x.current_error=error;
        x(i + parameters.freeze)=new_x;
        if(fabs(error) > parameters.dconv_3D) converged=false;
      }
      if(converged) output("response singles converged");
      else output("No convergence in response singles");
      return converged;
    }
    bool iterate_cc2_doubles( Pairs<CC_Pair> &doubles, const CC_vecfunction &singles)const;
    /// Compute the pair correlation energy of an electron pair function at mp2/CCD level (no singles contributions)
    double compute_mp2_pair_energy(CC_Pair &u)const;
    CC_vecfunction initialize_cc2_singles()const;
    Pairs<CC_Pair> initialize_pairs(const pairtype type, const double omega=0.0)const;
    /// Initialize an electron pair
    void initialize_electron_pair(CC_Pair &u)const;
    /// Calculate the current CC2 correlation energy
    double get_correlation_energy(const Pairs<CC_Pair> &doubles)const;
    /// update the pair energies of cc2
    std::vector<double> update_cc2_pair_energies(const Pairs<CC_Pair> &doubles, const CC_vecfunction &singles)const;
    /// Iterates the CC2 singles equations
    void iterate_singles(vecfuncT &singles, const Pairs<real_function_6d> &doubles)const;
    /// Iterates the CC2 doubles equations
    void iterate_doubles(const vecfuncT &singles, Pairs<real_function_6d> &doubles)const;
    /// Iterates a pair of the CC2 doubles equations
    bool iterate_pair(CC_Pair & pair, const CC_vecfunction &singles)const;
    bool iterate_pair(CC_Pair &pair,const CC_vecfunction &singles, const CC_vecfunction &response_singles,const calctype ctype) const;
    bool iterate_nonorthogonal_pair(CC_Pair &pair);
    /// Create formated output, std output with world rank 0
    void output(const std::string &msg)const{
      if(world.rank()==0) std::cout << msg << "\n";
    }
    /// Create formated output, New programm section
    void output_section(const std::string&msg)const{
      if(world.rank()==0){
	std::cout << std::setw(100) << std::setfill('#') << std::endl;
	std::cout << "\n" << msg << "\n";
	std::cout << std::setw(100) << std::setfill('#') << "\n" << std::endl;
      }
    }
    /// Create formated output, New programm subsection
    void output_subsection(const std::string&msg)const{
      if(world.rank()==0){
	std::cout << std::setw(50) << std::setfill('*') << std::endl;
	std::cout << "\n" << msg << "\n";
	std::cout << std::setw(50) << std::setfill('*') << "\n" << std::endl;
      }
    }
    void decompose_constant_part();

    void print_results(const Pairs<CC_Pair> &doubles, const CC_vecfunction &singles)const{
      const double Ecorr = get_correlation_energy(doubles);
      output("\n Results:\n");
      const size_t prec = std::max<double>(6,parameters.output_prec);
      std::cout <<std::setw(5)<<std::setfill(' ')<< "Pair" << "|"
	  <<std::setw(prec+1)<<std::setfill(' ')<< "omega" <<"|"
	  <<std::setw(prec+1)<<std::setfill(' ')<< "omega*2" << "|"
	  <<std::setw(7)<<std::setfill(' ')<< "error" << "|"
	  <<std::setw(7)<<std::setfill(' ')<< "deltaE" << "|"
	  <<std::setw(7)<<std::setfill(' ')<< "||uij||" << "|"
	  <<std::setw(7)<<std::setfill(' ')<< "||ti||" << "|"
	  <<std::setw(7)<<std::setfill(' ')<< "||tj||" <<"\n" ;
      for(const auto utmp:doubles.allpairs){
	const CC_Pair & u=utmp.second;
	double omega = u.current_energy;
	if(u.i!=u.j) omega = 2.0*u.current_energy;
	if(world.rank()==0){
	  std::cout << std::fixed << std::setprecision(prec)
	  <<std::setw(5)<<std::setfill(' ')<< u.name() << "|"
	  <<std::setw(prec+1)<<std::setfill(' ')<< u.current_energy <<"|"
	  <<std::setw(prec+1)<<std::setfill(' ')<< omega << "|"
	  << std::scientific << std::setprecision(3)
	  <<std::setw(7)<<std::setfill(' ')<< u.current_error << "|"
	  <<std::setw(7)<<std::setfill(' ')<< u.current_energy_difference << "|"
	  <<std::setw(7)<<std::setfill(' ')<< u.function.norm2() <<"|"
	  <<std::setw(7)<<std::setfill(' ')<< singles(u.i).function.norm2() <<"|"
	  <<std::setw(7)<<std::setfill(' ')<< singles(u.j).function.norm2() <<"\n";
	}
      }
      if(world.rank()==0) std::cout << "\n ---> overall correlation energy: " << std::fixed << std::setprecision(parameters.output_prec) << Ecorr << std::endl;
    }
  };

} /* namespace madness */

#endif /* CC2_H_ */
