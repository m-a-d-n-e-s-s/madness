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
	double solve_cispd();
	double solve_cispd(Pairs<CC_Pair> &doubles,const Pairs<CC_Pair> &mp2_pairs, const CC_vecfunction & cis_singles, const double cis_omega);
	bool iterate_cc2_singles(const Pairs<CC_Pair> &doubles, CC_vecfunction &singles);
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
	bool iterate_pair(CC_Pair &pair,const CC_vecfunction &singles, const CC_vecfunction response_singles,const calctype ctype) const;
	bool iterate_nonorthogonal_pair(CC_Pair &pair) const;
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
	  calctype ctype = CC2_;
	  if(CCOPS.make_norm(singles)==0.0) ctype = MP2_;
	  for(auto x : doubles.allpairs){
	    const std::string ij=std::to_string(x.second.i) + std::to_string(x.second.j);
	    const real_function_6d full_pair=CCOPS.make_full_pair_function(x.second,singles(x.second.i),singles(x.second.j));
	    const real_function_6d tmp=multiply(copy(full_pair),nemo.nuclear_correlation->function(),1);
	    const real_function_6d R12_full_pair=multiply(tmp,nemo.nuclear_correlation->function(),2);
	    double single_i=0.0;
	    if(ctype==CC2_) single_i=singles(x.second.i).function.norm2();
	    double single_j=0.0;
	    if(ctype==CC2_) single_j=singles(x.second.j).function.norm2();
	    double R_single_i= 0.0;
	    if(ctype==CC2_) R_single_i=(singles(x.second.i).function * nemo.nuclear_correlation->function()).norm2();
	    double R_single_j=0.0;
	    if(ctype==CC2_) R_single_j=(singles(x.second.j).function * nemo.nuclear_correlation->function()).norm2();
	    const double norm_pair=x.second.function.norm2();
	    const double norm_full_pair=full_pair.norm2();
	    const double norm_R12_full_pair=R12_full_pair.norm2();
	    if(world.rank() == 0){
	      std::cout.precision(parameters.output_prec);
	      std::cout << std::fixed;
	      std::cout << "Pair " << x.second.name() << ": Correlation Energy=" << x.second.current_energy  << std::endl;
	      std::cout << std::setw(15) << "|| u" + ij + "||=" << norm_pair << std::setw(15) << ", ||tau" + ij + "||=" << norm_full_pair << std::setw(15) << ", ||R12tau" + ij + "||="
		  << norm_R12_full_pair << "\n";
	      std::cout << std::setw(15) << "||tau" + std::to_string(x.second.i) + " ||=" << single_i << std::setw(15) << ", ||tau" + std::to_string(x.second.j) + "||=" << single_j << std::setw(15)
		  << ", ||Rtau" + std::to_string(x.second.i) + "||=" << R_single_i << std::setw(15) << ", ||Rtau" + std::to_string(x.second.j) + "||=" << R_single_j << "\n" << std::endl;
	    }
	  }
	}
};

} /* namespace madness */

#endif /* CC2_H_ */
