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
//#include <string>

// to debug
//#include<chem/mp2.h>

namespace madness {

typedef std::vector<Function<double, 3> > vecfuncT;

static double dipole_x(const coord_3d &r){return r[0];}
static double dipole_y(const coord_3d &r){return r[1];}
static double dipole_z(const coord_3d &r){return r[2];}
static double dipole_r(const coord_3d &r){return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);}


/// Structure that holds the CC intermediates and is able to refresh them
struct CC_Intermediates {
public:
	CC_Intermediates(World&world, const vecfuncT &bra, const vecfuncT &ket,
			const Nemo&nemo, const CC_Parameters &param) :
				world(world),
				parameters(param),
				mo_bra_(bra),
				mo_ket_(ket),
				poisson(std::shared_ptr < real_convolution_3d> (CoulombOperatorPtr(world,parameters.lo,parameters.thresh_poisson_3D))),
				f12op(std::shared_ptr< real_convolution_3d> (SlaterF12OperatorPtr(world,parameters.gamma(),parameters.lo,parameters.thresh_poisson_3D))),
				density_(make_density(bra, CC_vecfunction(ket,HOLE,parameters.freeze))),
				exchange_intermediate_(make_exchange_intermediate(bra, CC_vecfunction(ket,HOLE,parameters.freeze))),
				f12_exchange_intermediate_(make_f12_exchange_intermediate(bra,CC_vecfunction(ket,HOLE,parameters.freeze))),
				hartree_potential_(make_hartree_potential(density_)),
				integrals_hf_(make_two_electron_integrals_hf()) {
		sanity_check();

	}

	void sanity_check()const{
		if(parameters.debug){
			if(world.rank()==0) std::cout << "\n\nCC INTERMEDIATES SANITY CHECK\n\n";
			real_function_3d test_J = real_factory_3d(world);
			real_function_3d test_rho = real_factory_3d(world);
			for(size_t i=0;i<mo_ket_.size();i++){
				test_rho += (mo_bra_[i]*mo_ket_[i]).truncate();
			}
			test_J = ((*poisson)(test_rho)).truncate();
			double diff_J = (hartree_potential_ - test_J).norm2();
			if(diff_J > FunctionDefaults<3>::get_thresh()){
				warning("Hartree Potential Inaccurate " + stringify(diff_J));
				hartree_potential_.print_size("hartree_potential stored");
				test_J.print_size("hartree_potential recalc");
			}else if(world.rank()==0) std::cout << "Hartree Potantial is sane\n";


			real_function_3d tmp00 = (mo_ket_[0]*mo_bra_[0]).truncate();
			real_function_3d testK00 = ((*poisson)(tmp00)).truncate();
			double diff_K00 = (get_EX(0,0) - testK00).norm2();
			if(diff_K00 > FunctionDefaults<3>::get_thresh()){
				warning("Exchange Potential Inaccurate <0|g|0> " + stringify(diff_K00));
			}else if(world.rank()==0) std::cout << "<0|g|0> is sane \n";

			if(mo_ket_.size()>1){
				size_t x=mo_ket_.size()-1;
				real_function_3d tmp0x = (mo_bra_[0]*mo_ket_[x]).truncate();
				real_function_3d testK0x = ((*poisson)(tmp0x));
				double diff_K0x = (get_EX(0,x) - testK0x).norm2();
				if(diff_K0x > FunctionDefaults<3>::get_thresh()){
					warning("Exchange Potential Inaccurate <0|g|"+ stringify(x) +"> " +stringify(diff_K0x));
				}else if(world.rank()==0) std::cout << "<0|g|" << x << "> is sane\n";
			}


			if(world.rank()==0) std::cout << "\n\nCC INTERMEDIATES SANITY CHECK ENDED\n\n";
		}
	}

	/// Get the intermediates
	real_function_3d get_density() {
		return density_;
	}
	real_function_3d get_perturbed_density() const {
		return perturbed_density_;
	}
	real_function_3d get_hartree_potential() const {
		return hartree_potential_;
	}
	real_function_3d get_J() {
		return hartree_potential_;
	}
	real_function_3d get_perturbed_hartree_potential() const {
		return perturbed_hartree_potential_;
	}
	real_function_3d get_pJ() {
		return perturbed_hartree_potential_;
	}
	intermediateT get_EX()const {
		return exchange_intermediate_;
	}

	/// returns <k|g|l>
	real_function_3d get_EX(const size_t &k,const size_t &l)const{return exchange_intermediate_(k,l);}
	real_function_3d get_fEX(const size_t &k,const size_t &l)const{return f12_exchange_intermediate_(k,l);}
	intermediateT get_pEX()const{
		return perturbed_exchange_intermediate_;
	}
	/// returns <k|g|\tau_l>
	real_function_3d get_pEX(const size_t &k, const size_t &l)const{return perturbed_exchange_intermediate_(k,l);}
	real_function_3d get_pfEX(const size_t &k,const size_t &l)const{return perturbed_f12_exchange_intermediate_(k,l);}
	Tensor<double> get_intergrals_hf() const {
		return integrals_hf_;
	}
	Tensor<double> get_integrals_mixed_t1() const {
		return integrals_mixed_t1_;
	}
	Tensor<double> get_integrals_t1() const {
		return integrals_t1_;
	}

	/// refresh the intermediates that depend on the \tau functions
	void update(const CC_vecfunction &tau) {
		if (world.rank() == 0)
			std::cout << "Update Intermediates:\n";
		perturbed_density_ = make_density(mo_bra_, tau);
		perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
		perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,
				tau);
		perturbed_f12_exchange_intermediate_ = make_f12_exchange_intermediate(mo_bra_,
				tau);
		{
			if(world.rank()==0)std::cout << "\n---Updated Intermediates---\n";
			hartree_potential_.print_size("0-Hartree-Potential ");
			perturbed_hartree_potential_.print_size("T1-Hartree-Potential");
			if(world.rank()==0)std::cout << "Exchange Intermediates:\n";
			double size0=0.0;
			double sizet=0.0;
			for(auto x:exchange_intermediate_.allpairs) size0 += get_size(world,vecfuncT(1,x.second));
			for(auto x:perturbed_exchange_intermediate_.allpairs) sizet += get_size(world,vecfuncT(1,x.second));
		}

	}

	/// make a density from two input functions
	/// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
	/// @param[in] vecfuncT_bra
	/// @param[in] vecfuncT_ket
	/// @param[out] \sum_i bra_i * ket_i
	//real_function_3d make_density(const vecfuncT &bra,const vecfuncT &ket) const;
	real_function_3d make_density(const vecfuncT &bra,const CC_vecfunction &ket) const;
	/// Poisson operator
	std::shared_ptr<real_convolution_3d> get_poisson()const{return poisson;}

private:
	World &world;
	const CC_Parameters &parameters;
	const vecfuncT &mo_bra_;
	const vecfuncT &mo_ket_;
	const std::shared_ptr<real_convolution_3d> poisson;
	const std::shared_ptr<real_convolution_3d> f12op;
	/// const intermediates
	const real_function_3d density_;
	/// Exchange intermediate: EX(i,j) = <i|g|j>
	intermediateT exchange_intermediate_;
	/// The f12 exchange intermediate fEX(i,j) = <i|f12|j>
	intermediateT f12_exchange_intermediate_;
	/// Hartree_Potential  = J = \sum_k <k|g|k> = Poisson(density)
	const real_function_3d hartree_potential_;
	/// intermediates that need to be recalculated before every iteration
	/// Perturbed Density = \sum_k |k><\tau_k|
	real_function_3d perturbed_density_;
	/// Perturbed Hartree Poptential PJ = \sum_k <k|g|\tau_k> = Poisson(perturbed_density)
	real_function_3d perturbed_hartree_potential_;
	/// Perturbed Exchange Intermediate: PEX(i,j) = <i|g|\tau_j>
	intermediateT perturbed_exchange_intermediate_;
	/// Perturbed f12-exchange-intermediate: pfEX(i,j) = <i|f12|tau_j>
	intermediateT perturbed_f12_exchange_intermediate_;

	/// Two electron integrals
	/// The Integrals which consist of the hartree-fock ground state
	/// <ij|g|kl> = <ji|g|lk>
	const Tensor<double> integrals_hf_;
	/// The Integrals which consist of ground state and t1 amplitudes
	/// <ij|g|k\tau_l> = <ji|g|\tau_lk>
	Tensor<double> integrals_mixed_t1_;
	/// The Integrals from the t1 functions and the hf orbitals
	/// <ij|g|\tau_k\tau_l> = <ji|g|\tau_l\tau_k>
	Tensor<double> integrals_t1_;

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_INTERMEDIATES!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
				1);
	}
	void warning(const std::string &msg) const {
		std::cout << "\n\n\nWARNING IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
	}
public:
	/// Make the exchange intermediate: EX[j][i] <bra[i](r2)|1/r12|ket[j](r2)>
	intermediateT make_exchange_intermediate(const vecfuncT &bra,
			const CC_vecfunction &ket)const;
	intermediateT make_f12_exchange_intermediate(const vecfuncT &bra,
			const CC_vecfunction &ket)const;
	/// Calculates the hartree potential Poisson(density)
	/// @param[in] density: a 3d function on which the poisson operator is applied (can be the occupied density and the perturbed density)
	/// @param[out] poisson(density) = \int 1/r12 density(r2) dr2
	real_function_3d make_hartree_potential(
			const real_function_3d &density) const {
		real_function_3d hartree = (*poisson)(density);
		hartree.truncate();
		return hartree;
	}

	/// Calculates two electron integrals
	/// <ij|g|kl>
	Tensor<double> make_two_electron_integrals_hf() const;
	/// <ij|g|k\tau_l>
	Tensor<double> make_two_electron_integrals_mixed_t1(
			const vecfuncT &tau) const;
	// <ij|g|\tau_k \tau_l>
	Tensor<double> make_two_electron_integrals_t1(const vecfuncT &tau) const;
};

/// Coupled Cluster Operators (all closed shell)
class CC_Operators {
public:
	/// Constructor
	CC_Operators(World& world, const Nemo &nemo,
			const CorrelationFactor &correlationfactor, const CC_Parameters &param) : Q12(world),
			world(world), nemo(nemo), corrfac(correlationfactor),parameters(param), mo_bra_(
					make_mo_bra(nemo)), mo_ket_(make_mo_ket(nemo)),active_mo_(make_active_mo()),orbital_energies(init_orbital_energies(nemo)), intermediates_(
							world, mo_bra_, mo_ket_, nemo, param){
		// make operators

		// make the active mo vector (ket nemos, bra is not needed for that)
		MADNESS_ASSERT(active_mo_.size()==mo_ket_.size()-parameters.freeze);
		MADNESS_ASSERT(mo_ket_.size()==mo_bra_.size());
		// initialize the Q12 projector
		Q12.set_spaces(mo_bra_,mo_ket_,mo_bra_,mo_ket_);
		performance_S.current_iteration = 0;
		performance_D.current_iteration = 0;
	}



	// collect all the data for every function and every iteration
	mutable CC_performance performance_S;
	mutable CC_performance performance_D;

	// collect all the warnings that are put out over a calculation
	mutable std::vector<std::string> warnings;

	/// save a function
	template<typename T, size_t NDIM>
	void save_function(const Function<T,NDIM>& f, const std::string name) const;

	void plot(const real_function_3d &f, const std::string &msg)const{
		CC_Timer plot_time(world,"plotting " + msg);
		plot_plane(world,f,msg);
		plot_time.info();
	}

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_OPERATORS!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
				1);
	}
	void warning(const std::string &msg,CC_data &data) const{
		std::cout << "\n\n\nWARNING IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		warnings.push_back(msg);
		data.warnings.push_back(msg);
	}
	void warning(const std::string &msg) const{
		std::cout << "\n\n\nWARNING IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		warnings.push_back(msg);
	}

	void output_section(const std::string&msg)const{
		if(world.rank()==0){
			std::cout << "\n\n--------------\n";
			std::cout << msg << std::endl;
			std::cout << "\n";
		}
	}
	void output(const std::string &msg)const{
		if(world.rank()==0){
			std::cout << msg << std::endl;
		}
	}

	void update_intermediates(const CC_vecfunction &singles){
		CC_Timer update(world,"Update Intermediates");
		intermediates_.update(singles);
		update.info();
	}

	StrongOrthogonalityProjector<double,3> Q12;

	real_function_3d mo_ket(const size_t &i)const{
		return mo_ket_[i];
	}
	vecfuncT mo_ket()const{return mo_ket_;}
	real_function_3d mo_bra(const size_t &i)const{
		return mo_bra_[i];
	}
	vecfuncT mo_bra()const{return mo_bra_;}


	/// makes the t intermediate which is defined as: |t_i> = |\tau_i> + |i>
	CC_vecfunction make_t_intermediate(const CC_vecfunction &tau)const{
		CC_vecfunction result;
		for(auto x:tau.functions){
			real_function_3d tmpx = x.second.function + mo_ket_[x.second.i];
			CC_function tmpi(tmpx,x.second.i,MIXED);
			result.insert(x.second.i,tmpi);
		}
		return result;
	}

	vecfuncT get_CCS_potential(const CC_vecfunction &singles)const{

		// make a dummy doubles with no content
		Pairs<CC_Pair> doubles;

		vecfuncT result = potential_singles(doubles,singles,_reF3D_);

		result = add(world,result,potential_singles(doubles,singles,_S1_));  // brillouin term
		result = add(world,result,potential_singles(doubles,singles,_S5a_)); // brillouin term
		result = add(world,result,potential_singles(doubles,singles,_S3c_));
		result = add(world,result,potential_singles(doubles,singles,_S5b_));
		result = add(world,result,potential_singles(doubles,singles,_S5c_));
		result = add(world,result,potential_singles(doubles,singles,_S6_));

		Q(result);
		truncate(world,result);
		performance_S.current_iteration++;
		return result;
	}

	vecfuncT get_CC2_singles_potential(const CC_vecfunction &singles, const Pairs<CC_Pair> &doubles){
		vecfuncT fock_residue = potential_singles(doubles,singles,_reF3D_);

		vecfuncT result = potential_singles(doubles,singles,_S3c_);
		result = add(world,result,potential_singles(doubles,singles,_S5b_));
		result = add(world,result,potential_singles(doubles,singles,_S5c_));
		result = add(world,result,potential_singles(doubles,singles,_S6_));
		result = add(world,result,potential_singles(doubles,singles,_S2b_));
		result = add(world,result,potential_singles(doubles,singles,_S2c_));
		result = add(world,result,potential_singles(doubles,singles,_S4a_));
		result = add(world,result,potential_singles(doubles,singles,_S4b_));
		result = add(world,result,potential_singles(doubles,singles,_S4c_));

		Q(result);
		truncate(world,result);
		// need to store this for application of Fock oerator on singles ( F|taui> = Singles_potential[i] + \epsilon_i|taui>)
		current_singles_potential = result;
		result = add(world,result,fock_residue);
		Q(result);
		truncate(world,result);
		performance_S.current_iteration++;
		return result;
	}

	// only get the part of the singles that is produced exclusively by the doulbes in order to make a first guess for the singles
	vecfuncT get_CC2_singles_initial_potential(const Pairs<CC_Pair> &doubles)const{
		// make_zero guess
//		real_function_3d zeroguess = real_factory_3d(world);
//		vecfuncT tmp(mo_ket_.size(),zeroguess);
//		CC_vecfunction singles(tmp,PARTICLE,parameters.freeze,tmp.size());
//		MADNESS_ASSERT(singles.size()==mo_ket_.size()-parameters.freeze);

		vecfuncT result = zero_functions<double,3>(world,mo_ket_.size()-parameters.freeze);
//		{CC_Timer timer_S2b(world,"Singles Potential: S2b+X");
//		result =potential_singles(doubles,singles,_S2b_);
//		for(size_t i=0;i<result.size();i++) result[i].print_size("S2b_"+stringify(i));
//		timer_S2b.info();}
//		{CC_Timer timer_S2c(world,"Singles Potential: S2c+X");
//		vecfuncT s2c = potential_singles(doubles,singles,_S2c_);
//		for(size_t i=0;i<result.size();i++) s2c[i].print_size("S2c_"+stringify(i));
//		result = add(world,s2c,result);
//		timer_S2c.info();}
		return result;
	}

	real_function_6d get_CC2_doubles_from_singles_potential(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		real_function_6d result = real_factory_6d(world);
		const std::string output = "DoPo:";
		result = potential_doubles(taui,tauj,singles,_D6b_D8b_D9_);
		result += potential_doubles(taui,tauj,singles,_D4b_D6c_D8a_);
		if(world.rank()==0) performance_D.info(performance_D.current_iteration);
		return result;
	}

	real_function_6d get_CC2_doubles_potential(const vecfuncT &singles, const CC_Pair &u)const{
		MADNESS_EXCEPTION("CC2 doubles potential not yet implemented",1);
		return real_factory_6d(world);
	}

	real_function_6d get_MP2_potential_constant_part(const CC_Pair &u)const{
		CC_function mo_i(mo_ket_[u.i],u.i,HOLE);
		CC_function mo_j(mo_ket_[u.j],u.j,HOLE);
		CC_Timer timer_U(world,"Ue(R)|ij>");
		real_function_6d UePart = apply_transformed_Ue(mo_ket_[u.i],mo_ket_[u.j],u.i,u.j);
		UePart.print_size("Ue|"+stringify(u.i)+stringify(u.j)+">");
		timer_U.info();

		CC_Timer timer_KffK(world,"Kf|ij>");
		real_function_6d KffKPart = apply_exchange_commutator(mo_i,mo_j);
		KffKPart.print_size("[K,f]|"+stringify(u.i)+stringify(u.j)+">");
		timer_KffK.info();

		real_function_6d unprojected_result = (UePart - KffKPart).truncate();
		unprojected_result.print_size("Ue - [K,f]|"+stringify(u.i)+stringify(u.j)+">");
		CC_Timer timer_Q(world,"Apply Q12");
		real_function_6d result = Q12(unprojected_result);
		result.print_size("Q12(Ue - [K,f]|"+stringify(u.i)+stringify(u.j)+">)");
		timer_Q.info();
		return result;

	}

	/// returns the non constant part of the MP2 potential which is
	/// (2J-K+Un)|uij>
	real_function_6d get_MP2_potential_residue(const CC_Pair &u)const{
		CC_Timer timer(world,"(2J-K(R)+Un)|uij>");
		CC_data data("mp2_residue");
		real_function_6d result = fock_residue_6d(u);
		data.result_size=get_size(result);
		data.result_norm=result.norm2();
		data.time = timer.current_time();
		performance_D.insert(data.name,data);
		timer.info();
		return result;
	}

	/// reconstructs the full pair function from the regularized pair functions
	/// used to compute norms of the doubles to compare with LCAO codes
	/// used to debug the singles potential
	/// @param[in]: u the regularized function
	/// @param[out]: tau = u + Q12f12(|ij> + |taui,j> + |i,tauj> + |taui,tauj>) = u + Q12f12|titj> with ti = taui + i
	real_function_6d make_full_pair_function(const CC_Pair &u,const  CC_function &taui, const CC_function &tauj)const{
		const size_t i=u.i;
		const size_t j=u.j;
		MADNESS_ASSERT(i==taui.i);
		MADNESS_ASSERT(j==tauj.i);
		real_function_3d ti = mo_ket_[i] + taui.function;
		real_function_3d tj = mo_ket_[j] + tauj.function;

		real_function_6d Q12f12titj = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(ti)).particle2(copy(tj));
		Q12f12titj.fill_tree().truncate().reduce_rank();
		apply_Q12(Q12f12titj);

		real_function_6d result = u.function + Q12f12titj;
		return result;
	}

	// right now this is all copied from mp2.cc
	double compute_mp2_pair_energy(CC_Pair &pair)const;



	/// Projectors to project out the occupied space
	// 3D on vector of functions
	void Q(vecfuncT &f) const {
		for (size_t i = 0; i < f.size(); i++)
			Q(f[i]);
	}
	// 3D on single function
	// use the projector class, like in Q12
	void Q(real_function_3d &f) const {
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			f -= mo_bra_[i].inner(f) * mo_ket_[i];
		}
	}

	/// CCSD/CC2 singles potential parts

	/// Genereal function which evaluates a CC_singles potential
	vecfuncT potential_singles(const Pairs<CC_Pair> u, const CC_vecfunction & singles , const potentialtype_s &name) const {
		CC_Timer timer(world,assign_name(name));
		CC_data data(name);
		vecfuncT result;

		switch(name){
		case _reF3D_ :
			result =fock_residue_closed_shell(singles);
			break;
		case _S3c_ :
			result = add(world,S3c(singles),S3c_X(singles));
			break;
		case _S5b_ :
			result = add(world,S5b(singles),S5b_X(singles));
			break;
		case _S5c_ :
			result = add(world,S5c(singles),S5c_X(singles));
			break;
		case _S6_  :
			result = add(world,S6(singles),S6(singles));
			break;
		case _S2b_ :
			result = add(world,S2b_3D_part(u,singles,data),S2b_6D_part(u,singles,data));
			break;
		case _S2c_ :
			result = add(world,S2b_3D_part(u,singles,data),S2b_6D_part(u,singles,data));
			break;
		case _S4a_ :
			result = add(world,S4a_3D_part(singles,data),S4a_6D_part(u,singles,data));
			break;
		case _S4b_ :
			result = add(world,S4b_3D_part(u,singles,data),S4b_6D_part(u,singles,data));
			break;
		case _S4c_ :
			result = add(world,S4c_3D_part(u,singles,data),S4c_6D_part(u,singles,data));
			break;
		case _S1_ :
			result = S1(singles);
			break;
		case _S5a_ :
			result = S5a(singles);
			break;
		}
		truncate(world,result);
		Q(result);
		size_t num=0;
		if(world.rank()==0) std::cout << "result for singles potential " << assign_name(name) << std::endl;
		for(auto i:result){
			if(world.rank()==0) std::cout << "||tau" << num << "||=" << i.norm2() << std::endl;
			num++;
		}
		data.result_size=get_size(world,result);
		data.result_norm=norm2(world,result);
		data.time = timer.current_time();
		performance_S.insert(data.name,data);
		return result;
	}

	real_function_6d potential_doubles(const CC_function &taui, const CC_function &tauj,const  CC_vecfunction &singles, const potentialtype_d &name)const{
		CC_Timer timer(world,assign_name(name));
		CC_data data(assign_name(name));
		output("Now Doing "+assign_name(name)+ " \n\n");

		real_function_6d result = real_factory_6d(world);

		switch(name){
		case _D4b_ :
			result = G_D4b_decomposed(taui,tauj,singles);
			break;
		case _D6b_ :
			result = G_D6b_decomposed(taui,tauj,singles);
			break;
		case _D6c_ :
			result = G_D6c(taui,tauj,singles);
			break;
		case _D8a_ :
			result = G_D8a(taui,tauj,singles);
			break;
		case _D8b_ :
			result = G_D8b(taui,tauj,singles);
			break;
		case _D9_  :
			result = G_D9(taui,tauj,singles);
			break;
		case _D6b_D8b_D9_ :
			result = D6b_D8b_D9(taui,tauj,singles);
			break;
		case _D4b_D6c_D8a_ :
			result = D4b_D6c_D8a(taui,tauj,singles);
			break;

		case _reF6D_ :
			error("reF6D not yet supported by general doubles potential function");
			break;
		case _reCC2_ :
			error("reCC2 not yet supported by general doubles potential function");
			break;
		}
		result.print_size(assign_name(name));
		output("Finished with "+assign_name(name));
		data.result_norm = result.norm2();
		data.result_size = get_size(result);
		data.time=(timer.current_time());
		performance_D.insert(data.name,data);
		return result;
	}


	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_i = \sum_k <k|r12|k> |tau_i>
	// K_i = \sum_k <k|r12|tau_i> |k>
	vecfuncT fock_residue_closed_shell(const CC_vecfunction &tau) const;

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3c(const CC_vecfunction &tau) const;

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = -Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3c_X(const CC_vecfunction &tau) const;

	/// The S5b term
	//[i]    [Q]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// 2\sum_k <k|g|\tau_k> |\tau_i>
	// No Q is applied yet !
	vecfuncT S5b(const CC_vecfunction &tau) const;

	/// The S5b Exchange Term
	//[i]         [Q]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_k <k|g|\tau_i> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5b_X(const CC_vecfunction &tau) const;

	/// The S5c term
	//[Q]    [i]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// -2\sum_kl <kl|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	// May use alteriative algorithm with perturbed density intermediate
	vecfuncT S5c(const CC_vecfunction &tau) const;

	/// The S5c_X echange term
	//[Q]         [i]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_kl <lk|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5c_X(const CC_vecfunction &tau) const;

	/// The S6+X Term
	// \    /\    /...
	//  \  /  \  /   /\
	//  _\/_  _\/_  _\/_
	// -Q \sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>
	// Q is not applied yet!
	vecfuncT S6(const CC_vecfunction  &tau) const;


	// The two brillouin terms S1 and S5a of the singles potential
	vecfuncT S1(const CC_vecfunction &tau)const{
		vecfuncT result;
		for(auto tmpi:tau.functions){
			CC_function& i = tmpi.second;
			real_function_3d resulti = real_factory_3d(world);
			resulti = apply_F(CC_function(mo_ket_[i.i],i.i,UNDEFINED)); // undefined for the testing case where the mos are not converged
			Q(resulti);
			if(parameters.debug){
				real_convolution_3d G = BSHOperator<3>(world, sqrt(-2.0 * get_orbital_energies()[i.i]), parameters.lo, parameters.thresh_bsh_3D);
				resulti.print_size("QF"+stringify(i.i));
				plot_plane(world,resulti,"QF"+stringify(i.i));
				real_function_3d GQFi = G(resulti);
				GQFi.print_size("GQF"+stringify(i.i));
				plot_plane(world,GQFi,"GQF"+stringify(i.i));

			}
			result.push_back(resulti);
		}
		return result;
	}

	vecfuncT S5a(const CC_vecfunction &tau)const{
		vecfuncT result;
		for(auto tmpi:tau.functions){
			CC_function& i=tmpi.second;
			real_function_3d resulti = real_factory_3d(world);
			for(auto tmpk:tau.functions){
				CC_function& k=tmpk.second;
				real_function_3d tmp = apply_F(CC_function(i.function,i.i,UNDEFINED)); // undefined for the test case where the moi are not converged yet
				const double a = mo_bra_[k.i].inner(tmp);
				resulti += a*k.function;
			}
			result.push_back(resulti);
		}
		Q(result);
		return result;
	}


	/// CC2 singles diagrams with 6d functions as input
	/// Use GFInterface in function_interface.h as kernel (f*g) and do not reconstruct \tau = f12u(1,2) if possible
	/// Since the correlation factor of CC2 has Slater form like in MP2: g12f12 = g12(1-exp(-mu*r12)/r12) = g12 - exp(-mu*r12)/r12 = Coulomb_Operator - BSH_Operator(mu)

	/// S2b + X Term
	// [i]   [Q]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// @param[in] Structure which holds all current CC pair functions
	/// @param[in] Structure which holds all current CC single excitations
	/// @param[out] -Q\sum_k \left( 2<k|g|u_ik> - <k|g|u_ki> + 2<k|gQf|t_it_k> - <k|gQf|t_kt_i> \right), with t_i = i + \tau_i
	/// notation: <k|g|u_ik> = <k(2)|g12|u_ik(1,2)> (Integration over second particle)
	vecfuncT S2b(const Pairs<CC_Pair> u, const CC_vecfunction &singles) const;
	vecfuncT S2b_3D_part(const Pairs<CC_Pair> u, const CC_vecfunction &singles,CC_data &data) const;
	vecfuncT S2b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction &singles,CC_data &data) const;

	/// S2c + X Term
	// [Q]   [i]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// = -Q\sum_{kl}\left( 2<k|lgi|ulk> - <l|kgi|u_{lk}> + 2<k|lgi|t_k>*|t_l> - 2<l|kgi|t_k>*|t_l? \right)
	/// Notation: 6D Integration over second particle, intermediates: lgi = <l|g|i> is the exchange intermediate
	/// Notation: t are the t-intermediates: |t_i> = |i> + |\tau_i>
	/// @param[in] All the current coupled cluster Pairs
	/// @param[in] The coupled cluster singles
	/// @param[out] the S2c+X Potential
	vecfuncT S2c(const Pairs<CC_Pair> &u, const CC_vecfunction &singles) const;
	vecfuncT S2c_3D_part(const Pairs<CC_Pair> &u, const CC_vecfunction &singles, CC_data &data) const;
	vecfuncT S2c_6D_part(const Pairs<CC_Pair> &u, const CC_vecfunction &singles, CC_data &data) const;
	/// The S4a + X diagram
	//[Q]       [i]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k>
	vecfuncT S4a(const Pairs<CC_Pair> u, const CC_vecfunction & tau) const;
	vecfuncT S4a_3D_part(const CC_vecfunction & tau,CC_data &data) const;
	vecfuncT S4a_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,CC_data &data) const;

	/// The S4b
	//[i]       [Q]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum_{kl} (2<k(3)l(4)|g34f14|\tau_{i}(3)u_{kl}(1,4)>  // exchange part - <k(4)l(3)|g34f14|\tau_i(3)u_{lk}(1,4)>)
	// 1. make exchange intermedaite X(4) = <k(3)|g34|\tau_i(3)>_3 *  l(4)			Exchange part : Xx(4) = <l(3)|g34|\tau_i(3)>(4) * k(4)
	// 2. make 6d intermediate Y(1,4) = X(4)* u_{kl}(1,4)							Exchange part : Yx(1,4) = X(4)*u_{lk}(1,4)
	// 3. make f14 integration via delta function trick: result(1) = \int f14 Y(1,4) d4 = \int delta(5-1) (\int f54 Y(1,4) d4)d5
	// 3.1 do the convolution Z(1,5) = \int f54 Y(1,4) d4							Exchange part: Zx(1,5) = int f54 Yx(1,4)d4
	// 3.2 project out the unit function: result(1) = <I(5)|Z(1,5)>_5				Exchange part: resultx(1) = <I(5)|Zx(1,5>_5
	vecfuncT S4b(const Pairs<CC_Pair> u, const CC_vecfunction & tau) const;
	vecfuncT S4b_3D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,CC_data &data) const;
	vecfuncT S4b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,CC_data &data) const;

	/// The S4c + X + X + X + X Diagrams
	//            [i]   [Q]
	//   .......   \    /
	//  /\     /\   \  /
	// _\/_   _\/____\/_
	/// Q\sum_{kl}[ 4*<k(3)l(4)|g34 f14| \tau_k(3) u_{il}(1,4)> - 2* <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)>
	/// - 2* <k(3)l(4)|g34 f14| \tau_k(3) U_{li}(1,4)> + <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)>  ]
	// First and third Terms are solved like this:
	// 1. X(4) = \sum_k (<k(3)|g34|\tau_k(3)>_3(4)) * l(4) = perturbed_hartree_potential(4) * l(4)
	// 2. Y(1,4) = X(4) u_{il}(1,4)			Exchange Part: Yx(4,1) = X(4) u_{li}(4,1)
	// 3.1 Z(1,5) = \int f54 Y(1,4) d4		Exchange Part: Zx(5,1) = \int f54 Yx(4,1) d4
	// 3.2 result(1) = -4 <I(5)|Z(1,5)>_5 -2 <I(5)|Zx(1,5)>_5
	// Second and fourth terms can not use the perturbed hartree potential
	vecfuncT S4c(const Pairs<CC_Pair> u, const CC_vecfunction & tau) const;
	vecfuncT S4c_3D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,CC_data &data) const;
	vecfuncT S4c_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,CC_data &data) const;

	// CC2 Doubles Potential

	/// smalll helper function to track the time for the projetor
	void apply_Q12(real_function_6d &f,const std::string &msg = "6D-function")const{
		CC_Timer Q12_time(world,"Applying Q12 to "+ msg);
		f=Q12(f);
		Q12_time.info();
	}


	/// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
	/// @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
	/// @param[in] \tau_j
	/// @param[in] u, the uij pair structure which holds the consant part of MP2
	/// @param[out] Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
	/// Right now Calculated in the decomposed form: |titj> = |i,j> + |\taui,\tauj> + |i,\tauj> + |\taui,j>
	/// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
	real_function_6d make_cc2_residue(const CC_function &taui, const CC_function &tauj, const CC_Pair &u)const;



	// apply the kinetic energy operator to a decomposed 6D function
	/// @param[in] a 3d function x (will be particle 1 in the decomposed 6d function)
	/// @param[in] a 3d function y (will be particle 2 in the decomposed 6d function)
	/// @param[out] a 6d function: G(f12*T*|xy>)
	real_function_6d make_GQfT_xy(const real_function_3d &x, const real_function_3d &y, const size_t &i, const size_t &j)const;


	// MP2 Terms are
	// fock_residue_6d = (2J - Kn + Un) |u_{ij}> , KR = R12^{-1}*K*R12 (nuclear tranformed K)
	// Uen|ij> = R12{-1}*U_e*R12 |ij>

	/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
	real_function_6d fock_residue_6d(const CC_Pair &u) const;

	real_function_6d G_fock_residue_xy(const CC_function &taui, const CC_function &tauj)const{
		error("G_fock_residue_xy .... this function should not be used, if so check it");
		const size_t i=taui.i;
		const size_t j=tauj.i;
		const real_function_3d & x = taui.function;
		const real_function_3d & y = tauj.function;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		CC_Timer local_time(world,"Fock-residue-xy-local-part");
		// make x2 = (2.0*J + U2)x, and the same for y
		real_function_3d x2 = ((2.0*intermediates_.get_hartree_potential() + nemo.nuclear_correlation->U2())*x).truncate();
		real_function_3d y2 = ((2.0*intermediates_.get_hartree_potential() + nemo.nuclear_correlation->U2())*y).truncate();

		// Apply the Greens Function on the local parts
		// G(|x2y> + |xy2>)
		real_function_6d local_part;
		{
			real_function_6d x2y = make_f_xy(CC_function(x2,i,UNDEFINED),tauj); //CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x2)).particle2(copy(y));
			real_function_6d xy2 = make_f_xy(taui,CC_function(y2,j,UNDEFINED)); //CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(y2));
			apply_Q12(x2y,"x2y");
			apply_Q12(xy2,"xy2");
			real_function_6d Gx2y = G(x2y);
			real_function_6d Gxy2 = G(xy2);
			local_part = Gx2y + Gxy2;
		}
		local_time.info();

		CC_Timer unuc_time(world,"Fock-residue-xy-U1-part");
		real_function_6d U1_part;
		{
			std::vector < std::shared_ptr<real_derivative_3d> > gradop;
			gradop = gradient_operator<double, 3>(world);
			for(size_t axis=0;axis<3; ++ axis){
				real_function_3d dx = (*gradop[axis])(x);
				real_function_3d dy = (*gradop[axis])(y);
				real_function_3d U1 = nemo.nuclear_correlation->U1(axis);
				real_function_3d Udx = (U1*dx).truncate();
				real_function_3d Udy = (U1*dy).truncate();

				real_function_6d fUdxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Udx)).particle2(copy(y));
				real_function_6d fxUdy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(Udy));
				CC_Timer fill_tree_timer_1(world,"fill_tree_1");
				fUdxy.fill_tree(G).truncate().reduce_rank();
				fill_tree_timer_1.info();
				CC_Timer fill_tree_timer_2(world,"fill_tree_2");
				fxUdy.fill_tree(G).truncate().reduce_rank();
				fill_tree_timer_2.info();
				apply_Q12(fUdxy,"fUdxy");
				apply_Q12(fxUdy,"fxUdy");
				real_function_6d GQfUdxy = G(fUdxy);
				real_function_6d GQfxUdy = G(fxUdy);
				U1_part += (GQfUdxy + GQfxUdy);
			}
		}
		unuc_time.info();

		CC_Timer unuc_faster(world,"Fock-residue-xy-U1-part-faster");
		real_function_3d U1_part_x = real_factory_3d(world);
		real_function_3d U1_part_y = real_factory_3d(world);
		{
			std::vector < std::shared_ptr<real_derivative_3d> > gradop;
			gradop = gradient_operator<double, 3>(world);
			for(size_t axis=0;axis<3; ++ axis){
				real_function_3d dx = (*gradop[axis])(x);
				real_function_3d dy = (*gradop[axis])(y);
				real_function_3d U1 = nemo.nuclear_correlation->U1(axis);
				real_function_3d Udx = (U1*dx).truncate();
				real_function_3d Udy = (U1*dy).truncate();

				U1_part_x += Udx;
				U1_part_y += Udy;
			}
		}
		real_function_6d fUdxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(U1_part_x)).particle2(copy(y));
		real_function_6d fxUdy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(U1_part_y));
		CC_Timer fill_tree_timer_1(world,"fill_tree_1");
		fUdxy.fill_tree(G).truncate().reduce_rank();
		fill_tree_timer_1.info();
		CC_Timer fill_tree_timer_2(world,"fill_tree_2");
		fxUdy.fill_tree(G).truncate().reduce_rank();
		fill_tree_timer_2.info();
		apply_Q12(fUdxy,"fUdxy");
		apply_Q12(fxUdy,"fxUdy");
		real_function_6d GQfUdxy = G(fUdxy);
		real_function_6d GQfxUdy = G(fxUdy);
		real_function_6d U1_part_faster = GQfUdxy + GQfxUdy;
		unuc_faster.info();

		// compare two U1 results
		{
			real_function_6d diff = U1_part - U1_part_faster;
			diff.print_size("U1_part - U1_part_faster");
		}

		CC_Timer ex_time(world,"Fock-residue-xy-K-part");
		real_function_6d K_part;
		{
			CC_Timer make_exim_time(world,"Make Exchange Intermedaites for fock residue");
			real_function_3d Kx=real_factory_3d(world);
			real_function_3d Ky=real_factory_3d(world);
			vecfuncT mo_bra_x = mul(world,x,mo_bra_);
			vecfuncT mo_bra_y = mul(world,y,mo_bra_);
			vecfuncT mo_bra_g_x = apply(world,*poisson,mo_bra_x);
			vecfuncT mo_bra_g_y = apply(world,*poisson,mo_bra_y);
			make_exim_time.info();
			for(size_t k=0;k<mo_bra_.size();k++){
				Kx += mo_bra_g_x[k]*mo_ket_[k];
				Ky += mo_bra_g_y[k]*mo_ket_[k];
			}
			real_function_6d fKxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Kx)).particle2(copy(y));
			real_function_6d fxKy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(Ky));
			CC_Timer fill_tree_timer_1(world,"fill_tree_1");
			fKxy.fill_tree(G).truncate().reduce_rank();
			fill_tree_timer_1.info();
			CC_Timer fill_tree_timer_2(world,"fill_tree_2");
			fxKy.fill_tree(G).truncate().reduce_rank();
			fill_tree_timer_2.info();
			apply_Q12(fKxy,"fKxy");
			apply_Q12(fxKy,"fxKy");
			real_function_6d GQfKxy = G(fKxy);
			real_function_6d GQfxKy = G(fxKy);
			K_part = (GQfKxy + GQfxKy).truncate();
		}
		ex_time.info();

		real_function_6d result = (local_part + U1_part - K_part).truncate();
		return result;

	}

	/// Echange Operator on 3D function
	/// !!!!Prefactor (-1) is not included
	real_function_3d K(const CC_function &f)const;



	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// if i==j in uij then the symmetry will be exploited
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d K(const real_function_6d &u,
			const bool symmetric = false,const double & thresh = FunctionDefaults<6>::get_thresh()) const;

	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
	/// 1. X(3,2) = bra_k(3)*u(3,2)
	/// 2. Y(1,2) = \int X(3,2) g13 d3
	/// 3. result = Y(1,2)*ket_k(1)
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d apply_K(const real_function_6d &u,
			const size_t &particle,const double & thresh = FunctionDefaults<6>::get_thresh()) const;

	/// returns \sum_k (<k|g|f> *|k>).truncate()
	real_function_3d apply_K(const CC_function &f)const{
		real_function_3d result = real_factory_3d(world);
		{
			if(f.type == HOLE){
				for(size_t k=0;k<mo_ket_.size();k++) result += (intermediates_.get_EX(k,f.i)*mo_ket_[k]).truncate();
			}
			else if(f.type == PARTICLE){
				for(size_t k=0;k<mo_ket_.size();k++) result += (intermediates_.get_pEX(k,f.i)*mo_ket_[k]).truncate();
			}
			else{
				for(size_t k=0;k<mo_ket_.size();k++) result += (((*poisson)(mo_bra_[k]*f.function).truncate())*mo_ket_[k]).truncate();
			}
		}
		return result;
	}

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
	/// @param[out]  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
	real_function_6d apply_transformed_Ue(const real_function_3d x,
			const real_function_3d y, const size_t &i, const size_t &j) const;

	/// Apply the Exchange Commutator [K,f]|xy>
	real_function_6d apply_exchange_commutator(const CC_function &x, const CC_function &y)const;

	/// Apply the Exchange operator on a tensor product multiplied with f12
	/// !!! Prefactor of (-1) is not inclued in K here !!!!
	real_function_6d apply_Kf(const CC_function &x, const CC_function &y) const;

	/// Apply fK on a tensor product of two 3D functions
	/// fK|xy> = fK_1|xy> + fK_2|xy>
	/// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
	/// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
	real_function_6d apply_fK(const CC_function &x, const CC_function &y) const;


	real_function_3d apply_F(const CC_function &x)const{

		if(x.type == HOLE){
			return get_orbital_energies()[x.i]*x.function;
		}else if(x.type == PARTICLE){
			real_function_3d singles_potential = current_singles_potential[x.i-parameters.freeze];
			return (get_orbital_energies()[x.i]*x.function - singles_potential);
		}else if(x.type == MIXED){
			real_function_3d singles_potential = current_singles_potential[x.i-parameters.freeze];
			return (get_orbital_energies()[x.i]*x.function - singles_potential); // for mixed: eps(i)*x.i = epsi*(moi + taui)
		}else if(x.type == UNDEFINED){
		real_function_3d refined_x = copy(x.function).refine();
		// kinetic part
		CC_Timer T_time(world,"apply_T");
		std::vector < std::shared_ptr<real_derivative_3d> > gradop;
		gradop = gradient_operator<double, 3>(world);
		vecfuncT gradx, laplacex;
		for(size_t axis=0;axis<3;axis++){
			real_function_3d gradxi = (*gradop[axis])(refined_x);
			gradxi.refine();
			gradx.push_back(gradxi);
			real_function_3d grad2xi = (*gradop[axis])(gradxi);
			laplacex.push_back(grad2xi);
		}
		real_function_3d laplace_x = laplacex[0]+laplacex[1]+laplacex[2];
		real_function_3d Tx = laplace_x.scale(-0.5).truncate();
		T_time.info();

		CC_Timer J_time(world,"apply_J");
		real_function_3d Jx = (intermediates_.get_hartree_potential()*x.function).truncate();
		J_time.info();

		CC_Timer K_time(world,"apply_K");
		real_function_3d Kx = K(x);

		CC_Timer U_time(world,"apply_U");
		real_function_3d U2x = (nemo.nuclear_correlation->U2()*x.function).truncate();
		real_function_3d U1x = real_factory_3d(world);
		for(size_t axis=0;axis<3;axis++){
			const real_function_3d U1_axis = nemo.nuclear_correlation->U1(axis);
			const real_function_3d dx = gradx[axis];
			U1x += (U1_axis*dx).truncate();
		}
		U_time.info();

		return (Tx + 2.0*Jx - Kx + U2x + U1x);
		}
		error("apply_F: should not end up here");
		return real_factory_3d(world);
	}

	/// little helper function to pack a vector of CC_3D_functions (just structures which hold the function the index and the type)
	std::vector<CC_function> make_CC_3D_function(const vecfuncT &f, const functype &type){
		std::vector<CC_function> result(f.size());
		for(size_t i=0;i<f.size();i++){
			CC_function tmp(f[i],i,type);
			result[i] = tmp;
		}
		return result;
	}

	// gives back \epsilon_{ij} = \epsilon_i + \epsilon_j
	double get_epsilon(const size_t &i, const size_t &j) const{return (orbital_energies[i]+orbital_energies[j]);}
	// gives back the orbital energies
	std::vector<double> get_orbital_energies()const{return orbital_energies;}
	/// swap particles 1 and 2

	/// param[in]	f	a function of 2 particles f(1,2)
	/// return	the input function with particles swapped g(1,2) = f(2,1)
	real_function_6d swap_particles(const real_function_6d& f) const;

	// Calculate the CC2 energy equation which is
	// \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
	// with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
	double get_CC2_correlation_energy() const;
	double compute_ccs_correlation_energy(const CC_function &taui, const CC_function &tauj)const;
	double compute_cc2_pair_energy(const CC_Pair &u,
			const CC_function &taui, const CC_function &tauj) const;
	/// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
	// the bra elements are always the R2orbitals
	// the ket elements can be \tau_i , or orbitals dependet n the type given
	double make_ij_gQf_ij(const size_t &i, const size_t &j,CC_Pair &u)const;
	double make_ijgQfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	double make_ijgfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	/// Make two electron integral (expensive without intermediates) use just for debugging
	double make_ijgxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	double make_integral(const size_t &i, const size_t &j, const CC_function &x, const CC_function&y)const{
		if(x.type == HOLE){
				real_function_3d igx_y = (intermediates_.get_EX(i,x.i)*y.function).truncate();
				return mo_bra_[j].inner(igx_y);
		}else if(x.type == PARTICLE){
			if(y.type == HOLE){
				real_function_3d jgy_x = (intermediates_.get_EX(j,y.i)*x.function).truncate();
				return mo_bra_[i].inner(jgy_x);
			}else if(y.type == PARTICLE){
				real_function_3d jgy_x = (intermediates_.get_pEX(j,y.i)*x.function).truncate();
				return mo_bra_[i].inner(jgy_x);
			}
		}
		else if (x.type == MIXED or y.type == MIXED){
			real_function_3d igx = (*poisson)((mo_bra_[i]*x.function).truncate());
			double result = mo_bra_[j].inner(igx*y.function);
			return result;
		}
		else if(x.type == UNDEFINED or y.type == UNDEFINED){
			real_function_3d igx = (*poisson)((mo_bra_[i]*x.function).truncate());
			double result = mo_bra_[j].inner(igx*y.function);
			return result;
		}
		else{
			error("ERROR in make_integrals ... should not end up here");
			return 0.0;
		}
		error("ERROR in make_integrals ... should not end up here");
		return 0.0;
	}
	/// Make two electron integral with the pair function
	double make_ijgu(const size_t &i, const size_t &j, const CC_Pair &u)const;
	double make_ijgu(const size_t &i, const size_t &j, const real_function_6d &u)const;
	/// Make two electron integral with BSH operator
	double make_ijGu(const size_t &i, const size_t &j, const CC_Pair &u)const;
	/// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
	/// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
	real_function_3d apply_gf(const real_function_3d &f)const;
	real_function_6d apply_gf(const real_function_6d &f,const size_t &particle)const;

	/// @param[in] x: Function which is convoluted with
	/// @param[in] y: function over which is not integrated
	/// @param[in] z: function which is correlated with y over Q12f12
	/// @param[out] <x(2)|Q12f12|y(1)z(2)>_2
	/// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
	/// 1. <x|f12|z>*|y>
	/// 2. -\sum_m <x|m> <m|f12|z>*|y>
	/// 3. -\sum_n <nx|f12|zy> * |n>
	/// 4. +\sum_{mn} <g|n> <mn|f12|yz> * |m>
	real_function_3d convolute_x_gQf_yz(const real_function_3d &x, const real_function_3d &y, const real_function_3d &z)const;
	/// Description: Similar to convolute x_gQf_yz just without coulomb operator
	real_function_3d convolute_x_Qf_yz(const CC_function &x, const CC_function &y, const CC_function &z)const{

		real_function_3d xz = (x.function*z.function).truncate();

		// the unprojected part <x(2)|f|z(2)> |y(1)>
		real_function_3d unprojected_part = ((*f12op)(xz)*y.function).truncate();

		// the O1 part : \sum_m <x(2)|m(1)><m(1)|f12|y(1)z(2)> = <x(2)|mfy(2)|z(2)> |m(1)>
		real_function_3d O1_part=real_factory_3d(world);
		for(size_t m=0;m<mo_ket_.size();m++){
			real_function_3d mfy;
			if(y.type == HOLE) mfy = intermediates_.get_fEX(m,y.i);
			else if(y.type == PARTICLE) mfy = intermediates_.get_pfEX(m,y.i);
			else mfy = (*f12op)(mo_bra_[m]*y.function).truncate();
			double a = mfy.inner(xz);
			O1_part += a*mo_ket_[m];
		}

		// the O2 part: \sum_n <x(2)|n(2)><n(2)|f12|z(2)> |y(1)>
		real_function_3d O2_part = real_factory_3d(world);
		for(size_t n=0;n<mo_ket_.size();n++){
			if(x.type == PARTICLE) break;
			real_function_3d nfz;
			if(z.type == HOLE) nfz = intermediates_.get_fEX(n,z.i);
			else if(z.type == PARTICLE) nfz = intermediates_.get_pfEX(n,z.i);
			else nfz = (*f12op)(mo_bra_[n]*z.function).truncate();
			double xn = x.function.inner(mo_ket_[n]);
			O2_part += xn*(nfz*y.function).truncate();
		}

		// the O1O2 part: \sum_mn <x(2)|n(2)> <m(1)n(2)|f12|y(1)z(2)> |m(1)>
		real_function_3d O1O2_part = real_factory_3d(world);
		for(size_t n=0;n<mo_ket_.size();n++){
			if(x.type==PARTICLE) break;
			double xn = x.function.inner(mo_ket_[n]);
			real_function_3d nfz;
			if(z.type == HOLE) nfz = intermediates_.get_fEX(n,z.i);
			else if(z.type == PARTICLE) nfz = intermediates_.get_pfEX(n,z.i);
			else nfz = (*f12op)(mo_bra_[n]*z.function).truncate();
			for(size_t m=0;m<mo_ket_.size();m++){
				double mnfyz = mo_bra_[m].inner(nfz*y.function);
				O1O2_part += xn*mnfyz*mo_ket_[m];
			}
		}

		real_function_3d result = unprojected_part - O1_part - O2_part + O1O2_part;
		return result.truncate();
	}
	real_function_6d test_fill_tree()const{
		return real_factory_6d(world);
	}

	std::pair<std::vector<double>,vecfuncT> decompose_u(const CC_Pair &u)const{
		std::vector<double> cresult;
		vecfuncT gi = get_higher_moments(mo_ket_[u.i]);
		vecfuncT gj = get_higher_moments(mo_ket_[u.j]);
		Q(gi);
		Q(gj);
		for(size_t i=0;i<gi.size();i++){
			real_function_3d xu1 = u.function.project_out(gi[i],0);
			const double xxu1 = xu1.inner(gj[i]);
			real_function_3d diff = xxu1*xu1 - gj[i];
			if(world.rank()==0) std::cout << "Difference between residue and second function=" << diff.norm2() << std::endl;
			real_function_3d xu2 = u.function.project_out(gj[i],1);
			const double xxu2 = xu2.inner(gi[i]);
			if(world.rank()==0) std::cout << "Decomposition of u" << u.i << u.j << " for moment " << i << " gives " << xxu1 << " and " << xxu2 << std::endl;
			cresult.push_back(xxu1);
		}
		std::pair<std::vector<double>,vecfuncT> result(cresult,gi);
		return result;
	}

	vecfuncT get_higher_moments(const real_function_3d &f)const{
		vecfuncT result;
		real_function_3d fx = real_factory_3d(world).f(dipole_x);
		real_function_3d fy = real_factory_3d(world).f(dipole_y);
		real_function_3d fz = real_factory_3d(world).f(dipole_z);
		real_function_3d fr = real_factory_3d(world).f(dipole_r);
		result.push_back(fx*f);
		result.push_back(fy*f);
		result.push_back(fz*f);
		result.push_back(fr*f);

		result.push_back(fx*fx*f);
		result.push_back(fx*fy*f);
		result.push_back(fx*fz*f);
		result.push_back(fx*fr*f);

		//result.push_back(fy*fx*f);
		result.push_back(fy*fy*f);
		result.push_back(fy*fz*f);
		result.push_back(fy*fr*f);

		//result.push_back(fz*fx*f);
		//result.push_back(fz*fy*f);
		result.push_back(fz*fz*f);
		result.push_back(fz*fr*f);

		//result.push_back(fr*fx*f);
		//result.push_back(fr*fy*f);
		//result.push_back(fr*fz*f);
		result.push_back(fr*fr*f);

		for(auto x:result) x.scale(1.0/(x.norm2()));

		return result;

	}

	/// Doubles potentials
	/// G_D4b = G(Q12\sum_k <k|g|j>(1) |i\tau_k> + <k|g|i>(2) |\tau_k,j>)
	/// use that Q12 = Q1*Q2
	/// need to apply G to every single term in order to avoid entangelment

	/// the doubles diagramms of the form:  <i|g|x>*|y\tauj> which are D4b, D6c and D8a
	real_function_6d D4b_D6c_D8a(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		output_section("Now doing G_D4b_D6c_D8a");
		const size_t i=taui.i;
		const size_t j=tauj.i;
		output("6D thresh for all new functions at least = " +stringify(parameters.thresh_Ue));
		// make t intermediate: ti = taui + moi
		real_function_3d ti = taui.function + mo_ket_[i];
		real_function_3d tj = tauj.function + mo_ket_[j];
		real_function_6d result = real_factory_6d(world);
		result.set_thresh(parameters.thresh_Ue);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			real_function_3d kgti= (intermediates_.get_EX(k.i,i))+(intermediates_.get_pEX(k.i,i));
			real_function_3d kgtj= (intermediates_.get_EX(k.i,j))+(intermediates_.get_pEX(k.i,j));

			real_function_3d tmp1 = (kgtj*ti).truncate();
			real_function_3d tmp2 = (kgti*tj).truncate();

			screening(tmp1,k.function);
			screening(k.function,tmp2);

			real_function_6d ik = make_xy(CC_function(tmp1,99,UNDEFINED),k);
			real_function_6d kj = make_xy(k,CC_function(tmp2,99,UNDEFINED));

			result += (ik + kj);

		}
		result.scale(-1.0);
		output("6D thresh for all new functions = " +stringify(parameters.thresh_6D));
		result.set_thresh(parameters.thresh_6D);
		return result;
	}


	/// the doubles diagramms of the form:  integral * |\tauk,\taul> together (D9,D6b,D8b)
	real_function_6d D6b_D8b_D9(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		output_section("Now doing G_D6b_D8b_D9");
		output("6D thresh for all new functions at least = " +stringify(parameters.thresh_Ue));
		const size_t i=taui.i;
		const size_t j=tauj.i;
		CC_function moi(mo_ket_[i],i,HOLE);
		CC_function moj(mo_ket_[j],j,HOLE);
		CC_function ti(mo_ket_[i]+taui.function,i,MIXED);
		CC_function tj(mo_ket_[j]+tauj.function,j,MIXED);
		real_function_6d result = real_factory_6d(world);
		result.set_thresh(parameters.thresh_Ue);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;
				CC_Timer integral_time_1(world,"Integrals decomposed");
				double integral_D6b  = make_integral(k.i,l.i,moi,moj);
				double integral_D8b  = make_integral(k.i,l.i,moi,tauj);
				       integral_D8b += make_integral(k.i,l.i,taui,moj);
				double integral_D9   = make_integral(k.i,l.i,taui,tauj);
				double integral1 = integral_D6b + integral_D8b + integral_D9;
				integral_time_1.info();
				CC_Timer integral_time_2(world,"Integrals with t-intermediates");
				double integral2     = make_integral(k.i,l.i,ti,tj);
				integral_time_2.info();
				if(world.rank()==0){
					std::cout << "Integrals of D6b, D8b and D9 are:\n"
							  << integral_D6b << ", " << integral_D8b <<", " << integral_D9 << "\n"
							  << "Together they give:\n" << integral1 << "\n"
							  << "Integral from t-intermediate is:" << integral2 << std::endl;
				}
				if(fabs(integral1-integral2)>FunctionDefaults<3>::get_thresh())warning("Integrals from t-intermediate has different size than decompose form, diff="+stringify(integral1-integral2));
				// Greens Function on |\tauk,\taul>
				real_function_6d tmp = make_xy(k,l);
//				real_convolution_6d G = Operator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
//				real_function_6d tmp= G(k.function,l.function);
				result += integral1*tmp;
			}
		}
		// no truncate, we will add up small fractions
		output("6D thresh for all new functions = " +stringify(parameters.thresh_6D));
		result.set_thresh(parameters.thresh_6D);
		return result;
	}

	real_function_6d G_D4b_explicit(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			real_function_3d kgj_i = (intermediates_.get_EX(k.i,j)*mo_ket_[i]).truncate();
			real_function_3d kgi_j = (intermediates_.get_EX(j,i)*mo_ket_[j]).truncate();
			Q(kgj_i);
			Q(kgi_j);
			screening(kgj_i,k.function);
			screening(k.function,kgi_j);

			real_function_6d tmp = CompositeFactory<double,6,3>(world).particle1(copy(kgj_i)).particle2(copy(k.function));
			real_function_6d tmpx= CompositeFactory<double,6,3>(world).particle1(copy(k.function)).particle2(copy(kgi_j));
			tmp.fill_tree(G).truncate().reduce_rank();
			tmpx.fill_tree(G).truncate().reduce_rank();
			real_function_6d ket = tmp + tmpx;
			real_function_6d G_ket = G(ket);
			result += G_ket;
		}
		result.truncate();
		result.scale(-1.0);
		return result;
	}
	real_function_6d G_D4b_decomposed(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			real_function_3d kgj_i = (intermediates_.get_EX(k.i,j)*mo_ket_[i]).truncate();
			real_function_3d kgi_j = (intermediates_.get_EX(j,i)*mo_ket_[j]).truncate();
			Q(kgj_i);
			Q(kgi_j);
			screening(kgj_i,k.function);
			screening(k.function,kgi_j);
			real_function_6d appliedG = apply<real_convolution_6d,double,3>(G,kgj_i,k.function);
			real_function_6d appliedGx= apply<real_convolution_6d,double,3>(G,k.function,kgi_j);
			real_function_6d result_k = appliedG + appliedGx;
			result += result_k;
		}
		result.truncate();
		result.scale(-1.0);
		return result;
	}

	/// @\param[out] result=Q12 \sum_{kl} <kl|g|ij> G(|\tau_k,\tau_l>) (Q12 can directly be absorbed into tau states)
	real_function_6d G_D6b_explicit(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		// integrals can also be stored in intermediates (depend only on hole states)
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;
				double integral = 0.0;//intermediates_.get_intergrals_hf()(k.i,l,i,j);
				{
					real_function_3d kgi_j = (intermediates_.get_EX(k.i,i)*mo_ket_[j]).truncate();
					double tmp = mo_bra_[j].inner(kgi_j);
					integral = tmp;
				}
				screening(k.function,l.function);
				real_function_6d tmp = apply<real_convolution_6d,double,3>(G,k.function,l.function);
				tmp.truncate();
				tmp.scale(integral);
				tmp.print_size("<"+stringify(k.i)+stringify(l.i)+"|g|"+stringify(i)+stringify(j)+"> G|\tau"+stringify(k.i)+"tau"+stringify(l.i)+">");
				result += tmp;
			}
		}

		result.truncate();
		return result;
	}
	real_function_6d G_D6b_decomposed(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		// integrals can also be stored in intermediates (depend only on hole states)
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;
				double integral = 0.0;//intermediates_.get_intergrals_hf()(k.i,l,i,j);
				{
					real_function_3d kgi_j = (intermediates_.get_EX(k.i,i)*mo_ket_[j]).truncate();
					double tmp = mo_bra_[j].inner(kgi_j);
					integral = tmp;
				}
				screening(k.function,l.function);
				real_function_6d tmp = apply<real_convolution_6d,double,3>(G,k.function,l.function);
				tmp.truncate();
				tmp.scale(integral);
				tmp.print_size("<"+stringify(k.i)+stringify(l.i)+"|g|"+stringify(i)+stringify(j)+"> G|\tau"+stringify(k.i)+"tau"+stringify(l.i)+">");
				result += tmp;
			}
		}

		result.truncate();
		return result;
	}

	/// @param[out] result = -( GQ12(\tauk,<k|g|i>|\tauj>) - GQ12(<k|g|j>|\taui>,|\tau_k>) + GQ12(<k|g|\tau_j>|i>,|\tauk>) - GQ12(\tauk,<k|g|\tauj>|j>)
	real_function_6d G_D6c(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);

		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			real_function_3d kgi_tauj = (intermediates_.get_EX(k.i,i)*singles(j).function).truncate();
			real_function_3d kgj_taui = (intermediates_.get_EX(k.i,j)*singles(i).function).truncate();
			real_function_3d kgtaui_j = (intermediates_.get_pEX(k.i,i)*mo_ket_[j]).truncate();
			real_function_3d kgtauj_i = (intermediates_.get_pEX(k.i,j)*mo_ket_[i]).truncate();

			Q(kgi_tauj);
			Q(kgj_taui);
			Q(kgtaui_j);
			Q(kgtauj_i);

			screening(k.function,kgi_tauj);
			screening(kgj_taui,k.function);

			real_function_6d part1 = (apply<real_convolution_6d,double,3>(G,k.function,kgi_tauj)); // G(k.function,kgi_tauj)
			real_function_6d part2 = (apply<real_convolution_6d,double,3>(G,kgj_taui,k.function));
			real_function_6d part3 = (apply<real_convolution_6d,double,3>(G,kgtauj_i,k.function));
			real_function_6d part4 = (apply<real_convolution_6d,double,3>(G,kgtaui_j,k.function));
			real_function_6d result_k = part1 - part2 + part3 - part4;
			result += result_k;
			output("D6c details:\n");
			part1.print_size("part1");
			part2.print_size("part2");
			part3.print_size("part3");
			part4.print_size("part4");
		}
		result.truncate();
		result.scale(-1.0);
		return result;
	}

	/// may use particle swap in the future
	/// @paramp[out] result = GQ12( 2.0 <k|g|\tauj>(1) |\taui\tauk> - <k|g|\taui>(1) |\tauj,tauk> + 2.0 <k|g|\taui>(2) |\tauk\tauj> - <k|g|tauj>(2) |\tauk\taui>
	real_function_6d G_D8a(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		output("Now Doing G_D8a, particle swap may be exploited in the future\n\n");
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);

		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			real_function_3d kgtauj_taui = (intermediates_.get_pEX(k.i,j)*singles(i).function).truncate();
			real_function_3d kgtaui_tauj = (intermediates_.get_pEX(k.i,i)*singles(j).function).truncate();
			Q(kgtauj_taui);
			Q(kgtaui_tauj);

			// screenin_
			screening(kgtauj_taui,k.function);
			screening(kgtaui_tauj,k.function);


			// end screening

			real_function_6d part1 = (apply<real_convolution_6d,double,3>(G,kgtauj_taui,k.function)).truncate();
			real_function_6d part2 = (apply<real_convolution_6d,double,3>(G,kgtaui_tauj,k.function)).truncate();
			real_function_6d part3 = (apply<real_convolution_6d,double,3>(G,k.function,kgtaui_tauj)).truncate();
			real_function_6d part4 = (apply<real_convolution_6d,double,3>(G,k.function,kgtauj_taui)).truncate();
			real_function_6d result_k = 2.0*part1 - part2 + 2.0*part3 - part4;
			result += result_k;

			// set back the threshold
			FunctionDefaults<6>::set_thresh(parameters.thresh_6D);

			if(parameters.debug){
				real_function_6d test13 = swap_particles(part1);
				real_function_6d diff13 = test13 - part3;
				real_function_6d test24 = swap_particles(part2);
				real_function_6d diff24 = test24 - part4;
				double norm13 = diff13.norm2();
				double norm24 = diff24.norm2();
				if(world.rank()==0) std::cout << "DEBUG:D8a, testing particle swap exploitation, differences are " << norm13 <<" and " << norm24 << std::endl;

				kgtauj_taui.print_size("<k|tauj>*|taui>");
				k.function.print_size("|tauk>");
				real_function_6d part1_test = CompositeFactory<double,6,3>(world).particle1(copy(kgtauj_taui)).particle2(copy(k.function));
				if(world.rank()==0) std::cout << "G is detructive ? " << G.destructive() << std::endl;
				part1_test.fill_tree(G);
				part1_test.print_size("6D test function after fill_tree and before G application");
				real_function_6d G_part1_test = G(part1_test);
				double diff1 =( part1 - G_part1_test).norm2();
				if(diff1 > FunctionDefaults<6>::get_thresh()) warning("ERROR in part1 of G_D8a, diff="+stringify(diff1));

				kgtauj_taui.print_size("<k|tauj>*|taui>");
				k.function.print_size("|tauk>");
				real_function_6d part4_test = CompositeFactory<double,6,3>(world).particle1(copy(k.function)).particle2(copy(kgtauj_taui));
				part4_test.fill_tree(G);
				part1_test.print_size("6D test function after fill_tree and before G application");
				real_function_6d G_part4_test = G(part4_test);
				double diff4 =( part4 - G_part4_test).norm2();
				if(diff4 > FunctionDefaults<6>::get_thresh()) warning("ERROR in part4 of G_D8a, diff="+stringify(diff4));
			}

		}


		result.truncate();
		return result;
	}

	/// @param[out] result = \sum_{kl} (<kl|g|i,\tauj> + <kl|g|\taui,j>)GQ12|tauk,taul>, Q12 absorbed into tauk and taul
	real_function_6d G_D8b(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		output("Now Doing G_D8b\n\n");
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;
				real_function_3d kgi_tauj = (intermediates_.get_EX(k.i,i)*singles(j).function).truncate();
				real_function_3d lgj_taui = (intermediates_.get_EX(l.i,j)*singles(i).function).truncate();
				double klgitj = mo_bra_[l.i].inner(kgi_tauj);
				double klgtij = mo_bra_[k.i].inner(lgj_taui);

				screening(k.function,l.function);

				real_function_6d tmp = (apply<real_convolution_6d,double,3>(G,k.function,l.function)).truncate();
				result += (klgitj + klgtij)*tmp;
			}
		}

		result.truncate();
		return result;
	}

	///@param[out] result = \sum_{kl} <kl|g|taui,tauj> GQ|tauk,taul>
	real_function_6d G_D9(const CC_function &taui, const CC_function &tauj,const CC_vecfunction &singles)const{
		const size_t i=taui.i;
		const size_t j=tauj.i;
		output("Now Doing G_D9\n\n");
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d result = real_factory_6d(world);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;
				real_function_3d kgtaui_tauj = (intermediates_.get_pEX(k.i,i)*singles(j).function).truncate();
				double integral = mo_bra_[l.i].inner(kgtaui_tauj);

				screening(k.function,l.function);

				real_function_6d tmp = (apply<real_convolution_6d,double,3>(G,k.function,l.function)).truncate();
				result += integral*tmp;
			}
		}

		result.truncate();
		return result;
	}


	real_function_6d make_xy(const CC_function &x, const CC_function &y)const{
		double thresh = guess_thresh(x,y);
		output("Making |"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
		CC_Timer timer(world,"Making |"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
		real_function_6d xy = CompositeFactory<double,6,3>(world).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
		xy.fill_tree().truncate().reduce_rank();
		timer.info();
		return xy;
	}


//	real_function_6d make_screened_hartree_product(const CC_function &x, const CC_function &y)const{
//		double thresh = guess_thresh(x,y);
//		output("Making screened |"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
//		CC_Timer timer(world,"Making screened |"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
//		real_convolution_6d screen_G = BSHOperator<6>(world, sqrt(-2*get_epsilon(x.i,y.i)),parameters.lo, parameters.thresh_bsh_6D);
//		screen_G.modified()=true;
//		screen_G.destructive()=true;
//		real_function_6d xy = CompositeFactory<double,6,3>(world).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
//		xy.fill_tree(screen_G).truncate().reduce_rank();
//		timer.info();
//		return xy;
//	}

	real_function_6d make_f_xy(const CC_function &x, const CC_function &y)const{
		double thresh = guess_thresh(x,y);
		CC_Timer timer(world,"Making f|"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
		output("Making f|"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
		real_function_6d fxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
		fxy.fill_tree().truncate().reduce_rank();
		timer.info();
		return fxy;
	}

//	real_function_6d make_screened_f_xy(const CC_function &x, const CC_function &y)const{
//		double thresh = guess_thresh(x,y);
//		CC_Timer timer(world,"Making screened f|"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
//		output("Making screened f|"+x.name()+","+y.name()+"> with 6D thresh="+stringify(thresh));
//		real_convolution_6d screen_G = BSHOperator<6>(world, sqrt(-2*get_epsilon(x.i,y.i)),parameters.lo, parameters.thresh_bsh_6D);
//		screen_G.modified()=true;
//		screen_G.destructive()=true;
//		real_function_6d fxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x.function)).particle2(copy(y.function)).thresh(thresh);
//		fxy.fill_tree(screen_G).truncate().reduce_rank();
//		timer.info();
//		return fxy;
//	}



	real_function_6d apply_G(const real_function_6d &f, const size_t &i, const size_t &j)const{
		const double eps = get_epsilon(i,j);
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*eps),parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d Gf = G(f);
		Gf.truncate();
		return Gf;
	}

	void testing_mp2_const_part(const CC_Pair &u, const CC_vecfunction &singles)const{

		// get the stored const_mp2_part
		real_function_6d st_const = copy(u.constant_term);
		double st_norm = st_const.norm2();
		st_const.print_size("copy(stored_const_term)");
		u.constant_term.print_size("      stored_const_term");

		// recalculate const_mp2_part
		real_function_6d re_const = real_factory_6d(world);
		{
			real_function_6d kffk = apply_exchange_commutator(CC_function(mo_ket_[0],0,HOLE),CC_function(mo_ket_[0],0,HOLE));
			real_function_6d Ue = apply_transformed_Ue(mo_ket_[0],mo_ket_[0],0,0);
			real_function_6d tmp = Ue - kffk;
			tmp.scale(-2.0);
			apply_Q12(tmp);
			real_function_6d G_tmp = apply_G(tmp,0,0);
			apply_Q12(G_tmp);
			re_const = G_tmp;
		}
		double re_norm = re_const.norm2();
		real_function_6d diff = st_const - re_const;
		double diff_norm = diff.norm2();
		re_const.print_size("recalc_const_term");

		if(world.rank()==0){
			std::cout << "\n\n\n";
			std::cout << "End of Testing Section\n";
			std::cout << "||stored const part||=" << st_norm << "\n";
			std::cout << "||recalc const part||=" << re_norm << "\n";
			std::cout << "|| stored - recalc ||=" << diff_norm << "\n\n\n";
		}
		st_const.print_size("stored_const_term");
		re_const.print_size("recalc_const_term");
		diff.print_size("difference");


		MADNESS_EXCEPTION("TESTING MP2 CONST PART ENDED",1);
	}

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
	const vecfuncT mo_bra_;
	const vecfuncT mo_ket_;
	const vecfuncT active_mo_;
	/// The orbital energies
	const std::vector<double> orbital_energies;
	std::vector<double> init_orbital_energies(const Nemo &nemo)const{
		std::vector<double> eps;
		if(world.rank()==0) std::cout << "SCF Orbital Energies are:\n";
		for(size_t i=0;i<mo_ket_.size();i++){
			eps.push_back(nemo.get_calc()->aeps(i));
			if(world.rank()==0) std::cout << nemo.get_calc()->aeps(i);
		}
		if(world.rank()==0) std::cout <<"\n"<< std::endl;
		return eps;
	}
	/// Helper function to initialize the const mo_bra and ket elements
	vecfuncT make_mo_bra(const Nemo &nemo) const {
		vecfuncT tmp = mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		set_thresh(world,tmp,parameters.thresh_3D);
		return tmp;
	}

	vecfuncT make_mo_ket(const Nemo&nemo)const{
		vecfuncT tmp = nemo.get_calc()->amo;
		set_thresh(world,tmp,parameters.thresh_3D);
		return tmp;
	}
	vecfuncT make_active_mo()const{
		MADNESS_ASSERT(not mo_ket_.empty());
		vecfuncT tmp;
		for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
			tmp.push_back(mo_ket_[i]);
		}
		set_thresh(world,tmp,parameters.thresh_3D);
		return tmp;
	}
	/// The poisson operator (Coulomb Operator)
	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr
			< real_convolution_3d
			> (CoulombOperatorPtr(world,parameters.lo,
					parameters.thresh_poisson_3D));
	/// The BSH Operator for the f12g12 convolution which is with f12= 1/(2gamma)[1-exp(-gamma*r12)], f12g12 = 1/(2gamma) [CoulombOp - BSHOp(gamma)]
	std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr
			< real_convolution_3d
			> (BSHOperatorPtr3D(world, corrfac.gamma(),
					parameters.lo,
					parameters.thresh_poisson_3D));
	/// The f12 convolution operator
	std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr
			< real_convolution_3d
			> (SlaterF12OperatorPtr(world, corrfac.gamma(),
					parameters.lo,
					parameters.thresh_poisson_3D));
	/// Intermediates (some need to be refreshed after every iteration)
	CC_Intermediates intermediates_;
	/// The current singles potential (Q\sum singles_diagrams) , needed for application of the fock opeerator on a singles function
	vecfuncT current_singles_potential;




public:

	//	// Screen Potential
	//	screeningtype screen_potential(const CC_vecfunction &singles,const CC_function &ti, const CC_function &tj, const potentialtype_d &name)const{
	//		size_t i = ti.i;
	//		size_t j = tj.j;
	//		CC_function moi(mo_ket_[i],i,HOLE);
	//		CC_function moj(mo_ket_[j],j,HOLE);
	//		double estimate=0.0;
	//
	//		switch(name);
	//		{
	//		case _D6b_ :
	//			for(auto k:singles.functions){
	//				for(auto l:singles.functions){
	//					estimate += fabs(intermediates_.get_integral(mo_ket_[k.i],mo_ket_[l.i],moi,moj))*k.function.norm2()*l.function.norm2();
	//				}
	//			}
	//		case _D6c_:
	//			for(auto k:singles.functions){
	//				estimate +=  (intermediates_.get_EX(k.i,i).norm2()*k.function.norm2()*tj.function.norm2()
	//		                     +intermediates_.get_EX(k.i,j).norm2()*k.function.norm2()*ti.function.norm2()
	//							 +intermediates_.get_pEX(k.i,i).norm2()*k.function.norm2()*mo_ket_[i].norm2()
	//							 +intermediates_.get_pEX(k.i,j).norm2()*k.function.norm2()*mo_ket_[j].norm2());
	//			}
	//		case _D8a_:
	//
	//		case _D8b_:
	//
	//		case _D9_:
	//
	//		}
	//		estimate = fabs(estimate);
	//		if(estimate>thresh_6D) return _calculate_;
	//		else if(estimate>thresh_6D_tight) return _refine_;
	//		else return _neglect_;
	//	}

	void screening(const real_function_3d &x, const real_function_3d &y)const{
		double normx = x.norm2();
		double normy = y.norm2();
		double norm_xy = normx*normy;
		if(world.rank()==0) std::cout << "Screening |xy> 6D function, norm is: " << norm_xy << std::endl;
		//return norm_xy;
	}

	double guess_thresh(const CC_function &x, const CC_function &y)const{
		double norm = x.function.norm2() * y.function.norm2();
		double thresh = parameters.thresh_6D;
		if(norm > parameters.thresh_6D) thresh= parameters.thresh_6D;
		else if(norm > 0.5*parameters.thresh_6D) thresh= 0.5*parameters.thresh_6D;
		else if(norm > 0.1*parameters.thresh_6D) thresh= 0.1*parameters.thresh_6D;
		else if(norm > 0.05*parameters.thresh_6D) thresh= 0.05*parameters.thresh_6D;
		else if(norm > 0.01*parameters.thresh_6D) thresh= 0.01*parameters.thresh_6D;
		else{
			if(world.rank()==0) std::cout << "Norm of 6D function from 3D function will be " << norm << "... far under the given accuracy ... trying with most precise thresh " << 0.01*parameters.thresh_6D << std::endl;
			return 0.01*parameters.thresh_6D;
		}
		if(world.rank()==0) std::cout << "6D thresh of " << thresh << " is needed to make |xy> with estimated norm of |||xy>||=" << norm << std::endl;
		return thresh;
	}

	// Debug function, content changes from time to time
	void test_potentials()const{

	}
};

} /* namespace madness */

#endif /* CCOPERATORS_H_ */
