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
//#include <string>o

// to debug
//#include<chem/mp2.h>

namespace madness {

template<size_t NDIM>
static double rsquare(const Vector<double,NDIM> &r){
	double result = 0.0;
	for(size_t i=0;i<r.size();i++){
		result += r[i]*r[i];
	}
	return result;
}
template<size_t NDIM>
static double gauss_ND(const Vector<double,NDIM> &r){
	const double r2=rsquare<NDIM>(r);
	const double double_N = (double) NDIM;
	const double c = pow(1.0/(sqrt(2.0*M_PI)),double_N);
	return c*exp(-0.5*r2);
}
template<size_t NDIM>
static double unitfunction(const Vector<double,NDIM>&r){
	return 1.0;
}
static double test_function_3d(const coord_3d &r){
	return 0.1*exp(-rsquare<3>(r));
}
static double unit_function_6d(const coord_6d &r){
	return 1.0;
}

// functors for gauss function
static double f_gauss(const coord_3d &r) {
	return exp(-((r[0]) * (r[0]) + (r[1]) * (r[1]) + (r[2]) * (r[2])));
}
static double f_r2(const coord_3d &r) {
	return (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
static double f_r(const coord_3d &r){
	return sqrt(f_r2(r));
}
static double f_laplace_gauss(const coord_3d&r) {
	return -6.0 * f_gauss(r) + 4.0 * f_r2(r) * f_gauss(r);
}

typedef std::vector<Function<double, 3> > vecfuncT;

/// Structure that holds the CC intermediates and is able to refresh them
struct CC_Intermediates {
public:
	CC_Intermediates(World&world, const CC_vecfunction &bra,
			const CC_vecfunction &ket, const Nemo&nemo,
			const CC_Parameters &param) :
				world(world), parameters(param), mo_bra_(bra), mo_ket_(ket), poisson(
						std::shared_ptr < real_convolution_3d
						> (CoulombOperatorPtr(world, parameters.lo,
								parameters.thresh_poisson_3D))), f12op(
										std::shared_ptr < real_convolution_3d
										> (SlaterF12OperatorPtr(world, parameters.gamma(),
												parameters.lo, parameters.thresh_poisson_3D))), density_(
														make_density(bra, ket)), exchange_intermediate_(
																make_exchange_intermediate(bra, ket)), f12_exchange_intermediate_(
																		make_f12_exchange_intermediate(bra, ket)), hartree_potential_(
																				make_hartree_potential(density_)) {
		sanity_check();

	}

	void sanity_check() const {

	}

	/// Get the intermediates
	const real_function_3d get_density() const {
		return density_;
	}
	const real_function_3d get_perturbed_density() const {
		return perturbed_density_;
	}
	const real_function_3d get_hartree_potential() const {
		return hartree_potential_;
	}
	const real_function_3d get_perturbed_hartree_potential() const {
		return perturbed_hartree_potential_;
	}
	/// returns <k|g|l>
	const real_function_3d get_EX(const size_t &k, const size_t &l) const {
		return exchange_intermediate_(k, l);
	}
	const real_function_3d get_EX(const CC_function &k,
			const CC_function &l) const {
		return exchange_intermediate_(k.i, l.i);
	}
	const real_function_3d get_fEX(const size_t &k, const size_t &l) const {
		return f12_exchange_intermediate_(k, l);
	}
	const real_function_3d get_fEX(const CC_function &k,
			const CC_function &l) const {
		return f12_exchange_intermediate_(k.i, l.i);
	}

	/// returns <k|g|\tau_l>
	const real_function_3d get_pEX(const CC_function &k,
			const CC_function &l) const {
		return perturbed_exchange_intermediate_(k.i, l.i);
	}
	const real_function_3d get_pEX(const size_t &k, const size_t &l) const {
		return perturbed_exchange_intermediate_(k, l);
	}
	/// returns <k|f|\tau_l>
	const real_function_3d get_pfEX(const CC_function &k,
			const CC_function &l) const {
		return perturbed_f12_exchange_intermediate_(k.i, l.i);
	}
	const real_function_3d get_pfEX(const size_t &k, const size_t &l) const {
		return perturbed_f12_exchange_intermediate_(k, l);
	}

	/// refresh the intermediates that depend on the \tau functions
	void update(const CC_vecfunction &tau) {
		if (world.rank() == 0)
			std::cout << "Update Intermediates:\n";
		perturbed_density_ = make_density(mo_bra_, tau);
		perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
		perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,
				tau);
		perturbed_f12_exchange_intermediate_ = make_f12_exchange_intermediate(
				mo_bra_, tau);
		{
			if (world.rank() == 0)
				std::cout << "\n---Updated Intermediates---\n";
			hartree_potential_.print_size("0-Hartree-Potential ");
			perturbed_hartree_potential_.print_size("T1-Hartree-Potential");
			if (world.rank() == 0)
				std::cout << "Exchange Intermediates:\n";
			double size0 = 0.0;
			double sizet = 0.0;
			for (auto x : exchange_intermediate_.allpairs)
				size0 += get_size(world, vecfuncT(1, x.second));
			for (auto x : perturbed_exchange_intermediate_.allpairs)
				sizet += get_size(world, vecfuncT(1, x.second));
		}

	}

	/// make a density from two input functions
	/// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
	/// @param[in] vecfuncT_bra
	/// @param[in] vecfuncT_ket
	/// @param[out] \sum_i bra_i * ket_i
	//real_function_3d make_density(const vecfuncT &bra,const vecfuncT &ket) const;
	real_function_3d make_density(const CC_vecfunction &bra,
			const CC_vecfunction &ket) const;
	/// Poisson operator
	std::shared_ptr<real_convolution_3d> get_poisson() const {
		return poisson;
	}

private:
	World &world;
	const CC_Parameters &parameters;
	const CC_vecfunction &mo_bra_;
	const CC_vecfunction &mo_ket_;
	const std::shared_ptr<real_convolution_3d> poisson;
	const std::shared_ptr<real_convolution_3d> f12op;
	/// const intermediates
	const real_function_3d density_;
	/// Exchange intermediate: \f$EX(i,j) = <i|g|j>\f$
	intermediateT exchange_intermediate_;
	/// The f12 exchange intermediate \f$fEX(i,j) = <i|f12|j>\f$
	intermediateT f12_exchange_intermediate_;
	/// Hartree_Potential  \f$ = J = \sum_k <k|g|k> = \f$ Poisson(density)
	const real_function_3d hartree_potential_;
	/// intermediates that need to be recalculated before every iteration
	/// Perturbed Density \f$= \sum_k |k><\tau_k| \f$
	real_function_3d perturbed_density_;
	/// Perturbed Hartree Potential PJ \f$ = \sum_k <k|g|\tau_k> = \f$ Poisson(perturbed_density)
	real_function_3d perturbed_hartree_potential_;
	/// Perturbed Exchange Intermediate: \f$ PEX(i,j) = <i|g|\tau_j> \f$
	intermediateT perturbed_exchange_intermediate_;
	/// Perturbed f12-exchange-intermediate: \f$ pfEX(i,j) = <i|f12|tau_j> \f$
	intermediateT perturbed_f12_exchange_intermediate_;

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_INTERMEDIATES!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
				1);
	}
	void warning(const std::string &msg) const {
		std::cout << "\n\n\nWARNING IN CC_INTERMEDIATES:\n" << msg
				<< "\n\n\n!!!";
	}
public:
	/// Make the exchange intermediate: EX[j][i] \f$ <bra[i](r2)|1/r12|ket[j](r2)> \f$
	intermediateT make_exchange_intermediate(const CC_vecfunction &bra,
			const CC_vecfunction &ket) const;
	intermediateT make_f12_exchange_intermediate(const CC_vecfunction &bra,
			const CC_vecfunction &ket) const;
	/// Calculates the hartree potential Poisson(density)
	/// @param[in] density A 3d function on which the poisson operator is applied (can be the occupied density and the perturbed density)
	/// @return poisson(density) \f$ = \int 1/r12 density(r2) dr2 \f$
	real_function_3d make_hartree_potential(
			const real_function_3d &density) const {
		real_function_3d hartree = (*poisson)(density);
		hartree.truncate();
		return hartree;
	}
};

/// Coupled Cluster Operators (all closed shell)
class CC_Operators {
public:
	/// Constructor
	CC_Operators(World& world, const Nemo &nemo,
			const CorrelationFactor &correlationfactor,
			const CC_Parameters &param) :
				world(world), nemo(nemo), corrfac(correlationfactor), parameters(
						param), mo_bra_(make_mo_bra(nemo)), mo_ket_(
								make_mo_ket(nemo)), orbital_energies(
										init_orbital_energies(nemo)), intermediates_(world, mo_bra_,
												mo_ket_, nemo, param), Q12(world) {
		// make operators

		// make the active mo vector (ket nemos, bra is not needed for that)
		MADNESS_ASSERT(mo_ket_.size() == mo_bra_.size());
		// initialize the Q12 projector
		Q12.set_spaces(mo_bra_.get_vecfunction(), mo_ket_.get_vecfunction(),
				mo_bra_.get_vecfunction(), mo_ket_.get_vecfunction());
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
	void save_function(const Function<T, NDIM>& f,
			const std::string name) const;

	void plot(const real_function_3d &f, const std::string &msg) const {
		CC_Timer plot_time(world, "plotting " + msg);
		plot_plane(world, f, msg);
		plot_time.info();
	}

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_OPERATORS!!!!\n\n\n\n\n\n\n\n\n\n\n\n", 1);
	}
	void warning(const std::string &msg, CC_data &data) const {
		std::cout << "\n\n\nWARNING IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		warnings.push_back(msg);
		data.warnings.push_back(msg);
	}
	void warning(const std::string &msg) const {
		std::cout << "\n\n\nWARNING IN CC_OPERATORS:\n" << msg << "!!!\n\n\n";
		warnings.push_back(msg);
	}

	void output_section(const std::string&msg) const {
		if (world.rank() == 0) {
			std::cout << "\n\n--------------\n";
			std::cout << msg << std::endl;
			std::cout << "\n";
		}
	}
	void output(const std::string &msg) const {
		if (world.rank() == 0) {
			std::cout << msg << std::endl;
		}
	}

	void update_intermediates(const CC_vecfunction &singles) {
		CC_Timer update(world, "Update Intermediates");
		intermediates_.update(singles);
		update.info();
	}

	CC_function mo_ket(const size_t &i) const {
		return mo_ket_(i);
	}
	CC_vecfunction mo_ket() const {
		return mo_ket_;
	}
	CC_function mo_bra(const size_t &i) const {
		return mo_bra_(i);
	}
	CC_vecfunction mo_bra() const {
		return mo_bra_;
	}

	/// makes the t intermediate which is defined as: \f$ |t_i> = |\tau_i> + |i> \f$
	CC_function make_t_intermediate(const CC_function &tau) const {
		CC_function t(mo_ket_(tau.i).function + tau.function, tau.i, MIXED);
		return t;
	}
	CC_vecfunction make_t_intermediate(const CC_vecfunction &tau) const {
		CC_vecfunction result;
		for (auto x : tau.functions) {
			CC_function tmpi = make_t_intermediate(x.second);
			result.insert(tmpi.i, tmpi);
		}
		return result;
	}

	vecfuncT get_CCS_potential(const CC_vecfunction &singles) const {

		// make a dummy doubles with no content
		Pairs<CC_Pair> doubles;

		vecfuncT result = potential_singles(doubles, singles, pot_F3D_);

		result = add(world, result,
				potential_singles(doubles, singles, pot_S1_)); // brillouin term
		result = add(world, result,
				potential_singles(doubles, singles, pot_S5a_)); // brillouin term
		result = add(world, result,
				potential_singles(doubles, singles, pot_ccs_));
		Q(result);
		truncate(world, result);
		performance_S.current_iteration++;
		return result;
	}

	double make_norm(const CC_function &f)const{return make_norm(f.function);}
	double make_norm(const real_function_3d &f)const{
		const real_function_3d bra = f*nemo.nuclear_correlation -> square();
		const double norm2 = bra.inner(f);
		return sqrt(norm2);
	}

	void test_singles_potential(){

		output_section("Singles Potential Consistency Check with r*|i> singles and Q12f12|ij> doubles");
		// make test singles from mos: |taui> = r*|i>
		// make test doubles from mos: |uij>  = Q12f12|titj>
		real_function_3d r = real_factory_3d(world).f(f_r);
		vecfuncT singles_tmp;
		for(size_t i=parameters.freeze;i<mo_ket_.size();i++){
			real_function_3d tmp = r*mo_ket_(i).function;
			Q(tmp);
			double norm = tmp.norm2();
			tmp.scale(1.0/norm);
			tmp.scale(0.5);
			tmp.print_size("TestSingle: r|" + stringify(i) + ">" );
			singles_tmp.push_back(tmp);
		}
		CC_vecfunction singles(singles_tmp,PARTICLE,parameters.freeze);
		Pairs<CC_Pair> doubles = make_reg_residues(singles);

		update_intermediates(singles);

		CC_data dummy;

		output("\n\n Checking u-parts and r-parts of singles potentials with doubles\n\n");
		const potentialtype_s u_parts_tmp[] = {pot_S4a_u_, pot_S4b_u_, pot_S4c_u_, pot_S2b_u_, pot_S2c_u_};
		const potentialtype_s r_parts_tmp[] = {pot_S4a_r_, pot_S4b_r_, pot_S4c_r_, pot_S2b_r_, pot_S2c_r_};
		std::vector<std::pair<std::string,double> > results;
		for(size_t pot=0;pot<5;pot++){
			const potentialtype_s current_u = u_parts_tmp[pot];
			const potentialtype_s current_r = r_parts_tmp[pot];
			const std::string name = assign_name(current_u);
			double largest_error=0.0;
			output("\n\nConsistency Check of Singles Potential " + assign_name(current_u) + " with " + assign_name(current_r));
			const vecfuncT u = potential_singles(doubles,singles,current_u);
			const vecfuncT r = potential_singles(doubles,singles,current_r);
			const vecfuncT diff = sub(world,u,r);
			const double normdiff = norm2(world,u)-norm2(world,r);
			if(world.rank()==0) std::cout<< std::setw(20) << "||"+assign_name(current_u)+"||-||"+assign_name(current_r)+"||" << std::setfill(' ') << "=" << normdiff << std::endl;
			for(const auto d:diff){
				const double norm = d.norm2();
				if(norm>largest_error) largest_error=norm;
				if(world.rank()==0) std::cout<< std::setw(20) << "||"+assign_name(current_u)+"-"+assign_name(current_r)+"||" << std::setfill(' ') << "=" << norm << std::endl;
			}
			results.push_back(std::make_pair(name,largest_error));
			if(current_u == pot_S2b_u_){
				output("Making Integration Test for S2b potential:");
				// integrate the s2b potential against a function which not in the hole space = \sum_k 2<X,k|g|uik> - <k,X|g|uik>, with X=QX
				real_function_3d X = real_factory_3d(world);
				for(const auto&s: singles.functions){
					X+=s.second.function;
				}
				Q(X);
				X=X*nemo.nuclear_correlation -> square();
				Tensor<double> xs2b = inner(world,X,u);
				Tensor<double> xs2b_reg = inner(world,X,r);
				std::vector<double> control_6D;
				for(auto& itmp:singles.functions){
					const size_t i=itmp.first;
					double controli_6D=0.0;
					for(auto& ktmp:singles.functions){
						const size_t k=ktmp.first;
						real_function_6d g = TwoElectronFactory(world).dcut(parameters.lo);
						real_function_6d Xk_g =CompositeFactory<double, 6, 3>(world).particle1(copy(X)).particle2(copy(mo_bra_(k).function)).g12(g);
						real_function_6d g2 = TwoElectronFactory(world).dcut(parameters.lo);
						real_function_6d kX_g =CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(k).function)).particle2(copy(X)).g12(g2);
						const double tmp_6D = 2.0*Xk_g.inner(get_pair_function(doubles,i,k)) - kX_g.inner(get_pair_function(doubles,i,k));
						controli_6D += tmp_6D;
					}
					control_6D.push_back(controli_6D);
				}
				for(size_t i=0;i<control_6D.size();i++){
					const double diff = xs2b(i) - control_6D[i];
					const double diff2= xs2b_reg(i) - control_6D[i];
					std::cout << "||diffu||_" << i << "=" << fabs(diff) << std::endl;
					std::cout << "||diffr||_" << i << "=" << fabs(diff2) << std::endl;
					if(fabs(diff)>FunctionDefaults<6>::get_thresh()) warning("Integration Test of S2b failed !!!!!");
				}
			}
		}

		bool problems=false;
		for(const auto res:results){
			std::string status = "... passed!";
			if(res.second > FunctionDefaults<6>::get_thresh()){
				status="... failed!";
				problems=true;
			}
			if(world.rank()==0) std::cout << std::setw(10) << res.first << status << " largest error was " << res.second <<std::endl;
		}
		if(problems) warning("Inconsistencies in Singles Potential detected!!!!");
		else output("Singles Potentials seem to be consistent");
		output("\n\n Ending Testing Section\n\n");
	}

	vecfuncT get_CC2_singles_potential(const CC_vecfunction &singles,
			const Pairs<CC_Pair> &doubles) {
		if(parameters.debug){
			std::cout <<"Checking current singles\n";
			for(const auto itmp:singles.functions){
				std::cout << itmp.second.name() << ": type=" << assign_name(itmp.second.type) << ", i=" << itmp.second.i <<", norm=" << itmp.second.function.norm2() <<std::endl;
			}
		}

		Pairs<CC_Pair> doubles_wrapper;
		if (parameters.pair_function_in_singles_potential == FULL) {
			doubles_wrapper = make_full_pairs(doubles, singles);
		} else {
			doubles_wrapper = doubles;
		}

		vecfuncT fock_residue = potential_singles(doubles_wrapper, singles,pot_F3D_);
		vecfuncT result = potential_singles(doubles_wrapper, singles, pot_ccs_);
		result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S2b_u_));
		result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S2c_u_));
		result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4a_u_));
		result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4b_u_));
		result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4c_u_));
		if (parameters.pair_function_in_singles_potential != FULL){
			result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S2b_r_));
			result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S2c_r_));
			result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4a_r_));
			result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4b_r_));
			result = add(world, result,potential_singles(doubles_wrapper, singles, pot_S4c_r_));
		}
		// the fock residue does not get projected, but all the rest
		Q(result);
		truncate(world, result);
		// need to store this for application of Fock oerator on singles ( F|taui> = -Singles_potential[i] + \epsilon_i|taui>)
		current_singles_potential = copy(world,result);
		result = add(world, result, fock_residue);
		performance_S.current_iteration++;
		return result;
	}

	// only get the part of the singles that is produced exclusively by the doulbes in order to make a first guess for the singles
	vecfuncT get_CC2_singles_initial_potential(
			const Pairs<CC_Pair> &doubles) const {
		// make_zero guess
		//		real_function_3d zeroguess = real_factory_3d(world);
		//		vecfuncT tmp(mo_ket_.size(),zeroguess);
		//		CC_vecfunction singles(tmp,PARTICLE,parameters.freeze,tmp.size());
		//		MADNESS_ASSERT(singles.size()==mo_ket_.size()-parameters.freeze);

		vecfuncT result = zero_functions<double, 3>(world,
				mo_ket_.size() - parameters.freeze);
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

	real_function_6d get_CC2_doubles_potential(const CC_Pair &u,const CC_vecfunction &singles) const {
		const real_function_6d coulomb_part = potential_doubles(u, singles, pot_cc2_coulomb_);
		const real_function_6d cc2_residue = potential_doubles(u, singles, pot_cc2_residue_);
		const real_function_6d fock_residue = potential_doubles(u, singles, pot_F6D_);

		real_function_6d potential = coulomb_part + cc2_residue;
		apply_Q12(potential,"coulomb-part+cc2_residue");
		real_function_6d result = fock_residue+potential;
		result.truncate().reduce_rank();
		result.print_size("doubles potential");
		if (world.rank() == 0)performance_D.info(performance_D.current_iteration);
		return result;
	}

	real_function_6d make_cc2_coulomb_parts(const CC_function &taui, const CC_function &tauj, const CC_vecfunction &singles) const {
		const CC_function ti = make_t_intermediate(taui);
		const CC_function tj = make_t_intermediate(tauj);
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * get_epsilon(taui.i,tauj.i)),parameters.lo, parameters.thresh_bsh_6D);
		G.destructive()=true;
		// first do the O1 and O2 parts which are
		// Otau1(g|titj) = |tauk><k|(1)g|titj> = kgti(2)|tauktj>
		// same for Otau2 = kgtj(1)|titauk>
		real_function_6d G_O1tau_part = real_factory_6d(world);
		real_function_6d G_O2tau_part = real_factory_6d(world);
		for(const auto& ktmp:singles.functions){
			const size_t k=ktmp.first;
			const CC_function &tauk=ktmp.second;

			real_function_3d kgti_tj = apply_g12(mo_bra_(k),ti)*tj.function;
			real_function_3d kgtj_ti = apply_g12(mo_bra_(k),tj)*ti.function;
			Q(kgti_tj);
			Q(kgtj_ti);

			real_function_3d tauk_tmp = copy(tauk.function);
			G_O1tau_part +=-2.0*G(tauk_tmp,kgti_tj);
			tauk_tmp = copy(tauk.function);
			G_O2tau_part +=-2.0*G(kgtj_ti,tauk_tmp);
		}

		// works only if the greens operator is not directly applied
		//		// now the O12 part (negative sign)
		//		// Otau12(g|titj>) = Otau1(Otau2g|titj>) = Otau1(O2_part) = |tauk><k|(2)(O2_part)(1,2)
		//		real_function_6d O12tau_part = real_factory_6d(world);
		//		for(const auto& ktmp:singles.functions){
		//			const size_t k=ktmp.first;
		//			const CC_function &tauk=ktmp.second;
		//
		//			const real_function_3d k_O2tau = O2tau_part.project_out(mo_bra_(k).function,1); // 1 is particle 2
		//			O12tau_part = O12tau_part + make_xy(k_O2tau,tauk);
		//		}

		// GOtau12_part
		// make <kl|g|titj>*G(tauk,taul)
		real_function_6d G_O12tau_part = real_factory_6d(world);
		for(const auto& ktmp:singles.functions){
			const size_t k=ktmp.first;
			const CC_function brak = mo_bra_(k);
			for(const auto& ltmp:singles.functions){
				const size_t l=ltmp.first;
				const CC_function bral = mo_bra_(l);
				const real_function_3d taul = copy(ltmp.second.function); // copy because greens is destructive
				const real_function_3d tauk = copy(ktmp.second.function); // copy because greens is destructive
				const double kgftitj=make_integral(k,l,ti,tj);

				G_O12tau_part += -2.0*kgftitj*G(tauk,taul);

			}
		}

		G_O1tau_part.print_size( "G(|tauk><k|g|titj>_2)");
		G_O2tau_part.print_size( "G(|tauk><k|g|titj>_1)");
		G_O12tau_part.print_size("G(|tauk,taul><kl|g|titj>)");
		return G_O1tau_part + G_O2tau_part - G_O12tau_part;
	}

	// computes: G(f(F-eij)|titj> + Ue|titj> - [K,f]|titj>) and uses G-operator screening
	real_function_6d make_cc2_residue_sepparated(const CC_function &taui, const CC_function &tauj){
		const CC_function ti = make_t_intermediate(taui);
		const CC_function tj = make_t_intermediate(tauj);
		const double epsij =  get_epsilon(taui.i,tauj.i);
		const double epsi =  get_orbital_energies()[taui.i];
		const double epsj =  get_orbital_energies()[tauj.i];
		if((epsi+epsj!=epsij)) warning("Error in epsilon values: (epsi+epsj-epsij)="+stringify(epsi+epsj-epsij));
		// Greens operator to apply later:
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*epsij),parameters.lo, parameters.thresh_bsh_6D);
		G.destructive()=true;
		// Greens operator to screen
		real_convolution_6d Gscreen = BSHOperator<6>(world,sqrt(-2.0*epsij),parameters.lo, parameters.thresh_bsh_6D);
		Gscreen.modified() = true;

		const real_function_3d F_K_ti = apply_F(ti)-epsi*ti.function + apply_K(ti);
		const real_function_3d F_K_tj = apply_F(tj)-epsj*tj.function + apply_K(tj);

		output_section("CC2-Residue-Unprojected-Part");
		CC_Timer time_unprojected(world,"CC2-Residue:Unprojected-Partr");
		real_function_6d unprojected_result;
		{
			const real_function_6d fF_fK_part = make_f_xy_screened(F_K_ti,tj,Gscreen) + make_f_xy_screened(ti,F_K_tj,Gscreen);
			const real_function_6d U_part = apply_transformed_Ue(ti,tj);

			const real_function_6d ftitj = make_f_xy_screened(ti,tj,Gscreen);
			const real_function_6d Kf_part = K(ftitj,ti.i==tj.i);

			const real_function_6d V = -2.0*(fF_fK_part + U_part - Kf_part).truncate().reduce_rank();
			V.print_size("-2.0(F-eij+Ue-[K,f])"+ti.name()+tj.name());
			const real_function_6d tmp = G(V);
			unprojected_result = tmp;
			unprojected_result.print_size("G(-2.0(F-eij+Ue-[K,f]))"+ti.name()+tj.name());
		}
		time_unprojected.info();

		output_section("CC2-Residue-Projected-Part");
		CC_Timer time_projected(world,"CC2-Residue:Projected-Part");
		const double tight_thresh = parameters.thresh_6D;
		real_function_6d projected_result=real_factory_6d(world);
		projected_result.set_thresh(tight_thresh);
		output("Tighten thresh to "+stringify(tight_thresh));
		FunctionDefaults<6>::set_thresh(tight_thresh);
		{
			// the f(F-eij+K) operator is of type A12 = f12(A1+A2)
			// (O1+O1-O12)(A12) = k(1)*[(<k|A|x>(2)*y(2) - 1/2 <kl|A|xy> l(2)] + []*l(2)
			// 					= |k> (x) (kAxy_1 - 1/2 im_k) + (kAxy_2 - 1/2 im_k)(x)|k>
			// im_k = \sum_l <kl|A|xy> |l>
			//
			vecfuncT kAxy_1;
			vecfuncT kAxy_2;
			vecfuncT im_k;
			Tensor<double> h(mo_bra_.size(),mo_bra_.size());
			for(const auto & ktmp:mo_bra_.functions){
				const CC_function & k=ktmp.second;
				const real_function_3d kAxy1 = unprojected_result.project_out(k.function,0);
				const real_function_3d kAxy2 = unprojected_result.project_out(k.function,1);
				real_function_3d imk = real_factory_3d(world);
					for(const auto & ltmp:mo_bra_.functions){
					const CC_function & l=ltmp.second;
					h(k.i,l.i) = l.inner(kAxy1);
					imk += h(k.i,l.i)*mo_ket_(l).function;
				}
				kAxy_1.push_back(kAxy1);
				kAxy_2.push_back(kAxy2);
				imk.truncate();
				im_k.push_back(imk);
			}

			for(const auto & ktmp:mo_ket_.functions){
				const CC_function & k=ktmp.second;
				const real_function_3d tmp1  = kAxy_1[k.i] - 0.5*im_k[k.i];
				const real_function_6d part1 = G(-2.0*k.function,tmp1);
				const real_function_3d tmp2  = kAxy_2[k.i] - 0.5*im_k[k.i];
				const real_function_6d part2 = G(tmp2,-2.0*k.function);
				projected_result += (part1+part2).truncate(tight_thresh);
			}
		}
		time_projected.info();
		output("Lowering thresh back to "+stringify(parameters.thresh_6D));
		FunctionDefaults<6>::set_thresh(parameters.thresh_6D);
		return unprojected_result - projected_result;

	}

	// make G(O1+O2-O12)A12|xy>
	//  \left(\O{1}+\O{2}-\O{12}\right)\optwo{A}\ket{xy} = \ket{k}\otimes\left(A^{k}_{x}\ket{y} - \frac{1}{2} A^{kl}_{xy}\ket{l}\right)
    //												     + \left(A^{k}_{y}\ket{x} - \frac{1}{2} A^{kl}_{xy}\ket{l}\right)\otimes\ket{k}
	// operator A should act as: A(k,x,y) = <k|A12|xy>_1 so that <k|A12|xy>_2 = A(k,y,x)
//	real_function_3d apply_G_on_screened_operator(const CC_function &x, const CC_function &y, real_function_3d (*A)(const CC_function&,const CC_function&, const CC_function&))const{
//		const double tight_thresh = parameters.thresh_6D*0.1;
//		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*epsij),parameters.lo, parameters.thresh_bsh_6D);
//		G.destructive()=true;
//		Tensor<double> klAxy = make_matrix_elements(x,y,A);
//		real_function_6d result = real_factory_3d(world);
//		for(const auto & ktmp: mo_ket_.functions){
//			const CC_function & k = ktmp.second;
//			const real_function_3d kAx_y = A(mo_bra_(k),x,y); // = <k|A|xy>_1
//			const real_function_3d kAy_x = A(mo_bra_(k),y,x); // = <k|A|xy>_2
//			real_function_3d klAxy_l = real_factory_3d(world);
//			for(const auto & ltmp: mo_ket_.functions){
//				const CC_function & l = ltmp.second;
//				klAxy_l += klAxy(k.i,l.i)*l;
//			}
//			const real_function_3d part1_tmp = (kax_y - 0.5*klAxy_l).truncate();
//			const real_function_6d part1 = G(-2.0*k.function,part1_tmp);
//			const real_function_3d part2_tmp = (kay_x - 0.5*klAxy_l).truncate();
//			const real_function_6d part2 = G(part1_tmp,-2.0*k.function);
//			result += (part1+part2).truncate(tight_thresh);
//		}
//		return result;
//	}

	// returns \sum_k <k|operator|xy>_1
	real_function_3d screen_operator(const CC_vecfunction &bra,const CC_function &x, const CC_function &y,real_function_3d (*A)(const CC_function&,const CC_function&, const CC_function&))const{
		real_function_3d result = real_factory_3d(world);
		for(const auto & ktmp:bra.functions){
			const CC_function &k=ktmp.second;
			result += A(k,x,y);
		}
		return result;
	}

	// return <k|g|xy>_1 = \sum_k <k|g|x>(2)|y(2)>
	real_function_3d screened_g(const CC_function&k,const CC_function&x, const CC_function&y)const{
		return (apply_g12(k,x)*y.function).truncate();
	}
	// return <k|f|xy>_1 = <k|f|x>(2)|y(2)>
	real_function_3d screened_f(const CC_function&k,const CC_function&x, const CC_function&y)const{
		return (apply_f12(k,x)*y.function).truncate();
	}
	// return <k|F|xy>_1 = <k|F|x> |y(2)> + <k|x> |Fy(2)>
	real_function_3d screened_fF(const CC_function&k,const CC_function&x, const CC_function&y)const{
		// <k|F|x> can be calculated more effectively but since the full Fock operator is applied to y we should keep it consistent
		const real_function_3d Fx = apply_F(x);
		const real_function_3d Fy = apply_F(y);
		const double kFx = k.inner(Fx);
		const double kx =k.inner(x);
		return (kFx*y.function + kx*Fy).truncate();
	}
	// return <k|Kf|xy>_1 = <k|K1f|xy>_1 + <k|K2f|xy>
	// <k|K1f|xy> = <l|kgl*f|x>*|y>
	// <k|K2f|xy> = <l|kfx*g|y>*|l>
	real_function_3d screened_Kf(const CC_function&k,const CC_function&x, const CC_function&y)const{
		const real_function_3d kfx_y = screened_f(k,x,y);
		real_function_3d l_kgl = real_factory_3d(world);
		real_function_3d k_K2f_xy = real_factory_3d(world);
		for(const auto & ltmp:mo_bra_.functions){
			const CC_function &l=ltmp.second;
			l_kgl += l.function*apply_g12(k,mo_ket_(l));
			k_K2f_xy += apply_g12(l,kfx_y)*mo_ket_(l).function;
		}
		l_kgl.truncate();
		k_K2f_xy.truncate();
		const real_function_3d k_K1f_xy = (apply_f12(l_kgl,x)*y.function).truncate();


		return (k_K1f_xy + k_K2f_xy).truncate();

	}

	// return <k|Ff|xy>_1 = <k|F1f|xy>_1 + <k|F2f|xy>
	real_function_3d screened_Ff(const CC_function&k,const CC_function&x, const CC_function&y)const{

	}



	// return <kl|A|xy>
	// the given operator a should act as: A(x,y,z) = <x|A|yz>_1
	Tensor<double> make_matrix_elements(const CC_function &x, const CC_function &y,real_function_3d (*A)(const CC_function&,const CC_function&, const CC_function&))const{
		Tensor<double> result(mo_bra_.size(),mo_bra_.size());
		for(const auto & ktmp:mo_bra_.functions){
			const CC_function & k=ktmp.second;
			const real_function_3d kfxy = A(k,x,y);
			for(const auto & ltmp:mo_bra_.functions){
				const CC_function & l=ltmp.second;
				result(k.i,l.i) = l.inner(kfxy);
			}
		}
		return result;
	}


	real_function_6d get_MP2_potential_constant_part(const CC_Pair &u) const {
		CC_function mo_i = mo_ket_(u.i);
		CC_function mo_j = mo_ket_(u.j);

		CC_Timer timer_U(world, "U|ij>");
		real_function_6d UePart = apply_transformed_Ue(mo_i, mo_j);
		UePart.print_size("Ue|" + mo_i.name() + mo_j.name());
		timer_U.info();

		CC_Timer timer_KffK(world, "Kf|ij>");
		real_function_6d KffKPart = apply_exchange_commutator(mo_i, mo_j);
		KffKPart.print_size("[K,f]|" + mo_i.name() + mo_j.name());
		timer_KffK.info();

		real_function_6d unprojected_result = (UePart - KffKPart).truncate();
		unprojected_result.print_size(
				"Ue - [K,f]|" + mo_i.name() + mo_j.name());
		real_function_6d result = Q12(unprojected_result);
		result.print_size("Q12(Ue - [K,f]" + u.name());
		return result;

	}

	/// returns the non constant part of the MP2 potential which is
	/// \f$ (2J-K+Un)|uij> \f$
	real_function_6d get_MP2_potential_residue(const CC_Pair &u) const {
		CC_Timer timer(world, "(2J-K(R)+Un)|uij>");
		CC_data data("mp2_residue");
		real_function_6d result = fock_residue_6d(u);
		data.result_size = get_size(result);
		data.result_norm = result.norm2();
		data.time = timer.current_time();
		performance_D.insert(data.name, data);
		timer.info();
		return result;
	}

	/// reconstructs the full pair function from the regularized pair functions
	/// used to compute norms of the doubles to compare with LCAO codes
	/// used to debug the singles potential
	/// @param[in] u the regularized function
	/// @return Equation: \f$ \tau = u + Q12f12(|ij> + |taui,j> + |i,tauj> + |taui,tauj>) = u + Q12f12|titj> \f$ with \f$ ti = taui + i \f$
	real_function_6d make_full_pair_function(const CC_Pair &u,
			const CC_function &taui, const CC_function &tauj) const {
		const size_t i = u.i;
		const size_t j = u.j;
		MADNESS_ASSERT(i == taui.i);
		MADNESS_ASSERT(j == tauj.i);
		real_function_3d ti = mo_ket_(i).function + taui.function;
		real_function_3d tj = mo_ket_(j).function + tauj.function;
		real_function_6d Q12f12titj = make_f_xy(ti, tj);
		apply_Q12(Q12f12titj);
		real_function_6d result = u.function + Q12f12titj;
		return result;
	}
	CC_Pair make_full_pair(const CC_Pair &u, const CC_function &taui,
			const CC_function &tauj) const {
		real_function_6d full_pair_function = make_full_pair_function(u, taui,
				tauj);
		CC_Pair result(full_pair_function, u.i, u.j);
		return result;
	}
	Pairs<CC_Pair> make_full_pairs(const Pairs<CC_Pair> &pairs,
			const CC_vecfunction &singles) const {
		Pairs<CC_Pair> result;
		for (auto utmp : pairs.allpairs) {
			CC_Pair u = utmp.second;
			CC_Pair ufull = make_full_pair(u, singles(u.i), singles(u.j));
			result.insert(ufull.i, ufull.j, ufull);
		}
		return result;
	}

	Pairs<CC_Pair> make_reg_residues(const CC_vecfunction &singles) const {
		CC_Timer time(world,"Making Regularization-Tails of Pair-Functions");
		Pairs<CC_Pair> result;
		CC_vecfunction t = make_t_intermediate(singles);
		for (size_t i=parameters.freeze;i<mo_ket_.size();i++) {
			for(size_t j=i;j<mo_ket_.size();j++){
				real_function_6d Qftitj = make_f_xy(t(i),t(j));
				apply_Q12(Qftitj);
				Qftitj.print_size("Q12|t"+stringify(i)+"t"+stringify(j)+">");
				CC_Pair pair_tmp(Qftitj,i,j);
				result.insert(pair_tmp.i, pair_tmp.j, pair_tmp);
			}
		}
		time.info();
		return result;
	}

	// right now this is all copied from mp2.cc
	double compute_mp2_pair_energy(CC_Pair &pair) const;

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
			f -= mo_bra_(i).function.inner(f) * mo_ket_(i).function;
		}
	}

	/// CCSD/CC2 singles potential parts

	/// Genereal function which evaluates a CC_singles potential
	vecfuncT potential_singles(const Pairs<CC_Pair> u,
			const CC_vecfunction & singles, const potentialtype_s &name) const {
		//output_section("Now doing Singles Potential " + assign_name(name));
		if (singles.functions.size() != mo_ket_.size() - parameters.freeze)
			warning(
					"Somethings wrong: Size of singles unequal to size of orbitals minus freeze parameter");
		CC_Timer timer(world, assign_name(name));
		CC_data data(name);
		vecfuncT result;

		switch (name) {
		case pot_F3D_:
			result = fock_residue_closed_shell(singles);
			break;
		case pot_ccs_:
			result = ccs_potential(singles);
			break;
		case pot_S2b_u_:
			result = S2b_u_part(u, singles);
			break;
		case pot_S2c_u_:
			result = S2c_u_part(u, singles);
			break;
		case pot_S4a_u_:
			result = S4a_u_part(u, singles);
			break;
		case pot_S4b_u_:
			result = S4b_u_part(u, singles);
			break;
		case pot_S4c_u_:
			result = S4c_u_part(u, singles);
			break;
		case pot_S2b_r_:
			result = S2b_reg_part(singles);
			break;
		case pot_S2c_r_:
			result = S2c_reg_part(singles);
			break;
		case pot_S4a_r_:
			result = S4a_reg_part(singles);
			break;
		case pot_S4b_r_:
			result = S4b_reg_part(singles);
			break;
		case pot_S4c_r_:
			result = S4c_reg_part(singles);
			break;
		case pot_S1_:
			result = S1(singles);
			break;
		case pot_S5a_:
			result = S5a(singles);
			break;
		}

		truncate(world, result);
		data.result_size = get_size(world, result);
		data.result_norm = norm2(world, result);
		data.time = timer.current_time();
		data.info();
		performance_S.insert(data.name, data);
		if (result.size() != mo_ket_.size() - parameters.freeze)
			warning(
					"Somethings wrong: Size of singles-potential unequal to size of orbitals minus freeze parameter");
		return result;
	}

	real_function_6d potential_doubles(const CC_Pair &u, const CC_vecfunction &singles,
			const potentialtype_d &name) const {
		CC_Timer timer(world, assign_name(name));
		CC_data data(assign_name(name));
		output("Now Doing " + assign_name(name) + " \n\n");

		real_function_6d result = real_factory_6d(world);

		switch (name) {
		//		case pot_D6b_D8b_D9_:
		//			result = D6b_D8b_D9(taui, tauj, singles);
		//			break;
		//		case pot_D4b_D6c_D8a_:
		//			result = D4b_D6c_D8a(taui, tauj, singles);
		//			break;
		case pot_F6D_:
			result = fock_residue_6d(u);
			break;
		case pot_cc2_coulomb_:
			result = make_cc2_coulomb_parts(singles(u.i),singles(u.j),singles);
			break;
		case pot_cc2_residue_:
			result = make_cc2_residue(singles(u.i),singles(u.j));
			break;
		default:
			error(
					"unknown or unsupported key for doubles potential: "
					+ assign_name(name));
			break;
		}

		result.print_size(assign_name(name));
		output("Finished with " + assign_name(name));
		data.result_norm = result.norm2();
		data.result_size = get_size(result);
		data.time = (timer.current_time());
		performance_D.insert(data.name, data);
		return result;
	}

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_i = \sum_k <k|r12|k> |tau_i>
	// K_i = \sum_k <k|r12|tau_i> |k>
	vecfuncT fock_residue_closed_shell(const CC_vecfunction &tau) const;

	/// The CCS Potential without Brillouin terms and Fock residue
	vecfuncT ccs_potential(const CC_vecfunction &tau) const {
		// first form the intermediate t-functions: ti = i + taui
		const CC_vecfunction tfunctions = make_t_intermediate(tau);
		vecfuncT result;

		// get the perturbed hartree_potential: kgtk = sum_k <k|g|\tau_k>
		const real_function_3d kgtauk =
				intermediates_.get_perturbed_hartree_potential();

		for (const auto& ttmp : tfunctions.functions) {
			const CC_function& ti = ttmp.second;
			real_function_3d resulti = real_factory_3d(world);

			const real_function_3d kgtauk_ti = kgtauk * ti.function;
			real_function_3d kgti_tauk = real_factory_3d(world);
			for (const auto &ktmp : tau.functions) {
				const CC_function& tauk = ktmp.second;
				const real_function_3d kgti = (intermediates_.get_pEX(tauk, ti)
						+ intermediates_.get_EX(tauk, ti));
				kgti_tauk += kgti * tauk.function;
			}

			real_function_3d l_kgtauk_ti_taul = real_function_3d(world);
			real_function_3d l_kgti_tauk_taul = real_function_3d(world);
			for (const auto &ltmp : tau.functions) {
				const CC_function& taul = ltmp.second;
				l_kgtauk_ti_taul += mo_bra_(taul).inner(kgtauk_ti)
										* taul.function;
				l_kgti_tauk_taul += mo_bra_(taul).inner(kgti_tauk)
										* taul.function;
			}

			resulti = 2.0 * kgtauk_ti - kgti_tauk - 2.0 * l_kgtauk_ti_taul
					+ l_kgti_tauk_taul;
			result.push_back(resulti);
		}
		return result;
	}

	// result: \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 )
	vecfuncT S2b_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const {
		vecfuncT result;
		if (current_s2b_u_part.empty()) {
			for (const auto& itmp : singles.functions) {
				const size_t i = itmp.first;
				real_function_3d resulti = real_factory_3d(world);
				for (const auto& ktmp : singles.functions) {
					const size_t k = ktmp.first;
					const real_function_6d uik = get_pair_function(doubles, i,
							k);
					// S2b u-part
					{
						const real_function_6d kuik = multiply(copy(uik),
								copy(mo_bra_(k).function), 2);
						poisson->particle() = 2;
						const real_function_6d kguik = (*poisson)(kuik);
						resulti += 2.0 * kguik.dirac_convolution<3>();
					}
					// S2b u-part-exchange
					{
						const real_function_6d kuik = multiply(copy(uik),
								copy(mo_bra_(k).function), 1);
						poisson->particle() = 1;
						const real_function_6d kguik = (*poisson)(kuik);
						resulti -= kguik.dirac_convolution<3>();
					}
				}
				//DEBUG
				{
					const double test = mo_bra_(i).function.inner(resulti);
					std::cout << "<" << mo_bra_(i).name() << "|s2b_" << i << "> = " << test << std::endl;
				}//DEBUG END
				result.push_back(resulti);
			}
			current_s2b_u_part = copy(world,result);
		} else {
			output("found previously calculated S2b-u-part");
			result = copy(world,current_s2b_u_part);
		}
		return result;
	}

	// result: -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1)
	vecfuncT S2c_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const {
		vecfuncT result;
		if (current_s2c_u_part.empty()) {
			for (const auto& itmp : singles.functions) {
				const size_t i = itmp.first;
				real_function_3d resulti = real_factory_3d(world);
				for (const auto& ktmp : singles.functions) {
					const size_t k = ktmp.first;
					const real_function_3d kgi = intermediates_.get_EX(k, i);

					for (const auto& ltmp : singles.functions) {
						const size_t l = ltmp.first;
						const real_function_6d ukl = get_pair_function(doubles,
								k, l);
						const real_function_3d l_kgi = mo_bra_(l).function
								* kgi;
						resulti += -2.0 * ukl.project_out(l_kgi, 1); // 1 means second particle
						resulti += ukl.project_out(l_kgi, 0);
					}
				}
				result.push_back(resulti);
			}
			current_s2c_u_part = copy(world,result);
		} else {
			output("found previously calculated S2c-u-part");
			result = copy(world,current_s2c_u_part);
		}
		return result;
	}

	/// The Part of the CC2 singles potential which depends on singles and doubles (S4a, S4b, S4c)
	vecfuncT S4a_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const {
		// S4a can be computed from the S2b potential
		// (-2<lk|g|uik> + <kl|g|uik>)|tau_l> =( <l( (-2)*<k|g|uik>_2) + <l| (<k|g|uik>_1) )|tau_l> = <l|s2b_u_part>*|tau_l> = - \sum_l <l|s2b_i> |l> important: minus sign and the fact that the s2b potential needs to be unprojected
		vecfuncT s4a;
		for (const auto& itmp : singles.functions) {
			const size_t i = itmp.first;
			real_function_3d s4ai = real_factory_3d(world);
			real_function_3d s4ai_consistency = real_factory_3d(world); // here the unprojected s2b result will be used to check consistency since this is not expensive this will be used everytime the s2b part was stored
			for (const auto& ltmp : singles.functions) {
				const size_t l = ltmp.first;
				const CC_function& taul = ltmp.second;
				for (const auto& ktmp : singles.functions) {
					const size_t k = ktmp.first;
					s4ai += (-2.0
							* make_ijgu(l, k,
									get_pair_function(doubles, i, k))
					+ make_ijgu(k, l,
							get_pair_function(doubles, i, k)))
												* taul.function;
				}
				if(not current_s2b_u_part.empty()){
					s4ai_consistency -= (mo_bra_(l).function.inner(current_s2b_u_part[i-parameters.freeze]))*taul.function;
					std::cout << "||current_s2b_u_part[" << i-parameters.freeze << "]||=" << current_s2b_u_part[i-parameters.freeze].norm2() << std::endl;
					std::cout << "<l|current_s2b_u_part[" << i-parameters.freeze << "]="<< mo_bra_(l).function.inner(current_s2b_u_part[i-parameters.freeze]) << std::endl;
					std::cout << "||taul||=||" << taul.name() << "||=" << taul.function.norm2() << std::endl;
				}
			}
			if(not current_s2b_u_part.empty()){
				const double consistency = (s4ai - s4ai_consistency).norm2();
				if(world.rank()==0){
					std::cout << "||s4a||_" << i << " = " << s4ai.norm2() << std::endl;
					std::cout << "||-sum_l <l|s2b>|taul>||_" << i << " = " << s4ai_consistency.norm2() << std::endl;
					std::cout << "||s4a + sum_l <l|s2b>|taul>||_" << i << " = " << consistency << std::endl;
				}
				if(consistency>FunctionDefaults<6>::get_thresh()) warning("S4a Consistency Check above the 6D thresh");
			}
			s4a.push_back(s4ai);
		}
		return s4a;
	}

	// result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
	vecfuncT S4b_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const {
		vecfuncT result;
		for (const auto& itmp : singles.functions) {
			const size_t i = itmp.first;
			real_function_3d resulti = real_factory_3d(world);
			for (const auto& ktmp : singles.functions) {
				const size_t k = ktmp.first;
				const real_function_3d kgi = intermediates_.get_pEX(k, i);

				for (const auto& ltmp : singles.functions) {
					const size_t l = ltmp.first;
					const real_function_6d ukl = get_pair_function(doubles, k,
							l);
					const real_function_3d l_kgi = mo_bra_(l).function * kgi;
					resulti += -2.0 * ukl.project_out(l_kgi, 1); // 1 means second particle
					resulti += ukl.project_out(l_kgi, 0);
				}
			}
			result.push_back(resulti);
		}
		return result;
	}

	vecfuncT S4c_u_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles) const {
		vecfuncT result;
		// DEBUG
		const CC_vecfunction t = make_t_intermediate(singles);
		// DEBUG END
		for (const auto& itmp : singles.functions) {
			const size_t i = itmp.first;
			real_function_3d resulti = real_factory_3d(world);
			real_function_3d part1 = real_factory_3d(world);
			real_function_3d part2 = real_factory_3d(world);
			real_function_3d part3 = real_factory_3d(world);
			real_function_3d part4 = real_factory_3d(world);
			const real_function_3d kgtauk =intermediates_.get_perturbed_hartree_potential();

			for (const auto& ltmp : singles.functions) {
				const size_t l = ltmp.first;
				const real_function_3d l_kgtauk = mo_bra_(l).function * kgtauk;
				const real_function_6d uil = get_pair_function(doubles, i, l);
				part1 += uil.project_out(l_kgtauk, 1);
				part2 += uil.project_out(l_kgtauk, 0);

				for (const auto& ktmp : singles.functions) {
					const size_t k = ktmp.first;
					const real_function_3d k_lgtauk = mo_bra_(k).function
							* intermediates_.get_pEX(l, k);
					part3 += uil.project_out(k_lgtauk, 1);
					part4 += uil.project_out(k_lgtauk, 0);
				}
			}
			resulti = 4.0*part1-2.0*part2-2.0*part3+part4;
			result.push_back(resulti);
		}
		return result;
	}

	vecfuncT S2b_reg_part(const CC_vecfunction &singles) const {
		vecfuncT result;
		const CC_vecfunction tfunction = make_t_intermediate(singles);
		const real_function_3d ktk = intermediates_.make_density(mo_bra_,tfunction); // the case that tfunction is smaller than mo_bra_ (freeze!=0) is considered
		const real_function_3d kgftk = apply_gf(ktk);
		for (const auto& itmp : tfunction.functions) {			// convenience
			const size_t i = itmp.first;						// convenience
			const CC_function& ti = itmp.second;				// convenience
			real_function_3d resulti = real_factory_3d(world);// this will be the result
			real_function_3d Ipart    = real_factory_3d(world);
			real_function_3d Ipartx   = real_factory_3d(world);
			real_function_3d O1part   = real_factory_3d(world);
			real_function_3d O1partx  = real_factory_3d(world);
			real_function_3d O2part   = real_factory_3d(world);
			real_function_3d O2partx  = real_factory_3d(world);
			real_function_3d O12part  = real_factory_3d(world);
			real_function_3d O12partx = real_factory_3d(world);
			Ipart += 2.0 * kgftk * ti.function; // part1
			for (const auto& ktmp : tfunction.functions) {
				const size_t k = ktmp.first;
				const CC_function& tk = ktmp.second;
				const real_function_3d kti = mo_bra_(k).function * ti.function;
				const real_function_3d kgfti = apply_gf(kti);
				Ipartx += -1.0 * kgfti * tk.function; // part1x

				for (const auto& mtmp : mo_ket_.functions) {
					const size_t m = mtmp.first;
					const CC_function& mom = mtmp.second;
					const real_function_3d mftk = intermediates_.get_fEX(m, k)
											+ intermediates_.get_pfEX(m, k);
					const real_function_3d mfti = intermediates_.get_fEX(m, i)
											+ intermediates_.get_pfEX(m, i);
					const real_function_3d kgm = intermediates_.get_EX(k, m);
					const real_function_3d mfti_tk = mfti * tk.function;
					const real_function_3d mftk_ti = mftk * ti.function;
					O2part -= (2.0 * kgm * mftk_ti); //part3
					O2partx-= (-1.0*kgm * mfti_tk);
					const real_function_3d k_mfti_tk = mo_bra_(k).function
							* mfti_tk;
					const real_function_3d k_gmfti_tk = (*poisson)(k_mfti_tk);
					const real_function_3d k_mftk_ti = mo_bra_(k).function
							* mftk_ti;
					const real_function_3d k_gmftk_ti = (*poisson)(k_mftk_ti);
					O1part -= (2.0 * k_gmfti_tk * mom.function); //part2
					O1partx-=(-1.0*k_gmftk_ti * mom.function);
					for (const auto& ntmp : mo_ket_.functions) {
						const CC_function& mon = ntmp.second;
						const double nmftitk = mo_bra_(mon).inner(mftk_ti);
						const double nmftkti = mo_bra_(mon).inner(mfti_tk);
						O12part  += (2.0 * nmftitk * kgm * mon.function);
						O12partx += (-1.0*nmftkti * kgm * mon.function);
					} // end n
				} // end m
			} // end k
			resulti = Ipart + Ipartx + O1part + O1partx + O2part + O2partx + O12part + O12partx;
			//			Ipart.print_size("unitpart");
			//			O1part.print_size("O1part");
			//			O2part.print_size("O2part");
			//			O12part.print_size("O12part");
			//			Ipartx.print_size("unitpartx");
			//			O1partx.print_size("O1partx");
			//			O2partx.print_size("O2partx");
			//			O12partx.print_size("O12partx");
			result.push_back(resulti);
		} // end i
		return result;
	}

	vecfuncT S2c_reg_part(const CC_vecfunction &singles) const {
		vecfuncT result;
		const CC_vecfunction tfunctions = make_t_intermediate(singles);
		for (const auto& itmp : singles.functions) {
			const CC_function taui = itmp.second;
			real_function_3d resulti = real_factory_3d(world);

			for (const auto& ktmp : tfunctions.functions) {
				const CC_function tk = ktmp.second;
				for (const auto& ltmp : tfunctions.functions) {
					const CC_function tl = ltmp.second;
					const real_function_3d l_kgi_tmp = mo_bra_(tl).function
							* intermediates_.get_EX(tk, taui);
					const CC_function l_kgi(l_kgi_tmp, 99, UNDEFINED);
					resulti -= (2.0 * convolute_x_Qf_yz(l_kgi, tk, tl)
					- convolute_x_Qf_yz(l_kgi, tl, tk));
				}
			}
			result.push_back(resulti);
		}
		return result;
	}

	vecfuncT S4a_reg_part(const CC_vecfunction &singles) const {
		vecfuncT result;
		const CC_vecfunction tfunctions = make_t_intermediate(singles);
		for (const auto& itmp : tfunctions.functions) {
			const CC_function& ti = itmp.second;
			real_function_3d resulti = real_factory_3d(world);

			for (const auto& ktmp : tfunctions.functions) {
				const CC_function& tk = ktmp.second;
				const size_t k = ktmp.first;

				for (const auto& ltmp : singles.functions) {
					const CC_function& taul = ltmp.second;
					const size_t l = ltmp.first;

					const double lkgQftitk = make_ijgQfxy(l, k, ti, tk);
					const double klgQftitk = make_ijgQfxy(k, l, ti, tk);
					resulti -= (2.0 * lkgQftitk - klgQftitk) * taul.function;
				}
			}
			result.push_back(resulti);
		}
		return result;
	}

	/// result: -\sum_{kl}( 2 <l|kgtaui|Qftktl> - <l|kgtaui|Qftltk>
	/// this is the same as S2c with taui instead of i
	vecfuncT S4b_reg_part(const CC_vecfunction &singles) const {
		vecfuncT result;
		const CC_vecfunction tfunctions = make_t_intermediate(singles);
		for (const auto& itmp : singles.functions) {
			const CC_function taui = itmp.second;
			real_function_3d resulti = real_factory_3d(world);

			for (const auto& ktmp : tfunctions.functions) {
				const CC_function tk = ktmp.second;
				for (const auto& ltmp : tfunctions.functions) {
					const CC_function tl = ltmp.second;
					const real_function_3d l_kgi_tmp = mo_bra_(tl).function
							* intermediates_.get_pEX(tk, taui);
					const CC_function l_kgi(l_kgi_tmp, 99, UNDEFINED);
					resulti -= (2.0 * convolute_x_Qf_yz(l_kgi, tk, tl)
					- convolute_x_Qf_yz(l_kgi, tl, tk));
				}
			}
			result.push_back(resulti);
		}
		return result;
	}

	/// result: 4<l|kgtauk|Qftitl> - 2<l|kgtauk|Qftlti> - 2<k|lgtauk|Qftitl> + <k|lgtauk|Qftlti>
	vecfuncT S4c_reg_part(const CC_vecfunction &singles) const {
		vecfuncT result;
		const CC_vecfunction tfunctions = make_t_intermediate(singles);
		for (const auto& itmp : tfunctions.functions) {
			const CC_function& ti = itmp.second;
			real_function_3d resulti = real_factory_3d(world);

			const real_function_3d kgtauk =
					intermediates_.get_perturbed_hartree_potential();

			// first two parts
			real_function_3d part1 = real_factory_3d(world);
			real_function_3d part2 = real_factory_3d(world);
			for (const auto& ltmp : tfunctions.functions) {
				const CC_function& tl = ltmp.second;
				const size_t l = ltmp.first;
				const real_function_3d l_kgtauk =(mo_bra_(l).function * kgtauk);
				part1 += convolute_x_Qf_yz(CC_function(l_kgtauk, 99, UNDEFINED),
						ti, tl);
				part2 += convolute_x_Qf_yz(CC_function(l_kgtauk, 99, UNDEFINED),
						tl, ti);
			}

			// second two parts
			real_function_3d part3 = real_factory_3d(world);
			real_function_3d part4 = real_factory_3d(world);
			for (const auto& ktmp : singles.functions) {
				const CC_function& tauk = ktmp.second;
				const size_t k = ktmp.first;

				for (const auto& ltmp : tfunctions.functions) {
					const CC_function& tl = ltmp.second;
					const size_t l = ltmp.first;

					const real_function_3d k_lgtauk = (mo_bra_(k).function
							* apply_g12(mo_bra_(l), tauk));
					part3 += convolute_x_Qf_yz(
							CC_function(k_lgtauk, 99, UNDEFINED), ti, tl);
					part4 += convolute_x_Qf_yz(
							CC_function(k_lgtauk, 99, UNDEFINED), tl, ti);
				}
			}
			resulti = 4.0 * part1 - 2.0 * part2 - 2.0 * part3 + part4;
			result.push_back(resulti);
		}
		return result;
	}

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
	// -Q \sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\tau i\tau_k> |\tau_l>
	// Q is not applied yet!
	vecfuncT S6(const CC_vecfunction &tau) const;

	// The two brillouin terms S1 and S5a of the singles potential
	vecfuncT S1(const CC_vecfunction &tau) const {
		vecfuncT result;
		for (auto tmpi : tau.functions) {
			CC_function& i = tmpi.second;
			real_function_3d resulti = real_factory_3d(world);
			resulti = apply_F(CC_function(mo_ket_(i.i).function,i.i,UNDEFINED)); // undefined for the testing case where the mos are not converged
			result.push_back(resulti);
		}
		return result;
	}

	vecfuncT S5a(const CC_vecfunction &tau) const {
		vecfuncT result;
		for (auto tmpi : tau.functions) {
			CC_function& i = tmpi.second;
			real_function_3d resulti = real_factory_3d(world);
			for (auto tmpk : tau.functions) {
				CC_function& k = tmpk.second;
				real_function_3d tmp = apply_F(
						CC_function(i.function, i.i, UNDEFINED)); // undefined for the test case where the moi are not converged yet
				const double a = mo_bra_(k.i).function.inner(tmp);
				resulti -= a * k.function;
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
	// @param[in] Structure which holds all current CC pair functions
	// @param[in] Structure which holds all current CC single excitations
	/// \todo Parameter descriptions.
	/// @return Equation: \f$ -Q\sum_k \left( 2<k|g|u_ik> - <k|g|u_ki> + 2<k|gQf|t_it_k> - <k|gQf|t_kt_i> \right), with t_i = i + \tau_i \f$
	/// \note notation: \f$ <k|g|u_ik> = <k(2)|g12|u_ik(1,2)> \f$ (Integration over second particle)
	vecfuncT S2b(const Pairs<CC_Pair> u, const CC_vecfunction &singles,
			CC_data &data) const;
	vecfuncT S2b_3D_part(const CC_vecfunction &singles, CC_data &data) const;
	vecfuncT S2b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction &singles,
			CC_data &data) const;

	/// The s2c and s4b singles potential of CC2 merged together
	vecfuncT S2c4b(const Pairs<CC_Pair> &doubles, const CC_vecfunction &singles,
			CC_data &data) const;
	vecfuncT S2c4b_6D_part(const Pairs<CC_Pair> &doubles,
			const CC_vecfunction &singles, CC_data &data) const;
	vecfuncT S2c4b_3D_part(const CC_vecfunction &singles, CC_data &data) const;

	/// S2c + X Term
	// [Q]   [i]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// \f$ = -Q\sum_{kl}\left( 2<k|lgi|ulk> - <l|kgi|u_{lk}> + 2<k|lgi|t_k>*|t_l> - 2<l|kgi|t_k>*|t_l? \right) \f$
	/// Notation: 6D Integration over second particle, intermediates: lgi = <l|g|i> is the exchange intermediate
	/// Notation: t are the t-intermediates: \f$ |t_i> = |i> + |\tau_i> \f$
	// @param[in] All the current coupled cluster Pairs
	// @param[in] The coupled cluster singles
	/// \todo Parameter descriptions. Unknown stuff commented out by Matt.
	/// @return the S2c+X Potential
	vecfuncT S2c(const Pairs<CC_Pair> &u, const CC_vecfunction &singles,
			CC_data &data) const;
	vecfuncT S2c_3D_part(const CC_vecfunction &singles, CC_data &data) const;
	vecfuncT S2c_6D_part(const Pairs<CC_Pair> &u, const CC_vecfunction &singles,
			CC_data &data) const;
	/// The S4a + X diagram
	//[Q]       [i]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// \f$ -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k> \f$
	vecfuncT S4a(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;
	vecfuncT S4a_3D_part(const CC_vecfunction & tau, CC_data &data) const;
	vecfuncT S4a_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;

	/// The S4b
	//[i]       [Q]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// \f$ -Q\sum_{kl} (2<k(3)l(4)|g34f14|\tau_{i}(3)u_{kl}(1,4)>  // exchange part - <k(4)l(3)|g34f14|\tau_i(3)u_{lk}(1,4)>) \f$
	// 1. make exchange intermedaite X(4) = <k(3)|g34|\tau_i(3)>_3 *  l(4)			Exchange part : Xx(4) = <l(3)|g34|\tau_i(3)>(4) * k(4)
	// 2. make 6d intermediate Y(1,4) = X(4)* u_{kl}(1,4)							Exchange part : Yx(1,4) = X(4)*u_{lk}(1,4)
	// 3. make f14 integration via delta function trick: result(1) = \int f14 Y(1,4) d4 = \int delta(5-1) (\int f54 Y(1,4) d4)d5
	// 3.1 do the convolution Z(1,5) = \int f54 Y(1,4) d4							Exchange part: Zx(1,5) = int f54 Yx(1,4)d4
	// 3.2 project out the unit function: result(1) = <I(5)|Z(1,5)>_5				Exchange part: resultx(1) = <I(5)|Zx(1,5>_5
	vecfuncT S4b(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;
	vecfuncT S4b_3D_part(const CC_vecfunction & tau, CC_data &data) const;
	vecfuncT S4b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;

	/// The S4c + X + X + X + X Diagrams
	//            [i]   [Q]
	//   .......   \    /
	//  /\     /\   \  /
	// _\/_   _\/____\/_
	/// \f$ Q\sum_{kl}[ 4*<k(3)l(4)|g34 f14| \tau_k(3) u_{il}(1,4)> - 2* <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)> \f$
	/// \f$ - 2* <k(3)l(4)|g34 f14| \tau_k(3) U_{li}(1,4)> + <k(3)l(4)|g34 f14|\tau_k(4) u_{li}(1,3)>  ] \f$
	// First and third Terms are solved like this:
	// 1. X(4) = \sum_k (<k(3)|g34|\tau_k(3)>_3(4)) * l(4) = perturbed_hartree_potential(4) * l(4)
	// 2. Y(1,4) = X(4) u_{il}(1,4)			Exchange Part: Yx(4,1) = X(4) u_{li}(4,1)
	// 3.1 Z(1,5) = \int f54 Y(1,4) d4		Exchange Part: Zx(5,1) = \int f54 Yx(4,1) d4
	// 3.2 result(1) = -4 <I(5)|Z(1,5)>_5 -2 <I(5)|Zx(1,5)>_5
	// Second and fourth terms can not use the perturbed hartree potential
	vecfuncT S4c(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;
	vecfuncT S4c_3D_part(const CC_vecfunction & tau, CC_data &data) const;
	vecfuncT S4c_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & tau,
			CC_data &data) const;

	// CC2 Doubles Potential

	/// smalll helper function to track the time for the projetor
	void apply_Q12(real_function_6d &f,
			const std::string &msg = "6D-function") const {
		CC_Timer Q12_time(world, "Applying Q12 to " + msg);
		f = Q12(f);
		Q12_time.info();
	}

	/// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\tau i>+|i>
	// @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
	// @param[in] \tau_j
	/// \todo Parameter descriptions.
	/// @return Equation: \f$ Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj> \f$  with \f$ |ti> = |\tau i>+|i> \f$
	/// Right now Calculated in the decomposed form: \f$ |titj> = |i,j> + |\tau i,\tau j> + |i,\tau j> + |\tau i,j> \f$
	/// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
	real_function_6d make_cc2_residue(const CC_function &taui,
			const CC_function &tauj) const;

	// apply the kinetic energy operator to a decomposed 6D function
	/// @param[in] y a 3d function x (will be particle 1 in the decomposed 6d function)
	/// @param[in] x a 3d function y (will be particle 2 in the decomposed 6d function)
	/// @return a 6d function: G(f12*T*|xy>)
	real_function_6d make_GQfT_xy(const real_function_3d &x,
			const real_function_3d &y, const size_t &i, const size_t &j) const;

	// MP2 Terms are
	// fock_residue_6d = (2J - Kn + Un) |u_{ij}> , KR = R12^{-1}*K*R12 (nuclear tranformed K)
	// Uen|ij> = R12{-1}*U_e*R12 |ij>

	/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
	real_function_6d fock_residue_6d(const CC_Pair &u) const;

	real_function_6d G_fock_residue_xy(const CC_function &taui,
			const CC_function &tauj) const {
		error(
				"G_fock_residue_xy .... this function should not be used, if so check it");
		const size_t i = taui.i;
		const size_t j = tauj.i;
		const real_function_3d & x = taui.function;
		const real_function_3d & y = tauj.function;
		real_convolution_6d G = BSHOperator<6>(world,
				sqrt(-2.0 * get_epsilon(i, j)), parameters.lo,
				parameters.thresh_bsh_6D);
		CC_Timer local_time(world, "Fock-residue-xy-local-part");
		// make x2 = (2.0*J + U2)x, and the same for y
		real_function_3d x2 = ((2.0 * intermediates_.get_hartree_potential()
				+ nemo.nuclear_correlation->U2()) * x).truncate();
		real_function_3d y2 = ((2.0 * intermediates_.get_hartree_potential()
				+ nemo.nuclear_correlation->U2()) * y).truncate();

		// Apply the Greens Function on the local parts
		// G(|x2y> + |xy2>)
		real_function_6d local_part;
		{
			real_function_6d x2y = make_f_xy(CC_function(x2, i, UNDEFINED),
					tauj); //CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x2)).particle2(copy(y));
			real_function_6d xy2 = make_f_xy(taui,
					CC_function(y2, j, UNDEFINED)); //CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(y2));
			apply_Q12(x2y, "x2y");
			apply_Q12(xy2, "xy2");
			real_function_6d Gx2y = G(x2y);
			real_function_6d Gxy2 = G(xy2);
			local_part = Gx2y + Gxy2;
		}
		local_time.info();

		CC_Timer unuc_time(world, "Fock-residue-xy-U1-part");
		real_function_6d U1_part;
		{
			std::vector < std::shared_ptr<real_derivative_3d> > gradop;
			gradop = gradient_operator<double, 3>(world);
			for (size_t axis = 0; axis < 3; ++axis) {
				real_function_3d dx = (*gradop[axis])(x);
				real_function_3d dy = (*gradop[axis])(y);
				real_function_3d U1 = nemo.nuclear_correlation->U1(axis);
				real_function_3d Udx = (U1 * dx).truncate();
				real_function_3d Udy = (U1 * dy).truncate();

				real_function_6d fUdxy =
						CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(
								copy(Udx)).particle2(copy(y));
				real_function_6d fxUdy =
						CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(
								copy(x)).particle2(copy(Udy));
				CC_Timer fill_tree_timer_1(world, "fill_tree_1");
				fUdxy.fill_tree(G).truncate().reduce_rank();
				fill_tree_timer_1.info();
				CC_Timer fill_tree_timer_2(world, "fill_tree_2");
				fxUdy.fill_tree(G).truncate().reduce_rank();
				fill_tree_timer_2.info();
				apply_Q12(fUdxy, "fUdxy");
				apply_Q12(fxUdy, "fxUdy");
				real_function_6d GQfUdxy = G(fUdxy);
				real_function_6d GQfxUdy = G(fxUdy);
				U1_part += (GQfUdxy + GQfxUdy);
			}
		}
		unuc_time.info();

		CC_Timer unuc_faster(world, "Fock-residue-xy-U1-part-faster");
		real_function_3d U1_part_x = real_factory_3d(world);
		real_function_3d U1_part_y = real_factory_3d(world);
		{
			std::vector < std::shared_ptr<real_derivative_3d> > gradop;
			gradop = gradient_operator<double, 3>(world);
			for (size_t axis = 0; axis < 3; ++axis) {
				real_function_3d dx = (*gradop[axis])(x);
				real_function_3d dy = (*gradop[axis])(y);
				real_function_3d U1 = nemo.nuclear_correlation->U1(axis);
				real_function_3d Udx = (U1 * dx).truncate();
				real_function_3d Udy = (U1 * dy).truncate();

				U1_part_x += Udx;
				U1_part_y += Udy;
			}
		}
		real_function_6d fUdxy = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(U1_part_x)).particle2(copy(y));
		real_function_6d fxUdy = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(x)).particle2(copy(U1_part_y));
		CC_Timer fill_tree_timer_1(world, "fill_tree_1");
		fUdxy.fill_tree(G).truncate().reduce_rank();
		fill_tree_timer_1.info();
		CC_Timer fill_tree_timer_2(world, "fill_tree_2");
		fxUdy.fill_tree(G).truncate().reduce_rank();
		fill_tree_timer_2.info();
		apply_Q12(fUdxy, "fUdxy");
		apply_Q12(fxUdy, "fxUdy");
		real_function_6d GQfUdxy = G(fUdxy);
		real_function_6d GQfxUdy = G(fxUdy);
		real_function_6d U1_part_faster = GQfUdxy + GQfxUdy;
		unuc_faster.info();

		// compare two U1 results
		{
			real_function_6d diff = U1_part - U1_part_faster;
			diff.print_size("U1_part - U1_part_faster");
		}

		CC_Timer ex_time(world, "Fock-residue-xy-K-part");
		real_function_6d K_part;
		{
			CC_Timer make_exim_time(world,
					"Make Exchange Intermedaites for fock residue");
			real_function_3d Kx = real_factory_3d(world);
			real_function_3d Ky = real_factory_3d(world);
			vecfuncT mo_bra_x = mul(world, x, mo_bra_.get_vecfunction());
			vecfuncT mo_bra_y = mul(world, y, mo_bra_.get_vecfunction());
			vecfuncT mo_bra_g_x = apply(world, *poisson, mo_bra_x);
			vecfuncT mo_bra_g_y = apply(world, *poisson, mo_bra_y);
			make_exim_time.info();
			for (size_t k = 0; k < mo_bra_.size(); k++) {
				Kx += mo_bra_g_x[k] * mo_ket_(k).function;
				Ky += mo_bra_g_y[k] * mo_ket_(k).function;
			}
			real_function_6d fKxy = CompositeFactory<double, 6, 3>(world).g12(
					corrfac.f()).particle1(copy(Kx)).particle2(copy(y));
			real_function_6d fxKy = CompositeFactory<double, 6, 3>(world).g12(
					corrfac.f()).particle1(copy(x)).particle2(copy(Ky));
			CC_Timer fill_tree_timer_1(world, "fill_tree_1");
			fKxy.fill_tree(G).truncate().reduce_rank();
			fill_tree_timer_1.info();
			CC_Timer fill_tree_timer_2(world, "fill_tree_2");
			fxKy.fill_tree(G).truncate().reduce_rank();
			fill_tree_timer_2.info();
			apply_Q12(fKxy, "fKxy");
			apply_Q12(fxKy, "fxKy");
			real_function_6d GQfKxy = G(fKxy);
			real_function_6d GQfxKy = G(fxKy);
			K_part = (GQfKxy + GQfxKy).truncate();
		}
		ex_time.info();

		real_function_6d result = (local_part + U1_part - K_part).truncate();
		return result;

	}

	/// Exchange Operator on 3D function
	/// !!!!Prefactor (-1) is not included
	real_function_3d K(const CC_function &f) const;
	real_function_3d K(const real_function_3d &f) const;

	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// if i==j in uij then the symmetry will be exploited
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d K(const real_function_6d &u,
			const bool symmetric = false) const;

	/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
	/// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
	/// 1. X(3,2) = bra_k(3)*u(3,2)
	/// 2. Y(1,2) = \int X(3,2) g13 d3
	/// 3. result = Y(1,2)*ket_k(1)
	/// !!!!Prefactor (-1) is not included here!!!!
	real_function_6d apply_K(const real_function_6d &u,
			const size_t &particle) const;

	/// returns \sum_k (<k|g|f> *|k>).truncate()
	real_function_3d apply_K(const CC_function &f) const;

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
	/// @return  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
	real_function_6d apply_transformed_Ue(const CC_function &x,
			const CC_function &y) const;

	/// Apply the Exchange Commutator [K,f]|xy>
	real_function_6d apply_exchange_commutator(const CC_function &x,
			const CC_function &y) const;

	/// Apply the Exchange operator on a tensor product multiplied with f12
	/// !!! Prefactor of (-1) is not inclued in K here !!!!
	real_function_6d apply_Kf(const CC_function &x, const CC_function &y) const;

	/// Apply fK on a tensor product of two 3D functions
	/// fK|xy> = fK_1|xy> + fK_2|xy>
	/// @param[in] x the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
	/// @param[in] y the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
	real_function_6d apply_fK(const CC_function &x, const CC_function &y) const;

	real_function_3d apply_laplacian(const real_function_3d &x) const {
		// make gradient operator for new k and with new thresh
		size_t high_k = 8;
		double high_thresh = 1.e-6;
		std::vector < std::shared_ptr<Derivative<double, 3> > > gradop(3);
		for (std::size_t d = 0; d < 3; ++d) {
			gradop[d].reset(
					new Derivative<double, 3>(world, d,
							FunctionDefaults<3>::get_bc(),
							Function<double, 3>(), Function<double, 3>(),
							high_k));
		}

		// project the function to higher k grid
		real_function_3d f = project(x, high_k);
		f.set_thresh(high_thresh);
		f.refine();

		// apply laplacian
		real_function_3d empty = real_factory_3d(world);
		real_function_3d laplace_f = project(empty, high_k);
		laplace_f.set_thresh(high_thresh);
		for (size_t i = 0; i < gradop.size(); i++) {
			real_function_3d tmp = (*gradop[i])(f);
			real_function_3d tmp2 = (*gradop[i])(tmp);
			laplace_f += tmp2;
		}

		// project laplace_f back to the normal grid
		real_function_3d result = project(laplace_f,
				FunctionDefaults<3>::get_k());
		result.set_thresh(FunctionDefaults<3>::get_thresh());

		// debug and failsafe: make inverse of laplacian and apply
		real_convolution_3d G = BSHOperator<3>(world, 0.0, parameters.lo,
				parameters.thresh_bsh_3D);
		real_function_3d Gresult = -1.0 * G(result);
		real_function_3d difference = x - Gresult;
		double diff = difference.norm2();
		plot_plane(world, difference,
				"Laplacian_error_iteration_"
				+ stringify(performance_D.current_iteration));
		if (world.rank() == 0)
			std::cout << "Apply Laplace:\n" << "||x - G(Laplace(x))||=" << diff
			<< std::endl;
		if (diff > FunctionDefaults<6>::get_thresh())
			warning("Laplacian Error above 6D thresh");

		return result;
	}

	vecfuncT apply_F(const CC_vecfunction &x) const {
		vecfuncT result;
		for(const auto itmp:x.functions){
			const CC_function& xi = itmp.second;
			result.push_back(apply_F(xi));
		}
		return result;
	}

	real_function_3d apply_F(const CC_function &x) const {

		//		if (x.type == HOLE) {
		//			return get_orbital_energies()[x.i] * x.function;
		//		} else if (x.type == PARTICLE) {
		//			const real_function_3d singles_potential = current_singles_potential[x.i-parameters.freeze];
		//			return (get_orbital_energies()[x.i] * x.function - singles_potential);
		//		} else if (x.type == MIXED) {
		//			const real_function_3d singles_potential = current_singles_potential[x.i-parameters.freeze];
		//			return (get_orbital_energies()[x.i] * x.function - singles_potential); // for mixed: eps(i)*x.i = epsi*(moi + taui)
		//		} else if (x.type == UNDEFINED) {
		real_function_3d refined_x = copy(x.function).refine();
		// kinetic part
		CC_Timer T_time(world, "apply_T");
		std::vector < std::shared_ptr<real_derivative_3d> > gradop;
		gradop = gradient_operator<double, 3>(world);
		real_function_3d laplace_x = apply_laplacian(x.function);
		real_function_3d Tx = laplace_x.scale(-0.5).truncate();
		T_time.info();

		CC_Timer J_time(world, "apply_J");
		real_function_3d Jx = (intermediates_.get_hartree_potential()
				* x.function).truncate();
		J_time.info();

		CC_Timer K_time(world, "apply_K");
		real_function_3d Kx = K(x);

		CC_Timer U_time(world, "apply_U");
		real_function_3d U2x =
				(nemo.nuclear_correlation->U2() * x.function).truncate();
		real_function_3d U1x = real_factory_3d(world);
		for (size_t axis = 0; axis < 3; axis++) {
			const real_function_3d U1_axis = nemo.nuclear_correlation->U1(
					axis);
			const real_function_3d dx = (*gradop[axis])(x.function);
			U1x += (U1_axis * dx).truncate();
		}
		U_time.info();

		return (Tx + 2.0 * Jx - Kx + U2x + U1x);
		//}
	//error("apply_F: should not end up here");
		//return real_factory_3d(world);
	}

	/// little helper function to pack a vector of CC_3D_functions (just structures which hold the function the index and the type)
	std::vector<CC_function> make_CC_3D_function(const vecfuncT &f,
			const functype &type) {
		std::vector<CC_function> result(f.size());
		for (size_t i = 0; i < f.size(); i++) {
			CC_function tmp(f[i], i, type);
			result[i] = tmp;
		}
		return result;
	}

	// gives back \epsilon_{ij} = \epsilon_i + \epsilon_j
	double get_epsilon(const size_t &i, const size_t &j) const {
		return (orbital_energies[i] + orbital_energies[j]);
	}
	// gives back the orbital energies
	std::vector<double> get_orbital_energies() const {
		return orbital_energies;
	}
	/// swap particles 1 and 2

	/// param[in] all CC_Pairs
	/// param[in] the i index
	/// param[in] the j index
	/// param[out] a 6d function correspoding to electron pair ij
	/// if i>j the pair will be created via: fij(1,2) = fji(2,1)
	real_function_6d get_pair_function(const Pairs<CC_Pair> &pairs,
			const size_t i, const size_t j) const {
		if (i > j) {
			const real_function_6d & function = pairs(j, i).function;
			const real_function_6d & swapped_function = swap_particles(
					function);
			return swapped_function;

		} else {
			return pairs(i, j).function;
		}
	}

	/// param[in]	f	a function of 2 particles f(1,2)
	/// return	the input function with particles swapped g(1,2) = f(2,1)
	real_function_6d swap_particles(const real_function_6d& f) const;

	// Calculate the CC2 energy equation which is
	// \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
	// with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
	double get_CC2_correlation_energy() const;
	double compute_ccs_correlation_energy(const CC_function &taui,
			const CC_function &tauj) const;
	double compute_cc2_pair_energy(const CC_Pair &u, const CC_function &taui,
			const CC_function &tauj) const;
	/// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
	// the bra elements are always the R2orbitals
	// the ket elements can be \tau_i , or orbitals dependet n the type given
	double make_ij_gQf_ij(const size_t &i, const size_t &j, CC_Pair &u) const;
	double make_ijgQfxy(const size_t &i, const size_t &j, const CC_function &x,
			const CC_function &y) const;
	double make_ijgfxy(const size_t &i, const size_t &j,
			const real_function_3d &x, const real_function_3d &y) const;
	/// Make two electron integral (expensive without intermediates) use just for debugging
	double make_ijgxy(const size_t &i, const size_t &j,
			const real_function_3d &x, const real_function_3d &y) const;
	double make_integral(const size_t &i, const size_t &j, const CC_function &x,
			const CC_function&y) const {
		if (x.type == HOLE) {
			real_function_3d igx_y =
					(intermediates_.get_EX(i, x.i) * y.function).truncate();
			return mo_bra_(j).function.inner(igx_y);
		} else if (x.type == PARTICLE) {
			if (y.type == HOLE) {
				real_function_3d jgy_x = (intermediates_.get_EX(j, y.i)
						* x.function);
				return mo_bra_(i).function.inner(jgy_x);
			} else if (y.type == PARTICLE) {
				real_function_3d jgy_x = (intermediates_.get_pEX(j, y.i)
						* x.function);
				return mo_bra_(i).function.inner(jgy_x);
			}
		} else if (x.type == MIXED or y.type == MIXED) {
			real_function_3d igx =
					((*poisson)(mo_bra_(i).function * x.function)).truncate();
			double result = mo_bra_(j).function.inner(igx * y.function);
			return result;
		} else if (x.type == UNDEFINED or y.type == UNDEFINED) {
			real_function_3d igx =
					((*poisson)(mo_bra_(i).function * x.function)).truncate();
			double result = mo_bra_(j).function.inner(igx * y.function);
			return result;
		} else {
			error("ERROR in make_integrals ... should not end up here");
			return 0.0;
		}
		error("ERROR in make_integrals ... should not end up here");
		return 0.0;
	}
	/// Make two electron integral with the pair function
	double make_ijgu(const size_t &i, const size_t &j, const CC_Pair &u) const;
	double make_ijgu(const size_t &i, const size_t &j,
			const real_function_6d &u) const;
	/// Make two electron integral with BSH operator
	double make_ijGu(const size_t &i, const size_t &j, const CC_Pair &u) const;
	/// apply the operator \f$ gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma) \f$
	/// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
	real_function_3d apply_gf(const real_function_3d &f) const;

	real_function_3d apply_f12(const CC_function & bra,
			const CC_function &ket) const {
		if (bra.type != HOLE) {
			output("Apply_f12, bra state is no hole state");
			return (*f12op)(bra.function * ket.function);
		}
		if (ket.type == HOLE) {
			return intermediates_.get_fEX(bra, ket);
		} else if (ket.type == PARTICLE) {
			return intermediates_.get_pfEX(bra, ket);
		} else if (ket.type == MIXED) {
			return intermediates_.get_fEX(bra, ket)
					+ intermediates_.get_pfEX(bra, ket);
		} else {
			output("Apply_f12, ket state undefined");
			return (*f12op)(bra.function * ket.function);
		}
	}

	real_function_3d apply_g12(const CC_function & bra,
			const CC_function &ket) const {
		if (bra.type != HOLE) {
			output("Apply_g12, bra state is no hole state");
			return (*poisson)(bra.function * ket.function);
		}
		if (ket.type == HOLE) {
			return intermediates_.get_EX(bra, ket);
		} else if (ket.type == PARTICLE) {
			return intermediates_.get_pEX(bra, ket);
		} else if (ket.type == MIXED) {
			return intermediates_.get_EX(bra, ket)
					+ intermediates_.get_pEX(bra, ket);
		} else {
			output("Apply_g12, ket state undefined");
			return (*poisson)(bra.function * ket.function);
		}
	}

	/// @param[in] x Function which is convoluted with (it is assumed that x is already multiplied with R2)
	/// @param[in] y function over which is not integrated
	/// @param[in] z function which is correlated with y over Q12f12
	/// @return <x(2)|Q12f12|y(1)z(2)>_2
	/// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
	/// 1. <x|f12|z>*|y>
	/// 2. - \sum_m mxfyz |m>
	/// 3. - \sum_m <x|m>*mfz*|y>
	/// 4. +\sum_{mn} <x|m> nmfyz |n>
	/// Description: Similar to convolute x_gQf_yz just without coulomb operator
	real_function_3d convolute_x_Qf_yz(const CC_function &x,
			const CC_function &y, const CC_function &z) const {

		const real_function_3d xfz = (*f12op)(x.function * z.function);
		const real_function_3d xfz_y = (xfz * y.function).truncate();
		const real_function_3d part1 = xfz * y.function;

		real_function_3d part2 = real_factory_3d(world);
		real_function_3d part3tmp = real_factory_3d(world);
		real_function_3d part4 = real_factory_3d(world);
		for (const auto& mtmp : mo_ket_.functions) {
			const CC_function& mom = mtmp.second;
			const double mxfyz = mo_bra_(mom).function.inner(xfz_y);
			part2 -= mxfyz * mom.function;

			const double xm = x.function.inner(mom.function);

			const real_function_3d mfz = apply_f12(mo_bra_(mom), z);
			const real_function_3d mfz_y = mfz * y.function;

			part3tmp -= xm * mfz;

			for (const auto& ntmp : mo_ket_.functions) {
				const CC_function& mon = ntmp.second;
				const double nmfyz = mo_bra_(mon).function.inner(mfz_y);
				part4 += xm * nmfyz * mon.function;
			}

		}
		const real_function_3d part3 = part3tmp * y.function;
		real_function_3d result = part1 + part2 + part3 + part4;
		result.truncate();
		return result;

		//
		//		real_function_3d xz = (x.function * z.function).truncate();
		//
		//		// the unprojected part <x(2)|f|z(2)> |y(1)>
		//		real_function_3d unprojected_part =
		//				((*f12op)(xz) * y.function).truncate();
		//
		//		// the O1 part : \sum_m <x(2)|m(1)><m(1)|f12|y(1)z(2)> = <x(2)|mfy(2)|z(2)> |m(1)>
		//		real_function_3d O1_part = real_factory_3d(world);
		//		for (auto mo_bra_m : mo_bra_.functions) {
		//			const CC_function& m = mo_bra_m.second;
		//			real_function_3d mfy;
		//			if (y.type == HOLE)
		//				mfy = intermediates_.get_fEX(m, y);
		//			else if (y.type == PARTICLE)
		//				mfy = intermediates_.get_pfEX(m, y);
		//			else
		//				mfy = ((*f12op)(mo_bra_(m).function * y.function)).truncate();
		//			double a = mfy.inner(xz);
		//			O1_part += a * mo_ket_(m.i).function;
		//		}
		//
		//		// the O2 part: \sum_n <x(2)|n(2)><n(2)|f12|z(2)> |y(1)>
		//		real_function_3d O2_part = real_factory_3d(world);
		//		for (auto mo_bra_n : mo_bra_.functions) {
		//			const CC_function & n = mo_bra_n.second;
		//			if (x.type == PARTICLE)
		//				break;
		//			real_function_3d nfz;
		//			if (z.type == HOLE)
		//				nfz = intermediates_.get_fEX(n, z);
		//			else if (z.type == PARTICLE)
		//				nfz = intermediates_.get_pfEX(n, z);
		//			else
		//				nfz = (*f12op)(mo_bra_(n).function * z.function).truncate();
		//			double xn = x.inner(mo_ket_(n));
		//			O2_part += xn * (nfz * y.function).truncate();
		//		}
		//
		//		// the O1O2 part: \sum_mn <x(2)|n(2)> <m(1)n(2)|f12|y(1)z(2)> |m(1)>
		//		real_function_3d O1O2_part = real_factory_3d(world);
		//		for (auto mo_bra_n : mo_bra_.functions) {
		//			const CC_function & n = mo_bra_n.second;
		//			if (x.type == PARTICLE)
		//				break;
		//			double xn = x.function.inner(mo_ket_(n).function);
		//			real_function_3d nfz;
		//			if (z.type == HOLE)
		//				nfz = intermediates_.get_fEX(n, z);
		//			else if (z.type == PARTICLE)
		//				nfz = intermediates_.get_pfEX(n, z);
		//			else
		//				nfz = (*f12op)(mo_bra_(n).function * z.function);
		//			for (auto mo_bra_m : mo_bra_.functions) {
		//				const CC_function& m = mo_bra_m.second;
		//				double mnfyz = mo_bra_(m).function.inner(nfz * y.function);
		//				O1O2_part += xn * mnfyz * mo_ket_(m).function;
		//			}
		//		}
		//
		//		real_function_3d result = unprojected_part - O1_part - O2_part
		//				+ O1O2_part;
		//		return result.truncate();
	}

	/// Doubles potentials
	/// G_D4b = G(Q12\sum_k <k|g|j>(1) |i\tau_k> + <k|g|i>(2) |\tau_k,j>)
	/// use that Q12 = Q1*Q2
	/// need to apply G to every single term in order to avoid entangelment

	/// the doubles diagramms of the form: \f$ <i|g|x>*|y\tau j> \f$ which are D4b, D6c and D8a
	real_function_6d D4b_D6c_D8a(const CC_function &taui,
			const CC_function &tauj, const CC_vecfunction &singles) const {
		output_section("Now doing combined_D4b_D6c_D8a");
		const size_t i = taui.i;
		const size_t j = tauj.i;
		// make t intermediate: ti = taui + moi
		real_function_3d ti = taui.function + mo_ket_(i).function;
		real_function_3d tj = tauj.function + mo_ket_(j).function;
		real_function_6d result = real_factory_6d(world);
		for (auto tmpk : singles.functions) {
			CC_function& k = tmpk.second;
			real_function_3d kgti = (intermediates_.get_EX(k.i, i))
									+ (intermediates_.get_pEX(k.i, i));
			real_function_3d kgtj = (intermediates_.get_EX(k.i, j))
									+ (intermediates_.get_pEX(k.i, j));

			real_function_3d tmp1 = (kgtj * ti);
			real_function_3d tmp2 = (kgti * tj);

			real_function_6d ik = make_xy(CC_function(tmp1, 99, UNDEFINED), k);
			real_function_6d kj = make_xy(k, CC_function(tmp2, 99, UNDEFINED));

			result += (ik + kj);

		}
		result.scale(-1.0);
		return result;
	}

	/// the doubles diagramms of the form:  integral * |\tau k,\tau l> together (D9,D6b,D8b)
	real_function_6d D6b_D8b_D9(const CC_function &taui,
			const CC_function &tauj, const CC_vecfunction &singles) const {
		output_section("Now doing D6b_D8b_D9");
		const size_t i = taui.i;
		const size_t j = tauj.i;
		const CC_function& moi = mo_ket_(i);
		const CC_function& moj = mo_ket_(j);
		const CC_function ti = make_t_intermediate(taui);
		const CC_function tj = make_t_intermediate(tauj);
		real_function_6d result = real_factory_6d(world);
		for (auto tmpk : singles.functions) {
			CC_function& k = tmpk.second;
			for (auto tmpl : singles.functions) {
				CC_function& l = tmpl.second;
				CC_Timer integral_time_1(world, "Integrals decomposed");
				double integral_D6b = make_integral(k.i, l.i, moi, moj);
				double integral_D8b = make_integral(k.i, l.i, moi, tauj);
				integral_D8b += make_integral(k.i, l.i, taui, moj);
				double integral_D9 = make_integral(k.i, l.i, taui, tauj);
				double integral1 = integral_D6b + integral_D8b + integral_D9;
				integral_time_1.info();
				CC_Timer integral_time_2(world,
						"Integrals with t-intermediates");
				double integral2 = make_integral(k.i, l.i, ti, tj);
				integral_time_2.info();
				if (world.rank() == 0) {
					std::cout << "Integrals of D6b, D8b and D9 are:\n"
							<< integral_D6b << ", " << integral_D8b << ", "
							<< integral_D9 << "\n" << "Together they give:\n"
							<< integral1 << "\n"
							<< "Integral from t-intermediate is:" << integral2
							<< std::endl;
				}
				if (fabs(integral1 - integral2)
						> FunctionDefaults<3>::get_thresh())
					warning(
							"Integrals from t-intermediate has different size than decompose form, diff="
							+ stringify(integral1 - integral2));
				// Greens Function on |\tau k,\tau l>
				real_function_6d tmp = make_xy(k, l);
				//				real_convolution_6d G = Operator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
				//				real_function_6d tmp= G(k.function,l.function);
				result += integral1 * tmp;
			}
		}
		// no truncate, we will add up small fractions
		return result;
	}

	real_function_6d make_xy(const CC_function &x, const CC_function &y) const {
		double thresh = guess_thresh(x, y);
		output(
				"Making |" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		CC_Timer timer(world,
				"Making |" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		real_function_6d xy = CompositeFactory<double, 6, 3>(world).particle1(
				copy(x.function)).particle2(copy(y.function)).thresh(thresh);
		xy.fill_tree().truncate().reduce_rank();
		timer.info();
		return xy;
	}

	real_function_6d make_f_xy(const CC_function &x,
			const CC_function &y) const {
		double thresh = guess_thresh(x, y);
		CC_Timer timer(world,
				"Making f|" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		output(
				"Making f|" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		real_function_6d fxy = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(x.function)).particle2(
						copy(y.function)).thresh(thresh);
		fxy.fill_tree().truncate().reduce_rank();
		timer.info();
		return fxy;
	}

	real_function_6d make_f_xy_screened(const CC_function &x,
			const CC_function &y, const real_convolution_6d &screenG) const {
		double thresh = guess_thresh(x, y);
		CC_Timer timer(world,
				"Making f|" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		output(
				"Making f|" + x.name() + "," + y.name() + "> with 6D thresh="
				+ stringify(thresh));
		real_function_6d fxy = CompositeFactory<double, 6, 3>(world).g12(
				corrfac.f()).particle1(copy(x.function)).particle2(
						copy(y.function)).thresh(thresh);
		fxy.fill_tree(screenG).truncate().reduce_rank();
		timer.info();
		return fxy;
	}

	real_function_6d apply_G(const real_function_6d &f, const size_t &i,
			const size_t &j) const {
		const double eps = get_epsilon(i, j);
		real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * eps),
				parameters.lo, parameters.thresh_bsh_6D);
		real_function_6d Gf = G(f);
		Gf.truncate();
		return Gf;
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
	const CC_vecfunction mo_bra_;
	const CC_vecfunction mo_ket_;
	/// The orbital energies
	const std::vector<double> orbital_energies;
	std::vector<double> init_orbital_energies(const Nemo &nemo) const {
		std::vector<double> eps;
		if (world.rank() == 0)
			std::cout << "SCF Orbital Energies are:\n";
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			eps.push_back(nemo.get_calc()->aeps(i));
			if (world.rank() == 0)
				std::cout << nemo.get_calc()->aeps(i);
		}
		if (world.rank() == 0)
			std::cout << "\n" << std::endl;
		return eps;
	}
	/// Helper function to initialize the const mo_bra and ket elements
	CC_vecfunction make_mo_bra(const Nemo &nemo) const {
		vecfuncT tmp = mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		set_thresh(world, tmp, parameters.thresh_3D);
		CC_vecfunction mo_bra(tmp, HOLE);
		return mo_bra;
	}

	CC_vecfunction make_mo_ket(const Nemo&nemo) const {
		vecfuncT tmp = nemo.get_calc()->amo;
		set_thresh(world, tmp, parameters.thresh_3D);
		CC_vecfunction mo_ket(tmp, HOLE);
		return mo_ket;
	}
	/// The poisson operator (Coulomb Operator)
	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr
			< real_convolution_3d
			> (CoulombOperatorPtr(world, parameters.lo,
					parameters.thresh_poisson_3D));
	/// The BSH Operator for the f12g12 convolution which is with f12= 1/(2gamma)[1-exp(-gamma*r12)], f12g12 = 1/(2gamma) [CoulombOp - BSHOp(gamma)]
	std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr
			< real_convolution_3d
			> (BSHOperatorPtr3D(world, corrfac.gamma(), parameters.lo,
					parameters.thresh_poisson_3D));
	/// The f12 convolution operator
	std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr
			< real_convolution_3d
			> (SlaterF12OperatorPtr(world, corrfac.gamma(), parameters.lo,
					parameters.thresh_poisson_3D));
	/// Intermediates (some need to be refreshed after every iteration)
	CC_Intermediates intermediates_;
	/// The current singles potential (Q\sum singles_diagrams) , needed for application of the fock opeerator on a singles function
	vecfuncT current_singles_potential;
	/// The 6D part of S2b depends only on doubles and HOLE states (if the singles are iterated this does not change, and it is expensive to calculate)
	mutable vecfuncT current_S2b_6D_part_;
	/// The 6D part of S2c depends only on doubles and HOLE states (if the singles are iterated this does not change, and it is expensive to calculate)
	mutable vecfuncT current_s2c_6D_part_;
	/// The 6D part of S2c and S2b (only the parts which depends on the u-function, not the regularization tail (since the singles change)
	mutable vecfuncT current_s2b_u_part;
	mutable vecfuncT current_s2c_u_part;
	StrongOrthogonalityProjector<double, 3> Q12;

public:
	void check_stored_singles_potentials() {
		output("Removing stored singles potentials\n");
		if(not current_s2b_u_part.empty())  warning("S2b_u_part is not empty before singles iteration ... if this is the first macro-iteration this is OK");
		if(not current_s2c_u_part.empty())  warning("S2b_u_part is not empty before singles iteration ... if this is the first macro-iteration this is OK");
		if(not current_singles_potential.empty()) warning("Singles_potential is not empty before singles iteration ... if this is the first macro-iteration this is OK");
	}
	void remove_stored_singles_potentials() {
		output("Removing stored singles potentials\n");
		current_s2b_u_part.clear();
		current_s2c_u_part.clear();
		current_singles_potential.clear();
		// not used anyway
		current_S2b_6D_part_.clear();
		current_s2c_6D_part_.clear();
	}

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

	void screening(const real_function_3d &x, const real_function_3d &y) const {
		double normx = x.norm2();
		double normy = y.norm2();
		double norm_xy = normx * normy;
		if (world.rank() == 0)
			std::cout << "Screening |xy> 6D function, norm is: " << norm_xy
			<< std::endl;
		//return norm_xy;
	}

	double guess_thresh(const CC_function &x, const CC_function &y) const {
		double norm = x.function.norm2() * y.function.norm2();
		double thresh = parameters.thresh_6D;
		if (norm > parameters.thresh_6D)
			thresh = parameters.thresh_6D;
		else if (norm > 0.5 * parameters.thresh_6D)
			thresh = 0.5 * parameters.thresh_6D;
		else if (norm > 0.1 * parameters.thresh_6D)
			thresh = 0.1 * parameters.thresh_6D;
		else if (norm > 0.05 * parameters.thresh_6D)
			thresh = 0.05 * parameters.thresh_6D;
		else if (norm > 0.01 * parameters.thresh_6D)
			thresh = 0.01 * parameters.thresh_6D;
		else {
			if (world.rank() == 0)
				std::cout << "Norm of 6D function from 3D function will be "
				<< norm
				<< "... far under the given accuracy ... trying with most precise thresh "
				<< 0.01 * parameters.thresh_6D << std::endl;
			return 0.01 * parameters.thresh_6D;
		}
		if (world.rank() == 0)
			std::cout << "6D thresh of " << thresh
			<< " is needed to make |xy> with estimated norm of |||xy>||="
			<< norm << std::endl;
		return thresh;
	}

	// Debug function, content changes from time to time
	void test_potentials(const int k, const double thresh =
			FunctionDefaults<3>::get_thresh(), const bool refine = true) const {
		output_section("Little Debug and Testing Session");
		// testing laplace operator with different thresholds
		const size_t old_k = FunctionDefaults<3>::get_k();
		{
			CC_Timer time(world, "time");
			FunctionDefaults<3>::set_thresh(thresh);
			FunctionDefaults<3>::set_k(k);
			std::vector < std::shared_ptr<real_derivative_3d> > gradop;
			gradop = gradient_operator<double, 3>(world);
			std::string name = "_" + stringify(k) + "_" + stringify(thresh);
			if (world.rank() == 0)
				std::cout
				<< "Testing Laplace operator with threshold  "
				+ stringify(FunctionDefaults<3>::get_thresh())
				+ " and k="
				+ stringify(FunctionDefaults<3>::get_k())
				+ " and refinement=" + stringify(refine) + "\n";

			real_function_3d gauss = real_factory_3d(world).f(f_gauss);
			real_function_3d laplace_gauss_analytical =
					real_factory_3d(world).f(f_laplace_gauss);
			real_function_3d laplace_gauss_analytical_old_k = project(
					laplace_gauss_analytical, old_k);
			plot_plane(world, gauss, "gauss" + name);
			plot_plane(world, laplace_gauss_analytical,
					"laplace_gauss_analytical_old_k" + name);

			real_function_3d laplace_gauss_numerical = real_factory_3d(world);
			for (size_t i = 0; i < 3; i++) {
				real_function_3d tmp = (*gradop[i])(gauss);
				real_function_3d tmp2 = (*gradop[i])(tmp);
				laplace_gauss_numerical += tmp2;
			}

			real_function_3d laplace_gauss_diff = laplace_gauss_analytical
					- laplace_gauss_numerical;
			plot_plane(world, laplace_gauss_diff, "laplace_gauss_diff" + name);
			laplace_gauss_diff.print_size(
					"||laplace on gauss num and ana      ||");

			real_function_3d projected_numerical = project(
					laplace_gauss_numerical, old_k);
			real_function_3d laplace_gauss_diff2 =
					laplace_gauss_analytical_old_k - projected_numerical;
			plot_plane(world, laplace_gauss_diff2,
					"laplace_gauss_diff_old_k" + name);
			laplace_gauss_diff.print_size(
					"||laplace on gauss num and ana old k||");

			FunctionDefaults<3>::set_thresh(parameters.thresh_3D);
			FunctionDefaults<3>::set_k(old_k);
			world.gop.fence();
			time.info();
		}
	}

	real_function_3d smooth_function(const real_function_3d &f,
			const size_t mode) const {
		size_t k = f.get_impl()->get_k();
		real_function_3d fproj = project(f, k - 1);
		real_function_3d freproj = project(fproj, k);
		real_function_3d smoothed2 = 0.5 * (f + freproj);
		// double diff = (freproj - f).norm2();
		// double diff2 = (smoothed2 - f).norm2();
		// if(world.rank()==0) std::cout << "||f - f_smoothed|| =" << diff << std::endl;
		// if(world.rank()==0) std::cout << "||f - f_smoothed2||=" << diff2 << std::endl;
		if (mode == 1)
			return freproj;
		else if (mode == 2)
			return smoothed2;
		else {
			std::cout << "Unknown smoothing mode, returning unsmoothed function"
					<< std::endl;
			return f;
		}
	}

	// return <f|F-eps|f> matrix
	Tensor<double> make_reduced_fock_matrix(const CC_vecfunction &f, const double eps)const{
		vecfuncT Ff = apply_F(f);
		vecfuncT fbra = mul(world,nemo.nuclear_correlation->square(),f.get_vecfunction());
		vecfuncT epsf = f.get_vecfunction();
		scale(world,epsf,eps);
		Tensor<double> result(f.size(),f.size());
		for(size_t i=0;i<Ff.size();i++){
			for(size_t j=0;j<fbra.size();j++){
				const double eji = fbra[j].inner(Ff[i]) - fbra[j].inner(epsf[i]);
				result(j,i) = eji;
			}
		}
		return result;
	}


	template<size_t NDIM>
	void test_greens_operators(const double thresh, const size_t k, const double eps, const double bsh_thresh = 1.e-6)const{
		std::cout << "\n\nTesting " << NDIM << "-dimensional Greens Operator with thresh=" << thresh << " BSH-thresh=" << bsh_thresh << " and epsilon=" << eps<< std::endl;
		FunctionDefaults<NDIM>::set_k(k);
		FunctionDefaults<NDIM>::set_thresh(thresh);
		Function<double,NDIM> f = FunctionFactory<double,NDIM>(world).f(gauss_ND<NDIM>);
		f.print_size("TestGaussFunction");
		//Function<double,NDIM> one = FunctionFactory<double,NDIM>(world).f(unitfunction<NDIM>);
		SeparatedConvolution<double,NDIM> G = BSHOperator<NDIM>(world, sqrt(-2.0 * eps),1.e-6, 1.e-5);
		Function<double,NDIM> Lf = general_apply_laplacian<NDIM>(f);
		Lf.print_size("Laplace(f)");
		// Helmholtz eq: (Delta + 2eps)f = ...
		Lf += 2.0*eps*f;
		Lf.truncate();
		Lf.print_size("(Laplace +2eps)f");
		Function<double,NDIM> GLf= G(Lf);

		//const double integral_f = one.inner(f);
		//const double integral_GLf = one.inner(GLf);

		//std::cout << "integral(f)  =" << integral_f   << std::endl;
		//std::cout << "integral(GLf)=" << integral_GLf << std::endl;
		std::cout << "<f|f>     = " << f.inner(f)     << std::endl;
		std::cout << "<GLf|f>   = " << GLf.inner(f)   << std::endl;
		std::cout << "<f|GLf>   = " << f.inner(GLf)   << std::endl;
		std::cout << "<GLf|GLf> = " << GLf.inner(GLf) << std::endl;
		std::cout << "\n\n";


	}



	template<size_t NDIM>
	Function<double,NDIM> general_apply_laplacian(const Function<double,NDIM> &f)const{
		Function<double,NDIM> laplace_f = FunctionFactory<double,NDIM>(world);
		for (int axis = 0; axis < NDIM; ++axis) {
			Derivative<double,NDIM> D = free_space_derivative<double, NDIM>(world,axis);
			const Function<double,NDIM> Df = D(f);
			const Function<double,NDIM> D2f= D(Df).truncate();
			laplace_f += D2f;
		}
		return laplace_f;
	}

};



} /* namespace madness */

#endif /* CCOPERATORS_H_ */
