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



/// Structure that holds the CC intermediates and is able to refresh them
struct CC_Intermediates {
public:
	CC_Intermediates(World&world, const vecfuncT &bra, const vecfuncT &ket,
			const Nemo&nemo, const CC_Parameters &param) :
			world(world), parameters(param), mo_bra_(bra), mo_ket_(ket), poisson(
					std::shared_ptr < real_convolution_3d
							> (CoulombOperatorPtr(world,
									parameters.lo,
									parameters.thresh_poisson_3D))), density_(
					make_density(bra, ket)), exchange_intermediate_(
					make_exchange_intermediate(bra, ket)), hartree_potential_(
							make_hartree_potential(density_)), integrals_hf_(
					make_two_electron_integrals_hf()) {

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
	std::vector<vecfuncT> get_exchange_intermediate() const {
		return exchange_intermediate_;
	}
	std::vector<vecfuncT> get_EX() {
		return exchange_intermediate_;
	}
	std::vector<vecfuncT> get_perturbed_exchange_intermediate() const {
		return perturbed_exchange_intermediate_;
	}
	std::vector<vecfuncT> get_pEX() const {
		return perturbed_exchange_intermediate_;
	}
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
	void update(const vecfuncT &tau) {
		if (world.rank() == 0)
			std::cout << "Update Intermediates:\n";
		perturbed_density_ = make_density(mo_bra_, tau);
		perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
		perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,
				tau);
		integrals_mixed_t1_ = make_two_electron_integrals_mixed_t1(tau);
		integrals_t1_ = make_two_electron_integrals_t1(tau);
	}

	/// make a density from two input functions
	/// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
	/// @param[in] vecfuncT_bra
	/// @param[in] vecfuncT_ket
	/// @param[out] \sum_i bra_i * ket_i
	real_function_3d make_density(const vecfuncT &bra,
			const vecfuncT &ket) const;
	/// Poisson operator
	std::shared_ptr<real_convolution_3d> get_poisson()const{return poisson;}

private:
	World &world;
	const CC_Parameters &parameters;
	const vecfuncT &mo_bra_;
	const vecfuncT &mo_ket_;
	const std::shared_ptr<real_convolution_3d> poisson;
	/// const intermediates
	const real_function_3d density_;
	/// Exchange intermediate: EX(i,j) = <i|g|j>
	const std::vector<vecfuncT> exchange_intermediate_;
	/// Hartree_Potential  = J = \sum_k <k|g|k> = Poisson(density)
	const real_function_3d hartree_potential_;
	/// intermediates that need to be recalculated before every iteration
	/// Perturbed Density = \sum_k |k><\tau_k|
	real_function_3d perturbed_density_;
	/// Perturbed Hartree Poptential PJ = \sum_k <k|g|\tau_k> = Poisson(perturbed_density)
	real_function_3d perturbed_hartree_potential_;
	/// Perturbed Exchange Intermediate: PEX(i,j) = <i|g|\tau_j>
	std::vector<vecfuncT> perturbed_exchange_intermediate_;

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
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT &bra,
			const vecfuncT &ket)const;
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
					make_mo_bra(nemo)), mo_ket_(nemo.get_calc()->amo),orbital_energies(init_orbital_energies(nemo)), intermediates_(
					world, mo_bra_, mo_ket_, nemo, param){

		// initialize the Q12 projector
		Q12.set_spaces(mo_bra_,mo_ket_,mo_bra_,mo_ket_);
	}

	void error(const std::string &msg) const {
		std::cout << "\n\n\nERROR IN CC_OPERATORS:\n" << msg << "\n\n\n!!!";
		MADNESS_EXCEPTION(
				"\n\n!!!!ERROR IN CC_OPERATORS!!!!\n\n\n\n\n\n\n\n\n\n\n\n",
				1);
	}
	void warning(const std::string &msg) const {
		std::cout << "\n\n\nWARNING IN CC_OPERATORS:\n" << msg << "\n\n\n!!!";
	}

	void update_intermediates(const CC_Singles &singles)const{
		CC_Timer update(world,"Update Intermediates");
		vecfuncT tmp=vectorize_singles(singles);
		intermediates_.update(tmp);
		update.info();
	}

	vecfuncT vectorize_singles(const CC_Singles &singles)const{
		vecfuncT tau(singles.size());
		for(size_t i=0;i<singles.size();i++) tau[i] = singles[i].function();
		return tau;
	}

	StrongOrthogonalityProjector<double,3> Q12;

	real_function_3d mo_ket(const size_t &i)const{
		return mo_ket_[i];
	}
	real_function_3d mo_bra(const size_t &i)const{
		return mo_bra_[i];
	}

	vecfuncT get_CIS_potential(const vecfuncT &tau) {
		CC_Singles singles(tau.size());
		for(size_t x=0;x<tau.size();x++) singles[x].function()=tau[x];
		CC_Timer timer_update(world,"update intermediates");
		intermediates_.update(tau);
		timer_update.info();
		vecfuncT result = add(world, S3c(singles), S3c_X(singles));
		Q(result);
		return add(world, result, fock_residue_closed_shell(singles));
	}

	/// makes the t intermediate which is defined as: |t_i> = |\tau_i> + |i>
	vecfuncT make_t_intermediate(const CC_Singles &tau)const{
		vecfuncT result(tau.size());
		vecfuncT mos = mo_ket_;
		for(size_t i=0;i<tau.size();i++){
			result[i] = tau[i].function() + mos[i];
		}
		truncate(world,result);
		return result;
	}

	vecfuncT get_CC2_singles_potential(const CC_Singles &singles, const Pairs<CC_Pair> &doubles)const{
		vecfuncT result = zero_functions<double,3>(world,mo_ket_.size());
		{CC_Timer timer_FR(world,"Singles Potential: Fock Residue");
		result = fock_residue_closed_shell(singles);
		timer_FR.info();}
		{CC_Timer timer_S3c(world,"Singles Potential: S3c");
		result =add(world, S3c(singles),result);
		timer_S3c.info();}
		{CC_Timer timer_S3CX(world,"Singles Potential: S3cX");
		result =add(world, S3c_X(singles),result);
		timer_S3CX.info();}
		{CC_Timer timer_S5b(world,"Singles Potential: S5b");
		result =add(world, S5b(singles),result);
		timer_S5b.info();}
		{CC_Timer timer_S5bx(world,"Singles Potential: S5bX");
		result =add(world, S5b_X(singles),result);
		timer_S5bx.info();}
		{CC_Timer timer_S5c(world,"Singles Potential: S5c");
		result =add(world, S5c(singles),result);
		timer_S5c.info();}
		{CC_Timer timer_S5cx(world,"Singles Potential: S5cX");
		result =add(world, S5c_X(singles),result);
		timer_S5cx.info();}
		{CC_Timer timer_S6(world,"Singles Potential: S6");
		result =add(world, S6(singles),result);
		timer_S6.info();}
		{CC_Timer timer_S2b(world,"Singles Potential: S2b+X");
		result =add(world, S2b(doubles,singles),result);
		timer_S2b.info();}
		{CC_Timer timer_S2c(world,"Singles Potential: S2c+X");
		result =add(world, S2c(doubles),result);
		timer_S2c.info();}
		{CC_Timer timer_S4a(world,"Singles Potential: S4a+X");
		result =add(world, S4a(doubles,singles),result);
		timer_S4a.info();}
		{CC_Timer timer_S4b(world,"Singles Potential: S4b+X");
		result =add(world, S4b(doubles,singles),result);
		timer_S4b.info();}
		{CC_Timer timer_S4c(world,"Singles Potential: S4c+X");
		result =add(world, S4c(doubles,singles),result);
		timer_S4c.info();}
		Q(result);
		truncate(world,result);
		return result;
	}

	// only get the part of the singles that is produced exclusively by the doulbes in order to make a first guess for the singles
	vecfuncT get_CC2_singles_initial_potential(const Pairs<CC_Pair> &doubles)const{
		CC_Singles singles;
		// make_zero guess
		real_function_3d zeroguess = real_factory_3d(world);
		for(size_t i=0;i<mo_ket_.size();i++){
			CC_Single tmp(i,zeroguess);
			singles.push_back(tmp);
		}
		vecfuncT result = zero_functions<double,3>(world,mo_ket_.size());
		{CC_Timer timer_S2b(world,"Singles Potential: S2b+X");
		result =S2b(doubles,singles);
		for(size_t i=0;i<result.size();i++) result[i].print_size("S2b_"+stringify(i));
		timer_S2b.info();}
		{CC_Timer timer_S2c(world,"Singles Potential: S2c+X");
		vecfuncT s2c = S2c(doubles);
		for(size_t i=0;i<result.size();i++) s2c[i].print_size("S2c_"+stringify(i));
		result = add(world,s2c,result);
		timer_S2c.info();}
		return result;
	}

	real_function_6d get_CC2_doubles_potential(const vecfuncT &singles, const CC_Pair &u)const{
	MADNESS_EXCEPTION("CC2 doubles potential not yet implemented",1);
	return real_factory_6d(world);
	}

	real_function_6d get_MP2_potential_constant_part(CC_Pair &u)const{
		CC_Timer timer_U(world,"Ue(R)|ij>");
		real_function_6d UePart = apply_transformed_Ue(mo_ket_[u.i],mo_ket_[u.j],u.i,u.j,u);
		UePart.print_size("Ue|"+stringify(u.i)+stringify(u.j)+">");
		timer_U.info();

		CC_Timer timer_KffK(world,"Kf|ij>");
		real_function_6d KffKPart = apply_exchange_commutator(mo_ket_[u.i],mo_ket_[u.j],"occupied",u.i,u.j);
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
		CC_Timer timer_mp2res(world,"(2J-K(R)+Un)|uij>");
		real_function_6d result = fock_residue_6d(u);
		timer_mp2res.info();
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
	void Q(real_function_3d &f) const {
		for (size_t i = 0; i < mo_ket_.size(); i++) {
			f -= mo_bra_[i].inner(f) * mo_ket_[i];
		}
	}

	/// CCSD/CC2 singles potential parts

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_i = \sum_k <k|r12|k> |tau_i>
	// K_i = \sum_k <k|r12|tau_i> |k>
	vecfuncT fock_residue_closed_shell(const CC_Singles &tau) const;

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3c(const CC_Singles &tau) const;

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = -Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3c_X(const CC_Singles &tau) const;

	/// The S5b term
	//[i]    [Q]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// 2\sum_k <k|g|\tau_k> |\tau_i>
	// No Q is applied yet !
	vecfuncT S5b(const CC_Singles &tau) const;

	/// The S5b Exchange Term
	//[i]         [Q]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_k <k|g|\tau_i> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5b_X(const CC_Singles &tau) const;

	/// The S5c term
	//[Q]    [i]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// -2\sum_kl <kl|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	// May use alteriative algorithm with perturbed density intermediate
	vecfuncT S5c(const CC_Singles &tau) const;

	/// The S5c_X echange term
	//[Q]         [i]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_kl <lk|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5c_X(const CC_Singles &tau) const;

	/// The S6+X Term
	// \    /\    /...
	//  \  /  \  /   /\
	//  _\/_  _\/_  _\/_
	// -Q \sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>
	// Q is not applied yet!
	vecfuncT S6(const CC_Singles  &tau) const;

	/// CC2 singles diagrams with 6d functions as input
	/// Use GFInterface in function_interface.h as kernel (f*g) and do not reconstruct \tau = f12u(1,2) if possible
	/// Since the correlation factor of CC2 has Slater form like in MP2: g12f12 = g12(1-exp(-mu*r12)/r12) = g12 - exp(-mu*r12)/r12 = Coulomb_Operator - BSH_Operator(mu)

	/// S2b + X Term
	// [i]   [Q]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// Current procedure:
	/// use g12 = \int \delta(1-3) g32 d3
	/// <k(2)|g12|u(1,2)> = \int d2[ g12x(1,2 ] with x(1,2) = k(2)u(1,2)
	/// = int d2 [ int d3[ \delta(1-3) g32 ] x(1,2) ]
	/// = \int d3[\delta(1-3) \int d2 [ g32 x(1,2 ] ]
	/// = \int d3[\delta(1-3) h(1,3)] with h(1,3) = \int d2 g23 x(1,2)
	vecfuncT S2b(const Pairs<CC_Pair> u, const CC_Singles &singles) const;

	/// S2c + X Term
	// [Q]   [i]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// = \sum <k(3)l(4)|g34 f31| u_{lk}(1,3) i(4)> = \sum <k(3)|f13| X_{lk,li}(1,3) > with X_{lk,li}(1,3) = u_lk(1,3) * (<l(4)|g34>|i(4)>_4)(3)
	/// = \sum \int d5 \delta(5-1) \int d3 f53 k(3)*X_{lk,li}(5,3)
	vecfuncT S2c(const Pairs<CC_Pair> u) const;

	/// The S4a + X diagram
	//[Q]       [i]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k>
	vecfuncT S4a(const Pairs<CC_Pair> u, const CC_Singles & tau) const;

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
	vecfuncT S4b(const Pairs<CC_Pair> u, const CC_Singles & tau) const;

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
	vecfuncT S4c(const Pairs<CC_Pair> u, const CC_Singles & tau) const;

	// CC2 Doubles Potential

	// MP2 Terms are
	// fock_residue_6d = (2J - Kn + Un) |u_{ij}> , KR = R12^{-1}*K*R12 (nuclear tranformed K)
	// Uen|ij> = R12{-1}*U_e*R12 |ij>

	/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
	real_function_6d fock_residue_6d(const CC_Pair &u) const;

	/// Echange Operator on 3D function
	/// !!!!Prefactor (-1) is not included
	real_function_3d K(const real_function_3d &f,const size_t &i, const bool hc=false)const;

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
			const real_function_3d y, const size_t &i, const size_t &j, CC_Pair &u) const;

	/// Apply the Exchange Commutator [K,f]|xy>
	real_function_6d apply_exchange_commutator(const real_function_3d &x, const real_function_3d &y,const std::string &type, const size_t &i, const size_t &j)const;

	/// Apply the Exchange operator on a tensor product multiplied with f12
	/// !!! Prefactor of (-1) is not inclued in K here !!!!
	real_function_6d apply_Kf(const real_function_3d &x,
			const real_function_3d &y, const std::string &type, const size_t &i, const size_t &j) const;

	/// Apply fK on a tensor product of two 3D functions
	/// fK|xy> = fK_1|xy> + fK_2|xy>
	/// @param[in] x, the first 3D function in |xy>
	/// @param[in] y, the second 3D function in |xy>
	/// @param[in] type, specifies if |xy> = |ij> (occupied), |xy> = |\tau_i,j> (mixed) or |xy> = |\tau_i\tau_j> (virtual)
	/// @param[in] i, the number of the function: bsp if occupied then x_i = |i>, if virtual then x_i = \tau_i etc
	/// @param[in] j , index of the second function
	real_function_6d apply_fK(const real_function_3d &x,
			const real_function_3d &y, const std::string &type, const size_t &i,
			const size_t &j) const;

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
	double compute_cc2_pair_energy(const CC_Pair &u,
			const real_function_3d &taui, const real_function_3d &tauj) const;
	/// Calculate the integral <bra1,bra2|gQf|ket1,ket2>
	// the bra elements are always the R2orbitals
	// the ket elements can be \tau_i , or orbitals dependet n the type given
	double make_ij_gQf_ij(const size_t &i, const size_t &j,CC_Pair &u)const;
	double make_ijgQfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	double make_ijgfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	/// Make two electron integral (expensive without intermediates) use just for debugging
	double make_ijgxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const;
	/// Make two electron integral with the pair function
	double make_ijgu(const size_t &i, const size_t &j, const CC_Pair &u)const;
	/// Make two electron integral with BSH operator
	double make_ijGu(const size_t &i, const size_t &j, const CC_Pair &u)const;
	/// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
	/// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
	real_function_3d apply_gf(const real_function_3d &f)const;
	real_function_6d apply_gf(const real_function_6d &f,const size_t &particle)const;

	real_function_6d test_fill_tree()const{
		return real_factory_6d(world);
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
		CC_Timer init_bra(world,"Initialized molecular orbital bra_elements");
		return mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		init_bra.info();
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
					parameters.thresh_f12));
	/// The f12 convolution operator
	std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr
			< real_convolution_3d
			> (SlaterF12OperatorPtr(world, corrfac.gamma(),
					parameters.lo,
					parameters.thresh_f12));
	/// Intermediates (some need to be refreshed after every iteration)
	mutable CC_Intermediates intermediates_;

};

} /* namespace madness */

#endif /* CCOPERATORS_H_ */
