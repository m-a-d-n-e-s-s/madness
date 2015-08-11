/*
 * CCOperators.h
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */

////// TODO LIST
// Make two_electron integrals calculation parallel
// S5c diagram: Test alternative with perturbed density (not 3 for loops anymore) -> see time difference (but keep in mind that integral intermediate has to be calculated anyway for the exchange part
// Or do S5c and S5c_X together
#ifndef CCOPERATORS_H_
#define CCOPERATORS_H_

// Operators for coupled cluster and CIS
#include <chem/SCFOperators.h>
#include <chem/electronic_correlation_factor.h>
namespace madness {

// TAKEN FROM MP2.h
/// POD holding all electron pairs with easy access
template<typename T>
struct Pairs {

	typedef std::map<std::pair<int,int>, T> pairmapT;
	pairmapT allpairs;

	/// getter
	const T& operator()(int i, int j) const {
		return allpairs.find(std::make_pair(i, j))->second;
	}

	/// getter
	T& operator()(int i, int j) {
		return allpairs[std::make_pair(i, j)];
	}

	/// setter
	void insert(int i, int j, T pair) {
		std::pair<int, int> key = std::make_pair(i, j);
		allpairs.insert(std::make_pair(key, pair));
	}
};

static double unitfunction(const coord_3d &r) {return 1.0;}

// forward declaration
//class SCF;
//class Nemo;
//class NuclearCorrelationFactor;
//class XCfunctional;
//class Nuclear;
typedef std::vector<Function<double,3> > vecfuncT;


/// Structure that holds the CC intermediates and is able to refresh them
struct CC_Intermediates{
public:
	CC_Intermediates(World&world,const vecfuncT &bra, const vecfuncT &ket, const Nemo&nemo):
		world(world),
		mo_bra_(bra),
		mo_ket_(ket),
		poisson(std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, nemo.get_calc() -> param.lo,FunctionDefaults<3>::get_thresh()))),
		density_(make_density(bra,ket)),
		hartree_potential_(make_hartree_potential(density_)),
		exchange_intermediate_(make_exchange_intermediate(bra,ket)),
		integrals_hf_(make_two_electron_integrals_hf())
{

}


	/// Get the intermediates
	real_function_3d get_density(){return density_;}
	real_function_3d get_perturbed_density()const{return perturbed_density_;}
	real_function_3d get_hartree_potential()const{return hartree_potential_;}
	real_function_3d get_J(){return hartree_potential_;}
	real_function_3d get_perturbed_hartree_potential()const{return perturbed_hartree_potential_;}
	real_function_3d get_pJ(){return perturbed_hartree_potential_;}
	std::vector<vecfuncT> get_exchange_intermediate()const{return exchange_intermediate_;}
	std::vector<vecfuncT> get_EX(){return exchange_intermediate_;}
	std::vector<vecfuncT> get_perturbed_exchange_intermediate()const{return perturbed_exchange_intermediate_;}
	std::vector<vecfuncT> get_pEX()const{return perturbed_exchange_intermediate_;}
	Tensor<double> get_intergrals_hf()const{return integrals_hf_;}
	Tensor<double> get_integrals_mixed_t1()const{return integrals_mixed_t1_;}
	Tensor<double> get_integrals_t1()const{return integrals_t1_;}
	/// refresh the intermediates that depend on the \tau functions
	void update(const vecfuncT &tau){
		if(world.rank()==0) std::cout << "Update Intermediates:\n";
		perturbed_density_ = make_density(mo_bra_,tau);
		perturbed_hartree_potential_ = (*poisson)(perturbed_density_);
		perturbed_exchange_intermediate_ = make_exchange_intermediate(mo_bra_,tau);
		integrals_mixed_t1_ = make_two_electron_integrals_mixed_t1(tau);
		integrals_t1_ = make_two_electron_integrals_t1(tau);
	}

	/// make a density from two input functions
	/// For closed shell the density has to be scaled with 2 in most cases (this is not done here!)
	/// @param[in] vecfuncT_bra
	/// @param[in] vecfuncT_ket
	/// @param[out] \sum_i bra_i * ket_i
	real_function_3d make_density(const vecfuncT &bra,const vecfuncT &ket)const{
		if(bra.size()!=ket.size()) error("error in make density: unequal sizes (" + stringify(bra.size()) + " and " + stringify(ket.size()) + ")");
		if(bra.empty()) error("error in make_density: bra_element is empty");
		// make the density
		real_function_3d density = real_factory_3d(world);
		for(size_t i=0;i<bra.size();i++) density += bra[i]*ket[i];
		density.truncate();
		return density;
	}
	/// Poisson operator
	const std::shared_ptr<real_convolution_3d> poisson;



private:
	World &world;
	const vecfuncT &mo_bra_;
	const vecfuncT &mo_ket_;
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

	void error(const std::string &msg)const{
		std::cout << "\n\n\nERROR IN CC_INTERMEDIATES:\n" << msg << "\n\n\n!!!";
		MADNESS_EXCEPTION("\n\n!!!!ERROR IN CC_INTERMEDIATES!!!!\n\n\n\n\n\n\n\n\n\n\n\n",1);
	}

public:
	/// Make the exchange intermediate: EX[j][i] <bra[i](r2)|1/r12|ket[j](r2)>
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT &bra, const vecfuncT &ket)const{
		if(bra.size()!=ket.size() or bra.empty()) error("in make_exchange_intermediate, bra and ket empty or unequal sizes:\n bra_size: " + stringify(bra.size()) + ", ket_size: " + stringify(ket.size()));
		std::vector<vecfuncT> EX;
		EX.resize(bra.size());
		for(size_t i=0;i<bra.size();i++){
			EX[i].resize(ket.size());
			for(size_t j=0;j<ket.size();j++){
				EX[i][j] = (*poisson)(bra[j]*ket[i]);
			}
			truncate(world,EX[i]);
		}
		return EX;
	}
	/// Calculates the hartree potential Poisson(density)
	/// @param[in] density: a 3d function on which the poisson operator is applied (can be the occupied density and the perturbed density)
	/// @param[out] poisson(density) = \int 1/r12 density(r2) dr2
	real_function_3d make_hartree_potential(const real_function_3d &density)const{
		real_function_3d hartree = (*poisson)(density);
		hartree.truncate();
		return hartree;
	}

	/// Calculates two electron integrals
	/// <ij|g|kl>
	Tensor<double> make_two_electron_integrals_hf()const{
		Tensor<double> result(mo_bra_.size(),mo_bra_.size(),mo_ket_.size(),mo_ket_.size());
		for(size_t i=0;i<mo_bra_.size();i++){
			for(size_t j=0;j<mo_bra_.size();j++){
				for(size_t k=0;k<mo_ket_.size();k++){
					for(size_t l=0;l<mo_ket_.size();l++){
						result(i,j,k,l)=(mo_bra_[i]*mo_ket_[k]).inner(exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}
	/// <ij|g|k\tau_l>
	Tensor<double> make_two_electron_integrals_mixed_t1(const vecfuncT &tau)const{
		Tensor<double> result(mo_bra_.size(),mo_bra_.size(),mo_ket_.size(),tau.size());
		for(size_t i=0;i<mo_bra_.size();i++){
			for(size_t j=0;j<mo_bra_.size();j++){
				for(size_t k=0;k<mo_ket_.size();k++){
					for(size_t l=0;l<tau.size();l++){
						result(i,j,k,l)=(mo_bra_[i]*mo_ket_[k]).inner(perturbed_exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}
	// <ij|g|\tau_k \tau_l>
	Tensor<double> make_two_electron_integrals_t1(const vecfuncT &tau)const{
		Tensor<double> result(mo_bra_.size(),mo_bra_.size(),tau.size(),tau.size());
		for(size_t i=0;i<mo_bra_.size();i++){
			for(size_t j=0;j<mo_bra_.size();j++){
				for(size_t k=0;k<tau.size();k++){
					for(size_t l=0;l<tau.size();l++){
						result(i,j,k,l)=(mo_bra_[i]*tau[k]).inner(perturbed_exchange_intermediate_[l][j]);
					}
				}
			}
		}
		return result;
	}

	bool use_timer_;
	mutable double ttt, sss;
	void START_TIMER() const {
		if(use_timer_)world.gop.fence(); ttt=wall_time(); sss=cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if(use_timer_)END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if(use_timer_){
			ttt=wall_time()-ttt; sss=cpu_time()-sss;
			if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}
};


/// Coupled Cluster Operators (all closed shell)
class CC_Operators{
public:
	/// Constructor
	CC_Operators(World& world, const Nemo &nemo,const CorrelationFactor &correlationfactor): world(world),nemo(nemo),corrfac(correlationfactor), mo_bra_(make_mo_bra(nemo)), mo_ket_(nemo.get_calc()->amo), intermediates_(world,mo_bra_,mo_ket_,nemo),use_timer_(true)  {
	}

	vecfuncT get_CIS_potential(const vecfuncT &tau){
		START_TIMER();
		intermediates_.update(tau);
		END_TIMER("update intermediates");
		vecfuncT result = add(world,S3c(tau),S3c_X(tau));
		Q(result);
		return add(world,result,fock_residue_closed_shell(tau));
	}

	/// Projectors to project out the occupied space
	// 3D on vector of functions
	void Q(vecfuncT &f)const{
		for(size_t i=0;i<f.size();i++) Q(f[i]);
	}
	// 3D on single function
	void Q(real_function_3d &f)const{
		for(size_t i=0;i<mo_ket_.size();i++){
			f -= mo_bra_[i].inner(f)*mo_ket_[i];
		}
	}

	/// CCSD/CC2 singles potential parts

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_i = \sum_k <k|r12|k> |tau_i>
	// K_i = \sum_k <k|r12|tau_i> |k>
	vecfuncT fock_residue_closed_shell(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT J = mul(world,intermediates_.get_hartree_potential(),tau);
		truncate(world,J);
		scale(world,J,2.0);
		END_TIMER("J");
		START_TIMER();
		vecfuncT K;
		for(size_t i=0;i<tau.size();i++){
			real_function_3d tmp= real_factory_3d(world);
			vecfuncT vectmp = mul(world, intermediates_.get_perturbed_exchange_intermediate()[i],mo_ket_);
			for(size_t j=0;j<tau.size();j++) tmp += vectmp[j];
			tmp.truncate();
			K.push_back(tmp);
		}
		truncate(world,K);
		scale(world,K,-1);
		END_TIMER("K");
		return add(world,J,K);
	}

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3c(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result = mul(world,intermediates_.get_perturbed_hartree_potential(),mo_ket_);
		Q(result);
		truncate(world,result);
		scale(world,result,2.0);
		END_TIMER("S3c");
		return result;
	}
	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	//	vecfuncT S3c(const vecfuncT &tau)const{
	//		START_TIMER();
	//		vecfuncT result = mul(world,(*intermediates_.poisson)(intermediates_.make_density(mo_bra_,tau)),mo_ket_);
	//		Q(result);
	//		truncate(world,result);
	//		scale(world,result,2.0);
	//		END_TIMER("S3C_C");
	//		return result;
	//	}

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = -Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3c_X(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result;
		for(size_t i=0;i<tau.size();i++){
			real_function_3d tmp= real_factory_3d(world);
			vecfuncT vectmp = mul(world, intermediates_.get_exchange_intermediate()[i],tau);
			for(size_t j=0;j<tau.size();j++) tmp += vectmp[j];
			tmp.truncate();
			result.push_back(tmp);
		}
		Q(result);
		truncate(world,result);
		scale(world,result,-1.0);
		END_TIMER("S3c_X");
		return result;
	}

	/// The S5b term
	//[i]    [Q]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// 2\sum_k <k|g|\tau_k> |\tau_i>
	// No Q is applied yet !
	vecfuncT S5b(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result = mul(world,intermediates_.get_perturbed_hartree_potential(),mo_ket_);
		truncate(world,result);
		scale(world,result,2.0);
		return result;
		END_TIMER("S5b");
	}

	/// The S5b Exchange Term
	//[i]         [Q]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_k <k|g|\tau_i> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5b_X(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT tmp;
		vecfuncT result= zero_functions_compressed<double,3>(world,(tau.size()));
		for(size_t i=0;i<tau.size();i++){
			tmp = mul(world,intermediates_.get_perturbed_exchange_intermediate()[i],tau);
			for(size_t k=0;k<tau.size();k++){
				result[i] += tmp[k];
			}
		}
		truncate(world,result);
		scale(world,result,-1);
		END_TIMER("S5b_X");
		return result;
	}

	/// The S5c term
	//[Q]    [i]
	// \     /....
	//  \   /   / \
	//  _\_/_  _\_/_
	// -2\sum_kl <kl|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	// May use alteriative algorithm with perturbed density intermediate
	vecfuncT S5c(const vecfuncT&tau)const{
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double,3>(world,tau.size());
		for(size_t i=0;i<mo_bra_.size();i++){
			for(size_t k=0;k<mo_bra_.size();k++){
				for(size_t l=0;l<mo_bra_.size();l++){
					result[i]+= intermediates_.get_integrals_mixed_t1()(k,l,i,l)*tau[k];
				}
			}
		}
		truncate(world,result);
		scale(world,result,-2.0);
		END_TIMER("S5c");
		return result;
	}

	/// The S5c_X echange term
	//[Q]         [i]
	// \     ...../
	//  \   /\   /
	//  _\_/  \_/_
	// -\sum_kl <lk|g|i\tau_l> |\tau_k>
	// No Q is applied yet !
	vecfuncT S5c_X(const vecfuncT&tau)const{
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double,3>(world,tau.size());
		for(size_t i=0;i<mo_bra_.size();i++){
			for(size_t k=0;k<mo_bra_.size();k++){
				for(size_t l=0;l<mo_bra_.size();l++){
					result[i]+= intermediates_.get_integrals_mixed_t1()(l,k,i,l)*tau[k];
				}
			}
		}
		truncate(world,result);
		scale(world,result,-1.0);
		END_TIMER("S5c_X");
		return result;
	}

	/// The S6+X Term
	// \    /\    /...
	//  \  /  \  /   /\
	//  _\/_  _\/_  _\/_
	// -Q \sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>
	// Q is not applied yet!
	vecfuncT S6(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result = zero_functions_compressed<double,3>(world,tau.size());
		for(size_t i=0;i<tau.size();i++){
			for(size_t k=0;k<mo_bra_.size();k++){
				for(size_t l=0;l<mo_bra_.size();l++){
					result[i]+=(-2*intermediates_.get_integrals_t1()(k,l,k,i) - intermediates_.get_integrals_t1()(k,l,i,k))*tau[l];
				}
			}
		}
		truncate(world,result);
		END_TIMER("S6+X");
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
	/// Current procedure:
	/// use g12 = \int \delta(1-3) g32 d3
	/// <k(2)|g12|u(1,2)> = \int d2[ g12x(1,2 ] with x(1,2) = k(2)u(1,2)
	/// = int d2 [ int d3[ \delta(1-3) g32 ] x(1,2) ]
	/// = \int d3[\delta(1-3) \int d2 [ g32 x(1,2 ] ]
	/// = \int d3[\delta(1-3) h(1,3)] with h(1,3) = \int d2 g23 x(1,2)
	vecfuncT S2b(const Pairs<real_function_6d> u)const{
		double prefactor = 1.0/(2.0*corrfac.gamma());
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(mo_ket_.size());
		for(size_t i=0;i<mo_ket_.size();i++){
			real_function_3d resulti = real_factory_3d(world);
			for(size_t k=0;k<mo_ket_.size();k++){
				// calculate x(1,2) from u(1,2) and k(2), --> F.A.B uses multiply(copy(f),copy(bra) ...) deep copy of functions (dont know why)
				real_function_6d xik = multiply(u(i,k),mo_bra_[k],2);
				real_function_6d xki = multiply(u(k,i),mo_bra_[k],2);
				// calculate the convolution with fg = 1/(2gamma)*(Coulomb - BSH(gamma))
				real_function_6d hik;
				real_function_6d hki;
				{
					real_function_6d CoulombTerm = ((*poisson)(xik)).truncate();
					real_function_6d fBSHTerm = ((*fBSH)(xik)).truncate();
					hik = prefactor*(CoulombTerm + fBSHTerm);
				}{
					real_function_6d CoulombTerm = ((*poisson)(xki)).truncate();
					real_function_6d fBSHTerm = ((*fBSH)(xki)).truncate();
					hki = prefactor*(CoulombTerm + fBSHTerm);
				}
				// Make the projection to 3D with the unit function
				real_function_3d resultik = hik.project_out(unity,3);
				real_function_3d resultki = hki.project_out(unity,3);
				resultik.truncate();
				resulti += (2.0*resultik - resultki);
			}
			result[i]=resulti;
		}
		Q(result);
		truncate(world,result);
		return result;
	}

	/// S2c + X Term
	// [Q]   [i]
	//  \    /....
	//   \  /    /\
	//  __\/_____\/__
	/// = \sum <k(3)l(4)|g34 f31| u_{lk}(1,3) i(4)> = \sum <k(3)|f13| X_{lk,li}(1,3) > with X_{lk,li}(1,3) = u_lk(1,3) * (<l(4)|g34>|i(4)>_4)(3)
	/// = \sum \int d5 \delta(5-1) \int d3 f53 k(3)*X_{lk,li}(5,3)
	vecfuncT S2c(const Pairs<real_function_6d> u)const{
		real_function_3d unity = real_factory_3d(world).f(unitfunction);
		vecfuncT result(mo_ket_.size());
		for(size_t i=0;i<mo_ket_.size();i++){
			real_function_3d resulti = real_factory_3d(world);
			for(size_t k=0;k<mo_ket_.size();i++){
				for(size_t l=0;l<mo_ket_.size();l++){
				// make X_lkli(5,3) = ulk(5,3) * (<l(4)|g34>|i(4)>_4)(3) and X
				real_function_6d xlkli = multiply(u(l,k),intermediates_.get_exchange_intermediate()[i][l],2);
				real_function_6d xlkki = multiply(u(l,k),intermediates_.get_exchange_intermediate()[i][k],2);
				// make Atmp =  k(3) * X_lkli(5,3) and X
				real_function_6d Aklkli = multiply(xlkli,mo_bra_[k],2);
				real_function_6d Allkki = multiply(xlkki,mo_bra_[l],2);
				// make f12 convolution
				real_function_6d tmpklkli= (*f12op)(Aklkli);
				real_function_6d tmpllkki= (*f12op)(Allkki);
				// Project to 3d
				resulti -= (2.0*tmpklkli.project_out(unity,0)- tmpllkki.project_out(unity,0));
				}
			}
			resulti.truncate();
			result[i]=resulti;
		}
		Q(result);
		return result;
	}

	/// The S4a + X diagram
	//[Q]       [i]
	// \    ..../.....
	//  \  /\  /     /\
	//  _\/_ \/______\/_
	/// -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k>
	vecfuncT s4a(const Pairs<real_function_6d> u, const vecfuncT & tau)const{
		vecfuncT result(mo_ket_.size());
		for(size_t i=0;i<mo_ket_.size();i++){
			real_function_3d resulti = real_factory_3d(world);
			for(size_t k=0;k<mo_ket_.size();k++){
				for(size_t l=0;l<mo_ket_.size();l++){
					// Coulomb Part of f12g12 = g12 - BSH
					real_function_6d eri = TwoElectronFactory(world).dcut(FunctionDefaults<3>::get_thresh());
					real_function_6d kl_g =
							CompositeFactory<double, 6, 3>(world).particle1(
									copy(mo_bra_[k])).particle2(
									copy(mo_bra_[l])).g12(eri);
					resulti -= (2.0*inner(u(i,l),kl_g)*tau[k] - inner(u(i,k),kl_g)*tau[l]);
					resulti.truncate();
					// BSH part of f12g12
					// -> need to add this to the two-electron-factories
				}
			}
			result[i]=resulti;
		}
		Q(result);
		return result;
	}
private:
	/// The World
	World &world;
	/// Nemo
	const Nemo &nemo;
	/// Electronic correlation factor
	CorrelationFactor corrfac;
	/// The ket and the bra element of the occupied space
	/// if a  nuclear correlation factor is used the bra elements are the MOs multiplied by the squared nuclear correlation factor (done in the constructor)
	const vecfuncT mo_bra_;
	const vecfuncT mo_ket_;
	/// Helper function to initialize the const mo_bra and ket elements
	vecfuncT make_mo_bra(const Nemo &nemo)const{
		START_TIMER();
		return mul(world,nemo.nuclear_correlation -> square(),nemo.get_calc() -> amo);
		END_TIMER("Initialized molecular orbital bra_elements");
	}
	/// The poisson operator (Coulomb Operator)
	std::shared_ptr<real_convolution_3d> poisson = std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, nemo.get_calc() -> param.lo,FunctionDefaults<3>::get_thresh()));
	/// The BSH Operator for the f12g12 convolution which is with f12= 1/(2gamma)[1-exp(-gamma*r12)], f12g12 = 1/(2gamma) [CoulombOp - BSHOp(gamma)]
	std::shared_ptr<real_convolution_3d> fBSH = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world,corrfac.gamma(),nemo.get_calc()->param.lo,FunctionDefaults<3>::get_thresh()));
	/// The f12 convolution operator
	std::shared_ptr<real_convolution_3d> f12op = std::shared_ptr<real_convolution_3d>(SlaterF12OperatorPtr(world,corrfac.gamma(),nemo.get_calc()->param.lo,FunctionDefaults<3>::get_thresh()));
	/// Intermediates (some need to be refreshed after every iteration)
	CC_Intermediates intermediates_;

	/// Take the time of operations
	bool use_timer_;
	mutable double ttt, sss;
	void START_TIMER() const {
		if(use_timer_)world.gop.fence(); ttt=wall_time(); sss=cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if(use_timer_)END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if(use_timer_){
			ttt=wall_time()-ttt; sss=cpu_time()-sss;
			if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}
};

/// Coupled Cluster Operators that use and result in only 3D functions (so CCS and CIS)
class CC_3D_Operator{
public:
	//	CC_3D_Operator(World&world, const Nemo &nemo): world(world), mo_ket_(nemo.get_calc() -> amo),R2(init_R2(nemo.nuclear_correlation -> square())){
	//		poisson = std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, nemo.get_calc() -> param.lo, nemo.get_calc() ->param.econv));
	//		mo_bra_ = mul(world,R2,mo_ket_);
	//		exchange_intermediate_ = make_exchange_intermediate();
	//		if(mo_bra_.empty()) std::cout << "\n\n!!!!!WARNING: mo_bra_ vector is empty!!!!!\n\n";
	//		if(mo_ket_.empty()) std::cout << "\n\n!!!!!WARNING: mo_ket_ vector is empty!!!!!\n\n";
	//		//madness::Nuclear U(world,nemo);
	//		//nuclear_potential_ =U;
	//	}
	// If other mos than the one in the nemo struct are needed (e.g. if lower thresh is demanded -> guess calculations)
	CC_3D_Operator(World&world, const Nemo &nemo,const vecfuncT &mos): world(world),use_nuclear_correlation_factor_(true), mo_ket_(mos),R2(init_R2(nemo)){
		if(nemo.nuclear_correlation -> type() == NuclearCorrelationFactor::None){
			std::cout << "No nuclear correlation factor used" << std::endl;
			use_nuclear_correlation_factor_ = false;
		}
		poisson = std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, nemo.get_calc() -> param.lo,FunctionDefaults<3>::get_thresh()));
		if(use_nuclear_correlation_factor_)mo_bra_ = mul(world,R2,mo_ket_);
		else mo_bra_ = mo_ket_;
		dR2 = get_gradient(R2);
		plot_plane(world,dR2[0],"dxR2");
		set_thresh(world,mo_bra_,FunctionDefaults<3>::get_thresh());
		set_thresh(world,mo_ket_,FunctionDefaults<3>::get_thresh());
		truncate(world,mo_ket_);
		truncate(world,mo_bra_);
		exchange_intermediate_ = make_exchange_intermediate();
		sanitycheck();
		//		START_TIMER();
		//		CC_Operators tmp(world,nemo);
		//		END_TIMER("CC TEST INIT OPERATORS");
	}


	/// Make shure that R2 gets the right thresh and is constant
	real_function_3d init_R2(const Nemo &nemo )const{
		if(nemo.nuclear_correlation){
			real_function_3d tmp = copy(nemo.nuclear_correlation -> square());
			tmp.set_thresh(FunctionDefaults<3>::get_thresh());
			tmp.truncate();
			tmp.verify();
			return tmp;
		}
		real_function_3d constant = real_factory_3d(world);
		return (constant +1.0);
	}

	// Make the derivative of R2
	vecfuncT get_gradient(const real_function_3d f)const{
		std::vector < std::shared_ptr<real_derivative_3d> > gradop=gradient_operator<double, 3>(world);
		f.verify();
		vecfuncT gradf;
		for(size_t i=0;i<3;i++){
			real_function_3d dfi = (*gradop[i])(f);
			gradf.push_back(dfi);
		}
		set_thresh(world,gradf,FunctionDefaults<3>::get_thresh());
		truncate(world,gradf);
		return gradf;
	}


	void sanitycheck()const{
		if(mo_ket_.empty()) error("mo_ket_ is empty");
		if(mo_bra_.empty()) error("mo_bra_ is empty");
		if(exchange_intermediate_.empty()) error("exchange intermediate is empty");
		for(auto x: mo_ket_){
			if(x.thresh() != FunctionDefaults<3>::get_thresh()) error("Wrong thresh in mo_ket_ functions");
		}
		for(auto x: mo_bra_){
			if(x.thresh() != FunctionDefaults<3>::get_thresh()) error("Wrong thresh in mo_bra_ functions");
		}
		for(auto tmp: exchange_intermediate_){
			if(tmp.empty()) error("Exchange Intermediate contains empty vectors");
			for(auto x: tmp){
				if(x.thresh() != FunctionDefaults<3>::get_thresh()) error("Wrong thresh in Exchange Intermediate");
			}
		}
		R2.verify();
	}

	void memory_information(const vecfuncT &v, const std::string &msg = "vectorfunction size is: ")const{
		const double x = get_size(world,v);
		if(world.rank()==0) std::cout << msg << "("<< x <<" GB)" << " for " << v.size() <<" functions\n";
	}

	std::vector<vecfuncT> make_exchange_intermediate()const{
		std::vector<vecfuncT> intermediate;
		double memory = 0.0;
		for(size_t i=0;i<mo_bra_.size();i++){
			const vecfuncT integrant = mul(world,mo_bra_[i],mo_ket_);
			const vecfuncT intermediate_i = apply(world,*poisson,integrant);
			intermediate.push_back(intermediate_i);
			memory += get_size(world,intermediate_i);
		}
		if(world.rank()==0) std::cout << "Created exchange intermediate of dimension " << intermediate.size() << "x" << intermediate.front().size() << " and size (" << memory <<   " GB)\n";
		return intermediate;

	}


public:
	// The nuclear potential is missing (or the U potential for the regularized approach)

	// Closed Shell Triplet CIS potential without the nuclear potential
	// returns (2J - K)x + S3C_X
	// which is in components VCIS_j =  2*\sum_i <i|r12|i> |x_j> - \sum_i <i|r12|x_j> |i> - Q\sum_i <i|r12|j> |x_i>
	vecfuncT get_CIS_potential_triplet(const vecfuncT &x)const{
		return add(world,fock_residue_closed_shell(x),S3C_X(x));
	}
	// Closed Shell Singlet CIS potential without the nuclear potential
	// returns (2J - K)x + 2*S3C_C + S3C_X
	// which is in components VCIS_j =  2*\sum_i <i|r12|i> |x_j> - \sum_i <i|r12|x_j> |i> + 2*Q\sum_i <i|r12|x_i> |j> - Q\sum_i <i|r12|j> |x_i>
	vecfuncT get_CIS_potential_singlet(const vecfuncT &x,const Nemo &nemo)const{
		// test the CC2 singles
		//		START_TIMER();
		//		CC_Operators tmpops(world,nemo);
		//		END_TIMER("INitialize CCOPS ... not necessary ...");
		//		return tmpops.get_CIS_potential(x);

		vecfuncT S3CC = S3C_C(x);
		scale(world,S3CC,2.0);
		vecfuncT S3CX = S3C_X(x);
		return add(world,fock_residue_closed_shell(x),add(world,S3CX,S3CC));
	}

	// Closed Shell potential for Virtuals (or SCF MOs)
	vecfuncT get_SCF_potential(const vecfuncT &x)const{
		return fock_residue_closed_shell(x);
	}


	// get the ground state density
	real_function_3d make_density()const{
		return make_density(mo_bra_,mo_ket_);
	}

	// Make a density out of two vectorfunctions f and g
	// density = \sum_i |f_i><g_i|
	real_function_3d make_density(const vecfuncT &f, const vecfuncT &g)const{
		if(f.size() != g.size()) error("make_density: sizes of vectors are not equal");
		real_function_3d density = real_factory_3d(world);
		for(size_t i=0;i<f.size();i++) density += f[i]*g[i];
		density.truncate();
		return density;
	}

	// The Fock operator is partitioned into F = T + Vn + R
	// the fock residue R= 2J-K for closed shell is computed here
	// J_j = \sum_i <i|r12|i> |tau>
	// K_j = \sum_i <i|r12|tau_j> |i>
	vecfuncT fock_residue_closed_shell(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT J = mul(world,(*poisson)(make_density()),tau);
		truncate(world,J);
		scale(world,J,2.0);
		END_TIMER("J");
		START_TIMER();
		vecfuncT K;
		for(size_t j=0;j<tau.size();j++){
			real_function_3d Kj = real_factory_3d(world);
			for(size_t i=0;i<mo_bra_.size();i++){
				Kj += (*poisson)(mo_bra_[i]*tau[j])*mo_ket_[i];
			}
			K.push_back(Kj);
		}
		truncate(world,K);
		scale(world,K,-1);
		END_TIMER("K");
		return add(world,J,K);
	}

	// The same residue for the case that the Fock operator is the Kohn-Sham Operator
	vecfuncT KS_residue_closed_shell (const std::shared_ptr<SCF> scf,const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT J = mul(world,(*poisson)(make_density()),tau);
		truncate(world,J);
		scale(world,J,2.0);
		END_TIMER("J");
		START_TIMER();
		XCOperator xcoperator(world,scf.get(),0);
		real_function_3d vxc=xcoperator.make_xc_potential();
		vxc.truncate();
		vecfuncT applied_vxc = mul(world,vxc,tau);
		truncate(world,applied_vxc);
		END_TIMER("Vxc");
		return add(world,J,applied_vxc);
	}

	// Kinetik energy
	// -1/2 <x|R2Nabla2|x> = +1/2 <Nabla R2 x | Nabla x> = grad(R2x)*grad(x)
	// grad(R2x) = GradR2*x + R2*gradx
	// grad(R2x)*grad(y) = GradR2*x*Grady + R2*Gradx*Grady
	double get_matrix_element_kinetic_energy(const vecfuncT &ket, const vecfuncT &bra)const{
		//std::cout << " Making Kintic Energy Matrix Element 1:\n";
		double value=0.0;
		vecfuncT R2bra = mul(world,R2,bra);
		truncate(world,R2bra);
		std::vector<vecfuncT> dket, dbra;
		std::vector < std::shared_ptr<real_derivative_3d> > gradop;
		gradop = gradient_operator<double, 3>(world);
		for(size_t axis=0;axis<3;axis++){
			START_TIMER();
			const vecfuncT gradbra = apply(world,*gradop[axis],R2bra);
			END_TIMER("Gradient of R2Bra");
			START_TIMER();
			const vecfuncT gradket = apply(world,*gradop[axis],ket);
			END_TIMER("Gradient of Ket");
			START_TIMER();
			value += 0.5*inner(world,gradbra,gradket).sum();
			END_TIMER("Inner Product");
		}
		return value;
	}
	double get_matrix_element_kinetic_2(const vecfuncT &bra,const vecfuncT &ket)const{
		std::cout << "Making Kinetic 2 Element\n";
		double value =0.0;
		std::vector < std::shared_ptr<real_derivative_3d> > gradop = gradient_operator<double, 3>(world);
		for(size_t axis=0;axis<3;axis++){
			START_TIMER();
			vecfuncT gradbra = apply(world,*gradop[axis],bra);
			truncate(world,gradbra);
			END_TIMER("make gradbra");
			START_TIMER();
			vecfuncT gradket = apply(world,*gradop[axis],ket);
			truncate(world,gradket);
			END_TIMER("make gradket");
			START_TIMER();
			vecfuncT gradR2bra = mul(world,dR2[axis],bra);
			truncate(world,gradR2bra);
			END_TIMER("multiply dR2 and bra");
			START_TIMER();
			value += (inner(world,gradR2bra,gradket).sum());
			END_TIMER("Inner product1");
			START_TIMER();
			value += make_inner_product(gradbra,gradket);
			END_TIMER("Inner product2");

		}
		return 0.5*value;
	}

	// Kinetic part of the CIS perturbed fock matrix
	Tensor<double> get_matrix_kinetic(const std::vector<vecfuncT> &x)const{
		Tensor<double> result(x.size(),x.size());
		// make the x,y and z parts of the matrix and add them
		std::vector < std::shared_ptr<real_derivative_3d> > gradop= gradient_operator<double, 3>(world);
		for(size_t axis=0;axis<3;axis++){
			// make all gradients
			std::vector<vecfuncT> dx,dR2x;
			size_t i=0;
			for(auto xi:x){
				const vecfuncT dxi = apply(world, *(gradop[axis]), xi);
				dx.push_back(dxi);
				if(use_nuclear_correlation_factor_){
					const vecfuncT dR2xi = apply(world, *(gradop[axis]), mul(world,R2,xi));
					plot_plane(world,dR2xi.back(),"R2xi"+stringify(i)+"_"+stringify(axis));
					dR2x.push_back(dR2xi);
				}i++;
			}
			for(size_t i=0;i<x.size();i++){
				for(size_t j=0;j<x.size();j++){
					if(use_nuclear_correlation_factor_) result(i,j) += 0.5*inner(world,dR2x[j],dx[i]).sum();
					else result(i,j) += 0.5*inner(world,dx[j],dx[i]).sum();
				}
			}

		}
		return result;
	}

	// Diagrammatic Potentials:

	// The coulomb Term of the S3C diagram: Positive sign
	// \     /
	//  \---/  = Q\sum_j(<j|g12|tau_j>)|i>
	//  _\_/_
	vecfuncT S3C_C(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result = mul(world,(*poisson)(make_density(mo_bra_,tau)),mo_ket_);
		Q(result);
		truncate(world,result);
		END_TIMER("S3C_C");
		return result;
	}

	// The Exchange Term of the S3C diagram: Negative sign
	// \  /
	//  \/...   = Q\sum_j(<j|g12|i>|tau_j>)
	//     / \
	//    _\_/_
	vecfuncT S3C_X(const vecfuncT &tau)const{
		START_TIMER();
		vecfuncT result;
		for(size_t i=0;i<tau.size();i++){
			real_function_3d tmp= real_factory_3d(world);
			vecfuncT vectmp = mul(world, exchange_intermediate_[i],tau);
			for(size_t j=0;j<tau.size();j++) tmp += vectmp[j];
			tmp.truncate();
			result.push_back(tmp);
		}
		Q(result);
		truncate(world,result);
		scale(world,result,-1.0);
		END_TIMER("S3C_X");
		return result;
	}

	// Project out the occupied space
	void Q(vecfuncT &f)const{
		for(size_t i=0;i<f.size();i++) Q(f[i]);
	}
	void Q(real_function_3d &f)const{
		for(size_t i=0;i<mo_ket_.size();i++){
			f -= mo_bra_[i].inner(f)*mo_ket_[i];
		}
	}

	// Make an inner product between vecfunctions
	double make_inner_product(const vecfuncT &bra, const vecfuncT &ket)const{
		if(use_nuclear_correlation_factor_)return inner(world,mul(world,R2,bra),ket).sum();
		else return inner(world,bra,ket).sum();
	}
	// inner product between functions
	double make_inner_product(const real_function_3d &bra, const real_function_3d &ket)const{
		if(use_nuclear_correlation_factor_)return (bra*R2).inner(ket);
		else return bra.inner(ket);
	}
	// inner product between function and vecfunction
	double make_inner_product(const real_function_3d &bra, const vecfuncT &ket)const{
		if(use_nuclear_correlation_factor_)return inner(world,bra*R2,ket).sum();
		else return inner(world,bra,ket).sum();
	}

private:
	bool use_timer_=true;
	World &world;
	bool use_nuclear_correlation_factor_;
	vecfuncT mo_bra_,mo_ket_;
	/// The squared nuclear correlation factor and its derivative;
	const real_function_3d R2;
	vecfuncT dR2;
	std::vector<vecfuncT> exchange_intermediate_;
	std::shared_ptr<real_convolution_3d> poisson;
	//	Nuclear nuclear_potential_;
	void error(const std::string &msg)const{
		std::cout << "\n\n\n !!!! ERROR IN CC_3D_OPERATOR CLASS:\n ERROR MESSAGE IS: " << msg <<"\n";
		MADNESS_EXCEPTION("!!!!ERROR IN CC_3D_OPERATOR CLASS!!!!",1);
	}
	// Timer
	mutable double ttt, sss;
	void START_TIMER() const {
		if(use_timer_)world.gop.fence(); ttt=wall_time(); sss=cpu_time();
	}

	void END_TIMER(const std::string msg) const {
		if(use_timer_)END_TIMER(msg.c_str());
	}

	void END_TIMER(const char* msg) const {
		if(use_timer_){
			ttt=wall_time()-ttt; sss=cpu_time()-sss;
			if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
		}
	}

};


} /* namespace madness */

#endif /* CCOPERATORS_H_ */
