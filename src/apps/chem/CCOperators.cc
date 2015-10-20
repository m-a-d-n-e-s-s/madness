/*
 * CCOperators.cc
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */

#include "CCOperators.h"

namespace madness{

/// save a function
template<typename T, size_t NDIM>
void CC_Operators::save_function(const Function<T, NDIM>& f,
		const std::string name) const {
	if (world.rank() == 0)
		print("saving function", name);
	f.print_size(name);
	archive::ParallelOutputArchive ar(world, name.c_str(), 1);
	ar & f;
}

real_function_3d CC_Intermediates::make_density(const vecfuncT &bra,
		const CC_vecfunction &ket) const {
	if (bra.empty()) error("error in make_density: bra_element is empty");
	if (ket.empty()) error("error in make_density: bra_element is empty");
	// make the density
	real_function_3d density = real_factory_3d(world);
	for (auto x:ket.functions) density += (bra[x.i] * x.function).truncate();
	density.truncate();
	return density;
}

intermediateT CC_Intermediates::make_exchange_intermediate(const vecfuncT &bra,
		const CC_vecfunction &ket) const {
//	if (bra.size() != ket.size() or bra.empty())
//		error("in make_exchange_intermediate, bra and ket empty or unequal sizes:\n bra_size: "
//				+ stringify(bra.size()) + ", ket_size: "
//				+ stringify(ket.size()));
//	std::vector<vecfuncT> EX;
//	EX.resize(bra.size());
//	for (size_t i = 0; i < bra.size(); i++) {
//		EX[i].resize(ket.size());
//		for (size_t j = 0; j < ket.size(); j++) {
//			EX[i][j] = (*poisson)(bra[j] * ket[i]);
//		}
//		truncate(world, EX[i]);
//	}
//	return EX;
	intermediateT xim;
	for(size_t k=0;k<bra.size();k++){
		for(auto l:ket.functions){
			real_function_3d kl = (bra[k]*l.function).truncate();
			real_function_3d result = ((*poisson)(kl)).truncate();
			xim.insert(k,l.i,result);
		}
	}
	return xim;
}

/// Calculates two electron integrals
/// <ij|g|kl>
Tensor<double> CC_Intermediates::make_two_electron_integrals_hf() const {
	Tensor<double> result(mo_bra_.size(), mo_bra_.size(), mo_ket_.size(),
			mo_ket_.size());
	return result;
}
/// <ij|g|k\tau_l>
Tensor<double> CC_Intermediates::make_two_electron_integrals_mixed_t1(
		const vecfuncT &tau) const {
	Tensor<double> result(mo_bra_.size(), mo_bra_.size(), mo_ket_.size(),
			tau.size());
	return result;
}
// <ij|g|\tau_k \tau_l>
Tensor<double> CC_Intermediates::make_two_electron_integrals_t1(const vecfuncT &tau) const {
	Tensor<double> result(mo_bra_.size(), mo_bra_.size(), tau.size(),
			tau.size());
	return result;
}

double CC_Operators::compute_mp2_pair_energy(CC_Pair &pair)const{

	// this will be the bra space
	real_function_6d eri = TwoElectronFactory(world).dcut(parameters.lo);
	real_function_6d ij_g =
			CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[pair.i])).particle2(
							copy(mo_bra_[pair.j])).g12(eri);
	real_function_6d ji_g =
			CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[pair.i])).particle2(
							copy(mo_bra_[pair.j])).g12(eri);

	// compute < ij | g12 | psi >
	const double ij_g_uij = inner(pair.function, ij_g);
	if (world.rank() == 0)
		printf("<ij | g12       | psi^1>  %12.8f\n", ij_g_uij);

	if(parameters.debug){
		if(world.rank()==0){
			std::cout << "Debugging make_ijgu function with mp2 pair energy\n";
		}
		const double ijguij = make_ijgu(pair.i,pair.j,pair);
		if(fabs(ijguij-ij_g_uij)>FunctionDefaults<6>::get_thresh()) warning("make_ijgu and mp2 pair energy function give not the same value "+ stringify(ijguij) + " vs " + stringify(ij_g_uij));
		else if(world.rank()==0) std::cout << "make_ijgu function seems to be fine values are: " << ijguij << " and " << ij_g_uij << std::endl;
	}

	// compute < ji | g12 | psi > if (i/=j)
	const double ji_g_uij = (pair.i == pair.j) ? 0 : inner(pair.function, ji_g);
	if (world.rank() == 0)
		printf("<ji | g12       | psi^1>  %12.8f\n", ji_g_uij);

	// the singlet and triplet triplet pair energies
	if (pair.i == pair.j) {
		pair.e_singlet = ij_g_uij + pair.ij_gQf_ij;
		pair.e_triplet = 0.0;
	} else {
		pair.e_singlet = (ij_g_uij + pair.ij_gQf_ij)
																								+ (ji_g_uij + pair.ji_gQf_ij);
		pair.e_triplet = 3.0
				* ((ij_g_uij - ji_g_uij) + (pair.ij_gQf_ij - pair.ji_gQf_ij));
	}

	// print the pair energies
	if (world.rank() == 0) {
		printf("current energy %2d %2d %12.8f %12.8f\n", pair.i, pair.j,
				pair.e_singlet, pair.e_triplet);
	}

	// return the total energy of this pair
	return pair.e_singlet + pair.e_triplet;
}

// The Fock operator is partitioned into F = T + Vn + R
// the fock residue R= 2J-K for closed shell is computed here
// J_i = \sum_k <k|r12|k> |tau_i>
// K_i = \sum_k <k|r12|tau_i> |k>
vecfuncT CC_Operators::fock_residue_closed_shell(const CC_vecfunction &singles) const {
	vecfuncT tau = singles.vec();
	CC_Timer timer_J(world,"J");
	vecfuncT J = mul(world, intermediates_.get_hartree_potential(), tau);
	truncate(world, J);
	scale(world, J, 2.0);
	timer_J.info();
	CC_Timer timer_K(world,"K");
	vecfuncT vK;
	for(auto i:singles.functions){
		real_function_3d Ki = K(i.function);
		vK.push_back(Ki);
	}
	truncate(world, vK);
	scale(world, vK, -1);
	timer_K.info();
	return add(world, J, vK);
}

// The coulomb Term of the S3C diagram: Positive sign
// \     /
//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
//  _\_/_
vecfuncT CC_Operators::S3c(const CC_vecfunction &singles, const std::string &msg) const {
	MADNESS_ASSERT(active_mo_.size()==singles.size());
	CC_Timer timer(world,"time");
	CC_data data("S3c");
	vecfuncT result;
	for(auto i:singles.functions){
		real_function_3d tmp = (intermediates_.get_perturbed_hartree_potential()*mo_ket_[i.i]).truncate();
		result.push_back(tmp);
	}
	Q(result);
	truncate(world, result);
	scale(world, result, 2.0);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

// The Exchange Term of the S3C diagram: Negative sign
// \  /
//  \/...   = -Q\sum_k(<k|g12|i>|tau_k>)
//     / \
//    _\_/_
vecfuncT CC_Operators::S3c_X(const CC_vecfunction &singles, const std::string &msg) const {
	MADNESS_ASSERT(active_mo_.size()==singles.size());
	CC_Timer timer(world,"time");
	CC_data data("S3c_X");
	vecfuncT result;
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
	for(auto k:singles.functions){
		resulti += (intermediates_.get_EX(i.i,k.i)*k.function).truncate();
	}
	result.push_back(resulti);
	}
	Q(result);
	truncate(world, result);
	scale(world, result, -1.0);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S5b term
//[i]    [Q]
// \     /....
//  \   /   / \
//  _\_/_  _\_/_
// 2\sum_k <k|g|\tau_k> |\tau_i>
// No Q is applied yet !
vecfuncT CC_Operators::S5b(const CC_vecfunction &singles, const std::string &msg) const {
	MADNESS_ASSERT(active_mo_.size()==singles.size());
	CC_Timer timer(world,"time");
	CC_data data("S5b");
	vecfuncT result = mul(world,
			intermediates_.get_perturbed_hartree_potential(), active_mo_);
	truncate(world, result);
	scale(world, result, 2.0);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S5b Exchange Term
//[i]         [Q]
// \     ...../
//  \   /\   /
//  _\_/  \_/_
// -\sum_k <k|g|\tau_i> |\tau_k>
// No Q is applied yet !
vecfuncT CC_Operators::S5b_X(const CC_vecfunction &singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S5b_X");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			resulti += (intermediates_.get_pEX(k.i,i.i)*k.function).truncate();
		}
		result.push_back(resulti);
	}

	truncate(world, result);
	scale(world, result, -1);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
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
vecfuncT CC_Operators::S5c(const CC_vecfunction&singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S5c");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			for(auto l:singles.functions){
				real_function_3d tmp = (intermediates_.get_EX(k.i,i.i)*l.function).truncate();
				double integral = mo_bra_[l.i].inner(tmp);
				resulti += integral*k.function;
			}
		}
		result.push_back(resulti);
	}
	truncate(world, result);
	scale(world, result, -2.0);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S5c_X echange term
//[Q]         [i]
// \     ...../
//  \   /\   /
//  _\_/  \_/_
// -\sum_kl <lk|g|i\tau_l> |\tau_k>
// No Q is applied yet !
vecfuncT CC_Operators::S5c_X(const CC_vecfunction &singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S5c_X");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			for(auto l:singles.functions){
				real_function_3d tmp = (intermediates_.get_EX(l.i,i.i)*l.function).truncate();
				double integral = mo_bra_[k.i].inner(tmp);
				resulti += integral*k.function;
			}
		}
		result.push_back(resulti);
	}
	truncate(world, result);
	scale(world, result, -1.0);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S6+X Term
// \    /\    /...
//  \  /  \  /   /\
//  _\/_  _\/_  _\/_
// -Q (\sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>)
// Q is not applied yet!
vecfuncT CC_Operators::S6(const CC_vecfunction &singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S6");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			for(auto l:singles.functions){
				real_function_3d tmp_J = (intermediates_.get_pEX(k.i,k.i)*i.function).truncate();
				real_function_3d tmp_K = (intermediates_.get_pEX(k.i,i.i)*k.function).truncate();
				double integral_J = l.function.inner(tmp_J);
				double integral_K = l.function.inner(tmp_K);
				resulti += (2.0*integral_J - integral_K)*l.function;
			}
		}
		result.push_back(resulti);
	}
	truncate(world, result);
	scale(world,result,-1.0);
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
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

vecfuncT CC_Operators::S2b(const Pairs<CC_Pair> u, const CC_vecfunction & singles , const std::string &msg) const {
	CC_Timer timer(world,"time");
	CC_data data("S2b");
	CC_vecfunction t = make_t_intermediate(singles);
	real_function_3d tdensity = intermediates_.make_density(mo_bra_,t);
	vecfuncT result;
	for(auto i:t.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			{
				CC_Timer S2b2(world,msg+"<k|g|u> with dirac convolution");
				real_function_6d tmp = multiply(copy(u(i.i,k.i).function),copy(mo_bra_[k.i]),2);
				poisson -> particle() = 2;
				real_function_6d result = (*poisson)(tmp);
				real_function_3d r = result.dirac_convolution<3>();

				// Exchange Part
				real_function_6d tmpx = multiply(copy(u(i.i,k.i).function),copy(mo_bra_[k.i]),1);
				poisson -> particle() = 1;
				real_function_6d resultx= (*poisson)(tmpx);
				real_function_3d rx = resultx.dirac_convolution<3>();

				resulti += 2.0*r -rx;

				S2b2.info();
				if(parameters.debug){
					CC_Timer s2b_debug_time(world,msg+"s2b-debug-6D-part");
					double debug = make_ijgu(i.i,k.i,u(i.i,k.i));
					double debug2= r.inner(mo_bra_[i.i]);
					if(fabs(debug-debug2)>FunctionDefaults<6>::get_thresh()) warning("S2b potential: Different result for |u> part " + stringify(debug) + " , " + stringify(debug2),data);
					else if(world.rank()==0) std::cout << "6d part of S2b potential seems to be fine (debug integrals are:" << std::setprecision(parameters.output_prec) << debug << " and " << debug2 << std::endl;
					s2b_debug_time.info();
				}

			}
		}
		// The Part which operators on the t intermediate OPTIMIZE THE LOOPS (for helium debug not necessary)
		real_function_3d t1part;
		real_function_3d t1partx;
		{

			// 2<k|gf|t_it_k> - <k|gf|t_kt_i> // factors are added in the end
			real_function_3d unitpart = (apply_gf(tdensity)*i.function).truncate();
			real_function_3d unitpartx = real_factory_3d(world);
			for(auto k:t.functions) unitpartx = (apply_gf(mo_bra_[k.i]*i.function)*k.function).truncate();

			// make intermediates
			// <n|f|ti>
			vecfuncT nfti;
			for(auto mo:mo_bra_) nfti.push_back( (*f12op)((mo*i.function).truncate()) );


			// 2 \sum_n <k(2)|g12|n(1)><n(1)|f12|ti(1)tk(2)> = <k(2)|g12 nfti(2)|n(1)tk(2)> = poisson(k*nfti*tk)*n
			//   \sum_n <k(2)|g12|n(1)><n(1)|f12|tk(1)ti(2)> = <k(2)|g12 nftk(2)|n(1)ti(2)> = poisson(k*nftk*ti)*n
			// The sum over n also goes over the frozen orbitals
			real_function_3d O1part = real_factory_3d(world);
			real_function_3d O1partx = real_factory_3d(world);
			for(auto k:t.functions){
				for(size_t n=0;n<mo_ket_.size();n++){
					O1part += ((*poisson)(mo_bra_[k.i]*k.function*nfti[n])*mo_ket_[n]).truncate();
					real_function_3d nftk = ((*f12op)(mo_bra_[n]*k.function)).truncate();
					real_function_3d poisson_integrant = (mo_bra_[k.i]*i.function*nftk).truncate();
					O1partx+= ((*poisson)(poisson_integrant)*mo_ket_[n]).truncate();
				}
			}

			// 2 \sum_n <k(2)|g12|n(2)><n(2)|f12|ti(1)tk(2)> = <k(2)|g12 nftk(1) |ti(1)n(2)> = poisson(k*n)*(nftk*ti)
			//   \sum_n <k(2)|g12|n(2)><n(2)|f12|tk(1)ti(2)> = <k(2)|g12 nfti(1) |tk(1)n(2)> = poisson(k,n)*(nfti*tk)
			real_function_3d O2part = real_factory_3d(world);
			real_function_3d O2partx = real_factory_3d(world);
			for(auto k:t.functions){
				for(size_t n=0;n<mo_ket_.size();n++){
					real_function_3d kgn = intermediates_.get_EX(k.i,n);
					real_function_3d nftk = (*f12op)(mo_bra_[n]*k.function);
					real_function_3d nfti = (*f12op)(mo_bra_[n]*i.function);
					O2part   += (kgn*nftk*i.function).truncate();
					O2partx  += (kgn*nfti*k.function).truncate();
				}
			}

			// 2 \sum_nm <k(2)|g12|mn><mn|f12|ti(1)tk(2)> = <k|g|n>(1)|m(1)> *<mn|f|titk>
			//   \sum_nm <k(2)|g12|mn><mn|f12|tk(1)ti(2)> = <k|g|n>(1)|m(1)> *<mn|f|tkti>
			real_function_3d O12part = real_factory_3d(world);
			real_function_3d O12partx = real_factory_3d(world);
			for(auto k:t.functions){
				for(size_t n=0;n<mo_bra_.size();n++){
					for(size_t m=0;m<mo_ket_.size();m++){
						real_function_3d kgn = intermediates_.get_EX(k.i,n);
						real_function_3d nftk = (*f12op)(mo_bra_[n]*k.function);
						real_function_3d nfti = (*f12op)(mo_bra_[n]*i.function);
						double f = nftk.inner(mo_bra_[m]*i.function);
						double fx= nfti.inner(mo_bra_[m]*k.function);
						O12part += (f*kgn*mo_ket_[m]).truncate();
						O12partx +=(fx*kgn*mo_ket_[m]).truncate();
					}
				}
			}

			t1part = unitpart - O1part - O2part +O12part;
			t1partx = unitpartx - O1partx - O2partx + O12partx;

			if(parameters.debug){
				CC_Timer s2b_debug_time(world,msg+"s2b-debug-3D-part");
				// make \sum_k <ik|g12|titk> and compare
				double debug = 0.0;
				for(auto k:t.functions) debug += make_ijgQfxy(i.i,k.i,i.function,k.function);
				double debug2 = t1part.inner(mo_bra_[i.i]);
				if(fabs(debug-debug2)>FunctionDefaults<3>::get_thresh()) warning("S2b potential: Different result for t1 part " + stringify(debug) + " , " + stringify(debug2),data);
				else if(world.rank()==0) std::cout << "3d part of S2b potential seems to be fine (debug integrals are:" << debug << " and " << debug2 << std::endl;
				s2b_debug_time.info();
			}


			t1part.truncate();
		}

		resulti += 2.0*t1part-t1partx;
		result.push_back(resulti);
	}
	scale(world,result,-1.0);
	truncate(world,result);
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// S2c + X Term
// [Q]   [i]
//  \    /....
//   \  /    /\
//  __\/_____\/__
/// = -Q\sum_{kl}\left( 2<k|lgi|ulk> - <l|kgi|u_{lk}> + 2<k|lgiQ12f12|t_lt_k> - <l|kgiQ12f12|t_lt_k> \right)
/// Notation: 6D Integration over second particle, intermediates: lgi = <l|g|i> is the exchange intermediate
/// Notation: t are the t-intermediates: |t_i> = |i> + |\tau_i>
/// @param[in] All the current coupled cluster Pairs
/// @param[in] The coupled cluster singles
/// @param[out] the S2c+X Potential

vecfuncT CC_Operators::S2c(const Pairs<CC_Pair> &u, const CC_vecfunction &singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S2c");
	CC_vecfunction t = make_t_intermediate(singles);
	for(auto i:t.functions){
		real_function_3d resulti = real_factory_3d(world);
		// The 6D Part
		real_function_3d resulti_6d = real_factory_3d(world);
		CC_Timer s2c_6d_timer(world,msg+"s2c-6D-part");
		{
			for(auto k:singles.functions){
				for(auto l:singles.functions){
					real_function_3d klgi = (mo_bra_[k.i]*intermediates_.get_EX(l.i,i.i)).truncate();
					real_function_3d lkgi = (mo_bra_[l.i]*intermediates_.get_EX(k.i,i.i)).truncate();
					real_function_3d tmp = u(k.i,l.i).function.project_out(klgi,1); // 1 means second particle
					real_function_3d tmpx= u(k.i,l.i).function.project_out(lkgi,1); // 1 means second particle
					resulti_6d += 2.0*tmp - tmpx;

					if(parameters.debug){
						CC_Timer s2c_6d_debug(world,msg+"s2c-6D-debug");
						real_function_6d test = multiply(copy(u(k.i,l.i).function),copy(intermediates_.get_EX(l.i,i.i)),2);
						real_function_3d test_tmp = test.project_out(mo_bra_[k.i],1);
						double test_norm = (test_tmp - tmp).norm2();
						if(world.rank()==0) std::cout << "s2c-6d-debug: ||method1 - method2 ||=" << test_norm << " ... should be 0"<< std::endl;
						if(test_norm > FunctionDefaults<6>::get_thresh()) warning("S2c-6D-Test-failed!",data);
						s2c_6d_debug.info();
					}
				}
			}
		}
		s2c_6d_timer.info();
		// The 3D Part (maybe merge this later with 6D part to avoid recalculating of  mo_bra_[k]*exchange_intermediate[i][l]
		real_function_3d resulti_3d = real_factory_3d(world);
		CC_Timer s2c_3d_timer(world,msg+"s2c_3D_part");
		{
			for(auto k:t.functions){
				for(auto l:t.functions){
					real_function_3d klgi = (mo_bra_[k.i]*intermediates_.get_EX(l.i,i.i)).truncate();
					real_function_3d lkgi = (mo_bra_[l.i]*intermediates_.get_EX(k.i,i.i)).truncate();
					real_function_3d k_lgi_Q12f12_tltk = convolute_x_gQf_yz(klgi,l.function,k.function);
					real_function_3d l_kgi_Q12f12_tltk = convolute_x_gQf_yz(lkgi,l.function,k.function);
					if(parameters.debug){
						// cant think of any inexpensive test right now
						CC_Timer s2c_debug_3d(world,msg+"S2c-3D-debug");
						real_function_6d test = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(k.function)).particle2(copy(l.function));
						test.fill_tree().truncate().reduce_rank();
						test = Q12(test);
						real_function_3d check = test.project_out(klgi,0); // 0 means particle 1
						real_function_3d diff = k_lgi_Q12f12_tltk - check;
						double diffn = diff.norm2();
						if(diffn > FunctionDefaults<6>::get_thresh()) warning("S2c-3D-Test failed with difference: " + stringify(diffn),data);
						else if(world.rank()==0) std::cout << "S2c-3D-Test passed: seems to be fine, difference is: " << diffn << std::endl;
						s2c_debug_3d.info();
					}
					resulti_3d += 2.0*k_lgi_Q12f12_tltk - l_kgi_Q12f12_tltk;
				}
			}
		}
		s2c_3d_timer.info();
		resulti = (resulti_6d + resulti_3d).truncate();
		result.push_back(resulti);
	}
	scale(world,result,-1.0);
	truncate(world,result);
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}


/// The S4a + X diagram
//[Q]       [i]
// \    ..../.....
//  \  /\  /     /\
//  _\/_ \/______\/_
/// -Q\sum (2<kl|g|\tau_il>|\tau_k> - <kl|g|\tau_ik>|\tau_l>)  : <kl|g|\tau_il>|\tau_k> = <k>
vecfuncT CC_Operators::S4a(const Pairs<CC_Pair> u, const CC_vecfunction & singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S4a");
	double prefactor = 1.0/(2.0*corrfac.gamma());
	for (auto i:singles.functions) {
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			for(auto l:singles.functions){
				// Coulomb Part of f12g12 = g12 - BSH
				{
					double a  =  make_ijgu(k.i,l.i,u(i.i,l.i));
					double ax =  make_ijgu(k.i,l.i,u(i.i,k.i));
					resulti -= (2.0 * a * k.function - ax * l.function);
					resulti.truncate();
				}
				// BSH part of f12g12
				{
					double b  = make_ijGu(k.i,l.i,u(i.i,l.i));
					double bx = make_ijGu(k.i,l.i,u(i.i,k.i));
					resulti -= (2.0 * b * k.function- bx * l.function);
					resulti.truncate();
				}
			}
		}
		resulti.scale(prefactor);
		resulti.print_size("S4a_"+stringify(i.i)+" potential");
		result.push_back(resulti);
	}
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S4b -> Merge this with S2c later to save time
//[i]       [Q]
// \    ..../.....
//  \  /\  /     /\
//  _\/_ \/______\/_
/// -Q\sum_{kl} (2<k(3)l(4)|g34f14|\tau_{i}(3)u_{kl}(1,4)>  // exchange part - <k(4)l(3)|g34f14|\tau_i(3)u_{lk}(1,4)>)
/// =-Q\sum_{kl} (2<k(3)l(4)|g34f14|\tau_{i}(3)u_{kl}(1,4)>  // exchange part - <k(4)l(3)|g34f14|\tau_i(3)u_{kl}(4,1)>)
// current procedure
// make intermediates: <k|\tau_i> and <l|g|\tau_i> -> get the perturbed_exchange_intermediate
// process intermedaites: lkgi = l*<k|g|tau_i> amd klgi = k*<l|g|\tau_i> (k and l are from the bra space)
// Project out:  <lkgi|tkl>_2 and <klgi|tkl>_1
// return -2.0*Q(<lkgi|tkl>_2) + Q(<klgi|tkl>_1)
vecfuncT CC_Operators::S4b(const Pairs<CC_Pair> u, const CC_vecfunction & singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S4b");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){
			for(auto l:singles.functions){ // index l is also pair function index therefore it runs only over non frozen orbs
				// make the cuspy pair
				CC_Timer make_tau(world,msg+"Creating cuspy pair: f12*u(1,2)");
				real_function_6d tkl =  CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).ket(copy(u(k.i,l.i).function));
				tkl.fill_tree().truncate().reduce_rank();
				make_tau.info();
				tkl.print_size("f12*u(12)");

				// make intermediates
				real_function_3d kgi = intermediates_.get_pEX(k.i,i.i);//(*poisson)(mo_bra_[k]*tau[i]);
				real_function_3d lgi = intermediates_.get_pEX(l.i,i.i);//(*poisson)(mo_bra_[k]*tau[i]);
				real_function_3d lkgi = (mo_bra_[l.i]*kgi).truncate();
				real_function_3d klgi = (mo_bra_[k.i]*lgi).truncate();

				// project out (1 is particle 2 and 0 is particle 1)
				real_function_3d lkgitkl = tkl.project_out(lkgi,1);
				real_function_3d klgitkl = tkl.project_out(klgi,0);

				resulti += (-2.0*lkgitkl+klgitkl).truncate();
				Q(resulti);
			}
		}
		resulti.print_size("S4b_"+stringify(i.i)+" potential");
		result.push_back(resulti);
	}
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// The S4c + X + X + X + X Diagrams -> merge this with S2c later to save time
//            [i]   [Q]
//   .......   \    /
//  /\     /\   \  /
// _\/_   _\/____\/_
/// Q\sum_{kl}[ 4*<k(3)l(4)|g34 f14| \tau_k(3) u_{il}(1,4)> - 2* <k(3)l(4)|g34 f14|\tau_k(4) u_{il}(1,3)>
/// -           2*<k(3)l(4)|g34 f14| \tau_k(3) u_{il}(4,1)>	    +<k(3)l(4)|g34 f14|\tau_k(4) u_{il}(3,1)>  ]

// current procedure:
// make cuspy function: til = uil*f12
// make intermediates: lkgk= l*<k|g|\tau_k> and klgk= k*<l|g|\tau_k>
// Project out:
// 1st term: <lkgk|til>_2
// 2nd term: <klgk|til>_2
// 3rd term: <lkgk|til>_1
// 4th term: <klgk|til>_1

vecfuncT CC_Operators::S4c(const Pairs<CC_Pair> u, const CC_vecfunction & singles, const std::string &msg) const {
	vecfuncT result;
	CC_Timer timer(world,"time");
	CC_data data("S4c");
	for(auto i:singles.functions){
		real_function_3d resulti = real_factory_3d(world);
		for(auto k:singles.functions){	// maybe loop over all pairs
			for(auto l:singles.functions){
				// make the cuspy pair
				CC_Timer make_tau(world,msg+"Creating cuspy pair: f12*u(1,2)");
				real_function_6d til =  CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).ket(copy(u(i.i,l.i).function));
				til.fill_tree().truncate().reduce_rank();
				make_tau.info();
				til.print_size("f12*u(12)");

				// get intermediates
				real_function_3d lkgk = (mo_bra_[l.i]*intermediates_.get_pEX(k.i,k.i)).truncate();
				real_function_3d klgk = (mo_bra_[l.i]*intermediates_.get_pEX(k.i,l.i)).truncate();

				// project out
				real_function_3d term1 = til.project_out(lkgk,1);
				real_function_3d term2 = til.project_out(klgk,1);
				real_function_3d term3 = til.project_out(lkgk,0);
				real_function_3d term4 = til.project_out(klgk,0);

				resulti += (4.0*term1 -2.0*term2 - 2.0*term3 + term4).truncate();
				Q(resulti);
			}
		}
		resulti.print_size("S4c_"+stringify(i.i)+" potential");
		result.push_back(resulti);
	}
	Q(result);
	data.result_size=get_size(world,result);
	data.result_norm=norm2(world,result);
	data.time = timer.current_time();
	performance_S.insert(data.name,data);
	return result;
}

/// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
/// @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
/// @param[in] \tau_j
/// @param[in] u, the uij pair structure which holds the consant part of MP2
/// @param[out] Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
/// Right now Calculated in the decomposed form: |titj> = |i,j> + |\taui,\tauj> + |i,\tauj> + |\taui,j>
/// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
real_function_6d CC_Operators::make_cc2_residue(const CC_function &taui, const CC_function &tauj, const CC_Pair &u)const{
	output("Now doing the CC2-Regularization-Residue");
	// make intermediates
	MADNESS_ASSERT(u.i == taui.i);
	MADNESS_ASSERT(u.j == tauj.i);
	const size_t i=taui.i;
	const size_t j=tauj.i;

	// DEBUG
	if(parameters.debug){
		if(world.rank()==0) std::cout << "FOCK OPERATOR CONSISTENCY CHECK\n";
		real_function_3d Fi = apply_F(CC_function(mo_ket_.front(),0,HOLE));
		real_function_3d ei = get_orbital_energies()[0]*mo_ket_[0];
		real_function_3d diff = Fi - ei;
		double ndiff = diff.norm2();
		if(ndiff < FunctionDefaults<3>::get_thresh()) std::cout << "... Passed";
		else std::cout << "... Failed";
		std::cout << " ||F|i>-ei|i>||=" << ndiff << std::endl;
	}
	// DEBUG END

	// f12 parts first
	real_function_6d GQf_parts = real_function_6d(world);
	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);

	// this intermediates we need for all the blocks
	CC_Timer Ftaui_time(world,"(F-ei)|\tau"+stringify(taui.i)+">");
	real_function_3d Ftaui = apply_F(taui);
	Ftaui -= get_orbital_energies()[i]*Ftaui;
	Ftaui_time.info();

	if(parameters.debug){
		plot(taui.function,"tau"+stringify(i));
		plot(Ftaui,"Ftaui");
	}

	CC_Timer Ftauj_time(world,"(F-ej)|\tau"+stringify(tauj.i)+">");
	real_function_3d Ftauj = apply_F(tauj);
	Ftauj -= get_orbital_energies()[j]*Ftauj;
	Ftauj_time.info();
	real_function_6d part1,part2,part3;
	{

		// Qf(F -eij)|ij> is zero

		// Qf(F -eij)|\taui,j> = Qf(F1 - ei)|\taui,j> = Qf( |(F1-ei)\taui>(x)|j>
		output("Qf(F1-ei)|taui,j>");
		CC_Timer GQfFtaui_time(world,"GQf(F-ei)|\tau"+stringify(taui.i)+">");
		real_function_6d Qfxj = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Ftaui)).particle2(copy(mo_ket_[j]));
		CC_Timer fill_tree_1(world,"fill_tree(G):f|Fx,j>");
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified() = true;
		Qfxj.fill_tree(screenG).truncate().reduce_rank();
		fill_tree_1.info();
		apply_Q12(Qfxj,"f|Fx,j>");
		real_function_6d tmp = G(Qfxj);
		GQf_parts += tmp;
		part1 = tmp;
		GQfFtaui_time.info();
	}{
		// Qf(F-eij)|i\tauj> = Qf(F2-ej)|i\tauj>
		output("Qf(F2-ej)|i,tauj>");
		CC_Timer GQfFtauj_time(world,"GQf(F-ej)|\tau"+stringify(tauj.i)+">");
		real_function_6d Qfix = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(mo_ket_[i])).particle2(copy(Ftauj));
		CC_Timer fill_tree_2(world,"fill_tree(G):f|i,Fx>");
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified() = true;
		Qfix.fill_tree(screenG).truncate().reduce_rank();
		fill_tree_2.info();
		apply_Q12(Qfix,"f|i,Fx>");
		real_function_6d tmp = G(Qfix);
		GQf_parts += tmp;
		part2 = tmp;
		GQfFtauj_time.info();
	}{
		// Qf(F-eij)|\taui\tauj> = Qf(|(F1-ei)\taui>(x)|(F2-ej)|\tauj>)
		output("Qf(F-eij)|taui,tauj>");
		CC_Timer GQfFxFy_time(world,"GQf(F-eij)|\tau"+stringify(taui.i)+"\tau" + stringify(j) + ">");
		real_function_6d Qfxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Ftaui)).particle2(copy(Ftauj));
		CC_Timer fill_tree_3(world,"fill_tree(G):f|Fx,Fy>");
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified()=true;
		Qfxy.fill_tree(screenG).truncate().reduce_rank();
		fill_tree_3.info();
		apply_Q12(Qfxy,"f|Fx,Fy>");
		real_function_6d tmp =  G(Qfxy);
		GQf_parts += tmp;
		part3 = tmp;
		GQfFxFy_time.info();
	}


	if(world.rank()==0)std::cout << "\n\nFinished with GQf(F-eij)|xy> Parts:\n";
	part1.print_size("GQf(F1-ei)|taui,j>   :");
	part2.print_size("GQf(F2-ej)|i,tauj>   :");
	part3.print_size("GQf(F-eij)|taui,tauj>:");
	GQf_parts.print_size("GQf(F-eij)|xy>       :");
	if(world.rank()==0) std::cout << "\n\n" << std::endl;


	// Ue and KffK Parts
	// pack some CC_3D_functions which hold all necessary informations
	CC_function moi(mo_ket_[i],i,HOLE);
	CC_function moj(mo_ket_[j],j,HOLE);
	taui.function.print_size("original taui");
	tauj.function.print_size("original tauj");
	// |ij> part is already there (constant part of MP2)

	// |taui,j> part
	real_function_6d taui_j_part = real_factory_6d(world);
	{
		output("QUe|taui,j>");
		CC_Timer Ue_time(world,"GQUe|taui,j>");
		real_function_6d Uef = apply_transformed_Ue(taui.function,mo_ket_[j],i,j);
		//apply_Q12(QUef,"Ue|taui,j>");
		//GQ_Ue_taui_j = G(QUef);
		Ue_time.info();
		output("Q[K,f]|taui,j>");
		CC_Timer KffK_time(world,"GQ[K,f]|taui,j>");
		real_function_6d KffK_tmp = apply_exchange_commutator(taui,moj);
		//apply_Q12(KffK_tmp,"[K,f]|taui,j>");
		KffK_time.info();
		taui_j_part = Uef - KffK_tmp;
	}
	// |i,tauj> part
	real_function_6d i_tauj_part = real_factory_6d(world);
	if(i==j){
		if(world.rank()==0) std::cout << "Exploit Symmetry: P(1,2)|taui,j> = |j,taui> for i=j case\n";
		CC_Timer symex_timer(world,"GQUe|i,tauj> (i=j, symmetry exploit)");
		i_tauj_part = swap_particles(taui_j_part);
		symex_timer.info();
	}
	if(i!=j or parameters.debug){
		output("QUe|i,tauj>");
		CC_Timer Ue_time(world,"GQUe|i,tauj>");
		real_function_6d Uef = apply_transformed_Ue(mo_ket_[i],tauj.function,i,j);
		Ue_time.info();

		CC_Timer KffK_time(world,"GQ[K,f]|taui,j>");
		real_function_6d KffK_tmp = apply_exchange_commutator(moi,tauj);
		KffK_time.info();

		if(i!=j)i_tauj_part = Uef - KffK_tmp;
		else if(parameters.debug){
			real_function_6d diff = (Uef - KffK_tmp)-i_tauj_part;
			diff.print_size("DEBUG: U-KffK difference, calc and swap");
			if(diff.norm2()>FunctionDefaults<6>::get_thresh()) warning("ERROR in comparison between swapped and calced particles");
		}
	}

	// |taui,tauj> part
	// ADD SCREENING HERE (for small tau) or/and flexible thresh
	real_function_6d taui_tauj_part = real_factory_6d(world);
	{
		output("Ue|taui,tauj>");
		CC_Timer Ue_time(world,"QUe|taui,tauj>");
		real_function_6d Uef = apply_transformed_Ue(taui.function,tauj.function,i,j);
		Ue_time.info();
		output("[K,f]|taui,tauj>");
		CC_Timer KffK_time(world,"[K,f]|taui,tauj>");
		real_function_6d KffK_tmp = apply_exchange_commutator(taui,tauj);
		KffK_time.info();

		taui_tauj_part = Uef - KffK_tmp;
	}

	real_function_6d UeK_parts = (taui_j_part + i_tauj_part + taui_tauj_part);
	apply_Q12(UeK_parts,"UeK_Parts");
	real_function_6d G_UeK_parts = apply_G(UeK_parts,i,j);



	// make the result;
	real_function_6d result = GQf_parts + G_UeK_parts;
	result.truncate();

	output("Finished with QUe|xy> and Q[K,f]|xy> Parts:");
	taui_tauj_part.print_size("Q12(Ue-[K,f])|taui,tauj>");
	i_tauj_part.print_size("Q12(Ue-[K,f])|i,tauj>   ");
	taui_j_part.print_size("Q12(Ue-[K,f])|taui,j>   ");
	G_UeK_parts.print_size("G(UeKffKparts)");
	GQf_parts.print_size("GQf(F-e)|titj>_parts");
    output("\n\n");


	// this is expensive !!!!
	if(parameters.debug){
		output("Expensive Debug Calculation");

		const CC_function ti(taui.function + mo_ket_[i],i,MIXED);
		const CC_function tj(tauj.function + mo_ket_[j],j,MIXED);

		real_function_6d GQUK_test = (-2.0)*(G_UeK_parts) + u.constant_term;
		 GQUK_test.print_size("Ue and [K,f] + mp2_residue");
		apply_Q12(GQUK_test,"GQUK_test");
		 GQUK_test.print_size("Q12(Ue and [K,f] + mp2_residue)");

		CC_Timer Ue_time(world,"GQUe|titj>");
		real_function_6d QUe_titj = apply_transformed_Ue(ti.function,tj.function,i,j);
		apply_Q12(QUe_titj,"Uetitj");
		real_function_6d GQUe_titj = G(QUe_titj);
		Ue_time.info();

		CC_Timer KffK_time(world,"GQ[K,f]|titj>");
		real_function_6d QKffK_titj = apply_exchange_commutator(ti,tj);
		apply_Q12(QKffK_titj,"Uetitj");
		real_function_6d GQKffK_titj = G(QKffK_titj);
		KffK_time.info();

		real_function_6d test = (GQUe_titj - GQKffK_titj).truncate();
		test.scale(-2.0);
		apply_Q12(test,"test");
		real_function_6d diff = test - GQUK_test;
		diff.print_size("difference between sepparated + constantMP2 and titj interemdiates");
		double norm = diff.norm2();
		if(norm>FunctionDefaults<6>::get_thresh()) warning("ERROR in CC2 residue at the end, different results for ti intermediate calculations");
		else output("\n\n |titj> and decomposed+constant_part methods are the same!\n\n");
		output("Debug Calculation ended");
	}

	output("CC2-Regularization-Residue finished");
	result.print_size("result");
	return result;
}


// apply the kinetic energy operator with cusp to a decomposed 6D function
/// @param[in] a 3d function x (will be particle 1 in the decomposed 6d function)
/// @param[in] a 3d function y (will be particle 2 in the decomposed 6d function)
/// @param[out] a 6d function: G(f12*T*|xy>)
real_function_6d CC_Operators::make_GQfT_xy(const real_function_3d &x, const real_function_3d &y, const size_t &i, const size_t &j)const{
	// construct the greens operator
	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);

	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	vecfuncT gradx,grady;
	vecfuncT laplacex, laplacey;
	for(size_t axis=0;axis<3;axis++){
		real_function_3d gradxi = (*gradop[axis])(x);
		real_function_3d gradyi = (*gradop[axis])(y);
		gradx.push_back(gradxi);
		grady.push_back(gradyi);
		real_function_3d grad2xi = (*gradop[axis])(gradxi);
		real_function_3d grad2yi = (*gradop[axis])(gradyi);
		laplacex.push_back(grad2xi);
		laplacey.push_back(grad2yi);
	}
	real_function_3d laplace_x = laplacex[0]+laplacex[1]+laplacex[2];
	real_function_3d laplace_y = laplacey[0]+laplacey[1]+laplacey[2];
	real_function_3d Tx = laplace_x.scale(-0.5);
	real_function_3d Ty = laplace_y.scale(-0.5);
	// make the two screened 6D functions
	// fTxy = f12 |(\Delta x)y> , fxTy = f12 |x\Delta y> (delta = laplace_operator)
	real_function_6d fTxy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(Tx)).particle2(copy(y));
	real_function_6d fxTy = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(x)).particle2(copy(Ty));
	// for now construct explicitly and project out Q12 later: use BSH operator to screen
	if(world.rank()==0) std::cout << "Constructing fTxy with G as screening operator\n";
	CC_Timer fTxy_construction_time(world,"Screened 6D construction of fTxy");
	{
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified()=true;
		fTxy.fill_tree(screenG).truncate().reduce_rank();
	}{
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified()=true;
		fxTy.fill_tree(screenG).truncate().reduce_rank();
	}
	fTxy_construction_time.info();
	CC_Timer addition_time(world,"f(Tx)y + fxTy");
	real_function_6d result = (fTxy + fxTy).truncate();
	apply_Q12(result,"fT|xy>");

	if(parameters.debug){
		CC_Timer lapladebugtime(world,"Laplace-Debug");
		if(world.rank()==0) std::cout << "Debugging Laplace Operator:\n";
		double xd2x_0 = x.inner(laplacex[0]);
		double xd2x_1 = x.inner(laplacex[1]);
		double xd2x_2 = x.inner(laplacex[2]);
		double dxdx_0 = gradx[0].inner(gradx[0]);
		double dxdx_1 = gradx[1].inner(gradx[1]);
		double dxdx_2 = gradx[2].inner(gradx[2]);
		double diff0 = xd2x_0 + dxdx_0; // + because <f|f2x|f> = - <dxf|dxf>
		double diff1 = xd2x_1 + dxdx_1;
		double diff2 = xd2x_2 + dxdx_2;
		double thresh = FunctionDefaults<3>::get_thresh();
		if(diff0 >thresh or diff1>thresh or diff2>thresh) warning("Laplace operator was not precise");
		if(world.rank()==0){
			std::cout << "<f|d2x|f>, -<dxf|dxf>, diff\n";
			std::cout << xd2x_0 << ", " << dxdx_0 << ", " << diff0 << std::endl;
			std::cout << xd2x_1 << ", " << dxdx_1 << ", " << diff1 << std::endl;
			std::cout << xd2x_2 << ", " << dxdx_2 << ", " << diff2 << std::endl;
		}
		lapladebugtime.info();
	}

	CC_Timer apply_G(world,"G(fTxy)");
	real_function_6d G_result = G(result);
	G_result.truncate();
	apply_G.info();
	return G_result;
}



/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
real_function_6d CC_Operators::fock_residue_6d(const CC_Pair &u) const {
	const double eps = get_epsilon(u.i, u.j);
	// make the coulomb and local Un part with the composite factory
	real_function_3d local_part = (2.0
			* intermediates_.get_hartree_potential()
			+ nemo.nuclear_correlation->U2());
	local_part.print_size("vlocal");
	u.function.print_size("u");

	// Contruct the BSH operator in order to screen

	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps),
			parameters.lo, parameters.thresh_bsh_6D);
	// apparently the modified_NS form is necessary for the screening procedure
	op_mod.modified() = true;
	// Make the CompositeFactory
	real_function_6d vphi =
			CompositeFactory<double, 6, 3>(world).ket(copy(u.function)).V_for_particle1(
					copy(local_part)).V_for_particle2(copy(local_part));
	// Screening procedure
	vphi.fill_tree(op_mod);

	vphi.print_size("vlocal|u>");

	// the part with the derivative operators: U1
	for (int axis = 0; axis < 6; ++axis) {
		real_derivative_6d D = free_space_derivative<double, 6>(world,
				axis);
		// Partial derivative of the pari function
		const real_function_6d Du = D(u.function).truncate();

		// % operator gives division rest (modulo operator)
		if (world.rank() == 0)
			print("axis, axis^%3, axis/3+1", axis, axis % 3, axis / 3 + 1);
		const real_function_3d U1_axis = nemo.nuclear_correlation->U1(
				axis % 3);

		double tight_thresh = parameters.thresh_Ue;
		real_function_6d x;
		if (axis / 3 + 1 == 1) {
			x =
					CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle1(
							copy(U1_axis)).thresh(tight_thresh);

		} else if (axis / 3 + 1 == 2) {
			x =
					CompositeFactory<double, 6, 3>(world).ket(copy(Du)).V_for_particle2(
							copy(U1_axis)).thresh(tight_thresh);
		}
		x.fill_tree(op_mod);
		x.set_thresh(FunctionDefaults<6>::get_thresh());
		vphi += x;
		vphi.truncate().reduce_rank();
	}

	vphi.print_size("(Un + J1 + J2)|u>");

	// Exchange Part
	vphi = (vphi - K(u.function, u.i == u.j)).truncate().reduce_rank();
	vphi.print_size("(Un + J1 + J2 - K1 - K2)|U>");
	vphi.truncate();
	vphi.print_size("truncated: (Un + J1 + J2 - K1 - K2)|U>");
	return vphi;

}

/// Echange Operator on 3D function
/// !!!!Prefactor (-1) is not included
real_function_3d CC_Operators::K(const CC_function &x)const{
	const real_function_3d &f = x.function;
	real_function_3d result = apply_K(x);

	// Sanity Check (expensive when not helium)
	if(mo_ket_.size()<3){
		real_function_3d result2 = real_factory_3d(world);
		// multiply rhs with R2orbitals (the bra space)
		vecfuncT R2rhs = mul(world, f, mo_bra_);
		for (std::size_t k = 0; k < mo_ket_.size(); ++k) {
			result2 += mo_ket_[k] * (*poisson)(R2rhs[k]);
		}
		double sanity = (result-result2).norm2();
		if(sanity < FunctionDefaults<3>::get_thresh()){
			if(world.rank()==0) std::cout << "Sanity Check of K passed\n";
		}else{
			if(world.rank()==0) std::cout << "Sanity Check of K NOT passed (" << sanity <<  ")\n";
		}
	}

	return result;
}

/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
/// if i==j in uij then the symmetry will be exploited
/// !!!!Prefactor (-1) is not included here!!!!
real_function_6d CC_Operators::K(const real_function_6d &u,
		const bool symmetric) const {
	/// DEBUG
	if(world.rank()==0)std::cout << "Entering K" << std::endl;
	/// DEBUG END

	/// TEST IF THIS WILL WORK FOR THE += Operator
	real_function_6d result = real_factory_6d(world).compressed();
	// K(1) Part
	result += apply_K(u, 1);
	// K(2) Part
	if (symmetric)
		result += swap_particles(result);
	else
		result += apply_K(u, 2);

	return (result.truncate());
}

/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
/// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
/// 1. X(3,2) = bra_k(3)*u(3,2)
/// 2. Y(1,2) = \int X(3,2) g13 d3
/// 3. result = Y(1,2)*ket_k(1)
/// !!!!Prefactor (-1) is not included here!!!!
real_function_6d CC_Operators::apply_K(const real_function_6d &u,
		const size_t &particle) const {
	/// DEBUG
	if(world.rank()==0)std::cout << "Entering apply_K" << std::endl;
	/// DEBUG END
	MADNESS_ASSERT(particle == 1 or particle == 2);
	poisson->particle() = particle;
	/// WARNING: CHECK IF THIS WORKS -> bc of the += operator later
	real_function_6d result = real_factory_6d(world).compressed();
	for (size_t k = 0; k < mo_ket_.size(); k++) {
		real_function_6d X = (multiply(copy(u), copy(mo_bra_[k]), particle)).truncate();
		real_function_6d Y = (*poisson)(X);
		result += multiply(copy(Y), copy(mo_ket_[k]), particle).truncate();
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
real_function_6d CC_Operators::apply_transformed_Ue(const real_function_3d x,
		const real_function_3d y, const size_t &i, const size_t &j) const {
	real_function_6d Uxy = real_factory_6d(world);
	// Apply the untransformed U Potential
	const double eps = get_epsilon(i, j);
	Uxy = corrfac.apply_U(x, y, eps);

	// Get the 6D BSH operator in modified-NS form for screening
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps),
			parameters.lo,
			parameters.thresh_Ue);
	op_mod.modified() = true;

	// make shure the thresh is high enough
	double tight_thresh = parameters.thresh_Ue;

	// Apply the double commutator R^{-1}[[T,f,R]
	for (size_t axis = 0; axis < 3; axis++) {
		// Make the local parts of the Nuclear and electronic U potentials
		const real_function_3d Un_local = nemo.nuclear_correlation->U1(
				axis);
		const real_function_3d Un_local_x = (Un_local * x).truncate();
		const real_function_3d Un_local_y = (Un_local * y).truncate();
		const real_function_6d Ue_local = corrfac.U1(axis);
		// Now add the Un_local_x part to the first particle of the Ue_local potential
		real_function_6d UeUnx = CompositeFactory<double, 6, 3>(world).g12(
				Ue_local).particle1(Un_local_x).particle2(copy(y)).thresh(
						tight_thresh);
		// Fill the Tree were it will be necessary
		UeUnx.fill_tree(op_mod);
		// Set back the thresh
		UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());

		UeUnx.print_size("UeUnx");

		// Now add the Un_local_y part to the second particle of the Ue_local potential
		real_function_6d UeUny = CompositeFactory<double, 6, 3>(world).g12(
				Ue_local).particle1(copy(x)).particle2(Un_local_y).thresh(
						tight_thresh);
		// Fill the Tree were it will be necessary
		UeUny.fill_tree(op_mod);
		// Set back the thresh
		UeUny.set_thresh(FunctionDefaults<6>::get_thresh());

		UeUny.print_size("UeUny");

		// Construct the double commutator part and add it to the Ue part
		real_function_6d diff = (UeUnx - UeUny).scale(-1.0);
		diff.truncate();
		Uxy = (Uxy+diff).truncate();
	}

	// sanity check: <xy|R2 [T,g12] |xy> = <xy |R2 U |xy> - <xy|R2 g12 | xy> = 0
	real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
			copy(x*nemo.nuclear_correlation -> square())).particle2(copy(y*nemo.nuclear_correlation -> square()));

	const double a = inner(Uxy, tmp);
	const real_function_3d xx = (x * x*nemo.nuclear_correlation -> square());
	const real_function_3d yy = (y * y*nemo.nuclear_correlation -> square());
	const real_function_3d gxx = (*poisson)(xx);
	const double aa = inner(yy, gxx);
	const double error = std::fabs(a - aa);
	if (world.rank() == 0) {
		printf("<xy| U_R |xy>  %12.8f\n", a);
		printf("<xy|1/r12|xy>  %12.8f\n", aa);
		if (error > FunctionDefaults<6>::get_thresh())
			print("WARNING : Kutzelnigg's potential inaccurate (box size, thresh ?)");
		//if (error > FunctionDefaults<6>::get_thresh() * 10.0)
		//	MADNESS_EXCEPTION("Kutzelnigg's potential plain wrong (box size, thresh ?)", 1);
	}
	Uxy.print_size("Uphi0");

	return Uxy;
}

/// Apply the Exchange Commutator [K,f]|xy>
real_function_6d CC_Operators::apply_exchange_commutator(const CC_function &x, const CC_function &y)const{
	// make first part of commutator
	CC_Timer part1_time(world,"Kf");
	real_function_6d Kfxy = apply_Kf(x,y).truncate();
	part1_time.info();
	// make the second part of the commutator
	CC_Timer part2_time(world,"fK");
	real_function_6d fKxy = apply_fK(x,y).truncate();
	part2_time.info();
	real_function_6d result = (Kfxy - fKxy);


	// for sanity check:
	CC_Timer sanity_time(world,"[K,f]|xy> sanity check");
	if(parameters.debug){
		double expv_first_part = 0.0;
		double expv_second_part = 0.0;

		{
			real_function_3d bra_x, bra_y;
			if(x.type == HOLE) bra_x = mo_bra_[x.i];
			else bra_x = x.function*nemo.nuclear_correlation -> square();
			if(y.type == HOLE) bra_y = mo_bra_[y.i];
			else bra_y = y.function*nemo.nuclear_correlation -> square();
			{
				real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
						copy(bra_x)).particle2(copy(bra_y));
				expv_first_part = inner(Kfxy,tmp);
			}{
				real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
						copy(bra_x)).particle2(copy(bra_y));
				expv_second_part = inner(fKxy,tmp);
			}
			if(world.rank()==0){
				std::cout << "Apply [K,f]|x,y> sanity check:";
				std::cout <<  "\n<xy|Kf|xy> =" << expv_first_part;
				std::cout <<  "\n<xy|fK|xy> =" << expv_second_part <<std::endl;
			}
			{
				// sanity check: The Expectation value of the Kommutator must vanish (symmetry)
				// <xy|[A,B]|xy> = <xy|AB|xy> - <xy|BA|xy> since A=f and B=K -> both are hermitian and symmetric under particle exchange
				real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
						copy(bra_x)).particle2(copy(bra_y));
				const double a = inner(result, tmp);
				if (world.rank() == 0) {
					printf("<xy| R^2 R-1 [K,f] R |xy >  %12.8f\n", a);
					if (std::fabs(a) > FunctionDefaults<6>::get_thresh())
						warning("Exchange commutator inaccurate");
					if (std::fabs(a) > FunctionDefaults<6>::get_thresh() * 10.0)
						if(parameters.debug) warning("exchange commutator plain wrong");
						else{
							MADNESS_EXCEPTION("exchange commutator plain wrong", 1);
						}
				}
			}
		}

	}
	sanity_time.info();




	return result;
}

/// Apply the Exchange operator on a tensor product multiplied with f12
/// !!! Prefactor of (-1) is not inclued in K here !!!!
real_function_6d CC_Operators::apply_Kf(const CC_function &x, const CC_function &y) const {
	bool symmetric = false;
	if((x.type == y.type) and (x.i == y.i)) symmetric = true;

	CC_Timer timer_f12xy(world,"Constructed f12|xy>");
	// First make the 6D function f12|x,y>
	real_function_6d f12xy = CompositeFactory<double, 6, 3>(world).g12(
			corrfac.f()).particle1(copy(x.function)).particle2(copy(y.function));
	f12xy.fill_tree().truncate().reduce_rank();
	timer_f12xy.info();
	// Apply the Exchange Operator
	real_function_6d result = K(f12xy, symmetric);
	return result.truncate();
}

/// Apply fK on a tensor product of two 3D functions
/// fK|xy> = fK_1|xy> + fK_2|xy>
/// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
/// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
real_function_6d CC_Operators::apply_fK(const CC_function &x, const CC_function &y) const {
	const real_function_3d& phi_i = x.function;
	const real_function_3d& phi_j = y.function;

	const real_function_3d Kphi_i = K(x.function);
	const real_function_3d Kphi_j = K(y.function);

	real_function_6d fKphi0a = CompositeFactory<double, 6, 3>(world).g12(
			corrfac.f()).particle1(copy(phi_i)).particle2(copy(Kphi_j));
	fKphi0a.fill_tree().truncate();
	real_function_6d fKphi0b = CompositeFactory<double, 6, 3>(world).g12(
			corrfac.f()).particle1(copy(Kphi_i)).particle2(copy(phi_j));
	fKphi0b.fill_tree().truncate();

	real_function_6d fKphi0 = (fKphi0a + fKphi0b).truncate();
	return fKphi0;

}

/// swap particles 1 and 2

/// param[in]	f	a function of 2 particles f(1,2)
/// return	the input function with particles swapped g(1,2) = f(2,1)
real_function_6d CC_Operators::swap_particles(const real_function_6d& f) const {
	CC_Timer timer_swap(world,"swap particles");
	// this could be done more efficiently for SVD, but it works decently
	std::vector<long> map(6);
	map[0] = 3;
	map[1] = 4;
	map[2] = 5;	// 2 -> 1
	map[3] = 0;
	map[4] = 1;
	map[5] = 2;	// 1 -> 2
	timer_swap.info();
	return mapdim(f, map);
}

// Calculate the CC2 energy equation which is
// \omega = \sum_{ij} 2<ij|g|\tau_{ij}> - <ij|g|\tau_{ji}> + 2 <ij|g|\tau_i\tau_j> - <ij|g|\tau_j\tau_i>
// with \tau_{ij} = u_{ij} + Q12f12|ij> + Q12f12|\tau_i,j> + Q12f12|i,\tau_j> + Q12f12|\tau_i\tau_j>
double CC_Operators::get_CC2_correlation_energy() const {
	MADNESS_EXCEPTION("get_cc2_correlation_energy not implemented yet",1);
	return 0.0;
}

double CC_Operators::compute_ccs_correlation_energy(const CC_function &taui, const CC_function &tauj)const{
	double omega = 0.0;
	//omega += 2.0*intermediates_.get_integrals_t1()(u.i,u.j,u.i,u.j); //<ij|g|\taui\tauj>
	omega += 2.0*make_ijgxy(taui.i,tauj.i,taui.function,tauj.function);
	//omega -= intermediates_.get_integrals_t1()(u.i,u.j,u.j,u.i);     //<ij|g|\tauj\taui>
	omega -= make_ijgxy(taui.i,tauj.i,tauj.function,taui.function);
	return omega;
}
double CC_Operators::compute_cc2_pair_energy(const CC_Pair &u,
		const real_function_3d &taui, const real_function_3d &tauj) const {
	double omega = 0.0;
	const size_t i = u.i;
	const size_t j = u.j;
	// Contribution from u itself, we will calculate <uij|g|ij> instead of <ij|g|uij> and then just make the inner product (see also mp2.cc)
	{
		real_function_6d coulomb = TwoElectronFactory(world).dcut(
				FunctionDefaults<6>::get_thresh());
		real_function_6d g_ij =
				CompositeFactory<double, 6, 3>(world).particle1(
						copy(mo_bra_[i])).particle2(copy(mo_bra_[j])).g12(
								coulomb);
		real_function_6d g_ji =
				CompositeFactory<double, 6, 3>(world).particle1(
						copy(mo_bra_[j])).particle2(copy(mo_bra_[i])).g12(
								coulomb);
		const double uij_g_ij = inner(u.function, g_ij);
		const double uij_g_ji = inner(u.function, g_ji); // =uji_g_ij
		omega += 2.0 * uij_g_ij - uij_g_ji;
	}
	// Contribution from the mixed f12(|\tau_i,j>+|i,\tau_j>) part
	{
		omega += 2.0*make_ijgQfxy(u.i,u.j,mo_ket_[i],tauj);
		omega += 2.0*make_ijgQfxy(u.i,u.j,taui,mo_ket_[j]);
		omega -= make_ijgQfxy(u.j,u.i,mo_ket_[i],tauj);
		omega -= make_ijgQfxy(u.j,u.i,taui,mo_ket_[j]);
	}
	// Contribution from the f12|ij> part, this should be calculated in the beginning
	{
		omega += (2.0*u.ij_gQf_ij - u.ji_gQf_ij );
	}
	// Contribution from the f12|\tau_i\tau_j> part
	{
		omega += 2.0*make_ijgQfxy(u.i,u.j,taui,tauj);
		omega -= make_ijgQfxy(u.i,u.j,tauj,taui);
	}
	// Singles Contribution
	{
		// I should use intermediates later because the t1 integrals are also needed for the CC2 potential
		//omega += 2.0*intermediates_.get_integrals_t1()(u.i,u.j,u.i,u.j); //<ij|g|\taui\tauj>
		omega += 2.0*make_ijgxy(u.i,u.j,taui,tauj);
		//omega -= intermediates_.get_integrals_t1()(u.i,u.j,u.j,u.i);     //<ij|g|\tauj\taui>
		omega -= make_ijgxy(u.i,u.j,tauj,taui);
	}
	return omega;
}

/// General Function to make the intergral <ij|gQf|xy>
double CC_Operators::make_ijgQfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const{
	// Q12 = I12 - O1 - O2 + O12
	real_function_3d jy = mo_bra_[j]*y;
	real_function_3d ix = mo_bra_[i]*x;
	// I12 Part:
	double ijgfxy = (ix).inner(apply_gf(jy));
	// O1 Part
	double ijgO1fxy =0.0;
	for(size_t k=0;k<mo_ket_.size();k++){
		real_function_3d igk = intermediates_.get_EX(i,k);
		real_function_3d kfx = (*f12op)(mo_bra_[k]*x);
		real_function_3d igkkfx = (igk*kfx).truncate();
		ijgO1fxy += jy.inner(igkkfx);
	}
	// O2 Part
	double ijgO2fxy =0.0;
	for(size_t k=0;k<mo_ket_.size();k++){
		real_function_3d jgk = intermediates_.get_EX(j,k);
		real_function_3d kfy = (*f12op)(mo_bra_[k]*y);
		real_function_3d jgkkfy = (jgk*kfy).truncate();
		ijgO2fxy += ix.inner(jgkkfy);
	}
	// O12 Part
	double ijgO12fxy = 0.0;
	for(size_t k=0;k<mo_ket_.size();k++){
		real_function_3d igk = intermediates_.get_EX(i,k);
		real_function_3d kfx = (*f12op)(mo_bra_[k]*x);
		for(size_t l=0;l<mo_ket_.size();l++){
			double ijgkl = igk.inner(mo_bra_[j]*mo_ket_[l]);
			double klfxy = kfx.inner(mo_bra_[l]*y);
			ijgO12fxy += ijgkl*klfxy;
		}
	}

	return (ijgfxy - ijgO1fxy - ijgO2fxy + ijgO12fxy);
}

double CC_Operators::make_ijgfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const{
	real_function_3d jy = mo_bra_[j]*y;
	real_function_3d ix = mo_bra_[i]*x;
	// I12 Part:
	double ijgfxy = (ix).inner(apply_gf(jy));
	return ijgfxy;
}

/// General Function to make the two electron integral <ij|g|xy>
/// For Debugging -> Expensive without intermediates
double CC_Operators::make_ijgxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const{
	real_function_3d igx = (*poisson)(mo_bra_[i]*x).truncate();
	real_function_3d jy = (mo_bra_[j]*y).truncate();
	return jy.inner(igx);
}

/// General Function to make two electron integrals with pair functions (needed for energy)
double CC_Operators::make_ijgu(const size_t &i, const size_t &j, const CC_Pair &u)const{
	return make_ijgu(i,j,u.function);
}
double CC_Operators::make_ijgu(const size_t &i, const size_t &j, const real_function_6d &u)const{
	real_function_6d g = TwoElectronFactory(world).dcut(parameters.lo);
	real_function_6d ij_g =
			CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[i])).particle2(
							copy(mo_bra_[j])).g12(g);

	// compute < ij | g12 | u >
	const double ij_g_u = inner(u, ij_g);
	return ij_g_u;
}

/// General Function to make two electorn integrals with pair function and orbitals and the BSH Operator (needed for gf = g - bsh)
double CC_Operators::make_ijGu(const size_t &i, const size_t &j, const CC_Pair &u)const{
	real_function_6d G = TwoElectronFactory(world).BSH().gamma(corrfac.gamma()).dcut(parameters.lo);
	double bsh_prefactor = 4.0 * constants::pi;
	real_function_6d ij_G =
			CompositeFactory<double, 6, 3>(world).particle1(
					copy(mo_bra_[i])).particle2(
							copy(mo_bra_[j])).g12(G);

	// compute < ij | g12 | u >
	const double ij_G_u = inner(u.function, ij_G);
	return bsh_prefactor*ij_G_u;
}




/// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
/// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
real_function_3d CC_Operators::apply_gf(const real_function_3d &f)const{
	double bsh_prefactor = 4.0 * constants::pi;
	double prefactor = 1.0/(2.0*corrfac.gamma());
	return prefactor*((*poisson)(f) - bsh_prefactor*(*fBSH)(f)).truncate();
}
real_function_6d CC_Operators::apply_gf(const real_function_6d &f,const size_t &particle)const{
	poisson->particle()=particle;
	fBSH->particle()=particle;
	double bsh_prefactor = 4.0 * constants::pi;
	double prefactor = 1.0/(2.0*corrfac.gamma());
	return prefactor*((*poisson)(f) - bsh_prefactor*(*fBSH)(f)).truncate();
}

/// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
/// 1. <x|f12|z>*|y>
/// 2. -\sum_m <x|m> <m|f12|z>*|y> = -(\sum_m <x|m> <m|f12|z>) *|y>
/// 3. -\sum_n <nx|f12|zy> * |n>
/// 4. +\sum_{mn} <x|n> <mn|f12|yz> * |m>
/// Metric from nuclear cusp is not taken into account -> give the right bra elements to the function
real_function_3d CC_Operators::convolute_x_gQf_yz(const real_function_3d &x, const real_function_3d &y, const real_function_3d &z)const{
	// make intermediates
	vecfuncT moz = mul(world,z,mo_bra_);
	vecfuncT moy = mul(world,y,mo_bra_);
	// Do <x|f12|z>*|y>
	real_function_3d part1_tmp = (*f12op)(x*z);
	real_function_3d part1 = (part1_tmp*y).truncate();
	// Do -\sum_m <x|m> <m|f12|z>*|y>
	real_function_3d part2 = real_factory_3d(world);
	{
		Tensor<double> xm = inner(world,x,mo_ket_);
		vecfuncT f12mz = apply(world,*f12op,moz);
		real_function_3d tmp = real_factory_3d(world);
		for(size_t m=0;m<mo_bra_.size();m++) tmp += xm[m]*f12mz[m];
		part2 = tmp*y;
	}
	// Do -\sum_n <nx|f12|zy> * |n> |  <nx|f12|zy> = <n| xf12y |z>
	real_function_3d part3 = real_factory_3d(world);
	{
		real_function_3d xf12y = (*f12op)((x*y).truncate());
		for(size_t n=0;n<mo_bra_.size();n++){
			double nxfzy = xf12y.inner(mo_bra_[n]*z);
			part3 += nxfzy*mo_ket_[n];
		}
	}
	// Do +\sum_{mn} <x|n> <mn|f12|yz> * |m>
	real_function_3d part4 = real_factory_3d(world);
	{
		Tensor<double> xn = inner(world,x,mo_ket_);
		vecfuncT nf12z = apply(world,*f12op,moz);
		Tensor<double> mnfyz = inner(world,moy,nf12z);
		for(size_t m=0;m<mo_bra_.size();m++){
			for(size_t n=0;n<mo_ket_.size();n++){
				part4 += xn(n)*mnfyz(m,n)*mo_ket_[m];
			}
		}
	}
	real_function_3d result = part1 - part2 - part3 + part4;
	result.truncate();

	if(parameters.debug){
		CC_Timer function_debug(world,"Debug-Time for <k|Qf|xy>");
		real_function_6d test_tmp = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(y)).particle2(copy(z));
		test_tmp.fill_tree().truncate().reduce_rank();
		real_function_6d test_1 = copy(test_tmp);
		real_function_6d test_4 = copy(Q12(test_tmp));
		real_function_3d test1 = test_1.project_out(x,1);
		real_function_3d test4 = test_4.project_out(x,1);
		double check1 = (test1 - part1).norm2();
		double check4 = (test4 - result).norm2();
		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check1 value is: " << check1 << std::endl;
		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check4 value is: " << check4 << std::endl;
		if(check1 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check1 failed");
		if(check4 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check4 failed");
		function_debug.info();
	}

	return result;
}






}
