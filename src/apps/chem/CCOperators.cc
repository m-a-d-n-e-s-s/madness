/*
 * CCOperators.cc
 *
 *  Created on: Jul 6, 2015
 *      Author: kottmanj
 */


#include "CCOperators.h"

#include "/usr/include/math.h"
#include "/usr/include/math.h"
#include "../../madness/constants.h"
#include "../../madness/mra/derivative.h"
#include "../../madness/mra/funcdefaults.h"
#include "../../madness/mra/funcimpl.h"
#include "../../madness/mra/funcplot.h"
#include "../../madness/mra/function_factory.h"
#include "../../madness/mra/functypedefs.h"
#include "../../madness/mra/mra.h"
#include "../../madness/mra/operator.h"
#include "../../madness/mra/vmra.h"
#include "../../madness/tensor/srconf.h"
#include "../../madness/tensor/tensor.h"
#include "../../madness/world/madness_exception.h"
#include "../../madness/world/parallel_archive.h"
#include "../../madness/world/print.h"
#include "../../madness/world/world.h"
#include "electronic_correlation_factor.h"
#include "TDA.h"

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

real_function_3d CC_Intermediates::make_density(const CC_vecfunction &bra,
		const CC_vecfunction &ket) const {
	if (bra.size()==0) error("error in make_density: bra_element is empty");
	if (ket.size()==0) error("error in make_density: ket_element is empty");
	// make the density
	real_function_3d density = real_factory_3d(world);
	for (auto x:ket.functions){

		density += (bra(x.first).function * x.second.function);
	}
	density.truncate(FunctionDefaults<3>::get_thresh()*0.01);
	return density;
}

intermediateT CC_Intermediates::make_exchange_intermediate(const CC_vecfunction &bra,
		const CC_vecfunction &ket) const {
	intermediateT xim;
	for(auto tmpk:bra.functions){
		const CC_function & k = tmpk.second;
		for(auto tmpl:ket.functions){
			const CC_function& l=tmpl.second;
			real_function_3d kl = (bra(k).function*l.function);
			real_function_3d result = ((*poisson)(kl)).truncate();
			xim.insert(k.i,l.i,result);
		}
	}
	return xim;
}

intermediateT CC_Intermediates::make_f12_exchange_intermediate(const CC_vecfunction &bra,
		const CC_vecfunction &ket)const{
	intermediateT xim;
	for(auto tmpk:bra.functions){
		const CC_function & k = tmpk.second;
		for(auto tmpl:ket.functions){
			CC_function l=tmpl.second;
			real_function_3d kl = (bra(k).function*l.function);
			real_function_3d result = ((*f12op)(kl)).truncate();
			xim.insert(k.i,l.i,result);
		}
	}
	return xim;
}

double CC_Operators::compute_mp2_pair_energy(CC_Pair &pair)const{

	const size_t i = pair.i;
	const size_t j = pair.j;

	// this will be the bra space
	real_function_6d eri = TwoElectronFactory(world).dcut(parameters.lo);
	real_function_6d ij_g =CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(i).function)).particle2(copy(mo_bra_(j).function)).g12(eri);
	real_function_6d ji_g =CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(j).function)).particle2(copy(mo_bra_(i).function)).g12(eri);

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
		pair.e_singlet = (ij_g_uij + pair.ij_gQf_ij)+ (ji_g_uij + pair.ji_gQf_ij);
		pair.e_triplet = 3.0* ((ij_g_uij - ji_g_uij) + (pair.ij_gQf_ij - pair.ji_gQf_ij));
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
// the fock residue R= 2J-K+Un for closed shell is computed here
// J_i = \sum_k <k|r12|k> |tau_i>
// K_i = \sum_k <k|r12|tau_i> |k>
vecfuncT CC_Operators::fock_residue_closed_shell(const CC_vecfunction &singles) const {
	//	vecfuncT tau = singles.get_vecfunction();
	//CC_Timer timer_J(world,"J");
	//	vecfuncT J = mul(world, intermediates_.get_hartree_potential(), tau);
	vecfuncT J;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second;
		const real_function_3d Ji = intermediates_.get_hartree_potential()*i.function;
		if(parameters.debug){
			real_function_3d Jpot = real_factory_3d(world);
			for(const auto ktmp:mo_ket_.functions){
				const size_t k=ktmp.first;
				output("J Debug: k="+stringify(k));
				Jpot += mo_bra_(k).function*mo_ket_(k).function;
			}
			const real_function_3d J_check = (*poisson)(Jpot)*i.function;
			const double diff = (Ji-J_check).norm2();
			std::cout << "J-debug-difference-is:" << diff << std::endl;
		}
		J.push_back(Ji);
	}
	truncate(world, J);
	scale(world, J, 2.0);
	//timer_J.info();
	//CC_Timer timer_K(world,"K");
	vecfuncT vK;
	for(const auto& tmpi:singles.functions){
		const CC_function& taui=tmpi.second;
		const real_function_3d Ki = K(taui);
		if(parameters.debug){
			const real_function_3d Ki_check = K(CC_function(taui.function,taui.i,UNDEFINED));
			const double diff=(Ki-Ki_check).norm2();
			std::cout << "K-debug-difference-is:" << diff << std::endl;
		}
		vK.push_back(Ki);
	}
	scale(world, vK, -1.0);
	//timer_K.info();

	// apply nuclear potential
	Nuclear Uop(world,&nemo);
	vecfuncT Upot = Uop(singles.get_vecfunction());
	vecfuncT KU = add(world,vK,Upot);

	return add(world, J, KU);
}

// The coulomb Term of the S3C diagram: Positive sign
// \     /
//  \---/  = 2Q\sum_j(<j|g12|tau_j>)|i>
//  _\_/_
vecfuncT CC_Operators::S3c(const CC_vecfunction &singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d tmp = intermediates_.get_perturbed_hartree_potential()*mo_ket_(i).function;
		result.push_back(tmp);
	}
	scale(world, result, 2.0);
	return result;
}

// The Exchange Term of the S3C diagram: Negative sign
// \  /
//  \/...   = -Q\sum_k(<k|g12|i>|tau_k>)
//     / \
//    _\_/_
vecfuncT CC_Operators::S3c_X(const CC_vecfunction &singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto& tmpk:singles.functions){
			const CC_function& k=tmpk.second;
			resulti += (intermediates_.get_EX(k,i)*k.function);
		}
		result.push_back(resulti);
	}
	scale(world, result, -1.0);
	return result;
}

/// The S5b term
//[i]    [Q]
// \     /....
//  \   /   / \
//  _\_/_  _\_/_
// 2\sum_k <k|g|\tau_k> |\tau_i>
// No Q is applied yet !
vecfuncT CC_Operators::S5b(const CC_vecfunction &singles) const {
	real_function_3d kgtauk = intermediates_.get_perturbed_hartree_potential();
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i = tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = (kgtauk*i.function);
		result.push_back(resulti);
	}
	scale(world,result,2.0);
	return result;
}

/// The S5b Exchange Term
//[i]         [Q]
// \     ...../
//  \   /\   /
//  _\_/  \_/_
// -\sum_k <k|g|\tau_i> |\tau_k>
// No Q is applied yet !
vecfuncT CC_Operators::S5b_X(const CC_vecfunction &singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto& tmpk:singles.functions){
			const CC_function& k=tmpk.second;
			resulti += (intermediates_.get_pEX(k,i)*k.function);
		}
		result.push_back(resulti);
	}
	scale(world, result, -1.0);
	return result;
}

/// The S5c term
//[Q]    [i]
// \     /....
//  \   /   / \
//  _\_/_  _\_/_
// 2\sum_kl <kl|g|i\tau_l> |\tau_k>
// No Q is applied yet !
// May use alteriative algorithm with perturbed density intermediate
vecfuncT CC_Operators::S5c(const CC_vecfunction&singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto& tmpk:singles.functions){
			const CC_function& k=tmpk.second;
			for(const auto& tmpl:singles.functions){
				const CC_function& l=tmpl.second;
				real_function_3d tmp = (intermediates_.get_EX(k,i)*l.function);
				double integral = mo_bra_(l).inner(tmp);
				resulti += integral*k.function;
			}
		}
		result.push_back(resulti);
	}
	scale(world, result, 2.0);
	return result;
}

/// The S5c_X echange term
//[Q]         [i]
// \     ...../
//  \   /\   /
//  _\_/  \_/_
// -\sum_kl <lk|g|i\tau_l> |\tau_k>
// No Q is applied yet !
vecfuncT CC_Operators::S5c_X(const CC_vecfunction &singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto& tmpk:singles.functions){
			const CC_function& k=tmpk.second;
			for(const auto& tmpl:singles.functions){
				const CC_function& l=tmpl.second;
				real_function_3d tmp = (intermediates_.get_EX(l,i)*l.function);
				double integral = mo_bra_(k).inner(tmp);
				resulti += integral*k.function;
			}
		}
		result.push_back(resulti);
	}
	scale(world, result, -1.0);
	return result;
}

/// The S6+X Term
// \    /\    /...
//  \  /  \  /   /\
//  _\/_  _\/_  _\/_
// Q (\sum_kl 2<kl|g|\tau_k\tau_i> |\tau_l> - \sum_kl <kl|g|\taui\tau_k> |\tau_l>)
// Q is not applied yet!
// overall minus sign is added at the potential_singles function
vecfuncT CC_Operators::S6(const CC_vecfunction &singles) const {
	vecfuncT result;
	for(const auto& tmpi:singles.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto& tmpk:singles.functions){
			const CC_function& k=tmpk.second;
			for(const auto& tmpl:singles.functions){
				const CC_function& l=tmpl.second;
				real_function_3d tmp_J = (intermediates_.get_pEX(k,k)*i.function);
				real_function_3d tmp_K = (intermediates_.get_pEX(k,i)*k.function);
				double integral_J = mo_bra_(l).inner(tmp_J);
				double integral_K = mo_bra_(l).inner(tmp_K);
				resulti += (2.0*integral_J - integral_K)*l.function;
			}
		}
		result.push_back(resulti);
	}
	return result;
}

/// CC2 singles diagrams with 6d functions as input
/// Use GFInterface in function_interface.h as kernel (f*g) and do not reconstruct \tau = f12u(1,2) if possible
/// Since the correlation factor of CC2 has Slater form like in MP2: g12f12 = g12(1-exp(-mu*r12)/r12) = g12 - exp(-mu*r12)/r12 = Coulomb_Operator - BSH_Operator(mu)

/// the S2c and S4b part together
vecfuncT CC_Operators::S2c4b(const Pairs<CC_Pair> &doubles, const CC_vecfunction &singles, CC_data &data)const{
	if(parameters.pair_function_in_singles_potential==FULL){
		return S2c4b_6D_part(doubles,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S2c4b_6D_part(doubles,singles,data),S2c4b_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT();
	}
}

vecfuncT CC_Operators::S2c4b_6D_part(const Pairs<CC_Pair> &doubles, const CC_vecfunction &singles, CC_data &data)const{

	vecfuncT result;
	for(const auto& itmp:singles.functions){
		const CC_function& i = itmp.second;
		real_function_3d resulti = real_factory_3d(world);
		for(const auto & ltmp:singles.functions){
			const CC_function & l = ltmp.second;

			// make the <l|g|ti+i> intermediate
			const real_function_3d lgti = intermediates_.get_EX(l,i) + intermediates_.get_pEX(l,i);

			for(const auto & ktmp:singles.functions){
				const CC_function & k = ktmp.second;

				const real_function_3d klgti = mo_bra_(k).function*lgti;
				const real_function_6d ukl = get_pair_function(doubles,k.i,l.i);

				const real_function_3d tmp = 2.0*ukl.project_out(klgti,0) - ukl.project_out(klgti,1);
				resulti += tmp;
			}
		}
		result.push_back(resulti);
	}
	return result;
}

// part0: the part without projectors:
// 2.0*<klgti|f|k> |l> - <klgti|f|l> |k>
// part1: the part with projector on particle 1:
// sum_m <klgti|m> <m|f|k> |l> - <klgti|m> <m|f|l> |k>
// part2: the part with projector on particle 2:
// sum_n <klgti|nfl|k> |n> bzw <klgti|mfl|k> |m>
// part3: the part with both projectors:
// sum_mn <klgti|m> <mn|f|kl> |n>
// result = part0 - part1 - part2 + part3
vecfuncT CC_Operators::S2c4b_3D_part(const CC_vecfunction &singles, CC_data &data)const{

	vecfuncT result;
	for(const auto& itmp:singles.functions){
		const CC_function& i = itmp.second;
		real_function_3d part0 = real_factory_3d(world);
		real_function_3d part1 = real_factory_3d(world);
		real_function_3d part2 = real_factory_3d(world);
		real_function_3d part3 = real_factory_3d(world);

		for(const auto & ltmp:singles.functions){
			const CC_function & l = ltmp.second;
			const real_function_3d mol = mo_ket_(l).function;

			const real_function_3d lgti = intermediates_.get_EX(l,i) + intermediates_.get_pEX(l,i);

			for(const auto & ktmp:singles.functions){
				const CC_function & k = ktmp.second;
				const real_function_3d mok = mo_ket_(k).function;
				const real_function_3d klgti = mo_bra_(k).function*lgti;

				const real_function_3d klgtik = klgti*mok;
				const real_function_3d klgtil = klgti*mol;
				const real_function_3d part0_tmp  = (*f12op)(klgtik);
				const real_function_3d part0_tmpx = (*f12op)(klgtil);
				part0 += (2.0*part0_tmp*mol - part0_tmpx*mok);

				for(const auto & mtmp:mo_ket_.functions){
					const CC_function & m = mtmp.second;
					const real_function_3d & mom = m.function;
					const double klgtim = klgti.inner(mom);
					const real_function_3d mfk = intermediates_.get_fEX(m,k);
					const real_function_3d mfl = intermediates_.get_fEX(m,l);
					part1 = 2.0*klgtim*mfk*mol - klgtim*mfl*mok;


					const real_function_3d part2_tmp  = klgti*mfl*mol;
					const real_function_3d part2_tmpx = klgti*mfk*mok;
					part2 += (2.0*part2_tmp - part2_tmpx)*mom ;


					for(const auto & ntmp:mo_ket_.functions){
						const CC_function & n = ntmp.second;
						const double mnfkl = mo_bra_(n).inner(mfk*mol);
						const double mnflk = mo_bra_(n).inner(mfl*mok);
						part3 += (2.0*mnfkl-mnfkl)*klgtim*n.function;
					}
				}
			}
		}
		real_function_3d resulti = part0 + part1 + part2 + part3;
		result.push_back(resulti);
	}
	return result;
}


/// S2b + X Term
// [i]   [Q]
//  \    /....
//   \  /    /\
//  __\/_____\/__
/// @param[in] Structure which holds all current CC pair functions
/// @param[in] Structure which holds all current CC single excitations
/// @param[out] Q\sum_k \left( 2<k|g|u_ik> - <k|g|u_ki> + 2<k|gQf|t_it_k> - <k|gQf|t_kt_i> \right), with t_i = i + \tau_i
/// notation: <k|g|u_ik> = <k(2)|g12|u_ik(1,2)> (Integration over second particle)
vecfuncT CC_Operators::S2b(const Pairs<CC_Pair> u, const CC_vecfunction & singles , CC_data &data) const {
	if(parameters.pair_function_in_singles_potential==FULL){
		return S2b_6D_part(u,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S2b_6D_part(u,singles,data),S2b_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT();
	}
}
vecfuncT CC_Operators::S2b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & singles , CC_data &data) const {
	// check if the 6D part has been calculated before
	if(not current_S2b_6D_part_.empty()){
		output("loading previously calculated S2b_6D_part");
		return current_S2b_6D_part_;
	}else{

		vecfuncT result;
		for(const auto& tmpi:singles.functions){
			const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
			real_function_3d resulti = real_factory_3d(world);
			CC_Timer S2b_6D_time(world,"S2b: 6D-Part internal");
			for(const auto& tmpk:singles.functions){
				const CC_function& k=tmpk.second;
				{
					const real_function_6d uik = get_pair_function(u,i.i,k.i);

					const real_function_6d tmp = multiply(copy(uik),copy(mo_bra_(k).function),2);
					poisson -> particle() = 2;
					const real_function_6d result = (*poisson)(tmp);
					const real_function_3d r = result.dirac_convolution<3>();

					// Exchange Part
					const real_function_6d tmpx = multiply(copy(uik),copy(mo_bra_(k).function),1);
					poisson -> particle() = 1;
					const real_function_6d resultx= (*poisson)(tmpx);
					const real_function_3d rx = resultx.dirac_convolution<3>();

					const real_function_3d resulti_increment  = (2.0*r -rx);

					resulti_increment.print_size("s2b_i="+stringify(i.i)+"_k="+stringify(k.i));

					resulti += resulti_increment;
				}
			}
			S2b_6D_time.info();
			result.push_back(resulti);
		}
		current_S2b_6D_part_ = result;
		return result;
	}// end else
}
vecfuncT CC_Operators::S2b_3D_part(const CC_vecfunction & singles , CC_data & data) const {
	CC_vecfunction t = make_t_intermediate(singles);
	real_function_3d tdensity = intermediates_.make_density(mo_bra_,t);
	vecfuncT result;
	CC_Timer S2b_3D_time(world,"S2b: 3D-Part");
	for(const auto& tmpi:t.functions){
		const CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		// The Part which operators on the t intermediate OPTIMIZE THE LOOPS (for helium debug not necessary)
		real_function_3d t1part;
		real_function_3d t1partx;

		{

			// 2<k|gf|t_it_k> - <k|gf|t_kt_i> // factors are added in the end
			real_function_3d unitpart = (apply_gf(tdensity)*i.function);
			real_function_3d unitpartx = real_factory_3d(world);
			for(const auto& tmpk:t.functions){
				const CC_function& k=tmpk.second;
				unitpartx += (apply_gf(mo_bra_(k).function*i.function)*k.function);
			}

			// make intermediates
			// <n|f|ti>
			vecfuncT nfti;
			for(const auto& motmp:mo_bra_.functions){
				const CC_function & mo = motmp.second;
				nfti.push_back( ((*f12op)(mo.function*i.function)));
			}


			// 2 \sum_n <k(2)|g12|n(1)><n(1)|f12|ti(1)tk(2)> = <k(2)|g12 nfti(2)|n(1)tk(2)> = poisson(k*nfti*tk)*n
			//   \sum_n <k(2)|g12|n(1)><n(1)|f12|tk(1)ti(2)> = <k(2)|g12 nftk(2)|n(1)ti(2)> = poisson(k*nftk*ti)*n
			// The sum over n also goes over the frozen orbitals
			real_function_3d O1part = real_factory_3d(world);
			real_function_3d O1partx = real_factory_3d(world);
			for(const auto& tmpk:t.functions){
				const CC_function& k=tmpk.second;
				for(size_t n=0;n<nfti.size();n++){
					O1part += ((*poisson)(mo_bra_(k).function*k.function*nfti[n])*mo_ket_(n).function);
					real_function_3d nftk = ((*f12op)(mo_bra_(n).function*k.function));
					real_function_3d poisson_integrant = (mo_bra_(k).function*i.function*nftk);
					O1partx+= ((*poisson)(poisson_integrant)*mo_ket_(n).function);
				}
			}

			// 2 \sum_n <k(2)|g12|n(2)><n(2)|f12|ti(1)tk(2)> = <k(2)|g12 nftk(1) |ti(1)n(2)> = poisson(k*n)*(nftk*ti)
			//   \sum_n <k(2)|g12|n(2)><n(2)|f12|tk(1)ti(2)> = <k(2)|g12 nfti(1) |tk(1)n(2)> = poisson(k,n)*(nfti*tk)
			real_function_3d O2part = real_factory_3d(world);
			real_function_3d O2partx = real_factory_3d(world);
			for(auto tmpk:t.functions){
				CC_function& k=tmpk.second;
				for(size_t n=0;n<mo_ket_.size();n++){
					real_function_3d kgn = intermediates_.get_EX(k.i,n);
					real_function_3d nftk = (*f12op)(mo_bra_(n).function*k.function);
					real_function_3d nfti = (*f12op)(mo_bra_(n).function*i.function);
					O2part   += (kgn*nftk*i.function);
					O2partx  += (kgn*nfti*k.function);
				}
			}

			// 2 \sum_nm <k(2)|g12|mn><mn|f12|ti(1)tk(2)> = <k|g|n>(1)|m(1)> *<mn|f|titk>
			//   \sum_nm <k(2)|g12|mn><mn|f12|tk(1)ti(2)> = <k|g|n>(1)|m(1)> *<mn|f|tkti>
			real_function_3d O12part = real_factory_3d(world);
			real_function_3d O12partx = real_factory_3d(world);
			for(auto tmpk:t.functions){
				CC_function& k=tmpk.second;
				for(size_t n=0;n<mo_bra_.size();n++){
					for(size_t m=0;m<mo_ket_.size();m++){
						real_function_3d kgn = intermediates_.get_EX(k.i,n);
						real_function_3d nftk = (*f12op)(mo_bra_(n).function*k.function);
						real_function_3d nfti = (*f12op)(mo_bra_(n).function*i.function);
						double f = nftk.inner(mo_bra_(m).function*i.function);
						double fx= nfti.inner(mo_bra_(m).function*k.function);
						O12part += (f*kgn*mo_ket_(m).function);
						O12partx +=(fx*kgn*mo_ket_(m).function);
					}
				}
			}

			t1part = unitpart - O1part - O2part +O12part;
			t1partx = unitpartx - O1partx - O2partx + O12partx;
//			(2.0*unitpart).print_size("unitpart");
//			(2.0*O1part).print_size("O1part");
//			(2.0*O2part).print_size("O2part");
//			(2.0*O12part).print_size("O12part");
//			unitpartx.print_size("unitpartx");
//			O1partx.print_size("O1partx");
//			O2partx.print_size("O2partx");
//			O12partx.print_size("O12partx");
		}
		resulti += 2.0*t1part-t1partx;
		result.push_back(resulti);
	}
	S2b_3D_time.info();
	return result;
}

/// S2c + X Term
// [Q]   [i]
//  \    /....
//   \  /    /\
//  __\/_____\/__
/// = Q\sum_{kl}\left( 2<k|lgi|ulk> - <l|kgi|u_{lk}> + 2<k|lgiQ12f12|t_lt_k> - <l|kgiQ12f12|t_lt_k> \right)
/// Notation: 6D Integration over second particle, intermediates: lgi = <l|g|i> is the exchange intermediate
/// Notation: t are the t-intermediates: |t_i> = |i> + |\tau_i>
/// @param[in] All the current coupled cluster Pairs
/// @param[in] The coupled cluster singles
/// @param[out] the S2c+X Potential
vecfuncT CC_Operators::S2c(const Pairs<CC_Pair> &u, const CC_vecfunction &singles,CC_data &data) const {
	if(parameters.pair_function_in_singles_potential==FULL){
		return S2c_6D_part(u,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S2c_6D_part(u,singles,data),S2c_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT(); // should not ever get here
	}
}
vecfuncT CC_Operators::S2c_3D_part(const CC_vecfunction &singles, CC_data &data) const {
	vecfuncT result;
	CC_vecfunction t = make_t_intermediate(singles);
	for(auto tmpi:singles.functions){
		CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		// The 3D Part (maybe merge this later with 6D part to avoid recalculating of  mo_bra_[k]*exchange_intermediate[i][l]
		real_function_3d resulti_3d = real_factory_3d(world);
		CC_Timer s2c_3d_timer(world,"s2c_3D_part");
		{
			for(auto tmpk:t.functions){
				CC_function& k=tmpk.second;
				for(auto tmpl:t.functions){
					CC_function& l=tmpl.second;
					real_function_3d klgi = (mo_bra_(k)*intermediates_.get_EX(l,i)).function;
					real_function_3d lkgi = (mo_bra_(l)*intermediates_.get_EX(k,i)).function;
					real_function_3d k_lgi_Q12f12_tltk = convolute_x_Qf_yz(CC_function(klgi,99,UNDEFINED),l,k);
					real_function_3d l_kgi_Q12f12_tltk = convolute_x_Qf_yz(CC_function(lkgi,99,UNDEFINED),l,k);
					resulti_3d += 2.0*k_lgi_Q12f12_tltk - l_kgi_Q12f12_tltk;
				}
			}
		}
		s2c_3d_timer.info();
		result.push_back(resulti_3d);
	}
	return result;
}
vecfuncT CC_Operators::S2c_6D_part(const Pairs<CC_Pair> &u, const CC_vecfunction &singles, CC_data &data) const {
	vecfuncT result;
	for(auto tmpi:singles.functions){
		CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		{
			for(auto tmpk:singles.functions){
				CC_function& k=tmpk.second;
				for(auto tmpl:singles.functions){
					CC_function& l=tmpl.second;
					real_function_6d ukl = get_pair_function(u,k.i,l.i); // needed to avoid seg.faults when k>l

					const real_function_3d klgi = (mo_bra_(k).function*intermediates_.get_EX(l,i));
					const real_function_3d lkgi = (mo_bra_(l).function*intermediates_.get_EX(k,i));
					const real_function_3d tmp = ukl.project_out(lkgi,1); // 1 means second particle
					const real_function_3d tmpx= ukl.project_out(klgi,1); // 1 means second particle
					const real_function_3d resulti_increment = 2.0*tmp - tmpx;
					resulti += resulti_increment;
				}
			}
		}

		result.push_back(resulti);
	}
	return result;
}


/// The S4a + X diagram
//[Q]       [i]
// \    ..../.....
//  \  /\  /     /\
//  _\/_ \/______\/_
/// Q\sum (2<kl|g|\tau_il>|\tau_k> - <lk|g|\tau_il>|\tau_k>)
/// 6D Part -Q\sum_k ( 2<kl|g|uil> - <lk|g|uil> ) tau_k
/// 3D Part t1 intermediate integras: <kl|g|uil> --> <kl|gQf|titl>
/// 3D Part decomposed integrals <kl|g|uil> --> <kl|gQf|xy> , |xy> = |il> + |ti,l> + |i,tl> + |ti,tl>
vecfuncT CC_Operators::S4a(const Pairs<CC_Pair> u, const CC_vecfunction & singles , CC_data &data) const {
	if(parameters.pair_function_in_singles_potential==FULL){
		return S4a_6D_part(u,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S4a_6D_part(u,singles,data),S4a_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT(); // should not ever get here
	}
}
vecfuncT CC_Operators::S4a_3D_part(const CC_vecfunction & singles,  CC_data &data) const {
	vecfuncT result;
	CC_vecfunction t = make_t_intermediate(singles);
	for (auto tmpi:t.functions) {
		CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for (auto tmpk:singles.functions) {
			CC_function& k=tmpk.second;
			for (auto tmpl:t.functions) {
				CC_function& l=tmpl.second;
				double a = 2.0*make_ijgQfxy(k.i,l.i,i.function,l.function) - make_ijgQfxy(l.i,k.i,i.function,l.function);
				real_function_3d tmp_result =a*k.function;
				resulti += tmp_result;
			}
		}
		result.push_back(resulti);
	}
	return result;
}
vecfuncT CC_Operators::S4a_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & singles, CC_data &data) const {
	vecfuncT result;
	for (auto tmpi:singles.functions) {
		CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for (auto tmpk:singles.functions) {
			CC_function& k=tmpk.second;
			for (auto tmpl:singles.functions) {
				CC_function& l=tmpl.second;

				real_function_6d uil = get_pair_function(u,i.i,k.i); // needed to avoid seg.faults when k>l

				double integral = make_ijgu(k.i,l.i,uil);
				double integralx= make_ijgu(l.i,k.i,uil);
				if(k.i==l.i and (integral-integralx)>FunctionDefaults<3>::get_thresh()) warning("Warnings in S4a_6D_part: J and K integral differ for k==l :" + stringify(integral) + " and " + stringify(integralx));
				resulti += (2.0*integral - integralx)*k.function;
			}
		}
		result.push_back(resulti);
	}
	return result;
}

/// The S4b -> Merge this with S2c later to save time
//[i]       [Q]
// \    ..../.....
//  \  /\  /     /\
//  _\/_ \/______\/_
/// -Q\sum_kl ( 2<k(3)l(4)|g34|\taui(3)\taukl(1,4)> - <k(3)l(4)|g34|taui(4)\taukl(3,1)>  )
/// procedure for 6D part:
/// make intermediate kgtaui(4) = <k|g|\taui> and then  l_kgtaui = <l| * kgtaui
/// project out:  <l_kgtaui(4)|taukl(1,4)>
/// do the same for exchange part  <k_lgtaui(3)|taukl(3,1)>
/// procedure for 3D part:
/// make intermediate
vecfuncT CC_Operators::S4b(const Pairs<CC_Pair> u, const CC_vecfunction & singles , CC_data &data) const {
	if(parameters.pair_function_in_singles_potential==FULL){
		return S4b_6D_part(u,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S4b_6D_part(u,singles,data),S4b_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT(); // should not ever get here
	}
}
vecfuncT CC_Operators::S4b_3D_part(const CC_vecfunction & singles,  CC_data &data) const {
	vecfuncT result;
	CC_vecfunction t = make_t_intermediate(singles);
	for (auto tmpi:singles.functions){ CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
	real_function_3d resulti = real_factory_3d(world);
	for(auto tmpk:t.functions){ CC_function& k=tmpk.second;
	for(auto tmpl:t.functions){ CC_function& l=tmpl.second;
	const CC_function l_kgti = mo_bra_(l)*intermediates_.get_pEX(k,i);
	const CC_function k_lgti = mo_bra_(k)*intermediates_.get_pEX(l,i);

	real_function_3d all_parts = convolute_x_Qf_yz(l_kgti,k,l);
	real_function_3d all_xparts= convolute_x_Qf_yz(l_kgti,l,k);
	resulti = 2.0*all_parts - all_xparts;
	}
	}
	result.push_back(resulti);
	}
	truncate(world,result);
	return result;
}
/// for explanation see above
vecfuncT CC_Operators::S4b_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & singles, CC_data &data) const {
	vecfuncT result;
	for(auto tmpi:singles.functions){ CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
	real_function_3d resulti = real_factory_3d(world);
	for(auto tmpk:singles.functions){ CC_function& k=tmpk.second;
	for(auto tmpl:singles.functions){ CC_function& l=tmpl.second; // index l is also pair function index therefore it runs only over non frozen orbs

	real_function_6d ukl = get_pair_function(u,k.i,l.i); // needed to avoid seg.fault when k>l

	real_function_3d k_lgtaui = (mo_bra_(k)*intermediates_.get_pEX(l,i)).function;
	real_function_3d l_kgtaui = (mo_bra_(l)*intermediates_.get_pEX(k,i)).function;
	real_function_3d tmp  = ukl.project_out(k_lgtaui,1); // 1 is particle 2
	real_function_3d tmpx = ukl.project_out(l_kgtaui,0); // 0 is particle 1
	resulti += (2.0*tmp - tmpx).truncate();
	}
	}
	result.push_back(resulti);
	}
	return result;
}

/// The S4c + X + X + X + X Diagrams -> merge this with S4b later to save time
//            [i]   [Q]
//   .......   \    /
//  /\     /\   \  /
// _\/_   _\/____\/_
/// Q\sum_{kl}[ 4*<k(3)l(4)|g34| \tau_k(3) \tau_{il}(1,4)> - 2* <k(3)l(4)|g34|\tau_k(4) \tau_{il}(1,3)>
/// -           2*<k(3)l(4)|g34| \tau_k(3) \tau_{il}(4,1)>	    +<k(3)l(4)|g34|\tau_k(4) \tau_{il}(3,1)>  ]

/// 3D part <l*kgtk(2)|Q12f12|xi(1)yl(2)> would be the 4.0* part
/// Similar than S4b
vecfuncT CC_Operators::S4c(const Pairs<CC_Pair> u, const CC_vecfunction & singles , CC_data &data) const {
	if(parameters.pair_function_in_singles_potential==FULL){
		return S4c_6D_part(u,singles,data);
	}else if(parameters.pair_function_in_singles_potential==DECOMPOSED){
		return add(world,S4c_6D_part(u,singles,data),S4c_3D_part(singles,data));
	}else{
		error("parameters.pair_function_in_singles_potential definition not valid : " + stringify(parameters.pair_function_in_singles_potential));
		return vecfuncT(); // should not ever get here
	}
}
vecfuncT CC_Operators::S4c_3D_part(const CC_vecfunction & singles,  CC_data &data) const {
	vecfuncT result;
	CC_vecfunction t = make_t_intermediate(singles);
	for(auto tmpi:t.functions){ CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
	real_function_3d resulti = real_factory_3d(world);
	for(auto tmpk:singles.functions){ CC_function& k=tmpk.second;
	for(auto tmpl:t.functions){ CC_function& l=tmpl.second;
	CC_function l_kgtk(mo_bra_(l)*intermediates_.get_pEX(k,k));
	CC_function k_lgtk(mo_bra_(k)*intermediates_.get_pEX(l,k));

	real_function_3d part1 = convolute_x_Qf_yz(l_kgtk,i,l);
	real_function_3d part2 = convolute_x_Qf_yz(k_lgtk,i,l);
	real_function_3d part3 = convolute_x_Qf_yz(l_kgtk,l,i);
	real_function_3d part4 = convolute_x_Qf_yz(k_lgtk,l,i);
	real_function_3d all = 4.0*part1 - 2.0*part2 - 2.0*part3 + part4;
	resulti+=all;
	}
	}
	result.push_back(resulti);
	}
	truncate(world,result);
	return result;
}
vecfuncT CC_Operators::S4c_6D_part(const Pairs<CC_Pair> u, const CC_vecfunction & singles, CC_data &data) const {
	vecfuncT result;
	for(auto tmpi:singles.functions){
		CC_function& i=tmpi.second; std::cout << "now doing single " << i.i << std::endl;
		real_function_3d resulti = real_factory_3d(world);
		for(auto tmpk:singles.functions){
			CC_function& k=tmpk.second;	// maybe loop over all pairs
			real_function_3d kgtk = intermediates_.get_pEX(k.i,k.i);
			for(auto tmpl:singles.functions){
				CC_function& l=tmpl.second;

				real_function_6d uil = get_pair_function(u,i.i,l.i); // needed to avoid seg.faults when k>l

				real_function_3d lgtk = intermediates_.get_pEX(l,k);
				real_function_3d l_kgtk = (mo_bra_(l)*kgtk).function;
				real_function_3d k_lgtk = (mo_bra_(k)*lgtk).function;
				real_function_3d part1 = uil.project_out(l_kgtk,1); // 1 is particle 2 and 0 is particle 1
				real_function_3d part2 = uil.project_out(k_lgtk,1);
				real_function_3d part3 = swap_particles(uil).project_out(l_kgtk,0);
				real_function_3d part4 = swap_particles(uil).project_out(k_lgtk,0);
				resulti += (4.0*part1 - 2.0*part2 - 2.0*part3 + part4);
			}
		}
		result.push_back(resulti);
	}
	return result;
}

/// Make the CC2 Residue which is:  Q12f12(T-eij + 2J -K +Un )|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
/// @param[in] \tau_i which will create the |t_i> = |\tau_i>+|i> intermediate
/// @param[in] \tau_j
/// @param[in] u, the uij pair structure which holds the consant part of MP2
/// @param[out] Q12f12(F-eij)|titj> + Q12Ue|titj> - [K,f]|titj>  with |ti> = |\taui>+|i>
/// Right now Calculated in the decomposed form: |titj> = |i,j> + |\taui,\tauj> + |i,\tauj> + |\taui,j>
/// The G_Q_Ue and G_Q_KffK part which act on |ij> are already calculated and stored as constant_term in u (same as for MP2 calculations) -> this should be the biggerst (faster than |titj> form)
real_function_6d CC_Operators::make_cc2_residue(const CC_function &taui, const CC_function &tauj)const{
	const CC_function ti = make_t_intermediate(taui);
	const CC_function tj = make_t_intermediate(tauj);
	const real_function_3d Fti = apply_F(ti)-get_orbital_energies()[ti.i]*ti.function;
	const real_function_3d Ftj = apply_F(tj)-get_orbital_energies()[tj.i]*tj.function;

	// make the Fock operator part:  f(F-eij)|titj> = (F1+F2-ei-ej)|titj> = (F1-ei)|ti>|tj> + |ti>(F2-ei)|tj>
	const real_function_6d fF_titj = make_f_xy(Fti,tj) + make_f_xy(ti,Ftj);


	output("Making the CC2 Residue");

	// make the (U-[K,f])|titj> part
	// first the U Part
	const real_function_6d U_titj = apply_transformed_Ue(ti,tj);
	// then the [K,f] part
	const real_function_6d KffK_titj = apply_exchange_commutator(ti,tj);

	real_function_6d V = (fF_titj + U_titj - KffK_titj);
	V.scale(-2.0);
	V.truncate().reduce_rank();

	fF_titj.print_size(   "CC2-Residue: f12(F-eij)|"+ti.name()+tj.name()+">");
	U_titj.print_size(    "CC2-Residue:          U|"+ti.name()+tj.name()+">");
	KffK_titj.print_size( "CC2-Residue:      [K,f]|"+ti.name()+tj.name()+">");

	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * get_epsilon(taui.i,tauj.i)),parameters.lo, parameters.thresh_bsh_6D);
	G.destructive()=true;
	const real_function_6d GV = G(V);
	return GV;
}

// apply the kinetic energy operator with cusp to a decomposed 6D function
/// @param[in] a 3d function x (will be particle 1 in the decomposed 6d function)
/// @param[in] a 3d function y (will be particle 2 in the decomposed 6d function)
/// @param[out] a 6d function: G(f12*T*|xy>)
real_function_6d CC_Operators::make_GQfT_xy(const real_function_3d &x, const real_function_3d &y, const size_t &i, const size_t &j)const{
	error("make_GQfT should not be used");
	// construct the greens operator
	real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);

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
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified()=true;
		fTxy.fill_tree(screenG).truncate().reduce_rank();
	}{
		real_convolution_6d screenG = BSHOperator<6>(world, sqrt(-2.0*get_epsilon(i,j)),parameters.lo, parameters.thresh_bsh_6D);
		screenG.modified()=true;
		fxTy.fill_tree(screenG).truncate().reduce_rank();
	}
	fTxy_construction_time.info();
	CC_Timer addition_time(world,"f(Tx)y + fxTy");
	real_function_6d result = (fTxy + fxTy).truncate();
	apply_Q12(result,"fT|xy>");

	CC_Timer apply_G(world,"G(fTxy)");
	real_function_6d G_result = G(result);
	G_result.truncate();
	apply_G.info();
	return G_result;
}



/// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
real_function_6d CC_Operators::fock_residue_6d(const CC_Pair &u) const {
	//const double eps = get_epsilon(u.i, u.j);
	// make the coulomb and local Un part with the composite factory
	real_function_3d local_part = (2.0
			* intermediates_.get_hartree_potential()
			+ nemo.nuclear_correlation->U2());
	local_part.print_size("vlocal");
	u.function.print_size("u");

	// Contruct the BSH operator in order to screen

	//	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0 * eps),
	//			parameters.lo, parameters.thresh_bsh_6D);
	//	// apparently the modified_NS form is necessary for the screening procedure
	//	op_mod.modified() = true;
	// Make the CompositeFactory
	real_function_6d vphi =
			CompositeFactory<double, 6, 3>(world).ket(copy(u.function)).V_for_particle1(
					copy(local_part)).V_for_particle2(copy(local_part));
	// Screening procedure
	vphi.fill_tree();

	vphi.print_size("vlocal|u>");

	// the part with the derivative operators: U1
	for (int axis = 0; axis < 6; ++axis) {
		real_derivative_6d D = free_space_derivative<double, 6>(world,axis);
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
		x.fill_tree();
		x.set_thresh(FunctionDefaults<6>::get_thresh());
		vphi += x;
	}

	vphi.print_size("(Un + J1 + J2)|u>");

	// Exchange Part
	vphi = (vphi - K(u.function, u.i == u.j)).truncate().reduce_rank();
	return vphi;

}

/// Echange Operator on 3D function
/// !!!!Prefactor (-1) is not included
real_function_3d CC_Operators::K(const CC_function &x)const{
	return apply_K(x);
}
real_function_3d CC_Operators::K(const real_function_3d &x)const{
	const CC_function tmp(x,99,UNDEFINED);
	return apply_K(tmp);
}

/// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
/// if i==j in uij then the symmetry will be exploited
/// !!!!Prefactor (-1) is not included here!!!!
real_function_6d CC_Operators::K(const real_function_6d &u,
		const bool symmetric) const {
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
	MADNESS_ASSERT(particle == 1 or particle == 2);
	poisson->particle() = particle;
	real_function_6d result = real_factory_6d(world).compressed();
	for (size_t k = 0; k < mo_ket_.size(); k++) {
		real_function_6d X = (multiply(copy(u), copy(mo_bra_(k).function), particle)).truncate();
		real_function_6d Y = (*poisson)(X);
		result += (multiply(copy(Y), copy(mo_ket_(k).function), particle)).truncate();
	}
	return result;
}

// the K operator runs over ALL orbitals (also the frozen ones)
real_function_3d CC_Operators::apply_K(const CC_function &f)const{
	if(parameters.debug and world.rank()==0)  std::cout << "apply K on " << assign_name(f.type) << " function" << std::endl;
	if(parameters.debug and world.rank()==0) std::cout << "K" << f.name() << "=";
	real_function_3d result = real_factory_3d(world);
	switch(f.type){
	case HOLE:
		for(auto k_iterator:mo_ket_.functions){
			const CC_function& k = k_iterator.second;
			const real_function_3d tmp=intermediates_.get_EX(k,f);
			result += (tmp*k.function);
			if(parameters.debug and world.rank()==0) std::cout << "+ <" << k.name() << "|g|" <<f.name() <<">*"<<k.name();
		}
		break;
	case PARTICLE:
		for(auto k_iterator:mo_ket_.functions){
			const CC_function& k = k_iterator.second;
			result += (intermediates_.get_pEX(k,f)*k.function);
			if(parameters.debug and world.rank()==0) std::cout << "+ <" << k.name() << "|g|" <<f.name() <<">*"<<k.name();
		}
		break;
	case MIXED:
		for(auto k_iterator:mo_ket_.functions){
			const CC_function& k = k_iterator.second;
			result += (intermediates_.get_EX(k,f)+intermediates_.get_pEX(k,f))*k.function;
			if(parameters.debug and world.rank()==0) std::cout << "+ <" << k.name() << "|g|t" <<f.i <<">*"<<k.name();
		}
		break;
	default:
		for(auto k_iterator:mo_ket_.functions){
			const CC_function& k = k_iterator.second;
			real_function_3d tmp = ((*poisson)(mo_bra_(k).function*f.function)).truncate();
			result += tmp*k.function;
			if(parameters.debug and world.rank()==0) std::cout << "+ poisson(mo_bra_"<<k.i<<"*"<< f.name() <<")|mo_ket_"<<k.i<<">" << std::endl;
		}
		break;
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
real_function_6d CC_Operators::apply_transformed_Ue(const CC_function &x, const CC_function &y) const {
	output_section("Now applying the transformed Ue operator");
	// make shure the thresh is high enough
	const size_t i = x.i;
	const size_t j = y.i;
	double tight_thresh = guess_thresh(x,y);
	if(tight_thresh > parameters.thresh_Ue){
		tight_thresh = parameters.thresh_Ue;
		output("Thresh to low for Ue setting it back to "+stringify(tight_thresh));
	}
	output("Applying transformed Ue with 6D thresh = " +stringify(tight_thresh));

	real_function_6d Uxy = real_factory_6d(world);
	Uxy.set_thresh(tight_thresh);
	// Apply the untransformed U Potential
	const double eps = get_epsilon(i, j);
	Uxy = corrfac.apply_U(x.function, y.function, eps);
	Uxy.set_thresh(tight_thresh);

	// Get the 6D BSH operator in modified-NS form for screening
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0 * eps),
			parameters.lo,
			parameters.thresh_bsh_6D);
	op_mod.modified() = true;



	// Apply the double commutator R^{-1}[[T,f,R]
	for (size_t axis = 0; axis < 3; axis++) {
		// Make the local parts of the Nuclear and electronic U potentials
		const real_function_3d Un_local = nemo.nuclear_correlation->U1(
				axis);
		const real_function_3d Un_local_x = (Un_local * x.function).truncate();
		const real_function_3d Un_local_y = (Un_local * y.function).truncate();
		const real_function_6d Ue_local = corrfac.U1(axis);
		// Now add the Un_local_x part to the first particle of the Ue_local potential
		real_function_6d UeUnx = CompositeFactory<double, 6, 3>(world).g12(
				Ue_local).particle1(Un_local_x).particle2(copy(y.function)).thresh(
						tight_thresh);
		// Fill the Tree were it will be necessary
		UeUnx.fill_tree(op_mod);
		// Set back the thresh
		UeUnx.set_thresh(FunctionDefaults<6>::get_thresh());

		UeUnx.print_size("UeUnx");

		// Now add the Un_local_y part to the second particle of the Ue_local potential
		real_function_6d UeUny = CompositeFactory<double, 6, 3>(world).g12(
				Ue_local).particle1(copy(x.function)).particle2(Un_local_y).thresh(
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
	real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(copy(x.function*nemo.nuclear_correlation -> square())).particle2(copy(y.function*nemo.nuclear_correlation -> square()));

	const double a = inner(Uxy, tmp);
	const real_function_3d xx = (x.function * x.function*nemo.nuclear_correlation -> square());
	const real_function_3d yy = (y.function * y.function*nemo.nuclear_correlation -> square());
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
	Uxy.print_size("Result of applying Ue");
	return Uxy;
}

/// Apply the Exchange Commutator [K,f]|xy>
real_function_6d CC_Operators::apply_exchange_commutator(const CC_function &x, const CC_function &y)const{

	output_section("Now applying [K,f]");

	// make first part of commutator
	CC_Timer part1_time(world,"Kf");
	real_function_6d Kfxy = apply_Kf(x,y);
	part1_time.info();
	// make the second part of the commutator
	CC_Timer part2_time(world,"fK");
	real_function_6d fKxy = apply_fK(x,y).truncate();
	part2_time.info();
	real_function_6d result = (Kfxy - fKxy);

	// sanity check
	// <psi|[A,B]|psi> = <psi|AB|psi> - <psi|BA|psi> = <Apsi|Bpsi> - <Bpsi|Apsi> = 0 (if A,B hermitian)
	{
		CC_Timer sanity(world,"[K,f] sanity check");
		output("\n\n[K,f] Sanity Check");
		// make the <xy| bra state which is <xy|R2
		const real_function_3d brax = x.function * nemo.nuclear_correlation->square();
		const real_function_3d bray = y.function * nemo.nuclear_correlation->square();
		const real_function_6d xy = make_xy(CC_function(brax,x.i,x.type),CC_function(bray,y.i,y.type));
		const double xyfKxy = xy.inner(fKxy);
		const double xyKfxy = xy.inner(Kfxy);
		const double diff = xyfKxy - xyKfxy;
		if(world.rank()==0){
			std::cout << std::setprecision(parameters.output_prec);
			std::cout << "<" << x.name() << y.name() << "|fK|" << x.name() << y.name() << "> =" << xyfKxy << std::endl;
			std::cout << "<" << x.name() << y.name() << "|Kf|" << x.name() << y.name() << "> =" << xyKfxy << std::endl;
			std::cout << "difference = " << diff << std::endl;
		}
		if(fabs(diff)>FunctionDefaults<6>::get_thresh()) warning("Exchange Commutator Plain Wrong");
		else if(fabs(diff)>FunctionDefaults<6>::get_thresh()*0.1) warning("Exchange Commutator critical");
		output("[K,f] Sanity Check ended\n\n");

		sanity.info();

	}
	result.print_size("[K,f]|"+x.name()+y.name()+">");
	return result;
}

/// Apply the Exchange operator on a tensor product multiplied with f12
/// !!! Prefactor of (-1) is not inclued in K here !!!!
real_function_6d CC_Operators::apply_Kf(const CC_function &x, const CC_function &y) const {
	bool symmetric = false;
	if((x.type == y.type) and (x.i == y.i)) symmetric = true;

	// First make the 6D function f12|x,y>
	real_function_6d f12xy = make_f_xy(x,y);
	// Apply the Exchange Operator
	real_function_6d result = K(f12xy, symmetric);
	return result;
}

/// Apply fK on a tensor product of two 3D functions
/// fK|xy> = fK_1|xy> + fK_2|xy>
/// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
/// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
real_function_6d CC_Operators::apply_fK(const CC_function &x, const CC_function &y) const {

	const real_function_3d Kx = K(x);
	const real_function_3d Ky = K(y);
	const real_function_6d fKphi0a = make_f_xy(x,CC_function(Ky,y.i,UNDEFINED));
	const real_function_6d fKphi0b = make_f_xy(CC_function(Kx,x.i,UNDEFINED),y);
	const real_function_6d fKphi0 = (fKphi0a + fKphi0b);
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
	if(taui.i!=tauj.i)warning("ccs energy fock parts only defined for one orbital molecules");
	double omega = 0.0;
	// fock operator parts (zero when HF converged)
	double omega_f = 2.0*mo_bra_(taui.i).inner(apply_F(CC_function(taui.function,taui.i,UNDEFINED)));
	output("CCS Energy Fock part: 2.0*<i|F|taui>="+stringify(omega_f));
	//omega += 2.0*intermediates_.get_integrals_t1()(u.i,u.j,u.i,u.j); //<ij|g|\taui\tauj>
	omega += 2.0*make_ijgxy(taui.i,tauj.i,taui.function,tauj.function);
	//omega -= intermediates_.get_integrals_t1()(u.i,u.j,u.j,u.i);     //<ij|g|\tauj\taui>
	omega -= make_ijgxy(taui.i,tauj.i,tauj.function,taui.function);
	output("CCS Energy Coulomb part: 2.0<ij|g|\taui\tauj> - <ji|g|\taui\tauj>="+stringify(omega));
	return omega+omega_f;
}
double CC_Operators::compute_cc2_pair_energy(const CC_Pair &u,
		const CC_function &taui, const CC_function &tauj) const {
	double omega = 0.0;
	const size_t i = u.i;
	const size_t j = u.j;
	MADNESS_ASSERT(i==taui.i);
	MADNESS_ASSERT(j==tauj.i);
	double u_part = 0.0;
	double ij_part = 0.0;
	double mixed_part = 0.0;
	double titj_part = 0.0;
	double singles_part = 0.0;
	double tight_thresh = parameters.thresh_Ue;
	// Contribution from u itself, we will calculate <uij|g|ij> instead of <ij|g|uij> and then just make the inner product (see also mp2.cc)
	{
		real_function_6d coulomb = TwoElectronFactory(world).dcut(tight_thresh);
		real_function_6d g_ij = CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(i).function)).particle2(copy(mo_bra_(j).function)).g12(coulomb).thresh(tight_thresh);

		real_function_6d g_ji = CompositeFactory<double, 6, 3>(world).particle1(copy(mo_bra_(j).function)).particle2(copy(mo_bra_(i).function)).g12(coulomb).thresh(tight_thresh);
		const double uij_g_ij = inner(u.function, g_ij);
		const double uij_g_ji = inner(u.function, g_ji); // =uji_g_ij
		u_part = 2.0 * uij_g_ij - uij_g_ji;
	}
	// Contribution from the mixed f12(|\tau_i,j>+|i,\tau_j>) part
	{
		mixed_part += 2.0*make_ijgQfxy(u.i,u.j,mo_ket_(i),tauj);
		mixed_part += 2.0*make_ijgQfxy(u.i,u.j,taui,mo_ket_(j));
		mixed_part -= make_ijgQfxy(u.j,u.i,mo_ket_(i),tauj);
		mixed_part -= make_ijgQfxy(u.j,u.i,taui.function,mo_ket_(j));
	}
	// Contribution from the f12|ij> part, this should be calculated in the beginning
	{
		ij_part += (2.0*u.ij_gQf_ij - u.ji_gQf_ij );
	}
	// Contribution from the f12|\tau_i\tau_j> part
	{
		titj_part += 2.0*make_ijgQfxy(u.i,u.j,taui,tauj);
		titj_part -= make_ijgQfxy(u.i,u.j,tauj,taui);
	}
	// Singles Contribution
	{
		// I should use intermediates later because the t1 integrals are also needed for the CC2 potential
		//omega += 2.0*intermediates_.get_integrals_t1()(u.i,u.j,u.i,u.j); //<ij|g|\taui\tauj>
		singles_part += 2.0*make_ijgxy(u.i,u.j,taui.function,tauj.function);
		//omega -= intermediates_.get_integrals_t1()(u.i,u.j,u.j,u.i);     //<ij|g|\tauj\taui>
		singles_part -= make_ijgxy(u.i,u.j,tauj.function,taui.function);
	}

	omega = u_part + ij_part + mixed_part + titj_part + singles_part;
	if(world.rank()==0){
		std::cout << "\n\nEnergy Contributions to the correlation energy of pair " << i << j << "\n";
		std::cout << std::fixed << std::setprecision(parameters.output_prec);
		std::cout << "from   |u" << i << j << "            |: " << u_part << "\n";
		std::cout << "from Qf|HH" << i << j << "           |: " << ij_part << "\n";
		std::cout << "from Qf|HP" << i << j << "           |: " << mixed_part << "\n";
		std::cout << "from Qf|PPu" << i << j << "          |: " << titj_part << "\n";
		std::cout << "from   |tau" << i << ",tau" << j << "|: " << singles_part << "\n";
		std::cout << "all together = " << omega << "\n";
		std::cout << "\n\n";
	}
	return omega;
}

/// General Function to make the intergral <ij|gQf|xy>
double CC_Operators::make_ijgQfxy(const size_t &i, const size_t &j, const CC_function &x, const CC_function &y)const{
	// convenience
	const real_function_3d brai = mo_bra_(i).function;
	const real_function_3d braj = mo_bra_(j).function;
	// part 1, no projector: <ij|gf|xy>
	const real_function_3d jy = (braj*y.function).truncate();
	const real_function_3d ix = (brai*x.function).truncate();
	const real_function_3d jgfy = apply_gf(jy);
	const double part1 = ix.inner(jgfy);
	// part 2, projector on particle 1 <j|igm*mfx|y> = jy.inner(igm*mfx)
	double part2 = 0.0;
	for(const auto& mtmp:mo_ket_.functions){
		const CC_function& mom = mtmp.second;
		const size_t m = mtmp.first;
		const real_function_3d igm = apply_g12(mo_bra_(i),mom);
		const real_function_3d mfx = apply_f12(mo_bra_(m),x);
		part2 -= jy.inner(igm*mfx);
	}
	// part3, projector on particle 2 <i|jgn*nfy|x>
	double part3 = 0.0;
	for(const auto& ntmp:mo_ket_.functions){
		const CC_function& mon = ntmp.second;
		const size_t n = ntmp.first;
		const real_function_3d jgn = apply_g12(mo_bra_(j),mon);
		const real_function_3d nfy = apply_f12(mo_bra_(n),y);
		part2 -= ix.inner(jgn*nfy);
	}
	// part4, projector on both particles <ij|g|mn><mn|f|xy>
	double part4 = 0.0;
	for(const auto& mtmp:mo_ket_.functions){
		const CC_function& mom = mtmp.second;
		const size_t m = mtmp.first;
		const real_function_3d igm = apply_g12(mo_bra_(i),mom);
		const real_function_3d mfx = apply_f12(mo_bra_(m),x);
		for(const auto& ntmp:mo_ket_.functions){
			const CC_function& mon = ntmp.second;
			const size_t n = ntmp.first;
			const real_function_3d jn = braj*mon.function;
			const real_function_3d ny = mo_bra_(n).function*y.function;
			const double ijgmn = jn.inner(igm);
			const double mnfxy = ny.inner(mfx);
			part4 += ijgmn*mnfxy;
		}
	}

	return part1+part2+part3+part4;

	//	// Q12 = I12 - O1 - O2 + O12
	//	real_function_3d jy = mo_bra_(j).function*y;
	//	real_function_3d ix = mo_bra_(i).function*x;
	//	// I12 Part:
	//	double ijgfxy = (ix).inner(apply_gf(jy));
	//	// O1 Part
	//	double ijgO1fxy =0.0;
	//	for(size_t k=0;k<mo_ket_.size();k++){
	//		real_function_3d igk = intermediates_.get_EX(i,k);
	//		real_function_3d kfx = (*f12op)(mo_bra_(k).function*x);
	//		real_function_3d igkkfx = (igk*kfx).truncate();
	//		ijgO1fxy += jy.inner(igkkfx);
	//	}
	//	// O2 Part
	//	double ijgO2fxy =0.0;
	//	for(size_t k=0;k<mo_ket_.size();k++){
	//		real_function_3d jgk = intermediates_.get_EX(j,k);
	//		real_function_3d kfy = (*f12op)(mo_bra_(k).function*y);
	//		real_function_3d jgkkfy = (jgk*kfy).truncate();
	//		ijgO2fxy += ix.inner(jgkkfy);
	//	}
	//	// O12 Part
	//	double ijgO12fxy = 0.0;
	//	for(size_t k=0;k<mo_ket_.size();k++){
	//		real_function_3d igk = intermediates_.get_EX(i,k);
	//		real_function_3d kfx = (*f12op)(mo_bra_(k).function*x);
	//		for(size_t l=0;l<mo_ket_.size();l++){
	//			double ijgkl = igk.inner(mo_bra_(j).function*mo_ket_(l).function);
	//			double klfxy = kfx.inner(mo_bra_(l).function*y);
	//			ijgO12fxy += ijgkl*klfxy;
	//		}
	//	}
	//
	//	return (ijgfxy - ijgO1fxy - ijgO2fxy + ijgO12fxy);
}

double CC_Operators::make_ijgfxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const{
	real_function_3d jy = mo_bra_(j).function*y;
	real_function_3d ix = mo_bra_(j).function*x;
	// I12 Part:
	double ijgfxy = (ix).inner(apply_gf(jy));
	return ijgfxy;
}

/// General Function to make the two electron integral <ij|g|xy>
/// For Debugging -> Expensive without intermediates
double CC_Operators::make_ijgxy(const size_t &i, const size_t &j, const real_function_3d &x, const real_function_3d &y)const{
	real_function_3d igx = (*poisson)(mo_bra_(i).function*x).truncate();
	real_function_3d jy = (mo_bra_(j).function*y).truncate();
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
					copy(mo_bra_(i).function)).particle2(
							copy(mo_bra_(j).function)).g12(g);

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
					copy(mo_bra_(i).function)).particle2(
							copy(mo_bra_(j).function)).g12(G);

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

/// Calculation is done in 4 steps over: Q12 = 1 - O1 - O2 + O12
/// 1. <x|f12|z>*|y>
/// 2. -\sum_m <x|m> <m|f12|z>*|y> = -(\sum_m <x|m> <m|f12|z>) *|y>
/// 3. -\sum_n <nx|f12|zy> * |n>
/// 4. +\sum_{mn} <x|n> <mn|f12|yz> * |m>
/// Metric from nuclear cusp is not taken into account -> give the right bra elements to the function
//real_function_3d CC_Operators::convolute_x_Qf_yz(const real_function_3d &x, const real_function_3d &y, const real_function_3d &z)const{
//	// make intermediates
//	vecfuncT moz = mul(world,z,mo_bra_);
//	vecfuncT moy = mul(world,y,mo_bra_);
//	// Do <x|f12|z>*|y>
//	real_function_3d part1_tmp = (*f12op)(x*z);
//	real_function_3d part1 = (part1_tmp*y).truncate();
//	// Do -\sum_m <x|m> <m|f12|z>*|y>
//	real_function_3d part2 = real_factory_3d(world);
//	{
//		Tensor<double> xm = inner(world,x,mo_ket_);
//		vecfuncT f12mz = apply(world,*f12op,moz);
//		real_function_3d tmp = real_factory_3d(world);
//		for(size_t m=0;m<mo_bra_.size();m++) tmp += xm[m]*f12mz[m];
//		part2 = tmp*y;
//	}
//	// Do -\sum_n <nx|f12|zy> * |n> |  <nx|f12|zy> = <n| xf12y |z>
//	real_function_3d part3 = real_factory_3d(world);
//	{
//		real_function_3d xf12y = (*f12op)((x*y).truncate());
//		for(size_t n=0;n<mo_bra_.size();n++){
//			double nxfzy = xf12y.inner(mo_bra_[n]*z);
//			part3 += nxfzy*mo_ket_[n];
//		}
//	}
//	// Do +\sum_{mn} <x|n> <mn|f12|yz> * |m>
//	real_function_3d part4 = real_factory_3d(world);
//	{
//		Tensor<double> xn = inner(world,x,mo_ket_);
//		vecfuncT nf12z = apply(world,*f12op,moz);
//		Tensor<double> mnfyz = inner(world,moy,nf12z);
//		for(size_t m=0;m<mo_bra_.size();m++){
//			for(size_t n=0;n<mo_ket_.size();n++){
//				part4 += xn(n)*mnfyz(m,n)*mo_ket_[m];
//			}
//		}
//	}
//	real_function_3d result = part1 - part2 - part3 + part4;
//	result.truncate();
//
//	if(parameters.debug){
//		CC_Timer function_debug(world,"Debug-Time for <k|Qf|xy>");
//		real_function_6d test_tmp = CompositeFactory<double,6,3>(world).g12(corrfac.f()).particle1(copy(y)).particle2(copy(z));
//		test_tmp.fill_tree().truncate().reduce_rank();
//		real_function_6d test_1 = copy(test_tmp);
//		real_function_6d test_4 = copy(Q12(test_tmp));
//		real_function_3d test1 = test_1.project_out(x,1);
//		real_function_3d test4 = test_4.project_out(x,1);
//		double check1 = (test1 - part1).norm2();
//		double check4 = (test4 - result).norm2();
//		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check1 value is: " << check1 << std::endl;
//		if(world.rank()==0) std::cout << std::setprecision(parameters.output_prec) << "<k|Qf|xy> debug, difference to check4 value is: " << check4 << std::endl;
//		if(check1 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check1 failed");
//		if(check4 > FunctionDefaults<6>::get_thresh()) warning("<k|Qf|xy> check4 failed");
//		function_debug.info();
//	}
//
//	return result;
//}






}
