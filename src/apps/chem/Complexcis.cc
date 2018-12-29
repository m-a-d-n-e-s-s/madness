/*
 * Complexcis.cpp
 *
 *  Created on: 21 Nov 2018
 *      Author: fbischoff
 */

#include <chem/Complexcis.h>
#include <chem/GuessFactory.h>
#include <chem/SCFOperators.h>


namespace madness {

double Complex_cis::value() {
	if (cis_param.swap_ab()) {
		std::swap(nemo.aeps,nemo.beps);
		std::swap(nemo.amo,nemo.bmo);
		std::swap(Qa,Qb);
	}
	std::vector<root> roots=make_guess();
	iterate(roots);
	return 0.0;
}

/// iterate the roots
void Complex_cis::iterate(std::vector<root>& roots) const {

	// some numbers
	const int anoct=noct(nemo.aeps).size();
	const int bnoct=noct(nemo.beps).size();
	const int noct=anoct+bnoct;

	// timings
	double wall0=wall_time(), wall1=wall_time();

	// compute ground state densities
	const real_function_3d adens=real(dot(world,conj(world,nemo.amo),nemo.amo)).truncate();
	const real_function_3d bdens=real(dot(world,conj(world,nemo.bmo),nemo.bmo)).truncate();
	real_function_3d totdens=(adens+bdens).truncate();
	if (nemo.cparam.spin_restricted) totdens.scale(2.0);
	double nelectron=totdens.trace();
	print("nelectron from the total density",nelectron);
	totdens.print_size("totdens");

	for (int iter=0; iter<cis_param.maxiter(); ++iter) {
		print("iteration ",iter);
		print("root     energy         wf delta    energy change    elapsed time");
		for (int iroot=0; iroot<cis_param.guess_excitations(); ++iroot) {
			root& thisroot=roots[iroot];

			std::vector<complex_function_3d> td_a=mul(world,thisroot.afunction,active_mo(nemo.amo));
			std::vector<complex_function_3d> td_b=mul(world,thisroot.bfunction,active_mo(nemo.bmo));
			int i=0;
			for (auto& a : td_a) save(real(a),"transition_density_a"+stringify(i++));

			wall1=wall_time();
			printf("  %2d   %12.8f     %4.2e      %4.2e      %4.1fs\n",iroot, thisroot.omega, thisroot.delta,
					thisroot.energy_change,wall1-wall0);
			wall0=wall1;

			std::vector<complex_function_3d> apot=zero_functions_compressed<double_complex,3>(world,anoct);
			std::vector<complex_function_3d> bpot=zero_functions_compressed<double_complex,3>(world,bnoct);

			// compute the perturbed density
			real_function_3d denspt=real((dot(world,conj(world,thisroot.afunction),active_mo(nemo.amo)) +
					dot(world,conj(world,thisroot.bfunction),active_mo(nemo.bmo))));
			if (nemo.cparam.spin_restricted) denspt.scale(2.0);

			denspt.truncate();
			double bla=denspt.trace();
			print("trace over the perturbed density",bla);

			denspt.print_size("denspt");

			for (std::string spin : {"alpha","beta"}) {
				if (nemo.cparam.spin_restricted and (spin=="beta")) continue;

				const std::vector<complex_function_3d>& mo=(spin=="alpha") ? nemo.amo : nemo.bmo;
				if (mo.size()==0) continue;

				std::vector<complex_function_3d>& pot=(spin=="alpha") ? apot : bpot;
				std::vector<complex_function_3d> x= (spin=="alpha") ? thisroot.afunction : thisroot.bfunction;
				const QProjector<double_complex,3>& Q=(spin=="alpha") ? Qa : Qb;
				Tensor<double> occ = (spin=="alpha") ? Tensor<double>(nemo.amo.size()) : Tensor<double>(nemo.bmo.size());
				occ=1.0;

				/// zeroth order Fock operator acting on the x functions
				std::vector<complex_function_3d> vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo;
				nemo.compute_potentials(mo,totdens,x,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);
				if (spin=="beta") scale(world,spin_zeeman_nemo,-1.0);

				pot+=(vnemo+vlznemo+dianemo+spin_zeeman_nemo-knemo+jnemo);	// need += b/c of reference
				truncate(world,pot);

				// perturbed Fock operator acting on the reference orbitals
				Coulomb Jp(world);
				Jp.potential() = Jp.compute_potential(denspt);

				Exchange<double_complex,3> Kp(world);
				Kp.set_parameters(conj(world,mo),x,occ);
				pot+=Q(Jp(mo) - Kp(mo));
				truncate(world,pot);
			}

			std::vector<complex_function_3d> pot=append(apot,bpot);
			std::vector<complex_function_3d> residuals=compute_residuals(pot,thisroot);
			thisroot.delta=norm2(world,residuals);
			if (cis_param.omega()!=0) {
				thisroot.omega=cis_param.omega();
				printf("\n\nset excitation energy manually to %8.4f\n\n",thisroot.omega);
			}


			auto [ares, bres] = split(residuals,thisroot.afunction.size());
//			print_size(world,ares,"ares");
//			print_size(world,bres,"bres");

			if (ares.size()>0) thisroot.afunction-=ares;
			if (bres.size()>0) thisroot.bfunction-=bres;
			print_size(world,thisroot.afunction,"afunction");
			print_size(world,thisroot.bfunction,"bfunction");

		}
		normalize(roots);

	}

}

Tensor<double> Complex_cis::compute_expectation_values(const std::vector<root>& roots,
		const std::vector<complex_function_3d>& mo, const real_function_3d& density) const {

	int thisroot=0;

	std::vector<complex_function_3d> vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo;
	std::vector<complex_function_3d> arg=append(roots[thisroot].afunction,roots[thisroot].bfunction);
	nemo.compute_potentials(mo,density,arg,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);
	Tensor<double_complex> fockpt=nemo.compute_vmat(arg,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);
	Kinetic<double_complex,3> T(world);
	fockpt+=T(arg,arg);

}


std::vector<complex_function_3d> Complex_cis::compute_residuals(std::vector<complex_function_3d>& pot,
		root& root) const {

	// make the BSH operator
	std::vector<std::shared_ptr<SeparatedConvolution<double,3> > > bsh(pot.size());

	std::vector<complex_function_3d> x=append(root.afunction,root.bfunction);
	Tensor<double> eps=copy(concatenate(noct(nemo.aeps),noct(nemo.beps))) + root.omega;
	for (std::size_t p = 0; p < bsh.size(); p++) {

		// if eps is above zero we have an unbound state (or are early in the iteration) however this needs a shift of the potential
		// we shift to -0.05 (same as moldft does, no particular reason)
		if(eps[p]>0.0){
			print("potential shift needed for V" + std::to_string(p+cis_param.freeze()));
			double shift = eps[p]+0.05;
			eps[p] -= shift;
			pot[p] -= (-2.0*shift*x[p]);
		}
		MADNESS_ASSERT(not(eps[p]>0.0));
//		print("assigning eps to bsh ",p,eps[p]);
		bsh[p] = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world, sqrt(-2.0 * eps[p]), 1.e-4, cis_param.thresh()));
	}

	// apply the BSH operator
	std::vector<complex_function_3d> tmp0=(apply(world,bsh,-2.0*pot));

	auto [atmp, btmp] = split(tmp0,root.afunction.size());
	std::vector<complex_function_3d> tmp=append(Qa(atmp),Qb(btmp));

	std::vector<complex_function_3d> residual=x-tmp;
	truncate(world,residual);

	// Calculate Second Order Energy Update
	double tmp1 = real(inner(world,conj(world,residual),pot).sum());
	// squared norm of GVpsi (Psi_tilde)
	double tmp2 = real(inner(world,conj(world,tmp),tmp).sum());

	const double sou= -tmp1/tmp2;
	root.omega+=sou;
	root.energy_change=sou;

	return residual;

}

void Complex_cis::update_roots(std::vector<root>& aroot, std::vector<root>& broot, std::vector<root>& troot) const {

}


std::vector<Complex_cis::root> Complex_cis::read_guess(const std::string spin) const {
	std::vector<root> guess(cis_param.guess_excitations());
	MADNESS_EXCEPTION("implement Complex_cis::read_guess()",1);
	return guess;
}


std::vector<Complex_cis::root> Complex_cis::make_guess() const {

	std::vector<complex_function_3d> avirtuals, bvirtuals;
	Tensor<double> aveps,bveps;		// virtual orbital energies for alpha and beta
	real_function_3d adens=real(dot(world,conj(world,active_mo(nemo.amo)),active_mo(nemo.amo))).truncate();
	real_function_3d bdens=real(dot(world,conj(world,active_mo(nemo.bmo)),active_mo(nemo.bmo))).truncate();
	real_function_3d density=adens + bdens;
	if (nemo.cparam.spin_restricted) density.scale(2.0);

	density.print_size("density in make_guess");

	std::vector<root> guess(cis_param.guess_excitations());

	for (std::string spin : {"alpha","beta"}) {
		print("spin " ,spin);
		std::vector<complex_function_3d> virtuals;
		std::vector<complex_function_3d> active_mo= (spin=="alpha") ? nemo.amo : nemo.bmo;
		if (active_mo.size()==0) continue;

		// prepare the list of excitation operators and copied seeds
		std::vector<coord_3d> centers = guessfactory::compute_centroids(active_mo);
		std::vector<std::pair<std::vector<complex_function_3d>, std::string> > exlist;
		{
	//		std::vector<std::string> exop_strings=cis_param.exops();
			std::vector<std::string> exop_strings=(guessfactory::make_predefined_exop_strings(cis_param.guess_excitation_operators()));
			for(const auto ex: exop_strings){
				std::vector<complex_function_3d> cseed=copy(world,active_mo,false);
				exlist.push_back(std::make_pair(cseed,ex));
			}
			print(exop_strings);
		}
		world.gop.fence();
		std::cout << "will create " << exlist.size()*centers.size() << " virtuals, from " << centers.size()
				<< " seeds and " << exlist.size() << " excitation operators"   << std::endl;

		// create the virtuals by unary operations: multiply excitation operators with seeds
		for(auto it:exlist){
			virtuals=append(virtuals,guessfactory::apply_trigonometric_exop(it.first,it.second,centers,false));
		}
		world.gop.fence();


		if (spin=="alpha") virtuals=Qa(virtuals);
		else if (spin=="beta") virtuals=Qb(virtuals);

		// remove linear dependencies
		const size_t spre=virtuals.size();
		virtuals=orthonormalize_canonical(virtuals,1.e-6);
		if(virtuals.size()!=spre) std::cout << "removed " << spre-virtuals.size() << " virtuals due to linear dependencies" << std::endl;

		// canonicalize virtuals and set up the CIS matrix
		if (spin=="alpha") {
			canonicalize(active_mo,density,virtuals,aveps);
			avirtuals=virtuals;
		}
		if (spin=="beta") {
			canonicalize(active_mo,density,virtuals,bveps);
			bvirtuals=virtuals;
		}

	}

	// active orbital energies only
	int anoct=noct(nemo.aeps).size();
	int bnoct=noct(nemo.beps).size();
	Tensor<double> oeps=concatenate(noct(nemo.aeps),noct(nemo.beps));
	Tensor<double> veps=concatenate(aveps,bveps);
//	auto [ae,be] = split(veps,nemo.aeps.size());

	Tensor<double_complex> MCISa=make_CIS_matrix(aveps,noct(nemo.aeps));
	Tensor<double_complex> MCISb=make_CIS_matrix(bveps,noct(nemo.beps));
	Tensor<double_complex> MCIS(MCISa.dim(0)+MCISb.dim(0),MCISa.dim(1)+MCISb.dim(1));
	if (MCISa.size()>0) MCIS(Slice(0,MCISa.dim(0)-1,1),Slice(0,MCISa.dim(0)-1,1))=MCISa;
	if (MCISb.size()>0) MCIS(Slice(MCISa.dim(0),-1,1),Slice(MCISa.dim(0),-1,1))=MCISb;

	// compute matrix CIS excitation energies
	Tensor<double_complex> U,UU;
	Tensor<double> evals;
	syev(MCIS, U, evals);
	print("MCIS excitation energies ",evals);
//	print(real(U));
//	print("MCIS excitation vectors");
//	for (int i=0; i<cis_param.guess_excitations(); ++i) print(real(U(_,i)));


	// assemble the guess vectors

	std::vector<complex_function_3d> abvirtuals=append(avirtuals,bvirtuals);
	int nvirt = abvirtuals.size();
	int nva=avirtuals.size();
	int nvb=bvirtuals.size();
	auto get_vir_idx_a = [nva](int I) {return I%nva;};
	auto get_occ_idx_a = [nva](int I) {return I/nva;};
	auto get_vir_idx_b = [nvb](int I) {return I%nvb;};
	auto get_occ_idx_b = [nvb](int I) {return I/nvb;};

//	print_size(world,avirtuals,"avirtuals");
//	print_size(world,bvirtuals,"bvirtuals");

	// split U into alpha and beta parts
	Tensor<double_complex> Ua,Ub;
	if (nva>0) Ua=U(Slice(0,MCISa.dim(0)-1),_);	// alpha part for all excitations
	if (nvb>0) Ub=U(Slice(MCISa.dim(0),-1),_);	// beta part for all excitations
//	print("Ua, Ub");
//	print(Ua);
//	print(Ub);


	for (size_t I = 0, iexcitation=0; I < MCIS.dim(1); ++I, iexcitation++) {
		if (iexcitation >= cis_param.guess_excitations()) break;

		guess[iexcitation].afunction=zero_functions_compressed<double_complex,3>(world,anoct);
		guess[iexcitation].bfunction=zero_functions_compressed<double_complex,3>(world,bnoct);

		if (evals(I) < 0.0 && world.rank() == 0)
			MADNESS_EXCEPTION("NEGATIVE EIGENVALUE IN INITIAL DIAGONALIZATION: CHECK YOUR REFERENCE!\n",0);

		guess[iexcitation].omega = evals(I);
		guess[iexcitation].excitation = iexcitation;

		// alpha part
		if (Ua.size()>0) {
			for (size_t J = 0; J < Ua.dim(0); ++J) {
				const int b = get_vir_idx_a(J);
				const int j = get_occ_idx_a(J);
//				print("I, J, b, j",I, J, b, j);
				const double_complex xjb = Ua(J, I);
				guess[iexcitation].afunction[j] += xjb * avirtuals[b];
			}
		}

		// beta part
		if (Ub.size()>0) {
			for (size_t J = 0; J < Ub.dim(0); ++J) {
				const int b = get_vir_idx_b(J);
				const int j = get_occ_idx_b(J);
//				print("I, J, b, j",I, J, b, j);
				const double_complex xjb = Ub(J, I);
				guess[iexcitation].bfunction[j] += xjb * bvirtuals[b];
			}
		}
		print_size(world,guess[iexcitation].afunction,"afunction");
		print_size(world,guess[iexcitation].bfunction,"bfunction");
	}
	for (auto& x : guess) truncate(world, x.afunction, cis_param.thresh());
	for (auto& x : guess) truncate(world, x.bfunction, cis_param.thresh());
	normalize(guess);

	return guess;

}

void Complex_cis::canonicalize(const std::vector<complex_function_3d>& mo, const real_function_3d& density,
		std::vector<complex_function_3d>& virtuals, Tensor<double>& veps) const {

	std::vector<complex_function_3d> vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo;
	nemo.compute_potentials(mo,density,virtuals,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);
	Tensor<double_complex> fock=nemo.compute_vmat(virtuals,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);
	Kinetic<double_complex,3> T(world);
	fock+=T(virtuals,virtuals);
//	print("virtual Fock");
//	print(fock);

	// get the fock transformation matrix
	Tensor<double_complex> U;
	Tensor<double_complex> overlap=matrix_inner(world,virtuals,virtuals);
    sygvp(world, fock, overlap, 1, U, veps);
//    print("virtual orbital energies",veps);
    virtuals = madness::transform(world, virtuals, U);


}

Tensor<double_complex> Complex_cis::make_CIS_matrix(const Tensor<double>& veps, const Tensor<double>& act) const {

	// the cis matrix is indexed by ij and ab
	// we will use the combined indixes from ia and jb named I and J
	// in order to not be confused we use the following helper functions
	const int nact = act.size();
	// determines for which orbitals (couting from the HOMO downwards) the off-diagonal elements will be computed
	// this simplifies the guess
	const int nvirt = veps.size();

	const int dim=(nvirt*nact);
	std::cout << "CIS-Matrix for guess calculation will be of size " << dim << "x" << dim <<std::endl;
	// the number of the matrix where elements which are not determined by orbital energies and the fock matrix are computed (controlled over active_guess_orbitals parameter)
	Tensor<double> MCIS(dim,dim);

	auto get_vir_idx = [nvirt](int I) {return I%nvirt;};
	auto get_occ_idx = [nvirt](int I) {return I/nvirt;};

	for(int I=0;I<dim;++I){
		const int a=get_vir_idx(I);
		const int i=get_occ_idx(I);
		MCIS(I,I) = veps(a)-act(i);
	}
	return MCIS;

}


void Complex_cis::normalize(std::vector<root>& roots) const {

	MADNESS_ASSERT(roots.size()==1);
	for (auto root : roots) {
		double na=0.0, nb=0.0;
        if (root.afunction.size()>0) na = norm2(world, root.afunction);
        if (root.bfunction.size()>0) nb = norm2(world, root.bfunction);

        if (na>1.e-8) scale(world,root.afunction,1.0/na);
        if (nb>1.e-8) scale(world,root.bfunction,1.0/nb);
	}
}



} /* namespace madness */
