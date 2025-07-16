/*
 * Complexcis.cpp
 *
 *  Created on: 21 Nov 2018
 *      Author: fbischoff
 */

#include<madness/chem/GuessFactory.h>
#include<madness/chem/SCFOperators.h>
#include "zcis.h"

namespace madness {

double Zcis::value() {
	if (cis_param.swap_ab()) {
		std::swap(nemo->aeps,nemo->beps);
		std::swap(nemo->amo,nemo->bmo);
		std::swap(Qa,Qb);
	}
	std::vector<root> roots;
	try {
		roots=read_guess();
	} catch (...) {
		roots=make_guess();
	}
	iterate(roots);
	save_guess(roots);
	return 0.0;
}

/// iterate the roots
void Zcis::iterate(std::vector<root>& roots) const {

	// timings
	double wall0=wall_time(), wall1=wall_time();

	// compute ground state densities
	const real_function_3d adens=real(dot(world,conj(world,nemo->amo),nemo->amo)).truncate();
	const real_function_3d bdens=real(dot(world,conj(world,nemo->bmo),nemo->bmo)).truncate();
	real_function_3d totdens=(adens+bdens).truncate();
	if (nemo->get_calc_param().spin_restricted()) totdens.scale(2.0);
	double nelectron=totdens.trace();
	print("nelectron from the total density",nelectron);
	totdens.print_size("totdens");

	//const double shift=nemo->param.shift();
	const bool use_kain=true;

	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, vector_function_allocator<double_complex,3> >
			allsolver(vector_function_allocator<double_complex,3> (world,(active_mo(nemo->amo).size()+active_mo(nemo->bmo).size())*roots.size()));

	for (int iter=0; iter<cis_param.maxiter(); ++iter) {
		wall1=wall_time();
		printf("iteration %2d  %4.1fs\n",iter,wall1-wall0);
		wall0=wall1;
		print("root     energy         wf delta      energy change ");

		for (size_t i=0; i<roots.size(); ++i)
			printf("  %2lu   %12.8f     %4.2e      %5.2e\n",i, roots[i].omega, roots[i].delta, roots[i].energy_change);

		bool do_kain=use_kain and (iter>cis_param.guess_maxiter());
		if (not do_kain) allsolver.clear_subspace();	// fock_pt diag reshuffles the roots and confuses solver

		// compute the perturbed fock matrix and the excitation energy expectation values
		compute_potentials(roots, totdens);

		Tensor<double> omega;
		if (not do_kain) {
			Tensor<double_complex> fock_pt = compute_fock_pt(roots);

			// compute the expectation value of the excitation energy Eq. (34) of Kottmann2015
			Tensor<double_complex> ovlp(roots.size()),eovlp(roots.size());
			for (std::size_t i=0; i<roots.size(); ++i) {
				ovlp(i)=inner(roots[i],roots[i]);
				eovlp(i)=(inner(world,roots[i].afunction,roots[i].afunction)).trace(noct(nemo->aeps))
						+(inner(world,roots[i].bfunction,roots[i].bfunction)).trace(noct(nemo->beps));
			}

			if (cis_param.printlevel()>2) {
				print("fock_pt");
				print(fock_pt);

				print("ovlp (ab), eovlp (ab)");
				print(ovlp);
				print(eovlp);
			}

			omega=Tensor<double>(roots.size());
			for (std::size_t i=0; i<roots.size(); ++i) {
				omega(i)=real((fock_pt(i,i)-eovlp(i,i))/(ovlp(i,i)));
			}
			print("omega from the perturbed fock matrix approach");
			print(omega);

			orthonormalize(roots, fock_pt);
//			orthonormalize(roots);
		} else {		// do_kain

			// update the residuals
			std::vector<complex_function_3d> allres, oldx;
			for (int iroot=0; iroot<cis_param.guess_excitations(); ++iroot) {
				root& thisroot=roots[iroot];

				std::vector<complex_function_3d> residuals=compute_residuals(thisroot);
				allres=append(allres,residuals);
				oldx=append(oldx,thisroot.afunction);
				oldx=append(oldx,thisroot.bfunction);

				thisroot.delta=norm2(world,residuals);

				if (omega.size()>0) thisroot.omega=omega(iroot);
				if (cis_param.omega()!=0) {
					thisroot.omega=cis_param.omega();
					printf("\n\nset excitation energy manually to %8.4f\n\n",thisroot.omega);
				}

				print_size(world,thisroot.afunction,"afunction1");
				print_size(world,thisroot.bfunction,"bfunction1");
				int i=0;
				for (auto& f : thisroot.afunction) save(abs_square(f),"afunction_root"+stringify(iroot)+"mo"+stringify(i++));

			}

			truncate(world,allres);
			truncate(world,oldx);
			std::vector<complex_function_3d> newx=allsolver.update(oldx,allres,0.01,3);
			truncate(world,newx);

			// distribute new x functions
			for (int iroot=0; iroot<cis_param.guess_excitations(); ++iroot) {
				root& thisroot=roots[iroot];
				auto [atmp, rest] = split(newx,thisroot.afunction.size());
				auto [btmp, rest2] = split(rest,thisroot.bfunction.size());
				newx=rest2;
				thisroot.afunction=atmp;
				thisroot.bfunction=btmp;

			}
			MADNESS_ASSERT(newx.size()==0);


			orthonormalize(roots);
			normalize(roots);

			std::vector<double> rnorm=norm2s(world,allres);
			double rmax=*(std::max_element(rnorm.begin(),rnorm.end()));
			if (rmax<cis_param.dconv()) break;
		}	// if do_kain
	}	// iter
}

void Zcis::compute_potentials(std::vector<root>& roots, const real_function_3d& totdens) const {
	// some numbers
	const int anoct=noct(nemo->aeps).size();
	const int bnoct=noct(nemo->beps).size();
	//const int noct=anoct+bnoct;

	for (int iroot=0; iroot<cis_param.guess_excitations(); ++iroot) {
		root& thisroot=roots[iroot];
		if (thisroot.delta<cis_param.dconv()) continue;

		std::vector<complex_function_3d> td_a=mul(world,thisroot.afunction,active_mo(nemo->amo));
		std::vector<complex_function_3d> td_b=mul(world,thisroot.bfunction,active_mo(nemo->bmo));
		int i=cis_param.freeze();
		for (auto& a : td_a) save(abs_square(a),"transition_density_a"+stringify(i++)+"root"+stringify(iroot));

		std::vector<complex_function_3d> apot=zero_functions_compressed<double_complex,3>(world,anoct);
		std::vector<complex_function_3d> bpot=zero_functions_compressed<double_complex,3>(world,bnoct);

		// compute the perturbed density
		complex_function_3d denspt=((dot(world,conj(world,thisroot.afunction),active_mo(nemo->amo)) +
				dot(world,conj(world,thisroot.bfunction),active_mo(nemo->bmo))));
		if (nemo->get_calc_param().spin_restricted()) denspt.scale(2.0);

		denspt=conj(denspt);
		denspt.truncate();
		denspt.print_size("denspt");
		save(abs_square(denspt),"denspt_root"+stringify(iroot));

		for (std::string spin : {"alpha","beta"}) {
			if (nemo->get_calc_param().spin_restricted() and (spin=="beta")) continue;

			const std::vector<complex_function_3d>& mo=(spin=="alpha") ? nemo->amo : nemo->bmo;
			const std::vector<complex_function_3d>& act_mo=(spin=="alpha") ? active_mo(nemo->amo) : active_mo(nemo->bmo);
			if (act_mo.size()==0) continue;

			std::vector<complex_function_3d>& pot=(spin=="alpha") ? apot : bpot;
			std::vector<complex_function_3d> x= (spin=="alpha") ? thisroot.afunction : thisroot.bfunction;
			const QProjector<double_complex,3>& Q=(spin=="alpha") ? Qa : Qb;
			Tensor<double> occ = (spin=="alpha") ? Tensor<double>(nemo->amo.size()) : Tensor<double>(nemo->bmo.size());
			occ=1.0;

			/// zeroth order Fock operator acting on the x functions
			Znemo::potentials pot_mo=nemo->compute_potentials(mo,totdens,x);
			if (spin=="beta") scale(world,pot_mo.spin_zeeman_mo,-1.0);

//			nemo->compute_vmat(mo,vnemo,vlznemo,dianemo,spin_zeeman_nemo,knemo,jnemo);


			pot+=(pot_mo.vnuc_mo+pot_mo.lz_mo+pot_mo.diamagnetic_mo+pot_mo.spin_zeeman_mo
					-pot_mo.K_mo+pot_mo.J_mo);	// need += b/c of reference
			truncate(world,pot);

			// perturbed Fock operator acting on the reference orbitals
			Coulomb<double_complex,3> Jp(world);
			complex_function_3d Jp_pot = Jp.compute_potential(denspt);

			Exchange<double_complex,3> Kp(world,nemo->get_calc_param().lo());
            Kp.set_bra_and_ket(conj(world, act_mo), x);
			pot+=Q(Jp_pot*act_mo - Kp(act_mo));
			truncate(world,pot);
		}
		thisroot.apot=copy(world,apot);
		thisroot.bpot=copy(world,bpot);
	}

}

Tensor<double_complex> Zcis::compute_fock_pt(const std::vector<root>& roots) const {

	std::size_t dim=roots.size();
	Kinetic<double_complex,3> T(world);
	Tensor<double_complex> fock_pt_a(dim,dim), fock_pt_b(dim,dim);
	Tensor<double_complex> Tmat;
	for (std::size_t r=0; r<dim; ++r) {

		for (std::size_t p=0; p<dim; ++p) {
			fock_pt_a(r,p) = inner(roots[r].afunction,roots[p].apot);
			Tmat=T(roots[r].afunction,roots[p].afunction);
			for (int i=0; i<Tmat.dim(0); ++i) fock_pt_a(r,p)+=Tmat(i,i);

			if (roots[p].bpot.size()>0) {
				fock_pt_b(r,p) = inner(roots[r].bfunction,roots[p].bpot);
				Tmat=T(roots[r].bfunction,roots[p].bfunction);
				for (int i=0; i<Tmat.dim(0); ++i) fock_pt_b(r,p)+=Tmat(i,i);
			}
		}
	}
	return fock_pt_a+fock_pt_b;

}


std::vector<complex_function_3d> Zcis::compute_residuals(root& root) const {


	std::vector<complex_function_3d> x=append(root.afunction,root.bfunction);
	std::vector<complex_function_3d> pot = append(root.apot,root.bpot);
	Tensor<double> eps=copy(concatenate(noct(nemo->aeps),noct(nemo->beps))) + root.omega;

	double global_shift=nemo->param.shift();
	eps+=global_shift;
//	print("compute_residuals, eps",eps);

	// make the BSH operator
	std::vector<std::shared_ptr<SeparatedConvolution<double,3> > > bsh(pot.size());

	for (std::size_t p = 0; p < bsh.size(); p++) {

		// if eps is above zero we have an unbound state (or are early in the iteration)
		// however this needs a shift of the potential
		// we shift to -0.05 (same as moldft does, no particular reason)
		if(eps[p]>0.0){
			print("potential shift needed for V" + std::to_string(p+cis_param.freeze()));
			double shift = eps[p]+0.05;
			eps[p] -= shift;
			pot[p] -= (-2.0*shift*x[p]);
		}
		MADNESS_ASSERT(not(eps[p]>0.0));
//		print("assigning eps to bsh ",p,eps[p]);
		bsh[p] = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world,
				sqrt(-2.0 * (eps[p])), 1.e-4, cis_param.thresh()));
	}

	// apply the BSH operator
	std::vector<complex_function_3d> tmp0=(apply(world,bsh,-2.0*pot-2.0*global_shift*x));

	auto [atmp, btmp] = split(tmp0,root.afunction.size());
	std::vector<complex_function_3d> tmp=append(Qa(atmp),Qb(btmp));

	std::vector<complex_function_3d> residual=x-tmp;
	truncate(world,residual);

	// Calculate Second Order Energy Update	(complex conjugate implicit in inner)
	double_complex tmp1 = inner(world,residual,pot+global_shift*x).sum();
//	double_complex tmp11= inner(world,pot,x).sum();
	// squared norm of GVpsi (Psi_tilde)
	double_complex tmp2 = inner(world,tmp,tmp).sum();

//	print(" sou = Vp * res / Vp * tmp");
//	double_complex tmp1=inner(pot,residual);
//	double_complex tmp2=inner(x,tmp);

	print("tmp1,tmp2,tmp11",tmp1,tmp2);
	const double sou=-real(tmp1)/real(tmp2);
	root.omega+=sou;
	root.energy_change=sou;

	return residual*nemo->sbox;

}

std::vector<Zcis::root> Zcis::read_guess() const {
	std::vector<root> guess(cis_param.guess_excitations());
	std::string name="cis_guess";
	print("reading cis guess from file",name);

	archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str(), 1);
	std::size_t size1=0, size2=0, size3=0, size4=0; // zeroed to silence clang++ uninit warnings

	for (root& g : guess) {
		ar & g.omega & g.delta & g.energy_change ;
		ar & size1 & size2 & size3 & size4;
        g.afunction.resize(size1);
        g.bfunction.resize(size2);
        g.apot.resize(size3);
        g.bpot.resize(size4);
        for (auto& a: g.afunction) ar & a;
        for (auto& a: g.bfunction) ar & a;
        for (auto& a: g.apot) ar & a;
        for (auto& a: g.bpot) ar & a;
	}
	return guess;
}

void Zcis::save_guess(const std::vector<root>& roots) const {
	std::string name="cis_guess";
	print("saving cis guess to file",name);
	archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(), 1);
	for (const root& g : roots) {
		ar & g.omega & g.delta & g.energy_change ;
		std::size_t size1=g.afunction.size();
		std::size_t size2=g.bfunction.size();
		std::size_t size3=g.apot.size();
		std::size_t size4=g.bpot.size();
		ar & size1 & size2 & size3 & size4;
        for (const auto& a: g.afunction) ar & a;
        for (const auto& a: g.bfunction) ar & a;
        for (const auto& a: g.apot) ar & a;
        for (const auto& a: g.bpot) ar & a;
	}
}


std::vector<Zcis::root> Zcis::make_guess() const {

	std::vector<complex_function_3d> avirtuals, bvirtuals;
	Tensor<double> aveps,bveps;		// virtual orbital energies for alpha and beta
	real_function_3d adens=real(dot(world,conj(world,nemo->amo),nemo->amo)).truncate();
	real_function_3d bdens=real(dot(world,conj(world,nemo->bmo),nemo->bmo)).truncate();
	real_function_3d density=adens + bdens;
	if (nemo->get_calc_param().spin_restricted()) density.scale(2.0);

	density.print_size("density in make_guess");

	std::vector<root> guess(cis_param.guess_excitations());

	for (std::string spin : {"alpha","beta"}) {
		print("spin " ,spin);
		std::vector<complex_function_3d>& virtuals =(spin=="alpha") ? avirtuals : bvirtuals;
		std::vector<complex_function_3d> mo= (spin=="alpha") ? active_mo(nemo->amo) : active_mo(nemo->bmo);
		if (mo.size()==0) continue;

		// prepare the list of excitation operators and copied seeds
		std::vector<coord_3d> centers = guessfactory::compute_centroids(mo);
		std::vector<std::pair<std::vector<complex_function_3d>, std::string> > exlist;
		{
	//		std::vector<std::string> exop_strings=cis_param.exops();
			std::vector<std::string> exop_strings=(guessfactory::make_predefined_exop_strings(cis_param.guess_excitation_operators()));
			for(const auto& ex: exop_strings){
				std::vector<complex_function_3d> cseed=copy(world,mo,false);
				exlist.push_back(std::make_pair(cseed,ex));
			}
			print(exop_strings);
		}
		world.gop.fence();
		std::cout << "will create " << exlist.size()*centers.size() << " virtuals, from " << centers.size()
				<< " seeds and " << exlist.size() << " excitation operators"   << std::endl;

		// create the virtuals by unary operations: multiply excitation operators with seeds
		for(auto it : exlist){
			virtuals=append(virtuals,guessfactory::apply_trigonometric_exop(it.first,it.second,centers,false));
		}
		world.gop.fence();


		if (spin=="alpha") virtuals=Qa(virtuals);
		else if (spin=="beta") virtuals=Qb(virtuals);

		// remove linear dependencies
		const size_t spre=virtuals.size();
		virtuals=orthonormalize_canonical(virtuals,1.e-6);
		if(virtuals.size()!=spre) std::cout << "removed " << spre-virtuals.size() << " virtuals due to linear dependencies" << std::endl;
	}

	// canonicalize virtuals and set up the CIS matrix
	if (avirtuals.size()>0) canonicalize(nemo->amo,density,avirtuals,aveps);
	if (bvirtuals.size()>0) canonicalize(nemo->bmo,density,bvirtuals,bveps);

	// active orbital energies only
	int anoct=noct(nemo->aeps).size();
	int bnoct=noct(nemo->beps).size();
	Tensor<double> oeps=concatenate(noct(nemo->aeps),noct(nemo->beps));
	Tensor<double> veps=concatenate(aveps,bveps);
//	auto [ae,be] = split(veps,nemo->aeps.size());

	Tensor<double_complex> MCISa=make_CIS_matrix(aveps,noct(nemo->aeps));
	Tensor<double_complex> MCISb=make_CIS_matrix(bveps,noct(nemo->beps));
	Tensor<double_complex> MCIS(MCISa.dim(0)+MCISb.dim(0),MCISa.dim(1)+MCISb.dim(1));
	if (MCISa.size()>0) MCIS(Slice(0,MCISa.dim(0)-1,1),Slice(0,MCISa.dim(0)-1,1))=MCISa;
	if (MCISb.size()>0) MCIS(Slice(MCISa.dim(0),-1,1),Slice(MCISa.dim(0),-1,1))=MCISb;

	// compute matrix CIS excitation energies
	Tensor<double_complex> U,UU;
	Tensor<double> evals;
	syev(MCIS, U, evals);
	print("MCIS excitation energies ",evals);
	if (evals.size()<cis_param.guess_excitations()) {
		print("number of requested excitations larger than the number of MCIS energies");
		MADNESS_EXCEPTION("increase the guess",1);
	}
//	print(real(U));
//	print("MCIS excitation vectors");
//	for (int i=0; i<cis_param.guess_excitations(); ++i) print(real(U(_,i)));


	// assemble the guess vectors

	std::vector<complex_function_3d> abvirtuals=append(avirtuals,bvirtuals);
	//int nvirt = abvirtuals.size();
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


	for (size_t I = 0, iexcitation=0; I < size_t(MCIS.dim(1)); ++I, iexcitation++) {
	        if (iexcitation >= size_t(cis_param.guess_excitations())) break;

		guess[iexcitation].afunction=zero_functions_compressed<double_complex,3>(world,anoct);
		guess[iexcitation].bfunction=zero_functions_compressed<double_complex,3>(world,bnoct);

		if (evals(I) < 0.0 && world.rank() == 0)
			MADNESS_EXCEPTION("NEGATIVE EIGENVALUE IN INITIAL DIAGONALIZATION: CHECK YOUR REFERENCE!\n",0);

		guess[iexcitation].omega = evals(I);

		// alpha part
		if (Ua.size()>0) {
		        for (size_t J = 0; J < size_t(Ua.dim(0)); ++J) {
				const int b = get_vir_idx_a(J);
				const int j = get_occ_idx_a(J);
				const double_complex xjb = Ua(J, I);
				guess[iexcitation].afunction[j] += xjb * avirtuals[b];
			}
		}

		// beta part
		if (Ub.size()>0) {
		        for (size_t J = 0; J < size_t(Ub.dim(0)); ++J) {
				const int b = get_vir_idx_b(J);
				const int j = get_occ_idx_b(J);
				const double_complex xjb = Ub(J, I);
				guess[iexcitation].bfunction[j] += xjb * bvirtuals[b];
			}
		}
		print_size(world,guess[iexcitation].afunction,"afunction");
		print_size(world,guess[iexcitation].bfunction,"bfunction");
	}
	for (auto& x : guess) truncate(world, x.afunction, cis_param.thresh());
	for (auto& x : guess) truncate(world, x.bfunction, cis_param.thresh());
	orthonormalize(guess);

	return guess;

}

void Zcis::canonicalize(const std::vector<complex_function_3d>& mo, const real_function_3d& density,
		std::vector<complex_function_3d>& virtuals, Tensor<double>& veps) const {

	Znemo::potentials pot=nemo->compute_potentials(mo,density,virtuals);
	Tensor<double_complex> fock=nemo->compute_vmat(virtuals,pot);
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

Tensor<double_complex> Zcis::make_CIS_matrix(const Tensor<double>& veps, const Tensor<double>& act) const {

	// the cis matrix is indexed by ij and ab
	// we will use the combined indixes from ia and jb named I and J
	// in order to not be confused we use the following helper functions
	const int nact = act.size();
	// determines for which orbitals (couting from the HOMO downwards) the off-diagonal elements will be computed
	// this simplifies the guess
	const int nvirt = veps.size();

//	print("eps(occ), eps(vir)");
//	print(act);
//	print(veps);
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

void Zcis::orthonormalize(std::vector<root>& roots, const Tensor<double_complex>& fock_pt) const {
	normalize(roots);
	Tensor<double_complex> ovlp(roots.size(),roots.size());
	for (size_t i=0; i<roots.size(); ++i) {
		for (size_t j=0; j<roots.size(); ++j) {
			ovlp(i,j)=inner(roots[i],roots[j]);
		}
	}

	// orthonormalize alpha and beta part
	Tensor<double_complex> U;
	Tensor<double> e;
	sygv(fock_pt, ovlp, 1, U, e);
	std::vector<root> tmp = transform(world,roots,U);
	for (size_t i=0; i<tmp.size(); ++i) {
		roots[i].afunction=tmp[i].afunction;
		roots[i].bfunction=tmp[i].bfunction;
		roots[i].apot=tmp[i].apot;
		roots[i].bpot=tmp[i].bpot;
	}
	normalize(roots);
}


void Zcis::orthonormalize(std::vector<root>& roots) const {
	normalize(roots);

	double maxq;
	do {

		Tensor<double_complex> ovlp(roots.size(),roots.size());
		for (size_t i=0; i<roots.size(); ++i) {
			for (size_t j=0; j<roots.size(); ++j) {
				ovlp(i,j)=inner(roots[i].afunction,roots[j].afunction)+inner(roots[i].bfunction,roots[j].bfunction);
			}
		}

		Tensor<double_complex> Q = NemoBase::Q2(ovlp);
		maxq = 0.0;
		for (int j=1; j<Q.dim(0); j++)
			for (int i=0; i<j; i++) {
				maxq = std::max(std::abs(Q(j,i)),maxq);
			}
//		amo = transform(world, amo, Q, 0.0, true);
		std::vector<root> tmp = transform(world,roots,Q);
		for (size_t i=0; i<tmp.size(); ++i) {
			roots[i].afunction=tmp[i].afunction;
			roots[i].bfunction=tmp[i].bfunction;
			roots[i].apot=tmp[i].apot;
			roots[i].bpot=tmp[i].bpot;
		}
//		print("maxq",maxq);
	} while (maxq>0.01);

	normalize(roots);
}


void Zcis::normalize(std::vector<root>& roots) const {

	for (auto root : roots) {
		double na=0.0, nb=0.0;
        if (root.afunction.size()>0) na = norm2(world, root.afunction);
        if (root.bfunction.size()>0) nb = norm2(world, root.bfunction);
        double n=sqrt(na*na + nb*nb);
        scale(world,root.afunction,1.0/n);
        scale(world,root.bfunction,1.0/n);
	}
}



} /* namespace madness */
