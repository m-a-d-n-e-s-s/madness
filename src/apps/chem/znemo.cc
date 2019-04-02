/*
 * Nemocomplex.cc
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#include <madness/mra/mra.h>
#include "znemo.h"


namespace madness {

/// compute the molecular energy
double Znemo::value() {

	// compute the molecular potential
	const Molecule& m=molecule;
	auto molecular_potential = [m] (const coord_3d& r) {
		return m.nuclear_attraction_potential(r[0],r[1],r[2]);};
	vnuclear=real_factory_3d(world).functor(molecular_potential).thresh(FunctionDefaults<3>::get_thresh()*0.1);
	vnuclear.set_thresh(FunctionDefaults<3>::get_thresh());


	// read the guess orbitals
	try {
		read_orbitals();

	} catch(...) {
		amo=read_guess("alpha");
		if (have_beta()) bmo=read_guess("beta");
		aeps=Tensor<double>(amo.size());
		beps=Tensor<double>(bmo.size());
	}

	double thresh=FunctionDefaults<3>::get_thresh();
	double energy=1.e10;
	double oldenergy=0.0;

	// the diamagnetic box

	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solvera(allocator(world,amo.size()));
	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solverb(allocator(world,bmo.size()));
	solvera.set_maxsub(cparam.maxsub);
	solvera.do_print=(param.printlevel()>2);
	solverb.set_maxsub(cparam.maxsub);
	solverb.do_print=(param.printlevel()>2);

	// increase the magnetic field
	for (int i=0; i<param.B().size(); ++i) {
		B=param.B()[i];
		print("solving for magnetic field B=",B);

		diamagnetic_boxed=make_diamagnetic_boxed();

		// set end of iteration cycles for intermediate calculations
//		int maxiter = (i==param.B().size()-1) ? cparam.maxiter : 20;
		int maxiter = cparam.maxiter;
		double dconv = (i==param.B().size()-1) ? cparam.dconv : 1.e-2;
		double na=1.0,nb=1.0;	// residual norms

		solvera.clear_subspace();
		solverb.clear_subspace();
		bool converged=false;

		// iterate the SCF cycles
		for (int iter=0; iter<maxiter; ++iter) {

			// compute the density
			real_function_3d density=sum(world,abssq(world,amo));
			if (have_beta()) density+=sum(world,abssq(world,bmo));
			if (cparam.spin_restricted) density*=2.0;
			density.truncate();

			// compute the fock matrix
			std::vector<complex_function_3d> Vnemoa, Vnemob;
			Tensor<double_complex> focka, fockb(0l,0l);
			potentials apot(world,amo.size()), bpot(world,bmo.size());

			apot=compute_potentials(amo, density, amo);
			Vnemoa=apot.vnuc_mo+apot.lz_mo+apot.diamagnetic_mo+apot.spin_zeeman_mo-apot.K_mo+apot.J_mo;
			truncate(world,Vnemoa,thresh*0.1);
			Kinetic<double_complex,3> T(world);
			focka=T(amo,amo) + compute_vmat(amo,apot);

			if (have_beta()) {
				bpot=compute_potentials(bmo, density, bmo);
				scale(world,bpot.spin_zeeman_mo,-1.0);
				Vnemob=bpot.vnuc_mo+bpot.lz_mo+bpot.diamagnetic_mo+bpot.spin_zeeman_mo-bpot.K_mo+bpot.J_mo;
				truncate(world,Vnemob,thresh*0.1);
				fockb=T(bmo,bmo) + compute_vmat(bmo,bpot);
			}

			if (world.rank()==0 and (param.printlevel()>1)) {
				print("Fock matrix");
				print(focka);
				print(fockb);
			}

			oldenergy=energy;
			energy=compute_energy(amo,apot,bmo,bpot,param.printlevel()>1);


			Tensor<double_complex> ovlp=matrix_inner(world,amo,amo);

			canonicalize(amo,Vnemoa,solvera,focka,ovlp);

			if (have_beta()) {
				Tensor<double_complex> ovlp=matrix_inner(world,bmo,bmo);
				canonicalize(bmo,Vnemob,solverb,fockb,ovlp);
			}

			if (param.printlevel()>2) print("using fock matrix for the orbital energies");
			for (int i=0; i<focka.dim(0); ++i) aeps(i)=real(focka(i,i));
			for (int i=0; i<fockb.dim(0); ++i) beps(i)=real(fockb(i,i));


			if (world.rank()==0 and (param.printlevel()>1)) {
				print("orbital energies alpha",aeps);
				print("orbital energies beta ",beps);
			}
			if (world.rank() == 0) {
				printf("finished iteration %2d at time %8.1fs with energy, norms %12.8f %12.8f %12.8f\n",
						iter, wall_time(), energy, na, nb);
			}

			if (std::abs(oldenergy-energy)<cparam.econv and (sqrt(na*na+nb*nb)<dconv)) {
				print("energy converged");
				save_orbitals("converged");
				converged=true;
			}
			if (converged) break;

			// compute the residual of the Greens' function
			std::vector<complex_function_3d> resa=compute_residuals(Vnemoa,amo,aeps);
			truncate(world,resa,thresh*0.1);
			std::vector<double> normsa=norm2s(world,resa);

			na=0.0; nb=0.0;
			for (auto nn : normsa) {na+=nn*nn;}
			na=sqrt(na);
			std::vector<complex_function_3d> amo_new=solvera.update(amo,resa,0.01,3);
			amo_new=sbox*amo_new;
			do_step_restriction(amo,amo_new);
			amo=amo_new;
//			amo=orthonormalize_symmetric(amo);
			orthonormalize(amo);
			truncate(world,amo);
			orthonormalize(amo);


			if (have_beta()) {
				std::vector<complex_function_3d> resb=compute_residuals(Vnemob,bmo,beps);
				truncate(world,resb);
				std::vector<double> normsb=norm2s(world,resb);
				for (auto nn : normsb) {nb+=nn*nn;}
				nb=sqrt(nb);
				std::vector<complex_function_3d> bmo_new=solverb.update(bmo,resb,0.01,3);
				bmo_new=sbox*bmo_new;
				do_step_restriction(bmo,bmo_new);
				bmo=bmo_new;
				orthonormalize(bmo);
//				bmo=orthonormalize_symmetric(bmo);
				truncate(world,bmo);
			}
		}

		// final energy computation
		real_function_3d density=sum(world,abssq(world,amo));
		if (have_beta()) density+=sum(world,abssq(world,bmo));
		if (cparam.spin_restricted) density*=2.0;
		density.truncate();

		potentials apot=compute_potentials(amo,density,amo);
		potentials bpot=compute_potentials(bmo,density,bmo);
		scale(world,bpot.spin_zeeman_mo,-1.0);
		double energy=compute_energy(amo,apot,bmo,bpot,true);


		if (world.rank()==0) {
			print("orbital energies alpha",aeps);
			print("orbital energies beta ",beps);

		}

	}
	save_orbitals("final");
	save_orbitals();

	analyze();

	return energy;
}

void Znemo::analyze() const {

	// compute the current density
	std::vector<real_function_3d> j=compute_current_density(amo,bmo);
	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

	// compute the expectation values of the Lz operator
	std::vector<complex_function_3d> lzamo=Lz(amo);
	std::vector<complex_function_3d> lzbmo=Lz(bmo);
	Tensor<double_complex> lza_exp=inner(world,amo,lzamo);
	Tensor<double_complex> lzb_exp=inner(world,bmo,lzbmo);
	print("< amo | lz | amo >",lza_exp);
	print("< bmo | lz | bmo >",lzb_exp);

	// compute magnetic vector potential
	Tensor<double> Bvec(3);
    Bvec(2)=B;
    std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,Bvec);
    print("B",Bvec);
    std::vector<real_function_3d> Btest=rot(A);
    Btest[2]=Btest[2]-B;
//    Btest[2]-=bla;
    double n=norm2(world,Btest);
    print("n(Btest-B)",n);
    Tensor<double> p_exp=compute_kinetic_momentum();
    Tensor<double> A_exp=compute_magnetic_potential_expectation(A);

    print("<p>       ",p_exp);
    print("<A>       ",A_exp);
    print("(p-eA)    ",p_exp+A_exp);

    // compute the standard kinetic gauge origin, defined as the gauge origin, where the
    // expectation value of the kinetic momentum p vanishes
    const long nmo=amo.size()+bmo.size();
    Tensor<double> v=-1.0/nmo*p_exp;
    Tensor<double> S=compute_standard_gauge_shift(p_exp);
    print("standard gauge shift S",S);

    // compute the standardized components of the canonical momentum square
    const double v2=v.trace(v);	// term B^2(Sx^2+Sy^2)
    const double vp=v.trace(p_exp);
    const double vA=v.trace(A_exp);
    print("v2, vp, vA", v2, vp, vA);

    print("expectation values in standard gauge");
    print("Delta 1/2 <p^2>  ",0.5*(nmo*v2 + 2.0*vp));


}


double Znemo::compute_energy(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print) const {

    double fac= cparam.spin_restricted ? 2.0 : 1.0;

    double_complex kinetic=0.0;
    for (int axis = 0; axis < 3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex, 3>(world, axis);
        const std::vector<complex_function_3d> damo = apply(world, D, amo);
        const std::vector<complex_function_3d> dbmo = apply(world, D, bmo);
        kinetic += fac* 0.5 * (inner(world, damo, damo).sum() + inner(world, dbmo, dbmo).sum());
    }

    double_complex nuclear_potential=fac*(inner(world,amo,apot.vnuc_mo).sum()+inner(world,bmo,bpot.vnuc_mo).sum());
    double_complex diamagnetic=fac*(inner(world,amo,apot.diamagnetic_mo).sum()+inner(world,bmo,bpot.diamagnetic_mo).sum());
    double_complex lz=fac*(inner(world,amo,apot.lz_mo).sum()+inner(world,bmo,bpot.lz_mo).sum());
    double_complex spin_zeeman=fac*(inner(world,amo,apot.spin_zeeman_mo).sum()+inner(world,bmo,bpot.spin_zeeman_mo).sum());
    double_complex coulomb=fac*0.5*(inner(world,amo,apot.J_mo).sum()+inner(world,bmo,bpot.J_mo).sum());
    double_complex exchange=fac*0.5*(inner(world,amo,apot.K_mo).sum()+inner(world,bmo,bpot.K_mo).sum());

    double_complex energy=kinetic + nuclear_potential + molecule.nuclear_repulsion_energy() +
    		diamagnetic + lz + spin_zeeman + coulomb - exchange;

	if (world.rank()==0 and do_print) {
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", molecule.nuclear_repulsion_energy());
		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", real(lz));
		printf("  spin zeeman term    %12.8f \n", real(spin_zeeman));
		printf("  Coulomb             %12.8f \n", real(coulomb));
		printf("  exchange            %12.8f \n", real(exchange));
		printf("  total               %12.8f \n", real(energy));
	}

	if(fabs(imag(energy))>1.e-8) {


		print("real part");
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", molecule.nuclear_repulsion_energy());
		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", real(lz));
		printf("  spin zeeman term    %12.8f \n", real(spin_zeeman));
		printf("  Coulomb             %12.8f \n", real(coulomb));
		printf("  exchange            %12.8f \n", real(exchange));
		printf("  total               %12.8f \n", real(energy));

		print("imaginary part");
		printf("  kinetic energy      %12.8f \n", imag(kinetic));
		printf("  nuclear potential   %12.8f \n", imag(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", 0.0);
		printf("  diamagnetic term    %12.8f \n", imag(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", imag(lz));
		printf("  spin zeeman term    %12.8f \n", imag(spin_zeeman));
		printf("  Coulomb             %12.8f \n", imag(coulomb));
		printf("  exchange            %12.8f \n", imag(exchange));
		printf("  total               %12.8f \n", imag(energy));

		print("imaginary part of the energy",energy.imag());

		MADNESS_EXCEPTION("complex energy computation.. ",1);
	}

	return real(energy);
}

/// following Lazeretti, J. Mol. Struct, 313 (1994)
std::vector<real_function_3d> Znemo::compute_current_density(
		const std::vector<complex_function_3d>& alpha_mo,
		const std::vector<complex_function_3d>& beta_mo) const {

	// compute vec r
	std::vector<real_function_3d> r(3);
    r[0]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
    r[1]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});
    r[2]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[2];});
    real_function_3d one=real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});

	// the vector potential A=1/2 B x r
	std::vector<real_function_3d> Bvec=zero_functions_compressed<double,3>(world,3);
	Bvec[2]=B*one;
	reconstruct(world,Bvec);
	reconstruct(world,r);
	std::vector<real_function_3d> A=0.5*cross(Bvec,r);

	// test consistency
	Bvec=rot(A);
	double bnorm2=(Bvec[2]-B*one).norm2();
	MADNESS_ASSERT(bnorm2<1.e-8);

	// compute density and spin density
	real_function_3d adens=real_factory_3d(world);
	real_function_3d bdens=real_factory_3d(world);
	for (auto& mo : alpha_mo) adens+=abs_square(mo);
	for (auto& mo : beta_mo) bdens+=abs_square(mo);
	real_function_3d density=adens+bdens;
	real_function_3d spin_density=adens-bdens;

//	SeparatedConvolution<double,3> smooth=SmoothingOperator3D(world,1.e-3);
//	spin_density=smooth(spin_density);

	std::vector<real_function_3d> vspin_density=zero_functions_compressed<double,3>(world,3);;
	vspin_density[2]=spin_density;


	// compute first contribution to current density from orbitals:
	// psi^* p psi = i psi^* del psi
	std::vector<complex_function_3d> j=zero_functions_compressed<double_complex,3>(world,3);
	for (auto& mo : alpha_mo) j+=double_complex(0.0,1.0)*(conj(mo)*grad(mo));
	for (auto& mo : beta_mo) j+=double_complex(0.0,1.0)*(conj(mo)*grad(mo));

	// compute density contribution and spin contribution
	j-=convert<double,double_complex,3>(world,A*density);
	j+=convert<double,double_complex,3>(world,0.5*rot(vspin_density));

	std::vector<real_function_3d>  realj=real(j);

//	std::vector<double> n1=norm2s(world,real(j));
//	std::vector<double> n2=norm2s(world,imag(j));
//	print("norm(re(j))",n1);
//	print("norm(im(j))",n2);

	// sanity check

	real_function_3d null=div(realj);
	double n3=null.norm2();
	print("div(j)",n3);
//	MADNESS_ASSERT(n3<FunctionDefaults<3>::get_thresh());


	return realj;
}


void Znemo::test_compute_current_density() const {

	complex_function_3d pp=complex_factory_3d(world).f(p_plus);
	double norm=pp.norm2();
	pp.scale(1/norm);
	double_complex l=inner(amo[0],Lz(amo[0]));
	print("<p+|Lz|p+>",l);

//	std::vector<complex_function_3d> vpp(1,amo);
	std::vector<real_function_3d> j=compute_current_density(amo,std::vector<complex_function_3d>());

	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

}



void Znemo::do_step_restriction(const std::vector<complex_function_3d>& mo,
		std::vector<complex_function_3d>& mo_new) const {
    PROFILE_MEMBER_FUNC(SCF);
    std::vector<double> anorm = norm2s(world, sub(world, mo, mo_new));
    int nres = 0;
    for (unsigned int i = 0; i < mo.size(); ++i) {
        if (anorm[i] > cparam.maxrotn) {
            double s = cparam.maxrotn / anorm[i];
            ++nres;
            if (world.rank() == 0) {
                if (nres == 1)
                    printf("  restricting step ");
                printf(" %d", i);
            }
            mo_new[i].gaxpy(s, mo[i], 1.0 - s, false);
        }
    }
    if (nres > 0 && world.rank() == 0)
        printf("\n");

    world.gop.fence();
}


/// compute the action of the Lz =i r x del operator on rhs
std::vector<complex_function_3d> Znemo::Lz(const std::vector<complex_function_3d>& rhs) const {
	// the operator in cartesian components as
	// L_z =  - i (x del_y - y del_x)

	if (rhs.size()==0) return std::vector<complex_function_3d>(0);
	auto monomial_x = [] (const coord_3d& r) {return r[0];};
	auto monomial_y = [] (const coord_3d& r) {return r[1];};

	World& world=rhs.front().world();
    complex_derivative_3d Dx = free_space_derivative<double_complex,3>(world, 0);
    complex_derivative_3d Dy = free_space_derivative<double_complex,3>(world, 1);
    real_function_3d x=real_factory_3d(world).functor(monomial_x);
    real_function_3d y=real_factory_3d(world).functor(monomial_y);

    std::vector<complex_function_3d> delx=apply(world,Dx,rhs);
    std::vector<complex_function_3d> dely=apply(world,Dy,rhs);

    std::vector<complex_function_3d> result1=x*dely - y*delx;
    std::vector<complex_function_3d> result=double_complex(0.0,-1.0)*result1;
	return result;

}

/// read the guess orbitals from a previous nemo or moldft calculation
std::vector<complex_function_3d> Znemo::read_guess(const std::string& spin) const {

	int nmo= (spin=="alpha") ? cparam.nalpha : cparam.nbeta;
	std::vector<real_function_3d> real_mo=zero_functions<double,3>(world,nmo);

	// load the converged orbitals
    for (std::size_t imo = 0; imo < nmo; ++imo) {
    	print("loading mos ",spin,imo);
    	load(real_mo[imo], "nemo_"+spin + stringify(imo));
    }
    return convert<double,double_complex,3>(world,real_mo);
}

/// compute the potential operators applied on the orbitals
Znemo::potentials Znemo::compute_potentials(const std::vector<complex_function_3d>& mo,
		const real_function_3d& density,
		std::vector<complex_function_3d>& rhs) const {

	// prepare exchange operator
	Exchange<double_complex,3> K=Exchange<double_complex,3>(world);
	Tensor<double> occ(mo.size());
	occ=1.0;
	K.set_parameters(conj(world,mo),mo,occ,cparam.lo,cparam.econv);

	potentials pot(world,rhs.size());

	pot.J_mo=(*coulop)(density)*rhs;
	pot.K_mo=K(rhs);
	pot.vnuc_mo=vnuclear*rhs;
	pot.lz_mo=0.5*B*Lz(rhs);
	pot.diamagnetic_mo=0.125*B*B*diamagnetic(rhs);
	pot.spin_zeeman_mo=B*0.5*rhs;

	truncate(world,pot.J_mo);
	truncate(world,pot.K_mo);
	truncate(world,pot.vnuc_mo);
	truncate(world,pot.lz_mo);
	truncate(world,pot.diamagnetic_mo);
	truncate(world,pot.spin_zeeman_mo);
	return pot;
};


Tensor<double_complex> Znemo::compute_vmat(const std::vector<complex_function_3d>& mo,
		const potentials& pot) const {

	Tensor<double_complex> Vnucmat=matrix_inner(world,mo,pot.vnuc_mo);
	Tensor<double_complex> lzmat=matrix_inner(world,mo,pot.lz_mo);
	Tensor<double_complex> diamat=matrix_inner(world,mo,pot.diamagnetic_mo);
	Tensor<double_complex> spin_zeeman_mat=matrix_inner(world,mo,pot.spin_zeeman_mo);
	Tensor<double_complex> Kmat=matrix_inner(world,mo,pot.K_mo);
	Tensor<double_complex> Jmat=matrix_inner(world,mo,pot.J_mo);
//	print("vnuc, lz, dia, spin-zeeman, kmat, jmat");
//	print(Vnucmat);
//	print(lzmat);
//	print(diamat);
//	print(spin_zeeman_mat);
//	print(Kmat);
//	print(Jmat);

	Tensor<double_complex> vmat=Vnucmat+lzmat+diamat+spin_zeeman_mat-Kmat+Jmat;
	return vmat;
};


std::vector<complex_function_3d>
Znemo::compute_residuals(
		const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi,
		Tensor<double>& eps) const {

	double tol = FunctionDefaults < 3 > ::get_thresh();

    std::vector < std::shared_ptr<real_convolution_3d> > ops(psi.size());
    for (int i=0; i<eps.size(); ++i)
    		ops[i]=std::shared_ptr<real_convolution_3d>(
    				BSHOperatorPtr3D(world, sqrt(-2.*std::min(-0.05,eps(i)+param.shift())), cparam.lo, tol*0.1));

    std::vector<complex_function_3d> tmp = apply(world,ops,-2.0*Vpsi-2.0*param.shift()*psi);
    std::vector<complex_function_3d> res=psi-tmp;

    // update eps
    std::vector<double> norms=norm2s(world,tmp);
    std::vector<double> rnorms=norm2s(world,res);
    if ((world.rank()==0) and (param.printlevel()>1)) {
    	print("norm2(tmp)",norms);
    	print("norm2(res)",rnorms);
    }
    Tensor<double_complex> in=inner(world,Vpsi,res);	// no shift here!
    Tensor<double> delta_eps(eps.size());
    for (int i=0; i<eps.size(); ++i) delta_eps(i)=real(in(i))/(norms[i]*norms[i]);
    eps-=delta_eps;

    if ((world.rank()==0) and (param.printlevel()>1)) {
    	print("orbital energy update",delta_eps);
    }
    truncate(world,res);
    return res;

}


void
Znemo::canonicalize(std::vector<complex_function_3d>& amo,
		std::vector<complex_function_3d>& vnemo,
		XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator>& solver,
		Tensor<double_complex> fock, Tensor<double_complex> ovlp) const {

    Tensor<double_complex> U;
    Tensor<double> evals;
    sygv(fock, ovlp, 1, U, evals);
    // Fix phases.
    for (long i = 0; i < amo.size(); ++i) if (real(U(i, i)) < 0.0) U(_, i).scale(-1.0);

    fock = 0.0;
    for (unsigned int i = 0; i < amo.size(); ++i) fock(i, i) = evals(i);
    amo = transform(world, amo, U);
    vnemo = transform(world, vnemo, U);
    rotate_subspace(U,solver);

}


/// orthonormalize the vectors

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void
Znemo::orthonormalize(std::vector<complex_function_3d>& amo) const {
    normalize(amo);
    double maxq;
    do {
        Tensor<double_complex> Q = Q2(matrix_inner(world, amo, amo)); // Q3(matrix_inner(world, amo_new, amo_new))
        maxq = 0.0;
        for (int j=1; j<Q.dim(0); j++)
            for (int i=0; i<j; i++)
                maxq = std::max(std::abs(Q(j,i)),maxq);
        amo = transform(world, amo, Q, 0.0, true);
        truncate(world, amo);
    } while (maxq>0.01);
    normalize(amo);
}

} // namespace madness
