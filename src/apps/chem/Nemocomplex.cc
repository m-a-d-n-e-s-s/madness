/*
 * Nemocomplex.cc
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#include "Nemocomplex.h"
#include <madness/mra/mra.h>


namespace madness {

/// compute the molecular energy
double Nemo_complex::value() {

	// compute the molecular potential
	const Molecule& m=molecule;
	auto molecular_potential = [m] (const coord_3d& r) {
		return m.nuclear_attraction_potential(r[0],r[1],r[2]);};
	vnuclear=real_factory_3d(world).functor(molecular_potential).thresh(FunctionDefaults<3>::get_thresh()*0.1);
	vnuclear.set_thresh(FunctionDefaults<3>::get_thresh());

	// read the guess orbitals
	amo=read_guess("alpha");
	if (have_beta()) bmo=read_guess("beta");
	aeps=Tensor<double>(amo.size());
	beps=Tensor<double>(bmo.size());

	double thresh=FunctionDefaults<3>::get_thresh();
	double energy=1.e10;
	double oldenergy=0.0;
	double two_electron_alpha, two_electron_beta;

	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solvera(allocator(world,amo.size()));
	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solverb(allocator(world,bmo.size()));
	solvera.set_maxsub(10);
	solvera.do_print=(param.printlevel()>2);
	solverb.set_maxsub(10);
	solverb.do_print=(param.printlevel()>2);

	// increase the magnetic field
	for (int i=0; i<param.B().size(); ++i) {
		B=param.B()[i];
		print("solving for magnetic field B=",B);

		// set end of iteration cycles for intermediate calculations
		int maxiter = (i==param.B().size()-1) ? cparam.maxiter : 20;
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
			std::vector<complex_function_3d> Vnemo, lznemo, dianemo, Knemo, Jnemo, Vnemoa, Vnemob,
				spin_zeeman_nemo;
			Tensor<double_complex> focka, fockb(0l,0l), tmata, tmatb, vmata, vmatb;

			compute_potentials(amo, density, Vnemo, lznemo, dianemo,spin_zeeman_nemo, Knemo, Jnemo);
			vmata=compute_vmat(amo,Vnemo,lznemo,dianemo,spin_zeeman_nemo,Knemo,Jnemo);
			Vnemoa=Vnemo+lznemo+dianemo+spin_zeeman_nemo-Knemo+Jnemo;
			truncate(world,Vnemoa,thresh*0.1);
			two_electron_alpha= real(-inner(world,Knemo,amo).sum() + inner(world,Jnemo,amo).sum());

			if (have_beta()) {
				compute_potentials(bmo, density, Vnemo, lznemo, dianemo,spin_zeeman_nemo, Knemo, Jnemo);
				scale(world,spin_zeeman_nemo,-1.0);
				vmatb=compute_vmat(bmo,Vnemo,lznemo,dianemo,spin_zeeman_nemo,Knemo,Jnemo);
				Vnemob=Vnemo+lznemo+dianemo+spin_zeeman_nemo-Knemo+Jnemo;
				truncate(world,Vnemob);
				two_electron_beta= real(-inner(world,Knemo,bmo).sum() + inner(world,Jnemo,bmo).sum());
			}

			Kinetic<double_complex,3> T(world);
			focka=T(amo,amo) + vmata;
			if (have_beta()) fockb=T(bmo,bmo) + vmatb;

			if (world.rank()==0 and (param.printlevel()>1)) {
				print("Fock matrix");
				print(focka);
				print(fockb);
			}

			Tensor<double_complex> ovlp=matrix_inner(world,amo,amo);
			canonicalize(amo,Vnemoa,focka,ovlp);

			if (have_beta()) {
				Tensor<double_complex> ovlp=matrix_inner(world,bmo,bmo);
				canonicalize(bmo,Vnemob,fockb,ovlp);
			}

			// compute orbital and total energies
			oldenergy=energy;
			energy=aeps.sum() + beps.sum();
			energy=energy-0.5*(two_electron_alpha + two_electron_beta) + molecule.nuclear_repulsion_energy();
//			if (iter<5) {
				for (int i=0; i<focka.dim(0); ++i) aeps(i)=real(focka(i,i));
				for (int i=0; i<fockb.dim(0); ++i) beps(i)=real(fockb(i,i));
//			}
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
			amo=solvera.update(amo,resa,0.01,3);
			amo=sbox*amo;
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
				bmo=solverb.update(bmo,resb,0.01,3);
				bmo=sbox*bmo;
//				bmo=orthonormalize_symmetric(bmo);
				truncate(world,bmo);
			}
			save_orbitals(iter);
		}
		if (world.rank()==0) {
			print("orbital energies alpha",aeps);
			print("orbital energies beta ",beps);

			Tensor<double> oza=0.5*B*real(inner(world,amo,Lz(amo)));
			print("Orbital Zeeman term alpha ",oza);
			Tensor<double> bza=0.5*B*real(inner(world,bmo,Lz(bmo)));
			print("Orbital Zeeman term beta  ",bza);
			Tensor<double> diaa=0.125*B*B*real(inner(world,amo,diamagnetic()*amo));
			print("diamagnetic term alpha    ",diaa);
			Tensor<double> diab= (have_beta()) ? 0.125*B*B*real(inner(world,bmo,diamagnetic()*bmo)) : Tensor<double>();
			print("diamagnetic term beta     ",diab);
		}

	}

	return energy;
}

/// compute the action of the Lz =i r x del operator on rhs
std::vector<complex_function_3d> Nemo_complex::Lz(const std::vector<complex_function_3d>& rhs) const {
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
std::vector<complex_function_3d> Nemo_complex::read_guess(const std::string& spin) const {

	int nmo= (spin=="alpha") ? cparam.nalpha : cparam.nbeta;
	std::vector<real_function_3d> real_mo=zero_functions<double,3>(world,nmo);

	// load the converged orbitals
    for (std::size_t imo = 0; imo < nmo; ++imo) {
    	print("loading mos ",spin,imo);
    	load(real_mo[imo], "nemo" + stringify(imo));
    }
    return convert<double,double_complex,3>(world,real_mo);
}

/// compute the potential operators applied on the orbitals
void Nemo_complex::compute_potentials(const std::vector<complex_function_3d>& mo,
		real_function_3d& density,
		std::vector<complex_function_3d>& Vnemo,
		std::vector<complex_function_3d>& lznemo,
		std::vector<complex_function_3d>& dianemo,
		std::vector<complex_function_3d>& spin_zeeman_nemo,
		std::vector<complex_function_3d>& Knemo,
		std::vector<complex_function_3d>& Jnemo) const {
	Vnemo=vnuclear*mo;
	lznemo=0.5*B*Lz(mo);
	dianemo=0.125*B*B*diamagnetic()*mo;
	spin_zeeman_nemo=B*0.5*mo;
	Exchange<double_complex,3> K=Exchange<double_complex,3>(world);
	Tensor<double> occ(mo.size());
	occ=1.0;
	K.set_parameters(conj(world,mo),mo,occ,cparam.lo,cparam.econv);
	Knemo=K(mo);
	Jnemo=(*coulop)(density)*mo;

	truncate(world,Vnemo);
	truncate(world,lznemo);
	truncate(world,dianemo);
	truncate(world,Knemo);
	truncate(world,Jnemo);


};


Tensor<double_complex>
Nemo_complex::compute_vmat(
		const std::vector<complex_function_3d>& mo,
		const std::vector<complex_function_3d>& Vnemo,
		const std::vector<complex_function_3d>& lznemo,
		const std::vector<complex_function_3d>& dianemo,
		const std::vector<complex_function_3d>& spin_zeeman_nemo,
		const std::vector<complex_function_3d>& Knemo,
		const std::vector<complex_function_3d>& Jnemo) const {

	Tensor<double_complex> Vnucmat=matrix_inner(world,mo,Vnemo);
	Tensor<double_complex> lzmat=matrix_inner(world,mo,lznemo);
	Tensor<double_complex> diamat=matrix_inner(world,mo,dianemo);
	Tensor<double_complex> spin_zeeman_mat=matrix_inner(world,mo,spin_zeeman_nemo);
	Tensor<double_complex> Kmat=matrix_inner(world,mo,Knemo);
	Tensor<double_complex> Jmat=matrix_inner(world,mo,Jnemo);

	Tensor<double_complex> vmat=Vnucmat+lzmat+diamat+spin_zeeman_mat-Kmat+Jmat;
	return vmat;
};


std::vector<complex_function_3d>
Nemo_complex::compute_residuals(
		const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi,
		Tensor<double>& eps) const {

	double tol = FunctionDefaults < 3 > ::get_thresh();

    std::vector < std::shared_ptr<real_convolution_3d> > ops(psi.size());
    for (int i=0; i<eps.size(); ++i)
    		ops[i]=std::shared_ptr<real_convolution_3d>(
    				BSHOperatorPtr3D(world, sqrt(-2.*std::min(-0.05,eps(i))), cparam.lo, tol*0.1));
    std::vector<complex_function_3d> tmp = apply(world,ops,-2.0*Vpsi);
    std::vector<complex_function_3d> res=psi-tmp;

    // update eps
    std::vector<double> norms=norm2s(world,tmp);
    std::vector<double> rnorms=norm2s(world,res);
    if ((world.rank()==0) and (param.printlevel()>1)) {
    	print("norm2(tmp)",norms);
    	print("norm2(res)",rnorms);
    }
    Tensor<double_complex> in=inner(world,Vpsi,res);
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
Nemo_complex::canonicalize(std::vector<complex_function_3d>& amo,
		std::vector<complex_function_3d>& vnemo,
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

}


/// orthonormalize the vectors

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void
Nemo_complex::orthonormalize(std::vector<complex_function_3d>& amo) const {
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
