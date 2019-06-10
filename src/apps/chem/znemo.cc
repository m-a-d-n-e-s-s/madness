/*
 * Nemocomplex.cc
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#include <madness/mra/mra.h>
#include "znemo.h"
#include <chem/diamagneticpotentialfactor.h>
#include <chem/test_utilities.h>
#include <chem/masks_and_boxes.h>


namespace madness {

Znemo::Znemo(World& w) : world(w), param(world), molecule("input"), cparam() {
	cparam.read_file("input");

    FunctionDefaults<3>::set_k(cparam.k);
    FunctionDefaults<3>::set_thresh(cparam.econv);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-cparam.L, cparam.L);

    aobasis.read_file(cparam.aobasis);
    cparam.set_molecular_info(molecule, aobasis, 0);
	param.set_derived_values();

	print_info=printleveler(param.printlevel());

	if (world.rank()==0 and print_info.print_setup()) {
		param.print(param.params,"complex","end");
		cparam.print(world);
		molecule.print();
	}

	B={0,0,param.physical_B()};

    // compute the nuclei's positions in the "A" space
	diafac.reset(new Diamagnetic_potential_factor(world,param,molecule.get_all_coords_vec()));

	coord_3d box_offset={0.0,0.0,0.0};
	if (param.box().size()==5) box_offset={param.box()[2],param.box()[3],param.box()[4]};
	spherical_box<3> sbox2(param.box()[0],param.box()[1],box_offset);
	sbox=real_factory_3d(world).functor(sbox2);

	// the guess is read from a previous nemo calculation
	// make sure the molecule was not reoriented there
	if (not cparam.no_orient) {
		MADNESS_EXCEPTION("the molecule of the reference calculation was reoriented\n\n",1);
	}

	coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,cparam.lo,cparam.econv));
};

bool Znemo::test() const {
	bool success=true;
	try {
		success=success and diafac->test_me(1);
		success=success and test_U_potentials();
	} catch (...) {
		print("tests failed..");
//		MADNESS_EXCEPTION("tests failed",1);
	}
	return success;
}

/// test the identity <F| f (T + Vdia ) f |F> = <F|f^2 (T + Udia) |F>
bool Znemo::test_U_potentials() const {

	test_output t("entering Znemo::test_U_potentials .. ");

	std::vector<complex_function_3d> ket=diafac->make_fake_orbitals(2,{0.0,1.0,0.5});
	std::vector<complex_function_3d> bra=diafac->make_fake_orbitals(2,{0.1,-0.1,0.3});
	std::vector<complex_function_3d> ket_R2=diafac->factor_square()*ket;
	std::vector<complex_function_3d> bra_R2=diafac->factor_square()*bra;
	std::vector<complex_function_3d> ket_R=diafac->factor()*ket;
	std::vector<complex_function_3d> bra_R=diafac->factor()*bra;

	double thresh=FunctionDefaults<3>::get_thresh();
	Kinetic<double_complex,3> T(world);

	double_complex t1,t2;
	for (int i=0; i<ket_R.size(); ++i) {
		t1+=T(bra_R[i],ket_R[i]);
		t2+=T(bra_R2[i],ket[i]);
	}
	double_complex vdia=inner(bra_R,diafac->bare_diamagnetic_potential()*ket_R);
	double_complex udia=inner(bra_R2,diafac->apply_potential(ket));

	t.logger << "<amo | R T R | amo>    "<< t1 << std::endl;
	t.logger << "<amo | R2 T  | amo>    "<< t2 << std::endl;
	t.logger << "<amo | R2 V_dia | amo> "<< vdia<< std::endl;
	t.logger << "<amo | R2 U_dia | amo> "<< udia<< std::endl;

	double error=std::abs(t1 + vdia - t2 - udia);
	t.logger << "error  "<< error<< std::endl;
	t.logger << "Note: if all individual potentials are correct and precise, check \n"
			<< "diamagnetic_height and box_softness as the most probable source of imprecision " << std::endl;

	return t.end(error<thresh);
}


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

		orthonormalize(amo);
		orthonormalize(bmo);
	}

	test();

	iterate();

	// final energy computation
	real_function_3d density=sum(world,abssq(world,amo));
	if (have_beta()) density+=sum(world,abssq(world,bmo));
	if (cparam.spin_restricted) density*=2.0;
	density=density*diafac->factor_square();
	density.truncate();

	potentials apot(world,amo.size()), bpot(world,bmo.size());
	apot=compute_potentials(amo,density,amo);
	if (have_beta()) {
		bpot=compute_potentials(bmo,density,bmo);
		scale(world,bpot.spin_zeeman_mo,-1.0);
	}

	double energy=compute_energy(amo,apot,bmo,bpot,true);

	if (world.rank()==0) {
		print("orbital energies alpha",aeps);
		print("orbital energies beta ",beps);
	}

	save_orbitals("final");
	save_orbitals();

	analyze();

	return energy;
}


void Znemo::iterate() {

	double thresh=FunctionDefaults<3>::get_thresh();
	double energy=1.e10;
	double oldenergy=0.0;

	// the diamagnetic box

	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solvera(allocator(world,amo.size()));
	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solverb(allocator(world,bmo.size()));
	solvera.set_maxsub(cparam.maxsub); // @suppress("Method cannot be resolved")
	solvera.do_print=(param.printlevel()>2);
	solverb.set_maxsub(cparam.maxsub);
	solverb.do_print=(param.printlevel()>2);

	// set end of iteration cycles for intermediate calculations
	int maxiter = cparam.maxiter;
	double dconv=cparam.dconv;
	double na=1.0,nb=1.0;	// residual norms

	solvera.clear_subspace();
	solverb.clear_subspace();
	bool converged=false;

	// iterate the SCF cycles
	for (int iter=0; iter<maxiter; ++iter) {
//		if (iter%10==0)
			save_orbitals("iteration"+stringify(iter));

		// compute the density
		real_function_3d density=sum(world,abssq(world,amo));
		if (have_beta()) density+=sum(world,abssq(world,bmo));
		if (cparam.spin_restricted) density*=2.0;
		density=density*diafac->factor_square();
		density.truncate();

		// compute the dual of the orbitals
		std::vector<complex_function_3d> R2amo=diafac->factor_square()*amo;
		std::vector<complex_function_3d> R2bmo=diafac->factor_square()*bmo;
		truncate(world,R2amo);
		truncate(world,R2bmo);

		// compute the fock matrix
		std::vector<complex_function_3d> Vnemoa, Vnemob;
		Tensor<double_complex> focka, fockb(0l,0l);
		potentials apot(world,amo.size()), bpot(world,bmo.size());

		apot=compute_potentials(amo, density, amo);
		Vnemoa=apot.vnuc_mo+apot.lz_mo+apot.lz_commutator+
				apot.diamagnetic_mo+apot.spin_zeeman_mo-apot.K_mo+apot.J_mo;
		truncate(world,Vnemoa,thresh*0.1);
		Kinetic<double_complex,3> T(world);
		focka=T(R2amo,amo) + compute_vmat(amo,apot);

		if (have_beta()) {
			bpot=compute_potentials(bmo, density, bmo);
			scale(world,bpot.spin_zeeman_mo,-1.0);
			Vnemob=bpot.vnuc_mo+bpot.lz_mo+bpot.lz_commutator+
					bpot.diamagnetic_mo+bpot.spin_zeeman_mo-bpot.K_mo+bpot.J_mo;
			truncate(world,Vnemob,thresh*0.1);
			fockb=T(R2bmo,bmo) + compute_vmat(bmo,bpot);
		}

		if (world.rank()==0 and (print_info.print_convergence())) {
			print("Fock matrix");
			print(focka);
			print(fockb);
		}

		oldenergy=energy;
		energy=compute_energy(amo,apot,bmo,bpot,print_info.print_convergence());
		energy=compute_energy_no_confinement(amo,apot,bmo,bpot,print_info.print_convergence());


		Tensor<double_complex> ovlp=matrix_inner(world,R2amo,amo);

		canonicalize(amo,Vnemoa,apot,solvera,focka,ovlp);
		truncate(world,Vnemoa);
		if (have_beta()) {
			Tensor<double_complex> ovlp=matrix_inner(world,R2bmo,bmo);
			canonicalize(bmo,Vnemob,bpot,solverb,fockb,ovlp);
		}

		if (world.rank()==0 and (print_info.print_convergence())) print("using fock matrix for the orbital energies");
		for (int i=0; i<focka.dim(0); ++i) aeps(i)=real(focka(i,i));
		for (int i=0; i<fockb.dim(0); ++i) beps(i)=real(fockb(i,i));


		if (world.rank()==0 and (print_info.print_convergence())) {
			print("orbital energies alpha",aeps);
			print("orbital energies beta ",beps);
		}
		if ((world.rank() == 0) and print_info.print_energies()) {
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

		Tensor<double> normsa=real(inner(world,resa*diafac->factor_square(),resa));
		na=sqrt(normsa.sumsq());

		std::vector<complex_function_3d> amo_new=solvera.update(amo,resa,0.01,3);
		amo_new=sbox*amo_new;

		do_step_restriction(amo,amo_new);

		amo=amo_new;
		truncate(world,amo);
		orthonormalize(amo);

		if (have_beta()) {
			std::vector<complex_function_3d> resb=compute_residuals(Vnemob,bmo,beps);
			truncate(world,resb);
			Tensor<double> normsb=real(inner(world,diafac->factor_square()*resb,resb));
			nb=sqrt(normsb.sumsq());

			std::vector<complex_function_3d> bmo_new=solverb.update(bmo,resb,0.01,3);
			bmo_new=sbox*bmo_new;
			do_step_restriction(bmo,bmo_new);
			bmo=bmo_new;
			truncate(world,bmo);
			orthonormalize(bmo);
		} else {
			nb=0.0;
		}
	}
}

void Znemo::analyze() const {

	// compute the current density
	std::vector<real_function_3d> j=compute_current_density(amo,bmo);
	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

	Lz lz(world);
	// compute the expectation values of the Lz operator
	std::vector<complex_function_3d> lzamo=0.5*B[2]*lz(amo);
	std::vector<complex_function_3d> lzbmo=0.5*B[2]*lz(bmo);
	Tensor<double_complex> lza_exp=inner(world,amo,lzamo);
	Tensor<double_complex> lzb_exp=inner(world,bmo,lzbmo);
	print("< amo | lz | amo >",lza_exp);
	print("< bmo | lz | bmo >",lzb_exp);

    std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,B);

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


double Znemo::compute_energy_no_confinement(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print) const {

	timer tenergy(world,print_info.print_timings());
    double fac= cparam.spin_restricted ? 2.0 : 1.0;
    std::vector<complex_function_3d> dia2amo=amo*diafac->factor_square();
    std::vector<complex_function_3d> dia2bmo=bmo*diafac->factor_square();
    std::vector<complex_function_3d> diaamo=amo*diafac->factor();
    std::vector<complex_function_3d> diabmo=bmo*diafac->factor();
    truncate(world,dia2amo,FunctionDefaults<3>::get_thresh());
    truncate(world,dia2bmo,FunctionDefaults<3>::get_thresh());


    double_complex kinetic=0.0;
    for (int axis = 0; axis < 3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex, 3>(world, axis);
        const std::vector<complex_function_3d> damo = apply(world, D, diaamo);
        const std::vector<complex_function_3d> dbmo = apply(world, D, diabmo);
        kinetic += fac* 0.5 * (inner(world, damo, damo).sum() + inner(world, dbmo, dbmo).sum());
    }

    real_function_3d diapot=diafac->bare_diamagnetic_potential();
    complex_function_3d lzcomm=diafac->compute_lz_commutator();

    Lz lz(world,false);
    double_complex nuclear_potential=fac*(inner(world,dia2amo,apot.vnuc_mo).sum()+inner(world,dia2bmo,bpot.vnuc_mo).sum());
    double_complex diamagnetic=fac*(inner(world,diaamo,diapot*diaamo).sum()+inner(world,diabmo,diapot*diabmo).sum());
    double_complex lzval=fac*(inner(world,diaamo,0.5*B[2]*lz(diaamo)).sum()+inner(world,diabmo,0.5*B[2]*lz(diabmo)).sum())
    		+fac*(inner(world,diaamo,lzcomm*diaamo).sum()+inner(world,diabmo,lzcomm*diabmo).sum());
    double_complex spin_zeeman=fac*(inner(world,dia2amo,apot.spin_zeeman_mo).sum()+inner(world,dia2bmo,bpot.spin_zeeman_mo).sum());
    double_complex coulomb=fac*0.5*(inner(world,dia2amo,apot.J_mo).sum()+inner(world,dia2bmo,bpot.J_mo).sum());
    double_complex exchange=fac*0.5*(inner(world,dia2amo,apot.K_mo).sum()+inner(world,dia2bmo,bpot.K_mo).sum());

    double_complex energy=kinetic + nuclear_potential + molecule.nuclear_repulsion_energy() +
    		diamagnetic + lzval + spin_zeeman + coulomb - exchange;

	if (world.rank()==0 and do_print) {
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", molecule.nuclear_repulsion_energy());
		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", real(lzval));
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
		printf("  orbital zeeman term %12.8f \n", real(lzval));
		printf("  spin zeeman term    %12.8f \n", real(spin_zeeman));
		printf("  Coulomb             %12.8f \n", real(coulomb));
		printf("  exchange            %12.8f \n", real(exchange));
		printf("  total               %12.8f \n", real(energy));

		print("imaginary part");
		printf("  kinetic energy      %12.8f \n", imag(kinetic));
		printf("  nuclear potential   %12.8f \n", imag(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", 0.0);
		printf("  diamagnetic term    %12.8f \n", imag(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", imag(lzval));
		printf("  spin zeeman term    %12.8f \n", imag(spin_zeeman));
		printf("  Coulomb             %12.8f \n", imag(coulomb));
		printf("  exchange            %12.8f \n", imag(exchange));
		printf("  total               %12.8f \n", imag(energy));

		print("imaginary part of the energy",energy.imag());

		MADNESS_EXCEPTION("complex energy computation.. ",1);
	}

	tenergy.end("compute_energy no confinement");
	return real(energy);
}



double Znemo::compute_energy(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print) const {

	timer tenergy(world,print_info.print_timings());
    double fac= cparam.spin_restricted ? 2.0 : 1.0;
    std::vector<complex_function_3d> dia2amo=amo*diafac->factor_square();
    std::vector<complex_function_3d> dia2bmo=bmo*diafac->factor_square();
    truncate(world,dia2amo);
    truncate(world,dia2bmo);


    double_complex kinetic=0.0;
    for (int axis = 0; axis < 3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex, 3>(world, axis);
        const std::vector<complex_function_3d> damo = apply(world, D, amo);
        const std::vector<complex_function_3d> ddia2amo = apply(world, D, dia2amo);
        const std::vector<complex_function_3d> dbmo = apply(world, D, bmo);
        const std::vector<complex_function_3d> ddia2bmo = apply(world, D, dia2bmo);
        kinetic += fac* 0.5 * (inner(world, ddia2amo, damo).sum() + inner(world, ddia2bmo, dbmo).sum());
    }

    double_complex nuclear_potential=fac*(inner(world,dia2amo,apot.vnuc_mo).sum()+inner(world,dia2bmo,bpot.vnuc_mo).sum());
    double_complex diamagnetic=fac*(inner(world,dia2amo,apot.diamagnetic_mo).sum()+inner(world,dia2bmo,bpot.diamagnetic_mo).sum());
    double_complex lz=fac*(inner(world,dia2amo,apot.lz_mo).sum()+inner(world,dia2bmo,bpot.lz_mo).sum());
    double_complex lz_comm=fac*(inner(world,dia2amo,apot.lz_commutator).sum()+inner(world,dia2bmo,bpot.lz_commutator).sum());
    double_complex spin_zeeman=fac*(inner(world,dia2amo,apot.spin_zeeman_mo).sum()+inner(world,dia2bmo,bpot.spin_zeeman_mo).sum());
    double_complex coulomb=fac*0.5*(inner(world,dia2amo,apot.J_mo).sum()+inner(world,dia2bmo,bpot.J_mo).sum());
    double_complex exchange=fac*0.5*(inner(world,dia2amo,apot.K_mo).sum()+inner(world,dia2bmo,bpot.K_mo).sum());

    double_complex energy=kinetic + nuclear_potential + molecule.nuclear_repulsion_energy() +
    		diamagnetic + lz + lz_comm + spin_zeeman + coulomb - exchange;

	if (world.rank()==0 and do_print) {
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", molecule.nuclear_repulsion_energy());
		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", real(lz));
		printf("  orbital zeeman comm %12.8f \n", real(lz_comm));
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
		printf("  orbital zeeman comm %12.8f \n", real(lz_comm));
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
		printf("  orbital zeeman comm %12.8f \n", imag(lz_comm));
		printf("  spin zeeman term    %12.8f \n", imag(spin_zeeman));
		printf("  Coulomb             %12.8f \n", imag(coulomb));
		printf("  exchange            %12.8f \n", imag(exchange));
		printf("  total               %12.8f \n", imag(energy));

		print("imaginary part of the energy",energy.imag());

		MADNESS_EXCEPTION("complex energy computation.. ",1);
	}

	tenergy.end("compute_energy");
	return real(energy);
}

/// following Lazeretti, J. Mol. Struct, 313 (1994)
std::vector<real_function_3d> Znemo::compute_current_density(
		const std::vector<complex_function_3d>& alpha_mo,
		const std::vector<complex_function_3d>& beta_mo) const {

	// compute density and spin density
	real_function_3d adens=real_factory_3d(world);
	real_function_3d bdens=real_factory_3d(world);
	for (auto& mo : alpha_mo) adens+=abs_square(mo);
	for (auto& mo : beta_mo) bdens+=abs_square(mo);
	real_function_3d density=adens+bdens;
	real_function_3d spin_density=adens-bdens;

	std::vector<real_function_3d> vspin_density=zero_functions_compressed<double,3>(world,3);;
	vspin_density[2]=spin_density;

	// compute first contribution to current density from orbitals:
	// psi^* p psi = i psi^* del psi
	std::vector<complex_function_3d> j=zero_functions_compressed<double_complex,3>(world,3);
	for (auto& mo : alpha_mo) j+=double_complex(0.0,1.0)*(conj(mo)*grad(mo));
	for (auto& mo : beta_mo) j+=double_complex(0.0,1.0)*(conj(mo)*grad(mo));


    // compute the magnetic potential
	std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,B);

	// compute density contribution and spin contribution
	j-=convert<double,double_complex,3>(world,A*density);
	j+=convert<double,double_complex,3>(world,0.5*rot(vspin_density));

	std::vector<real_function_3d>  realj=real(j);

	// sanity check

	real_function_3d null=div(realj);
	double n3=null.norm2();
	print("div(j)",n3);
//	MADNESS_ASSERT(n3<FunctionDefaults<3>::get_thresh());


	return realj;
}


void Znemo::test_compute_current_density() const {

	Lz lz(world);
	complex_function_3d pp=complex_factory_3d(world).f(p_plus);
	double norm=pp.norm2();
	pp.scale(1/norm);
	double_complex l=inner(amo[0],0.5*B[2]*lz(amo[0]));
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
            if (world.rank() == 0 and print_info.print_convergence()) {
                if (nres == 1)
                    printf("  restricting step ");
                printf(" %d", i);
            }
            mo_new[i].gaxpy(s, mo[i], 1.0 - s, false);
        }
    }
    if (nres > 0 && (world.rank() == 0) and print_info.print_convergence())
        printf("\n");

    world.gop.fence();
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

    // confine the orbitals to an approximate Gaussian form corresponding to the
    // diamagnetic (harmonic) potential
    coord_3d remaining_B=B-coord_3d{0,0,param.explicit_B()};
    real_function_3d gauss=diafac->custom_factor(remaining_B,diafac->get_v(),1.0);
//    complex_function_3d pp=complex_factory_3d(world).f(p_plus);
//    save(real(pp),"p_plus");
//    complex_function_3d gauss=diafac->factor_with_phase(remaining_B,diafac->get_v());
//    save(real(gauss),"gauss");
//    return gauss*real_mo;
    return convert<double,double_complex,3>(world,real_mo*gauss);
}


/// compute the potential operators applied on the orbitals
Znemo::potentials Znemo::compute_potentials(const std::vector<complex_function_3d>& mo,
		const real_function_3d& density,
		std::vector<complex_function_3d>& rhs) const {

	timer potential_timer(world,print_info.print_timings());
	std::vector<complex_function_3d> dia2mo=mo*diafac->factor_square();
	truncate(world,dia2mo);

	// prepare exchange operator
	Exchange<double_complex,3> K=Exchange<double_complex,3>(world);
	Tensor<double> occ(mo.size());
	occ=1.0;
	K.set_parameters(conj(world,dia2mo),mo,occ,cparam.lo,cparam.econv);

	potentials pot(world,rhs.size());

	pot.J_mo=(*coulop)(density)*rhs;
	pot.K_mo=K(rhs);
	pot.vnuc_mo=vnuclear*rhs;

	pot.lz_mo=0.5*B[2]*Lz(world,false)(rhs);
	pot.lz_commutator=diafac->compute_lz_commutator()*rhs;

	pot.diamagnetic_mo=diafac->apply_potential(rhs);

	MADNESS_ASSERT(diafac->B_along_z(B));
	pot.spin_zeeman_mo=B[2]*0.5*rhs;

	truncate(world,pot.J_mo);
	truncate(world,pot.K_mo);
	truncate(world,pot.vnuc_mo);
	truncate(world,pot.lz_mo);
	truncate(world,pot.diamagnetic_mo);
	truncate(world,pot.spin_zeeman_mo);
	potential_timer.end("compute_potentials");
	return pot;
};


Tensor<double_complex> Znemo::compute_vmat(const std::vector<complex_function_3d>& mo,
		const potentials& pot) const {

	std::vector<complex_function_3d> dia2mo=mo*diafac->factor_square();
	truncate(world,dia2mo);
	Tensor<double_complex> Vnucmat=matrix_inner(world,dia2mo,pot.vnuc_mo);
	Tensor<double_complex> lzmat=matrix_inner(world,dia2mo,pot.lz_mo);
	Tensor<double_complex> lz_comm_mat=matrix_inner(world,dia2mo,pot.lz_commutator);
	Tensor<double_complex> diamat=matrix_inner(world,dia2mo,pot.diamagnetic_mo);
	Tensor<double_complex> spin_zeeman_mat=matrix_inner(world,dia2mo,pot.spin_zeeman_mo);
	Tensor<double_complex> Kmat=matrix_inner(world,dia2mo,pot.K_mo);
	Tensor<double_complex> Jmat=matrix_inner(world,dia2mo,pot.J_mo);

	Tensor<double_complex> vmat=Vnucmat+lzmat+lz_comm_mat+diamat+spin_zeeman_mat-Kmat+Jmat;
	return vmat;
};


std::vector<complex_function_3d>
Znemo::compute_residuals(
		const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi,
		Tensor<double>& eps) const {
	timer residual_timer(world,print_info.print_timings());

	double tol = FunctionDefaults < 3 > ::get_thresh();

    std::vector < std::shared_ptr<real_convolution_3d> > ops(psi.size());
    for (int i=0; i<eps.size(); ++i) {
    	ops[i]=std::shared_ptr<real_convolution_3d>(
    			BSHOperatorPtr3D(world, sqrt(-2.*std::min(-0.05,eps(i)+param.shift())),
    					cparam.lo*0.001, tol*0.001));
    }

    std::vector<complex_function_3d> tmp = apply(world,ops,-2.0*Vpsi-2.0*param.shift()*psi);

    std::vector<complex_function_3d> res=psi-tmp;

    // update eps
    Tensor<double> norms=real(inner(world,tmp*diafac->factor_square(),tmp));
    Tensor<double_complex> rnorms=inner(world,res*diafac->factor_square(),res);
    for (int i=0; i<norms.size(); ++i) {
    	norms(i)=sqrt(norms(i));
    	rnorms(i)=sqrt(rnorms(i));
    }
    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("norm2(tmp)",norms);
    	print("norm2(res)",rnorms);
    }
    Tensor<double_complex> in=inner(world,Vpsi*diafac->factor_square(),res);	// no shift here!
    Tensor<double> delta_eps(eps.size());
    for (int i=0; i<eps.size(); ++i) delta_eps(i)=real(in(i))/(norms[i]*norms[i]);
    eps-=delta_eps;

    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("orbital energy update",delta_eps);
    }
    truncate(world,res);
    residual_timer.end("compute_residuals");
    return res;

}


void
Znemo::canonicalize(std::vector<complex_function_3d>& amo,
		std::vector<complex_function_3d>& vnemo,
		potentials& pot,
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
    pot.transform(U);

}


void Znemo::normalize(std::vector<complex_function_3d>& mo) const {
	std::vector<complex_function_3d> dia2mo=mo*diafac->factor_square();
	Tensor<double_complex> n=inner(world,mo,dia2mo);
	for (int i=0; i<mo.size(); ++i) mo[i].scale(1.0/sqrt(real(n[i])));
}


/// orthonormalize the vectors

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void
Znemo::orthonormalize(std::vector<complex_function_3d>& amo) const {
	if (amo.size()==0) return;
    normalize(amo);
    double maxq;
    do {
    	std::vector<complex_function_3d> dia2mo=amo*diafac->factor_square();
        Tensor<double_complex> Q = Q2(matrix_inner(world, dia2mo, amo)); // Q3(matrix_inner(world, amo_new, amo_new))
        maxq = 0.0;
        for (int j=1; j<Q.dim(0); j++)
            for (int i=0; i<j; i++)
                maxq = std::max(std::abs(Q(j,i)),maxq);
        amo = transform(world, amo, Q, 0.0, true);
        truncate(world, amo);
    } while (maxq>0.01);
    normalize(amo);
}

void Znemo::save_orbitals(std::string suffix) const {
	suffix="_"+suffix;
	const real_function_3d& dia=diafac->factor();
	for (int i=0; i<amo.size(); ++i) save(amo[i],"amo"+stringify(i)+suffix);
	for (int i=0; i<bmo.size(); ++i) save(bmo[i],"bmo"+stringify(i)+suffix);
	for (int i=0; i<amo.size(); ++i) save(madness::abs(amo[i]),"absamo"+stringify(i)+suffix);
	for (int i=0; i<bmo.size(); ++i) save(madness::abs(bmo[i]),"absbmo"+stringify(i)+suffix);
	for (int i=0; i<amo.size(); ++i) save(madness::abs(amo[i]*dia),"diaamo"+stringify(i)+suffix);
	for (int i=0; i<bmo.size(); ++i) save(madness::abs(bmo[i]*dia),"diabmo"+stringify(i)+suffix);

}

} // namespace madness
