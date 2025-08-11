/*
 * Nemocomplex.cc
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#include <madness/mra/mra.h>
#include "znemo.h"
#include<madness/chem/diamagneticpotentialfactor.h>
#include <madness/world/test_utilities.h>
#include<madness/chem/masks_and_boxes.h>
#include<madness/chem/MolecularOrbitals.h>



namespace madness {

Znemo::Znemo(World& world, const commandlineparser& parser) : NemoBase(world), mol(world,parser), param(world,parser),
                    calc_param(world,parser) {
	calc_param.set_derived_value("k",9);
	calc_param.read_input_and_commandline_options(world,parser,"dft");

    FunctionDefaults<3>::set_k(calc_param.k());
    FunctionDefaults<3>::set_thresh(calc_param.econv());
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);

    aobasis.read_file(calc_param.aobasis());
//    cparam.set_molecular_info(mol, aobasis, 0);
    calc_param.set_derived_values(mol);

    FunctionDefaults<3>::set_cubic_cell(-calc_param.L(), calc_param.L());
    calc_param.set_derived_value("spin_restricted",false);

	param.set_derived_values();

	print_info=printleveler(param.printlevel());

	if (world.rank()==0 and print_info.print_setup()) {
		param.print("complex","end");
		calc_param.print("dft","end");
		mol.print();
	}

	B={0,0,param.physical_B()};


    // set the linear moments
    for (int i=0; i<3; ++i) rvec.push_back(real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];}));

	coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,calc_param.lo(),calc_param.econv()));
};

bool Znemo::need_recompute_factors_and_potentials(const double thresh) const {
	bool need=NemoBase::need_recompute_factors_and_potentials(thresh);
	if ((not (potentialmanager.get() and potentialmanager->vnuclear().is_initialized()))
	        		or (potentialmanager->vnuclear().thresh()>thresh)) need=true;
    if (not diafac.get()) need=true;
    return need;
}

void Znemo::invalidate_factors_and_potentials() {
    NemoBase::invalidate_factors_and_potentials();
	diafac.reset();
	potentialmanager.reset();
}

void Znemo::recompute_factors_and_potentials(const double thresh) {

	print("recomputing factors and potentials in znemo");
    // compute the nuclei's positions in the "A" space
	diafac.reset(new Diamagnetic_potential_factor(world,param,mol.get_all_coords_vec()));

	coord_3d box_offset={0.0,0.0,0.0};
	if (param.box().size()==5) box_offset={param.box()[2],param.box()[3],param.box()[4]};
	spherical_box<3> sbox2(param.box()[0],param.box()[1],box_offset);
	sbox=real_factory_3d(world).functor(sbox2);

	// the guess is read from a previous nemo calculation
	// make sure the molecule was not reoriented there
	if (not mol.parameters.no_orient()) {
		MADNESS_EXCEPTION("the molecule of the reference calculation was reoriented\n\n",1);
	}

	potentialmanager=std::shared_ptr<PotentialManager>(new PotentialManager(mol, molecule().parameters.core_type()));
	potentialmanager->make_nuclear_potential(world);
    construct_nuclear_correlation_factor(mol, potentialmanager, nemo_param.ncf());

	// currently we cannot have the diamagnetic factor and the nuclear correlation factor
	// at the same time
	if (ncf->type()!=NuclearCorrelationFactor::None and (param.explicit_B()!=0.0)) {
		MADNESS_EXCEPTION("no use of diamagnetic factor and nuclear correlation factor at the same time",1);
	}

    // (re) construct the Poisson solver
	coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,get_calc_param().lo(),get_calc_param().econv()));

}

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
	for (size_t i=0; i<ket_R.size(); ++i) {
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
double Znemo::value(const Tensor<double>& x) {

	if ((x-mol.get_all_coords()).normf()>1.e-12) invalidate_factors_and_potentials();

	// compute the molecular potential
	mol.set_all_coords(x.reshape(mol.natom(),3));
	if (need_recompute_factors_and_potentials(get_calc_param().econv()))
		recompute_factors_and_potentials(get_calc_param().econv());

	get_initial_orbitals();

	iterate();

	save_orbitals("final");
	save_orbitals();

	double energy=analyze();
    output_calc_info_schema(energy);

	return energy;
}

void Znemo::output_calc_info_schema(const double& energy) const {
    nlohmann::json j;
    j["scf_eigenvalues_a"]=tensor_to_json(aeps);
    j["scf_eigenvalues_b"]=tensor_to_json(beps);
    j["model"]="UHF";
    j["B"]=this->B;
    j["driver"]="energy";
    j["return_energy"]=energy;
    update_schema(get_calc_param().prefix()+".calc_info", j);
}

Tensor<double> Znemo::gradient(const Tensor<double>& x) {

	MADNESS_ASSERT(param.explicit_B()==0.0);
	const real_function_3d rhonemo=compute_nemo_density(amo,bmo);
	const Tensor<double> grad=NemoBase::compute_gradient(rhonemo,mol);

    if (world.rank() == 0) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print("  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(" ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (size_t i = 0; i < mol.natom(); ++i) {
            const Atom& atom = mol.get_atom(i);
            printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", int(i),
                   atom.x, atom.y, atom.z, grad[i * 3 + 0], grad[i * 3 + 1],
                   grad[i * 3 + 2]);
        }
    }
    return grad;
}


void Znemo::iterate() {

	double thresh=FunctionDefaults<3>::get_thresh();
	double energy=1.e10;
	double oldenergy=0.0;

	// the diamagnetic box

//	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator<double_complex,3> > solvera(allocator<double_complex,3> (world,amo.size()));
//	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator<double_complex,3> > solverb(allocator<double_complex,3> (world,bmo.size()));
    auto solvera= nonlinear_vector_solver<double_complex,3>(world,amo.size());
    auto solverb= nonlinear_vector_solver<double_complex,3>(world,bmo.size());
	solvera.set_maxsub(get_calc_param().maxsub()); // @suppress("Method cannot be resolved")
	solvera.do_print=(param.printlevel()>2);
	solverb.set_maxsub(get_calc_param().maxsub());
	solverb.do_print=(param.printlevel()>2);

	// set end of iteration cycles for intermediate calculations
	int maxiter = get_calc_param().maxiter();
	double dconv=get_calc_param().dconv();
	double na=1.0,nb=1.0;	// residual norms

	solvera.clear_subspace();
	solverb.clear_subspace();
	bool converged=false;

	// iterate the SCF cycles
	for (int iter=0; iter<maxiter; ++iter) {
//		if (iter%10==0) save_orbitals("iteration"+stringify(iter));
		save_orbitals();

		// compute the density
		real_function_3d density=compute_density(amo,bmo);

		// compute the dual of the orbitals
		std::vector<complex_function_3d> R2amo=make_bra(amo);
		std::vector<complex_function_3d> R2bmo=make_bra(bmo);

		// compute the fock matrix
		std::vector<complex_function_3d> Vnemoa, Vnemob;
		Tensor<double_complex> focka, fockb(0l,0l);
		potentials apot(world,amo.size()), bpot(world,bmo.size());

		apot=compute_potentials(amo, density, amo);		// 4.53
		Vnemoa=apot.vnuc_mo+apot.lz_mo+apot.lz_commutator+
				apot.diamagnetic_mo+apot.spin_zeeman_mo-apot.K_mo+apot.J_mo+apot.zeeman_R_comm;
		truncate(world,Vnemoa,thresh*0.1);
		Kinetic<double_complex,3> T(world);
		focka=T(R2amo,amo) + compute_vmat(amo,apot);		// 0.5??

		if (have_beta()) {
			bpot=compute_potentials(bmo, density, bmo);
			scale(world,bpot.spin_zeeman_mo,-1.0);
			Vnemob=bpot.vnuc_mo+bpot.lz_mo+bpot.lz_commutator+
					bpot.diamagnetic_mo+bpot.spin_zeeman_mo-bpot.K_mo+bpot.J_mo+bpot.zeeman_R_comm;
			truncate(world,Vnemob,thresh*0.1);
			fockb=T(R2bmo,bmo) + compute_vmat(bmo,bpot);
		}

		if (world.rank()==0 and (print_info.print_convergence())) {
			print("Fock matrix");
			print(focka);
			print(fockb);
		}

		oldenergy=energy;
		energy=compute_energy(amo,apot,bmo,bpot,print_info.print_convergence(),false);
//		energy=compute_energy_no_confinement(amo,apot,bmo,bpot,print_info.print_convergence());	// 2.48

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

		if (std::abs(oldenergy-energy)<get_calc_param().econv() and (sqrt(na*na+nb*nb)<dconv)) {
			print("energy converged");
			save_orbitals("converged");
			converged=true;
		}
		if (converged) break;

		// compute the residual of the Greens' function
		double levelshifta=param.shift()-param.physical_B()*0.5;	// account for spin zeeman term
		std::vector<complex_function_3d> resa=compute_residuals(Vnemoa,amo,aeps,levelshifta);	// 2.16

		Tensor<double> normsa=real(inner(world,make_bra(resa),resa));
		na=sqrt(normsa.sumsq());

		std::vector<complex_function_3d> amo_new=solvera.update(amo,resa,0.01,3);
		amo_new=truncate(sbox*amo_new);

		do_step_restriction(amo,amo_new);

		amo=truncate(orthonormalize(amo_new));

		if (have_beta()) {
			double levelshiftb=param.shift()+param.physical_B()*0.5;	// account for spin zeeman term
			std::vector<complex_function_3d> resb=compute_residuals(Vnemob,bmo,beps,levelshiftb);
			truncate(world,resb);
			Tensor<double> normsb=real(inner(world,make_bra(resb),resb));
			nb=sqrt(normsb.sumsq());

			std::vector<complex_function_3d> bmo_new=solverb.update(bmo,resb,0.01,3);
			bmo_new=truncate(sbox*bmo_new);
			do_step_restriction(bmo,bmo_new);
			bmo=truncate(orthonormalize(bmo_new));
		} else {
			nb=0.0;
		}
	}
}


std::vector<complex_function_3d> Znemo::orthonormalize(const std::vector<complex_function_3d>& mo) const {
	std::vector<complex_function_3d> result=copy(world,mo);
	NemoBase::orthonormalize(result,diafac->factor()*R);
	return result;
}

std::vector<complex_function_3d> Znemo::normalize(const std::vector<complex_function_3d>& mo) const {
	std::vector<complex_function_3d> result=copy(world,mo);
	NemoBase::normalize(result,diafac->factor()*R);
	return result;
}


real_function_3d Znemo::compute_density(const std::vector<complex_function_3d>& amo,
		const std::vector<complex_function_3d>& bmo) const {
	return (compute_nemo_density(amo,bmo)*diafac->factor_square()*R_square).truncate();
}

real_function_3d Znemo::compute_spin_density(const std::vector<complex_function_3d>& amo,
		const std::vector<complex_function_3d>& bmo) const {
	return (compute_nemo_spin_density(amo,bmo)*diafac->factor_square()*R_square).truncate();
}

std::vector<complex_function_3d> Znemo::make_bra(const std::vector<complex_function_3d>& rhs) const {
	return truncate((diafac->factor_square()*R_square)*rhs);
}


double Znemo::analyze() {

	if (need_recompute_factors_and_potentials(get_calc_param().econv()))
		recompute_factors_and_potentials(get_calc_param().econv());

	real_function_3d density=compute_density(amo,bmo);
	real_function_3d spindensity=compute_spin_density(amo,bmo);

	potentials apot(world,amo.size()), bpot(world,bmo.size());
	apot=compute_potentials(amo,density,amo);
	if (have_beta()) {
		bpot=compute_potentials(bmo,density,bmo);
		scale(world,bpot.spin_zeeman_mo,-1.0);
	}

	double energy=compute_energy(amo,apot,bmo,bpot,true,true);

	if (world.rank()==0) {
		print("orbital energies alpha",aeps);
		print("orbital energies beta ",beps);
	}

	// compute the current density
	std::vector<real_function_3d> j=compute_current_density(amo,bmo);
	save(j[0],"j0");
	save(j[1],"j1");
	save(j[2],"j2");

	Lz<double_complex,3> lz(world);
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

    print("saving orbitals for plotting");
	save_orbitals("plot");
	save(density,"density");
	save(spindensity,"spindensity");

	auto components=std::vector<std::vector<real_function_3d>>({real(amo),imag(amo),real(bmo),imag(bmo)});
	auto component_names=std::vector<std::string>({"real_amo","imag_amo","real_bmo","imag_bmo"});
	std::vector<real_function_3d> real_aos=SCF::project_ao_basis_only(world, aobasis, mol);
	for (size_t i=0; i<components.size(); ++i) {
		if (world.rank()==0) print("analysis of MO component ",component_names[i]);
		if (components[i].size()>0) SCF::analyze_vectors(world, components[i], real_aos,
			FunctionDefaults<3>::get_thresh()*0.1, molecule(), get_calc_param().print_level(), aobasis);

	}

    return energy;
}


//double Znemo::compute_energy_no_confinement(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
//		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print) const {
//
//	MADNESS_EXCEPTION("Znemo::compute_energy_no_confinement is inaccurate",1);
//	timer tenergy(world,print_info.print_timings());
//    double fac= cparam.spin_restricted() ? 2.0 : 1.0;
//    std::vector<complex_function_3d> dia2amo=make_bra(amo);
//    std::vector<complex_function_3d> dia2bmo=make_bra(bmo);
//    std::vector<complex_function_3d> diaamo=truncate(amo*diafac->factor()*R);
//    std::vector<complex_function_3d> diabmo=truncate(bmo*diafac->factor()*R);
//
//    double_complex kinetic=0.0;
//    for (int axis = 0; axis < 3; axis++) {
//        complex_derivative_3d D = free_space_derivative<double_complex, 3>(world, axis);
//        const std::vector<complex_function_3d> damo = apply(world, D, amo, false);
//        const std::vector<complex_function_3d> dbmo = apply(world, D, bmo);
//        double k1a=norm2(world,R*damo-R*ncf->U1(axis)*amo);
//        double k1b=norm2(world,R*dbmo-R*ncf->U1(axis)*bmo);
//        kinetic+=fac * 0.5* (k1a*k1a + k1b*k1b);
//    }
//
//    const real_function_3d diapot=diafac->bare_diamagnetic_potential();
//    const real_function_3d vnuc=potentialmanager->vnuclear();
//    const complex_function_3d lzcomm=diafac->compute_lz_commutator();
//    const real_function_3d density=compute_density(amo,bmo);
//
//    Lz lz(world,false);
//    double_complex nuclear_potential=inner(density,vnuc);
//    double_complex diamagnetic=inner(density,diapot);
//    double_complex lzval=fac*(inner(world,diaamo,0.5*B[2]*lz(diaamo)).sum()+inner(world,diabmo,0.5*B[2]*lz(diabmo)).sum())
//    		+fac*(inner(world,diaamo,lzcomm*diaamo).sum()+inner(world,diabmo,lzcomm*diabmo).sum());
//    double_complex spin_zeeman=fac*(inner(world,dia2amo,apot.spin_zeeman_mo).sum()+inner(world,dia2bmo,bpot.spin_zeeman_mo).sum());
//    double_complex coulomb=fac*0.5*(inner(world,dia2amo,apot.J_mo).sum()+inner(world,dia2bmo,bpot.J_mo).sum());
//    double_complex exchange=fac*0.5*(inner(world,dia2amo,apot.K_mo).sum()+inner(world,dia2bmo,bpot.K_mo).sum());
//    double_complex zeeman_R_comm=0.0;//fac*(inner(world,dia2amo,apot.zeeman_R_comm).sum()+inner(world,dia2bmo,bpot.zeeman_R_comm).sum());
//
//    double_complex energy=kinetic + nuclear_potential + mol.nuclear_repulsion_energy() +
//    		diamagnetic + lzval + spin_zeeman + coulomb - exchange + zeeman_R_comm;
//
//	if (world.rank()==0 and do_print) {
//		printf("  kinetic energy        %12.8f \n", real(kinetic));
//		printf("  nuclear potential     %12.8f \n", real(nuclear_potential));
//		printf("  nuclear repulsion     %12.8f \n", mol.nuclear_repulsion_energy());
//		printf("  diamagnetic term      %12.8f \n", real(diamagnetic));
//		printf("  orbital zeeman term   %12.8f \n", real(lzval));
//		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
//		printf("  spin zeeman term      %12.8f \n", real(spin_zeeman));
//		printf("  Coulomb               %12.8f \n", real(coulomb));
//		printf("  exchange              %12.8f \n", real(exchange));
//		printf("  total                 %12.8f \n", real(energy));
//	}
//
//	if(fabs(imag(energy))>1.e-7) {
//
//
//		print("real part");
//		printf("  kinetic energy      %12.8f \n", real(kinetic));
//		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
//		printf("  nuclear repulsion   %12.8f \n", mol.nuclear_repulsion_energy());
//		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
//		printf("  orbital zeeman term %12.8f \n", real(lzval));
//		printf("  spin zeeman term    %12.8f \n", real(spin_zeeman));
//		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
//		printf("  Coulomb             %12.8f \n", real(coulomb));
//		printf("  exchange            %12.8f \n", real(exchange));
//		printf("  total               %12.8f \n", real(energy));
//
//		print("imaginary part");
//		printf("  kinetic energy      %12.8f \n", imag(kinetic));
//		printf("  nuclear potential   %12.8f \n", imag(nuclear_potential));
//		printf("  nuclear repulsion   %12.8f \n", 0.0);
//		printf("  diamagnetic term    %12.8f \n", imag(diamagnetic));
//		printf("  orbital zeeman term %12.8f \n", imag(lzval));
//		printf("  orbital zeeman R comm %12.8f \n", imag(zeeman_R_comm));
//		printf("  spin zeeman term    %12.8f \n", imag(spin_zeeman));
//		printf("  Coulomb             %12.8f \n", imag(coulomb));
//		printf("  exchange            %12.8f \n", imag(exchange));
//		printf("  total               %12.8f \n", imag(energy));
//
//		print("imaginary part of the energy",energy.imag());
//
//		MADNESS_EXCEPTION("complex energy computation.. ",1);
//	}
//
//	tenergy.end("compute_energy no confinement");
//	return real(energy);
//}



double Znemo::compute_energy(const std::vector<complex_function_3d>& amo, const Znemo::potentials& apot,
		const std::vector<complex_function_3d>& bmo, const Znemo::potentials& bpot, const bool do_print,
		const bool no_confinement) const {

	timer tenergy(world,print_info.print_timings());
    double fac= get_calc_param().spin_restricted() ? 2.0 : 1.0;
    std::vector<complex_function_3d> dia2amo=make_bra(amo);
    std::vector<complex_function_3d> dia2bmo=make_bra(bmo);


    double_complex kinetic=0.0;
    for (int axis = 0; axis < 3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex, 3>(world, axis);
        const std::vector<complex_function_3d> damo = apply(world, D, amo);
        const std::vector<complex_function_3d> ddia2amo = apply(world, D, dia2amo);
        const std::vector<complex_function_3d> dbmo = apply(world, D, bmo);
        const std::vector<complex_function_3d> ddia2bmo = apply(world, D, dia2bmo);
        kinetic += fac* 0.5 * (inner(world, ddia2amo, damo).sum() + inner(world, ddia2bmo, dbmo).sum());
    }

    double_complex diamagnetic=double_complex(0.0);
    if (no_confinement) {
        const real_function_3d diapot=diafac->bare_diamagnetic_potential();
        const real_function_3d density=compute_density(amo,bmo);
    	diamagnetic=inner(density,diapot);
    } else {
    	diamagnetic=fac*(inner(world,dia2amo,apot.diamagnetic_mo).sum()+inner(world,dia2bmo,bpot.diamagnetic_mo).sum());
    }

    double_complex nuclear_potential=fac*(inner(world,dia2amo,apot.vnuc_mo).sum()+inner(world,dia2bmo,bpot.vnuc_mo).sum());
    double_complex lz=fac*(inner(world,dia2amo,apot.lz_mo).sum()+inner(world,dia2bmo,bpot.lz_mo).sum());
    double_complex lz_comm=fac*(inner(world,dia2amo,apot.lz_commutator).sum()+inner(world,dia2bmo,bpot.lz_commutator).sum());
    double_complex spin_zeeman=fac*(inner(world,dia2amo,apot.spin_zeeman_mo).sum()+inner(world,dia2bmo,bpot.spin_zeeman_mo).sum());
    double_complex coulomb=fac*0.5*(inner(world,dia2amo,apot.J_mo).sum()+inner(world,dia2bmo,bpot.J_mo).sum());
    double_complex exchange=fac*0.5*(inner(world,dia2amo,apot.K_mo).sum()+inner(world,dia2bmo,bpot.K_mo).sum());
    double_complex zeeman_R_comm=fac*(inner(world,dia2amo,apot.zeeman_R_comm).sum()+inner(world,dia2bmo,bpot.zeeman_R_comm).sum());

    double_complex energy=kinetic + nuclear_potential + mol.nuclear_repulsion_energy() +
    		diamagnetic + lz + lz_comm + spin_zeeman + coulomb - exchange + zeeman_R_comm;

	if (world.rank()==0 and do_print) {
		printf("  kinetic energy        %12.8f \n", real(kinetic));
		printf("  nuclear potential     %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion     %12.8f \n", mol.nuclear_repulsion_energy());
		printf("  diamagnetic term      %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term   %12.8f \n", real(lz));
		printf("  orbital zeeman comm   %12.8f \n", real(lz_comm));
		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
		printf("  spin zeeman term      %12.8f \n", real(spin_zeeman));
		printf("  Coulomb               %12.8f \n", real(coulomb));
		printf("  exchange              %12.8f \n", real(exchange));
		printf("  total                 %12.8f \n", real(energy));
	}

	if(fabs(imag(energy))>1.e-7) {


		print("real part");
		printf("  kinetic energy      %12.8f \n", real(kinetic));
		printf("  nuclear potential   %12.8f \n", real(nuclear_potential));
		printf("  nuclear repulsion   %12.8f \n", mol.nuclear_repulsion_energy());
		printf("  diamagnetic term    %12.8f \n", real(diamagnetic));
		printf("  orbital zeeman term %12.8f \n", real(lz));
		printf("  orbital zeeman comm %12.8f \n", real(lz_comm));
		printf("  orbital zeeman R comm %12.8f \n", real(zeeman_R_comm));
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
		printf("  orbital zeeman R comm %12.8f \n", imag(zeeman_R_comm));
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

/// following Lazzeretti, J. Mol. Struct, 313 (1994)

/// see also the first three terms of Eq. (47) of Hodge, Am. J. Phys., 82 (2014)
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

//	double ndiff=(spin_density*R_square).trace();
//	print("nalpha - nbeta",ndiff);
//	double nel=(density*R_square).trace();
//	print("nelectron",nel);

	std::vector<real_function_3d> vspin_density=zero_functions_compressed<double,3>(world,3);;
	vspin_density[2]=spin_density;

	// compute first contribution to current density from orbitals: Re (\Psi \hat p \Psi)
	// p = -i del
	// Re(-i z) + Re( -i a - i^2 b) = b  =Im (z)
	std::vector<real_function_3d> j=zero_functions_compressed<double,3>(world,3);
	for (auto& mo : alpha_mo) j+=imag(conj(mo)*grad(mo));
	for (auto& mo : beta_mo) j+=imag(conj(mo)*grad(mo));
//
//	real_function_3d null1=div(j*R_square);
//	double n31=null1.norm2();
//	print("div(j)",n31);

    // compute the magnetic potential and density contribution: A \rho
	std::vector<real_function_3d> A=compute_magnetic_vector_potential(world,B);
//	auto vAdens=A*density;
	j+=A*density;
//	real_function_3d null2=div(vAdens*R_square);
//	double n32=null2.norm2();
//	print("div(A*dens)",n32);

	// spin density contribution: rot(\Psi \hat S \Psi)
	j+=0.5*rot(vspin_density);
	j[0]+=0.5*2.0*(-1.0)*ncf->U1(1)*spin_density;
	j[1]-=0.5*2.0*(-1.0)*ncf->U1(0)*spin_density;

//	real_function_3d null3=div(j*R_square);
//	double n33=null3.norm2();
//	print("div(j)",n33);

//	save(j[0],"j0reg");
//	save(j[1],"j1reg");
//	save(j[2],"j2reg");

//	auto divj=div(j);
//	save(divj,"div_reg_j");
//	auto term1=2.0*dot(world,ncf->U1vec(),j);
//	save(term1,"U1term");

	// sanity check
	real_function_3d null=R_square*(-2.0*dot(world,ncf->U1vec(),j) + div(j,true));
	save(null,"null");
	double n3=null.norm2();
	print("final div(j)",n3);

	j=truncate(j*R_square);
//	auto jnorm=norm2s(world,j);
//	print("norm(j)",jnorm);

	return j;
}


void Znemo::test_compute_current_density() const {

	Lz<double_complex,3> lz(world);
	complex_function_3d pp=complex_factory_3d(world).functor(p_orbital(1,1.0,{0,0,0}));
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

void Znemo::test_landau_wave_function() {
	recompute_factors_and_potentials(get_calc_param().econv());

	int m=1;
	double B=param.physical_B();
	complex_function_3d f=complex_factory_3d(world).functor(landau_wave_function(m,-B));
	double fnorm=f.norm2();
	f.scale(1.0/fnorm);

	Lz<double_complex,3> lz(world);
	double_complex expval=0.5*B*lz(f,f);
	print("<f(m) | Lz | f(m)> ", m, expval);
	double_complex diaexpval=inner(f,diafac->bare_diamagnetic_potential()*f);
	print("<f(m) | dia | f(m)>", m, diaexpval);
	Kinetic<double_complex,3> T(world);
	double_complex kinexpval=T(f,f);
	print("<f(m) | T | f(m)>  ", m, kinexpval);

	print("total energy",expval + diaexpval +kinexpval);
	print("expect some contribution from the confinement in z-direction");



	MADNESS_EXCEPTION("end of test",1);
}



void Znemo::do_step_restriction(const std::vector<complex_function_3d>& mo,
		std::vector<complex_function_3d>& mo_new) const {
    PROFILE_MEMBER_FUNC(SCF);
    std::vector<double> anorm = norm2s(world, sub(world, mo, mo_new));
    int nres = 0;
    for (unsigned int i = 0; i < mo.size(); ++i) {
        if (anorm[i] > get_calc_param().maxrotn()) {
            double s = get_calc_param().maxrotn() / anorm[i];
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


void Znemo::get_initial_orbitals() {

	/* possible guesses
	 *
	 * 1) guess_functions as specified in the user input
	 * 2) custom guess (l,ml,exponenent) for atoms as specified in the user input
	 * 3) reference file from a previous znemo calculation
	 * 4) ao restart data from a previous znemo calculations
	 * 5) restartdata file from a previous moldft calculation
	 * 6) hcore guess
	 */

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	bool gotit=false;

	if (param.guess_functions().size()>0) {
		zmos=read_explicit_guess_functions();
		gotit=true;
	} else if (param.guess().size()>0) {
		zmos=custom_guess();
		gotit=true;
	} else {

		try {
			zmos=read_reference();
			gotit=true;
		} catch (...) {
			print("could not read reference orbitals");
		}

//		zmos=read_complex_guess();
		if (not gotit) {
			try {
				zmos=read_real_guess();
				print("successful read real guess orbitals");
				gotit=true;
			} catch (...) {
				print("could not read real guess orbitals");
			}
		}
		if (not gotit) {
			try {
				zmos=read_restartaodata();
				gotit=true;
			} catch (...) {
				print("could not read ao restart data");
			}
		}
		if (not gotit) {
			try {
				zmos=hcore_guess();
				gotit=true;
			} catch (...) {
				print("could not do an hcore guess");
			}
		}
	}
	if (not gotit) {
		MADNESS_EXCEPTION("unable to generate guess functions",1);
	}

	amo.clear();
	bmo.clear();
	for (int i=0; i<get_calc_param().nalpha(); ++i) amo.push_back(zmos.first.get_mos()[i]);
	for (int i=0; i<get_calc_param().nbeta(); ++i) bmo.push_back(zmos.second.get_mos()[i]);
	aeps=Tensor<double>(amo.size());
	beps=Tensor<double>(bmo.size());

	amo=orthonormalize(amo);
	bmo=orthonormalize(bmo);
	set_thresh(world,amo,FunctionDefaults<3>::get_thresh());
	set_thresh(world,bmo,FunctionDefaults<3>::get_thresh());

}


/// read the guess orbitals from a previous nemo or moldft calculation
std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::read_real_guess() const {

	print("reading mos from restartdata");
	auto mos=MolecularOrbitals<double,3>::read_restartdata(world, "restartdata", molecule(), get_calc_param().nalpha(), get_calc_param().nbeta());
	MolecularOrbitals<double,3> amo=mos.first,bmo=mos.second;

	// check for correct polynomial order
	auto k1=amo.get_mos()[0].get_impl()->get_k();
	auto k2=FunctionDefaults<3>::get_k();
	MADNESS_CHECK_THROW(k1==k2,"reference polynomial order is not equal to the current polynomial order");

	// check if beta orbitals are present if they should be
	if (get_calc_param().have_beta() and bmo.size()==0) {
		print("reading from a closed-shell reference -> replicating alpha to beta orbitals");
		std::vector<real_function_3d> tmp;
		for (int i=0; i<get_calc_param().nmo_beta(); ++i) tmp.push_back(amo.get_mos()[i]);
		bmo.set_mos(tmp);
	}

    // confine the orbitals to an approximate Gaussian form corresponding to the
    // diamagnetic (harmonic) potential
    coord_3d remaining_B=B-coord_3d{0,0,param.explicit_B()};
    real_function_3d gauss=diafac->custom_factor(remaining_B,diafac->get_v(),1.0);

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	auto inv=ncf->inverse();
	zmos.first.set_mos(truncate(convert<double,double_complex,3>(world,amo.get_mos()*inv)));		// alpha
	zmos.second.set_mos(truncate(convert<double,double_complex,3>(world,bmo.get_mos()*inv)));		// beta

	return zmos;
}


/// read the guess orbitals from a previous nemo or moldft calculation
std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::read_complex_guess() const {
	print("reading mos from complex restartdata");
	auto zmos=MolecularOrbitals<double_complex,3>::read_restartdata(world, "restartdata", molecule(), get_calc_param().nalpha(), get_calc_param().nbeta());
	return zmos;
}

/// read a list of functions for the guess
std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::read_explicit_guess_functions() const {
	// loop over all provided functions
        if (param.guess_functions().size()<size_t(std::max(get_calc_param().nalpha(),get_calc_param().nbeta()))) {
		MADNESS_EXCEPTION("you did not provide enough guess functions for this calculations",1);
	}
	std::vector<complex_function_3d> guess_vector;
	for (auto filename : param.guess_functions()) {
		complex_function_3d mo=complex_factory_3d(world);
		print("loading file",filename,"as guess function");
		load(mo,filename);
		guess_vector.push_back(mo);
	}

	std::vector<complex_function_3d> amos,bmos;
	for (int i=0; i<get_calc_param().nalpha(); ++i) amos.push_back(guess_vector[i]);
	for (int i=0; i<get_calc_param().nbeta(); ++i) bmos.push_back(guess_vector[i]);

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	zmos.first.set_mos(amos);
	zmos.second.set_mos(bmos);

	return zmos;
}

/// read a list of functions for the guess
std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::read_restartaodata() const {
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;

	print("reading orbitals ao projection from complex calculation");
	zmos=MolecularOrbitals<double_complex,3>::read_restartaodata(world, mol, get_calc_param().have_beta());
//		print("  .. failed ");
//		try {
//			print("reading orbitals ao projection from real calculation");
//			auto mos=MolecularOrbitals<double,3>::read_restartaodata(world, mol, get_calc_param().have_beta());
//			zmos.first.set_mos(convert<double,double_complex,3>(world,mos.first.get_mos()));
//			zmos.second.set_mos(convert<double,double_complex,3>(world,mos.second.get_mos()));
//
//		} catch (...) {
//			print("  .. failed ");
//		}
	return zmos;
}

/// read a list of functions for the guess
std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::read_reference() const {
	std::string name="reference";
	print("reading orbitals from file",name);

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;

	archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str(), 1);
	std::size_t namo=0, nbmo=0;

	std::vector<complex_function_3d> amos,bmos;
	Tensor<double> aeps,beps;

	ar & namo & nbmo & aeps & beps;
	amos.resize(namo);
	bmos.resize(nbmo);
	for (auto& a: amos) ar & a;
	for (auto& a: bmos) ar & a;
	zmos.first.set_mos(amos);
	zmos.second.set_mos(bmos);

	return zmos;
}


std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::custom_guess() const {
	MADNESS_ASSERT(molecule().natom()==1);

	Lz<double_complex,3> lz(world);

	std::vector<std::string> guess=param.guess();
	std::vector<complex_function_3d> guess_vector;
	for (size_t i=0; i<guess.size(); i+=3) {
		long l=atol(guess[i].c_str());
		long ml=atol(guess[i+1].c_str());
		double exponent=atof(guess[i+2].c_str());
		print("making guess with l,ml,e",l,ml,exponent);

		complex_function_3d g;
		coord_3d origin=molecule().get_atoms()[0].get_coords();
		if (l==0) g=complex_factory_3d(world).functor(s_orbital(exponent,origin));
		if (l==1) g=complex_factory_3d(world).functor(p_orbital(ml,exponent,origin));
		double norm=g.norm2();
		g.scale(1/norm);
		double_complex ii=inner(g,lz(g));
		print("<p|Lz|p>, ml",ii, ml);

		guess_vector.push_back(g);
	}

	std::vector<complex_function_3d> amos,bmos;
	amos.resize(get_calc_param().nalpha());
	bmos.resize(get_calc_param().nbeta());
	std::copy(guess_vector.begin(),guess_vector.begin()+get_calc_param().nalpha(),amos.begin());
	std::copy(guess_vector.begin(),guess_vector.begin()+get_calc_param().nbeta(),bmos.begin());

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	zmos.first.set_mos(amos);
	zmos.second.set_mos(bmos);

	return zmos;
}

std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
Znemo::hcore_guess() const {

	print("computing mos from hcore guess");

	// get the AOs from the basis set ..
	std::vector<real_function_3d> real_aos=SCF::project_ao_basis_only(world, aobasis, mol);
	std::vector<complex_function_3d> aos = convert<double,double_complex,3>(world,real_aos);

	// .. and from the Landau wave function (the eigenfunction of the free particle in the B-field)
    if (param.physical_B()>0.0) {
		for (int n=0; n<2; ++n) {
//			coord_3d origin=mol.get_atom(i).get_coords();
			coord_3d origin;
			double B=param.physical_B();
			complex_function_3d f=complex_factory_3d(world).functor(landau_wave_function(n,-B,origin));
			double fnorm=f.norm2();
			f.scale(1.0/fnorm);

			aos.push_back(conj(f));
		}
    }

	// compute the Hamilton matrix

	// nuclear and diamagnetic potential
	real_function_3d vlocal=(potentialmanager->vnuclear()+diafac->bare_diamagnetic_potential()).truncate();
	Tensor<double_complex> potential=matrix_inner(world,aos,vlocal*aos);

	// zeeman term
	Lz<double_complex,3> lz_operator(world,false);
	Tensor<double_complex> lzmat=0.5*B[2]*lz_operator(aos,aos);
	potential+=lzmat;

	// kinetic term
	Kinetic<double_complex,3> T(world);
	Tensor<double_complex>  kinetic=T(aos,aos);

	Tensor<double_complex> fock = convert<double_complex>(kinetic) + potential;
	fock = 0.5 * (fock + conj_transpose(fock));

	// diagonalize the Fock matrix
	Tensor<double_complex> c;
	Tensor<double> e;
	Tensor<double_complex> overlap =matrix_inner(world, aos, aos, true);
	sygvp(world, fock, overlap, 1, c, e);

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> > zmos;
	std::vector<complex_function_3d> amo = transform(world, aos, c(_, Slice(0, get_calc_param().nmo_alpha() - 1)));
	amo=truncate(normalize(amo));
	zmos.first.set_mos(truncate(amo*ncf->inverse()));		// alpha

	if (get_calc_param().have_beta()) {
		std::vector<complex_function_3d>
		bmo = transform(world, aos, c(_, Slice(0, get_calc_param().nmo_beta() - 1)));
		bmo=truncate(normalize(bmo));
		zmos.second.set_mos(truncate(bmo*ncf->inverse()));		// beta
	}

	return zmos;
}

/// compute the potential operators applied on the orbitals
Znemo::potentials Znemo::compute_potentials(const std::vector<complex_function_3d>& mo,
		const real_function_3d& density,
		const std::vector<complex_function_3d>& rhs) const {

	timer potential_timer(world,print_info.print_timings());
	std::vector<complex_function_3d> dia2mo=make_bra(mo);

	// prepare exchange operator
	Exchange<double_complex,3> K(world,get_calc_param().lo());
	Tensor<double> occ(mo.size());
	occ=1.0;
    K.set_bra_and_ket(conj(world, dia2mo), mo);

	Nuclear<double_complex,3> nuc(world,ncf);

	potentials pot(world,rhs.size());

	pot.J_mo=(*coulop)(density)*rhs;
	pot.K_mo=K(rhs);
	pot.vnuc_mo=nuc(rhs);

	pot.lz_mo=0.5*B[2]*Lz<double_complex,3>(world,false)(rhs);
	pot.lz_commutator=diafac->compute_lz_commutator()*rhs;

	pot.diamagnetic_mo=diafac->apply_potential(rhs);

	MADNESS_ASSERT(diafac->B_along_z(B));
	pot.spin_zeeman_mo=B[2]*0.5*rhs;
	pot.zeeman_R_comm = double_complex(0.0,0.5)*B[2]*(rvec[0]* ncf->U1(1)-rvec[1]* ncf->U1(0))*rhs;

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

	std::vector<complex_function_3d> dia2mo=make_bra(mo);
	truncate(world,dia2mo);
	Tensor<double_complex> Vnucmat=matrix_inner(world,dia2mo,pot.vnuc_mo);
	Tensor<double_complex> lzmat=matrix_inner(world,dia2mo,pot.lz_mo);
	Tensor<double_complex> lz_comm_mat=matrix_inner(world,dia2mo,pot.lz_commutator);
	Tensor<double_complex> diamat=matrix_inner(world,dia2mo,pot.diamagnetic_mo);
	Tensor<double_complex> spin_zeeman_mat=matrix_inner(world,dia2mo,pot.spin_zeeman_mo);
	Tensor<double_complex> Kmat=matrix_inner(world,dia2mo,pot.K_mo);
	Tensor<double_complex> Jmat=matrix_inner(world,dia2mo,pot.J_mo);
	Tensor<double_complex> ZRcom_mat=matrix_inner(world,dia2mo,pot.zeeman_R_comm);

	Tensor<double_complex> vmat=Vnucmat+lzmat+lz_comm_mat+diamat+spin_zeeman_mat-Kmat+Jmat+ZRcom_mat;
	return vmat;
};


std::vector<complex_function_3d>
Znemo::compute_residuals(
		const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi,
		Tensor<double>& eps,
		double levelshift) const {
	timer residual_timer(world,print_info.print_timings());

	double tol = FunctionDefaults < 3 > ::get_thresh();

    std::vector < std::shared_ptr<real_convolution_3d> > ops(psi.size());
    for (int i=0; i<eps.size(); ++i) {
    	ops[i]=std::shared_ptr<real_convolution_3d>(
    			BSHOperatorPtr3D(world, sqrt(-2.*std::min(-0.05,eps(i)+levelshift)),
    					get_calc_param().lo()*0.001, tol*0.001));
    }

    std::vector<complex_function_3d> tmp = apply(world,ops,-2.0*Vpsi-2.0*levelshift*psi);

    const std::vector<complex_function_3d> res=truncate(psi-tmp);
    const std::vector<complex_function_3d> bra_res=make_bra(res);

    // update eps
    Tensor<double> norms=real(inner(world,make_bra(tmp),tmp));
    Tensor<double_complex> rnorms=inner(world,bra_res,res);
    for (int i=0; i<norms.size(); ++i) {
    	norms(i)=sqrt(norms(i));
    	rnorms(i)=sqrt(rnorms(i));
    }
    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("norm2(tmp)",norms);
    	print("norm2(res)",rnorms);
    }
    Tensor<double_complex> in=inner(world,Vpsi,bra_res);	// no shift here!
    Tensor<double> delta_eps(eps.size());
    for (int i=0; i<eps.size(); ++i) delta_eps(i)=real(in(i))/(norms[i]*norms[i]);
    eps-=delta_eps;

    if ((world.rank()==0) and (print_info.print_convergence())) {
    	print("orbital energy update",delta_eps);
    }
    residual_timer.end("compute_residuals");
    return res;

}


void
Znemo::canonicalize(std::vector<complex_function_3d>& amo,
		std::vector<complex_function_3d>& vnemo,
		potentials& pot,
		XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, vector_function_allocator<double_complex,3> >& solver,
		Tensor<double_complex> fock, Tensor<double_complex> ovlp) const {

    Tensor<double_complex> U;
    Tensor<double> evals;
    sygv(fock, ovlp, 1, U, evals);
    // Fix phases.
    for (size_t i = 0; i < amo.size(); ++i) if (real(U(i, i)) < 0.0) U(_, i).scale(-1.0);

    fock = 0.0;
    for (unsigned int i = 0; i < amo.size(); ++i) fock(i, i) = evals(i);
    amo = transform(world, amo, U);
    vnemo = transform(world, vnemo, U);
    rotate_subspace(U,solver);
    pot.transform(U);

}

void Znemo::save_orbitals(std::string suffix) const
{
	suffix = "_" + suffix;
	const real_function_3d &dia = diafac->factor();
	for (size_t i = 0; i < amo.size(); ++i)
		save(amo[i], "amo" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(bmo[i], "bmo" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(real(amo[i]), "amo_real" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(imag(amo[i]), "amo_imag" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(real(bmo[i]), "bmo_real" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(imag(bmo[i]), "bmo_imag" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(madness::abs(amo[i]), "absamo" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(madness::abs(bmo[i]), "absbmo" + stringify(i) + suffix);
	for (size_t i = 0; i < amo.size(); ++i)
		save(madness::abs(amo[i] * dia), "diaamo" + stringify(i) + suffix);
	for (size_t i = 0; i < bmo.size(); ++i)
		save(madness::abs(bmo[i] * dia), "diabmo" + stringify(i) + suffix);
}

} // namespace madness
