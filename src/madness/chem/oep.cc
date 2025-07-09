/*
 * oep.cpp
 *
 *  Created on: Nov 6, 2019
 *      Author: fbischoff
 */

#include<madness/chem/oep.h>
#include<madness/chem/BSHApply.h>
#include <madness/world/test_utilities.h>



namespace madness {

/// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
/// for other functionals, slater potential must be modified
/// HF orbitals and eigenvalues are used as the guess here
/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps
/// converged if norm, total energy difference and orbital energy differences (if not OAEP) are converged
double OEP::solve(const vecfuncT& HF_nemo1) {

	// recompute HF Fock matrix and orbitals
	auto [HF_Fock, HF_nemo] = recompute_HF(HF_nemo1);

	// compute Slater potential Vs and average IHF from HF orbitals and eigenvalues
	const real_function_3d Vs = compute_slater_potential(HF_nemo);
	if (oep_param.saving_amount() >= 1) save(Vs, "Slaterpotential");
	save(Vs, "Slaterpotential");

	tensorT KS_Fock=copy(HF_Fock);
    real_function_3d Voep = copy(Vs);
	if (oep_param.restart()) {
		load_restartdata(KS_Fock);
	}

	// deep copy KS_nemo MOs
	vecfuncT KS_nemo = copy(world,calc->amo);

	for (std::string model : oep_param.model())
		results=iterate(model,HF_nemo,HF_Fock,KS_nemo,KS_Fock,Voep,Vs);

	print("KS_Fock after convergence");
	print(KS_Fock);
	auto [eval, evec] = syev(KS_Fock);
	calc->aeps=eval;
	calc->amo=KS_nemo;

	Vfinal=copy(Voep);
    save_restartdata(KS_Fock);

	printf("      +++ FINAL TOTAL %s ENERGY = %15.8f  Eh +++\n\n\n", oep_param.model().back().c_str(), results.Econv);
	return results.Econv;
}

void OEP::output_calc_info_schema(const double& energy) const {
    nlohmann::json j;
    j["scf_eigenvalues_a"]=tensor_to_json(calc->aeps);
    j["model"]=oep_param.model().back();
    j["driver"]="energy";
    j["return_energy"]=energy;
    update_schema(get_calc_param().prefix()+".oep_calc_info", j);
}


void OEP::save_restartdata(const Tensor<double>& fock) const {
	if (world.rank()==0) print("saving OEP orbitals to file restartdata_OEP");
    MolecularOrbitals<double,3> mo;
    mo.update_mos_and_eps(get_calc()->get_amo(),get_calc()->aeps);
    mo.update_occ(get_calc()->aocc);
    mo.recompute_irreps(get_calc()->param.pointgroup(),R_square);

    projector_irrep p=projector_irrep(calc->param.pointgroup())
            .set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(false);
    auto Vfinal1=p(Vfinal)[0];

    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, "restartdata_OEP");
    ar & mo & fock & Vfinal1;
//    MolecularOrbitals<double,3> amo=to_MO();
//    archive::ParallelOutputArchive ar(world, "restartdata_OEP");
//	ar & amo & fock & Vfinal;
}

void OEP::load_restartdata(Tensor<double>& fock) {
	if (world.rank()==0) print("loading OEP orbitals from file restartdata_OEP");
	archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, "restartdata_OEP");
	MolecularOrbitals<double,3> mo;
	ar & mo & fock & Vfinal;
	mo.pretty_print("OEP MOs from file");
	get_calc()->amo=mo.get_mos();
    get_calc()->aeps=mo.get_eps();
    get_calc()->aset=mo.get_localize_sets();
    get_calc()->aocc=mo.get_occ();
//    MolecularOrbitals<double,3> amo;
//    ar & amo & fock & Vfinal;
    projector_irrep p=projector_irrep(calc->param.pointgroup())
            .set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(false);
    Vfinal=p(Vfinal)[0];
//    calc->amo=amo.get_mos();
}

std::tuple<Tensor<double>, vecfuncT> OEP::recompute_HF(const vecfuncT& HF_nemo) const {

	timer timer1(world);
    const vecfuncT R2nemo=truncate(R_square*HF_nemo);
    Fock<double,3> F(world,get_reference().get());
    Tensor<double> HF_Fock=F(R2nemo,HF_nemo);
	vecfuncT nemo_new=HF_nemo;
	timer1.end("recompute HF");
	return std::make_tuple(HF_Fock, nemo_new);
}

/// the OEP Fock operator is the HF Fock operator without exchange but with the OEP
std::shared_ptr<Fock<double,3>> OEP::make_fock_operator() const {
    Fock<double,3> fock(world,get_reference().get());
    MADNESS_CHECK(fock.remove_operator(("K")));
    LocalPotentialOperator<double,3> Voep(world,"Voep",Vfinal);
    fock.add_operator("Voep",std::make_shared<LocalPotentialOperator<double,3> >(Voep));
    return std::make_shared<Fock<double,3>>(fock);
}

OEPResults OEP::compute_and_print_final_energies(const std::string model, const real_function_3d& Voep,
		const vecfuncT& KS_nemo, const tensorT& KS_Fock,
		const vecfuncT& HF_nemo, const tensorT& HF_Fock) const {

	// print final orbital energies
	auto [KS_eigvals, evec1] = syev(KS_Fock);

	print("final", model, "canonical orbital energies (no level alignment included):");
	print_orbens(KS_eigvals);

    auto [eval1, evec2] = syev(HF_Fock);
    double homoHF = eval1.max();
    double homoKS = KS_eigvals.max();
    print("canonical HF HOMO energy",homoHF);
    print("canonical KS HOMO energy",homoKS);

    // final Jnemo and Knemo have to be computed again in order to calculate final energy
	vecfuncT Jnemo, Knemo, Knemo_HF;
	compute_coulomb_potential(KS_nemo, Jnemo);
	compute_exchange_potential(KS_nemo, Knemo);

	Exchange<double,3> K(world,get_calc_param().lo());
    K.set_bra_and_ket(R_square * HF_nemo, HF_nemo);
	double Ex_HF=-inner(R_square*HF_nemo,K(HF_nemo));

	// compute final exchange energy using different methods and final kinetic energy
	double Ex_vir = compute_exchange_energy_vir(KS_nemo, Voep);
	double Ex_conv = compute_exchange_energy_conv(R_square*KS_nemo, Knemo);
//	double Ex_HF = compute_exchange_energy_conv(R_square*HF_nemo, Knemo_HF);
	double Ekin_KS = compute_kinetic_energy(KS_nemo); // like Ts in Ospadov_2017, equation (22)
	double Ekin_HF = compute_kinetic_energy(HF_nemo); // like T in Ospadov_2017, equation (22)
	double Tc = Ekin_HF - Ekin_KS; // like Tc = T - Ts in Ospadov_2017, equation (24)

	real_function_3d rho_KS = compute_density(KS_nemo);
	real_function_3d rho_HF = compute_density(HF_nemo);
	double Drho = compute_delta_rho(rho_HF, rho_KS);


	print("FINAL", model, "ENERGY Evir:");
	//double Evir = compute_energy(KS_nemo, Ex_vir)[0];

	print("FINAL", model, "ENERGY Econv:");
	double Econv = compute_energy(KS_nemo, Ex_conv)[0];


	printf("     Ex_vir       = %15.8f  Eh", Ex_vir);
	printf("\n     Ex_conv      = %15.8f  Eh", Ex_conv);
	printf("\n     Ex_HF        = %15.8f  Eh\n", Ex_HF);

	printf("\n     Ekin_HF (T)  = %15.8f  Eh", Ekin_HF);
	printf("\n     Ekin_KS (Ts) = %15.8f  Eh\n", Ekin_KS);

	printf("\n     DEvir_14     = %15.8f mEh", (Ex_vir - Ex_conv)*1000.0); // like in Kohut_2014, equation (45)
	printf("\n     DEvir_17     = %15.8f mEh\n", (Ex_vir - Ex_HF - 2.0*Tc)*1000.0); // like in Ospadov_2017, equation (28)
	print("     Drho         =     ", Drho, "e\n\n");

	if (not get_calc_param().do_localize()) {
		print("---------------------------------------------------------------------------");
		double E_0 = compute_E_zeroth(KS_eigvals);
		double E_1 = compute_E_first(R*KS_nemo, R*Jnemo, R*Knemo, Voep);

		printf("  E^(0)               = %15.8f  Eh", E_0);
		printf("\n  E^(1)               = %15.8f  Eh", E_1);
		printf("\n  E^(0) + E^(1)       = %15.8f  Eh", E_0 + E_1);
		printf("\n  difference to Econv = %15.8f mEh\n\n", (E_0 + E_1 - Econv)*1000.0);
	}
	Tensor<double> f_pp = compute_fock_diagonal_elements(calc->aeps, KS_nemo, Knemo, Voep);

//	print("KS Fock matrix elements using Voep/virial ", calc->aeps);
	print("KS Fock matrix elements using K operator  ", f_pp);

	OEPResults results;
	results.devir14= (Ex_vir - Ex_conv); // like in Kohut_2014, equation (45)
	results.devir17= (Ex_vir - Ex_HF - 2.0*Tc); // like in Ospadov_2017, equation (28)
	results.drho=Drho;
	results.Ex_conv=Econv;
	results.Ex_vir=Ex_vir;
	results.Ex_HF=Ex_HF;
	results.E_kin_HF=Ekin_HF;
	results.E_kin_KS=Ekin_KS;
	results.model=model;
	results.Econv=Econv;

	return results;
}

OEPResults OEP::iterate(const std::string model, const vecfuncT& HF_nemo, const tensorT& HF_Fock,
		vecfuncT& KS_nemo, tensorT& KS_Fock, real_function_3d& Voep, const real_function_3d Vs) const {

	bool print_debug=(get_calc_param().print_level()>=10) and (world.rank()==0);

	// compute the constant HF contributions to the OEP hierarchy
	const real_function_3d ocep_numerator_HF=-1.0*compute_energy_weighted_density_local(HF_nemo,HF_Fock);
	const real_function_3d dcep_numerator_HF=compute_total_kinetic_density(HF_nemo);
	const real_function_3d mrks_numerator_HF=compute_Pauli_kinetic_density(HF_nemo);
	if (oep_param.saving_amount() >= 1) {
		save(ocep_numerator_HF,"ocep_numerator_HF");
		save(dcep_numerator_HF,"dcep_numerator_HF");
		save(mrks_numerator_HF,"mrks_numerator_HF");
	}


//	typedef allocator<double, 3> allocT;
//	typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
//	allocT alloc(world, KS_nemo.size());
//	solverT solver(allocT(world, KS_nemo.size()),param.print_level()>4);
    auto solver= nonlinear_vector_solver<double,3>(world,KS_nemo.size());
	solver.set_maxsub(get_calc_param().maxsub());

	double energy=0.0;
	std::vector<double> energies(1,0.0);
	bool converged=false;

	timer timer1(world,get_calc_param().print_level()>=3);
	for (size_t iter = 0; iter < oep_param.maxiter(); ++iter) {

	    if (get_calc_param().do_localize()) {
	    	for (size_t i=0; i<KS_nemo.size(); ++i) calc->aeps(i)=KS_Fock(i,i);
	    	KS_nemo=localize(KS_nemo,get_calc_param().econv(),iter==0);
	    	if (get_calc_param().print_level()>=10)
	    		SCF::analyze_vectors(world, KS_nemo, calc->ao, calc->vtol, calc->molecule, get_calc_param().print_level(),
	    			calc->aobasis, calc->aocc, tensorT(), calc->aset );
	    	// calc->analyze_vectors(world,KS_nemo,calc->aocc,tensorT(),calc->aset);
	    }
	    if (do_symmetry()) {
		    std::vector<std::string> str_irreps;
	    	KS_nemo=symmetry_projector(KS_nemo,R_square,str_irreps);
            if (world.rank()==0) print("orbital irreps",str_irreps);
	    }

	    // compute parts of the KS Fock matrix J, Unuc and Voep
		vecfuncT Jnemo, Unemo, Fnemo;
		compute_nemo_potentials(KS_nemo, Jnemo, Unemo);
		vecfuncT R2KS_nemo = truncate(R_square*KS_nemo);

		timer1.tag("compute potentials");

		Voep=copy(Vs);

		// treat ocep correction separately as it depends on the KS Fock matrix
		tensorT Fock_ocep;
		if (need_ocep_correction(model)) {
			real_function_3d ocep_correction = compute_ocep_correction(ocep_numerator_HF, HF_nemo,KS_nemo,HF_Fock,KS_Fock);
			Fock_ocep=matrix_inner(world,R2KS_nemo,ocep_correction*KS_nemo);
			if (oep_param.save_iter_corrections()>0 and (iter%oep_param.save_iter_corrections()==0))
				save(ocep_correction,"ocep_correction"+std::to_string(iter));
		}

		if (need_dcep_correction(model)) {
			real_function_3d dcep_correction=compute_dcep_correction(dcep_numerator_HF, HF_nemo,KS_nemo);
			Voep += dcep_correction;
			if (oep_param.save_iter_corrections()>0 and (iter%oep_param.save_iter_corrections()==0))
				save(dcep_correction,"dcep_correction"+std::to_string(iter));
		}

		if (need_mrks_correction(model)) {
			real_function_3d mrks_correction=compute_mrks_correction(mrks_numerator_HF, HF_nemo,KS_nemo);
			Voep += mrks_correction;
			if (oep_param.save_iter_corrections()>0 and (iter%oep_param.save_iter_corrections()==0))
				save(mrks_correction,"mrks_correction"+std::to_string(iter));
		}

		Fnemo = (Jnemo + Unemo + Voep*KS_nemo);

		tensorT Fock_no_ocep = matrix_inner(world, R2KS_nemo, Fnemo, false);
		Kinetic<double,3> T(world);
		Fock_no_ocep += T(R2KS_nemo, KS_nemo);

		// recompute the OCEP correction with the updated Fock matrix
		if (need_ocep_correction(model)) {
			real_function_3d ocep_correction;
			for (int i=0; i<5; ++i) {
				ocep_correction= compute_ocep_correction(ocep_numerator_HF, HF_nemo,KS_nemo,HF_Fock,KS_Fock);
				Fock_ocep=matrix_inner(world,R2KS_nemo,ocep_correction*KS_nemo);
				tensorT KS_Fock_old=copy(KS_Fock);
				KS_Fock=Fock_no_ocep+Fock_ocep;
				if ((KS_Fock_old-KS_Fock).normf()<get_calc_param().econv()) break;
				if (print_debug) {
					print("fock with updated ocep Fock matrix, ocep micro",i);
					print(KS_Fock);
				}
			}
			Voep+=ocep_correction;
			Fnemo+=(ocep_correction*KS_nemo);
		} else {
			KS_Fock=copy(Fock_no_ocep);
		}

		timer1.tag("compute oep");

		if (not get_calc_param().do_localize()) {

			// report the off-diagonal Fock matrix elements because canonical orbitals are used
			tensorT F_offdiag = copy(KS_Fock);
			for (int i = 0; i < KS_Fock.dim(0); ++i) F_offdiag(i, i) = 0.0;
			double max_F_offidag = F_offdiag.absmax();
			if (print_debug) print("F max off-diagonal ", max_F_offidag);

			// diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
			// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
			tensorT overlap = matrix_inner(world, R*KS_nemo, R*KS_nemo, true);
			tensorT KSeig;
			tensorT X = calc->get_fock_transformation(world, overlap, KS_Fock, KSeig, calc->aocc,
					FunctionDefaults<3>::get_thresh());
			KS_nemo = truncate(transform(world, KS_nemo, X));
			normalize(KS_nemo,R);
			// rotate_subspace(world, X, solver, 0, KS_nemo.size());
			Fnemo = transform(world, Fnemo, X);

			timer1.tag("canonicalization");
		}
		Fnemo=truncate(Fnemo);

        projector_irrep p=projector_irrep(calc->param.pointgroup())
                .set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(false);
        Voep=p(Voep)[0];

        // compute new (current) energy
        //double old_energy = energy;
        double Ex_KS = compute_exchange_energy_vir(KS_nemo, Voep);

		std::vector<double> oldenergies=energies;
        energies=compute_energy(KS_nemo, Ex_KS);
        energy =energies[0];

		timer1.tag("compute energy");

		BSHApply<double,3> bsh_apply(world);
		bsh_apply.metric=R_square;
		bsh_apply.levelshift=oep_param.levelshift();
		bsh_apply.lo=get_calc()->param.lo();
		auto [residual,eps_update] =bsh_apply(KS_nemo,KS_Fock,Fnemo);
		timer1.tag("apply BSH");

		double bshnorm=norm2(world,residual)/sqrt(KS_nemo.size());

		// KAIN solver (helps to converge)
		vecfuncT nemo_new = (solver.update(KS_nemo, residual, oep_param.kain_param()[0], oep_param.kain_param()[1]));
		truncate(world, nemo_new);
		normalize(nemo_new,R);

		// estimate the orbital energies, as they enter the potential
		if (get_calc_param().print_level()>=10)  print(KS_Fock);

		// What is step restriction?
		calc->do_step_restriction(world, KS_nemo, nemo_new, "ab spin case");
		orthonormalize(nemo_new,R);
		KS_nemo = nemo_new;
		timer1.tag("post-process");
		if (oep_param.saving_amount() >= 3)
			for (size_t n=0; n<KS_nemo.size(); ++n) save(KS_nemo[n], "KS_nemo"+stringify(n)+"iter"+stringify(iter));

		double deltadensity=0.0;
		converged=check_convergence(energies,oldenergies,bshnorm,deltadensity,get_calc_param(),
				get_calc_param().econv(),get_calc_param().dconv());
		if (world.rank() == 0) {
			printf("\nfinished %s iteration %2lu at time %8.1fs with energy %12.8f; residual %12.8f\n\n",
					model.c_str(), iter, wall_time(), energy,bshnorm);
		}

		if (converged) break;
	}

	if (world.rank()==0) {
		if (converged) print("\n     +++ Iterations converged +++\n");
		else print("\n     --- Iterations failed ---\n\n");
	}
	auto results=compute_and_print_final_energies(model,Voep,KS_nemo,KS_Fock,HF_nemo,HF_Fock);
	if (not converged) results.Econv=0.0;
	return results;
}


/// The following function tests all essential parts of the OEP program qualitatively and some also quantitatively
bool OEP::selftest() {

    printf("\n   +++ starting test of the OEP program +++\n\n");

	// hack away constness for the refence here
    const_cast<Nemo*>(reference.get())->value();
    calc->copy_data(world,*(reference->get_calc()));

    print("HF Fock operator ", reference->make_fock_operator()->info());
    print("OEP Fock operator", make_fock_operator()->info());

    const vecfuncT& HF_nemo1=reference->get_calc()->get_amo();
	int ierr=0;
	set_protocol(get_calc_param().econv());

    // Start by testing all important functions. If something severe fails, they throw an exception
    print("\n   >> test calculation of all important functions - severe errors will cause an exception\n");

    test_output hfdens("test construction of HF density");
    const real_function_3d rho_HF = compute_density(HF_nemo1);
    hfdens.end(true);

    test_output hf_recompute("test recomputation of HF orbitals and fock matrix");
	auto [HF_Fock, HF_nemo] = recompute_HF(HF_nemo1);
	double err_hf=norm2(world,HF_nemo1-HF_nemo);
	hf_recompute.end(err_hf<get_calc_param().econv());

    test_output slater("test computation of the Slater potential");
    const real_function_3d Vs = compute_slater_potential(HF_nemo);
    double refn1=1.70581413e+01;
    double n1=Vs.norm2();
    slater.logger << "  norm of the Slater potential " << n1 << std::endl;
    ierr+=slater.end(fabs(n1-refn1) < get_calc_param().econv());

    test_output ihf("test computation of the energy_weighted density IHF");
    const real_function_3d IHF = compute_energy_weighted_density_local(HF_nemo,HF_Fock);
    double refn2=4.48800120e+00;
    double n2=IHF.norm2();
    ihf.logger << "  norm of the IHF potential " << n2 << std::endl;
    ierr+=ihf.end(fabs(n2-refn2) < get_calc_param().econv());

    test_output kinhf("test computation of the kinetic density for dcep");
    const real_function_3d kin_tot_HF = compute_total_kinetic_density(HF_nemo);
    double refn3=7.02281570e+00;
    double n3=kin_tot_HF.norm2();
    kinhf.logger << "  norm of the kin_tot_HF potential " << n3 << std::endl;
    ierr+=kinhf.end(fabs(n3-refn3) < get_calc_param().econv());

    test_output paulihf("test computation of the pauli kinetic density for mrks");
    const real_function_3d kin_P_HF = compute_Pauli_kinetic_density(HF_nemo);
    double refn4= 7.59945929e-02;
    double n4=kin_P_HF.norm2();
    paulihf.logger << "  norm of the kin_P_HF potential " << n4 << std::endl;
    ierr+=paulihf.end(fabs(n4-refn4) < get_calc_param().econv());

    print("\n   >> test some quantities based on the reference HF calculation\n");

    vecfuncT Knemo;
	compute_exchange_potential(HF_nemo, Knemo);

    test_output conv_hf("test conventional HF exchange energy");
    const double Exconv_HF_correct = -2.66691504; // exchange energy from nemo calculation
    conv_hf.logger << "HF exchange energy of the system should be " <<  Exconv_HF_correct << " Eh" << std::endl;
    const double Exconv_HF = compute_exchange_energy_conv(R*HF_nemo, R*Knemo);
    const double Exconv_HF_diff = fabs(Exconv_HF_correct - Exconv_HF);
    conv_hf.logger << "the HF exchange energy of the system is ... " <<  Exconv_HF <<" Eh" << std::endl;
    conv_hf.logger << "the error is ... " <<  std::scientific << Exconv_HF_diff <<" Eh" << std::endl;
    ierr+=conv_hf.end(Exconv_HF_diff < get_calc_param().econv());


    test_output slater_x("test HF exchange energy via Slater potential");
    const double ExVs_HF_correct = -2.666912; // exchange energy from nemo calculation
    slater_x.logger << "  HF exchange energy of the system should be "
    		<< std::scientific << std::setprecision(8) <<  ExVs_HF_correct << " Eh"<<std::endl;
    const double ExVs_HF = 0.5*inner(rho_HF, Vs);
    const double ExVs_HF_diff = fabs(ExVs_HF_correct - ExVs_HF);
    slater_x.logger <<"  the HF exchange energy of the system is ... " << ExVs_HF << " Eh" << std::endl;
    slater_x.logger <<"  error: "  << ExVs_HF_diff <<  " Eh" << std::endl;;
    ierr+=slater_x.end(ExVs_HF_diff < get_calc_param().econv());


    test_output slater_x_vir("test virial HF exchange energy (with Slater potential)");
    const double Exvir_HF_correct = -3.00661935; // exchange energy from a test calculation with HF reference
    slater_x_vir.logger << "  HF virial exchange energy of the system should be "
    		<< std::scientific << std::setprecision(8) << Exvir_HF_correct << " Eh" << std::endl;
    const double Exvir_HF = compute_exchange_energy_vir(HF_nemo, Vs);
    const double Exvir_HF_diff = fabs(Exvir_HF_correct - Exvir_HF);
    slater_x_vir.logger << "  the virial HF exchange energy of the system is ... " << Exvir_HF <<" Eh" << std::endl;
    slater_x_vir.logger << "  error: " << Exvir_HF_diff << " Eh" << std::endl;
    ierr+=slater_x_vir.end(Exvir_HF_diff < get_calc_param().econv());


    print("\n   >> test solve_oep function with oaep model for 2 iterations\n");

    tensorT KS_Fock=copy(HF_Fock);
	vecfuncT KS_nemo = copy(world,calc->amo);
	real_function_3d Voep = copy(Vs);
	//double energy=iterate("oaep",HF_nemo,HF_Fock,KS_nemo,KS_Fock,Voep,Vs);


    test_output ihf_vir("test virial ocep exchange energy");
    save(IHF,"IHF");
    vecfuncT empty; tensorT fock0(2,2);
    real_function_3d V_HFocep = compute_ocep_correction(IHF, HF_nemo, KS_nemo, HF_Fock, KS_Fock);
    const double Exvir_HFocep_correct = -5.01251194e+00; // exchange energy from a test calculation with HF reference
    ihf_vir.logger  << "HF OCEP virial exchange energy of the system should be "
    		<< std::scientific << std::setprecision(8) << Exvir_HFocep_correct << " Eh" << std::endl;
    const double Exvir_HFocep = compute_exchange_energy_vir(HF_nemo, V_HFocep);
    const double Exvir_HFocep_diff = fabs(Exvir_HFocep_correct - Exvir_HFocep);
    ihf_vir.logger << "  the virial HF OCEP exchange energy of the system is ... "<< Exvir_HFocep<<" Eh"<<std::endl;
    ihf_vir.logger << " error: "<< Exvir_HFocep_diff << " Eh" << std::endl;
    ierr+=ihf_vir.end(Exvir_HFocep_diff < get_calc_param().econv());


    test_output ihf_kin_vir("test virial dcep exchange energy");
    real_function_3d V_HFdcep = compute_dcep_correction(kin_tot_HF, HF_nemo, KS_nemo);
    const double Exvir_HFdcep_correct = 7.83817133e-02; // exchange energy from a test calculation with HF reference
    ihf_kin_vir.logger << "HF DCEP virial exchange energy of the system should be "
    		<< std::scientific << std::setprecision(8) << Exvir_HFdcep_correct << " Eh" << std::endl;
    const double Exvir_HFdcep = compute_exchange_energy_vir(HF_nemo, V_HFdcep);
    const double Exvir_HFdcep_diff = fabs(Exvir_HFdcep_correct - Exvir_HFdcep);
    ihf_kin_vir.logger << "  the virial HF DCEP exchange energy of the system is ... " << Exvir_HFdcep << " Eh" << std::endl;
    ihf_kin_vir.logger << "  error: " << Exvir_HFdcep_diff << " Eh" << std::endl;
    ierr+=ihf_kin_vir.end(Exvir_HFdcep_diff < get_calc_param().econv());


    test_output ihf_pauli_vir("test virial mrks exchange energy ");
    real_function_3d V_HFmrks = compute_mrks_correction(kin_P_HF,HF_nemo, KS_nemo);
    const double Exvir_HFmrks_correct = 2.97175117e-01; // exchange energy from a test calculation with HF reference
    ihf_pauli_vir.logger << "HF mRKS virial exchange energy of the system should be"
    		<< std::scientific << std::setprecision(8) << Exvir_HFmrks_correct << " Eh" << std::endl;
    const double Exvir_HFmrks = compute_exchange_energy_vir(HF_nemo, V_HFmrks);
    const double Exvir_HFmrks_diff = fabs(Exvir_HFmrks_correct - Exvir_HFmrks);
    ihf_pauli_vir.logger << "  the virial HF mRKS exchange energy of the system is ... " << Exvir_HFmrks << " Eh" << std::endl;;
    ihf_pauli_vir.logger << "  error:" << Exvir_HFmrks_diff << " Eh" << std::endl;
    ierr+=ihf_pauli_vir.end(Exvir_HFmrks_diff < get_calc_param().econv());


    test_output hf_kinetic("test HF kinetic energy");
    const double Ekin_HF_correct = 14.57304144; // HF kinetic energy from a test calculation with HF reference (OEP: maxiter = 2)
    hf_kinetic.logger << "HF kinetic energy of the system should be"
    		<< std::scientific << std::setprecision(8) << Ekin_HF_correct << " Eh" << std::endl;
    const double Ekin_HF = compute_kinetic_energy(HF_nemo);
    const double Ekin_HF_diff = fabs(Ekin_HF_correct - Ekin_HF);
    hf_kinetic.logger << "  the HF kinetic energy of the system is ... " << Ekin_HF << " Eh" << std::endl;
    hf_kinetic.logger << "  error:" << Ekin_HF_diff << " Eh" << std::endl;
    ierr+=hf_kinetic.end(Ekin_HF_diff < get_calc_param().econv());





    print("\n   >> test solve_oep function with mRKS model for 2 iterations\n");
    solve(HF_nemo);
    print("\n   >> solve_oep test finished, calculating test quantities based on the new KS orbitals and eigenvalues\n");

    KS_nemo = calc->amo;
    //tensorT& KS_eigvals = calc->aeps;

    const real_function_3d rho_KS = compute_density(KS_nemo);

    print("loading final potential ...");
    load(Voep, "OEPapprox_final");
    print(   "... done\n");

    vecfuncT Jnemo;
    compute_coulomb_potential(KS_nemo, Jnemo);
    compute_exchange_potential(KS_nemo, Knemo);

    test_output conv_ks("test conventional KS exchange energy");
    const double Ex_conv_correct = -2.68048325e+00;
    conv_ks.logger << "KS conventional exchange energy of the system should be "
    		<< std::scientific << std::setprecision(8) << Ex_conv_correct << " Eh" << std::endl;
    const double Ex_conv = compute_exchange_energy_conv(R*KS_nemo, R*Knemo);
    const double Ex_conv_diff = fabs(Ex_conv_correct - Ex_conv);
    conv_ks.logger << "  the conventional KS exchange energy of the system is ... " << Ex_conv << " Eh" << std::endl;;
    conv_ks.logger << "  error: "<< Ex_conv_diff << " Eh" << std::endl;
    ierr+=conv_ks.end(Ex_conv_diff < get_calc_param().econv());


    test_output vir_ks("test virial KS exchange energy");
    const double Ex_vir_correct = -2.67728104e+00; // exchange energy from a test calculation with HF reference (OEP: maxiter = 2)
    vir_ks.logger << "KS virial exchange energy of the system should be"
    		<< std::scientific << std::setprecision(8)  << Ex_vir_correct << " Eh" << std::endl;
    const double Ex_vir = compute_exchange_energy_vir(KS_nemo, Voep);
    const double Ex_vir_diff = fabs(Ex_vir_correct - Ex_vir);
    vir_ks.logger << "  the virial KS exchange energy of the system is ... " << Ex_vir << " Eh" << std::endl;
    vir_ks.logger << "  error: " << Ex_vir_diff << " Eh" << std::endl;
    ierr+=vir_ks.end(Ex_vir_diff < get_calc_param().econv());


    test_output total_energy("test final total energy");
    const double Etot_correct = -1.45708107e+01; // total energy (conv) from a test calculation with HF reference (OEP: maxiter = 2)
    total_energy.logger << "final total energy of the system should be "
    		<< Etot_correct << " Eh" << std::endl;
    const double Etot = compute_energy(KS_nemo, Ex_conv)[0];
    const double Etot_diff = fabs(Etot_correct - Etot);
    total_energy.logger << "  the final total energy of the system is ... " << Etot << " Eh" << std::endl;
    total_energy.logger << "  error: " << Etot_diff << "Eh" << std::endl;
    ierr+=total_energy.end(Etot_diff < get_calc_param().econv());


    test_output ks_kinetic_energy("test KS kinetic energy");
    const double Ekin_KS_correct = 1.46848100e+01; // KS kinetic energy from a test calculation with HF reference (OEP: maxiter = 2)
    ks_kinetic_energy.logger << "KS kinetic energy of the system should be " << Ekin_KS_correct << " Eh" << std::endl;
    const double Ekin_KS = compute_kinetic_energy(KS_nemo);
    const double Ekin_KS_diff = fabs(Ekin_KS_correct - Ekin_KS);
    ks_kinetic_energy.logger << "  the KS kinetic energy of the system is ... " << Ekin_KS << " Eh" << std::endl;
    ks_kinetic_energy.logger << "  error: " << Ekin_KS_diff << " Eh" << std::endl;
    ierr+=ks_kinetic_energy.end(Ekin_KS_diff < get_calc_param().econv());


    const double Tc = Ekin_HF - Ekin_KS;

    test_output corr_energy("test quantity Delta Evir (14) after 2 iterations");
    const double DEvir_14_correct = 3.20220918e+00; // DEvir_14 (in mEh) from a test calculation with HF reference (OEP: maxiter = 2)
    corr_energy.logger << "Delta Evir (14) of the system should be " << DEvir_14_correct << " mEh" << std::endl;
    const double DEvir_14 = (Ex_vir - Ex_conv)*1000.0;
    const double DEvir_14_diff = fabs(DEvir_14_correct - DEvir_14);
    corr_energy.logger << "  Delta Evir (14) of the system is ... " << DEvir_14 << " mEh" << std::endl;
    corr_energy.logger << "  error: "<< DEvir_14_diff << " mEh" << std::endl;
    ierr+=corr_energy.end(DEvir_14_diff*0.001 < get_calc_param().econv());


    test_output devir17("test quantity Delta Evir (17) after 2 iterations");
    const double DEvir_17_correct = 2.13166331e+02; // DEvir_17 (in mEh) from a test calculation with HF reference (OEP: maxiter = 2)
    devir17.logger << "Delta Evir (17) of the system should be " << DEvir_17_correct <<  " mEh" << std::endl;
    const double DEvir_17 = (Ex_vir - ExVs_HF - 2.0*Tc)*1000.0;
    const double DEvir_17_diff = fabs(DEvir_17_correct - DEvir_17);
    devir17.logger << "  Delta Evir (17) of the system is ... " << DEvir_17<< " mEh" << std::endl;
    devir17.logger << "  error: " << DEvir_17_diff << " mEh" << std::endl;
    ierr+=devir17.end(DEvir_17_diff*0.001 < get_calc_param().econv());


    test_output delta_rho("test quantity Delta rho after 2 iterations");
    const double Drho_correct = 3.66404755e-02; // Drho from a test calculation with HF reference (OEP: maxiter = 2)
    delta_rho.logger << "Delta rho of the system should be "<< Drho_correct << std::endl;
    const double Drho = compute_delta_rho(rho_HF, rho_KS);
    const double Drho_diff = fabs(Drho_correct - Drho);
    delta_rho.logger << "  Delta rho of the system is ... " << Drho <<std::endl;
    delta_rho.logger << "  error: " <<  Drho_diff << std::endl;
    ierr+=delta_rho.end(Drho_diff < get_calc_param().econv());


    // TODO: What else can be checked?

    print("+++ OEP test finished +++\n");

    if (ierr==0) print("\n  All calculated results are correct, everything ok!\n");
    else print("  ATTENTION! There are errors in the results, see above!\n");
    return ierr;

}

} /* namespace madness */
