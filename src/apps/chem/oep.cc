/*
 * oep.cpp
 *
 *  Created on: Nov 6, 2019
 *      Author: fbischoff
 */

#include <chem/oep.h>
#include <chem/BSHApply.h>



namespace madness {

/// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
/// for other functionals, slater potential must be modified
/// HF orbitals and eigenvalues are used as the guess here
/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps
/// converged if norm, total energy difference and orbital energy differences (if not OAEP) are converged
void OEP::solve(const vecfuncT& HF_nemo1, const tensorT& HF_eigvals) {

	// recompute HF Fock matrix and orbitals
	auto [HF_Fock, HF_nemo] = recompute_HF(HF_nemo1);

	// compute Slater potential Vs and average IHF from HF orbitals and eigenvalues
	const real_function_3d Vs = compute_slater_potential(HF_nemo);
	if (oep_param.saving_amount() >= 1) save(Vs, "Slaterpotential");

	// set KS_nemo as reference to MOs
	vecfuncT& KS_nemo = calc->amo;
	tensorT& KS_eigvals = calc->aeps; // 1d tensor of same length as KS_nemo

	real_function_3d Voep = copy(Vs);

	double energy=iterate(oep_param.model(),HF_nemo,HF_Fock,KS_nemo,KS_eigvals,Voep,Vs);

	save(Voep,"OEP_final");
//	if (calc->param.save()) calc->save_mos(world);

	printf("      +++ FINAL TOTAL %s ENERGY = %15.8f  Eh +++\n\n\n", oep_param.model().c_str(), energy);

}

std::tuple<Tensor<double>, vecfuncT> OEP::recompute_HF(const vecfuncT& HF_nemo) const {

	timer timer1(world);
    const vecfuncT R2nemo=truncate(R_square*HF_nemo);
	vecfuncT psi, Jnemo, Knemo, pcmnemo, Unemo;
	Nemo::compute_nemo_potentials(HF_nemo, psi, Jnemo, Knemo, pcmnemo, Unemo);
	vecfuncT Vnemo=Unemo+Jnemo-Knemo;
	if (do_pcm()) Vnemo+=pcmnemo;
	tensorT HF_Fock=matrix_inner(world,R2nemo,Vnemo,false);   // not symmetric actually
	Kinetic<double,3> T(world);
	HF_Fock+=T(R2nemo,HF_nemo);


//	BSHApply<double,3> bsh_apply(world);
//	bsh_apply.metric=R_square;
//	bsh_apply.lo=get_calc()->param.lo();
//	bsh_apply.do_coupling=calc->param.do_localize();
//	auto [residual,eps_update] =bsh_apply(HF_nemo,HF_Fock,Vnemo);
//
//	vecfuncT nemo_new=HF_nemo-residual;
//	if (param.print_level()>=10) {
//		print("HF Fock matrix");
//		print(HF_Fock);
//		double rnorm=norm2(world,residual);
//		print("residual for HF nemo recomputation",rnorm);
//	}
	vecfuncT nemo_new=HF_nemo;
	timer1.end("recompute HF");
	return std::make_tuple(HF_Fock, nemo_new);
}

double OEP::compute_and_print_final_energies(const std::string model, const real_function_3d& Voep,
		const vecfuncT& KS_nemo, const tensorT& KS_eigvals,
		const vecfuncT& HF_nemo, const tensorT& HF_eigvals) const {

	// print final orbital energies
	print("final shifted", model, "orbital energies:");
	print_orbens(KS_eigvals, homo_diff(HF_eigvals, KS_eigvals));
	print("HF/KS HOMO energy difference of", homo_diff(HF_eigvals, KS_eigvals), "Eh is already included\n");

	// final Jnemo and Knemo have to be computed again in order to calculate final energy
	vecfuncT Jnemo, Knemo, Knemo_HF;
	compute_coulomb_potential(KS_nemo, Jnemo);
	compute_exchange_potential(KS_nemo, Knemo);

	Exchange<double,3> K(world);
	K.set_parameters(R_square*HF_nemo,HF_nemo,calc->aocc);
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
	double Evir = compute_energy(KS_nemo, Ex_vir)[0];

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

	if (not calc->param.do_localize()) {
		print("---------------------------------------------------------------------------");
		double E_0 = compute_E_zeroth(KS_eigvals);
		double E_1 = compute_E_first(R*KS_nemo, R*Jnemo, R*Knemo, Voep);

		printf("  E^(0)               = %15.8f  Eh", E_0);
		printf("\n  E^(1)               = %15.8f  Eh", E_1);
		printf("\n  E^(0) + E^(1)       = %15.8f  Eh", E_0 + E_1);
		printf("\n  difference to Econv = %15.8f mEh\n\n", (E_0 + E_1 - Econv)*1000.0);
	}
	Tensor<double> f_pp = compute_fock_diagonal_elements(calc->aeps, KS_nemo, Knemo, Voep);

	print("KS Fock matrix elements using Voep/virial ", calc->aeps);
	print("KS Fock matrix elements using K operator  ", f_pp);

	return Econv;
}

double OEP::iterate(const std::string model, const vecfuncT& HF_nemo, const tensorT& HF_Fock,
		vecfuncT& KS_nemo, tensorT& KS_eigvals, real_function_3d& Voep, const real_function_3d Vs) const {

	bool print_debug=(calc->param.print_level()>=10) and (world.rank()==0);

	// get the HF orbital energies from its Fock matrix
	auto [HF_eigvals, evec] = syev(HF_Fock);
	KS_eigvals=copy(HF_eigvals);

	// compute the constant HF contributions to the OEP hierarchy
	const real_function_3d ocep_numerator_HF=-1.0*compute_energy_weighted_density_local(HF_nemo,HF_Fock);
	const real_function_3d dcep_numerator_HF=compute_total_kinetic_density(HF_nemo);
	const real_function_3d mrks_numerator_HF=compute_Pauli_kinetic_density(HF_nemo);

	typedef allocator<double, 3> allocT;
	typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
	allocT alloc(world, KS_nemo.size());
	solverT solver(allocT(world, KS_nemo.size()),calc->param.print_level()>4);
	solver.set_maxsub(calc->param.maxsub());

	std::deque<tensorT> eps_history;
	double energy=0.0;
	std::vector<double> energies(1,0.0);
	bool converged=false;
	tensorT old_eigvals = copy(KS_eigvals);

	tensorT F=copy(HF_Fock);

	timer timer1(world,calc->param.print_level()>=3);
	for (int iter = 0; iter < oep_param.maxiter(); ++iter) {
//		print("\n     ***", model, "iteration", iter, "***\n");

	    if (calc->param.do_localize()) KS_nemo=localize(KS_nemo,calc->param.econv(),iter==0);

	    // compute parts of the KS Fock matrix J, Unuc and Voep
		vecfuncT Jnemo, Unemo, Fnemo;
		compute_nemo_potentials(KS_nemo, Jnemo, Unemo);
		vecfuncT R2KS_nemo = truncate(R_square*KS_nemo);

		timer1.tag("compute potentials");

		Voep=copy(Vs);

		// treat ocep correction separately as it depends on the KS Fock matrix
		tensorT Fock_ocep;
		if (need_ocep_correction(model)) {
			real_function_3d ocep_correction = compute_ocep_correction(ocep_numerator_HF, HF_nemo,KS_nemo,HF_Fock,F);
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
		F=Fock_no_ocep+Fock_ocep;
		if (print_debug) {
			print("fock");
			print(F);
		}

		// recompute the OCEP correction with the updated Fock matrix
		if (need_ocep_correction(model)) {
			real_function_3d ocep_correction = compute_ocep_correction(ocep_numerator_HF, HF_nemo,KS_nemo,HF_Fock,F);
			Voep+=ocep_correction;
			Fnemo+=(ocep_correction*KS_nemo);

			Fock_ocep=matrix_inner(world,R2KS_nemo,ocep_correction*KS_nemo);
			F=Fock_no_ocep+Fock_ocep;
		}

		if (print_debug) {
			print("fock with updated ocep Fock matrix");
			print(Fock_no_ocep+Fock_ocep);
		}
		timer1.tag("compute oep");

		if (not calc->param.do_localize()) {

			// report the off-diagonal Fock matrix elements because canonical orbitals are used
			tensorT F_offdiag = copy(F);
			for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
			double max_F_offidag = F_offdiag.absmax();
			if (print_debug) print("F max off-diagonal ", max_F_offidag);

			// diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
			// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
			tensorT overlap = matrix_inner(world, R*KS_nemo, R*KS_nemo, true);
			tensorT KSeig;
			tensorT X = calc->get_fock_transformation(world, overlap, F, KSeig, calc->aocc,
					FunctionDefaults<3>::get_thresh());
			KS_nemo = truncate(transform(world, KS_nemo, X));
			normalize(KS_nemo,R);
			rotate_subspace(world, X, solver, 0, KS_nemo.size());
			Fnemo = transform(world, Fnemo, X);

			double delta_eig=(KSeig-KS_eigvals).normf();
			KS_eigvals=KSeig;

			if (print_debug) {
				print("delta eigenvalues",delta_eig);
			}
			timer1.tag("canonicalization");
		}
		Fnemo=truncate(Fnemo);

		// compute new (current) energy
        double old_energy = energy;
        double Ex_KS = compute_exchange_energy_vir(KS_nemo, Voep);

		std::vector<double> oldenergies=energies;
        energies=compute_energy(KS_nemo, Ex_KS);
        energy =energies[0];

		// print orbital energies:
		if (print_debug) {
			print("current orbital energies");
			print_orbens(KS_eigvals);
			print("HF/KS HOMO energy difference of", homo_diff(HF_eigvals, KS_eigvals), "Eh is not yet included");
		}

		timer1.tag("compute energy");

		BSHApply<double,3> bsh_apply(world);
		bsh_apply.metric=R_square;
		bsh_apply.lo=get_calc()->param.lo();
		bsh_apply.do_coupling=calc->param.do_localize();
		auto [residual,eps_update] =bsh_apply(KS_nemo,F,Fnemo);
		timer1.tag("apply BSH");

		double bshnorm=norm2(world,residual)/sqrt(KS_nemo.size());

		// KAIN solver (helps to converge)
		vecfuncT nemo_new = (solver.update(KS_nemo, residual, oep_param.kain_param()[0], oep_param.kain_param()[1]));
		truncate(world, nemo_new);
		normalize(nemo_new,R);

		// estimate the orbital energies, as they enter the potential
		eps_history.push_back(copy(F));
		if (eps_history.size()>solver.get_c().size()) eps_history.pop_front();
		tensorT fock1(F.dim(0),F.dim(1));
		int i=0;
		for (auto eps : eps_history) {
			if (calc->param.print_level()>10) print("c[i], eps",solver.get_c()[i],eps);
			fock1 += eps*solver.get_c()[i++];
		}
		F=copy(fock1);
		if (calc->param.print_level()>=10) print("KS_eigvals projection",F);

		// What is step restriction?
		calc->do_step_restriction(world, KS_nemo, nemo_new, "ab spin case");
		orthonormalize(nemo_new,R);
		KS_nemo = nemo_new;
		timer1.tag("post-process");

		double deltadensity=0.0;
		converged=check_convergence(energies,oldenergies,bshnorm,deltadensity,calc->param,
				calc->param.econv(),calc->param.dconv());
		if (world.rank() == 0) {
			printf("\nfinished %s iteration %2d at time %8.1fs with energy %12.8f; residual %12.8f\n\n",
					model.c_str(), iter, wall_time(), energy,bshnorm);
		}

		// evaluate convergence via norm error and energy difference
		if (converged) {
			if (oep_param.is_oaep()) converged = true;  // if OAEP, the following evaluation is not necessary
			else {
				// build vector of convergence information of every orbital energy
    			std::vector<bool> conv(KS_eigvals.size());
    			for (long i = 0; i < KS_eigvals.size(); i++) {
    				if (fabs(KS_eigvals(i) - old_eigvals(i)) < calc->param.dconv()) conv[i] = true;
    				else conv[i] = false;
    			}
    			if (IsAlltrue(conv)) converged = true; // converged if all are converged
			}
		}

		if (converged) break;
	}

	if (world.rank()==0) {
		if (converged) print("\n     +++ Iterations converged +++\n");
		else print("\n     --- Iterations failed ---\n\n");
	}
	energy=compute_and_print_final_energies(model,Voep,KS_nemo,KS_eigvals,HF_nemo,HF_eigvals);
	if (not converged) energy=0.0;
	return energy;
}


/// The following function tests all essential parts of the OEP program qualitatively and some also quantitatively
void OEP::test_oep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

    bool everything_ok = true;

    // Start by testing all important functions. If something severe fails, they throw an exception
    print("\n   >> test calculation of all important functions - severe errors will cause an exception\n");

    print("test calculation of HOMO index from HF calculation");
    print("     the HOMO index is ...", homo_ind(HF_eigvals));
    print("  HOMO index computed successfully\n");

    print("test construction of HF density");
    const real_function_3d rho_HF = compute_density(HF_nemo);
    print("  HF density computed successfully\n");

    print("test construction of Slater potential");
    const real_function_3d Vs = compute_slater_potential(HF_nemo);
    print("  Slater potential computed successfully\n");

    print("test construction of IHF");
    const real_function_3d IHF = compute_energy_weighted_density(HF_nemo, HF_eigvals);
    print("  compute_energy_weighted_density computed successfully\n");

    print("test construction of kin_tot_HF (tau/rho HF)");
    const real_function_3d kin_tot_HF = compute_total_kinetic_density(HF_nemo);
    print("  kin_tot_HF computed successfully\n");

    print("test construction of kin_P_HF (tau_P/rho HF)");
    const real_function_3d kin_P_HF = compute_Pauli_kinetic_density(HF_nemo);
    print("  kin_P_HF computed successfully\n");

    print("\n   >> test some quantities based on the reference HF calculation\n");

    vecfuncT Knemo;
	compute_exchange_potential(HF_nemo, Knemo);

    print("test conventional HF exchange energy");
    const double Exconv_HF_correct = -2.66691504; // exchange energy from nemo calculation
    print("HF exchange energy of the system should be", Exconv_HF_correct, "Eh");
    const double Exconv_HF = compute_exchange_energy_conv(R*HF_nemo, R*Knemo);
    const double Exconv_HF_diff = fabs(Exconv_HF_correct - Exconv_HF);
    print("     the HF exchange energy of the system is ...", Exconv_HF, "Eh");
    print("     error:", Exconv_HF_diff, "Eh");
    if (Exconv_HF_diff <= param.econv()) print("  conventional HF exchange energy is correct\n");
    else {
    	print("  ATTENTION: conventional HF exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test HF exchange energy via Slater potential");
    const double ExVs_HF_correct = -2.66691504; // exchange energy from nemo calculation
    print("HF exchange energy of the system should be", ExVs_HF_correct, "Eh");
    const double ExVs_HF = 0.5*inner(rho_HF, Vs);
    const double ExVs_HF_diff = fabs(ExVs_HF_correct - ExVs_HF);
    print("     the HF exchange energy of the system is ...", ExVs_HF, "Eh");
    print("     error:", ExVs_HF_diff, "Eh");
    if (ExVs_HF_diff <= param.econv()) print("  HF exchange energy via Slater potential is correct\n");
    else {
    	print("  ATTENTION: HF exchange energy (via Slater potential) error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test virial HF exchange energy (with Slater potential)");
    const double Exvir_HF_correct = -3.00658754; // exchange energy from a test calculation with HF reference
    print("HF virial exchange energy of the system should be", Exvir_HF_correct, "Eh");
    const double Exvir_HF = compute_exchange_energy_vir(HF_nemo, Vs);
    const double Exvir_HF_diff = fabs(Exvir_HF_correct - Exvir_HF);
    print("     the virial HF exchange energy of the system is ...", Exvir_HF, "Eh");
    print("     error:", Exvir_HF_diff, "Eh");
    if (Exvir_HF_diff <= param.econv()) print("  virial HF exchange energy is correct\n");
    else {
    	print("  ATTENTION: virial HF exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test virial HF exchange energy (with only HF OCEP: IKS = 0)");
    real_function_3d V_HFocep = Vs + IHF;
    const double Exvir_HFocep_correct = -0.47639487; // exchange energy from a test calculation with HF reference
    print("HF OCEP virial exchange energy of the system should be", Exvir_HFocep_correct, "Eh");
    const double Exvir_HFocep = compute_exchange_energy_vir(HF_nemo, V_HFocep);
    const double Exvir_HFocep_diff = fabs(Exvir_HFocep_correct - Exvir_HFocep);
    print("     the virial HF OCEP exchange energy of the system is ...", Exvir_HFocep, "Eh");
    print("     error:", Exvir_HFocep_diff, "Eh");
    if (Exvir_HFocep_diff <= param.econv()) print("  virial HF OCEP exchange energy is correct\n");
    else {
    	print("  ATTENTION: virial HF OCEP exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test virial HF exchange energy (with only HF DCEP: IKS = 0, kin_tot_KS = 0)");
    real_function_3d V_HFdcep = Vs + IHF + kin_tot_HF;
    const double Exvir_HFdcep_correct = 4.03650400; // exchange energy from a test calculation with HF reference
    print("HF DCEP virial exchange energy of the system should be", Exvir_HFdcep_correct, "Eh");
    const double Exvir_HFdcep = compute_exchange_energy_vir(HF_nemo, V_HFdcep);
    const double Exvir_HFdcep_diff = fabs(Exvir_HFdcep_correct - Exvir_HFdcep);
    print("     the virial HF DCEP exchange energy of the system is ...", Exvir_HFdcep, "Eh");
    print("     error:", Exvir_HFdcep_diff, "Eh");
    if (Exvir_HFdcep_diff <= param.econv()) print("  virial HF DCEP exchange energy is correct\n");
    else {
    	print("  ATTENTION: virial HF DCEP exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test virial HF exchange energy (with only HF mRKS: IKS = 0, kin_P_KS = 0)");
    real_function_3d V_HFmrks = Vs + IHF + kin_P_HF;
    const double Exvir_HFmrks_correct = -0.84506060; // exchange energy from a test calculation with HF reference
    print("HF mRKS virial exchange energy of the system should be", Exvir_HFmrks_correct, "Eh");
    const double Exvir_HFmrks = compute_exchange_energy_vir(HF_nemo, V_HFmrks);
    const double Exvir_HFmrks_diff = fabs(Exvir_HFmrks_correct - Exvir_HFmrks);
    print("     the virial HF mRKS exchange energy of the system is ...", Exvir_HFmrks, "Eh");
    print("     error:", Exvir_HFmrks_diff, "Eh");
    if (Exvir_HFmrks_diff <= param.econv()) print("  virial HF mRKS exchange energy is correct\n");
    else {
    	print("  ATTENTION: virial HF mRKS exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test HF kinetic energy");
    const double Ekin_HF_correct = 14.57304144; // HF kinetic energy from a test calculation with HF reference (OEP: maxiter = 2)
    print("HF kinetic energy of the system should be", Ekin_HF_correct, "Eh");
    const double Ekin_HF = compute_kinetic_energy(HF_nemo);
    const double Ekin_HF_diff = fabs(Ekin_HF_correct - Ekin_HF);
    print("     the HF kinetic energy of the system is ...", Ekin_HF, "Eh");
    print("     error:", Ekin_HF_diff, "Eh");
    if (Ekin_HF_diff <= param.econv()) print("  HF kinetic energy is correct\n");
    else {
    	print("  ATTENTION: HF kinetic energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("\n   >> test solve_oep function with mRKS model for 2 iterations\n");
    solve(HF_nemo, HF_eigvals);
    print("\n   >> solve_oep test finished, calculating test quantities based on the new KS orbitals and eigenvalues\n");

    vecfuncT& KS_nemo = calc->amo;
    tensorT& KS_eigvals = calc->aeps;

    const real_function_3d rho_KS = compute_density(KS_nemo);

    real_function_3d Voep = Vs;
    print("loading final potential ...");
    load(Voep, "OEPapprox_final");
    print(   "... done\n");

    vecfuncT Jnemo;
    compute_coulomb_potential(KS_nemo, Jnemo);
    compute_exchange_potential(KS_nemo, Knemo);

    print("test conventional KS exchange energy");
    const double Ex_conv_correct = -2.68888478; // exchange energy from a test calculation with HF reference (OEP: maxiter = 2)
    print("KS conventional exchange energy of the system should be", Ex_conv_correct, "Eh");
    const double Ex_conv = compute_exchange_energy_conv(R*KS_nemo, R*Knemo);
    const double Ex_conv_diff = fabs(Ex_conv_correct - Ex_conv);
    print("     the conventional KS exchange energy of the system is ...", Ex_conv, "Eh");
    print("     error:", Ex_conv_diff, "Eh");
    if (Ex_conv_diff <= param.econv()) print("  conventional KS exchange energy is correct\n");
    else {
    	print("  ATTENTION: conventional KS exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test virial KS exchange energy");
    const double Ex_vir_correct = -2.81673416; // exchange energy from a test calculation with HF reference (OEP: maxiter = 2)
    print("KS virial exchange energy of the system should be", Ex_vir_correct, "Eh");
    const double Ex_vir = compute_exchange_energy_vir(KS_nemo, Voep);
    const double Ex_vir_diff = fabs(Ex_vir_correct - Ex_vir);
    print("     the virial KS exchange energy of the system is ...", Ex_vir, "Eh");
    print("     error:", Ex_vir_diff, "Eh");
    if (Ex_vir_diff <= param.econv()) print("  virial KS exchange energy is correct\n");
    else {
    	print("  ATTENTION: virial KS exchange energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test final total energy");
    const double Etot_correct = -14.56855740; // total energy (conv) from a test calculation with HF reference (OEP: maxiter = 2)
    print("final total energy of the system should be", Etot_correct, "Eh");
    const double Etot = compute_energy(R*KS_nemo, Ex_conv)[0];
    const double Etot_diff = fabs(Etot_correct - Etot);
    print("     the final total energy of the system is ...", Etot, "Eh");
    print("     error:", Etot_diff, "Eh");
    if (Etot_diff <= param.econv()) print("  final total energy is correct\n");
    else {
    	print("  ATTENTION: final total energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test KS kinetic energy");
    const double Ekin_KS_correct = 14.75193175; // KS kinetic energy from a test calculation with HF reference (OEP: maxiter = 2)
    print("KS kinetic energy of the system should be", Ekin_KS_correct, "Eh");
    const double Ekin_KS = compute_kinetic_energy(KS_nemo);
    const double Ekin_KS_diff = fabs(Ekin_KS_correct - Ekin_KS);
    print("     the KS kinetic energy of the system is ...", Ekin_KS, "Eh");
    print("     error:", Ekin_KS_diff, "Eh");
    if (Ekin_KS_diff <= param.econv()) print("  KS kinetic energy is correct\n");
    else {
    	print("  ATTENTION: KS kinetic energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    const double Tc = Ekin_HF - Ekin_KS;

    print("test quantity Delta Evir (14) after 2 iterations");
    const double DEvir_14_correct = -127.84938639; // DEvir_14 (in mEh) from a test calculation with HF reference (OEP: maxiter = 2)
    print("Delta Evir (14) of the system should be", DEvir_14_correct, "mEh");
    const double DEvir_14 = (Ex_vir - Ex_conv)*1000.0;
    const double DEvir_14_diff = fabs(DEvir_14_correct - DEvir_14);
    print("     Delta Evir (14) of the system is ...", DEvir_14, "mEh");
    print("     error:", DEvir_14_diff, "mEh");
    if (DEvir_14_diff*0.001 <= param.econv()) print("  quantity Delta Evir (14) is correct\n"); // mind the units
    else {
    	print("  ATTENTION: error of quantity Delta Evir (14) is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test quantity Delta Evir (17) after 2 iterations");
    const double DEvir_17_correct = 207.96150201; // DEvir_17 (in mEh) from a test calculation with HF reference (OEP: maxiter = 2)
    print("Delta Evir (17) of the system should be", DEvir_17_correct, "mEh");
    const double DEvir_17 = (Ex_vir - ExVs_HF - 2.0*Tc)*1000.0;
    const double DEvir_17_diff = fabs(DEvir_17_correct - DEvir_17);
    print("     Delta Evir (17) of the system is ...", DEvir_17, "mEh");
    print("     error:", DEvir_17_diff, "mEh");
    if (DEvir_17_diff*0.001 <= param.econv()) print("  quantity Delta Evir (17) is correct\n"); // mind the units
    else {
    	print("  ATTENTION: error of quantity Delta Evir (17) is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("test quantity Delta rho after 2 iterations");
    const double Drho_correct = 0.05266930; // Drho from a test calculation with HF reference (OEP: maxiter = 2)
    print("Delta rho of the system should be", Drho_correct);
    const double Drho = compute_delta_rho(rho_HF, rho_KS);
    const double Drho_diff = fabs(Drho_correct - Drho);
    print("     Delta rho of the system is ...", Drho);
    print("     error:", Drho_diff);
    if (Drho_diff <= param.econv()) print("  quantity Delta rho is correct\n"); // unitless
    else {
    	print("  ATTENTION: error of quantity Delta rho is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    // TODO: What else can be checked?

    print("+++ OEP test finished +++\n");

    if (everything_ok) print("\n  All calculated results are correct, everything ok!\n");
    else print("  ATTENTION! There are errors in the results, see above!\n");

}

} /* namespace madness */
