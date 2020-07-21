/*
 * oep.cpp
 *
 *  Created on: Nov 6, 2019
 *      Author: fbischoff
 */

#include <chem/oep.h>

namespace madness {

/// Iterative energy calculation for approximate OEP with EXACT EXCHANGE functional
/// for other functionals, slater potential must be modified
/// HF orbitals and eigenvalues are used as the guess here
/// note that KS_nemo is a reference and changes oep->get_calc()->amo orbitals
/// same for orbital energies (eigenvalues) KS_eigvals which is oep->get_calc()->aeps
/// converged if norm, total energy difference and orbital energy differences (if not OAEP) are converged
void OEP::solve_oep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

	double energy = 0.0;
	bool converged = false;
	unsigned int iter_counter = 0;

	// compute Slater potential Vs and average IHF from HF orbitals and eigenvalues
	const real_function_3d Vs = compute_slater_potential(HF_nemo, homo_ind(HF_eigvals));
	const real_function_3d kin_tot_HF = compute_total_kinetic_density(HF_nemo, HF_eigvals);
	const real_function_3d kin_P_HF = compute_Pauli_kinetic_density(HF_nemo, HF_eigvals);
	const real_function_3d rho_HF = compute_density(HF_nemo);
	if (oep_param.saving_amount() >= 1) save(Vs, "Slaterpotential");
	if (oep_param.saving_amount() >= 2) {
		save(rho_HF, "density_HF");
        if (oep_param.is_dcep()) save(kin_tot_HF, "kin_tot_HF");
        if (oep_param.is_mrks()) save(kin_P_HF, "kin_P_HF");
	}

	// compute ab initio HF exchange energy using equation (21) from Ospadov_2017 and HF kinetic energy
	// edit: this is redundant because it is the same as -0.5*<phi|K|phi>
	const double Ex_HF = 0.5*inner(rho_HF, Vs);
	const double Ekin_HF = compute_kinetic_energy(R*HF_nemo); // like T in Ospadov_2017, equation (20)

	// set KS_nemo as reference to MOs
	vecfuncT& KS_nemo = calc->amo;
	tensorT& KS_eigvals = calc->aeps; // 1d tensor of same length as KS_nemo
	if (oep_param.saving_amount() >= 3) save(compute_density(KS_nemo), "density_start");

	// if desired: save HF orbitals and orbital contributions to total density (orbital squares)
	if (oep_param.saving_amount() >= 3) {
		vecfuncT HF_nemo_square = square(world, HF_nemo);
    	for (size_t i = 0; i < HF_nemo.size(); i++) {
    		save(R*HF_nemo[i], "HF_orb_" + stringify(i));
    		save(2.0*R_square*HF_nemo_square[i], "HF_orb_square_" + stringify(i)); // 2 because closed shell
    	}
	}

	// all necessary operators applied on nemos
	vecfuncT Jnemo, Unemo, Vnemo, Knemo;
	real_function_3d Voep = Vs;

	// copy Vs to all old potentials for damping
	// attention: Voep_old is only used if damping is used, so if oep_param.damp_num() > 1
	std::vector<real_function_3d> Voep_old;
	if (oep_param.do_damping()) {
		for (unsigned int i = 1; i < oep_param.damp_num(); i++) {
			Voep_old.push_back(Vs);
		}
	}

	// define the solver
	typedef allocator<double, 3> allocT;
	typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
//    	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
	allocT alloc(world, KS_nemo.size());
	solverT solver(allocT(world, KS_nemo.size()));

	// iterate until self-consistency
	for (size_t iter = 0; iter < oep_param.maxiter(); ++iter) {
		iter_counter++;
		print("\n     ***", oep_param.model(), "iteration", iter_counter, "***\n");

		if (oep_param.is_ocep() or oep_param.is_dcep() or oep_param.is_mrks()) {

			// damping for better convergence of Voep
			if (oep_param.do_damping()) {
    			for (unsigned int i = 1; i < oep_param.damp_num() - 1; i++) {
    				Voep_old[i] = Voep_old[i - 1];
    			}
    			Voep_old[0] = Voep;
			}

    		// compute OCEP potential from current nemos and eigenvalues
			real_function_3d corr_ocep, corr_dcep, corr_mrks;
			corr_ocep = compute_oep_correction("ocep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
			if (oep_param.is_dcep()) corr_dcep = compute_oep_correction("dcep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
			if (oep_param.is_mrks()) corr_mrks = compute_oep_correction("mrks", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);

			// and shift potential so that HOMO_HF = HOMO_KS, so potential += (HOMO_HF - HOMO_KS)
			double shift = homo_diff(HF_eigvals, KS_eigvals);
			print("building new Veop: orbital shift is", shift, "Eh");

			// damping
			Voep = oep_param.damp_coeff()[0]*(Vs + corr_ocep + shift);
			if (oep_param.is_dcep()) Voep += oep_param.damp_coeff()[0]*corr_dcep;
			else if (oep_param.is_mrks()) Voep += oep_param.damp_coeff()[0]*corr_mrks;
			if (oep_param.do_damping()) {
				for (unsigned int i = 1; i < oep_param.damp_num(); i++) {
					Voep += oep_param.damp_coeff()[i]*Voep_old[i - 1];
				}
			}

			// save certain functions if desired
			if (oep_param.save_iter_orbs() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_orbs() == 0) {
			    	for (size_t i = 0; i < KS_nemo.size(); i++) {
			    		save(R*KS_nemo[i], "KS_orb_" + stringify(i) + "_iter_" + stringify(iter_counter));
			    	}
				}
			}
			if (oep_param.save_iter_density() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_density() == 0) {
					save(compute_density(KS_nemo), "density_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_kin_tot_KS() > 0 and oep_param.is_dcep()) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_kin_tot_KS() == 0) {
					save(compute_total_kinetic_density(KS_nemo, KS_eigvals), "kin_tot_KS_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_kin_P_KS() > 0 and oep_param.is_mrks()) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_kin_P_KS() == 0) {
					save(compute_Pauli_kinetic_density(KS_nemo, KS_eigvals), "kin_P_KS_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_ocep_correction() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_ocep_correction() == 0) {
					save(corr_ocep + shift, "OCEP_correction_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_dcep_correction() > 0 and oep_param.is_dcep()) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_dcep_correction() == 0) {
					save(corr_dcep + shift, "DCEP_correction_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_mrks_correction() > 0 and oep_param.is_mrks()) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_mrks_correction() == 0) {
					save(corr_mrks + shift, "mRKS_correction_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_total_correction() > 0) {
				if (oep_param.is_dcep() and (iter_counter == 2 or iter_counter % oep_param.save_iter_total_correction() == 0)) {
					save(corr_ocep + corr_dcep + shift, "total_correction_iter_" + stringify(iter_counter));
				}
				else if (oep_param.is_mrks() and (iter_counter == 2 or iter_counter % oep_param.save_iter_total_correction() == 0)) {
					save(corr_ocep + corr_mrks + shift, "total_correction_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_effective_potential() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_effective_potential() == 0) {
					save(Voep, "effective_potential_iter_" + stringify(iter_counter));
				}
			}

		}

		vecfuncT R2KS_nemo = R_square*KS_nemo;
		truncate(world, R2KS_nemo);

		// compute parts of the Fock matrix J, Unuc and Voep
		compute_nemo_potentials(KS_nemo, Jnemo, Unemo, Voep, Vnemo);

		// compute Fock matrix F = J + Voep + Vnuc and kinetic energy
		vecfuncT Fnemo = Jnemo + Vnemo + Unemo;
		truncate(world, Fnemo);
		tensorT F = matrix_inner(world, R2KS_nemo, Fnemo, false); // matrix_inner gives 2d tensor
		Kinetic<double,3> T(world);
		F += T(R2KS_nemo, KS_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

		// report the off-diagonal Fock matrix elements because canonical orbitals are used
        tensorT F_offdiag = copy(F);
        for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
        double max_F_offidag = F_offdiag.absmax();
        if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

		// compute new (current) energy
        double old_energy = energy;
        print("energy contributions of iteration", iter_counter);
        double Ex_KS = compute_exchange_energy_vir(R*KS_nemo, Voep);
		energy = compute_energy(R*KS_nemo, R*Jnemo, Ex_KS);
		// compute_exchange_potential(KS_nemo, Knemo);
		// double Ex_KS = compute_exchange_energy_conv(R*KS_nemo, R*Knemo);
		// there should be no difference between these two methods, because energy is only needed
		// for checking convergence threshold; but: Evir should be much faster because K is expensive

		// copy old orbital energies for convergence criterium at the end
		tensorT old_eigvals = copy(KS_eigvals);

        // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
        tensorT X; // must be formed from R*nemos but can then be used for nemos also
        tensorT overlap = matrix_inner(world, R*KS_nemo, R*KS_nemo, true);
        X = calc->get_fock_transformation(world, overlap, F, KS_eigvals, calc->aocc,
        		FunctionDefaults<3>::get_thresh());
        KS_nemo = transform(world, KS_nemo, X, trantol(), true);
        rotate_subspace(world, X, solver, 0, KS_nemo.size());

        truncate(world, KS_nemo);
        normalize(KS_nemo,R);

		// calculate new orbital energies (current eigenvalues from Fock-matrix)
		for (size_t i = 0; i < KS_nemo.size(); ++i) {
			KS_eigvals(i) = std::min(-0.05, F(i, i)); // orbital energy is set to -0.05 if it was above
		}

		/// TODO: Question: is this necessary in our programme or even bad?
		// if requested: subtract orbital shift from orbital energies
		if (calc->param.orbitalshift() > 0.0) {
			if (world.rank() == 0) print("shifting orbitals by ",
					calc->param.orbitalshift(), " to lower energies");
			KS_eigvals -= calc->param.orbitalshift();
		}

		// print orbital energies:
		print("orbital energies of iteration", iter_counter);
		print_orbens(KS_eigvals);
		print("HF/KS HOMO energy difference of", homo_diff(HF_eigvals, KS_eigvals), "Eh is not yet included");

		// construct the BSH operators ops
		std::vector<poperatorT> G = calc->make_bsh_operators(world, KS_eigvals);

		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
		Fnemo = transform(world, Fnemo, X, trantol(), true);
		truncate(world, Fnemo);

		// apply the BSH operators G (here ops) on the wave function
		scale(world, Fnemo, -2.0);
		vecfuncT GFnemo = apply(world, G, Fnemo);
		truncate(world, GFnemo);

		double n1 = norm2(world, KS_nemo);
		double n2 = norm2(world, GFnemo);
		print("\nnorm of nemo and GFnemo, ratio ", n1, n2, n1/n2);

		// compute the residuals for KAIN
		vecfuncT residual = KS_nemo - GFnemo;
		const double norm = norm2(world, residual) / sqrt(KS_nemo.size());

		// KAIN solver (helps to converge)
		vecfuncT nemo_new;
		if (norm < 5.0e-1) {
			nemo_new = (solver.update(KS_nemo, residual, oep_param.kain_param()[0], oep_param.kain_param()[1]));
		} else {
			nemo_new = GFnemo;
		}
		truncate(world, nemo_new);
		normalize(nemo_new,R);

		// What is step restriction?
		calc->do_step_restriction(world, KS_nemo, nemo_new, "ab spin case");
		orthonormalize(nemo_new,R);
		KS_nemo = nemo_new;

		// evaluate convergence via norm error and energy difference
		if ((norm < calc->param.dconv()) and (fabs(energy - old_energy) < oep_param.conv_thresh())) {

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

		if (calc->param.save()) calc->save_mos(world);

		if (world.rank() == 0) {
			printf("\nfinished iteration %2d at time %8.1fs with energy %12.8f\n", iter_counter, wall_time(), energy);
			print("current residual norm", norm, "\n");
		}

		if (converged) break;

	}

	if (converged) {
		if (world.rank() == 0) {
			print("\n     +++ Iterations converged +++\n");
			print(oep_param.model(), "converged after", iter_counter, "iterations\n\n");
		}
	}
	else {
		if (world.rank() == 0) print("\n     --- Iterations failed ---\n\n");
		energy = 0.0;
	}

	// calculate and print all final numbers
	print("\n  computing final orbitals, IKS and density");

	double shift_final = homo_diff(HF_eigvals, KS_eigvals);
	real_function_3d kin_tot_KS = compute_total_kinetic_density(KS_nemo, KS_eigvals);
	real_function_3d kin_P_KS = compute_Pauli_kinetic_density(KS_nemo, KS_eigvals);
	real_function_3d rho_KS = compute_density(KS_nemo);
	double Drho = compute_delta_rho(rho_HF, rho_KS);
	if (oep_param.saving_amount() >= 1) save(rho_KS, "density_final");
	if (oep_param.saving_amount() >= 2) {
        if (oep_param.is_dcep()) save(kin_tot_KS, "kin_tot_KS_final");
        if (oep_param.is_mrks()) save(kin_P_KS, "kin_P_KS_final");
    	for (size_t i = 0; i < KS_nemo.size(); i++) {
    		save(R*KS_nemo[i], "KS_orb_" + stringify(i) + "_final");
    	}
	}

	// if desired: print final KS orbital contributions to total density (nemo squares)
	if (oep_param.saving_amount() >= 3) {
    	vecfuncT KS_nemo_square = square(world, KS_nemo);
    	for (size_t i = 0; i < KS_nemo_square.size(); i++) {
    		save(2.0*R_square*KS_nemo_square[i], "KS_orb_square_" + stringify(i)); // 2 because closed shell
    	}
	}

	print("     done");

	if (oep_param.is_oaep()) {
		print("\n  computing final OAEP with converged OAEP orbitals and eigenvalues");
    	Voep = Vs + shift_final;
    	if (oep_param.saving_amount() >= 1) save(Voep, "OEPapprox_final");
	}
	if (oep_param.is_ocep()) {
		print("\n  computing final OCEP with converged OCEP orbitals and eigenvalues");
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
    	Voep = Vs + ocep_correction_final + shift_final;
    	if (oep_param.saving_amount() >= 1) {
    		save(ocep_correction_final + shift_final, "OCEP_correction_final");
    		save(Voep, "OEPapprox_final");
    	}
	}
	if (oep_param.is_dcep()) {
		print("\n  computing final DCEP with converged DCEP orbitals and eigenvalues");
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
    	real_function_3d dcep_correction_final = compute_oep_correction("dcep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
    	Voep = Vs + ocep_correction_final + dcep_correction_final + shift_final;
    	if (oep_param.saving_amount() >= 2) {
        	save(ocep_correction_final + shift_final, "OCEP_correction_final");
        	save(dcep_correction_final + shift_final, "DCEP_correction_final");
    	}
    	if (oep_param.saving_amount() >= 1) {
    		save(ocep_correction_final + dcep_correction_final + shift_final, "total_correction_final");
    		save(Voep, "OEPapprox_final");
    	}
	}
	if (oep_param.is_mrks()) {
		print("\n  computing final mRKS potential with converged mRKS orbitals and eigenvalues");
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
    	real_function_3d mrks_correction_final = compute_oep_correction("mrks", HF_nemo,HF_eigvals, KS_nemo, KS_eigvals);
    	Voep = Vs + ocep_correction_final + mrks_correction_final + shift_final;
    	if (oep_param.saving_amount() >= 2) {
        	save(ocep_correction_final + shift_final, "OCEP_correction_final");
        	save(mrks_correction_final + shift_final, "mRKS_correction_final");
    	}
    	if (oep_param.saving_amount() >= 1) {
    		save(ocep_correction_final + mrks_correction_final + shift_final, "total_correction_final");
    		save(Voep, "OEPapprox_final");
    	}
	}
	print("     done\n");

	// print final orbital energies
	print("final shifted", oep_param.model(), "orbital energies:");
	print_orbens(KS_eigvals, homo_diff(HF_eigvals, KS_eigvals));
	print("HF/KS HOMO energy difference of", homo_diff(HF_eigvals, KS_eigvals), "Eh is already included\n");

	// final Jnemo and Knemo have to be computed again in order to calculate final energy
	compute_coulomb_potential(KS_nemo, Jnemo);
	compute_exchange_potential(KS_nemo, Knemo);

	// compute final exchange energy using different methods and final kinetic energy
	double Ex_vir = compute_exchange_energy_vir(R*KS_nemo, Voep);
	double Ex_conv = compute_exchange_energy_conv(R*KS_nemo, R*Knemo);
	double Ekin_KS = compute_kinetic_energy(R*KS_nemo); // like Ts in Ospadov_2017, equation (22)
	double Tc = Ekin_HF - Ekin_KS; // like Tc = T - Ts in Ospadov_2017, equation (24)

	print("FINAL", oep_param.model(), "ENERGY Evir:");
	//double Evir = compute_energy(R*KS_nemo, R*Jnemo, Ex_vir);

	print("FINAL", oep_param.model(), "ENERGY Econv:");
	double Econv = compute_energy(R*KS_nemo, R*Jnemo, Ex_conv);

	printf("      +++ FINAL TOTAL ENERGY = %15.8f  Eh +++\n\n\n", Econv);

	printf("     Ex_vir       = %15.8f  Eh", Ex_vir);
	printf("\n     Ex_conv      = %15.8f  Eh", Ex_conv);
	printf("\n     Ex_HF        = %15.8f  Eh\n", Ex_HF);

	printf("\n     Ekin_HF (T)  = %15.8f  Eh", Ekin_HF);
	printf("\n     Ekin_KS (Ts) = %15.8f  Eh\n", Ekin_KS);

	printf("\n     DEvir_14     = %15.8f mEh", (Ex_vir - Ex_conv)*1000.0); // like in Kohut_2014, equation (45)
	printf("\n     DEvir_17     = %15.8f mEh\n", (Ex_vir - Ex_HF - 2.0*Tc)*1000.0); // like in Ospadov_2017, equation (28)
	print("     Drho         =     ", Drho, "e\n\n");

	print("---------------------------------------------------------------------------");
	double E_0 = compute_E_zeroth(KS_eigvals);
	double E_1 = compute_E_first(R*KS_nemo, R*Jnemo, R*Knemo, Voep);

	printf("  E^(0)               = %15.8f  Eh", E_0);
	printf("\n  E^(1)               = %15.8f  Eh", E_1);
	printf("\n  E^(0) + E^(1)       = %15.8f  Eh", E_0 + E_1);
	printf("\n  difference to Econv = %15.8f mEh\n\n", (E_0 + E_1 - Econv)*1000.0);

	print("saving orbitals to restartdata");
	Tensor<double> f_pp = compute_fock_diagonal_elements(calc->aeps, KS_nemo, Knemo, Voep);

	print("KS Fock matrix elements ", calc->aeps);
	print("HF Fock matrix elements ", f_pp);

	calc->aeps = f_pp;
	if (calc->param.save()) calc->save_mos(world);

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
    const real_function_3d Vs = compute_slater_potential(HF_nemo, homo_ind(HF_eigvals));
    print("  Slater potential computed successfully\n");

    print("test construction of IHF");
    const real_function_3d IHF = compute_energy_weighted_density(HF_nemo, HF_eigvals);
    print("  compute_energy_weighted_density computed successfully\n");

    print("test construction of kin_tot_HF (tau/rho HF)");
    const real_function_3d kin_tot_HF = compute_total_kinetic_density(HF_nemo, HF_eigvals);
    print("  kin_tot_HF computed successfully\n");

    print("test construction of kin_P_HF (tau_P/rho HF)");
    const real_function_3d kin_P_HF = compute_Pauli_kinetic_density(HF_nemo, HF_eigvals);
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
    const double Exvir_HF = compute_exchange_energy_vir(R*HF_nemo, Vs);
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
    const double Exvir_HFocep = compute_exchange_energy_vir(R*HF_nemo, V_HFocep);
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
    const double Exvir_HFdcep = compute_exchange_energy_vir(R*HF_nemo, V_HFdcep);
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
    const double Exvir_HFmrks = compute_exchange_energy_vir(R*HF_nemo, V_HFmrks);
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
    const double Ekin_HF = compute_kinetic_energy(R*HF_nemo);
    const double Ekin_HF_diff = fabs(Ekin_HF_correct - Ekin_HF);
    print("     the HF kinetic energy of the system is ...", Ekin_HF, "Eh");
    print("     error:", Ekin_HF_diff, "Eh");
    if (Ekin_HF_diff <= param.econv()) print("  HF kinetic energy is correct\n");
    else {
    	print("  ATTENTION: HF kinetic energy error is larger than energy convergence threshold (econv)!\n");
    	everything_ok = false;
    }

    print("\n   >> test solve_oep function with mRKS model for 2 iterations\n");
    solve_oep(HF_nemo, HF_eigvals);
    print("\n   >> solve_oep test finished, calculating test quantities based on the new KS orbitals and eigenvalues\n");

    vecfuncT& KS_nemo = calc->amo;
    //tensorT& KS_eigvals = calc->aeps;

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
    const double Ex_vir = compute_exchange_energy_vir(R*KS_nemo, Voep);
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
    const double Etot = compute_energy(R*KS_nemo, R*Jnemo, Ex_conv);
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
    const double Ekin_KS = compute_kinetic_energy(R*KS_nemo);
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
