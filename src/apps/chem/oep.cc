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
	const real_function_3d IHF = compute_average_I(HF_nemo, HF_eigvals);
	const real_function_3d kin_tot_HF = compute_total_kinetic_density(HF_nemo, HF_eigvals);
	const real_function_3d kin_P_HF = compute_Pauli_kinetic_density(HF_nemo, HF_eigvals);
	const real_function_3d rho_HF = compute_density(HF_nemo);
	if (oep_param.saving_amount() >= 1) save(Vs, "Slaterpotential");
	if (oep_param.saving_amount() >= 2) {
		save(rho_HF, "density_HF");
        if (oep_param.is_ocep() or oep_param.is_dcep() or oep_param.is_mrks()) save(IHF, "IHF");
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
    	for (long i = 0; i < HF_nemo.size(); i++) {
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
	for (int iter = 0; iter < calc->param.maxiter(); ++iter) {
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
			corr_ocep = compute_oep_correction("ocep", IHF, KS_nemo, KS_eigvals);
			if (oep_param.is_dcep()) corr_dcep = compute_oep_correction("dcep", kin_tot_HF, KS_nemo, KS_eigvals);
			if (oep_param.is_mrks()) corr_mrks = compute_oep_correction("mrks", kin_P_HF, KS_nemo, KS_eigvals);

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
			    	for (long i = 0; i < KS_nemo.size(); i++) {
			    		save(R*KS_nemo[i], "KS_orb_" + stringify(i) + "_iter_" + stringify(iter_counter));
			    	}
				}
			}
			if (oep_param.save_iter_density() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_density() == 0) {
					save(compute_density(KS_nemo), "density_iter_" + stringify(iter_counter));
				}
			}
			if (oep_param.save_iter_IKS() > 0) {
				if (iter_counter == 2 or iter_counter % oep_param.save_iter_IKS() == 0) {
					save(compute_average_I(KS_nemo, KS_eigvals), "IKS_iter_" + stringify(iter_counter));
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
        normalize(KS_nemo);

		// calculate new orbital energies (current eigenvalues from Fock-matrix)
		for (int i = 0; i < KS_nemo.size(); ++i) {
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
		normalize(nemo_new);

		// What is step restriction?
		calc->do_step_restriction(world, KS_nemo, nemo_new, "ab spin case");
		orthonormalize(nemo_new);
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
	real_function_3d IKS = compute_average_I(KS_nemo, KS_eigvals);
	real_function_3d kin_tot_KS = compute_total_kinetic_density(KS_nemo, KS_eigvals);
	real_function_3d kin_P_KS = compute_Pauli_kinetic_density(KS_nemo, KS_eigvals);
	real_function_3d rho_KS = compute_density(KS_nemo);
	double Drho = compute_delta_rho(rho_HF, rho_KS);
	if (oep_param.saving_amount() >= 1) save(rho_KS, "density_final");
	if (oep_param.saving_amount() >= 2) {
        if (oep_param.is_ocep() or oep_param.is_dcep() or oep_param.is_mrks()) save(IKS, "IKS_final");
        if (oep_param.is_dcep()) save(kin_tot_KS, "kin_tot_KS_final");
        if (oep_param.is_mrks()) save(kin_P_KS, "kin_P_KS_final");
    	for (long i = 0; i < KS_nemo.size(); i++) {
    		save(R*KS_nemo[i], "KS_orb_" + stringify(i) + "_final");
    	}
	}

	// if desired: print final KS orbital contributions to total density (nemo squares)
	if (oep_param.saving_amount() >= 3) {
    	vecfuncT KS_nemo_square = square(world, KS_nemo);
    	for (long i = 0; i < KS_nemo_square.size(); i++) {
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
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", IHF, KS_nemo, KS_eigvals);
    	Voep = Vs + ocep_correction_final + shift_final;
    	if (oep_param.saving_amount() >= 1) {
    		save(ocep_correction_final + shift_final, "OCEP_correction_final");
    		save(Voep, "OEPapprox_final");
    	}
	}
	if (oep_param.is_dcep()) {
		print("\n  computing final DCEP with converged DCEP orbitals and eigenvalues");
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", IHF, KS_nemo, KS_eigvals);
    	real_function_3d dcep_correction_final = compute_oep_correction("dcep", kin_tot_HF, KS_nemo, KS_eigvals);
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
    	real_function_3d ocep_correction_final = compute_oep_correction("ocep", IHF, KS_nemo, KS_eigvals);
    	real_function_3d mrks_correction_final = compute_oep_correction("mrks", kin_P_HF, KS_nemo, KS_eigvals);
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
	double Evir = compute_energy(R*KS_nemo, R*Jnemo, Ex_vir);

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



} /* namespace madness */
