/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
  \file examples/oep.cc
  \brief optimized effective potentials for DFT
*/

#include <chem/nemo.h>
#include <chem/cheminfo.h>
#include <chem/SCFOperators.h>
#include <chem/projector.h>

using namespace madness;

struct dens_inv{

    /// @param[out] U   result
    /// @param[in]  t   numerator
    /// @param[in]  inv density to be inverted >0.0
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& inv) const {
        ITERATOR(
            U, double d = t(IND);
            double p = std::max(inv(IND), 1.e-8);
            U(IND) = d/p;
        );
   }
    template <typename Archive>
    void serialize(Archive& ar) {}

};

struct binary_munge{

	double longrangevalue;
	double threshold;

	binary_munge(const double lrv = 0.0, const double thresh = 1.0e-8) {
		longrangevalue = lrv;
		threshold = thresh;
	}

    /// @param[out] U   result
    /// @param[in]  f   function to be munged
    /// @param[in]  r   refdensity
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& f,
            const Tensor<double>& refdens) const {
        ITERATOR(
            U, double r = refdens(IND);
            double ff = f(IND);
            U(IND) = (r > threshold) ? ff : longrangevalue;
        );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}

};

/*
struct unary_munge{

    /// @param[out] U   result
    /// @param[in]  f   function to be munged
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& f) const {
        ITERATOR(
            U, double ff = f(IND);
            U(IND) = (ff > 1.e-8) ? ff : 0.0;
        );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}

};
*/

/// simple structure to take the pointwise logarithm of a function, shifted by +14
struct logme{
    typedef double resultT;
    struct logme1 {
        double operator()(const double& val) {return log(std::max(1.e-14,val))+14.0;}
    };
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& val) const {
        Tensor<double> result=copy(val);
        logme1 op;
        return result.unaryop(op);
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
};

class OEP : public Nemo {
	typedef std::shared_ptr<real_convolution_3d> poperatorT;

public:
    OEP(World& world, const std::shared_ptr<SCF> calc) : Nemo(world, calc) {}

    void solve_oaep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

        // Iterative energy calculation for OAEP with EXACT EXCHANGE functional
    	// for other functions, slater potential must be modified
    	// HF orbitals and eigenvectors are used as the guess here
    	// note that oaep_nemo is a reference and changes oep->get_calc()->amo orbitals
    	// same for orbital energies (eigenvalues) oaep_eigvals which is oep->get_calc()->aeps

    	double energy = 0.0;
    	bool converged = false;

		// compute OAEP Slater potential Vs and average IHF from HF orbitals and eigenvalues
    	const real_function_3d Vs = compute_slater_potential(HF_nemo);
    	const real_function_3d IHF = compute_average_I(HF_nemo, HF_eigvals);
    	save(IHF, "IHF");

    	// set oaep_nemo as reference to MOs
    	vecfuncT& oaep_nemo = calc->amo;
    	tensorT& oaep_eigvals = calc->aeps; // 1d tensor of same length as oaep_nemo

    	// all necessary operators applied on nemos (Knemo is used after the cycle):
    	vecfuncT Jnemo, Unemo, Slaternemo, Knemo;

    	//define the solver
    	typedef allocator<double, 3> allocT;
    	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
    	allocT alloc(world, oaep_nemo.size());
    	solverT solver(allocT(world, oaep_nemo.size()));


    	// iterations until energy is self-consistent:
    	for (int iter = 0; iter < calc->param.maxiter; ++iter) {
    		save_function(oaep_nemo, "oaep_nemo_it"+stringify(iter));
    		vecfuncT R2oaep_nemo = R_square*oaep_nemo;
    		truncate(world, R2oaep_nemo);

    		// compute parts of the Fock matrix J, Unuc and Vs
    		compute_oaep_nemo_potentials(oaep_nemo, Jnemo, Unemo, Vs, Slaternemo);

    		// compute Fock matrix F = J + Vs + Vnuc and kinetic energy
    		vecfuncT Fnemo = Jnemo + Slaternemo + Unemo;
    		truncate(world, Fnemo);
    		tensorT F = matrix_inner(world, R2oaep_nemo, Fnemo, false); // matrix_inner gives 2d tensor
    		Kinetic<double,3> T(world);
    		F += T(R2oaep_nemo, oaep_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

    		// report the off-diagonal Fock matrix elements because canonical orbitals are used
            tensorT F_offdiag = copy(F);
            for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
            double max_F_offidag = F_offdiag.absmax();
            if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

    		// compute new (current) energy
            double old_energy = energy;
    		energy = compute_energy(R*oaep_nemo, R*Jnemo, Vs, Knemo, true); // Knemo is not used here
    		// compute_exchange_potential(oaep_nemo, Knemo);
    		// energy = compute_energy(R*oaep_nemo, R*Jnemo, Vs, R*Knemo, false);
    		// there should be no difference between these two methods, because energy is only needed
    		// for checking convergence threshold; but: Evir should be much faster because K takes time

            // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
    		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
            tensorT X; // must be formed from R*nemos but can then be used for nemos also
            tensorT overlap = matrix_inner(world, R*oaep_nemo, R*oaep_nemo, true);
            X = calc->get_fock_transformation(world, overlap, F, oaep_eigvals, calc->aocc,
            		FunctionDefaults<3>::get_thresh());
            oaep_nemo = transform(world, oaep_nemo, X, trantol(), true);
            rotate_subspace(world, X, solver, 0, oaep_nemo.size());

            truncate(world, oaep_nemo);
            normalize(oaep_nemo);

    		// update the nemos
    		// set and modify orbital energies (current eigenvalues from Fock-matrix)
    		for (int i = 0; i < oaep_nemo.size(); ++i) {
    			oaep_eigvals(i) = std::min(-0.05, F(i, i));
    			// orbital energy is set to -0.05 if it was above
    		}

    		// if requested: subtract orbital shift from orbital energies
    		if (calc->param.orbitalshift > 0.0) {
    			if (world.rank() == 0) print("shifting orbitals by ",
    					calc->param.orbitalshift, " to lower energies");
    			oaep_eigvals -= calc->param.orbitalshift;
    		}

    		// print orbital energies:
    		printf("\norbital energies of iteration %3u\n", iter);
    		for (long i = oaep_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, oaep_eigvals(i));
    		}

    		// construct the BSH operators ops
    		std::vector<poperatorT> G = calc->make_bsh_operators(world, oaep_eigvals);

    		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
    		Fnemo = transform(world, Fnemo, X, trantol(), true);
    		truncate(world, Fnemo);

    		// apply the BSH operators G (here ops) on the wave function
    		scale(world, Fnemo, -2.0);
    		vecfuncT GFnemo = apply(world, G, Fnemo);
    		truncate(world, GFnemo);

    		double n1 = norm2(world, oaep_nemo);
    		double n2 = norm2(world, GFnemo);
    		print("\nnorm of oaep_nemo and GFnemo, ratio ", n1, n2, n1/n2);

    		// compute the residuals
    		vecfuncT residual = oaep_nemo - GFnemo;
    		const double norm = norm2(world, residual) / sqrt(oaep_nemo.size());

    		// What is happening here? (KAIN)
    		vecfuncT nemo_new;
    		if (norm < 5.0e-1) {
    			nemo_new = (solver.update(oaep_nemo, residual)).x;
    		} else {
    			nemo_new = GFnemo;
    		}
    		truncate(world, nemo_new);
    		normalize(nemo_new);

    		// What is step restriction?
    		calc->do_step_restriction(world, oaep_nemo, nemo_new, "ab spin case");
    		orthonormalize(nemo_new);
    		oaep_nemo = nemo_new;

    		// evaluate convergence via norm error and energy difference
    		if ((norm < calc->param.dconv) and (fabs(energy - old_energy) < calc->param.econv))
    			converged = true;

    		if (calc->param.save) calc->save_mos(world);

    		if (world.rank() == 0) {
    			printf("finished iteration %2d at time %8.1fs with energy %12.8f\n", iter, wall_time(), energy);
    			print("current residual norm ", norm, "\n");
    		}

    		if (converged) break;

    	}

    	if (converged) {
    		if (world.rank() == 0) print("\nIterations converged\n");
    	}
    	else {
    		if (world.rank() == 0) print("\nIterations failed\n");
    		energy = 0.0;
    	}

    	print("\nFINAL OAEP ENERGY Evir:");
    	double Evir = compute_energy(R*oaep_nemo, R*Jnemo, Vs, Knemo, true); // Knemo is not used here

    	print("FINAL OAEP ENERGY Econv:");
    	compute_exchange_potential(oaep_nemo, Knemo);
    	double Econv = compute_energy(R*oaep_nemo, R*Jnemo, Vs, R*Knemo, false); // Vs is not used here

    	printf("OAEP  Evir = %15.8f  Eh", Evir);
    	printf("\nOAEP Econv = %15.8f  Eh", Econv);
    	printf("\nOAEP DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

    	printf("\n  computing V_OCEP with converged OAEP orbitals and eigenvalues\n");
    	real_function_3d IKS = compute_average_I(oaep_nemo, oaep_eigvals);
    	save(IKS, "IKS_OAEP");
    	real_function_3d Vocep = Vs + compute_OCEP_correction(HF_nemo, HF_eigvals, oaep_nemo, oaep_eigvals);
    	printf("     done\n");

    }


    void solve_oaep_2(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

        // Iterative energy calculation for OAEP with EXACT EXCHANGE functional
    	// for other functions, slater potential must be modified
    	// HF orbitals and eigenvectors are used as the guess here
    	// note that oaep_nemo is a reference and changes oep->get_calc()->amo orbitals
    	// same for orbital energies (eigenvalues) oaep_eigvals which is oep->get_calc()->aeps

    	double energy = 0.0;
    	bool converged = false;

		// compute OAEP Slater potential Vs and average IHF from HF orbitals and eigenvalues
    	const real_function_3d Vs = compute_slater_potential(HF_nemo);
    	const real_function_3d IHF = compute_average_I(HF_nemo, HF_eigvals);
    	save(IHF, "IHF");

    	// set oaep_nemo as reference to MOs
    	vecfuncT& oaep_nemo = calc->amo;
    	tensorT& oaep_eigvals = calc->aeps; // 1d tensor of same length as oaep_nemo

    	// all necessary operators applied on nemos (Knemo is used after the cycle):
    	vecfuncT Jnemo, Unemo, Slaternemo, Knemo;

    	//define the solver
    	typedef allocator<double, 3> allocT;
    	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
    	allocT alloc(world, oaep_nemo.size());
    	solverT solver(allocT(world, oaep_nemo.size()));


    	// iterations until energy is self-consistent:
    	for (int iter = 0; iter < calc->param.maxiter; ++iter) {
    		for (int i = 0; i < oaep_nemo.size(); i++) {
    			save(oaep_nemo[i], "oaep_nemo_"+stringify(i)+"_it_"+stringify(iter));
    		}
    		vecfuncT R2oaep_nemo = R_square*oaep_nemo;
    		truncate(world, R2oaep_nemo);

    		// compute parts of the Fock matrix J, Unuc and Vs
    		compute_oaep_nemo_potentials(oaep_nemo, Jnemo, Unemo, Vs, Slaternemo);

    		// compute Fock matrix F = J + Vs + Vnuc and kinetic energy
    		vecfuncT Fnemo = Jnemo + Slaternemo + Unemo;
    		truncate(world, Fnemo);
    		tensorT F = matrix_inner(world, R2oaep_nemo, Fnemo, false); // matrix_inner gives 2d tensor
    		Kinetic<double,3> T(world);
    		F += T(R2oaep_nemo, oaep_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

    		// report the off-diagonal Fock matrix elements because canonical orbitals are used
            tensorT F_offdiag = copy(F);
            for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
            double max_F_offidag = F_offdiag.absmax();
            if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

    		// compute new (current) energy
            double old_energy = energy;
    		energy = compute_energy(R*oaep_nemo, R*Jnemo, Vs, Knemo, true); // Knemo is not used here
    		// compute_exchange_potential(oaep_nemo, Knemo);
    		// energy = compute_energy(R*oaep_nemo, R*Jnemo, Vs, R*Knemo, false);
    		// there should be no difference between these two methods, because energy is only needed
    		// for checking convergence threshold; but: Evir should be much faster because K takes time

            // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
    		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
            tensorT X; // must be formed from R*nemos but can then be used for nemos also
            tensorT overlap = matrix_inner(world, R*oaep_nemo, R*oaep_nemo, true);
            X = calc->get_fock_transformation(world, overlap, F, oaep_eigvals, calc->aocc,
            		FunctionDefaults<3>::get_thresh());
            oaep_nemo = transform(world, oaep_nemo, X, trantol(), true);
            rotate_subspace(world, X, solver, 0, oaep_nemo.size());

            truncate(world, oaep_nemo);
            normalize(oaep_nemo);

    		// update the nemos
    		// set and modify orbital energies (current eigenvalues from Fock-matrix)
    		for (int i = 0; i < oaep_nemo.size(); ++i) {
    			oaep_eigvals(i) = std::min(-0.05, F(i, i));
    			// orbital energy is set to -0.05 if it was above
    		}

    		// if requested: subtract orbital shift from orbital energies
    		if (calc->param.orbitalshift > 0.0) {
    			if (world.rank() == 0) print("shifting orbitals by ",
    					calc->param.orbitalshift, " to lower energies");
    			oaep_eigvals -= calc->param.orbitalshift;
    		}

    		// print orbital energies:
    		printf("\norbital energies of iteration %3u\n", iter);
    		for (long i = oaep_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, oaep_eigvals(i));
    		}

    		// construct the BSH operators ops
    		std::vector<poperatorT> G = calc->make_bsh_operators(world, oaep_eigvals);

    		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
    		Fnemo = transform(world, Fnemo, X, trantol(), true);
    		truncate(world, Fnemo);

    		// apply the BSH operators G (here ops) on the wave function
    		scale(world, Fnemo, -2.0);
    		vecfuncT GFnemo = apply(world, G, Fnemo);
    		truncate(world, GFnemo);

    		double n1 = norm2(world, oaep_nemo);
    		double n2 = norm2(world, GFnemo);
    		print("\nnorm of oaep_nemo and GFnemo, ratio ", n1, n2, n1/n2);

    		// compute the residuals
    		vecfuncT residual = oaep_nemo - GFnemo;
    		const double norm = norm2(world, residual) / sqrt(oaep_nemo.size());

    		// What is happening here? (KAIN)
    		vecfuncT nemo_new;
    		if (norm < 5.0e-1) {
    			nemo_new = (solver.update(oaep_nemo, residual)).x;
    		} else {
    			nemo_new = GFnemo;
    		}
    		truncate(world, nemo_new);
    		normalize(nemo_new);

    		// What is step restriction?
    		calc->do_step_restriction(world, oaep_nemo, nemo_new, "ab spin case");
    		orthonormalize(nemo_new);
    		oaep_nemo = nemo_new;

    		// evaluate convergence via norm error and energy difference
    		if ((norm < calc->param.dconv) and (fabs(energy - old_energy) < calc->param.econv))
    			converged = true;

    		if (calc->param.save) calc->save_mos(world);

    		if (world.rank() == 0) {
    			printf("finished iteration %2d at time %8.1fs with energy %12.8f\n", iter, wall_time(), energy);
    			print("current residual norm ", norm, "\n");
    		}

    		if (converged) break;

    	}

    	if (converged) {
    		if (world.rank() == 0) print("\nIterations converged\n");
    	}
    	else {
    		if (world.rank() == 0) print("\nIterations failed\n");
    		energy = 0.0;
    	}

    	print("\nFINAL OAEP ENERGY Evir:");
    	double Evir = compute_energy(R*oaep_nemo, R*Jnemo, Vs, Knemo, true); // Knemo is not used here

    	print("FINAL OAEP ENERGY Econv:");
    	compute_exchange_potential(oaep_nemo, Knemo);
    	double Econv = compute_energy(R*oaep_nemo, R*Jnemo, Vs, R*Knemo, false); // Vs is not used here

    	printf("OAEP  Evir = %15.8f  Eh", Evir);
    	printf("\nOAEP Econv = %15.8f  Eh", Econv);
    	printf("\nOAEP DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

    	printf("\n  computing V_OCEP with converged OAEP orbitals and eigenvalues\n");
    	real_function_3d IKS = compute_average_I(oaep_nemo, oaep_eigvals);
    	save(IKS, "IKS_OAEP");
    	real_function_3d correction = IHF - IKS;
    	save(correction, "correction_old");

    	real_function_3d rho = 2.0*R_square*dot(world, oaep_nemo, oaep_nemo);
    	save(rho, "density");
    	double homo_diff = oaep_eigvals(oaep_eigvals.size() - 1) - HF_eigvals(HF_eigvals.size() - 1);
    	printf("   homo_diff = %11.8f\n", homo_diff);
    	real_function_3d correction_munged = binary_op(correction, rho, binary_munge(homo_diff, 1.0e-5));
    	save(correction_munged, "correction_old_munged");

    	const real_function_3d Vocep = Vs + compute_OCEP_correction(HF_nemo, HF_eigvals, oaep_nemo, oaep_eigvals);
    	save(Vocep, "OCEPpotential_old");
    	printf("     done\n\n");



    	// new iteration: same as above, but Vocep instead of Vs and current nemos as start

    	converged = false;
    	vecfuncT OCEPnemo;

    	// iterations until energy is self-consistent:
    	for (int iter = 0; iter < calc->param.maxiter; ++iter) {
    		for (int i = 0; i < oaep_nemo.size(); i++) {
    			save(oaep_nemo[i], "ocep_nemo_"+stringify(i)+"_it_"+stringify(iter));
    		}
    		vecfuncT R2oaep_nemo = R_square*oaep_nemo;
    		truncate(world, R2oaep_nemo);

    		// compute parts of the Fock matrix J, Unuc and Vocep
    		compute_oaep_nemo_potentials(oaep_nemo, Jnemo, Unemo, Vocep, OCEPnemo);

    		// compute Fock matrix F = J + Vs + Vnuc and kinetic energy
    		vecfuncT Fnemo = Jnemo + OCEPnemo + Unemo;
    		truncate(world, Fnemo);
    		tensorT F = matrix_inner(world, R2oaep_nemo, Fnemo, false); // matrix_inner gives 2d tensor
    		Kinetic<double,3> T(world);
    		F += T(R2oaep_nemo, oaep_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

    		// report the off-diagonal Fock matrix elements because canonical orbitals are used
            tensorT F_offdiag = copy(F);
            for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
            double max_F_offidag = F_offdiag.absmax();
            if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

    		// compute new (current) energy
            double old_energy = energy;
    		energy = compute_energy(R*oaep_nemo, R*Jnemo, Vocep, Knemo, true); // Knemo is not used here
    		// compute_exchange_potential(oaep_nemo, Knemo);
    		// energy = compute_energy(R*oaep_nemo, R*Jnemo, Vs, R*Knemo, false);
    		// there should be no difference between these two methods, because energy is only needed
    		// for checking convergence threshold; but: Evir should be much faster because K takes time

            // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
    		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
            tensorT X; // must be formed from R*nemos but can then be used for nemos also
            tensorT overlap = matrix_inner(world, R*oaep_nemo, R*oaep_nemo, true);
            X = calc->get_fock_transformation(world, overlap, F, oaep_eigvals, calc->aocc,
            		FunctionDefaults<3>::get_thresh());
            oaep_nemo = transform(world, oaep_nemo, X, trantol(), true);
            rotate_subspace(world, X, solver, 0, oaep_nemo.size());

            truncate(world, oaep_nemo);
            normalize(oaep_nemo);

    		// update the nemos
    		// set and modify orbital energies (current eigenvalues from Fock-matrix)
    		for (int i = 0; i < oaep_nemo.size(); ++i) {
    			oaep_eigvals(i) = std::min(-0.05, F(i, i));
    			// orbital energy is set to -0.05 if it was above
    		}

    		// if requested: subtract orbital shift from orbital energies
    		if (calc->param.orbitalshift > 0.0) {
    			if (world.rank() == 0) print("shifting orbitals by ",
    					calc->param.orbitalshift, " to lower energies");
    			oaep_eigvals -= calc->param.orbitalshift;
    		}

    		// print orbital energies:
    		printf("\norbital energies of iteration %3u\n", iter);
    		for (long i = oaep_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, oaep_eigvals(i));
    		}

    		// construct the BSH operators ops
    		std::vector<poperatorT> G = calc->make_bsh_operators(world, oaep_eigvals);

    		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
    		Fnemo = transform(world, Fnemo, X, trantol(), true);
    		truncate(world, Fnemo);

    		// apply the BSH operators G (here ops) on the wave function
    		scale(world, Fnemo, -2.0);
    		vecfuncT GFnemo = apply(world, G, Fnemo);
    		truncate(world, GFnemo);

    		double n1 = norm2(world, oaep_nemo);
    		double n2 = norm2(world, GFnemo);
    		print("\nnorm of oaep_nemo and GFnemo, ratio ", n1, n2, n1/n2);

    		// compute the residuals
    		vecfuncT residual = oaep_nemo - GFnemo;
    		const double norm = norm2(world, residual) / sqrt(oaep_nemo.size());

    		// What is happening here? (KAIN)
    		vecfuncT nemo_new;
    		if (norm < 5.0e-1) {
    			nemo_new = (solver.update(oaep_nemo, residual)).x;
    		} else {
    			nemo_new = GFnemo;
    		}
    		truncate(world, nemo_new);
    		normalize(nemo_new);

    		// What is step restriction?
    		calc->do_step_restriction(world, oaep_nemo, nemo_new, "ab spin case");
    		orthonormalize(nemo_new);
    		oaep_nemo = nemo_new;

    		// evaluate convergence via norm error and energy difference
    		if ((norm < calc->param.dconv) and (fabs(energy - old_energy) < calc->param.econv))
    			converged = true;

    		if (calc->param.save) calc->save_mos(world);

    		if (world.rank() == 0) {
    			printf("finished iteration %2d at time %8.1fs with energy %12.8f\n", iter, wall_time(), energy);
    			print("current residual norm ", norm, "\n");
    		}

    		if (converged) break;

    	}

    	if (converged) {
    		if (world.rank() == 0) print("\nIterations converged\n");
    	}
    	else {
    		if (world.rank() == 0) print("\nIterations failed\n");
    		energy = 0.0;
    	}


    	print("\nFINAL OAEP ENERGY Evir:");
    	Evir = compute_energy(R*oaep_nemo, R*Jnemo, Vocep, Knemo, true); // Knemo is not used here

    	print("FINAL OAEP ENERGY Econv:");
    	compute_exchange_potential(oaep_nemo, Knemo);
    	Econv = compute_energy(R*oaep_nemo, R*Jnemo, Vocep, R*Knemo, false); // Vs is not used here

    	printf("OAEP  Evir = %15.8f  Eh", Evir);
    	printf("\nOAEP Econv = %15.8f  Eh", Econv);
    	printf("\nOAEP DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

    	printf("\n  computing V_OCEP with new converged orbitals and eigenvalues\n");
    	real_function_3d IKS_new = compute_average_I(oaep_nemo, oaep_eigvals);
    	save(IKS_new, "IKS_OAEP_new");
    	correction = IHF - IKS_new;
    	save(correction, "correction_new");
    	real_function_3d Vocep_new = Vs + compute_OCEP_correction(HF_nemo, HF_eigvals, oaep_nemo, oaep_eigvals);
    	save(Vocep_new, "OCEPpotential_new");
    	printf("     done\n");


    }


    void solve_ocep(const vecfuncT& HF_nemo, const tensorT& HF_eigvals) {

        // Iterative energy calculation for OCEP with EXACT EXCHANGE functional
    	// for other functions, slater potential must be modified
    	// HF orbitals and eigenvectors are used as the guess here
    	// note that ocep_nemo is a reference and changes oep->get_calc()->amo orbitals
    	// same for orbital energies (eigenvalues) ocep_eigvals which is oep->get_calc()->aeps

    	double energy = 0.0;
    	bool converged = false;

		// compute OAEP Slater potential Vs and average IHF from HF orbitals and eigenvalues
    	const real_function_3d Vs = compute_slater_potential(HF_nemo);
    	const real_function_3d IHF = compute_average_I(HF_nemo, HF_eigvals);
    	save(IHF, "IHF");

    	// set ocep_nemo as reference to MOs
    	vecfuncT& ocep_nemo = calc->amo;
    	tensorT& ocep_eigvals = calc->aeps; // 1d tensor of same length as ocep_nemo

    	// all necessary operators applied on nemos (Knemo is used after the cycle):
    	vecfuncT Jnemo, Unemo, OCEPnemo, Knemo;
    	real_function_3d Vocep = Vs;

    	//define the solver
    	typedef allocator<double, 3> allocT;
    	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
    	allocT alloc(world, ocep_nemo.size());
    	solverT solver(allocT(world, ocep_nemo.size()));


    	// iterations until energy is self-consistent:
    	for (int iter = 0; iter < calc->param.maxiter; ++iter) {
    		// save_function(ocep_nemo, "ocep_nemo_it_"+stringify(iter));
    		vecfuncT R2ocep_nemo = R_square*ocep_nemo;
    		truncate(world, R2ocep_nemo);

    		real_function_3d rho = 2.0*R_square*dot(world, ocep_nemo, ocep_nemo);  // 2 -> closed shell
    		save(rho, "density_it_"+stringify(iter));

    		// compute OCEP potential from current nemos and eigenvalues
    		// like Kohut, 2014, equation (26) with correction = IHF - IKS
    		// only in every 3rd iteration
    		if ((iter % 3 == 0) and (iter > 2)) {
    			real_function_3d corr = compute_OCEP_correction(HF_nemo, HF_eigvals, ocep_nemo, ocep_eigvals);
    			Vocep = Vs + corr;
    			save(corr, "OCEP_correction_it_"+stringify(iter));
    		}
    		save(Vocep, "OCEP_potential_it_"+stringify(iter));

    		// compute parts of the Fock matrix J, Unuc and Vs
    		compute_ocep_nemo_potentials(ocep_nemo, Jnemo, Unemo, Vocep, OCEPnemo);

    		// compute Fock matrix F = J + Vs + Vnuc and kinetic energy
    		vecfuncT Fnemo = Jnemo + OCEPnemo + Unemo;
    		truncate(world, Fnemo);
    		tensorT F = matrix_inner(world, R2ocep_nemo, Fnemo, false); // matrix_inner gives 2d tensor
    		Kinetic<double,3> T(world);
    		F += T(R2ocep_nemo, ocep_nemo); // 2d tensor = Fock-matrix  // R_square in bra, no R in ket

    		// report the off-diagonal Fock matrix elements because canonical orbitals are used
            tensorT F_offdiag = copy(F);
            for (int i = 0; i < F.dim(0); ++i) F_offdiag(i, i) = 0.0;
            double max_F_offidag = F_offdiag.absmax();
            if (world.rank() == 0) print("F max off-diagonal ", max_F_offidag);

    		// compute new (current) energy
            double old_energy = energy;
    		energy = compute_energy(R*ocep_nemo, R*Jnemo, Vocep, Knemo, true); // Knemo is not used here
    		// compute_exchange_potential(ocep_nemo, Knemo);
    		// energy = compute_energy(R*ocep_nemo, R*Jnemo, Vs, R*Knemo, false);
    		// there should be no difference between these two methods, because energy is only needed
    		// for checking convergence threshold; but: Evir should be much faster because K takes time

            // diagonalize the Fock matrix to get the eigenvalues and eigenvectors (canonical)
    		// FC = epsilonSC and X^dSX with transform matrix X, see Szabo/Ostlund (3.159) and (3.165)
            tensorT X; // must be formed from R*nemos but can then be used for nemos also
            tensorT overlap = matrix_inner(world, R*ocep_nemo, R*ocep_nemo, true);
            X = calc->get_fock_transformation(world, overlap, F, ocep_eigvals, calc->aocc,
            		FunctionDefaults<3>::get_thresh());
            ocep_nemo = transform(world, ocep_nemo, X, trantol(), true);
            rotate_subspace(world, X, solver, 0, ocep_nemo.size());

            truncate(world, ocep_nemo);
            normalize(ocep_nemo);

    		// update the nemos
    		// set and modify orbital energies (current eigenvalues from Fock-matrix)
    		for (int i = 0; i < ocep_nemo.size(); ++i) {
    			ocep_eigvals(i) = std::min(-0.05, F(i, i));
    			// orbital energy is set to -0.05 if it was above
    		}

    		// if requested: subtract orbital shift from orbital energies
    		if (calc->param.orbitalshift > 0.0) {
    			if (world.rank() == 0) print("shifting orbitals by ",
    					calc->param.orbitalshift, " to lower energies");
    			ocep_eigvals -= calc->param.orbitalshift;
    		}

    		// print orbital energies:
    		printf("\norbital energies of iteration %3u\n", iter);
    		for (long i = ocep_eigvals.size() - 1; i >= 0; i--) {
    			printf(" e%2.2lu = %12.8f\n", i, ocep_eigvals(i));
    		}

    		// construct the BSH operators ops
    		std::vector<poperatorT> G = calc->make_bsh_operators(world, ocep_eigvals);

    		// remember Fock matrix * nemos from above; make sure it's in phase with nemo (transform)
    		Fnemo = transform(world, Fnemo, X, trantol(), true);
    		truncate(world, Fnemo);

    		// apply the BSH operators G (here ops) on the wave function
    		scale(world, Fnemo, -2.0);
    		vecfuncT GFnemo = apply(world, G, Fnemo);
    		truncate(world, GFnemo);

    		double n1 = norm2(world, ocep_nemo);
    		double n2 = norm2(world, GFnemo);
    		print("\nnorm of ocep_nemo and GFnemo, ratio ", n1, n2, n1/n2);

    		// compute the residuals
    		vecfuncT residual = ocep_nemo - GFnemo;
    		const double norm = norm2(world, residual) / sqrt(ocep_nemo.size());

    		// What is happening here? (KAIN)
    		vecfuncT nemo_new;
    		if (norm < 5.0e-1) {
    			nemo_new = (solver.update(ocep_nemo, residual)).x;
    		} else {
    			nemo_new = GFnemo;
    		}
    		truncate(world, nemo_new);
    		normalize(nemo_new);

    		// What is step restriction?
    		calc->do_step_restriction(world, ocep_nemo, nemo_new, "ab spin case");
    		orthonormalize(nemo_new);
    		ocep_nemo = nemo_new;

    		// evaluate convergence via norm error and energy difference
    		if ((norm < calc->param.dconv) and (fabs(energy - old_energy) < calc->param.econv))
    			converged = true;

    		if (calc->param.save) calc->save_mos(world);

    		if (world.rank() == 0) {
    			printf("finished iteration %2d at time %8.1fs with energy %12.8f\n", iter, wall_time(), energy);
    			print("current residual norm ", norm, "\n");
    		}

    		if (converged) break;

    	}

    	if (converged) {
    		if (world.rank() == 0) print("\nIterations converged\n");
    	}
    	else {
    		if (world.rank() == 0) print("\nIterations failed\n");
    		energy = 0.0;
    	}

    	printf("\n  computing final IKS with converged OCEP orbitals and eigenvalues\n");
    	real_function_3d IKS = compute_average_I(ocep_nemo, ocep_eigvals);
    	save(IKS, "IKS");
    	printf("     done\n");
    	printf("\n  computing final V_OCEP with converged OCEP orbitals and eigenvalues\n");
    	Vocep = Vs + compute_OCEP_correction(HF_nemo, HF_eigvals, ocep_nemo, ocep_eigvals);
    	printf("     done\n");

    	print("\nFINAL OCEP ENERGY Evir:");
    	double Evir = compute_energy(R*ocep_nemo, R*Jnemo, Vocep, Knemo, true); // Knemo is not used here

    	print("FINAL OCEP ENERGY Econv:");
    	compute_exchange_potential(ocep_nemo, Knemo);
    	double Econv = compute_energy(R*ocep_nemo, R*Jnemo, Vocep, R*Knemo, false); // Vocep is not used here

    	printf("OCEP  Evir = %15.8f  Eh", Evir);
    	printf("\nOCEP Econv = %15.8f  Eh", Econv);
    	printf("\nOCEP DEvir = %15.8f mEh\n", (Evir - Econv)*1000);

    }


    // compute Slater potential as OAEP potential (Kohut, 2014, equation (15))
    real_function_3d compute_slater_potential(const vecfuncT& nemo) const {

        Exchange K(world, this, 0); // no - in K here, so factor -1 must be included at the end
        vecfuncT Knemo = K(nemo);
        real_function_3d numerator = dot(world, nemo, Knemo);
        real_function_3d rho = dot(world, nemo, nemo);
        // note that there is factor R_square and factor 2 (close shell) in numerator and rho,
        // but both cancel when dividing

        real_function_3d Vs = -1.0*binary_op(numerator, rho, dens_inv());
        save(Vs, "Slaterpotential");
        return Vs;

    }

    // compute average ionization energy I for OCEP potential (eigenvalues as 1d tensor)
    real_function_3d compute_average_I(const vecfuncT& nemo, const tensorT eigvals) const {

    	// transform tensor eigvals to vector epsilon
		std::vector<double> epsilon(eigvals.size());
		for (int i = 0; i < eigvals.size(); ++i) epsilon[i] = eigvals(i);

		vecfuncT nemo_square = square(world, nemo); // |nemo|^2
		scale(world, nemo_square, epsilon); // epsilon*|nemo|^2
		real_function_3d numerator = sum(world, nemo_square);
        real_function_3d rho = dot(world, nemo, nemo); // density
        // note that there is factor R_square and factor 2 (close shell) in numerator and rho,
        // but both cancel when dividing

        // like Kohut, 2014, equations (21) and (25)
        real_function_3d I = -1.0*binary_op(numerator, rho, dens_inv());
        return I;

    }

    // compute OCEP potential from OAEP Slater potential Vs
    real_function_3d compute_OCEP_correction(const vecfuncT& nemoHF, const tensorT eigvalsHF,
    		const vecfuncT& nemoKS, const tensorT eigvalsKS) const {

    	// Hartree-Fock and Kohn-Sham average ionization energie
    	real_function_3d IHF = compute_average_I(nemoHF, eigvalsHF);
    	real_function_3d IKS = compute_average_I(nemoKS, eigvalsKS);

    	// density with KS orbitals and HF/KS HOMO difference
    	real_function_3d rho = 2.0*R_square*dot(world, nemoKS, nemoKS); // 2 because closed shell
    	double homo_diff = eigvalsKS(eigvalsKS.size() - 1) - eigvalsHF(eigvalsHF.size() - 1);

    	real_function_3d correction = binary_op(IHF - IKS, rho, binary_munge(homo_diff, 1.0e-5));
    	// save(correction, "OCEP_correction");
    	return correction;

    }

    // compute all potentials from given nemos except kinetic energy
    void compute_oaep_nemo_potentials(const vecfuncT& nemo, vecfuncT& Jnemo, vecfuncT& Unemo,
    		const real_function_3d V, vecfuncT& OAEPnemo) const {

    	// compute Coulomb part
    	Coulomb J = Coulomb(world, this);
    	Jnemo = J(nemo);     // as in working equation for MRA
    	truncate(world, Jnemo);

    	// compute nuclear potential part
    	Nuclear Unuc(world, this->nuclear_correlation);
    	Unemo = Unuc(nemo);

    	// compute OAEP Slater potential part
    	OAEPnemo = V*nemo;

    }

    // compute all potentials from given nemos except kinetic energy
    void compute_ocep_nemo_potentials(const vecfuncT& nemo, vecfuncT& Jnemo, vecfuncT& Unemo,
    		const real_function_3d V, vecfuncT& OCEPnemo) const {

    	// compute Coulomb part
    	Coulomb J = Coulomb(world, this);
    	Jnemo = J(nemo);     // as in working equation for MRA
    	truncate(world, Jnemo);

    	// compute nuclear potential part
    	Nuclear Unuc(world, this->nuclear_correlation);
    	Unemo = Unuc(nemo);

    	// compute OCEP potential part
    	OCEPnemo = V*nemo;

    }

    // compute exchange potential (needed for compute_energy_conv)
    void compute_exchange_potential(const vecfuncT& nemo, vecfuncT& Knemo) const {

    	Exchange K = Exchange(world, this, 0);
    	Knemo = K(nemo);
    	truncate(world, Knemo);

    }

    // compute energy from given nemos and given OEP approximation for exchange
    // for example Slater potential for OAEP
    double compute_energy(const vecfuncT phi, const vecfuncT Jphi, const real_function_3d Vx,
    		const vecfuncT Kphi, const bool vir) const {

    	// compute kinetic energy
    	// it's ok to use phi here, no regularization necessary for this eigenvalue
    	double E_kin = 0.0;
    	for (int axis = 0; axis < 3; ++axis) {
    		real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    		const vecfuncT dphi = apply(world, D, phi);
    		E_kin += 0.5 * (inner(world, dphi, dphi)).sum();
    		// -1/2 sum <Psi|Nabla^2|Psi> = 1/2 sum <NablaPsi|NablaPsi>   (integration by parts)
    	}
    	E_kin *= 2.0; // because a = b closed shell

    	// check vir: if true -> compute Evir using Levy-Perdew virial relation (Kohut_2014, (43))
    	// if false -> compute exchange energy using the expectation value of the exchange operator
    	double E_X;
    	if (vir) {
        	// make vector of functions r = (x, y, z)
        	auto monomial_x = [] (const coord_3d& r) {return r[0];};
        	auto monomial_y = [] (const coord_3d& r) {return r[1];};
        	auto monomial_z = [] (const coord_3d& r) {return r[2];};
        	vecfuncT r(3);
        	r[0]=real_factory_3d(world).functor(monomial_x);
        	r[1]=real_factory_3d(world).functor(monomial_y);
        	r[2]=real_factory_3d(world).functor(monomial_z);
        	// compute density rho = dot(R*nemo, R*nemo), note that phi should be R*nemo
        	real_function_3d rho = 2.0*dot(world, phi, phi); // because a = b closed shell
        	real_function_3d rhoterm = 3*rho + dot(world, r, grad(rho));
        	E_X = inner(Vx, rhoterm);
    	}
    	else {
    		// this uses the oep->get_calc()->amo orbitals in the current form
    		E_X = -1.0*inner(world, phi, Kphi).sum();
    	}

    	// compute external potential (nuclear attraction)
    	real_function_3d Vext = calc->potentialmanager->vnuclear();
    	const vecfuncT Vextphi = Vext*phi;

    	// compute remaining energies: nuclear attraction, Coulomb, nuclear repulsion
    	// computed as ecpectation values (see Szabo, Ostlund (3.81))
    	const double E_ext = 2.0 * inner(world, phi, Vextphi).sum(); // because a = b closed shell
    	const double E_J = inner(world, phi, Jphi).sum();
    	const double E_nuc = calc->molecule.nuclear_repulsion_energy();
    	double energy = E_kin + E_ext + E_J + E_X + E_nuc;

    	if (world.rank() == 0) {
    		printf("\n                       kinetic energy %15.8f\n", E_kin);
    		printf("   electron-nuclear attraction energy %15.8f\n", E_ext);
    		printf("                       Coulomb energy %15.8f\n", E_J);
    		if (vir) printf("  exchange energy (OEP approximation) %15.8f\n", E_X);
    		else printf("  exchange energy (exchange operator) %15.8f\n", E_X);
    		printf("     nuclear-nuclear repulsion energy %15.8f\n", E_nuc);
    		printf("                         total energy %15.8f\n\n", energy);
            // printf("    works for exact exchange functional only...\n");
    	}

    	return energy;

    }


    // check if this is still useful (for DCEP)
    /// compute the kinetic energy potential using the Kohut trick Eq. (30)
    real_function_3d kinetic_energy_potential2(const vecfuncT& nemo) const {

        // the density
        real_function_3d rhonemo=dot(world,nemo,nemo);
        real_function_3d rho=R_square*rhonemo;

        // compute tauL: Eq. (20) of Kohut
        vecfuncT Rnemo=R*nemo;
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_Rnemo=Laplace(Rnemo);
        real_function_3d tauL=-0.5*dot(world,laplace_Rnemo,Rnemo);
        save(tauL,"tauL");
        real_function_3d tauL_div_rho=binary_op(tauL,rho,dens_inv());
        save(tauL_div_rho,"tauL_div_rho");

        // compute tau = | grad(mo) |^2
        //  = dR.dR * |nemo|^2 + 2*dR.grad(nemo) + R^2*|grad(nemo)|^2
        real_function_3d tau=real_factory_3d(world).compressed();
        vecfuncT dR=this->nuclear_correlation->U1vec();

        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

        for (const real_function_3d& n : nemo) {
            vecfuncT gradnemo=grad(n);
            tau+=U1dot*n*n;
            tau-=2.0*n*dot(world,dR,gradnemo);  // grad(R) = -U1
            tau+=dot(world,gradnemo,gradnemo);
        }
        tau=0.5*tau*R_square;
        real_function_3d tau_div_rho=binary_op(tau,rho,dens_inv());
        save(tau_div_rho,"tau_div_rho");

        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rho,logme());
        save(logdensa,"logdensa");
        vecfuncT gradzeta=grad(logdensa);
        madness::save_function(gradzeta,"gradzeta");
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        save(d2rho_div_rho,"d2rho_div_rho");
        real_function_3d result=tau_div_rho-0.25*d2rho_div_rho;
        save(result,"kinetic2");
        return result;
    }

    // check if this is still useful (for DCEP)
    real_function_3d kinetic_energy_potential(const vecfuncT& nemo) const {

        const Nuclear U_op(world,this->nuclear_correlation);
        const Nuclear V_op(world,this->get_calc().get());

        const vecfuncT Vnemo=V_op(nemo);  // eprec is important here!
        const vecfuncT Unemo=U_op(nemo);

        // nabla^2 nemo
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_nemo=Laplace(nemo);

        vecfuncT tmp=Unemo-this->nuclear_correlation->U2()*nemo;
        laplace_nemo-=2.0*tmp;
        vecfuncT D2Rnemo=R*laplace_nemo;

//        // double check result: recompute the density from its laplacian
//        vecfuncT nemo_rec=apply(world,*poisson,D2Rnemo);
//        scale(world,nemo_rec,-1./(4.*constants::pi));
        vecfuncT Rnemo=mul(world,R,nemo);
//        vecfuncT diff=sub(world,Rnemo,nemo_rec);
//        double dnorm=norm2(world,diff);
//        print("dnorm of laplacian phi ",dnorm);

        // compute \sum_i \phi_i \Delta \phi_i
        real_function_3d phiD2phi=dot(world,Rnemo,D2Rnemo);
        save(phiD2phi,"phiD2phi");

        // compute \sum_i \phi_i \epsilon_i \phi_i
        vecfuncT R2nemo=R_square*nemo;
        const real_function_3d rho=2.0*dot(world,nemo,R2nemo);
        const real_function_3d rhonemo=2.0*dot(world,nemo,nemo);

        // compute T1
        double T1 = 0.0;
        for (int axis = 0; axis < 3; axis++) {
            real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
            const vecfuncT dnemo = apply(world, D, Rnemo);
            T1 += 2.0* 0.5 * (inner(world, dnemo, dnemo)).sum();
        }
        printf("T1 %16.8f \n",T1);

        std::vector<double> eps(nemo.size());
        for (std::size_t i=0; i<eps.size(); ++i) eps[i]=calc->aeps(i);
        scale(world,R2nemo,eps);
        real_function_3d phiepsilonphi=dot(world,R2nemo,nemo);

        // divide by the density
        real_function_3d numerator=2.0*(-0.5*phiD2phi);   // fac 2 for sum over spin orbitals
//        real_function_3d numerator=2.0*(-0.5*phiD2phi-phiepsilonphi);   // fac 2 for sum over spin orbitals
        real_function_3d kinetic1=binary_op(numerator,rho,dens_inv());
        save(kinetic1,"kinetic1");
        real_function_3d nu_bar=kinetic1 + 2.0*(-0.5)*binary_op(phiepsilonphi,rho,dens_inv());

        // reintroduce the nuclear potential *after* smoothing
        real_function_3d uvnuc=calc->potentialmanager->vnuclear()-nuclear_correlation->U2();
        nu_bar=nu_bar-uvnuc;

        // compute T2
        vecfuncT dipole(3), drho(3);
        dipole[0]=real_factory_3d(world).functor(MomentFunctor(1,0,0));
        dipole[1]=real_factory_3d(world).functor(MomentFunctor(0,1,0));
        dipole[2]=real_factory_3d(world).functor(MomentFunctor(0,0,1));
        drho[0]=make_ddensity(rhonemo,0);
        drho[1]=make_ddensity(rhonemo,1);
        drho[2]=make_ddensity(rhonemo,2);
        real_function_3d one=real_factory_3d(world).functor(MomentFunctor(0,0,0));
        real_function_3d arg=3.0*rho+dot(world,dipole,drho);
        double T2=0.5*inner(nu_bar,arg);
        printf("T2 %16.8f \n",T2);

        Coulomb J(world,this);
        XCOperator xc(world,this);

        real_function_3d vne=calc->potentialmanager->vnuclear();
        real_function_3d vj=J.potential();
        real_function_3d vxc=xc.make_xc_potential();

        real_function_3d euler=nu_bar + vne + vj + vxc;
        double eulernorm=euler.norm2();
        print("eulernorm ",eulernorm);
        save(euler,"euler");
        save(vne,"vne");
        save(vj,"vj");
        save(vxc,"vxc");

        return nu_bar;
    }

    // check if this is still useful
    /// compute the laplacian of the density using the log trick
    real_function_3d make_laplacian_density_oep(const real_function_3d& rhonemo) const {


        // U1^2 operator
        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        real_function_3d result=(2.0*U1dot).truncate();

        // U2 operator
        const real_function_3d V=calc->potentialmanager->vnuclear();
        const real_function_3d U2=nuclear_correlation->U2();

        real_function_3d term2=2.0*(U2-V).truncate();
        result-=term2;


        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rhonemo,logme());
        vecfuncT gradzeta=grad(logdensa);

        // zeta contribution
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        result+=d2rho_div_rho;

        real_function_3d term3=4.0*dot(world,nuclear_correlation->U1vec(),gradzeta);
        result-=term3;

        result=(R_square*rhonemo*result).truncate();
        save(result,"d2rho");

        // double check result: recompute the density from its laplacian
        real_function_3d rho_rec=-1./(4.*constants::pi)*(*poisson)(result);
        save(rho_rec,"rho_reconstructed");

        real_function_3d rho=rhonemo*R_square;
        save(rho,"rho");
        real_function_3d laplace_rho=div(grad(rho));
        save(laplace_rho,"d2rho_direct");


        return result;
    }

};


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world, argc, argv);
    std::cout.precision(6);

    const std::string input = "input";
    std::shared_ptr<SCF> calc(new SCF(world, input.c_str()));
    if (world.rank()==0) {
        calc->molecule.print();
        print("\n");
        calc->param.print(world);
    }

    std::shared_ptr<OEP> oep(new OEP(world, calc));
    const double energy = oep->value();

    if (world.rank()==0) {
        printf("final energy   %12.8f\n", energy);
        printf("finished at time %.1f\n", wall_time());
    }

    // save converged HF MOs and orbital energies
    vecfuncT HF_MOs = copy(world, oep->get_calc()->amo);
    tensorT HF_orbens = copy(oep->get_calc()->aeps);

//    oep->kinetic_energy_potential(oep->get_calc()->amo);
//    oep->kinetic_energy_potential2(oep->get_calc()->amo);
//    real_function_3d rhonemo = dot(world,oep->get_calc()->amo,oep->get_calc()->amo);
//    oep->make_laplacian_density_oep(rhonemo);


    // OAEP final energy
    printf("   +++ starting OCEP iterative calculation +++\n");
    oep->solve_ocep(HF_MOs, HF_orbens);






    finalize();
    return 0;
}
