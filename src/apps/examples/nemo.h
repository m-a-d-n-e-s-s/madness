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

 $Id$
 */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
 \file examples/nemo.h
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#ifndef NEMO_H_
#define NEMO_H_

#include <mra/mra.h>
#include <mra/funcplot.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>
#include <moldft/moldft.h>
#include <examples/correlationfactor.h>
#include <examples/nonlinsol.h>
#include <mra/vmra.h>

namespace madness {

// this class needs to be moved to vmra.h !!

// This class is used to store information for the non-linear solver
template<typename T, std::size_t NDIM>
class vecfunc {
public:
	World& world;
	std::vector<Function<T, NDIM> > x;

	vecfunc(World& world, const std::vector<Function<T, NDIM> >& x1) :
			world(world), x(x1) {
	}

	vecfunc(const std::vector<Function<T, NDIM> >& x1) :
			world(x1[0].world()), x(x1) {
	}

	vecfunc(const vecfunc& other) :
			world(other.world), x(other.x) {
	}

	vecfunc& operator=(const vecfunc& other) {
		x = other.x;
		return *this;
	}

	vecfunc operator-(const vecfunc& b) {
		return vecfunc(world, sub(world, x, b.x));
	}

	vecfunc operator+=(const vecfunc& b) { // Operator+= necessary
		x = add(world, x, b.x);
		return *this;
	}

	vecfunc operator*(double a) { // Scale by a constant necessary
		scale(world, x, a);
		return *this;
	}
};

/// the non-linear solver requires an inner product
template<typename T, std::size_t NDIM>
T inner(const vecfunc<T, NDIM>& a, const vecfunc<T, NDIM>& b) {
	Tensor<T> i = inner(a.world, a.x, b.x);
	return i.sum();
}

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
template<typename T, std::size_t NDIM>
struct allocator {
	World& world;
	const int n;

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(World& world, const int nn) :
			world(world), n(nn) {
	}

	/// allocate a vector of n empty functions
	vecfunc<T, NDIM> operator()() {
		return vecfunc<T, NDIM>(world, zero_functions<T, NDIM>(world, n));
	}
};

/// The Nemo class
class Nemo: public OptimizationTargetInterface {
	typedef std::shared_ptr<real_convolution_3d> poperatorT;

public:

	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	calc	the Calculation
	Nemo(World& world1, std::shared_ptr<Calculation> calc) :
			world(world1), calc(calc), coords_sum(-1.0) {
		initialize();
	}

	double value() {
		return value(calc->molecule.get_all_coords());
	}

	double value(const Tensor<double>& x) {

		// fast return if the reference is already solved at this geometry
		double xsq = x.sumsq();
		if (xsq == coords_sum)
			return calc->current_energy;

		calc->molecule.set_all_coords(x.reshape(calc->molecule.natom(), 3));
		coords_sum = xsq;

		// Make the nuclear potential, initial orbitals, etc.
		calc->make_nuclear_potential(world);
		calc->potentialmanager->vnuclear().print_size("vnuc");
		calc->project_ao_basis(world);
		save_function(calc->potentialmanager->vnuclear(),"vnuc");

		// read converged wave function from disk if there is one
		if (calc->param.no_compute) {
			calc->load_mos(world);
			return calc->current_energy;
		}

		if (calc->param.restart) {
			calc->load_mos(world);
		} else {
			calc->initial_guess(world);
			calc->param.restart = true;
		}

		double energy = solve();
		if (calc->param.save) calc->save_mos(world);

		// save the converged orbitals and nemos
		vecfuncT psi = mul(world, R, calc->amo);
		for (std::size_t imo = 0; imo < calc->amo.size(); ++imo) {
			save_function(calc->amo[imo], "nemo" + stringify(imo));
			save_function(psi[imo], "psi" + stringify(imo));
		}

		return energy;
	}

	Tensor<double> gradient(const Tensor<double>& x) {

		MADNESS_EXCEPTION("no nemo gradients", 1);
//        value(x); // Ensures DFT equations are solved at this geometry
//        return calc->derivatives(world);
	}

	std::shared_ptr<Calculation> get_calc() const {return calc;}


private:

	/// the world
	World& world;

	std::shared_ptr<Calculation> calc;

public:

	/// the nuclear correlation factor
	std::shared_ptr<NuclearCorrelationFactor> nuclear_correlation;

	/// the nuclear correlation factor
	real_function_3d R;

	/// the inverse nuclear correlation factor
	real_function_3d R_inverse;

private:

	/// sum of square of coords at last solved geometry
	mutable double coords_sum;

	/// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson;

	void initialize() {
		// construct the nuclear correlation factor:
		// use a Gauss-Slater factor, or the identity
		if (calc->param.nuclear_corrfac == "GaussSlater") {
			nuclear_correlation = std::shared_ptr<NuclearCorrelationFactor>(
					new GaussSlater(world, *calc));
		} else if (calc->param.nuclear_corrfac == "none") {
			nuclear_correlation = std::shared_ptr<NuclearCorrelationFactor>(
					new PseudoNuclearCorrelationFactor(world, *calc));
		} else if (calc->param.nuclear_corrfac == "two") {
			nuclear_correlation = std::shared_ptr<NuclearCorrelationFactor>(
					new TwoNuclearCorrelationFactor(world, *calc));
		} else {
			if (world.rank()==0) print(calc->param.nuclear_corrfac);
			MADNESS_EXCEPTION("unknown nuclear correlation factor", 1);
		}

		// the nuclear correlation function
		R = nuclear_correlation->function();
		R_inverse = nuclear_correlation->inverse();

		// construct the Poisson solver
		poisson = std::shared_ptr<real_convolution_3d>(
				CoulombOperatorPtr(world, calc->param.lo, calc->param.econv));

		// print some output for the user
		if (world.rank() == 0) calc->molecule.print();
	}

	/// solve the HF equations
	double solve() {

		// guess: multiply the guess orbitals with the inverse R
		calc->amo = mul(world, R_inverse, calc->amo);

		vecfuncT& nemo = calc->amo;
		long nmo = nemo.size();

		// NOTE that nemos are somewhat sensitive to sparse operations (why??)
		// Therefore set all tolerance thresholds to zero, also in the mul_sparse
		const double trantol = 0.0;

		normalize(nemo);

		// apply all potentials (J, K, Vnuc) on the nemos
		vecfuncT psi, Jnemo, Knemo, Vnemo, JKVpsi, Unemo;
		tensorT fock;

		double energy = 0.0;
		bool converged = false;

		typedef allocator<double, 3> allocT;
		typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
		allocT alloc(world, nemo.size());
		solverT solver(allocT(world, nemo.size()));

		// iterate the residual equations
		for (int iter = 0; iter < calc->param.maxiter; ++iter) {

			// compute potentials the Fock matrix: J - K + Vnuc
			compute_nemo_potentials(nemo, psi, Jnemo, Knemo, Vnemo, Unemo);

			// compute the fock matrix
			double ekinetic = 0.0;
			JKVpsi = mul(world, R, add(world, sub(world, Jnemo, Knemo), Vnemo));
			fock = calc->make_fock_matrix(world, psi, JKVpsi, calc->aocc,
					ekinetic);

			// report the off-diagonal fock matrix elements
			tensorT fock_offdiag=copy(fock);
			for (int i=0; i<fock.dim(0); ++i) fock_offdiag(i,i)=0.0;
			double max_fock_offidag=fock_offdiag.absmax();
			if (world.rank()==0) print("F max off-diagonal  ",max_fock_offidag);
//			print(fock);

			double oldenergy=energy;
			energy = compute_energy(psi, mul(world, R, Jnemo),
					mul(world, R, Knemo));

			// Diagonalize overlap to get the eigenvalues and eigenvectors
			tensorT overlap = matrix_inner(world, psi, psi, true);
			const tensorT U=calc->get_fock_transformation(world,overlap,
		    		fock,calc->aeps,calc->aocc,FunctionDefaults<3>::get_thresh());

			nemo = transform(world, nemo, U, trantol, true);
			rotate_subspace(world, U, solver, 0, nemo.size(), trantol);

			truncate(world, nemo);
			normalize(nemo);

			// update the nemos

			// construct the BSH operator
			tensorT eps(nmo);
			for (int i = 0; i < nmo; ++i) eps(i) = std::min(-0.05, fock(i, i));
			std::vector<poperatorT> ops = calc->make_bsh_operators(world, eps);

			// make the potential * nemos term; make sure it's in phase with nemo
			START_TIMER(world);
			vecfuncT Vpsi = add(world, sub(world, Jnemo, Knemo), Unemo);
			Vpsi = transform(world, Vpsi, U, trantol, true);
			END_TIMER(world, "make Vpsi");

			// apply the BSH operator on the wave function
			START_TIMER(world);
			scale(world, Vpsi, -2.0);
			vecfuncT tmp = apply(world, ops, Vpsi);
			truncate(world, tmp);
			END_TIMER(world, "apply BSH");

			// compute the residuals
			vecfuncT residual = sub(world, nemo, tmp);
			const double norm = norm2(world, residual) / sqrt(nemo.size());

			// kain works best in the quadratic region
			vecfuncT nemo_new;
			if (norm < 5.e-1) {
				nemo_new = (solver.update(nemo, residual)).x;
			} else {
				nemo_new = tmp;
			}
			normalize(nemo_new);

			calc->do_step_restriction(world,nemo,nemo_new,"ab spin case");
			orthonormalize(nemo_new);
			nemo=nemo_new;

			if ((norm < calc->param.dconv) and
					(fabs(energy-oldenergy)<calc->param.econv))
				converged = true;

			if (world.rank() == 0) {
				printf(
					"finished iteration %2d at time %8.1fs with energy %12.8f\n",
						iter, wall_time(), energy);
				print("current residual norm  ", norm, "\n");
			}

			if (converged)
				break;
		}

		if (converged) {
			print("\nIterations converged\n");
		} else {
			print("\nIterations failed\n");
			energy = 0.0;
		}

		return energy;
	}

	/// given nemos, compute the HF energy
	double compute_energy(const vecfuncT& psi, const vecfuncT& Jpsi,
			const vecfuncT& Kpsi) const {

		const vecfuncT Vpsi = mul(world, calc->potentialmanager->vnuclear(),
				psi);
		const tensorT V = inner(world, Vpsi, psi);
		const double pe = 2.0 * V.sum();  // closed shell

		double ke = 0.0;
		for (int axis = 0; axis < 3; axis++) {
			real_derivative_3d D = free_space_derivative<double, 3>(world,
					axis);
			const vecfuncT dpsi = apply(world, D, psi);
			ke += 0.5 * (inner(world, dpsi, dpsi)).sum();
		}
		ke *= 2.0; // closed shell

		const double J = inner(world, psi, Jpsi).sum();
		const double K = inner(world, psi, Kpsi).sum();
		const double nucrep = calc->molecule.nuclear_repulsion_energy();

		const double energy = ke + J - K + pe + nucrep;
		if (world.rank() == 0) {
			printf("\n              kinetic %16.8f\n", ke);
			printf("   nuclear attraction %16.8f\n", pe);
			printf("              coulomb %16.8f\n", J);
			printf(" exchange-correlation %16.8f\n", -K);
			printf("    nuclear-repulsion %16.8f\n", nucrep);
			printf("                total %16.8f\n\n", energy);
		}
		return energy;
	}

	/// compute the reconstructed orbitals, and all potentials applied on nemo

	/// to use these potentials in the fock matrix computation they must
	/// be multiplied by the nuclear correlation factor
	/// @param[in]	nemo	the nemo orbitals
	/// @param[out]	psi		the reconstructed, full orbitals
	/// @param[out]	Jnemo	Coulomb operator applied on the nemos
	/// @param[out]	Knemo	exchange operator applied on the nemos
	/// @param[out]	Vnemo	nuclear potential applied on the nemos
	/// @param[out]	Unemo	regularized nuclear potential applied on the nemos
	void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& psi,
			vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& Vnemo,
			vecfuncT& Unemo) const {

		// reconstruct the orbitals
		START_TIMER(world);
		psi = mul(world, R, nemo);
		truncate(world, psi);
		END_TIMER(world, "reconstruct psi");

		// compute the density and the coulomb potential
		START_TIMER(world);
		real_function_3d J = get_coulomb_potential(psi);
		Jnemo = mul(world, J, nemo);
		truncate(world, Jnemo);
		END_TIMER(world, "compute Jnemo");

		// compute the exchange potential
		START_TIMER(world);
		Knemo = apply_exchange(nemo, psi);
		truncate(world, Knemo);
		END_TIMER(world, "compute Knemo");

		START_TIMER(world);
		const real_function_3d& Vnuc = calc->potentialmanager->vnuclear();
		Vnemo = mul(world, Vnuc, nemo);
		truncate(world, Vnemo);
		END_TIMER(world, "compute Vnemo");

		START_TIMER(world);
		Unemo.clear();
		for (std::size_t i = 0; i < nemo.size(); ++i) {
			Unemo.push_back(nuclear_correlation->apply_U(nemo[i]));
		}
		END_TIMER(world, "compute Unemo");

	}

	/// normalize the nemos
	void normalize(vecfuncT& nemo) const {

		// compute the norm of the reconstructed orbitals, includes the factor
		vecfuncT mos = mul(world, R, nemo);
		std::vector<double> norms = norm2s(world, mos);

		// scale the nemos, excludes the nuclear correlation factor
		std::vector<double> invnorm(norms.size());
		for (std::size_t i = 0; i < norms.size(); ++i)
			invnorm[i] = 1.0 / norms[i];
		scale(world, nemo, invnorm);
		truncate(world, nemo);
	}

    /// orthonormalize the vectors

    /// @param[inout]	amo_new	the vectors to be orthonormalized
    void orthonormalize(vecfuncT& nemo) const {
        TAU_START("Orthonormalize");
        START_TIMER(world);
		vecfuncT mos = mul(world, R, nemo);
        double trantol = 0.0;
        madness::normalize(world, mos);
        nemo = transform(world, nemo, Q3(matrix_inner(world, mos, mos)), trantol, true);
        truncate(world, nemo);
        normalize(nemo);
        END_TIMER(world, "Orthonormalize");
        TAU_STOP("Orthonormalize");
    }

	/// return the Coulomb potential
	real_function_3d get_coulomb_potential(const vecfuncT& psi) const {
		MADNESS_ASSERT(calc->param.spin_restricted);
		functionT rho = calc->make_density(world, calc->aocc, psi).scale(2.0);
		return calc->make_coulomb_potential(rho);
	}

	vecfuncT apply_exchange(const vecfuncT& nemo, const vecfuncT& psi) const {

		// IMPORTANT NOTE:
		// The mul_sparse in apply_hf_exchange uses a tolerance that is
		// too loose. Fails even for H2O, eprec=1.e-5
//    	return calc->apply_hf_exchange(world,calc->aocc,psi,nemo);
		vecfuncT result = zero_functions<double, 3>(world, int(nemo.size()));
		compress(world, result);
		for (std::size_t i = 0; i < nemo.size(); ++i) {
			for (std::size_t k = 0; k < psi.size(); ++k) {
				const real_function_3d ik = psi[i] * psi[k];
				result[i] += nemo[k] * (*poisson)(ik);
			}
		}
		return result;
	}

	template<typename solverT>
	void rotate_subspace(World& world, const tensorT& U, solverT& solver,
			int lo, int nfunc, double trantol) const {
		std::vector < vecfunc<double, 3> > &ulist = solver.get_ulist();
		std::vector < vecfunc<double, 3> > &rlist = solver.get_rlist();
		for (unsigned int iter = 0; iter < ulist.size(); ++iter) {
			vecfuncT& v = ulist[iter].x;
			vecfuncT& r = rlist[iter].x;
			transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), U, trantol,
					false);
			transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), U, trantol,
					true);
		}
	}

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const Function<T,NDIM>& f, const std::string name) const {
        if (world.rank()==0) print("saving function",name);
        f.print_size(name);
        archive::ParallelOutputArchive ar(world, name.c_str(), 1);
        ar & f;
    }


};

}

#endif /* NEMO_H_ */

