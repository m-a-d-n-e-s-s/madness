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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
 \file examples/nemo.cc
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#include <chem/nemo.h>
#include <chem/molecular_optimizer.h>
namespace madness {

static double ttt, sss;
void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss;
    if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}


extern Tensor<double> Q3(const Tensor<double>& s);

double Nemo::value(const Tensor<double>& x) {

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

    // construct the Poisson solver
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, calc->param.lo, calc->param.econv));

    // construct the nuclear correlation factor:
    nuclear_correlation=create_nuclear_correlation_factor(world,*calc);
    R = nuclear_correlation->function();
    R_inverse = nuclear_correlation->inverse();

	print_nuclear_corrfac();

	// read converged wave function from disk if there is one
	if (calc->param.no_compute) {
		calc->load_mos(world);
		return calc->current_energy;
	}

	if (calc->param.restart) {
		calc->load_mos(world);
	} else {
		calc->initial_guess(world);

		// guess: multiply the guess orbitals with the inverse R
		calc->amo = mul(world, R_inverse, calc->amo);
//		calc->param.restart = true;
	}

	double energy = solve();

	calc->current_energy=energy;

	// localize the orbitals
	if (calc->param.localize) {
		calc->amo=localize(calc->amo);
		tensorT fock=compute_fock_matrix(calc->amo,calc->aocc);
		if (world.rank()==0) print("localized Fock matrix \n",fock);
		for (std::size_t i=0; i<calc->amo.size(); ++i) {
			calc->aeps(i)=fock(i,i);
			if (world.rank()==0) print("orbital energy ",i,calc->aeps(i));
		}
	}

	if (calc->param.save) calc->save_mos(world);

	// save the converged orbitals and nemos
	vecfuncT psi = mul(world, R, calc->amo);
	truncate(world,psi);
	for (std::size_t imo = 0; imo < calc->amo.size(); ++imo) {
		save_function(calc->amo[imo], "nemo" + stringify(imo));
		save_function(psi[imo], "psi" + stringify(imo));
	}

	// compute the dipole moment
    functionT rho = calc->make_density(world, calc->aocc, psi).scale(2.0);
    calc->dipole(world,rho);

    // compute stuff
    functionT rhonemo = calc->make_density(world, calc->aocc, calc->amo).scale(2.0);
    real_derivative_3d Dz = free_space_derivative<double, 3>(world,2);
    real_function_3d rhonemoz=Dz(rhonemo);
    std::string filename="plot_rhonemoz";
    Vector<double,3> lo=vec<double>(0,0,-10);
    Vector<double,3> hi=vec<double>(0,0,10);
    plot_line(filename.c_str(),500, lo, hi, rhonemoz);

    filename="plot_rhonemo";
    plot_line(filename.c_str(),500, lo, hi, rhonemo);

	return energy;
}

void Nemo::print_nuclear_corrfac() const {

	// the nuclear correlation function
	save_function(R,"R");
	save_function(nuclear_correlation->U2(),"U2");
	save_function(nuclear_correlation->U1(0),"U1x");
	save_function(nuclear_correlation->U1(1),"U1y");
	save_function(nuclear_correlation->U1(2),"U1z");

	// FIXME: plot_plane doesn't work for more that 1 rank
	if (world.size()==0) {
		plot_plane(world,R,"R");
		plot_plane(world,nuclear_correlation->U2(),"U2");
		plot_plane(world,nuclear_correlation->U1(0),"U1x");
	}
}

/// localize the nemo orbitals according to Pipek-Mezey
vecfuncT Nemo::localize(const vecfuncT& nemo) const {
	DistributedMatrix<double> dUT;
	const double tolloc = 1e-3;
	double trantol = calc->vtol / std::min(30.0, double(nemo.size()));

	std::vector<int> aset=calc->group_orbital_sets(world,calc->aeps,
			calc->aocc, nemo.size());
	// localize using the reconstructed orbitals
	vecfuncT psi = mul(world, R, nemo);
	dUT = calc->localize_PM(world, psi, aset, tolloc, 0.25, true, true);
	dUT.data().screen(trantol);

	tensorT UT(calc->amo.size(),calc->amo.size());
	dUT.copy_to_replicated(UT); // for debugging
	tensorT U = transpose(UT);

	vecfuncT localnemo = transform(world, nemo, U, true);
	truncate(world, localnemo);
	normalize(localnemo);
	return localnemo;
}

/// compute the Fock matrix from scratch
tensorT Nemo::compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const {
	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT psi, Jnemo, Knemo, Vnemo, JKVpsi, Unemo;

	// compute potentials the Fock matrix: J - K + Vnuc
	compute_nemo_potentials(nemo, psi, Jnemo, Knemo, Vnemo, Unemo);

	// compute the fock matrix
	double ekinetic = 0.0;
	JKVpsi = mul(world, R, add(world, sub(world, Jnemo, Knemo), Vnemo));
	tensorT fock = calc->make_fock_matrix(world, psi, JKVpsi, occ,
			ekinetic);
	return fock;
}

/// solve the HF equations
double Nemo::solve() {

	// guess has already been performed in value()
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
		if (world.rank()==0) {
			print("F max off-diagonal  ",max_fock_offidag);
			print(fock);
		}

		double oldenergy=energy;
		energy = compute_energy(psi, mul(world, R, Jnemo),
				mul(world, R, Knemo));
        energy = compute_energy_regularized(nemo, Jnemo, Knemo, Unemo);

		// Diagonalize overlap to get the eigenvalues and eigenvectors
		tensorT overlap = matrix_inner(world, psi, psi, true);
		const tensorT U=calc->get_fock_transformation(world,overlap,
	    		fock,calc->aeps,calc->aocc,FunctionDefaults<3>::get_thresh());

		START_TIMER(world);
		nemo = transform(world, nemo, U, trantol, true);
		rotate_subspace(world, U, solver, 0, nemo.size(), trantol);

		truncate(world, nemo);
		normalize(nemo);
		END_TIMER(world, "transform orbitals");

		// update the nemos

		// construct the BSH operator
		tensorT eps(nmo);
		for (int i = 0; i < nmo; ++i) {
			eps(i) = std::min(-0.05, fock(i, i));
		}
		if (calc->param.orbitalshift>0.0) {
			if (world.rank()==0) print("shifting orbitals by "
					,calc->param.orbitalshift," to lower energies");
			eps-=calc->param.orbitalshift;
		}
		std::vector<poperatorT> ops = calc->make_bsh_operators(world, eps);

		// make the potential * nemos term; make sure it's in phase with nemo
		START_TIMER(world);
		vecfuncT Vpsi = add(world, sub(world, Jnemo, Knemo), Unemo);
		truncate(world,Vpsi);
		END_TIMER(world, "make Vpsi");

		START_TIMER(world);
		Vpsi = transform(world, Vpsi, U, trantol, true);
		truncate(world,Vpsi);
		END_TIMER(world, "transform Vpsi");

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

		if (calc->param.save) calc->save_mos(world);

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
		if (world.rank()==0) print("\nIterations converged\n");
	} else {
		if (world.rank()==0) print("\nIterations failed\n");
		energy = 0.0;
	}

	return energy;
}


/// given nemos, compute the HF energy
double Nemo::compute_energy(const vecfuncT& psi, const vecfuncT& Jpsi,
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

/// given nemos, compute the HF energy using the regularized expressions for T and V
double Nemo::compute_energy_regularized(const vecfuncT& nemo, const vecfuncT& Jnemo,
        const vecfuncT& Knemo, const vecfuncT& Unemo) const {

    const vecfuncT R2nemo=mul(world,nuclear_correlation->square(),nemo);

    const tensorT U = inner(world, R2nemo, Unemo);
    const double pe = 2.0 * U.sum();  // closed shell

    double ke = 0.0;
    for (int axis = 0; axis < 3; axis++) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        const vecfuncT dnemo = apply(world, D, nemo);
        const vecfuncT dr2nemo = apply(world, D, R2nemo);
        ke += 0.5 * (inner(world, dnemo, dr2nemo)).sum();
    }
    ke *= 2.0; // closed shell

    const double J = inner(world, R2nemo, Jnemo).sum();
    const double K = inner(world, R2nemo, Knemo).sum();
    const double nucrep = calc->molecule.nuclear_repulsion_energy();

    const double energy = ke + J - K + pe + nucrep;
    if (world.rank() == 0) {
        printf("\n  nuclear and kinetic %16.8f\n", ke + pe);
        printf("              coulomb %16.8f\n", J);
        printf(" exchange-correlation %16.8f\n", -K);
        printf("    nuclear-repulsion %16.8f\n", nucrep);
        printf("   regularized energy %16.8f\n", energy);
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
void Nemo::compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& psi,
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
void Nemo::normalize(vecfuncT& nemo) const {

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
void Nemo::orthonormalize(vecfuncT& nemo) const {
    START_TIMER(world);
	vecfuncT mos = mul(world, R, nemo);
    double trantol = 0.0;
    madness::normalize(world, mos);
    nemo = transform(world, nemo, Q3(matrix_inner(world, mos, mos)), trantol, true);
    truncate(world, nemo);
    normalize(nemo);
    END_TIMER(world, "Orthonormalize");
}

/// return the Coulomb potential
real_function_3d Nemo::get_coulomb_potential(const vecfuncT& psi) const {
	MADNESS_ASSERT(calc->param.spin_restricted);
	functionT rho = calc->make_density(world, calc->aocc, psi).scale(2.0);
	return calc->make_coulomb_potential(rho);
}

vecfuncT Nemo::apply_exchange(const vecfuncT& nemo, const vecfuncT& psi) const {

	// IMPORTANT NOTE:
	// The mul_sparse in apply_hf_exchange uses a tolerance that is
	// too loose. Fails even for H2O, eprec=1.e-5
//    	return calc->apply_hf_exchange(world,calc->aocc,psi,nemo);
	vecfuncT result = zero_functions_compressed<double, 3>(world, int(nemo.size()));
	for (std::size_t i = 0; i < nemo.size(); ++i) {
		for (std::size_t k = 0; k < psi.size(); ++k) {
			const real_function_3d ik = psi[i] * psi[k];
			result[i] += nemo[k] * (*poisson)(ik);
		}
	}
	return result;
}

/// rotate the KAIN subspace (cf. SCF.cc)
template<typename solverT>
void Nemo::rotate_subspace(World& world, const tensorT& U, solverT& solver,
        int lo, int nfunc, double trantol) const {
    std::vector < vecfunc<double, 3> > &ulist = solver.get_ulist();
    std::vector < vecfunc<double, 3> > &rlist = solver.get_rlist();
    for (unsigned int iter = 0; iter < ulist.size(); ++iter) {
        vecfuncT& v = ulist[iter].x;
        vecfuncT& r = rlist[iter].x;
        vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), U,
                trantol, false);
        vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), U,
                trantol, true);

        world.gop.fence();
        for (int i=0; i<nfunc; i++) {
            v[i] = vnew[i];
            r[i] = rnew[i];
        }
    }
    world.gop.fence();
}

/// compute the nuclear gradients
Tensor<double> Nemo::gradient(const Tensor<double>& x) {
    START_TIMER(world);

#if 0
    const double thresh=FunctionDefaults<3>::get_thresh();
    FunctionDefaults<3>::set_thresh(thresh*0.1);
    R.set_thresh(thresh*0.1);
    vecfuncT psi = mul(world, R, calc->amo);
    functionT rho = calc->make_density(world, calc->aocc, psi).scale(2.0);
    R.set_thresh(thresh);
    Tensor<double> grad=calc->derivatives(world,rho);
    FunctionDefaults<3>::set_thresh(thresh);
#else

    Tensor<double> grad(3*calc->molecule.natom());
    functionT rhonemo = calc->make_density(world, calc->aocc, calc->amo).scale(2.0);

    vecfuncT bra(3);
    for (int axis=0; axis<3; ++axis) {

        // compute \frac{\partial \rho}{\partial x_i}
        real_derivative_3d D = free_space_derivative<double, 3>(world,axis);
        real_function_3d Drhonemo=D(rhonemo);

        // compute the second term of the bra
        real_function_3d tmp=rhonemo*nuclear_correlation->U1(axis);
        tmp.scale(2.0);
        bra[axis]=(Drhonemo-tmp).truncate();
    }

//    for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
//        const Atom& atom=calc->molecule.get_atom(iatom);
//        std::shared_ptr< FunctionFunctorInterface<double,3> > r2v(new
//        NuclearCorrelationFactor::square_times_V_functor(nuclear_correlation.get(),atom));
//
//        for (int axis=0; axis<3; axis++) {
//            grad(3*iatom + axis)=-bra[axis].inner_adaptive(r2v);
//        }
//    }

        for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
            for (int axis=0; axis<3; ++axis) {
                NuclearCorrelationFactor::square_times_V_derivative_functor r2v(
                        nuclear_correlation.get(),this->molecule(),iatom,axis);
                grad(3*iatom + axis)=inner(rhonemo,r2v);

            }
        }

    // add the nuclear contribution
    for (int atom = 0; atom < calc->molecule.natom(); ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            grad[atom * 3 + axis] +=
                    calc->molecule.nuclear_repulsion_derivative(atom,axis);
        }
    }

#endif
    END_TIMER(world, "compute gradients");

    if (world.rank() == 0) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print(
              "  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(
              " ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (int i = 0; i < calc->molecule.natom(); ++i) {
            const Atom& atom = calc->molecule.get_atom(i);
            printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", i,
                   atom.x, atom.y, atom.z, grad[i * 3 + 0], grad[i * 3 + 1],
                   grad[i * 3 + 2]);
        }
    }
    return grad;
}


/// save a function
template<typename T, size_t NDIM>
void Nemo::save_function(const Function<T,NDIM>& f, const std::string name) const {
    if (world.rank()==0) print("saving function",name);
    f.print_size(name);
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    ar & f;
}


} // namespace madness
