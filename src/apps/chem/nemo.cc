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
 \file examples/nemo.cc
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#include <chem/nemo.h>
#include <chem/projector.h>
#include <chem/molecular_optimizer.h>
#include <chem/SCFOperators.h>
#include <madness/constants.h>


namespace madness {

double Nemo::value(const Tensor<double>& x) {

    // fast return if the reference is already solved at this geometry
	double xsq = x.sumsq();
	if (xsq == coords_sum)
		return calc->current_energy;

	calc->molecule.set_all_coords(x.reshape(calc->molecule.natom(), 3));
	coords_sum = xsq;

	construct_nuclear_correlation_factor();

    // construct the Poisson solver
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, calc->param.lo, FunctionDefaults<3>::get_thresh()));

	print_nuclear_corrfac();

	double energy=0.0;

	// read converged wave function from disk if there is one
	if (calc->param.no_compute) {
		calc->load_mos(world);

	} else {

        // guess: read from file or multiply the guess orbitals with the inverse R
	    if (calc->param.restart) {
	        calc->load_mos(world);

	    } else {

            calc->initial_guess(world);
	        real_function_3d R_inverse = nuclear_correlation->inverse();
	        calc->amo = mul(world, R_inverse, calc->amo);
	    }

	    for (protocol p(*this); not p.finished(); ++p) {
	        set_protocol(p.current_prec);
	        energy=solve(p);
	    }
	    set_protocol(get_calc()->param.econv);

	    calc->current_energy=energy;
	    if (calc->param.save) calc->save_mos(world);

	    // save the converged orbitals and nemos
	    vecfuncT psi = mul(world, R, calc->amo);
	    truncate(world,psi);
	    for (std::size_t imo = 0; imo < calc->amo.size(); ++imo) {
	        save_function(calc->amo[imo], "nemo" + stringify(imo));
	        save_function(psi[imo], "psi" + stringify(imo));
	    }
	}

	// compute the dipole moment
	functionT rho = 2.0*(R_square*make_density(calc->aocc, calc->amo)).truncate();
	calc->dipole(world,rho);

	// compute the hessian
	if (calc->param.hessian) hessian(x);

	return energy;
}


void Nemo::print_nuclear_corrfac() const {

	// the nuclear correlation function
	save_function(R,"R");
	save_function(nuclear_correlation->U2(),"U2");
	save_function(nuclear_correlation->U1(0),"U1x");
	save_function(nuclear_correlation->U1(1),"U1y");
	save_function(nuclear_correlation->U1(2),"U1z");

}

/// localize the nemo orbitals according to Pipek-Mezey
vecfuncT Nemo::localize(const vecfuncT& nemo) const {
	DistributedMatrix<double> dUT;
	const double tolloc = 1e-3;

	std::vector<int> aset=calc->group_orbital_sets(world,calc->aeps,
			calc->aocc, nemo.size());
	// localize using the reconstructed orbitals
	vecfuncT psi = mul(world, R, nemo);
	dUT = calc->localize_PM(world, psi, aset, tolloc, 0.25, true, true);
	dUT.data().screen(trantol());

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
double Nemo::solve(const protocol& proto) {

	// guess has already been performed in value()
	vecfuncT& nemo = calc->amo;
	long nmo = nemo.size();

	// NOTE that nemos are somewhat sensitive to sparse operations (why??)
	// Therefore set all tolerance thresholds to zero, also in the mul_sparse


	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT psi, Jnemo, Knemo, Vnemo, JKVpsi, Unemo;
	tensorT fock;

	double energy = 0.0;
	bool converged = false;
	bool localized=calc->param.localize;

	typedef allocator<double, 3> allocT;
	typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
	allocT alloc(world, nemo.size());
	solverT solver(allocT(world, nemo.size()));

	// iterate the residual equations
	for (int iter = 0; iter < calc->param.maxiter; ++iter) {

	    if (localized) nemo=localize(nemo);

		// compute potentials the Fock matrix: J - K + Vnuc
		compute_nemo_potentials(nemo, psi, Jnemo, Knemo, Vnemo, Unemo);

		// compute the fock matrix
		double ekinetic = 0.0;
		JKVpsi = mul(world, R, add(world, sub(world, Jnemo, Knemo), Vnemo));
		fock = calc->make_fock_matrix(world, psi, JKVpsi, calc->aocc,
				ekinetic);
		JKVpsi.clear();

		// report the off-diagonal fock matrix elements
		if (not localized) {
            tensorT fock_offdiag=copy(fock);
            for (int i=0; i<fock.dim(0); ++i) fock_offdiag(i,i)=0.0;
            double max_fock_offidag=fock_offdiag.absmax();
            if (world.rank()==0) print("F max off-diagonal  ",max_fock_offidag);
		}

		double oldenergy=energy;
//		energy = compute_energy(psi, mul(world, R, Jnemo),
//				mul(world, R, Knemo));
        energy = compute_energy_regularized(nemo, Jnemo, Knemo, Unemo);

        // Diagonalize the Fock matrix to get the eigenvalues and eigenvectors
        tensorT U;
        if (not localized) {
            tensorT overlap = matrix_inner(world, psi, psi, true);
            U=calc->get_fock_transformation(world,overlap,
                    fock,calc->aeps,calc->aocc,FunctionDefaults<3>::get_thresh());

            START_TIMER(world);
            nemo = transform(world, nemo, U, trantol(), true);
            rotate_subspace(world, U, solver, 0, nemo.size());

            truncate(world, nemo);
            normalize(nemo);
            END_TIMER(world, "transform orbitals");
        }

		// update the nemos

		// construct the BSH operator and add off-diagonal elements
		// (in case of non-canonical HF)
		tensorT eps(nmo);
		for (int i = 0; i < nmo; ++i) {
			eps(i) = std::min(-0.05, fock(i, i));
            fock(i, i) -= eps(i);
		}
        vecfuncT fnemo;
        if (localized) fnemo= transform(world, nemo, fock, trantol(), true);

        // Undo the damage
        for (int i = 0; i < nmo; ++i) fock(i, i) += eps(i);

		if (calc->param.orbitalshift>0.0) {
			if (world.rank()==0) print("shifting orbitals by "
					,calc->param.orbitalshift," to lower energies");
			eps-=calc->param.orbitalshift;
		}
		std::vector<poperatorT> ops = calc->make_bsh_operators(world, eps);

		// make the potential * nemos term; make sure it's in phase with nemo
		START_TIMER(world);
		vecfuncT Vpsi = add(world, sub(world, Jnemo, Knemo), Unemo);
        if (localized) gaxpy(world, 1.0, Vpsi, -1.0, fnemo);
		truncate(world,Vpsi);
		END_TIMER(world, "make Vpsi");

		START_TIMER(world);
		if (not localized) Vpsi = transform(world, Vpsi, U, trantol(), true);
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

		if ((norm < proto.dconv) and
				(fabs(energy-oldenergy)<proto.econv))
			converged = true;

		if (calc->param.save) calc->save_mos(world);

		if (world.rank() == 0) {
			printf("finished iteration %2d at time %8.1fs with energy %12.8f\n",
					iter, wall_time(), energy);
			print("current residual norm  ", norm, "\n");
		}

		if (converged) break;
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

	int ispin=0;
	double exc=0.0;
	if (calc->xc.is_dft()) {
	    XCOperator xcoperator(world,this,ispin);
	    exc=xcoperator.compute_xc_energy();
	}

	const double nucrep = calc->molecule.nuclear_repulsion_energy();
	double energy = ke + J + pe + nucrep;
	if (is_dft()) energy+=exc;
	else energy-=K;

	if (world.rank() == 0) {
		printf("\n              kinetic %16.8f\n", ke);
		printf("   nuclear attraction %16.8f\n", pe);
		printf("              coulomb %16.8f\n", J);
        if (is_dft()) {
            printf(" exchange-correlation %16.8f\n", exc);
        } else {
            printf("             exchange %16.8f\n", -K);
        }
		printf("    nuclear-repulsion %16.8f\n", nucrep);
		printf("                total %16.8f\n\n", energy);
        printf("  buggy if hybrid functionals are used..\n");
	}
	return energy;
}

/// given nemos, compute the HF energy using the regularized expressions for T and V
double Nemo::compute_energy_regularized(const vecfuncT& nemo, const vecfuncT& Jnemo,
        const vecfuncT& Knemo, const vecfuncT& Unemo) const {

    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

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

    int ispin=0;
    double exc=0.0;
    if (calc->xc.is_dft()) {
        XCOperator xcoperator(world,this,ispin);
        exc=xcoperator.compute_xc_energy();
    }

    const double nucrep = calc->molecule.nuclear_repulsion_energy();

    double energy = ke + J + pe + nucrep;
    if (is_dft()) energy+=exc;
    else energy-=K;

    if (world.rank() == 0) {
        printf("\n  nuclear and kinetic %16.8f\n", ke + pe);
        printf("              coulomb %16.8f\n", J);
        if (is_dft()) {
            printf(" exchange-correlation %16.8f\n", exc);
        } else {
            printf("             exchange %16.8f\n", -K);
        }
        printf("    nuclear-repulsion %16.8f\n", nucrep);
        printf("   regularized energy %16.8f\n", energy);
        printf("  buggy if hybrid functionals are used..\n");
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
	Coulomb J=Coulomb(world,this);
	Jnemo = J(nemo);
	truncate(world, Jnemo);
	END_TIMER(world, "compute Jnemo");

	// compute the exchange potential
    int ispin=0;
    Knemo=zero_functions_compressed<double,3>(world,nemo.size());
    if (calc->xc.hf_exchange_coefficient()>0.0) {
        START_TIMER(world);
        Exchange K=Exchange(world,this,ispin).same(true).small_memory(false);
        Knemo=K(nemo);
        scale(world,Knemo,calc->xc.hf_exchange_coefficient());
        truncate(world, Knemo);
        END_TIMER(world, "compute Knemo");
    }

	// compute the exchange-correlation potential
    if (calc->xc.is_dft()) {
        START_TIMER(world);
        XCOperator xcoperator(world,this,ispin);
        double exc=0.0;
        if (ispin==0) exc=xcoperator.compute_xc_energy();
        print("exc",exc);
        Knemo=sub(world,Knemo,xcoperator(nemo));   // minus times minus gives plus
        truncate(world,Knemo);
        double size=get_size(world,Knemo);
        END_TIMER(world, "compute XCnemo "+stringify(size));
    }

	START_TIMER(world);
	const real_function_3d& Vnuc = calc->potentialmanager->vnuclear();
	Vnemo = mul(world, Vnuc, nemo);
	truncate(world, Vnemo);
    double size=get_size(world,Vnemo);
	END_TIMER(world, "compute Vnemo "+stringify(size));

	START_TIMER(world);
	Nuclear Unuc(world,this->nuclear_correlation);
	Unemo=Unuc(nemo);
    size=get_size(world,Unemo);
	END_TIMER(world, "compute Unemo "+stringify(size));

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
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    normalize(nemo);
    double maxq;
    do {
        vecfuncT R2nemo=mul(world,R_square,nemo);
        tensorT Q = Q2(matrix_inner(world, R2nemo, nemo));
        maxq=0.0;
        for (int i=0; i<Q.dim(0); ++i)
            for (int j=0; j<i; ++j)
                maxq = std::max(maxq,std::abs(Q(i,j)));

        Q.screen(trantol()); // ???? Is this really needed?
        nemo = transform(world, nemo, Q, trantol(), true);
        truncate(world, nemo);
        if (world.rank() == 0) print("ORTHOG2: maxq trantol", maxq, trantol());

    } while (maxq>0.01);
    normalize(nemo);
    END_TIMER(world, "Orthonormalize");
}

/// return the Coulomb potential
real_function_3d Nemo::get_coulomb_potential(const vecfuncT& psi) const {
	MADNESS_ASSERT(calc->param.spin_restricted);
	functionT rho = make_density(calc->aocc, psi).scale(2.0);
	return calc->make_coulomb_potential(rho);
}

real_function_3d Nemo::make_density(const Tensor<double>& occ,
        const vecfuncT& nemo) const {
    return calc->make_density(world,occ,nemo).truncate();
}

real_function_3d Nemo::make_density(const tensorT & occ,
        const vecfuncT& bra, const vecfuncT& ket) const {

    vecfuncT vsq = mul(world, bra, ket);
    compress(world, vsq);
    functionT rho = factoryT(world).compressed();
    for (unsigned int i = 0; i < vsq.size(); ++i) {
        if (calc->get_aocc()[i])
            rho.gaxpy(1.0, vsq[i], calc->get_aocc()[i], false);
    }
    world.gop.fence();
    return rho.truncate();
}

/// compute the nuclear gradients
Tensor<double> Nemo::gradient(const Tensor<double>& x) {
    START_TIMER(world);

    vecfuncT nemo=calc->amo;
    {
        Tensor<double> grad2(3*calc->molecule.natom());
        vecfuncT R2nemo=mul(world,R_square,nemo);
        for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
            for (int iaxis=0; iaxis<3; iaxis++) {
                DNuclear Dnuc(world,this,iatom,iaxis);
                vecfuncT dnucnemo=Dnuc(nemo);
                grad2(3*iatom + iaxis)=2.0*(inner(world,R2nemo,dnucnemo)).sum();
            }
        }
        print("grad2");
        print(grad2);
    }


    // the pseudo-density made up of the square of the nemo orbitals
    functionT rhonemo = make_density(calc->aocc, nemo).scale(2.0);

    // the following block computes the gradients more precisely than the
    // direct evaluation of the derivative of the nuclear potential
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

    Tensor<double> grad(3*calc->molecule.natom());

    for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
        const Atom& atom=calc->molecule.get_atom(iatom);
        NuclearCorrelationFactor::square_times_V_functor r2v(nuclear_correlation.get(),atom);

        for (int axis=0; axis<3; axis++) {
            grad(3*iatom + axis)=-inner(bra[axis],r2v);
        }
    }
    print("grad");
    print(grad);

//    // this block is less precise
//    for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
//        for (int axis=0; axis<3; ++axis) {
//            NuclearCorrelationFactor::square_times_V_derivative_functor r2v(
//                    nuclear_correlation.get(),this->molecule(),iatom,axis);
//            grad(3*iatom + axis)=inner(rhonemo,r2v);
//
//        }
//    }

    // add the nuclear contribution
    for (int atom = 0; atom < calc->molecule.natom(); ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            grad[atom * 3 + axis] +=
                    calc->molecule.nuclear_repulsion_derivative(atom,axis);
        }
    }

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


/// compute the nuclear hessian
Tensor<double> Nemo::hessian(const Tensor<double>& x) {
    START_TIMER(world);

    const int natom=molecule().natom();
    const vecfuncT& nemo=get_calc()->amo;
    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    Tensor<double> hessian(3*natom,3*natom);

    // the perturbed MOs determined by the CPHF equations
    std::vector<vecfuncT> xi=compute_all_cphf();


    // compute the derivative of the density d/dx rho; 2: closed shell
    const real_function_3d dens=2.0*make_density(get_calc()->get_aocc(),nemo,R2nemo);
    vecfuncT drho(3);
    for (int i=0; i<3; ++i) {
        real_derivative_3d D = free_space_derivative<double, 3>(world,i);
        drho[i]=D(dens).truncate();
    }


    // add the electronic contribution to the hessian
    int i=0;
    for (int iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            i=iatom*3 + iaxis;

            // compute the full perturbed density
            // \rho_pt = R2 F_i F_i^X + R^X R2 F_i F_i
            real_function_3d dens_pt=4.0*make_density(calc->get_aocc(),R2nemo,xi[i]);
            NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,-1);
            const real_function_3d RX_div_R=real_factory_3d(world).functor2(rxr_func).truncate_on_project();
            dens_pt=(dens_pt+2.0*RX_div_R*dens).truncate();
            save_function(dens_pt,"dens_pt"+stringify(i));

            int j=0;
            for (int jatom=0; jatom<natom; ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    j=jatom*3 + jaxis;

                    MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                    double result=inner(dens_pt,mdf);

                    // integration by parts
                    if (iatom==jatom) result+=inner(drho[iaxis],mdf);

                    // skip diagonal elements because they are extremely noisy!
                    // use translational symmetry to reconstruct them from other
                    // hessian matrix elements (see below)
                    if (i==j) result=0.0;
                    hessian(i,j)=result;
                    ++j;
                }
            }
            ++i;
        }
    }
    if (world.rank() == 0) {
        print("\n raw electronic Hessian (a.u.)\n");
        print(hessian);
    }

    Tensor<double> asymmetric=0.5*(hessian-transpose(hessian));
    const double max_asymmetric=asymmetric.absmax();
    if (world.rank() == 0) {
        print("\n asymmetry in the electronic Hessian (a.u.)\n");
        print(asymmetric);
        print("max asymmetric element: ",max_asymmetric);
    }


    // symmetrize hessian
    hessian+=transpose(hessian);
    hessian.scale(0.5);

    // exploit translational symmetry to compute the diagonal elements:
    // translating all atoms in the same direction will make no energy change,
    // therefore the respective sum of hessian matrix elements will be zero:
    for (int i=0; i<3*natom; ++i) {
        double sum=0.0;
        for (int j=0; j<3*natom; j+=3) sum+=hessian(i,j+(i%3));
        hessian(i,i)=-sum;
    }

    if (world.rank() == 0) {
        print("\n electronic Hessian (a.u.)\n");
        print(hessian);
    }

    // add the nuclear-nuclear contribution
    hessian+=molecule().nuclear_repulsion_hessian();

    if (world.rank() == 0) {
        print("\n Hessian (a.u.)\n");
        print(hessian);
    }
    END_TIMER(world, "compute hessian");

    Tensor<double> frequencies=compute_frequencies(hessian,false,true);

    if (world.rank() == 0) {
        print("\n vibrational frequencies (a.u.)\n");
        print(frequencies);
        print("\n vibrational frequencies (cm-1)\n");
        print(constants::au2invcm*frequencies);
    }

    frequencies=compute_frequencies(hessian,true,true);

    if (world.rank() == 0) {
        print("\n vibrational frequencies projected (a.u.)\n");
        print(frequencies);
        print("\n vibrational frequencies projected (cm-1)\n");
        print(constants::au2invcm*frequencies);
    }

    return hessian;

}

/// solve the CPHF equation for the nuclear displacements

/// @param[in]  iatom   the atom A to be moved
/// @param[in]  iaxis   the coordinate X of iatom to be moved
/// @return     \frac{\partial}{\partial X_A} \varphi
vecfuncT Nemo::cphf(const int iatom, const int iaxis, const Tensor<double> fock,
        const vecfuncT& guess, const protocol& proto) const {

    print("\nsolving nemo cphf equations for atom, axis",iatom,iaxis);

    // guess for the perturbed MOs
    const vecfuncT nemo=calc->amo;
    vecfuncT xi=copy(world,calc->amo);
    scale(world,xi,0.5);
    if (guess.size()>0) {
        if (world.rank()==0) print("using guess for the CPHF vectors");
        xi=copy(world,guess);
    }
    const int nmo=nemo.size();

    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);
    QProjector<double,3> Q(world,R2nemo,nemo);

    // construct some intermediates
    const Tensor<double> occ=get_calc()->get_aocc();
    const real_function_3d rhonemo=2.0*make_density(occ,nemo); // closed shell
    NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor2(rxr_func).truncate_on_project();

    // construct quantities that are independent of xi

    // derivative of the (regularized) nuclear potential
    DNuclear Dunuc(world,this,iatom,iaxis);
    vecfuncT Vpsi2b=Dunuc(nemo);
    truncate(world,Vpsi2b);

    // part of the Coulomb operator with the derivative of the NCF
    // J <- \int dr' 1/|r-r'| \sum_i R^XR F_iF_i
    Coulomb Jconst(world);
    Jconst.potential()=Jconst.compute_potential(2.0*RXR*rhonemo);        // factor 2 for cphf
    vecfuncT Jconstnemo=Jconst(nemo);
    truncate(world,Jconstnemo);

    // part of the exchange operator with the derivative of the NCF
    // K <- \sum_k |F_k> \int dr' 1/|r-r'| 2R^XR F_k F_i
    // there is no constant term for DFT, since the potentials are not
    // linear in the density
    vecfuncT Kconstnemo=zero_functions_compressed<double,3>(world,nmo);
    if (not is_dft()) {
        Exchange Kconst(world);
        vecfuncT kbra=mul(world,RXR,nemo);
        scale(world,kbra,2.0);
        truncate(world,kbra);
        Kconst.set_parameters(kbra,nemo,occ);
        Kconstnemo=Kconst(nemo);
        truncate(world,Kconstnemo);
    }

    vecfuncT rhsconst=add(world,Vpsi2b,sub(world,Jconstnemo,Kconstnemo));
    truncate(world,rhsconst);

    // the part of the response that is contained in the occupied space
    const vecfuncT parallel=parallel_CPHF(nemo,iatom,iaxis);

    // construct unperturbed operators
    const Coulomb J(world,this);
    const Exchange K(world,this,0);
    const XCOperator xc(world,this,0);
    const Nuclear V(world,this);

    for (int iter=0; iter<25; ++iter) {

        const vecfuncT xi_complete=sub(world,xi,parallel);

        // make the rhs
        START_TIMER(world);
        vecfuncT Kxi;
        if (is_dft()) {
            Kxi=xc(xi);
            scale(world,Kxi,-1.0);
        } else {
            Kxi=K(xi);
        }
        vecfuncT Vpsi1=add(world,V(xi),sub(world,J(xi),Kxi));
        truncate(world,Vpsi1);
        END_TIMER(world, "CPHF: make rhs1");

        START_TIMER(world);

        // construct perturbed operators
        Coulomb Jp(world);
        // 2: closed shell, 2: cphf
        real_function_3d density_pert=4.0*make_density(occ,R2nemo,xi_complete);
        Jp.potential()=Jp.compute_potential(density_pert);

        vecfuncT Kp;
        if (is_dft()) {
            // reconstruct the full perturbed density
            real_function_3d full_dens_pt=density_pert + 2.0*RXR*rhonemo;
            real_function_3d gamma=-1.0*xc.make_xc_kernel()*full_dens_pt;
            Kp=mul(world,gamma,nemo);
        } else {
            Exchange Kp1=Exchange(world).small_memory(false).same(true);
            Kp1.set_parameters(R2nemo,xi_complete,occ);
            vecfuncT R2xi=mul(world,R_square,xi_complete);
            truncate(world,R2xi);
            Exchange Kp2=Exchange(world).small_memory(false);
            Kp2.set_parameters(R2xi,nemo,occ);

            Kp=add(world,Kp1(nemo),Kp2(nemo));
        }
        vecfuncT Vpsi2a=sub(world,Jp(nemo),Kp);
        vecfuncT Vpsi2=add(world,Vpsi2a,rhsconst);
        truncate(world,Vpsi2);
        Vpsi2=Q(Vpsi2);
        truncate(world,Vpsi2);

        vecfuncT Vpsi=add(world,Vpsi1,Vpsi2);
        truncate(world,Vpsi);
        END_TIMER(world, "CPHF make rhs2");

        // apply the BSH
        tensorT eps(nmo);
        for (int i = 0; i < nmo; ++i) eps(i) = std::min(-0.05, fock(i, i));
        std::vector<poperatorT> ops = calc->make_bsh_operators(world, eps);

        // add the coupling elements in case of localized orbitals
        if (get_calc()->param.localize) {
            Tensor<double> fcopy=copy(fock);
            for (int i = 0; i < nmo; ++i) fcopy(i, i) -= eps(i);
            vecfuncT fnemo= transform(world, xi, fcopy, trantol(), true);
            gaxpy(world, 1.0, Vpsi, -1.0, fnemo);
        }

        // apply the BSH operator on the wave function
        START_TIMER(world);
        scale(world, Vpsi, -2.0);
        vecfuncT tmp = apply(world, ops, Vpsi);
        truncate(world, tmp);
        END_TIMER(world, "apply BSH");

        tmp=Q(tmp);
        truncate(world,tmp);

        vecfuncT residual = sub(world, xi, tmp);

        std::vector<double> rnorm = norm2s(world, residual);
        double rms, maxval;
        calc->vector_stats(rnorm, rms, maxval);
        if (world.rank() == 0)
            print("CPHF BSH residual: rms", rms, "   max", maxval);

        const double norm = norm2(world,xi);
        xi=tmp;

        if (rms/norm<proto.dconv) break;
    }
    return xi;

}

std::vector<vecfuncT> Nemo::compute_all_cphf() {

    const Tensor<double> fock=this->compute_fock_matrix(get_calc()->amo,
            get_calc()->get_aocc());
    print("fock");
    print(fock);

    const int natom=molecule().natom();
    std::vector<vecfuncT> xi(3*natom);

    // read CPHF vectors from file if possible
    if (get_calc()->param.read_cphf) {
        xi.resize(3*natom);
        for (std::size_t i=0; i<xi.size(); ++i) {
            load_function(xi[i],"xi_"+stringify(i));
        }
        return xi;
    }

    for (protocol p(*this); not p.finished(); ++p) {
        set_protocol(p.current_prec);

        if (world.rank()==0) {
            printf("\nstarting CPHF equations at time %8.1fs \n",wall_time());
        }
        // double loop over all nuclear displacements
        for (int i=0, iatom=0; iatom<natom; ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis) {
                if (xi[i].size()>0) {
                    for (real_function_3d& xij : xi[i]) xij.set_thresh(p.current_prec);
                }
                xi[i]=cphf(iatom,iaxis,fock,xi[i],p);
                save_function(xi[i],"xi_"+stringify(i));
                ++i;
            }
        }
        printf("\nfinished CPHF equations at time %8.1fs \n",wall_time());
    }

    // reset the initial thresholds
    set_protocol(get_calc()->param.econv);

    if (world.rank()==0) print("\nadding the inhomogeneous part to xi\n");

    // double loop over all nuclear displacements
    for (int i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {

            const vecfuncT& nemo=calc->amo;
            load_function(xi[i],"xi_"+stringify(i));
            const vecfuncT parallel=parallel_CPHF(nemo,iatom,iaxis);
            xi[i]=sub(world,xi[i],parallel);
            truncate(world,xi[i]);
            save_function(xi[i],"xi_"+stringify(i));
            ++i;
        }
    }

    if (world.rank()==0) {
        printf("finished solving the CPHF equations at time %8.1fs \n", wall_time());
    }

    return xi;

}

vecfuncT Nemo::parallel_CPHF(const vecfuncT& nemo, const int iatom,
        const int iaxis) const {

    NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor2(rxr_func).truncate_on_project();
    vecfuncT RXRnemo=mul(world,RXR,nemo);   // skipping factor 2
    truncate(world,RXRnemo);

    Tensor<double> FRXRF=matrix_inner(world,nemo,RXRnemo);  // skipping factor 0.5
    vecfuncT parallel=transform(world,nemo,FRXRF);
    truncate(world,parallel);
    return parallel;

}

/// returns the vibrational frequencies

/// @param[in]  hessian the hessian matrix (not mass-weighted)
/// @param[in]  project_tr whether to project out translation and rotation
/// @param[in]  print_hessian   whether to print the hessian matrix
/// @return the frequencies in atomic units
Tensor<double> Nemo::compute_frequencies(const Tensor<double>& hessian,
        const bool project_tr=true, const bool print_hessian=false) const {

    // compute mass-weighing matrices
    Tensor<double> M=massweights(molecule());
    Tensor<double> Minv(3*molecule().natom(),3*molecule().natom());
    for (int i=0; i<3*molecule().natom(); ++i) Minv(i,i)=1.0/M(i,i);

    // mass-weight the hessian
    Tensor<double> mwhessian=inner(M,inner(hessian,M));

    // remove translation and rotation
    if (project_tr) MolecularOptimizer::remove_external_dof(mwhessian,molecule());

    if (print_hessian) {
        if (project_tr) {
            print("mass-weighted hessian with translation and rotation projected out");
        } else {
            print("mass-weighted unprojected hessian");
        }
        Tensor<double> mmhessian=inner(Minv,inner(mwhessian,Minv));
        print(mwhessian);
        print("mass-weighted unprojected hessian; mass-weighing undone");
        print(mmhessian);
    }

    Tensor<double> normalmodes,freq;
    syev(mwhessian,normalmodes,freq);
    for (long i=0; i<freq.size(); ++i) {
        if (freq(i)>0.0) freq(i)=sqrt(freq(i)); // real frequencies
        else freq(i)=-sqrt(-freq(i));           // imaginary frequencies
    }

    // compute the reduced mass
    Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule());
    Tensor<double> L=copy(normalmodes);
    Tensor<double> DL=inner(D,L);
    Tensor<double> MDL=inner(M,DL);

    for (int i=0; i<3*molecule().natom(); ++i) {
        double mu=0.0;
        for (int j=0; j<3*molecule().natom(); ++j) {
            mu+=MDL(j,i)*MDL(j,i);
        }
        mu=1.0/(mu*constants::atomic_mass_in_au);
        print("reduced mass for mode in amu",i,mu);
    }


    return freq;
}

/// compute the mass-weighting matrix for the hessian
Tensor<double> Nemo::massweights(const Molecule& molecule) const {

    Tensor<double> M(3*molecule.natom(),3*molecule.natom());
    for (int i=0; i<molecule.natom(); i++) {
        const double sqrtmass=1.0/sqrt(molecule.get_atom(i).get_mass_in_au());
        M(3*i  ,3*i  )=sqrtmass;
        M(3*i+1,3*i+1)=sqrtmass;
        M(3*i+2,3*i+2)=sqrtmass;
    }
    return M;
}

} // namespace madness
