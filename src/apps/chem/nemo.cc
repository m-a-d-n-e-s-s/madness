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

struct dens_inv{

    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& inv) const {
//        real_tensor U,
//        const real_tensor& rho) const {
ITERATOR(U,
     double d = t(IND);
     double p = inv(IND);
         U(IND) = d/p;
     );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}

};


double Nemo::value(const Tensor<double>& x) {

    // fast return if the reference is already solved at this geometry
	double xsq = x.sumsq();
	if (xsq == coords_sum)
		return calc->current_energy;

	calc->molecule.set_all_coords(x.reshape(calc->molecule.natom(), 3));
	coords_sum = xsq;

	if (world.rank()==0) {
	    print("\n");
	    calc->molecule.print();
	}

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

        protocol p(*this);

        // guess: read from file or multiply the guess orbitals with the inverse R
	    if (calc->param.restart) {
	        calc->load_mos(world);
	        p.start_prec=calc->amo[0].thresh();

	    } else {
            calc->initial_guess(world);
	        real_function_3d R_inverse = nuclear_correlation->inverse();
	        calc->amo = mul(world, R_inverse, calc->amo);
	    }


	    for (p.initialize() ; not p.finished(); ++p) {
	        set_protocol(p.current_prec);
	        energy=solve(p);
	    }
	    set_protocol(get_calc()->param.econv);

	    calc->current_energy=energy;
	    if (calc->param.save) calc->save_mos(world);

	    // save the converged orbitals and nemos
	    for (std::size_t imo = 0; imo < calc->amo.size(); ++imo) {
	        save(calc->amo[imo], "nemo" + stringify(imo));
	    }
	}

	// compute the dipole moment
	const real_function_3d rhonemo=2.0*make_density(calc->aocc, calc->amo);
	const real_function_3d rho = (R_square*rhonemo);
	save(rho,"rho");
	save(rhonemo,"rhonemo");
	calc->dipole(world,rho);

	return energy;
}


void Nemo::print_nuclear_corrfac() const {

	// the nuclear correlation function
//	save_function(R,"R");
//	save_function(nuclear_correlation->U2(),"U2");
//	save_function(nuclear_correlation->U1(0),"U1x");
//	save_function(nuclear_correlation->U1(1),"U1y");
//	save_function(nuclear_correlation->U1(2),"U1z");

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

    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    // compute potentials the Fock matrix: J - K + Vnuc
	compute_nemo_potentials(nemo, psi, Jnemo, Knemo, Vnemo, Unemo);

    vecfuncT JKUpsi=add(world, sub(world, Jnemo, Knemo), Unemo);
    tensorT fock=matrix_inner(world,R2nemo,JKUpsi,false);   // not symmetric actually
    Kinetic<double,3> T(world);
    fock+=T(R2nemo,nemo);
    JKUpsi.clear();

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
	vecfuncT psi, Jnemo, Knemo, Vnemo, Unemo;

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
	    vecfuncT R2nemo=mul(world,R_square,nemo);
	    truncate(world,R2nemo);

		// compute potentials the Fock matrix: J - K + Vnuc
		compute_nemo_potentials(nemo, psi, Jnemo, Knemo, Vnemo, Unemo);

		// compute the fock matrix
		vecfuncT JKUpsi=add(world, sub(world, Jnemo, Knemo), Unemo);
		tensorT fock=matrix_inner(world,R2nemo,JKUpsi,false);   // not symmetric actually
		Kinetic<double,3> T(world);
		fock+=T(R2nemo,nemo);
		JKUpsi.clear();


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

		double n1=norm2(world,nemo);
		double n2=norm2(world,tmp);
		print("norm of nemo and GVnemo; ratio ",n1,n2,n1/n2);

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

//	START_TIMER(world);
//	const real_function_3d& Vnuc = calc->potentialmanager->vnuclear();
//	Vnemo = mul(world, Vnuc, nemo);
//	truncate(world, Vnemo);
//    double size=get_size(world,Vnemo);
//	END_TIMER(world, "compute Vnemo "+stringify(size));

	START_TIMER(world);
	Nuclear Unuc(world,this->nuclear_correlation);
	Unemo=Unuc(nemo);
    double size1=get_size(world,Unemo);
	END_TIMER(world, "compute Unemo "+stringify(size1));

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
    return calc->make_density(world,occ,nemo);
}

real_function_3d Nemo::make_density(const tensorT & occ,
        const vecfuncT& bra, const vecfuncT& ket, const bool do_refine) const {

    // density may be twice as precise as the orbitals, if you refine
    if (do_refine) {
        refine(world,bra,false);
        if (&bra!=&ket) refine(world,ket,true);
    }

    vecfuncT vsq = mul(world, bra, ket);
    compress(world, vsq);
    functionT rho = factoryT(world).compressed();
    for (unsigned int i = 0; i < vsq.size(); ++i) {
        if (calc->get_aocc()[i])
            rho.gaxpy(1.0, vsq[i], calc->get_aocc()[i], false);
    }
    world.gop.fence();
    return rho;
}


real_function_3d Nemo::make_ddensity(const real_function_3d& rhonemo,
        const int axis) const {

    // 2 RXR * rhonemo
    NuclearCorrelationFactor::U1_functor U1_func(nuclear_correlation.get(),axis);
    real_function_3d RXR=real_factory_3d(world).functor(U1_func).truncate_on_project();
    real_function_3d term1=-2.0*RXR*rhonemo;

    // R^2 * \nabla \rho
    real_derivative_3d D = free_space_derivative<double, 3>(world,axis);
    real_function_3d rhonemo_copy=copy(rhonemo).refine();
    real_function_3d Drhonemo=D(rhonemo_copy);
    return R_square*(term1+Drhonemo);
}


real_function_3d Nemo::make_laplacian_density(const real_function_3d& rhonemo) const {

    // U1^2 operator
    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
    const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

    real_function_3d result=(2.0*U1dot*rhonemo).truncate();

    // U2 operator
    const Nuclear U_op(world,this->nuclear_correlation);
    const Nuclear V_op(world,this->get_calc().get());

    const real_function_3d Vrho=V_op(rhonemo);  // eprec is important here!
    const real_function_3d Urho=U_op(rhonemo);

    real_function_3d term2=4.0*(Urho-Vrho).truncate();
    result-=term2;

    // derivative contribution: R2 \Delta rhonemo
    real_function_3d laplace_rhonemo=real_factory_3d(world).compressed();
    real_function_3d rhonemo_refined=copy(rhonemo).refine();
    for (int axis=0; axis<3; ++axis) {
        real_derivative_3d D = free_space_derivative<double, 3>(world,axis);
        real_function_3d drhonemo=D(rhonemo_refined).refine();
        smoothen(drhonemo);
        real_function_3d d2rhonemo=D(drhonemo);
        laplace_rhonemo+=d2rhonemo;
    }
    save(laplace_rhonemo,"laplace_rhonemo");

    result+=(laplace_rhonemo).truncate();
    result=(R_square*result).truncate();
    save(result,"d2rho");

    // double check result: recompute the density from its laplacian
    real_function_3d rho_rec=-1./(4.*constants::pi)*(*poisson)(result);
    save(rho_rec,"rho_reconstructed");

    return result;
}



real_function_3d Nemo::kinetic_energy_potential(const vecfuncT& nemo) const {

    const Nuclear U_op(world,this->nuclear_correlation);
    const Nuclear V_op(world,this->get_calc().get());

    const vecfuncT Vnemo=V_op(nemo);  // eprec is important here!
    const vecfuncT Unemo=U_op(nemo);

    // nabla^2 nemo
    Laplacian<double,3> Laplace(world,0.0);
    vecfuncT laplace_nemo=Laplace(nemo);
    real_function_3d laplace_sum=sum(world,laplace_nemo);
    save(laplace_sum,"laplace_sum");

    // result=-2.0*(Unemo-Vnemo)  + laplace_nemo;
//    vecfuncT tmp=sub(world,Unemo,Vnemo);
//    vecfuncT tmp=Unemo;
    vecfuncT tmp=sub(world,Unemo,mul(world,this->nuclear_correlation->U2(),nemo));
    gaxpy(world,1.0,laplace_nemo,-2.0,tmp);
    vecfuncT D2Rnemo=mul(world,R,laplace_nemo);

    // double check result: recompute the density from its laplacian
    vecfuncT nemo_rec=apply(world,*poisson,D2Rnemo);
    scale(world,nemo_rec,-1./(4.*constants::pi));
    vecfuncT Rnemo=mul(world,R,nemo);
    vecfuncT diff=sub(world,Rnemo,nemo_rec);
    double dnorm=norm2(world,diff);
    print("dnorm of laplacian phi ",dnorm);

    // compute \sum_i \phi_i \Delta \phi_i
    real_function_3d phiD2phi=dot(world,Rnemo,D2Rnemo);
    save(phiD2phi,"phiD2phi");

    // compute \sum_i \phi_i \epsilon_i \phi_i
    vecfuncT R2nemo=mul(world,R_square,nemo);
    real_function_3d rho=2.0*dot(world,nemo,R2nemo);

    std::vector<double> eps(nemo.size());
    for (int i=0; i<eps.size(); ++i) eps[i]=calc->aeps(i);
    scale(world,R2nemo,eps);
    real_function_3d phiepsilonphi=dot(world,R2nemo,nemo);

    // divide by the density
    real_function_3d numerator=-0.5*phiD2phi-phiepsilonphi;
    real_function_3d nu_bar=binary_op(numerator,rho,dens_inv());

    // smooth the smooth part of the potential
    SeparatedConvolution<double,3> smooth=SmoothingOperator3D(world,0.001);
    nu_bar=smooth(nu_bar);
    save(nu_bar,"nu_bar_bare_smoothed");

    // reintroduce the nuclear potential *after* smoothing
    real_function_3d uvnuc=0.5*(calc->potentialmanager->vnuclear()-nuclear_correlation->U2());
    nu_bar=nu_bar-uvnuc;

    return nu_bar;
}


/// compute the reduced densities sigma (gamma) for GGA functionals

/// the reduced density is given by
/// \f[
///   \sigma = \nabla\rho_1 \nabla\rho_2
///          = (\nabla R^2 \rho_{R,1} + R^2 \nabla \rho_{R,1}) \cdot
///              (\nabla R^2 \rho_{R,2} + R^2 \nabla \rho_{R,2})
///          = 4R^4 U_1^2 \rho_{R,1}\rho_{R,2}
///             + 2R^4 \vec U_1 \left(\rho_{R,1}\nabla \rho_{R,2} + \nabla\rho_{R,1} \rho_{R,2}\right)
///             + R^4 \nabla \rho_{R,1}\cdot \nabla\rho_{R,2}
/// \f]
///
real_function_3d Nemo::make_sigma(const real_function_3d& rho1,
        const real_function_3d& rho2) const {

    const double tight=FunctionDefaults<3>::get_thresh()*0.001;
    const double thresh=FunctionDefaults<3>::get_thresh();
    FunctionDefaults<3>::set_thresh(tight);

    // do refine to have sigma more precise
    std::vector<real_function_3d> drho1=calc->nabla(rho1,true);
    std::vector<real_function_3d> drho2=calc->nabla(rho2,true);

    // first term
    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
    const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1);
    real_function_3d result=(4.0*U1dot*rho1*rho2);

    std::vector<real_function_3d> uvec=nuclear_correlation->U1vec();
    real_function_3d term2=-2.0*(rho1*dot(world,uvec,drho2) + rho2*dot(world,uvec,drho1));

    real_function_3d term3=dot(world,drho1,drho2);
    FunctionDefaults<3>::set_thresh(thresh);
    result+=term2+term3;
    return result*R_square*R_square;

}



/// compute the nuclear gradients
Tensor<double> Nemo::gradient(const Tensor<double>& x) {
    START_TIMER(world);

    const vecfuncT& nemo=calc->amo;

    // the pseudo-density made up of the square of the nemo orbitals
    functionT rhonemo = make_density(calc->aocc, nemo).scale(2.0);
    rhonemo=rhonemo.refine();

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
        bra[axis]=(Drhonemo-tmp);
    }

    Tensor<double> grad(3*calc->molecule.natom());

    for (int iatom=0; iatom<calc->molecule.natom(); ++iatom) {
        NuclearCorrelationFactor::square_times_V_functor r2v(nuclear_correlation.get(),
                calc->molecule,iatom);

        for (int axis=0; axis<3; axis++) {
            grad(3*iatom + axis)=-inner(bra[axis],r2v);
        }
    }

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

    const bool hessdebug=(false and (world.rank()==0));

    const int natom=molecule().natom();
    const vecfuncT& nemo=get_calc()->amo;
    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    Tensor<double> hessian(3*natom,3*natom);

    // the perturbed MOs determined by the CPHF equations
    std::vector<vecfuncT> xi=compute_all_cphf();

    START_TIMER(world);

    // compute the derivative of the density d/dx rho; 2: closed shell
    const real_function_3d rhonemo=2.0*make_density(get_calc()->get_aocc(),nemo);
    const real_function_3d dens=R_square*rhonemo;
    vecfuncT drho(3);
    drho[0]=make_ddensity(rhonemo,0);
    drho[1]=make_ddensity(rhonemo,1);
    drho[2]=make_ddensity(rhonemo,2);

    // compute the perturbed densities
    // \rho_pt = R2 F_i F_i^X + R^X R2 F_i F_i
    vecfuncT dens_pt(3*natom);
    for (int iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            dens_pt[i]=4.0*make_density(calc->get_aocc(),nemo,xi[i]);
            NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,-1);
            const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();
            dens_pt[i]=R_square*(dens_pt[i]+2.0*RXR*rhonemo);//.truncate();
            save(dens_pt[i],"fulldens_pt"+stringify(i));
        }
    }

    // add the electronic contribution to the hessian
    for (int iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            for (int jatom=0; jatom<natom; ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    int j=jatom*3 + jaxis;

                    MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                    double result=inner(dens_pt[i],mdf);

                    // integration by parts
                    if (iatom==jatom) result+=inner(drho[iaxis],mdf);

                    // skip diagonal elements because they are extremely noisy!
                    // use translational symmetry to reconstruct them from other
                    // hessian matrix elements (see below)
                    if (i==j) result=0.0;
                    hessian(i,j)=result;
                }
            }
        }
    }
//    if (hessdebug) {
        print("\n raw electronic Hessian (a.u.)\n");
        print(hessian);
//    }

    Tensor<double> asymmetric=0.5*(hessian-transpose(hessian));
    const double max_asymmetric=asymmetric.absmax();
    if (hessdebug) {
        print("\n asymmetry in the electronic Hessian (a.u.)\n");
        print(asymmetric);
    }
    if (world.rank()==0) print("max asymmetric element in the Hessian matrix: ",max_asymmetric);

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

//    if (hessdebug) {
        print("\n electronic Hessian (a.u.)\n");
        print(hessian);
//    }

    // add the nuclear-nuclear contribution
    hessian+=molecule().nuclear_repulsion_hessian();
        print("\n nuclear Hessian (a.u.)\n");
        print(molecule().nuclear_repulsion_hessian());

//    if (hessdebug) {
        print("\n Hessian (a.u.)\n");
        print(hessian);
//    }
    END_TIMER(world, "compute hessian");

    Tensor<double> normalmodes;
    Tensor<double> frequencies=MolecularOptimizer::compute_frequencies(molecule(),
            hessian,normalmodes,false,hessdebug);

    if (hessdebug) {
        print("\n vibrational frequencies (unprojected) (a.u.)\n");
        print(frequencies);
        print("\n vibrational frequencies (unprojected) (cm-1)\n");
        print(constants::au2invcm*frequencies);
    }

    frequencies=MolecularOptimizer::compute_frequencies(molecule(),hessian,
            normalmodes,true,hessdebug);
    Tensor<double> intensities=compute_IR_intensities(normalmodes,dens_pt);
    Tensor<double> reducedmass=MolecularOptimizer::compute_reduced_mass(
            molecule(),normalmodes);

    if (world.rank() == 0) {
        print("\nprojected vibrational frequencies (cm-1)\n");
        printf("frequency in cm-1   ");
        for (int i=0; i<frequencies.size(); ++i) {
            printf("%10.3f",constants::au2invcm*frequencies(i));
        }
        printf("\n");
        printf("intensity in km/mol ");
        for (int i=0; i<intensities.size(); ++i) {
            printf("%10.3f",intensities(i));
        }
        printf("\n");
        printf("reduced mass in amu ");
        for (int i=0; i<intensities.size(); ++i) {
            printf("%10.3f",reducedmass(i));
        }
        printf("\n\n");
    }

    return hessian;

}

Tensor<double> Nemo::make_incomplete_hessian() const {

    const int natom=molecule().natom();
    vecfuncT& nemo=get_calc()->amo;
    refine(world,nemo);
    real_function_3d rhonemo=2.0*make_density(get_calc()->get_aocc(),nemo);
    real_function_3d rho=R_square*rhonemo;

    Tensor<double> incomplete_hessian=molecule().nuclear_repulsion_hessian();

    // compute the perturbed densities (partial only!)
    // \rho_pt = R2 F_i F_i^X + R^X R2 F_i F_i
    vecfuncT dens_pt(3*natom);
    for (int iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;
            NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
            const real_function_3d RXR=real_factory_3d(world).functor(rxr_func);//.truncate_on_project();
            dens_pt[i]=2.0*RXR*rhonemo;//.truncate();
        }
    }

    // add the electronic contribution to the hessian
    for (int iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            for (int jatom=0; jatom<natom; ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    int j=jatom*3 + jaxis;

                    MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                    double result=inner(dens_pt[i],mdf);

                    // no integration by parts
                    MolecularSecondDerivativeFunctor m2df(molecule(), jatom, jaxis,iaxis);
                    if (iatom==jatom) result+=inner(rho,m2df);

                    // skip diagonal elements because they are extremely noisy!
                    // use translational symmetry to reconstruct them from other
                    // hessian matrix elements (see below)
                    incomplete_hessian(i,j)+=result;
                    if (i==j) incomplete_hessian(i,j)=0.0;
                }
            }
        }
    }

    return incomplete_hessian;
}

/// compute the complementary incomplete hessian

/// @param[in]  xi the response functions including the parallel part
Tensor<double> Nemo::make_incomplete_hessian_response_part(
        const std::vector<vecfuncT>& xi) const {

    int natom=calc->molecule.natom();
    const vecfuncT& nemo=calc->amo;

    Tensor<double> complementary_hessian(3*natom,3*natom);
    for (int i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis, ++i) {

            real_function_3d dens_pt=dot(world,xi[i],nemo);
            dens_pt=4.0*R_square*dens_pt;
            Tensor<double> h(3*molecule().natom());

            for (int jatom=0, j=0; jatom<molecule().natom(); ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis, ++j) {
                    if ((iatom==jatom) and (iaxis==jaxis)) continue;
                    MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                    h(j)=inner(dens_pt,mdf);
                }
            }
            complementary_hessian(i,_)=h;
        }
    }
    return complementary_hessian;
}


vecfuncT Nemo::make_cphf_constant_term(const int iatom, const int iaxis,
        const vecfuncT& R2nemo, const real_function_3d& rhonemo) const {
    // guess for the perturbed MOs
    const vecfuncT nemo=calc->amo;
    const int nmo=nemo.size();

    const Tensor<double> occ=get_calc()->get_aocc();
    QProjector<double,3> Q(world,R2nemo,nemo);

    START_TIMER(world);
    DNuclear Dunuc(world,this,iatom,iaxis);
    vecfuncT Vpsi2b=Dunuc(nemo);
    truncate(world,Vpsi2b);
    END_TIMER(world,"tag3");

    // construct some intermediates
    NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();

    // part of the Coulomb operator with the derivative of the NCF
    // J <- \int dr' 1/|r-r'| \sum_i R^XR F_iF_i
    START_TIMER(world);
    Coulomb Jconst(world);
    Jconst.potential()=Jconst.compute_potential(2.0*RXR*rhonemo);        // factor 2 for cphf
    vecfuncT Jconstnemo=Jconst(nemo);
    truncate(world,Jconstnemo);
    END_TIMER(world,"tag4");
    START_TIMER(world);

    // part of the exchange operator with the derivative of the NCF
    // K <- \sum_k |F_k> \int dr' 1/|r-r'| 2R^XR F_k F_i
    // there is no constant term for DFT, since the potentials are not
    // linear in the density
    vecfuncT Kconstnemo=zero_functions_compressed<double,3>(world,nmo);
    if (not is_dft()) {
        Exchange Kconst=Exchange(world).small_memory(false);
        vecfuncT kbra=mul(world,RXR,nemo);
        scale(world,kbra,2.0);
        truncate(world,kbra);
        Kconst.set_parameters(kbra,nemo,occ);
        Kconstnemo=Kconst(nemo);
        truncate(world,Kconstnemo);
    }
    END_TIMER(world,"tag5");
    START_TIMER(world);

    vecfuncT rhsconst=Vpsi2b+Jconstnemo-Kconstnemo;
    truncate(world,rhsconst);
    rhsconst=Q(rhsconst);
    END_TIMER(world,"tag6");
    return rhsconst;

}


/// solve the CPHF equation for the nuclear displacements

/// @param[in]  iatom   the atom A to be moved
/// @param[in]  iaxis   the coordinate X of iatom to be moved
/// @return     \frac{\partial}{\partial X_A} \varphi
vecfuncT Nemo::solve_cphf(const int iatom, const int iaxis, const Tensor<double> fock,
        const vecfuncT& guess, const vecfuncT& rhsconst,
        const Tensor<double> incomplete_hessian, const vecfuncT& parallel,
        const protocol& proto) const {

    print("\nsolving nemo cphf equations for atom, axis",iatom,iaxis);
    START_TIMER(world);

    vecfuncT xi=guess;
    // guess for the perturbed MOs
    const vecfuncT nemo=calc->amo;
    const int nmo=nemo.size();
    const Tensor<double> occ=get_calc()->get_aocc();
    const real_function_3d rhonemo=2.0*make_density(occ,nemo); // closed shell
    NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();

    int ii=3*iatom+iaxis;
    Tensor<double> ihr=incomplete_hessian(ii,_);
    Tensor<double> old_h(3*molecule().natom());


    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);
    QProjector<double,3> Q(world,R2nemo,nemo);

    END_TIMER(world,"tag1");
    START_TIMER(world);
    // construct quantities that are independent of xi

    // construct the BSH operator
    tensorT eps(nmo);
    for (int i = 0; i < nmo; ++i) eps(i) = fock(i, i);
    std::vector<poperatorT> bsh = calc->make_bsh_operators(world, eps);
    END_TIMER(world,"tag2");

    // derivative of the (regularized) nuclear potential
    START_TIMER(world);

    // construct the KAIN solver
    typedef allocator<double, 3> allocT;
    typedef XNonlinearSolver<vecfunc<double, 3>, double, allocT> solverT;
    allocT alloc(world, nemo.size());
    solverT solver(allocT(world, nemo.size()));
    solver.set_maxsub(5);

    // construct unperturbed operators
    const Coulomb J(world,this);
    const Exchange K=Exchange(world,this,0).small_memory(false);
    const XCOperator xc(world,this,0);
    const Nuclear V(world,this);
    END_TIMER(world,"tag8");
    Tensor<double> h_diff(3l);

    for (int iter=0; iter<10; ++iter) {

        const vecfuncT xi_complete=xi-parallel;

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
        // factor 4 from: closed shell (2) and cphf (2)
        real_function_3d density_pert=4.0*make_density(occ,R2nemo,xi_complete);
        Jp.potential()=Jp.compute_potential(density_pert);
        vecfuncT Kp;
        if (is_dft()) {
            // reconstruct the full perturbed density: do not truncate!
            const real_function_3d full_dens_pt=(density_pert + 2.0*RXR*rhonemo);
            const XCOperator xc1(world,this,0);
            real_function_3d gamma=-1.0*xc1.apply_xc_kernel(full_dens_pt);
            Kp=mul(world,gamma,nemo);
            truncate(world,Kp);
        } else {
            Exchange Kp1=Exchange(world).small_memory(false).same(true);
            Kp1.set_parameters(R2nemo,xi_complete,occ);
            vecfuncT R2xi=mul(world,R_square,xi_complete);
            truncate(world,R2xi);
            Exchange Kp2=Exchange(world).small_memory(false);
            Kp2.set_parameters(R2xi,nemo,occ);

            Kp=add(world,Kp1(nemo),Kp2(nemo));
        }
        vecfuncT Vpsi2=Jp(nemo)-Kp+rhsconst;
        truncate(world,Vpsi2);
        Vpsi2=Q(Vpsi2);
        truncate(world,Vpsi2);

        vecfuncT Vpsi=Vpsi1+Vpsi2;
        truncate(world,Vpsi);
        END_TIMER(world, "CPHF make rhs2");


        // add the coupling elements in case of localized orbitals
        if (get_calc()->param.localize) {
            Tensor<double> fcopy=copy(fock);
            for (int i = 0; i < nmo; ++i) fcopy(i, i) -= eps(i);
            vecfuncT fnemo= transform(world, xi, fcopy, trantol(), true);
            gaxpy(world, 1.0, Vpsi, -1.0, fnemo);
        }

        // apply the BSH operator on the wave function
        START_TIMER(world);
        vecfuncT tmp = apply(world, bsh, -2.0*Vpsi);
        truncate(world, tmp);
        END_TIMER(world, "apply BSH");

        tmp=Q(tmp);
        truncate(world,tmp);

        vecfuncT residual = xi-tmp;

        std::vector<double> rnorm = norm2s(world, residual);
        double rms, maxval;
        calc->vector_stats(rnorm, rms, maxval);
        const double norm = norm2(world,xi);

        if (rms < 1.0) {
            xi = (solver.update(xi, residual)).x;
        } else {
            xi = tmp;
        }

        // measure for hessian matrix elements
        real_function_3d dens_pt=dot(world,xi-parallel,nemo);
        dens_pt=4.0*R_square*dens_pt;
        Tensor<double> h(3*molecule().natom());
        for (int jatom=0, j=0; jatom<molecule().natom(); ++jatom) {
            for (int jaxis=0; jaxis<3; ++jaxis, ++j) {
                if ((iatom==jatom) and (iaxis==jaxis)) continue;
                MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                h(j)=inner(dens_pt,mdf);
            }
        }
        h_diff=h-old_h;

        if (world.rank() == 0)
            print("xi_"+stringify(ii),"CPHF BSH residual: rms", rms,
                    "   max", maxval, "H, Delta H", h.normf(), h_diff.normf());
        print("h+ihr");
        print(h+ihr);
        old_h=h;

        if ((proto.dconv<5.e-4) and iter==2) break;
        if (rms/norm<proto.dconv) break;
    }
    return xi;

}


std::vector<vecfuncT> Nemo::compute_all_cphf() {

    const int natom=molecule().natom();
    std::vector<vecfuncT> xi(3*natom);
    const vecfuncT& nemo=calc->amo;

    // read CPHF vectors from file if possible
    if (get_calc()->param.read_cphf) {
        xi.resize(3*natom);
        for (int i=0; i<3*natom; ++i) {

            load_function(xi[i],"xi_"+stringify(i));
            real_function_3d dens_pt=dot(world,xi[i],nemo);
            dens_pt=4.0*R_square*dens_pt;
            save(dens_pt,"dens_pt"+stringify(i));
        }


        return xi;
    }


    timer t1(world);
    vecfuncT R2nemo=mul(world,R_square,nemo);
    const real_function_3d rhonemo=2.0*make_density(calc->aocc, nemo);
    t1.tag("make density");

    // compute some intermediates

    // compute those contributions to the hessian that do not depend on the response
    Tensor<double> incomplete_hessian=make_incomplete_hessian();
    t1.tag("make incomplete hessian");
    print(incomplete_hessian);

    // compute the initial guess for the response
    const Tensor<double> fock=compute_fock_matrix(nemo,get_calc()->get_aocc());
    const int nmo=nemo.size();
    tensorT eps(nmo);
    for (int i = 0; i < nmo; ++i) eps(i) = fock(i, i);
    std::vector<poperatorT> bsh = calc->make_bsh_operators(world, eps);
    t1.tag("make fock matrix");

    // construct the leading and constant term rhs involving the derivative
    // of the nuclear potential
    std::vector<vecfuncT> parallel(3*natom);
    std::vector<vecfuncT> rhsconst(3*natom);

    for (int i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis, ++i) {
            parallel[i]=zero_functions_compressed<double,3>(world,3*natom);
            rhsconst[i]=zero_functions_compressed<double,3>(world,3*natom);

            timer t2(world);
            parallel[i]=compute_cphf_parallel_term(iatom,iaxis);
            rhsconst[i]=make_cphf_constant_term(iatom,iaxis,R2nemo,rhonemo);
            t2.end("tag1 xi_"+stringify(i));
        }
    }
    t1.tag("make constant and parallel terms");


    // initial guess from the constant rhs
    for (std::size_t i=0; i<rhsconst.size(); ++i) {
        xi[i]=zero_functions_compressed<double,3>(world,3*natom);

        START_TIMER(world);
        xi[i]=apply(world, bsh, -2.0*rhsconst[i]);
        truncate(world,xi[i]);
        END_TIMER(world,"tag2 xi_"+stringify(i));
    }
    t1.tag("make initial guess");

    // compute a first guess for the hessian
    std::vector<vecfuncT> full_xi(3*natom);
    for (int i=0; i<3*natom; ++i) full_xi[i]=xi[i]-parallel[i];
    Tensor<double> complementary_hessian=make_incomplete_hessian_response_part(full_xi);
    t1.tag("make complementary hessian");
    print("incomplete, complementary, and full hessian matrix");
    print(incomplete_hessian);
    print(complementary_hessian);
    print(incomplete_hessian+complementary_hessian);


    // solve the response equations
    protocol p(*this);

    for (p.initialize() ;not p.finished(); ++p) {
        set_protocol(p.current_prec);

        if (world.rank()==0) {
            printf("\nstarting CPHF equations at time %8.1fs \n",wall_time());
        }

        // double loop over all nuclear displacements
        for (int i=0, iatom=0; iatom<natom; ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis, ++i) {
                if (xi[i].size()>0) {
                    for (real_function_3d& xij : xi[i]) xij.set_thresh(p.current_prec);
                }
                xi[i]=solve_cphf(iatom,iaxis,fock,xi[i],rhsconst[i],
                        incomplete_hessian,parallel[i],p);
                save_function(xi[i],"xi_"+stringify(i));
            }
        }
        if (world.rank()==0) {
            printf("\nfinished CPHF equations at time %8.1fs \n",wall_time());
        }

        for (int i=0; i<3*natom; ++i) full_xi[i]=xi[i]-parallel[i];
        Tensor<double> complementary_hessian=make_incomplete_hessian_response_part(full_xi);
        print("full hessian matrix");
        print(incomplete_hessian+complementary_hessian);

    }

    // reset the initial thresholds
    set_protocol(get_calc()->param.econv);

    if (world.rank()==0) print("\nadding the inhomogeneous part to xi\n");

    // double loop over all nuclear displacements
    for (int i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis, ++i) {

            load_function(xi[i],"xi_"+stringify(i));
            xi[i]-=parallel[i];
            truncate(world,xi[i]);
            save_function(xi[i],"xi_"+stringify(i));
        }
    }

    if (world.rank()==0) {
        printf("finished solving the CPHF equations at time %8.1fs \n", wall_time());
    }

    return xi;

}


vecfuncT Nemo::compute_cphf_parallel_term(const int iatom, const int iaxis) const {

    const vecfuncT& nemo=calc->amo;
    int natom=calc->molecule.natom();
    vecfuncT parallel(nemo.size());
    NuclearCorrelationFactor::RX_functor rxr_func(nuclear_correlation.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();
    vecfuncT RXRnemo=mul(world,RXR,nemo);   // skipping factor 2
    truncate(world,RXRnemo);

    Tensor<double> FRXRF=matrix_inner(world,nemo,RXRnemo);  // skipping factor 0.5
    parallel=transform(world,nemo,FRXRF);
    truncate(world,parallel);
    return parallel;

}



Tensor<double> Nemo::compute_IR_intensities(const Tensor<double>& normalmodes,
        const vecfuncT& dens_pt) const {

    // compute the matrix of the normal modes: x -> q
    Tensor<double> M=molecule().massweights();
    Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule());
    Tensor<double> DL=inner(D,normalmodes);
    Tensor<double> nm=inner(M,DL);

    // square of the dipole derivative wrt the normal modes Q
    Tensor<double> mu_Qxyz2(dens_pt.size());

    // compute the derivative of the dipole moment wrt nuclear displacements X
    for (std::size_t component=0; component<3; ++component) {
        // electronic and nuclear dipole derivative wrt nucl. displacements X
        Tensor<double> mu_X(dens_pt.size()), mu_X_nuc(dens_pt.size());

        for (int iatom=0; iatom<molecule().natom(); ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis) {
                int i=iatom*3 + iaxis;
                mu_X(i)=-inner(dens_pt[i],DipoleFunctor(component));
                Tensor<double> dnucdipole=molecule().nuclear_dipole_derivative(iatom,iaxis);
                mu_X_nuc(i)=dnucdipole(component);
            }
        }
        Tensor<double> mu_X_total=mu_X+mu_X_nuc;
        Tensor<double> mu_Q=inner(nm,mu_X_total,0,0);
        mu_Qxyz2+=mu_Q.emul(mu_Q);
    }

    // have fun with constants
    // unit of the dipole moment:
    const double mu_au_to_SI=constants::atomic_unit_of_electric_dipole_moment;
    // unit of the dipole derivative
    const double dmu_au_to_SI=mu_au_to_SI/constants::atomic_unit_of_length;
    // unit of the mass-weighted dipole derivative
    const double dmuq_au_to_SI=dmu_au_to_SI*sqrt(constants::atomic_mass_in_au/
            constants::atomic_mass_constant);
    // some more fundamental constants
    const double pi=constants::pi;
    const double N=constants::Avogadro_constant;
    const double e=constants::dielectric_constant;
    const double c=constants::speed_of_light_in_vacuum;
    // final conversion: Eqs (13) and (14) of
    // J. Neugebauer, M. Reiher, C. Kind, and B. A. Hess,
    // J. Comp. Chem., vol. 23, no. 9, pp. 895–910, Apr. 2002.
    const double conversion=pi*N*dmuq_au_to_SI*dmuq_au_to_SI/(3.0*4.0*pi*e*c*c)/1.e3;
    return mu_Qxyz2.scale(conversion);
}




} // namespace madness
