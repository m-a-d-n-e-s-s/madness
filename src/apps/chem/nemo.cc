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
#include <chem/vibanal.h>
#include <chem/pcm.h>
#include <chem/pointgroupsymmetry.h>


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


class atomic_attraction : public FunctionFunctorInterface<double,3> {
    const Molecule& molecule;
    const size_t iatom;
public:
    atomic_attraction(const Molecule& mol, const size_t iatom1)
        : molecule(mol), iatom(iatom1) {}

    double operator()(const coord_3d& xyz) const {
        return -molecule.atomic_attraction_potential(iatom, xyz[0], xyz[1], xyz[2]);
    }

    std::vector<coord_3d> special_points() const {
        return molecule.get_all_coords_vec();
    }
};


/// compute the nuclear gradients
Tensor<double> NemoBase::compute_gradient(const real_function_3d& rhonemo, const Molecule& molecule) const {

    // the following block computes the gradients more precisely than the
    // direct evaluation of the derivative of the nuclear potential
    vecfuncT bra(3);
    for (int axis=0; axis<3; ++axis) {

        // compute \frac{\partial \rho}{\partial x_i}
        real_derivative_3d D = free_space_derivative<double, 3>(world,axis);
        const real_function_3d Drhonemo=D(rhonemo);

        // compute the second term of the bra
        const real_function_3d tmp=2.0*rhonemo*ncf->U1(axis);
        bra[axis]=(Drhonemo-tmp);
    }

    Tensor<double> grad(3*molecule.natom());

    // linearly scaling code
    bra=bra*R_square;
    compress(world,bra);
    for (size_t iatom=0; iatom<molecule.natom(); ++iatom) {
        atomic_attraction aa(molecule,iatom);
        for (int iaxis=0; iaxis<3; iaxis++) {
            grad(3*iatom + iaxis)=-inner(bra[iaxis],aa);
        }
    }

    // add the nuclear contribution
    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            grad[atom * 3 + axis] +=
                    molecule.nuclear_repulsion_derivative(atom,axis);
        }
    }
    return grad;
}


/// ctor

/// @param[in]	world1	the world
/// @param[in]	calc	the SCF
Nemo::Nemo(World& world, std::shared_ptr<SCF> calc, const std::string inputfile) :
		NemoBase(world), calc(calc), param(calc->param),
		ttt(0.0), sss(0.0), coords_sum(-1.0), ac(world,calc) {

    if (do_pcm()) pcm=PCM(world,this->molecule(),calc->param.pcm_data(),true);

    // reading will not overwrite the derived and defined values
    if (world.rank()==0) param.read(world,inputfile,"dft");
    world.gop.broadcast_serializable(param, 0);


    symmetry_projector=projector_irrep(calc->param.pointgroup())
    		.set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(true);;
    if (world.rank()==0) print("constructed symmetry operator for point group",
    		symmetry_projector.get_pointgroup());
	if (symmetry_projector.get_verbosity()>1) symmetry_projector.print_character_table();

	param.print("dft","end");
}


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

	SCFProtocol p(world,calc->param,"nemo_iterations",calc->param.restart());

	// read (pre-) converged wave function from disk if there is one
	if (calc->param.no_compute() or calc->param.restart()) {
	    calc->load_mos(world);

	    set_protocol(calc->amo[0].thresh());	// set thresh to current value
	    calc->ao=calc->project_ao_basis(world,calc->aobasis);

	} else {		// we need a guess

		FunctionDefaults<3>::set_thresh(p.start_prec);
		set_protocol(p.start_prec);	// set thresh to initial value

	    calc->ao=calc->project_ao_basis(world,calc->aobasis);

		calc->initial_guess(world);
		real_function_3d R_inverse = ncf->inverse();
		calc->amo = R_inverse*calc->amo;
		truncate(world,calc->amo);

	}

	if (not calc->param.no_compute()) {

		p.start_prec=calc->amo[0].thresh();
		p.current_prec=calc->amo[0].thresh();

		for (p.initialize() ; not p.finished(); ++p) {

			set_protocol(p.current_prec);
			calc->current_energy=solve(p);

		}
    }


    // save the converged orbitals and nemos
    for (std::size_t imo = 0; imo < calc->amo.size(); ++imo) {
        save(calc->amo[imo], "nemo" + stringify(imo));
    }

	// compute the dipole moment
	const real_function_3d rhonemo=2.0*make_density(calc->aocc, calc->amo);
	const real_function_3d rho = (R_square*rhonemo);
	save(rho,"rho");
	save(rhonemo,"rhonemo");
	calc->dipole(world,rho);

	if(world.rank()==0) std::cout << "Nemo Orbital Energies: " << calc->aeps << "\n";

	return calc->current_energy;
}



/// localize the nemo orbitals according to Pipek-Mezey or Foster-Boys
vecfuncT Nemo::localize(const vecfuncT& nemo, const double dconv, const bool randomize) const {
        DistributedMatrix<double> dUT;

        const double tolloc = std::min(1e-6,0.01*dconv);
        std::vector<int> aset=calc->group_orbital_sets(world,calc->aeps,
                        calc->aocc, nemo.size());
        // localize using the reconstructed orbitals
        vecfuncT psi = mul(world, R, nemo);
        if(calc->param.localize_method()=="pm")dUT = calc->localize_PM(world, psi, aset, tolloc, 0.1, randomize, true);
        else dUT = calc->localize_boys(world, psi, aset, tolloc, 0.1, randomize);
        dUT.data().screen(trantol());

        vecfuncT localnemo = transform(world, nemo, dUT);
        truncate(world, localnemo);
        normalize(localnemo,R);
        return localnemo;
}

/// compute the Fock matrix from scratch
tensorT Nemo::compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const {
	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT psi, Jnemo, Knemo, pcmnemo, JKVpsi, Unemo;

    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    // compute potentials the Fock matrix: J - K + Vnuc
	compute_nemo_potentials(nemo, psi, Jnemo, Knemo, pcmnemo, Unemo);

//    vecfuncT JKUpsi=add(world, sub(world, Jnemo, Knemo), Unemo);
    vecfuncT JKUpsi=Unemo+Jnemo-Knemo;
    if (do_pcm()) JKUpsi+=pcmnemo;
    tensorT fock=matrix_inner(world,R2nemo,JKUpsi,false);   // not symmetric actually
    Kinetic<double,3> T(world);
    fock+=T(R2nemo,nemo);
    JKUpsi.clear();

	return fock;
}

/// solve the HF equations
double Nemo::solve(const SCFProtocol& proto) {

	// guess has already been performed in value()
	vecfuncT& nemo = calc->amo;
	long nmo = nemo.size();

	// NOTE that nemos are somewhat sensitive to sparse operations (why??)
	// Therefore set all tolerance thresholds to zero, also in the mul_sparse

	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT psi, Jnemo, Knemo, pcmnemo, Unemo;

	double energy = 0.0;
	bool converged = false;
	bool localized=calc->param.do_localize();

	typedef allocator<double, 3> allocT;
	typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
	allocT alloc(world, nemo.size());
	solverT solver(allocT(world, nemo.size()));


	// iterate the residual equations
	for (int iter = 0; iter < calc->param.maxiter(); ++iter) {

	    if (localized) nemo=localize(nemo,proto.dconv,iter==0);
	    std::vector<std::string> str_irreps;
	    if (do_symmetry()) nemo=symmetry_projector(nemo,R_square,str_irreps);
	    if (world.rank()==0) print("orbital irreps",str_irreps);
	    save_function(nemo,"nemo_it"+stringify(iter));
	    vecfuncT R2nemo=mul(world,R_square,nemo);
	    truncate(world,R2nemo);

		// compute potentials the Fock matrix: J - K + Vnuc
		compute_nemo_potentials(nemo, psi, Jnemo, Knemo, pcmnemo, Unemo);

		// compute the fock matrix
		vecfuncT JKUpsi=Unemo+Jnemo-Knemo;
		if (do_pcm()) JKUpsi+=pcmnemo;
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
            normalize(nemo,R);
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
         
        if(localized) calc->aeps=eps;
        vecfuncT fnemo;
        if (localized) fnemo= transform(world, nemo, fock, trantol(), true);

        // Undo the damage
        for (int i = 0; i < nmo; ++i) fock(i, i) += eps(i);

		if (calc->param.orbitalshift()>0.0) {
			if (world.rank()==0) print("shifting orbitals by "
					,calc->param.orbitalshift()," to lower energies");
			eps-=calc->param.orbitalshift();
		}
		std::vector<poperatorT> ops = calc->make_bsh_operators(world, eps);

		// make the potential * nemos term; make sure it's in phase with nemo
		START_TIMER(world);
        vecfuncT Vpsi=Unemo+Jnemo-Knemo;
        if (do_pcm()) Vpsi+=pcmnemo;

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
			nemo_new = solver.update(nemo, residual);
		} else {
			nemo_new = tmp;
		}
		truncate(world,nemo_new);
		normalize(nemo_new,R);

		calc->do_step_restriction(world,nemo,nemo_new,"ab spin case");
		orthonormalize(nemo_new,R);
		nemo=nemo_new;

		if ((norm < proto.dconv) and
				(fabs(energy-oldenergy)<proto.econv))
			converged = true;

		if (calc->param.save()) calc->save_mos(world);

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
    START_TIMER(world);

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

    // possible advantageous to do it this way, the U1 term should cancel
    // with the nuclear potential, might be faster and more precise.
//    // do the kinetic energy again
//    double ke1=0.0;
//    for (int axis = 0; axis < 3; axis++) {
//        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
//        const vecfuncT dnemo = apply(world, D, nemo);
//
//        real_function_3d term1=dot(world,dnemo,dnemo);
//        ke1 += 0.5*inner(term1,R_square);
//
//        vecfuncT term2=R_square*nuclear_correlation->U1(axis)*nemo;
//        ke1 +=-2.0* 0.5 * (inner(world, dnemo, term2)).sum();
//    }
//    ke1 *= 2.0; // closed shell


    const double J = inner(world, R2nemo, Jnemo).sum();
    const double K = inner(world, R2nemo, Knemo).sum();

    int ispin=0;
    double exc=0.0;
    if (calc->xc.is_dft()) {
        XCOperator xcoperator(world,this,ispin);
        exc=xcoperator.compute_xc_energy();
    }

    double pcm_energy=0.0;
    if (do_pcm()) pcm_energy=pcm.compute_pcm_energy();

    const double nucrep = calc->molecule.nuclear_repulsion_energy();

    double energy = ke + J + pe + nucrep + pcm_energy;
    if (is_dft()) energy+=exc;
    else energy-=K;

    if (world.rank() == 0) {
        printf("\n  nuclear and kinetic %16.8f\n", ke + pe);
//        printf("\n  kinetic again %16.8f\n",  ke1);
        printf("              coulomb %16.8f\n", J);
        if (is_dft()) {
            printf(" exchange-correlation %16.8f\n", exc);
        } else {
            printf("             exchange %16.8f\n", -K);
        }
        if (do_pcm()) {
            printf("   polarization (PCM) %16.8f\n", pcm_energy);
        }
        printf("    nuclear-repulsion %16.8f\n", nucrep);
        printf("   regularized energy %16.8f\n", energy);
        printf("  buggy if hybrid functionals are used..\n");
    }
    END_TIMER(world, "compute energy");

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
		vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& pcmnemo,
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
        Exchange<double,3> K=Exchange<double,3>(world,this,ispin).same(false).small_memory(false);
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
        // copy???
        real_function_3d xc_pot = xcoperator.make_xc_potential();

        // compute the asymptotic correction of exchange-correlation potential
        if(do_ac()) {
        	std::cout << "Computing asymtotic correction!\n";
        	double charge = double(molecule().total_nuclear_charge())-calc->param.charge();
        	real_function_3d scaledJ = -1.0/charge*J.potential()*(1.0-calc->xc.hf_exchange_coefficient());
        	xc_pot = ac.apply(xc_pot, scaledJ);
        }

        Knemo=sub(world,Knemo,mul(world,xc_pot,nemo));   // minus times minus gives plus
        truncate(world,Knemo);
        double size=get_size(world,Knemo);
        END_TIMER(world, "compute XCnemo "+stringify(size));
    }


    // compute the solvent (PCM) contribution to the potential
    if (do_pcm()) {
        START_TIMER(world);
        const real_function_3d vpcm = pcm.compute_pcm_potential(J.potential());
        pcmnemo=vpcm*nemo;
        double size=get_size(world,pcmnemo);
        END_TIMER(world, "compute PCMnemo "+stringify(size));
    }

	START_TIMER(world);
	Nuclear Unuc(world,this->ncf);
	Unemo=Unuc(nemo);
    double size1=get_size(world,Unemo);
	END_TIMER(world, "compute Unemo "+stringify(size1));

}



/// return the Coulomb potential
real_function_3d Nemo::get_coulomb_potential(const vecfuncT& psi) const {
	MADNESS_ASSERT(calc->param.spin_restricted());
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
    NuclearCorrelationFactor::U1_functor U1_func(ncf.get(),axis);
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
    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(ncf.get());
    const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

    real_function_3d result=(2.0*U1dot*rhonemo).truncate();

    // U2 operator
    const Nuclear U_op(world,this->ncf);
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

    const Nuclear U_op(world,this->ncf);
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
    vecfuncT tmp=sub(world,Unemo,mul(world,this->ncf->U2(),nemo));
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
    for (std::size_t i=0; i<eps.size(); ++i) eps[i]=calc->aeps(i);
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
    real_function_3d uvnuc=0.5*(calc->potentialmanager->vnuclear()-ncf->U2());
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
    std::vector<real_function_3d> drho1=grad(rho1,true);
    std::vector<real_function_3d> drho2=grad(rho2,true);

    // first term
    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(ncf.get());
    const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1);
    real_function_3d result=(4.0*U1dot*rho1*rho2);

    std::vector<real_function_3d> uvec=ncf->U1vec();
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

    Tensor<double> grad=NemoBase::compute_gradient(rhonemo,calc->molecule);

//    // the following block computes the gradients more precisely than the
//    // direct evaluation of the derivative of the nuclear potential
//    vecfuncT bra(3);
//    for (int axis=0; axis<3; ++axis) {
//
//        // compute \frac{\partial \rho}{\partial x_i}
//        real_derivative_3d D = free_space_derivative<double, 3>(world,axis);
//        real_function_3d Drhonemo=D(rhonemo);
//
//        // compute the second term of the bra
//        real_function_3d tmp=rhonemo*ncf->U1(axis);
//        tmp.scale(2.0);
//        bra[axis]=(Drhonemo-tmp);
//    }
//
//    Tensor<double> grad(3*calc->molecule.natom());
//
//    // linearly scaling code
//    bra=bra*R_square;
//    compress(world,bra);
//    calc->potentialmanager->vnuclear().compress();
//    for (size_t iatom=0; iatom<calc->molecule.natom(); ++iatom) {
//        atomic_attraction aa(calc->molecule,iatom);
//        for (int iaxis=0; iaxis<3; iaxis++) {
//            grad(3*iatom + iaxis)=-inner(bra[iaxis],aa);
//        }
//    }
//
//
////  // quadratically scaling code..
////    for (size_t iatom=0; iatom<calc->molecule.natom(); ++iatom) {
////        NuclearCorrelationFactor::square_times_V_functor r2v(nuclear_correlation.get(),
////                calc->molecule,iatom);
////
////        for (int axis=0; axis<3; axis++) {
////            grad(3*iatom + axis)=-inner(bra[axis],r2v);
////        }
////    }
//
////    // this block is less precise
////    for (size_t iatom=0; iatom<calc->molecule.natom(); ++iatom) {
////        for (int axis=0; axis<3; ++axis) {
////            NuclearCorrelationFactor::square_times_V_derivative_functor r2v(
////                    nuclear_correlation.get(),this->molecule(),iatom,axis);
////            grad(3*iatom + axis)=inner(rhonemo,r2v);
////
////        }
////    }
//
//    // add the nuclear contribution
//    for (size_t atom = 0; atom < calc->molecule.natom(); ++atom) {
//        for (int axis = 0; axis < 3; ++axis) {
//            grad[atom * 3 + axis] +=
//                    calc->molecule.nuclear_repulsion_derivative(atom,axis);
//        }
//    }
//
//
//    END_TIMER(world, "compute gradients");

    if (world.rank() == 0) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print(
              "  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(
              " ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (size_t i = 0; i < calc->molecule.natom(); ++i) {
            const Atom& atom = calc->molecule.get_atom(i);
            printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", int(i),
                   atom.x, atom.y, atom.z, grad[i * 3 + 0], grad[i * 3 + 1],
                   grad[i * 3 + 2]);
        }
    }
    return grad;
}


/// compute the nuclear hessian
Tensor<double> Nemo::hessian(const Tensor<double>& x) {

    const bool hessdebug=(false and (world.rank()==0));

    const size_t natom=molecule().natom();
    const vecfuncT& nemo=get_calc()->amo;
    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    Tensor<double> hessian(3*natom,3*natom);

    // the perturbed MOs determined by the CPHF equations
    std::vector<vecfuncT> xi=compute_all_cphf();

    timer time_hessian(world);

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
    vecfuncT pre_dens_pt(3*natom);
    std::vector<vecfuncT> pre_mo_pt(3*natom);
    std::vector<vecfuncT> R2mo_pt(3*natom);

    for (size_t iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            real_function_3d FXF=4.0*make_density(calc->get_aocc(),nemo,xi[i]);
            NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,-1);
            const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();
            pre_dens_pt[i]=(FXF+2.0*RXR*rhonemo);//.truncate();
            pre_mo_pt[i]=(RXR*nemo + xi[i]);
            R2mo_pt[i]=R_square*pre_mo_pt[i];
        }
    }

    dens_pt=R_square * pre_dens_pt;
    save_function(dens_pt,"full_dens_pt");

    // add the electronic contribution to the hessian
    for (size_t iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            for (size_t jatom=0; jatom<natom; ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    int j=jatom*3 + jaxis;

                    // skip diagonal elements because they are extremely noisy!
                    // use translational symmetry to reconstruct them from other
                    // hessian matrix elements (see below)
                    if (i==j) continue;

                    MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                    hessian(i,j)=inner(dens_pt[i],mdf);

                    // integration by parts < 0 | H^{YX} | 0 >
                    if (iatom==jatom) hessian(i,j)+=inner(drho[iaxis],mdf);

                }
            }
        }
    }
//    if (hessdebug) {
        print("\n raw electronic Hessian (a.u.)\n");
        print(hessian);
//    }
    for (size_t i=0; i<3*natom; ++i) hessian(i,i)=0.0;
    if (calc->param.get<bool>("purify_hessian")) hessian=purify_hessian(hessian);

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
    for (size_t i=0; i<3*natom; ++i) {
        double sum=0.0;
        for (size_t j=0; j<3*natom; j+=3) sum+=hessian(i,j+(i%3));
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

    time_hessian.end("compute hessian");

    Tensor<double> normalmodes;
    Tensor<double> frequencies=compute_frequencies(molecule(),hessian,normalmodes,false,hessdebug);

    if (hessdebug) {
        print("\n vibrational frequencies (unprojected) (a.u.)\n");
        print(frequencies);
        print("\n vibrational frequencies (unprojected) (cm-1)\n");
        print(constants::au2invcm*frequencies);
    }

    frequencies=compute_frequencies(molecule(),hessian, normalmodes,true,hessdebug);
    Tensor<double> intensities=compute_IR_intensities(normalmodes,dens_pt);
    Tensor<double> reducedmass=compute_reduced_mass(molecule(),normalmodes);

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
        printf("done with computing the hessian matrix at time %8.1fs \n",wall_time());
        printf("final energy %16.8f", calc->current_energy);
    }

    return hessian;

}

/// purify and symmetrize the hessian

Tensor<double> Nemo::purify_hessian(const Tensor<double>& hessian) const {

    Tensor<double> purified=copy(hessian);
    double maxasymmetric=0.0;

    const size_t natom=calc->molecule.natom();

    for (size_t iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            for (size_t jatom=0; jatom<natom; ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    int j=jatom*3 + jaxis;

                    double mean=(purified(i,j)+purified(j,i))*0.5;
                    double diff=0.5*fabs(purified(i,j)-purified(j,i));
                    maxasymmetric=std::max(maxasymmetric,diff);

                    unsigned int ZA=calc->molecule.get_atom_number(iatom);
                    unsigned int ZB=calc->molecule.get_atom_number(jatom);
                    if (ZA<ZB) purified(i,j)=purified(j,i);
                    if (ZA>ZB) purified(j,i)=purified(i,j);
                    if (ZA==ZB) {
                        purified(i,j)=mean;
                        purified(j,i)=mean;
                    }
                }
            }
        }
    }
    print("purify: max asymmetric element ",maxasymmetric);
    print("purify: raw hessian ");
    print(hessian);
    print("purify: purified hessian ");
    print(purified);

    return purified;
}



Tensor<double> Nemo::make_incomplete_hessian() const {

    const size_t natom=molecule().natom();
    vecfuncT& nemo=get_calc()->amo;
    refine(world,nemo);
    real_function_3d rhonemo=2.0*make_density(get_calc()->get_aocc(),nemo);
    real_function_3d rho=R_square*rhonemo;

    Tensor<double> incomplete_hessian=molecule().nuclear_repulsion_hessian();

    // compute the perturbed densities (partial only!)
    // \rho_pt = R2 F_i F_i^X + R^X R2 F_i F_i
    vecfuncT dens_pt(3*natom);
    for (size_t iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;
            NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,2);
            const real_function_3d RXR=real_factory_3d(world).functor(rxr_func);//.truncate_on_project();
            dens_pt[i]=2.0*RXR*rhonemo;//.truncate();
        }
    }

    // add the electronic contribution to the hessian
    for (size_t iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            int i=iatom*3 + iaxis;

            for (size_t jatom=0; jatom<natom; ++jatom) {
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

    size_t natom=calc->molecule.natom();
    const vecfuncT& nemo=calc->amo;

    Tensor<double> complementary_hessian(3*natom,3*natom);
    for (size_t i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis, ++i) {

            real_function_3d dens_pt=dot(world,xi[i],nemo);
            dens_pt=4.0*R_square*dens_pt;
            Tensor<double> h(3*molecule().natom());

            for (size_t jatom=0, j=0; jatom<molecule().natom(); ++jatom) {
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


vecfuncT Nemo::make_cphf_constant_term(const size_t iatom, const int iaxis,
        const vecfuncT& R2nemo, const real_function_3d& rhonemo) const {
    // guess for the perturbed MOs
    const vecfuncT nemo=calc->amo;
    const int nmo=nemo.size();

    const Tensor<double> occ=get_calc()->get_aocc();
    QProjector<double,3> Q(world,R2nemo,nemo);

    DNuclear Dunuc(world,this,iatom,iaxis);
    vecfuncT Vpsi2b=Dunuc(nemo);
    truncate(world,Vpsi2b);

    // construct some intermediates
    NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();

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
        Exchange<double,3> Kconst=Exchange<double,3>(world).small_memory(false);
        vecfuncT kbra=2.0*RXR*nemo;
        truncate(world,kbra);
        Kconst.set_parameters(kbra,nemo,occ);
        Kconstnemo=Kconst(nemo);
        truncate(world,Kconstnemo);
    }

    vecfuncT rhsconst=Vpsi2b+Jconstnemo-Kconstnemo;
    truncate(world,rhsconst);
    rhsconst=Q(rhsconst);
    return rhsconst;

}


/// solve the CPHF equation for the nuclear displacements

/// @param[in]  iatom   the atom A to be moved
/// @param[in]  iaxis   the coordinate X of iatom to be moved
/// @return     \frac{\partial}{\partial X_A} \varphi
vecfuncT Nemo::solve_cphf(const size_t iatom, const int iaxis, const Tensor<double> fock,
        const vecfuncT& guess, const vecfuncT& rhsconst,
        const Tensor<double> incomplete_hessian, const vecfuncT& parallel,
        const SCFProtocol& proto, const std::string& xc_data) const {

    print("\nsolving nemo cphf equations for atom, axis",iatom,iaxis);

    vecfuncT xi=guess;
    // guess for the perturbed MOs
    const vecfuncT nemo=calc->amo;
    const int nmo=nemo.size();
    const Tensor<double> occ=get_calc()->get_aocc();
    const real_function_3d rhonemo=2.0*make_density(occ,nemo); // closed shell
    const real_function_3d arho=0.5*R_square*rhonemo;
    NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();

    int ii=3*iatom+iaxis;
    Tensor<double> ihr=incomplete_hessian(ii,_);
    Tensor<double> old_h(3*molecule().natom());


    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);
    QProjector<double,3> Q(world,R2nemo,nemo);

    // construct quantities that are independent of xi

    // construct the BSH operator
    tensorT eps(nmo);
    for (int i = 0; i < nmo; ++i) eps(i) = fock(i, i);
    std::vector<poperatorT> bsh = calc->make_bsh_operators(world, eps);
    for (poperatorT& b : bsh) b->destructive()=true;    // make it memory efficient

    // derivative of the (regularized) nuclear potential

    // construct the KAIN solver
    typedef allocator<double, 3> allocT;
    typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
    allocT alloc(world, nemo.size());
    solverT solver(allocT(world, nemo.size()));
    solver.set_maxsub(5);

    // construct unperturbed operators
    const Coulomb J(world,this);
    const Exchange<double,3> K=Exchange<double,3>(world,this,0).small_memory(false);
    const XCOperator xc(world, xc_data, not calc->param.spin_restricted(), arho, arho);
    const Nuclear V(world,this);

    Tensor<double> h_diff(3l);
    for (int iter=0; iter<10; ++iter) {


        // make the rhs
        START_TIMER(world);
        vecfuncT Vpsi1=V(xi) + J(xi);
        if (is_dft()) {
            Vpsi1+=(xc(xi));
        } else {
            Vpsi1-=(K(xi));
        }

        truncate(world,Vpsi1);
        END_TIMER(world, "CPHF: make rhs1");

        START_TIMER(world);

        // construct perturbed operators
        Coulomb Jp(world);
        const vecfuncT xi_complete=xi-parallel;

        // factor 4 from: closed shell (2) and cphf (2)
        real_function_3d density_pert=4.0*make_density(occ,R2nemo,xi_complete);
        Jp.potential()=Jp.compute_potential(density_pert);
        vecfuncT Kp;
        if (is_dft()) {
            // reconstruct the full perturbed density: do not truncate!
            const real_function_3d full_dens_pt=(density_pert + 2.0*RXR*rhonemo);
            real_function_3d gamma=-1.0*xc.apply_xc_kernel(full_dens_pt);
            Kp=truncate(gamma*nemo);
        } else {
            Exchange<double,3> Kp1=Exchange<double,3>(world).small_memory(false).same(true);
            Kp1.set_parameters(R2nemo,xi_complete,occ);
            vecfuncT R2xi=mul(world,R_square,xi_complete);
            truncate(world,R2xi);
            Exchange<double,3> Kp2=Exchange<double,3>(world).small_memory(false);
            Kp2.set_parameters(R2xi,nemo,occ);
            Kp=truncate(Kp1(nemo) + Kp2(nemo));
        }
        vecfuncT Vpsi2=truncate(Jp(nemo)-Kp+rhsconst);
        Vpsi2=Q(Vpsi2);
        truncate(world,Vpsi2);

        vecfuncT Vpsi=truncate(Vpsi1+Vpsi2);
        Vpsi1.clear();
        Vpsi2.clear();
        END_TIMER(world, "CPHF make rhs2");


        // add the coupling elements in case of localized orbitals
        if (get_calc()->param.do_localize()) {
            Tensor<double> fcopy=copy(fock);
            for (int i = 0; i < nmo; ++i) fcopy(i, i) -= eps(i);
            vecfuncT fnemo= transform(world, xi, fcopy, trantol(), true);
            gaxpy(world, 1.0, Vpsi, -1.0, fnemo);
            truncate(Vpsi);
        }

        // apply the BSH operator on the wave function
        START_TIMER(world);
        scale(world,Vpsi,-2.0);
        vecfuncT tmp = apply(world, bsh,Vpsi);
        Vpsi.clear();
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
            xi = solver.update(xi, residual);
        } else {
            xi = tmp;
        }
        truncate(xi);

        // measure for hessian matrix elements
        real_function_3d dens_pt=dot(world,xi-parallel,nemo);
        dens_pt=4.0*R_square*dens_pt;
        Tensor<double> h(3*molecule().natom());
        for (size_t jatom=0, j=0; jatom<molecule().natom(); ++jatom) {
            for (int jaxis=0; jaxis<3; ++jaxis, ++j) {
                if ((iatom==jatom) and (iaxis==jaxis)) continue;
                MolecularDerivativeFunctor mdf(molecule(), jatom, jaxis);
                h(j)=inner(dens_pt,mdf);
            }
        }
        h_diff=h-old_h;

        if (world.rank() == 0)
            print("xi_"+stringify(ii),"CPHF BSH residual: rms", rms,
                    "   max", maxval, "H, Delta H", h.normf(), h_diff.absmax());
        print("h+ihr");
        print(h+ihr);
        old_h=h;

        if (rms/norm<proto.dconv and (h_diff.absmax()<1.e-2)) break;
    }
    return xi;

}


std::vector<vecfuncT> Nemo::compute_all_cphf() {

    const int natom=molecule().natom();
    std::vector<vecfuncT> xi(3*natom);
    const vecfuncT& nemo=calc->amo;

    // read CPHF vectors from file if possible
    if (get_calc()->param.get<bool>("read_cphf")) {
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
            parallel[i]=compute_cphf_parallel_term(iatom,iaxis);
            rhsconst[i]=make_cphf_constant_term(iatom,iaxis,R2nemo,rhonemo);
        }
    }
    t1.tag("make constant and parallel terms");


    // initial guess from the constant rhs or from restart
    if (get_calc()->param.restart_cphf()) {
        for (int i=0; i<3*natom; ++i) {
            load_function(xi[i],"xi_guess"+stringify(i));
        }
        t1.tag("read guess from file");
    } else {
        for (std::size_t i=0; i<rhsconst.size(); ++i) {
            xi[i]=apply(world, bsh, -2.0*rhsconst[i]);
            truncate(world,xi[i]);
            save_function(xi[i],"xi_guess"+stringify(i));
        }
        t1.tag("make initial guess");
    }


    // do some pre-iterations using fast LDA xc functional
    if (world.rank()==0) {
        print("\ngenerating CPHF guess using the LDA functional\n");
    }

    SCFProtocol preiterations(world,calc->param,"cphf_preiterations",
            calc->param.restart_cphf());
    preiterations.end_prec*=10.0;
    preiterations.initialize();
    for (; not preiterations.finished(); ++preiterations) {
        set_protocol(preiterations.current_prec);

        if (world.rank()==0) {
            printf("\nstarting initial CPHF equations at time %8.1fs \n",wall_time());
        }

        // double loop over all nuclear displacements
        for (int i=0, iatom=0; iatom<natom; ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis, ++i) {
                if (xi[i].size()>0) {
                    for (real_function_3d& xij : xi[i]) xij.set_thresh(preiterations.current_prec);
                }
                xi[i]=solve_cphf(iatom,iaxis,fock,xi[i],rhsconst[i],
                        incomplete_hessian,parallel[i],preiterations,"LDA");
                save_function(xi[i],"xi_guess"+stringify(i));
            }
        }
        if (world.rank()==0) {
            printf("\nfinished CPHF equations at time %8.1fs \n",wall_time());
        }

        std::vector<vecfuncT> full_xi(3*natom);
        for (int i=0; i<3*natom; ++i) full_xi[i]=xi[i]-parallel[i];
        Tensor<double> complementary_hessian=make_incomplete_hessian_response_part(full_xi);
        print("full hessian matrix");
        print(incomplete_hessian+complementary_hessian);

    }

    // solve the response equations
    SCFProtocol p(world,calc->param,"cphf_final_iterations",calc->param.restart_cphf());
    p.start_prec=p.end_prec;
    p.initialize();

    for ( ; not p.finished(); ++p) {
        set_protocol(p.current_prec);

        if (world.rank()==0) {
            printf("\nstarting CPHF equations at time %8.1fs \n",wall_time());
            print("solving CPHF with the density functional",calc->param.xc());
        }

        // double loop over all nuclear displacements
        for (int i=0, iatom=0; iatom<natom; ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis, ++i) {
                if (xi[i].size()>0) {
                    for (real_function_3d& xij : xi[i]) xij.set_thresh(p.current_prec);
                }
                xi[i]=solve_cphf(iatom,iaxis,fock,xi[i],rhsconst[i],
                        incomplete_hessian,parallel[i],p,calc->param.xc());
                save_function(xi[i],"xi_guess"+stringify(i));
            }
        }
        if (world.rank()==0) {
            printf("\nfinished CPHF equations at time %8.1fs \n",wall_time());
        }

        std::vector<vecfuncT> full_xi(3*natom);
        for (int i=0; i<3*natom; ++i) full_xi[i]=xi[i]-parallel[i];
        Tensor<double> complementary_hessian=make_incomplete_hessian_response_part(full_xi);
        print("full hessian matrix");
        print(incomplete_hessian+complementary_hessian);

    }

    // reset the initial thresholds
    set_protocol(get_calc()->param.econv());

    if (world.rank()==0) print("\nadding the inhomogeneous part to xi\n");

    // double loop over all nuclear displacements
    for (int i=0, iatom=0; iatom<natom; ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis, ++i) {

            load_function(xi[i],"xi_guess"+stringify(i));
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


vecfuncT Nemo::compute_cphf_parallel_term(const size_t iatom, const int iaxis) const {

    const vecfuncT& nemo=calc->amo;
    vecfuncT parallel(nemo.size());
    NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,2);
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
    Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule(),{"Tx","Ty","Tz","Rx","Ry","Rz"});
    Tensor<double> DL=inner(D,normalmodes);
    Tensor<double> nm=inner(M,DL);

    // square of the dipole derivative wrt the normal modes Q
    Tensor<double> mu_Qxyz2(dens_pt.size());

    // compute the derivative of the dipole moment wrt nuclear displacements X
    for (std::size_t component=0; component<3; ++component) {
        // electronic and nuclear dipole derivative wrt nucl. displacements X
        Tensor<double> mu_X(dens_pt.size()), mu_X_nuc(dens_pt.size());

        for (size_t iatom=0; iatom<molecule().natom(); ++iatom) {
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
