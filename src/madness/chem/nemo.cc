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

/*!
 \file examples/nemo.cc
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#include<madness/chem/nemo.h>
#include<madness/chem/projector.h>
#include<madness/chem/molecular_optimizer.h>
#include<madness/chem/SCFOperators.h>
#include <madness/constants.h>
#include<madness/chem/vibanal.h>
#include<madness/chem/pcm.h>
#include<madness/chem/pointgroupsymmetry.h>
#include<madness/chem/BSHApply.h>
#include<madness/chem/localizer.h>
#include <madness/mra/macrotaskq.h>
#include <madness/mra/memory_measurement.h>


using namespace madchem;
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


Nemo::Nemo(World& world, const commandlineparser &parser) :
        NemoBase(world),
        calc(std::make_shared<SCF>(world, parser)),
        param(calc->param),
        coords_sum(-1.0),
        ac(world,calc) {
    if (do_pcm()) pcm=PCM(world,this->molecule(),param.pcm_data(),true);

    // reading will not overwrite the derived and defined values
    param.read_input_and_commandline_options(world,parser,"dft");

    symmetry_projector=projector_irrep(param.pointgroup())
            .set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(true);;
    if (symmetry_projector.get_verbosity()>1) symmetry_projector.print_character_table();
    calc->param=param;
};

Nemo::Nemo(World& world, const CalculationParameters& param, const Molecule& molecule) :
        NemoBase(world),
        calc(std::make_shared<SCF>(world, param, molecule)),
        param(calc->param),
        coords_sum(-1.0),
        ac(world,calc) {
    if (do_pcm()) pcm=PCM(world,this->molecule(),param.pcm_data(),true);

    symmetry_projector=projector_irrep(param.pointgroup())
            .set_ordering("keep").set_verbosity(0).set_orthonormalize_irreps(true);;
    if (symmetry_projector.get_verbosity()>1) symmetry_projector.print_character_table();
    calc->param=param;
};

double Nemo::value(const Tensor<double>& x) {

    // fast return if the reference is already solved at this geometry
    if (check_converged(x)) return calc->current_energy;
	double xsq = x.sumsq();

    if (world.rank()==0) print_header2("computing the nemo wave function");

    if ((xsq-calc->molecule.get_all_coords()).normf()>1.e-12) invalidate_factors_and_potentials();

    calc->molecule.set_all_coords(x.reshape(calc->molecule.natom(), 3));
	coords_sum = xsq;

	SCFProtocol p(world,param);

	// read (pre-) converged wave function from disk if there is one
	if (param.no_compute() or param.restart()) {
	    set_protocol(param.econv());	// set thresh to current value
        if (world.rank()==0 and param.print_level()>2) print("reading orbitals from disk");
	    calc->load_mos(world);
        if (world.rank()==0 and param.print_level()>2) {
            print("orbitals are converged to ",calc->converged_for_thresh);
        }
        p.start_prec=calc->converged_for_thresh;

	    calc->ao=calc->project_ao_basis(world,calc->aobasis);

	} else {		// we need a guess

		FunctionDefaults<3>::set_thresh(p.start_prec);
		set_protocol(p.start_prec);	// set thresh to initial value

	    calc->ao=calc->project_ao_basis(world,calc->aobasis);

	    if (not (calc->converged_for_thresh*0.999<p.end_prec)) {
	        // if the orbitals are not converged to the end precision, we need to recompute the ncf
	        calc->initial_guess(world);
	        real_function_3d R_inverse = ncf->inverse();
	        calc->amo = R_inverse*calc->amo;
	        truncate(world,calc->amo);
	    }

	}

    bool skip_solve=(param.no_compute()) or (calc->converged_for_thresh*0.9999<p.end_prec);
    if (skip_solve) {
        if (world.rank()==0) {
            print("skipping the solution of the SCF equations:");
            if (param.no_compute()) print(" -> the option no_compute =1");
            if (calc->converged_for_thresh*0.9999<p.end_prec) print(" -> orbitals are converged to the required threshold of",p.end_prec);
        }

    } else {

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
	Tensor<double> dipole=calc->dipole(world,rho);

	if(world.rank()==0) std::cout << "Nemo Orbital Energies: " << calc->aeps << "\n";

    calc->output_calc_info_schema();

    if (world.rank()==0) print_header2("end computing the nemo wave function");
    return calc->current_energy;
}


/// localize the nemo orbitals according to Pipek-Mezey or Foster-Boys
vecfuncT Nemo::localize(const vecfuncT& nemo, const double dconv, const bool randomize) const {

    Localizer localizer(world, get_calc()->aobasis, molecule(), get_calc()->ao);
    localizer.set_metric(ncf->function()).set_method(calc->param.localize_method());

    MolecularOrbitals<double, 3> mo(nemo, calc->aeps, {}, calc->aocc, calc->aset);
    Tensor<double> UT = localizer.compute_localization_matrix(world, mo, randomize);

    vecfuncT localnemo = transform(world, nemo, transpose(UT));
    truncate(world, localnemo);
    normalize(localnemo, R);
    return localnemo;
}

std::shared_ptr<Fock<double, 3>> Nemo::make_fock_operator() const {
    MADNESS_CHECK(param.spin_restricted());
    const int ispin = 0;

    std::shared_ptr<Fock<double,3> > fock(new Fock<double,3>(world));
    Coulomb<double,3> J(world,this);
    fock->add_operator("J",std::make_shared<Coulomb<double,3> >(J));
    fock->add_operator("V",std::make_shared<Nuclear<double,3> >(world,this));
    fock->add_operator("T",std::make_shared<Kinetic<double,3> >(world));
    if (calc->xc.hf_exchange_coefficient()>0.0) {
        Exchange<double,3> K=Exchange<double,3>(world,this,ispin).set_symmetric(false);
        fock->add_operator("K",{-1.0,std::make_shared<Exchange<double,3>>(K)});
    }
    if (calc->xc.is_dft()) {
        XCOperator<double,3> xcoperator(world,this,ispin);
        real_function_3d xc_pot = xcoperator.make_xc_potential();

        // compute the asymptotic correction of exchange-correlation potential
        if(do_ac()) {
            std::cout << "Computing asymtotic correction!\n";
            double charge = double(molecule().total_nuclear_charge())-param.charge();
            real_function_3d scaledJ = -1.0/charge*J.potential()*(1.0-calc->xc.hf_exchange_coefficient());
            xc_pot = ac.apply(xc_pot, scaledJ);
        }
        LocalPotentialOperator<double,3> xcpot(world);
        xcpot.set_potential(xc_pot);
        xcpot.set_info("xc potential");
        if (do_ac()) xcpot.set_info("xc potential with ac");
        fock->add_operator("Vxc",std::make_shared<LocalPotentialOperator<double,3>>(xcpot));
    }


    // compute the solvent (PCM) contribution to the potential
    if (do_pcm()) {
        const real_function_3d vpcm = pcm.compute_pcm_potential(J.potential());
        LocalPotentialOperator<double,3> pcmpot(world);
        pcmpot.set_potential(vpcm);
        pcmpot.set_info("pcm potential");
        fock->add_operator("Vpcm",std::make_shared<LocalPotentialOperator<double,3>>(pcmpot));
    }
    return fock;
}

/// compute the Fock matrix from scratch
tensorT Nemo::compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const {
	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT Jnemo, Knemo, xcnemo, pcmnemo, JKVpsi, Unemo;

    vecfuncT R2nemo=mul(world,R_square,nemo);
    truncate(world,R2nemo);

    // compute potentials the Fock matrix: J - K + Vnuc
	compute_nemo_potentials(nemo, Jnemo, Knemo, xcnemo, pcmnemo, Unemo);

//    vecfuncT JKUpsi=add(world, sub(world, Jnemo, Knemo), Unemo);
    vecfuncT JKUpsi=Unemo+Jnemo-Knemo;
    if (do_pcm()) JKUpsi+=pcmnemo;
    if (calc->xc.is_dft()) JKUpsi+=xcnemo;
    tensorT fock=matrix_inner(world,R2nemo,JKUpsi,false);   // not symmetric actually
    Kinetic<double,3> T(world);
    fock+=T(R2nemo,nemo);
    JKUpsi.clear();

	return 0.5*(fock+transpose(fock));
}

/// solve the HF equations
double Nemo::solve(const SCFProtocol& proto) {

	// guess has already been performed in value()
	vecfuncT& nemo = calc->amo;
	//long nmo = nemo.size();

	// NOTE that nemos are somewhat sensitive to sparse operations (why??)
	// Therefore set all tolerance thresholds to zero, also in the mul_sparse

	// apply all potentials (J, K, Vnuc) on the nemos
	vecfuncT Jnemo, Knemo, xcnemo, pcmnemo, Unemo;

	std::vector<double> energies(1,0.0);	// contains the total energy and all its contributions
	double energy=0.0;
	bool converged = false;
	bool localized=param.do_localize();
	real_function_3d density=real_factory_3d(world); 	// for testing convergence

    auto solver= nonlinear_vector_solver<double,3>(world,nemo.size());

	// iterate the residual equations
	for (int iter = 0; iter < param.maxiter(); ++iter) {

	    if (localized) nemo=localize(nemo,proto.dconv,iter==0);
	    std::vector<std::string> str_irreps;
	    if (do_symmetry()) nemo=symmetry_projector(nemo,R_square,str_irreps);
	    if (world.rank()==0 and param.print_level()>9) print("orbital irreps",str_irreps);
	    vecfuncT R2nemo=mul(world,R_square,nemo);
	    truncate(world,R2nemo);
        if (iter==0) solver.initialize(nemo);

		// compute potentials the Fock matrix: J - K + Vnuc
		compute_nemo_potentials(nemo, Jnemo, Knemo, xcnemo, pcmnemo, Unemo);

		// compute the energy
		std::vector<double> oldenergies=energies;
		energies = compute_energy_regularized(nemo, Jnemo, Knemo, Unemo);
        energy=energies[0];

		// compute the fock matrix
        timer t_fock(world,param.print_level()>2);
		vecfuncT Vnemo=Unemo+Jnemo-Knemo;
		if (do_pcm()) Vnemo+=pcmnemo;
        if (calc->xc.is_dft()) Vnemo+=xcnemo;
		tensorT fock=matrix_inner(world,R2nemo,Vnemo,false);   // not symmetric actually
		Kinetic<double,3> T(world);
		fock+=T(R2nemo,nemo);
		t_fock.end("compute fock matrix");


        // Diagonalize the Fock matrix to get the eigenvalues and eigenvectors
        if (not localized) {
            timer t(world,param.print_level()>2);
    		// report the off-diagonal fock matrix elements
            tensorT fock_offdiag=copy(fock);
            for (int i=0; i<fock.dim(0); ++i) fock_offdiag(i,i)=0.0;
            double max_fock_offidag=fock_offdiag.absmax();
            if (world.rank()==0 and param.print_level()>3) print("F max off-diagonal  ",max_fock_offidag);

            // canonicalize the orbitals, rotate subspace and potentials
            tensorT overlap = matrix_inner(world, R2nemo, nemo, true);
            tensorT U=calc->get_fock_transformation(world,overlap,
                    fock,calc->aeps,calc->aocc,FunctionDefaults<3>::get_thresh());

            nemo = transform(world, nemo, U, trantol(), true);
            Vnemo = transform(world, Vnemo, U, trantol(), true);
            // rotate_subspace(world, U, solver, 0, nemo.size());

            truncate(world, nemo);
            normalize(nemo,R);
            t.end("canonicalize orbitals");

        } else {
        	// if localized the orbital energies are the diagonal fock matrix elements
        	for (int i=0; i<calc->aeps.size(); ++i) calc->aeps[i]=fock(i,i);
        }

        timer t_bsh(world,param.print_level()>2);
		BSHApply<double,3> bsh_apply(world);
		bsh_apply.metric=R_square;
		bsh_apply.ret_value=BSHApply<double,3>::update;
		bsh_apply.lo=get_calc()->param.lo();
		bsh_apply.levelshift=param.orbitalshift();
		auto [update,eps_update] =bsh_apply(nemo,fock,Vnemo);
	    auto residual=nemo-update;
		t_bsh.tag("BSH apply");

		const double bsh_norm = norm2(world, residual) / sqrt(nemo.size());

		// vecfuncT nemo_new = truncate(solver.update(nemo, residual));
		vecfuncT nemo_new = truncate(solver.update(update));
        t_bsh.tag("solver.update");
		normalize(nemo_new,R);

		calc->do_step_restriction(world,nemo,nemo_new,"ab spin case");
		orthonormalize(nemo_new,R);
		nemo=nemo_new;

		real_function_3d olddensity=density;
		density=R_square*compute_density(nemo);
		double deltadens=(density-olddensity).norm2();
		converged=check_convergence(energies,oldenergies,bsh_norm,deltadens,param,
				proto.econv,proto.dconv);

		if (param.save()) calc->save_mos(world);
        t_bsh.tag("orbital update");

		if (world.rank() == 0 and param.print_level()>1) {
			printf("finished iteration %2d at time %8.1fs with energy  %12.8f\n",
					iter, wall_time(), energy);
		}

		if (converged) break;
	}

	if (converged) {
		if (world.rank()==0) print("\nIterations converged\n");
        calc->converged_for_thresh=param.econv();
        if (param.save()) calc->save_mos(world);
    } else {
		if (world.rank()==0) print("\nIterations failed\n");
		energy = 0.0;
	}

	return energy;
}

/// given nemos, compute the HF energy using the regularized expressions for T and V
std::vector<double> Nemo::compute_energy_regularized(const vecfuncT& nemo, const vecfuncT& Jnemo,
        const vecfuncT& Knemo, const vecfuncT& Unemo) const {
    timer t(world,param.print_level()>2);

    vecfuncT R2nemo=R_square*nemo;
    truncate(world,R2nemo);

    const tensorT U = inner(world, R2nemo, Unemo);
    const double pe = 2.0 * U.sum();  // closed shell

//    real_function_3d dens=dot(world,nemo,nemo)*R_square;
//    double pe1=2.0*inner(dens,calc->potentialmanager->vnuclear());

    // compute \sum_i <F_i | R^2 T | F_i>
    double ke = 0.0;
    for (int axis = 0; axis < 3; axis++) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        const vecfuncT dnemo = apply(world, D, nemo);
        const vecfuncT dr2nemo = apply(world, D, R2nemo);
        ke += 0.5 * (inner(world, dnemo, dr2nemo)).sum();
    }
    ke *= 2.0; // closed shell


    double ke0=compute_kinetic_energy(nemo);
//    double ke1=compute_kinetic_energy1(nemo);
//    double ke2=compute_kinetic_energy2(nemo);

    const double J = inner(world, R2nemo, Jnemo).sum();
    double K = inner(world, R2nemo, Knemo).sum();

    int ispin=0;
    double exc=0.0;
    if (calc->xc.is_dft()) {
        XCOperator<double,3> xcoperator(world,this,ispin);
        exc=xcoperator.compute_xc_energy();
    }

    double pcm_energy=0.0;
    if (do_pcm()) pcm_energy=pcm.compute_pcm_energy();
    const double nucrep = calc->molecule.nuclear_repulsion_energy();

    double energy = ke + J - K + exc + pe + nucrep + pcm_energy;

    if (world.rank() == 0 and param.print_level()>2) {
        printf("\n  nuclear and kinetic %16.8f\n", ke + pe);
        printf("         kinetic only %16.8f\n",  ke0);
//        printf("\n  kinetic only  %16.8f\n",  ke2);
//        printf("\n  kinetic plain %16.8f\n",  ke3);
//        printf("\n  nuclear only  %16.8f\n",  pe1);
//        printf("\n  nuclear and kinetic each  %16.8f\n",  pe1+ke1);
        printf("              coulomb %16.8f\n", J);
        if (is_dft()) printf(" exchange-correlation %16.8f\n", exc);
        if (calc->xc.hf_exchange_coefficient()!=0.0) printf("       exact exchange %16.8f\n", -K);
        if (do_pcm()) printf("   polarization (PCM) %16.8f\n", pcm_energy);
        printf("    nuclear-repulsion %16.8f\n", nucrep);
        printf("   regularized energy %16.8f\n", energy);
    }
    t.end( "compute energy");

    return std::vector<double>{energy,ke0,ke+pe,J,exc,-K,pcm_energy,nucrep};
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
void Nemo::compute_nemo_potentials(const vecfuncT& nemo,
		vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& xcnemo, vecfuncT& pcmnemo,
		vecfuncT& Unemo) const {

    {
        timer t(world,param.print_level()>2);
        real_function_3d vcoul;
        int ispin = 0;
        auto taskq = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(world, world.size()));
        Coulomb<double, 3> J = Coulomb<double, 3>(world, this).set_taskq(taskq);
        {
            t.tag("initialize Coulomb operator");

            // compute the density and the coulomb potential
            Jnemo = J(nemo);

            // compute the exchange potential
            int ispin = 0;
            Knemo = zero_functions_compressed<double, 3>(world, nemo.size());
            if (calc->xc.hf_exchange_coefficient() > 0.0) {
                Exchange<double, 3> K = Exchange<double, 3>(world, this, ispin).set_symmetric(true).set_taskq(taskq);
	            K.set_algorithm(Exchange<double,3>::Algorithm::multiworld_efficient_row);
                Knemo = K(nemo);
            }
            t.tag("initialize K operator");
            taskq->set_printlevel(param.print_level());
            if (param.print_level()>9) taskq->print_taskq();
            taskq->run_all();

            t.tag("compute Knemo");
            scale(world, Knemo, calc->xc.hf_exchange_coefficient());
            truncate(world, Knemo);
        }

        // compute the exchange-correlation potential
        if (calc->xc.is_dft()) {
            XCOperator<double, 3> xcoperator(world, this, ispin);
            //double exc = 0.0;
            //if (ispin == 0) exc = xcoperator.compute_xc_energy();
            real_function_3d xc_pot = xcoperator.make_xc_potential();

            // compute the asymptotic correction of exchange-correlation potential
            if (do_ac()) {
                std::cout << "Computing asymtotic correction!\n";
                double charge = double(molecule().total_nuclear_charge()) - param.charge();
                real_function_3d scaledJ = -1.0 / charge * J.potential() * (1.0 - calc->xc.hf_exchange_coefficient());
                xc_pot = ac.apply(xc_pot, scaledJ);
            }

            xcnemo=truncate(xc_pot*nemo);
            t.tag("compute XCnemo");
        }


        // compute the solvent (PCM) contribution to the potential
        if (do_pcm()) {
            const real_function_3d vpcm = pcm.compute_pcm_potential(J.potential());
            pcmnemo = vpcm * nemo;
            double size = get_size(world, pcmnemo);
            t.tag("compute PCMnemo " + stringify(size));
        }

        Nuclear<double, 3> Unuc(world, this->ncf);
        Unemo = Unuc(nemo);
        double size1 = get_size(world, Unemo);
        t.tag("compute Unemo " + stringify(size1));
    }
    world.gop.fence();

}



/// return the Coulomb potential
real_function_3d Nemo::get_coulomb_potential(const vecfuncT& psi) const {
	MADNESS_ASSERT(param.spin_restricted());
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
    const Nuclear<double,3> U_op(world,this->ncf);
    const Nuclear<double,3> V_op(world,this->get_calc().get());

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

    const Nuclear<double,3> U_op(world,this->ncf);
    const Nuclear<double,3> V_op(world,this->get_calc().get());

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

    const vecfuncT& nemo=calc->amo;

    // the pseudo-density made up of the square of the nemo orbitals
    functionT rhonemo = make_density(calc->aocc, nemo).scale(2.0);
    rhonemo=rhonemo.refine();

    Tensor<double> grad=NemoBase::compute_gradient(rhonemo,calc->molecule);


    if (world.rank() == 0) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print(
              "  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(
              " ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (size_t i = 0; i < calc->molecule.natom(); ++i) {
            const Atom& atom = calc->molecule.get_atom(i);
            printf(" %3s %3d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", atomic_number_to_symbol(atom.atomic_number).c_str(),int(i),
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

    timer time_hessian(world,param.print_level()>2);

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
    if (param.get<bool>("purify_hessian")) hessian=purify_hessian(hessian);

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

                    unsigned int ZA=calc->molecule.get_atomic_number(iatom);
                    unsigned int ZB=calc->molecule.get_atomic_number(jatom);
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
    QProjector<double,3> Q(R2nemo,nemo);

    DNuclear<double,3> Dunuc(world,this,iatom,iaxis);
    vecfuncT Vpsi2b=Dunuc(nemo);
    truncate(world,Vpsi2b);

    // construct some intermediates
    NuclearCorrelationFactor::RX_functor rxr_func(ncf.get(),iatom,iaxis,2);
    const real_function_3d RXR=real_factory_3d(world).functor(rxr_func).truncate_on_project();

    // part of the Coulomb operator with the derivative of the NCF
    // J <- \int dr' 1/|r-r'| \sum_i R^XR F_iF_i
    Coulomb<double,3> Jconst(world);
    Jconst.potential()=Jconst.compute_potential(2.0*RXR*rhonemo);        // factor 2 for cphf
    vecfuncT Jconstnemo=Jconst(nemo);
    truncate(world,Jconstnemo);

    // part of the exchange operator with the derivative of the NCF
    // K <- \sum_k |F_k> \int dr' 1/|r-r'| 2R^XR F_k F_i
    // there is no constant term for DFT, since the potentials are not
    // linear in the density
    vecfuncT Kconstnemo=zero_functions_compressed<double,3>(world,nmo);
    if (not is_dft()) {
        Exchange<double,3> Kconst(world,param.lo());
        vecfuncT kbra=2.0*RXR*nemo;
        truncate(world,kbra);
        Kconst.set_bra_and_ket(kbra, nemo);
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
    QProjector<double,3> Q(R2nemo,nemo);

    // construct quantities that are independent of xi

    // construct the BSH operator
    tensorT eps(nmo);
    for (int i = 0; i < nmo; ++i) eps(i) = fock(i, i);
    std::vector<poperatorT> bsh = calc->make_bsh_operators(world, eps,param);
    for (poperatorT& b : bsh) b->destructive()=true;    // make it memory efficient

    // derivative of the (regularized) nuclear potential

    // construct the KAIN solver
    typedef vector_function_allocator<double, 3> allocT;
    typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
    allocT alloc(world, nemo.size());
    solverT solver(allocT(world, nemo.size()));
    solver.set_maxsub(5);

    // construct unperturbed operators
    const Coulomb<double,3> J(world,this);
    const Exchange<double,3> K=Exchange<double,3>(world,this,0);
    const XCOperator<double,3> xc(world, xc_data, not param.spin_restricted(), arho, arho);
    const Nuclear<double,3> V(world,this);

    Tensor<double> h_diff(3l);
    for (int iter=0; iter<10; ++iter) {


        // make the rhs
        vecfuncT Vpsi1=V(xi) + J(xi);
        if (is_dft()) {
            Vpsi1+=(xc(xi));
        } else {
            Vpsi1-=(K(xi));
        }

        truncate(world,Vpsi1);


        // construct perturbed operators
        Coulomb<double,3> Jp(world);
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
            Exchange<double,3> Kp1(world,param.lo());
            Kp1.set_bra_and_ket(R2nemo, xi_complete).set_symmetric(true);
            vecfuncT R2xi=mul(world,R_square,xi_complete);
            truncate(world,R2xi);
            Exchange<double,3> Kp2(world,param.lo());
            Kp2.set_bra_and_ket(R2xi, nemo);
            Kp=truncate(Kp1(nemo) + Kp2(nemo));
        }
        vecfuncT Vpsi2=truncate(Jp(nemo)-Kp+rhsconst);
        Vpsi2=Q(Vpsi2);
        truncate(world,Vpsi2);

        vecfuncT Vpsi=truncate(Vpsi1+Vpsi2);
        Vpsi1.clear();
        Vpsi2.clear();


        // add the coupling elements in case of localized orbitals
        if (get_calc()->param.do_localize()) {
            Tensor<double> fcopy=copy(fock);
            for (int i = 0; i < nmo; ++i) fcopy(i, i) -= eps(i);
            vecfuncT fnemo= transform(world, xi, fcopy, trantol(), true);
            gaxpy(world, 1.0, Vpsi, -1.0, fnemo);
            truncate(Vpsi);
        }

        // apply the BSH operator on the wave function
        scale(world,Vpsi,-2.0);
        vecfuncT tmp = apply(world, bsh,Vpsi);
        Vpsi.clear();
        truncate(world, tmp);

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


    timer t1(world,param.print_level()>2);
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
    std::vector<poperatorT> bsh = calc->make_bsh_operators(world, eps,param);
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

    SCFProtocol preiterations(world,param);
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
    SCFProtocol p(world,param);
    p.start_prec=p.end_prec;
    p.initialize();

    for ( ; not p.finished(); ++p) {
        set_protocol(p.current_prec);

        if (world.rank()==0) {
            printf("\nstarting CPHF equations at time %8.1fs \n",wall_time());
            print("solving CPHF with the density functional",param.xc());
        }

        // double loop over all nuclear displacements
        for (int i=0, iatom=0; iatom<natom; ++iatom) {
            for (int iaxis=0; iaxis<3; ++iaxis, ++i) {
                if (xi[i].size()>0) {
                    for (real_function_3d& xij : xi[i]) xij.set_thresh(p.current_prec);
                }
                xi[i]=solve_cphf(iatom,iaxis,fock,xi[i],rhsconst[i],
                        incomplete_hessian,parallel[i],p,param.xc());
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
    Tensor<double> D=MolecularOptimizer::projector_external_dof(molecule(),{"tx","ty","tz","rx","ry","rz"});
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
