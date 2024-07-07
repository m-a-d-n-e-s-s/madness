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
/*!
 \file helium_exact.cc
 \brief Solves the two-particle system exactly
 \defgroup helium_exact Solves the two-particle system exactly
 \ingroup examples


 */

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/lbdeux.h>

#include <iostream>

#include<madness/chem/SCF.h>
#include<madness/chem/nemo.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/electronic_correlation_factor.h>
#include"madness/mra/commandlineparser.h"


// switch the electronic interaction on or off
#undef WITH_G12
#define WITH_G12

static double lo = 1.e-8;
static double bsh_eps = 1.e-8;

using namespace madness;

typedef std::shared_ptr<NuclearCorrelationFactor> ncf_ptr;

// apply the mixed double commutator
real_function_6d apply_U_mix(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf, const ncf_ptr& ncf) {

	real_function_6d v = real_factory_6d(world);

	// for screening when the multiplication is performed
	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo,
			bsh_eps);
	op_mod.modified() = true;

	// the part with the derivative operators: U1
	for (int axis = 0; axis < 3; ++axis) {
		if (world.rank() == 0)
			print("Umix working on axis: ", axis);

		// particles 1 and 2
		// note ncf->U1 returns -S'/S  => swap signs
		const real_function_3d U1_axis_ncf1 = copy(ncf->U1(axis)).scale(-1.0);
		const real_function_3d U1_axis_ncf2 = copy((ncf->U1(axis)));

		// is \vec r12/r12
		real_function_6d U1_axis_ecf = ecf.U1(axis);

		real_function_6d x = CompositeFactory<double, 6, 3>(world)
				.g12(U1_axis_ecf).ket(copy(psi));
		x.fill_tree(op_mod);

		real_function_6d x2 = CompositeFactory<double, 6, 3>(world)
					.V_for_particle1(U1_axis_ncf1).V_for_particle2(U1_axis_ncf2)
					.ket(x);
		x2.fill_tree(op_mod);

		v += x2;
		v.truncate();
	}
	return v;
}

// apply the U1 potential of the nuclear correlation factor
real_function_6d apply_U_ncf(World& world, const real_function_6d& psi,
		const ncf_ptr& ncf) {

	real_function_6d result = real_factory_6d(world);

	// for screening when the multiplication is performed
	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo,
			bsh_eps);
	op_mod.modified() = true;

	// the part with the derivative operators: U1
	for (int axis = 0; axis < 6; ++axis) {
		real_derivative_6d D = free_space_derivative<double, 6>(world, axis);
		const real_function_6d Drhs = D(psi).truncate();

		// note integer arithmetic
		if (world.rank() == 0)
			print("axis, axis^%3, axis/3+1", axis, axis % 3, axis / 3 + 1);
		const real_function_3d U1_axis = ncf->U1(axis % 3);

		real_function_6d x;
		if (axis / 3 + 1 == 1) {
			x = CompositeFactory<double, 6, 3>(world).ket(Drhs).V_for_particle1(
					copy(U1_axis));
		} else if (axis / 3 + 1 == 2) {
			x = CompositeFactory<double, 6, 3>(world).ket(Drhs).V_for_particle2(
					copy(U1_axis));
		}
		x.fill_tree(op_mod);
		result += x;
		result.truncate().reduce_rank();
	}

	real_function_3d U2 = ncf->U2();
	real_function_6d r2 =
			CompositeFactory<double, 6, 3>(world).ket(copy(psi)).V_for_particle1(
					copy(U2)).V_for_particle2(copy(U2));
	r2.fill_tree(op_mod);
	result = (result + r2).truncate();

	return result;
}


// apply the U1 potential of the nuclear correlation factor
real_function_6d apply_V(World& world, const real_function_6d& psi, const SCF& calc) {

	real_function_3d vnuc=copy(calc.potentialmanager->vnuclear());
	// for screening when the multiplication is performed
	const double eps = -0.5;
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo,
			bsh_eps);
	op_mod.modified() = true;

	real_function_6d r2 =
			CompositeFactory<double, 6, 3>(world).ket(copy(psi)).V_for_particle1(
					copy(vnuc)).V_for_particle2(copy(vnuc));
	r2.fill_tree(op_mod);

	return r2;
}


// apply the U1 potential of the electronic correlation factor
real_function_6d apply_U_ecf(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf) {

	const double eps = -2.9;
	return ecf.apply_U(psi, eps);
}

/// reconstruct the full wave function

/// Psi from the regularized wave function and the correlation factors
real_function_6d reconstruct_psi(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf, const ncf_ptr& ncf) {
	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world,sqrt(-2*eps),lo,bsh_eps);
	op_mod.modified() = true;

	// reconstruct the full wave function
	real_function_6d R1psi=CompositeFactory<double, 6, 3>(world)
							.V_for_particle1(copy(ncf->function()))
							.ket(copy(psi));
	R1psi.fill_tree(op_mod);

	real_function_6d R12psi=CompositeFactory<double, 6, 3>(world)
							.V_for_particle2(copy(ncf->function()))
							.ket(copy(R1psi));
	R12psi.fill_tree(op_mod);


#ifdef WITH_G12
	real_function_6d Rfpsi = CompositeFactory<double, 6, 3>(world)
				.g12(ecf.function()).ket(copy(R12psi));
	Rfpsi.fill_tree(op_mod);
#else
	real_function_6d Rfpsi=copy(R12psi);
#endif

	double Rfnorm=inner(Rfpsi,Rfpsi);
	if (world.rank()==0) print("< Rf psi | Rf psi> ",Rfnorm);
	return Rfpsi;
}


/// compute R_12^2f_12^2|psi>
real_function_6d compute_R2f2_psi(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf, const ncf_ptr& ncf) {
	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world,sqrt(-2*eps),lo,bsh_eps);
	op_mod.modified() = true;

	// reconstruct the full wave function
	real_function_6d R1psi=CompositeFactory<double, 6, 3>(world)
							.V_for_particle1(copy(ncf->square()))
							.ket(copy(psi));
	R1psi.fill_tree(op_mod);

	real_function_6d Rsq_psi=CompositeFactory<double, 6, 3>(world)
							.V_for_particle2(copy(ncf->square()))
							.ket(copy(R1psi));
	Rsq_psi.fill_tree(op_mod);


#ifdef WITH_G12
	real_function_6d R2f2psi = CompositeFactory<double, 6, 3>(world)
				.g12(ecf.square()).ket(copy(Rsq_psi));
	R2f2psi.fill_tree(op_mod);
#else
	real_function_6d R2f2psi=copy(Rsq_psi);
#endif

	double Rfnorm=inner(R2f2psi,psi);
	if (world.rank()==0) print("< R2f2 psi | psi> ",Rfnorm);
	return R2f2psi;
}

/// compute the energy using the regularization

/// @return the energy expectation value <E> = < psi | R2f2 (T+U) | psi >
double compute_energy_with_U(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf, const ncf_ptr& ncf, const real_function_6d& twoVpsi) {

	// this is the bra function <psi | R^2 f^2
	real_function_6d R2f2psi=compute_R2f2_psi(world,psi,ecf,ncf);

	// the norm of the real wave function (=denominator in the energy expression)
	const double norm=inner(R2f2psi,psi);

	// the pseudo-kinetic energy
	double ke=0;
	for (int i=0; i<6; ++i) {
		real_derivative_6d D = free_space_derivative<double, 6>(world, i);
		const real_function_6d Drhs = D(psi).truncate();
		const real_function_6d Dlhs = D(R2f2psi).truncate();
		const double ke_axis=inner(Dlhs,Drhs);
		ke+=0.5*ke_axis/norm;
		if (world.rank()==0) print("<psi | R2f2 T(i)  | psi>/<Psi | Psi> ",0.5*ke_axis);
	}
	if (world.rank()==0) print("<psi | R2f2 T | psi>/<Psi | Psi> ",ke);

	// the pseudo-potential energy (note factor -2 on Vpsi!)
	const double pe=-0.5*inner(R2f2psi,twoVpsi)/norm;
	if (world.rank()==0) print("<psi | R2f2 U | psi>/<Psi | Psi> ",pe);

	const double energy=ke+pe;
	if (world.rank()==0) print("<psi | R2f2 (T+U) | psi>/<Psi | Psi> ",energy);
	return energy;
}

/// compute the energy using the reconstructed wave function

/// @return the energy expectation value <E> = < Psi | H | Psi >
double compute_energy(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf, const ncf_ptr& ncf, const real_function_3d& vnuc) {

	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world,sqrt(-2*eps),lo,bsh_eps);
	op_mod.modified() = true;


	real_function_6d Rfpsi=reconstruct_psi(world,psi,ecf,ncf);
	real_function_6d psi2=(psi*psi).truncate();

#ifdef WITH_G12
	real_function_6d psif2 = CompositeFactory<double, 6, 3>(world)
				.g12(ecf.square()).ket(copy(psi2));
	psif2.fill_tree(op_mod);
#else
	real_function_6d psif2=copy(psi2);
#endif

	real_function_6d R2 = CompositeFactory<double, 6, 3>(world)
				.particle1(copy(ncf->square()))
				.particle2(copy(ncf->square()));

	double norm=inner(psif2,R2);
	if (world.rank()==0) print("<psi | psi> ",norm);

	{
		real_function_6d R2f2 = CompositeFactory<double, 6, 3>(world)
#ifdef WITH_G12
					.g12(ecf.square())
#endif
					.particle1(copy(ncf->square()))
					.particle2(copy(ncf->square()));
		double norm1=inner(psi2,R2f2);
		if (world.rank()==0) print("<psi | psi> ",norm1);

	}

	real_function_6d eri=TwoElectronFactory(world).dcut(bsh_eps);
	real_function_6d R2V = CompositeFactory<double, 6, 3>(world)
				.g12(eri).particle1(copy(ncf->square()))
				.particle2(copy(ncf->square()));

	double pe=inner(psif2,R2V);
	if (world.rank()==0) print("<psi | V_el | psi>/<psi | psi> ",pe/norm);


	const real_function_3d R2Vnuc=vnuc*ncf->square();
	real_function_6d R2Vnuc6_1  = CompositeFactory<double, 6, 3>(world)
				.particle1(copy(ncf->square()))
				.particle2(copy(R2Vnuc));
	double pe1=inner(psif2,R2Vnuc6_1);
	real_function_6d R2Vnuc6_2  = CompositeFactory<double, 6, 3>(world)
				.particle1(copy(R2Vnuc))
				.particle2(copy(ncf->square()));
	double pe2=inner(psif2,R2Vnuc6_2);

	if (world.rank()==0) print("<psi | V_nuc | psi>/<psi | psi> ",(pe1+pe2)/norm);

	double ke=0;
	for (int i=0; i<6; ++i) {
		real_derivative_6d D = free_space_derivative<double, 6>(world, i);
		const real_function_6d Drhs = D(Rfpsi).truncate();
		const double n=Drhs.norm2();
		ke+=0.5*n*n;
		if (world.rank()==0) print("<psi | T(i)  | psi>/<psi | psi> ",0.5*n*n/norm);
	}

	if (world.rank()==0) print("<psi | T     | psi>/<psi | psi> ",ke/norm);

	real_function_6d VRfpsi=CompositeFactory<double, 6, 3>(world)
							.V_for_particle1(copy(vnuc))
							.V_for_particle2(copy(vnuc))
							.ket(copy(Rfpsi));
	double pe_nuc=inner(VRfpsi,Rfpsi);
	if (world.rank()==0) print("<psi | V_nuc2| psi>/<psi | psi> ",pe_nuc/norm);

	return (pe1 + pe2 + pe + ke)/norm;
}

void save(World& world, const real_function_6d& f, std::string filename) {
	archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, filename.c_str(), 1);
	ar & f;
}

void load(World& world, real_function_6d& f, std::string filename) {
	archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, filename.c_str(), 1);
	ar & f;
}

/// the electronic U potential is U -1/r12 = f^-1 [T, f] = f^-1 T f - T
void test_U_el(World& world, const real_function_6d& psi,
		const CorrelationFactor2& ecf) {

	// expectation value of U
	real_function_6d Vpsi = apply_U_ecf(world, psi, ecf);
	double u=inner(psi,Vpsi);
	if (world.rank()==0) print("<psi | U | psi> ",u);


	real_function_6d eri=TwoElectronFactory(world).dcut(bsh_eps);
	real_function_6d gpsi = CompositeFactory<double, 6, 3>(world)
				.g12(eri).ket(copy(psi));
	double g=inner(psi,gpsi);
	if (world.rank()==0) print("<psi | 1/r12 | psi> ",g);


	// expectation value of the first term rhs
	const double eps = -2.9;
	real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo,
			bsh_eps);
	op_mod.modified() = true;
	real_function_6d fpsi = CompositeFactory<double, 6, 3>(world)
				.g12(ecf.function()).ket(copy(psi));
	fpsi.fill_tree(op_mod);

	real_function_6d finvpsi = CompositeFactory<double, 6, 3>(world)
				.g12(ecf.inverse()).ket(copy(psi));
	finvpsi.fill_tree(op_mod);

	double ftf=0;
	for (int i=0; i<6; ++i) {
		real_derivative_6d D = free_space_derivative<double, 6>(world, i);
		const real_function_6d Dbra = D(fpsi).truncate();
		const real_function_6d Dket = D(finvpsi).truncate();
		const double ke1=0.5*inner(Dbra,Dket);
		ftf+=ke1;
		if (world.rank()==0) print("<psi |f-1 T(i) f| psi> ",ke1);
	}

	if (world.rank()==0) print("<psi | f-1 T f | psi> ",ftf);


	// expectation value of the second term rhs
	double ke=0;
	for (int i=0; i<6; ++i) {
		real_derivative_6d D = free_space_derivative<double, 6>(world, i);
		const real_function_6d Drhs = D(psi).truncate();
		const double n=Drhs.norm2();
		ke+=0.5*n*n;
		if (world.rank()==0) print("<psi | T(i)  | psi> ",0.5*n*n);
	}

	if (world.rank()==0) print("<psi | T     | psi> ",ke);

	if (world.rank()==0) {
		print(" U - f-1 T f - T ",(u - g) - (ftf - ke));
	}

}


int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	startup(world, argc, argv,true);
	std::cout.precision(6);

	if (world.rank()==0) {
          print("           git source description ...", info::git_source_description());
	}

	if (world.rank() == 0) {
		print("main() in helium_exact.cc compiled at ", __TIME__, " on ",
				__DATE__);
	}

	// read out input and stuff
	commandlineparser parser(argc,argv);
    parser.set_keyval("geometry","he"); // it's always going to be helium
	Nemo nemo(world, parser);
	auto calc=nemo.get_calc();

	TensorType tt = TT_2D;
	FunctionDefaults<6>::set_tensor_type(tt);
	FunctionDefaults<6>::set_thresh(FunctionDefaults<3>::get_thresh());
	FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
    FunctionDefaults<3>::set_cubic_cell(-calc->param.L(), calc->param.L());
    FunctionDefaults<6>::set_cubic_cell(-calc->param.L(), calc->param.L());
    calc->set_protocol<6>(world,FunctionDefaults<3>::get_thresh());


    // get command line parameters
    bool do_energy=false;
    std::string testfilename;
    for(int ii = 1; ii < argc; ii++) {
        const std::string arg=argv[ii];
        if (arg=="energy_only") do_energy=true;
    }

	if (world.rank() == 0) {
		print("helium exact");
		print("polynomial order:    ", FunctionDefaults<6>::get_k());
		print("threshold:           ", FunctionDefaults<6>::get_thresh());
		print("truncation mode:     ", FunctionDefaults<6>::get_truncate_mode());
		print("tensor type:         ", FunctionDefaults<6>::get_tensor_type());
		print("cell size:           ", FunctionDefaults<6>::get_cell());
		print("");
		print("facReduce            ", GenTensor<double>::fac_reduce());
		print("max displacement     ", Displacements<6>::bmax_default());
		print("apply randomize      ",
				FunctionDefaults<6>::get_apply_randomize());
		print("world.size()         ", world.size());
		print("compute energy only  ", do_energy);
		print("");
	}


	// solve the one-particle equations for having a guess for the wave function
	nemo.value();

	// create the electronic and nuclear correlation factors
	ncf_ptr ncf = nemo.ncf;
	CorrelationFactor2 ecf(world);

	// get a guess for the wave function
	double eps = -2.919474;
	real_function_6d psi;

	std::string filename = "he_psi_exact";
	bool exists = archive::ParallelInputArchive<madness::archive::BinaryFstreamInputArchive>::exists(world,
			filename.c_str());
	if (exists) {
		if (world.rank() == 0) print("reading wave function from file");
		archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> iar(world, filename.c_str(), 1);
		iar & psi & eps;
		psi.set_thresh(FunctionDefaults<6>::get_thresh());
		real_function_3d vnuc=copy(calc->potentialmanager->vnuclear());
		compute_energy(world,psi,ecf,ncf,vnuc);

	} else {
		psi = hartree_product(calc->amo[0], calc->amo[0]);
	}
	psi.print_size("guess psi");


	if (do_energy) {
		real_function_6d  Vpsi1;
		load(world,Vpsi1,"Vpsi");
		const double en1=compute_energy_with_U(world,psi,ecf,ncf,Vpsi1);
		Vpsi1.clear();
		if (world.rank()==0) print("asymmetric energy :  ",en1);
		real_function_3d vnuc1=copy(calc->potentialmanager->vnuclear());
		const double en2=compute_energy(world,psi,ecf,ncf,vnuc1);
		if (world.rank()==0) print("symmetric energy  :  ",en2);
	} else {

		test_U_el(world,psi,ecf);

		if (world.rank()==0) {
			print(" total energy  2nd order update      < H > ");
		}


		for (int iter = 0; iter < nemo.get_calc()->param.maxiter(); ++iter) {

			real_function_6d Vpsi = real_factory_6d(world);
	#ifdef WITH_G12
			Vpsi = apply_U_mix(world, psi, ecf, ncf);
			Vpsi.print_size("Umix");

			Vpsi += apply_U_ecf(world, psi, ecf);
			Vpsi.truncate();
			Vpsi.print_size("Umix+U_ecf");
	#endif
			Vpsi += apply_U_ncf(world, psi, ncf);
			Vpsi.truncate();
	//		Vpsi += apply_V(world, psi, *nemo.get_calc().get());
			Vpsi.print_size("Umix+U_ecf+U_ncf");

			Vpsi.scale(-2.0).truncate();
			Vpsi.print_size("Vpsi");


			Vpsi.scale(-2.0);
			real_convolution_6d bsh = BSHOperator<6>(world, sqrt(-2 * eps), lo,
					bsh_eps);
			bsh.destructive() = true;
			save(world,Vpsi,"Vpsi");
			real_function_6d tmp = bsh(Vpsi).truncate();// green is destructive
			tmp.print_size("GV Psi");
			save(world,tmp,"GVpsi");

			// second order update of the energy
			real_function_6d r = tmp - psi;
			double rnorm = r.norm2();
			double norm = tmp.norm2();
			load(world,Vpsi,"Vpsi");
			double eps_new = eps - 0.5 * inner(Vpsi, r) / (norm * norm);
			if (world.rank() == 0) {
				print("norm=", norm, " eps=", eps, " err(psi)=", rnorm,
						" err(eps)=", eps_new - eps);
			}
			Vpsi.clear();
			tmp.scale(1.0 / norm);
			psi = copy(tmp);
			if (iter >1) {
				if (world.rank() == 0)
					print("updating energy");
				eps = eps_new;
			}
			archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, filename.c_str(), 1);
			ar & psi & eps;

			real_function_3d vnuc=copy(calc->potentialmanager->vnuclear());
			double energy=compute_energy(world,psi,ecf,ncf,vnuc);

			if (world.rank()==0) {
				print(" total energy ",eps,energy);
			}

		}
	}

	if (world.rank() == 0)
		printf("\nfinished at time %.1fs\n\n", wall_time());
	world.gop.fence();
	finalize();

	return 0;
}

