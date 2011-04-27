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
  \file helium_mp2.cc
  \brief Solves the Hartree-Fock and MP2 equations for the helium atom
  \defgroup examplehehf Hartree-Fock and MP2 for the helium atom
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/helium_mp2.cc>here</a>.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <iostream>


using namespace madness;

static const double dcut=1.e-5;
static const double shift=0.0;

//static const double L = 32.0;   // box size
//static const long k = 6 ;        // wavelet order
//static const double thresh = 1e-3; // precision
//static const TensorType tt = TT_3D;
//static const long truncate_mode = 0;

//template<typename T>
//static std::string stringify(T arg) {
//	std::ostringstream o;
//	if (!(o << arg))
//		throw std::domain_error("stringify(double)");
//	return o.str();
//}

static double guess(const coord_3d& r) {
    const double x=r[0]+shift, y=r[1]+shift, z=r[2];
    return 6.0*exp(-sqrt(x*x+y*y+z*z+1e-8));
}

static double V(const coord_3d& r) {
    const double x=r[0]+shift, y=r[1]+shift, z=r[2];
    return -1.0/(sqrt(x*x+y*y+z*z+1e-8));
}

static double HO_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return x*x + y*y + z*z;
}

static double HO_6d(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double value=(x1*x1 + y1*y1 + z1*z1) + (x2*x2 + y2*y2 + z2*z2);
    return value;
}


static double gauss_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= x*x + y*y + z*z;
    const double norm=0.712705695388313;
    return norm*exp(-r2);
}

static double gauss_6d(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double r2=(x1*x1 + y1*y1 + z1*z1) + (x2*x2 + y2*y2 + z2*z2);
    const double norm=0.5;
    return norm*exp(-r2);
}


static double HO_vphi_3d(const coord_3d& r) {
    const double v=HO_3d(r);
    const double phi=gauss_3d(r);
    return v*phi;
}


static double Z2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+dcut*dcut));
}



static double V_1(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
	const double y1=r[1], y2=r[4];
	const double z1=r[2], z2=r[5];


	const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
	const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);

	const double value=-2.0/r1 - 2.0/r2;
	return value;


}



static double g12(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    const double value=1.0/r12;
    return value;
}


static double helium_pot(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    const double value=-2.0/r1 - 2.0/r2 + 1.0/r12;
    return value;
}

// return the helium potential times the hylleraas function
static double V_times_phi(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    const double pot=-2.0/r1 - 2.0/r2;// + 1.0/r12;
    const double phi=exp(-1.8*(r1 + r2))*(1.0 + 0.5*r12);
    const double value=pot*phi;

    return value;
}


// according to Ed
static double he_orbital_3d(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);

    const double val=(0.00995312870402086*(31.05166416452748 - 7.405311261369526*r1 +
           r1*r1)*(4.335496673568937 + 0.24243181498262073*r1 +
           r1*r1)) * exp(-1.81607*r1);

    return val;
}

// according to McQuarrie
static double he_orbital_McQuarrie(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);

    const double val=exp(-(27.0/16.0)*r1);

    return val;
}


// according to Ed / McQuarrie
static double he_orbitals(const coord_6d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
    coord_3d r1;
    coord_3d r2;
    r1[0]=r[0];
    r1[1]=r[1];
    r1[2]=r[2];
    r2[0]=r[3];
    r2[1]=r[4];
    r2[2]=r[5];

//    const double val=he_orbital_3d(r1) * he_orbital_3d(r2);
    const double val=he_orbital_McQuarrie(r1) * he_orbital_McQuarrie(r2);

    return val;
}



static double f6d_svd(const coord_6d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

//    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
//    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
//    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);
//    const double value=exp(-1.8*(r1 + r2))*(1.0 + 0.5*r12);
    const double r1 = (x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = (x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= (xx*xx + yy*yy + zz*zz + dcut*dcut);
    const double value=exp(-1.8*(r1 + r2))*(1.0 + 0.5*r12);
//    const double value=exp(-1.8*(r1 + r2));
    return value;
}

static double f6d_sr(const coord_6d& r) {

	// separation for 3-way decomposition (x12 -- y12 -- z12)
	const double x1=r[0], x2=r[1];
    const double y1=r[2], y2=r[3];
    const double z1=r[4], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

	const double dcut=0.01;
    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);
    const double value=exp(-1.8*(r1 + r2))*(1.0 + 0.5*r12);
    return value;
}


static double hylleraas_3term(const coord_6d& r) {


    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];
    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    return -exp(-1.81607*(r1 + r2)) * (
    -1.33083943395992
    -0.388320016632985 * r12
    -0.174093511691879 *  ( r1*r1  + r2*r2  -2 * r1 * r2 )
    );
}


// Hylleraas 3-term minus He orbitals
static double he_correlation(const coord_6d& r) {
	return hylleraas_3term(r) - he_orbitals(r);
}



void iterate(World& world, const real_function_6d& Vpsi, real_function_6d& psi, double& eps) {

    real_convolution_6d op = BSHOperator<6>(world, sqrt(-2*eps), 0.001, 1e-6);

    print("starting convolution");
   	real_function_6d tmp = op(Vpsi).truncate();
   	tmp.scale(-2.0);
   	print("finished convolution");

    double norm = tmp.norm2();
    print("finished norm");
    real_function_6d r = tmp-psi;
    print("finished difference");
    double rnorm = inner(r,r);
    print("finished rnorm");
    double eps_new = eps + inner(r,Vpsi)/(norm*norm);
    print("finished inner(Vpair,r)");
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("eps_new",eps_new);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

void iterate(World& world, const real_function_3d& V, real_function_3d& psi, double& eps) {

	real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
	real_function_3d copy_of_Vpsi = copy(Vpsi);
	real_function_3d copy_of_Vpsi2 = copy(Vpsi);
	real_function_3d copy_of_Vpsi3 = copy(Vpsi);
	real_function_3d copy_of_V = copy(V);
	real_function_3d copy_of_V2 = copy(V);
	real_function_3d copy_of_psi = copy(psi);
	real_function_3d copy_of_psi2 = copy(psi);
	real_function_3d copy_of_psi3 = copy(psi);
	copy_of_psi2.scale(-2.0);

    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d tmp;

    // direct Vpsi
    if (0) {
    	// set the impl of arg to get the structure of the target tree
    	real_function_3d arg=CompositeFactory<double,3,3>(world).ket(copy_of_Vpsi.get_impl())
    			.muster(copy_of_Vpsi3.get_impl());

    	tmp = op(arg).truncate();
    }

    // composite Vpsi
    else {
//    	real_function_3d arg=CompositeFactory<double,3,3>(world).ket(copy_of_psi2.get_impl())
//					.g12(copy_of_V.get_impl()).muster(copy_of_psi.get_impl());
    	real_function_3d arg=CompositeFactory<double,3,3>(world).ket(copy_of_Vpsi2.get_impl())
					.muster(copy_of_psi.get_impl());

    	tmp = op(arg).truncate();
    }

//    // conventional
//    else {
//    	tmp=op(Vpsi).truncate();
//    }

    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("eps_new",eps_new);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

void compute_energy(World& world, const real_function_3d& psi, const real_function_3d& pot, double& ke, double& pe) {


    print("working on the kinetic energy");
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
    	real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    	real_function_3d dpsi = D(psi);
    	kinetic_energy += 0.5*inner(dpsi,dpsi);
    	print("done with axis",axis);
    }
    ke=kinetic_energy;
    if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());

    pe=inner(psi,pot*psi);
    if(world.rank() == 0) printf("\ncompute V at time %.1fs\n\n", wall_time());

    print("kinetic energy:  ", ke);
    print("potential energy:", pe);
    print("total energy:    ", pe+ke);


}



void compute_energy(World& world, const real_function_6d& pair,
		const real_function_3d& pot1, const real_function_3d& pot2, double& ke, double& pe) {

	// compute kinetic energy
	ke=0.0;
	if (1) {
		for (int axis=0; axis<6; axis++) {
			real_derivative_6d D = free_space_derivative<double,6>(world, axis);
			real_function_6d dpsi = D(pair);
			double a=0.5*inner(dpsi,dpsi);
			ke += a;
			print("done with axis",axis, a);
		}
	}
	print("kinetic energy:", ke);

	if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());


	// compute potential energy
	pe=0.0;
	if (0) {
		// doomed copy of pair, to save pair
		real_function_6d copy_of_pair=copy(pair);

		// two-electron interaction potential
		real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-6);

		real_function_6d v11=CompositeFactory<double,6,3>(world)
				.ket(copy_of_pair.get_impl())
				.g12(eri.get_impl())
				.V_for_particle1(copy(pot1).get_impl())
				.V_for_particle2(copy(pot2).get_impl())
				;


		double a=inner(pair,v11);
		print("<phi|V_tot|phi> ", a);
		pe=a;
	}

	if(world.rank() == 0) printf("\npotential at time %.1fs\n\n", wall_time());

}



int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);



    double L = 16;   // box size
    long k = 4 ;        // wavelet order
    double thresh = 1.e-2; // precision
    TensorType tt = TT_2D;
    long truncate_mode = 0;

    if (argc==6) {
		L = atof(argv[1]);   // box size
		k = atoi(argv[2]) ;        // wavelet order
		thresh = atof(argv[3]); // precision
		tt = TensorType(atoi(argv[4]));
		truncate_mode = atoi(argv[5]);
    }

    if(world.rank() == 0) printf("\nstarting at time %.1fs\n\n", wall_time());


    // hydrogen
#if 0

    tt=TT_FULL;

    FunctionDefaults<3>::set_tensor_type(tt);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(truncate_mode);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

    print("polynomial order:  ", FunctionDefaults<3>::get_k());
    print("threshold:         ", FunctionDefaults<3>::get_thresh());
    print("cell size:         ", L);
    print("truncation mode:   ", FunctionDefaults<3>::get_truncate_mode());
    print("tensor type:       ", FunctionDefaults<3>::get_tensor_type());



//    real_function_3d Vnuc = real_factory_3d(world).f(V).truncate_mode(0);
//    print("helium potential ",Vnuc.tree_size());
//    real_function_3d psi  = real_factory_3d(world).f(guess);
    real_function_3d psi  = real_factory_3d(world).f(gauss_3d);
    if(world.rank() == 0) printf("\nguess at time %.1fs\n\n", wall_time());
    print("helium guess tree size    ", psi.tree_size());
    print("helium guess number coeff ", psi.size());

    print("normalizing");
    double normsq=inner(psi,psi);
    psi.scale(1.0/sqrt(normsq));
    print("<psi | psi>", normsq);
    if(world.rank() == 0) printf("\ninner at time %.1fs\n\n", wall_time());
    real_function_3d pot = real_factory_3d(world).f(V);


    double ke=0.0;
    double pe=0.0;

    compute_energy(world,psi,pot,ke,pe);
    double te=ke+pe;

	// compute the potential energy on the fly
    if (1) {
		real_function_3d copy_of_psi=copy(psi);
		real_function_3d copy_of_pot=copy(pot);
		real_function_3d Vpsi2=pot*psi;

		real_function_3d Vpsi=CompositeFactory<double,3,3>(world)
							.ket(Vpsi2.get_impl())
							;
		print("computing <psi | V | psi> on the fly");
		pe=inner(psi,Vpsi);
		print("potential energy, on the fly",pe);
    }


    for (int i=0; i<10; i++) {
    	iterate(world,pot,psi,te);
        compute_energy(world,psi,pot,ke,pe);
        psi.reconstruct();
		trajectory<3> traj(-L/2,L/2,201);
		std::string filename="iteration"+stringify(i);
		plot_along<3>(world,traj,psi,filename);

    }

#endif

    // helium
#if 1

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_tensor_type(TT_2D);


    print("polynomial order:  ", FunctionDefaults<6>::get_k());
    print("threshold:         ", FunctionDefaults<6>::get_thresh());
    print("cell size:         ", L);
    print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
    print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());


    // one orbital at a time
	real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
//	real_function_3d orbital=real_factory_3d(world).f(gauss_3d);
    print("orbital.tree_size()",orbital.tree_size());
    print("orbital.size()     ",orbital.size());
    {
//    	real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
    	double norm=inner(orbital,orbital);
    	print("norm(orbital)",norm);
    	orbital.scale(1.0/sqrt(norm));

    	// compute kinetic energy
    	double kinetic_energy = 0.0;
    	for (int axis=0; axis<3; axis++) {
    		real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    		real_function_3d dpsi = D(orbital);
    		kinetic_energy += 0.5*inner(dpsi,dpsi);
    		print("done with axis",axis);
    	}
    	print("kinetic energy/electron  :", kinetic_energy);

    	// compute potential energy
    	real_function_3d one_el_pot=real_factory_3d(world).f(Z2);
    	double pe=inner(orbital,one_el_pot*orbital);
    	print("potential energy/electron:",pe);
    }

    real_function_6d pair=hartree_product(orbital,orbital);
//    real_function_6d pair=real_factory_6d(world).f(he_orbitals);
//	pair.get_impl()->print_stats();
    print("pair.tree_size()",pair.tree_size());
    print("pair.size()     ",pair.size());

    // normalize pair function
    double norm=inner(pair,pair);
    print("norm(ij_pair)",norm);
    pair.scale(1.0/sqrt(norm));
    if(world.rank() == 0) printf("\npair function at time %.1fs\n\n", wall_time());

    // one-electron potential
    real_function_3d pot1=real_factory_3d(world).f(Z2);
    real_function_3d pot2=real_factory_3d(world).f(Z2);
    if(world.rank() == 0) printf("\nproject at time %.1fs\n\n", wall_time());

    double ke=0.0;
    double pe=0.0;
    compute_energy(world,pair,pot1,pot2,ke,pe);
	double eps=ke+pe;
    return 0;

    // iterate
	for (unsigned int i=0; i<10; i++) {
		// doomed copy of pair, to save pair
		real_function_6d copy_of_pair=copy(pair);
		real_function_6d copy2_of_pair=copy(pair);

		// two-electron interaction potential
		real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-6);

		real_function_6d v11=CompositeFactory<double,6,3>(world)
							.ket(copy_of_pair.get_impl())
							.g12(eri.get_impl())
							.V_for_particle1(pot1.get_impl())
							.V_for_particle2(pot2.get_impl())
							.muster(copy2_of_pair.get_impl())
							;

		iterate(world,v11,pair,eps);
		compute_energy(world,pair,pot1,pot2,ke,pe);
//		pair.get_impl()->print_stats();
	    print("pair.tree_size()",pair.tree_size());
	    print("pair.size()     ",pair.size());

	}

    print("for he orbitals");
    print("total energy   :",ke+pe);
    print("expected energy:",-2.8477);

#endif

    // HO
#if 0
    // one orbital at a time
     {
     	real_function_3d orbital=real_factory_3d(world).f(gauss_3d);
     	double norm=inner(orbital,orbital);
     	print("norm(orbital)",norm);
     	orbital.scale(1.0/sqrt(norm));

     	// compute kinetic energy
     	double kinetic_energy = 0.0;
     	for (int axis=0; axis<3; axis++) {
     		real_derivative_3d D = free_space_derivative<double,3>(world, axis);
     		real_function_3d dpsi = D(orbital);
     		kinetic_energy += 0.5*inner(dpsi,dpsi);
     		print("done with axis",axis);
     	}
     	print("kinetic energy/electron  :", kinetic_energy);

     	// compute potential energy
     	real_function_3d one_el_pot=real_factory_3d(world).f(HO_3d);
     	double pe=inner(orbital,one_el_pot*orbital);
     	print("potential energy/electron:",pe);
     	double pe4=inner(orbital,one_el_pot);
     	double pe5=inner(orbital,orbital);
     	print("<phi|phi> /electron:",pe5);
     	print("<phi|V> /electron:  ",pe4);


     	// compute potential energy on demand
     	real_function_3d v_phi=real_factory_3d(world).f(HO_vphi_3d);

     	double pe3=inner(orbital,v_phi);
     	print("potential energy 2/electron:",pe3);


     	real_function_3d copy_of_orbital=copy(orbital);
     	real_factory_3d comp_factory=real_factory_3d(world).empty();
        std::shared_ptr<CompositeFunctorInterface<double,3,3> >
        	comp(new CompositeFunctorInterface<double,3,3>(comp_factory,
        			copy_of_orbital.get_impl(),
        			one_el_pot.get_impl(),
        			v_phi.get_impl(),
        			v_phi.get_impl()
        			));

     	real_function_3d v_phi_od=real_factory_3d(world).functor(comp).is_on_demand();
     	double pe2=inner(orbital,v_phi_od);
     	print("potential energy od/electron:",pe2);



     }



     // pair function
     real_function_6d pair=real_factory_6d(world).f(gauss_6d);

     // two-electron potential
     real_function_6d pot_r12=real_factory_6d(world).f(HO_6d).is_on_demand();
     real_function_6d pot_ho=real_factory_6d(world).f(HO_6d);

     // one-electron potential
     real_function_3d pot1=real_factory_3d(world).f(HO_3d);
     real_function_3d pot2=real_factory_3d(world).f(HO_3d);

     // normalize pair function
     double norm=inner(pair,pair);
     print("norm(ij_pair)",norm);
     pair.scale(1.0/sqrt(norm));

     // compute kinetic energy
     double kinetic_energy = 0.0;
     for (int axis=0; axis<6; axis++) {
     	real_derivative_6d D = free_space_derivative<double,6>(world, axis);
     	real_function_6d dpsi = D(pair);
     	kinetic_energy += 0.5*inner(dpsi,dpsi);
     	print("done with axis",axis);
     }
     print("kinetic energy:", kinetic_energy);

     real_function_6d copy_of_pair=copy(pair);
     real_factory_6d comp_factory=real_factory_6d(world).empty();
     std::shared_ptr<CompositeFunctorInterface<double,6,3> >
     	comp(new CompositeFunctorInterface<double,6,3>(comp_factory,
     			copy_of_pair.get_impl(),
     			pot_r12.get_impl(),
     			pot1.get_impl(),
     			pot2.get_impl()
     			));

     // this will project comp to the MRA function V_phi
     real_function_6d V_pair=real_factory_6d(world).functor(comp).is_on_demand();
     double potential_energy=inner(pair,V_pair);
     print("<phi|V_ftor>", potential_energy);
     print("for HO orbitals");
     print("total energy   :",kinetic_energy+potential_energy);

     // direct
     double direct=inner(pair,pot_ho*pair);
     print("direct 2e-pot*pair:",direct);


#endif

     // hylleraas for ACS
#if 0

     FunctionDefaults<3>::set_k(k);
     FunctionDefaults<3>::set_thresh(thresh);
     FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

     FunctionDefaults<6>::set_k(k);
     FunctionDefaults<6>::set_thresh(thresh);
     FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
     FunctionDefaults<6>::set_tensor_type(TT_2D);


     print("polynomial order:  ", FunctionDefaults<6>::get_k());
     print("threshold:         ", FunctionDefaults<6>::get_thresh());
     print("cell size:         ", L);
     print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
     print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());

	 // plot circle with radius 0.5, with electron 2 at (0.0, 0.5, 0.0)
     if (0) {
	     coord_3d el2(0.0);
	     const double phi=1.0;
	     const double radius=0.5;
	     el2[0]=radius * sin(phi);
	     el2[1]=radius * cos(phi);
 		 std::string filename="hylleraas_plotfile_k"+stringify(k)+"_eps"+stringify(thresh);
		 trajectory<6> traj(0.5,el2,201);
		 plot_along<6>(world,traj,hylleraas_3term,filename);
		 print("plotting done");
     }

     // pair function
//     real_function_6d pair=real_factory_6d(world).f(he_correlation);
     real_function_6d pair=real_factory_6d(world).f(hylleraas_3term);
//     print("Hylleraas correlation part: 3-term minus Ed's HF");
     print("Hylleraas-3term;  tree size",pair.tree_size());
     print("Hylleraas-3term;  size     ",pair.size());

     // normalize pair function
     double norm=inner(pair,pair);
     print("norm(ij_pair)",norm);
     pair.scale(1.0/sqrt(norm));

     print("Hylleraas-3term;  tree size",pair.tree_size());
     print("Hylleraas-3term;  size     ",pair.size());

	 // plot circle with radius 0.5, with electron 2 at (0.0, 0.5, 0.0)
     {
		 coord_3d el2(0.0);
		 const double phi=1.0;
		 const double radius=0.5;
		 el2[0]=radius * sin(phi);
		 el2[1]=radius * cos(phi);
		 std::string filename="plotfile_k"+stringify(k)+"_eps"+stringify(thresh);
		 trajectory<6> traj(0.5,el2,201);
		 plot_along<6>(world,traj,pair,filename);
     }

     // plot xy plane containing the origin
     for (int i=0; i<20; i++) {
    	 coord_6d fix_coord(0.0);
    	 // electron 2:
    	 fix_coord[3]=0.0;
    	 fix_coord[4]=0.1+0.1*i;
    	 fix_coord[5]=0.0;

    	 std::shared_ptr<FunctionImpl<double,6> > pair_impl=pair.get_impl();
    	 Tensor<double> cell(6,2);
    	 cell(Slice(_),0)=-2.0;
    	 cell(Slice(_),1)= 2.0;
    	 std::string filename="plot_plane_d"+stringify(fix_coord[4])+"_k"+stringify(k)+"_eps"+stringify(thresh);
    	 pair_impl->print_plane(filename,"xy",fix_coord,cell);

     }
     if (1) {
		 // plot function thru xy plane containing the origin
		 std::shared_ptr<FunctionImpl<double,6> > pair_impl=pair.get_impl();
		 Tensor<double> cell(6,2);
		 cell(Slice(_),0)=-2.0;
		 cell(Slice(_),1)= 2.0;
		 std::string filename="plot_function_k"+stringify(k)+"_eps"+stringify(thresh);
		 FILE * f = fopen(filename.c_str(), "w");
		 if(!f) MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

         fprintf(f,"\\psset{xunit=4cm}\n");
         fprintf(f,"\\psset{yunit=4cm}\n");
         fprintf(f,"\\begin{pspicture}(0,0)(0,4)\n");
         fprintf(f,"\\pslinewidth=0.05pt\n");


		 for (int i=-100; i<100; i++) {
			 for (int j=-100; j<100; j++) {
				 coord_6d fix_coord(0.0);
				 // electron 2:
				 fix_coord[0]=0.0+i*0.05;
				 fix_coord[1]=0.0+j*0.05;
				 fix_coord[2]=0.0;
				 fix_coord[3]=0.0;
				 fix_coord[4]=0.5;
				 fix_coord[5]=0.0;
				 print("coord",fix_coord);
				 fprintf(f,"%12.8f   %12.8f   %12.8f\n",fix_coord[0],fix_coord[1],pair(fix_coord));
			 }
		 }

         fprintf(f,"\\end{pspicture}\n");
		 fclose(f);

     }

#if 0

     // compute kinetic energy
     double kinetic_energy = 0.0;
     for (int axis=0; axis<6; axis++) {
     	real_derivative_6d D = free_space_derivative<double,6>(world, axis);
     	real_function_6d dpsi = D(pair);
     	kinetic_energy += 0.5*inner(dpsi,dpsi);
     	print("done with axis",axis);
     }
     print("kinetic energy:", kinetic_energy);

     if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());


     // compute potential energy
     double potential_energy=0.0;
     {
 		// doomed copy of pair, to save pair
 		real_function_6d copy_of_pair=copy(pair);

 		// two-electron interaction potential
 		real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

 	    // one-electron potential
 	    real_function_3d pot1=real_factory_3d(world).f(Z2);
 	    real_function_3d pot2=real_factory_3d(world).f(Z2);


 		real_function_6d v11=CompositeFactory<double,6,3>(world)
 				.ket(copy_of_pair.get_impl())
 				.g12(eri.get_impl())
 				.V_for_particle1(pot1.get_impl())
 				.V_for_particle2(pot2.get_impl())
 				;


 		double a=inner(pair,v11);
 		print("<phi|V_tot|phi> ", a);
 		potential_energy=a;
     }


     print("for Hylleraas 3-term");
     print("total energy   :",kinetic_energy+potential_energy);
     print("expected energy:",-2.902432);
#endif



#endif



    // start the MP2 bit
#if 0

//    if (world.rank() == 0) {
//        print("            Kinetic energy ", kinetic_energy);
//        print(" Nuclear attraction energy ", nuclear_attraction_energy);
//        print("       Two-electron energy ", two_electron_energy);
//        print("              Total energy ", total_energy);
//        print("                    Virial ", (nuclear_attraction_energy + two_electron_energy) / kinetic_energy);
//    }

#endif
//     print("Hylleraas correlation part: 3-term minus Ed's HF");

     if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
     world.gop.fence();


     finalize();
     return 0;


}


