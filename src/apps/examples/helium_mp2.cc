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
#include <iostream>


using namespace madness;

static const double dcut=1.e-3;
//static const double L = 32.0;   // box size
//static const long k = 6 ;        // wavelet order
//static const double thresh = 1e-3; // precision
//static const TensorType tt = TT_3D;
//static const long truncate_mode = 0;

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return 6.0*exp(-1.0*sqrt(x*x+y*y+z*z+1e-8));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
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


// according to McQuarrie eq. (7.20)
static double he_orbitals(const coord_6d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

//    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double value=exp(-27.0/16.0*(r1 + r2));
    return value;
}

// according to McQuarrie eq. (7.20)
static double he_orbital_3d(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double value=exp(-27.0/16.0*r1);
    return value;
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



void iterate(World& world, real_function_3d& V, real_function_3d& psi, double& eps) {
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
//    Vpsi.ftr2sr();
    real_function_3d tmp = op(Vpsi).truncate();
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);




    double L = 4;   // box size
    long k = 5 ;        // wavelet order
    double thresh = 1.e-3; // precision
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

    tt=TT_3D;

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
    real_function_3d psi  = real_factory_3d(world).f(guess);
    if(world.rank() == 0) printf("\nguess at time %.1fs\n\n", wall_time());
    print("helium guess tree size    ", psi.tree_size());
    print("helium guess number coeff ", psi.size());

    print("normalizing");
    print("<psi | psi>", inner(psi,psi));
    if(world.rank() == 0) printf("\ninner at time %.1fs\n\n", wall_time());
    double normsq=inner(psi,psi);
    psi.scale(1.0/sqrt(normsq));
    print("<psi | psi>", inner(psi,psi));
    if(world.rank() == 0) printf("\ninner3 at time %.1fs\n\n", wall_time());

    print("working on the kinetic energy");
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
    	real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    	real_function_3d dpsi = D(psi);
    	kinetic_energy += 0.5*inner(dpsi,dpsi);
    	print("done with axis",axis);
    }
    print("kinetic energy:", kinetic_energy);
    double ke=kinetic_energy;
    if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());

    print("projecting the potential");
    real_function_3d pot = real_factory_3d(world).f(V);
    if(world.rank() == 0) printf("\nproject V at time %.1fs\n\n", wall_time());
    double pe=inner(psi,pot*psi);
    if(world.rank() == 0) printf("\ncompute V at time %.1fs\n\n", wall_time());

    print("kinetic energy:  ", ke);
    print("potential energy:", pe);
    print("total energy:    ", pe+ke);


#endif

    // random test
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

    // helium
#if 1

    // one orbital at a time
    {
    	real_function_3d orbital=real_factory_3d(world).f(he_orbital_3d);
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


    // pair function
    real_function_6d pair=real_factory_6d(world).f(he_orbitals);

    // one-electron potential
    real_function_3d pot1=real_factory_3d(world).f(Z2).thresh(1.e-5);
    real_function_3d pot2=real_factory_3d(world).f(Z2).thresh(1.e-5);
    if(world.rank() == 0) printf("\nproject at time %.1fs\n\n", wall_time());

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

    if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());

    // compute potential energy
    double potential_energy=0.0;
    {
		// one-electron potential by direct product
		real_function_6d copy_of_pair=copy(pair);
		real_factory_6d comp_factory=real_factory_6d(world).empty();

		std::shared_ptr<CompositeFunctorInterface<double,6,3> >
			comp(new CompositeFunctorInterface<double,6,3>(comp_factory,
					copy_of_pair.get_impl(),
					pot1.get_impl(),
					pot2.get_impl()
					));

		// this will project comp to the MRA function V_phi
		real_function_6d v11=real_factory_6d(world).functor(comp).is_on_demand();
		double a=inner(pair,v11);
		print("<phi|V_tot|phi> ", a);
		potential_energy=a;
    }




    print("for he orbitals");
    print("total energy   :",kinetic_energy+potential_energy);
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
    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();


    finalize();
    return 0;
}
