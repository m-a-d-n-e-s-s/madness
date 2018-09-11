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
  \file complex_h2.cc
  \brief Solves the Hartree-Fock equations for the hydrogen molecule
  \defgroup examplesh2hf Hartree-Fock equations for the hydrogen molecule
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/functionio.cc>h2.cc</a>.

  The Hartree-Fock wave function is computed for the hydrogen molecule
  in three dimensions without using symmetry.

  Since all of the details except for the nuclear potential are the
  same, please refer to the \ref examplehehf helium atom HF example.

*/


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

static const double phi_denom=1.5;

static double_complex guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    double rho=(exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
//    double_complex result2=std::polar(rho,0.0);
    double_complex result2=std::polar(rho,constants::pi/phi_denom);
    return result2;
}

double monomial_x(const coord_3d& r) {return r[0];}
double monomial_y(const coord_3d& r) {return r[1];}

static double Vz(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
}
static double Vx(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(z*z+y*y+(x-R/2)*(x-R/2)+1e-8)+
           -1.0/sqrt(z*z+y*y+(x+R/2)*(x+R/2)+1e-8);
}


static double_complex V(const coord_3d& r) {
	return double_complex(Vx(r),0);
}

static double_complex p_plus(const coord_3d& xyz) {
	double r=xyz.normf();
	double theta=acos(xyz[2]/r);
	double phi=atan2(xyz[1],xyz[0]);
	return r*exp(-r/2.0)*sin(theta)*exp(double_complex(0.0,1.0)*phi);
}

double compute_kinetic_energy(const complex_function_3d& rhs) {
    double ke = 0.0;
    World& world=rhs.world();
    for (int axis=0; axis<3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex,3>(world, axis);
        complex_function_3d dpsi = D(rhs);
        ke += 0.5*real(inner(dpsi,dpsi));
    }
    return ke;
}

// compute the diamagnetic local potential
real_function_3d diamagnetic(World& world) {
    real_function_3d x=real_factory_3d(world).f(monomial_x);
    real_function_3d y=real_factory_3d(world).f(monomial_y);

	return 0.125*(x*x+y*y);
}

// compute the action of the Lz =i r x del operator on rhs
complex_function_3d Lz(const complex_function_3d& rhs) {
	// the operator in cartesian components as
	// L_z =  - i (x del_y - y del_x)

	World& world=rhs.world();
    complex_derivative_3d Dx = free_space_derivative<double_complex,3>(world, 0);
    complex_derivative_3d Dy = free_space_derivative<double_complex,3>(world, 1);
    real_function_3d x=real_factory_3d(world).f(monomial_x);
    real_function_3d y=real_factory_3d(world).f(monomial_y);

	complex_function_3d delx=Dx(rhs);
	complex_function_3d dely=Dy(rhs);

	complex_function_3d result1=x*dely - y*delx;
	complex_function_3d result=result1.scale(-double_complex(0.0,1.0));
	return result;

}



void iterate(World& world, complex_function_3d& V, complex_function_3d& psi, double& eps) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    complex_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    complex_function_3d tmp = op(Vpsi).truncate();
    double norm = tmp.norm2();
    complex_function_3d r = tmp-psi;
    psi = tmp.scale(1.0/norm);

    double rnorm = r.norm2();
    double_complex in1=inner(Vpsi,r);
    double_complex in2=inner(r,Vpsi);
    double_complex in3=inner(conj(r),Vpsi);
    print("  inner(Vpsi,r)",in1,in2,in3);
    double eps_new = eps - 0.5*real(in1)/(norm*norm);

    double ke=compute_kinetic_energy(psi);
    double pe=real(inner(V*psi,psi));

    double realnorm=real(psi).norm2();
    double imagnorm=imag(psi).norm2();
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("    real, imag norm", realnorm, imagnorm);
        print("     <T>, <V>, <H> ",ke, pe, ke+pe);
    }
    eps = eps_new;
    eps=ke+pe;

    complex_function_3d Lzpsi=Lz(psi);
    double_complex lz=inner(Lzpsi,psi);
    print("  <psi | L_z | psi> ",real(lz));
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);


    print("guess rotation phi: pi/",phi_denom,"\n");

    complex_function_3d pp = complex_factory_3d(world).f(p_plus);
    double n=pp.norm2();
    pp.scale(1.0/n);
    complex_function_3d lzp=Lz(pp);
    double_complex ev1=inner(pp,lzp);
    print("lz expectation value for the p+ orbital", ev1);


    complex_function_3d Vnuc = complex_factory_3d(world).f(V);
    complex_function_3d psi  = complex_factory_3d(world).f(guess);
    psi.truncate();
    double norm=psi.norm2();
    psi.scale(1.0/norm);

    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

    double eps = -0.6;
    for (int iter=0; iter<100; iter++) {
        real_function_3d rho = abssq(psi).truncate();
        complex_function_3d oprho=convert<double,double_complex,3>(op(rho).truncate());
        complex_function_3d potential = Vnuc + oprho;
        potential+=-0.5*Lz(psi);
        potential+=diamagnetic(world)*psi;
        iterate(world, potential, psi, eps);
    }

    double kinetic_energy=2.0*compute_kinetic_energy(psi);

    real_function_3d rho = abssq(psi);
    double two_electron_energy = inner(op(rho),rho);
    double nuclear_attraction_energy = 2.0*inner(real(Vnuc),rho);
    double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy +
        nuclear_attraction_energy + nuclear_repulsion_energy;
    double virial = (nuclear_attraction_energy + two_electron_energy + nuclear_repulsion_energy) / kinetic_energy;
    double_complex lz1=inner(psi,Lz(psi));

    if (world.rank() == 0) {
    	print("");
        print("                        Lz ", lz1);
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
    }


    real_function_3d r=real(psi);
    real_function_3d i=imag(psi);
    plot_plane(world,std::vector<real_function_3d>{r,i},"complex_function");


    world.gop.fence();

    finalize();
    return 0;
}
