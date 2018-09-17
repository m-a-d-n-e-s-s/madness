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
static const long k = 8;        // wavelet order
static const double thresh = 1e-5; // precision

static const double phi_denom=0.0;
static const bool singlet=false;
static const double B_init=0.1;

//static const coord_3d nuc1={0.0,0.0, R/2.0};
//static const coord_3d nuc2={0.0,0.0,-R/2.0};

static const coord_3d nuc1={ R/2.0,0.0,0.0};
static const coord_3d nuc2={-R/2.0,0.0,0.0};

struct spherical_box : public FunctionFunctorInterface<double,3> {
	const double radius;
	const double height;
	const double tightness;
	spherical_box(const double r, const double h, const double t) :
		radius(r), height(h), tightness(t) {}

	double operator()(const coord_3d& xyz) const {
		double r=xyz.normf();
		double v1=height/(1.0+exp(-tightness*height*(r-radius)));
		return 1.0-v1;
	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};

static double_complex guess(const coord_3d& r) {
	const double r1=(r-nuc1).normf();
	const double r2=(r-nuc2).normf();
    const double x=r[0], y=r[1], z=r[2];
    double rho=(exp(-sqrt(r1*r1+1e-8)) + exp(-sqrt(r2*r2+1e-8)));
    double_complex result2=std::polar(rho,0.0);
//    double_complex result2=std::polar(rho,constants::pi/phi_denom);
    return result2;
}

static double_complex guess_sigma_u(const coord_3d& r) {
	const double r1=(r-nuc1).normf();
	const double r2=(r-nuc2).normf();
    const double x=r[0], y=r[1], z=r[2];
    double rho=(exp(-sqrt(r1*r1+1e-8)) - exp(-sqrt(r2*r2+1e-8)));
//    double_complex result2=std::polar(rho,0.0);
//    double_complex result2=std::polar(rho,constants::pi/phi_denom);
//    return result2;
    return double_complex(rho,0);
}


double monomial_x(const coord_3d& r) {return r[0];}
double monomial_y(const coord_3d& r) {return r[1];}

static double_complex V(const coord_3d& r) {
	const double r1=(r-nuc1).normf();
	const double r2=(r-nuc2).normf();
    const double v1=-1.0/sqrt(r1*r1+1e-8) - 1.0/sqrt(r2*r2+1e-8);
    return double_complex(v1,0);

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

	return (x*x+y*y);
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


double compute_energy(const complex_function_3d& Vpsi, const complex_function_3d& psi,
		const std::string text="") {
	World& world=Vpsi.world();
	double ke=compute_kinetic_energy(psi);
	double_complex pe=(inner(Vpsi,psi));
	print(text,"ke, pe : ",ke,pe, "   ",ke+real(pe));
	return ke+real(pe);
}

void iterate_triplet(World& world, std::vector<complex_function_3d>& Vpsi,
		std::vector<complex_function_3d>& psi, std::vector<double>& eps) {

    print("norm(Vpsi)", norm2s(world,Vpsi));
//    eps[0]=compute_energy(Vpsi[0],psi[0],"  sigma_g");
//    eps[1]=compute_energy(Vpsi[1],psi[1],"  sigma_u");

	print("current eps ",eps);
	real_convolution_3d op0 = BSHOperator3D(world, sqrt(-2.*std::min(-0.05,eps[0])), 0.001, 1e-6);
	real_convolution_3d op1 = BSHOperator3D(world, sqrt(-2.*std::min(-0.05,eps[1])), 0.001, 1e-6);
	scale(world,Vpsi,-2.0);
    truncate(world,Vpsi);
    std::vector<complex_function_3d> tmp(2);
    tmp[0] = op0(Vpsi[0]).truncate();
    tmp[1] = op1(Vpsi[1]).truncate();

    std::vector<double> norms = norm2s(world,tmp);
    print("norm(tmp)",norms);
    print("norm(psi)",norm2s(world,psi));
    std::vector<complex_function_3d> r = tmp-psi;
    psi = copy(world,tmp);
    std::vector<double> invnorm={1.0/norms[0],1.0/norms[1]};
    scale(world,psi,invnorm);

    Tensor<double_complex> in=inner(world,Vpsi,r);
    print("  inner(Vpsi,r)",in);
    eps[0]-=0.5*real(in(0l))/(norms[0]*norms[0]);
    eps[1]-=0.5*real(in(1l))/(norms[1]*norms[1]);

    std::vector<double> rnorms = norm2s(world,r);

    if (world.rank() == 0) {
        print("norm=",norms," eps=",eps," err(psi)=",rnorms);
    }

    for (int i : {0,1}) {
		complex_function_3d Lzpsi=Lz(psi[i]);
		double_complex lz=inner(Lzpsi,psi[i]);
		print("  <psi | L_z | psi> ",i,real(lz));
    }
}


// return \sum_k \phi_k \int dr' 1/|r-r'| phi_k*phi_i
std::vector<complex_function_3d> apply_K(World& world, const std::vector<complex_function_3d> rhs) {
	real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

	std::vector<complex_function_3d> result=
			zero_functions_compressed<double_complex,3>(world,rhs.size());
	for (int i=0; i<result.size(); ++i) {
		for (int k=0; k<rhs.size(); ++k) {
			result[i]+=rhs[k]*op(rhs[i]*conj(rhs[k]));
		}
	}
	return result;
}

void orthonormalize(std::vector<complex_function_3d>& rhs) {
	World& world=rhs[0].world();
	std::vector<double> norms=norm2s(world,rhs);
	rhs[0].scale(1.0/norms[0]);
	double_complex ovlp=inner(rhs[1],rhs[0]);
	rhs[1]-=ovlp*rhs[0];
	double n=rhs[1].norm2();
	rhs[1].scale(1.0/n);
	double_complex ovlp1=inner(rhs[1],rhs[0]);

	print("orthonormalize: ovlp, norm2(1), ovlp(new)",ovlp,n, ovlp1);
	double inorm0=imag(rhs[0]).norm2();
	double inorm1=imag(rhs[1]).norm2();

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
    FunctionDefaults<3>::set_cubic_cell(-30.0, 30.0);

    double B=B_init;

    print("\n\nH2 in a strong magnetic field\n\n");
    print("B = ",B);
    print("nuclear configuration");
    print(nuc1);
    print(nuc2);

    print("guess rotation phi: pi/",phi_denom,"\n");

//    complex_function_3d pp = complex_factory_3d(world).f(p_plus);
//    double n=pp.norm2();
//    pp.scale(1.0/n);
//    complex_function_3d lzp=Lz(pp);
//    double_complex ev1=inner(pp,lzp);
//    print("lz expectation value for the p+ orbital", ev1);

    complex_function_3d Vnuc = complex_factory_3d(world).f(V);

	spherical_box sbox2(15,1,4);
	real_function_3d sbox=real_factory_3d(world).functor(sbox2);
	print("spherical damping box radius: ",sbox2.radius);


    if (singlet) {
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
			potential+=-0.5*B*Lz(psi);
			potential+=0.125*B*B*diamagnetic(world)*psi;
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


    } else { // triplet

    	std::vector<complex_function_3d> psi(2);
    	psi[0]=complex_factory_3d(world).f(guess);
    	psi[1]=complex_factory_3d(world).f(guess_sigma_u);

    	orthonormalize(psi);

		real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

		for (double B : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

			print("computing with B= ",B);
			std::vector<double> eps = {-0.6,-0.3};
			for (int iter=0; iter<20; iter++) {
				print("\n\n ---- Iteration ",iter,"\n\n");
				real_function_3d rho = (abssq(psi[0])+abssq(psi[1])).truncate();
				complex_function_3d J=convert<double,double_complex,3>(op(rho).truncate());
				std::vector<complex_function_3d> potential = Vnuc*psi + J*psi;
				potential[0]+=B*0.5*Lz(psi[0]);					// angular momentum
				potential[1]+=B*0.5*Lz(psi[1]);					// angular momentum
				potential+=0.125*B*B*diamagnetic(world)*psi;	// diamagnetic part
				potential+=0.5*B*psi;								// spin part
				potential-=apply_K(world,psi);

				iterate_triplet(world, potential, psi, eps);
				psi=psi*sbox;									// dampen at the boundaries
				orthonormalize(psi);
				plot_plane(world,real(psi[0]),real(psi[1]),"real_parts");
				save_function(psi,"psi");
				std::vector<real_function_3d> realpsi={real(psi[0]),real(psi[1])};
				save(realpsi[0],"realpsi0");
				save(realpsi[1],"realpsi1");
				std::vector<real_function_3d> abspsi={abssq(psi[0]),abssq(psi[1])};
				save(abspsi[0],"abssqpsi0_B"+stringify(B));
				save(abspsi[1],"abssqpsi1_B"+stringify(B));
			}
			print("final orbital energies for B=",B, eps);
		}

    }

//    real_function_3d r=real(psi);
//    real_function_3d i=imag(psi);
//    plot_plane(world,std::vector<real_function_3d>{r,i},"complex_function");


    world.gop.fence();

    finalize();
    return 0;
}
