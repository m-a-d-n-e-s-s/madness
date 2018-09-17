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
#include <madness/mra/nonlinsol.h>
#include <chem/SCFOperators.h>


using namespace madness;

static const double R = 1.4;    // bond length
static const long k = 8;        // wavelet order
static const double thresh = 1e-5; // precision

static const double phi_denom=1.0;
static const bool singlet=false;

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

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct allocator {
	World& world;
	const int n;

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(World& world, const int nn) :
			world(world), n(nn) {
	}

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(const allocator& other) :
			world(other.world), n(other.n) {
	}

	/// allocate a vector of n empty functions
	std::vector<complex_function_3d> operator()() {
		return zero_functions_compressed<double_complex, 3>(world, n);
	}
};

static double_complex guess(const coord_3d& r) {
	const double r1=(r-nuc1).normf();
	const double r2=(r-nuc2).normf();
    const double x=r[0], y=r[1], z=r[2];
    double rho=(exp(-sqrt(r1*r1+1e-8)) + exp(-sqrt(r2*r2+1e-8)));
//    double_complex result2=std::polar(rho,0.0);
    double_complex result2=std::polar(rho,constants::pi/phi_denom);
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

Tensor<double_complex> compute_fock_matrix(const std::vector<complex_function_3d>& Vpsi,
		const std::vector<complex_function_3d>& psi) {
	World& world=Vpsi[0].world();
	Kinetic<double_complex,3> T(world);
	Tensor<double_complex> Tmat=T(psi,psi);
	Tensor<double_complex> Vmat=matrix_inner(world,Vpsi,psi);
	print("tmat, vmat, fock");
	print(Tmat);
	print(Vmat);
	print(Tmat+Vmat);
	return Tmat+Vmat;
}


double compute_energy(const complex_function_3d& Vpsi, const complex_function_3d& psi,
		const std::string text="") {
	World& world=Vpsi.world();
	double ke=compute_kinetic_energy(psi);
	double_complex pe=(inner(Vpsi,psi));
	print(text,"ke, pe : ",ke,pe, "   ",ke+real(pe));
	return ke+real(pe);
}

complex_function_3d compute_residual(complex_function_3d Vpsi, const complex_function_3d psi,
		double& eps) {
	World& world=Vpsi.world();
	real_convolution_3d op0 = BSHOperator3D(world, sqrt(-2.*std::min(-0.05,eps)), 0.001, 1e-6);
    complex_function_3d tmp = op0(-2.0*Vpsi).truncate();
    complex_function_3d res=psi-tmp;

    // update eps
    double norms=tmp.norm2();
    double_complex in=inner(Vpsi,res);
    double eps1=real(in)/(norms*norms);
    eps-=eps1;
    double rnorm=res.norm2();
    print("rnorm, eps, delta(eps)",rnorm, eps, eps1);
    return res;

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


void orthonormalize_symmetric(std::vector<complex_function_3d>& rhs) {
	World& world=rhs[0].world();

	Tensor<double_complex> ovlp=matrix_inner(world,rhs,rhs);
	print("ovlp in symmetric orthogonalization",ovlp);
	Tensor<double> eval;
	Tensor<double_complex> evec;
	syev(ovlp,evec,eval);
	Tensor<double_complex> sqrtsmat(eval.size(),eval.size());
	for (int i=0; i<eval.size(); ++i) sqrtsmat(i,i)=1.0/sqrt(eval(i));
	Tensor<double_complex> evecT=conj_transpose(evec);
	Tensor<double_complex> sqrts=inner(evec,inner(sqrtsmat,evecT,1,0),1,0);
	rhs=transform(world,rhs,sqrts);


}

void orthonormalize_gs(std::vector<complex_function_3d>& rhs) {
	World& world=rhs[0].world();

	std::vector<double> norms=norm2s(world,rhs);
	rhs[0].scale(1.0/norms[0]);
	double_complex ovlp=inner(rhs[0],rhs[1]);
	rhs[1]-=ovlp*rhs[0];
	double n=rhs[1].norm2();
	rhs[1].scale(1.0/n);
	double_complex ovlp1=inner(rhs[1],rhs[0]);
	double inorm0=imag(rhs[0]).norm2();
	double inorm1=imag(rhs[1]).norm2();
	print("imaginary norm in GS orthogonalization ",inorm0, inorm1);
}

void orthonormalize(std::vector<complex_function_3d>& rhs) {
	orthonormalize_symmetric(rhs);
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

    double B_init=0.0;
    double global_shift_user=0.0;
    double damping_radius_user=0.0;

    print("\n\nH2 in a strong magnetic field\n\n");
    print("nuclear configuration");
    print(nuc1);
    print(nuc2);

    // get the command line options
    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);
        if (key=="-B") B_init=std::stod(val);
        if (key=="--shift") global_shift_user=std::stod(val);
        if (key=="--damping_radius") damping_radius_user=std::stod(val);

    }

    // test lz operator
    complex_function_3d pp=complex_factory_3d(world).f(p_plus);
    double_complex a=inner(pp,Lz(pp));
    double_complex n=pp.norm2();
    printf("testing lz with p+ orbital %12.8f\n",real(a)/(n*n));


    print("guess rotation phi: pi/",phi_denom,"\n");
    complex_function_3d Vnuc = complex_factory_3d(world).f(V);

	print("global shift and radius of the spherical damping box are still to be");
	print("determined by the user -- need to automate this");

    if (singlet) {
    	double B=B_init;
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


    	XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator> solver(allocator(world,2));
    	solver.set_maxsub(10);
    	solver.do_print=true;

		real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

		// set global shift such that the diamagnetic potential is zero
		// at the damping radius
		print("computing with B= ",B_init);

		std::vector<double> eps = {-0.6,-0.3};
//		for (double B : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
		for (double B_percent : {10,20, 30, 40, 50, 60, 70, 80, 90, 100}) {
			double B=B_init/100.0*B_percent;
			solver.clear_subspace();


			double global_shift=20.0*sqrt(2)*B;		// contain the harmonic wave function once
//			double global_shift=5.0*B*B;		// contain the harmonic wave function once
			if (global_shift_user!=0.0) global_shift=global_shift_user;
			print("global shift      ",global_shift);

			double damping_radius=B>0.0 ? sqrt(global_shift/(B*B)) : 100;
			if (damping_radius_user!=0.0) damping_radius=damping_radius_user;
			spherical_box sbox2(damping_radius,1,4); // set radius such that the potential is negative everywhere
			real_function_3d sbox=real_factory_3d(world).functor(sbox2);
			print("spherical damping box radius for pot: ",sbox2.radius);

			spherical_box sbox3(10,1,4); // set radius such that the potential is negative everywhere
			real_function_3d sboxpsi=real_factory_3d(world).functor(sbox3);
			print("spherical damping box radius for psi: ",sbox3.radius);


			int iter_end = (std::abs(B-B_init)<1.e-10) ? 25 : 15;
			for (int iter=0; iter<iter_end; iter++) {
				print("\n\n ---- Iteration ",iter," with B = ", B,"\n\n");
				real_function_3d rho = (abssq(psi[0])+abssq(psi[1])).truncate();
				complex_function_3d J=convert<double,double_complex,3>(op(rho).truncate());
				std::vector<complex_function_3d> potential = Vnuc*psi;
				potential+= J*psi;

				std::vector<complex_function_3d> Vlz(2);

				Vlz[0]=B*0.5*Lz(psi[0]);					// angular momentum
				Vlz[1]=B*0.5*Lz(psi[1]);					// angular momentum
				Tensor<double_complex> lzmat=matrix_inner(world,psi,Vlz);
				print("< Lz > ");
				print(lzmat);
				double inorm0=imag(psi[0]).norm2();
				double inorm1=imag(psi[1]).norm2();
				print("imaginary norm ",inorm0, inorm1);


				potential+=Vlz;

				std::vector<complex_function_3d> diapsi=0.125*B*B*diamagnetic(world)*psi;
				std::string name="potential_raw";
				plot_line(world,name.c_str(),200,coord_3d{-10,0,0},coord_3d{10,0,0},real(diapsi[1]));
				diapsi=diapsi*sbox;						// dampen at the boundaries
				name="potential_boxed";
				plot_line(world,name.c_str(),200,coord_3d{-10,0,0},coord_3d{10,0,0},real(diapsi[1]));


				potential+=diapsi;								// diamagnetic part
				potential+=0.5*B*psi;							// spin part
				potential-=global_shift*psi;							// global shift
				potential-=apply_K(world,psi);


				eps[0]=compute_energy(potential[0],psi[0],"energy(0)");
				eps[1]=compute_energy(potential[1],psi[1],"energy(1)");

				Tensor<double_complex> Fock=compute_fock_matrix(potential,psi);

				std::vector<complex_function_3d> res(2);
				res[0]=compute_residual(potential[0], psi[0], eps[0]);
				res[1]=compute_residual(potential[1], psi[1], eps[1]);
				psi=solver.update(psi,res,0.01,3);
				printf("current eps at B= %4.2f %12.8f %12.8f   --   %12.8f\n ",B,eps[0]+global_shift,eps[1]+global_shift,eps[0]-eps[1]);

				psi=psi*sboxpsi;
				orthonormalize(psi);
				std::string filename="B"+stringify(B);
				plot_plane(world,real(psi[0]),real(psi[1]),"realpsi_"+filename);
				save_function(psi,"psi_"+filename);
				std::vector<real_function_3d> abspsi={abssq(psi[0]),abssq(psi[1])};
				save(abspsi[0],"absqpsi0_"+filename);
				save(abspsi[1],"absqpsi1_"+filename);
			}
			printf("final orbital energies for B= %4.2f %12.8f %12.8f %12.8f\n",B, eps[0]+global_shift,eps[1]+global_shift,eps[0]-eps[1]);
		}

    }



    world.gop.fence();

    finalize();
    return 0;
}
