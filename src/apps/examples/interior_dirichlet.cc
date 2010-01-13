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
/// This file demonstrates solving a problem with interior (embedded) Dirichlet
/// conditions.
/// 
/// The SDF_SHAPE library is used to specify a sphere of radius 1 centered
/// at the origin.  A Dirichlet condition (dir_cond) is imposed on this sphere.
///
/// After constructing the mask and the imposed boundary condition, the
/// following routine is used to solve the equation (see the Lowengrub paper).
///
/// Suppose phi is the mask function (1 on the inside, 0 on the outside, blurry
/// on the boundary), u is the desired function, d is the imposed Dirichlet
/// condition on the boundary, f is the inhomogeneity, L is the differential
/// operator, and G is its free-space Green's function.
///
/// The DE is Lu == f in the domain, and
///
///     Lu - b(phi) (u - d) == phi f,
///
/// where b(phi) = 36 eps**(-3) phi**2 (1 - phi)**2 and eps is the thickness
/// of the surface layer.
///
/// Applying the Green's function:
///
///     u - G*( b(phi) u) == G*(phi f) - G*( b(phi) d)
///
/// Thus, solving this problem involves a linear solve, as provided by GMRES.
///
/// To run this code, at least two command-line arguments are needed:
///   1) the width of the surface layer (0.1 at the largest, 0.01 seems ok)
///   2) 0 for the Poisson equation, 1 for the (BS) Helmholtz equation
///   3) if Helmholtz, the frequency

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include <mra/sdf_shape_3D.h>

double eps;
double inveps;
bool is_helmholtz;
double helmholtz_k;

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef std::vector<functionT> vecFuncT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Tensor<double> tensorT;

const double L = 2.0;

static const double thresh = 1e-4;
static const double thresh1 = thresh*0.1;

/// the Dirichlet condition on the sphere
static double dir_cond(const coordT3d &pt) {
	// Y_1^0
	const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

	return pt[2] / r;
}

/// the inhomogeneity -- right-hand-side of the pde
static double f_rhs(const coordT3d &pt) {
	return 0.0;
}

// the exact solution, for comparison (only true inside the sphere)
static double exact_sol(const coordT3d & pt) {
	const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

	if(r < 1.0e-3)
		return 0.0;

	if(is_helmholtz) {
		return (sinh(helmholtz_k*r) - helmholtz_k*r*cosh(helmholtz_k*r)) /
			(sinh(helmholtz_k) - helmholtz_k*cosh(helmholtz_k)) * pt[2] /
			(r*r*r);
	}
	else {
		return pt[2];
	}
}

// gives the B(phi(x)) function == phi^2 (1-phi)^2
// this is a unary op, given the phi function, computes b pointwise
template <typename T>
inline static void b_phi(const Key<3> &key, Tensor<T> & t) {
	UNARY_OPTIMIZED_ITERATOR(T, t,
		*_p0 = (*_p0) * (*_p0) * (1.0-(*_p0)) * (1.0-(*_p0)));
}

/// the operator needed for solving for u with GMRES
class DirichletCondIntOp : public Operator<functionT> {
	protected:
		// the Green's function
		const SeparatedConvolution<double,3> &G;
		// the surface
		const functionT &b;

		void action(const functionT &invec, functionT &outvec) const {
                    outvec = invec - G(b*invec);
                    outvec.truncate();
		}

	public:
		DirichletCondIntOp(const SeparatedConvolution<double, 3> &gin,
			const functionT &bin)
			: G(gin), b(bin) {}
};


//*****************************************************************************
// C++ function to solve an embedded Dirichlet problem.
int main(int argc, char **argv) {
	coordT3d pt, axis;

	// get the helmholtz_k value from the command-line
	if(argc < 3) {
		std::cerr << "Usage error: ./app_name eps (0 Poisson, 1 Helmholtz) " \
			"[helmholtz_k]" << std::endl;
		return -1;
	}

	eps = atof(argv[1]);
	if(eps <= 0.0) {
		std::cerr << "eps must be positive, and hopefully small." << std::endl;
		return -1;
	}
	inveps = 1.0 / eps;

	switch(atoi(argv[2])) {
		case 0:
			is_helmholtz = false;
			helmholtz_k = 0.0;
			break;
		case 1:
			is_helmholtz = true;
			if(argc < 4 || (helmholtz_k = atof(argv[3])) == 0.0) {
				std::cerr << "Must specify a helmholtz_k != 0 for a Helmholtz problem"
					<< std::endl;
				return -1;
			}
			break;
		default:
			std::cerr << "Only 0 (Poisson) and 1 (Helmholtz) are accepted." << std::endl;
			return -1;
	}

        initialize(argc,argv);
	World world(MPI::COMM_WORLD);
	startup(world,argc,argv);

	// Function defaults
	int k = 6;
	FunctionDefaults<3>::set_k(k);
	FunctionDefaults<3>::set_cubic_cell(-L, L);
	FunctionDefaults<3>::set_thresh(thresh);
	//FunctionDefaults<3>::set_initial_level(4);
	/// the following line can be commented out if memory is not an issue
	FunctionDefaults<3>::set_max_refine_level(6);

	// create the forcing function inhomogeneity
	functionT f = factoryT(world).f(f_rhs);

	// create the Dirichlet boundary condition, expanded throughout the domain
	functionT d = factoryT(world).f(dir_cond);
	d.truncate();

	// create the domain mask, phi, and the surface function, b
	pt[0] = pt[1] = pt[2] = 0.0;
	functionT phi = factoryT(world).functor(shape_mask(eps, new SDFSphere(1.0, pt)));

	functionT b = copy(phi);
	phi.truncate();
	b.unaryop(&b_phi<double>);

	// apply the scale factor to b
	// scale is (36 / epsilon), which normalizes b, -1 from the Green's function
	b.scale(-36.0 * inveps);
	// add in more scaling (essentially a penalty)
	b.scale(inveps * inveps);
	b.truncate();

	// setup the Green's function
	// NOTE that CoulombOperator essentially makes the BSH w/ k == 0.0,
	// and then rescales by 4 pi.  This is more consistent.
	real_convolution_3d G = BSHOperator<3>(world, helmholtz_k, eps*0.1, thresh);

	// compute the inhomogeneous portion
	functionT usol = phi*f + b*d; // should be -b*d, but b accounts for -G
	functionT uinhomog = G(usol).truncate();
	uinhomog.scale(-1.0); // add the -1 from the Green's function
	uinhomog.truncate();
	world.gop.fence();
	usol.clear();

	// solve the linear system
	// make an initial guess
	usol = copy(uinhomog);
	DirichletCondIntOp dcio(G, b);
	FunctionSpace<double, 3> space;
	int maxiters = 10;
	double solver_thresh = 1.0e-4;
	GMRES(space, dcio, uinhomog, usol, maxiters, solver_thresh, true);

	functionT exact = factoryT(world).f(exact_sol);
	double error = ((usol - exact)*phi).norm2();
	printf("   u error: %.10e\n", error);

	// set up file output
	char filename[100];
	sprintf(filename, "interior.vts");
	Vector<double, 3> plotlo, plothi;
	Vector<long, 3> npts;
	for(int i = 0; i < 3; ++i) {
		plotlo[i] = -1.1;
		plothi[i] = 1.1;
		npts[i] = 71;
	}
	plotvtk_begin(world, filename, plotlo, plothi, npts);
	plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
	plotvtk_data(exact, "exact", world, filename, plotlo, plothi, npts);
	plotvtk_end<3>(world, filename);

        finalize();

	return 0;
}
