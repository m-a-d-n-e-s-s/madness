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

/** \file interior_dirichlet.cc
    \brief This file demonstrates solving a problem with interior (embedded)
    Dirichlet conditions.

    The signed_distance_functions shapes (mra/sdf_shape.h and
    mra/sdf_shape_3D.h) are used to specify a sphere of radius 1 centered
    at the origin.  A Dirichlet condition (dir_cond) is imposed on this sphere.

    After constructing the mask and the imposed boundary condition, the
    following routine is used to solve the equation (see the Lowengrub paper).

    Suppose \f$\varphi\f$ is the mask function (1 on the inside, 0 on the
    outside, blurry on the boundary), \f$u\f$ is the desired function, \f$d\f$
    is the imposed Dirichlet condition on the boundary, \f$f\f$ is the
    inhomogeneity, \f$\mathcal{L}\f$ is the differential operator, and \f$G\f$
    is its free-space Green's function.

    The DE is \f$ \mathcal{L} u = f\f$ in the domain, and
        \f[ \mathcal{L}u - b(\varphi) (u - d) = \varphi f, \f]
    where \f$b(\varphi) = 36 \varepsilon^{-3} \varphi^2 (1 - \varphi)^2\f$ and
    \f$\varepsilon\f$ is the thickness of the surface layer.

    Applying the Green's function:
        \f[ u - G*( b(\varphi) u) == G*(\varphi f) - G*( b(\varphi) d). \f]
    Thus, solving this problem involves a linear solve, as provided by GMRES.

    To run this code, at least two command-line arguments are needed:
      1) the width of the surface layer (0.1 at the largest, 0.01 seems ok)
      2) 0 for the Poisson equation, 1 for the (BS) Helmholtz equation
      3) if Helmholtz, the frequency
*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include <mra/sdf_shape_3D.h>

/// \brief The width of the sphere's surface, \f$\varepsilon\f$.
double eps;

/// \brief \f$1/\varepsilon\f$
double inveps;

/// \brief Is this a Helmholtz problem (\f$k\neq0\f$?)
bool is_helmholtz;

/// \brief \f$k\f$ for a Helmholtz problem
double helmholtz_k;

using namespace madness;

/// \brief Threshold for MADNESS.
static const double thresh = 1e-4;

/// \brief One-tenth of threshold
static const double thresh1 = thresh*0.1;

/** \brief The Dirichlet condition on the sphere.

    @param pt The point at which to evaluate; \f$|pt|=1\f$.
    @return The Dirichlet condition. */
static double dir_cond(const coord_3d &pt) {
	// Y_1^0
	const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

	return pt[2] / r;
}

/** \brief The inhomogeneity -- right-hand-side of the PDE.

    \param pt The point at which to evaluate.
    \return The inhomgeneity. */
static double f_rhs(const coord_3d &pt) {
	return 0.0;
}

/** \brief The exact solution, for comparison (only valid inside the sphere).

     \param pt The point at which to evaluate.
     \return The exact solution. */
static double exact_sol(const coord_3d & pt) {
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

/** \brief Gives the surface function, \f$B(\varphi(x)) = \varphi^2 (1-\varphi)^2\f$.

     This is a unary op, given the \f$\varphi(x)\f$ function, it computes
     \f$b\f$ pointwise.

     This function should probably be deprecated with the new sdf_shape
     library that accounts for surfaces... */
template <typename T>
inline static void b_phi(const Key<3> &key, Tensor<T> & t) {
	UNARY_OPTIMIZED_ITERATOR(T, t,
		*_p0 = (*_p0) * (*_p0) * (1.0-(*_p0)) * (1.0-(*_p0)));
}

/// \brief The operator needed for solving for \f$u\f$ with GMRES
class DirichletCondIntOp : public Operator<real_function_3d> {
	protected:
		/// \brief The Green's function
		const SeparatedConvolution<double,3> &G;
		/// \brief The surface function, \f$b\f$
		const real_function_3d &b;

      /** \brief Applies the operator to \c invec

		    \param[in] invec The input vector
		    \param[out] outvec The action of the Green's function on \c invec */
		void action(const real_function_3d &invec, real_function_3d &outvec)
			const {
                    outvec = invec - G(b*invec);
                    outvec.truncate();
		}

	public:
		DirichletCondIntOp(const SeparatedConvolution<double, 3> &gin,
			const real_function_3d &bin)
			: G(gin), b(bin) {}
};


//*****************************************************************************
// C++ function to solve an embedded Dirichlet problem.
int main(int argc, char **argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        // get the helmholtz_k value from the command-line
        if(argc < 3) {
            std::cerr << "Usage error: ./app_name eps (0 Poisson, 1 Helmholtz) " \
                "[helmholtz_k]" << std::endl;
            error("bad number of arguments");
        }
        
        eps = atof(argv[1]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");
    
        switch(atoi(argv[2])) {
        case 0:
            is_helmholtz = false;
            helmholtz_k = 0.0;
            break;
        case 1:
            is_helmholtz = true;
            if(argc < 4 || (helmholtz_k = atof(argv[3])) == 0.0)
                error("Must specify a helmholtz_k != 0 for a Helmholtz problem");
            break;
        default:
            error("Only 0 (Poisson) and 1 (Helmholtz) are accepted.");
	}
    }
    world.gop.broadcast(eps);
    world.gop.broadcast(is_helmholtz);
    world.gop.broadcast(helmholtz_k);
    inveps = 1.0 / eps;

    
    // Function defaults
    int k = 6;
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    //FunctionDefaults<3>::set_initial_level(4);
    /// the following line can be commented out if memory is not an issue
    FunctionDefaults<3>::set_max_refine_level(6);
    
    // create the forcing function inhomogeneity
    real_function_3d f = real_factory_3d(world).f(f_rhs);
    
    // create the Dirichlet boundary condition, expanded throughout the domain
    real_function_3d d = real_factory_3d(world).f(dir_cond);
    d.truncate();
    
    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    SharedPtr< SignedDFInterface<3> > sphere = SharedPtr<
        SignedDFInterface<3> >(new SDFSphere(1.0, pt));

    // use LLRV domain masking
    SharedPtr< DomainMaskInterface > llrvmask = SharedPtr<
        DomainMaskInterface >(new LLRVDomainMask(eps));

    // make the functor
    SharedPtr< DomainMaskSDFFunctor<3> > functor = SharedPtr<
        DomainMaskSDFFunctor<3> >(new DomainMaskSDFFunctor<3>(llrvmask, sphere));
    
    real_function_3d phi = real_factory_3d(world).functor(functor);

    // THIS FORMULATION OF B IS OLD SCHOOL AND SHOULD BE REPLACED BY
    // THE SURFACE OPTION ON THE FUNCTOR
    real_function_3d b = copy(phi);
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
    real_function_3d usol = phi*f + b*d; // should be -b*d, but b accounts for -G
    real_function_3d uinhomog = G(usol).truncate();
    uinhomog.scale(-1.0); // add the -1 from the Green's function
    uinhomog.truncate();
    world.gop.fence();
    usol.clear();
    
    // solve the linear system
    // make an initial guess
    usol = copy(uinhomog);
    DirichletCondIntOp dcio(G, b);
    FunctionSpace<double, 3> space(world);
    int maxiters = 20;
    double resid_thresh = 1.0e-2;
    double update_thresh = 1.0e-2;
    GMRES(space, dcio, uinhomog, usol, maxiters, resid_thresh, update_thresh,
        true);
    
    real_function_3d exact = real_factory_3d(world).f(exact_sol);
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
