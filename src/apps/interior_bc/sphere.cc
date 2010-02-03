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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>
#include <linalg/gmres.h>

using namespace madness;

/** \brief The Dirichlet condition */
static double dirichlet_cond(const coord_3d &pt);

/** \brief The exact solution, for comparison */
static double exact_sol(const coord_3d &pt);

/** \brief The operator needed for solving for \f$u\f$ with GMRES */
class DirichletCondIntOp : public Operator<real_function_3d> {
    protected:
        /// \brief The Green's function
        const real_convolution_3d &G;
        /// \brief The surface function (normalized)
        const real_function_3d &b;
        /// \brief The reciprocal square surface width \f$\varepsilon^{-2}\f$
        const double inveps2;

        /** \brief Applies the operator to \c invec

            \note \c G is actually \f$-G\f$.

            \param[in] invec The input vector
            \param[out] outvec The action of the operator on \c invec */
        void action(const real_function_3d &invec, real_function_3d &outvec)
            const {

                outvec = invec - G(b*invec);
                outvec.truncate();
        }

    public:
        DirichletCondIntOp(const real_convolution_3d &gin,
            const real_function_3d &bin, const double eps)
            : G(gin), b(bin), inveps2(1.0 / (eps*eps)) {}
};

int main(int argc, char **argv) {
    double eps;
    int k, maxrefine, maxiter;
    double thresh;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 6) {
            std::cerr << "Usage error: ./app_name k thresh max_refine eps " \
                "maxiter" << std::endl;
            error("bad number of arguments");
        }
        
        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");
        maxrefine = atoi(argv[3]);
        if(maxrefine < 6) error("cheapskate");
        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");
        maxiter = atoi(argv[5]);
        if(maxiter < 1) error("maxiter >= 1");
    }
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(maxrefine);
    world.gop.broadcast(k);
    
    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_max_refine_level(maxrefine);
    
    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    SharedPtr<SignedDFInterface<3> > sphere = SharedPtr<SignedDFInterface<3> >
        (new SDFSphere(1.0, pt));

    SharedPtr<DomainMaskInterface> llrv = SharedPtr<DomainMaskInterface>
        (new LLRVDomainMask(eps));

    SharedPtr<DomainMaskSDFFunctor<3> > functor = SharedPtr<
        DomainMaskSDFFunctor<3> >(new DomainMaskSDFFunctor<3>(llrv, sphere));

    real_function_3d phi = real_factory_3d(world).functor(functor)
                               .truncate_on_project();
    functor->setMaskFunction(DomainMaskSDFFunctor<3>::SURFACE);
    real_function_3d surf = real_factory_3d(world).functor(functor)
                                .truncate_on_project();

    // print out the errors in volume and surface area
    // these are checks of the diffuse domain approximation
    double vol = phi.trace();
    double surfarea = surf.trace();
    if(world.rank() == 0) {
        printf("Error in Volume    = %.3e\n",
            fabs(vol-4.0*constants::pi/3.0));
        printf("Error in Surf Area = %.3e\n",
            fabs(surfarea-4.0*constants::pi));
    }

    // scale the surface function by -eps^-2
    // eps^-2 for the surface penalty (see Dirichlet approx. 2 in Li et al.
    // paper; -1 since MADNESS G is really -G mathematically
    surf.scale(-1.0 / (eps * eps));

    // set up the solution to the auxiliary differential equation
    // get the boundary condition
    real_function_3d d = real_factory_3d(world).f(dirichlet_cond);
    d.truncate();

    // green's function
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // boundary condition
    // should be -surf*d, but surf already has a -1 that needs to be cancelled
    real_function_3d usol = surf*d; 
    // integrate this boundary function
    double intbound = usol.trace();
    if(world.rank() == 0)
        printf("Boundary Integral = %.3e\n", intbound);
    real_function_3d uinhomog = G(usol);
    uinhomog.scale(-1.0);
    uinhomog.truncate();
    world.gop.fence();
    usol.clear();

    // solve the linear system
    // make the initial guess -- make its norm, after one application of dcio,
    //                           close to that of uinhomog
    usol = copy(uinhomog);
    usol.scale(eps*eps);
    DirichletCondIntOp dcio(G, surf, eps);
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-4;
    double update_thresh = 1.0e-4;
    GMRES(space, dcio, uinhomog, usol, maxiter, resid_thresh, update_thresh,
        true);

    // compare to the exact solution
    real_function_3d uexact = real_factory_3d(world).f(exact_sol);
    double error = ((usol - uexact)*phi).norm2();
    if(world.rank() == 0)
        printf("u error %.10e\n", error);

    // file output
    char filename[100];
    sprintf(filename, "sphere.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -1.1;
        plothi[i] = 1.1;
        npts[i] = 71;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
    plotvtk_data(uexact, "uexact", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);

    finalize();
    
    return 0;
}

static double dirichlet_cond(const coord_3d &pt) {
    // use Y_1^0
    const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

    if(r < 1.0e-3)
        return 0.0;
    else
        return pt[2]/r;
}

static double exact_sol(const coord_3d &pt) {
    return pt[2];
}
