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

/** \file embedded_dirichlet.cc
    \brief Provides test problems for examining the convergence of
           embedded (Dirichlet) boundary conditions.

    \note Full details of the mathematics of this routine can be found in
          M.G. Reuter et al., Comput. Phys. Commun. 183, pp. 1-7 (2012).
          This code can generate the data in Figures 1, 2, and 4 of that
          article.

    The auxiliary PDE being solved is
    \f[ \nabla^2 u - p(\varepsilon) S (u-g) = \varphi f, \f]
    where
       - \f$u\f$ is the solution function
       - \f$\varepsilon\f$ is the thickness of the boundary layer
       - \f$p(\varepsilon)\f$ is the penalty prefactor, \f$2/\varepsilon\f$
         seems to work well.
       - \f$S\f$ is the surface function
       - \f$g\f$ is the Dirichlet condition to be enforced on the surface
       - \f$\varphi\f$ is the domain mask (1 inside, 0 outside, blurry on the
         border)
       - \f$f\f$ is the inhomogeneity.

    The available test problems are
       -# A sphere of radius \f$R\f$ with \f$g = Y_0^0\f$, homogeneous
          (ConstantSphere)
       -# A sphere of radius \f$R\f$ with \f$g = Y_1^0\f$, homogeneous
          (CosineSphere)
       -# A sphere of radius \f$R\f$ with \f$g = Y_2^0\f$, homogeneous
          (Y20Sphere)
       -# A sphere of radius \f$R\f$ with \f$g = Y_0^0\f$, inhomogeneous
          \f$ f = 1 \f$ (InhomoConstantSphere)

    This program allows testing of various parameters,
       -# The surface thickness
       -# The penalty prefactor
       -# The type of domain masking (LLRV or Gaussian)
       -# The curvature / shape of the domain
       .
    for their effect on convergence of the solution. */

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/tensor/gmres.h>
#include <madness/external/muParser/muParser.h>
#include "test_problems.h"

using namespace madness;

int main(int argc, char **argv) {
    double eps, penalty_prefact;
    int k, prob;
    double thresh, radius;
    Mask mask;

    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    // the structure for the problem
    std::shared_ptr<EmbeddedDirichlet > functor;

    if (world.rank() == 0) {
        if(argc != 7) {
            std::cerr << "Usage error: ./app_name k thresh prob eps penalty" \
                " mask" << std::endl;
            std::cerr << "    Where prob = 1 for constant sphere,\n" \
                         "                 2 for cosine theta sphere,\n" \
                         "                 3 for Y20 sphere\n" \
                         "                 4 for inhomogeneous const. sphere" \
                         "\n" << std::endl;
            std::cerr << "    Where mask = 1 for LLRV, 2 for Gaussian\n"
                << std::endl;
            std::cerr << "    Where penalty is the penalty prefactor, " \
                "specified as a function\n    of eps, e.g., 2/eps" << std::endl;
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        prob = atoi(argv[3]);
        if(prob < 1 || prob > 4) error("bad problem number");

        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        mu::Parser parser;
        try {
            parser.DefineVar("eps", &eps);
            parser.SetExpr(std::string(argv[5]));
            penalty_prefact = parser.Eval();
        }
        catch(mu::Parser::exception_type &e) {
            error(e.GetMsg().c_str());
        }
        if(penalty_prefact <= 0.0) error("penalty prefactor must be positive");

        switch(atoi(argv[6])) {
        case 1:
            mask = LLRV;
            break;
        case 2:
            mask = Gaussian;
            break;
        default:
            error("unknown domain mask type, should be 1 or 2");
            break;
        }

        radius = 1.0;
    }
    world.gop.broadcast(prob);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);
    world.gop.broadcast(penalty_prefact);
    world.gop.broadcast(radius);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_on_project(true);

    // do the final problem setup
    switch(prob) {
    case 1:
        functor.reset(new ConstantSphere(k,
                      thresh, eps, std::string(argv[5]), penalty_prefact,
                      radius, mask));
        break;
    case 2:
        functor.reset(new CosineSphere(k, thresh,
                      eps, std::string(argv[5]), penalty_prefact, radius,
                      mask));
        break;
    case 3:
        functor.reset(new Y20Sphere(k, thresh,
                      eps, std::string(argv[5]), penalty_prefact, radius,
                      mask));
        break;
    case 4:
        functor.reset(new InhomoConstantSphere(k,
                      thresh, eps, std::string(argv[5]), penalty_prefact,
                      radius, mask));
        break;
    default:
        error("shouldn't be here");
        break;
    }

    if(world.rank() == 0) {
        // print out the arguments
        functor->printout();
    }

    // project the surface function
    real_function_3d surf;
    functor->fop = SURFACE;

    // do some initial load balancing, if high k or low thresh
    if(k > 6 || thresh < 1.0e-4) {
        if(world.rank() == 0) {
            printf("Projecting the surface function (to low order)\n");
            fflush(stdout);
        }

        surf = real_factory_3d(world).k(6).thresh(1.0e-4).functor(functor);

        if(world.rank() == 0) {
            printf("Performing load balancing\n");
            fflush(stdout);
        }
        functor->load_balance(world, surf);
        surf.clear();
    }

    // reproject the surface function to the requested threshold & k
    if(world.rank() == 0) {
        printf("Projecting the surface function to requested order\n");
        fflush(stdout);
    }
    surf = real_factory_3d(world).functor(functor);

    // project the domain mask
    real_function_3d phi;
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    functor->fop = DOMAIN_MASK;
    phi = real_factory_3d(world).functor(functor);

    // print out the errors in surface area and volume
    // these are checks of the diffuse domain approximation
    double surf_integral, anals;
    double vol_integral, analv;

    surf_integral = surf.trace();
    anals = functor->SurfaceIntegral();
    vol_integral = phi.trace();
    analv = functor->VolumeIntegral();
    if(world.rank() == 0) {
        printf("Error in Surface Integral: %.6e\n",
            fabs(surf_integral/penalty_prefact - anals));
        printf("Error in Volume Integral: %.6e\n",
            fabs(vol_integral - analv));
    }

    // green's function
    // note that this is really -G...
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);
    //G.broaden();

    // project the r.h.s. function (phi*f - penalty*S*g)
    // and then convolute with G
    real_function_3d usol, rhs;
    if(world.rank() == 0) {
        printf("Projecting the r.h.s. function\n");
        fflush(stdout);
    }

    functor->fop = DIRICHLET_RHS;
    usol = real_factory_3d(world).functor(functor);
    rhs = G(usol);
    rhs.truncate();
    usol.clear();

    // load balance using the domain mask, the surface function, and the rhs
    if(world.rank() == 0){
        printf("Load Balancing\n");
        fflush(stdout);
    }

    LoadBalanceDeux<3> lb(world);
    lb.add_tree(phi, DirichletLBCost<3>(1.0, 1.0));
    lb.add_tree(surf, DirichletLBCost<3>(1.0, 1.0));
    lb.add_tree(rhs, DirichletLBCost<3>(1.0, 1.0));
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,
        false));

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to rhs in magnitude for
    //     starting in GMRES
    DirichletCondIntOp dcio(G, surf);
    usol = copy(rhs);
    usol.scale(1.0 / penalty_prefact);
    usol.compress();

    // make the operators and prepare GMRES
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    GMRES(space, dcio, rhs, usol, maxiter, resid_thresh, update_thresh,
              true);

    // compare to the exact solution
    real_function_3d uexact, uerror;
    double error, ratio, exactval;
    error = 0.0;
    ratio = 0.0;
    exactval = 0.0;

    std::vector<Vector<double, 3> > check_pts;

    functor->fop = EXACT;
    uexact = real_factory_3d(world).functor(functor);
    uerror = (usol - uexact)*phi; // only use interior solution
    error = uerror.norm2();

    if(world.rank() == 0) {
        printf("\nu interior error: %.10e\n", error);
        fflush(stdout);
    }

    // check the points prescribed by the problem
    check_pts = functor->check_pts();
    for(std::vector<Vector<double, 3> >::iterator iter =
        check_pts.begin(); iter != check_pts.end(); ++iter) {

        ratio = usol(*iter);
        exactval = functor->ExactSol(*iter);

        if(world.rank() == 0) {
            printf("u/uexact ratio at (%.2f, %.2f, %.2f): %.10e\n",
                (*iter)[0], (*iter)[1], (*iter)[2], ratio/exactval);
        }
    }

    // uncomment these lines for various plots
    // make a line plot along the positive z axis
    if(world.rank() == 0)
        printf("\n\n");
    double zmin = 0.0;
    double zmax = 2.0;
    int nz = 201;
    /*double zmin = 1.0 - 10.0*eps;
    double zmax = 1.0 + 10.0*eps;
    int nz = 51;*/
    coord_3d pt;
    pt[0] = pt[1] = 0.0;
    double dz = (zmax - zmin) / (nz - 1);
    for(int i = 0; i < nz; ++i) {
        pt[2] = zmin + i * dz;
        double uval = usol(pt);

        if(world.rank() == 0) {
            printf("%.4e %.4e\n", pt[2], uval);
        }
    }

    // print out the solution function
    /*char filename[100];
    sprintf(filename, "spheresol.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -2.0;
        plothi[i] = 2.0;
        npts[i] = 101;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);*/

    finalize();

    return 0;
}
