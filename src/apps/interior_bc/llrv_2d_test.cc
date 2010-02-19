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

/** \file llrv_2d_test.cc
    \brief Provides a 2-D test problem for examining the convergence of
           embedded boundary conditions.

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

    The test problem is the same as that used in LLRV:
    The unit circle with \f$g = 1/4\f$, inhomogeneous: \f$f=1\f$.

    This program allows testing of various parameters,
       -# The surface thickness
       -# The penalty prefactor
       -# The type of domain masking (LLRV or Gaussian)
       .
    for their effect on convergence of the solution. */

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_2D.h>
#include <linalg/gmres.h>
#include <muParser/muParser.h>

using namespace madness;

enum Mask { LLRV, Gaussian };

enum FunctorOutput {SURFACE, DIRICHLET_SURFACE, EXACT, INHOMO };

/** \brief Produces the surface Dirichlet condition */
class SurfaceProblem : public FunctionFunctorInterface<double, 2> {
    private:
        SurfaceProblem() {}

        SharedPtr<DomainMaskInterface> dmi;
        SharedPtr<SignedDFInterface<2> > sdfi;
        double penalty_prefact;

    public:
        /// which function to use when projecting:
        /// -# the weighted surface (SURFACE)
        /// -# the weighted Dirichlet condition (DIRICHLET_SURFACE)
        /// -# the exact solution (EXACT)
        /// -# the original r.h.s. inhomogeneity w/domain mask (INHOMO)
        FunctorOutput fop;

        SurfaceProblem(SharedPtr<DomainMaskInterface> dmi,
            SharedPtr<SignedDFInterface<2> > sdfi, double penalty_prefact)
            : dmi(dmi), sdfi(sdfi), penalty_prefact(penalty_prefact),
              fop(DIRICHLET_SURFACE)
        {}

        double operator() (const coord_2d &x) const {
            switch(fop) {
            case EXACT:
                return ExactSol(x);
                break;
            case INHOMO:
                return dmi->mask(sdfi->sdf(x)) * Inhomogeneity(x);
                break;
            case DIRICHLET_SURFACE:
                return DirichletCond(x) * dmi->surface(sdfi->sdf(x)) *
                    penalty_prefact;
                break;
            case SURFACE:
                return dmi->surface(sdfi->sdf(x)) * penalty_prefact;
                break;
            default:
                error("bad project_mode");
                return 0.0;
                break;
            }
        }

        double DirichletCond(const coord_2d &x) const {
            return 0.25;
        }

        double Inhomogeneity(const coord_2d &x) const {
            return 1.0;
        }

        double ExactSol(const coord_2d &x) const {
            double r2 = x[0]*x[0] + x[1]*x[1];

            if(r2 <= 1.0)
                return r2 * 0.25;
            else
                return 0.25;
        }
};

/** \brief The operator needed for solving for \f$u\f$ with GMRES */
class DirichletCondIntOp : public Operator<real_function_2d> {
    protected:
        /// \brief The Green's function
        const real_convolution_2d &G;
        /// \brief The surface function, weighted with eps^{-2}
        const real_function_2d &b;

        /** \brief Applies the operator to \c invec

            \note \c G is actually \f$-G\f$.

            \param[in] invec The input vector
            \param[out] outvec The action of the operator on \c invec */
        void action(const real_function_2d &invec, real_function_2d &outvec)
            const {

                outvec = invec + G(b*invec);
                outvec.scale(-1.0);
                outvec.truncate();
        }

    public:
        DirichletCondIntOp(const real_convolution_2d &gin,
            const real_function_2d &bin)
            : G(gin), b(bin) {}
};

int main(int argc, char **argv) {
    double eps, penalty_prefact;
    int i, k;
    double thresh;
    Mask mask;
    char maskname[20];

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 5) {
            std::cerr << "Usage error: ./app_name k thresh eps penalty mask "
                << std::endl;
            std::cerr << "    Where mask = 1 for LLRV, 2 for Gaussian\n"
                << std::endl;
            std::cerr << "    Where penalty is the penalty_prefact, " \
                "specified as a function\n    of eps, i.e. 2/eps" << std::endl;
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        eps = atof(argv[3]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        mu::Parser parser;
        try {
            parser.DefineVar("eps", &eps);
            parser.SetExpr(std::string(argv[4]));
            penalty_prefact = parser.Eval();
        }
        catch(mu::Parser::exception_type &e) {
            error(e.GetMsg().c_str());
        }
        if(penalty_prefact <= 0.0) error("penalty_prefact must be positive");

        switch(atoi(argv[5])) {
        case 1:
            mask = LLRV;
            sprintf(maskname, "LLRV");
            break;
        case 2:
            mask = Gaussian;
            sprintf(maskname, "Gaussian");
            break;
        default:
            error("unknown domain mask type, should be 1 or 2");
            break;
        }

        // print out the arguments
        printf("Wavelet order = %d\nThreshold = %.2e\nLayer Thickness = %.3e" \
            "\nPenalty Prefactor, %s = %.6e\nUsing %s Domain Masking\n\n",
            k, thresh, eps, argv[4], penalty_prefact, maskname);
        fflush(stdout);
    }
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);
    world.gop.broadcast(penalty_prefact);
    
    // Function defaults
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<2>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<2>::set_truncate_on_project(true);

    // calculate some nice initial projection level...
    // should be no lower than 6, but may need to be higher for small eps
    i = ceil(log(4.0 / eps) / log(2.0) - 4);
    if(i < 6)
        i = 6;
    if(world.rank() == 0) {
        printf("Initial projection level = %d\n", i);
        fflush(stdout);
    }
    FunctionDefaults<2>::set_initial_level(i);
    
    // create the domain mask, phi, and the surface function
    coord_2d pt(0.0); // Origin
    SharedPtr<SignedDFInterface<2> > sdf(new SDFCircle(1.0, pt));

    DomainMaskInterface *dmip = NULL;
    if(mask == LLRV)
        dmip = new LLRVDomainMask(eps);
    else if(mask == Gaussian)
        dmip = new GaussianDomainMask(eps);
    SharedPtr<DomainMaskInterface> dmask(dmip);

    // a functor for the domain mask
    SharedPtr<DomainMaskSDFFunctor<2> > phi_functor
        (new DomainMaskSDFFunctor<2>(dmask, sdf));

    // a functor for the domain surface...
    // the SURFACE option of phi_functor could be used here, however,
    // we really want (penalty_prefact*surface) for our calculations, and this
    // functor uses that factor when MADNESS projects the function.
    SharedPtr<SurfaceProblem> surf_functor
        (new SurfaceProblem(dmask, sdf, penalty_prefact));

    // project the surface function
    if(world.rank() == 0) {
        printf("Projecting the surface function\n");
        fflush(stdout);
    }
    surf_functor->fop = SURFACE;
    real_function_2d surf = real_factory_2d(world).functor(surf_functor);

    // phi_functor defaults to the domain mask
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    real_function_2d phi = real_factory_2d(world).functor(phi_functor);

    // print out the errors in volume and surface area
    // these are checks of the diffuse domain approximation
    double vol = phi.trace();
    double surfarea = surf.trace();
    if(world.rank() == 0) {
        printf("Error in Area      = %.3e\n",
            fabs(vol-constants::pi));
        printf("Error in Perimeter = %.3e\n",
            fabs(surfarea/penalty_prefact-2.0*constants::pi));
    }

    // green's function
    // note that this is really -G...
    real_convolution_2d G = BSHOperator<2>(world, 0.0, eps*0.1, thresh);

    // project the inhomogeneity
    if(world.rank() == 0) {
        printf("Projecting the r.h.s. inhomogeneity\n");
        fflush(stdout);
    }
    surf_functor->fop = INHOMO;
    real_function_2d d = real_factory_2d(world).functor(surf_functor);

    // make the Dirichlet condition
    // this will include the penalty prefactor
    if(world.rank() == 0) {
        printf("Projecting the Dirichlet condition\n");
        fflush(stdout);
    }
    surf_functor->fop = DIRICHLET_SURFACE;
    real_function_2d usol;
    real_function_2d uboundary = real_factory_2d(world).functor(surf_functor);
    
    // make a line plot
    /*{
    if(world.rank() == 0)
        printf("\n\n");
    double xmin = 0.0;
    double xmax = 2.0;
    int nx = 201;
    double dx = (xmax - xmin) / (nx - 1);
    for(int i = 0; i < nx; ++i) {
        pt[0] = xmin + i * dx;
        double uval = uboundary(pt);
        double down = d(pt);

        if(world.rank() == 0) {
            printf("%.4e %.4e %4e\n", pt[0], uval, down);
        }
    }
    if(world.rank() == 0)
        printf("\n\n");
    }*/

    // make the r.h.s. of the auxiliary differential equation
    d.compress();
    uboundary.compress();
d.gaxpy(1.0, uboundary, -1.0);
uboundary = d;
//    uboundary.gaxpy(-1.0, d, 1.0);
    //uboundary.scale(-1.0);
    //d.clear();
    
    // make a line plot
    /*{
    if(world.rank() == 0)
        printf("\n\n");
    double xmin = 0.0;
    double xmax = 2.0;
    int nx = 201;
    double dx = (xmax - xmin) / (nx - 1);
    for(int i = 0; i < nx; ++i) {
        pt[0] = xmin + i * dx;
        double uval = uboundary(pt);

        if(world.rank() == 0) {
            printf("%.4e %.4e\n", pt[0], uval);
        }
    }
    if(world.rank() == 0)
        printf("\n\n");
    }*/

    // transform Dirichlet condition into the right-hand side vector
    // rhs = -G*(condition)
    d = G(uboundary);
    d.truncate();

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to d in magnitude for
    //     starting in GMRES
    usol.clear();
    usol = copy(d);
    usol.scale(1.0 / penalty_prefact);
    usol.compress();

    // make the operators and prepare GMRES
    DirichletCondIntOp dcio(G, surf);
    FunctionSpace<double, 2> space(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    GMRES(space, dcio, d, usol, maxiter, resid_thresh, update_thresh, true);

    // compare to the exact solution
    surf_functor->fop = EXACT;
    real_function_2d uexact = real_factory_2d(world).functor(surf_functor);
    real_function_2d uerror = usol - uexact;
    double error = uerror.norm2();
    pt[0] = 0.0;
    pt[1] = 0.1;
    double cons = usol(pt);
    if(world.rank() == 0) {
        printf("\nu error = %.10e\n", error);
        printf("u/uexact ratio at %f = %.10e\n", 0.1,
            cons/surf_functor->ExactSol(pt));
    }
    pt[1] = 2.0;
    cons = usol(pt);
    if(world.rank() == 0) {
        printf("u/uexact ratio at %f = %.10e\n", 2.0,
            cons/surf_functor->ExactSol(pt));
    }
    pt[1] = 1.0;
    cons = usol(pt);
    if(world.rank() == 0) {
        printf("u/uexact ratio at %f = %.10e\n", 1.0,
            cons/surf_functor->ExactSol(pt));
    }

    // uncomment these lines for various plots
    
    // make a line plot
    {
    if(world.rank() == 0)
        printf("\n\n");
    double xmin = 0.0;
    double xmax = 2.0;
    int nx = 201;
    double dx = (xmax - xmin) / (nx - 1);
    for(int i = 0; i < nx; ++i) {
        pt[0] = xmin + i * dx;
        double uval = usol(pt);

        if(world.rank() == 0) {
            printf("%.4e %.4e\n", pt[0], uval);
        }
    }
    }

    finalize();
    
    return 0;
}
