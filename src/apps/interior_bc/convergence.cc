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
#include <mra/lbdeux.h>
#include <mra/sdf_shape_3D.h>
#include <linalg/gmres.h>
#include <muParser/muParser.h>
#include <string>
#include "llrv_gaussian.h"
#include "loadbalcost.h"

using namespace madness;

enum Problem { CONSTANT, COSTHETA, ELLIPSE };

enum Mask { LLRV, Gaussian };

enum FunctorOutput { SURFACE, DIRICHLET_SURFACE, EXACT };

/** \brief Produces the surface Dirichlet condition */
class SurfaceProblem : public FunctionFunctorInterface<double, 3> {
    private:
        SurfaceProblem() {}

        SharedPtr<DomainMaskInterface> dmi;
        SharedPtr<SignedDFInterface<3> > sdfi;
        double penalty_prefact;
        Problem prob;

    public:
        /// which function to use when projecting:
        /// -# the weighted surface (SURFACE)
        /// -# the weighted Dirichlet condition (DIRICHLET_SURFACE)
        /// -# the exact solution (EXACT)
        FunctorOutput fop;

        SurfaceProblem(SharedPtr<DomainMaskInterface> dmi,
            SharedPtr<SignedDFInterface<3> > sdfi, double penalty_prefact,
            Problem prob)
            : dmi(dmi), sdfi(sdfi), penalty_prefact(penalty_prefact),
              prob(prob), fop(DIRICHLET_SURFACE)
        {}

        double operator() (const coord_3d &x) const {
            switch(fop) {
            case EXACT:
                return ExactSol(x);
                break;
            case DIRICHLET_SURFACE:
                return DirichletCond(x) * dmi->surface(sdfi->sdf(x)) *
                    penalty_prefact;
                break;
            case SURFACE:
                return dmi->surface(sdfi->sdf(x)) * penalty_prefact;
                break;
            default:
                error("shouldn't be here...");
                return 0.0;
                break;
            }
        }

        double DirichletCond(const coord_3d &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            switch(prob) {
            case CONSTANT:
                // Y_0^0
                return 1.0;
                break;

            case COSTHETA:
                // Y_1^0
                if(r < 1.0e-3)
                    return 0.0;
                else
                    return x[2] / r;
                break;

            default:
                error("Unknown problem in SurfaceCondition::DirichletCond");
                return 0.0;
            }
        }

        double ExactSol(const coord_3d &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            switch(prob) {
            case CONSTANT:
                if(r <= 1.0)
                    return 1.0;
                else
                    return 1.0 / r;
                break;

            case COSTHETA:
                if(r <= 1.0)
                    return x[2];
                else
                    return x[2] / (r*r*r);
                break;

            default:
                error("Unknown problem in SurfaceCondition::DirichletCond");
                return 0.0;
            }
        }
};

/** \brief The operator needed for solving for \f$u\f$ with GMRES */
class DirichletCondIntOp : public Operator<real_function_3d> {
    protected:
        /// \brief The Green's function
        const real_convolution_3d &G;
        /// \brief The surface function (normalized)
        const real_function_3d &b;

        /** \brief Applies the operator to \c invec

            \note \c G is actually \f$-G\f$.

            \param[in] invec The input vector
            \param[out] outvec The action of the operator on \c invec */
        void action(const real_function_3d &invec, real_function_3d &outvec)
            const {

                outvec = invec + G(b*invec);
                outvec.scale(-1.0);
                outvec.truncate();
        }

    public:
        DirichletCondIntOp(const real_convolution_3d &gin,
            const real_function_3d &bin)
            : G(gin), b(bin) {}
};

int main(int argc, char **argv) {
    double eps, penalty_prefact;
    coord_3d radius;
    int k;
    double thresh;
    Problem prob;
    Mask mask;
    char probname[15], maskname[20];

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 6) {
            std::cerr << "Usage error: ./app_name k thresh prob eps penalty" \
                " mask [radius, prob = CONSTANT or COSTHETA]" << std::endl;
            std::cerr << "    Where prob = 1 for CONSTANT, 2 for COSTHETA," \
                " 3 for ELLIPSE\n" << std::endl;
            std::cerr << "    Where mask = 1 for LLRV, 2 for LLRV-Gaussian\n"
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

        switch(atoi(argv[3])) {
        case 1:
            prob = CONSTANT;
            sprintf(probname, "constant");
            break;
        case 2:
            prob = COSTHETA;
            sprintf(probname, "cos(theta)");
            break;
        case 3:
            prob = ELLIPSE;
            sprintf(probname, "ellipse");
            radius[0] = 0.5;
            radius[1] = 1.0;
            radius[2] = 1.5;
            break;
        default:
            error("unknown problem type, should be 1, 2, or 3");
            break;
        }

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
        if(penalty_prefact <= 0.0) error("penalty_prefact must be positive");

        switch(atoi(argv[6])) {
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

        if(argc > 7) {
            radius[0] = atof(argv[7]);
            if(radius[0] <= 0.0) error("radius must be positive");
        }
        else
            radius[0] = 1.0;

        // print out the arguments
        printf("Solving %s problem\nWavelet order = %d\nThreshold = %.2e\n" \
            "Layer Thickness = %.6e\nPenalty Prefactor, %s = %.6e\nUsing %s " \
            "Domain Masking\n",
            probname, k, thresh, eps, argv[5], penalty_prefact, maskname);
        if(prob == CONSTANT || prob == COSTHETA)
            printf("Sphere Radius = %.6e\n", radius[0]);
        else if(prob == ELLIPSE)
            printf("Ellipse Radii = %.6e, %.6e, %.6e\n", radius[0], radius[1],
                radius[2]);
        printf("\n");
        fflush(stdout);
    }
    world.gop.broadcast(prob);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);
    world.gop.broadcast(penalty_prefact);
    world.gop.broadcast(radius);

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_on_project(true);
    //FunctionDefaults<3>::set_initial_level(8);
    
    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    SignedDFInterface<3> *sdfi = NULL;
    if(prob == CONSTANT || prob == COSTHETA)
        sdfi = new SDFSphere(radius[0], pt);
    else if(prob == ELLIPSE)
        sdfi = new SDFEllipsoid(radius, pt);
    SharedPtr<SignedDFInterface<3> > sdf(sdfi);

    DomainMaskInterface *dmip = NULL;
    if(mask == LLRV)
        dmip = new LLRVDomainMask(eps);
    else if(mask == Gaussian)
        dmip = new LLRVGaussianDomainMask(eps);
    SharedPtr<DomainMaskInterface> dmask(dmip);

    // a functor for the domain mask
    SharedPtr<DomainMaskSDFFunctor<3> > phi_functor
        (new DomainMaskSDFFunctor<3>(dmask, sdf));

    // a functor for the domain surface...
    // the SURFACE option of phi_functor could be used here, however,
    // we really want (penalty_prefact*surface) for our calculations, and this
    // functor uses that factor when MADNESS projects the function.
    SharedPtr<SurfaceProblem> surf_functor
        (new SurfaceProblem(dmask, sdf, penalty_prefact, prob));

    // project the surface function
    surf_functor->fop = SURFACE;
    real_function_3d surf = real_factory_3d(world).k(6).thresh(1.0e-4)
        .functor(surf_functor);

    if(world.rank() == 0) {
        printf("Performing load rebalancing\n");
        fflush(stdout);
    }
    // make a load balancing map off of the surface
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(surf, LBCost(1.0, 1.0));
    // set this map as the default
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));

    // reproject the surface function to the requested threshold / k
    if(k > 6 || thresh < 1.0e-4) {
        surf.clear();
        surf = real_factory_3d(world).functor(surf_functor);
    }

    // phi_functor defaults to the domain mask
    real_function_3d phi = real_factory_3d(world).functor(phi_functor);

    // print out the errors in volume and surface area
    // these are checks of the diffuse domain approximation
    double vol = phi.trace();
    if(world.rank() == 0) {
        double analvol = 0.0;
        if(prob == CONSTANT || prob == COSTHETA)
            analvol = 4.0*constants::pi/3.0*radius[0]*radius[0]*radius[0];
        else if(prob == ELLIPSE)
            analvol = 4.0*constants::pi/3.0*radius[0]*radius[1]*radius[2];
        printf("Error in Volume    = %.3e\n", fabs(vol-analvol));
    }

    // surface area is only easy for a sphere
    if(prob == CONSTANT || prob == COSTHETA) {
        double surfarea = surf.trace();
        if(world.rank() == 0) {
            printf("Error in Surf Area = %.3e\n",
                fabs(surfarea/penalty_prefact
                     -4.0*constants::pi*radius[0]*radius[0]));
        }
    }

    // green's function
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // make the Dirichlet condition
    // this will include the penalty prefactor
    surf_functor->fop = DIRICHLET_SURFACE;
    real_function_3d usol;
    real_function_3d uboundary = real_factory_3d(world).functor(surf_functor);

    // transform Dirichlet condition into the right-hand side vector
    // rhs = -G*(condition)
    real_function_3d d = G(uboundary);
    d.scale(-1.0);
    d.truncate();

    // make an initial guess:
    // uguess = eps^2 rhs
    // the rescaling will make operator(uguess) close to d in magnitude for
    //     starting in GMRES
    usol.clear();
    usol = copy(d);
    usol.scale(-1.0 / penalty_prefact);
    usol.compress();

    // make the operators and prepare GMRES
    DirichletCondIntOp dcio(G, surf);
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    GMRES(space, dcio, d, usol, maxiter, resid_thresh, update_thresh, true);

    // compare to the exact solution
    surf_functor->fop = EXACT;
    real_function_3d uexact = real_factory_3d(world).functor(surf_functor);
    real_function_3d uerror = (usol - uexact);
    double error = uerror.norm2();
    pt[0] = pt[1] = 0.0;
    pt[2] = 0.1;
    double cons = usol(pt);
    if(world.rank() == 0) {
        printf("u error = %.10e\n", error);
        printf("u(0.1)/uexact(0.1) ratio = %.10e\n",
            cons/surf_functor->ExactSol(pt));
    }
    pt[2] = 2.0;
    cons = usol(pt);
    if(world.rank() == 0) {
        printf("u(2.0)/uexact(2.0) ratio = %.10e\n",
            cons/surf_functor->ExactSol(pt));
    }
    pt[2] = 1.0;
    cons = usol(pt);
    if(world.rank() == 0) {
        printf("u(1.0)/uxact(1.0) ratio = %.10e\n",
            cons/surf_functor->ExactSol(pt));
    }

    // uncomment these lines for various plots
    /*
    // make a line plot along the positive z axis
    {
    if(world.rank() == 0)
        printf("\n\n");
    double zmin = 0.0;
    double zmax = 2.0;
    int nz = 10001;
    pt[0] = pt[1] = 0.0;
    double dz = (zmax - zmin) / (nz - 1);
    for(int i = 0; i < nz; ++i) {
        pt[2] = zmin + i * dz;
        double uval = usol(pt);
        double ueval = uexact(pt);

        if(world.rank() == 0) {
            printf("%.4e %.4e %.4e\n", pt[2], uval, ueval);
        }
    }
    }

    // make a line plot along the positive z axis
    {
    if(world.rank() == 0)
        printf("\n\n");
    double zmin = -3.0;
    double zmax = 3.0;
    int nz = 1001;
    pt[0] = pt[1] = 0.0;
    double dz = (zmax - zmin) / (nz - 1);
    for(int i = 0; i < nz; ++i) {
        pt[2] = 1.0 + (zmin + i * dz) * eps;
        double uval = (usol(pt) - 1.0) * surf(pt) / penalty_prefact;

        if(world.rank() == 0) {
            printf("%.4e %.4e\n", (zmin + i *dz), uval);
        }
    }
    }*/

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
