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
#include "llrv_gaussian.h"
#include "loadbalcost.h"

using namespace madness;

enum Problem { CONSTANT, COSTHETA };

enum Mask { LLRV, LLRVGaussian };

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
        /// -# the weighted surface (false)
        /// -# the weighted Dirichlet condition (true)
        bool useDirichlet;

        /// use the exact solution when projecting?
        bool useExact;

        bool useDSurface;

        SurfaceProblem(SharedPtr<DomainMaskInterface> dmi,
            SharedPtr<SignedDFInterface<3> > sdfi, double penalty_prefact,
            Problem prob)
            : dmi(dmi), sdfi(sdfi), penalty_prefact(penalty_prefact),
              prob(prob), useDirichlet(true), useExact(false), useDSurface(false)
        {}

        double operator() (const coord_3d &x) const {
            if(useExact)
                return ExactSol(x);
            if(useDirichlet)
                return DirichletCond(x) * dmi->surface(sdfi->sdf(x)) *
                    penalty_prefact;
            else if(useDSurface)
                return dmi->dsurface(sdfi->sdf(x)) * penalty_prefact;
            else
                return dmi->surface(sdfi->sdf(x)) * penalty_prefact;
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
    int k, maxiter;
    double thresh;
    Problem prob;
    Mask mask;
    char probname[15], maskname[20];

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 7) {
            std::cerr << "Usage error: ./app_name prob k thresh eps mask "\
                "maxiter" << std::endl;
            std::cerr << "    Where prob = 1 for CONSTANT, 2 for COSTHETA"
                << std::endl;
            std::cerr << "    Where mask = 1 for LLRV, 2 for LLRV-Gaussian"
                << std::endl;
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        switch(atoi(argv[1])) {
        case 1:
            prob = CONSTANT;
            sprintf(probname, "constant");
            break;
        case 2:
            prob = COSTHETA;
            sprintf(probname, "cos(theta)");
            break;
        default:
            error("unknown problem type, should be 1 or 2");
            break;
        }

        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        thresh = atof(argv[3]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        k = atoi(argv[2]);
        if(k < 4) error("cheapskate");

        maxiter = atoi(argv[6]);
        if(maxiter < 1) error("maxiter >= 1");

        switch(atoi(argv[5])) {
        case 1:
            mask = LLRV;
            sprintf(maskname, "LLRV");
            break;
        case 2:
            mask = LLRVGaussian;
            sprintf(maskname, "LLRV-Gaussian");
            break;
        default:
            error("unknown domain mask type, should be 1 or 2");
            break;
        }

        // print out the arguments
        printf("Solving %s problem\nWavelet order = %d\nThreshold = %.2e\n" \
            "Layer Thickness = %.3e\nUsing %s Domain Masking\n\n",
            probname, k, thresh, eps, maskname);
    }
    world.gop.broadcast(prob);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(maxiter);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);

    penalty_prefact = 2.0/eps;
 
    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_on_project(true);
    FunctionDefaults<3>::set_initial_level(8);
    
    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    SharedPtr<SignedDFInterface<3> > sphere(new SDFSphere(1.0, pt));

    DomainMaskInterface *dmip = NULL;
    if(mask == LLRV)
        dmip = new LLRVDomainMask(eps);
    else if(mask == LLRVGaussian)
        dmip = new LLRVGaussianDomainMask(eps);
    SharedPtr<DomainMaskInterface> llrv(dmip);

    // a functor for the domain mask
    SharedPtr<DomainMaskSDFFunctor<3> > phi_functor
        (new DomainMaskSDFFunctor<3>(llrv, sphere));

    // a functor for the domain surface...
    // the SURFACE option of phi_functor could be used here, however,
    // we really want eps^{-2} surface for our calculations, and this functor
    // uses that factor when MADNESS projects the function.
    SharedPtr<SurfaceProblem> surf_functor
        (new SurfaceProblem(llrv, sphere, penalty_prefact, prob));

    // project the surface function
    surf_functor->useDirichlet = false;
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
    surf.clear();
    surf = real_factory_3d(world).functor(surf_functor);

    // phi_functor defaults to the domain mask
    real_function_3d phi = real_factory_3d(world).functor(phi_functor);

    // print out the errors in volume and surface area
    // these are checks of the diffuse domain approximation
    double vol = phi.trace();
    double surfarea = surf.trace();
    if(world.rank() == 0) {
        printf("Error in Volume    = %.3e\n",
            fabs(vol-4.0*constants::pi/3.0));
        printf("Error in Surf Area = %.3e\n",
            fabs(surfarea/penalty_prefact-4.0*constants::pi));
    }

    // green's function
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // make the Dirichlet condition
    // this will include the penalty prefactor
    surf_functor->useDirichlet = true;
    real_function_3d usol;
    real_function_3d uboundary = real_factory_3d(world).functor(surf_functor);

    // make a line plot along the positive z axis of eps^{-2} g S
    /*{
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

        if(world.rank() == 0) {
            printf("%.4e %.4e\n", pt[2], uval);
        }
    }
    }*/

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
    GMRES(space, dcio, d, usol, maxiter, resid_thresh, update_thresh, true);

    // compare to the exact solution
    surf_functor->useDirichlet = false;
    surf_functor->useDSurface = false;
    surf_functor->useExact = true;
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

    // uncomment these lines for line plots
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
