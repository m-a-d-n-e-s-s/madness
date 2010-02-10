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

/** \brief The exact solution, for comparison */
static double exact_sol(const coord_3d &pt);

/** \brief Li, Lowengrub, Ratz, Voight domain masking with a Gaussian for the
    surface function. */
class LLRVGaussianDomainMask : public DomainMaskInterface {
private:
    LLRVGaussianDomainMask() : epsilon(0.0) {} ///< Forbidden
        
protected:
    const double epsilon; ///< The width of the transition region
    
public:
    /** \brief Constructor for the domain mask

        \param[in] epsilon The effective width of the surface */
    LLRVGaussianDomainMask(double epsilon) 
        : epsilon(epsilon)
    {}

    /** \brief Value of characteristic function at normal distance d from
               the surface

        \param[in] d The signed distance.  Negative is ``inside,''
                     positive is ``outside.''
        \return The domain mask */
    double mask(double d) const {
        if (d > 8.0*epsilon) {
            return 0.0; // we're safely outside
        }
        else if (d < -8.0*epsilon) {
            return 1.0; // inside
        }
        else {
            return 0.5 * (1.0 - tanh(3.0 * d / epsilon));
        }
    }
    
    /** \brief Derivative of characteristic function with respect to the
               normal distance

        \param[in] d The signed distance
        \return The derivative */
    double dmask(double d) const {
        if (fabs(d) > 8.0*epsilon) {
            return 0.0; // we're safely outside or inside
        }
        else {
            double tanh3d = tanh(3.0*d/epsilon);
            return 1.5*(tanh3d*tanh3d - 1.0) / epsilon;
        }
    }
    
    /** \brief Value of surface function at distance d normal to surface

        \param[in] d The signed distance
        \return The value of the surface function */
    double surface(double d) const {
        return exp(-d*d*0.5/(epsilon*epsilon)) / (sqrt(2.0*constants::pi)
            * epsilon);
    }

    /** \brief Value of d(surface)/ddistance

        \param[in] d The signed distance
        \return The derivative of the surface function */
    double dsurface(double d) const {
        double phi = mask(d);
        double dphi = dmask(d);
        return 72.0*phi*(1.0-phi)*dphi*(1.0 - 2.0*phi)/epsilon;
    }
    
    virtual ~LLRVGaussianDomainMask() {}
};

/** \brief Produces the surface Dirichlet condition */
class SurfaceCondition : public FunctionFunctorInterface<double, 3> {
    private:
        SharedPtr<DomainMaskInterface> dmi;
        SharedPtr<SignedDFInterface<3> > sdfi;
        double inveps2;

        double DirichletCond(const coord_3d &x) const {
            // Y_0^0
            return 1.0;
        }

    public:
        /// which function to read:
        /// -# the weighted surface
        /// -# the weighted Dirichlet condition
        bool useDirichlet;

        SurfaceCondition(SharedPtr<DomainMaskInterface> dmi,
            SharedPtr<SignedDFInterface<3> > sdfi, double eps)
            : dmi(dmi), sdfi(sdfi), inveps2(1.0 / (eps*eps)),
              useDirichlet(true)
        {}

        double operator() (const coord_3d &x) const {
            if(useDirichlet)
                return DirichletCond(x) * dmi->surface(sdfi->sdf(x)) * inveps2;
            else
                return dmi->surface(sdfi->sdf(x)) * inveps2;
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
    double eps;
    int k, maxiter;
    double thresh;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 5) {
            std::cerr << "Usage error: ./app_name k thresh eps maxiter"
                << std::endl;
            error("bad number of arguments");
        }
        
        eps = atof(argv[3]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");
        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");
        maxiter = atoi(argv[4]);
        if(maxiter < 1) error("maxiter >= 1");
    }
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(maxiter);
    world.gop.broadcast(k);
    
    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-2.5, 2.5);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_on_project(true);
    
    // create the domain mask, phi, and the surface function, b
    coord_3d pt(0.0); // Origin
    SharedPtr<SignedDFInterface<3> > sphere(new SDFSphere(1.0, pt));

    SharedPtr<DomainMaskInterface> llrv(new LLRVGaussianDomainMask(eps));

    // a functor for the domain mask
    SharedPtr<DomainMaskSDFFunctor<3> > phi_functor
        (new DomainMaskSDFFunctor<3>(llrv, sphere));

    // a functor for the domain surface...
    // the SURFACE option of phi_functor could be used here, however,
    // we really want eps^{-2} surface for our calculations, and this functor
    // uses that factor when MADNESS projects the function.
    SharedPtr<SurfaceCondition> surf_functor
        (new SurfaceCondition(llrv, sphere, eps));

    // phi_functor defaults to the domain mask
    real_function_3d phi = real_factory_3d(world).functor(phi_functor);
    surf_functor->useDirichlet = false;
    real_function_3d surf = real_factory_3d(world).functor(surf_functor);

    // print out the errors in volume and surface area
    // these are checks of the diffuse domain approximation
    double vol = phi.trace();
    double surfarea = surf.trace();
    if(world.rank() == 0) {
        printf("Error in Volume    = %.3e\n",
            fabs(vol-4.0*constants::pi/3.0));
        printf("Error in Surf Area = %.3e\n",
            fabs(surfarea*eps*eps-4.0*constants::pi));
    }

    // green's function
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);
    //real_convolution_3d G = BSHOperator<3>(world, 0.0, 1.0e-10, 1.0e-10);

    // make the Dirichlet condition
    // this will include the eps^{-2} factor from the penalty
    surf_functor->useDirichlet = true;
    real_function_3d usol = real_factory_3d(world).functor(surf_functor);

    surfarea = usol.trace();
    if(world.rank() == 0) {
        printf("Surface Condition Integral = %.3e\n", surfarea * eps*eps);
    }

    // transform Dirichlet condition into the right-hand side vector
    // rhs = -G*(condition)
    real_function_3d d = G(usol);
    d.scale(-1.0);
    d.truncate();

    // make an initial guess:
    // uguess = eps^2 rhs
    // the rescaling will make operator(uguess) close to d in magnitude for
    //     starting in GMRES
    usol.clear();
    usol = copy(d);
    usol.scale(-eps*eps);
    usol.compress();

    // make the operators and prepare GMRES
    DirichletCondIntOp dcio(G, surf);
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-7;//1.0e-3;
    double update_thresh = 1.0e-7;//1.0e-4;
    GMRES(space, dcio, d, usol, maxiter, resid_thresh, update_thresh, true);

    // compare to the exact solution
    real_function_3d uexact = real_factory_3d(world).f(exact_sol);
    real_function_3d uerror = (usol - uexact)*phi;
    double error = uerror.norm2();
    double cons = usol(pt);
    if(world.rank() == 0) {
        printf("u error = %.10e\n", error);
        printf("u(0) error = %.10e\n", 1.0 - cons);
    }

    finalize();
    
    return 0;
}

static double exact_sol(const coord_3d &pt) {
    const double r = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);

    if(r <= 1.0)
        return 1.0;
    else
        return 1.0 / r;
}
