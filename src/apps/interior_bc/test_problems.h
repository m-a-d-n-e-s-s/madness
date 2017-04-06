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

/** \file test_problems.h
    \brief Provides test problems for examining the convergence of
           embedded (Dirichlet) boundary conditions.

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

    This file sets up the various details of the problems... the main program
    is found in embedded_dirichlet.cc. */

#ifndef MADNESS_INTERIOR_BC_TEST_PROBLEMS_H__INCLUDED
#define MADNESS_INTERIOR_BC_TEST_PROBLEMS_H__INCLUDED

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <madness/mra/sdf_shape_3D.h>
#include <string>

using namespace madness;

enum Mask { LLRV, Gaussian };

enum FunctorOutput { SURFACE, DIRICHLET_RHS, EXACT, DOMAIN_MASK };

// load balancing structure lifted from dataloadbal.cc
template <int NDIM>
struct DirichletLBCost {
    double leaf_value;
    double parent_value;
    DirichletLBCost(double leaf_value = 1.0, double parent_value = 1.0)
        : leaf_value(leaf_value), parent_value(parent_value) {}

    double operator() (const Key<NDIM> &key,
        const FunctionNode<double, NDIM> &node) const {

        if(key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if(node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

/** \brief Abstract base class for embedded Dirichlet problems. */
class EmbeddedDirichlet : public FunctionFunctorInterface<double, 3> {
    private:
        EmbeddedDirichlet() {}
        int initial_level;
        int k;
        double thresh;
        std::string penalty_name;

    protected:
        DomainMaskInterface *dmi;
        SignedDFInterface<3> *sdfi;
        double penalty_prefact, eps;
        std::string problem_name;
        std::string problem_specific_info;
        std::string domain_mask_name;

       // is the problem homogeneous?
       virtual bool isHomogeneous() const =  0;

    public:
        /// which function to use when projecting:
        /// -# the weighted surface (SURFACE)
        /// -# the rhs of the auxiliary DE (DIRICHLET_RHS)
        /// -# the exact solution (EXACT)
        /// -# the domain mask (DOMAIN_MASK)
        FunctorOutput fop;

        /// \brief Sets up the data for the problem-inspecific parts.
        ///
        /// Also sets the FunctionDefaults for the appropriate dimension.
        EmbeddedDirichlet(double penalty_prefact, std::string penalty_name,
            double eps, int k, double thresh, Mask mask)
            : k(k), thresh(thresh), penalty_name(penalty_name), dmi(nullptr),
              sdfi(nullptr), penalty_prefact(penalty_prefact), eps(eps),
              fop(DIRICHLET_RHS) {

            // calculate some nice initial projection level
            // should be no lower than 6, but may need to be higher for small
            // eps
            initial_level = ceil(log(4.0 / eps) / log(2.0) - 4);
            if(initial_level < 6)
                initial_level = 6;
            FunctionDefaults<3>::set_initial_level(initial_level);

            switch(mask) {
            case LLRV:
                domain_mask_name = "LLRV";
                dmi = new LLRVDomainMask(eps);
                break;
            case Gaussian:
                domain_mask_name = "Gaussian";
                dmi = new GaussianDomainMask(eps);
                break;
            default:
                error("Unknown mask");
                break;
            }
        }

        virtual ~EmbeddedDirichlet() {
            if(sdfi != nullptr)
                delete sdfi;
            if(dmi != nullptr)
                delete dmi;
        }

        /// \brief Load balances using the provided Function
        void load_balance(World &world, const Function<double, 3> &f)
            const {

            LoadBalanceDeux<3> lb(world);
            lb.add_tree(f, DirichletLBCost<3>(1.0, 1.0));
            // set this map as the default
            FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,
                false));
        }

        /// \brief Do a standard Printout of the problem details
        void printout() const {
            printf("Solving problem: %s\nWavelet Order: %d\nThreshold: %.6e" \
                "\nEpsilon: %.6e\nPenalty Prefactor, %s: %.6e\nUsing %s " \
                "domain masking\n%s\n",
                problem_name.c_str(), k, thresh, eps, penalty_name.c_str(),
                penalty_prefact, domain_mask_name.c_str(),
                problem_specific_info.c_str());
            fflush(stdout);
        }

        /// \brief The operator for projecting a MADNESS function.
        double operator() (const Vector<double, 3> &x) const {
            switch(fop) {
            case EXACT:
                return ExactSol(x);
                break;
            case DIRICHLET_RHS:
                if(isHomogeneous())
                    return -DirichletCond(x) * dmi->surface(sdfi->sdf(x)) *
                        penalty_prefact;
                else
                    return dmi->mask(sdfi->sdf(x)) * Inhomogeneity(x) -
                        DirichletCond(x) * dmi->surface(sdfi->sdf(x)) *
                        penalty_prefact;
                break;
            case SURFACE:
                return dmi->surface(sdfi->sdf(x)) * penalty_prefact;
                break;
            case DOMAIN_MASK:
                return dmi->mask(sdfi->sdf(x));
                break;
            default:
                error("shouldn't be here...");
                return 0.0;
                break;
            }
        }

        virtual double DirichletCond(const Vector<double, 3> &x) const = 0;

        virtual double ExactSol(const Vector<double, 3> &x) const = 0;

        virtual double Inhomogeneity(const Vector<double, 3> &x) const = 0;

        /// \brief The surface area of the domain
        virtual double SurfaceIntegral() const = 0;

        /// \brief The volume of the domain
        virtual double VolumeIntegral() const = 0;

        /// \brief A list of points where we should compare the computed
        ///        solution to the exact solution
        virtual std::vector< Vector<double, 3> > check_pts() const {
            return std::vector< Vector<double, 3> >();
        }
};

/** \brief The constant on a sphere problem */
class ConstantSphere : public EmbeddedDirichlet {
    protected:
        double radius;

        bool isHomogeneous() const { return true; }

    public:
        ConstantSphere(int k, double thresh, double eps, std::string penalty_name,
            double penalty_prefact, double radius, Mask mask)
            : EmbeddedDirichlet(penalty_prefact, penalty_name, eps, k,
              thresh, mask), radius(radius) {

            char str[80];
            sprintf(str, "Sphere radius: %.6e\n", radius);
            problem_specific_info = str;
            problem_name = "Constant Sphere";

            // set up the domain masks, etc.
            coord_3d pt(0.0); // origin
            sdfi = new SDFSphere(radius, pt);
        }

        double DirichletCond(const Vector<double, 3> &x) const {
            return 1.0;
        }

        double ExactSol(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r <= radius)
                return 1.0;
            else
                return radius / r;
        }

        double Inhomogeneity(const Vector<double, 3> &x) const {
            return 0.0;
        }

        double SurfaceIntegral() const {
            return 4.0*constants::pi*radius*radius;
        }

        double VolumeIntegral() const {
            return 4.0*constants::pi*radius*radius*radius / 3.0;
        }

        virtual std::vector< Vector<double, 3> > check_pts() const {
            std::vector< Vector<double, 3> > vec;
            Vector<double, 3> pt;

            pt[0] = pt[1] = 0.0;
            pt[2] = 0.1 * radius;
            vec.push_back(pt);

            pt[2] = radius;
            vec.push_back(pt);

            pt[2] = 2.0;
            vec.push_back(pt);

            return vec;
        }
};

/** \brief The constant on a sphere problem, with inhomogeneity */
class InhomoConstantSphere : public EmbeddedDirichlet {
    protected:
        double radius;

        bool isHomogeneous() const { return false; }

    public:
        InhomoConstantSphere(int k, double thresh, double eps, std::string penalty_name,
            double penalty_prefact, double radius, Mask mask)
            : EmbeddedDirichlet(penalty_prefact, penalty_name, eps, k,
              thresh, mask), radius(radius) {

            char str[80];
            sprintf(str, "Sphere radius: %.6e\n", radius);
            problem_specific_info = str;
            problem_name = "Inhomogeneous Constant Sphere";

            // set up the domain masks, etc.
            coord_3d pt(0.0); // origin
            sdfi = new SDFSphere(radius, pt);
        }

        double DirichletCond(const Vector<double, 3> &x) const {
            return 1.0;
        }

        double ExactSol(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r <= radius)
                return r*r / (radius*radius);
            else
                return radius / r;
        }

        double Inhomogeneity(const Vector<double, 3> &x) const {
            return 6.0 / (radius*radius);
        }

        double SurfaceIntegral() const {
            return 4.0*constants::pi*radius*radius;
        }

        double VolumeIntegral() const {
            return 4.0*constants::pi*radius*radius*radius / 3.0;
        }

        virtual std::vector< Vector<double, 3> > check_pts() const {
            std::vector< Vector<double, 3> > vec;
            Vector<double, 3> pt;

            pt[0] = pt[1] = 0.0;
            pt[2] = 0.1 * radius;
            vec.push_back(pt);

            pt[2] = radius;
            vec.push_back(pt);

            pt[2] = 2.0;
            vec.push_back(pt);

            return vec;
        }
};

/** \brief The cos(theta) on a sphere problem */
class CosineSphere : public EmbeddedDirichlet {
    protected:
        double radius;

        bool isHomogeneous() const { return true; }

    public:
        CosineSphere(int k, double thresh, double eps, std::string penalty_name,
            double penalty_prefact, double radius, Mask mask)
            : EmbeddedDirichlet(penalty_prefact, penalty_name, eps, k,
              thresh, mask), radius(radius) {

            char str[80];
            sprintf(str, "Sphere radius: %.6e\n", radius);
            problem_specific_info = str;
            problem_name = "Cosine Sphere";

            // set up the domain masks, etc.
            coord_3d pt(0.0); // origin
            sdfi = new SDFSphere(radius, pt);
        }

        double DirichletCond(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r < 1.0e-4)
                return 0.0;
            else
                return x[2] / r;
        }

        double ExactSol(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r <= radius)
                return x[2] / radius;
            else
                return x[2] * radius * radius / (r*r*r);
        }

        double Inhomogeneity(const Vector<double, 3> &x) const {
            return 0.0;
        }

        double SurfaceIntegral() const {
            return 4.0*constants::pi*radius*radius;
        }

        double VolumeIntegral() const {
            return 4.0*constants::pi*radius*radius*radius / 3.0;
        }

        virtual std::vector< Vector<double, 3> > check_pts() const {
            std::vector< Vector<double, 3> > vec;
            Vector<double, 3> pt;

            pt[0] = pt[1] = 0.0;
            pt[2] = 0.1 * radius;
            vec.push_back(pt);

            pt[2] = radius;
            vec.push_back(pt);

            pt[2] = 2.0;
            vec.push_back(pt);

            return vec;
        }
};

/** \brief The Y_2^0 on a sphere problem */
class Y20Sphere : public EmbeddedDirichlet {
    protected:
        double radius;

        bool isHomogeneous() const { return true; }

    public:
        Y20Sphere(int k, double thresh, double eps, std::string penalty_name,
            double penalty_prefact, double radius, Mask mask)
            : EmbeddedDirichlet(penalty_prefact, penalty_name, eps, k,
              thresh, mask), radius(radius) {

            char str[80];
            sprintf(str, "Sphere radius: %.6e\n", radius);
            problem_specific_info = str;
            problem_name = "Y20 Sphere";

            // set up the domain masks, etc.
            coord_3d pt(0.0); // origin
            sdfi = new SDFSphere(radius, pt);
        }

        double DirichletCond(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r < 1.0e-4)
                return 0.0;
            else
                return 3.0* x[2] * x[2] / (r*r) - 1.0;
        }

        double ExactSol(const Vector<double, 3> &x) const {
            double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

            if(r <= radius)
                //return x[2] / radius;
                return (3.0*x[2]*x[2] - r*r) / (radius * radius);   
            else
                return radius*radius*radius / (r*r*r) * (3.0*x[2]*x[2]/(r*r) - 1.0);
        }

        double Inhomogeneity(const Vector<double, 3> &x) const {
            return 0.0;
        }

        double SurfaceIntegral() const {
            return 4.0*constants::pi*radius*radius;
        }

        double VolumeIntegral() const {
            return 4.0*constants::pi*radius*radius*radius / 3.0;
        }

        virtual std::vector< Vector<double, 3> > check_pts() const {
            std::vector< Vector<double, 3> > vec;
            Vector<double, 3> pt;

            pt[0] = pt[1] = 0.0;
            pt[2] = 0.1 * radius;
            vec.push_back(pt);

            pt[2] = radius;
            vec.push_back(pt);

            pt[2] = 2.0;
            vec.push_back(pt);

            return vec;
        }
};

/** \brief The operator needed for solving for \f$u\f$ with GMRES */
class DirichletCondIntOp : public Operator<Function<double, 3> > {
    protected:
        /// \brief The Green's function
        const SeparatedConvolution<double, 3> &G;
        /// \brief The surface function (normalized)
        const Function<double, 3> &b;

        /** \brief Applies the operator to \c invec

            \note \c G is actually \f$-G\f$.

            \param[in] invec The input vector
            \param[out] outvec The action of the operator on \c invec */
        void action(const Function<double, 3> &invec,
                    Function<double, 3> &outvec) const {

                Function<double, 3> f = b*invec;
                f.broaden();
                f.broaden();
                outvec = invec + G(f);
                f.clear();
                outvec.scale(-1.0);
                outvec.truncate();
        }

    public:
        DirichletCondIntOp(const SeparatedConvolution<double, 3> &gin,
            const Function<double, 3> &bin)
            : G(gin), b(bin) {}
};

#endif // MADNESS_INTERIOR_BC_TEST_PROBLEMS_H__INCLUDED
