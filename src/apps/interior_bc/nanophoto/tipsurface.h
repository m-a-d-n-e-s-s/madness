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
  
  $Id: test_problems.h 1856 2010-04-06 14:03:52Z mgr522 $
*/

/** \file tipsurface.h
    \brief 

*/

#ifndef MADNESS_INTERIOR_BC_TIPSURFACE_H__INCLUDED
#define MADNESS_INTERIOR_BC_TIPSURFACE_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/lbdeux.h>
#include <mra/sdf_shape_3D.h>
#include <string>

using namespace madness;

enum FunctorOutput { SURFACE, DIRICHLET_RHS, DOMAIN_MASK };

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

/** \brief Setup the tip-surface problem. */
class TipSurface : public FunctionFunctorInterface<double, 3> {
    private:
        TipSurface() : dmi(1.0) {}

    protected:
        GaussianDomainMask dmi;
        SignedDFInterface<3> *tip, *solid;
        double penalty_prefact, eps;
        double phi, d;

    public:
        /// which function to use when projecting:
        /// -# the weighted surface (SURFACE)
        /// -# the rhs of the auxiliary DE (DIRICHLET_RHS)
        /// -# the domain mask (DOMAIN_MASK)
        FunctorOutput fop;

        /// \brief Sets up the data for the problem-inspecific parts.
        TipSurface(double eps, double penalty, double phi, double d)
            : dmi(eps), tip(NULL), solid(NULL), penalty_prefact(penalty),
              eps(eps), phi(0.5*phi), d(d), fop(DIRICHLET_RHS) {

            // note on phi: distribute the potential difference across the
            // two surfaces (only half on each)

            // make the sdfs
            coord_3d normal, point;

            // solid surface is the xy-plane
            normal[0] = normal[1] = 0.0;
            normal[2] = -1.0;
            point[0] = point[1] = point[2] = 0.0;
            solid = new SDFPlane(normal, point);

            // tip apex is at (0, 0, d)
            point[2] = d;
            normal[2] = 1.0;
            tip = new SDFParaboloid(25.0 / 0.052918, point, normal);
        }

        virtual ~TipSurface() {
            if(solid != NULL)
                delete solid;
            if(tip != NULL)
                delete tip;
        }

        /// \brief The operator for projecting a MADNESS function.
        double operator() (const Vector<double, 3> &x) const {
            switch(fop) {
            case DIRICHLET_RHS:
                return dmi.mask(solid->sdf(x)) * dmi.mask(-tip->sdf(x))
                    * Inhomogeneity(x) - DirichletCond(x) *
                    (dmi.surface(solid->sdf(x)) + dmi.surface(-tip->sdf(x))) *
                    penalty_prefact;
                break;
            case SURFACE:
                return (dmi.surface(solid->sdf(x)) + dmi.surface(-tip->sdf(x)))
                    * penalty_prefact;
                break;
            case DOMAIN_MASK:
                return dmi.mask(solid->sdf(x)) * dmi.mask(-tip->sdf(x));
                break;
            default:
                error("shouldn't be here...");
                return 0.0;
                break;
            }
        }

        virtual double DirichletCond(const Vector<double, 3> &x) const {
            if(x[2] > 0.5*d)
                return phi;
            else
                return -phi;
        }

        /** \brief The PDE's inhomogeneity.

            The Inhomogeneity is 0. */
        virtual double Inhomogeneity(const Vector<double, 3> &x) const {
            return 0.0;
        }
};

/** \brief The operator needed for solving for \f$u\f$ with GMRES */
template<int NDIM>
class DirichletCondIntOp : public Operator<Function<double, NDIM> > {
    protected:
        /// \brief The Green's function
        const SeparatedConvolution<double, NDIM> &G;
        /// \brief The surface function (normalized)
        const Function<double, NDIM> &b;

        /** \brief Applies the operator to \c invec

            \note \c G is actually \f$-G\f$.

            \param[in] invec The input vector
            \param[out] outvec The action of the operator on \c invec */
        void action(const Function<double, NDIM> &invec,
                    Function<double, NDIM> &outvec) const {

                outvec = invec + G(b*invec);
                outvec.scale(-1.0);
                outvec.truncate();
        }

    public:
        DirichletCondIntOp(const SeparatedConvolution<double, NDIM> &gin,
            const Function<double, NDIM> &bin)
            : G(gin), b(bin) {}
};

#endif // MADNESS_INTERIOR_BC_TIPSURFACE_H__INCLUDED
