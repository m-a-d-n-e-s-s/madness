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

/*!
  \file mra/sdf_shape.h
  \brief Defines abstract interfaces and concrete classes for use with diffuse boundary calculations
  \ingroup mrabcint
*/

#ifndef MADNESS_MRA_SDF_SHAPE_H__INCLUDED
#define MADNESS_MRA_SDF_SHAPE_H__INCLUDED

#include <mra/mra.h>

namespace madness {

    /// The interface to be provided by shapes defined by signed distance from surface

    ///\ingroup mrabcint
    template <int NDIM>
    class ShapeDFInterface {
    public:
        /// Should return the signed normal distance from the surface
        virtual double sdf(const Vector<double,NDIM>& x) const = 0;

        /// Should return the gradient of the signed normal distance from the surface (i.e., \c dsdf(x)/dx[i] )
        virtual Vector<double,NDIM> grad_sdf(const Vector<double,NDIM>& x) const = 0;

        virtual ~ShapeDFInterface() {}
    };

    /// The interface to be provided by masking functions for shapes defined by normal distance from surface

    ///\ingroup mrabcint
    template <int NDIM>
    class ShapeMaskInterface {
    public:
        /// Should return the characteristic function (1 in interior, 0 in exterior, 1/2 on boundary)

        /// @param d The signed normal distance from the surface
        /// @param d The mask or characteristic function
        virtual double mask(double d) const = 0;

        /// Should return the derivative of the characteristic function w.r.t. the distance

        /// @param d The signed normal distance from the surface
        /// @return Derivative of the mask or characteristic function
        virtual double dmask(double d) const = 0;

        /// Should return the value of the normalized surface layer function

        /// Normalized means the volume integral of this function should converge to the
        /// surface area in the limit of a either an infinitely thin surface or
        /// zero curvature.  The function thus acts as a "delta function" located at the boundary.
        /// @param d  The signed normal distance from the surface
        /// @return Normalized surface layer function
        virtual double surface(double d) const = 0;

        /// Should return the derivative of the normalized surface layer function

        /// @param d The signed normal distance from the surface
        /// @return Derivative of the normalized surface layer function
        virtual double dsurface(double d) const = 0;

        virtual ~ShapeMaskInterface() {}
    };

    /// Provides MADNESS interface to compute mask or characteristic function

    /// \ingroup mrabcint
    template <int NDIM>
    class ShapeMaskFunctor : public FunctionFunctorInterface<double,NDIM> {
        ShapeMaskInterface<NDIM>* mask;
        ShapeDFInterface<NDIM>* shape;
    public:
        /// Constructor for mask functor

        /// The constructor takes ownership of the provided mask and shape and deletes them when finished.
        /// @param mask Pointer to object providing maksing function (deleted when finished)
        /// @param shape Pointer to object providing shape distance function (deleted when finished)
        ShapeMaskFunctor(ShapeMaskInterface<NDIM>* mask, ShapeDFInterface<NDIM>* shape) 
            : mask(mask), shape(shape)
        {}

        /// For MADNESS FunctionFunctorInterface

        /// @param x Point to compute value
        /// @return Value of characteristic function
        double operator()(const Vector<double,NDIM>& x) const {
            return mask->mask(shape->sdf(x));
        }

        virtual ~ShapeMaskFunctor() {
            delete mask;
            delete shape;
        };
    };

    /// Provides MADNESS interface to compute surface

    /// \ingroup mrabcint
    template <int NDIM>
    class ShapeSurfaceFunctor : public FunctionFunctorInterface<double,NDIM> {
        ShapeMaskInterface<NDIM>* mask;
        ShapeDFInterface<NDIM>* shape;
    public:
        /// Constructor for surface functor

        /// The constructor takes ownership of the provided mask and shape and deletes them when finished.
        /// @param mask Pointer to object providing maksing function (deleted when finished)
        /// @param shape Pointer to object providing shape distance function (deleted when finished)
        ShapeSurfaceFunctor(ShapeMaskInterface<NDIM>* mask, ShapeDFInterface<NDIM>* shape) 
            : mask(mask), shape(shape)
        {}

        /// For MADNESS FunctionFunctorInterface

        /// @param x Point to compute value
        /// @return Value of surface 
        double operator()(const Vector<double,NDIM>& x) const {
            return mask->surface(shape->sdf(x));
        }

        virtual ~ShapeSurfaceFunctor() {
            delete mask;
            delete shape;
        };
    };

    /// Provides MADNESS interface to compute d(surface)/dx[i]

    /// \ingroup mrabcint
    /// Takes ownership of provided mask and shape and deletes them
    /// when finished.
    template <int NDIM>
    class ShapeDSurfaceDxFunctor : public FunctionFunctorInterface<double,NDIM> {
        ShapeMaskInterface<NDIM>* mask;
        ShapeDFInterface<NDIM>* shape;
        const int axis;
    public:
        /// Constructor for derivative surface functor

        /// The constructor takes ownership of the provided mask and shape and deletes them when finished.
        /// @param mask Pointer to object providing maksing function (deleted when finished)
        /// @param shape Pointer to object providing shape distance function (deleted when finished)
        /// @param axis Axis (0,...,NDIM-1) to differentiate w.r.t.
        ShapeDSurfaceDxFunctor(ShapeMaskInterface<NDIM>* mask, ShapeDFInterface<NDIM>* shape, int axis) 
            : mask(mask), shape(shape), axis(axis)
        {}

        /// For MADNESS FunctionFunctorInterface

        /// @param x Point to compute value
        /// @return Value of d(surface)/dx[axis]
        double operator()(const Vector<double,NDIM>& x) const {
            Vector<double,NDIM> g = shape->grad_sdf(x);
            return g[axis]*mask->dsurface(shape->sdf(x));
        }

        virtual ~ShapeDSurfaceDxFunctor() {
            delete mask;
            delete shape;
        };
    };

    /*!
      \brief Provides the Lowengrub characteristic and surface functions

      \ingroup mrabcint

      X. Li, J. Lowengrub, A. R&auml;tz, and A. Voight, "Solving PDEs in Complex
      Geometries: A Diffuse Domain Approach," Commun. Math. Sci., 7, p81-107, 2009.

      The derived class should implement the \c sdf() (signed distance
      function) interface to define the surface: \c sdf is 0 on the
      surface, positive outside, and negative inside.  It should be
      monotonic as points move across the surface from one side to the
      other, and within distance \f$ 8 \epsilon \f$ of the
      surface should be proportional to the normal distance from the
      surface (beyond this distance the switching function is one/zero to
      machine precision).

      The user-specified width (\f$ \epsilon \f$) of the surface is
      employed in the masking or characteristic or phase-field 
      function as follows
      \f[ 
      \phi(x) = \frac{1}{2} \left( 1 - \tanh \frac{ 3 \mbox{sdf}(x) }{\epsilon} \right) 
      \f] 
      where \f$ x \f$ is the point. 

      The normalized surface layer used by Lowengrub in the Dirichlet
      algorithm is 
      \f[ 
      B(s) = 36 \epsilon^{-1} \phi(s)^2 (1 - \phi(s))^2 
      \f]
      where the constant 36 is chosen so that 
      \f[ 
      \int_{-\infty}^{\infty} B(s) \, ds = 1 
      \f]
      For this function the parameter \f$ \epsilon \f$ is an effective measure of the
      full width of the surface layer since 
      \f[
      \int_{-\epsilon/2}^{\epsilon/2} B(s) \, ds = 0.987 
      \f]
      and 
      \f[
      \int_{-\epsilon}^{\epsilon} B(s) \, ds = 0.999963 
      \f]
      The \f$ \epsilon^{-1} \f$ is included here to normalize the volume integral of
      the resulting function to be the surface area.
    */
    template <int NDIM>
    class ShapeMask : public ShapeMaskInterface<NDIM> {
    private:
        ShapeMask() {} ///< Forbidden
        
    protected:
        const double epsilon; ///< The width
        const double sign;    ///< 1 for standard surface, -1 to complement
        
    public:
        /// Constructor for mask

        /// @param epsilon The effective width of the surface
        /// @param complement If true, compute complement the volume (i.e., invert the sense of the surface)
        ShapeMask(double epsilon, bool complement=false) 
            : epsilon(epsilon)
            , sign(complement?-1.0:1.0) 
        {}

        /// Value of characteristic function at normal distance d from surface
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
        
        /// Derivative of characteristic function w.r.t. normal distance
        double dmask(double d) const {
            if (d > 8.0*epsilon) {
                return 0.0; // we're safely outside
            }
            else if (d < -8.0*epsilon) {
                return 0.0; // inside
            }
            else {
                double tanh3d = tanh(3.0*d/epsilon);
                return -1.5*(1.0 - tanh3d*tanh3d) / epsilon;
            }
        }
        
        /// Value of surface function at distance d normal to surface
        double surface(double d) const {
            double phi = mask(d);
            double phic = 1.0 - phi;
            return 36.0*phi*phi*phic*phic/epsilon;
        }

        /// Value of d(surface)/ddistance
        double dsurface(double d) const {
            double phi = mask(d);
            double dphi = dmask(d);
            return 72.0*phi*(1.0-phi)*dphi*(1.0 - 2.0*phi)/epsilon;
        }
        
        virtual ~ShapeMask() {}
    };

    /// Convenience wrapper in 3D combining Lowengrub mask with arbitrary shape to compute characteristic function (mask)

    /// \ingroup mrabcint
    SharedPtr< FunctionFunctorInterface<double,3> > shape_mask(double epsilon, ShapeDFInterface<3>* sdf) {
        return SharedPtr< FunctionFunctorInterface<double,3> >(new ShapeMaskFunctor<3>(new ShapeMask<3>(epsilon), sdf));
    }

    /// Convenience wrapper in 3D combining Lowengrub mask with arbitrary shape to compute surface function

    /// \ingroup mrabcint
    SharedPtr< FunctionFunctorInterface<double,3> > shape_surface(double epsilon, ShapeDFInterface<3>* sdf) {
        return SharedPtr< FunctionFunctorInterface<double,3> >(new ShapeSurfaceFunctor<3>(new ShapeMask<3>(epsilon), sdf));
    }

    /// Convenience wrapper in 3D combining Lowengrub mask with arbitrary shape to compute derivative of surface function

    /// \ingroup mrabcint
    SharedPtr< FunctionFunctorInterface<double,3> > shape_surface_derivative(double epsilon, ShapeDFInterface<3>* sdf, int axis) {
        return SharedPtr< FunctionFunctorInterface<double,3> >(new ShapeDSurfaceDxFunctor<3>(new ShapeMask<3>(epsilon), sdf, axis));
    }

} // end of madness namespace

#endif // MADNESS_MRA_SDF_SHAPE_H__INCLUDED
