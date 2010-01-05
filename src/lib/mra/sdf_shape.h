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
  \brief An abstract signed distance function for use when forming subdomains.
  \ingroup mrabcint
*/

#ifndef MADNESS_MRA_SDF_SHAPE_H__INCLUDED
#define MADNESS_MRA_SDF_SHAPE_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>

namespace madness {

    /*!
      \brief Computes the characteristic function of an interior surface as a MADNESS functor.  

      \ingroup mrabcint

      The characteristic function is one in the interior, zero
      exterior, and one half on the surface.  It is computed from 
      the distance function implemented by the derived class.
      
      The derived class should implement the \c sdf() (signed distance
      function) interface to define the surface: \c sdf is 0 on the
      surface, positive outside, and negative inside.  It should be
      monotonic as points move across the surface from one side to the
      other, and within distance \f$ 8 \epsilon \f$ of the
      surface should be proportional to the normal distance from the
      surface (beyond this distance the switching function is one/zero to
      machine precision).

      The \c operator() required by the MADNESS FunctorInterface is
      implemented here to compute the characteristic function following
      the switching function of the Lowengrub paper.  To be concrete,
      the user-specified width (\f$ \epsilon \f$) of the surface is
      employed in the masking or characteristic or phase-field 
      function as follows
      \f[ 
      \phi(x) = \frac{1}{2} \left( 1 - \tanh \frac{ 3 \mbox{sdf}(x) }{\epsilon} \right) 
      \f] 
      where \f$ x \f$ is the point.  Roughly, at distance \f$ n \epsilon \f$ from
      the boundary, the value of the switching function is \f$ 10^{-2n} \f$.

      X. Li, J. Lowengrub, A. R&auml;tz, and A. Voight, "Solving PDEs in Complex
      Geometries: A Diffuse Domain Approach," Commun. Math. Sci., 7, p81-107, 2009.

    */
    template <typename Q, int dim>
    class SurfaceLayerInterface : public FunctionFunctorInterface<Q, dim> {
    private:
        SurfaceLayerInterface() {} ///< Forbidden
        
    protected:
        const double eps; ///< The width
        
    public:
        /// Constructor for surface layer interface

        /// @param width The effective width of the surface
        SurfaceLayerInterface(double width) : eps(width) {}
        
        /// Required overload for FunctionFunctorInterface
        virtual Q operator() (const Vector<double, dim> &pt) const {
            double val = sdf(pt);
            if (val > 8.0*eps) {
                return 0.0; // we're safely outside
            }
            else if (val < -8.0*eps) {
                return 1.0; // inside
            }
            else {
                return 0.5 * (1.0 - tanh(3.0 * val / eps));
            }
        }
        
        /// Derived class provides normal distance from surface for any point

        /// @param pt The point in user coordinates
        /// @return Signed distance from surface
        virtual double sdf(const Vector<double, dim> &pt) const = 0;
        
        virtual ~SurfaceLayerInterface() {}
    };
    
    /// Complement unary operation to reverse inside and outside of a mask.
    ///
    /// This is to be used on a madness function generated from the SDF object.
    template<typename Q, int dim>
    inline static void mask_complement(const Key<dim> &key, Tensor<Q> &t) {
	UNARY_OPTIMIZED_ITERATOR(Q, t,
                                 *_p0 = 1.0 - *_p0);
    }
    
} // end of madness namespace

#endif // MADNESS_MRA_SDF_SHAPE_H__INCLUDED
