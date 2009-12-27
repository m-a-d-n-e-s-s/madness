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
/// This file sets up an abstract signed distance function for use when
/// breaking up the total domain into subdomains.
///
/// SurfaceLayer uses a signed distance function to produce a mask function.
/// The surface function is 0 outside (signed distance function is positive)
/// and 1 inside (negative).

#ifndef MADNESS_MRA_SDF_SHAPE_H__INCLUDED
#define MADNESS_MRA_SDF_SHAPE_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>

namespace madness {

/// Represents a surface as a MADNESS functor.  The sdf function (signed
/// distance function) defines the surface: sdf is 0 on the surface, positive
/// outside, and negative inside.  It should be monotonic as points move
/// across the surface from one side to the other.
///
/// The sdf function should be overloaded in derived classes for the specific
/// surface.
///
/// The operator () is defined here, following the Lowengrub paper.  Other
/// variables needed by the constructor are the width of the surface layer (eps)
/// and the desired threshold of the MADNESS calculation (thresh).
///
/// NOTE: eps is the width of one side of the surface layer.
///
///   mask(x) = 0.5 * (1.0 - tanh(n * SDF(x) / eps))
///
///   such that the mask function is 0 or 1 (to the desired accuracy) when
///   |SDF(x)| > eps
///
///   n = 0.5 * ln ((2 - thresh) / thresh)
template <typename Q, int dim>
class SurfaceLayer : public FunctionFunctorInterface<Q, dim> {
	private:
		// bury the default constructor
		SurfaceLayer() {}

	protected:
		// the width
		Q eps;
		// the prefactor (includes epsilon)
		Q n;

	public:
		SurfaceLayer(const Q width, const Q thresh) : eps(width) {
			n = 0.5 * (log(2.0 - thresh) - log(thresh)) / width;
		}

		// the required overload from FunctionFunctorInterface
		virtual Q operator() (const Vector<Q, dim> &pt) const {
			double val = sdf(pt);
			if(val > eps)
				return 0.0; // we're safely outside
			else if(val < -eps)
				return 1.0; // inside
			else
				return 0.5 * (1.0 - tanh(n * val));
		}

		// the abstract signed distance function (sdf) function
		virtual Q sdf(const Vector<Q, dim> &pt) const = 0;
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
