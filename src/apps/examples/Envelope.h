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


/*!
  \file examples/Envelope.h
  \brief provides a few common envelope functions

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/Envelope.h>here</a>.

*/


#ifndef ENVELOPE_H_
#define ENVELOPE_H_


using namespace madness;

namespace madness {


/// a Gaussian envelope with width sigma
template<typename T, std::size_t NDIM>
class GaussianEnvelope {
	const double sigma;

public:
	/// constructor takes the width, use the box size as default
	GaussianEnvelope(const double sigma=FunctionDefaults<3>::get_cell_min_width())
		: sigma(sigma) {}

	double operator()(const coordT& r) const {
		double r2=0.0;
		for (std::size_t i=0; i<NDIM; ++i) r2+=r[i]*r[i];
		return exp(-2.0/(sigma*sigma)*r2);
	}
};

/// the envelope function of MRA, basic theory, 2004, Eq. (22)

/// designed for forcing a function to vanish at the box boundaries
template<typename T, std::size_t NDIM>
class CubicSoftBoxEnvelope {

	/// rim size as fraction of the total box length
	const double tau;

public:
	/// constructor takes the rim size as input

	/// @param[in]	tau rim size as a fraction of the total box length
	CubicSoftBoxEnvelope(const double tau) : tau(tau) {}

	double operator()(const coordT& r) const {
		coordT rsim;
		user_to_sim(r,rsim);
		double result=1.0;
		for (std::size_t i=0; i<NDIM; ++i) result*=m(rsim[i],tau);
		return result;
	}

private:
	/// @param[in]	the coordinate in simulation dimensions
	double m(const double x, const double tau) const {
		if (x<=tau) return s(x/tau);
		else if (x>(1.0-tau)) return s((1.0-x)/tau);
		else return 1.0;
	}

	double s(const double& x) const {return (x*x*(3.0-2.0*x));}

};

/// muffin tin function

/// designed for restricting a function into a given part of the box
template<typename T, std::size_t NDIM>
class SphericalSoftBoxEnvelope {

	/// cutoff radius in user coordinates
	double cutoff_;

public:
	/// constructor

	/// @param[in]	c 	cutoff radius in user coordinates
	SphericalSoftBoxEnvelope(const double c) : cutoff_(c) {}

	/// return the function value

	/// @param[in]	the distance to the center of the box
	double operator()(double d) const {
		d/=cutoff_;
		if (d>1.0) return 1.0;
		return d*d*(3.0-2.0*d);
	}

	/// return the function value

	/// @param[in]	relative coordinates r-location in user coords
	double operator()(const coordT& r) const {
		const double d=r.normf()/cutoff_;
		if (d>1.0) return 1.0;
		return d*d*(3.0-2.0*d);
	}

	/// return the gradient of the function \nabla m

	/// @param[in]	relative coordinates in user coord
	coordT grad(const coordT& r) const {
		const double d=r.normf()/cutoff_;
		if (d>1.0) return coordT(0.0);
		return n12(r)*6.0*(d-d*d);
	}

	/// return the scalar part of the gradient of the function \nabla m

	/// @param[in]	relative coordinates r-location in user coords
	double scalar_grad(double d) const {
		d/=cutoff_;
		if (d>1.0) return 0.0;
		return 6.0*(d-d*d);
	}

//	/// return the laplacian of the function
//
//	/// @param[in]	relative coordinates in user coords
//	double laplacian(const coordT& r) const {
//		const double d=r.normf()/cutoff_;
//		if (d>1.0) return 0.0;
//		return 18.0-24.0*d;
//	}

	/// return the laplacian of the function

	/// @param[in]	the distance to the center of the box
	double laplacian(double d) const {
		d/=cutoff_;
		if (d>1.0) return 0.0;
		return 18.0-24.0*d;
	}

	/// return the cutoff
	const double& cutoff() const {return cutoff_;}
private:

	/// helper function unit vector in direction r
	static coord_3d n12(const coord_3d& r) {
		const double norm=r.normf();
		if (norm<1.e-6) return coord_3d(0.0);
		return r*(1.0/norm);
	}

};


}

#endif /* ENVELOPE_H_ */


