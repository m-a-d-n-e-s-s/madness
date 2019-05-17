/*
 * masks_and_boxes.h
 *
 *  Created on: 16 May 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_MASKS_AND_BOXES_H_
#define SRC_APPS_CHEM_MASKS_AND_BOXES_H_


#include <madness/world/vector.h>
#include <math.h>

namespace madness {

/// construct a smooth approximation for the piecewise function f(x)={{x, x<1}, {1,x>1}}.
/// see https://doi.org/10.1186/s40064-016-3278-y
struct max_of_x_1_smooth {

	/// return a smooth transition of linear r to a constant at radius r

	/// @param[in]	r			norm of the n-dimensional vector
	/// @param[in]	tightness 	see the compute_tightness() function
	/// @param[in]	max_radius	where the curve flattens
	/// @return		factor		multiply your vector with the factor xyz->xyz*factor
	static double compute_factor(const double& r, const double& tightness, const double& rmax) {

		const double a=0.5;
		const double b=0.5;
		const double c=-0.5;
		const double beta=1.0;
		const double A=a-c*beta;
		const double B=b+c;
		const double C=2.0*c/tightness;

		const double f=rmax * (A +B*r/rmax + C*log(1.0+exp(-tightness*(r/rmax-beta))));
		const double factor=(r>1.e-10) ? f/r : 1.0;	// prevent division by zero
		return factor;
	}

	/// compute the tightness parameter given the deviation at 0.8 radius

	/// simple linear fit
	static double compute_tightness(const double deviation_at_dot_8_radius, const double rmax) {
		const double a=-1.0577335859862533;
		const double b=-1.0/0.09618399253086536;
		return (log10(deviation_at_dot_8_radius/rmax)-a)*b;
	}

};




}


#endif /* SRC_APPS_CHEM_MASKS_AND_BOXES_H_ */
