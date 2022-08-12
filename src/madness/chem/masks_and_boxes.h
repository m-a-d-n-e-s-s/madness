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



/// an 2-dimensional smooth mask that is 1 inside the radius and 0 outside
template<std::size_t NDIM>
struct spherical_box {

	typedef Vector<double,NDIM> coordT;

	const double radius;
	const double tightness=11.5;	// 99.9% accurate at 0.8 radius
	const coordT offset=coordT(0.0);
	const coordT B_direction=coordT(0.0);

	spherical_box(const double r, const double deviation,
			const coordT o={0.0,0.0,0.0}, const coordT B_dir={0.0,0.0,1.0}) :
		radius(r), tightness(compute_tightness(deviation)), offset(o), B_direction(B_dir) {}

	double operator()(const coordT& xyz) const {
		// project out all contributions from xyz along the direction of the B field
		coordT tmp=(xyz-offset)*B_direction;
		const double inner=tmp[0]+tmp[1]+tmp[2];
		coordT proj=(xyz-offset)-B_direction*inner;
		double r=proj.normf();
		double v1=1.0/(1.0+exp(-tightness*(r-radius)));
		return 1.0-v1;
	}

	/// compute the tightness parameter given the deviation at 0.8 radius
	static double compute_tightness(const double deviation_at_dot_8_radius) {
		return 5.0/3.0*log((1.0-deviation_at_dot_8_radius)/deviation_at_dot_8_radius);
	}

};


}


#endif /* SRC_APPS_CHEM_MASKS_AND_BOXES_H_ */
