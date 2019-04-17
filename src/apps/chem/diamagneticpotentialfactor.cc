/*
 * diamagneticpotentialfactor.cc
 *
 *  Created on: 16 Apr 2019
 *      Author: fbischoff
 */

#include <chem/diamagneticpotentialfactor.h>

namespace madness {

void Diamagnetic_potential_factor::recompute_functions(const double& pB, const double& eB,
		const std::vector<coord_3d>& v) {

	set_B_explicit_B({0,0,pB},{0,0,eB});

	diamagnetic_factor_=real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});
	diamagnetic_factor_square=real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});
	diamagnetic_boxed=real_factory_3d(world).functor([] (const coord_3d& r) {return 0.0;});

	if (physical_B.normf()!=0.0) {

		if (param.use_diamagnetic_factor()) {
			diamagnetic_factor_=real_factory_3d(world).functor([] (const coord_3d& r) {return 0.0;});
			// loop over all atoms
			for (auto& vv: v) {
				real_function_3d tmp1=custom_factor(explicit_B,vv,1.0);
				diamagnetic_factor_+=tmp1;
			}
			diamagnetic_factor_square=diamagnetic_factor_.square();
		}

		double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
		if (diapot_diffB!=0.0) diamagnetic_boxed=real_factory_3d(world)
				.functor(diamagnetic_boxed_functor(compute_radius(diapot_diffB)));

		save(0.125*diapot_diffB*diapot_diffB*diamagnetic_boxed,"diabox");
		save(diamagnetic_factor_,"diamagnetic_factor");
	}
}

/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<complex_function_3d> Diamagnetic_potential_factor::apply_potential(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,rhs.size());
	if (physical_B.normf()==0.0) return result;

	// the diamagnetic potential
	double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
	if (diapot_diffB!=0.0) {
		const double radius=this->compute_radius(diapot_diffB);
		result+= 0.125*diapot_diffB*diapot_diffB*(rhs*diamagnetic_boxed + radius*radius*rhs);
	}

	if (param.use_diamagnetic_factor() and explicit_B.normf()!=0.0) {
		Derivative<double_complex,3> dx(world,0);
		Derivative<double_complex,3> dy(world,1);
		real_function_3d x=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
		real_function_3d y=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});
		std::vector<complex_function_3d> dxrhs=apply(world,dx,rhs);
		std::vector<complex_function_3d> dyrhs=apply(world,dy,rhs);

		result+=0.5*fabs(explicit_B.normf())*(rhs + (x*dxrhs + y*dyrhs));
	}
	truncate(world,result);
	return result;
}


} /* namespace madness */
