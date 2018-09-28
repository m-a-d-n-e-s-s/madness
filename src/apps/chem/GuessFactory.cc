/*
 * GuessFactory.cc
 *
 *  Created on: Sep 27, 2018
 *      Author: kottmanj
 */

#include "GuessFactory.h"

namespace madness {

void ExopUnaryOpStructure::operator ()(const Key<3>& key, Tensor<double>& t) const {
	Tensor<double> exop(t.ndim(), t.dims());
	const Tensor<double>& qp = cdata.quad_x;
	fcube(key, (*exfunc), qp, exop);
	t.emul(exop);
}

/// compute the centroid of a function i.e. c[xi]=<f|xi|f>/<f|f> i.e. position expectation value
coord_3d compute_centroid(const real_function_3d& f) {
	coord_3d result(0.0);
	real_function_3d density = f * f;
	const double integral = density.trace();
	for (size_t x = 0; x < 3; ++x) {
		const auto mf = xyz(x);
		real_function_3d m = real_factory_3d(f.world()).functor(mf);
		result[x] = (m * density).trace() / integral;
	}
	return result;
}
std::vector<coord_3d> compute_centroids(const vecfuncT & vf){
	std::vector<coord_3d> results(vf.size());
	if(vf.empty()) return results;

	World& world=vf.front().world();
	const Tensor<double> denom=inner(world,vf,vf);
	std::vector<Tensor<double> >nums(3);
	for (size_t x = 0; x < 3; ++x) {
		const auto mf = xyz(x);
		real_function_3d m = real_factory_3d(world).functor(xyz(x));
		nums[x]=inner(world,vf,m*vf);
	}

	// create correct data structure
	for(size_t i=0;i<vf.size();++i){
		for(size_t x=0;x<3;++x){
			results[i][x]=(nums[x][i]/denom[i]);
		}
	}

	return results;
}

} /* namespace madness */
