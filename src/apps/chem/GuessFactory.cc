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

/// Makes an excitation operator string based on predefined keywords
std::vector<std::string> make_predefined_exop_strings(const std::string what){
	std::vector<std::string> exop_strings;
	if (what == "dipole") {
		exop_strings.resize(3);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
	} else if (what == "x") {
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0";
	} else if (what == "y") {
		exop_strings.resize(1);
		exop_strings[0] = "y 1.0";
	} else if (what == "z") {
		exop_strings.resize(1);
		exop_strings[0] = "z 1.0";
	} else if (what == "r2") {
		exop_strings.resize(1);
		exop_strings[0] = "x 2.0 , y 2.0 , z 2.0";
	} else if (what == "quadrupole") {
		exop_strings.resize(9);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0";
		exop_strings[4] = "y 2.0";
		exop_strings[5] = "z 2.0";
		exop_strings[6] = "x 1.0 y 1.0";
		exop_strings[7] = "x 1.0 z 1.0";
		exop_strings[8] = "y 1.0 z 1.0";
	} else if (what == "dipole+") {
		exop_strings.resize(4);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
	} else if (what == "dipole+diffuse") {
		exop_strings.resize(4);
		exop_strings[0] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[1] = "x 2.0 y 1.0 , y 3.0 , y 1.0 z 2.0";
		exop_strings[2] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[3] = "x 4.0 , y 4.0 , z 4.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	} else if (what == "dipole+diffuse_big") {
		exop_strings.resize(8);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
		exop_strings[4] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[5] = "x 2.0 y 1.0 , y 3.0 ,y 1.0 z 2.0";
		exop_strings[6] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[7] = "x 4.0 , 4 2.0 , 4 2.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	} else if (what == "c2v") {
		exop_strings.resize(4);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0 , x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[1] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0 , x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[2] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0 , x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[3] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0, x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0 ";
	} else if (what == "water_first") {
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0 y 1.0, x 1.0 y 1.0 z 1.0";
	} else if (what == "c2v_big") {
		exop_strings.resize(8);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0";
		exop_strings[1] = "x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[2] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0";
		exop_strings[3] = "x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[4] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0";
		exop_strings[5] = "x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[6] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0";
		exop_strings[7] = "x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0";
	} else if (what == "big_fock") {
		exop_strings = make_auto_polynom_strings(6);
	} else if (what == "small_fock") {
		exop_strings = make_auto_polynom_strings(4);
	} else if (what == "big_fock_2") {
		exop_strings = make_auto_polynom_strings(2);
	} else if (what == "big_fock_3") {
		exop_strings = make_auto_polynom_strings(3);
	} else if (what == "big_fock_4") {
		exop_strings = make_auto_polynom_strings(4);
	} else if (what == "big_fock_5") {
		exop_strings = make_auto_polynom_strings(5);
	} else if (what == "big_fock_6") {
		exop_strings = make_auto_polynom_strings(6);
	} else if (what == "big_fock_7") {
		exop_strings = make_auto_polynom_strings(7);
	} else if (what == "big_fock_8") {
		exop_strings = make_auto_polynom_strings(8);
	} else if (what == "big_fock_9") {
		exop_strings = make_auto_polynom_strings(9);
	} else
		std::cout << "Keyword " << what << " is not known" << std::endl;

	return exop_strings;
}

/// Makes an automated excitation operator string for the excitation operators needed to create virtuals from the reference orbitals
std::vector<std::string> make_auto_polynom_strings(const size_t order){
	std::vector<std::string> exop_strings;
	for (size_t i = 0; i < order + 1; i++) {
		for (size_t j = 0; j < order + 1; j++) {
			for (size_t k = 0; k < order + 1; k++) {
				if (i + j + k > order)
					MADNESS_ASSERT(i + j + k > order); // do nothing
				else if (i == 0 && j == 0 && k == 0)
					MADNESS_ASSERT(i == 0); // do nothing
				else {
					if (i == 0 && j != 0 && k != 0)
						exop_strings.push_back(" y " + madness::stringify(j) + " z " + madness::stringify(k));
					else if (j == 0 && i != 0 && k != 0)
						exop_strings.push_back("x " + madness::stringify(i) + " z " + madness::stringify(k));
					else if (k == 0 && i != 0 && j != 0)
						exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j));
					else if (i == 0 && j == 0)
						exop_strings.push_back(" z " + madness::stringify(k));
					else if (i == 0 && k == 0)
						exop_strings.push_back(" y " + madness::stringify(j));
					else if (j == 0 && k == 0)
						exop_strings.push_back("x " + madness::stringify(i));
					else
						exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j) + " z " + madness::stringify(k));
				}
			}
		}
	}
	return exop_strings;
}

} /* namespace madness */
