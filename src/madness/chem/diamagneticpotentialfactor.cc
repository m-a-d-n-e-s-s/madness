/*
 * diamagneticpotentialfactor.cc
 *
 *  Created on: 16 Apr 2019
 *      Author: fbischoff
 */

#include<madness/chem/diamagneticpotentialfactor.h>
#include <madness/world/test_utilities.h>
#include<madness/chem/masks_and_boxes.h>

namespace madness {



/// functor for the diamagnetic term in a box

/// The diamagnetic term is a 2D harmonic potential, which makes the iterations diverge.
/// Approximate it by making it constant at a specific radius: construct a smooth approximation
/// for the piecewise function f(x)={{x, x<1}, {1,x>1}}, which is squared. For this approximation
/// see https://doi.org/10.1186/s40064-016-3278-y
struct harmonic_potential_boxed : public FunctionFunctorInterface<double,3> {
	using FunctionFunctorInterface<double,3>::operator();

	double radius;		//
	double tightness;	// alpha in the article
	harmonic_potential_boxed(const double r, const double deviation=FunctionDefaults<3>::get_thresh()) :
		radius(r) {
		tightness=max_of_x_1_smooth::compute_tightness(deviation,radius);
	}

	double operator()(const coord_3d& xyz) const {

		double r=sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);	// radius in x and y directions only
		coord_3d xyz_boxed=xyz*max_of_x_1_smooth::compute_factor(r,tightness,radius);
		xyz_boxed[2]=0.0;
		double result=inner(xyz_boxed,xyz_boxed);
		MADNESS_ASSERT(!isnan(result));
		return result;
	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};

void Diamagnetic_potential_factor::recompute_factors_and_potentials() {
	diamagnetic_U2=compute_U2();
}


real_function_3d Diamagnetic_potential_factor::compute_U2() const {

	MADNESS_ASSERT(B_along_z(physical_B));
	const double pB=physical_B.normf();
	// cut vector potential: its diamagnetic potential; otherwise (B^2/8) rho^2, flattened beyond potential_radius
	real_function_3d U2_function;
	if (cut_A) {
		U2_function=cut_A_diamagnetic_potential();
	} else {
		real_function_3d harmonic_potential=real_factory_3d(world).functor(harmonic_potential_boxed(get_potential_radius()));
		U2_function=0.125*pB*pB*harmonic_potential;
	}
	return U2_function.truncate();
}

/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<complex_function_3d> Diamagnetic_potential_factor::apply_potential(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,rhs.size());
	if (physical_B.normf()==0.0) return result;

	double plateau=diamagnetic_U2(20.0, 20.0, 20.0);

	result=(diamagnetic_U2-plateau)*rhs + plateau*rhs;
	truncate(world,result);
	return result;
}


bool Diamagnetic_potential_factor::test_harmonic_potential() const {
	test_output t("entering Diamagnetic_potential_factor::test_harmonic_potential .. ");

	double radius=3;
	real_function_3d pot=real_factory_3d(world).functor(harmonic_potential_boxed(radius));
	double error1=std::abs(pot({0,0,0}) - 0.0);
	double error2=std::abs(pot({0,0,10}) - 0.0);	// orientied along z
	double error3=std::abs(pot({0,10,0}) - 9.0);	// flattened out at r=3
	double error4=std::abs(pot({10,0,0}) - 9.0);	// flattened out at r=3
	double error5=std::abs(pot({2.4,0,0}) - 5.76);	// deviation of thresh at 0.8*r

	double thresh=FunctionDefaults<3>::get_thresh();
	t.logger << "thresh " << thresh << std::endl;
	t.logger << "error1 " << error1 << std::endl;
	t.logger << "error2 " << error2 << std::endl;
	t.logger << "error3 " << error3 << std::endl;
	t.logger << "error4 " << error4 << std::endl;
	t.logger << "error5 " << error5 << std::endl;
	bool success=((error1<thresh) and (error2<thresh) and (error3<thresh) and (error4<thresh) and (error5<10*thresh));
	return t.end(success);

}


std::vector<complex_function_3d> Diamagnetic_potential_factor::make_fake_orbitals(const int n,
		const coord_3d& offset) const {

	std::vector<complex_function_3d> result;
	for (int i=0; i<n; ++i) {
		double exponent=2.0/double(i+1);
		coord_3d offset1=offset+coord_3d({0.1,0.2,0.3*i});
		double_complex phase(0.0,2.0*i);
		auto gaussian = [&exponent, &offset1, &phase](const coord_3d& r) {
			double arg=(r-offset1).normf();
			return exp(-exponent*arg*arg+phase);
		};
		complex_function_3d tmp=complex_factory_3d(world).functor(gaussian);
		result.push_back(tmp);
	}

	return result;
}


} /* namespace madness */
