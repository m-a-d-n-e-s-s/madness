/*
 * diamagneticpotentialfactor.cc
 *
 *  Created on: 16 Apr 2019
 *      Author: fbischoff
 */

#include <chem/diamagneticpotentialfactor.h>

namespace madness {


/// functor for the diamagnetic term in a box

/// The diamagnetic term is a 2D harmonic potential, which makes the iterations diverge.
/// Approximate it by making it constant at a specific radius: construct a smooth approximation
/// for the piecewise function f(x)={{x, x<1}, {1,x>1}}, which is squared. For this approximation
/// see https://doi.org/10.1186/s40064-016-3278-y
struct diamagnetic_boxed_functor : public FunctionFunctorInterface<double,3> {
	double radius;	//
	const double tightness;	// alpha in the article
	diamagnetic_boxed_functor(const double r, const double t=15.0) :
		radius(r), tightness(t) {}

	double operator()(const coord_3d& xyz) const {
		double r=sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
		const double a=0.5;
		const double b=0.5;
		const double c=-0.5;
		const double beta=1.0;
		const double A=a-c*beta;
		const double B=b+c;
		const double C=2.0*c/tightness;

		const double f=radius * (A +B*r/radius + C*log(1.0+exp(-tightness*(r/radius-beta))));
		return f*f-radius*radius;
	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};

void Diamagnetic_potential_factor::recompute_functions(const coord_3d& pB, const coord_3d& eB,
		const std::vector<coord_3d>& v1, const real_function_3d& sbox) {

	set_B_explicit_B(pB,eB);
	v=v1;
	v=std::vector<coord_3d> (1,{0.0,0.0,0.0});
	print("set v in diamagnetic_potential_factor to {0,0,0}");

	if (explicit_B.normf()>0.0) {

		diamagnetic_factor_=custom_factor(explicit_B,v,1.0);
		diamagnetic_factor_square=custom_factor(explicit_B,v,2.0);

		double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
		MADNESS_ASSERT(B_along_z(explicit_B));
		if (diapot_diffB!=0.0) diamagnetic_U2=real_factory_3d(world)
				.functor(diamagnetic_boxed_functor(compute_radius(diapot_diffB)));
		diamagnetic_U2.print_size("diamagnetic_U2");

		// accumulate the local terms
		print("diamagnetic(0,0,0)",diamagnetic_U2(coord_3d{0,0,0}));
		diamagnetic_U2.scale(0.125*diapot_diffB*diapot_diffB);
		print("diamagnetic(0,0,0) after scaling",diamagnetic_U2(coord_3d{0,0,0}));

		double scalar_terms=0.5*explicit_B.normf();
		print("B/2",0.5*explicit_B.normf());
		for (auto& vv : v) {
			scalar_terms-=0.5*inner(vv,vv);
			print("vv, ",0.5*inner(vv,vv));
		}
		diamagnetic_U2=diamagnetic_U2+scalar_terms;
		print("diamagnetic(0,0,0) with scalar terms",diamagnetic_U2(coord_3d{0,0,0}));


		// accumulate the vector terms
		MADNESS_ASSERT(B_along_z(explicit_B));
		for (int i=0; i<2; ++i) {
//			diamagnetic_U1[i]=real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];});
			real_function_3d ri=real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];});
			diamagnetic_U1[i]=ri*sbox;
		}
		scale(world,diamagnetic_U1,0.5*explicit_B.normf());

		coord_3d tmp;
		for (auto& vv : v) {
			print("vv, explicit_B, cross",vv,explicit_B,cross(vv,explicit_B));
			tmp+=1.0/explicit_B.normf()*cross(vv,explicit_B);
		}
		print("U1, tmp",tmp);
		for (int i=0; i<2; ++i) diamagnetic_U1[i]=diamagnetic_U1[i]-tmp[i];

	} else if (physical_B.normf()>0.0) {			// explicit_B==0
		diamagnetic_U2=real_factory_3d(world)
				.functor(diamagnetic_boxed_functor(compute_radius(physical_B.normf())));
		diamagnetic_U2.scale(0.125*physical_B.normf()*physical_B.normf());
	}


	// save the relevant terms
	save(diamagnetic_U2,"diamagnetic_U2");
	save(diamagnetic_U1[0],"diamagnetic_U1x");
	save(diamagnetic_U1[1],"diamagnetic_U1y");
	save(diamagnetic_U1[2],"diamagnetic_U1z");
	save(diamagnetic_factor_,"diamagnetic_factor");
	save(diamagnetic_factor_square,"diamagnetic_factor_square");

}

/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<complex_function_3d> Diamagnetic_potential_factor::apply_potential(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,rhs.size());
	if (physical_B.normf()==0.0) return result;

	// the diamagnetic potential
	double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
	double radius=0.0;
	if (diapot_diffB!=0.0) radius=this->compute_radius(diapot_diffB);
	result+=rhs*diamagnetic_U2 +  0.125*diapot_diffB*diapot_diffB*radius*radius*rhs;

	if (param.use_diamagnetic_factor() and explicit_B.normf()!=0.0) {
		Derivative<double_complex,3> dx(world,0);
		Derivative<double_complex,3> dy(world,1);
		std::vector<complex_function_3d> dxrhs=apply(world,dx,rhs);
		std::vector<complex_function_3d> dyrhs=apply(world,dy,rhs);

		std::vector<complex_function_3d> tmp=(diamagnetic_U1[0]*dxrhs + diamagnetic_U1[1]*dyrhs);
		save(abs(tmp[0]),"apply_potential_tmp0");
		save(abs(tmp[1]),"apply_potential_tmp1");
		result+=tmp;
//		result+=(diamagnetic_U1[0]*dxrhs + diamagnetic_U1[1]*dyrhs);
	}
	truncate(world,result);
	return result;
}


} /* namespace madness */
