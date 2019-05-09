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
struct harmonic_potential_boxed : public FunctionFunctorInterface<double,3> {
	double radius;	//
	const double tightness;	// alpha in the article
	harmonic_potential_boxed(const double r, const double t=15.0) :
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
		return f*f;
	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};

void Diamagnetic_potential_factor::print_info() const {

	print("\ninformation about the diamagnetic_potential_factor\n");
	print("coords in A space (v) ",v);
	print("coords in real space  ",coords);
	print("diamagnetic height    ",diamagnetic_height);
	print("physical magnetic field B      ",physical_B);
	print("explicit magnetic field B      ",explicit_B);
	print("");

}


void Diamagnetic_potential_factor::recompute_factors_and_potentials() {

	if (explicit_B.normf()>0.0) {

		diamagnetic_factor_=custom_factor(explicit_B,v,1.0);
		diamagnetic_factor_square=custom_factor(explicit_B,v,2.0);
		save(diamagnetic_factor_,"diamagnetic_factor");
		save(diamagnetic_factor_square,"diamagnetic_factor_square");


		double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
		MADNESS_ASSERT(B_along_z(explicit_B));
		if (diapot_diffB!=0.0) diamagnetic_U2=complex_factory_3d(world)
				.functor(harmonic_potential_boxed(compute_radius(diapot_diffB)));
		save(real(diamagnetic_U2),"harmonic_potential");
		diamagnetic_U2.scale(0.125*diapot_diffB*diapot_diffB);

		// accumulate the local terms
		double scalar_terms=0.5*explicit_B.normf();
		for (auto& vv : v) scalar_terms-=0.5*inner(vv,vv);
		diamagnetic_U2=diamagnetic_U2+scalar_terms;

		// accumulate the vector terms
		MADNESS_ASSERT(B_along_z(explicit_B));
		for (int i=0; i<2; ++i) {
			diamagnetic_U1[i]=real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];});
		}
		scale(world,diamagnetic_U1,0.5*explicit_B.normf());

		coord_3d tmp;
		for (auto& vv : v) {
			tmp+=1.0/explicit_B.normf()*cross(vv,explicit_B);
		}
		print("U1, tmp",tmp);
		for (int i=0; i<2; ++i) diamagnetic_U1[i]=diamagnetic_U1[i]-tmp[i];

	} else if (physical_B.normf()>0.0) {			// explicit_B==0
		diamagnetic_U2=complex_factory_3d(world)
				.functor(harmonic_potential_boxed(compute_radius(physical_B.normf())));
		diamagnetic_U2.scale(0.125*physical_B.normf()*physical_B.normf());
	}


//	// save the relevant terms
//	save(diamagnetic_U2,"diamagnetic_U2");
//	save(diamagnetic_U1[0],"diamagnetic_U1x");
//	save(diamagnetic_U1[1],"diamagnetic_U1y");
//	save(diamagnetic_U1[2],"diamagnetic_U1z");
//	save(diamagnetic_factor_,"diamagnetic_factor");
//	save(diamagnetic_factor_square,"diamagnetic_factor_square");

}


complex_function_3d Diamagnetic_potential_factor::compute_lz_commutator() const {
	std::vector<real_function_3d> r(3);
	for (int i=0; i<3; ++i) r[i]=real_factory_3d(world).functor([&i] (const coord_3d& r) {return r[i];});

	real_function_3d result=real_factory_3d(world);
	for (auto& vv : v) result+=vv[0]*r[0] + vv[1]*r[1] + vv[2]*r[2];
	return result*double_complex(0.0,1.0)*explicit_B.normf();
}


real_function_3d Diamagnetic_potential_factor::compute_T_commutator_scalar_term() const {

	struct U2b {
		const coord_3d B;
		const std::vector<coord_3d>& vv;
		const double epsilon;

		U2b(const coord_3d& B, const std::vector<coord_3d>& v, const double epsilon)
			: B(B), vv(v), epsilon(epsilon) {}

		double operator()(const coord_3d& r) const {
			double numerator=0.0;
			double denominator=0.0;
			for (auto& v : vv) {
				double v2=inner(v,v);
				coord_3d Bxv=cross(B,v);
				coord_3d A=0.5*cross(B,r);
				coord_3d Amv=A-v;
				double Ri=exp(-1.0/B.normf()*inner(Amv,Amv));
				double ui=(inner(r,Bxv) + v2);
				numerator+=Ri*ui;
				denominator+=Ri;
			}
			return numerator/(denominator + epsilon);
		}
	};



	MADNESS_ASSERT(B_along_z(explicit_B));
	MADNESS_ASSERT(B_along_z(physical_B));
	double eB=explicit_B.normf();
	double pB=physical_B.normf();

	double radius=compute_radius(sqrt(pB*pB - eB*eB));
	real_function_3d harmonic_potential=real_factory_3d(world).functor(harmonic_potential_boxed(radius));
	real_function_3d U2b_function=real_factory_3d(world).functor(U2b(explicit_B,get_v(),1.e-8));
	save(U2b_function,"U2b_function");

	real_function_3d U2_function=-0.5*(-eB + 0.25*(eB*eB -pB*pB)*harmonic_potential  + U2b_function);
	return U2_function;
}

real_function_3d Diamagnetic_potential_factor::compute_R_times_T_commutator_scalar_term_numerically() const {
	// the physical potential
	MADNESS_ASSERT(B_along_z(physical_B));
	double pB=physical_B.normf();
	double radius=compute_radius(pB);
	real_function_3d harmonic_potential=real_factory_3d(world).functor(harmonic_potential_boxed(radius));
	harmonic_potential.scale(0.125*pB*pB) ;

	return -0.5*div(grad_bspline_one(factor())) + harmonic_potential*factor();
}


std::vector<complex_function_3d> Diamagnetic_potential_factor::compute_T_commutator_vector_term() const {
	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,3);
	const coord_3d& B=explicit_B;
	const std::vector<coord_3d>& vv=this->v;

	for (int idim=0; idim<3; ++idim) {
		auto U1 = [&B, &vv, &idim](const coord_3d& r) {
			double B2=inner(B,B);
			double Br=inner(B,r);
			double r2=inner(r,r);
			double result=(0.25 * (B2*r2 - Br*Br) - sqrt(B2))*vv.size();
			for (auto& v : vv) {
				double v2=inner(v,v);
				coord_3d Bxv=cross(B,v);
				result+=(inner(r,Bxv) + v2);
			}
			return -0.5*result;
		};
		result[idim]=complex_factory_3d(world).functor(U1);
	}
	return result;
}


/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<std::vector<complex_function_3d> >
Diamagnetic_potential_factor::apply_potential_greensp(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> r=zero_functions_compressed<double_complex,3>(world,rhs.size());
	std::vector<std::vector<complex_function_3d> > result(3,r);
	if (physical_B.normf()==0.0) return result;

	if (explicit_B.normf()!=0.0) {

		result[0]=diamagnetic_U1[0]*rhs;
		result[1]=diamagnetic_U1[1]*rhs;
		result[2]=diamagnetic_U1[2]*rhs;
	}
	return result;
}


/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<complex_function_3d> Diamagnetic_potential_factor::apply_potential_scalar(
		const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,rhs.size());
	if (physical_B.normf()==0.0) return result;

	// the diamagnetic potential
	double diapot_diffB=sqrt(inner(physical_B,physical_B) - inner(explicit_B,explicit_B));
	double radius=0.0;
	if (diapot_diffB!=0.0) radius=this->compute_radius(diapot_diffB);
	result+=rhs*diamagnetic_U2 +  0.125*diapot_diffB*diapot_diffB*radius*radius*rhs;

	// add the integration by parts term if greensp is used
	// factor 2 for dU1x/dx + dU1y/dy = 2.0, since U1z=0
	result-=0.5*explicit_B.normf()*2.0*rhs;

	truncate(world,result);
	return result;
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
	save(abs(result[1]),"U2rhs1");

	if (explicit_B.normf()!=0.0) {
		Derivative<double_complex,3> dx(world,0);
		Derivative<double_complex,3> dy(world,1);
		dx.set_bspline1();
		dy.set_bspline1();
		std::vector<complex_function_3d> dxrhs=apply(world,dx,rhs);
		std::vector<complex_function_3d> dyrhs=apply(world,dy,rhs);

		std::vector<complex_function_3d> tmp=(diamagnetic_U1[0]*dxrhs + diamagnetic_U1[1]*dyrhs);
		result+=tmp;
	}
	truncate(world,result);
	return result;
}

void Diamagnetic_potential_factor::test_Gp_potential(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> nonGp=apply_potential(rhs);
	std::vector<complex_function_3d> Gpscalar=apply_potential_scalar(rhs);
	std::vector<std::vector<complex_function_3d> > Gpvector=apply_potential_greensp(rhs);

	Derivative<double_complex,3> dx(world,0);
	Derivative<double_complex,3> dy(world,1);
	Derivative<double_complex,3> dz(world,2);
	std::vector<complex_function_3d> divGpvector=apply(world,dx,Gpvector[0])
		+apply(world,dy,Gpvector[1])+apply(world,dz,Gpvector[2]);

	double n1=norm2(world,nonGp);
	double n2=norm2(world,Gpscalar);
	double n3=norm2(world,divGpvector);
	print(" -- test_Gp_potential, n1, n2, n3",n1,n2,n3);

	std::vector<complex_function_3d> diff=nonGp-(divGpvector+Gpscalar);
	double n4=norm2(world,diff);
	print(" -- test_Gp_potential, diffnorm ",n4);



}

/// compute a factor for comparison in coordinate space
bool Diamagnetic_potential_factor::test_factor() const {
	test_output t("entering Diamagnetic_potential_factor::test_factor .. ");

	// compute the factor in coordiate space
	const auto& B=explicit_B;
	const auto& c=coords;
	auto factor = [&B, &c](const coord_3d& r) {
		double result=0.0;
		double absB=B.normf();
		for (auto& cc : c) {
			coord_3d arg1=(r-cc);
			result+=exp(-0.25/absB*(inner(B,B)*inner(arg1,arg1) - inner(B,arg1)*inner(B,arg1)));
		}
		return result;
	};

	auto factor_square = [&B, &c](const coord_3d& r) {
		double result=0.0;
		double absB=B.normf();
		for (auto& cc : c) {
			coord_3d arg1=(r-cc);
			result+=exp(-0.25/absB*(inner(B,B)*inner(arg1,arg1) - inner(B,arg1)*inner(B,arg1)));
		}
		return result*result;
	};

	real_function_3d factor_in_coordinate_space=real_factory_3d(world).functor(factor);
	real_function_3d factor_in_A_space=custom_factor(explicit_B,v,1.0);

	real_function_3d factor_square_in_coordinate_space=real_factory_3d(world).functor(factor_square);
	real_function_3d factor_square_in_A_space=custom_factor(explicit_B,v,2.0);


	double diffnorm1=(factor_in_coordinate_space-factor_in_A_space).norm2();
	double diffnorm2=(factor_square_in_coordinate_space-factor_square_in_A_space).norm2();
	t.testos << "diffnorm1 " << diffnorm1 << std::endl;
	t.testos << "diffnorm2 " << diffnorm2 << std::endl;
	bool success=((diffnorm1<1.e-10) and (diffnorm2<1.e-10));
	return t.end(success);
}

bool Diamagnetic_potential_factor::test_scalar_potentials() const {

	test_output t("entering Diamagnetic_potential_factor::test_scalar_potentials .. ");
	const std::vector<complex_function_3d> mo=make_fake_orbitals(2);
	const std::vector<complex_function_3d> diamo=factor()*mo;
	const std::vector<complex_function_3d> dia2mo=factor_square()*mo;

	real_function_3d T_commutator_scalar_term=compute_T_commutator_scalar_term();

	double_complex t3=inner(dia2mo,T_commutator_scalar_term*mo);
	double_complex t3a=inner(diamo,compute_R_times_T_commutator_scalar_term_numerically()*mo);
//	double_complex t3=inner(dia2mo,mo);
//	double_complex t3a=inner(diamo,diamo);


	t.testos << "t3   " << t3 << std::endl;
	t.testos << "t3a  " << t3a << std::endl;
	t.testos << "error " << std::abs(t3-t3a) << std::endl;

	save(real(T_commutator_scalar_term*factor()),"t3factor");
	save(real(compute_R_times_T_commutator_scalar_term_numerically()),"t3factor_numerical");
	save(real(T_commutator_scalar_term),"t3");
	save(real(factor()),"factor");

	bool success=(std::abs(t3-t3a)<FunctionDefaults<3>::get_thresh());
	return t.end(success);
}

bool Diamagnetic_potential_factor::test_vector_potentials() const {

	test_output t("entering Diamagnetic_potential_factor::test_vector_potentials .. ");
	const std::vector<complex_function_3d> mo=make_fake_orbitals(2);
	const std::vector<complex_function_3d> diamo=factor()*mo;
	const std::vector<complex_function_3d> dia2mo=factor_square()*mo;

//	double_complex t3=inner(dia2mo[0],compute_T_commutator_vector_term()*mo[0]);
//	double_complex t3a=inner(diamo[0],grad(mo[0]));
//
//	bool success=(std::abs(t3-t3a)<1.e-10);
	bool success=false;
	return t.end(success);
}


std::vector<complex_function_3d> Diamagnetic_potential_factor::make_fake_orbitals(const int n) const {

	std::vector<complex_function_3d> result;
	for (int i=0; i<n; ++i) {
		double exponent=2.0/double(i+1);
		coord_3d offset={0.1,0.2,0.3*i};
		double_complex phase(0.0,2.0*i);
		auto gaussian = [&exponent, &offset, &phase](const coord_3d& r) {
			double arg=(r-offset).normf();
			return exp(-exponent*arg*arg+phase);
		};
		complex_function_3d tmp=complex_factory_3d(world).functor(gaussian);
		result.push_back(tmp);
	}

	return result;
}


} /* namespace madness */
