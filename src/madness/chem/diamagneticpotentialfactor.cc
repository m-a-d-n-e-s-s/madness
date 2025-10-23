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



/// compute the expression (\sum_i R_i)^(-1) \sum_i R_i arg(r,B,v_i)

/// see compute_T_commutator_vector_term() for an example
/// @param[in]	arg	the function argument arg(r,B,v)
/// @param[in]	B	magnetic field
/// @param[in]	v	the position of the nuclei (centers of the Gaussian) in A space
struct R_times_arg_div_R {

	const coord_3d& B;
	const std::vector<coord_3d>& vv;
	std::function<double(const coord_3d& r, const coord_3d& B, const coord_3d& vv)>& arg;
	const double radius;
	double tightness;

	R_times_arg_div_R(std::function<double(const coord_3d& r, const coord_3d& B, const coord_3d& v)>& arg,
			const coord_3d& B, const std::vector<coord_3d>& vv, const double radius) : B(B), vv(vv), arg(arg),
					radius(radius) {
		tightness=max_of_x_1_smooth::compute_tightness(FunctionDefaults<3>::get_thresh(),radius);
	}

	double operator()(const coord_3d& r) const {
		coord_3d r_boxed=r;
		r_boxed*=max_of_x_1_smooth::compute_factor(r.normf(),tightness,radius);
		int imax=0;
		double Riuimax=-1.0;
		//double Rmax=0.0;
		double umax=0.0;
		double Amv_square_max=0.0;
		for (size_t i=0; i<vv.size(); ++i) {
			const coord_3d& v = vv[i];
			coord_3d Amv=0.5*cross(B,r_boxed)-v;
			double Amv_square=inner(Amv,Amv);
			double Ri=exp(-1.0/B.normf()*Amv_square);
			double ui=arg(r_boxed,B,v);

			if (Ri*ui>Riuimax) {
			        //Rmax=Ri;
				umax=ui;
				Amv_square_max=Amv_square;
				Riuimax=Ri*ui;
				imax=i;
			}
		}
		double numerator=umax;
		double denominator=1.0;
		for (size_t i=0; i<vv.size(); ++i) {
		        if (i==size_t(imax)) continue;
			const coord_3d& v = vv[i];
			coord_3d Amv=0.5*cross(B,r_boxed)-v;
			double Amv_square=inner(Amv,Amv);

			double Ri_div_Rmax=exp(-1.0/B.normf()*(Amv_square-Amv_square_max));
			double ui=arg(r_boxed,B,v);

			numerator+=Ri_div_Rmax*ui;
			denominator+=Ri_div_Rmax;
		}
        MADNESS_ASSERT(!isnan(numerator/denominator));
        return numerator/denominator;
	}
};



/// functor for the diamagnetic term in a box

/// The diamagnetic term is a 2D harmonic potential, which makes the iterations diverge.
/// Approximate it by making it constant at a specific radius: construct a smooth approximation
/// for the piecewise function f(x)={{x, x<1}, {1,x>1}}, which is squared. For this approximation
/// see https://doi.org/10.1186/s40064-016-3278-y
struct harmonic_potential_boxed : public FunctionFunctorInterface<double,3> {
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

void Diamagnetic_potential_factor::print_info() const {

	print("\ninformation about the diamagnetic_potential_factor\n");
	print("coords in A space (v) ",v);
	print("coords in real space  ",coords);
	print("physical magnetic field B      ",physical_B);
	print("explicit magnetic field B      ",explicit_B);
	print("");

}


void Diamagnetic_potential_factor::recompute_factors_and_potentials() {
	diamagnetic_factor_=custom_factor(explicit_B,v,1.0);
	diamagnetic_factor_square=custom_factor(explicit_B,v,2.0);
	diamagnetic_U2=compute_U2();
	diamagnetic_U1=-1.0*compute_nabla_R_div_R();
//	save(diamagnetic_U2,"U2");
//	save(diamagnetic_factor_,"factor");
//	save(diamagnetic_U1[0],"U1x");
//	save(diamagnetic_U1[1],"U1y");
//	save(diamagnetic_U1[2],"U1z");
}


complex_function_3d Diamagnetic_potential_factor::compute_lz_commutator() const {
	std::vector<real_function_3d> r(3), A(3);
    r[0]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
    r[1]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});
    r[2]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[2];});

    A[0]=explicit_B[1]*r[2]-explicit_B[2]*r[1];
    A[1]=explicit_B[2]*r[0]-explicit_B[0l]*r[2];
    A[2]=explicit_B[0l]*r[1]-explicit_B[1]*r[0];

	scale(world,A,0.5);
	complex_function_3d result=double_complex(0.0,1.0)*dot(world,A,diamagnetic_U1);
	return result;
}


real_function_3d Diamagnetic_potential_factor::compute_U2() const {

	MADNESS_ASSERT(B_along_z(explicit_B));
	MADNESS_ASSERT(B_along_z(physical_B));
	double eB=explicit_B.normf();
	double pB=physical_B.normf();

	std::function<double(const coord_3d& r, const coord_3d& B, const coord_3d& v)>
	U2c = [] (const coord_3d& r, const coord_3d& B, const coord_3d& v) {
		coord_3d Amv=0.5*cross(B,r)-v;
		double Amv_square=inner(Amv,Amv);
		double ui=Amv_square;
		return ui;
	};

	real_function_3d U2c_function(world);
	double radius=get_potential_radius();
	if (eB>0.0) U2c_function=real_factory_3d(world).functor(R_times_arg_div_R(U2c,explicit_B,v,radius)).truncate_on_project();
	real_function_3d harmonic_potential=real_factory_3d(world).functor(harmonic_potential_boxed(radius));

	real_function_3d U2_function;
	U2_function=-0.5*(-eB - 0.25*pB*pB*harmonic_potential  + U2c_function);
	return U2_function.truncate();
}

real_function_3d Diamagnetic_potential_factor::compute_R_times_T_commutator_scalar_term_numerically() const {
	// the physical potential
	MADNESS_ASSERT(B_along_z(physical_B));
	double pB=physical_B.normf();
	double radius=get_potential_radius();
	real_function_3d harmonic_potential=real_factory_3d(world).functor(harmonic_potential_boxed(radius));
	harmonic_potential.scale(0.125*pB*pB) ;

	return (-0.5*div(grad_bspline_one(factor())) + harmonic_potential*factor()).truncate();
}


std::vector<real_function_3d> Diamagnetic_potential_factor::compute_nabla_R_div_R() const {

	std::vector<real_function_3d> result=zero_functions_compressed<double,3>(world,3);

	if (explicit_B.normf()>0.0) {
		MADNESS_ASSERT(B_along_z(physical_B));
		for (int idim=0; idim<3; ++idim) {

			std::function<double(const coord_3d& r, const coord_3d& B, const coord_3d& v)>
			U1 = [&idim] (const coord_3d& r, const coord_3d& B, const coord_3d& v) {
				const double Bnorm=B.normf();
				const double B2=inner(B,B);
				const double rB=inner(r,B);
				const coord_3d Bxv=cross(B,v);
				return (-1.0/(2.0*Bnorm))*((B2*r[idim] - B[idim]*rB) + 2.0*Bxv[idim]);
			};

			// no flattening of the U1 potential: use cell diameter as radius
			double radius=FunctionDefaults<3>::get_cell().normf();
			result[idim]=real_factory_3d(world).functor(R_times_arg_div_R(U1,explicit_B,v,radius));
		}
		truncate(world,result);
	}
	return result;
}


/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
std::vector<complex_function_3d> Diamagnetic_potential_factor::apply_potential(const std::vector<complex_function_3d>& rhs) const {

	std::vector<complex_function_3d> result=zero_functions_compressed<double_complex,3>(world,rhs.size());
	if (physical_B.normf()==0.0) return result;

	double plateau=diamagnetic_U2(20.0, 20.0, 20.0);

	result=(diamagnetic_U2-plateau)*rhs + plateau*rhs;
	if (explicit_B.normf()!=0.0) {
		Derivative<double_complex,3> dx(world,0);
		Derivative<double_complex,3> dy(world,1);
		std::vector<complex_function_3d> dxrhs=apply(world,dx,rhs);
		std::vector<complex_function_3d> dyrhs=apply(world,dy,rhs);

		std::vector<complex_function_3d> tmp=(diamagnetic_U1[0]*dxrhs + diamagnetic_U1[1]*dyrhs);
		result+=tmp;
	}
	truncate(world,result);
	return result;
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
		if (absB==0.0) return 1.0;
		for (auto& cc : c) {
			coord_3d arg1=(r-cc);
			result+=exp(-0.25/absB*(inner(B,B)*inner(arg1,arg1) - inner(B,arg1)*inner(B,arg1)));
		}
		return result;
	};

	auto factor_square = [&B, &c](const coord_3d& r) {
		double result=0.0;
		double absB=B.normf();
		if (absB==0.0) return 1.0;
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
	t.logger << "diffnorm1 " << diffnorm1 << std::endl;
	t.logger << "diffnorm2 " << diffnorm2 << std::endl;
	bool success=((diffnorm1<1.e-10) and (diffnorm2<1.e-10));
	return t.end(success);
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


bool Diamagnetic_potential_factor::test_scalar_potentials() const {

	test_output t("entering Diamagnetic_potential_factor::test_scalar_potentials .. ");
	const std::vector<complex_function_3d> mo=make_fake_orbitals(2);
	const std::vector<complex_function_3d> diamo=factor()*mo;
	const std::vector<complex_function_3d> dia2mo=factor_square()*mo;

	real_function_3d T_commutator_scalar_term=compute_U2();
	real_function_3d T_commutator_scalar_term_numerical=compute_R_times_T_commutator_scalar_term_numerically();

	double_complex t3=inner(dia2mo,T_commutator_scalar_term*mo);
	double_complex t3a=inner(diamo,T_commutator_scalar_term_numerical*mo);

	t.logger << "t3   " << t3 << std::endl;
	t.logger << "t3a  " << t3a << std::endl;
	t.logger << "error " << std::abs(t3-t3a) << std::endl;

	bool success=(std::abs(t3-t3a)<FunctionDefaults<3>::get_thresh());
	return t.end(success);
}

bool Diamagnetic_potential_factor::test_vector_potentials() const {

	test_output t("entering Diamagnetic_potential_factor::test_vector_potentials .. ");
	double thresh=FunctionDefaults<3>::get_thresh();

	const std::vector<real_function_3d> U1_analytical=compute_nabla_R_div_R();
	const std::vector<real_function_3d> U1_numerical_times_factor=grad(factor());

	double total_error=0.0;
	for (int i=0; i<3; ++i) {
		double norm=U1_analytical[i].norm2();
		double error=(U1_analytical[i]*factor()-U1_numerical_times_factor[i]).norm2();
		double relative_error= (error>thresh) ? error/norm : 0.0;	// norm might be zero
		t.logger << "error[" << i << "]" << relative_error << std::endl;
		total_error+=relative_error;
	}

	bool success=(total_error<thresh);
	return t.end(success);
}


bool Diamagnetic_potential_factor::test_lz_commutator() const {

	test_output t("entering Znemo::test_lz_commutator .. ");
	double thresh=FunctionDefaults<3>::get_thresh();

	Lz<double_complex,3> lz(world);

	complex_function_3d lzR=0.5*explicit_B[2]*lz(convert<double,double_complex,3>(factor()));
	complex_function_3d lz_commR=compute_lz_commutator()*factor();
	double norm1=lzR.norm2();
	double norm2=lz_commR.norm2();
	double error=(lzR-lz_commR).norm2();

	t.logger << "norm lz numerical  " << norm1 << std::endl;
	t.logger << "norm lz analytical " << norm2 << std::endl;
	t.logger << "error lz numerical vs analytical " << error << std::endl;

	return t.end(error<thresh);
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
