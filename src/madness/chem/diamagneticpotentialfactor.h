/*
 * diamagneticpotentialfactor.h
 *
 *  Created on: 16 Apr 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_DIAMAGNETICPOTENTIALFACTOR_H_
#define SRC_APPS_CHEM_DIAMAGNETICPOTENTIALFACTOR_H_


#include <madness/mra/mra.h>
#include<madness/chem/znemo.h>

namespace madness {


/// to be put in a separate file
class Diamagnetic_potential_factor {

public:

	/// constructor takes a world and the parameters for the calculation
	Diamagnetic_potential_factor(World& world, const Nemo_complex_Parameters& param,
			const std::vector<coord_3d>& coords)
		: world(world), coords(coords) {

		use_v_vector=param.use_v_vector();
		physical_B={0,0,param.physical_B()};
		explicit_B={0,0,param.explicit_B()};

		v=compute_v_vector(explicit_B,coords,use_v_vector);

		potential_radius=param.potential_radius();
		cut_A=param.cut_vector_potential();
		if (cut_A) {
			MADNESS_CHECK_THROW(param.explicit_B()==0.0,"cut_vector_potential needs explicit_B = 0");
			cut_width=param.potential_cut_width();
		}

		// initialize all functions to their default values
		diamagnetic_factor_=real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});
		diamagnetic_factor_square=real_factory_3d(world).functor([] (const coord_3d& r) {return 1.0;});
		diamagnetic_U2=real_factory_3d(world).functor([] (const coord_3d& r) {return 0.0;});
		diamagnetic_U1=zero_functions<double,3>(world,3);

		if (physical_B.normf()>0.0) recompute_factors_and_potentials();
		if (param.printlevel()>9) print_info();
	}

	void reset_explicit_B_and_v(const coord_3d& eB) {
		explicit_B=eB;
		v=compute_v_vector(explicit_B,coords,use_v_vector);
	}

	static std::vector<coord_3d> compute_v_vector(const coord_3d& B,
			const std::vector<coord_3d>& coords, bool use_v_vector) {
		std::vector<coord_3d> v;
		if (use_v_vector) {
			for (auto& c : coords) v.push_back(0.5*cross(B,c));
		} else {
			v=std::vector<coord_3d> (1,{0.0,0.0,0.0});
		}
		return v;
	}


	void print_info() const;

	/// return the diamagnetic factor
	real_function_3d factor() const {return diamagnetic_factor_;}

	/// return the square of the diamagnetic factor
	real_function_3d factor_square() const {return diamagnetic_factor_square;}

	/// return a custom factor for a given magnetic field
	real_function_3d custom_factor(const coord_3d& B, const std::vector<coord_3d>& vv,
			const double extra_exponent=1.0) const {
		const double absB=B.normf();
		if (absB==0.0) return real_factory_3d(world).functor([](const coord_3d& r){return 1.0;}).truncate_on_project();
		const double& ee=extra_exponent;
		auto diamagnetic_HO = [&absB, &B, &vv, &ee](const coord_3d& r) {
			double result=0.0;
			const coord_3d A=0.5*cross(B,r);
			for (auto& v : vv) {
				coord_3d arg=A-v;
				result+=exp(-1.0/absB*inner(arg,arg));
			}
			if (ee!=1.0) result=std::pow(result,ee);
			return result;
		};
		real_function_3d result=real_factory_3d(world).functor(diamagnetic_HO);
		return result;
	}


	/// return a custom factor for a given magnetic field
	complex_function_3d factor_with_phase(const coord_3d& B, const std::vector<coord_3d>& vv) const {
		const double absB=B.normf();
		if (absB==0.0) return complex_factory_3d(world).functor([](const coord_3d& r){return 1.0;}).truncate_on_project();
		auto diamagnetic_HO = [&absB, &B, &vv](const coord_3d& r) {
			double result=0.0;
			const coord_3d A=0.5*cross(B,r);
			for (auto& v : vv) {
				coord_3d arg=A-v;
				result+=exp(-1.0/absB*inner(arg,arg));
			}
			double theta=acos(r[2]/r.normf());
			double phi=atan2(r[1],r[0]);
			double_complex phase=r.normf()*sin(theta)*exp(absB*double_complex(0.0,1.0)*phi);
			return result*phase;
		};
		complex_function_3d result=complex_factory_3d(world).functor(diamagnetic_HO);
		return result;
	}

	/// compute the bare potential without confinement or factors
	real_function_3d bare_diamagnetic_potential() const {
		MADNESS_ASSERT(B_along_z(physical_B));
		const double b_square=inner(physical_B,physical_B);
		real_function_3d result=real_factory_3d(world)
				.functor([& b_square](const coord_3d& r){return 0.125*b_square*(r[0]*r[0] + r[1]*r[1]);});
		return result;
	}

	/// apply the diamagnetic potential on rhs
	std::vector<complex_function_3d> apply_potential(const std::vector<complex_function_3d>& rhs) const;

	/// field profile of the cut vector potential: B(rho) = B F(rho), F = erfc((rho - R_c)/w)/2
	static double cut_field_profile(const double rho, const double R_c, const double w) {
		return 0.5*std::erfc((rho-R_c)/w);
	}

	/// flux-conserving cut of the symmetric-gauge vector potential, A_c = S(rho) A

	/// A_c is the vector potential of the field B F(rho) about the gauge origin (rho = distance from the z axis):
	/// S(rho) = (2/rho^2) int_0^rho F(rho') rho' drho', -> 1 inside R_c, -> R_c^2/rho^2 outside (a pure-gauge
	/// Aharonov-Bohm potential).  In closed form, by parts, with x = rho' - R_c:
	///   int_0^rho F rho' drho' = rho^2 F(rho)/2 + 1/(2 w sqrt(pi)) int_{-R_c}^{rho-R_c} (x + R_c)^2 exp(-x^2/w^2) dx
	static double cut_A_factor(const double rho, const double R_c, const double w) {
		const double sqrtpi=std::sqrt(constants::pi);
		// small rho: the closed form cancels; S = F(0) + (2/3) F'(0) rho + O(rho^2)
		if (rho<0.01*w) return cut_field_profile(0.0,R_c,w)-(2.0/3.0)*rho*std::exp(-R_c*R_c/(w*w))/(w*sqrtpi);
		const double a=-R_c, b=rho-R_c;
		const double ea=std::exp(-a*a/(w*w)), eb=std::exp(-b*b/(w*w));
		const double E0=0.5*w*sqrtpi*(std::erf(b/w)-std::erf(a/w));	// int exp(-x^2/w^2)
		const double E1=-0.5*w*w*(eb-ea);								// int x exp(-x^2/w^2)
		const double E2=-0.5*w*w*(b*eb-a*ea)+0.5*w*w*E0;				// int x^2 exp(-x^2/w^2)
		const double I=0.5*rho*rho*cut_field_profile(rho,R_c,w)+(E2+2.0*R_c*E1+R_c*R_c*E0)/(2.0*w*sqrtpi);
		return 2.0*I/(rho*rho);
	}

	/// is the vector potential cut (cut_vector_potential)?
	bool cut_vector_potential() const {return cut_A;}

	/// S(rho) of the cut vector potential as a function (1 if the vector potential is not cut)
	real_function_3d cut_A_factor_function() const {
		if (not cut_A) return real_factory_3d(world).functor([](const coord_3d& r){return 1.0;});
		const double R_c=potential_radius, w=cut_width;
		return real_factory_3d(world).functor([R_c,w](const coord_3d& r) {
			return cut_A_factor(std::sqrt(r[0]*r[0]+r[1]*r[1]),R_c,w);
		});
	}

	/// diamagnetic potential of the cut vector potential, A_c^2/2 = (B^2/8) S^2 rho^2, plus the exterior
	/// offset (|B|/2)(1 - F), which lifts the field-free exterior to the Landau zero point (for the bsh kernel;
	/// without it the exterior continuum lies below states bound under the lowest Landau level)
	real_function_3d cut_A_diamagnetic_potential() const {
		MADNESS_CHECK_THROW(cut_A,"cut_A_diamagnetic_potential: the vector potential is not cut");
		const double R_c=potential_radius, w=cut_width, B=physical_B.normf();
		return real_factory_3d(world).functor([R_c,w,B](const coord_3d& r) {
			const double rho=std::sqrt(r[0]*r[0]+r[1]*r[1]);
			const double S=cut_A_factor(rho,R_c,w);
			return 0.125*B*B*S*S*rho*rho+0.5*B*(1.0-cut_field_profile(rho,R_c,w));
		});
	}

	/// make sure the magnetic field is oriented along the z axis
	static bool B_along_z(const coord_3d& B) {
		return((B[0]==0.0) && (B[1]==0.0));
	}

	std::vector<coord_3d> get_v() const {return v;}
	std::vector<coord_3d> get_coords() const {return coords;}

	coord_3d get_explicit_B() const {return explicit_B;}

	coord_3d get_physical_B() const {return physical_B;}

	/// given the explicit and the physical B, estimate the radius
	/// of the wave function
	double estimate_wavefunction_radius(const double eps=1.e-8) const {
		double remaining_B=(get_physical_B()-get_explicit_B()).normf();
		return 2.0*sqrt(-log(eps)/remaining_B);
	}

private:
	World& world;				///< the world
	coord_3d physical_B={0,0,0};		///< the actual magnetic field strength
	coord_3d explicit_B={0,0,0};		///< the magnetic field strength encoded in the diamagnetic factor

	/// the diamagnetic factor to cancel the diamagnetic potential
	real_function_3d diamagnetic_factor_;
	real_function_3d diamagnetic_factor_square;

	/// radius where the diamagnetic potential flattens out, or (cut_A) where the field is cut
	double potential_radius;

	/// flux-conserving cut of the vector potential at potential_radius, edge width cut_width
	bool cut_A=false;
	double cut_width=0.0;

	/// the boxed diamagnetic potential (for a given B)
	real_function_3d diamagnetic_U2;
	std::vector<real_function_3d> diamagnetic_U1;

	/// the position of the nuclei in the coordinate space:
	std::vector<coord_3d> coords;

	/// the position of the nuclei in the "A" space: v = 1/2 B cross R
	std::vector<coord_3d> v;

	bool use_v_vector=true;

public:
	/// recompute the factor and the potentials for given physical and explicit magnetic fields
	void recompute_factors_and_potentials();

	/// compute the commutator of the orbital-zeeman term with the diamagnetic factor

	/// @return	a local potential: iB \sum_i \vec r \cdot \vec v_i
	complex_function_3d compute_lz_commutator() const;

	real_function_3d compute_U2() const;
	real_function_3d compute_R_times_T_commutator_scalar_term_numerically() const;

	/// returns \f$R^{-1} \vec\nabla R\f$
	std::vector<real_function_3d> compute_nabla_R_div_R() const;

	/// run the tests

	/// @param[in]	level: 1: basic tests, .. , 3: all the tests
	bool test_me(const int level) const {
		bool success=true;
		if ((level<1) or (level>3)) {
			print("testing diamagneticpotentialfactor with an invalid level of ",level);
		}
		if (level<2) {
			success=test_factor() and success;
			success=test_lz_commutator() and success;
			success=test_harmonic_potential() and success;
			success=test_scalar_potentials() and success;
			success=test_vector_potentials() and success;
		}
		if (level<3) {
			;
		}
		if (level<4) {
			;
		}
		return success;
	}

private:


	/// compute a factor for comparison in coordinate space
	bool test_factor() const;

	/// test the harmonic potential
	bool test_harmonic_potential() const;

	/// test analytical vs numerical computation of the potentials
	bool test_scalar_potentials() const;

	/// test analytical vs numerical computation of the potentials
	bool test_vector_potentials() const;

	bool test_lz_commutator() const;

public:
	/// make a set orbitals for testing (not orthonormalized!)
	std::vector<complex_function_3d> make_fake_orbitals(const int n,
			const coord_3d& offset={0.0,0.0,0.0}) const;

	/// compute the radius for the diamagnetic potential
	double get_potential_radius() const {
		return potential_radius;
	}

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_DIAMAGNETICPOTENTIALFACTOR_H_ */
