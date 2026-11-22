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

		physical_B={0,0,param.physical_B()};
		potential_radius=param.potential_radius();
		cut_A=param.cut_vector_potential();
		if (cut_A) cut_width=param.potential_cut_width();

		diamagnetic_U2=real_factory_3d(world).functor([] (const coord_3d& r) {return 0.0;});
		if (physical_B.normf()>0.0) recompute_factors_and_potentials();
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

	std::vector<coord_3d> get_coords() const {return coords;}

	coord_3d get_physical_B() const {return physical_B;}

private:
	World& world;				///< the world
	coord_3d physical_B={0,0,0};		///< the actual magnetic field strength

	/// radius where the diamagnetic potential flattens out, or (cut_A) where the field is cut
	double potential_radius;

	/// flux-conserving cut of the vector potential at potential_radius, edge width cut_width
	bool cut_A=false;
	double cut_width=0.0;

	/// the boxed diamagnetic potential (for a given B)
	real_function_3d diamagnetic_U2;

	/// the position of the nuclei in the coordinate space:
	std::vector<coord_3d> coords;

public:
	/// recompute the diamagnetic potential for the physical magnetic field
	void recompute_factors_and_potentials();

	/// the (boxed or cut) diamagnetic potential of the bsh iteration
	real_function_3d compute_U2() const;

	/// run the tests

	/// @param[in]	level: 1: basic tests, .. , 3: all the tests
	bool test_me(const int level) const {
		bool success=true;
		if ((level<1) or (level>3)) {
			print("testing diamagneticpotentialfactor with an invalid level of ",level);
		}
		if (level<2) {
			success=test_harmonic_potential() and success;
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


	/// test the harmonic potential
	bool test_harmonic_potential() const;

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
