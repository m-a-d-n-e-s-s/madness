/*
 * BasisFunctions.cc
 *
 *  Created on: Apr 18, 2018
 *      Author: kottmanj
 */

#include <PNOGuessFunctions.h>
#include <cmath>
#ifdef MADNESS_HAS_BOOST
#include <boost/math/special_functions/legendre.hpp>
#endif

namespace madness {

 vector_real_function_3d BasisFunctions::guess_virtuals_from_file() const {
	MyTimer time_1 = MyTimer(world).start();
	// Read Exponents from file
	std::map<std::string, std::vector<std::vector<double> > > exponents = read_basis_from_file("bas", molecule.get_atoms());
	print("Exponents from file bas:");
	for (const auto& x : exponents) {
		if (world.rank() == 0) {
			std::cout << x.first << " Exponents\n";
			for (size_t l = 0; l < x.second.size(); ++l) {
				std::cout << "l=" << l << "\n";
                                using madness::operators::operator<<;
				std::cout << x.second[l] << "\n";
			}
		}
	}
	time_1.stop().print("Read Exponents From File");
	vector_real_function_3d virtuals;
	MyTimer time_2 = MyTimer(world).start();
	for (const madness::Atom& atom : molecule.get_atoms()) {
		std::vector<std::vector<double> > exp_atom = exponents.at(atomic_number_to_symbol(atom.atomic_number));
		for (size_t l = 0; l < exp_atom.size(); ++l) {
			for (const double e : exp_atom[l]) {
				virtuals=append(virtuals, guess_virtual_gaussian_shell(atom, l, e));
			}
		}
	}
	time_2.stop().print("Creating Guess Basis");
	return virtuals;
}

 vector_real_function_3d BasisFunctions::guess_virtuals_internal(const std::map<std::string, std::vector<int> > guess_map) const {
	vector_real_function_3d virtuals;
	for (size_t i = 0; i < molecule.natom(); ++i) {
		// get the partial wave basis for the atom
		const Atom& atom = molecule.get_atom(i);
		const std::string symbol = atomic_number_to_symbol(atom.atomic_number);
		std::vector<int> pw_guess;
		if (guess_map.find(symbol) != guess_map.end()) {
			pw_guess = guess_map.find(symbol)->second;
		} else {
			if (world.rank() == 0)
				print("symbol", symbol);

			//if (world.rank() == 0) print("pw_guess",param.pw_guess);
			MADNESS_EXCEPTION("could not find basis for atom", 1);
		}
		if (world.rank() == 0) {
			print("working on atom", symbol);
			print("pw_guess", pw_guess);
		}
		for (size_t l = 0; l < pw_guess.size(); ++l) {
			for (int i = 0; i < pw_guess[l]; ++i) {
				double e = double(pw_guess[l]) / (double(i + 1.0));
				virtuals=append(virtuals, guess_virtual_gaussian_shell(atom, l, e));
				if (world.rank() == 0)
					print("l,e", l, e);

				if (world.rank() == 0)
					print("virtuals.size() = ", virtuals.size());
			}
		}
	}
	if (world.rank() == 0)
		print("number of guess virtuals: ", virtuals.size());

	// plot guess
	for (std::size_t i = 0; i < virtuals.size(); ++i) {
		std::string name = "virtuals" + stringify(i);
		plot_plane(world, virtuals[i], name);
	}
	return virtuals;
}

 vector_real_function_3d BasisFunctions::guess_virtual_gaussian_shell(const Atom& atom, const int l, const double e) const {
	vector_real_function_3d virtuals;
	auto gg = make_solidharmonic_guess(atom, l, e);
	for (std::size_t m = 0; m < gg.size(); ++m) {
		virtuals.push_back(real_factory_3d(world).functor(gg[m]).truncate_on_project());
	}
	normalize(world, virtuals);
	return virtuals;
}

 std::vector<std::vector<double> > BasisFunctions::read_basis_from_file(const std::string& filename, const std::string& atom) const {
	std::ifstream f(filename.c_str());
	position_stream(f, atom);
	std::string s;
	std::string str_pw_guess, symbol;
	std::vector<std::vector<double> > result(8);
	std::size_t stars = 0;
	while (f >> s) {
		if (s == "*") {
			++stars;
			if (stars == 2)
				break;
		}
		if (s == "s" || s == "p" || s == "d" || s == "f" || s == "g" || s == "h" || s == "i" || s == "k") {
			double exponent;
			f >> exponent;
			result[lqtoint(s)].push_back(exponent);
		}
	}
	return result;
}

 std::map<std::string, std::vector<std::vector<double> > > BasisFunctions::read_basis_from_file(const std::string& filename, const std::vector<madness::Atom> atoms) const {
	std::map<std::string, std::vector<std::vector<double> > > cabs_exponents;
	for (const madness::Atom& atom : atoms) {
		const std::string symbol = atomic_number_to_symbol(atom.atomic_number);
		if (cabs_exponents.find(symbol) == cabs_exponents.end()) {
			std::vector<std::vector<double> > tmpbas = read_basis_from_file(filename, symbol);
			cabs_exponents[symbol] = tmpbas;
		}
	}
	return cabs_exponents;
}

 size_t BasisFunctions::lqtoint(const std::string& l) const {
	if (l == "s")
		return 0;
	else if (l == "p")
		return 1;
	else if (l == "d")
		return 2;
	else if (l == "f")
		return 3;
	else if (l == "g")
		return 4;
	else if (l == "h")
		return 5;
	else if (l == "i")
		return 6;
	else if (l == "k")
		return 7;
	else {
		MADNESS_EXCEPTION(("l-quantum number:" + l + " not supported").c_str(), 1);
		return 99999999;
	}

	MADNESS_EXCEPTION("Should not end up here", 1);
}



 double BasisFunctions::SolidHarmonicGaussian::operator ()(const coord_3d& xyz) const {
#ifdef MADNESS_HAS_BOOST
        static const int64_t fac[] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000 };
	static const auto sqrt_2 = 1.41421356237309504880168872421;
	double xx = xyz[0] - x;
	double yy = xyz[1] - y;
	double zz = xyz[2] - z;
	const double r2 = (xx * xx + yy * yy + zz * zz);
	const double r = std::sqrt(r2);
	const double cos_theta = zz / r;
	const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
	const double cos_phi = xx / (r * sin_theta);
	const double phi = (yy > 0 ? 1.0 : -1.0) * std::acos(cos_phi);
	const auto abs_m = std::abs(m);
	assert(size_t(l + abs_m) <= (sizeof(fac) / sizeof(int64_t) - 1));
	// wrong sign for m < 0
	//const auto P_l_m = boost::math::assoc_legendre(l, abs_m, cos_theta);
	const auto P_l_m = boost::math::legendre_p(l, abs_m, cos_theta);
	// this excludes sqrt((2l+1)/4pi) since that gets cancelled by its inverse in the definition of the solid harmonics
	// this also excludes (-1)^m since that is also included in the real spherical harmonics phase
	const auto Y_normconst = std::sqrt(static_cast<double>(fac[l - m]) / static_cast<double>(fac[l + m]));
	// (-1)^m was cancelled in Y_normconst
	const auto real_azimuthal_prefactor = (m != 0) ? sqrt_2 * (m > 0 ? std::cos(m * phi) : std::sin(abs_m * phi)) : 1.0;
	return pow(r, l) * (Y_normconst * real_azimuthal_prefactor * P_l_m) * exp(-exponent * r2);
#else
	MADNESS_EXCEPTION("can not create SolidHarmonicGaussian without boost.math library, compile MADNESS with boost by using cmake flag -D ENABLE_BOOST=ON ",1);
	return 0.0;
#endif


}

} /* namespace madness */
