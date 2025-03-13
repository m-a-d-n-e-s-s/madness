/*
 * BasisFunctions.h
 *
 *  Created on: Apr 18, 2018
 *      Author: kottmanj
 */

#ifndef PAPER_CODE_BASISFUNCTIONS_H_
#define PAPER_CODE_BASISFUNCTIONS_H_


#include <exception>
#include <stdlib.h>
#include <cassert>
#include <cmath>
#include <numeric>
#include <valarray>
#include <stdio.h>
#include <madness.h>
#include<madness/chem/PNOStructures.h>
#include<madness/chem/projector.h>
#include<madness/chem/GuessFactory.h>
#include<madness/chem/molecule.h>

#include <tuple>


namespace madness {


/// class which provides guess functions for pnos and cabs basis sets
/// the guess functions are not orthonormalized or Q-projected
class BasisFunctions{
public:
	typedef std::vector<std::tuple<int, std::vector<double>, std::vector<double> > > cbfT;
	BasisFunctions(World& world,const Molecule& mol, const size_t& l): world(world), molecule(mol), lmax(l) {}
	World& world;
	const Molecule& molecule;
	const size_t lmax;

	// print out a basis of contracted orbitals
	void print_contracted_basis(std::map<std::string, cbfT >& molbas)const{
		for(const auto& tmp:molbas){
		cbfT abf=tmp.second;
		for(const auto& bf:abf){
			const int type = std::get<0>(bf);
			const std::vector<double> ex=std::get<1>(bf);
			const std::vector<double> c=std::get<2>(bf);
			std::cout << ex.size() << " " << type << "\n";
			for(size_t k=0;k<ex.size();++k) std::cout << ex[k] << "  " << c[k] << "\n";

		}
		}
	}

	vector_real_function_3d guess_virtuals_from_file() const;
	vector_real_function_3d guess_contracted_virtuals_from_file() const {
		MyTimer time_1 = MyTimer(world).start();
		// Read Exponents from file
		std::map<std::string, cbfT > molbas = read_contracted_basis_from_file("bas", molecule.get_atoms());
		print("Exponents from file bas:");
		print_contracted_basis(molbas);
		time_1.stop().print("Read Exponents From File");
		vector_real_function_3d virtuals;

		MyTimer time_2 = MyTimer(world).start();
		for (const madness::Atom& atom : molecule.get_atoms()) {
			cbfT abf = molbas.at(atomic_number_to_symbol(atom.atomic_number));
			for(const auto& bf:abf){
			const int type = std::get<0>(bf);
			const std::vector<double> ex=std::get<1>(bf);
			const std::vector<double> c=std::get<2>(bf);
			MADNESS_ASSERT(ex.size()==c.size());
			vector_real_function_3d contracted_function;
			for(size_t k=0;k<ex.size();++k){
				vector_real_function_3d gshell = guess_virtual_gaussian_shell(atom, type, ex[k]);
				print_size(world,gshell,"ex="+std::to_string(ex[k]));
				normalize(world,gshell);
				truncate(world,gshell,FunctionDefaults<3>::get_thresh());
				print_size(world,gshell,"ex="+std::to_string(ex[k]));
				if(contracted_function.empty()) contracted_function=c[k]*gshell;
				else contracted_function += c[k]*gshell;
			}
			normalize(world,contracted_function);
			virtuals=append(virtuals,contracted_function);

			}
//			for (size_t l = 0; l < exp_atom.size(); ++l) {
//				for (const double e : exp_atom[l]) {
//					virtuals=append(virtuals, guess_virtual_gaussian_shell(atom, l, e));
//				}
//			}
		}
		time_2.stop().print("Creating Contracted Guess Basis");
		return virtuals;
	}

	vector_real_function_3d guess_with_psi4(const vector_real_function_3d& mos)const{
		MADNESS_EXCEPTION("Not there for new madness version",1);
	}

	/// make guess virtuals by exciting with polynomials: v = poly*f
	vector_real_function_3d guess_with_exop(const vector_real_function_3d& f, const std::string& type="dipole+", const bool trigo=true)const{
		if(world.rank()==0) std::cout << "Create guess functions by multiplying plane-waves to " << f.size() << " functions \n";
		MADNESS_ASSERT(not f.empty());

		std::vector<coord_3d> centers = guessfactory::compute_centroids(f);
		std::vector<std::string> exop_list =  guessfactory::make_predefined_exop_strings(type);



		// make the list with seeds and excitation operators
		std::vector<std::pair<vector_real_function_3d,std::string> > init_list;

		for(const auto& exop:exop_list){
			vector_real_function_3d cs=copy(world,f,false);
			init_list.push_back(std::make_pair(cs,exop));
		}
		world.gop.fence();

		vector_real_function_3d virtuals;
		if(trigo) for(auto& it:init_list) virtuals=append(virtuals, guessfactory::apply_trigonometric_exop(it.first,it.second,centers,false));
		else for(auto& it:init_list){
			if(world.rank()==0) std::cout << "polynomial guess!\n";
			virtuals=append(virtuals, guessfactory::apply_polynomial_exop(it.first,it.second,centers,false));
		}
		world.gop.fence();

		return (virtuals);
	}

	vector_real_function_3d guess_virtuals_internal(const std::map<std::string, std::vector<int> > guess_map) const;
	vector_real_function_3d predefined_guess(const std::string name)const{
		// determine start zeta value
		size_t start=0;
		bool empty=false;
		if(name=="pvdz") start=2;
		else if(name=="pvtz") start=3;
		else if(name=="pvqz") start=4;
		else if(name=="pv5z") start=5;
		else if(name=="pv6z") start=6;
		else if(name=="pv7z") start=7;
		else if(name=="zero" or name=="0" or name=="void" or name=="none" or name=="empty") empty=true;
		else MADNESS_EXCEPTION(("unknown predefined guess:"+name).c_str(),1);

		std::map<std::string, std::vector<int> > guess_map;
		for(const auto& atom :molecule.get_atoms()){// (int i = 0; i < nemo.molecule().natom(); ++i) {
			const std::string symbol = atomic_number_to_symbol(atom.atomic_number);
			const std::size_t n=atom.atomic_number;
			if(empty) guess_map[symbol]=std::vector<int>(lmax,0);
			else{
				if(n<3) guess_map[symbol]=fill_peterson(start); // first row
				else if(n<11) guess_map[symbol]=fill_peterson(start+1); // first row
				else if(n<19) guess_map[symbol]=fill_peterson(start+2); // third row
				else MADNESS_EXCEPTION("Predefined guesses only up to third row",1);
			}

		}
		return guess_virtuals_internal(guess_map);
	}

	// helper to fill up l numbers (Peterson Style)
	std::vector<int> fill_peterson(const size_t& s)const{
		std::vector<int> result(lmax,0);
		for(size_t i=0;i<lmax;i++){
			const int tmp=(s-i);
			if(tmp>0) result[i]=tmp;
			else result[i]=0;
		}
		return result;
	}

	// helper to get vectors of right size
	std::vector<int> fill_up(std::vector<int> v)const{
		std::vector<int> result(lmax,0);
		for(size_t i=0;i<(lmax>v.size() ? v.size() : lmax);++i) result[i]=v[i];
		return result;
	}

	/// return a shell of l-quantum l and exponent e, including all cartesian
	/// components
	vector_real_function_3d guess_virtual_gaussian_shell(const Atom& atom, const int l, const double e) const;

	/// read external CABS
	std::map<std::string, std::vector<std::vector<double> > > read_basis_from_file(const std::string& filename, const std::vector<madness::Atom> atoms) const;
	std::vector<std::vector<double> > read_basis_from_file(const std::string& filename, const std::string& atom) const;




	std::map<std::string, cbfT > read_contracted_basis_from_file(const std::string& filename, const std::vector<madness::Atom> atoms) const {
		std::map<std::string, cbfT > molbas;
		for (const madness::Atom& atom : atoms) {
			const std::string symbol = atomic_number_to_symbol(atom.atomic_number);
			if (molbas.find(symbol) == molbas.end()) {
				cbfT tmpbas = read_contracted_basis_from_file(filename, symbol);
				molbas[symbol] = tmpbas;
			}
		}
		return molbas;
	}

	cbfT read_contracted_basis_from_file(const std::string& filename, const std::string& atom) const {
		std::ifstream f(filename.c_str());
		position_stream(f, atom);
		std::string s;
		int n=-1;
		std::string str_pw_guess, symbol;
		cbfT result;
		std::size_t stars = 0;
		while (f >> s) {
			if (s == "*") {
				f>>n; // initial step
				++stars;
				if (stars == 2)
					break;
			}else{
				n=std::stoi(s);
				MADNESS_ASSERT(n>0);
			}
			if(n>0){
				f>>s;
				if (s == "s" || s == "p" || s == "d" || s == "f" || s == "g" || s == "h" || s == "i" || s == "k") {
					std::vector<double> exponents;
					std::vector<double> coeff;
					for(int i=0;i<n;++i){
						double ex,c;
						f>>ex;
						f>>c;
						exponents.push_back(ex);
						coeff.push_back(c);
					}

					result.push_back(std::make_tuple(lqtoint(s),exponents,coeff));
				} else MADNESS_EXCEPTION("UNKNOWN ANGULAR QUANTUM NUMBER",1);
			}
		}

		return result;
	}

	/// little helper function for l-quantum numbers
	size_t lqtoint(const std::string& l) const;

	class CartesianGaussian : public FunctionFunctorInterface<double, 3> {
	public:
		CartesianGaussian(const Atom& atom, const double e, const int i, const int j, const int k);

		double x, y, z;   ///< origin
		double exponent;  ///< exponent
		int i, j, k;      ///< cartesian exponents

		double operator()(const coord_3d& xyz) const {
			double xx = xyz[0] - x;
			double yy = xyz[1] - y;
			double zz = xyz[2] - z;
			const double e = exponent * (xx * xx + yy * yy + zz * zz);
			return pow(xx, i) * pow(yy, j) * pow(zz, k) * exp(-e);
		}
	};

	/// make a set of Cartesian Gaussian functions located on atom

	/// @param[in]  atom    where the Cartesian Gaussian function is located
	/// @param[in]  l       the l-quantum number
	/// @param[in]  exponent    exponent of the Cartesian Gaussian function
	std::vector<CartesianGaussian> make_cartesian_guess(const Atom& atom, int l,
			const double exponent) {
		std::vector<CartesianGaussian> gg;

		// loop over all Cartesian components of l
		for (int kx = l; kx >= 0; kx--) {
			for (int ky = l - kx; ky >= 0; ky--) {
				int kz = l - kx - ky;
				gg.push_back(CartesianGaussian(atom, exponent, kx, ky, kz));
			}
		}

		return gg;
	}

	/// real solid harmonic Gaussian
	/// \note see https://en.wikipedia.org/wiki/Spherical_harmonics and https://en.wikipedia.org/wiki/Solid_harmonics
	class SolidHarmonicGaussian : public FunctionFunctorInterface<double, 3> {
	public:
		SolidHarmonicGaussian(const Atom& atom, const double e, int l, int m)
	: x(atom.x), y(atom.y), z(atom.z), exponent(e), Z(atom.atomic_number), l(l), m(m) {}

		double x, y, z;   ///< origin
		double exponent;  ///< exponent
		double Z;
		int l, m;         ///< (real) solid harmonic quanta

		template <typename T> int sgn(T val) const {
			return (T(0) < val) - (val < T(0));
		}

		double operator ()(const coord_3d& xyz) const;

		std::vector<coord_3d> special_points() const final {
			coord_3d coord;
			coord[0]=x;
			coord[1]=y;
			coord[2]=z;
			return std::vector<coord_3d>(1,coord);
		}

		Level special_level() const final {
			const double width = 1.0/sqrt(2.0*exponent); //  with width of the gaussian
			const int level = FunctionDefaults<3>::set_length_scale(width); // length scale for special points (center)
			return level;
		}
	};

	/// make a set of Cartesian Gaussian functions located on atom

	/// @param[in]  atom    where the Cartesian Gaussian function is located
	/// @param[in]  l       the l-quantum number
	/// @param[in]  exponent    exponent of the Cartesian Gaussian function
	std::vector<SolidHarmonicGaussian> make_solidharmonic_guess(const Atom& atom, int l,
			const double exponent)const {
		std::vector<SolidHarmonicGaussian> gg;

		// loop over all angular components of l
		for(int m=-l; m<=l; ++m) {
			gg.push_back(SolidHarmonicGaussian(atom, exponent, l, m));
		}

		return gg;
	}

};

}// namespace madness

#endif /* PAPER_CODE_BASISFUNCTIONS_H_ */
