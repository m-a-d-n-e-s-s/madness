/*
 * Nemocomplex.h
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_NEMOCOMPLEX_H_
#define SRC_APPS_CHEM_NEMOCOMPLEX_H_


#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <chem/molecule.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>
#include <chem/SCFOperators.h>
#include <chem/CalculationParametersBaseMap.h>
#include <chem/CalculationParameters.h>
#include <chem/molecularbasis.h>
#include <madness/mra/operator.h>


namespace madness {


struct spherical_box : public FunctionFunctorInterface<double,3> {
	const double radius;
	const double height;
	const double tightness;
	spherical_box(const double r, const double h, const double t) :
		radius(r), height(h), tightness(t) {}

	double operator()(const coord_3d& xyz) const {
		double r=xyz.normf();
		double v1=height/(1.0+exp(-tightness*height*(r-radius)));
		return 1.0-v1;
	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};


// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct allocator {
	World& world;
	const int n;

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(World& world, const int nn) :
			world(world), n(nn) {
	}

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(const allocator& other) :
			world(other.world), n(other.n) {
	}

	/// allocate a vector of n empty functions
	std::vector<complex_function_3d> operator()() {
		return zero_functions_compressed<double_complex, 3>(world, n);
	}
};


class Nemo_complex_Parameters : public CalculationParametersBase {
public:
	enum parameterenum {B_, box_, shift_,printlevel_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<std::vector<double> >(B_,{"B",{0.0}}),
        		init<std::vector<double> >(box_,{"box",{15.0,1.0,4.0}}),
				init<double>(shift_,{"shift",0.0}),
				init<int>(printlevel_,{"printlevel",1})		// 0: energies, 1: fock matrix, 2: function sizes
    };

	/// ctor reading out the input file
	Nemo_complex_Parameters(World& world) {

		// read input file
		read(world,"input","complex",params);

		// set derived values
//		params[param2_].set_derived_value(this->get<int>(param1_)*10.0);

		// print final parameters
//		if (world.rank()==0) print(params,"Our parameters");
	}

	int printlevel() const {return get<int>(printlevel_);}
	std::vector<double> B() const {return get<std::vector<double> >(B_);}
	std::vector<double> box() const {return get<std::vector<double> >(box_);}

	/// return the value of the parameter
	template<typename T>
	T get(parameterenum k) const {
		if (params.find(int(k))!=params.end()) {
			return params.find(int(k))->second.get_parameter<T>().get();
		} else {
			MADNESS_EXCEPTION("could not fine parameter ",1);
		}
	}
};



class Nemo_complex {
public:
	Nemo_complex(World& w) : world(w), param(world), molecule("input"), cparam() {
		cparam.read_file("input");

	    FunctionDefaults<3>::set_k(cparam.k);
	    FunctionDefaults<3>::set_thresh(cparam.econv);
	    FunctionDefaults<3>::set_refine(true);
	    FunctionDefaults<3>::set_initial_level(5);
	    FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_cubic_cell(-cparam.L, cparam.L);

        aobasis.read_file(cparam.aobasis);
        cparam.set_molecular_info(molecule, aobasis, 0);

		if (world.rank()==0) {
			param.print(param.params,"complex","end");
			cparam.print(world);
			molecule.print();
		}

		// the guess is read from a previous nemo calculation
		// make sure the molecule was not reoriented there
		if (not cparam.no_orient) {
			MADNESS_EXCEPTION("the molecule of the reference calculation was reoriented\n\n",1);
		}

		coulop=std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world,cparam.lo,cparam.econv));
		spherical_box sbox2(param.box()[0],param.box()[1],param.box()[2]);
		sbox=real_factory_3d(world).functor(sbox2);
	};

	/// compute the molecular energy
	double value();

	/// solve the SCF iterations
	void solve_SCF();

	/// are there explicit beta orbitals
	bool have_beta() const {
		return ((not cparam.spin_restricted) and (cparam.nbeta>0));
	}

	void save_orbitals(int iter) const {
		for (int i=0; i<amo.size(); ++i) save(amo[i],"amo"+stringify(i)+"_iter"+stringify(iter));
		for (int i=0; i<bmo.size(); ++i) save(bmo[i],"bmo"+stringify(i)+"_iter"+stringify(iter));
		for (int i=0; i<amo.size(); ++i) save(madness::abssq(amo[i]),"absamo"+stringify(i)+"_iter"+stringify(iter));
		for (int i=0; i<bmo.size(); ++i) save(madness::abssq(bmo[i]),"absbmo"+stringify(i)+"_iter"+stringify(iter));
	}

	/// read the guess orbitals from a previous nemo or moldft calculation
	std::vector<complex_function_3d> read_guess(const std::string& spin) const;

	/// compute the action of the Lz =i r x del operator on rhs
	std::vector<complex_function_3d> Lz(const std::vector<complex_function_3d>& rhs) const;

	/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
	real_function_3d diamagnetic() const {
		auto dia = [](const coord_3d& r) {return r[0]*r[0] + r[1]*r[1];};
		real_function_3d result=real_factory_3d(world).functor(dia);
		return result;
	}

	/// compute the potential operators applied on the orbitals
	void compute_potentials(const std::vector<complex_function_3d>& mo,
			real_function_3d& density,
			std::vector<complex_function_3d>& Vnemo,
			std::vector<complex_function_3d>& lznemo,
			std::vector<complex_function_3d>& dianemo,
			std::vector<complex_function_3d>& spin_zeeman_nemo,
			std::vector<complex_function_3d>& Knemo,
			std::vector<complex_function_3d>& Jnemo) const;

	Tensor<double_complex> compute_vmat(
			const std::vector<complex_function_3d>& mo,
			const std::vector<complex_function_3d>& Vnemo,
			const std::vector<complex_function_3d>& lznemo,
			const std::vector<complex_function_3d>& dianemo,
			const std::vector<complex_function_3d>& spin_zeeman_nemo,
			const std::vector<complex_function_3d>& Knemo,
			const std::vector<complex_function_3d>& Jnemo) const;

	std::vector<complex_function_3d> compute_residuals(
			const std::vector<complex_function_3d>& Vpsi,
			const std::vector<complex_function_3d>& psi,
			Tensor<double>& eps) const;

	void canonicalize(std::vector<complex_function_3d>& amo,
			std::vector<complex_function_3d>& vnemo,
			Tensor<double_complex> fock, Tensor<double_complex> ovlp) const;

	void orthonormalize(std::vector<complex_function_3d>& amo) const;

	void normalize(std::vector<complex_function_3d>& mo) const {
		std::vector<double> n=norm2s(world,mo);
		for (int i=0; i<mo.size(); ++i) mo[i].scale(1.0/(n[i]));
	}

	Tensor<double_complex> Q2(const Tensor<double_complex> & s) const {
		Tensor<double_complex> Q = -0.5*s;
		for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
		return Q;
	}

private:

	World& world;
	Molecule molecule;
	Nemo_complex_Parameters param;
    AtomicBasisSet aobasis;

	/// standard calculation parameters
	CalculationParameters cparam;

	/// nuclear potential
	real_function_3d vnuclear;

	/// the molecular orbitals -- alpha
	std::vector<complex_function_3d> amo, bmo;

	/// the orbital energies
	Tensor<double> aeps, beps;

	/// the external magnetic field in z-direction
	double B=0.0;

	/// the spherical damping box
	real_function_3d sbox;

	std::shared_ptr<real_convolution_3d> coulop;

};




} // namespace madness
#endif /* SRC_APPS_CHEM_NEMOCOMPLEX_H_ */
