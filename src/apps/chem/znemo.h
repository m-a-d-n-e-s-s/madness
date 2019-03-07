/*
 * Nemocomplex.h
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_ZNEMO_H_
#define SRC_APPS_CHEM_ZNEMO_H_


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
	const coord_3d offset;
	const coord_3d B_direction;
	spherical_box(const double r, const double h, const double t,
			const coord_3d o={0.0,0.0,0.0}, const coord_3d B_dir={0.0,0.0,1.0}) :
		radius(r), height(h), tightness(t), offset(o), B_direction(B_dir) {}

	double operator()(const coord_3d& xyz) const {
		// project out all contributions from xyz along the direction of the B field
		coord_3d tmp=(xyz-offset)*B_direction;
		const double inner=tmp[0]+tmp[1]+tmp[2];
		coord_3d proj=(xyz-offset)-B_direction*inner;
		double r=proj.normf();
		double v1=height/(1.0+exp(-tightness*height*(r-radius)));
		return 1.0-v1;

	}

    std::vector<coord_3d> special_points() const {
    	return std::vector<coord_3d>();
    }

};


/// functor for the diamagnetic term in a box

/// The diamagnetic term is a 2D harmonic potential, which makes the iterations diverge.
/// Approximate it by making it constant at a specific radius: construct a smooth approximation
/// for the piecewise function f(x)={{x, x<1}, {1,x>1}}, which is squared. For this approximation
/// see https://doi.org/10.1186/s40064-016-3278-y
struct diamagnetic_boxed_functor : public FunctionFunctorInterface<double,3> {
	double radius;	//
	const double tightness;	// alpha in the article
	diamagnetic_boxed_functor(const double r, const double t=30.0) :
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
	enum parameterenum {B_, box_, shift_,printlevel_,diamagnetic_height_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<std::vector<double> >(B_,{"B",{0.0}}),
        		init<std::vector<double> >(box_,{"box",{15.0,1.0,4.0,0.0,0.0,0.0}}),
				init<double>(shift_,{"shift",0.0}),
				init<int>(printlevel_,{"printlevel",1}),		// 0: energies, 1: fock matrix, 2: function sizes
				init<double>(diamagnetic_height_,{"diamagnetic_height",30})
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
	double shift() const {return get<double>(shift_);}
	std::vector<double> B() const {return get<std::vector<double> >(B_);}
	std::vector<double> box() const {return get<std::vector<double> >(box_);}
	double diamagnetic_height() const {return get<double>(diamagnetic_height_);}


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



class Znemo {
	friend class Zcis;

	struct potentials {
		potentials(World& world, const std::size_t nmo) {
			vnuc_mo=zero_functions<double_complex,3>(world,nmo);
			diamagnetic_mo=zero_functions<double_complex,3>(world,nmo);
			lz_mo=zero_functions<double_complex,3>(world,nmo);
			J_mo=zero_functions<double_complex,3>(world,nmo);
			K_mo=zero_functions<double_complex,3>(world,nmo);
			spin_zeeman_mo=zero_functions<double_complex,3>(world,nmo);
		}
		std::vector<complex_function_3d> vnuc_mo;
		std::vector<complex_function_3d> diamagnetic_mo;
		std::vector<complex_function_3d> lz_mo;
		std::vector<complex_function_3d> J_mo;
		std::vector<complex_function_3d> K_mo;
		std::vector<complex_function_3d> spin_zeeman_mo;
	};

public:
	Znemo(World& w) : world(w), param(world), molecule("input"), cparam() {
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
		coord_3d box_offset{param.box()[3],param.box()[4],param.box()[5]};
		spherical_box sbox2(param.box()[0],param.box()[1],param.box()[2],box_offset);
		sbox=real_factory_3d(world).functor(sbox2);
		save(sbox,"sbox");
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

	void read_orbitals() {
		std::string name="reference";
		print("reading orbitals from file",name);

		archive::ParallelInputArchive ar(world, name.c_str(), 1);
		std::size_t namo, nbmo;

		ar & namo & nbmo & aeps & beps;
		amo.resize(namo);
		bmo.resize(nbmo);
		for (auto& a: amo) ar & a;
		for (auto& a: bmo) ar & a;
	}

	void save_orbitals() const {
		std::string name="reference";
		print("saving orbitals to file",name);

		archive::ParallelOutputArchive ar(world, name.c_str(), 1);

		ar & amo.size() & bmo.size() & aeps & beps;
		for (auto& a: amo) ar & a;
		for (auto& a: bmo) ar & a;
	}

	void do_step_restriction(const std::vector<complex_function_3d>& mo,
			std::vector<complex_function_3d>& mo_new) const;

	/// rotate the KAIN subspace (cf. SCF.cc)
	template<typename solverT>
	void rotate_subspace(const Tensor<double_complex>& U, solverT& solver) const {
		double trantol=FunctionDefaults<3>::get_thresh()*0.1;
	    std::vector < std::vector<Function<double_complex, 3> > >&ulist = solver.get_ulist();
	    std::vector < std::vector<Function<double_complex, 3> > >&rlist = solver.get_rlist();
	    for (unsigned int iter = 0; iter < ulist.size(); ++iter) {
	    	std::vector<Function<double_complex, 3> >& v = ulist[iter];
	    	std::vector<Function<double_complex, 3> >& r = rlist[iter];
	    	std::vector<Function<double_complex, 3> > vnew = transform(world, v, U, trantol, false);
	    	std::vector<Function<double_complex, 3> > rnew = transform(world, r, U, trantol, true);

	        world.gop.fence();
	        for (int i=0; i<v.size(); i++) {
	            v[i] = vnew[i];
	            r[i] = rnew[i];
	        }
	    }
	    world.gop.fence();
	}
	double compute_energy(const std::vector<complex_function_3d>& amo, const potentials& apot,
			const std::vector<complex_function_3d>& bmo, const potentials& bpot, const bool do_print) const;

	/// compute the action of the Lz =i r x del operator on rhs
	std::vector<complex_function_3d> Lz(const std::vector<complex_function_3d>& rhs) const;

	/// compute the diamagnetic local potential (B is in z direction -> dia = x^2 + y^2
	std::vector<complex_function_3d> diamagnetic(const std::vector<complex_function_3d>& rhs) const {
//		auto dia = [](const coord_3d& r) {return r[0]*r[0] + r[1]*r[1];};
		if (B!=0.0) {
			const double height=param.diamagnetic_height();
			const double radius=sqrt(8.0*height/(B*B));
			return rhs*diamagnetic_boxed + radius*radius*rhs;
		} else {
			return zero_functions<double_complex,3>(world,rhs.size());
		}
	}

	real_function_3d make_diamagnetic_boxed() const {
		if (B!=0.0) {
			const double height=param.diamagnetic_height();
			const double radius=sqrt(8.0*height/(B*B));
			if (world.rank()==0) print("computed radius for diamagnetic term from its cutoff height to",radius);
			real_function_3d diabox=real_factory_3d(world).functor(diamagnetic_boxed_functor(radius));
			save(diabox*(0.125*(B*B)),"diabox_B"+stringify(B));
			return diabox;
		} else {
			auto one = [](const coord_3d& r) {return 1.0;};
			real_function_3d nobox=real_factory_3d(world).functor(one);
			return nobox;
		}
	}

	/// compute the potential operators applied on the orbitals
	potentials compute_potentials(const std::vector<complex_function_3d>& mo,
			const real_function_3d& density,
			std::vector<complex_function_3d>& rhs) const;

	Tensor<double_complex> compute_vmat(
			const std::vector<complex_function_3d>& mo,
			const potentials& pot) const;

	std::vector<complex_function_3d> compute_residuals(
			const std::vector<complex_function_3d>& Vpsi,
			const std::vector<complex_function_3d>& psi,
			Tensor<double>& eps) const;

	void canonicalize(std::vector<complex_function_3d>& amo,
			std::vector<complex_function_3d>& vnemo,
			XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator>& solver,
			Tensor<double_complex> fock, Tensor<double_complex> ovlp) const;

	void orthonormalize(std::vector<complex_function_3d>& amo) const;

	void normalize(std::vector<complex_function_3d>& mo) const {
		std::vector<double> n=norm2s(world,mo);
		for (int i=0; i<mo.size(); ++i) mo[i].scale(1.0/(n[i]));
	}

	static Tensor<double_complex> Q2(const Tensor<double_complex> & s) {
		Tensor<double_complex> Q = -0.5*s;
		for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
		return Q;
	}

protected:

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

	/// the boxed diamagnetic potential (for a given B)
	real_function_3d diamagnetic_boxed;

	std::shared_ptr<real_convolution_3d> coulop;

};




} // namespace madness
#endif /* SRC_APPS_CHEM_ZNEMO_H_ */
