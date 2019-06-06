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

class Diamagnetic_potential_factor;


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


struct printleveler {
	bool print_nothing_=true;
	bool print_setup_=false;
	bool print_energies_=false;
	bool print_timings_=false;
	bool print_convergence_=false;
	bool print_debug_information_=false;

	printleveler(int printlevel) {
		print_debug_information_=(printlevel>=10);
		print_convergence_=(printlevel>=4);
		print_timings_=(printlevel>=3);
		print_energies_=(printlevel>=2);
		print_setup_=(printlevel>=1);
	}

	bool print_setup() const {return print_setup_;}
	bool print_energies() const {return print_energies_;}
	bool print_timings() const {return print_timings_;}
	bool print_convergence() const {return print_convergence_;}
	bool print_debug() const {return print_debug_information_;}

};


class Nemo_complex_Parameters : public CalculationParametersBase {
public:
	enum parameterenum {physical_B_, explicit_B_, box_, box_softness_, shift_, printlevel_,
		 potential_radius_, use_v_vector_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<double>(physical_B_,{"physical_B",0.0}),
        		init<double>(explicit_B_,{"explicit_B",0.0}),
        		init<std::vector<double> >(box_,{"box",{1.0, 0.01, 0.0, 0.0, 0.0}}),
				init<double>(shift_,{"shift",0.0}),
				init<int>(printlevel_,{"printlevel",2}),		// 0: energies, 1: fock matrix, 2: function sizes
				init<bool>(use_v_vector_,{"use_v_vector",false}),
				init<double>(potential_radius_,{"potential_radius",-1.0}),
    };

	/// ctor reading out the input file
	Nemo_complex_Parameters(World& world) {

		// read input file
		read(world,"input","complex",params);

	}

	void set_derived_values() {
		double pb=physical_B();
		double eb=explicit_B();
		double remaining_B=fabs(pb-eb);
		double thresh=FunctionDefaults<3>::get_thresh()*0.1;

		double wave_function_radius=2.0*sqrt(-log(thresh)/remaining_B);
		double potential_radius=wave_function_radius*1.4;
		double box_radius=wave_function_radius*1.25;

		// set the diamagnetic height unless explicitly given
		params[box_].set_derived_value(std::vector<double>({box_radius,0.01}));
		params[potential_radius_].set_derived_value(potential_radius);


		// set derived values
//		params[param2_].set_derived_value(this->get<int>(param1_)*10.0);

		// print final parameters
//		if (world.rank()==0) print(params,"Our parameters");
	}

	int printlevel() const {return get<int>(printlevel_);}
	double shift() const {return get<double>(shift_);}
	double physical_B() const {return get<double>(physical_B_);}
	double explicit_B() const {return get<double>(explicit_B_);}
	std::vector<double> box() const {return get<std::vector<double> >(box_);}
	double potential_radius() const {return get<double>(potential_radius_);}
	bool use_v_vector() const {return get<bool>(use_v_vector_);}


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
			for (int i=0; i<3; ++i) GpVmo.push_back(zero_functions<double_complex,3>(world,nmo));
			Gpscalar=zero_functions<double_complex,3>(world,nmo);

		}

		void transform(const Tensor<double_complex>& U) {
			World& world=vnuc_mo[0].world();
			vnuc_mo = ::madness::transform(world, vnuc_mo, U);
			diamagnetic_mo = ::madness::transform(world, diamagnetic_mo, U);
			lz_mo = ::madness::transform(world, lz_mo, U);
			J_mo = ::madness::transform(world, J_mo, U);
			K_mo = ::madness::transform(world, K_mo, U);
			spin_zeeman_mo = ::madness::transform(world, spin_zeeman_mo, U);
			for (auto& a : GpVmo) a=::madness::transform(world,a,U);
			Gpscalar = ::madness::transform(world, Gpscalar, U);
		}


		std::vector<complex_function_3d> vnuc_mo;
		std::vector<complex_function_3d> diamagnetic_mo;
		std::vector<complex_function_3d> lz_mo;
		std::vector<complex_function_3d> J_mo;
		std::vector<complex_function_3d> K_mo;
		std::vector<complex_function_3d> spin_zeeman_mo;
		std::vector<complex_function_3d> lz_commutator;
		std::vector<std::vector<complex_function_3d> > GpVmo;	// potentials for the derivative of the BSH operator
		std::vector<complex_function_3d> Gpscalar;				// scalar terms arising from the Gp treatment
	};
public:
	struct timer {
        World& world;
	    double ttt,sss;
	    bool do_print=true;

	    timer(World& world, bool do_print=true) : world(world), do_print(do_print) {
	        world.gop.fence();
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void tag(const std::string msg) {
            world.gop.fence();
	        double tt1=wall_time()-ttt;
	        double ss1=cpu_time()-sss;
	        if (world.rank()==0 and do_print) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void end(const std::string msg) {
            world.gop.fence();
            double tt1=wall_time()-ttt;
            double ss1=cpu_time()-sss;
            if (world.rank()==0 and do_print) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
        }
	};


public:
	Znemo(World& w);

	/// compute the molecular energy
	double value();

	void iterate();

	bool test() const;

	/// test the identity <F| f (T + Vdia ) f |F> = <F|f^2 (T + Udia) |F>
	bool test_U_potentials() const;

	// analyse the results only
	void analyze() const;

	/// compute the expectation value of the kinetic momentum p
	Tensor<double> compute_kinetic_momentum() const {
	    Tensor<double> p_exp(3);
	    for (auto& mo : amo) p_exp+=imag(inner(world,mo,grad(mo)));
	    for (auto& mo : bmo) p_exp+=imag(inner(world,mo,grad(mo)));
	    return p_exp;
	}

	/// compute the expectation value of the linear moment r
	Tensor<double> compute_linear_moment() const {
		std::vector<real_function_3d> r(3);
	    r[0]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
	    r[1]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});
	    r[2]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[2];});

	    Tensor<double> r_exp(3);
	    for (auto& mo : amo) r_exp+=real(inner(world,mo,r*mo));
	    for (auto& mo : bmo) r_exp+=real(inner(world,mo,r*mo));
	    return r_exp;
	}

	/// compute the expectation value of the magnetic vector potential A
	Tensor<double> compute_magnetic_potential_expectation(const std::vector<real_function_3d>& A) const {
	    Tensor<double> A_exp(3);
	    for (auto& mo : amo) A_exp+=real(inner(world,mo,A*mo));
	    for (auto& mo : bmo) A_exp+=real(inner(world,mo,A*mo));
	    return A_exp;
	}

	/// compute the shift of the molecule such that the kinetic momentum vanishes
	Tensor<double> compute_standard_gauge_shift(const Tensor<double>& p_exp) const {
		Tensor<double> S(3);
		const double B=param.physical_B();

	    S(0l)=-p_exp(1);
	    S(1) =p_exp(0l);
	    S(2) =p_exp(2);
	    S*=-2.0/(B*(amo.size()+bmo.size()));
	    return S;
	}

	/// compute the current density
	std::vector<real_function_3d> compute_current_density(
			const std::vector<complex_function_3d>& alpha_mo,
			const std::vector<complex_function_3d>& beta_mo) const;

	/// solve the SCF iterations
	void solve_SCF();

	/// compute the magnetic vector potential A
	static std::vector<real_function_3d> compute_magnetic_vector_potential(World& world,
			const coord_3d& Bvec) {
		std::vector<real_function_3d> r(3), A(3);
	    r[0]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
	    r[1]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});
	    r[2]=real_factory_3d(world).functor([] (const coord_3d& r) {return r[2];});

	    A[0]=Bvec[1]*r[2]-Bvec[2]*r[1];
	    A[1]=Bvec[2]*r[0]-Bvec[0l]*r[2];
	    A[2]=Bvec[0l]*r[1]-Bvec[1]*r[0];

		return 0.5*A;
	}

	/// are there explicit beta orbitals
	bool have_beta() const {
		return ((not cparam.spin_restricted) and (cparam.nbeta>0));
	}

	void save_orbitals(std::string suffix) const;

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

	double compute_energy_no_confinement(const std::vector<complex_function_3d>& amo, const potentials& apot,
			const std::vector<complex_function_3d>& bmo, const potentials& bpot, const bool do_print) const;

	/// compute the potential operators applied on the orbitals
	potentials compute_potentials(const std::vector<complex_function_3d>& mo,
			const real_function_3d& density,
			std::vector<complex_function_3d>& rhs) const;

	Tensor<double_complex> compute_vmat(
			const std::vector<complex_function_3d>& mo,
			const potentials& pot) const;

	std::vector<complex_function_3d> compute_residuals(
			const std::vector<complex_function_3d>& Vpsi,
			const std::vector<std::vector<complex_function_3d> > GpVpsi,
			const std::vector<complex_function_3d>& psi,
			Tensor<double>& eps) const;

	void canonicalize(std::vector<complex_function_3d>& amo,
			std::vector<complex_function_3d>& vnemo,
			potentials& pot,
			XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, allocator>& solver,
			Tensor<double_complex> fock, Tensor<double_complex> ovlp) const;

	void orthonormalize(std::vector<complex_function_3d>& amo) const;

	void normalize(std::vector<complex_function_3d>& mo) const;

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
    printleveler print_info=printleveler(2);

	/// standard calculation parameters
	CalculationParameters cparam;

	std::shared_ptr<Diamagnetic_potential_factor> diafac;

	/// the magnetic field B=rot(A)
	coord_3d B;

	/// nuclear potential
	real_function_3d vnuclear;

	/// the molecular orbitals -- alpha
	std::vector<complex_function_3d> amo, bmo;

	/// the orbital energies
	Tensor<double> aeps, beps;

	/// the spherical damping box
	real_function_3d sbox;

	std::shared_ptr<real_convolution_3d> coulop;

	static double_complex p_plus(const coord_3d& xyz1) {
		coord_3d xyz=xyz1-coord_3d({5.0,0,0});
		double r=xyz.normf();
		double theta=acos(xyz[2]/r);
		double phi=atan2(xyz[1],xyz[0]);
		return r*exp(-r/2.0)*sin(theta)*exp(double_complex(0.0,1.0)*phi);
//		return r*exp(-r/2.0)*exp(double_complex(0.0,1.0)*phi);
	}

	void test_compute_current_density() const;

};




} // namespace madness
#endif /* SRC_APPS_CHEM_ZNEMO_H_ */