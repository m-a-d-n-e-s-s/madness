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


class Nemo_complex_Parameters : public CalculationParametersBase {
public:
	enum parameterenum {physical_B_, explicit_B_, box_, box_softness_, shift_, printlevel_, diamagnetic_height_,
		use_greensp_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<double>(physical_B_,{"physical_B",0.0}),
        		init<double>(explicit_B_,{"explicit_B",0.0}),
        		init<std::vector<double> >(box_,{"box",{-1.0, 1.0, 4.0, 0.0, 0.0, 0.0}}),
        		init<double>(box_softness_,{"box_softness",1.0}),
				init<double>(shift_,{"shift",0.0}),
				init<int>(printlevel_,{"printlevel",1}),		// 0: energies, 1: fock matrix, 2: function sizes
				init<double>(diamagnetic_height_,{"diamagnetic_height",-1.0}),
				init<bool>(use_greensp_,{"greensp",false}),
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
		double differential_Bsquare=sqrt(pb*pb-eb*eb);
		double thresh=FunctionDefaults<3>::get_thresh()*0.1;
		::madness::print("thresh, rB",thresh,remaining_B);

		double wave_function_radius=2.0*sqrt(-log(thresh)/remaining_B);
		double potential_radius=wave_function_radius*1.4;
		double box_radius=wave_function_radius*1.1;

		::madness::print("wave function radius  ",wave_function_radius);
		::madness::print("potential_radius      ",potential_radius);
		::madness::print("differential_B_square ",differential_Bsquare);
		::madness::print("diamagnetic_height    ",params[diamagnetic_height_]);

		// set the diamagnetic height unless explicitly given
		params[box_].set_derived_value(std::vector<double>({box_radius,1.0,box_softness()}));
		params[diamagnetic_height_].set_derived_value(0.125*
				differential_Bsquare*differential_Bsquare*potential_radius*potential_radius);


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
	double box_softness() const {return get<double>(box_softness_);}
	double diamagnetic_height() const {return get<double>(diamagnetic_height_);}
	bool use_greensp() const {return get<bool>(use_greensp_);}


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
		std::vector<std::vector<complex_function_3d> > GpVmo;	// potentials for the derivative of the BSH operator
		std::vector<complex_function_3d> Gpscalar;				// scalar terms arising from the Gp treatment
	};

	struct timer {
        World& world;
	    double ttt,sss;
	    timer(World& world) : world(world) {
	        world.gop.fence();
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void tag(const std::string msg) {
            world.gop.fence();
	        double tt1=wall_time()-ttt;
	        double ss1=cpu_time()-sss;
	        if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void end(const std::string msg) {
            world.gop.fence();
            double tt1=wall_time()-ttt;
            double ss1=cpu_time()-sss;
            if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
        }

	};

public:
	Znemo(World& w);

	/// compute the molecular energy
	double value();

	void test();
	void test2();
	void test_gp(const std::vector<complex_function_3d>& arg,
			const std::vector<std::vector<complex_function_3d> > varg) const;
	void test_gp2(const std::vector<complex_function_3d>& arg,
			const std::vector<std::vector<complex_function_3d> > varg) const;

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

	/// turn the magnetic strength from the input file to a coord_3d
	coord_3d Bvec(const double B) {return coord_3d{0.0,0.0,B};}

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

	static std::vector<coord_3d> compute_v_vector(World& world, const coord_3d& B, const Molecule& mol) {
		std::vector<coord_3d> v;
		for (auto& c : mol.get_all_coords_vec()) v.push_back(0.5*cross(B,c));
		return v;
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

	/// compute the action of the Lz =i r x del operator on rhs
	std::vector<complex_function_3d> Lz(const std::vector<complex_function_3d>& rhs) const;

	/// compute the Lz operator, prepared with integration by parts for the derivative of the BSH operator
	std::vector<std::vector<complex_function_3d> > Lz_Gp(const std::vector<complex_function_3d>& rhs) const;

	/// compute the action of the Lz =i r x del operator on rhs
	complex_function_3d Lz(const complex_function_3d& rhs) const {
		std::vector<complex_function_3d> vrhs(1,rhs);
		return Lz(vrhs)[0];
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

	/// standard calculation parameters
	CalculationParameters cparam;

	std::shared_ptr<Diamagnetic_potential_factor> diafac;

	/// the magnetic potential A = 1/2 B cross r
	std::vector<real_function_3d> A;

	/// the magnetic field B=rot(A)
	coord_3d B;

	/// the position of the nuclei in the "A" space: v = 1/2 B cross R
	std::vector<coord_3d> v;

	/// nuclear potential
	real_function_3d vnuclear;

	/// the molecular orbitals -- alpha
	std::vector<complex_function_3d> amo, bmo;

	/// the orbital energies
	Tensor<double> aeps, beps;

	/// the spherical damping box
	real_function_3d sbox;

	std::shared_ptr<real_convolution_3d> coulop;

	static double_complex p_plus(const coord_3d& xyz) {
		double r=xyz.normf();
		double theta=acos(xyz[2]/r);
		double phi=atan2(xyz[1],xyz[0]);
		return r*exp(-r/2.0)*sin(theta)*exp(double_complex(0.0,1.0)*phi);
	}

	void test_compute_current_density() const;

};




} // namespace madness
#endif /* SRC_APPS_CHEM_ZNEMO_H_ */
