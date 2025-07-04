/*
 * Nemocomplex.h
 *
 *  Created on: 14 Nov 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_ZNEMO_H_
#define SRC_APPS_CHEM_ZNEMO_H_

#include <madness/chem/CalculationParameters.h>
#include <madness/chem/MolecularOrbitals.h>
#include "madness/mra/QCCalculationParametersBase.h"
#include <madness/chem/SCFOperators.h>
#include <madness/chem/molecular_optimizer.h>
#include <madness/chem/molecularbasis.h>
#include <madness/chem/molecule.h>
#include <madness/chem/nemo.h>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>

namespace madness {

class Diamagnetic_potential_factor;



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


class Nemo_complex_Parameters : public QCCalculationParametersBase {
public:

	/// ctor reading out the input file
    Nemo_complex_Parameters() {
        initialize<double>("physical_B",0.0);
        initialize<double>("explicit_B",0.0);
        initialize<std::vector<double> >("box",{1.0, 0.01, 0.0, 0.0, 0.0});
        initialize<double>("shift",0.0);
        initialize<int>("printlevel",2);		// 0: energies, 1: fock matrix, 2: function sizes
        initialize<bool>("use_v_vector",false);
        initialize<double>("potential_radius",-1.0);
        initialize<std::vector<std::string> >("guess",std::vector<std::string>(),"list of l,m,exponent");	// atomic guess functions l, ml, exponent
        initialize<std::vector<std::string> >("guess_functions",std::vector<std::string>(),"list function names");	// atomic guess functions l, ml, exponent
    }

	Nemo_complex_Parameters(World& world, const commandlineparser& parser) : Nemo_complex_Parameters() {
		// read input file
        read_input_and_commandline_options(world,parser,"complex");
	}

    std::string get_tag() const override {
		return std::string("complex");
	}

	void set_derived_values() {
		double pb=physical_B();
		double eb=explicit_B();
		double remaining_B=fabs(pb-eb);
		double thresh=FunctionDefaults<3>::get_thresh()*0.1;

		double wave_function_radius=2.0*sqrt(-log(thresh)/remaining_B);
		double potential_radius=wave_function_radius*1.6;
		double box_radius=wave_function_radius*1.33;

		set_derived_value("box",std::vector<double>{box_radius,0.01});
		set_derived_value("potential_radius",potential_radius);
		set_derived_value("shift",physical_B());

	}

	int printlevel() const {return get<int>("printlevel");}
	double shift() const {return get<double>("shift");}
	double physical_B() const {return get<double>("physical_b");}
	double explicit_B() const {return get<double>("explicit_b");}
	std::vector<double> box() const {return get<std::vector<double> >("box");}
	double potential_radius() const {return get<double>("potential_radius");}
	bool use_v_vector() const {return get<bool>("use_v_vector");}
	std::vector<std::string> guess() const {return get<std::vector<std::string>>("guess");}
	std::vector<std::string> guess_functions() const {return get<std::vector<std::string>>("guess_functions");}

};


class Znemo : public NemoBase, public QCPropertyInterface {
	friend class Zcis;

	struct potentials {
		potentials(World& world, const std::size_t nmo) {
			vnuc_mo=zero_functions<double_complex,3>(world,nmo);
			diamagnetic_mo=zero_functions<double_complex,3>(world,nmo);
			lz_mo=zero_functions<double_complex,3>(world,nmo);
			J_mo=zero_functions<double_complex,3>(world,nmo);
			K_mo=zero_functions<double_complex,3>(world,nmo);
			spin_zeeman_mo=zero_functions<double_complex,3>(world,nmo);
			zeeman_R_comm=zero_functions<double_complex,3>(world,nmo);

		}

		void transform(const Tensor<double_complex>& U) {
			World& world=vnuc_mo[0].world();
			vnuc_mo = ::madness::transform(world, vnuc_mo, U);
			diamagnetic_mo = ::madness::transform(world, diamagnetic_mo, U);
			lz_mo = ::madness::transform(world, lz_mo, U);
			J_mo = ::madness::transform(world, J_mo, U);
			K_mo = ::madness::transform(world, K_mo, U);
			spin_zeeman_mo = ::madness::transform(world, spin_zeeman_mo, U);
			zeeman_R_comm = ::madness::transform(world, zeeman_R_comm, U);
		}


		std::vector<complex_function_3d> vnuc_mo;
		std::vector<complex_function_3d> diamagnetic_mo;
		std::vector<complex_function_3d> lz_mo;
		std::vector<complex_function_3d> J_mo;
		std::vector<complex_function_3d> K_mo;
		std::vector<complex_function_3d> spin_zeeman_mo;
		std::vector<complex_function_3d> lz_commutator;
		std::vector<complex_function_3d> zeeman_R_comm;
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
	Znemo(World& w, const commandlineparser& parser);

	/// compute the molecular energy
	double value() {return value(mol.get_all_coords());}

	/// compute the molecular energy
	double value(const Tensor<double>& x) override;

    void output_calc_info_schema(const double& energy) const;

    std::string name() const override {return "znemo";};

    static void help() {
        print_header2("help page for ZNEMO ");
        print("The znemo code computes Hartree-Fock energies and gradients in the presence of a strong");
        print("magnetic field");
        print("");
        print("You can print all available calculation parameters by running\n");
        print("znemo --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("znemo --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        Nemo_complex_Parameters zparam;
        Nemo::NemoCalculationParameters param;
        print("default calculations parameters for the znemo program are the same as for the nemo program");
        print("default parameters for the magnetic field of the znemo program are");
        zparam.print("complex","end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    bool selftest() override {
        return true;
    };

	/// adapt the thresholds consistently to a common value
    void recompute_factors_and_potentials(const double thresh);

    bool need_recompute_factors_and_potentials(const double thresh) const override;
    void invalidate_factors_and_potentials() override;

	void iterate();

	Tensor<double> gradient(const Tensor<double>& x) override;
	Tensor<double> gradient() {
		return gradient(mol.get_all_coords());
	}

	bool test() const;

	const CalculationParameters& get_calc_param() const {return calc_param;};
	Molecule& molecule() override {return mol;};
	const Molecule& molecule() const {return mol;};

	/// test the identity <F| f (T + Vdia ) f |F> = <F|f^2 (T + Udia) |F>
	bool test_U_potentials() const;

	/// analyse the results only

	/// @return	the energy
	double analyze();

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
		return ((not get_calc_param().spin_restricted()) and (get_calc_param().nbeta()>0));
	}

	void save_orbitals(std::string suffix) const;

	/// get initial orbitals from guesses
	void get_initial_orbitals();

	/// read the guess orbitals from a previous nemo or moldft calculation
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	read_real_guess() const;

	/// read the guess orbitals from a previous nemo or moldft calculation
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	read_complex_guess() const;

	/// read a list of functions for the guess
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	read_explicit_guess_functions() const;

	/// read guess as projection of the orbitals onto an AO basis set (for geometry optimization)
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	read_restartaodata() const;

	/// read orbitals from a znemo calculation
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	read_reference() const;

	/// hcore guess including the magnetic field
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	hcore_guess() const;

	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	custom_guess() const;

	void save_orbitals() const {
		std::string name="reference";
		print("saving orbitals to file",name);

		archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(), 1);

		ar & amo.size() & bmo.size() & aeps & beps;
		for (auto& a: amo) ar & a;
		for (auto& a: bmo) ar & a;

		// save AO restart data
		MolecularOrbitals<double_complex,3> amo1, bmo1;
		amo1.update_mos_and_eps(amo,aeps);
		bmo1.update_mos_and_eps(bmo,beps);

		MolecularOrbitals<double_complex,3>::save_restartaodata(world, mol, amo1, bmo1, aobasis);
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
	        for (size_t i=0; i<v.size(); i++) {
	            v[i] = vnew[i];
	            r[i] = rnew[i];
	        }
	    }
	    world.gop.fence();
	}

	double compute_energy(const std::vector<complex_function_3d>& amo, const potentials& apot,
			const std::vector<complex_function_3d>& bmo, const potentials& bpot, const bool do_print, const bool no_confinement) const;

	/// compute the potential operators applied on the orbitals
	potentials compute_potentials(const std::vector<complex_function_3d>& mo,
			const real_function_3d& density,
			const std::vector<complex_function_3d>& rhs) const;

	Tensor<double_complex> compute_vmat(
			const std::vector<complex_function_3d>& mo,
			const potentials& pot) const;

	std::vector<complex_function_3d> compute_residuals(
			const std::vector<complex_function_3d>& Vpsi,
			const std::vector<complex_function_3d>& psi,
			Tensor<double>& eps,
			const double levelshift=0.0) const;

	void canonicalize(std::vector<complex_function_3d>& amo,
			std::vector<complex_function_3d>& vnemo,
			potentials& pot,
			XNonlinearSolver<std::vector<complex_function_3d> ,double_complex, vector_function_allocator<double_complex,3> >& solver,
			Tensor<double_complex> fock, Tensor<double_complex> ovlp) const;

	std::vector<complex_function_3d> orthonormalize(const std::vector<complex_function_3d>& mo) const;

	std::vector<complex_function_3d> normalize(const std::vector<complex_function_3d>& mo) const;


	real_function_3d compute_nemo_density(const std::vector<complex_function_3d>& amo,
			const std::vector<complex_function_3d>& bmo) const {
		real_function_3d density=NemoBase::compute_density(amo);
		if (have_beta()) density+=NemoBase::compute_density(bmo);
		if (get_calc_param().spin_restricted()) density=density.scale(2.0);
		return density;
	}


	real_function_3d compute_nemo_spin_density(const std::vector<complex_function_3d>& amo,
			const std::vector<complex_function_3d>& bmo) const {
		if (get_calc_param().spin_restricted()) return real_function_3d(world);

		real_function_3d density=NemoBase::compute_density(amo);
		if (have_beta()) density-=NemoBase::compute_density(bmo);
		return density;
	}


	real_function_3d compute_density(const std::vector<complex_function_3d>& amo,
			const std::vector<complex_function_3d>& bmo) const;


	real_function_3d compute_spin_density(const std::vector<complex_function_3d>& amo,
			const std::vector<complex_function_3d>& bmo) const;

	std::vector<complex_function_3d> make_bra(const std::vector<complex_function_3d>& mo) const;


protected:

	Molecule mol;
	Nemo_complex_Parameters param;
    AtomicBasisSet aobasis;
    printleveler print_info=printleveler(2);

	/// standard calculation parameters
    Nemo::NemoCalculationParameters nemo_param;
	CalculationParameters calc_param;

	std::shared_ptr<Diamagnetic_potential_factor> diafac;

	/// the magnetic field B=rot(A)
	coord_3d B;

	/// nuclear potential
	std::shared_ptr<PotentialManager> potentialmanager;
public:
	/// the molecular orbitals -- alpha
	std::vector<complex_function_3d> amo, bmo;
protected:
	/// the orbital energies
	Tensor<double> aeps, beps;

	/// the spherical damping box
	real_function_3d sbox;

	/// the linear moments r={x,y,z}
	std::vector<real_function_3d> rvec;

	std::shared_ptr<real_convolution_3d> coulop;


	struct s_orbital : public FunctionFunctorInterface<double_complex,3> {

		double exponent=1.0;
		coord_3d origin={0.0,0.0,0.0};
		s_orbital(const double& expo, const coord_3d& o) :
			exponent(expo), origin(o) {
			MADNESS_ASSERT(exponent>0.0);
		}

		double_complex operator()(const coord_3d& xyz1) const {
			coord_3d xyz=xyz1-origin;
			double r=xyz.normf();
			return exp(-exponent*r);
		}

	};

	struct p_orbital : public FunctionFunctorInterface<double_complex,3> {

		long m=0;
		double exponent=1.0;
		coord_3d origin={0.0,0.0,0.0};
		p_orbital(const long& mm, const double& expo, const coord_3d& o)
			: m(mm), exponent(expo), origin(o) {
			MADNESS_ASSERT(m>-2 and m<2);
			MADNESS_ASSERT(exponent>0.0);
		}

		double_complex operator()(const coord_3d& xyz1) const {
			coord_3d xyz=xyz1-origin;
			double r=xyz.normf();

			if (m==0) {
				double theta=acos(xyz[2]/r);
				//double phi=atan2(xyz[1],xyz[0]);
				return r*exp(-exponent*r)*cos(theta);
			} else {
				double theta=acos(xyz[2]/r);
				double phi=atan2(xyz[1],xyz[0]);
				return r*exp(-exponent*r)*sin(theta)*exp(double(m)*double_complex(0.0,1.0)*phi);
			}
		}

	};

	/// following Doucot Pascier 2005
	struct landau_wave_function {

		int m=0;
		double l=0.0;		///< radius Eq. (29)
		coord_3d origin={0.0,0.0,0.0};

		landau_wave_function(const int m, const double B,
				const coord_3d& origin={0.0,0.0,0.0})
			: m(m), l(1.0/sqrt(B)), origin(origin) {
			print("origin",origin);
			MADNESS_ASSERT(m>=0);
			MADNESS_ASSERT(B>0);
		}

		/// following Eq. (43)
		double_complex operator()(const coord_3d& xyz1) const {
			const coord_3d xyz=xyz1-origin;
			const double_complex z(xyz[0],xyz[1]);			// z = x + iy
			const double_complex zbar(xyz[0],-xyz[1]);

			double z_confinement=exp(-0.01*xyz[2]);		// make wf decay in z-direction
			return std::pow(zbar/l,double(m)) * exp(-0.25*zbar*z/(l*l)) * z_confinement;
		}
	};

	public:
	void test_compute_current_density() const;

public:
	void test_landau_wave_function();
};




} // namespace madness
#endif /* SRC_APPS_CHEM_ZNEMO_H_ */
