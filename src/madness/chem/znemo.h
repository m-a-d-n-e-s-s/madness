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
#include <madness/mra/magnetic_kernel.h>

namespace madness {




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

        // magnetic Green's-function kernels (madness/mra/magnetic_kernels.md); the kernel parameters are in the
        // magnetic_kernel block (MagneticKernelParameters)
        initialize<std::string>("kernel","chain_symmetric","magnetic Green's function: ft, chain_symmetric, "
                "chain_landau, mehler",{"ft","chain_symmetric","chain_landau","mehler"});
        initialize<std::vector<int> >("guess_lz",std::vector<int>(),"hcore guess: L_z of the occupied orbitals (alpha in "
                "order, beta the first nbeta; B along a linear molecule only); empty: aufbau");
    }

	Nemo_complex_Parameters(World& world, const commandlineparser& parser) : Nemo_complex_Parameters() {
		// read input file
        read_input_and_commandline_options(world,parser,"complex");
	}

    std::string get_tag() const override {
		return std::string("complex");
	}

	double physical_B() const {return get<double>("physical_b");}

	std::string kernel() const {return get<std::string>("kernel");}
	std::vector<int> guess_lz() const {return get<std::vector<int>>("guess_lz");}
};


class Znemo : public NemoBase, public QCPropertyInterface {
	friend class Zcis;

	struct potentials {
		potentials(World& world, const std::size_t nmo) {
			vnuc_mo=zero_functions<double_complex,3>(world,nmo);
			J_mo=zero_functions<double_complex,3>(world,nmo);
			K_mo=zero_functions<double_complex,3>(world,nmo);
			zeeman_R_comm=zero_functions<double_complex,3>(world,nmo);
			vbare_mo=zero_functions<double_complex,3>(world,nmo);

		}

		void transform(const Tensor<double_complex>& U) {
			World& world=vnuc_mo[0].world();
			vnuc_mo = ::madness::transform(world, vnuc_mo, U);
			J_mo = ::madness::transform(world, J_mo, U);
			K_mo = ::madness::transform(world, K_mo, U);
			zeeman_R_comm = ::madness::transform(world, zeeman_R_comm, U);
			vbare_mo = ::madness::transform(world, vbare_mo, U);
		}


		std::vector<complex_function_3d> vnuc_mo;
		std::vector<complex_function_3d> J_mo;
		std::vector<complex_function_3d> K_mo;
		std::vector<complex_function_3d> zeeman_R_comm;
		/// magnetic kernels with a correlation factor: the bare (smoothed) nuclear potential on the nemos, for the
		/// Fock matrix and the energy with the symmetric kinetic term (kinetic_pi_matrix)
		std::vector<complex_function_3d> vbare_mo;
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

	const CalculationParameters& get_calc_param() const {return calc_param;};
	Molecule& molecule() override {return mol;};
	const Molecule& molecule() const {return mol;};

	/// test the identity <F| f (T + Vdia ) f |F> = <F|f^2 (T + Udia) |F>

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

	/// initial orbitals: the restart sequence of moldft and nemo, with the generalized hcore as the initial guess
	/// @return false if the orbitals read are the answer and the SCF equations are not solved
	bool get_initial_orbitals();

	/// hcore guess: the one-electron Hamiltonian of the solver in magnetically translated, Landau-squeezed AOs
	std::pair<MolecularOrbitals<double_complex,3>, MolecularOrbitals<double_complex,3> >
	hcore_guess() const;

	/// write the orbitals to <prefix>.restartdata (RestartMetadata header, representation znemo)
	void save_orbitals() const;

	/// read the orbitals from <prefix>.restartdata
	void load_restartdata();

	/// the nuclear correlation factor for the restart header, as nemo writes it
	std::string restart_ncf() const {
		return nemo_param.ncf().first+":"+std::to_string(nemo_param.ncf().second);
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
			const std::vector<complex_function_3d>& bmo, const potentials& bpot, const bool do_print) const;

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
			Tensor<double>& eps) const;

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

	/// right-hand side V psi of the iteration: the potentials not contained in the Green's function
	std::vector<complex_function_3d> assemble_rhs(const potentials& pot, const std::vector<complex_function_3d>& mo) const;

	/// magnetic kernels: <pi phi_i|pi phi_j>/2 with pi = -i nabla + A in the gauge of the kernel
	Tensor<double_complex> kinetic_pi_matrix(const std::vector<complex_function_3d>& mo) const;

	/// magnetic kernels with a nuclear correlation factor: symmetric kinetic term, bare V in the Fock matrix
	bool kernel_nemo() const {return ncf and ncf->type()!=NuclearCorrelationFactor::None;}

	/// pi_d f = -i d_d f + A_d f in the gauge of the kernel
	std::vector<complex_function_3d> pi_apply(const std::vector<complex_function_3d>& f, const int axis) const;

	/// nemo: tmin from the right-hand sides g_i (see the definition)
	void adapt_kernel_tmin(const std::vector<complex_function_3d>& g, const Tensor<double>& E, const bool first);

	/// smoothing radii of the nuclear potentials: from eprec (error per atom), shrunk for the field-compressed density
	void update_smoothing_radii();

	/// derived parameters of the magnetic kernel from the molecule, the complex and the nemo block
	void set_kernel_derived_values(MagneticKernelParameters& kp, const double thresh) const;

	/// the one-electron potentials (nuclear, paramagnetic, diamagnetic) applied on rhs; J and K are zero
	potentials one_electron_potentials(const std::vector<complex_function_3d>& rhs) const;


protected:

	Molecule mol;
	Nemo_complex_Parameters param;
    AtomicBasisSet aobasis;
    printleveler print_info=printleveler(2);

	/// standard calculation parameters; nemo_param: the nuclear correlation factor (nemo block, ncf)
    Nemo::NemoCalculationParameters nemo_param;
	CalculationParameters calc_param;

	/// magnetic Green's-function kernel (complex block: kernel), holding its parameters (magnetic_kernel
	/// block); rebuilt when the derived parameters change (recompute_factors_and_potentials, adapt_kernel_tmin)
	std::shared_ptr<MagneticKernel> mkernel;
	/// magnetic kernels: sum of the kernel energies - (J - K) + E_nuc, the energy of the original drivers
	double kernel_energy=0.0;

	/// restart bookkeeping as in SCF: the energy and the precision the orbitals are converged to
	double current_energy=1.e10;
	double converged_for_thresh=1.e10;
	double converged_for_dconv=1.e10;

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

		using FunctionFunctorInterface<double_complex,3>::operator();

		double_complex operator()(const coord_3d& xyz1) const override {
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

		using FunctionFunctorInterface<double_complex,3>::operator();

		double_complex operator()(const coord_3d& xyz1) const override {
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


	public:
	void test_compute_current_density() const;
};




} // namespace madness
#endif /* SRC_APPS_CHEM_ZNEMO_H_ */
