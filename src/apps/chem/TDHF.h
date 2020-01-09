/*
 * TDHF.h
 *
 *  Created on: Aug 11, 2016
 *      Author: kottmanj
 */

#ifndef SRC_APPS_CHEM_TDHF_H_
#define SRC_APPS_CHEM_TDHF_H_

#include "CCStructures.h"
#include "nemo.h"
#include "projector.h"
#include "SCFOperators.h"
#include <math.h>
#include "GuessFactory.h"



namespace madness {

/// The TDHF class
/// solves CIS/TDA equations and hopefully soon the full TDHF/TDDFT equations
class TDHF{
public:

	/// the TDHF parameter class
	struct Parameters : public QCCalculationParametersBase {

		Parameters() {
			initialize_all();
		}

		Parameters(const Parameters& other) : QCCalculationParametersBase(other) {}

		/// todo: read_from_file compatible with dist. memory computation
		Parameters(World& world, const std::shared_ptr<SCF>& scf, const std::string& input) {
			initialize_all();
			read(world,input,"response");
			set_derived_values(scf);
		}

		void initialize_all() {

			// MRA stuff
			initialize<double>("thresh",1.e-5);
			initialize<double>("econv",1.e-5);
			initialize<double>("dconv",1.e-4);

			// physics
			initialize<std::string>("calculation","cis","currently only cis=tda possible, TBD: thdf",{"cis"});
			initialize<bool>("triplet",false,"calculate triplet excitation energies (only works for CIS)");
			initialize<bool>("do_oep",false,"use OEP potentials for the ground state exchange");
			initialize<std::size_t>("excitations",1);
			initialize<std::size_t>("freeze",0,"the number of frozen occupied orbitals");
			initialize<std::string>("irrep","all","compute only irreps of the respective point group");

			// solver
			initialize<size_t>("maxiter",25,"maximum number of iterations in the final iterations");
			initialize<std::size_t>("kain_subspace",8,"use kain (kain subspace<=0 == no kain, kain_subspace==1 should have the same effect)");

			// guess
			initialize<double>("guess_econv",1.e-4);
			initialize<double>("guess_dconv",1.e-3);
			initialize<std::size_t>("iterating_excitations",2);
			initialize<std::size_t>("guess_excitations",4);
			initialize<std::size_t>("guess_occ_to_virt",5);
			initialize<double>("damping_width",0.0,"every exop is multiplied with e^(-exponent*r2) to avoid noise at the boundaries");

			initialize<bool>("debug",false);
			initialize<bool>("plot",false);
			initialize<bool>("no_compute",false);

			initialize<std::vector<size_t> >("restart",std::vector<size_t>(),"excitations which will be read from disk");
			initialize<double> ("lo",1.e10,"smallest length scale we need to resolve");

			initialize<std::string>("guess_excitation_operators","dipole+","guess typ",{"dipole+","quadrupole","big_fock_2","big_fock_3","big_fock_4","custom"});

			/// add center of mass functions determined by the homo-energy
			/// will add s,px,py,pz functions in the center of mass with exponent: -(e_homo/c) and c=guess_cm is the value of this parameter
			initialize<double>("guess_cm",2.0,"center-of-mass functions, s/p shell with exponent -(e_homo/c)");

			/// use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)
			/// much faster
			initialize<bool>("guess_diag",true,"use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)");

			/// determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess
			/// guess_active_orbitals=0 is the same as guess_diag
			initialize<std::size_t>("guess_active_orbitals",0,"determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess");


			initialize<bool>("store_potential",true,"store the potential for orthogonalizations or recalculate it");

			initialize<size_t>("guess_maxiter",0,"maximum number of guess iterations ");

			//		/// determine how the virtuals for the guess are constructed: scf, external, custom, dipole, quadrupole
			//		/// scf: read in the ao set from scf (scales with system size)
			//		/// external: read in virtuals from disk
			//		/// custom or predefined strings like dipole, dipole+, ... : create virtuals form occupied orbitals by multiplying with polynomials
			//		/// |v> = |occ>*poly
			//		/// if custom is chosen:
			//		/// the polynomials can be determined by: exop x n1 y n2 z n3 c n4, x n11 ... x n44, ...  which will be n4*x^n1*y^n2*z^n3 + n44*x^n11* ...
			//		/// e.g. for a dipole guess enter the exop keyword 3 times as:
			//		/// exop x 1.0
			//		/// exop y 1.0
			//		/// exop z 1.0
			//		/// the options dipole, dipole+, dipole+diffuse and quadrupole give predefined exops without explicitly stating them
			//		/// see the end of TDHF.cc for predefined keys
			//		std::string guess_excitation_operators="dipole+";
			//
			/// Vector of strings which contains the polynomial excitation operators
			/// For this to be used the tda_guess key has to be "custom"
			/// The strings are given in a format like: "c c1 x x1 y y1 z z1, c c2 x x2 y y2 z z2, ..." which will be interpreted as: c1*x^x1*y^y1*z^z1 + c2*x^x2*y^y2*z^z2 + ....
			initialize<std::vector<std::string> >("exops",{""},"applies only if guess_excitation_operator is custom");



		}

		void set_derived_values(const std::shared_ptr<SCF>& scf);

		// physical part
		std::size_t excitations() const {return get<std::size_t>("excitations");}
		std::size_t freeze() const {return get<std::size_t>("freeze");}
		std::string irrep() const {return get<std::string>("irrep");}
		bool triplet() const {return get<bool>("triplet");}
		bool do_oep() const {return get<bool>("do_oep");}

		// precision
		double thresh() const {return get<double>("thresh");}
		double econv() const {return get<double>("econv");}
		double dconv() const {return get<double>("dconv");}
		double lo() const {return get<double>("lo");}

		// restart and plotting
		bool debug() const {return get<bool>("debug");}
		bool no_compute() const {return get<bool>("no_compute");}
		std::vector<size_t>  restart() const {return get<std::vector<size_t> >("restart");}
		bool plot() const {return get<bool>("plot");}

		// solver parameters
		std::size_t iterating_excitations() const {return get<std::size_t>("iterating_excitations");}
		std::size_t maxiter() const {return get<std::size_t>("maxiter");}
		std::size_t kain_subspace() const {return get<std::size_t>("kain_subspace");}
		bool store_potential() const {return get<bool>("store_potential");}

		// guess parameters
		std::size_t guess_occ_to_virt() const {return get<std::size_t>("guess_occ_to_virt");}
		std::vector<std::string> exops() const {return get<std::vector<std::string> >("exops");}
		std::size_t guess_active_orbitals() const {return get<std::size_t>("guess_active_orbitals");}
		bool guess_diag() const {return get<bool>("guess_diag");}
		std::size_t guess_excitations() const {return get<std::size_t>("guess_excitations");}
		std::string guess_excitation_operators() const {return get<std::string>("guess_excitation_operators");}
		double damping_width() const {return get<double>("damping_width");}
		double guess_cm() const {return get<double>("guess_cm");}
		double guess_econv() const {return get<double>("guess_econv");}
		double guess_dconv() const {return get<double>("guess_dconv");}
		std::size_t guess_maxiter() const {return get<std::size_t>("guess_maxiter");}

		/// make parameters for convolution operator
		typename CCConvolutionOperator::Parameters get_ccc_parameters()const{
			typename CCConvolutionOperator::Parameters result;
			result.freeze=freeze();
			result.lo=lo();
			result.thresh_op=thresh();
			result.gamma=1.0;
			return result;
		}
	}; // end of parameter class

	TDHF(World & world,const Nemo &nemo, const std::string& input="input");
	TDHF(World & world,const Nemo &nemo, const Parameters& param);
	virtual
	~TDHF() {};


	/// check consistency of the input parameters
	void check_consistency() const;

	/// plot planes and cubes
	void plot(const vector_real_function_3d& vf, const std::string& name)const;

	/// sort the xfunctions according to their excitation energy and name the excitation energies accordingly
	std::vector<CC_vecfunction> sort_xfunctions(std::vector<CC_vecfunction> x)const;

	/// print information
	void print_xfunctions(const std::vector<CC_vecfunction> & f, const bool& fullinfo=false)const;

	/// Initialize the CIS functions

	/// @param[in\out] on input the already obtained guess functions (or empty vector), on output new guess functions are added
	void initialize(std::vector<CC_vecfunction> &start)const;

	void symmetrize(std::vector<CC_vecfunction>& v) const;

	/// Solve the CIS equations

	/// @param[in/out] CC_vecfunction
	/// on input the guess functions (if empty or not enough the a guess will be generated)
	/// on output the solution
	std::vector<CC_vecfunction> solve_cis()const;

	std::vector<CC_vecfunction> solve_cis(std::vector<CC_vecfunction>& start) const;

	/// Solve TDHF equations (not ready)
	void solve_tdhf(std::vector<CC_vecfunction>& guess)const;
	/// iterate the CIS guess vectors
	/// @param[in,out] x: on input the guess, on output the iterated guess
	/// see CC_Structures.h CCParameters class for convergence criteria
	bool iterate_cis_guess_vectors(std::vector<CC_vecfunction> &x)const;
	/// iterate the final CIS vectors
	/// @param[in,out] x: on input the guess, on output the iterated guess
	/// see CC_Structures.h CCParameters class for convergence criteria
	bool iterate_cis_final_vectors(std::vector<CC_vecfunction> &x)const;
	/// General function to iterate vectors
	/// @param[in,out] x: the CIS (or TDHF x) functions
	/// @param[in,out] the TDHF y functions (empty for CIS)
	/// @param[in] iterate_y, if true the y equation for TDHF is iterated
	/// @param[in] dconv: wavefunction convergence (for the vector norm of the vectorfunction)
	/// @param[in] econv: Energy convergece
	/// @param[in] iter: maximum number of iterations
	/// @param[in] kain: use kain if true (kainsubspace is controlled over CCParameters class)
	bool iterate_vectors(std::vector<CC_vecfunction> &x,const std::vector<CC_vecfunction> &y,bool iterate_y,const double dconv, const double econv, const double iter, const bool kain)const;
	/// Apply the Greens function to a vector of vectorfunction with a given potential
	/// @param[in] x: the vector of vectorfunctions where G will be applied to
	/// @param[in] V: the vector of potentials to the vectorfunctions, will be cleared afterwards (potentials are all potentials excpet the nuclear: 2J - K + Q(2pJ - pK)
	/// @param[out] the vectorfunctions after G has been applied
	/// the energy is assumed to be stored in the CC_vecfunctions member omega
	/// the wavefunction error is stored in the CC_vecfunctions member current_error
	std::vector<vector_real_function_3d> apply_G(std::vector<CC_vecfunction> &x,std::vector<vector_real_function_3d> &V)const;
	/// Make the old CIS Guess
	/// the routine is now used to create virtuals
	std::vector<CC_vecfunction> make_old_guess(const vector_real_function_3d& f)const;

	/// Create a set of virtual orbitals for the initial guess
	vector_real_function_3d make_virtuals() const;

	/// multiply excitation operators defined in the parameters with the seed functions
	/// @param[in] the seeds, define the function which are multiplied by the excitation operators
	/// @param[in] use_trigo, if false polynomials are used for excitation operators, else trigonometric functions (i.e. x^2y vs sin^2(x)*sin(y))
	/// Trigonometric functions are prefered since they are bounded (no weird behaviour at the boundaries for large exponents)
	vector_real_function_3d apply_excitation_operators(const vector_real_function_3d& seed,const bool& use_trigo=true) const;

	/// make the initial guess by explicitly diagonalizing a CIS matrix with virtuals from the make_virtuals routine
	vector<CC_vecfunction> make_guess_from_initial_diagonalization() const;
	/// canonicalize a set of orbitals (here the virtuals for the guess)
	vector_real_function_3d canonicalize(const vector_real_function_3d& v, Tensor<double>& veps)const;

	/// compute the CIS matrix for a given set of virtuals

	/// @param[in]	virtuals	the virtual orbitals
	/// @param[in]	veps		the orbital energies of the virtuals
	Tensor<double> make_cis_matrix(const vector_real_function_3d virtuals, const Tensor<double>& veps)const;

	/// initialize the excitation functions
	bool
	initialize_singles(CC_vecfunction &singles,const FuncType type,const int ex) const;


	/// Make the potentials to a given vector of vecfunctions (excitations)
	/// @param[in] The vector of excitations
	/// @param[out] The potentials
	std::vector<vector_real_function_3d> make_potentials(const std::vector<CC_vecfunction> &x)const;
	//    /// Make the CIS potential for a single excitation vector
	//	vecfuncT get_cis_potential(const CC_vecfunction& x) const {
	//		return CCOPS.make_cis_potential(x);
	//	}
	/// Make the TDA potential for a single excitation vector
	vector_real_function_3d get_tda_potential(const CC_vecfunction &x)const;
	/// Make the TDHF potential (not ready)
	std::vector<vector_real_function_3d> make_tdhf_potentials(std::vector<CC_vecfunction> &x,const std::vector<CC_vecfunction> &y)const;
	/// orthonormalize a vector of excitations
	/// @param[in,out] input: the excitations, output: the orthonormalized excitations
	/// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
	/// output: the transformed potentials
	void orthonormalize(std::vector<CC_vecfunction> &x,std::vector<vector_real_function_3d> &V)const;
	/// Calculate the perturbed fock matrix for a given vector of excitations
	/// @param[in] input: the excitations
	/// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
	Tensor<double> make_perturbed_fock_matrix(const std::vector<CC_vecfunction> &x, const std::vector<vector_real_function_3d> &V)const;
	Tensor<double> make_overlap_matrix(const std::vector<CC_vecfunction> &x)const;
	std::vector<vector_real_function_3d> transform(const std::vector<vector_real_function_3d> &x,const madness::Tensor<double> U) const{
		std::vector<CC_vecfunction> tmp;
		for(const auto& xi:x) tmp.push_back(CC_vecfunction(xi));
		std::vector<CC_vecfunction> tmp2= transform(tmp,U);
		std::vector<vector_real_function_3d> result;
		for(const auto&xi:tmp2) result.push_back(xi.get_vecfunction());
		return result;
	}
	/// Interface to the SCF.h fock_transform function
	std::vector<CC_vecfunction> transform(const std::vector<CC_vecfunction> &x,const madness::Tensor<double> U) const;


	/// Helper function to initialize the const mo_bra and ket elements
	CC_vecfunction make_mo_bra(const Nemo &nemo) const {
		vector_real_function_3d tmp = mul(world, nemo.ncf->square(),
				nemo.get_calc()->amo);
		set_thresh(world, tmp, parameters.thresh());
		truncate(world,tmp);
		reconstruct(world,tmp);
		CC_vecfunction mo_bra(tmp, HOLE);
		return mo_bra;
	}

	CC_vecfunction make_mo_ket(const Nemo&nemo) const {
		vector_real_function_3d tmp = nemo.get_calc()->amo;
		set_thresh(world, tmp, parameters.thresh());
		truncate(world,tmp);
		reconstruct(world,tmp);
		CC_vecfunction mo_ket(tmp, HOLE);
		return mo_ket;
	}

	double get_orbital_energy(const size_t i)const{
		return nemo.get_calc()->aeps(i);
	}

	/// convenience
	vector_real_function_3d make_bra(const CC_vecfunction &ket)const{
		return make_bra(ket.get_vecfunction());
	}
	real_function_3d make_bra(const real_function_3d &ket)const{
		vector_real_function_3d v(1,ket);
		return make_bra(v).front();

	}
	/// maybe move this into nuclear_correlation class ?
	vector_real_function_3d make_bra(const vector_real_function_3d &ket)const{
		CCTimer time(world,"Make Bra");
		real_function_3d nucf = nemo.ncf->square();
		vector_real_function_3d result= mul(world,nucf,ket);
		time.info(parameters.debug());
		return result;
	}

	template<typename T, size_t NDIM>
	bool load_function(Function<T, NDIM>& f, const std::string name) const {
		bool exists = archive::ParallelInputArchive::exists(world,name.c_str());
		if(exists){
			if (world.rank() == 0) print("loading function", name);
			archive::ParallelInputArchive ar(world, name.c_str());
			ar & f;
			f.print_size(name);
			return true;
		}else return false;
	}

	const vector_real_function_3d get_active_mo_ket()const{
		vector_real_function_3d result;
		for(size_t i=parameters.freeze();i<mo_ket_.size();i++) result.push_back(mo_ket_(i).function);
		return result;
	}
	const vector_real_function_3d get_active_mo_bra()const{
		vector_real_function_3d result;
		for(size_t i=parameters.freeze();i<mo_ket_.size();i++) result.push_back(mo_bra_(i).function);
		return result;
	}

	/// compute the oscillator strength in the length representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]  root    a converged root
	double oscillator_strength_length(const CC_vecfunction& x) const;

	/// compute the oscillator strength in the velocity representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]  root    a converged root
	double oscillator_strength_velocity(const CC_vecfunction& x) const;

	/// analyze the root: oscillator strength and contributions from occupied orbitals
	void analyze(const std::vector<CC_vecfunction> &x) const;
	/// Fock matrix for occupied orbitals
	Tensor<double> F_occ;
	/// The MPI Communicator
	World& world;
	/// The Parameters for the Calculations
	const Parameters parameters;
	/// The Nemo structure (convenience)
	const Nemo& nemo;
	/// Operator Structure which can handle intermediates (use for exchange with GS orbitals)
	/// Can be replaced by another potential manager
	CCConvolutionOperator g12;
	/// MO bra and ket
	const CC_vecfunction mo_ket_;
	const CC_vecfunction mo_bra_;
	/// the Projector to the virtual space
	const QProjector<double,3> Q;
	/// the symmetry projector
	projector_irrep symmetry_projector;
	/// the messenger IO
	CCMessenger msg;
	/// converged roots
	mutable std::vector<CC_vecfunction> converged_roots;
	/// stored guess roots roots to feed into the cycle, sorted backwards for easier pop_back calling
	mutable std::vector<CC_vecfunction> guess_roots;
};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_TDHF_H_ */
