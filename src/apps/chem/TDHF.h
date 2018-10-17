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
	struct Parameters{

		// disable default constructor (need scf information)
		Parameters();
		Parameters(const std::shared_ptr<SCF>& scf);
		Parameters(const std::shared_ptr<SCF>& scf, const std::string& input);

		/// assign default for everything that has not been assigned from file:
		void complete_with_defaults(const std::shared_ptr<SCF>& scf);

		/// reads parameters in from file
		/// parameters are given between key and end keywords
		/// filename default is "input"
		void read_from_file(const std::string input, const std::string& key = "response");

		bool assign_bool(std::ifstream& s)const{
			std::string tmp;
			s>>tmp;
			return assign_bool(tmp);
		}
		bool assign_bool(std::string& s)const{
			if(s=="0" or s=="false" or s=="f") return false;
			else return true;
		}

		bool debug=false;
		bool plot=false;
		bool no_compute=false;

		/// full linear-response (TDHF/TDDFT) or Tamm-Dancoff (CIS/DFT-TDA)
		std::string calculation="cis"; // possible: cis and tdhf (dft part is controlled throught moldft parameters)

		/// the threshold for the CIS calculation
		double thresh=FunctionDefaults<3>::get_thresh();
		double thresh_op=-1.0; // i.e. uninitialized;
		///Smallest length scale to resolve
		double lo=1.e-10;

		/// make parameters for convolution operator
		typename CCConvolutionOperator::Parameters get_ccc_parameters()const{
			typename CCConvolutionOperator::Parameters result;
			result.freeze=freeze;
			result.lo=lo;
			result.thresh_op=thresh_op;
			result.gamma=1.0;
			return result;
		}

		/// the number of frozen occupied orbitals (not taken into account for response)
		int freeze=0;

		/// the number of frozen occupied orbitals (not taken into account for response)
		std::string irrep="all";

		/// excitations which will be read from disk
		std::vector<size_t> restart;

		/// The number of occupied orbitals which are used to create the virtuals for the guess
		/// <= 0 means that all orbitals will be used
		int guess_occ_to_virt=5;

		/// The number of excitation vectors for which the alorithm will solve
		int excitations=1;
		/// The number of guess_excitation vectors for the first iterations
		int guess_excitations=-1;
		/// The number of excitation vectors which will be iterated parallel
		int iterating_excitations=-1;

		/// determine how the virtuals for the guess are constructed: scf, external, custom, dipole, quadrupole
		/// scf: read in the ao set from scf (scales with system size)
		/// external: read in virtuals from disk
		/// custom or predefined strings like dipole, dipole+, ... : create virtuals form occupied orbitals by multiplying with polynomials
		/// |v> = |occ>*poly
		/// if custom is chosen:
		/// the polynomials can be determined by: exop x n1 y n2 z n3 c n4, x n11 ... x n44, ...  which will be n4*x^n1*y^n2*z^n3 + n44*x^n11* ...
		/// e.g. for a dipole guess enter the exop keyword 3 times as:
		/// exop x 1.0
		/// exop y 1.0
		/// exop z 1.0
		/// the options dipole, dipole+, dipole+diffuse and quadrupole give predefined exops without explicitly stating them
		/// see the end of TDHF.cc for predefined keys
		std::string guess_excitation_operators="dipole+";

		/// add center of mass functions determined by the homo-energy
		/// will add s,px,py,pz functions in the center of mass with exponent: -(e_homo/c) and c=guess_cm is the value of this parameter
		double guess_cm=2.0;

		/// use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)
		/// much faster
		bool guess_diag=true;

		/// determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess
		/// guess_active_orbitals=0 is the same as guess_diag
		int guess_active_orbitals=0;


		/// convergence for the excitation vectors
		double guess_dconv=-1.0;
		double dconv=-1.0;
		/// convergence for the excitation energy
		double guess_econv=-1.0;
		double econv=-1.0;

		/// store the potential for orthogonalizations or recalculate it
		bool store_potential=true;

		/// maximum number of iterations in the final iterations
		size_t maxiter=25;
		/// maximum number of guess iterations (pre iterate guess functions)
		size_t guess_maxiter=0;
		/// Vector of strings which contains the polynomial excitation operators
		/// For this to be used the tda_guess key has to be "custom"
		/// The strings are given in a format like: "c c1 x x1 y y1 z z1, c c2 x x2 y y2 z z2, ..." which will be interpreted as: c1*x^x1*y^y1*z^z1 + c2*x^x2*y^y2*z^z2 + ....
		std::vector<std::string> exops;
		/// smoothing exponent
		/// every exop is multiplied with e^(-exponent*r2) to avoid noise at the boundaries
		double damping_width=0.0;

		/// calculate triplet excitation energies (only works for CIS)
		bool triplet=false;

		/// use kain (kain subspace<=0 == no kain, kain_subspace==1 should have the same effect)
		int kain_subspace=8;

		/// general new key/value pair for convenience in development
		std::map<std::string,std::string> generalkeyval;

		void print(World& world) const;

	}; // end of parameter class

	TDHF(World & world,const Nemo &nemo, const std::string& input="input");
	TDHF(World & world,const Nemo &nemo, const Parameters& param);
	virtual
	~TDHF();


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

	void symmetrize(std::vector<CC_vecfunction>& v) const {

		std::vector<std::string> irreps, orbital_irreps, reduced;
		vector_real_function_3d bla=symmetry_projector(nemo.get_calc()->amo,nemo.R_square,orbital_irreps);
		if (nemo.do_symmetry()) {
			for (auto& f : v) {
				vector_real_function_3d tmp=symmetry_projector(f.get_vecfunction(),nemo.R_square,irreps);
				f.set_functions(tmp,RESPONSE,parameters.freeze);
				for (int i=0; i<irreps.size(); ++i) {
					reduced=symmetry_projector.reduce(irreps[i],orbital_irreps[i]);

					if (not ((reduced[0]==f.irrep) or (reduced[0]=="null"))) {
						print("reduced, irrep",reduced[0],f.irrep);
						MADNESS_ASSERT(0);
					}
				}
			}
		}
	}

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
	vector_real_function_3d canonicalize(const vector_real_function_3d& v)const;
	/// compute the CIS matrix for a given set of virtuals
	Tensor<double> make_cis_matrix(const vector_real_function_3d virtuals)const;

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
		vector_real_function_3d tmp = mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		set_thresh(world, tmp, parameters.thresh);
		truncate(world,tmp);
		reconstruct(world,tmp);
		CC_vecfunction mo_bra(tmp, HOLE);
		return mo_bra;
	}

	CC_vecfunction make_mo_ket(const Nemo&nemo) const {
		vector_real_function_3d tmp = nemo.get_calc()->amo;
		set_thresh(world, tmp, parameters.thresh);
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
		real_function_3d nucf = nemo.nuclear_correlation ->square();
		vector_real_function_3d result= mul(world,nucf,ket);
		time.info(parameters.debug);
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
		for(size_t i=parameters.freeze;i<mo_ket_.size();i++) result.push_back(mo_ket_(i).function);
		return result;
	}
	const vector_real_function_3d get_active_mo_bra()const{
		vector_real_function_3d result;
		for(size_t i=parameters.freeze;i<mo_ket_.size();i++) result.push_back(mo_bra_(i).function);
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
