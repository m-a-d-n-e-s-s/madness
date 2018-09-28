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

		// disable default constructor (need nemo information)
		Parameters();
		Parameters(const std::shared_ptr<SCF>& scf);
		Parameters(const std::shared_ptr<SCF>& scf, const std::string& input);

		/// assign default for everything that has not been assigned from file:
		void complete_with_defaults(const std::shared_ptr<SCF>& scf);

		/// reads parameters in from file
		/// parameters are given between key and end keywords
		/// filename default is "input"
		void read_from_file(const std::string input, const std::string& key = "response");

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
		std::string guess_virtuals="dipole+";

		/// add center of mass functions determined by the homo-energy
		/// will add s,px,py,pz functions in the center of mass with exponent: -(e_homo/c) and c=guess_cm is the value of this parameter
		double guess_cm=3.5;

		/// use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)
		/// much faster
		bool guess_diag=true;

		/// determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess
		/// guess_active_orbitals=0 is the same as guess_diag
		int guess_active_orbitals=0;


		/// convergence for the excitation vectors
		double guess_dconv=-1.0;
		double dconv=FunctionDefaults<3>::get_thresh()*10.0;
		/// convergence for the excitation energy
		double guess_econv=-1.0;
		double econv=FunctionDefaults<3>::get_thresh();

		/// store the potential for orthogonalizations or recalculate it
		bool store_potential=true;

		/// maximum number of iterations in the final iterations
		size_t maxiter=25;
		/// maximum number of guess iterations (mostly more than the final ones and always without KAIN)
		size_t guess_maxiter=10;
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

	/// plot planes and cubes
	void plot(const vecfuncT& vf, const std::string& name)const{
		if(parameters.plot){
		CCTimer timer(world,"plot planes and cubes");
		madness::plot(vf,name,nemo.get_calc()->molecule.cubefile_header());
		timer.print();
		}
	}

	/// print information
	void print_xfunctions(std::vector<CC_vecfunction> & f, const bool& fullinfo=false)const{
		for(const auto& x:f){
			const double mem=get_size(world,x.get_vecfunction());
			if(world.rank()==0){
				std::cout << "ex. vector " << x.excitation << " | " << x.omega << " | " << mem << "(Gbyte)" << "\n";
			}
			if(fullinfo){
				print_size(world,x.get_vecfunction(),"ex. vector "+std::to_string(x.excitation));
			}
		}
	}

	/// Initialize the CIS functions

	/// @param[in\out] on input the already obtained guess functions (or empty vector), on output new guess functions are added
	void initialize(std::vector<CC_vecfunction> &start)const;

	/// Solve the CIS equations

	/// @param[in/out] CC_vecfunction
	/// on input the guess functions (if empty or not enough the a guess will be generated)
	/// on output the solution
	std::vector<CC_vecfunction> solve_cis()const{
		std::vector<CC_vecfunction> ccs;
		// look for restart options
		for(size_t k=0;k<parameters.restart.size();k++){
			CC_vecfunction tmp;
			const bool found= initialize_singles(tmp,RESPONSE,parameters.restart[k]);
			if(found) ccs.push_back(tmp);
		}
		return solve_cis(ccs);
	}

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
	std::vector<vecfuncT> apply_G(std::vector<CC_vecfunction> &x,std::vector<vecfuncT> &V)const;
	/// Guess for TDHF y functions (not ready)
	std::vector<CC_vecfunction> make_y_guess(const std::vector<CC_vecfunction> & x, std::vector<CC_vecfunction> & y)const;
	/// Make the old CIS Guess
	/// the routine is now used to create virtuals
	std::vector<CC_vecfunction> make_old_guess(const vecfuncT& f)const;

	/// Create a set of virtual orbitals for the initial guess
	vecfuncT make_virtuals() const;

	/// multiply excitation operators defined in the parameters with the seed functions
	/// @param[in] the seeds, define the function which are multiplied by the excitation operators
	/// @param[in] use_trigo, if false polynomials are used for excitation operators, else trigonometric functions (i.e. x^2y vs sin^2(x)*sin(y))
	vecfuncT apply_excitation_operators(const vecfuncT& seed,const bool& use_trigo=true) const;

	/// make the initial guess by explicitly diagonalizing a CIS matrix with virtuals from the make_virtuals routine
	std::vector<CC_vecfunction> make_guess_from_initial_diagonalization() const;
	/// canonicalize a set of orbitals (here the virtuals for the guess)
	vecfuncT canonicalize(const vecfuncT& v)const{
		CCTimer time(world,"canonicalize");
		Fock F(world, &nemo);
		const vecfuncT vbra=make_bra(v);
		Tensor<double> Fmat = F(vbra,v);
		Tensor<double> S = matrix_inner(world, vbra, v);
		Tensor<double> occ(v.size());
		occ=1.0;
		Tensor<double> evals;
		if(parameters.debug and world.rank()==0) std::cout << "Canonicalize: Fock Matrix\n" << Fmat(Slice(0,std::min(10,int(v.size()))-1),Slice(0,std::min(10,int(v.size()))-1));
		if(parameters.debug and world.rank()==0) std::cout << "Canonicalize: Overlap Matrix\n" << S(Slice(0,std::min(10,int(v.size()))-1),Slice(0,std::min(10,int(v.size()))-1));
		Tensor<double> U = nemo.get_calc()->get_fock_transformation(world, S, Fmat, evals, occ, std::min(parameters.thresh,1.e-4));
		vecfuncT result = madness::transform(world, v, U);
		time.print();
		return result;
	}
	/// compute the CIS matrix for a given set of virtuals
	Tensor<double> make_cis_matrix(const vecfuncT virtuals)const{

		// make bra elements
		const vecfuncT virtuals_bra = make_bra(virtuals);
		// make Fock Matrix of virtuals for diagonal elements
		Fock F(world, &nemo);
		Tensor<double> Fmat = F(virtuals_bra, virtuals);


		if (parameters.debug) {
			const int dim = std::min(10,int(virtuals.size()));
			if (world.rank() == 0)
				std::cout << "Debug Part of Virtual Fock Matrix\n" << Fmat(Slice(0,dim-1),Slice(0,dim-1)) << "\n";

			Tensor<double> S = matrix_inner(world, virtuals_bra, virtuals);
			if (world.rank() == 0)
				std::cout << "Debug Overlap of virtuals\n" << S(Slice(0,dim-1),Slice(0,dim-1)) << "\n";
		}


		CCTimer time_cis(world, "make CIS matrix");

		// the cis matrix is indexed by ij and ab
		// we will use the combined indixes from ia and jb named I and J
		// in order to not be confused we use the following helper functions
		const int nocc = get_active_mo_ket().size();
		// determines for which orbitals (couting from the HOMO downwards) the off-diagonal elements will be computed
		// this simplifies the guess
		int active_guess_orbitals = parameters.guess_active_orbitals;
		const int nvirt = virtuals.size();
		auto get_com_idx = [nvirt](int i, int a) { return i*nvirt+a; };
		auto get_vir_idx = [nvirt](int I) {return I%nvirt;};
		auto get_occ_idx = [nvirt](int I) {return I/nvirt;};
		auto delta = [](int x, int y) {if (x==y) return 1; else return 0;};


		const int dim=(virtuals.size()*nocc);
		if(world.rank()==0) std::cout << "CIS-Matrix for guess calculation will be of size " << dim << "x" << dim << "\n";
		// the number of the matrix where elements which are not determined by orbital energies and the fock matrix are computed (controlled over active_guess_orbitals parameter)
		const int dim2=(virtuals.size()*active_guess_orbitals);
		if(dim2<dim and world.rank()==0) std::cout << "Effective size through neglect of some orbitals will be: " << dim2 << "x" << dim2 << "\n";
		const int start_ij = nocc-active_guess_orbitals;
		Tensor<double> MCIS(dim,dim);

		// make CIS matrix
		// first do the "diagonal" entries
		if(nemo.get_calc()->param.localize){
			Tensor<double> Focc = F(get_active_mo_bra(),get_active_mo_ket());
			for(int I=0;I<dim;++I){
				const int a=get_vir_idx(I);
				const int i=get_occ_idx(I);
				for(int J=0;J<dim;++J){
					const int b=get_vir_idx(J);
					const int j=get_occ_idx(J);
					MCIS(I,J) = Fmat(a,b)*delta(i,j)-Focc(i,j)*delta(a,b);
				}
			}
		}else{
		for(int I=0;I<dim;++I){
			const int a=get_vir_idx(I);
			const int i=get_occ_idx(I);
			MCIS(I,I) = Fmat(a,a)-get_orbital_energy(i+parameters.freeze);
		}
		}

		if(not parameters.guess_diag){
		int I = -1; // combined index from i and a, start is -1 so that initial value is 0 (not so important anymore since I dont use ++I)
		for (int i = start_ij; i < get_active_mo_ket().size(); ++i) {
			const real_function_3d brai = get_active_mo_bra()[i];
			const vecfuncT igv = g12(brai * virtuals);
			for (int a = 0; a < virtuals.size(); ++a) {
				I=get_com_idx(i,a);
				int J =-1;
				for (int j = start_ij; j < get_active_mo_ket().size(); ++j) {
					const real_function_3d braj =get_active_mo_bra()[j];
					for (int b = 0; b < virtuals.size(); ++b) {
						J=get_com_idx(j,b);
						if(J<=I){
							const real_function_3d igj = g12(mo_bra_(i+parameters.freeze),mo_ket_(j+parameters.freeze)); // use exchange intermediate
							const double rIJ = 2.0 * inner(braj * virtuals[b], igv[a]) - inner(virtuals_bra[a] * virtuals[b],igj);
							MCIS(J,I) += rIJ;
							MCIS(I,J) += rIJ;
						}
					}
				}
			}
		}
		}
		if(world.rank()==0){
			int sdim=std::min(int(MCIS.dim(0)),10);
			std::cout << "Part of the CIS Matrix:\n" << MCIS(Slice(dim-sdim,-1),Slice(dim-sdim,-1)) << "\n";
			if(parameters.debug) std::cout << "Debug: Full CIS Matrix:\n" << MCIS<< "\n";
		}

		// test if symmetric
		if (parameters.debug) {
			const double symm_norm = (MCIS - transpose(MCIS)).normf();
			if (world.rank() == 0)
				std::cout << "Hermiticity of CIS Matrix:\n" << "||MCIS-transpose(MCIS)||=" << symm_norm << "\n";

			if (symm_norm > 1.e-4) {
				int sliced_dim = 8;
				if (8 > MCIS.dim(0))
					sliced_dim = MCIS.dim(0);

				if (world.rank() == 0)
					std::cout << "first " << sliced_dim << "x" << sliced_dim << " block of MCIS Matrix\n" << MCIS(_, Slice(sliced_dim - 1, sliced_dim - 1));
			}
		}
		time_cis.info();

		return MCIS;
	}

	/// initialize the excitation functions
	bool
	initialize_singles(CC_vecfunction &singles,const FuncType type,const int ex) const;


	/// Make the potentials to a given vector of vecfunctions (excitations)
	/// @param[in] The vector of excitations
	/// @param[out] The potentials
	std::vector<vecfuncT> make_potentials(const std::vector<CC_vecfunction> &x)const;
	//    /// Make the CIS potential for a single excitation vector
	//	vecfuncT get_cis_potential(const CC_vecfunction& x) const {
	//		return CCOPS.make_cis_potential(x);
	//	}
	/// Make the TDA potential for a single excitation vector
	vecfuncT get_tda_potential(const CC_vecfunction &x)const;
	/// Make the TDHF potential (not ready)
	std::vector<vecfuncT> make_tdhf_potentials(std::vector<CC_vecfunction> &x,const std::vector<CC_vecfunction> &y)const;
	/// orthonormalize a vector of excitations
	/// @param[in,out] input: the excitations, output: the orthonormalized excitations
	/// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
	/// output: the transformed potentials
	void orthonormalize(std::vector<CC_vecfunction> &x,std::vector<vecfuncT> &V)const;
	/// Calculate the perturbed fock matrix for a given vector of excitations
	/// @param[in] input: the excitations
	/// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
	Tensor<double> make_perturbed_fock_matrix(const std::vector<CC_vecfunction> &x, const std::vector<vecfuncT> &V)const;
	Tensor<double> make_overlap_matrix(const std::vector<CC_vecfunction> &x)const;
	std::vector<vecfuncT> transform(const std::vector<vecfuncT> &x,const madness::Tensor<double> U) const{
		std::vector<CC_vecfunction> tmp;
		for(const auto& xi:x) tmp.push_back(CC_vecfunction(xi));
		std::vector<CC_vecfunction> tmp2= transform(tmp,U);
		std::vector<vecfuncT> result;
		for(const auto&xi:tmp2) result.push_back(xi.get_vecfunction());
		return result;
	}
	/// Interface to the SCF.h fock_transform function
	std::vector<CC_vecfunction> transform(const std::vector<CC_vecfunction> &x,const madness::Tensor<double> U) const;
	/// Make guess function strings for given key
	std::vector<std::string> make_predefined_guess_strings(const std::string what)const;
	/// Make guess functions strings  to given order
	std::vector<std::string> make_auto_polynom_guess(const size_t order)const;

	/// Helper function to initialize the const mo_bra and ket elements
	CC_vecfunction make_mo_bra(const Nemo &nemo) const {
		vecfuncT tmp = mul(world, nemo.nuclear_correlation->square(),
				nemo.get_calc()->amo);
		set_thresh(world, tmp, parameters.thresh);
		truncate(world,tmp);
		reconstruct(world,tmp);
		CC_vecfunction mo_bra(tmp, HOLE);
		return mo_bra;
	}

	CC_vecfunction make_mo_ket(const Nemo&nemo) const {
		vecfuncT tmp = nemo.get_calc()->amo;
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
	vecfuncT make_bra(const CC_vecfunction &ket)const{
		return make_bra(ket.get_vecfunction());
	}
	real_function_3d make_bra(const real_function_3d &ket)const{
		vecfuncT v(1,ket);
		return make_bra(v).front();

	}
	/// maybe move this into nuclear_correlation class ?
	vecfuncT make_bra(const vecfuncT &ket)const{
		CCTimer time(world,"Make Bra");
		real_function_3d nucf = nemo.nuclear_correlation ->square();
		vecfuncT result= mul(world,nucf,ket);
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

	const vecfuncT get_active_mo_ket()const{
		vecfuncT result;
		for(size_t i=parameters.freeze;i<mo_ket_.size();i++) result.push_back(mo_ket_(i).function);
		return result;
	}
	const vecfuncT get_active_mo_bra()const{
		vecfuncT result;
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
	/// the messenger IO
	CCMessenger msg;
};



} /* namespace madness */

#endif /* SRC_APPS_CHEM_TDHF_H_ */
