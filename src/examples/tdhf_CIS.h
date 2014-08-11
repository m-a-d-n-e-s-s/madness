/*
 * tdhfCIS.h
 *
 *  Created on: May 5, 2014
 *      Author: kottmanj
 */

#ifndef TDHFCIS_H_
#define TDHFCIS_H_

#include <chem/projector.h>
//#include "mp2.h"

#include "nonlinsol.h"
#include <chem/SCF.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>




using namespace madness;



typedef std::vector<Function<double,3> > vecfuncT;
struct noise : public FunctionFunctorInterface<double,3> {

typedef std::shared_ptr<FunctionFunctorInterface<double,3> > functorT;

public:
	noise(double size){noise(size,1.0);}
	noise(double size,double width)  {
	 Random random(wall_time()*1.e8);
	 x=size*(random.get()-random.get());
	 y=size*(random.get()-random.get());
	 z=size*(random.get()-random.get());
	 prefactor= x*y*z/(fabs(x*y*z));
	 sigma=width;
	// print("Noise constructor x,y,z,pre,sigma are: ",x,y,z,prefactor,sigma);
	}


private:
	real_function_3d gauss_function;
	double x,y,z;
	double sigma;
	double prefactor;

public:
	double operator()(const coord_3d &r) const{
		double r2 = (r[0]-x)*(r[0]-x)+(r[1]-y)*(r[1]-y)+(r[2]-z)*(r[2]-z);
		return prefactor*exp(-r2/(2.0*sigma*sigma));
	}

};


/// POD holding excitation energy and response vector for a single excitation
struct root {
	root(World& world) : world(world), omega(0.0),converged(false),err(10.0),delta(10.0) {}
	root(World& world,vecfuncT& x, double omega) :world(world), x(x), omega(omega),converged(false),err(10.0),delta(10.0) {}
	root(World& world, const vecfuncT& x1) : world(world), x(x1),converged(false),err(10.0),delta(10.0){}
    root(const root& other) : world(other.world), x(other.x), omega(other.omega), converged(other.converged),err(other.err),delta(other.delta){}

	World &world;
	vecfuncT x;
	double omega;
	bool converged;
	double err;
	double delta;
	//solverT solver; // Maybe later
	std::vector<double> amplitudes_;

	// Operators needed by the nonlinear solver
    root& operator=(const root& other) {
    	x=other.x;
    	// This is new, needed for sort function
    	// Not shure if this changes something when KAIN is used again
    	omega=other.omega;
    	converged=other.converged;
    	err=other.err;
    	delta=other.delta;
    	amplitudes_=other.amplitudes_;
    	return *this;
    }

    root operator-(const root& b) {
        return root(world,sub(world,x,b.x));
    }

    root operator+=(const root& b) { // Operator+= necessary
    	x=add(world,x,b.x);
    	return *this;
    }

    root operator*(double a) { // Scale by a constant necessary
    	scale(world,x,a);
        return *this;
    }
    // An operator for the sort function
    bool operator<(const root &other){
    	return (this->omega<other.omega);
    }
    bool operator>(const root &other){
    	return (this->omega>other.omega);
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
    allocator(World& world, const int nn) : world(world), n(nn) {}

    /// allocate a vector of n empty functions
    root operator()() {
        return root(world,zero_functions<double,3>(world,n));
    }
};

/// for convenience
static double inner(const root& a, const root& b) {
	if (a.x.size()==0) return 0.0;
	return madness::inner(a.x[0].world(),a.x,b.x).sum();
}

/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};

namespace madness {




/// The CIS class holds all machinery to compute excited state properties
class CIS {

	typedef SeparatedConvolution<double,3> operatorT;
	typedef std::shared_ptr<operatorT> poperatorT;
	typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

	struct prec{
		double thresh;
		double econv;
		double dconv;
		bool fock;
	};

private:
	#define TRUE  1
	#define FALSE 0

	static double rfunction(const coord_3d& r) {
	    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	}

	static double rfunction2(const coord_3d& r) {
	    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])+1.0;
	}

	static double x(const coord_3d &r){return r[0];}
	static double x2(const coord_3d &r){return r[0]-0.7;}
	static double x3(const coord_3d &r){return r[0]+0.7;}
	static double y(const coord_3d &r){return r[1];}
	static double z(const coord_3d &r){return r[2];}


	// Random number Generator from Madness lib misc/ran.h
	static double random_number(const coord_3d &r){
		Random random;
		return random.get();
	}

	// For the sorting of the roots
	static bool compare_roots(const root &a,const root &b){
		return a.omega<b.omega;
	}

public:

	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	hf		the HartreeFock reference state
	/// @param[in]	input	the input file name
	CIS(World& world, const SCF& calc, const std::string input)
: world(world),
  calc_(calc),
  guess_("physical"),
  nroot_(8),
  nfreeze_(0),
  econv_(calc.param.econv),
  dconv_(calc.param.econv),
  guess_thresh_(1.e-4),
  guess_econv_(0.001),
  guess_dconv_(0.03),
  guess_iter_(20),
  guess_iter_fock_(1),
  guess_roots_(8),
  guess_mode_("all_orbitals"),
  thresh_(dconv_*0.01),
  bsh_eps_(1.e-6),
  iter_max_(50),
  print_grid_(false),
  plot_(false),
  fixed_point_(false),

  guess_save_(false),
  guess_damp_(false),
  guess_damp_iter_(6),
  guess_damp_fock_(false),
  noise_(false),
  noise_box_(get_calc().molecule.bounding_cube()),
  noise_comp_(0.01),
  noise_width_(3.5),
  noise_gaussnumber_(25),
  hf_(false),
  triplet_(false),
  read_and_save_koala_(false),
  guess_exf_("dipole"){

		size_t nmo = get_calc().amo.size();
		guess_omega_=-0.9*get_calc().aeps(nmo-1);
		omega_=std::vector<double>(9,100.0);


		std::ifstream f(input.c_str());
		position_stream(f, "CIS");
		std::string s, tag;
		while (std::getline(f,s)) {
			std::istringstream ss(s);
			ss >> tag;
			if (tag == "end") break;
			else if (tag == "guess") ss >> guess_;
			else if (tag == "guess_preopt") ss >> preopt_;
			else if (tag == "nroot") ss >> nroot_;
			// dangerous right now else if (tag == "guess_roots") ss >> guess_roots_;
			else if (tag == "guess_thresh") ss >> guess_thresh_;
			else if (tag == "guess_econv") ss >> guess_econv_;
			else if (tag == "guess_dconv") ss >> guess_dconv_;
			else if (tag == "guess_iter") ss >> guess_iter_;
			else if (tag == "guess_iter_fock") ss >> guess_iter_fock_;
			else if (tag == "guess_roots") ss >> guess_roots_;
			else if (tag == "guess_omega") ss >> guess_omega_;
			else if (tag == "guess_mode") ss >> guess_mode_;
			else if (tag == "thresh") ss >> thresh_;
			else if (tag == "bsh_eps") ss >> bsh_eps_;
			else if (tag == "iter_max") ss >> iter_max_;
			else if (tag == "freeze") ss >> nfreeze_;
			else if (tag == "econv") ss >> econv_;
			else if (tag == "dconv") ss >> dconv_;
			else if (tag == "fixed_point") fixed_point_=true;
			else if (tag == "print_grid") print_grid_=true;
			else if (tag == "plot") plot_=true;
			else if (tag == "guess_save") guess_save_=true;
			else if (tag == "guess_damp") guess_damp_=true;
			else if (tag == "guess_damp_iter") ss >> guess_damp_iter_;
			else if (tag == "guess_damp_fock") guess_damp_fock_ = true;
			else if (tag == "noise") noise_=true;
			else if (tag == "noise_box") ss >> noise_box_;
			else if (tag == "noise_comp") ss >> noise_comp_;
			else if (tag == "noise_width") ss >> noise_width_;
			else if (tag == "noise_gaussnumber") ss >> noise_gaussnumber_;
			else if (tag == "hf") hf_=true;
			else if (tag == "triplet") triplet_=true;
			else if (tag == "koala_read_and_save") read_and_save_koala_=true;
			else if (tag == "guess_exf") ss >> guess_exf_;
			else if (tag == "omega0") ss >> omega_[0];
			else if (tag == "omega1") ss >> omega_[1];
			else if (tag == "omega2") ss >> omega_[2];
			else if (tag == "omega3") ss >> omega_[3];
			else if (tag == "omega4") ss >> omega_[4];
			else if (tag == "omega5") ss >> omega_[5];
			else if (tag == "omega6") ss >> omega_[6];
			else if (tag == "omega7") ss >> omega_[7];
			else if (tag == "omega8") ss >> omega_[8];
			else continue;
		}






		if (world.rank() == 0) {
			madness::print("\n ======= CIS info =======\n");
			if (nfreeze_==0) madness::print("   # frozen orbitals ","none");
			if (nfreeze_>0) madness::print("   # frozen orbitals ",0, " to ",nfreeze_-1);
			madness::print("     active orbitals ", nfreeze_," to ",get_calc().param.nalpha-1);

			madness::print("          guess from ", guess_);
			madness::print("        threshold 3D ", FunctionDefaults<3>::get_thresh());
			madness::print("  energy convergence ", econv_);
			madness::print("max residual (dconv) ", dconv_);
			madness::print("     number of roots ", nroot_);
			madness::print(" omega ", omega_[0],omega_[1],omega_[2]);

		}

		lo=get_calc().param.lo;
	}

	void set_prot(World &world,prec &a){
		FunctionDefaults<3>::set_thresh(a.thresh);
		econv_ = a.econv;
		dconv_ = a.dconv;
	}

	/// return the HF reference
	const SCF& get_calc() const {return calc_;}

	// print information of root or root vector
	void print_roots(const std::vector<root> &roots) const;
	void print_roots(const std::vector<root> &roots,const int iter) const;
	void print_root(const root &root) const;

	/// If roots[j] < roots[i] and i<j then the roots will switch places
	void sort_roots(std::vector<root> & roots);

	/// solve the CIS equations for n roots
	void solve();

	// Internal solver for parallel or sequential optimization
	bool solve_internal_par(const std::string mode,std::vector<root> &roots,const int iter_max);

	/// return the roots of the response equation
	std::vector<root>& roots();

	/// are we supposed to print the grid for an external guess
	bool print_grid() const;

private:

	/// the world
	World& world;

	/// the HartreeFock reference state
	const SCF& calc_;

	/// the excited states aka the roots of the response equation
	std::vector<root> roots_;

	/// intermediate for the two-electron interaction term, Eq. (8)

	/// the intermediate is the same for all roots:
	/// \[
	///   int[p,i] = \int 1/r12 \phi_i(1) * \phi_p(1)
	/// \]
	/// with p \in noct, i \in nocc
	std::vector<vecfuncT> exchange_intermediate_;

	/// the coulomb potential
	mutable real_function_3d coulomb_;


	/// where we get our guess from (all_virtual, koala)
	std::string guess_;

	/// Preoptimization mode
	std::string preopt_;

	/// number of roots we are supposed to solve
	int nroot_;

	/// number of frozen orbitals
	int nfreeze_;

	/// energy convergence threshold
	double econv_;

	/// density convergence threshold (=residual)
	double dconv_;

	/// the phases of the guess and ours might differ
	Tensor<double> guess_phases_;

	// Thresh for guess calculation
	double guess_thresh_;
	/// Energy convergence criterium for guess
	double guess_econv_;

	/// Wavefunction convergence for guess
	double guess_dconv_;

	/// how many iterations for guess roots
	int guess_iter_;
	int guess_iter_fock_;

	/// number of roots to guess
	int guess_roots_;

	/// Guess excitation energy
	double guess_omega_;

	/// Creating the guess with ... (all_orbitals, mo or homo) * (x,y,z,r)
	std::string guess_mode_;

	/// Thresh for CIS calculation
	double thresh_;

	/// Thresh for BSH Operator (std 1.e-6)
	double bsh_eps_;

	/// Maximal number of iterations
	int iter_max_;

	/// guess for the excitation energies
	std::vector<double> omega_;


	/// flag if the grid for the density should be printed

	/// external programs (e.g. Koala) need this grid to export the guess roots
	bool print_grid_;

	/// Plot functions in each iteration step
	bool plot_;

	/// perform a fixed-point iteration of the given excitation energies
	bool fixed_point_;

	/// Save the guess roots
	bool guess_save_;

	/// Damp the first guess_damp_iter guess interations
	bool guess_damp_;
	int guess_damp_iter_;
	bool guess_damp_fock_;
	/// Add noise to the guess, noise_comp_ is the prefactor to make it small
	/// noise_width is the width for the gauss functions
	/// noise_gaussnumber_ is the number of gauss functions computed for the noise
	bool noise_;
	double noise_box_;
	double noise_comp_;
	double noise_width_;
	double noise_gaussnumber_;

	/// Compute only Hartree Fock (for later restarts)
	bool hf_;

	/// Compute triplets (when false singlets are computed, default is false)
	bool triplet_;

	/// Just read and save the koala guess
	bool read_and_save_koala_;

	/// Dipole or quadrupole guess
	std::string guess_exf_;

	double lo;

	/// Initialize the roots
	/// Depending on the keyword guess_MO, guess_read, guess_koala or guess_physical will be called
	void initialize_roots(World &world,std::vector<root> &roots);

	void guess_MO(World &world,std::vector<root> &roots);
	void guess_read(World &world,std::vector<root> &roots);
	void guess_koala(World &world,std::vector<root> &roots);
	void guess_physical(World &world,std::vector<root> &roots);

	void add_noise(World & world,std::vector<root> &roots)const;

	std::vector<root> guess();

	// Guess Function
	std::vector<root> guess_big(const std::string exf);

	// Excitation functions for the guess (x,y,z,r,x^2,xy,xz ...)
	std::vector<real_function_3d> excitation_functions(const std::string exf) const;

	root guess_seq(int iroot) const;

	/// note that the orbitals from an external guess might be rotated wrt the
	/// MRA orbitals, which would spoil the guess. Therefore we rotate the
	/// guess with the overlap matrix so that the guess conforms with the MRA
	/// phases. The overlap to the external orbitals (if applicable) is computed
	/// only once
	/// @param[in]	iroot	guess for root iroot
	root guess_amplitudes(const int iroot);

	bool orthogonalize_fock(World& world,std::vector<root> &roots) const;
	bool orthogonalize_fock(World& world,std::vector<root> &roots,int iter) const;

	/// solve the CIS equations in a sequential way
	template<typename solverT>
	bool iterate_sequential(World& world, std::vector<solverT> & solver,std::vector<root> & guess_roots) const;

	/// solve the CIS equations for all roots

	/// @param[in]	world	the world
	/// @param[in]	solver	the KAIN solver (unused right now)
	/// @param[inout]	roots	on entry: guess for the roots
	///                         on successful exit: converged root
	/// @param[guess] for guess optimization with shifted HF
	/// @return	convergence reached or not
	template<typename solverT>
	bool iterate_all_CIS_roots(World& world, std::vector<solverT>& solver,
			std::vector<root>& roots,const std::string mode,const int iter_max) const;


	/// iterate the TDHF or CIS equations

	/// follow Eq (4) of
	/// T. Yanai, R. J. Harrison, and N. Handy,
	/// Multiresolution quantum chemistry in multiwavelet bases: time-dependent
	/// density functional theory with asymptotically corrected potentials in
	/// local density and generalized gradient approximations,
	/// Mol. Phys., vol. 103, no. 2, pp. 413-424, 2005.
	///
	/// The convergence criterion is that the excitation amplitudes don't change
	/// @param[in]		world	the world
	/// @param[in]		solver	the KAIN solver (not used right now..)
	/// @param[inout]	root	the current root that we solve
	/// @return			the residual error

	bool check_convergence(std::vector<root> &roots)const;

	template<typename solverT>
	double iterate_one_CIS_root(World& world, solverT& solver, root& thisroot,const std::string mode)const;

	/// iterate the TDHF or CIS equations -- deprecated !!

	/// follow Eq (4) of
	/// T. Yanai, R. J. Harrison, and N. Handy,
	/// Multiresolution quantum chemistry in multiwavelet bases: time-dependent
	/// density functional theory with asymptotically corrected potentials in
	/// local density and generalized gradient approximations,
	/// Mol. Phys., vol. 103, no. 2, pp. 413-424, 2005.
	///
	/// The convergence criterion is that the excitation amplitudes don't change
	/// @param[in]		world	the world
	/// @param[in]		solver	the KAIN solver
	/// @param[inout]	root	the current root that we solve
	/// @param[in]		excited_states all lower-lying excited states for orthog.
	/// @param[in]		iteration	the current iteration
	/// @return			if the iteration on this root has converged

	template<typename solverT>
	bool iterate_CIS(World& world, solverT& solver, root& thisroot,
			const std::vector<root >& excited_states, const int iteration) const;

	/// apply the gamma potential of Eq. (6) on the x vector

	/// @param[in]	x		the response amplitudes
	/// @param[in]	act		the active orbitals p
	/// @param[in]	rho0	the projector on all (frozen & active) MOs
	/// @return		(1-\rho^0) \Gamma \phi_p
	vecfuncT apply_gamma(const vecfuncT& x, const vecfuncT& act,
			const Projector<double,3>& rho0) const ;

	/// make the zeroth-order potential J-K+V_nuc

	/// note the use of all (frozen & active) orbitals in the computation
	/// @param[in]	x	the response vector
	/// @return		(J-K+V_nuc) x
	vecfuncT apply_fock_potential(const vecfuncT& x) const;

	/// return the Coulomb potential
	real_function_3d get_coulomb_potential() const;
	/// make the 2-electron interaction intermediate

	/// the intermediate is the same for all roots:
	/// \f[
	///   Int[i,p] = \int \frac{1}{r_{12}} \phi_i(1) * \phi_p(1)
	/// \f]
	/// both i and p are active MOs
	/// @param[in]	active_mo	active orbitals in the CIS computation
	/// @param[in]	amo			all MOs of the HF calculation
	/// @return		a vector of vectors of functions: [noct][nocc]
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT& active_mo,
			const vecfuncT& amo) const;

	/// compute the perturbed fock matrix

	/// the matrix is given by
	/// \f[
	///   F_{pq} = < x^p_i | F | x^q_i > + < x^p_i | \Gamma^q | \phi_i >
	/// \f]
	/// where an amplitude is given by its components x_i, such that and
	/// summation over the occupied orbitals i is implied
	/// \f[
	///   < \vec x^p | \vec x^q > = \sum_i x^p_i x^q_i
	/// \f]
	/// and similar for the Fock matrix elements.
	/// Note this is NOT the response matrix A
	Tensor<double> make_perturbed_fock_matrix(const std::vector<root>& roots,
			const vecfuncT& act, const Projector<double,3>& rho0) const;

	/// load a converged root from disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[out]	omega	the excitation energy
	/// @return	successfully loaded a root or not
	/// @param[in] filename_end end of filename (for guess roots)
	bool load_root(World& world, const int i, root& root)  const;
	bool load_root(World& world, const int i, root& root,const std::string filename_end)  const;

	/// save a converged root to disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[in]	omega	the excitation energy
	/// @param[in] filename for guess roots to be saved
	void save_root(World& world, const int i, const root& root) const;
	void save_root(World& world, const int i, const root& root, const std::string filename_end) const;

	/// normalize the excitation amplitudes

	/// normalize the set of excitation amplitudes, such that the sum of square
	/// of all amplitudes equals 1.
	/// @param[in]		world the world
	/// @param[inout]	x	the excitation vector
	void normalize(World& world, root& x) const;

	/// orthonormalize all roots using the perturbed fock-matrix (energy update included)
	void orthonormalize_fock(World &world,std::vector<root> &roots)const;

	/// orthonormalize all roots using Gram-Schmidt
	void orthonormalize(World& world, std::vector<root>& roots) const;

	/// compute the overlap between 2 sets of roots
	Tensor<double> overlap(const std::vector<root>& r1,
			const std::vector<root>& r2) const;


	/// compute the oscillator strength in the length representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_length(const root& root) const;

	/// compute the oscillator strength in the velocity representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_velocity(const root& root) const;

	/// analyze the root: oscillator strength and contributions from occ
	void analyze(const std::vector<root>& roots) const;

	/// return the active MOs only, note the shallow copy
	const vecfuncT active_mo() const;

};

}




#endif /* TDHFCIS_H_ */
