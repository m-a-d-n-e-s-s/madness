/*
 * PNO.h
 *
 *  Created on: Oct 22, 2018
 *      Author: kottmanj
 */

#ifndef PAPER_CODE_PNO_H_
#define PAPER_CODE_PNO_H_

// convenience macro
#define PAIRLOOP(it) for(ElectronPairIterator it=pit();it;++it)
#define TIMER(timer) MyTimer timer=MyTimer(world).start();

#include<madness/chem/PNOF12Potentials.h>
#include <madness/world/worldmem.h>

#include<madness/chem/CC2.h>
#include<madness/chem/QCPropertyInterface.h>
#include<madness/chem/molecule.h>
#include<madness/chem/PNOTensors.h>

namespace madness {
// needed to plot cubefiles with madness
extern std::vector<std::string> cubefile_header(std::string filename="input", const bool& no_orient=false);

class PNO : public QCPropertyInterface {
public:
	typedef std::shared_ptr<operatorT> poperatorT;
	PNO(World& world, const Nemo& nemo, const PNOParameters& parameters, const F12Parameters& paramf12)
	: world(world),
	  param(parameters),
	  nemo(nemo),
	  J(world, &nemo),
	  K(ParametrizedExchange(world, nemo, parameters.exchange())),
	  T(world),
	  V(world, nemo.ncf),
	  F(world, &nemo),
	  Q( nemo.get_calc()->amo),
	  basis(world,nemo.get_calc()->molecule,8),
	  f12(world,nemo,basis,paramf12),
	  msg(world)
	{
		poisson = std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(
				world, nemo.get_calc()->param.lo(),param.op_thresh()));
		MADNESS_ASSERT(param.freeze() == f12.param.freeze());
	}

    std::string name() const {return "PNO";};

    static void help() {
        print_header2("help page for PNO");
        print("The PNO code computes MP2 energies using pair natural orbitals");
        print("You can print all available calculation parameters by running\n");
        print("pno --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("pno --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        PNOParameters param;
        print("default parameters for the pno program are\n");
        param.print("pno", "end");

        F12Parameters f12param(param);
        print("\n\nadditional parameters for the correlation factor are\n");
        f12param.print("f12", "end");

        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    virtual bool selftest() {
        return true;
    };

	/// Compute the projected MP2 energies: 2<ij|g|uij> - <ji|g|uij> for all pairs
	/// gives back the PairEnergies structure with all singlet and triplet energies
	PairEnergies compute_projected_mp2_energies(PNOPairs& pairs)const;

	/// solves whatever was specified in input
	/// will call solve_mp2, solve_cispd etc
	void solve()const{
		std::vector<PNOPairs> dummy;
		solve(dummy);
	}
	/// Solve for the PairType that is given
	void solve(std::vector<PNOPairs>& all_pairs) const;

	/// interfaces the BasisFunctions class
	/// guesses a set of virtuals depending on the guesstype and what was specified in the parameters
	vector_real_function_3d guess_virtuals(const vector_real_function_3d& f = vector_real_function_3d(), const GuessType& inpgt = UNKNOWN_GUESSTYPE) const;

	/// convenience
	PNOPairs orthonormalize_cholesky(PNOPairs& pairs) const;

	/// only the f12 part of the excited state correction of CIS(D) (namely the terms s2b and s2c)
	PairEnergies compute_cispd_f12_correction_es(const vector_real_function_3d& xcis, PairEnergies& energies) const;
	/// excited state correction of CIS(D) (namely the terms s2b and s2c)
	/// f12 part is excluded (no need to recompute after every iteration)
	PairEnergies compute_cispd_correction_es(const vector_real_function_3d& xcis, PNOPairs& pairs) const;

	/// ground state correction of CIS(D) (s4a, s4b, s4c)
	/// f12 is excluded to avoid recomputation in iterations
	PairEnergies compute_cispd_correction_gs(const vector_real_function_3d& xcis,const PNOPairs& pairs) const;

	/// The f12 part of the ground state correction of CIS(D) (s4a, s4b, s4c)
	PairEnergies compute_cispd_f12_correction_gs(const vector_real_function_3d& xcis, PairEnergies& energies) const;

	/// solve PNO-MP2 equations
	void solve_mp2(std::vector<PNOPairs>& pairs)const{
		if(pairs.empty()){
			PNOPairs mp2(MP2_PAIRTYPE,f12.acmos.size());
			solve_mp2(mp2);
			pairs.push_back(mp2);
		}else{
			solve_mp2(pairs.front());
		}
	}
	PNOPairs solve_mp2(PNOPairs& pairs) const;
	/// solve the PNO-CIS(D) equations
	std::vector<PNOPairs> solve_cispd(std::vector<PNOPairs>& pairs) const;

	/// compute the average and the maximum rank of a set of PNOs
	std::pair<size_t, size_t> get_average_rank(const std::valarray<vector_real_function_3d>& va) const;

	/// change the state of insignificant pairs to frozen
	/// This is based on:
	/// 1. External Parameters (are some pairs demanded to be frozen)
	/// 2. Energies (if the current energy is significantly below the MRA threshold there is no point in further optimizing)
	/// 3. Spatial overlapp (if the norm of the product of the originating MOs is small: i.e. if ||(phi_i)*(phi_j)||_2 is small then pair_ij is insignificant)
	PNOPairs freeze_insignificant_pairs(PNOPairs& pairs)const;

	/// Do MRA truncation on all pairs
	/// pairs are extracted into one large vector and then compressed and truncated
	/// this saves time since there are less fences
	PNOPairs truncate_pairs(PNOPairs& pairs) const;

	/// Truncate the ranks of pairs to the given maxrank
	PNOPairs truncate_pair_ranks(PNOPairs& pairs) const;

	// small helper function to keep track of all keywords which mess with the guess
	bool is_guess_from_scratch(const PairType& ct)const{
		if(param.restart()==ALL_PAIRTYPE) return false;
		else if(param.restart()==ct) return false;
		else return true;
	}

	/// save PNOs on disc
	void save_pnos(const PNOPairs& pairs) const;

	/// load PNOs from disc
	PNOPairs load_pnos(PNOPairs& pairs) const;

	/// Initialize PNO pairs
	/// If the guesstype is the default (UNKNOWN_GUESSTYPE) this will take the guess type from the parameters
	PNOPairs initialize_pairs(PNOPairs& pairs, const GuessType& inpgt = UNKNOWN_GUESSTYPE) const;

	/// compute all fluctuation potentials and store them in the pair structure
	void update_fluctuation_potentials(PNOPairs& pairs) const;
	/// Compute the MP2 fluctuation potential of a speficif pair
	PNOPairs compute_fluctuation_potential(const ElectronPairIterator& it, PNOPairs& pairs) const;
	/// Compute the CIS(D) fluctuation potential of a specific pair
 	PNOPairs compute_cispd_fluctuation_potential(const ElectronPairIterator& it, PNOPairs& pairs) const;
 	/// Iterate the pairs
 	/// The type of iteration depends on the parameters and the type of PNOPairs structure
 	/// This will call iterate_pairs_internal at some point
 	/// Depending on the parameters the pairs will be iterated adaptively or just once with one initial guess
	PNOPairs iterate_pairs(PNOPairs& pairs) const;
	/// Solve with adaptive solver descriped in the JPC paper
	PNOPairs adaptive_solver(PNOPairs& pairs) const;
	/// The actual function which interates the pairs
	/// Guess functions have to be created before
	PNOPairs iterate_pairs_internal(PNOPairs& pairs, const int maxiter, const double econv) const;

	/// Grow the rank of pairs by creating guess functions and adding them
	/// The type of guess is controlled by the parameters
	PNOPairs grow_rank(PNOPairs& pairs, std::string exop) const;

	/// convenience
	EnergyType energytype()const{return ( (param.f12() and f12.param.energytype()==PROJECTED_ENERGYTYPE) ? PROJECTED_ENERGYTYPE : HYLLERAAS_ENERGYTYPE);}

	/// solve the MP2 and CIS(D) amplitude equations
	PairEnergies t_solve(PNOPairs& pairs, const Tensor<double>& F_occ,
			const double R_convergence_threshold = 1e-6, const size_t max_niter = 100) const;

	/// Compress the PNOs -> lose all PNOs which eigenvalue is below tpno
	std::valarray<Tensor<double> > pno_compress(PNOPairs& pairs, const double tpno) const;

	/// transform the pnos and all the intermediates/potentials and matrices stored in the PNOPair structure
	PNOPairs transform_pairs(PNOPairs& pairs, const std::valarray<Tensor<double> >& U_ij) const;

	/// Print information about the PNO rank of all pairs (i.e. the number of PNO functions)
	void print_ranks(const PNOPairs& pairs)const;

	/// Update all PNOs: i.e. compute Fock-residue potential (V+2J-K, K is reused if Kpno_ij intermeidate in pairs is initialized) and apply the Green's function
	/// The fluctuation potentials have to be precomputed and stored in the pairs structure
	bool update_pno(PNOPairs& pairs, const std::valarray<Tensor<double> >& rdm_evals_ij,const Tensor<double>& F_occ) const;

	/// Convenience function to initialize all bsh operators
	std::vector<poperatorT> make_bsh_operators(World& world, const tensorT& evals) const;

	/// get the CIS potentials without Fock residue, i.e Q(2tJ - 2tK)|i> , with transformed K and J
	vector_real_function_3d compute_CIS_potentials(const vector_real_function_3d& xcis) const;


	/// the terms are expanded as follows:
	/// Q (-J1 +K1) | i(1) >  < a(2) | j(2) >
	///  +  Q | i(1) > < a(2) | -J(2) + K(2) | j(2) >
	///  +  i(1) * \int \dr2 1/|r12| a(2) j(2)
	/// the first line is zero due to orthogonality of a and j, the second line is zero due to action of Q on i
	template<typename projector>
	vector_real_function_3d compute_V_aj_i(const real_function_3d& moi,
			const real_function_3d& moj,
			const vector_real_function_3d& virtuals,
			const projector& Qpr) const;

	// not used by CIS(D) -> can still use K intermediates
	vector_real_function_3d compute_Vreg_aj_i(const size_t& i, const size_t& j,
			const vector_real_function_3d& virtuals,
			const vector_real_function_3d& Kpno) const;

	// used by CIS(D) (could also be used by MP2)
	template<typename projector>
	vector_real_function_3d compute_Vreg_aj_i(const real_function_3d& moi,
			const real_function_3d& moj,
			const vector_real_function_3d& virtuals, const projector& Qpr,
			const vector_real_function_3d& Kpno =
					vector_real_function_3d()) const;

	// the Fock residue of the regularized potential for CIS(D) (vanishes for MP2)
	// the minus sign is not included here
	// this only evalues one part (has to be called twice: -compute_Vreg_aj_i_fock_residue(Vxi,moj) - compute_Vreg_aj_i_fock_residue(moi,Vxj)
	// returns <a|Qf12|ket1ket2> = Q(<a|f12|ket2>*|ket1>
	vector_real_function_3d compute_Vreg_aj_i_fock_residue(
			const real_function_3d& ket1, const real_function_3d& ket2,
			const vector_real_function_3d& virtuals) const;
	// returns <a|(OVxQ + QOVx)f|ij>_2 = OVx(<a|f12|j>*|i>) + Q(<a|OVxf12|j>*|i>), can apply Q global since QOVx=OVx
	// use also <ta| = <a|OVx = (OVx^t|a>)^t
	vector_real_function_3d compute_Vreg_aj_i_commutator_response(
			const real_function_3d& moi, const real_function_3d& moj,
			const vector_real_function_3d& virtuals,
			const vector_real_function_3d& Vx) const;



	// compute Fluctuation Matrix without fluctuation potential (save memory in first iteration)
	Tensor<double> compute_fluctuation_matrix(const ElectronPairIterator& it, const vector_real_function_3d& pnos, const vector_real_function_3d& Kpnos_in = vector_real_function_3d()) const;

	Tensor<double> compute_cispd_fluctuation_matrix(const ElectronPairIterator& it, PNOPairs& pairs) const;

	void check_orthonormality(const vector_real_function_3d& v) const;


	/// \param[in,out] A on input = real symmetric matrix, on output = diagonal
	/// matrix of eigenvalues (from smallest to largest) \param[in,out] v on input
	/// = basis in which input A is represented, on output = the eigenbasis of A
	void canonicalize(PNOPairs& v) const;

public:
	World& world;
	PNOParameters param;  ///< calculation parameters
	Nemo nemo;
	Coulomb<double,3> J;
	ParametrizedExchange K;
	Kinetic<double, 3> T;
	Nuclear<double,3> V;
	Fock<double,3> F;
	QProjector<double, 3> Q;
	std::shared_ptr<real_convolution_3d> poisson;
	BasisFunctions basis; ///< class which holds all methods to read or create guess functions for PNO or CABS
	F12Potentials f12;
	/// convenience
	size_t nocc()const{return nemo.get_calc()->amo.size();}
	size_t nact()const{return nemo.get_calc()->amo.size()-param.freeze();}
	ElectronPairIterator pit()const{ return f12.pit();} /// convenience
	OrbitalIterator oit()const{return OrbitalIterator(nemo.get_calc()->amo.size(),param.freeze());}
	CCMessenger msg;
};



} /* namespace madness */

#endif /* PAPER_CODE_PNO_H_ */
