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
#include "electronic_correlation_factor.h"

namespace madness {

  class TDHF{
  public:
    TDHF(World & world, const CCParameters& param,const Nemo &nemo);
    virtual
    ~TDHF();
    /// Initialize the CIS functions
    /// @param[in\out] on input the already obtained guess functions (or empty vector), on output new guess functions are added
    void initialize(std::vector<CC_vecfunction> &start)const;
    /// Solve the CIS equations
    /// @param[in/out] CC_vecfunction
    /// on input the guess functions (if empty or not enough the a guess will be generated)
    /// on output the solution
    void solve_cis()const{
      std::vector<CC_vecfunction> ccs;
      for(size_t k=0;k<parameters.excitations_.size();k++){
  	CC_vecfunction tmp;
  	const bool found= initialize_singles(tmp,RESPONSE,parameters.excitations_[k]);
  	if(found) ccs.push_back(tmp);
      }
      solve_cis(ccs);
    }
    void solve_cis(std::vector<CC_vecfunction>& start) const;
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
    /// Make the CIS Guess
    /// the type of guess is  ontrolled over the tda_guess keyword in the CCParameters class (CC_Structures.h)
    std::vector<CC_vecfunction> make_guess()const;
    /// Guess only takes the Homos into account
    /// Guess functions are created from the application of the excitation operators (see tda_guess keyword)
    /// applied to the Homo orbital and to the Homo-1 ... Homo-N orbital with N=tda_guess_orbitals parameter
    /// With this we can get more irreps without the need for large polynomials
    std::vector<CC_vecfunction> make_homo_guess()const;
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
    /// Make the TDHF potential
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
      set_thresh(world, tmp, parameters.thresh_3D);
      truncate(world,tmp);
      reconstruct(world,tmp);
      CC_vecfunction mo_bra(tmp, HOLE);
      return mo_bra;
    }

    CC_vecfunction make_mo_ket(const Nemo&nemo) const {
      vecfuncT tmp = nemo.get_calc()->amo;
      set_thresh(world, tmp, parameters.thresh_3D);
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

    /// The MPI Communicator
    World& world;
    /// The Parameters for the Calculations
    const CCParameters& parameters;
    /// The Nemo structure (convenience)
    const Nemo& nemo;
    /// Operator Structure which can handle intermediates (use for exchange with GS orbitals)
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
