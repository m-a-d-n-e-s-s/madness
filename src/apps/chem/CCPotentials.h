/*
 * CCPotentials.h
 *
 *  Created on: 4 Jan 2017
 *      Author: kottmanj
 */

#ifndef SRC_APPS_CHEM_CCPOTENTIALS_H_
#define SRC_APPS_CHEM_CCPOTENTIALS_H_

#include "CCStructures.h"
#include "nemo.h"
#include "projector.h"
#include "SCFOperators.h"
#include "electronic_correlation_factor.h"

namespace madness {

/// Class which calculates all types of CC2 Potentials
  class CCPotentials{
  public:
    CCPotentials(World& world_,const Nemo& nemo,const CCParameters& param);

    virtual
    ~CCPotentials() {};

    /// forms the regularized functions from Q and Qt Ansatz for CIS(D) where tau=0 and t=mo so that Qt=Q
    void test_pair_consistency(const CCPairFunction &u,const size_t i, const size_t j, const CC_vecfunction &x)const;
    bool test_compare_pairs(const CCPair &pair1, const CCPair &pair2)const;
    void test_pairs();
    void test_singles_potential()const;
    void test();
    /// Returns the full 6D function from a CCPair structure (only recommended for debug purposes)
    real_function_6d make_6D_pair(const CCPair &pair)const;

    /// Function to load a function from disc
    /// @param[in] the function which will be loaded
    /// @param[in] name of the file in which the function was stored
    /// @return true or false depending on if the data was found on disc
    template<typename T, size_t NDIM>
    bool load_function(Function<T, NDIM>& f, const std::string name) const{
      bool exists = archive::ParallelInputArchive::exists(world,name.c_str());
      if(exists){
  	if (world.rank() == 0) print("loading function", name);
  	archive::ParallelInputArchive ar(world, name.c_str());
  	ar & f;
  	f.print_size(name);
  	f.set_thresh(FunctionDefaults<NDIM>::get_thresh());
  	f.truncate();
  	f.print_size(name);
  	return true;
      }else return false;
    }

    /// Plotting (convenience)
    void plot(const vector_real_function_3d &f, const std::string &msg) const;
    /// Plotting (convenience)
    void plot(const real_function_3d &f, const std::string &msg, const bool doprint=true)const;

    /// print size of a function
    template<size_t NDIM>
    void print_size(const Function<double,NDIM> &f, const std::string &msg, const bool print=true)const{
      if(print) f.print_size(msg);
    }

    /// get a reference to the nemo structure
    const Nemo& get_nemo()const{return nemo_;}

    /// returns a vector of all orbital energies
    std::vector<double> get_orbital_energies()const{return orbital_energies_;}

    /// returns epsilon_i + epsilon_j (needed for bsh operator of pairs)
    double get_epsilon(const size_t i, const size_t j)const{
      return orbital_energies_[i] + orbital_energies_[j];
    }

    /// returns a vector of all active mos without nuclear correlation factor (nemos)
    vector_real_function_3d get_active_mo_ket()const{
      vector_real_function_3d result;
      for(size_t i=parameters.freeze;i<mo_ket_.size();i++) result.push_back(mo_ket_(i).function);
      return result;
    }

    /// returns a vector of all active mos multiplied with the squared nuclear currelation factor: mo_bra = R^2*mo_ket
    vector_real_function_3d get_active_mo_bra()const{
      vector_real_function_3d result;
      for(size_t i=parameters.freeze;i<mo_bra_.size();i++) result.push_back(mo_bra_(i).function);
      return result;
    }

    /// get the corresponding mo bra vectors to a ket vector
    vector_real_function_3d get_mo_bra(const CC_vecfunction& ket)const{
      vector_real_function_3d result;
      for(const auto ktmp:ket.functions){
	result.push_back(mo_bra_(ktmp.first).function);
      }
      return result;
    }

    /// apply the poisson operator and use intermediates if it has such
    real_function_3d apply_poisson(const real_function_3d &f)const{
      return g12(f);
    }

    /// returns a specific mo
    CCFunction mo_ket(const size_t &i) const {
      return mo_ket_(i);
    }

    // returns constant mo_ket_
    CC_vecfunction mo_ket() const {
      return mo_ket_;
    }

    /// returns a specific mo multiplied with the squared nuclear correlation factor
    CCFunction mo_bra(const size_t &i) const {
      return mo_bra_(i);
    }

    // returns constant mo_bra_
    CC_vecfunction mo_bra() const {
      return mo_bra_;
    }

    /// make bra element: R^2|t>
    vector_real_function_3d make_bra(const CC_vecfunction& t)const{
      return make_bra(t.get_vecfunction());
    }
    vector_real_function_3d make_bra(const vector_real_function_3d &t)const{
      vector_real_function_3d bra= mul(world,nemo_.nuclear_correlation->square(),t);
      truncate(world,bra);
      return bra;
    }

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + factor*tau(i)
    /// if factor!=1 then we can not use intermediates and set the type to UNDEFINED
    CC_vecfunction make_t_intermediate(const CC_vecfunction &tau, const double factor = 1.0)const;

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + factor*tau(i)
    /// if the core is frozen the core ti will just be mo_ket_
    CC_vecfunction make_full_t_intermediate(const CC_vecfunction &tau)const;

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + tau
    /// i = tau.i
    CCFunction make_t_intermediate(const CCFunction &tau)const;

  private:
    /// Helper function to initialize the const mo_bra and ket elements adn orbital energies
    CC_vecfunction
    make_mo_bra(const Nemo& nemo) const;
    /// Helper function to initialize the const mo_bra and ket elements
    CC_vecfunction
    make_mo_ket(const Nemo& nemo) const;
    /// Transfers Eigenvalues of HF MOs from Nemo (Tensor Format) to std::vector<double> format
    std::vector<double>
    init_orbital_energies(const Nemo& nemo) const;
  public:

    // Pair functions

    /// creates pair functions for ground states
    /// the full pair-functions is: pair = u + Qt12f12|titj>
    /// this is decomposed in different types of CCPairFunction
    /// the CC_Pair2 class holds an vector of CCPairFunction
    /// we use: Qt12 = 1 - Ot1 - Ot2 + Ot12 = 1 - Ot1(1-0.5*Ot2) - Ot2(1-0.5O*t1)
    /// depending on parameters.QtAnsatz we use the t functions also in the projector
    /// functions[0] = u , the pure 6D part
    /// functions[1] = f12|titj>
    /// functions[2] = -Ot1(1-0.5*Ot2)f12|titj>
    /// functions[3] = -(1-0.5*Ot1)Ot2f12|titj>
    CCPair
    make_pair_gs(const real_function_6d& u,const CC_vecfunction& tau,const size_t i,const size_t j) const;

    /// creates pair functions for Excited states
    /// the full pair-functions is: pair = u + Qt12f12|xitj> + Qt12f12|tixj> - OxQtf12|titj> - QtOxf12|titj>
    /// the old Ansatz (not QtAnsatz) was: pair = u + Q12f12|xitj> + Q12f12|tixj>
    /// this is decomposed in different types of CCPairFunction
    /// the CC_Pair2 class holds an vector of CCPairFunction
    /// we use: Qt12 = 1 - Ot1 - Ot2 + Ot12 = 1 - Ot1(1-0.5*Ot2) - Ot2(1-0.5O*t1)
    /// functions[0] = u , the pure 6D part
    /// functions[1] = f12|xitj>
    /// functions[2] = f12|tixj>
    /// functions[3] = -Ot1(1-0.5*Ot2)f12|xitj>
    /// functions[4] = -(1-0.5*Ot1)Ot2f12|xitj>
    /// functions[5] = -Ot1(1-0.5*Ot2)f12|tixj>
    /// functions[6] = -(1-0.5*Ot1)Ot2f12|tixj>
    /// functions[7] = -OxQtf12|titj> only for QtAnsatz
    /// functions[8] = -QtOxf12|titj> only for QtAnsatz
    CCPair
    make_pair_ex(const real_function_6d& u,const CC_vecfunction& tau,const CC_vecfunction& x,const size_t i,const size_t j,const CalcType ctype) const;

    // Energy operations

    /// Compute pair correlation energies of MP2 or CC2 Ground State
    // Off diagonal pair energies are multiplied with 2.0 to acount for their permuted partners
    /// @param[in] The Pair_function
    /// @param[in] The Singles (for MP2 give empty function) for the energy contribution over disconnected doubles
    /// @param[out] 2*<ij|g|u> - <ji|g|u> , where i and j are determined by u (see CC_Pair class)
    double
    compute_pair_correlation_energy(const CCPair& u,const CC_vecfunction& singles=CC_vecfunction(PARTICLE)) const;

    /// Compute CC2 correlation energy
    /// @param[in] The Pair_function
    /// @param[out] \sum_{ij} 2*<ij|g|u> - <ji|g|u> + 2*<ij|g|\tau_i\tau_j> - <ji|g|\tau_i\tau_j> , where i and j are determined by u (see CC_Pair class)
    /// since we do not compute all pairs (symmetry reasons) the off diagonal pair energies are conted twice
    /// the cc2 pair functions are dependent on the doubles (see CC_Pair structure, and make_pair function) so make shure they are updated
    double
    compute_cc2_correlation_energy(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;


    double
    compute_kinetic_energy(const vector_real_function_3d& xbra,const vector_real_function_3d& xket) const;

    /// returns \f$  <x|T|x> + <x|V|x>  \f$
    double
    compute_cis_expectation_value(const CC_vecfunction& x,const vector_real_function_3d& V,const bool print=true) const;

    /// Something like a pair energy for CIS(D)/LRCC2 to estimate energy convergence
    /// calculates the response part of s2b and s2c which are independent of the mp2 amplitudes
    double
    compute_excited_pair_energy(const CCPair& d,const CC_vecfunction& x) const;

    /// Compute the CIS(D) Energy Correction to CIS
    double
    compute_cispd_energy(const CC_vecfunction& x,const Pairs<CCPair> mp2,const Pairs<CCPair> cispd) const;

    /// Compute the CC2 Excitation Energy
    double
    compute_cc2_excitation_energy(const CC_vecfunction& stau,const CC_vecfunction& sx,const Pairs<CCPair> dtau,const Pairs<CCPair> dx) const;

    // 6D operations

    /// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
    real_function_6d
    fock_residue_6d(const CCPair& u) const;

    /// Function evaluates the consant part of the ground state for MP2
    /// @param[out]The result is \f$ Q12(G(Q12(Vreg|titj>))) \f$
    /// @param[in] ti, first particle -> should be HOLE state
    /// @param[in] tj, second particle -> should be HOLE state
    /// @param[in] Gscreen pointer to bsh operator (in order to screen), has to be in modified NS form
    real_function_6d
    make_constant_part_mp2(const CCFunction& ti,const CCFunction& tj,const real_convolution_6d* Gscreen=NULL) const;

    /// Function evaluates the consant part of the ground state for CC2
    /// @param[out]The result is \f$ Q12(G(Q12((Vreg+V_{coupling})|titj>))) \f$ with \f$ |t_k> = |tau_k> + |k> \f$
    /// @param[in] u, The Pair function
    /// @param[in] tau, The CC2 singles of the coupling potential -> should be PARTICLE state
    /// @param[in] Gscreen pointer to bsh operator (in order to screen), has to be in modified NS form
    /// for the coulomb coupling potential we use:
    /// \f$ V_{CC} = -QOtau -OtauQ + OtauOtau
    ///            = - Otau(Qt(1/2)) - Qt(1/2)(Otau) \f$
    /// where t(1/2) = |i> + 1/2|tau_i> , t(1/2) = th
    real_function_6d
    make_constant_part_cc2_gs(const CCPair& u,const CC_vecfunction& tau,const real_convolution_6d* Gscreen=NULL) const;

    /// Function evaluates the consant part of the ground state for CC2 if the Qt Ansatz is used
    /// @param[out]The result is \f$ Q12(G(Qt12((Vreg+V_{coupling})|titj> + [F,Qt]f12|titj>))) \f$ with \f$ |t_k> = |tau_k> + |k>  and Qt = Q - \sum_k |tau_k><k| \f$
    /// @param[in] u, The Pair function
    /// @param[in] tau, The CC2 singles of the coupling potential -> should be PARTICLE state
    /// @param[in] Gscreen pointer to bsh operator (in order to screen), has to be in modified NS form
    /// for the coulomb coupling potential we use:
    /// \f$ V_{CC} = -QOtau -OtauQ + OtauOtau
    ///            = - Otau(Qt(1/2)) - Qt(1/2)(Otau) \f$
    /// where t(1/2) = |i> + 1/2|tau_i> , t(1/2) = th
    real_function_6d
    make_constant_part_cc2_Qt_gs(const CCPair& u,const CC_vecfunction& tau,const real_convolution_6d* Gscreen=NULL) const;

    /// Function evaluates the consant part of the Excited state for CIS(D) if the Q Ansatz is used
    real_function_6d
    make_constant_part_cispd(const CCPair& u,const CC_vecfunction& x,const real_convolution_6d* Gscreen=NULL) const;
    /// Function evaluates the consant part of the Excited state for CIS(D) if the Qt Ansatz is used
    real_function_6d
    make_constant_part_cispd_Qt(const CCPair& u,const CC_vecfunction& x,const real_convolution_6d* Gscreen=NULL) const;
    /// Function evaluates the consant part of the Excited state for CC2 if the Q Ansatz is used
    real_function_6d
    make_constant_part_cc2_ex(const CCPair& u,const CC_vecfunction& tau,const CC_vecfunction& x,const real_convolution_6d* Gscreen=NULL);
    /// Function evaluates the consant part of the Excited state for CC2 if the Qt Ansatz is used
    real_function_6d
    make_constant_part_cc2_Qt_ex(const CCPair& u,const CC_vecfunction& tau,const CC_vecfunction& x,const real_convolution_6d* Gscreen=NULL);
    /// Apply the Regularization potential
    /// \f$ V_{reg} = [ U_e - [K,f12] + f12(F12-eij) ]|titj> \f$
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    /// @param[out] the regularization potential (unprojected), see equation above
    real_function_6d
    apply_Vreg(const CCFunction& ti,const CCFunction& tj,const real_convolution_6d* Gscreen=NULL) const;

    /// evaluates: \f$ (F(1)-ei)|ti> (x) |tj> + |ti> (x) (F(2)-ej)|tj> \f$ with the help of the singles potential
    /// singles equation is: (F-ei)|ti> = - V(ti)
    /// response singles equation: (F-ei-omega)|xi> = - V(xi)
    /// response:  \f$ (F12-ei-ej-omega)|xitj> = (F1 - ei - omega)|xi> (x) |tj> + |xi> (x) (F2-ej)|tj> \f$
    /// so in both cases the result will be: |V(ti),tj> + |ti,V(tj)>
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    real_function_6d
    apply_reduced_F(const CCFunction& ti,const CCFunction& tj,const real_convolution_6d* Gscreen=NULL) const;



    /// Apply Ue on a tensor product of two 3d functions: Ue(1,2) |x(1)y(2)> (will be either |ij> or |\tau_i\tau_j> or mixed forms)
    /// The Transformed electronic regularization potential (Kutzelnigg) is R_{12}^{-1} U_e R_{12} with R_{12} = R_1*R_2
    /// It is represented as: R_{12}^{-1} U_e R_{12} = U_e + R^-1[Ue,R]
    /// where R^-1[Ue,R] = R^-1 [[T,f],R] (see: Regularizing the molecular potential in electronic structure calculations. II. Many-body
    /// methods, F.A.Bischoff)
    /// The double commutator can be evaluated as follows:  R^-1[[T,f],R] = -Ue_{local}(1,2)*(Un_{local}(1) - Un_{local}(2))
    /// @param[in] x the 3D function for particle 1
    /// @param[in] y the 3D function for particle 2
    /// @param[in] The BSH operator to screen: Has to be in NS form, Gscreen->modified == true
    /// @return  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
    real_function_6d
    apply_transformed_Ue(const CCFunction& x,const CCFunction& y,const real_convolution_6d* Gscreen=NULL) const;

    /// Apply Ue on a tensor product of two 3d functions: Ue(1,2) |x(1)y(2)> (will be either |ij> or |\tau_i\tau_j> or mixed forms)
    /// The Transformed electronic regularization potential (Kutzelnigg) is R_{12}^{-1} U_e R_{12} with R_{12} = R_1*R_2
    /// It is represented as: R_{12}^{-1} U_e R_{12} = U_e + R^-1[Ue,R]
    /// where R^-1[Ue,R] = R^-1 [[T,f],R] (see: Regularizing the molecular potential in electronic structure calculations. II. Many-body
    /// methods, F.A.Bischoff)
    /// The double commutator can be evaluated as follows:  R^-1[[T,f],R] = -Ue_{local}(1,2)*(Un_{local}(1) - Un_{local}(2))
    /// @param[in] x the 3D function for particle 1
    /// @param[in] y the 3D function for particle 2
    /// @param[in] Gscreen pointer to the BSH operator in order to screen, has to be in modified NS form, Gscreen->modified==true
    /// the f12K|xy> part will be screened with the BSH while the Kf12|xy> can not be screened with the BSH operator but maybe with the coulomb
    /// @return  R^-1U_eR|x,y> the transformed electronic smoothing potential applied on |x,y> :
    real_function_6d
    apply_exchange_commutator(const CCFunction& x,const CCFunction& y,const real_convolution_6d* Gscreen=NULL) const;
    /// This applies the exchange commutator, see apply_exchange_commutator function for information
    real_function_6d
    apply_exchange_commutator1(const CCFunction& x,const CCFunction& y,const real_convolution_6d* Gscreen=NULL) const;

    /// Helper Function which performs the operation \f$ <xy|g12f12|ab> \f$
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function a,
    /// @param[in] function b,
    double
    make_xy_gf_ab(const CCFunction& x,const CCFunction& y,const CCFunction& a,const CCFunction& b) const;

    double make_xy_ff_ab(const CCFunction &x, const CCFunction &y, const CCFunction &a, const CCFunction &b)const{
      error("xy_ff_ab not yet implemented");
      return 0.0;
    }

    /// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
    /// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
    real_function_3d
    apply_gf(const real_function_3d& f) const;

    /// returns <xy|op|u>
    /// loops over every entry in the vector and accumulates results
    /// helper function for CIS(D) energy
    double
    make_xy_op_u(const CCFunction& x,const CCFunction& y,const CCConvolutionOperator& op,const std::vector<CCPairFunction>& u) const;

    /// returns <xy|u> for a vector of CCPairFunction
    /// the result is accumulated for every vercotr
    /// helper functions for CIS(D) energy
    double
    make_xy_u(const CCFunction& x,const CCFunction& y,const std::vector<CCPairFunction>& u) const;

    /// Functions which operate with the CCPairFunction structure
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] CCPairFunction u,
    double
    make_xy_op_u(const CCFunction& x,const CCFunction& y,const CCConvolutionOperator& op,const CCPairFunction& u) const;

    /// Helper Function which returns
    /// @return <xy|op|ab>
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function a,
    /// @param[in] function b,
    double
    make_xy_op_ab(const CCFunction& x,const CCFunction& y,const CCConvolutionOperator& op,const CCFunction& a,const CCFunction& b) const;

    /// get the correct pair function as vector of CCPairFunction functions
    /// @param[in] The pair functions
    /// @param[out] The demanded pair function as vector of CCPairFunction functions (includes regularization tails)
    std::vector<CCPairFunction>
    get_pair_function(const Pairs<CCPair>& pairs,const size_t i,const size_t j) const;


    /// returns <a|g12|u>_2
    real_function_3d
    apply_s2b_operation(const CCFunction& bra,const CCPairFunction& u,const size_t particle) const;

    /// dummy to avoid confusion and for convenience
    real_function_6d swap_particles(const real_function_6d &f)const{
      return madness::swap_particles<double>(f);
    }

    /// swap the particles of the CCPairFunction and return a new vector of swapped functions
    std::vector<CCPairFunction> swap_particles(const std::vector<CCPairFunction> & f) const{
      std::vector<CCPairFunction> swapped;
      for(size_t i=0;i<f.size();i++) swapped.push_back(f[i].swap_particles());
      return swapped;
    }

    /// Calculate the Overlap between two CCPairFunction <f1|f2>
    /// @param[in] 6D function 1
    /// @param[in] 6D function 2
    double
    overlap(const CCPairFunction& f1,const CCPairFunction& f2) const;

    /// Computes the squared norm of the pair function <x|x>
    /// including all regularization tails
    /// means: if  x = u + f12|mn> + |ab> --> <x|x> = <u|u> + <u|f12|mn> + <u|ab> + ....
    double
    overlap(const CCPair& x) const;

    // Projectors

    /// Apply a projector to a CC_vecfunction
    /// result: \f$ P(tau)f[i] = \sum_k tau[k]*<k|f[i]>   \f$
    /// ket state of the projector can vary
    /// bra state are always MOs
    vector_real_function_3d
    apply_projector(const CC_vecfunction& f,const CC_vecfunction& ket_) const;

    /// Apply Qt projector on 6D function
    real_function_6d
    apply_Q12t(const real_function_6d& f,const CC_vecfunction& t) const;


    /// Apply the Q projector to a CC_vecfunction
    /// result: \f$ 1-c*P(tau)f[i] = 1 - c*\sum_k tau[k]*<k|f[i]> \f$
    /// ket state of the projector can vary
    /// bra state are always MOs
    /// the factor c is usually 1 or 1/2
    vector_real_function_3d
    apply_Qt(const CC_vecfunction& f,const CC_vecfunction& ket_,const double c=1.0) const;

    /// Apply the Qt projector on a CCPairFunction
    /// works in principle like apply_Ot
    CCPairFunction
    apply_Qt(const CCPairFunction& f,const CC_vecfunction& t,const size_t particle,const double c=1.0) const;

    /// Apply Ot projector on decomposed or op_decomposed 6D function
    /// The function does not work with type==pure right now (not needed)
    /// for CCPairFunction type==decomposed_ : f = \sum_k |c_k,d_k> where c and d are stored as vectors
    /// then the result is \sum_k |a_k,b_k> with a,b stored as vectors and
    /// \f$ a_k = t_k \f$
    /// \f$ b_k = Ot(d_k) or vise versa for particle==2
    /// for CCPairFunction type == op_decomposd the function si f=op|xy> and we have for particle==1
    /// \f$ a_k = t_k \f$
    /// \f$ b_k = <mo_k|op|x>*y \f$
    CCPairFunction
    apply_Ot(const CCPairFunction& f,const CC_vecfunction& t,const size_t particle) const;

    /// Apply the Greens Operator to a CCPairFunction
    /// For CCPairFunction only type pure and type decomposed is supported
    /// for the op_decomposed type a pure function can be constructed (not needed therefore not implemented yet)
    real_function_6d
    apply_G(const CCPairFunction& u,const real_convolution_6d& G) const;

    /// Apply BSH Operator and count time
    real_function_6d apply_G(const real_function_6d &f, const real_convolution_6d &G)const{
      CCTimer time(world,"Apply Greens Operator");
      const real_function_6d result = G(f);
      time.info(true,result.norm2());
      return result;
    }


    // Coupled Cluster Singles potentials which work with CCPairFunction pairs
    // notation follows Shavit & Bartlett Many-Body Methods ...

    /// Calculates the CC2 singles potential for the ground state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    vector_real_function_3d
    get_CC2_singles_potential_gs(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;

    /// Calculates the CCS/CIS singles potential for the excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    /// the expectation value is calculated and updated
    vector_real_function_3d
    get_CCS_potential_ex(CC_vecfunction& x,const bool print=false) const;

    /// Calculates the CIS potential (no expectation value or storing)
    vector_real_function_3d
    get_CIS_potential(const CC_vecfunction& x) const;



    /// Calculates the CC2 singles potential for the Excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    vector_real_function_3d
    get_CC2_singles_potential_ex(const CC_vecfunction& gs_singles,const Pairs<CCPair>& gs_doubles,CC_vecfunction& ex_singles,const Pairs<CCPair>& response_doubles) const;

    /// Calculates the CC2 singles potential for the Excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    vector_real_function_3d
    get_ADC2_singles_potential(const Pairs<CCPair>& gs_doubles,CC_vecfunction& ex_singles,const Pairs<CCPair>& response_doubles) const;

    /// The potential manager for the ground state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] Bra vectorfunction to contract the potential
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    double
    potential_energy_gs(const CC_vecfunction& bra,const CC_vecfunction& singles,const Pairs<CCPair>& doubles,const PotentialType& name) const;

    /// The potential manager for the ground state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    vector_real_function_3d
    potential_singles_gs(const CC_vecfunction& singles,const Pairs<CCPair>& doubles,const PotentialType& name) const;

    /// The integra manager for the excited state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] The bra element
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Singles of the Excited State
    /// @param[in] Doubles of the Excited State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    double
    potential_energy_ex(const CC_vecfunction& bra,const CC_vecfunction& singles_gs,const Pairs<CCPair>& doubles_gs,const CC_vecfunction& singles_ex,const Pairs<CCPair>& doubles_ex,
			const PotentialType& name) const;

    /// The potential manager for the excited state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Singles of the Excited State
    /// @param[in] Doubles of the Excited State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    vector_real_function_3d
    potential_singles_ex(const CC_vecfunction& singles_gs,const Pairs<CCPair>& doubles_gs,const CC_vecfunction& singles_ex,const Pairs<CCPair>& doubles_ex,const PotentialType& name) const;

    /// The Fock operator is partitioned into F = T + Vn + R
    /// the fock residue R= 2J-K+Un for closed shell is computed here
    /// J_i = \sum_k <k|r12|k> |tau_i>
    /// K_i = \sum_k <k|r12|tau_i> |k>
    vector_real_function_3d
    fock_residue_closed_shell(const CC_vecfunction& singles) const;

    /// the K operator runs over ALL orbitals (also the frozen ones)
    real_function_3d
    K(const CCFunction& f) const;

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// if i==j in uij then the symmetry will be exploited
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d
    K(const real_function_6d& u,const bool symmetric) const;

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
    /// 1. X(3,2) = bra_k(3)*u(3,2)
    /// 2. Y(1,2) = \int X(3,2) g13 d3
    /// 3. result = Y(1,2)*ket_k(1)
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d
    apply_K(const real_function_6d& u,const size_t& particle) const;

    /// Apply the Exchange operator on a tensor product multiplied with f12
    /// !!! Prefactor of (-1) is not inclued in K here !!!!
    real_function_6d
    apply_Kf(const CCFunction& x,const CCFunction& y) const;

    /// Apply fK on a tensor product of two 3D functions
    /// fK|xy> = fK_1|xy> + fK_2|xy>
    /// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    /// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    /// @param[in] BSH operator to screen, has to be in modified NS form, Gscreen->modified()==true;
    real_function_6d
    apply_fK(const CCFunction& x,const CCFunction& y,const real_convolution_6d* Gscreen=NULL) const;

    /// Creates a 6D function with the correlation factor and two given CCFunctions
    real_function_6d
    make_f_xy(const CCFunction& x,const CCFunction& y,const real_convolution_6d* Gscreen=NULL) const;

    /// the ccs potential without terms from the fock operator
    /// returns: \f$ (1-|\tau_k><k|)(2 <k|g|tau_k> |t_i> - <k|g|t_i> |\tau_k>)  \f$
    vector_real_function_3d
    ccs_potential_gs(const CC_vecfunction& tau) const;

    /// the ccs potential for the response equations
    /// returns \f$ d\dtau Q(tau)(unprojected_ccs(t,tau) = Q(tau)(unprojected(x,tau)+unprojected(t,x)) - O(x)(unprojected(t,tau))\f$
    vector_real_function_3d
    ccs_potential_ex(const CC_vecfunction& singles_gs,const CC_vecfunction& singles_ex) const;


    /// unprojected ccs potential
    /// returns 2kgtk|ti> - kgti|tk>
    /// the ccs potential: ti = ti and tk = tauk
    vector_real_function_3d
    ccs_unprojected(const CC_vecfunction& ti,const CC_vecfunction& tk) const;


    real_function_3d
    make_density(const CC_vecfunction& x) const;

    vector_real_function_3d
    cis_potential_ex(const CC_vecfunction& x) const;

    // integrals from singles potentials

    /// <xi|T|ti> + <xi|Vn|ti> + 2.0*<xi,k|g|ti,k> - <xi,k|g|k,ti>
    double
    x_s3a(const CC_vecfunction& x,const CC_vecfunction& t) const;

    /// -<xi|e_i|ti>
    double
    x_s3b(const CC_vecfunction& x,const CC_vecfunction& t) const;

    /// 2<xik|g|itk> - <xik|g|tki>
    double
    x_s3c(const CC_vecfunction& x,const CC_vecfunction& t) const;
    /// 2<xik|g|t1it2k> - <xik|g|t2kt1i>
    double
    x_s5b(const CC_vecfunction& x,const CC_vecfunction& t1,const CC_vecfunction& t2) const;
    /// -(2<lk|g|it1k> - <lk|g|t1ki>)* <xi|t2l>
    double
    x_s5c(const CC_vecfunction& x,const CC_vecfunction& t1,const CC_vecfunction& t2) const;
    /// -(2<lk|g|t3i,t1k> - <lk|g|t1k,t3i>)* <xi|t2l>
    double
    x_s6(const CC_vecfunction& x,const CC_vecfunction& t1,const CC_vecfunction& t2,const CC_vecfunction& t3) const;
    /// 2.0 <xik|g|uik>- <kxi|g|uik>
    double
    x_s2b(const CC_vecfunction& x,const Pairs<CCPair>& u) const;
    /// -( 2.0 <xi,k|(l_g_i)(2)|ukl> - <k,xi,(l_g_i)(1)|ukl> )
    double
    x_s2c(const CC_vecfunction& x,const Pairs<CCPair>& u) const;
    /// -( <lk|g|uik> - <kl|g|uik> )*<xi|tl>
    double
    x_s4a(const CC_vecfunction& x,const CC_vecfunction& t,const Pairs<CCPair>& u) const;
    /// -( 2.0 <xi,k|(l_g_ti)(2)|ukl> - <k,xi,(l_g_ti)(1)|ukl> )
    double
    x_s4b(const CC_vecfunction& x,const CC_vecfunction& t,const Pairs<CCPair>& u) const;

    /// 4.0 <xil|(k_g_tk)|uil> - 2.0 <lxi|(k_g_tk)|uil> - 2.0 <xik|(l_g_tk)|uil> + <kxi|(l_g_tk)|uil>
    double
    x_s4c(const CC_vecfunction& x,const CC_vecfunction& t,const Pairs<CCPair>& u) const;

    // result: \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 )
    // singles are not needed explicitly but to determine if it is response or ground state
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[out] \f$ \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    /// if the s2b potential has already been calculated it will be loaded from the intermediate_potentials structure
    vector_real_function_3d
    s2b(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;

    // result: -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1)
    // singles are not needed explicitly but to determine if it is response or ground state
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[out] \f$ -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1) \f$
    /// Q-Projector is not applied, sign is correct
    vector_real_function_3d
    s2c(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;

    /// the S4a potential can be calcualted from the S2b potential
    /// result is \f$ s4a_i = - <l|s2b_i>*|tau_l> \f$
    vector_real_function_3d
    s4a_from_s2b(const vector_real_function_3d& s2b,const CC_vecfunction& singles) const;

    // result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[out] \f$ -( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied, sign is correct
    vector_real_function_3d
    s4b(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;


    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param[out] \f$ ( 4<l|kgtauk|uil>_2 - 2<l|kgtauk|uil>_1 - 2<k|lgtauk|uil>_2 + <k|lgtauk|uil>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    vector_real_function_3d
    s4c(const CC_vecfunction& singles,const Pairs<CCPair>& doubles) const;

    // update the intermediates
    void update_intermediates(const CC_vecfunction &t){
      g12.update_elements(mo_bra_,t);
      f12.update_elements(mo_bra_,t);
    }

    /// clear stored potentials
    /// if a response function is given only the response potentials are cleared (the GS potentials dont change anymore)
    void clear_potentials(const CC_vecfunction &t)const{

      if(t.type==RESPONSE){
	output("Clearing Response Singles-Potentials");
	get_potentials.clear_response();
      }
      else{
	output("Clearing all stored Singles-Potentials");
	get_potentials.clear_all();
      }
    }

  private:
    // member variables
    /// MPI World
    World & world;
    /// the nemo structure which holds the nuclear correlation factor and all other functions from the reference calculation
    const Nemo& nemo_;
    /// CC_Parameters: All parameters (see class description)
    const CCParameters& parameters;
    /// ket functions of the reference determinant
    CC_vecfunction mo_ket_;
    /// bra function of the reference determinant: bra = R2*ket, with R^2 beiing the squared nuclear correlation factor
    CC_vecfunction mo_bra_;
    /// orbital energies
    std::vector<double> orbital_energies_;
    /// the coulomb operator with all intermediates
    CCConvolutionOperator g12;
    /// the f12 operator with all intermediates
    CCConvolutionOperator f12;
    /// the correlation factor, holds necessary regularized potentials
    CorrelationFactor corrfac;
    /// Manager for stored intermediate potentials which are s2c, s2b and the whole singles potentials without fock-residue for GS and EX state
    mutable CCIntermediatePotentials get_potentials;
  public:
    /// Messenger structure for formated output and to store warnings
    CCMessenger output;

  };

} /* namespace madness */

#endif /* SRC_APPS_CHEM_CCPOTENTIALS_H_ */
