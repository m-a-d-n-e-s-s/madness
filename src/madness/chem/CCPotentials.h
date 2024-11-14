/*
 * CCPotentials.h
 *
 *  Created on: 4 Jan 2017
 *      Author: kottmanj
 */

#ifndef SRC_APPS_CHEM_CCPOTENTIALS_H_
#define SRC_APPS_CHEM_CCPOTENTIALS_H_

#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/nemo.h>
#include<madness/chem/projector.h>
#include<madness/chem/SCFOperators.h>
#include<madness/chem/electronic_correlation_factor.h>

namespace madness {
/// Class which calculates all types of CC2 Potentials
class CCPotentials {
public:
    CCPotentials(World& world_, const std::shared_ptr<Nemo> nemo, const CCParameters& param);

    void reset_nemo(const std::shared_ptr<Nemo> nemo) {
        nemo_ = nemo;
        mo_ket_ = (make_mo_ket(*nemo));
        mo_bra_ = (make_mo_bra(*nemo));
        orbital_energies_ = init_orbital_energies(*nemo);
    };

    Info update_info(const CCParameters& parameters, const std::shared_ptr<Nemo> nemo) const {
        Info info;
        info.mo_bra = mo_bra().get_vecfunction();
        info.mo_ket = mo_ket().get_vecfunction();
        info.molecular_coordinates = nemo->get_calc()->molecule.get_all_coords_vec();
        info.parameters = parameters;
        info.R_square = nemo->R_square;
        info.R = nemo->R;
        info.U1 = nemo->ncf->U1vec();
        info.U2 = nemo->ncf->U2();
        info.intermediate_potentials = get_potentials;
        info.orbital_energies = orbital_energies_;
        info.fock=nemo->compute_fock_matrix(nemo->get_calc()->amo, nemo->get_calc()->aocc);
        return info;
    }

    virtual
    ~CCPotentials() {};

    /// forms the regularized functions from Q and Qt Ansatz for CIS(D) where tau=0 and t=mo so that Qt=Q
    void test_pair_consistency(const CCPairFunction<double, 6>& u, const size_t i, const size_t j,
                               const CC_vecfunction& x) const;

    bool test_compare_pairs(const CCPair& pair1, const CCPair& pair2) const;

    void test_pairs();

    void test_singles_potential(Info& info) const;

    void test();

    /// Returns the full 6D function from a CCPair structure (only recommended for debug purposes)
    real_function_6d make_6D_pair(const CCPair& pair) const;

    /// Function to load a function from disc
    /// @param[in] f the function which will be loaded
    /// @param[in] name of the file in which the function was stored
    /// @param do_print
    /// @return true or false depending on if the data was found on disc
    template <typename T, size_t NDIM>
    bool load_function(Function<T, NDIM>& f, const std::string name, bool do_print) const {
        bool exists = archive::ParallelInputArchive<
            archive::BinaryFstreamInputArchive>::exists(world, name.c_str());
        if (exists) {
            if ((world.rank() == 0) and do_print) print("loading function", name);
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str());
            ar & f;
            if (do_print) f.print_size(name);
            if (f.is_compressed()) {
                if (world.rank()==0 and do_print) print("function is compressed -- reconstructing");
                f.change_tree_state(reconstructed);
                if (do_print) f.print_size(name+" reconstructed");
                save(f, name);
            }
            f.set_thresh(FunctionDefaults<NDIM>::get_thresh());
            f.truncate();
            f.print_size(name);
            return true;
        } else {
            if ((world.rank()==0) and do_print) print("could not find function",name);
        }
        return false;
    }

    /// Plotting (convenience)
    void plot(const vector_real_function_3d& f, const std::string& msg) const;

    /// Plotting (convenience)
    void plot(const real_function_3d& f, const std::string& msg, const bool doprint = true) const;

    /// print size of a function
    template <size_t NDIM>
    void print_size(const Function<double, NDIM>& f, const std::string& msg, const bool print = true) const {
        if (print) f.print_size(msg);
    }

    /// returns a vector of all orbital energies
    std::vector<double> get_orbital_energies() const { return orbital_energies_; }

    /// returns epsilon_i + epsilon_j (needed for bsh operator of pairs)
    double get_epsilon(const size_t i, const size_t j) const {
        return orbital_energies_[i] + orbital_energies_[j];
    }

    /// returns epsilon_i + epsilon_j (needed for bsh operator of pairs)
    static double get_epsilon(const size_t i, const size_t j, const Info& info) {
        return info.orbital_energies[i] + info.orbital_energies[j];
    }

    /// returns a vector of all active mos without nuclear correlation factor (nemos)
    vector_real_function_3d get_active_mo_ket() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_ket_.size(); i++) result.push_back(mo_ket_(i).function);
        return result;
    }

    /// returns a vector of all active mos multiplied with the squared nuclear currelation factor: mo_bra = R^2*mo_ket
    vector_real_function_3d get_active_mo_bra() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_bra_.size(); i++) result.push_back(mo_bra_(i).function);
        return result;
    }

    /// get the corresponding mo bra vectors to a ket vector
    vector_real_function_3d get_mo_bra(const CC_vecfunction& ket) const {
        vector_real_function_3d result;
        for (const auto& ktmp : ket.functions) {
            result.push_back(mo_bra_(ktmp.first).function);
        }
        return result;
    }

    /// returns a specific mo
    CCFunction<double, 3> mo_ket(const size_t& i) const {
        return mo_ket_(i);
    }

    // returns constant mo_ket_
    CC_vecfunction mo_ket() const {
        return mo_ket_;
    }

    /// returns a specific mo multiplied with the squared nuclear correlation factor
    CCFunction<double, 3> mo_bra(const size_t& i) const {
        return mo_bra_(i);
    }

    // returns constant mo_bra_
    CC_vecfunction mo_bra() const {
        return mo_bra_;
    }

    /// make bra element: R^2|t>
    vector_real_function_3d make_bra(const CC_vecfunction& t) const {
        return make_bra(t.get_vecfunction());
    }

    vector_real_function_3d make_bra(const vector_real_function_3d& t) const {
        vector_real_function_3d bra = mul(world, nemo_->ncf->square(), t);
        truncate(world, bra);
        return bra;
    }

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + factor*tau(i)
    /// if factor!=1 then we can not use intermediates and set the type to UNDEFINED
    CC_vecfunction make_t_intermediate(const CC_vecfunction& tau, const CCParameters& parameters) const;

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + factor*tau(i)
    /// if the core is frozen the core ti will just be mo_ket_
    CC_vecfunction make_full_t_intermediate(const CC_vecfunction& tau) const;

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + tau(i)
    /// if the core is frozen the core ti will just be mo_ket_
    static CC_vecfunction make_full_t_intermediate(const CC_vecfunction& tau, const Info& info);

    /// makes the t intermediates

    /// t_i = mo_ket_(i) + tau(i)
    /// skip frozen orbitals
    static CC_vecfunction make_active_t_intermediate(const CC_vecfunction& tau, const Info& info);

    /// makes the t intermediates
    /// t_i = mo_ket_(i) + tau
    /// i = tau.i
    CCFunction<double, 3> make_t_intermediate(const CCFunction<double, 3>& tau) const;

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
    /// return the regularized MP2 ansatz: |\tau_ij> = |u_ij> + Q12 f12 |ij>
    static CCPair make_pair_mp2(World& world, const real_function_6d& u, const size_t i, const size_t j,
                                const Info& info, bool compute_Q12_f12_ij);

    /// return the regularized CC2 ansatz: |\tau_ij> = |u_ij> + Q12t f12 |t_i t_j>
    static CCPair make_pair_cc2(World& world, const real_function_6d& u,
                                const CC_vecfunction& gs_singles, const size_t i, const size_t j, const Info& info, const bool compute_Q12_f12_ij);

    /// return the regularized CC2 ansatz: |x_ij> = |u_ij> + Q12t f12 |t_i t_j> + ?????
    static CCPair make_pair_lrcc2(World& world, const CalcType& ctype,
                                  const real_function_6d& u, const CC_vecfunction& gs_singles,
                                  const CC_vecfunction& ex_singles, const size_t i, const size_t j, const Info& info,
                                  const bool compute_Q12_f12);

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
    make_pair_gs(const real_function_6d& u, const CC_vecfunction& tau, const size_t i, const size_t j) const;

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
    make_pair_ex(const real_function_6d& u, const CC_vecfunction& tau, const CC_vecfunction& x, const size_t i,
                 const size_t j, const CalcType ctype) const;

    // Energy operations

    /// Compute pair correlation energies of MP2 or CC2 Ground State
    // Off diagonal pair energies are multiplied with 2.0 to acount for their permuted partners
    /// @param[in] u the Pair_function
    /// @param[in] singles the Singles (for MP2 give empty function) for the energy contribution over disconnected doubles
    /// @param[out] 2*<ij|g|u> - <ji|g|u> , where i and j are determined by u (see CC_Pair class)
    static double
    compute_pair_correlation_energy(World& world,
                                    const CCPair& u,
                                    const CC_vecfunction& singles,
                                    const Info& info);

    /// Compute CC2 correlation energy
    /// @param[in] The Pair_function
    /// @param[out] \sum_{ij} 2*<ij|g|u> - <ji|g|u> + 2*<ij|g|\tau_i\tau_j> - <ji|g|\tau_i\tau_j> , where i and j are determined by u (see CC_Pair class)
    /// since we do not compute all pairs (symmetry reasons) the off diagonal pair energies are conted twice
    /// the cc2 pair functions are dependent on the doubles (see CC_Pair structure, and make_pair function) so make shure they are updated
    /// @param world
    /// @param info
    static double
    compute_cc2_correlation_energy(World& world, const CC_vecfunction& singles, const Pairs<CCPair>& doubles,
                                   const Info& info, const std::string msg="");


    static double
    compute_kinetic_energy(World& world, const vector_real_function_3d& xbra, const vector_real_function_3d& xket);

    /// compute the expectation value excitation energy using the CIS/CCS/CC2 singles
    static double
    compute_cis_expectation_value(World& world, const CC_vecfunction& x,
                                  const vector_real_function_3d& V, const bool print, const Info& info);

    /// Something like a pair energy for CIS(D)/LRCC2 to estimate energy convergence
    /// calculates the response part of s2b and s2c which are independent of the mp2 amplitudes
    static double
    compute_excited_pair_energy(World& world, const CCPair& d, const CC_vecfunction& x, const Info& info);

    /// Compute the CIS(D) Energy Correction to CIS
    double
    compute_cispd_energy(const CC_vecfunction& x, const Pairs<CCPair> mp2, const Pairs<CCPair> cispd) const;

    /// Compute the CC2 Excitation Energy
    double
    compute_cc2_excitation_energy(const CC_vecfunction& stau, const CC_vecfunction& sx, const Pairs<CCPair> dtau,
                                  const Pairs<CCPair> dx) const;

    // 6D operations

    /// The 6D Fock residue on the cusp free pair function u_{ij}(1,2) is: (2J - Kn - Un)|u_{ij}>
    real_function_6d
    fock_residue_6d(const CCPair& u) const;

    /// Static function for the 6D Fock residue for use in macrotask
    static madness::real_function_6d
    fock_residue_6d_macrotask(World& world, const CCPair& u, const CCParameters& parameters,
                              const std::vector<madness::Vector<double, 3>>& all_coords_vec,
                              const std::vector<real_function_3d>& mo_ket,
                              const std::vector<real_function_3d>& mo_bra,
                              const std::vector<real_function_3d>& U1,
                              const real_function_3d& U2);

    /// Static version of make_constant_part_mp2 to be called from macrotask.
    static madness::real_function_6d
    make_constant_part_mp2_macrotask(World& world, const CCPair& pair, const std::vector<real_function_3d>& mo_ket,
                                     const std::vector<real_function_3d>& mo_bra,
                                     const CCParameters& parameters, const real_function_3d& Rsquare,
                                     const std::vector<real_function_3d>& U1,
                                     const std::vector<std::string> argument);

    /// Compute the constant part of MP2, CC2 or LR-CC2
    ///
    /// depending on pair.calc_type different terms are included in the constant part.
    /// @param[in] pair         the (empty) pair function, determines the terms in the constant part, contains some bookkeeping information (bsh_eps, i, j)
    /// @param[in] gs_singles   the ground-state singles for CC2 (used for the T1-transformed SO projector), may be left empty for MP2
    /// @param[in] ex_singles   the excited-state singles for CC2 (used for the T1-transformed SO projector), may be left empty for MP2 and GS-CC2
    /// @param[in] info         the Info object, containing the some basic quantities (MOs, parameters, etc)
    /// @return            the constant part of the MP2, CC2 or LR-CC2: G(Q12(g~|titj>))
    static madness::real_function_6d
    make_constant_part_macrotask(World& world, const CCPair& pair,
                                 const CC_vecfunction& gs_singles, const CC_vecfunction& ex_singles,
                                 const Info& info);


    /// Static function to iterate the mp2 pairs from macrotask
    static madness::real_function_6d
    update_pair_mp2_macrotask(World& world, const CCPair& pair, const CCParameters& parameters,
                              const std::vector<madness::Vector<double, 3>>& all_coords_vec,
                              const std::vector<real_function_3d>& mo_ket,
                              const std::vector<real_function_3d>& mo_bra,
                              const std::vector<real_function_3d>& U1,
                              const real_function_3d& U2, const real_function_6d& mp2_coupling);


    /// iterate a pair for MP2, CC2, LRCC2 on constant singles
    static CCPair iterate_pair_macrotask(World& world,
                                         const CCPair& pair, const CC_vecfunction& gs_singles,
                                         const CC_vecfunction& ex_singles,
                                         const real_function_6d& coupling, const Info& info, const long maxiter);


    /// Function evaluates the consant part of the Excited state for CIS(D) if the Q Ansatz is used
    real_function_6d
    make_constant_part_cispd(const CCPair& u, const CC_vecfunction& x,
                             const real_convolution_6d* Gscreen = NULL) const;

    /// Function evaluates the consant part of the Excited state for CIS(D) if the Qt Ansatz is used
    real_function_6d
    make_constant_part_cispd_Qt(const CCPair& u, const CC_vecfunction& x,
                                const real_convolution_6d* Gscreen = NULL) const;

    /// Apply the Regularization potential
    /// \f$ V_{reg} = [ U_e - [K,f12] + f12(F12-eij) ]|titj> \f$
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    /// @param[out] the regularization potential (unprojected), see equation above
    real_function_6d
    apply_Vreg(const CCFunction<double, 3>& ti, const CCFunction<double, 3>& tj,
               const real_convolution_6d* Gscreen = NULL) const;

    /// Apply the Regularization potential
    /// \f$ V_{reg} = [ U_e - [K,f12] + f12(F12-eij) + [F,Qt] ]|titj> \f$
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    /// @param[out] the regularization potential (unprojected), see equation above
    std::vector<CCPairFunction<double, 6>>
    static apply_Vreg(World& world, const CCFunction<double, 3>& ti, const CCFunction<double, 3>& tj,
                      const CC_vecfunction& gs_singles, const CC_vecfunction& ex_singles,
                      const Info& info, const std::vector<std::string>& argument, const double bsh_eps);

    /// Static version of apply_Vreg to be used from a macrotask. Will eventually replace former.
    madness::real_function_6d
    static
    apply_Vreg_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                         const std::vector<real_function_3d>& mo_bra,
                         const CCParameters& parameters, const real_function_3d& Rsquare,
                         const std::vector<real_function_3d>& U1, const size_t& i, const size_t& j,
                         const FuncType& x_type, const FuncType& y_type,
                         const std::vector<std::string> argument,
                         const real_convolution_6d* Gscreen = NULL);

    /// evaluates: \f$ (F(1)-ei)|ti> (x) |tj> + |ti> (x) (F(2)-ej)|tj> \f$ with the help of the singles potential
    /// singles equation is: (F-ei)|ti> = - V(ti)
    /// response singles equation: (F-ei-omega)|xi> = - V(xi)
    /// response:  \f$ (F12-ei-ej-omega)|xitj> = (F1 - ei - omega)|xi> (x) |tj> + |xi> (x) (F2-ej)|tj> \f$
    /// so in both cases the result will be: |V(ti),tj> + |ti,V(tj)>
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    real_function_6d
    apply_reduced_F1(const CCFunction<double, 3>& ti, const CCFunction<double, 3>& tj,
                     const real_convolution_6d* Gscreen = NULL) const;

    /// evaluates: \f$ (F(1)-ei)|ti> (x) |tj> + |ti> (x) (F(2)-ej)|tj> \f$ with the help of the singles potential
    /// singles equation is: (F-ei)|ti> = - V(ti)
    /// response singles equation: (F-ei-omega)|xi> = - V(xi)
    /// response:  \f$ (F12-ei-ej-omega)|xitj> = (F1 - ei - omega)|xi> (x) |tj> + |xi> (x) (F2-ej)|tj> \f$
    /// so in both cases the result will be: |V(ti),tj> + |ti,V(tj)>
    /// @param[in] ti, first function in the ket, for MP2 it is the Orbital, for CC2 the relaxed Orbital t_i=\phi_i + \tau_i
    /// @param[in] tj, second function in the ket ...
    /// @param[in] pointer to bsh operator (in order to screen)
    real_function_6d
    static apply_reduced_F(World& world, const CCFunction<double, 3>& ti, const CCFunction<double, 3>& tj,
                           const Info& info, const real_convolution_6d* Gscreen = NULL);

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
    apply_transformed_Ue(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                         const real_convolution_6d* Gscreen = NULL) const;

    /// Static version of apply_transformed_Ue for the use in a macrotask.
    /// Will eventually replace the former.
    real_function_6d
    static apply_transformed_Ue_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                          const CCParameters& parameters, const real_function_3d& Rsquare,
                                          const std::vector<real_function_3d>& U1, const size_t& i, const size_t& j,
                                          const FuncType& x_type, const FuncType& y_type,
                                          const real_convolution_6d* Gscreen = NULL);

    real_function_6d
    static apply_Ue(World& world, const CCFunction<double, 3>& phi_i, const CCFunction<double, 3>& phi_j,
                    const Info& info, const real_convolution_6d* Gscreen);


    static real_function_6d
    apply_KffK(World& world, const CCFunction<double, 3>& phi_i, const CCFunction<double, 3>& phi_j,
               const Info& info, const real_convolution_6d* Gscreen);
    static CCPairFunction<double, 6>
    apply_commutator_F_Qt_f12(World& world, const CCFunction<double, 3>& phi_i, const CCFunction<double, 3>& phi_j,
                              const CC_vecfunction& gs_singles, const CC_vecfunction& ex_singles,
                              const Info& info, const real_convolution_6d* Gscreen);

    static CCPairFunction<double, 6>
    apply_commutator_F_dQt_f12(World& world, const CCFunction<double, 3>& phi_i, const CCFunction<double, 3>& phi_j,
                               const CC_vecfunction& gs_singles, const CC_vecfunction& ex_singles,
                               const Info& info, const real_convolution_6d* Gscreen);

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
    apply_exchange_commutator(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                              const real_convolution_6d* Gscreen = NULL) const;

    real_function_6d
    static apply_exchange_commutator_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                                               const std::vector<real_function_3d>& mo_bra,
                                               const real_function_3d& Rsquare,
                                               const size_t& i, const size_t& j, const CCParameters& parameters,
                                               const FuncType& x_type, const FuncType& y_type,
                                               const real_convolution_6d* Gscreen = NULL);

    /// This applies the exchange commutator, see apply_exchange_commutator function for information
    real_function_6d
    apply_exchange_commutator1(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                               const real_convolution_6d* Gscreen = NULL) const;

    /// Helper Function which performs the operation \f$ <xy|g12f12|ab> \f$
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function a,
    /// @param[in] function b,
    double
    make_xy_gf_ab(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y, const CCFunction<double, 3>& a,
                  const CCFunction<double, 3>& b) const;

    double make_xy_ff_ab(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                         const CCFunction<double, 3>& a, const CCFunction<double, 3>& b) const {
        error("xy_ff_ab not yet implemented");
        return 0.0;
    }

    /// apply the operator gf = 1/(2\gamma)*(Coulomb - 4\pi*BSH_\gamma)
    /// works only if f = (1-exp(-\gamma*r12))/(2\gamma)
    static real_function_3d
    apply_gf(World& world, const real_function_3d& f, const Info& info);

    /// returns <xy|op|u>
    /// loops over every entry in the vector and accumulates results
    /// helper function for CIS(D) energy
    static double
    make_xy_op_u(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                 const CCConvolutionOperator<double, 3>& op,
                 const std::vector<CCPairFunction<double, 6>>& u);

    /// returns <xy|u> for a vector of CCPairFunction
    /// the result is accumulated for every vercotr
    /// helper functions for CIS(D) energy
    static double
    make_xy_u(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
              const std::vector<CCPairFunction<double, 6>>& u);

    /// Functions which operate with the CCPairFunction structure
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] CCPairFunction u,
    static double
    make_xy_op_u(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                 const CCConvolutionOperator<double, 3>& op,
                 const CCPairFunction<double, 6>& u);

    /// Helper Function which returns
    /// @return <xy|op|ab>
    /// @param[in] function x, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function y, if nuclear correlation is used make sure this is the correct bra function
    /// @param[in] function a,
    /// @param[in] function b,
    double
    make_xy_op_ab(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                  const CCConvolutionOperator<double, 3>& op, const CCFunction<double, 3>& a,
                  const CCFunction<double, 3>& b) const;

    /// get the correct pair function as vector of CCPairFunction functions
    /// @param[in] The pair functions
    /// @param[out] The demanded pair function as vector of CCPairFunction functions (includes regularization tails)
    static std::vector<CCPairFunction<double, 6>>
    get_pair_function(const Pairs<CCPair>& pairs, const size_t i, const size_t j);


    /// returns <a|g12|u>_2
    static real_function_3d
    apply_s2b_operation(World& world, const CCFunction<double, 3>& bra, const CCPairFunction<double, 6>& u,
                        const size_t particle, const Info& info);

    /// dummy to avoid confusion and for convenience
    real_function_6d swap_particles(const real_function_6d& f) const {
        return madness::swap_particles<double>(f);
    }

    /// swap the particles of the CCPairFunction and return a new vector of swapped functions
    static std::vector<CCPairFunction<double, 6>> swap_particles(const std::vector<CCPairFunction<double, 6>>& f) {
        std::vector<CCPairFunction<double, 6>> swapped;
        for (size_t i = 0; i < f.size(); i++) swapped.push_back(f[i].swap_particles());
        return swapped;
    }

    /// Calculate the Overlap between two CCPairFunction <f1|f2>
    /// @param[in] 6D function 1
    /// @param[in] 6D function 2
    double
    overlap(const CCPairFunction<double, 6>& f1, const CCPairFunction<double, 6>& f2) const {
        return inner(f1, f2, nemo_->ncf->square());
    };

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
    apply_projector(const CC_vecfunction& f, const CC_vecfunction& ket_) const;

    /// Apply Qt projector on 6D function
    real_function_6d
    apply_Q12t(const real_function_6d& f, const CC_vecfunction& t) const;


    /// Apply the Q projector to a CC_vecfunction
    /// result: \f$ 1-c*P(tau)f[i] = 1 - c*\sum_k tau[k]*<k|f[i]> \f$
    /// ket state of the projector can vary
    /// bra state are always MOs
    /// the factor c is usually 1 or 1/2
    vector_real_function_3d
    apply_Qt(const CC_vecfunction& f, const CC_vecfunction& ket_, const double c = 1.0) const;

    /// Apply the Qt projector on a CCPairFunction
    /// works in principle like apply_Ot
    CCPairFunction<double, 6>
    apply_Qt(const CCPairFunction<double, 6>& f, const CC_vecfunction& t, const size_t particle,
             const double c = 1.0) const;

    /// Apply Ot projector on decomposed or op_decomposed 6D function
    /// The function does not work with type==pure right now (not needed)
    /// for CCPairFunction type==decomposed_ : f = \sum_k |c_k,d_k> where c and d are stored as vectors
    /// then the result is \sum_k |a_k,b_k> with a,b stored as vectors and
    /// \f$ a_k = t_k \f$
    /// \f$ b_k = Ot(d_k) or vise versa for particle==2
    /// for CCPairFunction type == op_decomposd the function si f=op|xy> and we have for particle==1
    /// \f$ a_k = t_k \f$
    /// \f$ b_k = <mo_k|op|x>*y \f$
    CCPairFunction<double, 6>
    apply_Ot(const CCPairFunction<double, 6>& f, const CC_vecfunction& t, const size_t particle) const;

    /// Apply the Greens Operator to a CCPairFunction
    /// For CCPairFunction only type pure and type decomposed is supported
    /// for the op_decomposed type a pure function can be constructed (not needed therefore not implemented yet)
    real_function_6d
    apply_G(const CCPairFunction<double, 6>& u, const real_convolution_6d& G) const;

    /// Apply BSH Operator and count time
    real_function_6d apply_G(const real_function_6d& f, const real_convolution_6d& G) const {
        CCTimer time(world, "Apply Greens Operator");
        const real_function_6d result = G(f);
        time.info(true, result.norm2());
        return result;
    }


    // Coupled Cluster Singles potentials which work with CCPairFunction pairs
    // notation follows Shavit & Bartlett Many-Body Methods ...

    /// Calculates the CC2 singles potential for the ground state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    static vector_real_function_3d
    get_CC2_singles_potential_gs(World& world, const CC_vecfunction& singles, const Pairs<CCPair>& doubles,
                                 Info& info);

    /// Calculates the CCS/CIS singles potential for the excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    /// the expectation value is calculated and updated
    static vector_real_function_3d
    get_CCS_potential_ex(World& world, const CC_vecfunction& x, const bool print, Info& info);

    /// Calculates the CC2 singles potential for the Excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    static vector_real_function_3d
    get_CC2_singles_potential_ex(World& world, const CC_vecfunction& gs_singles,
                                 const Pairs<CCPair>& gs_doubles, const CC_vecfunction& ex_singles,
                                 const Pairs<CCPair>& response_doubles, Info& info);

    /// Calculates the CC2 singles potential for the Excited state: result = Fock_residue + V
    /// the V part is stored in the intermediate_potentials structure
    static vector_real_function_3d
    get_ADC2_singles_potential(World& world, const Pairs<CCPair>& gs_doubles,
                               CC_vecfunction& ex_singles, const Pairs<CCPair>& response_doubles, Info& info);

    /// The potential manager for the ground state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] Bra vectorfunction to contract the potential
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    /// @param world
    double
    potential_energy_gs(World& world, const CC_vecfunction& bra, const CC_vecfunction& singles,
                        const Pairs<CCPair>& doubles, const PotentialType& name) const;

    /// The potential manager for the ground state potential
    /// CC2 singles potential parts of the ground state
    /// General function which evaluates a CC_singles potential
    /// @param[in] result_index corresponds to indices of external lines in the diagram
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    /// @param world
    static std::tuple<madness::vector_real_function_3d, madness::vector_real_function_3d>
    potential_singles_gs(World& world, const std::vector<int>& result_index, const CC_vecfunction& singles,
                         const Pairs<CCPair>& doubles,
                         const PotentialType& name, const Info& info);

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
    /// @param world
    double
    potential_energy_ex(World& world, const CC_vecfunction& bra, const CC_vecfunction& singles_gs,
                        const Pairs<CCPair>& doubles_gs, const CC_vecfunction& singles_ex,
                        const Pairs<CCPair>& doubles_ex, const PotentialType& name) const;

    /// The potential manager for the excited state potential
    /// CC2 singles potential parts of the ground state
    /// Genereal function which evaluates a CC_singles potential
    /// @param[in] result_index dummy vector for partitioning the macrotasks, must have the length of the singles vector
    /// @param[in] Singles of the Ground State
    /// @param[in] Doubles of the Ground State
    /// @param[in] Singles of the Excited State
    /// @param[in] Doubles of the Excited State
    /// @param[in] Name of the potential
    /// @param[out] the potential (without Q application)
    /// @param world
    static std::tuple<madness::vector_real_function_3d, madness::vector_real_function_3d>
    potential_singles_ex(World& world, const std::vector<int> result_index, const CC_vecfunction& singles_gs,
                         const Pairs<CCPair>& doubles_gs, const CC_vecfunction& singles_ex,
                         const Pairs<CCPair>& doubles_ex, const PotentialType& name, const Info& info);

    /// The Fock operator is partitioned into F = T + Vn + R
    /// the fock residue R= 2J-K+Un for closed shell is computed here
    /// J_i = \sum_k <k|r12|k> |tau_i>
    /// K_i = \sum_k <k|r12|tau_i> |k>
    static vector_real_function_3d
    fock_residue_closed_shell(World& world, const CC_vecfunction& singles, const Info& info);

    /// the K operator runs over ALL orbitals (also the frozen ones)
    static real_function_3d
    K(World& world, const CCFunction<double, 3>& f, const Info& info);

    /// static version of k above for access from macrotask. will eventually replace former.
    real_function_3d
    static K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                       const std::vector<real_function_3d>& mo_bra, const real_function_3d& f,
                       const CCParameters& parameters);

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// if i==j in uij then the symmetry will be exploited
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d
    K(const real_function_6d& u, const bool symmetric) const;

    /// static version of k above for access from macrotask. will eventually replace former.
    real_function_6d
    static K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                       const std::vector<real_function_3d>& mo_bra,
                       const real_function_6d& u, const bool symmetric, const CCParameters& parameters);

    /// Exchange Operator on Pair function: -(K(1)+K(2))u(1,2)
    /// K(1)u(1,2) = \sum_k <k(3)|g13|u(3,2)> |k(1)>
    /// 1. X(3,2) = bra_k(3)*u(3,2)
    /// 2. Y(1,2) = \int X(3,2) g13 d3
    /// 3. result = Y(1,2)*ket_k(1)
    /// !!!!Prefactor (-1) is not included here!!!!
    real_function_6d
    apply_K(const real_function_6d& u, const size_t& particle) const;

    /// Static version of apply_K above for access from macrotask. Will eventually replace former.
    real_function_6d
    static apply_K_macrotask(World& world, const std::vector<real_function_3d>& mo_ket,
                             const std::vector<real_function_3d>& mo_bra,
                             const real_function_6d& u, const size_t& particle, const CCParameters& parameters);

    /// Apply the Exchange operator on a tensor product multiplied with f12
    /// !!! Prefactor of (-1) is not inclued in K here !!!!
    real_function_6d
    apply_Kf(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y) const;

    /// Apply fK on a tensor product of two 3D functions
    /// fK|xy> = fK_1|xy> + fK_2|xy>
    /// @param[in] x, the first 3D function in |xy>, structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    /// @param[in] y, the second 3D function in |xy>  structure holds index i and type (HOLE, PARTICLE, MIXED, UNDEFINED)
    /// @param[in] BSH operator to screen, has to be in modified NS form, Gscreen->modified()==true;
    real_function_6d
    apply_fK(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
             const real_convolution_6d* Gscreen = NULL) const;

    /// Creates a 6D function with the correlation factor and two given CCFunctions
    real_function_6d
    make_f_xy(const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
              const real_convolution_6d* Gscreen = NULL) const;

    /// Creates a 6D function with the correlation factor and two given CCFunctions
    real_function_6d
    static make_f_xy(World& world, const CCFunction<double, 3>& x, const CCFunction<double, 3>& y,
                     const Info& info, const real_convolution_6d* Gscreen = NULL);

    real_function_6d
    static make_f_xy_macrotask(World& world, const real_function_3d& x_ket, const real_function_3d& y_ket,
                               const real_function_3d& x_bra, const real_function_3d& y_bra,
                               const size_t& i, const size_t& j, const CCParameters& parameters,
                               const FuncType& x_type, const FuncType& y_type,
                               const real_convolution_6d* Gscreen = NULL);

    /// unprojected ccs (S3c) potential
    /// returns 2kgtk|ti> - kgti|tk>
    /// the ccs potential: ti = ti and tk = tauk
    static vector_real_function_3d
    ccs_unprojected(World& world, const CC_vecfunction& ti, const CC_vecfunction& tk, const Info& info);

    /// return RMS norm and max norm of residuals
    template <typename T, std::size_t NDIM>
    static std::pair<double, double> residual_stats(const std::vector<Function<T, NDIM>>& residual) {
        if (residual.size() == 0) return std::make_pair(0.0, 0.0);
        World& world = residual.front().world();
        auto errors = norm2s(world, residual);
        double rnorm = 0.0, maxrnorm = 0.0;
        for (double& e : errors) {
            maxrnorm = std::max(maxrnorm, e);
            rnorm += e * e;
        }
        rnorm = sqrt(rnorm / errors.size());
        return std::make_pair(rnorm, maxrnorm);
    }

    static void print_convergence(const std::string name, const double rmsresidual, const double maxresidual,
                                  const double energy_diff, const int iteration) {
        const std::size_t bufsize = 255;
        char msg[bufsize];
        std::snprintf(msg, bufsize,
                      "convergence of %s in iteration %2d at time %8.1fs: rms/max residual, energy change %.1e %.1e %.1e",
                      name.c_str(), iteration, wall_time(), rmsresidual, maxresidual,energy_diff);
        print(msg);
    }

    // integrals from singles potentials

    /// <xi|T|ti> + <xi|Vn|ti> + 2.0*<xi,k|g|ti,k> - <xi,k|g|k,ti>
    double
    x_s3a(const CC_vecfunction& x, const CC_vecfunction& t) const;

    /// -<xi|e_i|ti>
    double
    x_s3b(const CC_vecfunction& x, const CC_vecfunction& t) const;

    /// 2<xik|g|itk> - <xik|g|tki>
    double
    x_s3c(const CC_vecfunction& x, const CC_vecfunction& t) const;

    /// 2<xik|g|t1it2k> - <xik|g|t2kt1i>
    double
    x_s5b(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2) const;

    /// -(2<lk|g|it1k> - <lk|g|t1ki>)* <xi|t2l>
    double
    x_s5c(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2) const;

    /// -(2<lk|g|t3i,t1k> - <lk|g|t1k,t3i>)* <xi|t2l>
    double
    x_s6(const CC_vecfunction& x, const CC_vecfunction& t1, const CC_vecfunction& t2,
         const CC_vecfunction& t3) const;

    /// 2.0 <xik|g|uik>- <kxi|g|uik>
    double
    x_s2b(const CC_vecfunction& x, const Pairs<CCPair>& u) const;

    /// -( 2.0 <xi,k|(l_g_i)(2)|ukl> - <k,xi,(l_g_i)(1)|ukl> )
    double
    x_s2c(const CC_vecfunction& x, const Pairs<CCPair>& u) const;

    /// -( <lk|g|uik> - <kl|g|uik> )*<xi|tl>
    double
    x_s4a(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const;

    /// -( 2.0 <xi,k|(l_g_ti)(2)|ukl> - <k,xi,(l_g_ti)(1)|ukl> )
    double
    x_s4b(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const;

    /// 4.0 <xil|(k_g_tk)|uil> - 2.0 <lxi|(k_g_tk)|uil> - 2.0 <xik|(l_g_tk)|uil> + <kxi|(l_g_tk)|uil>
    double
    x_s4c(const CC_vecfunction& x, const CC_vecfunction& t, const Pairs<CCPair>& u) const;

    /// result: \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 )

    /// singles are not needed explicitly but to determine if it is response or ground state
    /// function return a vector of size 2*singles, the first half being the s2b potential, the second half
    /// an intermediate that will be stored by the caller in the info structure
    ///@param world
    ///@param external_indices
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param builder
    ///@param info
    ///@param[out] \f$ \sum_k( 2<k|g|uik>_2 - <k|g|uik>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    /// if the s2b potential has already been calculated it will be loaded from the intermediate_potentials structure
    static std::tuple<madness::vector_real_function_3d, madness::vector_real_function_3d>
    s2b(World& world, std::vector<int> external_indices, const CC_vecfunction& singles, const Pairs<CCPair>& doubles, const
        CCPairBuilder& builder, const Info& info);

    /// result: -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1)

    /// singles are not needed explicitly but to determine if it is response or ground state
    /// function return a vector of size 2*singles, the first half being the s2c potential, the second half
    /// an intermediate that will be stored by the caller in the info structure
    ///@param world
    ///@param external_index
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param builder
    ///@param info
    ///@param[out] \f$ -\sum_k( <l|kgi|ukl>_2 - <l|kgi|ukl>_1) \f$
    /// Q-Projector is not applied, sign is correct
    static std::tuple<vector<Function<double, 3>>, vector<Function<double, 3>>>
    s2c(World& world, std::vector<int> external_index, const CC_vecfunction& singles, const Pairs<CCPair>& doubles, const
        CCPairBuilder& builder, const Info& info);

    /// the S4a potential can be calcualted from the S2b potential
    /// result is \f$ s4a_i = - <l|s2b_i>*|tau_l> \f$
    vector_real_function_3d
    s4a_from_s2b(const vector_real_function_3d& s2b, const CC_vecfunction& singles) const;

    // result: -\sum_k( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui>
    ///@param world
    ///@param external_index
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param info
    ///@param[out] \f$ -( <l|kgtaui|ukl>_2 - <l|kgtaui|ukl>_1) | kgtaui = <k|g|taui> | taui=singles_i \f$
    /// Q-Projector is not applied, sign is correct
    static vector_real_function_3d
    s4b(World& world, std::vector<int> external_index, const CC_vecfunction& singles, const Pairs<CCPair>& doubles, const Info& info);


    ///@param world
    ///@param external_index
    ///@param[in] singles:CC_vecfunction fof type response or particle (depending on this the correct intermediates will be used) the functions themselves are not needed
    ///@param[in] doubles:Pairs of CC_Pairs (GS or Response)
    ///@param info
    ///@param[out] \f$ ( 4<l|kgtauk|uil>_2 - 2<l|kgtauk|uil>_1 - 2<k|lgtauk|uil>_2 + <k|lgtauk|uil>_1 ) \f$
    /// Q-Projector is not applied, sign is correct
    static vector_real_function_3d
    s4c(World& world, std::vector<int> external_index, const CC_vecfunction& singles, const Pairs<CCPair>& doubles, const Info& info);

    // update the intermediates
    void update_intermediates(const CC_vecfunction& t) {
        g12->update_elements(mo_bra_, t);
        //        g12.sanity();
        f12->update_elements(mo_bra_, t);
        //        f12.sanity();
    }

    /// clear stored potentials
    /// if a response function is given only the response potentials are cleared (the GS potentials dont change anymore)
    void clear_potentials(const CC_vecfunction& t) const {
        if (t.type == RESPONSE) {
            output("Clearing Response Singles-Potentials");
            get_potentials.clear_response();
        }
        else {
            output("Clearing all stored Singles-Potentials");
            get_potentials.clear_all();
        }
    }

public:
    // member variables
    /// MPI World
    World& world;
    /// the nemo structure which holds the nuclear correlation factor and all other functions from the reference calculation
    std::shared_ptr<Nemo> nemo_;
    /// CC_Parameters: All parameters (see class description)
    const CCParameters& parameters;
    /// ket functions of the reference determinant
    CC_vecfunction mo_ket_;
    /// bra function of the reference determinant: bra = R2*ket, with R^2 beiing the squared nuclear correlation factor
    CC_vecfunction mo_bra_;
    /// orbital energies
    std::vector<double> orbital_energies_;
    /// the coulomb operator with all intermediates
public:
    std::shared_ptr<CCConvolutionOperator<double, 3>> g12;
    /// the f12 operator with all intermediates
    std::shared_ptr<CCConvolutionOperator<double, 3>> f12;
    /// the correlation factor, holds necessary regularized potentials
    CorrelationFactor corrfac;
    /// Manager for stored intermediate potentials which are s2c, s2b and the whole singles potentials without fock-residue for GS and EX state
    mutable CCIntermediatePotentials get_potentials;
    /// POD for basis and intermediates
    Info info;

public:
    /// Messenger structure for formated output and to store warnings
    CCMessenger output;
};
} /* namespace madness */

#endif /* SRC_APPS_CHEM_CCPOTENTIALS_H_ */
