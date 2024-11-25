/*
 * CC2.h
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */

#ifndef CC2_H_
#define CC2_H_

#include<madness/chem/projector.h>
#include<madness/chem/SCF.h>
#include<madness/chem/nemo.h>
#include<madness/chem/CCPotentials.h>
#include<madness/chem/mp3.h>
#include <madness/mra/operator.h>
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include<madness/chem/TDHF.h>
#include <madness/mra/nonlinsol.h>

#include "BSHApply.h"

namespace madness {

class CC2 : public OptimizationTargetInterface, public QCPropertyInterface {
public:

    CC2(World& world_, const commandlineparser& parser, const std::shared_ptr<Nemo> nemo_)
            : world(world_),
              parameters(world_,parser),
              nemo(nemo_),
              CCOPS(world, nemo, parameters),
              output(CCOPS.output) {

        output.section("CC2 Class has been initialized with the following parameters");
        // set the threshholds
        // Set Protocoll
        output("Set Protocol 3D");
        nemo_->get_calc()->set_protocol<3>(world, parameters.thresh_3D());
        output("Set Protocol 6D");
        nemo_->get_calc()->set_protocol<6>(world, parameters.thresh_6D());

        FunctionDefaults<3>::set_thresh(parameters.thresh_3D());
        FunctionDefaults<6>::set_thresh(parameters.thresh_6D());
        // Make sure that k is the same in 3d and 6d functions
        FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
        // by default SCF sets the truncate_mode to 1
        FunctionDefaults<3>::set_truncate_mode(3);
        FunctionDefaults<6>::set_truncate_mode(3);
        // set initial level to 3 (important to avoid undersampling issues)
        // the "real initial level" is then 2 since the initial level gets substracted by one if refine is true (which is the default)
        FunctionDefaults<6>::set_initial_level(2);
        FunctionDefaults<3>::set_initial_level(2);
        FunctionDefaults<6>::set_special_level(
                FunctionDefaults<6>::set_length_scale(parameters.dmin(), FunctionDefaults<6>::get_k()));
        FunctionDefaults<3>::set_special_level(
                FunctionDefaults<3>::set_length_scale(parameters.dmin(), FunctionDefaults<3>::get_k()));
        parameters.information(world);
        parameters.sanity_check(world);

        tdhf.reset(new TDHF(world,parser,nemo));

    }

    virtual ~CC2() {}


    double value() {
        return value(nemo->molecule().get_all_coords());
    }

    double value(const Tensor<double>& x) {
        solve();
        return 0.0;
    }

    void output_calc_info_schema(const std::string model, const double& energy) const;


    std::string name() const {return "CC2";};

    static void help() {
        print_header2("help page for CC2");
        print("The CC2 code computes correlated ground and excited state energies:\n");
        print(" - MP2 ground state");
        print(" - CC2 ground and excited states");
        print(" - ADC(2) and CIS(D) excited states\n");
        print("You need a SCF reference calculation from the nemo program. If there no such calculation can");
        print("be found CC2 will perform its own. If excited states are requested also a CIS calculation is ");
        print("necessary.\n");
        print("Note that for correlated calculations the k parameter must be chosen small, typically k=5 or k=6 ");
        print("because the curse of dimensions make higher k extremely expensive\n");
        print("You can print all available calculation parameters by running\n");
        print("cc2 --print_parameters\n");
        print("You can perform a simple MP2 calculation by running\n");
        print("cc2 --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.\n");

    }

    static void print_parameters() {
        CCParameters param;
        print("\ndefault parameters for the cc2 program are\n");
        param.print("cc2","end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    virtual bool selftest() {
        return true;
    };
    void
    plot(const real_function_3d& f, const std::string& msg = "unspecified function") const {
        plot_plane(world, f, msg);
        output("Plotted " + msg);
    }

    /// The World
    World& world;
    /// Structure holds all the parameters used in the CC2 calculation
    CCParameters parameters;
    /// The SCF Calculation
    std::shared_ptr<Nemo> nemo;
    /// The excited state cis calculation
    std::shared_ptr<TDHF> tdhf;
    /// The CC Operator Class
    CCPotentials CCOPS;
    /// Formated Output (same as used in CC2Potentials structure)
    CCMessenger& output;
    /// map Pair struct to vector
    PairVectorMap triangular_map;

    /// solve the CC2 ground state equations, returns the correlation energy
    void solve();


    std::vector<CC_vecfunction>
    solve_ccs() const;

    double compute_mp3(const Pairs<CCPair>& mp2pairs) const {
        MP3 mp3(CCOPS);
        double mp3_contribution=mp3.mp3_energy_contribution_macrotask_driver(mp2pairs);
        return mp3_contribution;
    }

    double
    solve_cc2(CC_vecfunction& tau, Pairs<CCPair>& u, Info& info) const;

    /// solve the excited state LR-CC2 equations for a given excitation

    /// @param[in] gs_doubles: the ground state doubles
    /// @param[in] gs_singles: the ground state singles
    /// @param[in] cis: the CIS singles
    /// @param[in] excitation: the excitation number
    /// @return a tuple with the excited state doubles, the excited state singles and the excitation energy
    std::tuple<Pairs<CCPair>, CC_vecfunction, double>
    solve_lrcc2(Pairs<CCPair>& gs_doubles, const CC_vecfunction& gs_singles, const CC_vecfunction& cis,
        const std::size_t excitation, Info& info) const;

    double
    solve_cispd(Pairs<CCPair>& doubles, const Pairs<CCPair>& mp2_pairs, const CC_vecfunction& cis_singles);

    /// convencience function to iterate the CC2 ground state singles,
    /// makes the right call on the iterate_singles functions
    static bool
    iterate_cc2_singles(World& world, CC_vecfunction& singles, Pairs<CCPair>& doubles, Info& info) {
        // CCOPS.clear_potentials(singles);
        info.intermediate_potentials.clear_all();
        Pairs<CCPair> empty;
        return iterate_singles(world, singles, CC_vecfunction(RESPONSE), doubles,
            empty, CT_CC2, info.parameters.iter_max_3D(), info);
    }

    bool
    iterate_adc2_singles(Pairs<CCPair>& mp2, CC_vecfunction& singles, Pairs<CCPair>& x, Info& info) {
        MADNESS_ASSERT(singles.type == RESPONSE);
        // CCOPS.clear_potentials(singles);
        info.intermediate_potentials.clear_response();
        return iterate_singles(world, singles, CC_vecfunction(UNDEFINED), mp2, x, CT_ADC2, parameters.iter_max_3D(), info);
    }

    static bool
    iterate_lrcc2_singles(World& world, const CC_vecfunction& cc2_s, Pairs<CCPair>& cc2_d, CC_vecfunction& lrcc2_s, Pairs<CCPair> lrcc2_d, Info& info) {
        MADNESS_ASSERT(cc2_s.type == PARTICLE);
        MADNESS_ASSERT(lrcc2_s.type == RESPONSE);
        info.intermediate_potentials.clear_response();
        // CCOPS.clear_potentials(lrcc2_s);
        return iterate_singles(world, lrcc2_s, cc2_s, cc2_d, lrcc2_d,
            CT_LRCC2, info.parameters.iter_max_3D(), info);
    }

    /// convencience function to iterate the CCS Response singles,
    /// makes the right call on the iterate_singles functions
    bool
    iterate_ccs_singles(CC_vecfunction& x, Info& info) const {
        Pairs<CCPair> empty;
        // CCOPS.clear_potentials(x);
        info.intermediate_potentials.clear_response();
        return iterate_singles(world, x, CC_vecfunction(PARTICLE), empty, empty, CT_LRCCS, info.parameters.iter_max_3D(), info);
    }

    static bool
    /// Iterates the singles equations for CCS, CC2, LRCC2
    /// The corresponding regulairzation tails of the doubles are updated in every iteration (therefore doubles are not marked as const)
    /// @param[in] : singles, the singles that are iterated
    /// @param[in] : singles2, additional set of singles for potential (in LRCC2 this are the Ground State singles)
    /// @param[in] : gs_doubles, Ground State doubles (Needed for CC2 and LRCC2)
    /// @param[in] : ex_doubles, Excited State doubles (Needed for LRCC2)
    /// @param[in] : ctype: the calculation type: CCS, CC2, CC2_response_
    /// @param[in] : maxiter: maxmial number of iterations
    /// @param[out]: true if the overall change of the singles is below 10*donv_6D
    iterate_singles(World& world, CC_vecfunction& singles, const CC_vecfunction singles2, Pairs<CCPair>& gs_doubles,
                    Pairs<CCPair>& ex_doubles, const CalcType ctype, const std::size_t maxiter, Info& info) {
        CCMessenger output(world);
        output.subsection("Iterate " + assign_name(ctype) + "-Singles");
        CCTimer time_all(world, "Overall Iteration of " + assign_name(ctype) + "-Singles");
        bool converged = true;

        CC_vecfunction old_singles(singles);
        for (auto& tmp : singles.functions)
            old_singles(tmp.first).function = copy(tmp.second.function);
        double old_omega=0.0;

        // KAIN solver
        typedef vector_function_allocator<double, 3> allocT;
        typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
        solverT solver(allocT(world, singles.size()));
        solver.do_print = (world.rank() == 0);

        print_size(world, singles.get_vecfunction(), "singles before iteration");

        for (size_t iter = 0; iter < maxiter; iter++) {
            // output.subsection("Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
            // CCTimer time(world, "Microiteration " + std::to_string(iter) + " of " + assign_name(ctype) + "-Singles");
            double omega = 0.0;
            if (ctype == CT_LRCC2) omega = singles.omega;
            else if (ctype == CT_LRCCS) omega = singles.omega;
            else if (ctype == CT_ADC2) omega = singles.omega;
            print("omega 1" ,omega);

            // consistency check
            switch (ctype) {
            case CT_CC2:
                if (singles.type != PARTICLE)
                    output.warning("iterate_singles: CC2 demanded but singles are not of type PARTICLE");
                break;
            case CT_MP2: MADNESS_EXCEPTION("Demanded Singles Calculation for MP2 ????", 1);
                break;
            case CT_LRCC2:
                if (singles.type != RESPONSE or singles2.type != PARTICLE)
                    output.warning("iterate_singles: CC2_response_ singles have wrong types");
                break;
            case CT_LRCCS:
                if (singles.type != RESPONSE)
                    output.warning("iterate_singles: CCS_response_ singles have wrong types");
                break;
            case CT_CISPD: MADNESS_EXCEPTION("Demanded Singles Calculation for CIS(D)", 1);
                break;
            case CT_ADC2:
                MADNESS_ASSERT(singles.type == RESPONSE);
                break;
            case CT_TEST: MADNESS_EXCEPTION("Iterate Singles not implemented for Experimental calculation", 1);
                break;
            default: MADNESS_EXCEPTION(
                    ("Unknown calculation type in iterate singles: " + assign_name(ctype)).c_str(), 1);
            }

            // get potentials
            CCTimer time_V(world, assign_name(ctype) + "-Singles Potential");
            vector_real_function_3d V;
            if (ctype == CT_CC2) V = CCPotentials::get_CC2_singles_potential_gs(world, singles, gs_doubles, info);
            else if (ctype == CT_LRCC2)
                V = CCPotentials::get_CC2_singles_potential_ex(world, singles2, gs_doubles, singles, ex_doubles, info);
            else if (ctype == CT_LRCCS) V = CCPotentials::get_CCS_potential_ex(world,singles,false, info);
            //            else if (ctype == CT_ADC2) V = CCOPS.get_ADC2_singles_potential(world, gs_doubles, singles, ex_doubles, info);
            else MADNESS_EXCEPTION("iterate singles: unknown type", 1);

            // add local coupling
            V-=compute_local_coupling(singles.get_vecfunction(),info);
            truncate(world, V);
            time_V.info(true, norm2(world, V));

            // update excitation energy
            if (ctype==CT_LRCC2 or ctype==CT_LRCCS or ctype==CT_ADC2) {
                old_omega=omega;
                omega = CCPotentials::compute_cis_expectation_value(world, singles, V, true, info);
                singles.omega = omega;
            }
            if (world.rank()==0 and info.parameters.debug())
                print("omega entering the update in the singles" ,omega);

            // make bsh operators
            scale(world, V, -2.0); // moved to BSHApply
            std::vector<std::shared_ptr<SeparatedConvolution<double, 3> > > G(singles.size());
            for (size_t i = 0; i < G.size(); i++) {
                const double bsh_eps = info.orbital_energies[i + info.parameters.freeze()] + omega;
                G[i] = std::shared_ptr<SeparatedConvolution<double, 3> >(
                        BSHOperatorPtr3D(world, sqrt(-2.0 * bsh_eps), info.parameters.lo(), info.parameters.thresh_bsh_3D()));
            }
            world.gop.fence();

            // apply bsh operators
            CCTimer time_applyG(world, "Apply G-Operators");
            vector_real_function_3d GV = apply<SeparatedConvolution<double, 3>, double, 3>(world, G, V);
            world.gop.fence();
            time_applyG.info();

            // apply Q-Projector to result
            QProjector<double,3> Q(info.mo_bra,info.mo_ket);
            GV = Q(GV);

            // Normalize Singles if it is excited state
            if (ctype == CT_LRCCS or ctype == CT_LRCC2 or ctype == CT_ADC2) {
                output("Normalizing new singles");
                const double norm=inner(GV,info.R_square*GV);
                scale(world, GV, 1.0 / norm);
            } else output("Singles not normalized");

            // residual
            const vector_real_function_3d residual = sub(world, singles.get_vecfunction(), GV);

            // information
            const Tensor<double> R2xinnerx = inner(world, info.R_square*singles.get_vecfunction(),
                                                   singles.get_vecfunction());
            const Tensor<double> R2GVinnerGV = inner(world, info.R_square*GV, GV);
            const Tensor<double> R2rinnerr = inner(world, info.R_square*residual, residual);
            const double R2vector_error = sqrt(R2rinnerr.sum());
            auto [rmsresidual, maxresidual]=CCPotentials::residual_stats(residual);

            // print information
            if (world.rank() == 0) std::cout << "\n\n-----Results of current interation:-----\n";
            if (world.rank() == 0)
                std::cout << "\nName: ||" << singles.name(0) << "||, ||GV" << singles.name(0) << ", ||residual||" << "\n";
            if (world.rank() == 0)
                std::cout << singles.name(0) << ": " << std::scientific << std::setprecision(info.parameters.output_prec())
                          << sqrt(R2xinnerx.sum()) << ", " << sqrt(R2GVinnerGV.sum()) << ", " << sqrt(R2rinnerr.sum())
                          << "\n----------------------------------------\n";
            for (size_t i = 0; i < GV.size(); i++) {
                if (world.rank() == 0)
                    std::cout << singles(i + info.parameters.freeze()).name() << ": " << std::scientific
                              << std::setprecision(info.parameters.output_prec())
                              << sqrt(R2xinnerx(i)) << ", " << sqrt(R2GVinnerGV(i)) << ", " << sqrt(R2rinnerr(i))
                              << "\n";
            }
            if (world.rank() == 0) std::cout << "\n----------------------------------------\n\n";

            // make second order update (only for response)
            if (ctype == CT_LRCC2 or ctype == CT_LRCCS) {
                double Rtmp = inner(world, info.R_square*residual, V).sum();
                double Rtmp2 = inner(world, info.R_square*GV, GV).sum();
                const double Rdelta = (0.5 * Rtmp / Rtmp2);
                if (world.rank() == 0) std::cout << "omega, second-order update (FYI): " << std::fixed
                              << std::setprecision(info.parameters.output_prec() + 2) << omega << ", " << Rdelta << "\n\n";
            }

            // update singles
            singles.omega = omega;
            vector_real_function_3d new_singles = truncate(GV);
            if (info.parameters.kain()) new_singles = solver.update(singles.get_vecfunction(), residual);
            if (info.parameters.debug()) print_size(world, new_singles, "new_singles");
            // if (ctype == CT_LRCCS or ctype == CT_LRCC2 or ctype == CT_ADC2) Nemo::normalize(new_singles, info.R);
            // if (info.parameters.debug()) print_size(world, new_singles, "new_singles normalized");

            for (size_t i = 0; i < GV.size(); i++) {
                singles(i + info.parameters.freeze()).function = copy(new_singles[i]);
            }

            // update reg_residues of doubles
            if (ctype==CT_CC2) update_reg_residues_gs(world, singles,gs_doubles, info);
            else if(ctype==CT_LRCC2) update_reg_residues_ex(world, singles2,singles,ex_doubles, info);

            if (world.rank()==0) CCPotentials::print_convergence(singles.name(0),rmsresidual,
                rmsresidual,omega-old_omega,iter);
            converged = (R2vector_error < info.parameters.dconv_3D());

            // time.info();
            if (converged) break;
            if (ctype == CT_LRCCS) break; // for CCS just one iteration to check convergence
        }
        time_all.info();
        print_size(world, singles.get_vecfunction(), "singles after iteration");

        // Assign the overall changes
        bool no_change = true;
        if (world.rank() == 0)
            std::cout << "Change in Singles functions after all the Microiterations" << std::endl;
        for (auto& tmp : singles.functions) {
            const double change = (tmp.second.function - old_singles(tmp.first).function).norm2();
            tmp.second.current_error = change;
            if (change > info.parameters.dconv_3D()) no_change = false;
            if (world.rank() == 0)
                std::cout << "Change of " << tmp.second.name() << "=" << tmp.second.current_error << std::endl;
        }
        // update reg_residues of doubles
        if (ctype == CT_CC2) update_reg_residues_gs(world, singles, gs_doubles, info);
        else if (ctype == CT_LRCC2) update_reg_residues_ex(world, singles2, singles, ex_doubles, info);

        //CCOPS.plot(singles);
        if (no_change) output("Change of Singles was below  = " + std::to_string(info.parameters.dconv_3D()) + "!");
        return no_change;
    }

    /// store singles to file
    void store_singles(const CC_vecfunction& singles, const int ex = -1) const;

    /// read singles from file or initialize new ones
    CC_vecfunction initialize_singles(const FuncType type, const int ex = -1) const;

    /// read pairs from file or initialize new ones
    bool initialize_pairs(Pairs<CCPair>& pairs, const CCState ftype, const CalcType ctype, const CC_vecfunction& tau,
                          const CC_vecfunction& x, const size_t extitation, const Info& info) const;

    static void
    update_reg_residues_gs(World& world, const CC_vecfunction& singles, Pairs<CCPair>& doubles, const Info& info);

    static void
    update_reg_residues_ex(World& world, const CC_vecfunction& singles, const CC_vecfunction& response, Pairs<CCPair>& doubles,
        const Info& info);

    /// Iterates a pair of the CC2 doubles equations
    bool
    iterate_pair(CCPair& pair, const CC_vecfunction& singles = CC_vecfunction(UNDEFINED)) const;

    bool
    iterate_adc2_pairs(Pairs<CCPair>& cispd, const CC_vecfunction& ccs);

    static bool
    iterate_lrcc2_pairs(World& world, const CC_vecfunction& cc2_s, const CC_vecfunction lrcc2_s,
                        Pairs<CCPair>& lrcc2_d, const Info& info);

    bool update_constant_part_cc2_gs(const CC_vecfunction& tau, CCPair& pair) {
        MADNESS_ASSERT(pair.ctype == CT_CC2);
        MADNESS_ASSERT(pair.type == GROUND_STATE);
        // make screening Operator
        real_convolution_6d Gscreen = BSHOperator<6>(world, sqrt(-2.0 * pair.bsh_eps), parameters.lo(),
                                                     parameters.thresh_bsh_6D());
        Gscreen.modified() = true;

        if (parameters.QtAnsatz())pair.constant_part = CCOPS.make_constant_part_cc2_Qt_gs(pair, tau, &Gscreen);
        else pair.constant_part = CCOPS.make_constant_part_cc2_gs(pair, tau, &Gscreen);
        save(pair.constant_part, pair.name() + "_const");
        return true;
    }

    bool update_constant_part_cispd(const CC_vecfunction& ccs, CCPair& pair) {
        MADNESS_ASSERT(pair.ctype == CT_CISPD);
        MADNESS_ASSERT(pair.type == EXCITED_STATE);
        MADNESS_ASSERT(pair.bsh_eps == CCOPS.get_epsilon(pair.i, pair.j) + ccs.omega);
        if (pair.constant_part.is_initialized()) return false; // the CIS(D) constant part does not change because there is no singles iteration (like MP2)
        // make screening Operator
        real_convolution_6d Gscreen = BSHOperator<6>(world, sqrt(-2.0 * pair.bsh_eps), parameters.lo(),
                                                     parameters.thresh_bsh_6D());
        Gscreen.modified() = true;

        if (parameters.QtAnsatz()) pair.constant_part = CCOPS.make_constant_part_cispd_Qt(pair, ccs, &Gscreen);
        else pair.constant_part = CCOPS.make_constant_part_cispd(pair, ccs, &Gscreen);
        save(pair.constant_part, pair.name() + "_const");
        return true;

    }

    bool update_constant_part_adc2(const CC_vecfunction& ccs, CCPair& pair) {
        std::cout << assign_name(pair.ctype);
        MADNESS_ASSERT(pair.ctype == CT_ADC2);
        MADNESS_ASSERT(pair.type == EXCITED_STATE);
        MADNESS_ASSERT(pair.bsh_eps == CCOPS.get_epsilon(pair.i, pair.j) + ccs.omega);
        // make screening Operator
        real_convolution_6d Gscreen = BSHOperator<6>(world, sqrt(-2.0 * pair.bsh_eps), parameters.lo(),
                                                     parameters.thresh_bsh_6D());
        Gscreen.modified() = true;

        if (parameters.QtAnsatz()) pair.constant_part = CCOPS.make_constant_part_cispd_Qt(pair, ccs, &Gscreen);
        else pair.constant_part = CCOPS.make_constant_part_cispd(pair, ccs, &Gscreen);
        save(pair.constant_part, pair.name() + "_const");
        return true;

    }

    bool update_constant_part_lrcc2(CCPair& pair, const CC_vecfunction& tau, const CC_vecfunction& x) {
        MADNESS_ASSERT(pair.ctype == CT_LRCC2);
        MADNESS_ASSERT(tau.type == PARTICLE);
        MADNESS_ASSERT(x.type == RESPONSE);

        // make screening Operator
        real_convolution_6d Gscreen = BSHOperator<6>(world, sqrt(-2.0 * pair.bsh_eps), parameters.lo(),
                                                     parameters.thresh_bsh_6D());
        Gscreen.modified() = true;

        if (parameters.QtAnsatz())pair.constant_part = CCOPS.make_constant_part_cc2_Qt_ex(pair, tau, x, &Gscreen);
        else pair.constant_part = CCOPS.make_constant_part_cc2_ex(pair, tau, x, &Gscreen);
        save(pair.constant_part, pair.name() + "_const");
        return true;
    }

    /// forward to the other function (converting CCPair to real_function)
    static Pairs<real_function_6d> compute_local_coupling(const std::vector<CCPair> &vpairs, const Info& info) {
        // create new pairs structure
        Pairs<CCPair> pairs;
        for (auto& tmp_pair : vpairs) pairs.insert(tmp_pair.i, tmp_pair.j, tmp_pair);
        auto ccpair2function = [](const CCPair& a) {return a.function();};
        return compute_local_coupling(pairs.convert<real_function_6d>(pairs,ccpair2function), info);
    };

    /// compute the coupling of singles function if orbitals are localized

    /// @return the coupling terms c_i = -\sum_(j\neq i) f_ij |\phi_j>  (for whatever phi is)
    static std::vector<real_function_3d> compute_local_coupling(const std::vector<real_function_3d>& singles,
        const Info& info) {

        MADNESS_CHECK_THROW(singles.size()>0,"compute_local_coupling: singles vector is empty");
        World& world=singles.front().world();
        auto active=Slice(info.parameters.freeze(),-1);
        Tensor<double> Fact=info.fock(active,active);
        for (int i=0; i<Fact.dim(0); ++i) Fact(i,i)=0.0;
        vector_real_function_3d fock_coupling = madness::transform(world, singles, Fact);
        return fock_coupling;
    }

    /// add the coupling terms for local MP2

    /// \sum_{k\neq i} f_ki |u_kj> + \sum_{l\neq j} f_lj |u_il>
    static Pairs<real_function_6d> compute_local_coupling(const Pairs<real_function_6d>& pairs, const Info& info);


    double solve_mp2_coupled(Pairs<CCPair> &doubles, Info& info);

    bool check_core_valence_separation(const Tensor<double>& fmat) const;

    /// @return     the new fock matrix
    Tensor<double> enforce_core_valence_separation(const Tensor<double>& fmat);
};


} /* namespace madness */

#endif /* CC2_H_ */
