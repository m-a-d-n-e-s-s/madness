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

    /// compute the MP2 correlation energy
    static double compute_mp2_energy(const Pairs<CCPair>& pairs, const Info& info, const std::string msg="");
    static double compute_mp2_energy(const std::vector<CCPair>& pairs, const Info& info, const std::string msg="");

    /// compute the CC2 correlation energy
    static double compute_cc2_energy(const CC_vecfunction& singles, const Pairs<CCPair>& pairs,
        const Info& info, const std::string msg="");
    static double compute_cc2_energy(const CC_vecfunction& singles, const std::vector<CCPair>& pairs,
        const Info& info, const std::string msg="");

    double compute_mp3(const Pairs<CCPair>& mp2pairs, const Info& info) const {

        MP3 mp3(world, info);
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

    /// convenience function to iterate the CCS Response singles,
    /// makes the right call on the iterate_singles functions
    bool
    iterate_ccs_singles(CC_vecfunction& x, Info& info) const {
        Pairs<CCPair> empty;
        // CCOPS.clear_potentials(x);
        info.intermediate_potentials.clear_response();
        return iterate_singles(world, x, CC_vecfunction(PARTICLE), empty, empty, CT_LRCCS, info.parameters.iter_max_3D(), info);
    }

    /// Iterates the singles equations for CCS, CC2, LRCC2
    /// The corresponding regulairzation tails of the doubles are updated in every iteration (therefore doubles are not marked as const)
    /// @param[in] : singles, the singles that are iterated
    /// @param[in] : singles2, additional set of singles for potential (in LRCC2 this are the Ground State singles)
    /// @param[in] : gs_doubles, Ground State doubles (Needed for CC2 and LRCC2)
    /// @param[in] : ex_doubles, Excited State doubles (Needed for LRCC2)
    /// @param[in] : ctype: the calculation type: CCS, CC2, CC2_response_
    /// @param[in] : maxiter: maxmial number of iterations
    /// @param[out]: true if the overall change of the singles is below 10*donv_6D
    static bool
    iterate_singles(World& world, CC_vecfunction& singles, const CC_vecfunction singles2, const Pairs<CCPair>& gs_doubles,
                    const Pairs<CCPair>& ex_doubles, const CalcType ctype, const std::size_t maxiter, Info& info);

    /// return the file name for singles
    static std::string singles_name(const CalcType& ctype, const FuncType& type, int ex=-1) {
        std::string fname=assign_name(ctype)+"_"+madness::name(type,ex);
        return fname;
    }

    /// read singles from file or initialize new ones

    /// type: PARTICLE (cc2) or RESPONSE (lrcc2)
    /// default_to_zero: if true, initialize with zero functions, otherwise return empty vector
    /// ex: if type is RESPONSE, the excitation number
    static CC_vecfunction
    initialize_singles(World&, const CalcType& ctype, const FuncType type, int ex=-1);

    static CC_vecfunction
    initialize_singles_to_zero(World& world, const CalcType& ctype, const FuncType type, const Info& info);

    /// read pairs from file or initialize new ones
    bool initialize_pairs(Pairs<CCPair>& pairs, const CCState ftype, const CalcType ctype, const CC_vecfunction& tau,
                          const CC_vecfunction& x, const size_t extitation, const Info& info) const;

    /// Iterates a pair of the CC2 doubles equations
    bool
    iterate_pair(CCPair& pair, const CC_vecfunction& singles = CC_vecfunction(UNDEFINED)) const;

    bool
    iterate_adc2_pairs(Pairs<CCPair>& cispd, const CC_vecfunction& ccs);

    static bool
    iterate_lrcc2_pairs(World& world, const CC_vecfunction& cc2_s, const CC_vecfunction lrcc2_s,
                        Pairs<CCPair>& lrcc2_d, const Info& info);

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

    /// make sure the orbitals are block diagonalized

    /// changes the orbitals in member variable nemo
    /// will throw if the fock matrix is not block diagonal
    /// @return     the new fock matrix
    Tensor<double> enforce_core_valence_separation(const Tensor<double>& fmat);
};


} /* namespace madness */

#endif /* CC2_H_ */
