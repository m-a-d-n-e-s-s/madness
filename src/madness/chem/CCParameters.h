//
// Created by Florian Bischoff on 04.06.25.
//

#ifndef CCPARAMETER_H
#define CCPARAMETER_H

#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>
#include <madness/world/world.h>

namespace madness {
    /// Calculation Types used by CC2
    enum CalcType {
        CT_UNDEFINED, CT_MP2, CT_MP3, CT_CC2, CT_LRCCS, CT_LRCC2, CT_CISPD, CT_ADC2, CT_TDHF, CT_TEST
    };

    /// Calculation TDHFParameters for CC2 and TDA calculations
    /// Maybe merge this with calculation_parameters of SCF at some point, or split into TDA and CC
    struct CCParameters : public QCCalculationParametersBase {
        static constexpr char const* tag = "cc2";


        CCParameters() {
            initialize_parameters();
        };

        /// copy constructor
        CCParameters(const CCParameters& other) =default;

        /// ctor reading out the input file
        CCParameters(World& world, const commandlineparser& parser) {
            initialize_parameters();
            read_input_and_commandline_options(world,parser,"cc2");
            set_derived_values();
        };


        void initialize_parameters() {
            double thresh=1.e-3;
            double thresh_operators=1.e-6;
            initialize < std::string > ("calc_type", "mp2", "the calculation type", {"mp2", "mp3", "cc2", "cis", "lrcc2", "cispd", "adc2", "test"});
            initialize < double > ("lo", 1.e-7, "the finest length scale to be resolved by 6D operators");
            initialize < double > ("dmin", 1.0, "defines the depth of the special level");
            initialize < double > ("thresh_6d", thresh, "threshold for the 6D wave function");
            initialize < double > ("tight_thresh_6d", 0.1*thresh, "tight threshold for the 6D wave function");
            initialize < double > ("thresh_3d", 0.01*thresh, "threshold for the 3D reference wave function");
            initialize < double > ("tight_thresh_3d", 0.001*thresh, "tight threshold for the 3D reference wave function");
            initialize < double > ("thresh_bsh_3d", thresh_operators, "threshold for BSH operators");
            initialize < double > ("thresh_bsh_6d", thresh_operators, "threshold for BSH operators");
            initialize < double > ("thresh_poisson", thresh_operators, "threshold for Poisson operators");
            initialize < double > ("thresh_f12", thresh_operators, "threshold for Poisson operators");
            initialize < double > ("thresh_Ue", thresh_operators, "ue threshold");
            initialize < double > ("econv", thresh, "overal convergence threshold ");
            initialize < double > ("econv_pairs", 0.1*thresh, "convergence threshold for pairs");
            initialize < double > ("dconv_3d", 0.3*thresh, "convergence for cc singles");
            initialize < double > ("dconv_6d", 3.0*thresh, "convergence for cc doubles");
            initialize < std::size_t > ("iter_max", 10, "max iterations");
            initialize < std::size_t > ("iter_max_3d", 10, "max iterations for singles");
            initialize < std::size_t > ("iter_max_6d", 10, "max iterations for doubles");
            initialize < std::pair<int, int>> ("only_pair", {-1, -1}, "compute only a single pair");
            initialize < bool > ("restart", false, "restart");
            initialize < bool > ("no_compute", false, "no compute");
            initialize < bool > ("no_compute_gs", false, "no compute");
            initialize < bool > ("no_compute_mp2_constantpart", false, "no compute");
            initialize < bool > ("no_compute_response", false, "no compute");
            initialize < bool > ("no_compute_mp2", false, "no compute");
            initialize < bool > ("no_compute_cc2", false, "no compute");
            initialize < bool > ("no_compute_cispd", false, "no compute");
            initialize < bool > ("no_compute_lrcc2", false, "no compute");
            initialize < double > ("corrfac_gamma", 1.0, "exponent for the correlation factor");
            initialize < std::size_t > ("output_prec", 8, "for formatted output");
            initialize < bool > ("debug", false, "");
            initialize < bool > ("plot", false, "");
            initialize < bool > ("kain", true, "");
            initialize < std::size_t > ("kain_subspace", 3, "");
            initialize < long > ("freeze", -1, "number of frozen orbitals: -1: automatic");
            initialize < bool > ("test", false, "");
            // choose if Q for the constant part of MP2 and related calculations should be decomposed: GQV or GV - GO12V
            initialize < bool > ("decompose_Q", true, "always true",{true});
            // if true the ansatz for the CC2 ground state pairs is |tau_ij> = |u_ij> + Qtf12|titj>, with Qt = Q - |tau><phi|
            // if false the ansatz is the same with normal Q projector
            // the response ansatz is the corresponding response of the gs ansatz
            initialize < bool > ("QtAnsatz", true, "always true",{true});
            // a vector containing the excitations which shall be optizmized later (with CIS(D) or CC2)
            initialize < std::vector<size_t>>
            ("excitations", {}, "vector containing the excitations");
        }

        void set_derived_values();

        CalcType calc_type() const {
            std::string value = get<std::string>("calc_type");
            if (value == "mp2") return CT_MP2;
            if (value == "mp3") return CT_MP3;
            if (value == "cc2") return CT_CC2;
            if (value == "cis") return CT_LRCCS;
            if (value == "lrcc2") return CT_LRCC2;
            if (value == "cispd") return CT_CISPD;
            if (value == "adc2") return CT_ADC2;
            if (value == "test") return CT_TEST;
            MADNESS_EXCEPTION("faulty CalcType", 1);
        }

        bool response() const {return calc_type()==CT_ADC2 or calc_type()==CT_CISPD or calc_type()==CT_LRCC2 or calc_type()==CT_LRCCS;}
        double lo() const { return get<double>("lo"); }

        double dmin() const { return get<double>("dmin"); }

        double thresh_3D() const { return get<double>("thresh_3d"); }

        double tight_thresh_3D() const { return get<double>("tight_thresh_3d"); }

        double thresh_6D() const { return get<double>("thresh_6d"); }

        double tight_thresh_6D() const { return get<double>("tight_thresh_6d"); }

        double thresh_bsh_3D() const { return get<double>("thresh_bsh_3d"); }

        double thresh_bsh_6D() const { return get<double>("thresh_bsh_6d"); }

        double thresh_poisson() const { return get<double>("thresh_poisson"); }

        double thresh_f12() const { return get<double>("thresh_f12"); }

        double thresh_Ue() const { return get<double>("thresh_ue"); }

        double econv() const { return get<double>("econv"); }

        double econv_pairs() const { return get<double>("econv_pairs"); }

        double dconv_3D() const { return get<double>("dconv_3d"); }

        double dconv_6D() const { return get<double>("dconv_6d"); }

        std::size_t iter_max() const { return get<std::size_t>("iter_max"); }

        std::size_t iter_max_3D() const { return get<std::size_t>("iter_max_3d"); }

        std::size_t iter_max_6D() const { return get<std::size_t>("iter_max_6d"); }

        std::pair<int, int> only_pair() const { return get<std::pair<int, int>>("only_pair"); }

        bool restart() const { return get<bool>("restart"); }

        bool no_compute() const { return get<bool>("no_compute"); }

        bool no_compute_gs() const { return get<bool>("no_compute_gs"); }

        bool no_compute_mp2_constantpart() const { return get<bool>("no_compute_mp2_constantpart"); }

        bool no_compute_response() const { return get<bool>("no_compute_response"); }

        bool no_compute_mp2() const { return get<bool>("no_compute_mp2"); }

        bool no_compute_cc2() const { return get<bool>("no_compute_cc2"); }

        bool no_compute_cispd() const { return get<bool>("no_compute_cispd"); }

        bool no_compute_lrcc2() const { return get<bool>("no_compute_lrcc2"); }

        bool debug() const { return get<bool>("debug"); }

        bool plot() const { return get<bool>("plot"); }

        bool kain() const { return get<bool>("kain"); }

        bool test() const { return get<bool>("test"); }

        bool decompose_Q() const { return get<bool>("decompose_q"); }

        bool QtAnsatz() const { return get<bool>("qtansatz"); }

        std::size_t output_prec() const { return get<std::size_t>("output_prec"); }

        std::size_t kain_subspace() const { return get<std::size_t>("kain_subspace"); }

        long freeze() const { return get<long>("freeze"); }

        std::vector<std::size_t> excitations() const { return get<std::vector<std::size_t>>("excitations"); }

        double gamma() const {return get<double>("corrfac_gamma");}

        /// print out the parameters
        void information(World& world) const;

        /// check if parameters are set correct
        void sanity_check(World& world) const;

        void error(World& world, const std::string& msg) const {
            if (world.rank() == 0)
                std::cout << "\n\n\n\n\n!!!!!!!!!\n\nERROR IN CC_PARAMETERS:\n    ERROR MESSAGE IS: " << msg
                          << "\n\n\n!!!!!!!!" << std::endl;
            MADNESS_EXCEPTION("ERROR IN CC_PARAMETERS", 1);
        }

        size_t warning(World& world, const std::string& msg) const {
            if (world.rank() == 0) std::cout << "WARNING IN CC_PARAMETERS!: " << msg << std::endl;
            return 1;
        }

    };
} // namespace madness

#endif //CCPARAMETER_H
