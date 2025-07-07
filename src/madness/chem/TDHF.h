/*
 * TDHF.h
 *
 *  Created on: Aug 11, 2016
 *      Author: kottmanj
 */

#ifndef SRC_APPS_CHEM_TDHF_H_
#define SRC_APPS_CHEM_TDHF_H_

#include<madness/chem/CCStructures.h>
#include<madness/chem/nemo.h>
#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/projector.h>
#include<madness/chem/SCFOperators.h>
#include <math.h>
#include<madness/chem/GuessFactory.h>
#include"madness/mra/commandlineparser.h"


namespace madness {
    /// the TDHF parameter class
    struct TDHFParameters : public QCCalculationParametersBase {

        static constexpr char const *tag = "response";

        TDHFParameters() {
            initialize_all();
        }

        TDHFParameters(const TDHFParameters &other) = default;

        /// todo: read_from_file compatible with dist. memory computation
        TDHFParameters(World &world, const commandlineparser& parser) {
            initialize_all();
            read_input_and_commandline_options(world, parser, "response");
        }

        void initialize_all() {

            // MRA stuff
            initialize < double > ("thresh", 1.e-5);
            initialize < double > ("econv", 1.e-5);
            initialize < double > ("dconv", 1.e-4);

            // physics
            initialize < std::string >
            ("calculation", "cis", "currently only cis=tda possible, TBD: thdf", {"cis"});
            initialize < std::string >
            ("response_kernel", "default", "default: corresponds to the ground state, libxc-notation otherwise");
            initialize < bool > ("triplet", false, "calculate triplet excitation energies (only works for CIS)");
            initialize < bool > ("do_oep", false, "use OEP potentials for the ground state exchange");
            initialize < std::size_t > ("nexcitations", 1,"number of excitation to be computed");
            initialize < long > ("freeze", -1, "the number of frozen occupied orbitals (-1: automatic)");
            initialize < std::string > ("irrep", "all", "compute only irreps of the respective point group");

            // solver
            initialize < size_t > ("maxiter", 25, "maximum number of iterations in the final iterations");
            initialize < std::size_t >
            ("kain_subspace", 8, "use kain (kain subspace<=0 == no kain, kain_subspace==1 should have the same effect)");

            // guess
            initialize < double > ("guess_econv", 1.e-4);
            initialize < double > ("guess_dconv", 1.e-3);
            initialize < std::size_t > ("iterating_excitations", 2);
            initialize < std::size_t > ("guess_excitations", 4);
            initialize < std::size_t > ("guess_occ_to_virt", 5);
            initialize < double >
            ("damping_width", 0.0, "every exop is multiplied with e^(-exponent*r2) to avoid noise at the boundaries");

            initialize < bool > ("debug", false);
            initialize < bool > ("plot", false);
//            initialize < bool > ("no_compute", false);
            initialize<int>  ("print_level",3,"0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");

//            initialize <std::vector<size_t>> ("restart", std::vector<size_t>(), "excitations which will be read from disk");
            initialize <std::string> ("restart", "iterate", "restart excitations from disk", {"no_restart","iterate","no_compute"});
            initialize <std::vector<size_t>> ("excitations", std::vector<size_t>(), "ordering of the excitations read from disk");


            initialize < std::string >
            ("guess_excitation_operators", "quadrupole", "guess type", {"dipole+", "quadrupole", "octopole", "custom"});

            /// add center of mass functions determined by the homo-energy
            /// will add s,px,py,pz functions in the center of mass with exponent: -(e_homo/c) and c=guess_cm is the value of this parameter
            initialize < double >
            ("guess_cm", 2.0, "center-of-mass functions, s/p shell with exponent -(e_homo/c)");

            /// use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)
            /// much faster
            initialize < bool >
            ("guess_diag", true, "use the diagonal approximation for the guess (only e_a -e_i terms in CIS matrix)");

            /// determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess
            /// guess_active_orbitals=0 is the same as guess_diag
            initialize < std::size_t >
            ("guess_active_orbitals", 0, "determine active orbitals in guess (for all inactive orbitals only the diagonal  e_a-e_i term is computed in the guess");


            initialize < bool >
            ("store_potential", true, "store the potential for orthogonalizations or recalculate it");

            initialize < size_t > ("guess_maxiter", 5, "maximum number of guess iterations ");

            //		/// determine how the virtuals for the guess are constructed: scf, external, custom, dipole, quadrupole
            //		/// scf: read in the ao set from scf (scales with system size)
            //		/// external: read in virtuals from disk
            //		/// custom or predefined strings like dipole, dipole+, ... : create virtuals form occupied orbitals by multiplying with polynomials
            //		/// |v> = |occ>*poly
            //		/// if custom is chosen:
            //		/// the polynomials can be determined by: exop x n1 y n2 z n3 c n4, x n11 ... x n44, ...  which will be n4*x^n1*y^n2*z^n3 + n44*x^n11* ...
            //		/// e.g. for a dipole guess enter the exop keyword 3 times as:
            //		/// exop x 1.0
            //		/// exop y 1.0
            //		/// exop z 1.0
            //		/// the options dipole, dipole+, dipole+diffuse and quadrupole give predefined exops without explicitly stating them
            //		/// see the end of TDHF.cc for predefined keys
            //		std::string guess_excitation_operators="dipole+";
            //
            /// Vector of strings which contains the polynomial excitation operators
            /// For this to be used the tda_guess key has to be "custom"
            /// The strings are given in a format like: "c c1 x x1 y y1 z z1, c c2 x x2 y y2 z z2, ..." which will be interpreted as: c1*x^x1*y^y1*z^z1 + c2*x^x2*y^y2*z^z2 + ....
            initialize < std::vector<std::string> >
            ("exops", {""}, "applies only if guess_excitation_operator is custom");


        }

        void set_derived_values(const std::shared_ptr<SCF> &scf);

        // physical part
        std::size_t nexcitations() const { return get<std::size_t>("nexcitations"); }

        long freeze() const { return get<long>("freeze"); }

        std::string irrep() const { return get<std::string>("irrep"); }

        bool triplet() const { return get<bool>("triplet"); }

        bool do_oep() const { return get<bool>("do_oep"); }

        std::string response_kernel() const { return get<std::string>("response_kernel"); }

        // precision
        double thresh() const { return get<double>("thresh"); }

        double econv() const { return get<double>("econv"); }

        double dconv() const { return get<double>("dconv"); }

        // restart and plotting
        bool debug() const { return get<bool>("debug"); }

        std::string restart() const { return get<std::string>("restart"); }
        bool no_compute() const { return (restart()=="no_compute"); }

        int print_level() const {return get<int>("print_level");}

        std::vector<size_t> excitations() const { return get<std::vector<size_t> >("excitations"); }

        bool plot() const { return get<bool>("plot"); }

        // solver parameters
        std::size_t iterating_excitations() const { return get<std::size_t>("iterating_excitations"); }

        std::size_t maxiter() const { return get<std::size_t>("maxiter"); }

        std::size_t kain_subspace() const { return get<std::size_t>("kain_subspace"); }

        bool store_potential() const { return get<bool>("store_potential"); }

        // guess parameters
//        std::size_t guess_occ_to_virt() const { return get<std::size_t>("guess_occ_to_virt"); }

        std::vector<std::string> exops() const { return get<std::vector<std::string> >("exops"); }

//        std::size_t guess_active_orbitals() const { return get<std::size_t>("guess_active_orbitals"); }

        bool guess_diag() const { return get<bool>("guess_diag"); }

        std::size_t guess_excitations() const { return get<std::size_t>("guess_excitations"); }

        std::string guess_excitation_operators() const { return get<std::string>("guess_excitation_operators"); }

        double damping_width() const { return get<double>("damping_width"); }

        double guess_cm() const { return get<double>("guess_cm"); }

        double guess_econv() const { return get<double>("guess_econv"); }

        double guess_dconv() const { return get<double>("guess_dconv"); }

        std::size_t guess_maxiter() const { return get<std::size_t>("guess_maxiter"); }

        /// make parameters for convolution operator
        typename CCConvolutionOperator<double,3>::Parameters get_ccc_parameters(const double lo) const {
            typename CCConvolutionOperator<double,3>::Parameters result;
            result.freeze = freeze();
            result.lo = lo;
            result.thresh_op = thresh();
            result.gamma = 1.0;
            return result;
        }
    }; // end of parameter class


/// The TDHF class
/// solves CIS/TDA equations and hopefully soon the full TDHF/TDDFT equations
class TDHF : public QCPropertyInterface {
public:

    TDHF(World &world, const TDHFParameters &parameters, std::shared_ptr<const Nemo> reference);

    TDHF(World &world, const commandlineparser &parser);

    TDHF(World &world, const commandlineparser &parser, std::shared_ptr<const Nemo> nemo);

    virtual ~TDHF() {}

    void initialize();

    std::string name() const {return "TDHF";};

    static void help() {
        print_header2("help page for CIS ");
        print("The CIS code computes Hartree-Fock and DFT excitations.");
        print("A moldft or nemo calculation will be performed automatically unless there is already a ");
        print("wave function on file ('restartdata.00000'). Both local and canonical orbitals can be ");
        print("used, for the latter individual irreps may be chosen for computation");
        print("Relevant parameters are derived from the reference calculation (moldft or nemo),");
        print("such as k, nuclear_correlation_factor, charge, etc");
        print("You can print all available calculation parameters by running\n");
        print("cis --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("cis --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    void print_frozen_orbitals() const;

    MolecularOrbitals<double,3> enforce_core_valence_separation(const Tensor<double>& fmat) const;

    static void print_parameters() {
        TDHFParameters param;
        print("default parameters for the CIS program are\n");
        param.print("response", "end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }

    virtual bool selftest() {
        return true;
    };

    ///  sets the reference wave function (nemo or oep)
    void set_reference(std::shared_ptr<NemoBase> reference) {
        reference_=reference;
    }

    std::shared_ptr<const NemoBase> get_reference() const {
        return reference_;
    }

    std::shared_ptr<SCF> get_calc() const {
        auto n=std::dynamic_pointer_cast<const Nemo>(reference_);
        if (not n) MADNESS_EXCEPTION("could not cast NemoBase to Nemo",1);
        return n->get_calc();
    }

    std::shared_ptr<const Nemo> get_nemo() const {
        auto n=std::dynamic_pointer_cast<const Nemo>(reference_);
        if (not n) MADNESS_EXCEPTION("could not cast NemoBase to Nemo",1);
        return n;
    }

    void prepare_calculation();

    CalculationParameters get_calcparam() const {
        auto n=std::dynamic_pointer_cast<const Nemo>(reference_);
        if (not n) MADNESS_EXCEPTION("could not cast NemoBase to Nemo",1);
        return n->param;
    }

    projector_irrep get_symmetry_projector() const {
        return symmetry_projector;
    }
    static int test(World &world, commandlineparser& parser);

    /// check consistency of the input parameters
    void check_consistency() const;

    /// plot planes and cubes
    void plot(const vector_real_function_3d &vf, const std::string &name) const;

    /// sort the xfunctions according to their excitation energy and name the excitation energies accordingly
    std::vector<CC_vecfunction> sort_xfunctions(std::vector<CC_vecfunction> x) const;

    /// print information
    void print_xfunctions(const std::vector<CC_vecfunction>& f, const std::string message) const;

    /// Initialize the CIS functions

    /// @param[in\out] on input the already obtained guess functions (or empty vector), on output new guess functions are added
    void initialize(std::vector<CC_vecfunction> &start) const;

    void symmetrize(std::vector<CC_vecfunction> &v) const;

    /// Solve the CIS equations

    /// @param[in/out] CC_vecfunction
    /// on input the guess functions (if empty or not enough the a guess will be generated)
    /// on output the solution
    std::vector<CC_vecfunction> solve_cis() const;

    /// analyze the root: oscillator strength and contributions from occupied orbitals
    void analyze(const std::vector<CC_vecfunction> &x) const;

    TDHFParameters get_parameters() const {return parameters;};

    std::vector<CC_vecfunction> get_converged_roots() const {return converged_roots;}

private:
    std::vector<CC_vecfunction> solve_cis(std::vector<CC_vecfunction> &start) const;

    /// Solve TDHF equations (not ready)
    void solve_tdhf(std::vector<CC_vecfunction> &guess) const;

    /// iterate the CIS guess vectors
    /// @param[in,out] x: on input the guess, on output the iterated guess
    /// see CC_Structures.h CCParameters class for convergence criteria
    bool iterate_cis_guess_vectors(std::vector<CC_vecfunction> &x) const;

    /// iterate the final CIS vectors
    /// @param[in,out] x: on input the guess, on output the iterated guess
    /// see CC_Structures.h CCParameters class for convergence criteria
    bool iterate_cis_final_vectors(std::vector<CC_vecfunction> &x) const;

    /// General function to iterate vectors
    /// @param[in,out] x: the CIS (or TDHF x) functions
    /// @param[in,out] the TDHF y functions (empty for CIS)
    /// @param[in] iterate_y, if true the y equation for TDHF is iterated
    /// @param[in] dconv: wavefunction convergence (for the vector norm of the vectorfunction)
    /// @param[in] econv: Energy convergece
    /// @param[in] iter: maximum number of iterations
    /// @param[in] kain: use kain if true (kainsubspace is controlled over CCParameters class)
    bool iterate_vectors(std::vector<CC_vecfunction> &x, const std::vector<CC_vecfunction> &y, bool iterate_y,
                         const double dconv, const double econv, const double iter, const bool kain) const;

    /// Apply the Greens function to a vector of vectorfunction with a given potential
    /// @param[in] x: the vector of vectorfunctions where G will be applied to
    /// @param[in] V: the vector of potentials to the vectorfunctions, will be cleared afterwards (potentials are all potentials excpet the nuclear: 2J - K + Q(2pJ - pK)
    /// @param[out] the vectorfunctions after G has been applied
    /// the energy is assumed to be stored in the CC_vecfunctions member omega
    /// the wavefunction error is stored in the CC_vecfunctions member current_error
    std::vector<vector_real_function_3d>
    apply_G(std::vector<CC_vecfunction> &x, std::vector<vector_real_function_3d> &V) const;

    /// Make the old CIS Guess
    /// the routine is now used to create virtuals
//    std::vector<CC_vecfunction> make_old_guess(const vector_real_function_3d &f) const;

    /// Create a set of virtual orbitals for the initial guess
    vector_real_function_3d make_virtuals() const;

    /// multiply excitation operators defined in the parameters with the seed functions
    /// @param[in] the seeds, define the function which are multiplied by the excitation operators
    /// @param[in] use_trigo, if false polynomials are used for excitation operators, else trigonometric functions (i.e. x^2y vs sin^2(x)*sin(y))
    /// Trigonometric functions are prefered since they are bounded (no weird behaviour at the boundaries for large exponents)
    vector_real_function_3d
    apply_excitation_operators(const vector_real_function_3d &seed, const bool &use_trigo = true) const;

    /// make the initial guess by explicitly diagonalizing a CIS matrix with virtuals from the make_virtuals routine
    vector<CC_vecfunction> make_guess_from_initial_diagonalization() const;

    /// canonicalize a set of orbitals (here the virtuals for the guess)
    vector_real_function_3d canonicalize(const vector_real_function_3d &v, Tensor<double> &veps) const;

    /// compute the CIS matrix for a given set of virtuals

    /// @param[in]	virtuals	the virtual orbitals
    /// @param[in]	veps		the orbital energies of the virtuals
    Tensor<double> make_cis_matrix(const vector_real_function_3d& virtuals, const Tensor<double> &veps) const;

    std::string filename_for_roots(const int ex) const {
        return get_calcparam().prefix()+"_root_"+std::to_string(ex);
    }

    /// Make the potentials to a given vector of vecfunctions (excitations)
    /// @param[in] The vector of excitations
    /// @param[out] The potentials
    std::vector<vector_real_function_3d> make_potentials(const std::vector<CC_vecfunction> &x) const;

    //    /// Make the CIS potential for a single excitation vector
    //	vecfuncT get_cis_potential(const CC_vecfunction& x) const {
    //		return CCOPS.make_cis_potential(x);
    //	}
    /// Make the TDA potential for a single excitation vector
    vector_real_function_3d get_tda_potential(const CC_vecfunction &x) const;

    /// Make the TDHF potential (not ready)
    std::vector<vector_real_function_3d>
    make_tdhf_potentials(std::vector<CC_vecfunction> &x, const std::vector<CC_vecfunction> &y) const;

    /// orthonormalize a vector of excitations
    /// @param[in,out] input: the excitations, output: the orthonormalized excitations
    /// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
    /// output: the transformed potentials
    void orthonormalize(std::vector<CC_vecfunction> &x, std::vector<vector_real_function_3d> &V) const;

    /// Calculate the perturbed fock matrix for a given vector of excitations
    /// @param[in] input: the excitations
    /// @param[in] input: the potentials, if empty the potentials will be recalculated but NOT stored
    Tensor<double> make_perturbed_fock_matrix(const std::vector<CC_vecfunction> &x,
                                              const std::vector<vector_real_function_3d> &V) const;

    Tensor<double> make_overlap_matrix(const std::vector<CC_vecfunction> &x) const;

    std::vector<vector_real_function_3d>
    transform(const std::vector<vector_real_function_3d> &x, const madness::Tensor<double> U) const {
        std::vector<CC_vecfunction> tmp;
        for (const auto &xi:x) tmp.push_back(CC_vecfunction(xi));
        std::vector<CC_vecfunction> tmp2 = transform(tmp, U);
        std::vector<vector_real_function_3d> result;
        for (const auto &xi:tmp2) result.push_back(xi.get_vecfunction());
        return result;
    }

    /// Interface to the SCF.h fock_transform function
    std::vector<CC_vecfunction>
    transform(const std::vector<CC_vecfunction> &x, const madness::Tensor<double> U) const;


    /// Helper function to initialize the const mo_bra and ket elements
    CC_vecfunction make_mo_bra(const std::vector<Function<double,3>>& amo) const {
        vector_real_function_3d tmp = get_reference()->get_ncf_ptr()->square()* amo;
        set_thresh(world, tmp, parameters.thresh());
        truncate(world, tmp);
        reconstruct(world, tmp);
        CC_vecfunction mo_bra(tmp, HOLE);
        return mo_bra;
    }

    CC_vecfunction make_mo_ket(const std::vector<Function<double,3>>& amo) const {
        vector_real_function_3d tmp = copy(world,amo);
        set_thresh(world, tmp, parameters.thresh());
        truncate(world, tmp);
        reconstruct(world, tmp);
        CC_vecfunction mo_ket(tmp, HOLE);
        return mo_ket;
    }

    double get_orbital_energy(const size_t i) const {
        auto n=std::dynamic_pointer_cast<const Nemo>(reference_);
        if (not n) MADNESS_EXCEPTION("could not cast NemoBase to Nemo",1);
        return n->get_calc()->aeps(i);
    }

    /// convenience
    vector_real_function_3d make_bra(const CC_vecfunction &ket) const {
        return make_bra(ket.get_vecfunction());
    }

    real_function_3d make_bra(const real_function_3d &ket) const {
        vector_real_function_3d v(1, ket);
        return make_bra(v).front();

    }

    /// maybe move this into nuclear_correlation class ?
    vector_real_function_3d make_bra(const vector_real_function_3d &ket) const {
        CCTimer time(world, "Make Bra");
        real_function_3d nucf = reference_->ncf->square();
        vector_real_function_3d result = mul(world, nucf, ket);
        time.info(parameters.debug());
        return result;
    }

    const vector_real_function_3d get_active_mo_ket() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_ket_.size(); i++) result.push_back(mo_ket_(i).function);
        return result;
    }

    const vector_real_function_3d get_active_mo_bra() const {
        vector_real_function_3d result;
        for (size_t i = parameters.freeze(); i < mo_ket_.size(); i++) result.push_back(mo_bra_(i).function);
        return result;
    }

    /// compute the oscillator strength in the length representation

    /// the oscillator strength is given by
    /// \f[
    /// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
    /// \f]
    /// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
    /// @param[in]  root    a converged root
    double oscillator_strength_length(const CC_vecfunction &x) const;

    /// compute the oscillator strength in the velocity representation

    /// the oscillator strength is given by
    /// \f[
    /// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
    /// \f]
    /// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
    /// @param[in]  root    a converged root
    double oscillator_strength_velocity(const CC_vecfunction &x) const;

    /// Fock matrix for occupied orbitals
    Tensor<double> F_occ;
    /// The MPI Communicator
    World &world;
    /// The Nemo structure (convenience)
    std::shared_ptr<const NemoBase> reference_;
    /// The TDHFParameters for the Calculations
    TDHFParameters parameters;
    /// Operator Structure which can handle intermediates (use for exchange with GS orbitals)
    /// Can be replaced by another potential manager
    std::shared_ptr<CCConvolutionOperator<double,3>> g12;
    /// MO bra and ket
    CC_vecfunction mo_ket_;
    CC_vecfunction mo_bra_;
    /// the Projector to the virtual space
    QProjector<double, 3> Q;
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
