// Copyright 2021 Adrian Hurtado

/// \file ResponseParameters
/// \brief Input parameters for a response calculation.

#ifndef SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_
#define SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_

#include <string>
#include <vector>

#include <madness/chem/molecule.h>
#include <madness/chem/xcfunctional.h>
#include <madness/mra/mra.h>
#include <madness/world/parallel_archive.h>
#include <molresponse/ground_parameters.h>
#include "madness/mra/QCCalculationParametersBase.h"

namespace madness
{

    struct ResponseParameters : public QCCalculationParametersBase
    {
        ResponseParameters(const ResponseParameters &other) = default;
        ResponseParameters(World &world, const commandlineparser &parser) : ResponseParameters()
        {
            read_input_and_commandline_options(world, parser, "response");
            // convenience option -- needs to be moved to the MolecularOptimizer class
        }
        ResponseParameters()
        {
            initialize<std::string>("archive", "../moldft.restartdata", "file to read ground parameters from");
            initialize<bool>("nwchem", false, "Using nwchem files for intelligent starting guess");
            initialize<std::string>("nwchem_dir", "none", "Root name of nwchem files for intelligent starting guess");
            initialize<size_t>("states", 1, "Number of excited states requested");
            initialize<int>("print_level", 3, "0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");
            initialize<bool>("tda", false, "turn on Tam-Danchof approximation (excitations energy");
            initialize<bool>("first_run", true, "Are we on the default guess");
            initialize<bool>("plot", false, "turn on plotting of final orbitals. Output format is .vts");
            initialize<bool>("plot_range", false, "controls which orbitals will be plotted");
            initialize<std::vector<int>>("plot_data", std::vector<int>{0}, "Orbitals to plot");
            initialize<std::vector<double>>("plot_cell", std::vector<double>(), "lo-hi plot cell (default is all space)");
            initialize<double>("plot_l", -1.0, "Controls the plotting box size");
            initialize<size_t>("plot_pts", 81, "Controls number of points in plots");
            initialize<bool>("plot_all_orbitals", false, "Turn on 2D plotting of response orbitals ");
            initialize<size_t>("maxiter", 25, "maximum number of iterations");
            initialize<double>("dconv", 1.e-4, "recommended values: 1.e-4 < dconv < 1.e-8");
            initialize<bool>("conv_only_dens", false, "if true remove bsh_residual from convergence criteria (deprecated)");
            initialize<bool>("dconv_set", false, "Convergence flage for the orbtial density");
            initialize<bool>("guess_xyz", true, "ExcitedState intial guess functions ground MO * <x,y,z>");
            initialize<double>("lo", 1.e-10, "smallest length scale we need to resolve");
            initialize<std::vector<double>>("protocol", {1.e-4, 1.e-6}, "Defines convergence and truncation protocol");
            initialize<size_t>("larger_subspace", 0,
                               "Number of iterations to diagonalize in a subspace "
                               "consisting of old and new vectors");
            initialize<int>("k", -1, "polynomial order");
            initialize<std::string>("deriv", "abgv", "derivative method", {"abgv", "bspline", "ble"});
            initialize<std::string>("dft_deriv", "abgv", "derivative method for gga potentials", {"abgv", "bspline", "ble"});
            initialize<bool>("random", true, "Use random guess for initial response functions");
            initialize<bool>("store_potential", true, "Store the potential instead of computing each iteration");
            initialize<bool>("e_range", false, "Use an energy range to excite from");
            initialize<double>("e_range_lo", 0, "Energy range (lower end) for orbitals to excite from");
            initialize<double>("e_range_hi", 1, "Energy range (upper end) for orbitals to excite from");
            initialize<bool>("plot_initial", false, "Flag to plot the ground state orbitals read in from archivie");
            initialize<bool>("restart", false, "Flag to restart scf loop from file");
            initialize<std::string>("restart_file", "none", "file to read ground parameters from");
            initialize<bool>("kain", false, "Turn on Krylov Accelarated Inexact Newton Solver");
            initialize<double>("maxrotn", 100, "Max orbital rotation per iteration");
            initialize<double>("maxbsh", 10, "Max bsh residual");
            initialize<size_t>("maxsub", 10, "size of iterative subspace ... set to 0 or 1 to disable");
            initialize<std::string>("xc", "hf", "XC input line");
            initialize<std::string>("hfexalg", "multiworld",
                                    "hf exchange algorithm: choose from multiworld "
                                    "(default), multiworld_row, smallmem, largemem");
            initialize<bool>("save", false, "if true save orbitals to disk");
            initialize<std::string>("save_file", "none", "File name to save orbitals for restart");
            initialize<bool>("save_density", false, "Flag to save density at each iteration");
            initialize<int>("vnucextra", 2, "load balance parameter for nuclear pot");
            initialize<int>("loadbalparts", 2, "??");
            initialize<std::string>("save_density_file", "none", "File name to save density for restart");
            initialize<bool>("load_density", false, "Flag to load density for restart");
            initialize<std::string>("load_density_file", "none", "File name to load density for restart");
            initialize<size_t>("guess_max_iter", 5, "maximum number of guess iterations");
            initialize<std::string>("calc_type", "full", "full,static,tda");
            initialize<bool>("excited_state", false, "Flag to turn on excited state calc");
            initialize<bool>("first_order", false, "Flag to turn on first order response calc");
            initialize<bool>("second_order", false, "Flag to turn on first order response calc");
            initialize<std::string>("perturbation", "dipole", "dipole,nuclear");
            initialize<bool>("dipole", false, "Sets RHS to dipole operator 3 x num_orbitals");
            initialize<bool>("quadratic", false, "Compute quadratic response");
            initialize<std::vector<double>>("freq_range", {0.0}, "Frequency range for quadratic response");
            initialize<bool>("nuclear", false, "Sets RHS to nuclear derivative 3 x num_atoms x num_orbitals");
            initialize<double>("omega", 0.0, "Incident energy for dynamic response");
            initialize<double>("l", 20, "user coordinates box size");
            initialize<size_t>("num_orbitals", 0, "number of ground state orbtials");
            initialize<bool>("spinrestricted", true, "is spinrestricted calculation");
            initialize<std::string>("localize", "canon", "localization method", {"pm", "boys", "new", "canon"});
        }

    public:
        using QCCalculationParametersBase::read_input_and_commandline_options;

        [[nodiscard]] std::string perturbation() const { return get<std::string>("perturbation"); }
        [[nodiscard]] std::string localize() const { return get<std::string>("localize"); }
        [[nodiscard]] std::string archive() const { return get<std::string>("archive"); }
        [[nodiscard]] std::string calc_type() const { return get<std::string>("calc_type"); }
        [[nodiscard]] std::string nwchem_dir() const { return get<std::string>("nwchem_dir"); }
        [[nodiscard]] bool nwchem() const { return get<bool>("nwchem"); }
        [[nodiscard]] size_t num_states() const { return get<size_t>("states"); }
        [[nodiscard]] size_t num_orbitals() const { return get<size_t>("num_orbitals"); }
        [[nodiscard]] int print_level() const { return get<int>("print_level"); }
        [[nodiscard]] bool tda() const { return get<bool>("tda"); }
        [[nodiscard]] bool plot() const { return get<bool>("plot"); }
        [[nodiscard]] double plot_l() const { return get<double>("plot_l"); }
        [[nodiscard]] size_t plot_pts() const { return get<size_t>("plot_pts"); }
        [[nodiscard]] bool plot_all_orbitals() const { return get<bool>("plot_all_orbitals"); }
        [[nodiscard]] size_t maxiter() const { return get<size_t>("maxiter"); }
        [[nodiscard]] double dconv() const { return get<double>("dconv"); }
        [[nodiscard]] bool guess_xyz() const { return get<bool>("guess_xyz"); }
        [[nodiscard]] double lo() const { return get<double>("lo"); }
        [[nodiscard]] std::vector<double> protocol() const { return get<std::vector<double>>("protocol"); }
        [[nodiscard]] size_t larger_subspace() const { return get<size_t>("larger_subspace"); }
        [[nodiscard]] int k() const { return get<int>("k"); }
        [[nodiscard]] bool random() const { return get<bool>("random"); }
        [[nodiscard]] bool store_potential() const { return get<bool>("store_potential"); }
        [[nodiscard]] vector<double> freq_range() const { return get<vector<double>>("freq_range"); }
        [[nodiscard]] bool quadratic() const { return get<bool>("quadratic"); }
        [[nodiscard]] bool plot_initial() const { return get<bool>("plot_initial"); }
        [[nodiscard]] bool restart() const { return get<bool>("restart"); }
        [[nodiscard]] std::string restart_file() const { return get<std::string>("restart_file"); }
        [[nodiscard]] bool kain() const { return get<bool>("kain"); }
        [[nodiscard]] size_t maxsub() const { return get<size_t>("maxsub"); }
        [[nodiscard]] std::string deriv() const { return get<std::string>("deriv"); }
        [[nodiscard]] std::string dft_deriv() const { return get<std::string>("dft_deriv"); }
        [[nodiscard]] std::string xc() const { return get<std::string>("xc"); }
        [[nodiscard]] std::string hfexalg() const { return get<std::string>("hfexalg"); }
        [[nodiscard]] bool save() const { return get<bool>("save"); }
        [[nodiscard]] std::string save_file() const { return get<std::string>("save_file"); }
        [[nodiscard]] size_t guess_max_iter() const { return get<size_t>("guess_max_iter"); }
        [[nodiscard]] bool property() const { return get<bool>("property"); }
        [[nodiscard]] int loadbalparts() const { return get<int>("loadbalparts"); }
        [[nodiscard]] bool excited_state() const { return get<bool>("excited_state"); }
        [[nodiscard]] bool first_order() const { return get<bool>("first_order"); }
        [[nodiscard]] bool second_order() const { return get<bool>("second_order"); }
        [[nodiscard]] bool third_order() const { return get<bool>("third_order"); }
        [[nodiscard]] bool dipole() const { return get<bool>("dipole"); }
        [[nodiscard]] bool nuclear() const { return get<bool>("nuclear"); }
        [[nodiscard]] std::string d2_types() const { return get<std::string>("d2_types"); }
        [[nodiscard]] double omega() const { return get<double>("omega"); }
        [[nodiscard]] double L() const { return get<double>("l"); }
        [[nodiscard]] bool spinrestricted() const { return get<bool>("spinrestricted"); }
        void set_ground_state_calculation_data(const GroundStateCalculation &g_params)
        {

            set_derived_value<size_t>("num_orbitals", g_params.n_orbitals());
            set_derived_value<bool>("spinrestricted", g_params.is_spinrestricted());
            set_derived_value<double>("l", g_params.get_L());
            set_derived_value<double>("lo", g_params.molecule().smallest_length_scale());
            set_derived_value<std::string>("xc", g_params.get_xc());
            set_derived_value<std::string>("localize", g_params.get_localize_method());
        }

        void set_derived_values(World &world, const Molecule &molecule)
        {
            // read the parameters from file and brodcast
            // tag
            vector<std::string> calculation_type;
            vector<bool> calc_flags;
            if (excited_state())
            {
                if (tda())
                {
                    set_derived_value<std::string>("calc_type", "tda");
                }
                else
                {
                    set_derived_value<std::string>("calc_type", "full");
                }
            }
            else if (dipole())
            {
                set_derived_value<std::string>("perturbation", "dipole");
                set_derived_value<size_t>("states", 3);
            }
            else if (nuclear())
            {
                set_derived_value<std::string>("perturbation", "nuclear");
                set_derived_value<size_t>("states", 3 * molecule.natom());
            }
        }

        // convenience getters
        [[nodiscard]] double econv() const { return get<double>("econv"); }
        [[nodiscard]] bool first_run() const { return get<bool>("first_run"); }
        [[nodiscard]] std::string local() const { return get<std::string>("local"); }
    };

    void from_json(const nlohmann::json &, ResponseParameters &p);

    bool operator==(const ResponseParameters &p1, const ResponseParameters &p2);

    bool operator!=(const ResponseParameters &p1, const ResponseParameters &p2);

    // convenience getters
    void to_json(nlohmann::json &j);
}; // namespace madness

// namespace madness

#endif // SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_
