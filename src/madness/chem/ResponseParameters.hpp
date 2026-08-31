#pragma once
#include <madness/mra/QCCalculationParametersBase.h>

using namespace madness;

struct ResponseParameters : public QCCalculationParametersBase {
    static constexpr char const* tag = "response";
    ResponseParameters(const ResponseParameters& other) = default;
    ResponseParameters(World& world, const commandlineparser& parser) : ResponseParameters() {
        read_input_and_commandline_options(world, parser, tag);
        if (!is_user_defined("state_parallel_groups")) {
            const auto world_groups = static_cast<size_t>(world.size() > 0 ? world.size() : 1);
            set_derived_value("state_parallel_groups", world_groups);
        }
        set_derived_properties();
        validate_user_specified_properties();
    }
    ResponseParameters() {
        initialize<std::string>("prefix", "response", "prefixes your output/restart/json/plot/etc files");
        initialize<std::string>("fock_json_file", "moldft.fock.json", "data file for fock matrix");
        initialize<std::string>("archive", "../moldft.restartdata", "file to read ground parameters from");
        initialize<bool>("nwchem", false, "Using nwchem files for intelligent starting guess");
        initialize<std::string>("nwchem_dir", "none", "Root name of nwchem files for intelligent starting guess");
        initialize<int>("print_level", 3, "0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");
        initialize<bool>("kain", false, "Turn on Krylov Accelarated Inexact Newton Solver");
        initialize<double>("maxrotn", .50, "Max orbital rotation per iteration");
        initialize<size_t>("maxsub", 8, "size of iterative subspace ... set to 0 or 1 to disable");
        initialize<std::string>("xc", "hf", "XC input line");
        initialize<std::string>("hfexalg", "multiworld_row",
                                "hf exchange algorithm: choose from multiworld "
                                "(default), multiworld_row, smallmem, largemem");
        initialize<double>("dconv", 1e-6, "density convergence");
        initialize<bool>("step_restrict", true, "Toggles step restriction");
        initialize<std::vector<std::string>>("requested_properties", {"polarizability"},
                                             "properties to calculate: polarizability, "
                                             "hyperpolarizability, raman");
        initialize<std::string>("engine", "v3",
                                "response engine backend: v3 (molresponse_v3 — alpha, beta, "
                                "single-component raman; ES experimental). 'v2' was removed "
                                "(M1) and is rejected with a migration hint.",
                                {"v2", "v3"});
        initialize<bool>("hdf5", false,
                         "opt-in HDF5-backed response-state restart I/O (writes .h5 "
                         "archives; readers auto-detect). Requires a build with "
                         "-DMADNESS_ENABLE_HDF5=ON; closed-shell states only.");
        initialize<int>("subworlds", 0,
                        "fan independent FD response states across this many "
                        "node-aligned subworlds (F2 state-parallel path); 0 = "
                        "single-World reference path. Use <= nodes, and note "
                        "PMIx on some clusters caps tasks/node at 2.");
        initialize<bool>("beta.shg", true,
                         "compute only SHG beta triplets (omegaB=omegaC, "
                         "omegaA=-(omegaB+omegaC))");
        initialize<bool>("beta.or", false,
                         "compute only optical-rectification beta triplets "
                         "(omegaB=0, omegaA=-omegaC)");
        initialize<bool>("beta.all_triplets", false, "compute full beta triplet grid over all (omegaB, omegaC) pairs");
        initialize<std::string>("state_parallel", "off", "state-level subgroup scheduling mode (off, auto, on)", {"off", "auto", "on"});
        initialize<size_t>("state_parallel_groups", 1, "number of processor groups for state-level subgroup scheduling");
        initialize<size_t>("state_parallel_min_states", 4,
                           "minimum number of generated states before auto state-parallel mode "
                           "activates");
        initialize<size_t>("state_parallel_property_group", 0, "subgroup id used for property assembly in state-parallel mode");
        initialize<size_t>("state_parallel_point_start_protocol", 1,
                           "first protocol index where state-parallel mode may fan out by "
                           "state-frequency point ownership");
        initialize<bool>("force_retry_removed_frequencies", false,
                         "allow retry of frequencies previously marked "
                         "remove_from_frequency_set");
        initialize<bool>("excited.enable", false, "enable excited-state bundle planning metadata scaffolding");
        initialize<size_t>("excited.num_states", 1, "number of excited states to target when enabled");
        initialize<bool>("excited.tda", false, "use Tamm-Dancoff approximation in excited-state stage");
        initialize<size_t>("excited.guess_max_iter", 5, "maximum iterations for excited-state guess stage");
        initialize<size_t>("excited.maxiter", 20, "maximum iterations for excited-state solve stage");
        initialize<size_t>("excited.maxsub", 8, "subspace size for excited-state iterative solves");
        initialize<size_t>("excited.owner_group", 0, "subgroup lane reserved for excited-state bundle execution");
        //** if properites are requested, then one should specify directions,
        // frequencies, and atom_indices(for nuclear response) */
        initialize<bool>("property", false, "Compute properties");
        initialize<bool>("dipole", false, "Compute linear dipole response");
        initialize<std::vector<double>>("dipole.frequencies", {0.0}, "frequencies for dipole response");
        initialize<std::string>("dipole.directions", "xyz", "directions for dipole response");
        initialize<bool>("nuclear", false, "Compute nuclear response");
        initialize<std::string>("nuclear.directions", "xyz", "directions for nuclear response");
        initialize<std::vector<double>>("nuclear.frequencies", {0.0}, "frequencies for nuclear response");
        initialize<bool>("quadratic", false, "Compute quadratic response properties from defined perturbations");
        initialize<std::string>("dalton.dir", "",
                                "seed-from-directory import (showcase W3): path to a DALTON "
                                "linear-response output directory (molden.inp + RSPVEC, loose or "
                                "in the *.tar.gz). The planned dipole FD states warm-start from "
                                "the projected RSPVEC vectors. Import-only — madness never "
                                "invokes DALTON. Geometry fingerprint mismatch is a hard error; "
                                "frequencies must match exactly 1-to-1.");
        initialize<std::string>("localize", "canon", "localization method", {"pm", "boys", "new", "canon"});
        initialize<size_t>("maxiter", 25, "maximum number of response iterations");
        initialize<std::string>("deriv", "abgv", "derivative method", {"abgv", "bspline", "ble"});
        initialize<std::string>("dft_deriv", "abgv", "derivative method for gga potentials", {"abgv", "bspline", "ble"});
    }

    std::string get_tag() const override {
        return std::string(tag);
    }

public:
    using QCCalculationParametersBase::read_input_and_commandline_options;

    [[nodiscard]] std::string prefix() const {
        return get<std::string>("prefix");
    }
    [[nodiscard]] std::string fock_json_file() const {
        return get<std::string>("fock_json_file");
    }
    [[nodiscard]] std::string localize() const {
        return get<std::string>("localize");
    }
    [[nodiscard]] std::string archive() const {
        return get<std::string>("archive");
    }
    [[nodiscard]] std::string nwchem_dir() const {
        return get<std::string>("nwchem_dir");
    }
    [[nodiscard]] bool nwchem() const {
        return get<bool>("nwchem");
    }
    [[nodiscard]] int print_level() const {
        return get<int>("print_level");
    }
    [[nodiscard]] std::string engine() const {
        return get<std::string>("engine");
    }
    [[nodiscard]] bool hdf5() const {
        return get<bool>("hdf5");
    }

    [[nodiscard]] int subworlds() const {
        return get<int>("subworlds");
    }
    [[nodiscard]] bool step_restrict() const {
        return get<bool>("step_restrict");
    }
    [[nodiscard]] size_t maxiter() const {
        return get<size_t>("maxiter");
    }
    [[nodiscard]] double dconv() const {
        return get<double>("dconv");
    }
    [[nodiscard]] bool quadratic() const {
        return get<bool>("quadratic");
    }
    [[nodiscard]] std::string dalton_dir() const {
        return get<std::string>("dalton.dir");
    }
    [[nodiscard]] bool kain() const {
        return get<bool>("kain");
    }
    [[nodiscard]] size_t maxsub() const {
        return get<size_t>("maxsub");
    }
    [[nodiscard]] std::string deriv() const {
        return get<std::string>("deriv");
    }
    [[nodiscard]] std::string dft_deriv() const {
        return get<std::string>("dft_deriv");
    }
    [[nodiscard]] std::string xc() const {
        return get<std::string>("xc");
    }
    [[nodiscard]] std::string hfexalg() const {
        return get<std::string>("hfexalg");
    }
    [[nodiscard]] double maxrotn() const {
        return get<double>("maxrotn");
    }
    [[nodiscard]] bool property() const {
        return get<bool>("property");
    }
    [[nodiscard]] std::vector<std::string> requested_properties() const {
        return get<std::vector<std::string>>("requested_properties");
    }
    [[nodiscard]] bool beta_shg() const {
        return get<bool>("beta.shg");
    }
    [[nodiscard]] bool beta_or() const {
        return get<bool>("beta.or");
    }
    [[nodiscard]] bool beta_all_triplets() const {
        return get<bool>("beta.all_triplets");
    }
    [[nodiscard]] bool dipole() const {
        return get<bool>("dipole");
    }
    [[nodiscard]] std::string state_parallel() const {
        return get<std::string>("state_parallel");
    }
    [[nodiscard]] size_t state_parallel_groups() const {
        return get<size_t>("state_parallel_groups");
    }
    [[nodiscard]] size_t state_parallel_min_states() const {
        return get<size_t>("state_parallel_min_states");
    }
    [[nodiscard]] size_t state_parallel_property_group() const {
        return get<size_t>("state_parallel_property_group");
    }
    [[nodiscard]] size_t state_parallel_point_start_protocol() const {
        return get<size_t>("state_parallel_point_start_protocol");
    }
    [[nodiscard]] bool force_retry_removed_frequencies() const {
        return get<bool>("force_retry_removed_frequencies");
    }
    [[nodiscard]] bool excited_enable() const {
        return get<bool>("excited.enable");
    }
    [[nodiscard]] size_t excited_num_states() const {
        return get<size_t>("excited.num_states");
    }
    [[nodiscard]] bool excited_tda() const {
        return get<bool>("excited.tda");
    }
    [[nodiscard]] size_t excited_guess_max_iter() const {
        return get<size_t>("excited.guess_max_iter");
    }
    [[nodiscard]] size_t excited_maxiter() const {
        return get<size_t>("excited.maxiter");
    }
    [[nodiscard]] size_t excited_maxsub() const {
        return get<size_t>("excited.maxsub");
    }
    [[nodiscard]] size_t excited_owner_group() const {
        return get<size_t>("excited.owner_group");
    }
    [[nodiscard]] std::vector<double> dipole_frequencies() const {
        return get<std::vector<double>>("dipole.frequencies");
    }
    [[nodiscard]] std::string dipole_directions() const {
        return get<std::string>("dipole.directions");
    }
    [[nodiscard]] bool nuclear() const {
        return get<bool>("nuclear");
    }
    [[nodiscard]] std::vector<double> nuclear_frequencies() const {
        return get<std::vector<double>>("nuclear.frequencies");
    }
    [[nodiscard]] std::string nuclear_directions() const {
        return get<std::string>("nuclear.directions");
    }

private:
    void validate_user_specified_properties() {
        // only validate if the user explicitly set requested_properties
        if (property()) {
            auto props = requested_properties();
            for (auto const& prop : props) {
                if (prop == "polarizability" || prop == "hyperpolarizability") {
                    if (!is_user_defined("dipole.frequencies") || !is_user_defined("dipole.directions"))
                        throw std::runtime_error("When requesting '" + prop +
                                                 "', you must also set dipole.frequencies "
                                                 "and dipole.directions.");
                }
                if (prop == "raman") {
                    if (!is_user_defined("dipole.frequencies") || !is_user_defined("dipole.directions") || !is_user_defined("nuclear.frequencies") || !is_user_defined("nuclear.directions"))
                        throw std::runtime_error("When requesting 'raman', you must set both dipole.* and "
                                                 "nuclear.* parameters.");
                }
            }
        }
    }

    void set_derived_properties() {
        // only override if user did NOT explicitly set requested_properties

        std::vector<std::string> props;
        if (!property()) {
            bool dip = dipole();
            bool nuc = nuclear();
            bool quad = quadratic();

            if (quad) {
                // quadratic response
                if (dip and nuc) {
                    // both nuclear & dipole kicks → Raman
                    props = {"polarizability", "hyperpolarizability", "raman"};
                } else {
                    // any pure quadratic dipole perturbations → α & β
                    props = {"polarizability", "hyperpolarizability"};
                }
            } else {
                // linear response only
                if (dip) {
                    props.push_back("polarizability");
                }
                // you could add a nuclear‐only property here if desired:
                // if (nuc) props.push_back("nuclear_response");
            }
        }
        if (!props.empty()) {
            // set_derived_value will only apply if precedence < derived
            set_derived_value("requested_properties", props);
        }
    }
};
