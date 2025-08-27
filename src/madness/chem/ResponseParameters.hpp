#pragma once
#include <madness/mra/QCCalculationParametersBase.h>

using namespace madness;

struct ResponseParameters : public QCCalculationParametersBase {
  static constexpr char const *tag = "response";
  ResponseParameters(const ResponseParameters &other) = default;
  ResponseParameters(World &world, const commandlineparser &parser) : ResponseParameters() {
    read_input_and_commandline_options(world, parser, tag);
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
    initialize<size_t>("maxiter", 5, "maximum number of iterations");
    initialize<bool>("kain", false, "Turn on Krylov Accelarated Inexact Newton Solver");
    initialize<double>("maxrotn", .50, "Max orbital rotation per iteration");
    initialize<double>("maxbsh", 10, "Max bsh residual");
    initialize<size_t>("maxsub", 8, "size of iterative subspace ... set to 0 or 1 to disable");
    initialize<std::string>("xc", "hf", "XC input line");
    initialize<std::string>("hfexalg", "multiworld_row",
                            "hf exchange algorithm: choose from multiworld "
                            "(default), multiworld_row, smallmem, largemem");
    initialize<double>("dconv", 1e-6, "density convergence");
    initialize<bool>("step_restrict", true, "Toggles step restriction");
    initialize<std::vector<std::string>>("requested_properties", {"polarizability"},
                                         "properties to calculate (polarizability, hyperpolarizability, "
                                         "Raman.)");
    //** if properites are requested, then one should specify directions,
    // frequencies, and atom_indices(for nuclear response) */
    initialize<bool>("property", false, "Compute properties");
    initialize<bool>("dipole", false, "Compute linear dipole response");
    initialize<std::vector<double>>("dipole.frequencies", {0.0}, "frequencies for dipole response");
    initialize<std::string>("dipole.directions", "xyz", "directions for dipole response");
    initialize<bool>("nuclear", false, "Compute nuclear response");
    initialize<std::string>("nuclear.directions", "xyz", "directions for nuclear response");
    initialize<std::vector<double>>("nuclear.frequencies", {0.0}, "frequencies for nuclear response");
    initialize<std::vector<int>>("nuclear.atom_indices", {0}, "atom indices for nuclear response");
    initialize<bool>("quadratic", false, "Compute quadratic response properties from defined perturbations");
    initialize<std::string>("localize", "canon", "localization method", {"pm", "boys", "new", "canon"});
  }

  std::string get_tag() const override { return std::string(tag); }

public:
  using QCCalculationParametersBase::read_input_and_commandline_options;

  [[nodiscard]] std::string prefix() const { return get<std::string>("prefix"); }
  [[nodiscard]] std::string fock_json_file() const { return get<std::string>("fock_json_file"); }
  [[nodiscard]] std::string localize() const { return get<std::string>("localize"); }
  [[nodiscard]] std::string archive() const { return get<std::string>("archive"); }
  [[nodiscard]] std::string nwchem_dir() const { return get<std::string>("nwchem_dir"); }
  [[nodiscard]] bool nwchem() const { return get<bool>("nwchem"); }
  [[nodiscard]] int print_level() const { return get<int>("print_level"); }
  [[nodiscard]] bool step_restrict() const { return get<bool>("step_restrict"); }
  [[nodiscard]] size_t maxiter() const { return get<size_t>("maxiter"); }
  [[nodiscard]] double dconv() const { return get<double>("dconv"); }
  [[nodiscard]] bool quadratic() const { return get<bool>("quadratic"); }
  [[nodiscard]] bool kain() const { return get<bool>("kain"); }
  [[nodiscard]] size_t maxsub() const { return get<size_t>("maxsub"); }
  [[nodiscard]] std::string deriv() const { return get<std::string>("deriv"); }
  [[nodiscard]] std::string dft_deriv() const { return get<std::string>("dft_deriv"); }
  [[nodiscard]] std::string xc() const { return get<std::string>("xc"); }
  [[nodiscard]] std::string hfexalg() const { return get<std::string>("hfexalg"); }
  [[nodiscard]] double maxrotn() const { return get<double>("maxrotn"); }
  [[nodiscard]] bool property() const { return get<bool>("property"); }
  [[nodiscard]] std::vector<std::string> requested_properties() const {
    return get<std::vector<std::string>>("requested_properties");
  }
  [[nodiscard]] bool dipole() const { return get<bool>("dipole"); }
  [[nodiscard]] std::vector<double> dipole_frequencies() const {
    return get<std::vector<double>>("dipole.frequencies");
  }
  [[nodiscard]] std::string dipole_directions() const { return get<std::string>("dipole.directions"); }
  [[nodiscard]] bool nuclear() const { return get<bool>("nuclear"); }
  [[nodiscard]] std::vector<double> nuclear_frequencies() const {
    return get<std::vector<double>>("nuclear.frequencies");
  }
  [[nodiscard]] std::string nuclear_directions() const { return get<std::string>("nuclear.directions"); }
  [[nodiscard]] std::vector<int> nuclear_atom_indices() const { return get<std::vector<int>>("nuclear.atom_indices"); }
  [[nodiscard]] bool first_order() const { return get<bool>("first_order"); }
  [[nodiscard]] bool second_order() const { return get<bool>("second_order"); }
  [[nodiscard]] bool third_order() const { return get<bool>("third_order"); }

private:
  void validate_user_specified_properties() {
    // only validate if the user explicitly set requested_properties
    if (property()) {
      auto props = requested_properties();
      for (auto const &prop : props) {
        if (prop == "polarizability" || prop == "hyperpolarizability") {
          if (!is_user_defined("dipole.frequencies") || !is_user_defined("dipole.directions"))
            throw std::runtime_error("When requesting '" + prop +
                                     "', you must also set dipole.frequencies "
                                     "and dipole.directions.");
        }
        if (prop == "raman") {
          if (!is_user_defined("dipole.frequencies") || !is_user_defined("dipole.directions") ||
              !is_user_defined("nuclear.frequencies") || !is_user_defined("nuclear.directions") ||
              !is_user_defined("nuclear.atom_indices"))
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
