#pragma once
#include "MolecularProperty.hpp"
#include "ResponseState.hpp"
#include <molecule.h>
#include <set>
#include <string>

using namespace madness;

struct GeneratedStateData {
  std::vector<ResponseState> states;
  std::map<std::string, ResponseState> state_map;

  static void print_generated_state_map(
      const std::map<std::string, ResponseState> &state_map) {
    std::cout << "📋 Generated Response States:\n";
    std::cout << std::setw(5) << "#" << "  " << std::setw(40) << std::left
              << "State Description" << std::setw(20) << "Type" << std::setw(10)
              << "Static" << std::setw(10) << "Num Freqs" << "\n";

    std::cout << std::string(90, '-') << "\n";

    size_t count = 0;
    for (const auto &[key, state] : state_map) {
      std::string type = (state.type == PerturbationType::Dipole) ? "Dipole"
                         : (state.type == PerturbationType::NuclearDisplacement)
                             ? "Nuclear"
                             : "Other";

      std::cout << std::setw(5) << count++ << "  " << std::setw(40) << std::left
                << key << std::setw(20) << type << std::setw(10)
                << (state.is_static() ? "Yes" : "No") << std::setw(10)
                << state.frequencies.size() << "\n";
    }
  }
};

class StateGenerator {
public:
  StateGenerator(const Molecule &mol,
                 const std::vector<MolecularProperty> &requested_properties,
                 const std::vector<double> &thresholds, bool spinrestricted)
      : molecule_(mol), requested_properties_(requested_properties),
        thresholds_(thresholds), spin_restricted_(spinrestricted) {}

  [[nodiscard]] GeneratedStateData generateStates() const {
    std::set<std::string> seen_ids; // Prevent duplicates
    GeneratedStateData result;
    for (const auto &prop : requested_properties_) {
      switch (prop.type) {
      case MolecularPropertyType::Polarizability:
        for (char dir : prop.directions) {
          DipolePerturbation pert{dir};
          ResponseState state(pert, PerturbationType::Dipole, prop.frequencies,
                              thresholds_, spin_restricted_);
          auto state_pert_description = state.perturbationDescription();
          if (seen_ids.insert(state.description()).second) {

            result.states.push_back(state);
            result.state_map[state_pert_description] = state;
          }
        }
        break;

      case MolecularPropertyType::Raman:
        for (char dir : prop.directions) {
          DipolePerturbation dipole{dir};
          ResponseState dstate(dipole, PerturbationType::Dipole,
                               prop.frequencies, thresholds_, spin_restricted_);
          auto state_pert_description = dstate.perturbationDescription();
          if (seen_ids.insert(dstate.description()).second) {
            result.states.push_back(dstate);
            result.state_map[state_pert_description] = dstate;
          }
        }

        for (int i = 0; i < molecule_.natom(); ++i) {

          std::string atom_id = "n" + std::to_string(i);

          for (char dir : {'x', 'y', 'z'}) {
            NuclearDisplacementPerturbation nuc{i, dir};
            ResponseState nstate(nuc, PerturbationType::NuclearDisplacement,
                                 prop.frequencies, thresholds_,
                                 spin_restricted_);
            auto state_pert_description = nstate.perturbationDescription();
            if (seen_ids.insert(nstate.description()).second) {
              result.states.push_back(nstate);
              result.state_map[state_pert_description] = nstate;
            }
          }
        }
        break;

      case MolecularPropertyType::Hyperpolarizability:
        // TODO: Multi-frequency design
        auto input_frequencies = prop.frequencies;
        auto input_directions = prop.directions;

        std::set<double> unique_frequencies;

        for (int b = 0; b < input_frequencies.size(); ++b) {
          unique_frequencies.insert(input_frequencies[b]);
          for (int c = b; c < input_directions.size(); ++c) {
            unique_frequencies.insert(input_directions[c]);
            unique_frequencies.insert(input_directions[c] +
                                      input_frequencies[b]);
          }
        }
        auto freq_vector = std::vector<double>(unique_frequencies.begin(),
                                               unique_frequencies.end());

        for (char dir : input_directions) {
          DipolePerturbation dipole{dir};
          ResponseState dstate(dipole, PerturbationType::Dipole, freq_vector,
                               thresholds_, spin_restricted_);
          if (seen_ids.insert(dstate.description()).second) {
            result.states.push_back(dstate);
            result.state_map[dstate.description()] = dstate;
          }
        }

        break;
      }
    }

    return result;
  }

private:
  const Molecule &molecule_;
  std::vector<MolecularProperty> requested_properties_;
  std::vector<double> thresholds_;
  bool spin_restricted_;
};

class PropertyGenerator {
public:
  PropertyGenerator(const Molecule &mol,
                    const std::vector<MolecularProperty> &requested_properties,
                    const std::vector<double> &thresholds, bool spinrestricted)
      : molecule_(mol), requested_properties_(requested_properties),
        thresholds_(thresholds), spin_restricted_(spinrestricted) {}
  std::vector<ResponseState> generateStates() const {
    std::set<std::string> seen_ids; // Prevent duplicates
    std::vector<ResponseState> all_states;

    for (const auto &prop : requested_properties_) {
      switch (prop.type) {
      case MolecularPropertyType::Polarizability:
        for (char dir : prop.directions) {
          DipolePerturbation pert{dir};
          ResponseState state(pert, PerturbationType::Dipole, prop.frequencies,
                              thresholds_, spin_restricted_);
          if (seen_ids.insert(state.description()).second) {
            all_states.push_back(state);
          }
        }
        break;

      case MolecularPropertyType::Raman:
        for (char dir : prop.directions) {
          DipolePerturbation dipole{dir};
          ResponseState dstate(dipole, PerturbationType::Dipole,
                               prop.frequencies, thresholds_, spin_restricted_);
          if (seen_ids.insert(dstate.description()).second) {
            all_states.push_back(dstate);
          }
        }

        for (int i = 0; i < molecule_.natom(); ++i) {

          std::string atom_id = "n" + std::to_string(i);

          for (char dir : {'x', 'y', 'z'}) {
            NuclearDisplacementPerturbation nuc{i, dir};
            ResponseState nstate(nuc, PerturbationType::NuclearDisplacement,
                                 prop.frequencies, thresholds_,
                                 spin_restricted_);
            if (seen_ids.insert(nstate.description()).second) {
              all_states.push_back(nstate);
            }
          }
        }
        break;

      case MolecularPropertyType::Hyperpolarizability:
        // TODO: Multi-frequency design
        auto input_frequencies = prop.frequencies;
        auto input_directions = prop.directions;

        std::set<double> unique_frequencies;

        for (int b = 0; b < input_frequencies.size(); ++b) {
          unique_frequencies.insert(input_frequencies[b]);
          for (int c = b; c < input_directions.size(); ++c) {
            unique_frequencies.insert(input_directions[c]);
            unique_frequencies.insert(input_directions[c] +
                                      input_frequencies[b]);
          }
        }
        auto freq_vector = std::vector<double>(unique_frequencies.begin(),
                                               unique_frequencies.end());

        for (char dir : input_directions) {
          DipolePerturbation dipole{dir};
          ResponseState dstate(dipole, PerturbationType::Dipole, freq_vector,
                               thresholds_, spin_restricted_);
          if (seen_ids.insert(dstate.description()).second) {
            all_states.push_back(dstate);
          }
        }

        break;
      }
    }

    return all_states;
  }

private:
  const Molecule &molecule_;
  std::vector<MolecularProperty> requested_properties_;
  std::vector<double> thresholds_;
  bool spin_restricted_;
};

enum class PropertyTensorType { Alpha, Beta };

struct PropertyComponentPlan {
  PropertyTensorType type;
  std::string description; // e.g., "alpha_xx", "beta_xyz"

  std::vector<std::string>
      required_perturbation_ids;         // "dipole_x", "dipole_y", etc.
  std::vector<double> input_frequencies; // ω or [ω₁, ω₂]
  double output_frequency = 0.0;

  std::vector<std::string> output_component_ids; // e.g., {"X_x", "X_y"}

  // For future compute step
  // std::function<void(...args)> compute_function;
};
