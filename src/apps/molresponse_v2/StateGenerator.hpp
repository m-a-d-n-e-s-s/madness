#pragma once
#include <molecule.h>

#include <set>
#include <string>

#include "../../madness/chem/ResponseParameters.hpp"
#include "MolecularProperty.hpp"
#include "Perturbation.hpp"
#include "ResponseState.hpp"

using namespace madness;

struct GeneratedStateData {
  std::vector<LinearResponseDescriptor> states;
  std::map<std::string, LinearResponseDescriptor> state_map;

  static void print_generated_state_map(const std::map<std::string, LinearResponseDescriptor> &state_map) {
    std::cout << "📋 Generated Response States:\n";
    std::cout << std::setw(5) << "#" << "  " << std::setw(40) << std::left << "State Description" << std::setw(20)
              << "Type" << std::setw(10) << "Num Freqs" << "\n";

    std::cout << std::string(90, '-') << "\n";

    size_t count = 0;
    for (const auto &[key, state] : state_map) {
      std::string type = perturbation_type_string(state.perturbation);
      std::cout << std::setw(5) << count++ << "  " << std::setw(40) << std::left << key << std::setw(20) << type
                << std::setw(10) << state.frequencies.size() << "\n";
    }
  }
};

class StateGenerator {
public:
  StateGenerator(const Molecule &mol, const std::vector<double> &thresholds, bool spinrestricted,
                 const ResponseParameters &response_parameters)
      : molecule_(mol), thresholds_(thresholds), spin_restricted_(spinrestricted), rp(response_parameters) {
    requested_properties_ = response_parameters.requested_properties();
  }

  [[nodiscard]] GeneratedStateData generateStates() const {
    struct Entry {
      Perturbation pert;
      std::set<double> freqs;
    };
    std::map<std::string, Entry> table;

    // helper to insert/merge one perturbation+freqs into table
    auto addPerturbation = [&](const Perturbation &p, const std::vector<double> &f) {
      // build a throwaway ResponseState so we get the exact key
      LinearResponseDescriptor tmp{p, f, thresholds_, spin_restricted_};
      std::string key = describe_perturbation(p);
      auto &e = table[key];
      // on first visit fill in type+pert
      if (e.freqs.empty()) {
        e.pert = p;
      }
      // merge in all f into the set
      e.freqs.insert(f.begin(), f.end());
    };

    auto dipole_dirs = rp.dipole_directions();
    auto dipole_freqs = rp.dipole_frequencies();
    auto nuclear_atom_indices = rp.nuclear_atom_indices();
    auto nuclear_directions = rp.nuclear_directions();
    auto freqs = rp.dipole_frequencies();
    auto nuclear_freqs = rp.nuclear_frequencies();
    auto num_freqs = freqs.size();
    auto num_nuclear_freqs = nuclear_freqs.size();

    enum class PropertyType { Alpha, Hessian, Beta, Raman };

    PropertyType prop_type;
    for (const auto &prop : requested_properties_) {

      auto prop_string = std::string(prop);
      // get rid of first and last characters
      prop_string = prop_string.substr(1, prop_string.size() - 2);

      if (prop_string == "polarizability") {
        prop_type = PropertyType::Alpha;
      } else if (prop_string == "hessian") {
        prop_type = PropertyType::Hessian;
      } else if (prop_string == "hyperpolarizability") {
        prop_type = PropertyType::Beta;
      } else if (prop_string == "raman") {
        prop_type = PropertyType::Raman;
      } else {
        throw std::runtime_error("Unknown property requested: " + prop);
      }

      auto augmented_dipole_freqs = dipole_freqs;
      // check if we need to augment dipole frequencies
      if (prop_type == PropertyType::Beta) {
        auto num_freqs = freqs.size();

        for (size_t b = 0; b < num_freqs; ++b)
          for (size_t c = b; c < num_freqs; ++c)
            augmented_dipole_freqs.push_back(freqs[b] + freqs[c]);
        std::sort(augmented_dipole_freqs.begin(), augmented_dipole_freqs.end());
        augmented_dipole_freqs.erase(std::unique(augmented_dipole_freqs.begin(), augmented_dipole_freqs.end()),
                                     augmented_dipole_freqs.end());
      }

      if (prop_type == PropertyType::Hessian) {
        // Hessian = all nuclear displacements at all nuclear frequencies
        for (auto atom_index : nuclear_atom_indices) {
          for (char d : nuclear_directions) {
            addPerturbation(NuclearDisplacementPerturbation{atom_index, d}, nuclear_freqs);
          }
        }
      }

      if (prop_type == PropertyType::Alpha) {
        for (char d : dipole_dirs) {
          addPerturbation(DipolePerturbation{d}, augmented_dipole_freqs);
        }
      } else if (prop_type == PropertyType::Raman) {
        for (char d : dipole_dirs) {
          addPerturbation(DipolePerturbation{d}, nuclear_freqs);
        }
        for (auto atom_index : nuclear_atom_indices) {
          for (char d : nuclear_directions) {
            addPerturbation(NuclearDisplacementPerturbation{atom_index, d}, nuclear_freqs);
          }
        }
      }
    }
    // 2) finally, flatten into GeneratedStateData
    GeneratedStateData out;
    for (auto const &[key, e] : table) {
      std::vector<double> freqs(e.freqs.begin(), e.freqs.end());
      LinearResponseDescriptor st{e.pert, freqs, thresholds_, spin_restricted_};
      out.states.push_back(st);
      out.state_map[key] = st;
    }
    return out;
  }

private:
  const Molecule &molecule_;
  std::vector<std::string> requested_properties_;
  std::vector<double> thresholds_;
  bool spin_restricted_;
  const ResponseParameters &rp;
};

enum class PropertyTensorType { Alpha, Beta };

struct PropertyComponentPlan {
  PropertyTensorType type;
  std::string description; // e.g., "alpha_xx", "beta_xyz"

  std::vector<std::string> required_perturbation_ids; // "dipole_x", "dipole_y", etc.
  std::vector<double> input_frequencies;              // ω or [ω₁, ω₂]
  double output_frequency = 0.0;

  std::vector<std::string> output_component_ids; // e.g., {"X_x", "X_y"}

  // For future compute step
  // std::function<void(...args)> compute_function;
};
