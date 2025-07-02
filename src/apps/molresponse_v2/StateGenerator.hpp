#pragma once
#include <molecule.h>

#include <set>
#include <string>

#include "MolecularProperty.hpp"
#include "../../madness/chem/ResponseParameters.hpp"
#include "ResponseState.hpp"

using namespace madness;

struct GeneratedStateData {
  std::vector<LinearResponseDescriptor> states;
  std::map<std::string, LinearResponseDescriptor> state_map;

  static void print_generated_state_map(
      const std::map<std::string, LinearResponseDescriptor> &state_map) {
    std::cout << "📋 Generated Response States:\n";
    std::cout << std::setw(5) << "#" << "  " << std::setw(40) << std::left
              << "State Description" << std::setw(20) << "Type" << std::setw(10)
              << "Num Freqs" << "\n";

    std::cout << std::string(90, '-') << "\n";

    size_t count = 0;
    for (const auto &[key, state] : state_map) {
      std::string type = (state.type == PerturbationType::Dipole) ? "Dipole"
                         : (state.type == PerturbationType::NuclearDisplacement)
                             ? "Nuclear"
                             : "Other";

      std::cout << std::setw(5) << count++ << "  " << std::setw(40) << std::left
                << key << std::setw(20) << type << std::setw(10)
                << state.frequencies.size() << "\n";
    }
  }
};

class StateGenerator {
 public:
  StateGenerator(const Molecule &mol, const std::vector<double> &thresholds,
                 bool spinrestricted,
                 const ResponseParameters &response_parameters)
      : molecule_(mol),
        thresholds_(thresholds),
        spin_restricted_(spinrestricted),
        rp(response_parameters) {
    requested_properties_ = response_parameters.requested_properties();
  }

  [[nodiscard]] GeneratedStateData generateStates() const {
    struct Entry {
      PerturbationType type;
      Perturbation pert;
      std::set<double> freqs;
    };
    std::map<std::string, Entry> table;

    // helper to insert/merge one perturbation+freqs into table
    auto addPerturbation = [&](PerturbationType t, const Perturbation &p,
                               const std::vector<double> &f) {
      // build a throwaway ResponseState so we get the exact key
      LinearResponseDescriptor tmp{p, t, f, thresholds_, spin_restricted_};
      std::string key = describe_perturbation(p);
      auto &e = table[key];
      // on first visit fill in type+pert
      if (e.freqs.empty()) {
        e.type = t;
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

    for (const auto &prop : requested_properties_) {
      auto myfreqs = dipole_freqs;
      if (prop == "hyperpolarizability") {
        auto num_freqs = freqs.size();
        auto freqs = rp.dipole_frequencies();

        for (size_t b = 0; b < num_freqs; ++b)
          for (size_t c = b; c < num_freqs; ++c)
            myfreqs.push_back(freqs[b] + freqs[c]);
        std::sort(myfreqs.begin(), myfreqs.end());
        myfreqs.erase(std::unique(myfreqs.begin(), myfreqs.end()),
                      myfreqs.end());
      }

      if (prop == "polarizability") {
        for (char d : dipole_dirs) {
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          myfreqs);
        }
      } else if (prop == "raman") {
        for (char d : dipole_dirs) {
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          nuclear_freqs);
        }
        auto nuclear_atom_indices = rp.nuclear_atom_indices();
        // 2) finally, flatten into GeneratedStateData
        GeneratedStateData out;
        for (auto const &[key, e] : table) {
          std::vector<double> freqs(e.freqs.begin(), e.freqs.end());
          LinearResponseDescriptor st{e.pert, e.type, freqs, thresholds_,
                                      spin_restricted_};
          out.states.push_back(st);
          out.state_map[key] = st;
        }
        return out;
        auto nuclear_directions = rp.nuclear_directions();
        for (size_t i = 0; i < nuclear_atom_indices.size(); ++i) {
          addPerturbation(PerturbationType::NuclearDisplacement,
                          NuclearDisplacementPerturbation{
                              nuclear_atom_indices[i], nuclear_directions[i]},
                          nuclear_freqs);
        }
      } else if (prop == "hyperpolarizability") {
        for (char d : dipole_dirs) {
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          myfreqs);
        }
      }
    }
    // 2) finally, flatten into GeneratedStateData
    GeneratedStateData out;
    for (auto const &[key, e] : table) {
      std::vector<double> freqs(e.freqs.begin(), e.freqs.end());
      LinearResponseDescriptor st{e.pert, e.type, freqs, thresholds_,
                                  spin_restricted_};
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
  std::string description;  // e.g., "alpha_xx", "beta_xyz"

  std::vector<std::string>
      required_perturbation_ids;          // "dipole_x", "dipole_y", etc.
  std::vector<double> input_frequencies;  // ω or [ω₁, ω₂]
  double output_frequency = 0.0;

  std::vector<std::string> output_component_ids;  // e.g., {"X_x", "X_y"}

  // For future compute step
  // std::function<void(...args)> compute_function;
};
