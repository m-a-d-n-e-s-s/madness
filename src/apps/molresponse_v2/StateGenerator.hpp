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
  StateGenerator(const Molecule &mol,
                 const std::vector<MolecularProperty> &requested_properties,
                 const std::vector<double> &thresholds, bool spinrestricted)
      : molecule_(mol), requested_properties_(requested_properties),
        thresholds_(thresholds), spin_restricted_(spinrestricted) {}

  GeneratedStateData generateStates() const {
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
      ResponseState tmp{p, t, f, thresholds_, spin_restricted_};
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

    // 1) scan all requested_properties, accumulate into `table`
    for (auto const &prop : requested_properties_) {
      // start with exactly the user‐requested freqs
      std::vector<double> myfreqs = prop.frequencies;

      if (prop.type == MolecularPropertyType::Hyperpolarizability) {
        // also generate all pairwise sums ω_b + ω_c
        for (size_t b = 0; b < prop.frequencies.size(); ++b)
          for (size_t c = b; c < prop.frequencies.size(); ++c)
            myfreqs.push_back(prop.frequencies[b] + prop.frequencies[c]);
      }

      // dedupe & sort
      std::sort(myfreqs.begin(), myfreqs.end());
      myfreqs.erase(std::unique(myfreqs.begin(), myfreqs.end()), myfreqs.end());

      switch (prop.type) {
      case MolecularPropertyType::Polarizability:
        // just dipoles
        for (char d : prop.directions)
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          myfreqs);
        break;

      case MolecularPropertyType::Raman:
        // dipoles
        for (char d : prop.directions)
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          myfreqs);
        // **and** nuclear displacements
        for (int i = 0; i < molecule_.natom(); ++i)
          for (char d : {'x', 'y', 'z'})
            addPerturbation(PerturbationType::NuclearDisplacement,
                            NuclearDisplacementPerturbation{i, d}, myfreqs);
        break;

      case MolecularPropertyType::Hyperpolarizability:
        // same as polarizability but with the extended myfreqs
        for (char d : prop.directions)
          addPerturbation(PerturbationType::Dipole, DipolePerturbation{d},
                          myfreqs);
        break;
      }
    }

    // 2) finally, flatten into GeneratedStateData
    GeneratedStateData out;
    for (auto const &[key, e] : table) {
      std::vector<double> freqs(e.freqs.begin(), e.freqs.end());
      ResponseState st{e.pert, e.type, freqs, thresholds_, spin_restricted_};
      out.states.push_back(st);
      out.state_map[key] = st;
    }
    return out;
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
