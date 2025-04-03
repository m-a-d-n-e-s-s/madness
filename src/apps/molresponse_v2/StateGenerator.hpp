#pragma once
#include "MolecularProperty.hpp"
#include "ResponseState.hpp"
#include <molecule.h>
#include <set>

using namespace madness;

struct ResponseStateComparator {
  bool operator()(const ResponseState &lhs, const ResponseState &rhs) const {
    if (lhs.type != rhs.type)
      return lhs.type < rhs.type;
    if (lhs.frequency != rhs.frequency)
      return lhs.frequency < rhs.frequency;
    return lhs.perturbationDescription() < rhs.perturbationDescription();
  }
};

class StateGenerator {
public:
  StateGenerator(const Molecule &mol,
                 const std::vector<MolecularProperty> &requested_properties,
                 const std::vector<double> &thresholds)
      : molecule_(mol), requested_properties_(requested_properties),
        thresholds_(thresholds) {}

  std::vector<ResponseState> generateStates() const {
    std::set<std::string> seen_ids; // For deduplication
    std::vector<ResponseState> all_states;

    for (const auto &prop : requested_properties_) {
      switch (prop.type) {
      case MolecularPropertyType::Polarizability:
        for (double freq : prop.frequencies) {
          for (char dir : prop.directions) {
            DipolePerturbation pert{dir};
            ResponseState state(pert, PerturbationType::Dipole, freq,
                                thresholds_);
            if (seen_ids.insert(state.description()).second) {
              all_states.push_back(state);
            }
          }
        }
        break;

      case MolecularPropertyType::Raman:
        for (double freq : prop.frequencies) {
          for (char dir : prop.directions) {
            DipolePerturbation dipole{dir};
            ResponseState state(dipole, PerturbationType::Dipole, freq,
                                thresholds_);
            if (seen_ids.insert(state.description()).second) {
              all_states.push_back(state);
            }
          }
          for (int i = 0; i < molecule_.natom(); ++i) {
            std::string atom_id = "nuclear_" + std::to_string(i);
            for (char dir : {'x', 'y', 'z'}) {
              NuclearDisplacementPerturbation nd{i, dir};
              ResponseState state(nd, PerturbationType::NuclearDisplacement,
                                  freq, thresholds_);
              if (seen_ids.insert(state.description()).second) {
                all_states.push_back(state);
              }
            }
          }
        }
        break;

      case MolecularPropertyType::Hyperpolarizability:
        // TODO: Add logic for multi-frequency perturbation sets
        break;
      }
    }

    return all_states;
  }

private:
  const Molecule &molecule_;
  std::vector<MolecularProperty> requested_properties_;
  std::vector<double> thresholds_;
};
