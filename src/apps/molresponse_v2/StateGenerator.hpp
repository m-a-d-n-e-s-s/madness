#pragma once
#include "MolecularProperty.hpp"
#include "ResponseState.hpp"
#include <molecule.h>
#include <set>
#include <string>

using namespace madness;

class StateGenerator {
public:
  StateGenerator(const Molecule &mol,
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
