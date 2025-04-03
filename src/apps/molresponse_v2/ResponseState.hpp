#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include "Perturbation.hpp"
#include <string>
#include <vector>

struct ResponseState {
  PerturbationType type;
  Perturbation perturbation;
  double frequency;

  std::vector<double> thresholds;  // Accuracy levels to loop over
  size_t current_thresh_index = 0; // Track which threshold we're working on
  bool is_converged = false;

  ResponseState(Perturbation pert, PerturbationType ptype, double freq,
                const std::vector<double> &thresh)
      : type(ptype), perturbation(pert), frequency(freq), thresholds(thresh) {}

  double current_threshold() const { return thresholds[current_thresh_index]; }

  bool at_final_accuracy() const {
    return current_thresh_index == thresholds.size() - 1;
  }

  void advance_protocol() {
    if (!at_final_accuracy()) {
      ++current_thresh_index;
      is_converged = false; // reset convergence at new protocol
    }
  }

  std::string description() const {
    return perturbationDescription() + " at freq " + std::to_string(frequency) +
           " (thresh=" + std::to_string(current_threshold()) + ")";
  }
  // Clearly named helper functions to get human-readable perturbation info
  [[nodiscard]] std::string perturbationDescription() const {
    switch (type) {
    case PerturbationType::Dipole:
      return "Dipole " +
             std::string(1,
                         std::get<DipolePerturbation>(perturbation).direction);
    case PerturbationType::NuclearDisplacement: {
      auto nuc = std::get<NuclearDisplacementPerturbation>(perturbation);
      return "NuclearDisplacement atom " + std::to_string(nuc.atom_index) +
             " direction " + nuc.direction;
    }
    case PerturbationType::Magnetic:
      return "Magnetic " +
             std::string(
                 1, std::get<MagneticPerturbation>(perturbation).direction);
    }
    return "Unknown";
  }
};

#endif // RESPONSE_STATE_HPP
