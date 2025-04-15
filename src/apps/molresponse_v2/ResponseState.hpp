#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include "GroundStateData.hpp"
#include "Perturbation.hpp"
#include "vmra.h"
#include <SCF.h>
#include <filesystem>
#include <madness/chem/projector.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;
namespace fs = std::filesystem;

struct ResponseState {
  PerturbationType type;
  Perturbation perturbation;

  bool spin_restricted = false; // Is the system open shell?

  std::vector<double> frequencies;
  std::vector<double> thresholds; // Accuracy levels to loop over
  //
  size_t current_frequency_index = 0;
  size_t current_thresh_index = 0; // Track which threshold we're working on
  bool is_converged = false;

  ResponseState(Perturbation pert, PerturbationType ptype,
                const std::vector<double> &freq,
                const std::vector<double> &thresh, bool spin_restricted)
      : type(ptype), perturbation(pert), frequencies(freq), thresholds(thresh),
        current_frequency_index(0), current_thresh_index(0),
        spin_restricted(spin_restricted), is_converged(false) {}

  [[nodiscard]] double current_threshold() const {
    return thresholds[current_thresh_index];
  }
  [[nodiscard]] double current_frequency() const {
    return frequencies[current_frequency_index];
  }

  [[nodiscard]] bool at_final_frequency() const {
    return current_frequency_index == frequencies.size() - 1;
  }
  [[nodiscard]] bool at_final_threshold() const {
    return current_thresh_index == thresholds.size() - 1;
  }

  void set_frequency_index(size_t index) {
    if (index < frequencies.size()) {
      current_frequency_index = index;
    } else {
      throw std::out_of_range("Frequency index out of range");
    }
  }

  void advance_threshold() {
    if (!at_final_threshold()) {
      ++current_thresh_index;
      is_converged = false; // reset convergence at new protocol
    }
  }
  void advance_frequency() {
    if (!at_final_frequency()) {
      ++current_frequency_index;
    }
  }
  [[nodiscard]] bool is_static() const {
    return std::abs(current_frequency()) < 1e-8;
  }
  [[nodiscard]] bool is_dynamic() const { return !is_static(); }
  [[nodiscard]] bool is_spin_restricted() const { return spin_restricted; }

  [[nodiscard]] std::string response_filename() const {
    std::ostringstream oss;
    oss << "responses/" << perturbationDescription() << "_p"
        << current_threshold() << "_f" << current_frequency() << ".response";
    return oss.str();
  }
  [[nodiscard]] std::string response_filename(const size_t &thresh_index,
                                              const size_t &freq_index) const {
    std::ostringstream oss;
    oss << "responses/" << perturbationDescription() << "_p"
        << thresholds[thresh_index] << "_f" << frequencies[freq_index]
        << ".response";
    return oss.str();
  }

  [[nodiscard]] std::string description() const {
    return perturbationDescription() + " at freq " +
           std::to_string(current_frequency()) +
           " (thresh=" + std::to_string(current_threshold()) + ")";
  }

  // Clearly named helper functions to get human-readable perturbation info
  [[nodiscard]] std::string perturbationDescription() const {
    switch (type) {
    case PerturbationType::Dipole:
      return "Dipole_" +
             std::string(1,
                         std::get<DipolePerturbation>(perturbation).direction);
    case PerturbationType::NuclearDisplacement: {
      auto nuc = std::get<NuclearDisplacementPerturbation>(perturbation);
      return "NuclearDisplacement_atom_" + std::to_string(nuc.atom_index) +
             " direction " + nuc.direction;
    }
    case PerturbationType::Magnetic:
      return "Magnetic_" +
             std::string(
                 1, std::get<MagneticPerturbation>(perturbation).direction);
    }
    return "Unknown";
  }

  vector_real_function_3d
  perturbation_vector(World &world, const GroundStateData &ground_state) const {

    vector_real_function_3d Vp;

    switch (type) {
    case PerturbationType::Dipole: {
      auto dipole = std::get<DipolePerturbation>(perturbation);
      ;
      std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};
      std::vector<int> f(3, 0);
      f[dipole_map[dipole.direction]] = 1;
      real_function_3d d =
          real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));

      Vp = mul(world, d, ground_state.orbitals, true);
      Vp = ground_state.Qhat(Vp);
      truncate(world, Vp, FunctionDefaults<3>::get_thresh(), true);
      auto vp_norms = norm2s_T(world, Vp);

      return Vp;
    }
    }
    throw std::runtime_error("Unknown perturbation type");
  }
};

#endif // RESPONSE_STATE_HPP
