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

struct AbstractResponseDescriptor {
  virtual bool is_spin_restricted() const = 0;
  virtual bool is_static() const = 0;
  virtual std::string response_filename() const = 0;
  virtual std::string response_filename(const size_t &thresh_index,
                                        const size_t &freq_index) const = 0;
  virtual ~AbstractResponseDescriptor() = default;
};

struct ResponseState : public AbstractResponseDescriptor {
  PerturbationType type;
  Perturbation perturbation;

  bool spin_restricted = false; // Is the system open shell?

  std::vector<double> frequencies;
  std::map<double, int> frequency_map; // Frequency to index map
  std::vector<double> thresholds;      // Accuracy levels to loop over
  //
  size_t current_frequency_index = 0;
  size_t current_thresh_index = 0; // Track which threshold we're working on
  bool is_converged = false;
  ResponseState() = default;

  ResponseState(Perturbation pert, PerturbationType ptype,
                const std::vector<double> &freq,
                const std::vector<double> &thresh, bool spin_restricted)
      : type(ptype), perturbation(pert), frequencies(freq), thresholds(thresh),
        current_frequency_index(0), current_thresh_index(0),
        spin_restricted(spin_restricted), is_converged(false) {

    for (size_t i = 0; i < frequencies.size(); ++i) {
      frequency_map[frequencies[i]] = i;
    }
  }

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
    std::ostringstream oss;
    oss << perturbationDescription() << " at freq " << current_frequency()
        << " (thresh=" << std::scientific << current_threshold() << ")";
    return oss.str();
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

struct SecondOrderResponseState : public AbstractResponseDescriptor {
  PerturbationType type; // likely "Dipole" for VBC
  std::pair<DipolePerturbation, DipolePerturbation> perturbations;
  std::pair<double, double> frequencies;
  double threshold;
  bool spin_restricted = true;

  SecondOrderResponseState(DipolePerturbation p1, DipolePerturbation p2,
                           double f1, double f2, double thresh,
                           bool spin_restricted)
      : type(PerturbationType::Dipole), perturbations(p1, p2),
        frequencies(f1, f2), threshold(thresh),
        spin_restricted(spin_restricted) {}

  [[nodiscard]] std::string perturbationDescription() const {
    return "VBC_" + std::string(1, perturbations.first.direction) +
           std::string(1, perturbations.second.direction);
  }

  [[nodiscard]] std::string response_filename() const {

    auto pathname = fs::path("vbc");
    if (!fs::exists(pathname)) {
      fs::create_directory(pathname);
    }

    std::ostringstream oss;

    oss << "vbc/" << perturbationDescription() << "_p" << threshold << "_f"
        << frequencies.first << "_" << frequencies.second << ".response";
    return oss.str();
  }
  [[nodiscard]] std::string response_filename(const size_t &thresh_index,
                                              const size_t &freq_index) const {
    std::ostringstream oss;
    oss << "vbc/" << perturbationDescription() << "_p" << threshold << "_f"
        << frequencies.first << "_" << frequencies.second << ".response";
    return oss.str();
  }

  [[nodiscard]] bool is_spin_restricted() const { return spin_restricted; }

  [[nodiscard]] bool is_static() const {
    return std::abs(frequencies.first) < 1e-8 &&
           std::abs(frequencies.second) < 1e-8;
  }
};

#endif // RESPONSE_STATE_HPP
