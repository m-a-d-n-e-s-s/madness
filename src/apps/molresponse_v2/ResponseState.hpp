#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include "GroundStateData.hpp"
#include "Perturbation.hpp"
#include "molecular_functors.h"
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
    oss << describe_perturbation(pert) << " at freq " << current_frequency()
        << " (thresh=" << std::scientific << current_threshold() << ")";
    return oss.str();
  }

};

struct SecondOrderResponseState : public AbstractResponseDescriptor {
  std::pair<Perturbation, Perturbation> perturbations;
  std::pair<double, double> frequencies;
  double threshold;
  bool spin_restricted = true;

  SecondOrderResponseState(Perturbation p1, Perturbation p2, double f1,
                           double f2, double thresh, bool spin_restricted)
      : perturbations(p1, p2), frequencies(f1, f2), threshold(thresh),
        spin_restricted(spin_restricted) {}

  [[nodiscard]] std::string perturbationDescription() const {

    return describe_perturbation(perturbations.first) + "_" +
           describe_perturbation(perturbations.second);
  }

  [[nodiscard]] std::string response_filename() const {

    std::ostringstream oss;

    oss << "responses/VBC_" << perturbationDescription() << "_p"
        << std::scientific << threshold << "_f" << frequencies.first << "_"
        << frequencies.second << ".response";
    return oss.str();
  }
  [[nodiscard]] std::string response_filename(const size_t &thresh_index,
                                              const size_t &freq_index) const {
    return response_filename();
  }

  [[nodiscard]] bool is_spin_restricted() const { return spin_restricted; }

  [[nodiscard]] bool is_static() const {
    return false; // Second order response is always dynamic (x an y response
                  // functions)
  }
};

// -----------------------------------------------------------------------------
// 1) Raw operator in real space, before applying to orbitals:
//
//    e.g. for a dipole:  V(r) = x, y or z moment
// -----------------------------------------------------------------------------
inline real_function_3d raw_perturbation_operator(World &world,
                                                  const GroundStateData &gs,
                                                  const ResponseState &state) {
  using P = PerturbationType;
  switch (state.type) {
  case P::Dipole: {
    auto d = std::get<DipolePerturbation>(state.perturbation);
    // build the moment functor f = (1,0,0) or (0,1,0) or (0,0,1)
    std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};
    std::vector<int> dir(3, 0);
    dir[dipole_map.at(d.direction)] = 1;
    real_function_3d f =
        real_factory_3d(world).functor(real_functor_3d{new MomentFunctor(dir)});
    f.truncate(FunctionDefaults<3>::get_thresh());
    return f;
  }
  case P::NuclearDisplacement: {
    auto n = std::get<NuclearDisplacementPerturbation>(state.perturbation);
    // you’d have whatever MomentDisplacementFunctor exists:
    real_function_3d f = real_factory_3d(world).functor(
        real_functor_3d{new madchem::MolecularDerivativeFunctor(
            gs.molecule, n.atom_index, n.direction)});
    f.truncate(FunctionDefaults<3>::get_thresh());
    return f;
  }
  case P::Magnetic: {
    // Not implemented yet...
    //
    //
    throw std::runtime_error("Magnetic perturbation not implemented yet");
  }
  }
  throw std::runtime_error("Unknown perturbation type");
}

// -----------------------------------------------------------------------------
// 2) Apply it to the ground‐state orbitals to get your Vp basis functions.
//    (You already have this in ResponseState::perturbation_vector.)
// -----------------------------------------------------------------------------
inline vector_real_function_3d
project_perturbation_onto_orbitals(World &world, const GroundStateData &gs,
                                   const real_function_3d &raw_op) {
  auto vp = mul(world, raw_op, gs.orbitals, /*fence=*/true);
  vp = gs.Qhat(vp);
  truncate(world, vp, FunctionDefaults<3>::get_thresh(), /*fence=*/true);
  return vp;
}

#endif // RESPONSE_STATE_HPP
