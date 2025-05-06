#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include <SCF.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>

#include <madness/external/nlohmann_json/json.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "GroundStateData.hpp"
#include "Perturbation.hpp"
#include "ResponseVector.hpp"
#include "molecular_functors.h"
#include "vmra.h"

using json = nlohmann::json;
namespace fs = std::filesystem;

struct AbstractResponseDescriptor {
  [[nodiscard]] virtual bool is_spin_restricted() const = 0;
  [[nodiscard]] virtual bool is_static() const = 0;
  [[nodiscard]] virtual std::string response_filename() const = 0;
  [[nodiscard]] virtual std::string response_filename(const size_t &thresh_index, const size_t &freq_index) const = 0;
  virtual ~AbstractResponseDescriptor() = default;
};

struct ResponseState : public AbstractResponseDescriptor {
  PerturbationType type;
  Perturbation perturbation;

  bool spin_restricted = false;  // Is the system open shell?

  std::vector<double> frequencies;
  std::map<double, int> frequency_map;  // Frequency to index map
  std::vector<double> thresholds;       // Accuracy levels to loop over
  //
  size_t current_frequency_index = 0;
  size_t current_thresh_index = 0;  // Track which threshold we're working on
  bool is_converged = false;
  ResponseState() = default;

  ResponseState(Perturbation pert, PerturbationType ptype, const std::vector<double> &freq, const std::vector<double> &thresh, bool spin_restricted)
      : type(ptype),
        perturbation(pert),
        frequencies(freq),
        thresholds(thresh),
        current_frequency_index(0),
        current_thresh_index(0),
        spin_restricted(spin_restricted),
        is_converged(false) {
    for (size_t i = 0; i < frequencies.size(); ++i) {
      frequency_map[frequencies[i]] = i;
    }
  }

  [[nodiscard]] double current_threshold() const { return thresholds[current_thresh_index]; }
  [[nodiscard]] double current_frequency() const { return frequencies[current_frequency_index]; }

  [[nodiscard]] bool at_final_frequency() const { return current_frequency_index == frequencies.size() - 1; }
  [[nodiscard]] bool at_final_threshold() const { return current_thresh_index == thresholds.size() - 1; }

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
      is_converged = false;  // reset convergence at new protocol
    }
  }
  void advance_frequency() {
    if (!at_final_frequency()) {
      ++current_frequency_index;
    }
  }
  [[nodiscard]] bool is_static() const override { return std::abs(current_frequency()) < 1e-8; }
  [[nodiscard]] bool is_dynamic() const { return !is_static(); }
  [[nodiscard]] bool is_spin_restricted() const override { return spin_restricted; }

  [[nodiscard]] std::string response_filename() const override {
    std::ostringstream oss;
    oss << "responses/" << describe_perturbation(perturbation) << "_p" << current_threshold() << "_f" << current_frequency() << ".response";
    return oss.str();
  }
  [[nodiscard]] std::string response_filename(const size_t &thresh_index, const size_t &freq_index) const override {
    std::ostringstream oss;
    oss << "responses/" << describe_perturbation(perturbation) << "_p" << thresholds[thresh_index] << "_f" << frequencies[freq_index] << ".response";
    return oss.str();
  }

  [[nodiscard]] std::string perturbationDescription() const { return describe_perturbation(perturbation); }

  [[nodiscard]] std::string description() const {
    std::ostringstream oss;
    oss << describe_perturbation(perturbation) << " at freq " << current_frequency() << " (thresh=" << std::scientific << current_threshold() << ")";
    return oss.str();
  }
};

struct SecondOrderResponseDescriptor : public AbstractResponseDescriptor {
  std::pair<PerturbationType, PerturbationType> ptypes_;
  std::pair<Perturbation, Perturbation> perturbations_;
  std::pair<double, double> frequencies_;

  double thresh;
  bool spin_restricted_ = false;  // Is the system open shell?

  SecondOrderResponseDescriptor(PerturbationType t1, PerturbationType t2, Perturbation p1, Perturbation p2, double f1, double f2, double thresh,
                                bool spin_restricted)
      : ptypes_(t1, t2), perturbations_(p1, p2), frequencies_(f1, f2), thresh(thresh), spin_restricted_(spin_restricted) {}

  [[nodiscard]] ResponseState B_state() const { return ResponseState(perturbations_.first, ptypes_.first, {frequencies_.first}, {thresh}, spin_restricted_); }
  [[nodiscard]] ResponseState C_state() const {
    return ResponseState(perturbations_.second, ptypes_.second, {frequencies_.second}, {thresh}, spin_restricted_);
  }

  [[nodiscard]] std::pair<ResponseState, ResponseState> get_states() const { return {B_state(), C_state()}; }

  [[nodiscard]] double current_threshold() const { return thresh; }

  [[nodiscard]] double current_frequency() const { return frequencies_.first + frequencies_.second; }

  [[nodiscard]] virtual const char *prefix() const = 0;

  [[nodiscard]] std::string perturbationDescription() const {
    return describe_perturbation(perturbations_.first) + "_" + describe_perturbation(perturbations_.second);
  }

  [[nodiscard]] std::string response_filename() const override {
    std::ostringstream oss;

    oss << "responses/" << prefix() << perturbationDescription() << "_p" << std::scientific << thresh << "_f" << frequencies_.first << "_"
        << frequencies_.second << ".response";
    return oss.str();
  }
  [[nodiscard]] std::string response_filename(const size_t &thresh_index, const size_t &freq_index) const override { return response_filename(); }

  [[nodiscard]] bool is_spin_restricted() const override { return spin_restricted_; }

  [[nodiscard]] bool is_static() const override {
    return false;
    // Second order response is always dynamic (x an y response
    // functions)
  }
};

//-----------------------------------------------------------------------------
// Now two trivial subclasses for VBC vs. XBC
//-----------------------------------------------------------------------------

struct VBCResponseState : public SecondOrderResponseDescriptor {
  using SecondOrderResponseDescriptor::SecondOrderResponseDescriptor;
  [[nodiscard]] const char *prefix() const override { return "VBC"; }
};

struct XBCResponseState : public SecondOrderResponseDescriptor {
  using SecondOrderResponseDescriptor::SecondOrderResponseDescriptor;
  [[nodiscard]] const char *prefix() const override { return "XBC"; }
};

// -----------------------------------------------------------------------------
// 1) Raw operator in real space, before applying to orbitals:
//
//    e.g. for a dipole:  V(r) = x, y or z moment
// -----------------------------------------------------------------------------
inline real_function_3d raw_perturbation_operator(World &world, const GroundStateData &gs, const ResponseState &state) {
  using P = PerturbationType;
  switch (state.type) {
    case P::Dipole: {
      auto d = std::get<DipolePerturbation>(state.perturbation);
      // build the moment functor f = (1,0,0) or (0,1,0) or (0,0,1)
      std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};
      std::vector<int> dir(3, 0);
      dir[dipole_map.at(d.direction)] = 1;
      real_function_3d f = real_factory_3d(world).functor(real_functor_3d{new MomentFunctor(dir)});
      f.truncate(FunctionDefaults<3>::get_thresh());
      return f;
    }
    case P::NuclearDisplacement: {
      auto n = std::get<NuclearDisplacementPerturbation>(state.perturbation);
      // you’d have whatever MomentDisplacementFunctor exists:
      real_function_3d f = real_factory_3d(world).functor(real_functor_3d{new madchem::MolecularDerivativeFunctor(gs.molecule, n.atom_index, n.direction)});
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
inline vector_real_function_3d project_perturbation_onto_orbitals(World &world, const GroundStateData &gs, const real_function_3d &raw_op) {
  auto vp = mul(world, raw_op, gs.orbitals, /*fence=*/true);
  vp = gs.Qhat(vp);
  truncate(world, vp, FunctionDefaults<3>::get_thresh(), /*fence=*/true);
  return vp;
}

// Linear (first-order) response
inline madness::vector_real_function_3d perturbation_vector(madness::World &world, GroundStateData const &gs, ResponseState const &state) {
  auto raw_op = raw_perturbation_operator(world, gs, state);
  auto Vp = project_perturbation_onto_orbitals(world, gs, raw_op);
  if (state.is_dynamic()) {
    // duplicate for dynamic response
    Vp.insert(Vp.end(), Vp.begin(), Vp.end());
  }
  return Vp;
}

// Second-order (VBC) response
inline madness::vector_real_function_3d perturbation_vector(madness::World &world, GroundStateData const &gs, XBCResponseState const &sos) {
  // build (or reuse) the VBCComputer for this ground state
  // you’ll need to pass it the same directions & frequency list
  // that you used to set up your ResponseStates originally:
  /*static thread_local VBCComputer2 vbc(*/
  /**/
  /*// find the indices of sos.frequencies in that computer’s list:*/
  /*size_t bi = vbc.frequency_index(sos.frequencies.first);*/
  /*size_t ci = vbc.frequency_index(sos.frequencies.second);*/
  /*// and the BC-pair index from its perturbation characters:*/
  /*size_t bc = vbc.BC_pair_index(sos.perturbations.first.direction,*/
  /*                              sos.perturbations.second.direction);*/
  /**/
  /*// will load from disk if already there, otherwise compute & save:*/
  /*return get_flat(vbc.compute_and_save(bc, bi, ci));*/
  return {};
}

#endif  // RESPONSE_STATE_HPP
