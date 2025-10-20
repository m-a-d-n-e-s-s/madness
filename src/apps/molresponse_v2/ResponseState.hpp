#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include <SCF.h>
#include <madness/chem/projector.h>
#include <madness/mra/mra.h>

#include "ResponseVector.hpp"
#include <madness/external/nlohmann_json/json.hpp>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "GroundStateData.hpp"
#include "Perturbation.hpp"
#include "molecular_functors.h"
#include "vmra.h"

using json = nlohmann::json;
namespace fs = std::filesystem;

struct AbstractResponseDescriptor {
  [[nodiscard]] virtual bool is_spin_restricted() const = 0;

  [[nodiscard]] virtual bool is_static() const = 0;

  [[nodiscard]] virtual std::string response_filename() const = 0;

  [[nodiscard]] virtual std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const = 0;

  virtual ~AbstractResponseDescriptor() = default;
};

struct LinearResponseDescriptor : public AbstractResponseDescriptor {
  Perturbation perturbation;

  bool spin_restricted = false; // Is the system open shell?

  std::vector<double> frequencies;
  std::map<double, size_t> frequency_map; // Frequency to index map
  std::vector<double> thresholds;         // Accuracy levels to loop over
  //
  size_t current_frequency_index = 0;
  size_t current_thresh_index = 0; // Track which threshold we're working on
  bool is_converged = false;

  LinearResponseDescriptor() = default;

  // default copy constructor

  LinearResponseDescriptor(Perturbation pert, const std::vector<double> &freq,
                           const std::vector<double> &thresh,
                           bool spin_restricted)
      : perturbation(pert), frequencies(freq), thresholds(thresh),
        spin_restricted(spin_restricted) {
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
    return current_frequency_index + 1 == frequencies.size();
  }
  [[nodiscard]] bool at_final_threshold() const {
    return current_thresh_index + 1 == thresholds.size();
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

  [[nodiscard]] bool is_static() const override {
    return std::abs(current_frequency()) < 1e-8;
  }
  [[nodiscard]] bool is_static(size_t freq_index) const {
    return std::abs(frequencies[freq_index]) < 1e-8;
  }
  [[nodiscard]] bool is_spin_restricted() const override {
    return spin_restricted;
  }

  // helper that builds the core "<perturbation>_p<thresh>_f<freq>"
  [[nodiscard]] std::string make_key(double thresh, double freq) const {
    std::ostringstream oss;
    // describe_perturbation gives "alpha", "beta", etc.
    double f = std::clamp(freq, 0.0, 100.0);

    oss << describe_perturbation(perturbation)
        // force exactly 3 decimal places:
        << "_f" << std::fixed << std::setprecision(3)
        << f
        // then threshold in scientific (1 digit mantissa)
        << "_p" << std::scientific << std::setprecision(2) << thresh;

    return oss.str();
  }

  [[nodiscard]] std::string make_key(size_t ti, size_t fi) const {
    return make_key(thresholds[ti], frequencies[fi]);
  }

  [[nodiscard]] std::string response_filename() const override {
    return response_filename(current_thresh_index, current_frequency_index);
  }

  [[nodiscard]] std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const override {
    return make_key(thresh_index, freq_index);
    ;
  }

  [[nodiscard]] std::string perturbationDescription() const {
    return describe_perturbation(perturbation);
  }
  [[nodiscard]] std::string description() const {
    return make_key(current_threshold(), current_frequency());
  }
};

struct SecondOrderResponseDescriptor : public AbstractResponseDescriptor {
  std::pair<Perturbation, Perturbation> perturbations_;
  std::pair<double, double> frequencies_;

  double thresh;
  bool spin_restricted_ = false; // Is the system open shell?

  SecondOrderResponseDescriptor(Perturbation p1, Perturbation p2, double f1,
                                double f2, double thresh, bool spin_restricted)
      : perturbations_(p1, p2), frequencies_(f1, f2), thresh(thresh),
        spin_restricted_(spin_restricted) {}

  [[nodiscard]] LinearResponseDescriptor B_state() const {
    return LinearResponseDescriptor(perturbations_.first, {frequencies_.first},
                                    {thresh}, spin_restricted_);
  }

  [[nodiscard]] LinearResponseDescriptor C_state() const {
    return LinearResponseDescriptor(perturbations_.second,
                                    {frequencies_.second}, {thresh},
                                    spin_restricted_);
  }

  [[nodiscard]] std::pair<LinearResponseDescriptor, LinearResponseDescriptor>
  get_states() const {
    return {B_state(), C_state()};
  }

  [[nodiscard]] double current_threshold() const { return thresh; }

  [[nodiscard]] double current_frequency() const {
    return frequencies_.first + frequencies_.second;
  }

  [[nodiscard]] virtual const char *prefix() const = 0;

  [[nodiscard]] std::string perturbationDescription() const {
    return describe_perturbation(perturbations_.first) + "_" +
           describe_perturbation(perturbations_.second);
  }

  [[nodiscard]] bool is_spin_restricted() const override {
    return spin_restricted_;
  }

  [[nodiscard]] bool is_static() const override {
    return false;
    // Second order response is always dynamic (x an y response
    // functions)
  }
  [[nodiscard]] bool is_static(size_t freq_index) const { return false; }

  // Build the core "<prefix><pert1>_<pert2>_f<f1>_<f2>_p<thresh>"
  [[nodiscard]] std::string make_key(double f1, double f2,
                                     double thresh) const {
    std::ostringstream oss;

    // clamp both frequencies
    f1 = std::clamp(f1, 0.0, 100.0);
    f2 = std::clamp(f2, 0.0, 100.0);

    oss
        // your existing prefix (e.g. "sr_" or "dyn_")
        << prefix()
        // two perturbations joined by "_"
        << describe_perturbation(perturbations_.first) << "_"
        << describe_perturbation(perturbations_.second)
        // f1 with fixed + 3 decimals
        << "_f" << std::fixed << std::setprecision(3)
        << f1
        // underscore + f2 in the same format
        << "_" << std::fixed << std::setprecision(3)
        << f2
        // threshold in scientific (1 digit mantissa)
        << "_p" << std::scientific << std::setprecision(1) << thresh;

    return oss.str();
  }

  [[nodiscard]] std::string response_filename() const override {
    return make_key(frequencies_.first, frequencies_.second, thresh) +
           ".response";
  }

  [[nodiscard]] std::string
  response_filename(const size_t &thresh_index,
                    const size_t &freq_index) const override {
    return response_filename();
  }

  [[nodiscard]] ResponseVector make_vector(int num_orbitals, size_t fi) const {
    bool stat = is_static();
    bool urstr = !spin_restricted_;
    return make_response_vector(num_orbitals, stat, urstr);
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

inline int dir_index(char c) {
  switch (std::tolower(static_cast<unsigned char>(c))) {
  case 'x':
    return 0;
  case 'y':
    return 1;
  case 'z':
    return 2;
  }
  throw std::invalid_argument("Direction must be x/y/z");
}

inline real_function_3d
make_perturbation_operator(World &world, const GroundStateData &g_s,
                           const DipolePerturbation &d) {
  std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};
  std::vector<int> dir(3, 0);
  dir[dipole_map.at(d.direction)] = 1;
  real_function_3d f =
      real_factory_3d(world).functor(real_functor_3d{new MomentFunctor(dir)});
  f.truncate(FunctionDefaults<3>::get_thresh());
  return f;
}

inline real_function_3d
make_perturbation_operator(World &world, const GroundStateData &gs,
                           const NuclearDisplacementPerturbation &n) {
  std::map<char, int> dipole_map = {{'x', 0}, {'y', 1}, {'z', 2}};

  madchem::MolecularDerivativeFunctor mdfunctor(gs.molecule, n.atom_index,
                                                dipole_map.at(n.direction));
  // you’d have whatever MomentDisplacementFunctor exists:
  real_function_3d dvdx = real_factory_3d(world)
                              .functor(mdfunctor)
                              .nofence()
                              .truncate_on_project()
                              .truncate_mode(0);
  dvdx.truncate();
  return dvdx;
}

inline real_function_3d
make_perturbation_operator(World &world, const GroundStateData &gs,
                           const MagneticPerturbation &m) {
  // Not implemented yet...
  //
  //
  throw std::runtime_error("Magnetic perturbation not implemented yet");
}

inline real_function_3d raw_perturbation_operator(World &world,
                                                  const GroundStateData &gs,
                                                  const Perturbation &p) {
  return std::visit(
      [&](const auto &pp) -> real_function_3d {
        return make_perturbation_operator(world, gs, pp);
      },
      p);
}

// -----------------------------------------------------------------------------
// 1) Raw operator in real space, before applying to orbitals:
//
//    e.g. for a dipole:  V(r) = x, y or z moment
// -----------------------------------------------------------------------------

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

// Linear (first-order) response
inline madness::vector_real_function_3d
perturbation_vector(madness::World &world, GroundStateData const &gs,
                    LinearResponseDescriptor const &state) {
  auto raw_op = raw_perturbation_operator(world, gs, state.perturbation);
  auto Vp = project_perturbation_onto_orbitals(world, gs, raw_op);
  if (!state.is_static()) {
    // duplicate for dynamic response
    Vp.insert(Vp.end(), Vp.begin(), Vp.end());
  }
  return Vp;
}

// Second-order (VBC) response
inline madness::vector_real_function_3d
perturbation_vector(madness::World &world, GroundStateData const &gs,
                    XBCResponseState const &sos) {
  return {};
}

#endif // RESPONSE_STATE_HPP
