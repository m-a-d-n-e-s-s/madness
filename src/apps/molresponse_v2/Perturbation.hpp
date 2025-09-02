#pragma once
#include <string>
#include <variant>

template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct DipolePerturbation {
  char direction; // X, Y, or Z
  //
};

struct NuclearDisplacementPerturbation {
  int atom_index;
  char direction; // X, Y, or Z
};

struct MagneticPerturbation {
  char direction; // X, Y, or Z
};

using Perturbation = std::variant<DipolePerturbation, NuclearDisplacementPerturbation, MagneticPerturbation>;

// helper to stringify any single perturbation:
inline std::string describe_perturbation(const Perturbation &p) {

  auto describe = overloaded{
      [](const DipolePerturbation &d) { return std::string("Dipole_") + d.direction; },
      [](const NuclearDisplacementPerturbation &n) {
        return "NucA_" + std::to_string(n.atom_index) + std::string(1, n.direction);
      },
      [](const MagneticPerturbation &m) { return "Mag_" + std::string(1, m.direction); },
  };
  return std::visit(describe, p);
}

inline std::string perturbation_type_string(const Perturbation &p) {
  return std::visit(overloaded{[](const DipolePerturbation &) { return "Dipole"; },
                               [](const NuclearDisplacementPerturbation &) { return "Nuclear"; },
                               [](const MagneticPerturbation &) { return "Magnetic"; }},
                    p);
}
