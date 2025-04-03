#pragma once
#include <string>
#include <variant>

enum class PerturbationType {
  Dipole,
  NuclearDisplacement,
  Magnetic,
};

struct DipolePerturbation {
  char direction; // X, Y, or Z
};

struct NuclearDisplacementPerturbation {
  int atom_index;
  char direction; // X, Y, or Z
};

struct MagneticPerturbation {
  char direction; // X, Y, or Z
};


using Perturbation = std::variant<DipolePerturbation,
                                   NuclearDisplacementPerturbation,
                                   MagneticPerturbation>; 
