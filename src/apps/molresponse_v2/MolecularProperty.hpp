#pragma once
#include "ResponseState.hpp"
#include <string>
#include <vector>

// Enumerate clearly your supported molecular properties
enum class MolecularPropertyType { Polarizability, Raman, Hyperpolarizability };

// A molecular property request explicitly specifies frequencies and possibly
// directions
struct MolecularProperty {
  MolecularPropertyType type;
  std::vector<double> frequencies; // Frequencies explicitly required
  std::vector<char> directions;    // Directions explicitly required (optional)

  MolecularProperty(MolecularPropertyType t, std::vector<double> freqs,
                    std::vector<char> dirs = {'x', 'y', 'z'})
      : type(t), frequencies(freqs), directions(dirs) {}
};

