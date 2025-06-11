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




// Helper to convert enum to string
std::string property_type_to_string(MolecularPropertyType type) {
  switch (type) {
    case MolecularPropertyType::Polarizability:
      return "Polarizability";
    case MolecularPropertyType::Raman:
      return "Raman";
    case MolecularPropertyType::Hyperpolarizability:
      return "Hyperpolarizability";
    default:
      return "Unknown";
  }
}

// Function to print requested properties
void print_requested_properties(const std::vector<MolecularProperty> &properties) {
  std::cout << "\n📋 Requested Molecular Properties:\n";
  for (const auto &prop : properties) {
    std::cout << "- " << property_type_to_string(prop.type) << ":\n";
    for (double freq : prop.frequencies) {
      for (char dir : prop.directions) {
        std::cout << "  • Direction: " << dir << ", Frequency: "
                  << std::scientific << std::setprecision(3) << freq << "\n";
      }
    }
  }
}

