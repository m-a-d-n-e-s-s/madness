#ifndef ORBITALSET_HPP
#define ORBITALSET_HPP
#include "ResponseState.hpp"
// #include "madness/mra/functypedefs.h"
#include <madness/chem/molecule.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>
#include <vector>

using namespace madness;

class GroundStateData {
public:
  int original_k;
  bool spinrestricted;
  unsigned int num_orbitals{};
  Tensor<double> energies;
  Tensor<double> occ;
  double L;
  int k;
  Molecule molecule;
  std::vector<real_function_3d> orbitals;
  std::string xc;
  std::string localize_method;
  double converged_for_thresh;

  explicit GroundStateData(World &world, const std::string &archiveFile,
                           const Molecule &mol)
      : molecule(mol) {
    load(world, archiveFile);
  }

  void load(World &world, const std::string &archiveFile) {
    archive::ParallelInputArchive input(world, archiveFile.c_str());

    unsigned int version;
    double dummy_energy;
    std::vector<int> dummy_sets;

    input & version;
    input & dummy_energy;
    input & spinrestricted;
    input & L;
    input & k;
    input & molecule;
    input & xc;
    input & localize_method;

    if (version > 3) {
      input & converged_for_thresh;
    }

    input & num_orbitals;
    input & energies;
    input & occ;
    input & dummy_sets;

    if (k < 1 || k > 30) {
      if (world.rank() == 0)
        madness::print(
            "Invalid wavelet order read from archive, setting to default k=8.");
      k = 8;
    }

    // FunctionDefaults<3>::set_k(k);

    world.gop.fence();
    orbitals.clear();
    world.gop.fence();

    for (unsigned int i = 0; i < num_orbitals; ++i) {
      real_function_3d orbital;
      input & orbital;
      orbitals.push_back(orbital);
    }

    world.gop.fence();
    truncate(world, orbitals);

    if (world.rank() == 0) {
      print_info();
    }
    original_k = k;
  }

  void print_info() const {
    madness::print("\nGround State Orbital Information:");
    madness::print("------------------------");
    madness::print("XC Functional:", xc);
    madness::print("Localization Method:", localize_method);
    madness::print("Spin Restricted:", spinrestricted);
    madness::print("Number of Orbitals:", num_orbitals);
    madness::print("Box Size L:", L);
    madness::print("Wavelet Order k:", k);
    madness::print("Converged for Threshold:", converged_for_thresh);
    madness::print("Orbital Energies:", energies);
  }
};

#endif // ORBITALSET_HPP
