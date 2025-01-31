
/// \file GroundStateCalculation
/// \brief Input parameters for a response calculation, read from a specified
/// archive.

#ifndef MADNESS_APPS_GROUNDPARAMS_H_INCLUDED
#define MADNESS_APPS_GROUNDPARAMS_H_INCLUDED

#include <utility>

#include "../chem/molecule.h"
#include "Plot_VTK.h"
#include "basic_operators.h"
#include "madness/chem/NWChem.h"  // For nwchem interface
#include "madness/chem/SCFOperators.h"
#include "madness/chem/pointgroupsymmetry.h"
#include "madness/chem/potentialmanager.h"
#include "madness/chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "madness/mra/functypedefs.h"
#include "madness/tensor/tensor.h"

using namespace madness;

class GroundStateCalculation {
  // Ground state parameters that are read in from archive
  std::string inFile{
      "../moldft.restartdata"};  ///< Name of input archive to read in ground state
  bool spinrestricted{
      true};  ///< Indicates if ground state calc. was open or closed
  double converged_for_thresh{1.e10};
  ///< shell
  unsigned int num_orbitals{};  ///< Number of orbitals in ground state
  Tensor<double> energies{};    ///< Energy of ground state orbitals
  Tensor<double> occ{};         ///< Occupancy of ground state orbitals
  double
      L{};  ///< Box size of ground state - response calcluation is in same box
  int k{};  ///< Order of polynomial used in ground state
  Molecule molecule_in{};  ///< The molecule used in ground state calculation
  std::vector<real_function_3d> g_orbitals{};  ///< The ground state orbitals
  std::string xc{};  ///< Name of xc functional used in ground state
  std::string
      localize_method{};  ///< Name of xc functional used in ground state

  // Default constructor
 public:
  explicit GroundStateCalculation(World& world) { read(world); }

  explicit GroundStateCalculation(World& world, std::string input_file)
      : inFile{std::move(input_file)} {
    read(world);
  }

  GroundStateCalculation(const GroundStateCalculation& other) = default;

  bool is_spinrestricted() const { return spinrestricted; }

  unsigned int n_orbitals() const { return num_orbitals; }

  Tensor<double> get_energies() const { return energies; }

  Tensor<double> get_occ() const { return occ; }

  Molecule molecule() const { return molecule_in; }

  double get_L() const { return L; }

  int get_k() const { return k; }

  vector_real_function_3d& orbitals() { return g_orbitals; }

  std::string get_xc() const { return xc; }
  std::string get_localize_method() const { return localize_method; }

  std::string get_archive() const { return xc; }

  // Initializes ResponseParameters using the contents of file \c filename
  void read(World& world) {
    // Save the filename

    unsigned int dummyversion;
    double dummy1;
    std::vector<int> dummy2;

    archive::ParallelInputArchive input(world, inFile.c_str());
    input & dummyversion;
    input & dummy1;           // double
    input & spinrestricted;   // bool
    input & L;                // double            box size
    input & k;                // int               wavelet order
    input & molecule_in;      // Molecule
    input & xc;               // std:string        xc functional
    input & localize_method;  // std:string        localize  method
    if (dummyversion > 3) {
      input & converged_for_thresh;
    }
    input & num_orbitals;  // int
    input & energies;      // Tensor<double>    orbital energies
    input & occ;           // Tensor<double>    orbital occupations
    input & dummy2;        // std::vector<int>  sets of orbitals(?)

    // Check that order is positive and less than 30
    if (k < 1 or k > 30) {
      if (world.rank() == 0)
        print(
            "\n   ***PLEASE NOTE***\n   Invalid wavelet order read from "
            "archive, setting to 8.\n   This seems to happen when the default "
            "wavelet order is used in moldft.");
      k = 8;
    }
    // Set this so we can read in whats
    // written in the archive
    FunctionDefaults<3>::set_k(k);
    // Possible to call this function multiple times now
    // Do this to ensure everything works.
    world.gop.fence();
    g_orbitals.clear();
    world.gop.fence();
    // Read in ground state orbitals
    for (unsigned int i = 0; i < num_orbitals; i++) {
      real_function_3d reader;
      input & reader;
      g_orbitals.push_back(reader);
    }
    world.gop.fence();
    //projector_irrep c2v("c2v");
    //g_orbitals = c2v(g_orbitals);
    // Clean up
    truncate(world, g_orbitals);
  }

  // Prints all information
  void print_params() const {
    madness::print("\n     Ground State Parameters");
    madness::print("     -----------------------");
    madness::print("    Ground State Archive:", inFile);
    madness::print(" Ground State Functional:", xc);
    madness::print(" Localize Method  Functional:", localize_method);
    madness::print("         Spin Restricted:", spinrestricted);
    madness::print("      Number of orbitals:", num_orbitals);
    madness::print("                       L:", L);
    madness::print("           Wavelet Order:", k);
    madness::print("        Orbital Energies:", energies);
  }
};

#endif
