#ifndef GROUNDSTATEDATA_HPP
#define GROUNDSTATEDATA_HPP

#include "potentialmanager.h"
#include <madness/chem/SCF.h>
#include "broadcast_json.hpp"
#include <madness/chem/molecule.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>
#include <madness/chem/projector.h>
#include <vector>

using namespace madness;

class GroundStateData {

private:
  std::string archive_File;
  int original_k;
  bool spinrestricted;
  unsigned int num_orbitals{};
  Tensor<double> energies;
  Tensor<double> occ;
  double L;
  int k;

public:
  Molecule molecule;
  std::vector<real_function_3d> orbitals;
  QProjector<double, 3> Qhat;
  std::string xc;
  std::string localize_method;
  double converged_for_thresh;
  bool fock_loaded_from_file = false;
  [[nodiscard]] bool isFockLoadedFromFile() const { return fock_loaded_from_file; }

  [[nodiscard]] const std::vector<real_function_3d> &getOrbitals() const { return orbitals; }
  [[nodiscard]] const Tensor<double> &getEnergies() const { return energies; }
  [[nodiscard]] const Tensor<double> &getOcc() const { return occ; }
  [[nodiscard]] long getNumOrbitals() const { return num_orbitals; }
  [[nodiscard]] double getL() const { return L; }
  [[nodiscard]] int getK() const { return k; }
  [[nodiscard]] bool isSpinRestricted() const { return spinrestricted; }
  [[nodiscard]] const Molecule &getMolecule() const { return molecule; }
  // The following
  XCfunctional xcf_;
  std::shared_ptr<PotentialManager> potential_manager_;
  // The following variables are used for the response calculation
  real_function_3d density_;
  real_function_3d V_local;
  Tensor<double> Hamiltonian;
  Tensor<double> Hamiltonian_no_diag;

  explicit GroundStateData(World &world, std::string archiveFile, Molecule mol);
  explicit GroundStateData(World &world, std::shared_ptr<SCF> scf_calc);
  void load(World &world);
  void prepareOrbitals(World &world, int current_k, double thresh);
  void print_info() const;
  // These functions are used to compute the preliminaries
  real_function_3d computeDensity(World &world) const;
  [[nodiscard]] real_function_3d computeNuclearPotential(double vtol) const;
  [[nodiscard]] real_function_3d computeCoulombPotential(const operatorT &coulop, double vtol) const;
  void computePreliminaries(World &world, const operatorT &coulop, double vtol, const std::string &fock_json_file);
  [[nodiscard]] real_function_3d computeXCPotential(World &world) const;
  [[nodiscard]] Tensor<double> computeKineticEnergy(World &world) const;
  [[nodiscard]] vector_real_function_3d computeHFExchangeEnergy(World &world) const;
  Tensor<double> tryLoadHamiltonianFromJson(World &world, const json &fock_json, double thresh, int k);
};

#endif // GROUNDSTATEDATA_HPP
