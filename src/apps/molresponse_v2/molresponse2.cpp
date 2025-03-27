//
// Created by ahurtado on 3/25/25.
//
// main.cpp
#include "AccuracyProtocols.hpp"
#include "MolecularProperties.hpp"
#include "ResponseManager.hpp"
#include "ResponseState.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  // Initial Setup
  std::cout << "Initializing molecular response calculation...\n";

  // Calculation-wide parameters
  CalculationSettings calcSettings(
      /*L=*/10.0,
      /*k=*/8,
      /*xc=*/"PBE",
      /*spin_restricted=*/true,
      /*localize_method=*/"boys",
      /*molecule=*/loadMoleculeFromFile("molecule.dat"));

  // Ground-state orbitals would typically be loaded from your SCF calculation
  OrbitalSet groundOrbitals("ground_state_orbitals.dat");

  // Instantiate the ResponseManager
  ResponseManager responseManager(groundOrbitals, calcSettings);

  // Define frequencies and perturbations of interest
  std::vector<double> frequencies = {0.0, 0.0656, 0.10};
  std::vector<ResponseState::Perturbation> perturbations = {
      ResponseState::Perturbation::X, ResponseState::Perturbation::Y,
      ResponseState::Perturbation::Z};

  // Define an accuracy protocol (coarse to fine thresholds)
  IncrementalAccuracy accuracyProtocol({1e-4, 1e-6, 1e-8});

  // Compute all states incrementally by accuracy
  for (auto threshold : accuracyProtocol.thresholds()) {
    std::cout << "Starting calculations at threshold: " << threshold << "\n";
    for (auto perturbation : perturbations) {
      for (auto freq : frequencies) {
        ResponseState state{perturbation, freq, threshold, false};

        // Check convergence first, then compute
        if (!responseManager.isConverged(state)) {
          responseManager.computeState(state);
          std::cout << "Computed state: " << state.getIdentifier() << "\n";
        } else {
          std::cout << "State already converged: " << state.getIdentifier()
                    << "\n";
        }
      }
    }
  }

  // Example calculation: Compute Polarizability at freq = 0.0656
  std::cout << "Computing polarizability...\n";
  Polarizability polar(0.0656);
  polar.compute(responseManager);

  // Example calculation: Compute Hyperpolarizability at freq1 = 0.0656, freq2 =
  // 0.10
  std::cout << "Computing hyperpolarizability...\n";
  Hyperpolarizability hyperpolar(0.0656, 0.10);
  hyperpolar.compute(responseManager);

  std::cout << "Molecular response calculations complete!\n";

  return 0;
}
