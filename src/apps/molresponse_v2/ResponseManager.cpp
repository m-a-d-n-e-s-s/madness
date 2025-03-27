#ifndef RESPONSEMANAGER_HPP
#include "ResponseManager.hpp"
#include "ResponsePreliminariesCalculator.hpp"
#include "molresponse_v2/ResponsePreliminaries.hpp"

ResponseManager::ResponseManager(World &world,
                                 const std::string &ground_state_archive,
                                 const Molecule &mol)
    : world_(world), ground_state_archive_(ground_state_archive),
      ground_state_(world, ground_state_archive, mol) {
  // Initially loaded at maximum polynomial order (settings_.k)
}

void ResponseManager::setProtocol(double thresh, int override_k) {
  int k;
  if (thresh >= 0.9e-2)
    k = 4;
  else if (thresh >= 0.9e-4)
    k = 6;
  else if (thresh >= 0.9e-6)
    k = 8;
  else if (thresh >= 0.9e-8)
    k = 10;
  else
    k = 12;

  if (override_k != -1) {
    k = override_k;
  }

  if (override_k == -1) {
    FunctionDefaults<3>::set_k(k);
    //        	param.k=k;
  } else {
    FunctionDefaults<3>::set_k(override_k);
  }
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_apply_randomize(false);
  FunctionDefaults<3>::set_project_randomize(false);

  if (world_.rank() == 0) {
    print("\nMRA Protocol Set: thresh =", thresh, ", k =", k, "\n");
  }
}

void ResponseManager::prepareOrbitalsForAccuracyStep() {

  auto thresh = FunctionDefaults<3>::get_thresh();
  auto current_k = FunctionDefaults<3>::get_k();

  bool needs_reload = false;

  // Check if the current orbitals' polynomial order is higher (projectable)
  if (ground_state_.original_k < current_k) {
    if (world_.rank() == 0) {
      print("Cannot project orbitals: current k (", ground_state_.k,
            ") is lower than requested k (", current_k, ").");
    }
    MADNESS_EXCEPTION("Orbital polynomial order too low for projection.",
                      ground_state_.k);
  }

  // Check if orbitals need reloading due to polynomial order mismatch or
  // convergence threshold
  if (ground_state_.k != current_k ||
      ground_state_.converged_for_thresh > thresh) {
    needs_reload = true;
  }

  if (needs_reload) {
    if (world_.rank() == 0)
      print("Preparing orbitals: reloading from archive.");

    // Reload original orbitals from archive (assume these have highest accuracy
    // and polynomial order)
    ground_state_.load(world_, ground_state_archive_);

    if (world_.rank() == 0)
      print("Preparing orbitals: reconstructing.");

    reconstruct(world_, ground_state_.orbitals);

    if (world_.rank() == 0)
      print("Preparing orbitals: projecting to requested polynomial order and "
            "threshold.");

    for (auto &orbital : ground_state_.orbitals) {
      orbital = project(orbital, current_k, thresh, false);
    }
    world_.gop.fence();

    if (world_.rank() == 0)
      print("Preparing orbitals: truncating projected orbitals.");

    truncate(world_, ground_state_.orbitals, thresh);

    // Update orbital metadata
    ground_state_.k = current_k;
    ground_state_.converged_for_thresh = thresh;
  }

  if (world_.rank() == 0)
    print("Orbital preparation complete for k =", current_k,
          " thresh =", thresh);
}

void ResponseManager::computePreliminaries() {
  // Use previously set FunctionDefaults (k and thresh)
  double thresh = FunctionDefaults<3>::get_thresh();
  double coulomb_lo_thresh = thresh; // or adjust as per your logic
  double potential_thresh = thresh;

  ResponsePreliminariesCalculator calculator(
      world_, ground_state_, coulomb_lo_thresh, potential_thresh);

  response_preliminaries_ = calculator.computeAll();

  if (world_.rank() == 0) {
    print("Response preliminaries computed successfully at thresh =", thresh);
  }
}

[[nodiscard]] const ResponsePreliminaries
ResponseManager::currentPreliminaries() const {
  return response_preliminaries_;
}

[[nodiscard]] const GroundStateData
ResponseManager::currentGroundState() const {
  return ground_state_;
}

#endif // RESPONSEMANAGER_HPP
