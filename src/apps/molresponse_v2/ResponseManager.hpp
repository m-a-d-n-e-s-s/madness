#ifndef RESPONSEMANAGER_HPP
#define RESPONSEMANAGER_HPP
#include "GroundStateData.hpp"
#include "ResponseState.hpp"
#include "molresponse_v2/ResponsePreliminaries.hpp"
#include <madness/mra/funcdefaults.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>

#include <vector>

using namespace madness;

class ResponseManager {
public:
  ResponseManager(World &world, const std::string &ground_state_archive,
                  const Molecule &molecule);
  void setProtocol(double thresh, int override_k = -1);

  // Prepares the orbitals for a given polynomial order (k) at an accuracy step
  void prepareOrbitalsForAccuracyStep();

  void computePreliminaries();

  [[nodiscard]] const ResponsePreliminaries currentPreliminaries() const;

  // Access the current orbitals
  [[nodiscard]] const GroundStateData currentGroundState() const;

  // Compute and store a response state
  void computeState(const ResponseState &state);

  // Check if a response state is already converged
  bool isConverged(const ResponseState &state) const;

private:
  World &world_;
  std::string ground_state_archive_;

  GroundStateData ground_state_;
  ResponsePreliminaries
      response_preliminaries_; // Resuable data for response calculations
};

#endif // RESPONSEMANAGER_HPP
