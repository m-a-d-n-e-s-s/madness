#ifndef RESPONSE_PRELIMINARIES_CALCULATOR_HPP
#define RESPONSE_PRELIMINARIES_CALCULATOR_HPP
#include "GroundStateData.hpp"
#include "ResponsePreliminaries.hpp"
#include "functypedefs.h"
#include "madness/chem/SCF.h"
#include "molresponse_v2/GroundStateData.hpp"
#include <madness/world/world.h>

using namespace madness;

// This class computes all the reusuable data from  the ground state
// calculation that is needed for each iteration of any response calculation.
// Call this function at each change of protocol
class ResponsePreliminariesCalculator {
public:
  ResponsePreliminariesCalculator(World &world,
                                  const GroundStateData &ground_state,
                                  double coulomb_lo_thresh,
                                  double potential_thresh);

  Tensor<double> computeKineticEnergy();
  real_function_3d computeDensity();
  real_function_3d computeNuclearPotential();
  real_function_3d computeCoulombPotential();
  vector_real_function_3d computeHFExchangeEnergy();
  real_function_3d computeXCPotential();
  Tensor<double> computeHamiltonian();

  ResponsePreliminaries computeAll();

private:
  World &world_;
  const GroundStateData &ground_state_;
  XCfunctional xcf_;

  std::shared_ptr<PotentialManager> potential_manager_;
  std::shared_ptr<operatorT> coulomb_operator_;

  double coulomb_lo_thresh_;
  double potential_thresh_;
  real_function_3d density_;
};
#endif // RESPONSE_PRELIMINARIES_CALCULATOR_HPP
