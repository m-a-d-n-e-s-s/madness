#ifndef RESPONSEMANAGER_HPP
#define RESPONSEMANAGER_HPP
#include <CalculationParameters.h>
#include <madness/chem/SCF.h>

#include "ResponseState.hpp"

using namespace madness;

class ResponseManager {
 public:
  ResponseManager(World &world, const CalculationParameters &r_params);

  void setProtocol(World &world, double L, double thresh, int override_k = -1);

  // Getters
  [[nodiscard]] double getVtol() const { return vtol; }
  [[nodiscard]] poperatorT getCoulombOp() const { return coulop; }
  [[nodiscard]] std::vector<std::shared_ptr<real_derivative_3d>> getGradOp() const { return gradop; }
  // Compute and store a response state
  void computeState(const AbstractResponseDescriptor &state);
  // Check if a response state is already converged
  [[nodiscard]] bool isConverged(const AbstractResponseDescriptor &state) const;

  [[nodiscard]] CalculationParameters params() const { return calc_params; }

 private:
  CalculationParameters calc_params;
  double vtol;
  poperatorT coulop;
  std::vector<std::shared_ptr<real_derivative_3d>> gradop;
};

#endif  // RESPONSEMANAGER_HPP
