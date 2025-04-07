#ifndef RESPONSEMANAGER_HPP
#define RESPONSEMANAGER_HPP

#include "../molresponse/response_parameters.h"
#include "ResponseState.hpp"
#include <madness/chem/SCF.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>
#include <vector>

using namespace madness;

class ResponseManager {
public:
  ResponseManager(World &world, const ResponseParameters &r_params);

  void setProtocol(World &world, double L, double thresh, int override_k = -1);

  // Getters
  [[nodiscard]] double getVtol() const { return vtol; }
  [[nodiscard]] poperatorT getCoulombOp() const { return coulop; }
  [[nodiscard]] std::vector<std::shared_ptr<real_derivative_3d>>
  getGradOp() const {
    return gradop;
  }
  // Compute and store a response state
  void computeState(const ResponseState &state);
  // Check if a response state is already converged
  [[nodiscard]] bool isConverged(const ResponseState &state) const;

  [[nodiscard]] ResponseParameters params() const { return r_params; }

private:
  double vtol;
  poperatorT coulop;
  std::vector<std::shared_ptr<real_derivative_3d>> gradop;
  ResponseParameters r_params;
};

#endif // RESPONSEMANAGER_HPP
