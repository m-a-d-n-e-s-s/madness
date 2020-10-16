#ifndef SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <memory>
#include <string>
#include <vector>

namespace madness {
// base class for a density
// operator used to create it
// homogeneous sol----x and y functions
// particular sol --- depends on lower order functions used to create it
// it also needs an xc functional
// The Rparams and Gparmas used to create the density
//
class FirstOrderDensity {
 private:
  Tensor<double> omega;
  std::string property;

  int num_response_states;
  int num_ground_states;

  XCfunctional xcf;

  ResponseParameters Rparams;
  GroundParameters Gparams;

  ResponseFunction x;
  ResponseFunction y;

  // first order frequency response densities
  std::vector<real_function_3d> rho_omega;

 public:
  // Collective constructor
  FirstOrderDensity(World &world, const char *filename)
      : FirstOrderDensity(
            world,
            (world.rank() == 0 ? std::make_shared<std::ifstream>(filename)
                               : nullptr)) {}

  FirstOrderDensity(World &world, std::shared_ptr<std::istream> density_input) {
    TDHF calc(world, density_input);
    if (calc.Rparams.property) {
      calc.ComputeFrequencyResponse(world);
    } else {
      calc.solve(world);
    }
    // right now everything uses copy
    property = calc.Rparams.property;
    omega = calc.Rparams.omega;

    x = calc.GetResponseFunctions("x");
    y = calc.GetResponseFunctions("y");

    num_response_states = x.size();
    num_ground_states = x[0].size();
    // stuff i want to save
    Rparams = calc.GetResponseParameters();
    Gparams = calc.GetGroundParameters();
    // get the response densities for our states
    rho_omega = calc.transition_density(world, Gparams.orbitals, x, y);
  }
};

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
