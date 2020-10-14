#ifndef SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <string>

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
  double omega;
  std::string property;

  ResponseFunction x;
  ResponseFunction y;

 public:
  FirstOrderDensity(World &world, std::shared_ptr<std::istream> density_input) {
    TDHF calc(world, density_input);
    if (calc.Rparams.property) {
      calc.ComputeFrequencyResponse(world);
    }
    property = calc.Rparams.property;
    omega = calc.Rparams.omega;
    x = calc.GetResponseFunctions("x");
    y = calc.GetResponseFunctions("y");
  }
};

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
