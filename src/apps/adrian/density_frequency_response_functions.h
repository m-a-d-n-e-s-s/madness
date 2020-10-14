#ifndef SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_

#include <ResponseFunction2.h>

#include <string>

namespace madness {
// base class for a density
// operator used to create it
// homogeneous sol----x and y functions
// particular sol --- depends on lower order functions used to create it
//
class FirstOrderDensity {
 private:
  double omega;
  std::string property;

  ResponseFunction x;
  ResponseFunction y;

 public:
  FirstOrderDensity(double frequency, std::string operator_property) {
    omega = frequency;
    property = operator_property;
  }
};

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
