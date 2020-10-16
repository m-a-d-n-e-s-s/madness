#ifndef SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"

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
  void PlotResponseDensity(World &world) {
    // Doing line plots along each axis
    // Doing line plots along each axis
    if (world.rank() == 0) print("\n\nStarting plots");
    coord_3d lo, hi;
    char plotname[500];
    double Lp = std::min(Gparams.L, 24.0);
    if (world.rank() == 0) print("x:");
    // x axis
    lo[0] = 0.0;
    lo[1] = 0.0;
    lo[2] = 0.0;
    hi[0] = Lp;
    hi[1] = 0.0;
    hi[2] = 0.0;

    for (int i = 0; i < num_response_states; i++) {
      sprintf(plotname, "plot_transition_density_%d_%d_x.plt",
              FunctionDefaults<3>::get_k(), i);
      plot_line(plotname, 5001, lo, hi, rho_omega[i]);
    }
  }
};

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
