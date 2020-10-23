#ifndef SRC_APPS_ADRIAN_DENSITY_H_
#define SRC_APPS_ADRIAN_DENSITY_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "adrian/property_operators.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;

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
  Tensor<double> omega;  // frequency or fruquencies
  std::string property;  // excited state, nuclear,dipole
  // operator used create first order density
  Property property_operator;  // dipole, nuclear, or none

  int num_response_states;
  int num_ground_states;

  XCfunctional xcf;

  ResponseParameters Rparams;
  GroundParameters Gparams;

  ResponseFunction x;
  ResponseFunction y;

  // first order frequency response densities
  VectorFunction3DT rho_omega;

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
    property = calc.Rparams.response_type;
    omega = calc.Rparams.omega;
    if (property.compare("dipole") == 0) {
      property_operator = Property(world, "dipole");
    } else if (property.compare("nuclear") == 0) {
      property_operator = Property(world, "nuclear", Gparams.molecule);
    }

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

  int GetNumberResponseStates() { return num_response_states; }
  int GetNumberGroundStates() { return num_ground_states; }
  VectorFunction3DT GetDensityVector() { return rho_omega; }
  const Molecule GetMolecule() { return Gparams.molecule; }
  TensorT GetFrequencyOmega() { return omega; }

  void PrintDensityInformation() {
    // print
    //
    print("Response Density Information");
    print(property, " response at", omega(0, 0), "frequency using ", Rparams.xc,
          " exchange functional");
    print("Number of Response States : ", num_response_states);
    print("Number of Ground States : ", num_ground_states);
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
      std::snprintf(plotname, sizeof(plotname),
                    "plot_transition_density_%d_%d_x.plt",
                    FunctionDefaults<3>::get_k(), i);
      plot_line(plotname, 5001, lo, hi, rho_omega[i]);
    }
  }
  Tensor<double> ComputeSecondOrderPropertyTensor(World &world) {
    return matrix_inner(world, rho_omega, property_operator.operator_vector,
                        true);
  }
};

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_H_
