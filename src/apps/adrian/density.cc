
#include "adrian/density.h"

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "adrian/global_functions.h"
#include "adrian/property_functions.h"
#include "adrian/property_operators.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;

// base class for a density
// operator used to create it
// homogeneous sol----x and y functions
// particular sol --- depends on lower order functions used to create it
// it also needs an xc functional
// The Rparams and Gparmas used to create the density
//

FirstOrderDensity::FirstOrderDensity(World &world, ResponseParameters Rparams,
                                     GroundParameters Gparams) {
  this->Rparams = Rparams;
  this->Gparams = Gparams;

  TDHF calc(world, Rparams, Gparams);
  if (calc.Rparams.property) {
    calc.ComputeFrequencyResponse(world);
  } else {
    calc.solve(world);
  }
  // stuff i want to save
  Rparams = calc.GetResponseParameters();
  Gparams = calc.GetGroundParameters();
  // right now everything uses copy
  property = calc.Rparams.response_type;
  omega = calc.GetFrequencyOmega();

  if (property.compare("dipole") == 0) {
    if (world.rank() == 0) print("creating dipole property operator");
    property_operator = Property(world, "dipole");
  } else if (property.compare("nuclear") == 0) {
    if (world.rank() == 0) print("creating nuclear property operator");
    property_operator = Property(world, "nuclear", Gparams.molecule);
  }

  x = calc.GetResponseFunctions("x");
  y = calc.GetResponseFunctions("y");

  num_response_states = x.size();
  num_ground_states = x[0].size();
  // get the response densities for our states
  rho_omega = calc.transition_density(world, Gparams.orbitals, x, y);
}

int FirstOrderDensity::GetNumberResponseStates() { return num_response_states; }
int FirstOrderDensity::GetNumberGroundStates() { return num_ground_states; }
VectorFunction3DT FirstOrderDensity::GetDensityVector() { return rho_omega; }
const Molecule FirstOrderDensity::GetMolecule() { return Gparams.molecule; }
TensorT FirstOrderDensity::GetFrequencyOmega() { return omega; }
ResponseParameters FirstOrderDensity::GetResponseParameters() {
  return Rparams;
}

void FirstOrderDensity::PrintDensityInformation() {
  // print
  //
  print("Response Density Information");
  print(property, " response at", omega(0, 0), "frequency using ", Rparams.xc,
        " exchange functional");
  print("Number of Response States : ", num_response_states);
  print("Number of Ground States : ", num_ground_states);
}

void FirstOrderDensity::PlotResponseDensity(World &world) {
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
Tensor<double> FirstOrderDensity::ComputeSecondOrderPropertyTensor(
    World &world) {
  // do some printing before we compute so we know what we are working with

  for (size_t i = 0; i < property_operator.operator_vector.size(); i++) {
    if (world.rank() == 0) {
      print("property operator vector i = ", i,
            "norm = ", property_operator.operator_vector[i].norm2());
    }
  }

  return matrix_inner(world, rho_omega, property_operator.operator_vector,
                      true);
}
void FirstOrderDensity::SaveDensity(World &world, std::string name) {
  // Archive to write everything to
  archive::ParallelOutputArchive ar(world, name.c_str(), 1);
  // Just going to enforce 1 io server

  ar &property;
  ar &omega;
  ar &num_response_states;
  ar &num_ground_states;
  ar &xcf;
  // Save response functions x and y
  // x first
  for (int i = 0; i < num_response_states; i++) {
    for (int j = 0; j < num_ground_states; j++) {
      ar &x[i][j];
    }
  }

  // y second
  for (int i = 0; i < num_response_states; i++) {
    for (int j = 0; j < num_ground_states; j++) {
      ar &y[i][j];
    }
  }
}
// Load a response calculation
void FirstOrderDensity::LoadDensity(World &world, std::string name) {
  // The archive to read from
  archive::ParallelInputArchive ar(world, name.c_str());
  // Reading in, in this order;

  ar &property;
  if (property.compare("dipole") == 0) {
    if (property.compare("dipole") == 0) {
      if (world.rank() == 0) print("creating dipole property operator");
      this->property_operator = Property(world, "dipole");
    } else if (property.compare("nuclear") == 0) {
      if (world.rank() == 0) print("creating nuclear property operator");
      this->property_operator = Property(world, "nuclear", Gparams.molecule);
    }
  }

  ar &omega;
  ar &num_response_states;
  ar &num_ground_states;
  ar &xcf;

  x = ResponseFunction(world, num_response_states, num_ground_states);
  y = ResponseFunction(world, num_response_states, num_ground_states);

  for (int i = 0; i < Rparams.states; i++) {
    for (unsigned int j = 0; j < Gparams.num_orbitals; j++) {
      ar &x[i][j];
    }
  }
  world.gop.fence();

  for (int i = 0; i < Rparams.states; i++) {
    for (unsigned int j = 0; j < Gparams.num_orbitals; j++) {
      ar &y[i][j];
      world.gop.fence();
    }
  }
}
