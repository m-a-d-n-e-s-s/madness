#ifndef SRC_APPS_ADRIAN_DENSITY_H_
#define SRC_APPS_ADRIAN_DENSITY_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "adrian/global_functions.h"
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
class FirstOrderDensity {
 private:
  std::string property;  // excited state, nuclear,dipole
  // operator used create first order density
  Property property_operator;  // dipole, nuclear, or none
  Tensor<double> omega;        // frequency or fruquencies

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

  FirstOrderDensity(World &world, ResponseParameters Rparams,
                    GroundParameters Gparams);

  int GetNumberResponseStates();
  int GetNumberGroundStates();
  VectorFunction3DT GetDensityVector();
  const Molecule GetMolecule();
  TensorT GetFrequencyOmega();
  ResponseParameters GetResponseParameters();

  void PrintDensityInformation();

  void PlotResponseDensity(World &world);

  Tensor<double> ComputeSecondOrderPropertyTensor(World &world);
  void SaveDensity(World &world, std::string name) {
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
  void LoadDensity(World &world, std::string name) {
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
};
#endif  // SRC_APPS_ADRIAN_DENSITY_H_
