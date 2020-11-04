#ifndef SRC_APPS_ADRIAN_DENSITY_H_
#define SRC_APPS_ADRIAN_DENSITY_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "adrian/global_functions.h"
#include "adrian/property.h"

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
  FirstOrderDensity(ResponseParameters Rparams, GroundParameters Gparams);

  void ComputeDensity(World &world);

  int GetNumberResponseStates();
  int GetNumberGroundStates();
  VectorFunction3DT GetDensityVector();
  const Molecule GetMolecule();
  TensorT GetFrequencyOmega();
  ResponseParameters GetResponseParameters();

  void PrintDensityInformation();

  void PlotResponseDensity(World &world);

  Tensor<double> ComputeSecondOrderPropertyTensor(World &world);
  void PrintSecondOrderAnalysis(World &world,
                                const Tensor<double> alpha_tensor);
  void SaveDensity(World &world, std::string name);
  // Load a response calculation
  void LoadDensity(World &world, std::string name, ResponseParameters Rparams,
                   GroundParameters Gparams);
};
#endif  // SRC_APPS_ADRIAN_DENSITY_H_
