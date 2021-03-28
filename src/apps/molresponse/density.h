#ifndef SRC_APPS_molresponse_DENSITY_H_
#define SRC_APPS_molresponse_DENSITY_H_

#include <TDDFT.h>
#include <response_functions.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "molresponse/global_functions.h"
#include "molresponse/property.h"

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
 protected:
  std::string property;  // excited state, nuclear,dipole
  // operator used create first order density
  Property property_operator;  // dipole, nuclear, or none
  Tensor<double> omega;        // frequency or frequencies

  size_t num_states;         // number of response states
  size_t num_ground_states;  // number of ground state orbitals

  XCfunctional xcf;  // xc functional

  ResponseParameters Rparams;  // Response Parameters
  GroundParameters Gparams;    // Ground Parameters

  X_space Chi;
  response_space x;  // The x response functions virt/occ
  response_space y;  // occ/virt

  response_space P;  // rhs vector P
  response_space Q;  // rhs vector Q

  // first order frequency response densities
  VectorFunction3DT rho_omega;  // the response density vector

 public:
  // Collective constructor
  FirstOrderDensity(ResponseParameters Rparams, GroundParameters Gparams);

  virtual void ComputeResponse(World& world);

  size_t GetNumberResponseStates();
  VectorFunction3DT ComputeDensityVector(World& world, bool is_static);
  size_t GetNumberGroundStates();
  VectorFunction3DT GetDensityVector();
  const Molecule GetMolecule();
  TensorT GetFrequencyOmega();
  ResponseParameters GetResponseParameters();

  void PrintDensityInformation();

  void PlotResponseDensity(World& world);

  Tensor<double> ComputeSecondOrderPropertyTensor(World& world);
  void PrintSecondOrderAnalysis(World& world,
                                const Tensor<double> alpha_tensor);
  void SaveDensity(World& world, std::string name);
  // Load a response calculation
  void LoadDensity(World& world,
                   std::string name,
                   ResponseParameters Rparams,
                   GroundParameters Gparams);
};

class DipoleDensity : public FirstOrderDensity {
 public:
  DipoleDensity(World& world, ResponseParameters R, GroundParameters G)
      : FirstOrderDensity(R, G) {
    this->property = Rparams.response_type;
    this->num_states = 3;
    this->num_ground_states = Gparams.num_orbitals;
    this->x = response_space(world, num_states, num_ground_states);
    this->y = response_space(world, num_states, num_ground_states);
    this->P = response_space(world, num_states, num_ground_states);
    this->Q = response_space(world, num_states, num_ground_states);
  }
};

class NuclearResponseDensity : public FirstOrderDensity {
 public:
  NuclearResponseDensity(World& world, ResponseParameters R, GroundParameters G)
      : FirstOrderDensity(R, G) {
    this->property = Rparams.response_type;
    this->num_states = Rparams.states;
    this->num_ground_states = Gparams.num_orbitals;
    this->x = response_space(world, num_states, num_ground_states);
    this->y = response_space(world, num_states, num_ground_states);
    this->P = response_space(world, num_states, num_ground_states);
    this->Q = response_space(world, num_states, num_ground_states);
  }
};

class ExcitedStateDensity : public FirstOrderDensity {
 public:
  ExcitedStateDensity(World& world, ResponseParameters R, GroundParameters G)
      : FirstOrderDensity(R, G) {
    this->property = Rparams.response_type;
    this->num_states = Rparams.states;
    this->num_ground_states = Gparams.num_orbitals;
    this->x = response_space(world, num_states, num_ground_states);
    this->y = response_space(world, num_states, num_ground_states);
    this->P = response_space(world, num_states, num_ground_states);
    this->Q = response_space(world, num_states, num_ground_states);
  }
};
#endif  // SRC_APPS_molresponse_DENSITY_H_
