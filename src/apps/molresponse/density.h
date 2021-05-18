#ifndef SRC_APPS_molresponse_DENSITY_H_
#define SRC_APPS_molresponse_DENSITY_H_

#include <TDDFT.h>
#include <response_functions.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "molresponse/global_functions.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/property.h"
#include "molresponse/response_parameters.h"

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
// The r_params and Gparmas used to create the density
//
class density_vector {
 protected:
  std::string property;  // excited state, nuclear,dipole
  // operator used create first order density
  Property property_operator;  // dipole, nuclear, or none
  Tensor<double> omega;        // frequency or frequencies

  size_t num_states;         // number of response states
  size_t num_ground_states;  // number of ground state orbitals

  XCfunctional xcf;  // xc functional

  ResponseParameters r_params;  // Response Parameters
  GroundParameters g_params;

  X_space Chi;
  X_space PQ;

  // first order frequency response densities
  VectorFunction3DT rho_omega;  // the response density vector

 public:
  // Collective constructor
  density_vector(World& world, ResponseParameters r_params, GroundParameters g_params);

  virtual void compute_response(World& world);

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
  void PrintSecondOrderAnalysis(World& world, const Tensor<double> alpha_tensor);
  void SaveDensity(World& world, std::string name);
  // Load a response calculation
  void LoadDensity(World& world, std::string name, ResponseParameters r_params, GroundParameters g_params);
};

class dipole_density_vector : public density_vector {
 public:
  dipole_density_vector(World& world, ResponseParameters R, GroundParameters G) : density_vector(world, R, G) {}
};

class nuclear_density_vector : public density_vector {
 public:
  nuclear_density_vector(World& world, ResponseParameters R, GroundParameters G) : density_vector(world, R, G) {}
};

class excited_state_density_vector : public density_vector {
 public:
  excited_state_density_vector(World& world, ResponseParameters R, GroundParameters G) : density_vector(world, R, G) {}
};

density_vector set_density_type(World& world, ResponseParameters R, GroundParameters G) {
  if (R.response_type().compare("excited_state") == 0) {
    return excited_state_density_vector(world, R, G);
  } else if (R.response_type().compare("dipole") == 0) {
    return dipole_density_vector(world, R, G);

  } else if (R.response_type().compare("nuclear") == 0) {
    return nuclear_density_vector(world, R, G);
  } else if (R.response_type().compare("order2") == 0) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);
  } else if (R.response_type().compare("order3") == 0) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);

  } else {
    MADNESS_EXCEPTION("what is this????", 0);
    return density_vector(world, R, G);
  }
};
#endif  // SRC_APPS_molresponse_DENSITY_H_
