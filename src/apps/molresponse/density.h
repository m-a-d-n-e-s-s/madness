#ifndef SRC_APPS_molresponse_DENSITY_H_
#define SRC_APPS_molresponse_DENSITY_H_

#include <molresponse/response_functions.h>
#include <molresponse/x_space.h>

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
typedef std::vector<real_function_3d> VectorFunction3DT;

// base class for a density

// operator used to create it.
// homogeneous sol----x and y functions.
// particular sol --- depends on lower order functions used to create it.
// it also needs a xc functional
// The r_params and Gparmas used to create the density.
//
class density_vector {
 protected:
  // operator used create first order density
  Tensor<double> omega;        // frequency or frequencies
  const size_t num_states;     // number of response states
  const size_t num_orbitals;   // number of ground state orbitals
  const std::string property;  // excited state, nuclear,dipole

  const ResponseParameters r_params;  // Response Parameters
  const GroundParameters g_params;

  XCfunctional xcf;                // xc functional
  PropertyBase property_operator;  // dipole, nuclear, or none
  X_space Chi;
  X_space PQ;
  VectorFunction3DT orbitals;

  // first order frequency response densities
  VectorFunction3DT rho_omega;  // the response density vector
  Molecule molecule;

 public:
  friend class TDDFT;
  // Collective constructor
  density_vector(World& world, ResponseParameters r_params, GroundParameters g_params);
  density_vector(const density_vector& other) = default;

  ResponseParameters GetResponseParameters();
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

density_vector set_density_type(World& world, ResponseParameters R, GroundParameters G);
#endif  // SRC_APPS_molresponse_DENSITY_H_
