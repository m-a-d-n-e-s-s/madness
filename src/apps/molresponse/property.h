#ifndef SRC_APPS_molresponse_PROPERTY_H_
#define SRC_APPS_molresponse_PROPERTY_H_

#include <molresponse/response_functions.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include <madness/mra/function_interface.h>
#include <madness/mra/functypedefs.h>

// Type definitions
using namespace madness;


//
class MolecularDerivativeFunctor : public FunctionFunctorInterface<double, 3> {

 private:
  const Molecule &molecule;
  const size_t atom;
  const size_t axis;

 public:
  MolecularDerivativeFunctor(const Molecule &molecule, size_t atom, size_t axis);
  double operator()(const coord_3d &x) const override;
  vector_coord_3d special_points() const override;
};  // namespace madness
// A property class...creates a object with operator vector and property name
// Used to compute properties or compute rhs vectors
class PropertyBase {
 public:
  size_t num_operators;  // number of operators in vectors
  vector_real_function_3d operator_vector;

  // default constructor
  PropertyBase();
};

class DipoleVector : public PropertyBase {
 public:
  DipoleVector(World &world);
};

class NuclearVector : public PropertyBase {
 public:
  NuclearVector(World &world, Molecule &molecule);
};

#endif  // SRC_APPS_molresponse_PROPERTY_H_
