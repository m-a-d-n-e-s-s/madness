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

// Type definitions

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;
//
class MolecularDerivativeFunctor : public FunctionFunctorInterface<double, 3> {
  typedef Vector<double, 3> coordT;

 private:
  const Molecule &molecule;
  const size_t atom;
  const size_t axis;

 public:
  MolecularDerivativeFunctor(const Molecule &molecule, size_t atom, size_t axis);
  double operator()(const coordT &x) const;
  std::vector<coordT> special_points() const;
};  // namespace madness
// A proerty class...creates a object with operator vector and property name
// Used to compute proerties or compute rhs vectors
class PropertyBase {
 public:
  size_t num_operators;  // number of operators in vectors
  std::vector<real_function_3d> operator_vector;

  // default constructor
  PropertyBase();

  PropertyBase(World &world, std::string property_type, Molecule molecule);
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
