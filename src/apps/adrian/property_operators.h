
#ifndef SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_
#define SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"

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
  const int atom;
  const int axis;

 public:
  MolecularDerivativeFunctor(const Molecule &molecule, int atom, int axis);
  double operator()(const coordT &x) const;
  std::vector<coordT> special_points() const;
};  // namespace madness
// A proerty class...creates a object with operator vector and property name
// Used to compute proerties or compute rhs vectors
class Property {
 public:
  int num_operators;  // number of operators in vectors
  std::string property;
  std::vector<real_function_3d> operator_vector;

  // default constructor
  Property();
  Property(World &world, std::string property_type);

  Property(World &world, std::string property_type, Molecule molecule);
};
// returns tensor alpha_ij=inner(rho_b,op_c)
Tensor<double> ComputeSecondOrderPropertyTensor(World &world,
                                                const VectorFunction3DT &rho_b,
                                                const Property &op_c);

void PrintSecondOrderAnalysis(World &world, const Tensor<double> alpha_tensor,
                              const Tensor<double> omega,
                              const ResponseParameters Rparams);

#endif  // SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_
