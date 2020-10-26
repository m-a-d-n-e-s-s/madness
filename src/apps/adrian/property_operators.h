
#ifndef SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_
#define SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
namespace madness {
// Type definitions
typedef Tensor<double> tensorT;
typedef Function<double, 3> functionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> functorT;
typedef FunctionFactory<double, 3> factoryT;
typedef Vector<double, 3> coordT;

//
class MolecularDerivativeFunctor : public FunctionFunctorInterface<double, 3> {
  typedef Vector<double, 3> coordT;

 private:
  const Molecule &molecule;
  const int atom;
  const int axis;

 public:
  MolecularDerivativeFunctor(const Molecule &molecule, int atom, int axis)
      : molecule(molecule), atom(atom), axis(axis) {}

  double operator()(const coordT &x) const {
    return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0],
                                                            x[1], x[2]);
  }

  std::vector<coordT> special_points() const {
    return std::vector<coordT>(1, molecule.get_atom(atom).get_coords());
  }
};  // namespace madness
// A proerty class...creates a object with operator vector and property name
// Used to compute proerties or compute rhs vectors
class Property {
 public:
  int num_operators;  // number of operators in vectors
  std::string property;
  std::vector<real_function_3d> operator_vector;

  // default constructor
  Property() : num_operators(0), property(""), operator_vector() {}
  Property(World &world, std::string property_type)
      : num_operators(3), operator_vector(num_operators) {
    property = property_type;
    MADNESS_ASSERT(property.compare("dipole") == 0);
    for (int i = 0; i < 3; i++) {
      std::vector<int> f(3, 0);
      f[i] = 1;
      operator_vector.at(i) = real_factory_3d(world).functor(
          real_functor_3d(new BS_MomentFunctor(f)));
    }
  }

  Property(World &world, std::string property_type, Molecule molecule)
      : num_operators(molecule.natom() * 3), operator_vector(num_operators) {
    property = property_type;
    MADNESS_ASSERT(property.compare("nuclear") == 0);

    vecfuncT dv(molecule.natom() * 3);  // default constructor for vector?

    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
      for (int axis = 0; axis < 3; ++axis) {
        // question here....MolecularDerivativeFunctor takes derivative with
        // respect to axis atom and axis
        functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
        // here we save
        operator_vector.at(atom * 3 + axis) =
            functionT(factoryT(world)
                          .functor(func)
                          .nofence()
                          .truncate_on_project()
                          .truncate_mode(0));

        print("norm of vector ", operator_vector[atom * 3 + axis].norm2());
        // need to project
        //        operator_vector[atom * 3 + axis] = mul_sparse(
        //           world, dv[atom * 3 + axis], Gparams.orbitals,
        //           Rparams.small);

        // project rhs vectors for state

        // core projector contribution
      }
    }
  }
};

}  // namespace madness

#endif  // SRC_APPS_ADRIAN_PROPERTY_OPERATORS_H_
