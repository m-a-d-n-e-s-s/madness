#include "molresponse/property.h"

#include <molresponse/TDDFT.h>
#include <molresponse/response_functions.h>

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
typedef Vector<double, 3> coordT;
typedef std::vector<real_function_3d> VectorFunction3DT;

MolecularDerivativeFunctor::MolecularDerivativeFunctor(const Molecule &molecule, size_t atom, size_t axis)
    : molecule(molecule), atom(atom), axis(axis) {}

double MolecularDerivativeFunctor::operator()(const coordT &x) const {
  return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0], x[1], x[2]);
}

std::vector<coordT> MolecularDerivativeFunctor::special_points() const {
  return std::vector<coordT>(1, molecule.get_atom(atom).get_coords());
}

PropertyBase::PropertyBase() : num_operators(0), operator_vector() {}

DipoleVector::DipoleVector(World &world) : PropertyBase() {
  num_operators = 3;
  operator_vector = vecfuncT(num_operators);
  for (size_t i = 0; i < 3; i++) {
    std::vector<int> f(3, 0);
    f[i] = 1;
    operator_vector.at(i) = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
    // prsize_t k L truncation
    print("norm of dipole function ", operator_vector[i].norm2());
  }

  truncate(world, operator_vector, true);

  for (size_t i = 0; i < 3; i++) {
    print("norm of dipole function after truncate ", operator_vector[i].norm2());
  }
};

NuclearVector::NuclearVector(World &world, Molecule &molecule) : PropertyBase() {
  num_operators = size_t(molecule.natom() * 3);
  operator_vector = vecfuncT(num_operators);

  vecfuncT dv(molecule.natom() * 3);  // default constructor for vector?

  print("Creating Nuclear Derivative Operator");
  for (size_t atom = 0; atom < molecule.natom(); ++atom) {
    for (size_t axis = 0; axis < 3; ++axis) {
      // question here....MolecularDerivativeFunctor takes derivative with
      // respect to axis atom and axis
      FunctorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
      // here we save
      operator_vector.at(atom * 3 + axis) =
          FunctionT(FactoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
      // prsize_t k L truncation
      print("norm of derivative function ", operator_vector[atom * 3 + axis].norm2());
    }
  }

  truncate(world, operator_vector, true);

  for (size_t atom = 0; atom < molecule.natom(); ++atom) {
    for (size_t axis = 0; axis < 3; ++axis) {
      print("norm of derivative function after truncate ", operator_vector[atom * 3 + axis].norm2());
    }
  }
}
