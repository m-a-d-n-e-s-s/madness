#include "molresponse/property.h"
#include <madness/chem/SCF.h>

#include <molresponse/response_functions.h>

#include <algorithm>
#include <vector>

using namespace madchem;


PropertyBase::PropertyBase() : num_operators(0), operator_vector() {}

DipoleVector::DipoleVector(World &world) : PropertyBase() {
  num_operators = 3;
  operator_vector = vecfuncT(num_operators);
  for (size_t i = 0; i < 3; i++) {
    std::vector<int> f(3, 0);
    f[i] = 1;
    operator_vector.at(i) = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
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
      functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
      // here we save
      operator_vector.at(atom * 3 + axis) =
          functionT (factoryT (world).functor(func).nofence().truncate_on_project().truncate_mode(0));
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
