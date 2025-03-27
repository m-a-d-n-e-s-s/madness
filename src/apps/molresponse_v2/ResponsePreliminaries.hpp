#ifndef RESPONSEPRELIMINARIES_HPP
#define RESPONSEPRELIMINARIES_HPP

#include "ResponseState.hpp"
#include <madness/mra/funcdefaults.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>
#include <vector>

using namespace madness;

struct ResponsePreliminaries {

  real_function_3d V_local;  // Nuclear + Coulomb + optionally XC potential 
  Tensor<double> Hamiltonian; // Hamiltonian matrix
  Tensor<double>
      Hamiltonian_no_diag; // Hamiltonian matrix without diagonal elements

  ResponsePreliminaries() = default;
};

#endif // RESPONSEPRELIMINARIES_HPP
