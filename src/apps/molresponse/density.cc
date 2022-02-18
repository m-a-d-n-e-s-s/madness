// Copyright 2021 Adrian Hurtado
#include "molresponse/density.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "TDDFT.h"
#include "molresponse/global_functions.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef std::vector<real_function_3d> VectorFunction3DT;

density_vector::density_vector(World &world, ResponseParameters &other_rparams,
                               GroundStateCalculation &other_gparams)
    : num_states(other_rparams.num_states()),
      num_orbitals(other_rparams.num_orbitals()),
      property(),
      r_params(other_rparams),  // should be a copy
      g_params(other_gparams),
      Chi(world, num_states, num_orbitals),
      PQ(world, num_states, num_orbitals),
      orbitals(copy(world, other_gparams.orbitals())),
      molecule(other_gparams.molecule()) {

  if (r_params.excited_state()) {
    this->omega = Tensor<double>(r_params.num_states());
  } else {
    this->omega = Tensor<double>(1);
    this->omega(0, 0) = r_params.omega();
  }
}

// right now everything uses copy

ResponseParameters density_vector::GetResponseParameters() { return r_params; }

density_vector set_density_type(World &world, ResponseParameters & R, GroundStateCalculation & G) {
  if (R.excited_state()) {
    return excited_state_density_vector(world, R, G);
  } else if (R.dipole())
    return dipole_density_vector(world, R, G);
  else if (R.nuclear()) {
    return nuclear_density_vector(world, R, G);
  } else if (R.second_order()) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);
  } else if (R.third_order()) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);

  } else {
    MADNESS_EXCEPTION("what is this????", 0);
    return density_vector(world, R, G);
  }
};
