#include "ResponseManager.hpp"

#include <madness/chem/SCF.h>

#include "CalculationParameters.h"
#include "madness/world/world.h"

using namespace madness;

ResponseManager::ResponseManager(World &world, const CalculationParameters &params) : calc_params(params), vtol(0.0) {}

// Initially loaded at maximum polynomial order (settings_.k)

void ResponseManager::setProtocol(World &world, double L, double thresh, int override_k) {
  int k;
  if (thresh >= 0.9e-2)
    k = 4;
  else if (thresh >= 0.9e-4)
    k = 6;
  else if (thresh >= 0.9e-6)
    k = 8;
  else if (thresh >= 0.9e-8)
    k = 10;
  else
    k = 12;

  if (override_k == -1) {
    FunctionDefaults<3>::set_k(k);
    //        	param.k=k;
  } else {
    FunctionDefaults<3>::set_k(calc_params.k());
  }
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_apply_randomize(false);
  FunctionDefaults<3>::set_project_randomize(false);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
  GaussianConvolution1DCache<double>::map.clear();
  double safety = 0.1;
  vtol = FunctionDefaults<3>::get_thresh() * safety;
  coulop = poperatorT(CoulombOperatorPtr(world, calc_params.lo(), 0.001 * thresh));
  gradop = gradient_operator<double, 3>(world);

  if (world.rank() == 0) {
    print("\nMRA Protocol Set: thresh =", thresh, ", k =", k, "\n");
  }
  world.gop.fence();
}
