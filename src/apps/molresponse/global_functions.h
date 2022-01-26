#ifndef SRC_APPS_molresponse_GLOBAL_FUNCTIONS_H_
#define SRC_APPS_molresponse_GLOBAL_FUNCTIONS_H_

#include <memory>
#include <string>
#include <vector>

#include "molresponse/density.h"
#include "response_parameters.h"

struct CalcParams {
  GroundStateCalculation ground_calculation;
  Molecule molecule;
  ResponseParameters response_parameters;
};
void print_molecule(World &world, const GroundStateCalculation &g_params);
class plotCoords {
 public:
  coord_3d lo, hi;

  plotCoords() {
    lo[0] = 0.0;
    lo[1] = 0.0;
    lo[2] = 0.0;
    hi[0] = 0.0;
    hi[1] = 0.0;
    hi[2] = 0.0;
  }
  plotCoords(size_t direction, double Lp) {
    lo[0] = 0.0;
    lo[1] = 0.0;
    lo[2] = 0.0;

    hi[0] = 0.0;
    hi[1] = 0.0;
    hi[2] = 0.0;

    lo[direction] = -Lp;
    hi[direction] = Lp;
  }
};
plotCoords SetPlotCoord(size_t i, double Lp);

CalcParams initialize_calc_params(World &world, std::string input_file);
/// Mask function to switch from 0 to 1 smoothly at boundary
// Masking function to switch from 0 to 1 smoothly at boundary
// Pulled from SCF.h
inline double mask1(double x) {
  /* Iterated first beta function to switch smoothly
   from 0->1 in [0,1].  n iterations produce 2*n-1
   zero derivatives at the end points. Order of polyn
   is 3^n.

  Currently use one iteration so that first deriv.
      is zero at interior boundary and is exactly representable
          by low order multiwavelet without refinement */

  x = (x * x * (3. - 2. * x));
  return x;
}

static double mask3(const coord_3d &ruser) {
  coord_3d rsim;
  user_to_sim(ruser, rsim);
  double x = rsim[0], y = rsim[1], z = rsim[2];
  double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
  double rlo = 1.0 / lo;

  if (x < lo)
    result *= mask1(x * rlo);
  else if (x > hi)
    result *= mask1((1.0 - x) * rlo);
  if (y < lo)
    result *= mask1(y * rlo);
  else if (y > hi)
    result *= mask1((1.0 - y) * rlo);
  if (z < lo)
    result *= mask1(z * rlo);
  else if (z > hi)
    result *= mask1((1.0 - z) * rlo);

  return result;
}
#endif  // SRC_APPS_molresponse_GLOBAL_FUNCTIONS_H_
