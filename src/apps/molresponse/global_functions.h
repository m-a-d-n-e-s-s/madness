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
#endif  // SRC_APPS_molresponse_GLOBAL_FUNCTIONS_H_
