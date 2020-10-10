/* My function to save a vtk file that visit can use
 *
 *   Input params:
 *
 *       npt_plot  -  number of points in each direction for plot
 *       L         -  box size
 *       plotlo    -  lowest orbital number to plot
 *       plothi    -  highest orbital number to plot
 *       molecule  -  molecule object, for creating the .xyz file
 *       densities -  vector of densities to be ploted
 *       name      -  name you would like for orbital plots
 */

#ifndef SRC_APPS_ADRIAN_PLOT_VTK_H_
  #define SRC_APPS_ADRIAN_PLOT_VTK_H_

  #include <madness/mra/mra.h>

  #include <cstdint>
  #include <string>
  #include <vector>

  #include "../chem/molecule.h"
namespace madness {
void do_vtk_plots(World &world, int npt_plot, double L, int plotlo, int plothi,
                  Molecule molecule, std::vector<real_function_3d> densities,
                  std::string name);

}  // namespace madness
#endif  // SRC_APPS_ADRIAN_PLOT_VTK_H_

// Dueces
