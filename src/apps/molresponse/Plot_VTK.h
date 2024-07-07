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

#ifndef SRC_APPS_molresponse_PLOT_VTK_H_
#define SRC_APPS_molresponse_PLOT_VTK_H_

#include <madness/mra/mra.h>

#include <cstdint>
#include <string>
#include <vector>

#include "../chem/molecule.h"
#include "x_space.h"
namespace madness {

    void write_molecules_to_file(const Molecule &molecule, const std::string &geo_file);
    void do_response_orbital_vtk_plots(World &world, int npt_plot, double L, const Molecule &molecule, const vector_real_function_3d &ground_orbs, const response_matrix &responseMatrix);
    void do_response_density_vtk_plots(World &world, int npt_plot, double L, const Molecule &molecule, const real_function_3d &ground_density, const vector_real_function_3d &response_density);

    void do_vtk_plots(World &world,
                      int npt_plot,
                      double L,
                      int lowest_orbital,
                      int highest_orbital,
                      const Molecule &molecule,
                      std::vector<real_function_3d> densities,
                      const std::string &name);

    void do_vtk_plots(World &world,
                      int npt_plot,
                      double L,
                      Molecule molecule,
                      real_function_3d &rho_0,
                      std::vector<real_function_3d> &rho_omega,
                      std::vector<real_function_3d> &ground_orbitals,
                      X_space &Chi);

}// namespace madness
#endif// SRC_APPS_molresponse_PLOT_VTK_H_

// Dueces
