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

#include "Plot_VTK.h"

#include <madness/mra/mra.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "../chem/molecule.h"
namespace madness {
void do_vtk_plots(World &world,
                  int npt_plot,
                  double L,
                  int plotlo,
                  int plothi,
                  Molecule molecule,
                  std::vector<real_function_3d> densities,
                  std::string name) {
  // Stuff needed to plot
  //
  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);

  std::string geo_file;
  const char *filename;
  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};

  // Plot the whole box?
  Vector<double, 3> box_lo{-L, -L, -L};
  Vector<double, 3> box_hi{L, L, L};

  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE *f = 0;
  geo_file = vtk_dir + "/geometry.xyz";

  f = fopen(geo_file.c_str(), "w");

  // Write the header
  fprintf(f, "%zu", molecule.natom());
  fprintf(f, "\n\n");

  // Get the data
  std::vector<Vector<double, 3>> coords = molecule.get_all_coords_vec();

  // Write the data
  size_t Natoms = molecule.natom();
  for (size_t i = 0; i < Natoms; i++) {
    std::fprintf(f,
                 "%5s   %16.12f %16.12f %16.12f\n",
                 atomic_number_to_symbol(molecule.get_atom_number(i)).c_str(),
                 coords[i][0],
                 coords[i][1],
                 coords[i][2]);
  }

  // Clean up
  fclose(f);

  std::string response_file;
  // Needed to plot the full electron density
  real_function_3d rho = real_factory_3d(world);

  // Plot each orbital requested
  for (int i = plotlo; i < plothi; i++) {
    // Add to total
    rho += densities[i];

    // Create filename in such a way that visit associates them together
    response_file = vtk_dir + "/" + name + std::to_string(i) + ".vts";
    filename = response_file.c_str();

    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(densities[i],
                            "electrondensity",
                            world,
                            filename,
                            box_lo,
                            box_hi,
                            points,
                            true,
                            false);
    plotvtk_end<3>(world, filename, true);
  }
  std::string b;

  // Plot the full density
  b = vtk_dir + "/" + "total-electrondensity.vts";
  filename = b.c_str();
  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(rho,
                          "total-electrondensity",
                          world,
                          filename,
                          box_lo,
                          box_hi,
                          points,
                          true,
                          false);
  plotvtk_end<3>(world, filename, true);
}
void do_vtk_plots(World &world,
                  int npt_plot,
                  double L,
                  Molecule molecule,
                  real_function_3d &rho_0,
                  std::vector<real_function_3d> &rho_omega,
                  std::vector<real_function_3d> &ground_orbitals,
                  X_space &Chi) {
  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);
  std::string geo_file;
  const char *filename;
  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};
  // Plot the whole box?
  Vector<double, 3> box_lo{-L, -L, -L};
  Vector<double, 3> box_hi{L, L, L};
  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE *f = 0;
  geo_file = vtk_dir + "/geometry.xyz";

  f = fopen(geo_file.c_str(), "w");
  // Write the header
  fprintf(f, "%zu", molecule.natom());
  fprintf(f, "\n\n");
  // Get the data
  std::vector<Vector<double, 3>> coords = molecule.get_all_coords_vec();
  // Write the data
  size_t Natoms = molecule.natom();
  for (size_t i = 0; i < Natoms; i++) {
    std::fprintf(f,
                 "%5s   %16.12f %16.12f %16.12f\n",
                 atomic_number_to_symbol(molecule.get_atom_number(i)).c_str(),
                 coords[i][0],
                 coords[i][1],
                 coords[i][2]);
  }
  // Clean up
  fclose(f);
  std::string rho0_file = "rho_0.vts";
  filename = rho0_file.c_str();

  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(rho_0,
                          "ground_density",
                          world,
                          filename,
                          box_lo,
                          box_hi,
                          points,
                          true,
                          false);
  plotvtk_end<3>(world, filename, true);
  // ground orbitals
  std::string g_orb_file = "phi_";
  std::string fname;
  for (size_t i = 0; i < Chi.num_orbitals(); ++i) {
    fname = vtk_dir + "/" + g_orb_file + std::to_string(i) + ".vts";
    filename = fname.c_str();
    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(ground_orbitals[i],
                            "ground_orbital",
                            world,
                            filename,
                            box_lo,
                            box_hi,
                            points,
                            true,
                            false);
    plotvtk_end<3>(world, filename, true);
  }

  std::string x_orb_file = "x_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    for (size_t j = 0; j < Chi.num_orbitals(); ++j) {
      fname = vtk_dir + "/" + x_orb_file + std::to_string(i) + "_" +
              std::to_string(j) + ".vts";
      filename = fname.c_str();
      // VTK plotting stuff
      plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
      plotvtk_data<double, 3>(Chi.X[i][j],
                              "x_orbitals",
                              world,
                              filename,
                              box_lo,
                              box_hi,
                              points,
                              true,
                              false);
      plotvtk_end<3>(world, filename, true);
    }
  }
  std::string y_orb_file = "y_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    for (size_t j = 0; j < Chi.num_orbitals(); ++j) {
      fname = vtk_dir + "/" + y_orb_file + std::to_string(i) + "_" +
              std::to_string(j) + ".vts";
      filename = fname.c_str();
      // VTK plotting stuff
      plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
      plotvtk_data<double, 3>(Chi.Y[i][j],
                              "y_orbitals",
                              world,
                              filename,
                              box_lo,
                              box_hi,
                              points,
                              true,
                              false);
      plotvtk_end<3>(world, filename, true);
    }
  }

  std::string rho_1_file = "rho_1_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    fname = vtk_dir + "/" + rho_1_file + std::to_string(i) + ".vts";
    filename = fname.c_str();
    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(rho_omega[i],
                            "transition_density",
                            world,
                            filename,
                            box_lo,
                            box_hi,
                            points,
                            true,
                            false);
    plotvtk_end<3>(world, filename, true);
  }

  //
  //
  //
  //
  //
}
}  // namespace madness

// Dueces
