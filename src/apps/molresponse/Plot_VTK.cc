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
#include <string>
#include <vector>

#if defined(__has_include)
#if __has_include(<filesystem>)
#define MADCHEM_HAS_STD_FILESYSTEM
// <filesystem> is not reliably usable on Linux with gcc < 9
#if defined(__GNUC__)
#if __GNUC__ >= 7 && __GNUC__ < 9
#undef MADCHEM_HAS_STD_FILESYSTEM
#endif
#endif
#if defined(MADCHEM_HAS_STD_FILESYSTEM)

#include <filesystem>
namespace madness {
void write_molecules_to_file(const Molecule& molecule, const std::string& geo_file, const double& scale) {

  FILE* f = nullptr;
  f = fopen(geo_file.c_str(), "w");
  // Write the header
  fprintf(f, "%zu", molecule.natom());
  fprintf(f, "\n\n");
  // Get the data
  auto coords = molecule.get_all_coords_vec();

  // Write the data
  size_t Natoms = molecule.natom();
  for (size_t i = 0; i < Natoms; i++) {
    auto coords_i = coords[i];
    auto x = coords_i[0] / scale;
    auto y = coords_i[1] / scale;
    auto z = coords_i[2] / scale;
    std::fprintf(f, "%5s   %16.12f %16.12f %16.12f\n", atomic_number_to_symbol(molecule.get_atomic_number(i)).c_str(),
                 x, y, z);
  }

  // Clean up
  fclose(f);
}

void do_response_orbital_vtk_plots(World& world, int npt_plot, double L, const Molecule& molecule,
                                   const vector_real_function_3d& ground_orbs, const response_matrix& responseMatrix) {
  // Stuff needed to plot
  //
  double box_size = L / 2;
  Vector<double, 3> box_lo{-box_size, -box_size, -box_size};
  Vector<double, 3> box_hi{box_size, box_size, box_size};

  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);

  std::string geo_file;

  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};
  // Plot the whole box?

  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE* f = nullptr;
  geo_file = vtk_dir + "/geometry.xyz";
  write_molecules_to_file(molecule, geo_file);

  std::string response_file;
  // Needed to plot the full electron density
  real_function_3d rho = real_factory_3d(world);

  int orb_num = 0;

  auto orbital_file = vtk_dir + "/" + "orbitals.vts";
  plotvtk_begin<3>(world, orbital_file.c_str(), box_lo, box_hi, points, true);

  std::for_each(ground_orbs.begin(), ground_orbs.end(), [&](const auto& phi0_i) {
    auto orb_name = "/phi0_" + std::to_string(orb_num);
    plotvtk_data<double, 3>(phi0_i, orb_name.c_str(), world, orbital_file.c_str(), box_lo, box_hi, points, true, false);
    orb_num++;
  });

  auto state_number = 0;
  std::for_each(responseMatrix.begin(), responseMatrix.end(), [&](const auto& xy) {
    auto num_orbitals = xy.size() / 2;
    // plot the x first
    auto orb_num = 0;
    std::for_each(xy.begin(), xy.begin() + num_orbitals, [&](const auto& xi) {
      auto field_name = "x_orbital_" + std::to_string(state_number) + "_" + std::to_string(orb_num);
      plotvtk_data<double, 3>(xi, field_name.c_str(), world, orbital_file.c_str(), box_lo, box_hi, points, true, false);
      orb_num++;
    });
    orb_num = 0;
    std::for_each(xy.begin() + num_orbitals, xy.end(), [&](const auto& yi) {
      auto field_name = "y_orbital_" + std::to_string(state_number) + "_" + std::to_string(orb_num);
      plotvtk_data<double, 3>(yi, field_name.c_str(), world, orbital_file.c_str(), box_lo, box_hi, points, true, false);
      orb_num++;
    });
    state_number++;
  });
  plotvtk_end<3>(world, orbital_file.c_str(), true);
}
void do_response_density_vtk_plots(World& world, int npt_plot, double L, const Molecule& molecule,
                                   const real_function_3d& ground_density,
                                   const vector_real_function_3d& response_density) {
  // Stuff needed to plot
  //
  double box_size = L / 2;
  Vector<double, 3> box_lo{-box_size, -box_size, -box_size};
  Vector<double, 3> box_hi{box_size, box_size, box_size};

  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);

  std::string geo_file;
  const char* filename;

  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};

  std::string response_file;
  //***********************************ground density plot
  auto density_file = vtk_dir + "/" + "rho0.vts";
  filename = density_file.c_str();
  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(ground_density, "ground_density", world, filename, box_lo, box_hi, points, true, false);
  plotvtk_end<3>(world, filename, true);

  //***********************************ground density plot
  int state_number = 0;
  std::for_each(response_density.begin(), response_density.end(), [&](const auto& rho_i) {
    auto density_file = vtk_dir + "/" + "response_rho_" + std::to_string(state_number) + ".vts";
    filename = density_file.c_str();
    auto field_name = "r_density_" + std::to_string(state_number);
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(rho_i, field_name.c_str(), world, filename, box_lo, box_hi, points, true, false);
    plotvtk_end<3>(world, filename, true);
    state_number++;
  });
}
void do_response_density_vtk_plots_new(World& world, int npt_plot, double L, const Molecule& molecule,
                                       const real_function_3d& ground_density,
                                       const vector_real_function_3d& response_density,
                                       const std::string& density_name) {
  double box_size = L;
  Vector<double, 3> box_lo{-box_size, -box_size, -box_size};
  Vector<double, 3> box_hi{box_size, box_size, box_size};

  std::string vtk_dir = "../vtk_plots";
  std::filesystem::create_directories(vtk_dir);

  FILE* f = nullptr;
  std::string geo_file;
  geo_file = vtk_dir + "/geometry.xyz";
  double scale = 1.0;
  write_molecules_to_file(molecule, geo_file, scale);

  const char* filename;

  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};

  std::string response_file;
  //***********************************ground density plot
  auto density_file = vtk_dir + "/" + density_name + ".vts";
  filename = density_file.c_str();
  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(ground_density, "ground", world, filename, box_lo, box_hi, points, true, false);
  //***********************************ground density plot
  int state_number = 0;
  std::for_each(response_density.begin(), response_density.end(), [&](const auto& rho_i) {
    auto field_name = "response_" + std::to_string(state_number);
    plotvtk_data<double, 3>(rho_i, field_name.c_str(), world, filename, box_lo, box_hi, points, true, false);
    state_number++;
  });
  plotvtk_end<3>(world, filename, true);
}
void do_vtk_plots(World& world, int npt_plot, double L, int lowest_orbital, int highest_orbital,
                  const Molecule& molecule, std::vector<real_function_3d> densities, const std::string& name) {
  // Stuff needed to plot
  //
  Vector<double, 3> box_lo{-L, -L, -L};
  Vector<double, 3> box_hi{L, L, L};

  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);

  std::string geo_file;
  const char* filename;

  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};
  // Plot the whole box?

  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE* f = nullptr;

  std::string response_file;
  // Needed to plot the full electron density
  real_function_3d rho = real_factory_3d(world);

  // Plot each orbital requested
  for (int i = lowest_orbital; i < highest_orbital; i++) {
    // Add to total
    rho += densities[i];

    // Create filename in such a way that visit associates them together
    response_file = vtk_dir + "/" + name + std::to_string(i) + ".vts";
    filename = response_file.c_str();

    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);

    // function
    // field anem
    // world
    // file name
    // plot lo
    // plot hi
    // npts
    // binary
    // plot refile
    plotvtk_data<double, 3>(densities[i], "density", world, filename, box_lo, box_hi, points, true, false);
    plotvtk_end<3>(world, filename, true);
  }
  std::string b;

  // Plot the full density
  b = vtk_dir + "/" + "total-electrondensity.vts";
  filename = b.c_str();
  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(rho, "total-electrondensity", world, filename, box_lo, box_hi, points, true, false);
  plotvtk_end<3>(world, filename, true);
}
void do_vtk_plots(World& world, int npt_plot, double L, Molecule molecule, real_function_3d& rho_0,
                  std::vector<real_function_3d>& rho_omega, std::vector<real_function_3d>& ground_orbitals,
                  X_space& Chi) {
  std::string vtk_dir = "vtk_plots";
  std::filesystem::create_directories(vtk_dir);
  std::string geo_file;
  const char* filename;
  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};
  // Plot the whole box?
  Vector<double, 3> box_lo{-L, -L, -L};
  Vector<double, 3> box_hi{L, L, L};
  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE* f = 0;
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
    std::fprintf(f, "%5s   %16.12f %16.12f %16.12f\n", atomic_number_to_symbol(molecule.get_atomic_number(i)).c_str(),
                 coords[i][0], coords[i][1], coords[i][2]);
  }
  // Clean up
  fclose(f);
  std::string rho0_file = "rho_0.vts";
  filename = rho0_file.c_str();

  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(rho_0, "ground_density", world, filename, box_lo, box_hi, points, true, false);
  plotvtk_end<3>(world, filename, true);
  // ground orbitals
  std::string g_orb_file = "phi_";
  std::string fname;
  for (size_t i = 0; i < Chi.num_orbitals(); ++i) {
    fname = vtk_dir + "/" + g_orb_file + std::to_string(i) + ".vts";
    filename = fname.c_str();
    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(ground_orbitals[i], "ground_orbital", world, filename, box_lo, box_hi, points, true, false);
    plotvtk_end<3>(world, filename, true);
  }

  std::string x_orb_file = "x_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    for (size_t j = 0; j < Chi.num_orbitals(); ++j) {
      fname = vtk_dir + "/" + x_orb_file + std::to_string(i) + "_" + std::to_string(j) + ".vts";
      filename = fname.c_str();
      // VTK plotting stuff
      plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
      plotvtk_data<double, 3>(Chi.x[i][j], "x_orbitals", world, filename, box_lo, box_hi, points, true, false);
      plotvtk_end<3>(world, filename, true);
    }
  }
  std::string y_orb_file = "y_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    for (size_t j = 0; j < Chi.num_orbitals(); ++j) {
      fname = vtk_dir + "/" + y_orb_file + std::to_string(i) + "_" + std::to_string(j) + ".vts";
      filename = fname.c_str();
      // VTK plotting stuff
      plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
      plotvtk_data<double, 3>(Chi.y[i][j], "y_orbitals", world, filename, box_lo, box_hi, points, true, false);
      plotvtk_end<3>(world, filename, true);
    }
  }

  std::string rho_1_file = "rho_1_";
  for (size_t i = 0; i < Chi.num_states(); ++i) {
    fname = vtk_dir + "/" + rho_1_file + std::to_string(i) + ".vts";
    filename = fname.c_str();
    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(rho_omega[i], "transition_density", world, filename, box_lo, box_hi, points, true, false);
    plotvtk_end<3>(world, filename, true);
  }

  //
  //
  //
  //
  //
}

}  // namespace madness

#endif
#endif
#endif

#include "../chem/molecule.h"
