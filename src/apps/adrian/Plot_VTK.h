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

#ifndef MADNESS_APPS_TDA_PLOT_VTK
#define MADNESS_APPS_TDA_PLOT_VTK

void do_vtk_plots(World &world, int npt_plot, double L, int plotlo, int plothi,
                  Molecule molecule, std::vector<real_function_3d> densities,
                  std::string name) {
  // Stuff needed to plot
  std::string b;
  const char *filename;
  Vector<long, 3> points{npt_plot, npt_plot, npt_plot};

  // Plot the whole box?
  Vector<double, 3> box_lo{-L, -L, -L};
  Vector<double, 3> box_hi{L, L, L};

  // Write an .xyz file with current geometry (to deal with molecular
  // reorientations that might occur)
  FILE *f = 0;
  b = "geometry.xyz";
  f = fopen(b.c_str(), "w");

  // Write the header
  fprintf(f, "%zu", molecule.natom());
  fprintf(f, "\n\n");

  // Get the data
  std::vector<Vector<double, 3>> coords = molecule.get_all_coords_vec();

  // Write the data
  int Natoms = molecule.natom();
  for (int i = 0; i < Natoms; i++) {
    std::fprintf(f, "%5s   %16.12f %16.12f %16.12f\n",
                 atomic_number_to_symbol(molecule.get_atom_number(i)).c_str(),
                 coords[i][0], coords[i][1], coords[i][2]);
  }

  // Clean up
  fclose(f);

  // Needed to plot the full electron density
  real_function_3d rho = real_factory_3d(world);

  // Plot each orbital requested
  for (int i = plotlo; i < plothi; i++) {
    // Add to total
    rho += densities[i];

    // Create filename in such a way that visit associates them together
    b = name + std::to_string(i) + ".vts";
    filename = b.c_str();

    // VTK plotting stuff
    plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
    plotvtk_data<double, 3>(densities[i], "electrondensity", world, filename,
                            box_lo, box_hi, points, true, false);
    plotvtk_end<3>(world, filename, true);
  }

  // Plot the full density
  b = "total-electrondensity.vts";
  filename = b.c_str();
  plotvtk_begin<3>(world, filename, box_lo, box_hi, points, true);
  plotvtk_data<double, 3>(rho, "total-electrondensity", world, filename, box_lo,
                          box_hi, points, true, false);
  plotvtk_end<3>(world, filename, true);
}
#endif

// Dueces
