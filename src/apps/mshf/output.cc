#include "output.h"
#include "states.h"

char direct[] = "./";

// Make output for plotting
void output(World& world,  const real_functionT& rho_p,
                           const real_functionT& rho_n,
                           const real_functionT& tau,
                           const real_functionT& lap_p,
                           const real_functionT& lap_n,
                           const real_functionT& U,
                           const real_functionT& Uc,
                           const real_functionT& Uex,
                           const int iter,
                           const double L,
                           const int vtk_output,
                           const int txt_output)

{
    if (world.rank() == 0) {print(" "); print("Make output: ", iter);}
    double Lp = std::min(L,24.0);

    real_functionT lap = lap_n + lap_p;
    real_functionT rho = rho_p + rho_n;

    // Paraview Output
    if (vtk_output == 1) {
        Vector<double, 3> plotlo, plothi;
        Vector<long, 3> npts;
        char filename[100];
        sprintf(filename, "%s/data_%d.vts",direct, iter);
        world.gop.fence();

        for (int i = 0; i < 3; i++) {plotlo[i] = -Lp; plothi[i] = Lp; npts[i] = 51;}
        world.gop.fence();

        plotvtk_begin(world, filename, plotlo, plothi, npts);
        plotvtk_data(rho,   "rho", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_p, "rho_p", world, filename, plotlo, plothi, npts);
        plotvtk_data(rho_n, "rho_n", world, filename, plotlo, plothi, npts);
        plotvtk_data(tau,   "tau", world, filename, plotlo, plothi, npts);
        plotvtk_data(lap,   "lap", world, filename, plotlo, plothi, npts);
        plotvtk_data(U,     "U", world, filename, plotlo, plothi, npts);
        plotvtk_end<3>(world, filename);
    }

    // Simple txt output 
    if (txt_output == 1) {
      coord_3d lo,hi;

      // Ouput for profile along y
      lo[0] = 0.0; lo[1] = -Lp; lo[2] = 0.0;
      hi[0] = 0.0; hi[1] =  Lp; hi[2] = 0.0;
      char plotname[500];
      sprintf(plotname, "%s/data_den_y_%d.txt",direct, iter);
      plot_line(plotname, 501, lo, hi, rho, tau, lap);

      // Output for profile along x
      lo[0] = -Lp; lo[1] = 0.0; lo[2] = 0.0;
      hi[0] =  Lp; hi[1] = 0.0; hi[2] = 0.0;
      sprintf(plotname, "%s/data_den_x_%d.txt",direct, iter);
      plot_line(plotname, 5001, lo, hi, rho, tau, lap);

      // Output for profiles along z
      lo[0] = 0.0; lo[1] = 0.0; lo[2] = -Lp;
      hi[0] = 0.0; hi[1] = 0.0; hi[2] =  Lp;
      sprintf(plotname, "%s/data_den_z_%d.txt",direct, iter);
      plot_line(plotname, 501, lo, hi, rho, tau, lap);
    
      // Output of potential and states
      lo[0] = -Lp; lo[1] = 0.0; lo[2] = 0.0;
      hi[0] =  Lp; hi[1] = 0.0; hi[2] = 0.0;
      sprintf(plotname, "%s/data_pot_x_%d.txt",direct, iter);
      plot_line(plotname, 501, lo, hi, U);
    }
    world.gop.fence();
}











