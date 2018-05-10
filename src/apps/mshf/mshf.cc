/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  ------------------
  This code was written by I. Sagert and G. I. Fann
*/  


/*
// For google-code version of MADNESS
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/vmra.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <assert.h>
*/

#include "input.h"
#include "ground.h"


// Main Code
int main(int argc, char** argv)
{
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    
    Skyrme skyrme("skyrme.input");
    Setup setup("mshf.input");
    Nuclear nuclear("mshf.input");
    Mixing mixing("mshf.input");
    Output output("mshf.input");
    Additional additional("mshf.input");

    if( world.rank() == 0){std::cout << " " << std::endl;}
    if( world.rank() == 0){std::cout << "---------------------------" << std::endl;}
    if( world.rank() == 0){std::cout << "MADNESS Skyrme-Hartree-Fock" << std::endl;}
    if( world.rank() == 0){std::cout << "---------------------------" << std::endl;}
    if( world.rank() == 0){std::cout << " " << std::endl;}


    // basic setup parameters 
    double A          = setup.A;               // Mass number
    double Z          = setup.Z;               // Charge number
    double length     = setup.length;          // Lenght of one side of the sim. volume (cubic)
    int    initial    = setup.initial;         // Initialization parameter
    int    boundary   = setup.boundary;        // Boundary conditions
    double knumber    = setup.knumber;         // Wavelet order
    double thresh     = setup.thresh;          // Truncation threshold
    int    IO_nodes   = setup.IO_nodes;        // Number of nodes for checkpointing
    int    project    = setup.project;         // Switch for projection to higher wavelet numbers

    // nuclear matter parameters 
    int    spinorbit  = nuclear.spinorbit;       // Spin-orbit potential
    int    meff       = nuclear.meff;            // Effective mass potential
    int    jellium    = nuclear.jellium;         // Jellium approximation 
    int    screening  = nuclear.screening;       // Electron screening
    double screenl    = nuclear.screenl;         // Electron screening length
    int    lap_comp   = nuclear.lap_comp;        // Calculation of laplacian

    // mixing parameters 
    int    avg_pot    = mixing.avg_pot;         // Switch for potential mixing
    int    avg_lap    = mixing.avg_lap;         // Switch for laplacian mixing
    int    avg_wav    = mixing.avg_wav;         // Switch for state mixing

    // output parameters 
    int    vtk_output = output.vtk_output;      // Switch for vtk output
    int    txt_output = output.txt_output;      // Switch for txt output 

    // additional parameters 
    int    use_infile = additional.use_infile;      // Switch for using resolution parameters from inpu and ignoring checkpointed parameters
    double prec       = additional.prec;            // Additional precision factor
    double chi        = additional.chi;             // Fraction of old wavefunctions for mixing
    double tol        = additional.tol;             // Tolerance (used in e.g. BSH and Coulomb)
    double brad       = additional.brad;            // Gaussian smoothing radius
    int    timing     = additional.timing;          // Switch for timing 


    double t0         = skyrme.t0;
    double t1         = skyrme.t1;
    double t2         = skyrme.t2;
    double t3         = skyrme.t3;
    double t4         = skyrme.t4;
    double x0         = skyrme.x0;
    double x1         = skyrme.x1;
    double x2         = skyrme.x2;
    double x3         = skyrme.x3;
    double bp4        = skyrme.bp4;
    double alpha      = skyrme.alpha;
    double k_fn       = skyrme.k_fn;
    double k_fp       = skyrme.k_fp;

    real_tensorT b(5);
    real_tensorT bp(5);

    b[0]  = t0 * (1.e0 + 0.5e0 * x0);
    b[1]  = (t1 + 0.5e0 * x1 * t1 + t2 + 0.5e0 * x2 * t2)/4.e0;
    b[2]  = (3.e0 * t1*(1.e0 + 0.5e0 * x1) - t2 * (1.e0 + 0.5e0 * x2))/8.e0;
    b[3]  = t3 * (1.e0 + 0.5e0 * x3)/4.e0;
    b[4]  = 0.5e0 * t4;
    bp[0] = t0 * (0.5e0 + x0);
    bp[1] = (t1 * (0.5e0 + x1) - t2 * (0.5e0 + x2))/4.e0;
    bp[2] = (3.e0 * t1 * (0.5e0 + x1) + t2 * (0.5e0 + x2))/8.e0;
    bp[3] = t3 * (0.5e0 + x3)/4.e0;
    bp[4] = bp4;


    try {
        startup(world,argc,argv);
        //print_meminfo(world.rank(), "startup");
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap< Key<3> >(world)));
        std::cout.precision(6);

        ground_state(world, A, Z, length, initial, boundary, jellium, spinorbit, meff, screening, screenl, 
             avg_pot, avg_lap, avg_wav, lap_comp, prec, vtk_output, txt_output, use_infile, knumber, 
             thresh, chi, tol, brad, IO_nodes, b, bp, alpha, k_fn, k_fp, project, timing);
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    // Nearly all memory will be freed at this point
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
    finalize();

    return 0;
}


