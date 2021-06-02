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


  $Id$
*/

/// \file molresponse.cc
/// \brief Molecular Response DFT code
/// \defgroup molresponse The molecular density funcitonal response code
#include <chem/SCF.h>
#include <madness/world/worldmem.h>
#include <stdlib.h>

#include "TDDFT.h"  // All response functions/objects enter through this
#include "molresponse/density.h"
#include "molresponse/global_functions.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
static inline int file_exists(const char* inpname) {
  struct stat buffer;
  size_t rc = stat(inpname, &buffer);
  return (rc == 0);
}
#endif

using namespace madness;

int main(int argc, char** argv) {
  initialize(argc, argv);
  {  // limite lifetime of world so that finalize() can execute cleanly
    World world(SafeMPI::COMM_WORLD);
    molresponse::start_timer(world);
    try {
      startup(world, argc, argv, true);
      print_meminfo(world.rank(), "startup");
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3> >(world)));

      std::cout.precision(6);
      // This makes a default input file name of 'input'
      const char* inpname = "input";
      // Process 0 reads input information and broadcasts
      for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
          inpname = argv[i];
          break;
        }
      }

      if (world.rank() == 0) print("input filename: ", inpname);
      if (!file_exists(inpname)) throw "input file not found";
      // first step is to read the input for r_params and g_params
      GroundParameters g_params;
      ResponseParameters r_params;
      r_params.read_and_set_derived_values(world, inpname, "response");
      std::string ground_file = r_params.archive();
      g_params.read(world, ground_file);

      density_vector d1 = set_density_type(world, r_params, g_params);
      if (world.rank() == 0) {
        g_params.print_params();
        r_params.print();
        d1.PrintDensityInformation();
      }
      // Create the TDDFT object
      if (r_params.load_density()) {
        print("Loading Density");
        d1.LoadDensity(world, r_params.load_density_file(), r_params, g_params);
      } else {
        print("Computing Density");
        d1.compute_response(world);
      }
      //
      // densityTest.PlotResponseDensity(world);
      d1.PrintDensityInformation();

      if (r_params.response_type().compare("dipole") == 0) {  //
        print("Computing Alpha");
        Tensor<double> alpha = d1.ComputeSecondOrderPropertyTensor(world);
        print("Second Order Analysis");
        d1.PrintSecondOrderAnalysis(world, alpha);
      }
    } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
      print(e);
      error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
      print(e);
      error("caught a Tensor exception");
    } catch (const char* s) {
      print(s);
      error("caught a string exception");
    } catch (const std::string& s) {
      print(s);
      error("caught a string (class) exception");
    } catch (const std::exception& e) {
      print(e.what());
      error("caught an STL exception");
    } catch (...) {
      error("caught unhandled exception");
    }

    world.gop.fence();
    world.gop.fence();
  }

  finalize();

  return 0;
}
