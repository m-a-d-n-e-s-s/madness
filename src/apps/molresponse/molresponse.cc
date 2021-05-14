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

namespace madness {
int main(int argc, char** argv) {
  // Initialize MADNESS mpi
  initialize(argc, argv);
  {  // limite lifetime of world so that finalize() can execute cleanly
    World world(SafeMPI::COMM_WORLD);
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
    // first step is to read the input for Rparams and Gparams
    ResponseParameters Rparams;

    Rparams.read_and_set_derived_values(world, inpname, "response");

    TDDFT calc(world, inpname);

    // Read the ground parameters from the archive
    Gparams.read(world, Rparams.archive);
    if (world.rank() == 0) {
      Gparams.print_params();
      print_molecule(world, Gparams);
    }
    // if a proerty calculation set the number of states
    if (Rparams.property) {
      Rparams.SetNumberOfStates(Gparams.molecule);
    }

    // print params
    if (world.rank() == 0) {
      Rparams.print_params();
    }
    // Broadcast to all other nodes
    density_vector densityTest = SetDensityType(world, Rparams.response_type, Rparams, Gparams);
    // Create the TDDFT object
    if (Rparams.load_density) {
      print("Loading Density");
      densityTest.LoadDensity(world, Rparams.load_density_file, Rparams, Gparams);
    } else {
      print("Computing Density");
      densityTest.compute_response(world);
    }
    //
    // densityTest.PlotResponseDensity(world);
    densityTest.PrintDensityInformation();

    if (Rparams.response_type.compare("dipole") == 0) {  //
      print("Computing Alpha");
      Tensor<double> alpha = densityTest.ComputeSecondOrderPropertyTensor(world);
      print("Second Order Analysis");
      densityTest.PrintSecondOrderAnalysis(world, alpha);
    }

    world.gop.fence();
    world.gop.fence();
  }  // world is dead -- ready to finalize
  finalize();

  return 0;
}
}  // namespace madness
