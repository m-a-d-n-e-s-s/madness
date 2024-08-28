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
*/

#include "calc_manager.hpp"

#include <madness/chem/SCF.h>
#include <madness/chem/molopt.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include "calc_factory.hpp"
#include "parameter_manager.hpp"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// static inline int file_exists(const char *inpname) {
//     struct stat buffer;
//     int rc = stat(inpname, &buffer);
//     return (rc == 0);
// }

#endif

using namespace madness;

static double ttt, sss;

static void START_TIMER(World& world) {
  world.gop.fence();
  ttt = wall_time();
  sss = cpu_time();
}

static void END_TIMER(World& world, const char* msg) {
  ttt = wall_time() - ttt;
  sss = cpu_time() - sss;
  if (world.rank() == 0)
    printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

int main(int argc, char** argv) {

  World& world = initialize(argc, argv);
  if (world.rank() == 0) {
    print_header1("MADQC -- Multiresolution Quantum Chemsitry Code ");
  }

  {  // limit lifetime of world so that finalize() can execute cleanly
    START_TIMER(world);
    try {
      // Load info for MADNESS numerical routines
      startup(world, argc, argv, true);
      if (world.rank() == 0)
        print(info::print_revision_information());

      commandlineparser parser(argc, argv);
      ParameterManager params;
      // Initialize the necessary components
      Molecule molecule;  // Initialize your molecule here
      params = ParameterManager(world, {"input.json"});

      // Define the properties to be calculated for each model
      property_map properties;

      // Example: Setting up properties for Moldft
      properties["moldft"] = {{"energy", true}, {"gradient", false}, {"dipole", true}};

      // Example: Setting up properties for Response
      properties["response"] = {{"alpha", true}, {"beta", true}, {"shg", false}};

      // Choose the model name
      std::string model_name = "response";  // This could be "moldft", "MP2", "CIS", etc.

      // Create the CalcManager using the factory function
      auto calc_manager = createCalcManager(model_name, params, properties);

      // Run the calculations
      calc_manager->runCalculations(world);

      std::cout << "Calculations completed successfully." << std::endl;

      if (parser.key_exists("help")) {
        ParameterManager::help();
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

    // Nearly all memory will be freed at this point
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }  // world is dead -- ready to finalize
  finalize();

  return 0;
}
