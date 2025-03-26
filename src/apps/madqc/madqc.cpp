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

#include <fstream>
#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) &&                   \
    defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// static inline int file_exists(const char *inpname) {
//     struct stat buffer;
//     int rc = stat(inpname, &buffer);
//     return (rc == 0);
// }

#endif

#include "madqc/calc_factory.hpp"
#include "madqc/calc_manager.hpp"
#include "madqc/parameter_manager.hpp"
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>

using namespace madness;

static double ttt, sss;

static void START_TIMER(World &world) {
  world.gop.fence();
  ttt = wall_time();
  sss = cpu_time();
}

static void END_TIMER(World &world, const char *msg) {
  ttt = wall_time() - ttt;
  sss = cpu_time() - sss;
  if (world.rank() == 0)
    printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  if (world.rank() == 0) {
    print_header1("MADQC -- Multiresolution Quantum Chemistry Code ");
  }

  { // limit lifetime of world so that finalize() can execute cleanly
    START_TIMER(world);
    //  try {
    // Load info for MADNESS numerical routines
    startup(world, argc, argv, true);
    if (world.rank() == 0)
      print(info::print_revision_information());

    commandlineparser parser(argc, argv);
    ParameterManager params;
    // Initialize the necessary components
    path input_file(argv[1]);
    params = ParameterManager(world, input_file);

    auto task_params = params.get_task_params();

    auto method = task_params.method;
    auto driver = task_params.driver;
    auto properties = task_params.properties;
    auto molecule = params.get_molecule();

    /*if (world.rank() == 0) {*/
    /*  task_params.print();*/
    /*}*/

    // Create the CalcManager using the factory function
    // if driver == energy use createCalcManager
    // if driver == optimize use createOptimizeManager
    // if driver == custom use createCustomManager

    path cwd = std::filesystem::current_path();
    path root;
    if (driver == "energy") {
      root = cwd;

    } else if (driver == "optimize") {
      root = cwd / "optimize";
      std::filesystem::create_directory(root);
      std::filesystem::current_path(root);
    }

    auto calc_manager =
        createEnergyDriver(world, method, params, properties, root);

    if (driver == "optimize") {

      // calc_manager = createOptimizationDriver(method, params);
      auto &opt_params = params.get_optimization_params();
      MolOpt opt(opt_params.get_maxiter(),            // geometry max iter
                 0.1,                                 // geometry step size
                 opt_params.get_value_precision(),    // value precision
                 opt_params.get_geometry_tolerence(), // geometry tolerance
                 1e-3,                                // XTOL
                 1e-5,                                // EPREC
                 opt_params.get_gradient_precision(), // gradient precision
                 (world.rank() == 0) ? 1 : 0,         // print_level
                 opt_params.get_algopt());            // algorithm options

      auto new_molecule = opt.optimize(molecule, *calc_manager);
      // Get output directory
      std::filesystem::current_path(cwd);

      if (world.rank() == 0) {
        json final_output;
        print("Optimization completed successfully.");
        std::ifstream ifs(calc_manager->get_output_path());
        ifs >> final_output;
        print(final_output.dump(4));
        ifs.close();
        std::ofstream ofs(root / "output.json");
        ofs << final_output[method].dump(4);
        ofs.close();
      }
    } else if (driver == "energy") {
      calc_manager->runCalculations(molecule.get_all_coords().flat());
    } else {
      throw std::runtime_error("Invalid driver");
    }
    // Run the calculations


    if (parser.key_exists("help")) {
      ParameterManager::help();
    }
    /*} catch (const SafeMPI::Exception& e) {*/
    /*  print(e);*/
    /*  error("caught an MPI exception");*/
    /*} catch (const madness::MadnessException& e) {*/
    /*  print(e);*/
    /*  error("caught a MADNESS exception");*/
    /*} catch (const madness::TensorException& e) {*/
    /*  print(e);*/
    /*  error("caught a Tensor exception");*/
    /*} catch (const char* s) {*/
    /*  print(s);*/
    /*  error("caught a string exception");*/
    /*} catch (const std::string& s) {*/
    /*  print(s);*/
    /*  error("caught a string (class) exception");*/
    /*} catch (const std::exception& e) {*/
    /*  print(e.what());*/
    /*  error("caught an STL exception");*/
    /*} catch (...) {*/
    /*  error("caught unhandled exception");*/
    /*}*/

    // Nearly all memory will be freed at this point
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  } // world is dead -- ready to finalize
  finalize();

  return 0;
}
