

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
#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && \
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

#include <madness/misc/info.h>
#include <madness/world/worldmem.h>

#include "Drivers.hpp"
#include "MoldftLib.hpp"
#include "MolresponseLib.hpp"
#include "CCLib.hpp"
#include "ParameterManager.hpp"

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
    print_header1("MADQC -- Multiresolution Quantum Chemistry Code ");
  }

  {
    // limit lifetime of world so that finalize() can execute cleanly
    START_TIMER(world);
    /*try {*/
    // Load info for MADNESS numerical routines
    startup(world, argc, argv, true);
    if (world.rank() == 0) print(info::print_revision_information());

    commandlineparser parser(argc, argv);
    print("input file  :",parser.value("input"));
    print("workflow    :",parser.value("workflow"));


    // Note that the molecule is part of Params!
    Params pm(world, parser);
    pm.get<CalculationParameters>().print("dft");

    // Create workflow
    qcapp::Workflow wf;

    // always do single-point DFT, can be reset to SCF or Nemo
    std::shared_ptr<Application> reference;
    std::string user_workflow = parser.value("workflow");

    if (user_workflow=="response") {
      reference.reset(new SCFApplication<moldft_lib,SCF>(world, pm));
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
      // directory with the ground-state DFT outputs
      std::filesystem::path gsDir = "MyCalc/task_0/moldft";

      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
          std::make_unique<ResponseApplication<molresponse_lib>>(world, pm,
                                                                 gsDir)));
    } else if (user_workflow=="mp2" or user_workflow=="cc2") {

      // set the tensor type
      TensorType tt = TT_2D;
      FunctionDefaults<6>::set_tensor_type(tt);

      pm.get<CalculationParameters>().set_derived_value("k", 5);
      reference.reset(new SCFApplication<moldft_lib,Nemo>(world, pm));
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));

      auto moldft_app = std::dynamic_pointer_cast<SCFApplication<moldft_lib,Nemo>>(reference);
      if (!moldft_app) {
        MADNESS_EXCEPTION("Could not cast reference to SCFApplication<moldft_lib>", 1);
      }
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
          std::make_unique<CC2Application<cc_lib,SCFApplication<moldft_lib,Nemo>>>(world, pm, *moldft_app)));
      pm.get<CCParameters>().print("cc2");

    } else if (user_workflow=="cis") {
      reference.reset(new SCFApplication<moldft_lib,Nemo>(world, pm));
      auto moldft_app = std::dynamic_pointer_cast<SCFApplication<moldft_lib,Nemo>>(reference);
      if (!moldft_app) {
        MADNESS_EXCEPTION("Could not cast reference to SCFApplication<moldft_lib>", 1);
      }
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
        std::make_unique<TDHFApplication<SCFApplication<moldft_lib,Nemo>>>(world,pm,*moldft_app)));

    } else if (user_workflow=="oep") {
      // add tight convergence criteria
      CalculationParameters& cparam=pm.get<CalculationParameters>();
      auto convergence_crit=cparam.get<std::vector<std::string> >("convergence_criteria");
      if (std::find(convergence_crit.begin(),convergence_crit.end(),"each_energy")==convergence_crit.end()) {
        convergence_crit.push_back("each_energy");
      }
      cparam.set_derived_value("convergence_criteria",convergence_crit);

      reference.reset(new SCFApplication<moldft_lib,Nemo>(world, pm));
      auto moldft_app = std::dynamic_pointer_cast<SCFApplication<moldft_lib,Nemo>>(reference);
      if (!moldft_app) {
        MADNESS_EXCEPTION("Could not cast reference to SCFApplication<moldft_lib>", 1);
      }

      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
      wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
        std::make_unique<OEPApplication<SCFApplication<moldft_lib,Nemo>>>(world,pm,*moldft_app)));
    } else {
      std::string msg= "Unknown workflow: " + user_workflow +
          "\nAvailable workflows are: response, mp2, cc2, cis";
      MADNESS_EXCEPTION(msg.c_str(), 1);
    }

    // Execute both in "MyCalc" directory
    wf.run("MyCalc");

    if (true) {
      qcapp::Workflow opt_wf;

      std::function<std::unique_ptr<Application>(Params)> scfFactory =
          [&](Params p) {
            return std::make_unique<SCFApplication<moldft_lib>>(world, p);
          };

      opt_wf.addDriver(
          std::make_unique<qcapp::OptimizeDriver>(world, scfFactory, pm));
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
  }  // world is dead -- ready to finalize
  finalize();

  return 0;
}
