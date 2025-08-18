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
#include <madness_exception.h>

#include <CCLib.hpp>
#include <Drivers.hpp>
#include <MoldftLib.hpp>
#include <MolresponseLib.hpp>
#include <ParameterManager.hpp>

using namespace madness;

/// TODO:
///  - restart options
///  - nemo projections
///  - numerical parameters

void help(std::string wf) {
  print("Usage: madqc [options] [input_file]");
  print("\nOptions:");
  print("  --help=<workflow>           : show this help message");
  print("  --print_parameters=<group>  : print all parameters and exit");
  print(
      "  --workflow=<name>           : specify the workflow to run (default: "
      "scf)");
  print("\nAvailable workflows: scf, nemo, response, mp2, cc2, cis, oep");
  print("Available groups: dft, nemo, response, cc2, cis, oep, geometry");
  print("");
  if (wf == "scf") {
    print("madqc --wf=scf");
    print("molecular Hartree-Fock and DFT code");
    print("\nexamples: ");
    print("madqc --wf=scf --geometry=h2o.xyz");
    print("madqc --wf=scf --geometry=he --dft=\"k=8\",xc=lda");
    print("madqc --wf=scf --geometry=h2o --dft=\"k=8; econv=1.e-4; L=30\"");
  } else if (wf == "nemo") {
    print("madqc --wf=nemo");
    print("molecular Hartree-Fock and DFT code with regularized orbitals");
    print("\nexamples: ");
    print("madqc --wf=nemo --geometry=h2o ");
  } else if (wf == "response") {
    print("madqc --wf=response");
    print("Response theory for DFT and Hartree-Fock");
    print("\nexamples: ");
  } else if (wf == "mp2" or wf == "cc2") {
    print("madqc --wf=cc2");
    print("MP2/CC2 code for correlated wave functions");
    print("\nexamples: ");
    print(
        "madqc --wf=cc2 --geometry=h2o --dft=\"k=5; econv=1.e-4; L=30\" "
        "--cc2=\"freeze=1\"");
  } else if (wf == "cis") {
    print("madqc --wf=cis");
    print("CIS code for excited states");
  } else if (wf == "oep") {
    print("madqc --wf=oep");
    print("Optimized effective potential code for DFT");
  }
}

void print_parameters(World& world, const commandlineparser& parser,
                      const std::string group) {
  Params pm;
  if (group == "") {
    print("please specify a data group to print parameters for");
    print("\n  --print_parameters=<group>  : print all parameters and exit");
    print(
        "\nAvailable data groups: scf, nemo, response, cc2, cis, oep, "
        "geometry");
  } else if (group == "dft") {
    print("Available parameters for data group: dft");
    pm.get<CalculationParameters>().print();
  } else if (group == "nemo") {
    print("Available parameters for data group: nemo");
    pm.get<Nemo::NemoCalculationParameters>().print();
  } else if (group == "response") {
    print("Available parameters for data group: response");
    pm.get<ResponseParameters>().print();
  } else if (group == "cc2") {
    print("Available parameters for data group: cc2");
    pm.get<CCParameters>().print();
  } else if (group == "cis") {
    print("Available parameters for data group: cis");
    pm.get<TDHFParameters>().print();
  } else if (group == "oep") {
    print("Available parameters for data group: oep");
    pm.get<OEP_Parameters>().print();
  } else if (group == "geometry") {
    Molecule::GeometryParameters geometryparam;
    geometryparam.print("geometry", "end");
  } else {
    std::string msg = "Unknown data group: " + group +
                      "\nAvailable data group are: dft, nemo, response, cc2, "
                      "cis, oep, geometry\n";
    print(msg);
  }
}

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  if (world.rank() == 0) {
    print_header1("MADQC -- Multiresolution Quantum Chemistry Code ");
  }

  commandlineparser parser(argc, argv);

  if (parser.key_exists("help")) {
    help(parser.value("help"));
  } else if (parser.key_exists("print_parameters")) {
    print_parameters(world, parser, parser.value("print_parameters"));

  } else {
    // limit lifetime of world so that finalize() can execute cleanly
    try {
      // Load info for MADNESS numerical routines
      startup(world, argc, argv, true);
      if (world.rank() == 0) print(info::print_revision_information());

      // Create workflow
      qcapp::Workflow wf;
      std::string user_workflow = "scf";
      if (parser.key_exists("workflow"))
        user_workflow = parser.value("workflow");
      else if (parser.key_exists("wf"))
        user_workflow = parser.value("wf");
      else if (parser.key_exists("w"))
        user_workflow = parser.value("w");

      if (world.rank() == 0) {
        print("input file  :", parser.value("input"));
        print("workflow    :", user_workflow);
        print("");
      }

      // read in all parameters from the input file and the command line
      // logic and interdependent parameter follow later
      Params pm(world, parser);

      if (user_workflow == "scf") {
        auto reference = std::shared_ptr<Application>(
            new SCFApplication<moldft_lib, SCF>(world, pm));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));

      } else if (user_workflow == "nemo") {
        pm.get<CalculationParameters>().set_derived_value("k", 8);
        auto reference = std::shared_ptr<Application>(
            new SCFApplication<moldft_lib, Nemo>(world, pm));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));

      } else if (user_workflow == "response") {
        auto reference = std::shared_ptr<Application>(
            new SCFApplication<moldft_lib, SCF>(world, pm));
        //
        std::filesystem::path gsDir(pm.prefix() + "/task_0/moldft");

        // prefix/task0/moldft
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
        // prefix/task1/molresponse
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
            std::make_unique<ResponseApplication<molresponse_lib>>(world, pm,
                                                                   gsDir)));

      } else if (user_workflow == "mp2" or user_workflow == "cc2") {
        // set the tensor type
        TensorType tt = TT_2D;
        FunctionDefaults<6>::set_tensor_type(tt);

        // do the parameter logic and print parameters
        auto& calc_param = pm.get<CalculationParameters>();
        auto& cc_param = pm.get<CCParameters>();
        auto& molecule = pm.get<Molecule>();

        calc_param.set_derived_value("k", 5);
        calc_param.set_derived_value("print_level", 2);
        calc_param.set_derived_value("econv",
                                     cc_param.get<double>("thresh_6d") * 0.01);

        calc_param.set_derived_values(molecule);
        cc_param.set_derived_values();

        auto reference = std::shared_ptr<SCFApplication<moldft_lib, Nemo>>(
            new SCFApplication<moldft_lib, Nemo>(world, pm));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));

        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
            std::make_unique<
                CC2Application<cc_lib, SCFApplication<moldft_lib, Nemo>>>(
                world, pm, *reference)));

      } else if (user_workflow == "cis") {
        auto reference = std::shared_ptr<SCFApplication<moldft_lib, Nemo>>(
            new SCFApplication<moldft_lib, Nemo>(world, pm));

        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
            std::make_unique<TDHFApplication<SCFApplication<moldft_lib, Nemo>>>(
                world, pm, *reference)));

      } else if (user_workflow == "oep") {
        // add tight convergence criteria
        CalculationParameters& cparam = pm.get<CalculationParameters>();
        auto convergence_crit =
            cparam.get<std::vector<std::string>>("convergence_criteria");
        if (std::find(convergence_crit.begin(), convergence_crit.end(),
                      "each_energy") == convergence_crit.end()) {
          convergence_crit.push_back("each_energy");
        }
        cparam.set_derived_value("convergence_criteria", convergence_crit);

        auto reference = std::shared_ptr<Application>(
            new SCFApplication<moldft_lib, Nemo>(world, pm));
        auto moldft_app =
            std::dynamic_pointer_cast<SCFApplication<moldft_lib, Nemo>>(
                reference);
        if (!moldft_app) {
          MADNESS_EXCEPTION(
              "Could not cast reference to SCFApplication<moldft_lib>", 1);
        }

        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
            std::make_unique<OEPApplication<SCFApplication<moldft_lib, Nemo>>>(
                world, pm, *moldft_app)));
      } else {
        static std::string msg =
            "Unknown workflow: " + user_workflow +
            "\nAvailable workflows are: response, mp2, cc2, cis";
        MADNESS_EXCEPTION(msg.c_str(), 1);
      }

      // Execute both in "MyCalc" directory
      if (world.rank() == 0) {
        print_header1("Calculation parameters");
        pm.get<Molecule>().print();
        wf.print_parameters(world);
        print("");
      }

      if (world.rank() == 0) print_header1("Starting calculations");
      // TODO: if reading from json file, then we need to set the prefix of
      // CalculationParameter   , first attempt is to modify in ParameterManager
      // ctor
      std::string prefix = pm.prefix();

      wf.run(prefix);

      if (false) {
        qcapp::Workflow opt_wf;

        std::function<std::unique_ptr<Application>(Params)> scfFactory =
            [&](Params p) {
              return std::make_unique<SCFApplication<moldft_lib>>(world, p);
            };

        opt_wf.addDriver(
            std::make_unique<qcapp::OptimizeDriver>(world, scfFactory, pm));
      }
    } catch (std::exception& e) {
      if (world.rank() == 0) {
        print_header2("caught an exception in the main loop");
        print(e.what());
        print_header2("ending program run");
      }
    }

    // Nearly all memory will be freed at this point
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }  // world is dead -- ready to finalize
  finalize();

  return 0;
}
