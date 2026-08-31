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

#include <madness/misc/info.h>
#include <madness/world/worldmem.h>

#include <fstream>
#include <iostream>

#include <Drivers.hpp>
#include <ParameterManager.hpp>
#include <WorkflowBuilders.hpp>
#include <madness/chem/ResultsSummary.hpp>
#include <madness/chem/VizManifest.hpp>
#include <apps/molresponse/madqc_adapter.hpp>  // molresponse_v3_lib (R3 engine)
#include <madness_exception.h>

using namespace madness;

/// TODO:
///  - restart options
///  - nemo projections
///  - numerical parameters

void help(const std::string &wf) {
  print("Usage: madqc [options] [input_file]");
  print("\nOptions:");
  print("  --help=<workflow>           : show this help message");
  print("  --print_parameters=<group>  : print all parameters and exit");
  print("  --restart_info=<prefix>     : describe a restart archive and exit");
  print("  --workflow=<name>           : specify the workflow to run (default: "
        "scf)");
  print("  --optimize                  : optimize the geometry of the "
        "workflow's reference (scf|nemo)");
  print("\nAvailable workflows: " + workflow_builders::runnable_workflow_list());
  print("Parameter groups (for --print_parameters): dft, nemo, response, cc2, "
        "cis, oep, optimization, geometry");
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
    print("(runs the ground-state SCF first, then linear/quadratic response)");
    print("\nexamples: ");
    print("madqc --wf=response h2o_response.in");
    print("see  madqc --print_parameters=response  for all response knobs");
  } else if (wf == "mp2" or wf == "cc2") {
    print("madqc --wf=cc2");
    print("MP2/CC2 code for correlated wave functions");
    print("\nexamples: ");
    print("madqc --wf=cc2 --geometry=h2o --dft=\"k=5; econv=1.e-4; L=30\" "
          "--cc2=\"freeze=1\"");
  } else if (wf == "cis") {
    print("madqc --wf=cis");
    print("CIS code for excited states");
  } else if (wf == "oep") {
    print("madqc --wf=oep");
    print("Optimized effective potential code for DFT");
    print("\nexamples: ");
    print("madqc --wf=oep --geometry=h2o");
  } else if (wf == "optimize") {
    print("madqc --optimize --wf=<scf|nemo>");
    print("Geometry optimization: --wf names the reference method, --optimize "
          "asks for its");
    print("geometry to be optimized. (`--wf=optimize` is not a workflow.)");
    print("\nexamples: ");
    print("madqc --optimize --wf=scf --geometry=h2o");
    print("madqc --optimize --wf=nemo --geometry=lih");
    print("madqc --optimize --wf=scf --geometry=h2o --optimization=\"gtol=1.e-4; "
          "maxiter=10\"");
    print("\nsee  madqc --print_parameters=optimization  for all knobs.");
    print("the in-SCF form (--dft=\"gopt=1\") has been removed; its keyvals are");
    print("retired and error with a pointer to the `optimization` group.");
  }
}

void print_parameters(World &world, const commandlineparser &parser,
                      const std::string &group) {
  Params pm;
  if (group.empty()) {
    print("please specify a data group to print parameters for");
    print("\n  --print_parameters=<group>  : print all parameters and exit");
    print("\nAvailable data groups: dft, nemo, response, cc2, cis, oep, "
          "geometry");
  } else if (group == "dft" || group == "scf") {
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
  } else if (group == "optimization") {
    print("Available parameters for data group: optimization");
    pm.get<OptimizationParameters>().print(OptimizationParameters::tag, "end");
  } else if (group == "geometry") {
    Molecule::GeometryParameters geometryparam;
    geometryparam.print("geometry", "end");
  } else {
    std::string msg = "Unknown data group: " + group +
                      "\nAvailable data group are: dft, nemo, response, cc2, "
                      "cis, oep, optimization, geometry\n";
    print(msg);
  }
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  // exit code: a run that died in the main loop must not report success --
  // CI and driver scripts have no other way to tell the difference
  int rc = 0;
  if (world.rank() == 0) {
    print_header1("MADQC -- Multiresolution Quantum Chemistry Code ");
  }

  commandlineparser parser(argc, argv);

  if (parser.key_exists("help")) {
    help(parser.value("help"));
  } else if (parser.key_exists("print_parameters")) {
    print_parameters(world, parser, parser.value("print_parameters"));
  } else if (parser.key_exists("restart_info")) {
    // Report what a restart archive holds and stop. No startup() needed: the
    // header and the orbital bookkeeping are plain Molecule/Tensor data, and no
    // MRA function is constructed.
    std::string p = parser.value("restart_info");
    if (p.empty() or p == "restart_info") p = "mad";
    madness::print_restartdata_info(world, p);
  } else {
    // limit lifetime of world so that finalize() can execute cleanly
    try {
      // Load info for MADNESS numerical routines
      startup(world, argc, argv, true);
      if (world.rank() == 0)
        print(info::print_revision_information());

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

      // Deck-level `io` block (roadmap change 5): run-wide restart-backend
      // selection. `io backend hdf5` turns the HDF5 opt-in on for every task
      // that supports it; the response-only `response.hdf5 true` keeps
      // working as an alias (applied here and re-checked in the response
      // adapter). The effective config is stamped into calc_info for ALL
      // tasks via Workflow::set_provenance below.
      // HDF5 restart I/O is a molresponse-only path today, so only honor the
      // deck flag (and stamp the provenance below) for the response workflow —
      // a stray `response.hdf5`/`io backend hdf5` in a cc2/cis/oep/moldft deck
      // must NOT flip the global flag or mis-stamp backend=hdf5 (review LOW).
      // `--optimize` turns the workflow's reference into an optimization task
      // (WorkflowBuilders::add_optimize_workflow_drivers).
      const bool optimize_geometry = parser.key_exists("optimize");
      const bool is_response_wf =
          workflow_builders::workflow_kind_from_name(user_workflow) ==
          workflow_builders::WorkflowKind::Response;
      const bool deck_io_hdf5 =
          pm.get<IOParameters>().hdf5() || pm.get<ResponseParameters>().hdf5();
      if (deck_io_hdf5 && is_response_wf) {
#ifdef MADNESS_HAS_HDF5
        molresponse_v3::set_hdf5_io_enabled(true);
#else
        throw std::runtime_error(
            "io backend hdf5 (or response.hdf5) requested but this build has "
            "no HDF5 support — configure with -DMADNESS_ENABLE_HDF5=ON");
#endif
      } else if (deck_io_hdf5 && !is_response_wf && world.rank() == 0) {
        print("NOTE: `io backend hdf5` / `response.hdf5` only affects the "
              "response workflow's MRA restart I/O — ignored for --wf=",
              user_workflow);
      }

      // The response workflow is built HERE in the app — MADchem's
      // WorkflowBuilders cannot reference v3 (circular lib dependency), and the
      // v2 engine was removed (M1 decoupling Stage 2). molresponse_v3 is THE
      // response engine; `engine v2` in old decks fails loudly below.
      if (is_response_wf && optimize_geometry) {
        // The response step would have to run at the optimized geometry, which
        // needs the downstream consumer work of ARCHITECTURE_ROADMAP change 2.
        throw std::runtime_error(
            "--optimize is not supported for --wf=response; optimize the "
            "geometry first (--optimize --wf=scf) and run the response "
            "calculation on the result");
      }

      if (is_response_wf) {
        const std::string engine = pm.get<ResponseParameters>().engine();
        if (engine == "v2") {
          if (world.rank() == 0)
            print("ERROR: response engine 'v2' (MolresponseLib) was removed.\n"
                  "       Delete the `engine v2` line (v3 is the default) — the\n"
                  "       input deck is otherwise unchanged; see\n"
                  "       molresponse_v3/MIGRATION_FROM_V2.md.");
          throw std::runtime_error("response engine v2 removed — use v3");
        }
        if (world.rank() == 0)
          print("response engine : molresponse_v3");
        pm.get<CalculationParameters>().set_derived_value("save", true);
        auto reference =
            std::make_shared<SCFApplication<moldft_lib>>(world, pm);
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
        wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
            std::make_unique<ResponseApplication<molresponse_v3_lib>>(
                world, pm, reference->calc())));
      } else {
        workflow_builders::add_workflow_drivers(world, pm, user_workflow, wf,
                                                optimize_geometry);
      }

      // io provenance for ALL tasks (today only response writes MRA restart
      // state, but the stamp is run-level: an HDF5 run and a native run must
      // be distinguishable from any task's calc_info).
#ifdef MADNESS_HAS_HDF5
      wf.set_provenance(
          "io", {{"backend",
                  molresponse_v3::hdf5_io_enabled() ? "hdf5" : "native"},
                 {"hdf5_compiled", true}});
#else
      wf.set_provenance("io",
                        {{"backend", "native"}, {"hdf5_compiled", false}});
#endif

      // Echo the effective input of the whole run before any work starts: the
      // molecule, then each driver's own parameter group. This is the only
      // record of what was actually computed (defaults and derived values
      // included) for a run whose deck sets just a handful of keys.
      if (world.rank() == 0) {
        print_header1("Calculation parameters");
        pm.get<Molecule>().print();
        wf.print_parameters(world);
        print("");
        print_header1("Starting calculations");
      }

      std::string prefix = pm.prefix();
      wf.run(prefix);

      // Human-readable chemistry report: to stdout and <prefix>.out. The
      // machine-readable <prefix>.calc_info.json remains the source of truth.
      if (world.rank() == 0) {
        qcapp::write_results_summary(std::cout, wf.results());
        std::ofstream report(prefix + ".out");
        qcapp::write_results_summary(report, wf.results());
        print("Wrote results summary :", prefix + ".out");

        // Index any visualization artifacts (cube/dx/line plots) the run
        // produced into <prefix>.viz_manifest.json for gecko/ParaView/VMD.
        qcapp::write_viz_manifest(prefix, wf.results());
      }
    } catch (const MadnessException &e) {
      rc = 1;
      if (world.rank() == 0) {
        print_header2("caught a MADNESS exception in the main loop");
        print(e.what(), e.filename, e.msg, e.line);
        print_header2("ending program run");
      }
    } catch (const json::exception &e) {
      rc = 1;
      if (world.rank() == 0) {
        print_header2("caught a JSON exception in the main loop");
        print(e.what());
        print_header2("ending program run");
      }
    }

    catch (std::exception &e) {
      rc = 1;
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
  } // world is dead -- ready to finalize
  finalize();

  return rc;
}
