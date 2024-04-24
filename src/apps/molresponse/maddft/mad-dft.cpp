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

/// \file mad-dft.cc
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density functional and Hartree-Fock code

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && \
    defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <unistd.h>

static inline int file_exists(const char* input_name) {
  struct stat buffer {};
  size_t rc = stat(input_name, &buffer);
  return (rc == 0);
}

#endif
#include <madness/misc/info.h>
#include "maddft/response_manager.hpp"

using path = std::filesystem::path;

using madness::World;

auto main(int argc, char* argv[]) -> int {
  madness::initialize(argc, argv);
  std::cout.precision(6);

  {
    World world(SafeMPI::COMM_WORLD);
    ParameterManager params;
    if (world.rank() == 0) {
      print_header1(
          "MOLDFT -- molecular DFT and Hartree-Fock code with TDHF and TDDFT");
    }
    startup(world, argc, argv, true);
    if (world.rank() == 0)
      print(info::print_revision_information());
    commandlineparser parser(argc, argv);

    if (parser.key_exists("help")) {
      ParameterManager::help();

    } else if (parser.key_exists("print_parameters")) {
      params.print_params();
    } else {
      if (parser.key_exists("dft")) {
        print("Running DFT");
      }
      if (parser.key_exists("response")) {
        print("Running Response");
      }
      if (parser.key_exists("quadratic")) {
        print("Running Quadratic Response");
      }

      try {
        // I need to write a help and a print parameters function which will be
        // called by the commandlineparser
        print_meminfo(world.rank(), "startup");

        ParameterManager params;
        if (argc == 1) {
          if (world.rank() == 0) {
            print("No input file found");
            print("For help type: ./mad-dft --help");
            print("For print parameters type: ./mad-dft --print_parameters");
          }
          return 1;

        } else if (argc == 2) {
          path input_file(argv[1]);
          commandlineparser parser(argc, argv);
          if (world.rank() == 0) {
            print("Input file found");
            print("Parsing Command Line");
          }
          parser.print_map();
          params = ParameterManager(world, parser);
        } else if (argc == 3) {
          if (world.rank() == 0) {
            print("Input and mol file found");
          }
          path input_file(argv[1]);
          path mol_input(argv[2]);
          params = ParameterManager(world, input_file, mol_input);
        } else {
          error("Too many arguments");
        }

        // print params
        auto response_manager = ResponseCalcManager(world, params);

        if (world.rank() == 0) {
          print("Running MOLDFT");
          print("Calc Info Path: ", response_manager.calc_info_json_path);
          print("Moldft Path: ", response_manager.moldft_path);
        }
        response_manager.output_calc_path_json();

        if (params.get_run_moldft()) {
          if (std::filesystem::exists(response_manager.calc_info_json_path) &&
              std::filesystem::exists(response_manager.moldft_path)) {
            print(" MOLDFT and Calc Info Found... Skipping MOLDFT");
          } else {
            response_manager.run_molresponse(world);
          }
        }
        // if both exist, read the calc_info json

        //
        if (std::filesystem::exists(response_manager.calc_info_json_path) &&
            std::filesystem::exists(response_manager.moldft_path)) {
          response_manager.run_molresponse(world);
        } else {
          if (world.rank() == 0) {
            print("Running MOLDFT since no previous calculation was found");
          }
          response_manager.run_moldft(world, true);
          world.gop.fence();

          response_manager.run_molresponse(world);
          world.gop.fence();
        }

        // if quadratic response is requested
        //
        if (params.get_molresponse_params().quadratic()) {
          if (world.rank() == 0) {
            print("Compute Quadratic Response Properties ");
          }
          if (!std::filesystem::exists(response_manager.quadratic_json_path)) {
            response_manager.run_quadratic_response(world);
          } else {
            print(
                "Quadratic Response Properties Found... Skipping Quadratic "
                "Response");
          }
        }
      } catch (const SafeMPI::Exception& e) {
        print(e.what());
        error("caught an MPI exception");
      } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
      } catch (const madness::TensorException& e) {
        print(e.what());
        error("caught a Tensor exception");
      } catch (const nlohmann::detail::exception& e) {
        print(e.what());
        error("Caught JSON exception");
      } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << ex.what() << "\n";
      } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
      } catch (...) {
        error("caught unhandled exception");
      }
      // Nearly all memory will be freed at this point
      print_stats(world);
      if (world.rank() == 0) {
        print("Finished All Frequencies");
      }
    }
  }
  finalize();
  return 0;
}
