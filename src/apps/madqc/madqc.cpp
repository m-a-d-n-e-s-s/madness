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

#include "madchem.h"
#include "parameters_manager.hpp"
#include "tasks.hpp"

using namespace madness;

int main(int argc, char** argv) {

  World& world = initialize(argc, argv, false);
  if (world.rank() == 0) {
    print_header1("MADQC -- numerical quantum chemistry in MADNESS");
  }

  startup(world, argc, argv, true);
  std::cout.precision(6);
  if (world.rank() == 0)
    print(info::print_revision_information());
  commandlineparser parser(argc, argv);
  ParameterManager params;

  try {
    
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
      if (world.rank() == 0) {
        print("Input file found");
        print("Parsing Command Line");
      }
      params = ParameterManager(world, input_file);
    } else if (argc == 3) {
      if (world.rank() == 0) {
        print("Input and mol file found");
      }
      path input_file(argv[1]);
      path mol_input(argv[2]);
      params = ParameterManager(world, {input_file, mol_input});
    } else {
      error("Too many arguments");
    }

    std::ofstream out_file;
    out_file.open("madqc_input.json");
    params.write_input_file(out_file);
    out_file.close();

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
  // create parameter classes
  // 1. read in all input blocks independently
  // 2. set up parameter logic
  // 2a from the model downstream
  // 2b from the task downstream

  // read input file
  // read into parameter handler

  // create class corresponding to qc model

  // check for the existence of the input file

  finalize();
  return 0;
}
