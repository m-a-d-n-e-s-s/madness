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
#include<madness/chem/SCF.h>
#include <madness/world/worldmem.h>
#include "ResponseExceptions.hpp"

#include "TDDFT.h"  // All response functions/objects enter through this
#include "molresponse/density.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <unistd.h>

static inline int file_exists(const char *inpname) {
  struct stat buffer {};
  size_t rc = stat(inpname, &buffer);
  return (rc == 0);
}

#endif



using namespace madness;

int main(int argc, char **argv) {
  initialize(argc, argv);
  {  // limit lifetime of world so that finalize() can execute cleanly
    World world(SafeMPI::COMM_WORLD);
    molresponse::start_timer(world);
    // try catch would start here
    try {
      startup(world, argc, argv, true);
      print_meminfo(world.rank(), "startup");
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

      std::cout.precision(6);
      // This makes a default input file name of 'input'
      const char *input_file = "response.in";
      // Process 0 reads input information and broadcasts
      for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
          input_file = argv[i];
          break;
        }
      }

      if (world.rank() == 0) print("input filename: ", input_file);
      if (!file_exists(input_file)) throw Input_Error{};

      auto [ground_calculation, molecule, response_parameters] =
      initialize_calc_params(world, std::string(input_file));
      vecfuncT ground_orbitals = ground_calculation.orbitals();

      density_vector rho = set_density_type(world, response_parameters, ground_calculation);
      // first step is to read the input for r_params and g_params
      // Create the TDDFT object
      TDDFT calc = TDDFT(world, rho);

      // Warm and fuzzy for the user
      if (world.rank() == 0) {
        print("\n\n");
        print(
            " MADNESS Time-Dependent Density Functional Theory Response "
            "Program");
        print(" ----------------------------------------------------------\n");
        print("\n");
        calc.molecule.print();
        print("\n");
        calc.r_params.print("response");
        // put the response parameters in a j_molrespone json object
        calc.r_params.to_json(calc.j_molresponse);
      }
      molresponse::end_timer(world, "initialize");
      // Come up with an initial OK data map
      if (world.size() > 1) {
        calc.set_protocol<3>(world, 1e-4);
        calc.make_nuclear_potential(world);
        calc.initial_load_bal(world);
      }
      // set protocol to the first
      if (calc.r_params.excited_state()) {
        calc.solve_excited_states(world);
      } else if (calc.r_params.first_order()) {
        calc.solve_response_states(world);
      } else if (calc.r_params.second_order()) {
      } else {
        throw Response_Input_Error{};
      }
      calc.output_json();
      if (calc.r_params.dipole()) {  //
        print("Computing Alpha");
        Tensor<double> alpha = calc.polarizability();
        print("Second Order Analysis");
        calc.PrintPolarizabilityAnalysis(world, alpha);
      }
    } catch (const Input_Error &e) {
      print(e);
      error("Input File Error");
    } catch (const SafeMPI::Exception &e) {
      print(e);
      error("caught an MPI exception");
    } catch (const madness::MadnessException &e) {
      print(e);
      error("caught a MADNESS exception");
    } catch (const madness::TensorException &e) {
      print(e);
      error("caught a Tensor exception");
    } catch (const char *s) {
      print(s);
      error("caught a string exception");
    } catch (const std::string &s) {
      print(s);
      error("caught a string (class) exception");
    } catch (const std::exception &e) {
      print(e.what());
      error("caught an STL exception");
    } catch (...) {
      error("caught unhandled exception");
    }
    world.gop.fence();
    print_stats(world);
    finalize();
  }

  return 0;
}
