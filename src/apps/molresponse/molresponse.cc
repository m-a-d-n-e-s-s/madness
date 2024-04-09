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
#include <madness/chem/SCF.h>
#include <madness/world/worldmem.h>

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
// #include "ResponseExceptions.hpp"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>

static inline int file_exists(const char *inpname) {
    struct stat buffer {};
    size_t rc = stat(inpname, &buffer);
    return (rc == 0);
}

#endif

using namespace madness;

int main(int argc, char **argv) {
    World &world = madness::initialize(argc, argv);
    startup(world, argc, argv, true);
    if (world.rank() == 0) {
        print_header1("MOLRESPONSE -- MADNESS Time-Dependent Density Functional "
                      "Theory Excited-State Program ");
    }

    //    sleep(10);
    int result = 0;
    std::cout.precision(6);
    std::string filename = "response.in";

    commandlineparser parser(argc, argv);

    if (parser.key_exists("help")) {
        FrequencyResponse::help();
    } else if (parser.key_exists("print_parameters")) {
        FrequencyResponse::print_parameters();
    } else {

        molresponse::start_timer(world);
        // try catch would start here
        try {
            auto calc_params = initialize_calc_params(world, filename);
            if (calc_params.response_parameters.excited_state()) {

                ExcitedResponse calc(world, calc_params);
                if (world.rank() == 0) {
                    print("\n\n");
                    print(" MADNESS Time-Dependent Density Functional Theory Excited-State "
                          "Program");
                    print(" ----------------------------------------------------------\n");
                    print("\n");
                    calc_params.molecule.print();
                    print("\n");
                    calc_params.response_parameters.print("response");
                    // put the response parameters in a j_molrespone json object
                }
                calc_params.response_parameters.to_json(calc.j_molresponse);
                // set protocol to the first
                calc.solve(world);
                calc.output_json();
            } else if (calc_params.response_parameters.first_order()) {
                RHS_Generator rhsGenerator;
                if (calc_params.response_parameters.dipole()) {
                    rhsGenerator = dipole_generator;
                } else if (calc_params.response_parameters.nuclear()) {
                    rhsGenerator = nuclear_generator;
                }
                auto omega = calc_params.response_parameters.omega();
                FrequencyResponse calc(world, calc_params, omega, rhsGenerator);
                // Warm and fuzzy for the user
                if (world.rank() == 0) {
                    print("\n\n");
                    print(" MADNESS Time-Dependent Density Functional Theory Frequency "
                          "Response "
                          "Program");
                    print(" ----------------------------------------------------------\n");
                    print("\n");
                    calc_params.molecule.print();
                    print("\n");
                    calc_params.response_parameters.print("response");
                }
                calc_params.response_parameters.to_json(calc.j_molresponse);
                // set protocol to the first
                calc.solve(world);
                calc.output_json();
            } else {
                if (world.rank() == 0) { print("Response not implemented"); }
            }

        } catch (const SafeMPI::Exception &e) { print(e); } catch (const madness::MadnessException &e) {
            std::cout << e << std::endl;
        } catch (const madness::TensorException &e) { print(e); } catch (const char *s) {
            print(s);
        } catch (const std::string &s) { print(s); } catch (const std::exception &e) {
            print(e.what());
        } catch (...) { error("caught unhandled exception"); }
    }

    finalize();
    return result;
}
