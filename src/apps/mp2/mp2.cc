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
/*!
  \file mp2.cc
  \brief Solves molecular MP2 equations
  \defgroup Solves molecular MP2 equations
  \ingroup examples
*/


#include "madness/chem/mp2.h"
#include "madness/misc/info.h"

using namespace madness;

#ifdef USE_GENTENSOR

int main(int argc, char** argv) {
    World& world=initialize(argc, argv,false);
    if (world.rank() == 0) {
        print_header1("MP2 -- second order correlation energies");
        printf("starting at time %.1f\n", wall_time());
    }

    startup(world,argc,argv,true);
    std::cout.precision(6);
    if (world.rank()==0) print(info::print_revision_information());


    commandlineparser parser(argc,argv);
    if (parser.key_exists("help")) {
        MP2::help();

    } else if (parser.key_exists("print_parameters")) {
        MP2::print_parameters();

    } else {

        TensorType tt = TT_2D;
        if (parser.key_exists("TT")) {
            if (parser.value("TT") == "TT_2D") tt = TT_2D;
            if (parser.value("TT") == "TT_TENSORTRAIN") tt = TT_TENSORTRAIN;
        }

        FunctionDefaults<6>::set_tensor_type(tt);
        FunctionDefaults<6>::set_apply_randomize(true);

        try {
            MP2 mp2(world, parser);

            if (world.rank() == 0) printf("\nstarting at time %.1fs\n", wall_time());

			const double hf_energy=mp2.get_hf().value();
			const double mp2_energy=mp2.value();
//            const double mp2_energy=0.0;
			if(world.rank() == 0) {
				printf("final hf/mp2/total energy %12.8f %12.8f %12.8f\n",
						hf_energy,mp2_energy,hf_energy+mp2_energy);
			}
        } catch (std::exception& e) {

            if (world.rank() == 0) {
                print("\ncaught an exception: \n", e.what());
            }
        }
    }

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}

#else


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    if(world.rank() == 0) {

    	print("\nYou can't run mp2 because you have configured MADNESS ");
    	print("without the --enable-gentensor flag");
    	print("You need to reconfigure and recompile\n");

    }
    world.gop.fence();
    finalize();

    return 0;

}
#endif
