/*
 * cc2.cc
 *
 *  Created on: Aug 10, 2015
 *      Author: kottmanj
 */
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

#include "madness/chem/CC2.h"
#include "madness/misc/info.h"

using namespace madness;

#ifdef USE_GENTENSOR

int main(int argc, char **argv) {

    World& world=initialize(argc, argv);

    if (world.rank() == 0) {
        print_header1("CC2: Coupled Cluster approximate Doubles");
        printf("starting at time %.1f\n", wall_time());
    }

    printf_msg_energy_time("message %8.4f %3.2fs",0.0,wall_time());

    startup(world, argc, argv,true);
    std::cout.precision(6);
    FunctionDefaults<3>::set_truncate_mode(1);
    if (world.rank()==0) print(info::print_revision_information());

    // set the tensor type
    TensorType tt = TT_2D;
    FunctionDefaults<6>::set_tensor_type(tt);
//    FunctionDefaults<6>::set_apply_randomize(true);

    commandlineparser parser(argc, argv);
    if (parser.key_exists("help")) {
        CC2::help();

    } else if (parser.key_exists("print_parameters")) {
        CC2::print_parameters();

    } else {

        std::shared_ptr<Nemo> nemo(new Nemo(world, parser));
        nemo->param.set_derived_value("print_level", 2);
        nemo->get_calc()->param.set_derived_value("print_level", 2);
        nemo->param.set_derived_value("k", 5);
        nemo->get_calc()->param.set_derived_value("k", 5);
        // nemo->param.set_derived_value<std::string>("localize", "canon");
        // nemo->get_calc()->param.set_derived_value<std::string>("localize", "canon");
        nemo->param.set_derived_values(nemo->molecule());
        nemo->get_calc()->param.set_derived_values(nemo->molecule());
        CC2 cc2(world, parser, nemo);

        std::shared_ptr<SCF> calc = nemo->get_calc();
        if (world.rank() == 0) {
            print("\n");
            cc2.parameters.print("cc2","end");
            print("\n");
            calc->param.print("dft","end");
            print("\n");
            cc2.tdhf->get_parameters().print("response","end");
            print("\n");
            nemo->molecule().print();
        }
        double hf_energy = nemo->value();
        try {
            cc2.solve();
        } catch (std::exception& e) {
            print("Caught exception: ", e.what());
        }

        if (world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
        world.gop.fence();
    }
    finalize();

    return 0;
}// end main

#else
int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    if(world.rank() == 0) {

        print("\nYou can't run cc2 because you have configured MADNESS ");
        print("without the --enable-gentensor flag");
        print("You need to reconfigure and recompile\n");

    }
    world.gop.fence();
    finalize();

    return 0;

}
#endif





