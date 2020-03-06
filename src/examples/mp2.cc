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


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/info.h>
#include <chem/mp2.h>
#include <madness/misc/gitinfo.h>

using namespace madness;

#ifdef USE_GENTENSOR

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    if (world.rank()==0) print(info::print_revision_information());


    if (world.rank()==0) {
#ifdef MADNESS_HAS_GOOGLE_PERF_MINIMAL
    	print("using gperftools, clearing memory at each fence()");
#endif
    }

    TensorType tt=TT_2D;

    // get command line parameters (overrides input file)
    bool do_test=false;
    std::string testfilename;
    for(int ii = 1; ii < argc; ii++) {
        const std::string arg=argv[ii];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="test") {
        	do_test=true;
        	testfilename=val;
        }
        if (key=="TT") {
            if (val=="TT_2D") tt=TT_2D;
            if (val=="TT_TENSORTRAIN") tt=TT_TENSORTRAIN;
        }
    }

    FunctionDefaults<6>::set_tensor_type(tt);
    FunctionDefaults<6>::set_apply_randomize(true);

    try {
    	MP2 mp2(world,"input");

    	if(world.rank() == 0) printf("\nstarting at time %.1fs\n", wall_time());

		if (do_test) mp2.test(testfilename);
		else {
			const double hf_energy=mp2.get_hf().value();
			const double mp2_energy=mp2.value();
			if(world.rank() == 0) {
				printf("final hf/mp2/total energy %12.8f %12.8f %12.8f\n",
						hf_energy,mp2_energy,hf_energy+mp2_energy);
			}
		}
    } catch (std::exception& e) {

    	if (world.rank()==0) {
    		print("\ncaught an exception: \n",e.what());
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
