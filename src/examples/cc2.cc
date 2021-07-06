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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/tdhf.cc
  \brief compute the time-dependent HF equations (currently CIS approximation)

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

 */
#include <chem/CC2.h>

using namespace madness;

#ifdef USE_GENTENSOR

int main(int argc, char** argv) {

	initialize(argc, argv);

	World world(SafeMPI::COMM_WORLD);

	if (world.rank() == 0) {
		print("\n  CC2: Coupled Cluster approximate Doubles  \n");
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);

	}
	startup(world,argc,argv);
	std::cout.precision(6);
	FunctionDefaults<3>::set_truncate_mode(1);
	if(world.rank()==0) print("Truncate mode set to ",FunctionDefaults<3>::get_truncate_mode());

#ifdef GITREVISION
	const  char* gitrev =  GITREVISION;
	const std::string gitrevision(gitrev);
	if (world.rank()==0) {
		print("           git revision ...",gitrevision);
	}
#endif

typedef std::vector<real_function_3d> vecfuncT;

// set the tensor type
TensorType tt=TT_2D;
FunctionDefaults<6>::set_tensor_type(tt);
FunctionDefaults<6>::set_apply_randomize(true);

// Process 0 reads input information and broadcasts
const char * inpname = "input";
for (int i=1; i<argc; i++) {
    if (argv[i][0] != '-') {
        inpname = argv[i];
        break;
    }
}
if (world.rank() == 0) print("input filename: ", inpname);

//SCF calc(world,input.c_str());
std::shared_ptr<SCF> calc(new SCF(world,inpname));
std::shared_ptr<Nemo> nemo(new Nemo(world,calc,inpname));
if (world.rank()==0) {
    calc->molecule.print();
    print("\n");
    calc->param.print("cc2");
}
double hf_energy = nemo->value();
if(world.rank()==0) std::cout << "\n\n\n\n\n\n Reference Calclation Ended\n SCF Energy is: " << hf_energy
		<<"\n current wall-time: " << wall_time()
		<<"\n current cpu-time: " << cpu_time()<< "\n\n\n";


// Make CC2
CC2 cc2(world,inpname,nemo);

cc2.solve();

if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
world.gop.fence();
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





