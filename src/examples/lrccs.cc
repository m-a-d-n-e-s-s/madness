/*
 * lrccs.cc
 *
 *  Created on: 4 Jan 2017
 *      Author: kottmanj
 */

/*
/*
 * lrccs.cc
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
#include <chem/TDHF.h>

using namespace madness;

int main(int argc, char** argv) {

	initialize(argc, argv);

	World world(SafeMPI::COMM_WORLD);

	// read out command line arguments
	bool analyze_only=false;		// default: compute the excitations

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];
        if (arg=="--analyze") analyze_only=true;
    }


	if (world.rank() == 0) {
		print("\n  CC2 without Gentensor Flag : Only Linear Response For CCS/CIS  \n");
		print("!!!DO NOT USE 6D APPLICATIONS HERE!!!");
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);

	}
	startup(world,argc,argv);
	std::cout.precision(6);
	FunctionDefaults<3>::set_truncate_mode(1);
	print("Truncate mode set to ",FunctionDefaults<3>::get_truncate_mode());

#ifdef GITREVISION
	const  char* gitrev =  GITREVISION;
	const std::string gitrevision(gitrev);
	if (world.rank()==0) {
		print("           git revision ...",gitrevision);
	}
#endif

	typedef std::vector<real_function_3d> vecfuncT;

	// Make reference
	const std::string input = "input";
	//SCF calc(world,input.c_str());
	std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
	Nemo nemo(world,calc);
	if (world.rank()==0) {
		calc->molecule.print();
		print("\n");
		calc->param.print(world);
	}
	double hf_energy = nemo.value();
	if(world.rank()==0) std::cout << "\n\n\n\n\n\n Reference Calclation Ended\n SCF Energy is: " << hf_energy
			<<"\n current wall-time: " << wall_time()
			<<"\n current cpu-time: " << cpu_time()<< "\n\n\n";

	CCParameters parameters(input,nemo.get_calc()->param.lo);
	if (analyze_only) parameters.no_compute_response=true;

	if(world.rank()==0) std::cout << "Setting 3D Thresh to " << parameters.thresh_3D << "\n";
	FunctionDefaults<3>::set_thresh(parameters.thresh_3D);

	TDHF tdhf(world,parameters,nemo);

	// try to restart
	std::vector<CC_vecfunction> roots = tdhf.read_vectors();

	// solve the CIS equations
	if (not analyze_only) tdhf.solve_cis(roots);

	// analyze the results
	tdhf.analyze(roots);

	if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
	world.gop.fence();
	finalize();

	return 0;
}// end main


