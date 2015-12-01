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
  \file examples/nemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

 */
#include "smooth.h"


using namespace madness;


int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) {
		print("\n  Density Smoothing for XC kernel \n");
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);

	}
	startup(world,argc,argv);
	std::cout.precision(6);

#ifdef MADNESS_GITREVISION
	const  char* gitrev =  MADNESS_GITREVISION;
	const std::string gitrevision(gitrev);
	if (world.rank()==0) {
		print("           git revision ...",gitrevision);
	}
#endif

	try {

		if(world.rank()==0){
			std::cout << "Converge Orbitals with Nemo" << std::endl;
		}

		const std::string input="input";
		std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
		if (world.rank()==0) {
			calc->molecule.print();
			print("\n");
			calc->param.print(world);
		}

		// compute the energy to get converged orbitals
		Nemo nemo(world,calc);
		const double energy=nemo.value();
		if (world.rank()==0) {
			printf("final energy   %12.8f\n", energy);
			printf("finished at time %.1f\n", wall_time());
		}

		if(world.rank()==0){
			std::cout << "\n\nSmooth Density:" << std::endl;
			std::cout << "3D Threshold is " << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout << "\n\n";
		}

		smooth smoothing(world,nemo);
		//smoothing.test_1d();
		//smoothing.test();
		//smoothing.make_smooth_gradient();
		//smoothing.make_smooth_XC_kernel(nemo.get_calc()->param.xc_data);
		smoothing.apply_smooth_slater_kernel();




	} catch (const SafeMPI::Exception& e) {
		print(e);
		error("caught an MPI exception");
	} catch (const madness::MadnessException& e) {
		print(e);
		error("caught a MADNESS exception");
	} catch (const madness::TensorException& e) {
		print(e);
		error("caught a Tensor exception");
	} catch (char* s) {
		print(s);
		error("caught a string exception");
	} catch (const char* s) {
		print(s);
		error("caught a string exception");
	} catch (const std::string& s) {
		print(s);
		error("caught a string (class) exception");
	} catch (const std::exception& e) {
		print(e.what());
		error("caught an STL exception");
	} catch (...) {
		error("caught unhandled exception");
	}


	finalize();
	return 0;
}
