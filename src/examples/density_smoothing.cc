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
  \file examples/density_smoothing.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/density_smoothing.cc>here</a>.

 */
#include "smooth.h"
#include"madness/mra/commandlineparser.h"
#include<madness/misc/info.h>


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

	if (world.rank()==0) {
		print("           git source description ...", info::git_source_description());
	}

	try {

		if(world.rank()==0){
			std::cout << "Converge Orbitals with Nemo" << std::endl;
		}

		commandlineparser parser(argc,argv);
        Nemo nemo(world,parser);
        auto calc=nemo.get_calc();
		if (world.rank()==0) {
			calc->molecule.print();
			print("\n");
			calc->param.print("dft");
		}

		// compute the energy to get converged orbitals
		const double energy=nemo.value();
		if (world.rank()==0) {
			printf("final energy   %12.8f\n", energy);
			printf("finished at time %.1f\n", wall_time());
		}

		if(world.rank()==0){
			std::cout << "\n\nSmooth Density:" << std::endl;
			std::cout << "Orbital energie \n";
			std::cout << "3D Threshold is " << FunctionDefaults<3>::get_thresh() << std::endl;
			std::cout << "\n\n";
		}

		const size_t norb = nemo.get_calc()->amo.size();
		std::vector<double> eps;
		for(size_t i=0;i<norb;i++){
			eps.push_back(nemo.get_calc()->aeps(i));
		}
		using madness::operators::operator<<;
		std::cout << "Orbital energies\n " << eps << std::endl;

		real_function_3d density = real_factory_3d(world);
		double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
		real_function_3d R2 = nemo.ncf->square();
		for(auto x:nemo.get_calc()->amo){
			//x.truncate();
			x.refine();
			x.scale(-1.0);
			real_function_3d xsquare = x*x;
			plot_plane(world,xsquare,"xsquare");
			plot_plane(world,x,"xorbital");
			plot_plane(world,x,"R2");
			coord_3d start(0.0); start[0]=-width;
			coord_3d end(0.0); end[0]=width;
			const std::string line = "line";
			plot_line((line+"_xorbital").c_str(),1000,start,end,x);
			plot_line((line + "_xsquare").c_str(),1000,start,end,xsquare);
			plot_line((line+"_R2").c_str(),1000,start,end,R2);
			density+=xsquare;
//			real_function_3d amoR = x*nemo.nuclear_correlation->square();
//			density += amoR*x;
		}
		density=density*nemo.ncf->square();
		density.truncate(FunctionDefaults<3>::get_thresh()*0.01);
		smooth<double,3> smoothing(world);
		smoothing.smooth_density_from_orbitals(nemo.get_calc()->amo);
		//smoothing.set_molecule_mask(1.0,1.0,nemo.get_calc()->molecule.get_atoms());





	} catch (const SafeMPI::Exception& e) {
		print(e);
		error("caught an MPI exception");
	} catch (const madness::MadnessException& e) {
		print(e);
		error("caught a MADNESS exception");
	} catch (const madness::TensorException& e) {
		print(e);
		error("caught a Tensor exception");
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
