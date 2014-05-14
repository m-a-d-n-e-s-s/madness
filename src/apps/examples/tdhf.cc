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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/tdhf.cc
  \brief compute the time-dependent HF equations (currently CIS approximation)

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

*/

#include <examples/tdhf_CIS.h>
//#include <examples/TD.h>

#include<iomanip>
#include<iostream>

using namespace madness;

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  TDHF -- time-dependent Hartree-Fock in the CIS approximation  \n");
    	printf("starting at time %.1f\n", wall_time());
       	print("\nmain() compiled at ",__TIME__," on ",__DATE__);

    }
    startup(world,argc,argv);
    std::cout.precision(6);
    typedef std::vector<functionT> vecfuncT;

    // take the HF orbitals to start
    const std::string input="input";
	Calculation calc(world,input.c_str());
    calc.molecule.print();
    print("\n");
    calc.param.print(world);

    // solve the ground state energy; calc is a reference
    MolecularEnergy me(world,calc);
    double hf_energy=me.value(calc.molecule.get_all_coords());

    if (world.rank()==0) print("MRA hf energy", hf_energy);
    if (world.rank()==0) {
    	printf("\n\n starting TDHF section at time %.1f\n",wall_time());
    	print("nuclear repulsion: ", calc.molecule.nuclear_repulsion_energy());
    	print("hf energy:         ", hf_energy);
    	print("orbital energies:  ");
    	for (std::size_t i=0; i<calc.amo.size(); ++i)
    		print("     ",calc.aeps[i]);
    }



    // construct the CIS solver, it requires a converged HF reference
    CIS cis(world,calc,input);

    // print grid information to file to get a better guess from external
    if (cis.print_grid()) {
    	if (world.rank()==0) print("printing grid for koala\n");
    	real_function_3d density=cis.get_calc().make_density(world,calc.aocc,
    		calc.amo);
    	density.get_impl()->print_grid("grid");
    } else {
    	
    	
    	// solve the response equation
    	cis.solve();
    }

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
