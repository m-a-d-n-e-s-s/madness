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

// for moldft consistency check
//#include<examples/dft_solver.h>
//
#include <chem/TDA.h>
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
    FunctionDefaults<3>::set_truncate_mode(1);
    print("Truncate mode set to ",FunctionDefaults<3>::get_truncate_mode());

#ifdef GITREVISION
    const  char* gitrev =  GITREVISION;
    const std::string gitrevision(gitrev);
    if (world.rank()==0) {
    	print("           git revision ...",gitrevision);
    }
#endif

    typedef std::vector<functionT> vecfuncT;

    // take the HF orbitals to start
    const std::string input="input";
    const std::string high_input="high_input";
	SCF calc(world,input.c_str());
    calc.molecule.print();
    print("\n");
    calc.param.print(world);

    // solve the ground state energy; calc is a reference
    MolecularEnergy me(world,calc);
    double hf_energy=me.value(calc.molecule.get_all_coords());

    // Check convergence and functional use of MolDFT
//    dft_solver dft(world,calc);
//    dft.solve();

    if (world.rank()==0) print("MRA hf energy", hf_energy);
    if (world.rank()==0) {
    	printf("\n\n starting TDHF section at time %.1f\n",wall_time());
    	print("nuclear repulsion: ", calc.molecule.nuclear_repulsion_energy());
    	print("hf energy:         ", hf_energy);
    	print("orbital energies:  ");
    	for (std::size_t i=0; i<calc.amo.size(); ++i)
    		print("     ",calc.aeps[i]);
    }




	// Print the coordinates
	if (world.rank()==0){
		print("Coordinates after MolDFT:");
		Tensor<double> Coord = calc.molecule.get_all_coords();
		for(int i=0;i<calc.molecule.natom();i++){
			print(calc.molecule.get_atom(i).atomic_number,Coord(i,0),Coord(i,1),Coord(i,2));
		}
	}


	// default threshs
	const double high_thresh = FunctionDefaults<3>::get_thresh();
	double low_thresh = sqrt(high_thresh);// sqrt(high_thresh);

	// fetch low thresh from the input section
	// high thresh is always the same as in moldft
	std::ifstream f(input.c_str());
	position_stream(f, "excite");
	std::string s, tag;
	while (std::getline(f,s)) {
		std::stringstream ss(s);
		ss >> tag;
		if (tag == "end") break;
		else if(tag=="low_thresh") ss>>low_thresh;
		else if(tag=="lo") ss>>low_thresh;
		else if(tag=="lo_thresh") ss>>low_thresh;
		else continue;
	}

	// create the vector of xfunctions (now empty, will be filled by cis class if it is empty)
	xfunctionsT xfunctions;

	// solving with low Threshold

	// make a deep copy of the occupied orbitals and change the thresh

	print("\n\n-----------------------------");
	std::cout << "high thresh is:" << high_thresh << std::endl;
	std::cout << "low thresh is:" << low_thresh << std::endl;
	print("-----------------------------\n\n");
	std::vector<real_function_3d> mos;

	calc.amo[0].print_size("hi, ");
	FunctionDefaults<3>::set_thresh(low_thresh);
	for(size_t i=0;i<calc.amo.size();i++){
		real_function_3d tmp = copy(calc.amo[i]);
		tmp.set_thresh(low_thresh);
		tmp.truncate();
		mos.push_back(tmp);
	}


    // construct the CIS solver, it requires a converged HF reference
    TDA cis(world,calc,mos,input,true);
bool asd =false;
    // print grid information to file to get a better guess from external
    if (asd){//cis.print_grid_TDA()) {
    	if (world.rank()==0) print("printing grid for koala\n");
    	real_function_3d density=cis.get_calc().make_density(world,calc.aocc,
    		calc.amo);
    	density.get_impl()->print_grid("grid");
    } else {

    	// solve the response equation
    	cis.solve(xfunctions);
    }

    // Now solve with high thresh
    vecfuncT high_mos;
	FunctionDefaults<3>::set_thresh(high_thresh);
	for(size_t i=0;i<calc.amo.size();i++){
		real_function_3d tmp = copy(calc.amo[i]);
		tmp.set_thresh(high_thresh);
		high_mos.push_back(tmp);
	}

	// set the thresh of the xfunctions up
	xfunctions.clear();
	xfunctionsT high_xfunctions;
	for(size_t i=0; i<cis.get_converged_xfunctions().size();i++){
		xfunction tmp(world);
		tmp.x = copy(world,cis.get_converged_xfunctions()[i].x);
		tmp.omega = cis.get_converged_xfunctions()[i].omega;
		tmp.converged = false;
		tmp.iterations = cis.get_converged_xfunctions()[i].iterations;
		tmp.delta = cis.get_converged_xfunctions()[i].delta;
		tmp.expectation_value = cis.get_converged_xfunctions()[i].expectation_value;
		tmp.error = cis.get_converged_xfunctions()[i].error;
		tmp.number = i;
		high_xfunctions.push_back(tmp);
	}
	for(size_t i=0;i<xfunctions.size();i++)set_thresh(world,xfunctions[i].x,high_thresh);

	// sort after energy
	std::sort(high_xfunctions.begin(),high_xfunctions.end());

	TDA high_cis(world,calc,high_mos,input,false);
	high_cis.solve_sequential(high_xfunctions);



    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
