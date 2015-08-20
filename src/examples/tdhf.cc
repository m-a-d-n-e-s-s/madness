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

// for moldft consistency check
//#include<examples/dft_solver.h>
//
#include <chem/TDA.h>

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

#ifdef MADNESS_GITREVISION
	const  char* gitrev =  MADNESS_GITREVISION;
	const std::string gitrevision(gitrev);
	if (world.rank()==0) {
		print("           git revision ...",gitrevision);
	}
#endif

	typedef std::vector<functionT> vecfuncT;

	// Set the threshold defaults
	double ground_state_thresh = FunctionDefaults<3>::get_thresh();
	double guess_thresh = ground_state_thresh*50.0;
	double solve_thresh = ground_state_thresh*10.0;
	double solve_seq_thresh = ground_state_thresh*10.0;
	bool print_grid=false;
	bool no_compute=false;
	bool only_sequential=false;
	bool use_nemo=true;

	for(size_t i=0;i<argc;i++){
		if(strcmp(argv[i],"-tda_print_grid")==0) print_grid = true;
		if(strcmp(argv[i],"-tda_no_compute")==0) no_compute = true;
		if(strcmp(argv[i],"-tda_analyze")==0) no_compute = true;
		if(strcmp(argv[i],"-tda_sequential")==0) only_sequential = true;
		if(strcmp(argv[i],"-nonemo")==0) use_nemo = false;
	}

	// First solve the ground state
	const std::string input = "input";
	//SCF calc(world,input.c_str());
	std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
	Nemo nemo(world,calc);
    if (world.rank()==0) {
        calc->molecule.print();
        print("\n");
        calc->param.print(world);
    }

    //
    double hf_energy =0;
    if(not use_nemo){
    	std::cout << "\n\n\n\n !!!!!!!!!! -- nonemo input parameter found ... using moldft ...  -- !!!!!!!!!!!!\n\n\n";
    	if(nemo.get_calc()->param.nuclear_corrfac != "none"){
    		std::cout << "WARNING NUCLEAR CORRELATION FACTOR IS SET TO: " << nemo.get_calc() -> param.nuclear_corrfac << "\n\n\n"<< std::endl;
    	}
    	// constructing the correlation factor as constant function
    	nemo.construct_nuclear_correlation_factor();
    	MolecularEnergy E(world, *nemo.get_calc());
    	hf_energy = E.value(calc -> molecule.get_all_coords().flat());
    }else{
    	std::cout << "\nSTARTING GROUND STATE CALCULATION\n";
    	hf_energy=nemo.value();
    }

    if (world.rank()==0) {
        printf("final energy   %12.8f\n", hf_energy);
        printf("finished at time %.1f\n", wall_time());
    }


	// Get the custom thresholds from the input file
	std::ifstream f(input.c_str());
	position_stream(f, "excite");
	std::string s, tag;
	while (std::getline(f,s)) {
		std::stringstream ss(s);
		ss >> tag;
		if (tag == "end") break;
		else if(tag == "guess_thresh") ss >> guess_thresh;
		else if(tag == "solve_thresh") ss >> solve_thresh;
		else if(tag == "solve_seq_thresh") ss >> solve_seq_thresh;
		else if(tag == "print_grid") print_grid=true;
		else if(tag == "no_compute") no_compute = true;
		else continue;
	}

	if (world.rank()==0) print("MRA hf energy", hf_energy);
	if (world.rank()==0) {
		printf("\n\n starting TDHF section at time %.1f\n",wall_time());
		print("nuclear repulsion: ", calc -> molecule.nuclear_repulsion_energy());
		print("hf energy:         ", hf_energy);
		print("orbital energies:  ");
		for (std::size_t i=0; i<calc -> amo.size(); ++i)
			print("     ",calc -> aeps[i]);
		print("\n\n-----------------------------");
		std::cout << "guess_thresh is:" << guess_thresh << std::endl;
		std::cout << "solve_thresh is:" << solve_thresh << std::endl;
		std::cout << "solve_seq_thresh is:" << solve_seq_thresh << std::endl;
		print("-----------------------------\n\n");

	}

	// Print the coordinates
	if (world.rank()==0){
		print("Coordinates after MolDFT:");
		Tensor<double> Coord = calc -> molecule.get_all_coords();
		for(int i=0;i<calc -> molecule.natom();i++){
			print(calc -> molecule.get_atom(i).atomic_number,Coord(i,0),Coord(i,1),Coord(i,2));
		}
	}
	// Check if center of charge is 0,0,0
	Tensor<double> Coord = calc -> molecule.get_all_coords();
	coord_3d cm(0.0);
	for(size_t i=0;i<3;i++){
		for(int j=0;j<calc -> molecule.natom();j++){
			cm[i] = calc -> molecule.get_atom(j).atomic_number*Coord(i,j);
		}
	}
	std::cout << "Center of Charge is " << cm << std::endl;

	// Make the empty vectors for the molecular orbitals and xfunctions with different thresholds
	vecfuncT guess_mos;
	vecfuncT solve_mos;
	vecfuncT solve_seq_mos;
	xfunctionsT guess_xfunctions;
	xfunctionsT solve_xfunctions;
	xfunctionsT solve_seq_xfunctions;

	// Print the grid (for koala guess or other external applications)
	if(print_grid){
		if(world.rank()==0) std::cout << " Printing grid for external application ... then stop the calculation ... remeber to remove print_grid keyword" << std::endl;
		real_function_3d density=calc -> make_density(world,calc -> aocc,calc -> amo);
    	density.get_impl()->print_grid("grid");
    	if (world.rank() == 0) printf("\n\n-------------------------\nfinished at time %.1f\n", wall_time());
    	finalize();
    	return 0;
	}

	// just read and analyze the xfunctions
	if(no_compute){
		MADNESS_EXCEPTION("Can not read xfunction currently due to strange crash in read_xfunctions",1);
//		if(world.rank()==0) std::cout << "\n\n\n---no_compute keyword detected ... just read and analyze---\n\n\n " << std::endl;
//		TDA cis(world,calc,calc -> amo,input);
//		xfunctionsT read_roots;
//		cis.read_xfunctions(read_roots);
//		xfunctionsT read_converged_roots = cis.get_converged_xfunctions();
//
//		if(not read_roots.empty()){
//			if(world.rank()==0) std::cout << "\nAnalyze the given " << read_roots.size() << " excitation vectors ..." << std::endl;
//			cis.analyze(read_roots);
//		}
//		if(not read_converged_roots.empty()){
//			if(world.rank()==0) std::cout << "\nAnalyze the given " << read_converged_roots.size() << " excitation vectors ..." << std::endl;
//			cis.analyze(read_converged_roots);
//		}else if(world.rank()==0) std::cout << "\nno converged excitation vectors found ..." << std::endl;
//
//    	finalize();
		return 0;
	}

	if(only_sequential){
		MADNESS_EXCEPTION("Can not read xfunction currently due to strange crash in read_xfunctions",1);
		return 0;
	}

	// Initialize and pre converge the guess xfunctions
	{
		// Set the guess_thresh
		FunctionDefaults<3>::set_thresh(guess_thresh);
		if(world.rank()==0) std::cout << "\n\nMODULE 1: INITIALIZE AND PRE-CONVERGE GUESS EXCITATIONS \n used threshold"
				<< FunctionDefaults<3>::get_thresh() << "\n\n"  << std::endl;

		// Make deep copys of the ground state mos and reduce the thresh
		guess_mos = copy(world,calc -> amo);
		set_thresh(world,guess_mos,guess_thresh);

		// Init. and solve the guess
		TDA guess_cis(world,nemo,guess_mos,input);
		guess_cis.solve_guess(guess_xfunctions);

		// Get the pre-converged xfunctions for the next solve routine
		solve_xfunctions = guess_cis.get_converged_xfunctions();
		guess_mos.clear();
		guess_xfunctions.clear();
	}// Destroy the guess_cis object

	// Solve the pre converged guess_xfunctions all together
	{
		// set the solve_thresh
		FunctionDefaults<3>::set_thresh(solve_thresh);
		if(world.rank()==0) std::cout << "\n\nMODULE 2: ITERATE GUESS EXCITATIONS \n used threshold: "
				<< FunctionDefaults<3>::get_thresh() << "\n\n"  << std::endl;

		// Make deep copys of the ground state mos and reduce the thresh
		solve_mos = copy(world,calc -> amo);
		set_thresh(world,solve_mos,solve_thresh);

		// iterate the solve_xfunctions
		TDA solve_cis(world,nemo,solve_mos,input);
		solve_cis.solve(solve_xfunctions);

		// Get the converged xfunctions and prepare for the last sequential iterations with kain
		solve_seq_xfunctions = solve_cis.get_converged_xfunctions();
		solve_mos.clear();
		solve_xfunctions.clear();

	}// destroy solve_cis object

	// Final iterations, sequantial and with kain (keyword: kain_subspace 1, is without kain)
	{
		// set the sequential thresh
		FunctionDefaults<3>::set_thresh(solve_seq_thresh);
		if(world.rank()==0) std::cout << "\n\nMODULE 3: SEQUENTIALLY ITERATE GUESS EXCITATIONS \n used threshold: "
				<< FunctionDefaults<3>::get_thresh() << "\n\n"  << std::endl;

		// Make deep copys of the ground state mos and reduce the thresh
		solve_seq_mos = copy(world,calc -> amo);
		set_thresh(world,solve_seq_mos,solve_seq_thresh);

		// make last iterations
		TDA solve_seq_cis(world,nemo,solve_seq_mos,input);
		solve_seq_cis.solve_sequential(solve_seq_xfunctions);
		solve_seq_mos.clear();
	}

	if (world.rank() == 0) printf("\n\n-------------------------\nfinished at time %.1f\n", wall_time());
	finalize();
	return 0;
}// end main

