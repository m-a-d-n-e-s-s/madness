/*
 * cis.cpp
 *
 *  Created on: Aug 10, 2015
 *      Author: Jakob S. Kottmann
 */

#include <madness/world/worldmem.h>
#include <chem/TDHF.h>

using namespace madness;


int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0){
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);
#ifdef GITREVISION
			const  char* gitrev =  GITREVISION;
			const std::string gitrevision(gitrev);
			if (world.rank()==0) {
				print("           git revision ...",gitrevision);
			}
#endif
	}
	//const double time_start = wall_time();
	std::cout.precision(6);

	startup(world,argc,argv,true);
	print_meminfo(world.rank(), "startup");

	if(world.rank()==0){
		std::cout << "\n\n";
		std::cout << "-------------------------------------------------------------------------------------\n";
		std::cout << "SOLVING MRA-CIS as described in \n";
		std::cout << "J.S. Kottmann, S. Höfener ,F.A. Bischoff\n";
		std::cout << "Numerically accurate linear response-properties in the configuration-interaction singles (CIS) approximation \n";
		std::cout << "Phys. Chem. Chem. Phys., 2015, 17, 31453-31462\n";
		std::cout << "DOI: 10.1039/C5CP00345H\n";
		std::cout << "-------------------------------------------------------------------------------------\n";
		std::cout << "\n\n";
	}

	// Get the name of the input file (if given)
	const std::string input = (argc > 1) ? argv[1] : "input";

	// Compute the SCF Reference
	const double time_scf_start = wall_time();
	std::shared_ptr<SCF> calc(new SCF(world, input));
	Nemo nemo(world, calc, input);
	nemo.get_calc()->param.print();
	const double scf_energy = nemo.value();
	if (world.rank() == 0) print("nemo energy: ", scf_energy);
	if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
	const double time_scf_end = wall_time();

	// Compute MRA-CIS
	const double time_cis_start = wall_time();
	TDHF tdhf(world,nemo,input);
	// solve the CIS equations
	std::vector<CC_vecfunction> roots=tdhf.solve_cis();

	const double time_cis_end = wall_time();
	if(world.rank()==0){
		std::cout << std::setfill(' ');
		std::cout << "\n\n\n";
		std::cout << "--------------------------------------------------\n";
		std::cout << "MRA-CIS ended \n";
		std::cout << "--------------------------------------------------\n";
		std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
		std::cout << std::setw(25) << "energy scf" << " = " << scf_energy << "\n";
		std::cout << std::setw(25) << "time cis" << " = " << time_cis_end - time_cis_start << "\n";
		std::cout << "--------------------------------------------------\n";
	}
	tdhf.analyze(roots);

	world.gop.fence();
	if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
	print_stats(world);
	finalize();
	return 0;
}



