/*
 * cis.cpp
 *
 *  Created on: Aug 10, 2015
 *      Author: Jakob S. Kottmann
 */

#include <madness/world/worldmem.h>
#include <chem/TDHF.h>
#include <chem/commandlineparser.h>
#include <madness/misc/gitinfo.h>


using namespace madness;


int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0){
		printf("starting at time %.1f\n", wall_time());
		print("\nmain() compiled at ",__TIME__," on ",__DATE__);
        print(info::print_revision_information());
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

	commandlineparser parser(argc,argv);
	parser.print_map();

	int error=0;
	if (parser.key_exists("test")) {
	    print("entering test mode");
	    error=TDHF::test(world,parser);
	} else {

	    TDHF tdhf(world,parser);

        tdhf.get_calcparam().print("dft");
        tdhf.parameters.print("response");

        // solve the CIS equations
        const double time_scf_start = wall_time();
        tdhf.prepare_calculation();
        const double time_scf_end = wall_time();
        if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

        const double time_cis_start = wall_time();
        std::vector<CC_vecfunction> roots = tdhf.solve_cis();
        const double time_cis_end = wall_time();
        if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

        if (world.rank() == 0) {
            std::cout << std::setfill(' ');
            std::cout << "\n\n\n";
            std::cout << "--------------------------------------------------\n";
            std::cout << "MRA-CIS ended \n";
            std::cout << "--------------------------------------------------\n";
            std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
            std::cout << std::setw(25) << "time cis" << " = " << time_cis_end - time_cis_start << "\n";
            std::cout << "--------------------------------------------------\n";
        }
        tdhf.analyze(roots);
    }
	world.gop.fence();
	if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
	print_stats(world);
	finalize();
	return error;
}



