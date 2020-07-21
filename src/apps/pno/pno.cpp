//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <iomanip>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/PNO.h>

using namespace madness;

// DEFINE PARAMETER TAGS FOR THE INPUT FILE
const std::string TAG_PNO = "pno";
const std::string TAG_F12 = "f12";
const std::string TAG_CP = "computeprotocol";


int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
	if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
	//const double time_start = wall_time();
	std::cout.precision(6);

	startup(world,argc,argv,true);
	print_meminfo(world.rank(), "startup");

	if(world.rank()==0){
		std::cout << "\n\n";
		std::cout << "-------------------------------------------------------------------------------------\n";
		std::cout << "SOLVING MRA-PNO-F12 as described in \n";
		std::cout << "J.S. Kottmann, F.A. Bischoff, E.F. Valeev\n";
		std::cout << "Direct determination of optimal pair-natural orbitals in a real-space representation:\n";
		std::cout << "the second-order Møller-Plesset energy\n";
		std::cout << "Journal of Chemical Physics ... 2019\n";
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

	// Compute MRA-PNO-MP2-F12
	const double time_pno_start = wall_time();
	PNOParameters parameters(world,input,nemo.get_calc()->molecule,TAG_PNO);
	F12Parameters paramf12(world, input, parameters, TAG_F12);
	PNO pno(world, nemo, parameters, paramf12);
	std::vector<PNOPairs> all_pairs;
	pno.solve(all_pairs);
	const double time_pno_end = wall_time();


	if(world.rank()==0){
		std::cout << std::setfill(' ');
		std::cout << "\n\n\n";
		std::cout << "--------------------------------------------------\n";
		std::cout << "MRA-PNO-MP2-F12 ended \n";
		std::cout << "--------------------------------------------------\n";
		std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
		std::cout << std::setw(25) << "energy scf" << " = " << scf_energy << "\n";
		std::cout << "--------------------------------------------------\n";
	}
	double mp2_energy = 0.0;
	std::cout<< std::setw(25) << "time pno" << " = " << time_pno_end - time_pno_start << "\n";
	for(const auto& pairs: all_pairs){
		if(pairs.type == MP2_PAIRTYPE){
			mp2_energy = pairs.energies.total_energy();
		}
		std::pair<size_t, size_t> ranks= pno.get_average_rank(pairs.pno_ij);
		if(world.rank()==0){
			std::string name;
			std::stringstream ss;
			ss << pairs.type;
			ss >> name;
			std::cout<< std::setw(25) << "energy "+name << " = " << pairs.energies.total_energy() << "\n";
			std::cout<< std::setw(25) << "average pno rank " + name << " = " << ranks.first << "\n";
			std::cout<< std::setw(25) << "max pno rank " + name << " = " << ranks.second << "\n";
		}
	}
	if(world.rank()==0 and mp2_energy != 0.0){
		std::cout << "--------------------------------------------------\n";
			std::cout<< std::setw(25) << "energy(total)" << " = " << scf_energy + mp2_energy << "\n";
			std::cout << "--------------------------------------------------\n";
			std::cout << "\n\n\n";
	}

	world.gop.fence();
	if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
	print_stats(world);
	finalize();
	return 0;
}
