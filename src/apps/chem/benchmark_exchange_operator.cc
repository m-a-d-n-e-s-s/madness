
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness.h>
#include <chem/SCFOperators.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/test_utilities.h>

using namespace madness;

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    srand (time(NULL));

    std::string structure="water2";
    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="--structure") structure=val;
    }

    std::string input1=R"input(
dft
    			econv 1.e-6
				xc hf
				protocol []
				multiworld true
)input";
    std::string input2=R"input(
end
)input";

	test_inputfile ifile("input",input1+"molecular_structure "+structure+"\n"+input2);
	ifile.keepfile=true;
	world.gop.fence();

	double cpu0=cpu_time();

    SCF calc(world,"input");
    calc.param.print("","");
    calc.molecule.print();

    calc.set_protocol<3>(world,1.e-4);
	MolecularEnergy me(world, calc);
	me.value(calc.molecule.get_all_coords());
	Exchange<double,3> K=Exchange<double,3>(world,&calc,0).multiworld(true);

	if (world.size() > 1) {
		LoadBalanceDeux < 3 > lb(world);
		for (unsigned int i = 0; i < calc.amo.size(); ++i) {
			lb.add_tree(calc.amo[i], lbcost<double, 3>(1.0, 8.0), false);
		}
		world.gop.fence();
		FunctionDefaults < 3 > ::redistribute(world, lb.load_balance(calc.param.loadbalparts())); // 6.0 needs retuning after param.vnucextra
		world.gop.fence();
	}

	vecfuncT tmp, tmp1, reference;
	double err;

	double cpu1=cpu_time();
	if (world.rank()==0) printf("\ntimings for preparation   %8.2fs\n",cpu1-cpu0);

    cpu0=cpu1;
    K.ntask_per_subworld=100;
    K.multiworld_=false;
    K.small_memory(false);
    K.same(true);
    reference=K(calc.amo);
    cpu1=cpu_time();
    if (world.rank()==0) printf("\ntimings exchange operator no multiworld largemem   %8.2fs\n",cpu1-cpu0);
//    err=norm2(world,reference-tmp);
//    print("error wrt small-memory result",err);

    cpu0=cpu1;
    K.ntask_per_subworld=100;
    K.multiworld_=true;
    K.efficient_=true;
    tmp=K(calc.amo);
    cpu1=cpu_time();
    if (world.rank()==0) printf("\ntimings exchange operator efficient (ntask_per_subworld=100) %8.2fs\n",cpu1-cpu0);
    err=norm2(world,reference-tmp);
    print("error wrt small-memory result",err);

	cpu0=cpu1;
	K.ntask_per_subworld=100;
	K.multiworld_=false;
    K.efficient_=false;
	K.small_memory(true);
	reference=K(calc.amo);
	cpu1=cpu_time();
	if (world.rank()==0) printf("\ntimings exchange operator no multiworld smallmem   %8.2fs\n",cpu1-cpu0);
    err=norm2(world,reference-tmp);
    print("error wrt small-memory result",err);

	cpu0=cpu1;
	K.ntask_per_subworld=100;
	K.multiworld_=true;
	tmp=K(calc.amo);
	cpu1=cpu_time();
	if (world.rank()==0) printf("\ntimings exchange operator (ntask_per_subworld=100) %8.2fs\n",cpu1-cpu0);
	err=norm2(world,reference-tmp);
	print("error wrt small-memory result",err);

	cpu0=cpu1;
	K.ntask_per_subworld=10;
	tmp=K(calc.amo);
	cpu1=cpu_time();
	if (world.rank()==0) printf("\ntimings exchange operator (ntask_per_subworld=10)  %8.2fs\n",cpu1-cpu0);
	err=norm2(world,reference-tmp);
	print("error wrt small-memory result",err);

	cpu0=cpu1;
	K.ntask_per_subworld=3;
	tmp=K(calc.amo);
	cpu1=cpu_time();
	if (world.rank()==0) printf("\ntimings exchange operator (ntask_per_subworld=3)   %8.2fs\n",cpu1-cpu0);
	err=norm2(world,reference-tmp);
	print("error wrt small-memory result",err);

    madness::finalize();
    return 0;
}
