/// \file mra/test.cc

#include <iostream>
using std::cout;
using std::endl;
#include <mra/mraX.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
//#include <octtree/sendrecv.h>

using namespace madness;

double fred(double x, double y, double z) {
    const double PI = 3.14;
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5; y-=0.5; z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

int main(int argc, char* argv[]) {
    // The following should be used to setup all calculations
    // 1) Initialize parallel environment
    // 2) Setup communication information
    // 3) Redirect standard output+err for parallel processes
    // 4) Load coeffs and quadrature information from file
    // 5) Setup default OctTreeLayout
    // 6) Sanity check
    MPI::Init(argc, argv);
    Communicator comm;
    madness::comm_default = &comm;
    redirectio(comm);
    comm.print();
//    comm.set_debug(true);
    load_coeffs(comm);
    load_quadrature(comm);

    OctTree<FunctionNode> *t = new OctTree<FunctionNode>();

    std::vector<SharedPtr<OctTree<FunctionNode> > > treeList;
    std::cout << "created tree list" << std::endl;
    treeList.push_back(SharedPtr<OctTree<FunctionNode> >(t->create_default(comm,1)));
    std::cout << "pushed back tree list" << std::endl;
    FunctionDefaults::tree = SharedPtr<FunctionOctTree> (new FunctionOctTree(treeList));
    std::cout << "made tree" << std::endl;

if (FunctionDefaults::tree->tree(0))
{
    std::cout << "Depth First Traversal:" << std::endl;
    FunctionDefaults::tree->tree(0)->depthFirstTraverse();
    std::cout << "End Depth First Traversal" << std::endl;
}
else
{
    std::cout << "no depth first traversal; tree does not exist" << std::endl;
}

    if (!gauss_legendre_test()) comm.Abort();
    std::cout << "tested gauss-legendre" << std::endl;
    if (!test_two_scale_coefficients()) comm.Abort();
    std::cout << "tested two scale coeffs" << std::endl;

    FunctionDefaults::k=7;
    FunctionDefaults::initial_level=2;
    std::cout << "set function defaults k and initial_level" << std::endl;

    MPI::COMM_WORLD.Barrier();

//    xterm_debug(comm,0,0);

    try {
	std::cout << "about to create function" << std::endl;

	MPI::COMM_WORLD.Barrier();
	double t1 = MPI::Wtime();
	std::cout << "about to init function" << std::endl;
	Function<double> f = FunctionFactory<double>(fred).thresh(1e-5).compress(0);
//	Function<double> f = FunctionFactory<double>(fred).norefine().thresh(1e-2).compress(0);
//	Function<double> f = FunctionFactory<double>(fred).thresh(1e-2).compress(0);
	MPI::COMM_WORLD.Barrier();
	double t2 = MPI::Wtime();
	std::cout << "created function" << std::endl;

	print("Tree in scaling function form");
	f.pnorms();

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t3 = MPI::Wtime();
	f.compress();
	MPI::COMM_WORLD.Barrier();
	double t4 = MPI::Wtime();

	print("Tree in wavelet form");
	f.pnorms();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t5 = MPI::Wtime();
	f.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t6 = MPI::Wtime();
	//load bal
	std::cout << "before load balancing: " << std::endl;
    	FunctionDefaults::tree->depthFirstTraverse();

	std::cout << "about to load balance" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t7 = MPI::Wtime();
	balanceFunctionOctTree(FunctionDefaults::tree);
//	balanceFunctionOctTree(f.data->trees);
	MPI::COMM_WORLD.Barrier();
	double t8 = MPI::Wtime();

	f.data->trees->depthFirstTraverse();

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t9 = MPI::Wtime();
	f.compress();
	MPI::COMM_WORLD.Barrier();
	double t10 = MPI::Wtime();

	print("Tree in wavelet form");
	f.pnorms();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t11 = MPI::Wtime();
	f.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t12 = MPI::Wtime();

	std::cout << "GAME OVER!!!" << std::endl;

	std::cout << std::endl << "Timings: " << std::endl;
	std::cout << "           Initialize    Load Balance    Compress     Reconstruct" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "Default     " << t2-t1 << "     --------        " << t4-t3 << "       " << t6-t5 
		<< std::endl;
	std::cout << "Balanced    --------     " << t8-t7 << "        " << t10-t9 << "       " << t12-t11 
		<< std::endl;

    }
    catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.Abort();
    }
    catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
        comm.Abort();
    }
    catch (TensorException& e) {
        std::cerr << e << std::endl;
        comm.Abort();
    }
    catch (MPI::Exception& e) {
        std::cerr << "Exception (mpi): code=" << e.Get_error_code() 
                  << ", class=" << e.Get_error_class() 
                  << ", string=" << e.Get_error_string() << std::endl;
        comm.Abort();
    }
    catch (...) {
        std::cerr << "Exception (general)" << std::endl;
        comm.Abort();
    }
    
    comm.close(); 
    MPI::Finalize();
    return 0;
}

