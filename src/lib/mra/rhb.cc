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
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);
    std::cout << "about to make ft" << std::endl;
    SharedPtr<FunctionOctTree> ft = SharedPtr<FunctionOctTree> (new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2)));
    std::cout << "made ft" << std::endl;
    FunctionDefaults::tree = std::vector<SharedPtr<FunctionOctTree> >();
    FunctionDefaults::tree.push_back(SharedPtr<FunctionOctTree>(ft));

    std::cout << "pushed back ft" << std::endl;
    std::cout << "made tree" << std::endl;

if (FunctionDefaults::tree[0]->tree())
{
    std::cout << "Depth First Traversal:" << std::endl;
    FunctionDefaults::tree[0]->tree()->depthFirstTraverse();
    std::cout << "End Depth First Traversal" << std::endl;
}
else
{
    std::cout << "no depth first traversal; tree does not exist" << std::endl;
}

MPI::COMM_WORLD.Barrier();
    if (!gauss_legendre_test()) comm.Abort();
std::cout << "tested gauss-legendre" << std::endl;
MPI::COMM_WORLD.Barrier();
    if (!test_two_scale_coefficients()) comm.Abort();
std::cout << "tested two scale coeffs" << std::endl;
MPI::COMM_WORLD.Barrier();

    FunctionDefaults::k=7;
    FunctionDefaults::initial_level=2;
std::cout << "set function defaults k and initial_level" << std::endl;
/*
if (FunctionDefaults::tree[0]->tree())
{
    std::cout << "Depth First Traversal:" << std::endl;
    FunctionDefaults::tree[0]->tree()->depthFirstTraverse();
    std::cout << "End Depth First Traversal" << std::endl;
}
else
{
    std::cout << "no depth first traversal; tree does not exist" << std::endl;
}
*/
MPI::COMM_WORLD.Barrier();

    try {
std::cout << "about to create function" << std::endl;
/*
if (FunctionDefaults::tree[0]->tree())
{
    std::cout << "Depth First Traversal:" << std::endl;
    FunctionDefaults::tree[0]->tree()->depthFirstTraverse();
    std::cout << "End Depth First Traversal" << std::endl;
}
else
{
    std::cout << "no depth first traversal; tree does not exist" << std::endl;
}
*/
MPI::COMM_WORLD.Barrier();
std::cout << "about to init function" << std::endl;
//	Function<double> f = FunctionFactory<double>(fred).thresh(1e-5).compress(0);
	Function<double> f = FunctionFactory<double>(fred).thresh(1e-2).compress(0);
std::cout << "waiting for barrier" << std::endl;
MPI::COMM_WORLD.Barrier();
std::cout << "created function" << std::endl;
	print("Tree in scaling function form");
	f.pnorms();
std::cout << "about to compress" << std::endl;
	f.compress();
std::cout << "waiting for barrier" << std::endl;
MPI::COMM_WORLD.Barrier();
	print("Tree in wavelet form");
	f.pnorms();
std::cout << "about to reconstruct" << std::endl;
	f.reconstruct();
	//load bal
	int tlen = FunctionDefaults::tree.size();
	std::cout << "before load balancing: " << std::endl;
	for (int i = 0; i < tlen; i++)
	{
	    std::cout << "Subtree:" << std::endl;
    	    FunctionDefaults::tree[0]->tree()->depthFirstTraverse();
	}
std::cout << "about to load balance" << std::endl;
	balanceFunctionOctTree(&FunctionDefaults::tree);
std::cout << "after load balance" << std::endl;
//	int tlen = FunctionDefaults::tree.size();
	tlen = FunctionDefaults::tree.size();
	for (int i = 0; i < tlen; i++)
	{
	    std::cout << "Subtree:" << std::endl;
    	    FunctionDefaults::tree[0]->tree()->depthFirstTraverse();
	}
	std::cout << "GAME OVER!!!" << std::endl;
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

