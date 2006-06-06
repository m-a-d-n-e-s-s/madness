/// \file mra/test.cc

#include <iostream>
using std::cout;
using std::endl;
#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>

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
    input_data testdata;

    try {
	FunctionDefaults::k=7;
	FunctionDefaults::initial_level=2;
	Function<double> f = FunctionFactory<double>(V).thresh(1e-3).compress(0);
	print("Tree in scaling function form");
	f.pnorms();
	f.compress();
	print("Tree in wavelet form");
	f.pnorms();
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

