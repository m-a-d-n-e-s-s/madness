#include <iostream>
using std::cout;
using std::endl;

#include <cstring>
using std::strcmp;

/// \file mra/test.cc

#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
#include <misc/madexcept.h>

using namespace madness;

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

double mary(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.4;
    y-=0.6;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

double_complex cfred(double x, double y, double z) {
    return x*x+y*y*z*z;
}

int main(int argc, char* argv[]) {
    // The following should be used to setup all calculations
    // 1) Initialize parallel environment
    // 2) Setup communication information
    // 3) Redirect standard output+err for parallel processes
    // 4) Load coeffs and quadrature information from file
    // 5) Setup default OctTreeLayout
    // 6) Sanity check
    // 7) Top level catching of exceptions
    MPI::Init(argc, argv);
    Communicator comm;
    madness::comm_default = &comm;
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);
    FunctionDefaults::tree = new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2));
    if (!gauss_legendre_test()) comm.Abort();
    if (!test_two_scale_coefficients()) comm.Abort();

    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i],"-d") == 0) xterm_debug(comm,0,0);
    }

    // To ensure reliable cleanup catch all C++ exceptions here
    try {
        // Do useful stuff here
        FunctionDefaults::k=6;
        FunctionDefaults::initial_level=1;
        Function<double> f = FunctionFactory<double>(fred).thresh(1e-7).nocompress();
        print("Tree in scaling function basis");
        //f.pnorms();
        f.compress();
        print("Tree in wavelet basis");
        //f.pnorms();
        f.reconstruct();
        print("Tree in scaling function basis");
        //f.pnorms();
        print("values",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        print("truncating fred");
        f.truncate();
        print("values",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        print("truncating fred again");
        f.truncate(1e-3);
        print("values",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        print("truncating fred yet again");
        f.truncate(1e-3);
        print("values",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));              
        //f.autorefine();
        //print("Tree in scaling function basis after autorefine");
        //f.pnorms();
        print("values",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        f.square();
        print("values",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        Function<double> m = FunctionFactory<double>(mary).thresh(1e-3).nocompress();
		Function<double>fm = f*m;
		print("f",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
		print("m",mary(0.45,0.53,0.48),m(0.45,0.53,0.48));
		print("fm",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48)*mary(0.45,0.53,0.48),fm(0.45,0.53,0.48));       
    } catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.Abort();
    } catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
        comm.Abort();
    } catch (MadnessException& e) {
    	std::cerr << e << std::endl;
        comm.Abort();
    } catch (TensorException& e) {
        std::cerr << e << std::endl;
        comm.Abort();
    } catch (MPI::Exception& e) {
        std::cerr << "Exception (mpi): code=" << e.Get_error_code()
        << ", class=" << e.Get_error_class()
        << ", string=" << e.Get_error_string() << std::endl;
        comm.Abort();
    } catch (...) {
        std::cerr << "Exception (general)" << std::endl;
        comm.Abort();
    }

    // The follwing should be used for succesful termination
    comm.close();
    MPI::Finalize();
    return 0;
}

