#include <iostream>
using std::cout;
using std::endl;

/// \file mra/test.cc

#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>

using namespace madness;

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5; y-=0.5; z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}
//double fred(double x, double y, double z) {
//    return x*x+y*y*z*z;
//}

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
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);
    FunctionDefaults::tree = new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2));
    if (!gauss_legendre_test()) comm.Abort();
    if (!test_two_scale_coefficients()) comm.Abort();
    
    // To ensure reliable cleanup catch all C++ exceptions here
        // Do useful stuff below here
        FunctionDefaults::k=9;
        FunctionDefaults::initial_level=0;
        Function<double> f = FunctionFactory<double>(fred).refine(1).compress(0).initial_level(4).thresh(1e-9);
        Function<double> g = FunctionFactory<double>(fred).refine(1).compress(0).initial_level(4).thresh(1e-9);
        //print("normsq after projection    ",f.norm2sq_local());
        print("normsq after projection    ",f.norm2sq());
        
        for (int i=0; i<=16; i++) {
            double z = i/16.0;
            double value = 0;
            f.eval_local(z,z,z,&value); 
            print("test eval",z,value,fred(z,z,z));
        }
        f.compress();
        print("norm2sq after compression   ",f.norm2sq());
        f.reconstruct();
        print("norm2sq after reconstruction   ",f.norm2sq());
        f.truncate();
        print("norm2sq after truncation   ",f.norm2sq());
	double inner_test = f.inner(g);
        cout << " inner_test = " << inner_test << endl;
        print("norm2sq after inner   ",f.norm2sq());
/*
	print("Tree in scaling function form");
	f.pnorms();
	print("Tree in wavelet form");
	f.compress();
	f.pnorms();
	f.norm2sq_local();
	

        f.reconstruct();
        
        Function<double> g;
        g = copy(f);
        print("start of statement");
        g = 2.0 + g + 3.0*g - f*2.0 - 1.0;
        print("end of statement");
        
        print(f(0.5,0.5,0.5));
        f.square();
        print(f(0.5,0.5,0.5));
        
        //Function<double_complex> cf = FunctionFactory<double_complex>(cfred).refine(1).initial_level(1);
        //print("normsq",cf.norm2sq_local());
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
    
*/
    // The follwing should be used for succesful termination
    comm.close(); 
    MPI::Finalize();
    return 0;
}

