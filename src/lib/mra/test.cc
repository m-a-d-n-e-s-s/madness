#include <iostream>
using std::cout;
using std::endl;

//#include <papi.h>

#include <cstring>
using std::strcmp;

#include <vector>
using std::vector;

/// \file mra/test.cc

#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
#include <misc/madexcept.h>
#include <serialize/vecar.h>

#include <unistd.h>

using namespace madness;

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

double dfred_dx(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*x;
}

double dfred_dy(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*y;
}

double dfred_dz(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*z;
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

//    int retval, Events[3] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_INS};
//    long_long values[3];
//    if (PAPI_start_counters(Events,3) != PAPI_OK) return 1;
    

    Communicator& comm = startup(argc,argv);

    // To ensure reliable cleanup catch all C++ exceptions here
    try {
        //comm.set_debug(true);
        // Do useful stuff here
        FunctionDefaults::k=9;
        FunctionDefaults::initial_level=2;
        Function<double> f = FunctionFactory<double>(fred).thresh(1e-7).nocompress().refine();
        Function<double> df,dfexact;

        print("valuesX",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        f.compress();
        print("DONE WITH THE COMPRESS");
        f.reconstruct();
        print("DONE WITH THE RECONSTRUCT");
        f.pnorms();
        print("valuesX",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        //f.pnorms();
        
        //f.square();
        //print("valuesSQ",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        
        
//        Function<double> p = f*f;
//        print("valuesSQ",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48),p(0.45,0.53,0.48));
        
        
//        Function<double> y = f.copy().autorefine();
//        print("err in autoref",y.norm2sq(),f.norm2sq(),(f-y).norm2sq());
//        y.compress();
//        y.reconstruct();
//        print("err in autoref",y.norm2sq(),f.norm2sq(),(f-y).norm2sq());
//
//        df = y.diff(0);
//        dfexact = FunctionFactory<double>(dfred_dx).thresh(1e-9);
//        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq(),y.norm2sq(),(f-y).norm2sq());
//        print("diff x",df(0.45,0.53,0.48),dfred_dx(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

        print("ABOUT TO DIFF");
        df = f.diff(0);
        print("DONE THE DIFF");
        dfexact = FunctionFactory<double>(dfred_dx).thresh(1e-7).nocompress().norefine();
        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff x",df(0.45,0.53,0.48),dfred_dx(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

        goto done;

        df = f.diff(1);
        dfexact = FunctionFactory<double>(dfred_dy).thresh(1e-7).nocompress();
        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff y",df(0.45,0.53,0.48),dfred_dy(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

        df = f.diff(2);
        dfexact = FunctionFactory<double>(dfred_dz).thresh(1e-7).nocompress();
        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff z",df(0.45,0.53,0.48),dfred_dz(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

       
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


//    if (PAPI_stop_counters(values,3) != PAPI_OK) comm.Abort();
//    print("PAPI",values[0],values[1],values[2]);
//    comm.global_sum(values,3);
//    print("PAPI SUM",values[0],values[1],values[2]);
    


    // The following should be used for succesful termination
    done:
    comm.close();
    MPI::Finalize();
    return 0;
}

