#include <iostream>
using std::cout;
using std::endl;

//#define MADNESS_PAPI

#ifdef MADNESS_PAPI
#include <papi.h>
#endif

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

void vfred(long n, const double *vx, const double *vy, const double *vz, double *f) {
    double fac = pow(2.0*65.0/PI,0.75);
    for (long i=0; i<n; i++) {
        double x = vx[i]-0.5;
        double y = vy[i]-0.5;
        double z = vz[i]-0.5;
        f[i] = fac*exp(-65.0*(x*x+y*y+z*z));
    }
}

void vdfred_dx(long n, const double *vx, const double *vy, const double *vz, double *f) {
    double fac = pow(2.0*65.0/PI,0.75);
    for (long i=0; i<n; i++) {
        double x = vx[i]-0.5;
        double y = vy[i]-0.5;
        double z = vz[i]-0.5;
        f[i] = fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*x;
    }
}

void vdfred_dy(long n, const double *vx, const double *vy, const double *vz, double *f) {
    double fac = pow(2.0*65.0/PI,0.75);
    for (long i=0; i<n; i++) {
        double x = vx[i]-0.5;
        double y = vy[i]-0.5;
        double z = vz[i]-0.5;
        f[i] = fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*y;
    }
}

void vdfred_dz(long n, const double *vx, const double *vy, const double *vz, double *f) {
    double fac = pow(2.0*65.0/PI,0.75);
    for (long i=0; i<n; i++) {
        double x = vx[i]-0.5;
        double y = vy[i]-0.5;
        double z = vz[i]-0.5;
        f[i] = fac*exp(-65.0*(x*x+y*y+z*z))*-65.0*2.0*z;
    }
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

#ifdef MADNESS_PAPI
class papi {
public:
    static void start() {
        int Events[] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_INS};
        MADNESS_ASSERT(PAPI_start_counters(Events,3) == PAPI_OK);
    };
    static void stop(const char*msg) {
        long_long values[3];
        double d[3];
        MADNESS_ASSERT(PAPI_stop_counters(values,3) == PAPI_OK);
        d[0] = values[0]/2.6e9;
        d[1] = values[1]/d[0]/1e9;
        d[2] = values[2]/d[0]/1e9;
        print("PAPI:  local",msg,d[0],"s",d[1],"gop/s",d[2],"gflop/s");
        if (madness::comm_default->nproc() > 1) {
            madness::comm_default->global_sum(d,3);
            print("PAPI: global", msg, d[0],"s",d[1],"gop/s",d[2],"gflop/s");
        }
    };
};
#else
class papi {
    static double cpu,wall;
public:
    static void start(){cpu = cputime(); wall = walltime();};
    static void stop(const char* msg) {
        double d[2] = {cputime()-cpu,walltime()-wall};
        print("TIMES: local",msg,d[0],d[1]);
        if (madness::comm_default->nproc() > 1) {
            madness::comm_default->global_sum(d,2,MPI::MAX);
            print("TIMES:   max",msg,d[0],d[1]);
        }
     };
};
double papi::cpu, papi::wall;
#endif


int main(int argc, char* argv[]) {

    Communicator& comm = startup(argc,argv);
    FunctionDefaults::k=15;
    FunctionDefaults::initial_level=0;
    Function<double> dummy;

    print("INITIAL TENSOR COUNT",BaseTensor::get_instance_count());


    // To ensure reliable cleanup catch all C++ exceptions here
    try {
        //comm.set_debug(true);
        // Do useful stuff here
        comm.Barrier();
        papi::start();
        Function<double> f = FunctionFactory<double>().vf(vfred).thresh(1e-13).nocompress().refine();
        papi::stop("project + refine of f");
        Function<double> df,dfexact;

        print("valuesX",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        comm.Barrier();
        papi::start();
        f.compress();
        papi::stop("compress of f");

        comm.Barrier();
        papi::start();
        f.reconstruct();
        papi::stop("reconstruct of f");


        print("valuesX",fred(0.45,0.53,0.48),f(0.45,0.53,0.48));
        
        papi::start();
        Function<double> p = f*f;
        papi::stop("squaring");
        print("valuesSQ",fred(0.45,0.53,0.48)*fred(0.45,0.53,0.48),p(0.45,0.53,0.48));
        
        papi::start();
        Function<double> q = p.copy();
        papi::start();
        p.truncate(1e-9);
        papi::stop("truncate");
        print("err after truncate 1",(p-q).norm2sq());
        p.reconstruct();
        p.truncate(1e-5);
        print("err after truncate 2",(p-q).norm2sq());
        p.truncate(1e-3);
        print("err after truncate 3",(p-q).norm2sq());
        
        p = q.copy();       
        q.compress(true);
        q.nsclean(true);
        print("err after NS compress and clean",(p-q).norm2sq());
        
        papi::start();
        df = f.diff(0);
        papi::stop("df/dx");

        papi::start();
        dfexact = FunctionFactory<double>().vf(vdfred_dx).thresh(1e-11).nocompress();
        papi::stop("project and refine of df/dx");

        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff x",df(0.45,0.53,0.48),dfred_dx(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

        df = f.diff(1);
        dfexact = FunctionFactory<double>().vf(vdfred_dy).thresh(1e-11).nocompress();
        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff y",df(0.45,0.53,0.48),dfred_dy(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());

        df = f.diff(2);
        dfexact = FunctionFactory<double>().vf(vdfred_dz).thresh(1e-11).nocompress();
        print("diff norms",df.norm2sq(),dfexact.norm2sq(),f.norm2sq());
        print("diff z",df(0.45,0.53,0.48),dfred_dz(0.45,0.53,0.48),"normerrsq",(df-dfexact).norm2sq());


        print("TENSOR COUNT AT END OF SCOPE",BaseTensor::get_instance_count());

       
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

    // The following should be used for succesful termination
    done:

    print("TENSOR COUNT OUTSIDE OF SCOPE",BaseTensor::get_instance_count());

    comm.close();
    MPI::Finalize();
    return 0;
}

