#ifdef MADNESS_PAPI
#include <papi.h>

namespace madness {

class papi {
public:
    static void start() {
        static int Events[] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_INS};
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
#include <misc/misc.h>
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
}
