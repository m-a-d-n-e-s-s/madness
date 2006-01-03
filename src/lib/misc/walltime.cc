#include <time.h>
#include <sys/time.h>

/// \file walltime.cc
/// \brief Implements wall and cpu timers

namespace madness {

    /// Returns the wall time in seconds relative to arbitrary origin 
    double walltime() {
#ifdef USE_GETTIMEOFDAY
        struct timeval tv;
        gettimeofday(&tv,0);
        return tv.tv_sec + 1e-6*tv.tv_usec;
#else
        return 0.0;
#endif
    }

    /// Returns the cpu time in seconds relative to arbitrary origin (may wrap)
    double cputime() {
        return clock()*1e-6;
    }

}
