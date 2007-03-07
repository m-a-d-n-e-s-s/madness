#ifndef WORLD_TIME_H
#define WORLD_TIME_H

#include <time.h>
#include <sys/time.h>

/// \file worldtime.h
/// \brief Wrappers around platform dependent timers and performance info


namespace madness {

    /// Returns the wall time in seconds relative to arbitrary origin

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the gettime of day system call.
    double wall_time();

    /// Returns the cpu time in seconds relative to arbitrary origin

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the clock system call.
    double cpu_time();


    /// On some machines we have access to a cycle count
    
    /// For small intervals this is probably the most lightweight and accurate timer
    /// but may not be meaningful over long intervals due to O/S scheduling,
    /// migration to different cores, frequency shifts, etc.  On x86 uses rtdsc.
    /// Otherwise uses wall_time() in nanoseconds.
    static inline uint64_t cycle_count() {
        uint64_t x;
#ifdef X8632
        __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
#elif defined(X8664)
     unsigned int a,d;
     __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
      x = ((uint64_t)a) | (((uint64_t)d)<<32); 
#else
        x = wall_time()*1e9;
#endif
        return x;
    }

    /// Estimates frequency of the processor in Hz

    /// First call may take about 0.1s to execute.  Subsequent
    /// calls return value cached from the first call so does
    /// not respond to changing processor frequency.
    ///
    /// If cycle_count() is returning wall_time() in nanoseconds
    /// this will return 1GHz.
    ///
    /// If not available returns 0.
    double cpu_frequency();

}


#endif
