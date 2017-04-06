/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/** 
 \file timers.h
 \brief Wrappers around platform dependent timers and performance info.
 \ingroup parallel_runtime
*/

#ifndef MADNESS_WORLD_TIMERS_H__INCLUDED
#define MADNESS_WORLD_TIMERS_H__INCLUDED

#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <madness/madness_config.h>

#ifdef _CRAY
#include <catamount/dclock.h>
#endif

#ifdef HAVE_IBMBGP
#  define BG_CYCLES_PER_MICROSECOND 850
#  define BG_SECONDS_PER_CYCLE 1.176470588235294033e-09
#  include </bgsys/drivers/ppcfloor/arch/include/bpcore/ppc450_inlines.h>
#endif

#ifdef HAVE_IBMBGQ
#  define BG_CYCLES_PER_MICROSECOND 1600
#  define BG_SECONDS_PER_CYCLE 6.25e-10
#  include <hwi/include/bqc/A2_inlines.h>
#endif


namespace madness {

    /// Returns the wall time in seconds relative to an arbitrary origin.

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the gettime of day system call.
    /// \return The wall time (in seconds).
    double wall_time();

    /// On some machines we have access to a cycle count.

    /// For small intervals this is probably the most lightweight and accurate timer
    /// but may not be meaningful over long intervals due to O/S scheduling,
    /// migration to different cores, frequency shifts, etc.  On x86 uses rtdsc.
    /// Otherwise uses wall_time() in nanoseconds.
    /// \return Timing, in cycle count.
    static inline uint64_t cycle_count() {
        uint64_t x;
#if defined(HAVE_IBMBGP)
        unsigned int rx, ry, rz;
        do
        {
            asm volatile ( "mftbu %0" : "=r"(rx) );
            asm volatile ( "mftb %0" : "=r"(ry) );
            asm volatile ( "mftbu %0" : "=r"(rz) );
        }
        while ( rx != rz );
        x = rx;
        x = (x << 32) | ry;
#elif defined(HAVE_IBMBGQ)
/* Jeff could use the asm above but is pretending this is more portable */
        x = GetTimeBase();
#elif defined(X86_32)
        __asm__ volatile(".byte 0x0f, 0x31" : "=A"(x));
#elif defined(X86_64)
        unsigned int a,d;
        __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));
        x = ((uint64_t)a) | (((uint64_t)d)<<32);
#else
        x = wall_time()*1e9;
#endif
        return x;
    }

    /// Estimate the processor frequency, in Hz.

    /// First call may take about 0.1s to execute. Subsequent
    /// calls return the value, cached from the first call, so it does
    /// not respond to changing processor frequency.
    ///
    /// If \c cycle_count() returns \c wall_time() in nanoseconds,
    /// this will return 1GHz.
    ///
    /// If not available returns 0.
    /// \return CPU frequency, in Hz.
    double cpu_frequency();

    /// Returns the cpu time in seconds relative to an arbitrary origin.

    /// As accurate and lightweight as we can get it, but may not
    /// be any better than the clock system call.
    /// \return The cpu time, in second.
    static inline double cpu_time() {
#if defined(X86_32) || defined(X86_64) || defined(HAVE_IBMBGP)
        static const double rfreq = 1.0/cpu_frequency();
        return cycle_count()*rfreq;
#elif defined(_CRAY)
        return dclock();
#elif defined(HAVE_IBMBGP)
        return BG_SECONDS_PER_CYCLE * _bgp_GetTimeBase();
#elif defined(HAVE_IBMBGQ)
        return BG_SECONDS_PER_CYCLE * GetTimeBase();
#else
        return double(clock())/CLOCKS_PER_SEC;
#endif
    }


    /// Do nothing and especially do not touch memory.

    /// \todo Can we provide some context for this function?
    inline void cpu_relax() {
#if defined(X86_32) || defined(X86_64)
        asm volatile("rep;nop" : : : "memory");
#elif defined(HAVE_IBMBGP) || defined(HAVE_IBMBGQ)
        asm volatile ("nop\n");
#else
        /* Jeff has no idea if this is actually portable.
         * See https://en.wikipedia.org/wiki/NOP for details. */
        asm volatile ("nop\n");
#endif
    }

    /// Sleep or spin for specified number of microseconds.

    /// Wrapper to ensure desired behavior across various platforms.
    /// \param[in] us The number of microseconds.
    static inline void myusleep(unsigned int us) {
#if defined(HAVE_CRAYXT)
        double secs = us*1e-6;
        double start = cpu_time();
        while (cpu_time()-start < secs) {
            for (int i=0; i<100; ++i) cpu_relax();
        }
#elif defined(HAVE_IBMBGP) || defined(HAVE_IBMBGQ)
        int count = BG_CYCLES_PER_MICROSECOND*us; // ??????
        for (int i=0; i<count; i++) {
            asm volatile ("nop\n");
        }
#else
        usleep(us);
#endif
    }
}

#endif // MADNESS_WORLD_TIMERS_H__INCLUDED
