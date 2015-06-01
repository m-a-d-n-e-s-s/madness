/*
  This file is part of MADNESS.

  Copyright (C) 2013, Argonne National Laboratory
  All rights reserved.

  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the 
  following conditions are met:

  * Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, 
    this list of conditions and the following disclaimer in the documentation 
    and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY 
  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
  DAMAGE.

  Original author:

  Jeff R. Hammond
  Argonne National Laboratory
  9700 S. Cass Ave.
  Argonne, IL 60439

  email: jhammond@alcf.anl.gov

*/

#include <madness/world/hardware.h>

namespace madness {

    void Hardware::Initialize(void) {
#if defined(__bgp__) || defined(__bgq__)
        int init = 0;
        MPI_Is_initialized(&init);
        if (init<1)
            MADNESS_EXCEPTION("MPI is not initialized!");
#  if defined(__bgp__)
        /* We need MPI to initialize DCMF so DCMF_Hardware will work.
         * Alternatively, we must initialize DCMF ourselves, which is fine since the init call is not a singleton like MPI. */
        DCMF_Hardware(&bgphw);
#  elif defined(__bgq__)
        /* MPIX requires MPI, obviously. */
        MPIX_Hardware(&bgqhw);
#  endif
#endif
        Hardware::InitializeHWThreads();
        Hardware::InitializeCPUFrequency();
        Hardware::InitializeMemorySize();

        return;
    }

    void Hardware::Print(std::ostream & out) {
        out.flush();
        out << "\n    MADNESS Hardware Information \n";
        out << "    --------------------------------\n";
    }

    /* Taken from ThreadBase.
     * Jeff rewrote the parts he didn't contribute because of licensing issues. */
    void Hardware::InitializeHWThreads(void) {
        int n=1;
#if defined(__bgp__)
        n = 4/bghw.tSize;
#elif defined(__bgq__)
        n = 64/bghw.ppn;
#elif defined(HARDWARE_USE_SYSCTL)
        size_t len = sizeof(n);
#  if defined(__APPLE__) && defined(__MACH__)
        /* Apple has deprecated HW_NCPU */
        int rc = sysctlbyname("hw.logicalcpu", &n, &len, nullptr, 0);
        if (rc!=0) 
            MADNESS_EXCEPTION("sysctlbyname failed");
#  else
        int mib[2] = {CTL_HW, HW_NCPU};
        int rc = sysctl(mib, 2, &n, &len, nullptr, 0);
        if (rc!=0) 
            MADNESS_EXCEPTION("sysctl failed");
#  endif
#elif defined(HARDWARE_USE_SYSCONF)
        n = (int)sysconf(_SC_NPROCESSORS_CONF);
        if (n<1) 
            MADNESS_EXCEPTION("sysconf failed");
#endif
        hwthreads = n;
        return;
    }

    void Hardware::InitializeCPUFrequency(void) {
        double f=1.0e9;
#if defined(__bgp__) || defined(__bgq__)
        f = bghw.clockMHz*1.0e6;
#elif defined(HARDWARE_USE_SYSCTL)
#  if defined(__APPLE__) && defined(__MACH__)
        long freq = 0;
        size_t len = sizeof(freq);
        int rc = sysctlbyname("hw.cpufrequency", &freq, &len, nullptr, 0);
        f = (double)freq;
#  else
        /* this does not work but I don't have an alternative for BSD */
        struct clockinfo c;
        int mib[2] = {CTL_KERN, KERN_CLOCKRATE};
        size_t len = sizeof(c);
        int rc = sysctl(mib, 2, &c, &len, nullptr, 0);
        f = (double) c.hz;
#  endif
        if (rc!=0) 
            MADNESS_EXCEPTION("sysctl failed");
#elif defined(HARDWARE_USE_SYSCONF)
        long cps = sysconf(_SC_CLK_TCK);
        f = (double)cps;
        if (cps<1) 
            MADNESS_EXCEPTION("sysconf failed");
#endif
        cpufrequency = f;
    }

    void Hardware::InitializeMemorySize(void) {
        int64_t m=1;
#if defined(__bgp__)
        m = bghw.memSize;
#elif defined(__bgq__)
        m = bghw.memSize;
#elif defined(HARDWARE_USE_SYSCTL)
        size_t len = sizeof(m);
        int mib[2] = {CTL_HW, HW_MEMSIZE};
        int rc = sysctl(mib, 2, &m, &len, nullptr, 0);
        if (rc!=0) 
            MADNESS_EXCEPTION("sysctl failed");
#elif defined(HARDWARE_USE_SYSCONF)
        long np = sysconf(_SC_PHYS_PAGES);
        long ps = sysconf(_SC_PAGESIZE);
        m = np*ps;
        if (np<1 || ps<1) 
            MADNESS_EXCEPTION("sysconf failed");
#endif
        memorysize = m;
        return;
    }



} // namespace madness






