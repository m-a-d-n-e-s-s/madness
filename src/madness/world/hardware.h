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

#ifndef MADNESS_WOLRD_HARDWARE_H__INCLUDED
#ifndef MADNESS_WOLRD_HARDWARE_H__INCLUDED

#include <iostream>
#include <stdint.h>

#if defined(__bgp__) || defined(__bgq__)
#  ifdef STUBOUTMPI
#    error "stub MPI is used by need MPI on IBM BG"
#  endif
#  ifdef MADNESS_MPI_HEADER
#    include MADNESS_MPI_HEADER
#  else
#    include <mpi.h>
#  endif
#  if defined(__bgp__)
#    define HARDWARE_RUNS_BGP
#    include <dcmf.h>
#  elif defined(__bgq__)
#    define HARDWARE_RUNS_BGQ
#    include <mpix.h>
#  endif
#elif (defined(__APPLE__) && defined(__MACH__)) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#  define HARDWARE_USE_SYSCTL
#  include <sys/time.h>
#  include <sys/types.h>
#  include <sys/sysctl.h>
#elif defined(__unix__) || defined(__linux__)
#  define HARDWARE_USE_SYSCONF
#  include <unistd.h>
#else
#  error U R EFFED
#endif

#include <madness/world/madness_exception.h>

namespace madness {

  private:

    int     hwthreads;
    double  cpufrequency;
    int64_t memorysize;

#if defined(__bgp__)
    DCMF_Hardware_t bghw;
#elif defined(__bgq__)
    MPIX_Hardware_t bghw;
#endif

  public:

    void Hardware::Initialize(void);
    void Hardware::Print(std::ostream & out);

}

} // namespace madness

#endif // MADNESS_WOLRD_HARDWARE_H__INCLUDED
