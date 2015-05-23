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

  Original authors:

  Nichols A. Romero
  Jeff Hammond
  Argonne National Laboratory
  9700 S. Cass Ave.
  Argonne, IL 60439


  email: naromero@alcf.anl.gov, jhammond@alcf.anl.gov

*/

#define MEMINFO_CC

#include <madness/madness_config.h>
#include <iostream>
#include <fstream>
#include <sstream>

#if defined(HAVE_IBMBGQ)
#include <spi/include/kernel/memory.h>
#elif defined(ON_A_MAC) 
#include <malloc/malloc.h>
#elif defined(X86_32) 
#include <malloc.h>
#elif defined(X86_64)
#include <sys/sysinfo.h>
#endif


using namespace std;

/// \file meminfo.cc
/// \brief Memory stats per rank, usually platform specific calls. 
/// Invoke with rank and tag

namespace madness {

  void print_meminfo(int id, const std::string& tag) {
#if defined(MEMINFO)
    ofstream memoryfile;
    ostringstream filename;

    filename << "MEMORY." << id;

    memoryfile.open(filename.str().c_str(), ios::out | ios::app);
    memoryfile << tag << endl;

    const double fac = 1024.0*1024.0; /* Convert from bytes to megabytes */
#if defined(HAVE_IBMBGQ)
    uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

    Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

    memoryfile << "Allocated heap (MB): " << (heap/fac) << ", avail. heap: " << (heapavail/fac) << endl;
    memoryfile << "Allocated stack (MB): " << (stack/fac) << ", avail. stack: " << (stackavail/fac) << endl;
    memoryfile << "Memory: shared: " << (shared/fac) << ", persist: " << (persist/fac) << ", guard: " << (guard/fac) << ", mmap: " << (mmap/fac) << endl;
#elif defined(ON_A_MAC)
  /* Mac OS X specific hack - un-tested post Snow Leopard */
  struct malloc_statistics_t mi; /* structure in bytes */

  malloc_zone_statistics(nullptr, &mi);

  memoryfile << "Allocate heap (MB): " << (mi.size_in_use/fac) << endl; 
#elif defined(X86_32) 
  struct mallinfo mi; /* structure in bytes */
  
  mi = mallinfo();
  
  memoryfile << "Non-mmap (MB): " << (mi.arena/fac) << endl;
  memoryfile << "Mmap (MB): " << (mi.hblkhd/fac) << endl;
  memoryfile << "Total malloc chunks (MB): " << (mi.uordblks/fac) << endl;
#elif defined(X86_64)
  /* Better than nothing, mallinfo unreliable on 64-bit machine due to
     use of int in mallinfo data structure. Requires Linux
     kernel. Inaccurate if other processes besides MADNESS are
     running. Memory differences appear to be reliable. */
  struct sysinfo si; /* structure in bytes */

  sysinfo(&si);

  memoryfile << "Total RAM (MB): " << (si.totalram/fac) << endl;
  memoryfile << "Free RAM (MB): " << (si.freeram/fac) << endl;
  memoryfile << "Buffer (MB): " << (si.bufferram/fac) << endl;
  memoryfile << "RAM in use (MB): " << ((si.totalram - si.freeram + si.bufferram)/fac) << endl;
#endif // platform specific
    memoryfile.close();
#else // MEMINFO
    return;
#endif 
  }
}
