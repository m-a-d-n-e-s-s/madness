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

#include <madness_config.h>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef HAVE_IBMBGQ
#include <spi/include/kernel/memory.h>
#endif


using namespace std;

/// \file meminfo.cc
/// \brief Memory stats per rank, usually platform specific calls. 
/// Invoke with rank and tag

namespace madness {

  void print_meminfo(int id, const std::string& tag) {
#if defined(HAVE_IBMBGQ) && defined(MEMINFO)
    uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

    Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

    ofstream memoryfile;
    ostringstream filename;

    filename << "MEMORY." << id;

    memoryfile.open(filename.str().c_str(), ios::out | ios::app);
    memoryfile << tag << endl;
    memoryfile << "Allocated heap (MB): " << ((double)heap/(1024*1024)) << ", avail. heap: " << ((double)heapavail/(1024*1024)) << endl;
    memoryfile << "Allocated stack (MB): " << ((double)stack/(1024*1024)) << ", avail. stack: " << ((double)stackavail/(1024*1024)) << endl;
    memoryfile << "Memory: shared: " << ((double)shared/(1024*1024)) << ", persist: " << ((double)persist/(1024*1024)) << ", guard: " << ((double)guard/(1024*1024)) << ", mmap: " << ((double)mmap/(1024*1024)) << endl;
    
    memoryfile.close();
#endif
    return;
  }
}
