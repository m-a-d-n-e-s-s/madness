#ifndef MADNESS_WORLD_WORLDMEM_H__INCLUDED
#define MADNESS_WORLD_WORLDMEM_H__INCLUDED

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


  $Id: world.cc 454 2008-01-26 03:08:15Z rjharrison $
*/

#include <madness/madness_config.h>
#include <string>
#ifdef WORLD_GATHER_MEM_STATS
#include <new>
#endif // WORLD_GATHER_MEM_STATS
#include <cstddef>
#include <fstream>
#include <sstream>

#if defined(HAVE_IBMBGQ)
#include <spi/include/kernel/memory.h>
#elif defined(ON_A_MAC)
#include <malloc/malloc.h>
#elif defined(X86_32)
#include <malloc.h>
#elif defined(X86_64)
#include <sys/types.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <string>
#endif

namespace madness {

    /// Used to output memory statistics and control tracing, etc.

    /// There is currently only one instance of this class in worldmem.cc
    /// that is used by special versions of global new+delete enabled
    /// by compiling with WORLD_GATHER_MEM_STATS enabled.
    ///
    /// Default for max_mem_limit is unlimited.
    class WorldMemInfo {
#ifdef WORLD_GATHER_MEM_STATS
        friend void* ::operator new(size_t size) throw (std::bad_alloc);
        friend void ::operator delete(void *p) throw();
#endif
    private:
        /// Invoked when user pointer p is allocated with size bytes
        void do_new(void *p, std::size_t size);

        /// Invoked when user pointer p is deleted with size bytes
        void do_del(void *p, std::size_t size);

    public:
        static const unsigned long overhead = 4*sizeof(long) + 16;
        /// If you add new stats be sure that the initialization in worldmem.cc is OK
        unsigned long num_new_calls;   ///< Counts calls to new
        unsigned long num_del_calls;   ///< Counts calls to delete
        unsigned long cur_num_frags;   ///< Current number of allocated fragments
        unsigned long max_num_frags;   ///< Lifetime maximum number of allocated fragments
        unsigned long cur_num_bytes;   ///< Current amount of allocated memory in bytes
        unsigned long max_num_bytes;   ///< Lifetime maximum number of allocated bytes
        unsigned long max_mem_limit;   ///< if size+cur_num_bytes>max_mem_limit new will throw MadnessException
        bool trace;

        /// Prints memory use statistics to std::cout
        void print() const;

        /// Resets all counters to zero
        void reset();

        /// If trace is set true a message is printed for every new and delete
        void set_trace(bool trace) {
            this->trace = trace;
        }

        /// Set the maximum memory limit (trying to allocate more will throw MadnessException)
        void set_max_mem_limit(unsigned long max_mem_limit) {
            this->max_mem_limit = max_mem_limit;
        }
    };

    /// Returns pointer to internal structure
    WorldMemInfo* world_mem_info();

    namespace detail {
      template <typename Char> const Char* Vm_cstr();

      bool& print_meminfo_flag_accessor();
    }

    /// \name print_meminfo
    /// These functions are useful to instrument *aggregate* (not fine-grained, like WorldMemInfo) memory profiling.
    /// \note To make this functionality usable must set the `ENABLE_MEM_PROFILE` CMake cache variable to `ON`
    ///@{

    /// disables print_meminfo() profiling (i.e. calling it is a no-op)
    void print_meminfo_disable();
    /// enables print_meminfo() profiling
    /// \note print_meminfo() profiling is on by default.
    void print_meminfo_enable();
    /// @return true if print_meminfo() will generate profiling data.
    bool print_meminfo_enabled();

    /// \brief print aggregate memory stats to file \c filename_prefix.<rank> , tagged with \c tag
    /// \param[in] rank process rank
    /// \param[in] tag record tag as any string type, e.g. \c const char[] , \c std::string , or \c std::wstring
    /// \param[in] filename_prefix file name prefix; the default value is "MEMORY"
    /// \note To make print_meminfo() usable must set the ENABLE_MEM_PROFILE CMake cache variable to ON; this makes print_meminfo() enabled by default.
    ///       To disable/enable programmatically use print_meminfo_disable()/print_meminfo_enable() .
    /// \note must set global locale properly with \c std::locale::global() if \c tag has
    ///       nontrivial encoding
    /// \warning this does not fence, it is up to the user to ensure proper synchronization
    template <typename String> void print_meminfo(
        int rank, const String& tag,
        const std::string filename_prefix = std::string("MEMORY")) {
      using namespace std;
#if defined(WORLD_MEM_PROFILE_ENABLE)
      if (print_meminfo_enabled()) {
        using Char = typename std::iterator_traits<decltype(
            std::begin(tag))>::value_type;
        basic_ofstream<Char> memoryfile;
        ostringstream filename;

        filename << filename_prefix << "." << rank;

        memoryfile.open(filename.str().c_str(), ios::out | ios::app);
        memoryfile << tag << endl;

        const double to_MiB =
            1 / (1024.0 * 1024.0); /* Convert from bytes to MiB */
#if defined(HAVE_IBMBGQ)
        uint64_t shared, persist, heapavail, stackavail, stack, heap, guard,
            mmap;

        Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
        Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

        memoryfile << "Heap size (MiB): " << (heap * to_MiB)
                   << ", available: " << (heapavail * to_MiB) << endl;
        memoryfile << "Stack size (MiB): " << (stack * to_MiB)
                   << ", available: " << (stackavail * to_MiB) << endl;
        memoryfile << "Memory (MiB): shared: " << (shared * to_MiB)
                   << ", persist: " << (persist * to_MiB)
                   << ", guard: " << (guard * to_MiB) << ", mmap: " << (mmap * to_MiB)
                   << endl;
#elif defined(ON_A_MAC)
        /* Mac OS X specific hack - un-tested post Snow Leopard */
        struct malloc_statistics_t mi; /* structure in bytes */

        malloc_zone_statistics(nullptr, &mi);

        memoryfile << "Heap allocated  (MiB): " << (mi.size_allocated * to_MiB) << endl;
        memoryfile << "Heap used       (MiB): " << (mi.size_in_use * to_MiB) << endl;
        memoryfile << "Heap max used   (MiB): " << (mi.max_size_in_use * to_MiB) << endl;
#elif defined(X86_32) // 32-bit Linux
        struct mallinfo mi; /* structure in bytes */

        mi = mallinfo();

        memoryfile << "Non-mmap (MiB): " << (mi.arena * to_MiB) << endl;
        memoryfile << "Mmap (MiB): " << (mi.hblkhd * to_MiB) << endl;
        memoryfile << "Total malloc chunks (MiB): " << (mi.uordblks * to_MiB)
                   << endl;
#elif defined(X86_64) // 64-bit Linux
        // try parsing /proc/PID/status first, fallback on sysinfo
        string status_fname =
            string("/proc/") + to_string(getpid()) + string("/status");
        basic_ifstream<Char> status_stream(status_fname);
        if (status_stream.good()) {
          basic_string<Char> line;
          while (getline(status_stream, line)) {
            if (line.find(detail::Vm_cstr<Char>()) == 0)
              memoryfile << line << endl;
          }
          status_stream.close();
        } else {
          /* Better than nothing, mallinfo unreliable on 64-bit machine due to
             use of int in mallinfo data structure. Requires Linux
             kernel. Inaccurate if other processes besides MADNESS are
             running. Memory differences appear to be reliable. */
          struct sysinfo si; /* structure in bytes */

          sysinfo(&si);

          memoryfile << "Total RAM (MiB): " << (si.totalram * to_MiB) << endl;
          memoryfile << "Free RAM (MiB): " << (si.freeram * to_MiB) << endl;
          memoryfile << "Buffer (MiB): " << (si.bufferram * to_MiB) << endl;
          memoryfile << "RAM in use (MiB): "
                     << ((si.totalram - si.freeram + si.bufferram) * to_MiB)
                     << endl;
        }
#endif                // platform specific
        memoryfile.close();
      }
#else   // WORLD_MEM_PROFILE_ENABLE
      return;
#endif  // WORLD_MEM_PROFILE_ENABLE
    }

    ///@}

}  // namespace madness

#endif // MADNESS_WORLD_WORLDMEM_H__INCLUDED
