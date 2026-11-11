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
*/

/// \file test_foreign_threads.cc
/// \brief Tasks submitted concurrently from threads the task backend did not
///        create.
///
/// The MADNESS comm thread and other runtimes' callback threads (e.g. CUDA
/// host-function callbacks) submit tasks from outside the backend's thread
/// pool. With the PaRSEC backend each such thread needs an execution stream.
/// In builds with PaRSEC tracing (PARSEC_PROF_TRACE) the task profiler records
/// every schedule call in that stream, so this test turns the profiler on.

#include <madness/world/MADworld.h>
#ifdef HAVE_PARSEC
#include <madness/world/parsec.h>
#endif

#include <atomic>
#include <barrier>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

using namespace madness;

int main(int argc, char** argv) {
#if defined(HAVE_PARSEC) && defined(PARSEC_PROF_TRACE)
    // trace to a scratch file; overwrite=0 keeps a user's own setting
    const auto prof = std::filesystem::temp_directory_path() /
                      ("test_foreign_threads_" + std::to_string(getpid()));
    setenv("PARSEC_MCA_profile_filename", prof.c_str(), 0);
    setenv("PARSEC_MCA_mca_pins", "task_profiler", 0);
#endif
    World& world = initialize(argc, argv);

    const int nthreads = 4;
    const long ntasks = 20000;
    std::atomic<long> ndone{0};
    std::atomic<int> nleft_es{0};  // threads left with a PaRSEC execution stream installed
    std::barrier start(nthreads);
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        threads.emplace_back([&] {
            start.arrive_and_wait();  // all threads submit at the same time
            for (long t = 0; t < ntasks; ++t) world.taskq.add([&ndone] { ++ndone; });
#ifdef HAVE_PARSEC
            // the thread may belong to another runtime and outlive MADNESS:
            // scheduling must not leave PaRSEC state behind on it
            if (parsec_my_execution_stream() != nullptr) ++nleft_es;
#endif
        });
    }
    for (auto& t : threads) t.join();
    world.taskq.fence();

    MADNESS_CHECK(ndone == nthreads * ntasks);
    MADNESS_CHECK(nleft_es == 0);
    finalize();

#if defined(HAVE_PARSEC) && defined(PARSEC_PROF_TRACE)
    for (const auto& f : std::filesystem::directory_iterator(prof.parent_path()))
        if (f.path().filename().string().starts_with(prof.filename().string()))
            std::filesystem::remove(f.path());
#endif
    return 0;
}
