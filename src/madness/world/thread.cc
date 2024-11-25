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
 \file thread.cc
 \brief Implements Dqueue, Thread, ThreadBase and ThreadPool.
 \ingroup threads
*/

#include <madness/world/worldinit.h>
#include <madness/world/thread.h>
#include <madness/world/worldprofile.h>
#include <madness/world/madness_exception.h>
#include <madness/world/print.h>
#include <madness/world/worldpapi.h>
#include <madness/world/safempi.h>
#include <madness/world/atomicint.h>
#include <cstring>
#include <fstream>

#if defined(HAVE_IBMBGQ) and defined(HPM)
extern "C" unsigned int HPM_Prof_init_thread(void);
extern "C" void HPM_Prof_start(unsigned int);
extern "C" void HPM_Prof_stop(unsigned int);
#endif

#if defined(HAVE_IBMBGP)
// This header causes tinyxml.h to barf but we only need it in the implementation, not the header.
#  include <spi/kernel_interface.h>
#  include <spi/bgp_kernel_inlines.h>
#  include <common/bgp_personality.h>
#  include <common/bgp_personality_inlines.h>
#endif

#if defined(HAVE_IBMBGQ)
#  include <spi/include/kernel/location.h>
#  include <spi/include/kernel/process.h>
#endif

#include <thread>

namespace madness {
    ThreadBinder binder;
    thread_local bool ThreadBinder::bound = false;
    pthread_key_t ThreadBase::thread_key;

    ThreadPool* ThreadPool::instance_ptr = 0;
    double ThreadPool::await_timeout = 900.0;
#if HAVE_INTEL_TBB
    std::unique_ptr<tbb::global_control> ThreadPool::tbb_control    = nullptr;
    std::unique_ptr<tbb::task_arena> ThreadPool::tbb_arena          = nullptr;
#endif
#ifdef MADNESS_TASK_PROFILING
    Mutex profiling::TaskProfiler::output_mutex_;
    const char* profiling::TaskProfiler::output_file_name_;
#endif // MADNESS_TASK_PROFILING
#if defined(HAVE_IBMBGQ) and defined(HPM)
    unsigned int ThreadPool::main_hpmctx;
    bool ThreadBase::main_instrumented;
    bool ThreadBase::all_instrumented;
    int ThreadBase::hpm_thread_id;
#endif

    void* ThreadBase::main(void* self) {
#ifdef HAVE_PAPI
        begin_papi_measurement();
#endif
#if defined(HAVE_IBMBGQ) and defined(HPM)
	unsigned int slave_hpmctx; // HPM context for the slave threads
	int pool_num = static_cast<ThreadBase*>(self)->pool_num;
	// int all_instrumented = static_cast<ThreadBase*>(self)->all_instrumented;
	// int hpm_thread_id = static_cast<ThreadBase*>(self)->hpm_thread_id;
	bool this_slave_instrumented;

	if ((hpm_thread_id == pool_num) || all_instrumented) {
	  this_slave_instrumented = true;
	} else
	  this_slave_instrumented = false;

	if (this_slave_instrumented) {
	  slave_hpmctx = HPM_Prof_init_thread();
	  HPM_Prof_start(slave_hpmctx);
	}
#endif
	const int rc = pthread_setspecific(thread_key, self);
        if(rc != 0)
            MADNESS_EXCEPTION("pthread_setspecific failed", rc);

        // mark as a MADNESS thread
        set_thread_tag(ThreadTag_MADNESS);

        try {
            ((ThreadBase*)(self))->run();
        }
        catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        }
        catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        }
        catch (const char* s) {
            print(s);
            error("caught a string exception");
        }
        catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        }
        catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        }
        catch (...) {
            error("caught unhandled exception");
        }

#ifdef HAVE_PAPI
        end_papi_measurement();
#endif

#if defined(HAVE_IBMBGQ) and defined(HPM)
	if (this_slave_instrumented) HPM_Prof_stop(slave_hpmctx);
#endif
        return 0;
    }

    // Start the thread running
    void ThreadBase::start() {
        pthread_attr_t attr;
        // Want detached thread with kernel scheduling so can use multiple cpus
        // RJH ... why not checking success/failure????
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
#ifndef HAVE_IBMBGP
        pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif
        int result = pthread_create(&id, &attr, &ThreadBase::main, (void *) this);
        if (result) MADNESS_EXCEPTION("failed creating thread", result);

        pthread_attr_destroy(&attr);
    }

    // Get no. of actual hardware processors
    int ThreadBase::num_hw_processors() {
#if defined(HAVE_IBMBGP)
    #if 0 /* total overkill - what was i thinking? */
        int ncpu=0;
        _BGP_Personality_t pers;
        Kernel_GetPersonality(&pers, sizeof(pers));
        if      ( BGP_Personality_processConfig(&pers) == _BGP_PERS_PROCESSCONFIG_SMP ) ncpu = 4;
        else if ( BGP_Personality_processConfig(&pers) == _BGP_PERS_PROCESSCONFIG_2x2 ) ncpu = 2;
        else if ( BGP_Personality_processConfig(&pers) == _BGP_PERS_PROCESSCONFIG_VNM ) ncpu = 1;
        return ncpu;
    #else
        /* Returns the number of Processes (Virtual Nodes) running on this Physical Node. */
        return 4/Kernel_ProcessCount();
    #endif
#elif defined(HAVE_IBMBGQ)
        /* Return number of processors (hardware threads) within the current process. */
        return Kernel_ProcessorCount();
#elif defined(_SC_NPROCESSORS_CONF)
        int ncpu = sysconf(_SC_NPROCESSORS_CONF);
        if (ncpu <= 0)
           MADNESS_EXCEPTION("ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)", ncpu);
        return ncpu;
#elif defined(HC_NCPU)
        int mib[2]={CTL_HW,HW_NCPU};
        size_t len = sizeof(ncpu);
        if (sysctl(mib, 2, &ncpu, &len, nullptr, 0) != 0)
            MADNESS_EXCEPTION("ThreadBase: sysctl(CTL_HW,HW_NCPU) failed", 0);
        //std::cout << "NCPU " << ncpu << std::endl;
#else
        return std::thread::hardware_concurrency();
	//return get_nprocs();
#endif
    }


#ifdef MADNESS_TASK_PROFILING

    namespace profiling {

        void TaskProfiler::write_to_file() {
            // Get output filename: NAME_[rank]x[threads + 1]
            if(output_file_name_ != nullptr) {
                // Construct the actual output filename
                std::stringstream file_name;
                file_name << output_file_name_ << "_"
                        << SafeMPI::COMM_WORLD.Get_rank() << "x"
                        << ThreadPool::size() + 1;

                // Lock file for output
                ScopedMutex<Mutex> locker(TaskProfiler::output_mutex_);

                // Open the file for output
                std::ofstream file(file_name.str().c_str(), std::ios_base::out | std::ios_base::app);
                if(! file.fail()) {
                    // Print the task profile data
                    // and delete the data since it is not needed anymore
                    const TaskEventListBase* next = nullptr;
                    while(head_ != nullptr) {
                        next = head_->next();
                        file << *head_;
                        delete head_;
                        head_ = const_cast<TaskEventListBase*>(next);
                    }

                    tail_ = nullptr;
                } else {
                    std::cerr << "!!! ERROR: TaskProfiler cannot open file: "
                            << file_name.str() << "\n";
                }

                // close the file
                file.close();
            } else {
                // Nothing is written so just cleanup data
                const TaskEventListBase* next = nullptr;
                while(head_ != nullptr) {
                    next = head_->next();
                    delete head_;
                    head_ = const_cast<TaskEventListBase*>(next);
                }

                tail_ = nullptr;
            }
        }


    } // namespace profiling

#endif // MADNESS_TASK_PROFILING

#if HAVE_PARSEC
    ParsecRuntime *ThreadPool::parsec_runtime = nullptr;
#endif

    // The constructor is private to enforce the singleton model
    ThreadPool::ThreadPool(int nthread)
    : threads(nullptr)
    , main_thread()
    , nthreads(nthread)
    , finish(false)
    {
        nfinished = 0;
        instance_ptr = this;
        if (nthreads < 0) nthreads = default_nthread();
        MADNESS_ASSERT(nthreads >=0);

        const int rc = pthread_setspecific(ThreadBase::thread_key,
                static_cast<void*>(&main_thread));
        if(rc != 0)
            MADNESS_EXCEPTION("pthread_setspecific failed", rc);
#if HAVE_PARSEC
        MADNESS_ASSERT(nullptr == parsec_runtime);
        const int num_parsec_threads = nthreads + 1;
        parsec_runtime = new ParsecRuntime(num_parsec_threads);
#elif HAVE_INTEL_TBB
// #if HAVE_INTEL_TBB

        if(nthreads < 1)
            nthreads = 1;

        //  nranks > 1: there are nthreads+2 because the main and communicator thread
        //              are counted as part of tbb.
        // nranks == 1: there are nthreads+1 because the main
        //              is counted as part of tbb.
        const int num_tbb_threads = (SafeMPI::COMM_WORLD.Get_size() > 1) ? nthreads + 2 : nthreads + 1;
        tbb_control = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, num_tbb_threads);
        tbb_arena   = std::make_unique<tbb::task_arena>(num_tbb_threads);
#else

        if (nthreads>64)
            MADNESS_EXCEPTION("When configured with MADNESS_TASK_BACKEND=Pthreads MAD_NUM_THREADS cannot exceed 64",1);

        try {
            if (nthreads > 0)
                threads = new ThreadPoolThread[nthreads];
            else
                threads = 0;
        }
        catch (...) {
            MADNESS_EXCEPTION("memory allocation failed", 0);
        }

        for (int i=0; i<nthreads; ++i) {
            threads[i].set_pool_thread_index(i);
            threads[i].start(pool_thread_main, (void *)(threads+i));
        }
#endif
        /****************************/
    }

    // Get number of threads from the environment
    int ThreadPool::default_nthread() {
        int nthread;
        int shift = 0;
        char *cnthread = getenv("MAD_NUM_THREADS");
        // MAD_NUM_THREADS is total no. of application threads whereas
        // POOL_NTHREAD is just the number in the pool (one less)
        if (cnthread) shift = 1;
        if (cnthread == 0) cnthread = getenv("POOL_NTHREAD");

        if (cnthread) {
            int result = sscanf(cnthread, "%d", &nthread);
            if (result != 1)
                MADNESS_EXCEPTION("POOL_NTHREAD is not an integer", result);
            nthread -= shift;
        }
        else {
            nthread = ThreadBase::num_hw_processors();
            if (nthread < 2)
                nthread = 2;
            nthread = nthread - 1; // One less than # physical processors
        }
        return nthread;
    }

    void ThreadPool::thread_main(ThreadPoolThread* const thread) {
        PROFILE_MEMBER_FUNC(ThreadPool);
	binder.bind();

#if !HAVE_PARSEC
#define MULTITASK
#ifdef  MULTITASK
        while (!finish) {
            run_tasks(true, thread);
        }
#else
        while (!finish) {
            run_task(true, thread);
        }
#endif
#endif

#ifdef MADNESS_TASK_PROFILING
        thread->profiler().write_to_file();
#endif // MADNESS_TASK_PROFILING

        nfinished++;
    }

    // Forwards thread to bound member function
    void* ThreadPool::pool_thread_main(void *v) {
        instance()->thread_main((ThreadPoolThread*)(v));
        return 0;
    }

    void ThreadPool::begin(int nthread) {
        // Check that the singleton has not been previously initialized
        if(instance_ptr) return;

        ThreadBase::init_thread_key();

        // Construct the thread pool singleton
        instance_ptr = new ThreadPool(nthread);

        const char* mad_wait_timeout = getenv("MAD_WAIT_TIMEOUT");
        if(mad_wait_timeout) {
            std::stringstream ss(mad_wait_timeout);
            ss >> await_timeout;
            if(await_timeout < 0.0) {
                if(SafeMPI::COMM_WORLD.Get_rank() == 0 && !madness::quiet())
                    std::cout << "!!MADNESS WARNING: Invalid wait timeout.\n"
                              << "!!MADNESS WARNING: MAD_WAIT_TIMEOUT = " << mad_wait_timeout << "\n";
                await_timeout = 900.0;
            }
            if(SafeMPI::COMM_WORLD.Get_rank() == 0 && !madness::quiet()) {
                if(await_timeout >= 1.0) {
                    std::cout << "MADNESS wait timeout set to " << await_timeout << " seconds.\n";
                } else {
                    std::cout << "MADNESS wait timeout disabled.\n";
                }
            }
        }

#ifdef MADNESS_TASK_PROFILING
        // Initialize the output file name for the task profiler.
        profiling::TaskProfiler::output_file_name_ =
                getenv("MAD_TASKPROFILER_NAME");
        if(! profiling::TaskProfiler::output_file_name_) {
            if(SafeMPI::COMM_WORLD.Get_rank() == 0)
                std::cerr
                    << "!!! WARNING: MAD_TASKPROFILER_NAME not set.\n"
                    << "!!! WARNING: There will be no task profile output.\n";
        } else {
            // Construct the actual output filename
            std::stringstream file_name;
            file_name << profiling::TaskProfiler::output_file_name_ << "_"
                    << SafeMPI::COMM_WORLD.Get_rank() << "x"
                    << ThreadPool::size() + 1;

            // Erase the profiler output file
            std::ofstream file(file_name.str().c_str(), std::ios_base::out | std::ios_base::trunc);
            file.close();
        }
#endif  // MADNESS_TASK_PROFILING

#if defined(HAVE_IBMBGQ) and defined(HPM)
	if (ThreadBase::main_instrumented) {
	  main_hpmctx = HPM_Prof_init_thread();
	  HPM_Prof_start(main_hpmctx);
	}
#endif
    }

    void ThreadPool::end() {
#if !HAVE_INTEL_TBB
        if (!instance_ptr) return;
        instance()->finish = true;
#if !HAVE_PARSEC
        for (int i=0; i<instance()->nthreads; ++i) {
            add(new PoolTaskNull);
        }
	instance_ptr->flush_prebuf();
        while (instance_ptr->nfinished != instance_ptr->nthreads);
#else  /* HAVE_PARSEC */
        /* Remove the fake task we used to keep the engine up and running */
        parsec_runtime->wait();
#endif
#ifdef MADNESS_TASK_PROFILING
        instance_ptr->main_thread.profiler().write_to_file();
#endif // MADNESS_TASK_PROFILING

        ThreadBase::delete_thread_key();
#endif

#if defined(HAVE_IBMBGQ) and defined(HPM)
	if (ThreadBase::main_instrumented) HPM_Prof_stop(main_hpmctx);
#endif
        delete instance_ptr;
        instance_ptr = nullptr;
    }

    // Returns queue statistics
    const DQStats& ThreadPool::get_stats() {
        return instance()->queue.get_stats();
    }

} // namespace madness
