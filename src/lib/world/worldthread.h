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

  $Id$
*/
#ifndef MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
#define MADNESS_WORLD_WORLDTHREAD_H__INCLUDED

/// \file worldthread.h
/// \brief Implements Dqueue, Thread, ThreadBase and ThreadPool

#include <world/safempi.h>
#include <world/worldexc.h>
#include <world/print.h>
#include <world/worldmutex.h>
#include <world/worldpapi.h>
#include <world/worldprofile.h>
#include <world/atomicint.h>

#ifndef _SC_NPROCESSORS_CONF
// Old macs don't have necessary support thru sysconf to determine the
// no. of processors so must use sysctl
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace madness {

    void error(const char *msg);

    struct DQStats { // Dilly bar, blizzard, ...
        uint64_t npush_back;    //< #calls to push_back
        uint64_t npush_front;   //< #calls to push_front
        uint64_t npop_front;    //< #calls to pop_front
        uint64_t ngrow;         //< #calls to grow
        uint64_t nmax;          //< Lifetime max. entries in the queue

        DQStats()
                : npush_back(0), npush_front(0), npop_front(0), ngrow(0), nmax(0) {}
    };


    /// A thread safe, fast but simple doubled-ended queue.

    /// Since the point is speed, the implementation is a circular
    /// buffer rather than a linked list so as to avoid the new/del
    /// overhead.  It will grow as needed, but presently will not
    /// shrink.  Had to modify STL API to make things thread safe.
    ///
    /// It is now rather heavily specialized to its only use.
    template <typename T>
    class DQueue : private CONDITION_VARIABLE_TYPE {
        char pad[64]; // To put the lock and the data in separate cache lines
        volatile size_t n __attribute__((aligned(64)));        ///< Number of elements in the buffer
        volatile size_t sz;              ///< Current capacity
        volatile T* volatile buf;        ///< Actual buffer
        volatile int _front;  ///< Index of element at front of buffer
        volatile int _back;    ///< Index of element at back of buffer
        DQStats stats;

        void grow() {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            stats.ngrow++;
            if (sz != n) MADNESS_EXCEPTION("assertion failure in dqueue::grow", static_cast<int>(sz));
            size_t oldsz = sz;
            if (sz < 32768)
                sz = 65536;
            else if (sz <= 1048576)
                sz *= 2;
            else
                sz += 1048576;
            volatile T* volatile nbuf = new T[sz];
            int lo = sz/2 - oldsz/2;
            for (int i=_front; i<int(oldsz); i++,lo++) {
                nbuf[lo] = buf[i];
            }
            if (_front > 0) {
                for (int i=0; i<=_back; i++,lo++) {
                    nbuf[lo] = buf[i];
                }
            }
            _front = sz/2 - oldsz/2;
            _back = _front + n - 1;
            buf = nbuf;
            //sanity_check();
        }

        void sanity_check() const {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            int num = _back - _front + 1;
            if (num < 0) num += sz;
            if (num==int(sz) && n==0) num=0;
            if (num==0 && n==sz) num=sz;
            //if (long(n) != num) print("size",sz,"front",_front,"back",_back,"n",n,"num",num);
            MADNESS_ASSERT(long(n) == num);
        }

        void push_back_with_lock(const T& value) {
            size_t nn = n;
            size_t ss = sz;
            if (nn == ss) {
                grow();
                ss = sz;
            }
            nn++;
            if (nn > stats.nmax) stats.nmax = nn;
            n = nn;

            int b = _back + 1;
            if (b >= int(ss)) b = 0;
            buf[b] = value;
            _back = b;
            stats.npush_back++;

            signal();
        }


    public:
        DQueue(size_t hint=200000) // was 32768
                : n(0)
                , sz(hint>2 ? hint : 2)
                , buf(new T[sz])
                , _front(sz/2)
                , _back(_front-1) {}

        virtual ~DQueue() {
            delete buf;
        }

        /// Insert value at front of queue
        void push_front(const T& value) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
            //sanity_check();

            size_t nn = n;
            size_t ss = sz;
            if (nn == ss) {
                grow();
                ss = sz;
            }
            nn++;
            if (nn > stats.nmax) stats.nmax = nn;
            n = nn;

            int f = _front - 1;
            if (f < 0) f = ss - 1;
            buf[f] = value;
            _front = f;
            stats.npush_front++;

            //sanity_check();
            signal();
            //broadcast();
        }

        /// Insert element at back of queue (default is just one copy)
        void push_back(const T& value, int ncopy=1) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
            //sanity_check();
            while (ncopy--)
                push_back_with_lock(value);
            //sanity_check();
            //broadcast();
        }

        template <typename opT>
        void scan(opT& op) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);

            int f = _front;
            size_t nn = n;
            int size = int(sz);
            std::cout << "IN Q " << nn << std::endl;

            while (nn--) {
                T* p = const_cast<T*>(buf + f);
                if (!op(p)) break;
                f++;
                if (f >= size) f = 0;
            }
        }

        /// Pop multiple values off the front of queue ... returns number popped ... might be zero

        /// r must refer to an array of dimension at least nmax ... you are presently
        /// given no more than max(size()/64,1) values ... arbitrary choice.
        ///
        /// multi-threaded tasks might cause fewer tasks to be taken
        int pop_front(int nmax, T* r, bool wait) {
            madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);

            size_t nn = n;

            if (nn==0 && wait) {
                while (n == 0) // !!! Must be n (memory) not nn (local copy)
                    CONDITION_VARIABLE_TYPE::wait();

                nn = n;
            }

            stats.npop_front++;
            if (nn) {
                size_t thesize = sz;
                //sanity_check();

                nmax = std::min(nmax,std::max(int(nn>>6),1));
                int retval; // Will return the number of items taken


                int f = _front;

                // Original loop was this
                //retval = nmax;
                //while (nmax--) {
                //    *r++ = buf[f++];
                //    if (f >= int(sz)) f = 0;
                //}

                // New loop includes checking for replicated multi-threaded task
                // ... take one task and then check that subsequent tasks differ
                nmax--;
                *r++ = buf[f++];
                if (f >= int(thesize)) f = 0;
                retval=1;
                while (nmax--) {
                    T ptr = buf[f];
                    if (ptr == *r) {
                        break;
                    }
                    else if (ptr) { // Null pointer indicates stolen task
                        *r++ = ptr;
                        f++;
                        if (f >= int(thesize)) f = 0;
                        retval++;
                    }
                }

                n = nn - retval;
                _front = f;

                //sanity_check();
                return retval;
            }
            else {
                return 0;
            }
        }

        /// Pop value off the front of queue
        std::pair<T,bool> pop_front(bool wait) {
            T r;
            int ngot = pop_front(1, &r, wait);
            return std::pair<T,bool>(r,ngot==1);
        }

        size_t size() const {
            return n;
        }

        bool empty() const {
            return n==0;
        }

        const DQStats& get_stats() const {
            return stats;
        }
    };

    class ThreadPool;           // Forward decl.
    class WorldTaskQueue;

    /// Simplified thread wrapper to hide pthread complexity

    /// If the thread is using any of the object state you cannot
    /// delete the object until the thread has terminated.
    ///
    /// The cleanest solution is to put the object on the heap and
    /// have the run method "delete this" at its end.
    class ThreadBase {
        friend class ThreadPool;
        static bool bind[3];
        static int cpulo[3];
        static int cpuhi[3];

        static void* main(void* self) {
#ifdef HAVE_PAPI
            begin_papi_measurement();
#endif

            try {
                ((ThreadBase*)(self))->run();
            }
            catch (const MPI::Exception& e) {
                //        print(e);
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
            return 0;
        }

        int pool_num; ///< Stores index of thread in pool or -1
        pthread_t id;

        void set_pool_thread_index(int i) {
            pool_num = i;
        }

    public:

        /// Default constructor ... must invoke \c start() to actually begin the thread.
        ThreadBase() : pool_num(-1) {};

        virtual ~ThreadBase() {};

        /// You implement this to do useful work
        virtual void run() = 0;

        /// Start the thread running
        void start() {
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

        /// A thread can call this to terminate its execution
        static void exit() {
            pthread_exit(0);
        }

        /// Get the pthread id of this thread (if running)
        const pthread_t& get_id() const {
            return id;
        }

        /// Get index of thread in ThreadPool (0,...,nthread-1) or -1 if not in ThreadPool
        int get_pool_thread_index() const {
            return pool_num;
        }

        /// Cancel this thread
        int cancel() const {
            return pthread_cancel(get_id());
        }


	/// Get no. of actual hardware processors
	static int num_hw_processors() {
#ifdef _SC_NPROCESSORS_CONF
	  int ncpu = sysconf(_SC_NPROCESSORS_CONF);
	  if (ncpu <= 0) MADNESS_EXCEPTION("ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)", ncpu);
#elif defined(HC_NCPU)
	  int mib[2]={CTL_HW,HW_NCPU}, ncpu;
	  size_t len = sizeof(ncpu);
	  if (sysctl(mib, 2, &ncpu, &len, NULL, 0) != 0)
	    MADNESS_EXCEPTION("ThreadBase: sysctl(CTL_HW,HW_NCPU) failed", 0);
	  std::cout << "NCPU " << ncpu << std::endl;
#else
	  int ncpu=1;
#endif
	  return ncpu;
	}

        /// Specify the affinity pattern or how to bind threads to cpus
        static void set_affinity_pattern(const bool bind[3], const int cpu[3]) {
            memcpy(ThreadBase::bind, bind, 3*sizeof(bool));
            memcpy(ThreadBase::cpulo, cpu, 3*sizeof(int));

	    int ncpu = num_hw_processors();

            // impose sanity and compute cpuhi
            for (int i=0; i<3; i++) {
                if (cpulo[i] < 0) cpulo[i] = 0;
                if (cpulo[i] >= ncpu) cpulo[i] = ncpu-1;

                if (i<2 && bind[i]) cpuhi[i] = cpulo[i];
                else cpuhi[i] = ncpu-1;

                //std::cout << "PATTERN " << i << " " << bind[i] << " " << cpulo[i] << " " << cpuhi[i] << std::endl;
            }
        }

        static void set_affinity(int logical_id, int ind=-1) {
            if (logical_id < 0 || logical_id > 2) {
                std::cout << "ThreadBase: set_affinity: logical_id bad?" << std::endl;
                return;
            }

            if (!bind[logical_id]) return;

            // If binding the main or rmi threads the cpu id is a specific cpu.
            //
            // If binding a pool thread, the cpuid is the lowest cpu
            // to be used.
            //
            // If floating any thread it is floated from the cpuid to ncpu-1

            int lo=cpulo[logical_id], hi=cpuhi[logical_id];

            if (logical_id == 2) {
                if (ind < 0) {
                    std::cout << "ThreadBase: set_affinity: pool thread index bad?" << std::endl;
                    return;
                }
                if (bind[2]) {
                    int nnn = hi-lo+1;
                    lo += (ind % nnn);
                    hi = lo;
                }
            }

#ifdef ON_A_MAC
#else
            cpu_set_t mask;
            CPU_ZERO(&mask);
            for (int i=lo; i<=hi; i++) CPU_SET(i,&mask);
            if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
                perror("system error message");
                std::cout << "ThreadBase: set_affinity: Could not set cpu Affinity" << std::endl;
            }
#endif
        }
    };


    /// Simplified thread wrapper to hide pthread complexity
    class Thread : public ThreadBase {
        void* (*f)(void *);
        void* args;

        void run() {
            f(args);
        }
    public:
        /// Default constructor ... must invoke \c start() to actually begin the thread.
        Thread() : f(0), args(0) {};

        /// Create a thread and start it running f(args)
        Thread(void* (*f)(void *), void* args=0)
                : f(f), args(args) {
            ThreadBase::start();
        }

        void start(void* (*f)(void *), void* args=0) {
            this->f = f;
            this->args = args;
            ThreadBase::start();
        }

        virtual ~Thread() {}
    };


    /// Contains attributes of a task

    /// \c generator : Setting this hints that a task will produce
    /// additional tasks and is used by the scheduler to
    /// increase/throttle parallelism. The default is false.
    ///
    /// \c stealable : Setting this indicates that a task may be
    /// migrated to another process for dynamic load balancing.  The
    /// default value is false.
    ///
    /// \c highpriority : indicates a high priority task. The default value is false.
    ///
    /// \c nthread : indicates number of threads. 0 threads is interpreted as 1 thread
    /// for backward compatibility and ease of specifying defaults. The default value
    /// is 0 (==1).
    class TaskAttributes {
        unsigned long flags;
    public:
    	static const unsigned long NTHREAD      = 0xff;          // Mask for nthread byte
        static const unsigned long GENERATOR    = 1ul<<8;        // Mask for generator bit
        static const unsigned long STEALABLE    = GENERATOR<<1;  // Mask for stealable bit
        static const unsigned long HIGHPRIORITY = GENERATOR<<2;  // Mask for priority bit

        TaskAttributes(unsigned long flags = 0) : flags(flags) {}

        TaskAttributes(const TaskAttributes& attr) : flags(attr.flags) {}

        bool is_generator() const {
            return flags&GENERATOR;
        }

        bool is_stealable() const {
            return flags&STEALABLE;
        }

        bool is_high_priority() const {
            return flags&HIGHPRIORITY;
        }

        void set_generator(bool generator_hint) {
            if (generator_hint) flags |= GENERATOR;
            else flags &= ~GENERATOR;
        }

        void set_stealable(bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
        }

        void set_highpriority(bool hipri) {
            if (hipri) flags |= HIGHPRIORITY;
            else flags &= ~HIGHPRIORITY;
        }

        /// Are you sure this is what you want to call?

        /// Only call this for a \c TaskAttributes that is \em not a base class
        /// of a task object.
        ///
        /// If you are trying to set the number of threads in an \em existing
        /// task you should call \c TaskInterface::set_nthread() instead.
        /// No doubt there is some virtual/protected/something voodoo to prevent
        /// you from doing harm.
        void set_nthread(int nthread) {
            MADNESS_ASSERT(nthread>=0 && nthread<256);
            flags = (flags & (~NTHREAD)) | (nthread & NTHREAD);
        }

        int get_nthread() const {
        	int n = flags & NTHREAD;
        	if (n == 0) n = 1;
        	return n;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & flags;
        }

        static TaskAttributes generator() {
            return TaskAttributes(GENERATOR);
        }

        static TaskAttributes hipri() {
            return TaskAttributes(HIGHPRIORITY);
        }

        static TaskAttributes multi_threaded(int nthread) {
            TaskAttributes t;
            t.set_nthread(nthread);
            return t;
        }
    };

    /// Used to pass info about thread environment into users task
    class TaskThreadEnv {
        const int _nthread; //< No. of threads collaborating on task
        const int _id;      //< Id of this thread (0,...,nthread-1)
        Barrier* _barrier;  //< Pointer to shared barrier, null if single thread

    public:
        TaskThreadEnv(int nthread, int id, Barrier* barrier)
            : _nthread(nthread), _id(id), _barrier(barrier)
        {}

        int nthread() const {return _nthread;}

        int id() const {return _id;}

        bool barrier() const {
            if (_nthread == 1)
                return true;
            else {
                MADNESS_ASSERT(_barrier);
                return _barrier->enter(_id);
            }
        }
    };


    /// Lowest level task interface

    /// The pool invokes run_multi_threaded() that does any necessary
    /// setup for multiple threads and then invokes the users \c run method.
    class PoolTaskInterface : public TaskAttributes {
    	friend class ThreadPool;

    private:
        Barrier* barrier;     //< Barrier, only allocated for multithreaded tasks
    	AtomicInt count;  //< Used to count threads as they start

    	/// Returns true for the one thread that should invoke the destructor
    	bool run_multi_threaded() {
            // As a thread enters this routine it increments the shared counter
            // to generate a unique id without needing any thread-local storage.
            // A downside is this does not preserve any relationships between thread
            // numbering and the architecture ... more work ahead.
            int nthread = get_nthread();
            if (nthread == 1) {
                run(TaskThreadEnv(1,0,0));
                return true;
            }
            else {
                int id = count++;
                volatile bool barrier_flag;
                barrier->register_thread(id, &barrier_flag);

                run(TaskThreadEnv(nthread, id, barrier));

                return barrier->enter(id);
            }
    	}

    public:
        PoolTaskInterface()
            : TaskAttributes()
            , barrier(0)
        {
            count = 0;
        }

        explicit PoolTaskInterface(const TaskAttributes& attr)
            : TaskAttributes(attr)
            , barrier(attr.get_nthread()>1 ? new Barrier(attr.get_nthread()) : 0)

        {
            count = 0;
        }

        /// Call this to reset the number of threads before the task is submitted

        /// Once a task has been constructed /c TaskAttributes::set_nthread()
        /// is insufficient because a multithreaded task includes a
        /// barrier that needs to know the number of threads.
        void set_nthread(int nthread) {
            if (nthread != get_nthread()) {
                TaskAttributes::set_nthread(nthread);
                delete barrier;
                if (nthread > 1)
                    barrier = new Barrier(nthread);
                else
                    barrier = 0;

            }
        }

        /// Override this method to implement a multi-threaded task

        /// \c info.nthread() will be the number of threads collaborating on this task
        ///
        /// \c info.id() will be the index of the current thread \c id=0,...,nthread-1
        ///
        /// \c info.barrier() will be a barrier for all of the threads, and returns
        /// \c true for the last thread to enter the barrier (other threads get false)
        virtual void run(const TaskThreadEnv& info) = 0;

        virtual ~PoolTaskInterface() {
            delete barrier;
        }
    };

    /// A no-op task used for various purposes
    class PoolTaskNull : public PoolTaskInterface {
    public:
        void run(const TaskThreadEnv& /*info*/) {}
        virtual ~PoolTaskNull() {}
    };


    /// A singleton pool of threads for dynamic execution of tasks.

    /// YOU MUST INSTANTIATE THE POOL WHILE RUNNING WITH JUST ONE THREAD
    class ThreadPool {
    private:
        friend class WorldTaskQueue;
        Thread *threads;              ///< Array of threads
        DQueue<PoolTaskInterface*> queue; ///< Queue of tasks
        int nthreads;		  ///< No. of threads
        volatile bool finish;              ///< Set to true when time to stop
        AtomicInt nfinished;

        static ThreadPool* instance_ptr;

        /// The constructor is private to enforce the singleton model
        ThreadPool(int nthread=-1) : nthreads(nthread), finish(false) {
            nfinished = 0;
            instance_ptr = this;
            if (nthreads < 0) nthreads = default_nthread();
            //std::cout << "POOL " << nthreads << std::endl;

            try {
                if (nthreads > 0)
                	threads = new Thread[nthreads];
                else
                	threads = 0;
            }
            catch (...) {
                MADNESS_EXCEPTION("memory allocation failed", 0);
            }

            for (int i=0; i<nthreads; i++) {
                threads[i].set_pool_thread_index(i);
                threads[i].start(pool_thread_main, (void *)(threads+i));
            }
        }

        ThreadPool(const ThreadPool&);           // Verboten
        void operator=(const ThreadPool&);       // Verboten

        /// Get number of threads from the environment
        int default_nthread() {
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
	      if (nthread < 2) nthread = 2;
	      nthread = nthread - 1; // One less than # physical processors
            }
            return nthread;
        }

        /// Run next task ... returns true if one was run ... blocks if wait is true
        bool run_task(bool wait) {
            if (!wait && queue.empty()) return false;
            std::pair<PoolTaskInterface*,bool> t = queue.pop_front(wait);
            // Task pointer might be zero due to stealing
            if (t.second && t.first) {
                PROFILE_BLOCK(working);
                if (t.first->run_multi_threaded())         // What we are here to do
                    delete t.first;
            }
            return t.second;
        }

        bool run_tasks(bool wait) {
            static const int nmax=128; // WAS 100 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
            PoolTaskInterface* taskbuf[nmax];
            int ntask = queue.pop_front(nmax, taskbuf, wait);
            for (int i=0; i<ntask; i++) {
                PROFILE_BLOCK(working);
                if (taskbuf[i]) { // Task pointer might be zero due to stealing
                    if (taskbuf[i]->run_multi_threaded()) {
                	delete taskbuf[i];
                    }
                }
            }
            return (ntask>0);
        }

        void thread_main(Thread* thread) {
            PROFILE_MEMBER_FUNC(ThreadPool);
            thread->set_affinity(2, thread->get_pool_thread_index());

#define MULTITASK
#ifdef  MULTITASK
            while (!finish) {
                run_tasks(true);
            }
#else
            while (!finish) {
                run_task(true);
            }
#endif
            nfinished++;
        }

        /// Forwards thread to bound member function
        static void* pool_thread_main(void *v) {
            instance()->thread_main((Thread*)(v));
            return 0;
        }


        /// Return a pointer to the only instance constructing as necessary
        static ThreadPool* instance(int nthread=-1) {
            if (!instance_ptr) instance_ptr = new ThreadPool(nthread);
            return instance_ptr;
        }


    public:
        /// Please invoke while in single threaded environment
        static void begin(int nthread=-1) {
            instance(nthread);
        }

        static void end() {
            if (!instance_ptr) return;
            instance()->finish = true;
            for (int i=0; i<instance()->nthreads; i++) {
                add(new PoolTaskNull);
            }
            while (instance_ptr->nfinished != instance_ptr->nthreads);
        }

        /// Add a new task to the pool
        static void add(PoolTaskInterface* task) {
            if (!task) MADNESS_EXCEPTION("ThreadPool: inserting a NULL task pointer", 1);
            int nthread = task->get_nthread();
            // Currently multithreaded tasks must be shoved on the end of the q
            // to avoid a race condition as multithreaded task is starting up
            if (task->is_high_priority() && nthread==1) {
                instance()->queue.push_front(task);
            }
            else {
                instance()->queue.push_back(task, nthread);
            }
        }

        template <typename opT>
        void scan(opT& op) {
            queue.scan(op);
        }


        /// Add a vector of tasks to the pool
        static void add(const std::vector<PoolTaskInterface*>& tasks) {
            typedef std::vector<PoolTaskInterface*>::const_iterator iteratorT;
            for (iteratorT it=tasks.begin(); it!=tasks.end(); ++it) {
                add(*it);
            }
        }

        /// An otherwise idle thread can all this to run a task

        /// Returns true if one was run
        static bool run_task() {
            return instance()->run_tasks(false);
        }

        /// Returns number of threads in the pool
        static size_t size() {
            return instance()->nthreads;
        }

        /// Returns number of tasks in the queue
        static size_t queue_size() {
            return instance()->queue.size();
        }

        /// Returns queue statistics
        static const DQStats& get_stats() {
            return instance()->queue.get_stats();
        }

        ~ThreadPool() {};
    };

}

#endif // MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
