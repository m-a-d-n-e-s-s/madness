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

  $Id $
*/
#ifndef MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
#define MADNESS_WORLD_WORLDTHREAD_H__INCLUDED

/// \file worldthread.h
/// \brief Implements Dqueue, Thread, ThreadBase and ThreadPool

#include <world/dqueue.h>
#include <vector>
#include <cstddef>
#include <cstdio>
#include <pthread.h>
#include <world/typestuff.h>
#include <typeinfo>

#ifdef MADNESS_TASK_PROFILING
#include <execinfo.h> // for backtrace_symbols
#include <cxxabi.h> // for abi::__cxa_demangle
#include <sstream> // for std::istringstream
#include <cstring> // for strchr & strrchr
#endif // MADNESS_TASK_PROFILING

#ifdef HAVE_INTEL_TBB
#include "tbb/tbb.h"
#endif


#ifndef _SC_NPROCESSORS_CONF
// Old macs don't have necessary support thru sysconf to determine the
// no. of processors so must use sysctl
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace madness {

    // Forward decl.
    class Barrier;
    class ThreadPool;
    class WorldTaskQueue;
    class AtomicInt;
    void error(const char *msg);

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
        static pthread_key_t thread_key; ///< Thread id key

        static void* main(void* self);

        int pool_num; ///< Stores index of thread in pool or -1
        pthread_t id;

        static void init_thread_key() {
           const int rc = pthread_key_create(&thread_key, NULL);
           if(rc != 0)
               MADNESS_EXCEPTION("pthread_key_create failed", rc);
        }

        static void delete_thread_key() {
           pthread_key_delete(thread_key);
        }

        void set_pool_thread_index(int i) { pool_num = i; }

    public:

        /// Default constructor ... must invoke \c start() to actually begin the thread.
        ThreadBase() : pool_num(-1) { }

        virtual ~ThreadBase() { }

        /// You implement this to do useful work
        virtual void run() = 0;

        /// Start the thread running
        void start();

        /// A thread can call this to terminate its execution
        static void exit() { pthread_exit(0); }

        /// Get the pthread id of this thread (if running)
        const pthread_t& get_id() const { return id; }

        /// Get index of thread in ThreadPool (0,...,nthread-1) or -1 if not in ThreadPool
        int get_pool_thread_index() const { return pool_num; }

        /// Cancel this thread
        int cancel() const { return pthread_cancel(get_id()); }


        /// Get no. of actual hardware processors
        static int num_hw_processors();

        /// Specify the affinity pattern or how to bind threads to cpus
        static void set_affinity_pattern(const bool bind[3], const int cpu[3]);

        static void set_affinity(int logical_id, int ind=-1);

        static ThreadBase* this_thread() {
            return static_cast<ThreadBase*>(pthread_getspecific(thread_key));
        }
    }; // class ThreadBase

    /// Simplified thread wrapper to hide pthread complexity
    class Thread : public ThreadBase {
        void* (*f)(void *);
        void* args;

        void run() { f(args); }

    public:
        /// Default constructor ... must invoke \c start() to actually begin the thread.
        Thread() : f(0), args(0) { }

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

        explicit TaskAttributes(unsigned long flags = 0) : flags(flags) {}

        TaskAttributes(const TaskAttributes& attr) : flags(attr.flags) {}

        virtual ~TaskAttributes() {}

        bool is_generator() const { return flags&GENERATOR; }

        bool is_stealable() const { return flags&STEALABLE; }

        bool is_high_priority() const { return flags&HIGHPRIORITY; }

        void set_generator(bool generator_hint) {
            if (generator_hint)
                flags |= GENERATOR;
            else
                flags &= ~GENERATOR;
        }

        void set_stealable(bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
        }

        void set_highpriority(bool hipri) {
            if (hipri)
                flags |= HIGHPRIORITY;
            else
                flags &= ~HIGHPRIORITY;
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
        	if (n == 0)
        	    n = 1;
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

#if HAVE_INTEL_TBB
        // I can not get the TaskThreadEnv to work with Barrier
        // Need to figure out why
        TaskThreadEnv(int nthread, int id)
            : _nthread(1), _id(id), _barrier(NULL)
        {};
#endif

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


#ifdef MADNESS_TASK_PROFILING

    namespace profiling {

        /// Task event object

        /// This object is used to record the task trace information, including
        /// submit, start, and stop times as well as identification information.
        class TaskEvent {
        private:
            std::pair<void*, unsigned long> id_; ///< Task identification information
            double times_[3]; ///< Task trace times: { submit, start, stop }
            unsigned long threads_; ///< Number of threads used by task

            /// Print demangled symbol name

            /// Add the demangled symbol name to \c os . If demangling fails,
            /// the unmodified symbol name is used instead. If symbol is NULL,
            /// "UNKNOWN" is used instead. A tab character is added after the
            /// symbol name.
            /// \param os The output stream
            /// \param symbol The symbol to add to the stream
            static void print_demangled(std::ostream& os, const char* symbol) {
                // Get the demagled symbol name
                if(symbol) {
                    int status = 0;
                    const char* name = abi::__cxa_demangle(symbol, 0, 0, &status);

                    // Append the demangled symbol name to the output stream
                    if(status == 0) {
                        os << name << "\t";
                        free((void*)name);
                    } else {
                        os << symbol << "\t";
                    }
                } else {
                    os << "UNKNOWN\t";
                }
            }

            /// Get name of the function pointer

            /// \return The mangled function name
            std::string get_name() const {

                // Get the backtrace symbol for the function address,
                // which contains the function name.
                void* const * func_ptr = const_cast<void* const *>(& id_.first);
                char** bt_sym = backtrace_symbols(func_ptr, 1);

                // Extract the mangled function name from the backtrace
                // symbol.
                std::string mangled_name;

#ifdef ON_A_MAC
                // Format of bt_sym is:
                // <frame #> <file name> <address> <mangled name> + <function offset>
                std::istringstream iss(bt_sym[0]);
                long frame;
                std::string file, address;
                iss >> frame >> file >> address >> mangled_name;
#else // Assume Linux
                // Format of bt_sym is:
                // <file>(<mangled name>+<function offset>) [<address>]
                const char* first = strchr(bt_sym[0],'(');
                if(first) {
                    ++first;
                    const char* last = strrchr(first,'+');
                    if(last)
                        mangled_name.assign(first, (last - first) - 1);
                }
#endif // ON_A_MAC

                // Free the backtrace buffer
                free(bt_sym);

                return mangled_name;
            }

        public:

            /// No constructors are defined.

            /// Record the start time of the task and collect task information

            /// \param id The task identifier (a function pointer or const char*)
            /// and an integer to differentiate the different types
            /// \param threads The number of threads this task uses
            /// \param submit_time The time that the task was submitted to the
            /// task queue
            void start(const std::pair<void*, unsigned int>& id,
                    const unsigned long threads, const double submit_time)
            {
                id_ = id;
                threads_ = threads;
                times_[0] = submit_time;
                times_[1] = wall_time();
            }

            /// Record the stop time of the task
            void stop() { times_[2] = wall_time(); }

            /// Output the task data using a tab seperated list

            /// Output information includes id pointer; the function, member
            /// function, object type name; number of threads used by the task;
            /// submit time; start time; and stop time.
            /// \param os The output stream
            /// \param te The task event to be output
            /// \return The \c os reference.
            friend std::ostream& operator<<(std::ostream& os, const TaskEvent& te) {
                // Add address to output stream
                os << std::hex << std::showbase << te.id_.first <<
                        std::dec << std::noshowbase << "\t";

                // Print the name
                switch(te.id_.second) {
                    case 1:
                        {
                            const std::string mangled_name = te.get_name();

                            // Print the demangled name
                            if(! mangled_name.empty())
                                print_demangled(os, mangled_name.c_str());
                            else
                                os << "UNKNOWN\t";
                        }
                        break;
                    case 2:
                        print_demangled(os, static_cast<const char*>(te.id_.first));
                        break;
                    default:
                        os << "UNKNOWN\t";
                }

                // Print:
                // # of threads, submit time, start time, stop time
                os << te.threads_;
                const std::streamsize precision = os.precision();
                os.precision(6);
                os << std::fixed << "\t" << te.times_[0]
                        << "\t" << te.times_[1] << "\t" << te.times_[2];
                os.precision(precision);
                return os;
            }

        }; // class TaskEvent

        /// Task event list base class

        /// This base class allows the data to be stored in a linked list
        class TaskEventListBase {
        private:
            TaskEventListBase* next_;

            TaskEventListBase(const TaskEventListBase&);
            TaskEventListBase& operator=(const TaskEventListBase&);

        public:

            /// default constructor
            TaskEventListBase() : next_(NULL) { }

            /// virtual destructor
            virtual ~TaskEventListBase() { }

            /// Get the next event list in the linked list
            TaskEventListBase* next() const { return next_; }

            /// Insert \c list after this list

            /// \param list The list to be inserted
            void insert(TaskEventListBase* list) {
                if(next_)
                    list->next_ = next_;
                next_ = list;
            }

            /// Output task event list to an output stream

            /// \param os The ouptut stream
            /// \param tel The task event list to be output
            /// \return The modified output stream
            friend inline std::ostream& operator<<(std::ostream& os, const TaskEventListBase& tel) {
                return tel.print_events(os);
            }

        private:

            /// print the events
            virtual std::ostream& print_events(std::ostream&) const = 0;

        }; // class TaskEventList

        /// A list of task events

        /// \tparam N The maximum number of events held by the list
        /// This object is used by the thread pool to record task data
        template <std::size_t N>
        class TaskEventList : public TaskEventListBase {
        private:
            std::size_t n_; ///< The number of events recorded
            TaskEvent events_[N]; ///< The event list

            // Not allowed
            TaskEventList(const TaskEventList&);
            TaskEventList& operator=(const TaskEventList&);

        public:

            /// Default constructor
            TaskEventList() : TaskEventListBase(), n_(0ul) { }

            /// virtual destructor
            virtual ~TaskEventList() { }

            /// Get a new event from this list

            /// \warning This function can only be called \c N times. It is the
            /// callers resonsibility to ensure that it is not called too many
            /// times.
            TaskEvent* event() { return events_ + (n_++); }

        private:

            /// Print events recorded in this list.
            virtual std::ostream& print_events(std::ostream& os) const {
                const int thread_id = ThreadBase::this_thread()->get_pool_thread_index();
                for(std::size_t i = 0; i < n_; ++i)
                    os << thread_id << "\t" << events_[i] << std::endl;
                return os;
            }

        }; // class TaskEventList

        /// This object collects and prints task profiling data

        /// \note Each thread has its own \c TaskProfiler object, so only one
        /// thread will ever operate on this object at a time and all operations
        /// are inheirently thread safe.
        class TaskProfiler {
        private:
            TaskEventListBase* head_; ///< The head of the linked list of data
            TaskEventListBase* tail_; ///< The tail of the linked list of data

            static Mutex output_mutex_; ///< Mutex used to lock the output file

            // Not allowed
            TaskProfiler(const TaskProfiler&);
            TaskProfiler& operator=(const TaskProfiler&);

        public:
            static const char* output_file_name_; ///< The output file name
                    ///< This variable is initialized by \c ThreadPool::begin .
                    ///< This variable is assigned the value given by the
                    ///< environment variable MAD_TASKPROFILER_NAME.
        public:
            /// Default constructor
            TaskProfiler() : head_(NULL), tail_(NULL) { }

            /// TaskProfiler destructor
            ~TaskProfiler() {
                // Cleanup linked list
                TaskEventListBase* next = NULL;
                while(head_ != NULL) {
                    next = head_->next();
                    delete head_;
                    head_ = next;
                }
            }

            /// Create a new task event list

            /// \tparam N The maximum number of elements that the list can contain
            /// \return A new task event list
            template <std::size_t N>
            TaskEventList<N>* new_list() {
                // Create a new event list
                TaskEventList<N>* list = new TaskEventList<N>();

                // Append the list to the tail of the linked list
                if(head_ != NULL) {
                    tail_->insert(list);
                    tail_ = list;
                } else {
                    head_ = list;
                    tail_ = list;
                }
                return list;
            }

            /// Write the profile data to file

            /// The data is cleared after it is written to the file, so this
            /// function may be called more than once.
            /// \warning This function should only be called from the thread
            /// that owns it, otherwise data will likely be corrupted.
            /// \note This function is thread safe, in that it may be called by
            /// different objects in different threads simultaniously.
            void write_to_file();
        }; // class TaskProfiler

    } // namespace profiling

#endif // MADNESS_TASK_PROFILING


    /// Lowest level task interface

    /// The pool invokes run_multi_threaded() that does any necessary
    /// setup for multiple threads and then invokes the users \c run method.
    class PoolTaskInterface :
#if HAVE_INTEL_TBB
            public tbb::task,
#endif // HAVE_INTEL_TBB
            public TaskAttributes
    {
        friend class ThreadPool;

    private:
        Barrier* barrier;     //< Barrier, only allocated for multithreaded tasks
    	AtomicInt count;  //< Used to count threads as they start

#ifdef MADNESS_TASK_PROFILING
    	profiling::TaskEvent* task_event_;
    	double submit_time_;
        std::pair<void*, unsigned long> id_;

        void set_event(profiling::TaskEvent* task_event) {
            task_event_ = task_event;
        }

        /// Collect info on the task and record the sbumit time.
        void submit() {
            submit_time_ = wall_time();
            this->get_id(id_);
        }
#endif // MADNESS_TASK_PROFILING

        /// Object that is used to convert function and member function pointers into void*

        /// \note This is technically not supported by the C++ standard but
        /// it will likely not cause any issues here (Famous last words?).
        template <typename T>
        union FunctionPointerGrabber {
            T in;
            void* out;
        };

    protected:


        template <typename fnT>
        static typename enable_if_c<detail::function_traits<fnT>::value ||
                detail::memfunc_traits<fnT>::value>::type
        make_id(std::pair<void*,unsigned long>& id, fnT fn) {
            FunctionPointerGrabber<fnT> poop;
            poop.in = fn;
            id.first = poop.out;
            id.second = 1ul;
        }

        template <typename fnobjT>
        static typename disable_if_c<detail::function_traits<fnobjT>::value ||
                detail::memfunc_traits<fnobjT>::value>::type
        make_id(std::pair<void*,unsigned long>& id, const fnobjT&) {
            id.first = reinterpret_cast<void*>(const_cast<char*>(typeid(fnobjT).name()));
            id.second = 2ul;
        }

    private:

        virtual void get_id(std::pair<void*,unsigned long>& id) const {
            id.first = NULL;
            id.second = 0ul;
        }

    	/// Returns true for the one thread that should invoke the destructor
    	bool run_multi_threaded() {
#if HAVE_INTEL_TBB
            MADNESS_EXCEPTION("run_multi_threaded should not be called when using Intel TBB", 1);
#else
            // As a thread enters this routine it increments the shared counter
            // to generate a unique id without needing any thread-local storage.
            // A downside is this does not preserve any relationships between thread
            // numbering and the architecture ... more work ahead.
            int nthread = get_nthread();
            if (nthread == 1) {
#ifdef MADNESS_TASK_PROFILING
                task_event_->start(id_, nthread, submit_time_);
#endif // MADNESS_TASK_PROFILING
                run(TaskThreadEnv(1,0,0));
#ifdef MADNESS_TASK_PROFILING
                task_event_->stop();
#endif // MADNESS_TASK_PROFILING
                return true;
            }
            else {
                int id = count++;
                volatile bool barrier_flag;
                barrier->register_thread(id, &barrier_flag);

#ifdef MADNESS_TASK_PROFILING
                if(id == 0)
                    task_event_->start(id_, nthread, submit_time_);
#endif // MADNESS_TASK_PROFILING

                run(TaskThreadEnv(nthread, id, barrier));

#ifdef MADNESS_TASK_PROFILING
                const bool cleanup = barrier->enter(id);
                if(cleanup) task_event_->stop();
                return cleanup;
#else
                return barrier->enter(id);
#endif // MADNESS_TASK_PROFILING
            }
#endif // HAVE_INTEL_TBB
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

#if HAVE_INTEL_TBB
        tbb::task* execute() {
            int nthread = get_nthread();
            int id = count++;
//            volatile bool barrier_flag;
//            barrier->register_thread(id, &barrier_flag);

            run( TaskThreadEnv(nthread, id) );
            return NULL;
        }
#endif

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
        virtual void get_id(std::pair<void*,unsigned long>& id) const {
            PoolTaskInterface::make_id(id, &PoolTaskNull::run);
        }
    };

    /// ThreadPool thread object

    /// This object holds thread local data for thread pool threads. It can be
    /// accessed via \c ThreadBase::this_thread() .
    class ThreadPoolThread : public Thread {
    private:
        // Thread local data for thread pool
#ifdef MADNESS_TASK_PROFILING
        profiling::TaskProfiler profiler_;
#endif // MADNESS_TASK_PROFILING

    public:

        /// Default constructor
        ThreadPoolThread() : Thread() { }

        /// Virtual destructor
        virtual ~ThreadPoolThread() { }

#ifdef MADNESS_TASK_PROFILING
        /// Task profiler accessor
        profiling::TaskProfiler& profiler() { return profiler_; }
#endif // MADNESS_TASK_PROFILING
    };

    /// A singleton pool of threads for dynamic execution of tasks.

    /// YOU MUST INSTANTIATE THE POOL WHILE RUNNING WITH JUST ONE THREAD
    class ThreadPool {
    private:
        friend class WorldTaskQueue;

        // Thread pool data
        ThreadPoolThread *threads; ///< Array of threads
        ThreadPoolThread main_thread; ///< Placeholder for main thread tls
        DQueue<PoolTaskInterface*> queue; ///< Queue of tasks
        int nthreads; ///< No. of threads
        volatile bool finish; ///< Set to true when time to stop
        AtomicInt nfinished; ///< Thread pool exit counter

        // Static data
        static ThreadPool* instance_ptr; ///< Singleton pointer
#ifdef __bgq__
#warning WE NEED TO TUNE THE nmax PARAMETER
#endif
        static const int nmax=128; // WAS 100 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG

        /// The constructor is private to enforce the singleton model
        ThreadPool(int nthread=-1);

        ThreadPool(const ThreadPool&);           // Verboten
        void operator=(const ThreadPool&);       // Verboten

        /// Get number of threads from the environment
        int default_nthread();

        /// Run next task ... returns true if one was run ... blocks if wait is true
        bool run_task(bool wait, ThreadPoolThread* this_thread) {
#if HAVE_INTEL_TBB
            MADNESS_EXCEPTION("run_task should not be called when using Intel TBB", 1);
#else

#ifdef MADNESS_TASK_PROFILING
            profiling::TaskEventList<1>* event_list =
                    this_thread->profiler().new_list<1>();
#endif // MADNESS_TASK_PROFILING
            if (!wait && queue.empty()) return false;
            std::pair<PoolTaskInterface*,bool> t = queue.pop_front(wait);
            // Task pointer might be zero due to stealing
            if (t.second && t.first) {
#ifdef MADNESS_TASK_PROFILING
                t.first->set_event(event_list->event());
#endif // MADNESS_TASK_PROFILING
                if (t.first->run_multi_threaded())         // What we are here to do
                    delete t.first;
            }
            return t.second;
#endif
        }

        bool run_tasks(bool wait, ThreadPoolThread* this_thread) {
#if HAVE_INTEL_TBB
//            if (!wait && tbb_task_list->empty()) return false;
//            tbb::task* t = &tbb_task_list->pop_front();
//            if (t) {
//                tbb_parent_task->increment_ref_count();
//                tbb_parent_task->enqueue(*t);
//            }

//            wait = (tbb_parent_task->ref_count() >= 1) ? false : true;
//            return wait;

            MADNESS_EXCEPTION("run_tasks should not be called when using Intel TBB", 1);
#else

#ifdef MADNESS_TASK_PROFILING
            profiling::TaskEventList<nmax>* event_list =
                    this_thread->profiler().new_list<nmax>();
#endif // MADNESS_TASK_PROFILING
            PoolTaskInterface* taskbuf[nmax];
            int ntask = queue.pop_front(nmax, taskbuf, wait);
            for (int i=0; i<ntask; ++i) {
                if (taskbuf[i]) { // Task pointer might be zero due to stealing
#ifdef MADNESS_TASK_PROFILING
                    taskbuf[i]->set_event(event_list->event());
#endif // MADNESS_TASK_PROFILING
                    if (taskbuf[i]->run_multi_threaded()) {
                        delete taskbuf[i];
                    }
                }
            }
            return (ntask>0);
#endif
        }

        void thread_main(ThreadPoolThread* const thread);

        /// Forwards thread to bound member function
        static void* pool_thread_main(void *v);


        /// Return a pointer to the only instance constructing as necessary
        static ThreadPool* instance(int nthread=-1) {
            if (!instance_ptr) instance_ptr = new ThreadPool(nthread);
            return instance_ptr;
        }


    public:

#if HAVE_INTEL_TBB
        // all tasks run as children of tbb_parent_task
        // be sure to allocate tasks with tbb_parent_task->allocate_child()
        static tbb::empty_task* tbb_parent_task;
        static tbb::task_scheduler_init* tbb_scheduler;
#endif

        /// Please invoke while in single threaded environment
        static void begin(int nthread=-1);

        static void end();

        /// Add a new task to the pool
        static void add(PoolTaskInterface* task) {
#ifdef MADNESS_TASK_PROFILING
            task->submit();
#endif // MADNESS_TASK_PROFILING
#if HAVE_INTEL_TBB
            ThreadPool::tbb_parent_task->increment_ref_count();
            ThreadPool::tbb_parent_task->spawn(*t);
#else
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
#endif
        }

        template <typename opT>
        void scan(opT& op) {
            queue.scan(op);
        }


        /// Add a vector of tasks to the pool
        static void add(const std::vector<PoolTaskInterface*>& tasks) {
#if HAVE_INTEL_TBB
            MADNESS_EXCEPTION("Do not add tasks to the madness task queue when using Intel TBB.", 1);
#else
            typedef std::vector<PoolTaskInterface*>::const_iterator iteratorT;
            for (iteratorT it=tasks.begin(); it!=tasks.end(); ++it) {
                add(*it);
            }
#endif
        }

        /// An otherwise idle thread can all this to run a task

        /// Returns true if one was run
        static bool run_task() {
#ifdef MADNESS_TASK_PROFILING
            return instance()->run_tasks(false, static_cast<ThreadPoolThread*>(ThreadBase::this_thread()));
#else
            return instance()->run_tasks(false, NULL);
#endif // MADNESS_TASK_PROFILING
        }

        /// Returns number of threads in the pool
        static std::size_t size() {
            return instance()->nthreads;
        }

        /// Returns number of tasks in the queue
        static std::size_t queue_size() {
            return instance()->queue.size();
        }

        /// Returns queue statistics
        static const DQStats& get_stats();

        ~ThreadPool() {
#if HAVE_INTEL_TBB
            delete(tbb_parent_task);
            tbb_scheduler->terminate();
            delete(tbb_scheduler);
#endif
        }
    };

}

#endif // MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
