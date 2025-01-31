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

#ifndef MADNESS_WORLD_THREAD_H__INCLUDED
#define MADNESS_WORLD_THREAD_H__INCLUDED

/**
 \file thread.h
 \brief Implements Dqueue, Thread, ThreadBase and ThreadPool.
 \ingroup threads
*/

#include <madness/world/thread_info.h>
#include <madness/world/dqueue.h>
#include <madness/world/function_traits.h>
#include <vector>
#include <cstddef>
#include <cstdio>
#include <pthread.h>
#include <type_traits>
#include <typeinfo>
#include <new>

//////////// Parsec Related Begin ////////////////////
#ifdef HAVE_PARSEC
#include "parsec.h"
#endif
//////////// Parsec Related End ////////////////////

#ifdef MADNESS_TASK_PROFILING
#include <execinfo.h> // for backtrace_symbols
#ifndef USE_LIBIBERTY
#include <cxxabi.h> // for abi::__cxa_demangle
#else
extern "C" {
  extern char * cplus_demangle (const char *mangled, int options);
#define DMGL_NO_OPTS     0              /* For readability... */
}
#endif
#include <sstream> // for std::istringstream
#include <cstring> // for strchr & strrchr
#endif // MADNESS_TASK_PROFILING

#ifdef HAVE_INTEL_TBB
#include <tbb/task_arena.h>
#ifndef TBB_PREVIEW_GLOBAL_CONTROL
# define TBB_PREVIEW_GLOBAL_CONTROL 1
#endif
# include <tbb/global_control.h>
#endif


#ifndef _SC_NPROCESSORS_CONF
// Old macs don't have necessary support thru sysconf to determine the
// no. of processors so must use sysctl
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace madness {

    // Forward decls.
    class Barrier;
    class ThreadPool;
    class WorldTaskQueue;
    class AtomicInt;
    void error(const char *msg);

    class ThreadBinder {
      static const size_t maxncpu = 1024;
      bool print;
      size_t ncpu = 0;
      bool do_bind = true;
      size_t cpus[maxncpu];
      std::atomic<size_t> nextcpu = 0;
      static thread_local bool bound;
      
    public:
      
      ThreadBinder(bool print = false) : print(print) {
#ifndef ON_A_MAC
	ncpu = 0;
        cpu_set_t mask;
        sched_getaffinity(0, sizeof(mask), &mask);
        for (size_t i=0; i<maxncpu; i++) {
	  if (CPU_ISSET(int(i),&mask)) {
	    MADNESS_CHECK(ncpu <= maxncpu);
	    cpus[ncpu++] = i;
	  }
	}
	if (print) {
	  std::cout << "ncpu: " << get_ncpu() << std::endl;
	  for (size_t i=0; i<get_ncpu(); i++) {
	    std::cout << get_cpus()[i] << " " ;
	  }
	  std::cout << std::endl;
	}
	nextcpu = ncpu/2;
#endif
          if (this->print) { };
      }

      void set_do_bind(bool value) {do_bind = value;}
      
      const size_t* get_cpus() const { return cpus; }
      
      const size_t get_ncpu() const { return ncpu; }
      
      void bind() {
#ifndef ON_A_MAC
	if (do_bind && !bound) { // In TBB this is called by each task, so check do_bind first
	  bound = true;
	  cpu_set_t mask;
	  CPU_ZERO(&mask);
	  size_t cpu = cpus[nextcpu++ % ncpu];
	  CPU_SET(cpu, &mask);
	  sched_setaffinity(0, sizeof(mask), &mask);
	  if (print) std::cout << "bound thread to " << cpu << std::endl;
	}
#endif
      }
    };

    extern ThreadBinder binder;

    /// \addtogroup threads
    /// @{

    /// Simplified thread wrapper to hide pthread complexity.

    /// If the thread is using any of the object state, you cannot
    /// delete the object until the thread has terminated.
    ///
    /// The cleanest solution is to put the object on the heap and
    /// have the run method `delete this` at its end.
    class ThreadBase {
        friend class ThreadPool;

        static pthread_key_t thread_key; ///< Thread id key.

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in,out] self Description needed.
        /// \return Description needed.
        static void* main(void* self);

        int pool_num; ///< Stores index of thread in pool or -1.
        pthread_t id; ///< \todo Brief description needed.

        /// \todo Brief description needed.
        static void init_thread_key() {
           const int rc = pthread_key_create(&thread_key, nullptr);
           if(rc != 0)
               MADNESS_EXCEPTION("pthread_key_create failed", rc);
        }

        /// \todo Brief description needed.
        static void delete_thread_key() {
           pthread_key_delete(thread_key);
        }

        /// Sets the index of this thread within the pool.

        /// \todo Verify documentation.
        /// \param[in] i The index of this thread.
        void set_pool_thread_index(int i) {
            pool_num = i;
        }

#if defined(HAVE_IBMBGQ) and defined(HPM)
        static const int hpm_thread_id_all = -10; ///< \todo Brief description needed.
        static const int hpm_thread_id_main = -2; ///< \todo Brief description needed.
        static bool main_instrumented; ///< \todo Brief description needed.
        static bool all_instrumented; ///< \todo Brief description needed.
        static int hpm_thread_id; ///< \todo Brief description needed.
#endif

    public:

        /// Default constructor.

        /// Sets up the thread; however, \c start() must be invoked to
        /// actually begin the thread.
        ThreadBase() : pool_num(-1) { }

        virtual ~ThreadBase() { }

        /// Function to be executed by the thread.

        /// Override this to do work.
        virtual void run() = 0;

        /// Start the thread running.
        void start();

        /// A thread can call this to terminate its execution.
        static void exit() {
            pthread_exit(0);
        }

        /// Get the pthread id of this thread (if running).
        const pthread_t& get_id() const {
            return id;
        }

        /// Get index of this thread in \c ThreadPool.

        /// \return (0,...,nthread-1) or -1 if not in the \c ThreadPool.
        int get_pool_thread_index() const {
            return pool_num;
        }

        /// Cancel this thread.
        int cancel() const {
            return pthread_cancel(get_id());
        }


        /// Get number of actual hardware processors.

        /// \return The number of hardward processors.
        static int num_hw_processors();

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
        static ThreadBase* this_thread() {
            return static_cast<ThreadBase*>(pthread_getspecific(thread_key));
        }

#if defined(HAVE_IBMBGQ) and defined(HPM)
        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in] hpm_thread_id Description needed.
        static void set_hpm_thread_env(int hpm_thread_id);
#endif
    }; // class ThreadBase

    /// Simplified thread wrapper to hide pthread complexity.
    class Thread : public ThreadBase {
        void* (*f)(void *); ///< The function called for executing this thread. \todo should we replace this by a std::function?
        void* args; ///< The arguments passed to this thread for execution.

        /// Invokes the function for this thread.
        void run() {
            f(args);
        }

    public:
        /// Default constructor.

        /// \c start() must be invoked to actually execute the thread.
        Thread() : f(nullptr), args(nullptr) { }

        /// Create a thread and start it running `f(args)`.

        /// \param[in] f The function to be called.
        /// \param[in,out] args The arguments to the function.
        Thread(void* (*f)(void *), void* args=nullptr)
                : f(f), args(args) {
            ThreadBase::start();
        }

        /// Start the thread by running `f(args)`.

        /// \param[in] f The function to be called.
        /// \param[in,out] args The arguments to the function.
        void start(void* (*f)(void *), void* args=nullptr) {
            this->f = f;
            this->args = args;
            ThreadBase::start();
        }

        virtual ~Thread() = default;
    }; // class Thread


    /// Contains attributes of a task.

    /// The current attributes are:
    /// - \c generator : Setting this hints that a task will produce
    ///   additional tasks and is used by the scheduler to
    ///   increase/throttle parallelism. The default is false.
    /// - \c stealable : Setting this indicates that a task may be
    ///   migrated to another process for dynamic load balancing. The
    ///   default value is false.
    /// - \c highpriority : indicates a high priority task. The default
    ///   value is false.
    /// - \c nthread : indicates number of threads. 0 threads is interpreted
    ///   as 1 thread for backward compatibility and ease of specifying
    ///   defaults. The default value is 0 (==1).
    class TaskAttributes {
        unsigned long flags; ///< Byte-string storing the specified attributes.

    public:
    	static const unsigned long NTHREAD = 0xff; ///< Mask for nthread byte.
        static const unsigned long GENERATOR = 1ul<<8; ///< Mask for generator bit.
        static const unsigned long STEALABLE = GENERATOR<<1; ///< Mask for stealable bit.
        static const unsigned long HIGHPRIORITY = GENERATOR<<2; ///< Mask for priority bit.

        /// Sets the attributes to the desired values.

        /// `flags`, if unspecified sets all attributes to their default
        /// values.
        /// \param[in] flags The attribute values.
        explicit TaskAttributes(unsigned long flags = 0)
            : flags(flags) {}

        /// Copy constructor.

        /// \param[in] attr The attributes to copy.
        TaskAttributes(const TaskAttributes& attr)
            : flags(attr.flags) {}

        virtual ~TaskAttributes() {}

        /// Test if the generator attribute is true.

        /// \return True if this task is a generator, false otherwise.
        bool is_generator() const {
            return flags&GENERATOR;
        }

        /// Test if the stealable attribute is true.

        /// \return True if this task is stealable, false otherwise.
        bool is_stealable() const {
            return flags&STEALABLE;
        }

        /// Test if the high priority attribute is true.

        /// \return True if this task is a high priority, false otherwise.
        bool is_high_priority() const {
            return flags&HIGHPRIORITY;
        }

        /// Sets the generator attribute.

        /// \param[in] generator_hint The new value for the generator attribute.
        TaskAttributes& set_generator(bool generator_hint) {
            if (generator_hint)
                flags |= GENERATOR;
            else
                flags &= ~GENERATOR;
            return *this;
        }

        /// Sets the stealable attribute.

        /// \param[in] stealable The new value for the stealable attribute.
        TaskAttributes& set_stealable(bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
            return *this;
        }

        /// Sets the high priority attribute.

        /// \param[in] hipri The new value for the high priority attribute.
        TaskAttributes& set_highpriority(bool hipri) {
            if (hipri)
                flags |= HIGHPRIORITY;
            else
                flags &= ~HIGHPRIORITY;
            return *this;
        }

        /// Set the number of threads.

        /// \attention Are you sure this is what you want to call? Only call
        /// this for a \c TaskAttributes that is \em not a base class of a task
        /// object.
        /// \p If you are trying to set the number of threads in an \em existing
        /// task you should call \c TaskInterface::set_nthread() instead. No
        /// doubt there is some virtual/protected/something voodoo to prevent
        /// you from doing harm.
        ///
        /// \todo Perhaps investigate a way to make this function only accessible
        /// from the intended functions (using the so-called voodoo)?
        ///
        /// \param[in] nthread The new number of threads.
        void set_nthread(int nthread) {
            MADNESS_ASSERT(nthread>=0 && nthread<256);
            flags = (flags & (~NTHREAD)) | (nthread & NTHREAD);
        }

        /// Get the number of threads.

        /// \return The number of threads.
        int get_nthread() const {
        	int n = flags & NTHREAD;
        	if (n == 0)
        	    n = 1;
        	return n;
        }

        /// Serializes the attributes for I/O.

        /// tparam Archive The archive type.
        /// \param[in,out] ar The archive.
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & flags;
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
        static TaskAttributes generator() {
            return TaskAttributes(GENERATOR);
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
        static TaskAttributes hipri() {
            return TaskAttributes(HIGHPRIORITY);
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
        static TaskAttributes multi_threaded(int nthread) {
            TaskAttributes t;
            t.set_nthread(nthread);
            return t;
        }
    };

    /// Used to pass information about the thread environment to a user's task.
    class TaskThreadEnv {
        const int _nthread; ///< Number of threads collaborating on task.
        const int _id; ///< ID of this thread (0,...,nthread-1).
        Barrier* _barrier; ///< Pointer to the shared barrier, `null` if there is only a single thread.

    public:
        /// Constructor collecting necessary environmental information.

        /// \todo Verify this documentation.
        /// \param[in] nthread The number of threads collaborating on this task.
        /// \param[in] id The ID of this thread.
        /// \param[in] barrier Pointer to the shared barrier.
        TaskThreadEnv(int nthread, int id, Barrier* barrier)
            : _nthread(nthread), _id(id), _barrier(barrier)
        {::madness::binder.bind();}

#if HAVE_INTEL_TBB
        /// Constructor collecting necessary environmental information.

        /// \todo Verify this documentation.
        /// \param[in] nthread The number of threads collaborating on this task.
        /// \param[in] id The ID of this thread.
        ///
        /// \todo I cannot get the TaskThreadEnv to work with Barrier.
        /// Need to figure out why.
        TaskThreadEnv(int nthread, int id)
            : _nthread(nthread), _id(id), _barrier(nullptr)
      {::madness::binder.bind();};
#endif

        /// Get the number of threads collaborating on this task.

        /// \return The number of threads.
        int nthread() const {
            return _nthread;
        }

        /// Get the ID of this thread.

        /// \return The ID of this thread.
        int id() const {
            return _id;
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
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

        /// Task event class.

        /// This class is used to record the task trace information, including
        /// submit, start, and stop times, as well as identification information.
        class TaskEvent {
        private:
            double times_[3]; ///< Task trace times: { submit, start, stop }.
            std::pair<void*, unsigned short> id_; ///< Task identification information.
            unsigned short threads_; ///< Number of threads used by the task.

            /// Print demangled symbol name.

            /// Add the demangled symbol name to \c os. If demangling fails,
            /// the unmodified symbol name is used instead. If symbol is NULL,
            /// "UNKNOWN" is used instead. A tab character is added after the
            /// symbol name.
            /// \param[in,out] os The output stream.
            /// \param[in] symbol The symbol to add to the stream.
            static void print_demangled(std::ostream& os, const char* symbol) {
                // Get the demagled symbol name
                if(symbol) {
                    int status = 0;
#ifndef USE_LIBIBERTY
                    const char* name = abi::__cxa_demangle(symbol, 0, 0, &status);
#else
		    char* name = cplus_demangle(symbol, DMGL_NO_OPTS);
#endif
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

            /// Get name of the function pointer.

            /// \return The mangled function name.
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

            // Only default constructors are needed.

            /// Record the start time of the task and collect task information.

            /// \param[in,out] id The task identifier (a function pointer or const char*)
            ///     and an integer to differentiate the different types.
            /// \param[in] threads The number of threads this task uses.
            /// \param[in] submit_time The time that the task was submitted to the
            ///     task queue.
            void start(const std::pair<void*, unsigned short>& id,
                    const unsigned short threads, const double submit_time)
            {
                id_ = id;
                threads_ = threads;
                times_[0] = submit_time;
                times_[1] = wall_time();
            }

            /// Record the stop time of the task.
            void stop() {
                times_[2] = wall_time();
            }

            /// Output the task data using a tab-separated list.

            /// Output information includes
            /// - the ID pointer
            /// - the function, member function, and object type name
            /// - the number of threads used by the task
            /// - the submit time
            /// - the start time
            /// - the stop time.
            ///
            /// \param[in,out] os The output stream.
            /// \param[in] te The task event to be output.
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

        /// Task event list base class.

        /// This base class allows the data to be stored in a linked list.
        class TaskEventListBase {
        private:
            TaskEventListBase* next_; ///< The next task event in the list.

            TaskEventListBase(const TaskEventListBase&) = delete;
            TaskEventListBase& operator=(const TaskEventListBase&) = delete;

        public:

            /// Default constructor.
            TaskEventListBase()
                : next_(nullptr) { }

            /// Virtual destructor.
            virtual ~TaskEventListBase() = default;

            /// Get the next event list in the linked list.

            /// \return The next event list.
            TaskEventListBase* next() const {
                return next_;
            }

            /// Insert \c list after this list.

            /// \param[in] list The list to be inserted.
            void insert(TaskEventListBase* list) {
                if(next_)
                    list->next_ = next_;
                next_ = list;
            }

            /// Output a task event list to an output stream.

            /// \param[in,out] os The ouptut stream.
            /// \param[in] tel The task event list to be output.
            /// \return The modified output stream.
            friend inline std::ostream& operator<<(std::ostream& os, const TaskEventListBase& tel) {
                return tel.print_events(os);
            }

        private:

            /// Print the events.
            virtual std::ostream& print_events(std::ostream&) const = 0;

        }; // class TaskEventList

        /// A list of task events.

        /// This object is used by the thread pool to record task data.
        class TaskEventList : public TaskEventListBase {
        private:
            unsigned int n_; ///< The number of events recorded.
            std::unique_ptr<TaskEvent[]> events_; ///< The event array.

            TaskEventList(const TaskEventList&) = delete;
            TaskEventList& operator=(const TaskEventList&) = delete;

        public:

            /// Default constructor.

            /// \param[in] nmax The maximum number of task events.
            /// \todo Should nmax be stored? I think it used to be a template
            ///    parameter (N), which is no longer present.
            TaskEventList(const unsigned int nmax) :
                TaskEventListBase(), n_(0ul), events_(new TaskEvent[nmax])
            { }

            /// Virtual destructor.
            virtual ~TaskEventList() = default;

            /// Get a new event from this list.

            /// \warning This function can only be called \c nmax times. It is
            /// the caller's resonsibility to ensure that it is not called too
            /// many times.
            /// \return The new event from the list.
            TaskEvent* event() {
                return events_.get() + (n_++);
            }

        private:

            /// Print events recorded in this list.

            /// \param[in,out] os The output stream.
            /// \return The modified output stream.
            virtual std::ostream& print_events(std::ostream& os) const {
                const int thread_id = ThreadBase::this_thread()->get_pool_thread_index();
                for(std::size_t i = 0; i < n_; ++i)
                    os << thread_id << "\t" << events_[i] << std::endl;
                return os;
            }

        }; // class TaskEventList

        /// This class collects and prints task profiling data.

        /// \note Each thread has its own \c TaskProfiler object, so only one
        /// thread will ever operate on this object at a time and all operations
        /// are inheirently thread safe.
        class TaskProfiler {
        private:
            TaskEventListBase* head_; ///< The head of the linked list of data.
            TaskEventListBase* tail_; ///< The tail of the linked list of data.

            static Mutex output_mutex_; ///< Mutex used to lock the output file.

            TaskProfiler(const TaskProfiler&) = delete;
            TaskProfiler& operator=(const TaskProfiler&) = delete;

        public:
            /// The output file name.

            /// This variable is initialized by \c ThreadPool::begin and is
            /// assigned the value given by the environment variable
            /// `MAD_TASKPROFILER_NAME`.
            static const char* output_file_name_;

        public:
            /// Default constructor.
            TaskProfiler()
                : head_(nullptr), tail_(nullptr)
            { }

            /// Destructor.
            ~TaskProfiler() {
                // Cleanup linked list
                TaskEventListBase* next = nullptr;
                while(head_ != nullptr) {
                    next = head_->next();
                    delete head_;
                    head_ = next;
                }
            }

            /// Create a new task event list.

            /// \param[in] nmax The maximum number of elements that the list
            ///     can contain.
            /// \return A new task event list.
            TaskEventList* new_list(const std::size_t nmax) {
                // Create a new event list
                TaskEventList* list = new TaskEventList(nmax);

                // Append the list to the tail of the linked list
                if(head_ != nullptr) {
                    tail_->insert(list);
                    tail_ = list;
                } else {
                    head_ = list;
                    tail_ = list;
                }
                return list;
            }

            /// Write the profile data to file.

            /// The data is cleared after it is written to the file, so this
            /// function may be called more than once.
            ///
            /// \warning This function should only be called from the thread
            /// that owns it, otherwise data will likely be corrupted.
            ///
            /// \note This function is thread safe, in that it may be called by
            /// different objects in different threads simultaneously.
            void write_to_file();
        }; // class TaskProfiler

    } // namespace profiling

#endif // MADNESS_TASK_PROFILING


    /// Lowest level task interface.

    /// The pool invokes \c run_multi_threaded(), which does any necessary
    /// setup for multiple threads, and then invokes the user's \c run() method.
    class PoolTaskInterface : public TaskAttributes
    {
        friend class ThreadPool;

    private:

#ifdef MADNESS_TASK_PROFILING
    	profiling::TaskEvent* task_event_; ///< \todo Description needed.
    	double submit_time_; ///< \todo Description needed.
        std::pair<void*, unsigned short> id_; ///< \todo Description needed.

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in,out] task_event Description needed.
        void set_event(profiling::TaskEvent* task_event) {
            task_event_ = task_event;
        }

        /// Collect info on the task and record the submit time.
        void submit() {
            submit_time_ = wall_time();
            this->get_id(id_);
        }
#endif // MADNESS_TASK_PROFILING

        /// Object that is used to convert function and member function pointers into `void*`.

        /// \note This is technically not supported by the C++ standard but
        /// it will likely not cause any issues here (famous last words?).
        /// \todo Descriptions needed.
        /// \tparam T Description needed.
        template <typename T>
        union FunctionPointerGrabber {
            T in; ///< \todo Description needed.
            void* out; ///< \todo Description needed.
        };

    protected:

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam fnT Description needed.
        /// \param[in,out] id Description needed.
        /// \param[in] fn Description needed.
        /// \return Description needed.
        template <typename fnT>
        static typename std::enable_if<detail::function_traits<fnT>::value ||
                detail::memfunc_traits<fnT>::value>::type
        make_id(std::pair<void*,unsigned short>& id, fnT fn) {
            FunctionPointerGrabber<fnT> poop;
            poop.in = fn;
            id.first = poop.out;
            id.second = 1ul;
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed. What is the purpose of the second argument?
        /// \tparam fnobjT Description needed.
        /// \param[in,out] id Description needed.
        template <typename fnobjT>
        static typename std::enable_if<!(detail::function_traits<fnobjT>::value ||
                detail::memfunc_traits<fnobjT>::value) >::type
        make_id(std::pair<void*,unsigned short>& id, const fnobjT&) {
            id.first = reinterpret_cast<void*>(const_cast<char*>(typeid(fnobjT).name()));
            id.second = 2ul;
        }

    private:

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in,out] id Description needed.
        virtual void get_id(std::pair<void*,unsigned short>& id) const {
            id.first = nullptr;
            id.second = 0ul;
        }

#ifndef HAVE_INTEL_TBB

        Barrier* barrier; ///< Barrier, only allocated for multithreaded tasks.
        AtomicInt count; ///< Used to count threads as they start.

    	/// Returns true for the one thread that should invoke the destructor.

        /// \return True for the one thread that should invoke the destructor.
    	bool run_multi_threaded() {
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
        }

    public:

        /// Default constructor.
        PoolTaskInterface()
            : TaskAttributes()
            , barrier(nullptr)
#if HAVE_PARSEC
            , parsec_task(ParsecRuntime::task(is_high_priority(), this))
#endif
        {
    	    count = 0;
    	}

        /// Constructor setting the specified task attributes.

        /// \param[in] attr The task attributes.
        explicit PoolTaskInterface(const TaskAttributes& attr)
            : TaskAttributes(attr)
            , barrier(attr.get_nthread()>1 ? new Barrier(attr.get_nthread()) : 0)
#if HAVE_PARSEC
            , parsec_task(ParsecRuntime::task(is_high_priority(), this))
#endif
        {
            count = 0;
        }

        /// Destructor.
        /// \todo Should we either use a unique_ptr for barrier or check that barrier != nullptr here?
        virtual ~PoolTaskInterface() {
#if HAVE_PARSEC
          *(reinterpret_cast<PoolTaskInterface**>(&(parsec_task->locals[0]))) = nullptr;
          ParsecRuntime::delete_parsec_task(parsec_task);
          parsec_task = nullptr;
#endif
            delete barrier;
        }

        /// Call this to reset the number of threads before the task is submitted.

        /// Once a task has been constructed, /c TaskAttributes::set_nthread()
        /// is insufficient because a multithreaded task includes a barrier
        /// that needs to know the number of threads.
        ///
        /// \param[in] nthread The new number of threads.
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
#if HAVE_PARSEC
	    //////////// Parsec Related Begin ////////////////////
	    parsec_task_t                       *parsec_task;
	    //////////// Parsec Related End   ///////////////////
#endif

#else

    public:

        /// Default constructor.
        PoolTaskInterface() : TaskAttributes() {
	}

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in] attr Description needed.
        explicit PoolTaskInterface(const TaskAttributes& attr) :
            TaskAttributes(attr)
        {
	}

        /// Destructor.
        virtual ~PoolTaskInterface() = default;

        /// Call this to reset the number of threads before the task is submitted

        /// Once a task has been constructed /c TaskAttributes::set_nthread()
        /// is insufficient because a multithreaded task includes a
        /// barrier that needs to know the number of threads.
        void set_nthread(int nthread) {
            if (nthread != get_nthread())
                TaskAttributes::set_nthread(nthread);
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \return Description needed.
        void execute() {
            const int nthread = get_nthread();
            run( TaskThreadEnv(nthread, 0) );
        }

#endif // HAVE_INTEL_TBB

        /// Override this method to implement a multi-threaded task.

        /// \c info.nthread() will be the number of threads collaborating on this task.
        ///
        /// \c info.id() will be the index of the current thread \c id=0,...,nthread-1.
        ///
        /// \c info.barrier() will be a barrier for all of the threads, and returns
        ///     true for the last thread to enter the barrier (other threads get false).
        ///
        /// \todo Description needed.
        /// \param[in] info Description needed.


        virtual void run(const TaskThreadEnv& info) = 0;

        };

    /// A no-operation task used for various purposes.
    class PoolTaskNull : public PoolTaskInterface {
    public:
        /// Execution function that does nothing.
        void run(const TaskThreadEnv& /*info*/) {}

        /// Destructor.
        virtual ~PoolTaskNull() {}

    private:
        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \param[in,out] id Description needed.
        virtual void get_id(std::pair<void*,unsigned short>& id) const {
            PoolTaskInterface::make_id(id, &PoolTaskNull::run);
        }
    };

    /// \c ThreadPool thread object.

    /// This class holds thread local data for thread pool threads. It can be
    /// accessed via \c ThreadBase::this_thread().
    class ThreadPoolThread : public Thread {
    private:
        // Thread local data for thread pool
#ifdef MADNESS_TASK_PROFILING
        profiling::TaskProfiler profiler_; ///< \todo Description needed.
#endif // MADNESS_TASK_PROFILING

    public:
        ThreadPoolThread() : Thread() { }
        virtual ~ThreadPoolThread() = default;

#ifdef MADNESS_TASK_PROFILING
        /// Task profiler accessor.

        /// \todo Description needed.
        /// \return Description needed.
        profiling::TaskProfiler& profiler() {
            return profiler_;
        }
#endif // MADNESS_TASK_PROFILING
    };

    /// A singleton pool of threads for dynamic execution of tasks.

    /// \attention You must instantiate the pool while running with just one
    /// thread.
    class ThreadPool {
    public:
      // non-copyable and non-movable
      ThreadPool(const ThreadPool&) = delete;
      ThreadPool(ThreadPool&&) = delete;
      void operator=(const ThreadPool&) = delete;
      void operator=(ThreadPool&&) = delete;

      /// Get the number of threads from the environment.

      /// \return The number of threads.
      static int default_nthread();

    private:
        friend class WorldTaskQueue;

        // Thread pool data
        ThreadPoolThread *threads; ///< Array of threads.
        ThreadPoolThread main_thread; ///< Placeholder for main thread tls.
        DQueue<PoolTaskInterface*> queue; ///< Queue of tasks.
        int nthreads; ///< Number of threads.
        volatile bool finish; ///< Set to true when time to stop.
        AtomicInt nfinished; ///< Thread pool exit counter.

        // Static data
        static ThreadPool* instance_ptr; ///< Singleton pointer.
        static const int nmax = 128; ///< Number of task a worker thread will pop from the task queue
        static double await_timeout; ///< Waiter timeout.

#if defined(HAVE_IBMBGQ) and defined(HPM)
        static unsigned int main_hpmctx; ///< HPM context for main thread.
#endif
        /// The constructor is private to enforce the singleton model.

        /// \todo Description needed.
        /// \param[in] nthread Description needed.
        ThreadPool(int nthread=-1);

       /// Run the next task.

        /// \todo Verify and complete this documentation.
        /// \param[in] wait Block of true.
        /// \param[in,out] this_thread Description needed.
        /// \return True if a task was run.
        bool run_task(bool wait, ThreadPoolThread* this_thread) {
#if HAVE_INTEL_TBB
            MADNESS_EXCEPTION("run_task should not be called when using Intel TBB", 1);
#else

            if (!wait && queue.empty()) return false;
            std::pair<PoolTaskInterface*,bool> t = queue.pop_front(wait);
#ifdef MADNESS_TASK_PROFILING
            profiling::TaskEventList* event_list =
                    this_thread->profiler().new_list(1);
#endif // MADNESS_TASK_PROFILING
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

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in] wait Description needed.
        /// \param[in,out] this_thread Description needed.
        /// \return Description needed.
        bool run_tasks(bool wait, ThreadPoolThread* const this_thread) {
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

            PoolTaskInterface* taskbuf[nmax];
            int ntask = queue.pop_front(nmax, taskbuf, wait);
#ifdef MADNESS_TASK_PROFILING
            profiling::TaskEventList* event_list =
                    this_thread->profiler().new_list(ntask);
#endif // MADNESS_TASK_PROFILING
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
#if HAVE_PARSEC
            ////////////////// Parsec Related Begin //////////////////
            if(0 == ntask) {
                ntask = parsec_runtime->test();
            }
            ///////////////// Parsec Related End ////////////////////
#endif
            return (ntask>0);
#endif
        }

        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \param[in,out] thread Description needed.
        void thread_main(ThreadPoolThread* const thread);

        /// Forwards the thread to bound member function.

        /// \todo Descriptions needed.
        /// \param[in] v Description needed.
        /// \return Description needed.
        static void* pool_thread_main(void *v);

    public:
        /// Return a pointer to the only instance, constructing as necessary.

        /// \return A pointer to the only instance.
        static ThreadPool* instance() {
#ifndef MADNESS_ASSERTIONS_DISABLE
            if(! instance_ptr) {
                std::cerr << "!!! ERROR: The thread pool has not been initialized.\n"
                          << "!!! ERROR: Call madness::initialize before submitting tasks to the task queue.\n";
                MADNESS_EXCEPTION("ThreadPool::instance_ptr is NULL", 0);
            }
#endif
            return instance_ptr;
        }

	void flush_prebuf() {
#if !(defined(HAVE_INTEL_TBB) || defined(HAVE_PARSEC))
	  queue.lock_and_flush_prebuf();
#endif
	}

#if HAVE_PARSEC
	    ////////////////// Parsec Related Begin //////////////////
	    static ParsecRuntime *parsec_runtime;
	    ///////////////// Parsec Related End ////////////////////
#endif

#if HAVE_INTEL_TBB
        static std::unique_ptr<tbb::global_control> tbb_control; ///< \todo Description needed.
        static std::unique_ptr<tbb::task_arena>     tbb_arena;
#endif

        /// Please invoke while in a single-threaded environment.

        /// \todo Verify documentation.
        /// \param[in] nthread The number of threads.
        static void begin(int nthread=-1);

        /// \todo Description needed.
        static void end();

        /// Add a new task to the pool.

        /// \todo Description needed.
        /// \param[in,out] task Description needed.
        static void add(PoolTaskInterface* task) {
#ifdef MADNESS_TASK_PROFILING
            task->submit();
#endif // MADNESS_TASK_PROFILING

#if HAVE_PARSEC
            //////////// Parsec Related Begin ////////////////////
            parsec_runtime->schedule(task);
            //////////// Parsec Related End ////////////////////
#elif HAVE_INTEL_TBB
//#ifdef MADNESS_CAN_USE_TBB_PRIORITY
//            if(task->is_high_priority())
//                tbb::task::enqueue(*task, tbb::priority_high);
//            else
//#endif  // MADNESS_CAN_USE_TBB_PRIORITY
            tbb_arena->enqueue(
                    //use unique_ptr to automatically delete task ptr
                    [task_p = std::unique_ptr<PoolTaskInterface>(task)] () noexcept {
                        //exceptions are not expected here, as nobody will catch them for enqueued tasks
                        task_p->execute();
            });
#else
            if (!task) MADNESS_EXCEPTION("ThreadPool: inserting a NULL task pointer", 1);
            int task_threads = task->get_nthread();
            // Currently multithreaded tasks must be shoved on the end of the q
            // to avoid a race condition as multithreaded task is starting up
            if (task->is_high_priority() && (task_threads == 1)) {
                instance()->queue.push_front(task);
            }
            else {
                instance()->queue.push_back(task, task_threads);
            }
#endif // HAVE_INTEL_TBB
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam opT Description needed.
        /// \param[in,out] op Description needed.
        template <typename opT>
        void scan(opT& op) {
            queue.scan(op);
        }

        /// Add a vector of tasks to the pool.

        /// \param[in] tasks Vector of tasks to add to the pool.
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

        /// An otherwise idle thread can all this to run a task.

        /// \return True if a task was run.
        static bool run_task() {
#ifdef HAVE_INTEL_TBB
            return false;
#else

#ifdef MADNESS_TASK_PROFILING
            ThreadPoolThread* const thread = static_cast<ThreadPoolThread*>(ThreadBase::this_thread());
#else
            ThreadPoolThread* const thread = nullptr;
#endif // MADNESS_TASK_PROFILING

            return instance()->run_tasks(false, thread);
#endif // HAVE_INTEL_TBB
        }

        /// Returns the number of threads in the pool.

        /// \return The number of threads in the pool.
        static std::size_t size() {
            return instance()->nthreads;
        }

        /// Returns the number of tasks in the queue.

        /// \return The number of tasks in the queue.
        static std::size_t queue_size() {
            return instance()->queue.size();
        }

        /// Returns queue statistics.

        /// \return Queue statistics.
        static const DQStats& get_stats();

        /// Access the pool thread array
        /// \return ptr to the pool thread array, its size is given by \c size()
        static const ThreadPoolThread* get_threads() {
          return const_cast<const ThreadPoolThread*>(instance()->threads);
        }

        /// Gracefully wait for a condition to become true, executing any tasks in the queue.

        /// Probe should be an object that, when called, returns the status.
        /// \todo Descriptions needed/need verification.
        /// \tparam Probe Type of the probe.
        /// \param[in] probe The probe.
        /// \param[in] dowork Do work while waiting - default is true
        /// \param[in] sleep Sleep instead of spin while waiting (e.g., to avoid pounding on MPI) - default is false
        template <typename Probe>
	    static void await(const Probe& probe, bool dowork = true, bool sleep = false) {
          if (!probe()) {
            double start = cpu_time();
            const double timeout = await_timeout;
            int counter = 0;

            MutexWaiter waiter;
            while (!probe()) {

                const bool working = (dowork ? ThreadPool::run_task() : false);
                const double current_time = cpu_time();

                if (working) {	// Reset timeout logic
                    waiter.reset();
                    start = current_time;
                    counter = 0;
                } else {
                    if(((current_time - start) > timeout) && (timeout > 1.0)) { // Check for timeout
                      std::cerr << "!!MADNESS: Hung queue?" << std::endl;
                      if (counter++ > 3) {
                        const long bufsize=256;
                        char errstr[bufsize];
                        snprintf(errstr,bufsize, "ThreadPool::await() timed out after %.1lf seconds", timeout);
                        throw madness::MadnessException(errstr, 0, 1,
                                                        __LINE__, __FUNCTION__,
                                                        __FILE__);
                      }
                    }
		    if (sleep) {
		      // THIS NEEDS TO BECOME AN EXTERNAL PARAMETER
		      // Problem is exacerbated when running with many
		      // (e.g., 512 or more) send/recv buffers, and
		      // also with many threads.  More outstanding
		      // requests means each call into MPI takes
		      // longer and more threads means more calls in
		      // spots where all threads are messaging.  Old
		      // code was OK on dancer.icl.utk.edu with just
		      // 32 bufs and 20 threads, but 512 bufs caused
		      // intermittent hangs I think due to something
		      // not being able to make progress or general
		      // confusion (this with MPICH) ... maybe using a
		      // fair mutex somewhere would help.
		      //
		      // 100us is a long time ... will try 10us. mmm ... perhaps need 100 at least on dancer with 17 threads per node
		      myusleep(100);
		    }
		    else {
		      waiter.wait();
		    }
                }
            }
          }  // if !probe()
        }

        /// Destructor.
        ~ThreadPool() {
#if HAVE_PARSEC
            ////////////////// Parsec related Begin /////////////////
            delete parsec_runtime;
            ////////////////// Parsec related End /////////////////
#elif HAVE_INTEL_TBB
#else
            delete[] threads;
#endif
        }

        /// \sa madness::threadpool_wait_policy
        static void set_wait_policy(
          WaitPolicy policy,
          int sleep_duration_in_microseconds = 0) {
#if !HAVE_INTEL_TBB && !HAVE_PARSEC
          instance()->queue.set_wait_policy(policy,
                                            sleep_duration_in_microseconds);
#endif
        }

    };

    // clang-format off
    /// Controls how aggressively ThreadPool holds on to the OS threads
    /// while waiting for work. Currently useful only for Pthread pool when it's using spinlocks;
    /// NOT used for TBB or PaRSEC.
    /// \param policy specifies how to wait for work;
    ///        - WaitPolicy::Busy -- threads are kept busy (default); recommended when intensive work is only performed by MADNESS threads
    ///        - WaitPolicy::Yield -- thread yields; recommended when intensive work is performed primarily by non-MADNESS threads
    ///        - WaitPolicy::Sleep -- thread sleeps for \p sleep_duration_in_microseconds ; recommended when intensive work is performed by MADNESS nd non-MADNESS threads
    /// \param sleep_duration_in_microseconds if `policy==WaitPolicy::Sleep` this specifies the duration of sleep, in microseconds
    // clang-format on
    inline void threadpool_wait_policy(WaitPolicy policy,
                                       int sleep_duration_in_microseconds = 0) {
      ThreadPool::set_wait_policy(policy, sleep_duration_in_microseconds);
    }

    /// @}
}

#endif // MADNESS_WORLD_THREAD_H__INCLUDED
