#ifndef MAD_WORLDTHREAD_H
#define MAD_WORLDTHREAD_H

/// \file worldthread.h
/// \brief Implements Dqueue, Thread, ThreadBase and ThreadPool

#include <world/safempi.h>
#include <world/worldexc.h>
#include <world/print.h>
#include <world/worldmutex.h>
#include <world/worldpapi.h>

namespace madness {

    void error(const char *msg);
    
    struct DQStats { // Dilly bar, blizzard, ...
        uint64_t nmax;          //< Lifetime max. entries in the queue
        uint64_t npush_back;    //< #calls to push_back
        uint64_t npush_front;   //< #calls to push_front
        uint64_t npop_back;     //< #calls to pop_back
        uint64_t npop_front;    //< #calls to pop_front
        uint64_t ngrow;         //< #calls to grow

        DQStats() 
            : nmax(0), npush_back(0), npush_front(0), npop_back(0), npop_front(0), ngrow(0)
        {}
    };


    /// A thread safe, fast but simple doubled-ended queue.

    /// Since the point is speed, the implementation is a circular
    /// buffer rather than a linked list so as to avoid the new/del
    /// overhead.  It will grow as needed, but presently will not
    /// shrink.  Had to modify STL API to make things thread safe.
    ///
    /// This queue is part of the horror show ... again previous
    /// attempts to disentangle from the "CV" and the thread pool have
    /// lead to a big slow down on the XT.
    template <typename T>
    class DQueue : private ConditionVariable {
        volatile size_t sz;              ///< Current capacity
        volatile T* volatile buf;        ///< Actual buffer
        volatile int _front;  ///< Index of element at front of buffer
        volatile int _back;    ///< Index of element at back of buffer
        volatile size_t n;      ///< Number of elements in the buffer
        DQStats stats;

        void grow() {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            stats.ngrow++;
            if (sz != n) throw "assertion failure in dqueue::grow";
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
            sanity_check();
        }

        void sanity_check() const {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
            int num = _back - _front + 1;
            if (num < 0) num += sz;
            if (num==int(sz) && n==0) num=0;
            if (num==0 && n==sz) num=sz;
            //if (long(n) != num) print("size",sz,"front",_front,"back",_back,"n",n,"num",num);
            if (long(n) != num) throw "assertion failure in dqueue::sanity";
        }

    public:
        DQueue(size_t hint=32768) 
            : sz(hint>2 ? hint : 2)
            , buf(new T[sz])
            , _front(sz/2)
            , _back(_front-1)
            , n(0)
        {}

        virtual ~DQueue() {
            delete buf;
        }

        /// Insert value at front of queue
        void push_front(const T& value) {
            madness::ScopedMutex<ConditionVariable> obolus(this);
            stats.npush_front++;
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

            signal();
        }

        /// Insert element at back of queue
        void push_back(const T& value) {
            madness::ScopedMutex<ConditionVariable> obolus(this);
            stats.npush_back++;
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

            int b = _back + 1;
            if (b >= int(ss)) b = 0;
            buf[b] = value;
            _back = b;

            signal();
        }

        /// Pop value off the front of queue
        std::pair<T,bool> pop_front(bool wait) {
            madness::ScopedMutex<ConditionVariable> obolus(this);
            stats.npop_front++;
            if (wait) ConditionVariable::wait();
            size_t nn = n;
            if (nn) {
                //sanity_check();
                n = nn - 1;
                
                int f = _front;
                T result = buf[f];
                buf[f] = T();   // For better stupidity detection
                f++;
                if (f >= int(sz)) f = 0;
                _front = f;
                return std::pair<T,bool>(result,true);
            }
            else {
                return std::pair<T,bool>(T(),false);
            }
        }


        /// Pop value off the back of queue
        std::pair<T,bool> pop_back(bool wait=true) {
            madness::ScopedMutex<ConditionVariable> obolus(this);
            stats.npop_back++;
            if (wait) ConditionVariable::wait();
            bool status=true; 
            T result;
            size_t nn = n;
            if (nn) {
                //sanity_check();
                n = nn - 1;

                int b = _back;
                result = buf[b];
                b--;
                if (b < 0) b = sz-1;
                _back = b;
            }
            else {
                status = false;
            }
            return std::pair<T,bool>(result,status);
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
            } catch (const MPI::Exception& e) {
                //        print(e);
                error("caught an MPI exception");
            } catch (const madness::MadnessException& e) {
                print(e);
                error("caught a MADNESS exception");
            } catch (const char* s) {
                print(s);
                error("caught a string exception");
            } catch (char* s) {
                print(s);
                error("caught a string exception");
            } catch (const std::string& s) {
                print(s);
                error("caught a string (class) exception");
            } catch (const std::exception& e) {
                print(e.what());
                error("caught an STL exception");
            } catch (...) {
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
            pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); 
            if (pthread_create(&id, &attr, &ThreadBase::main, (void *) this))
                throw "failed creating thread";
            
            pthread_attr_destroy(&attr);
        }
        
        /// A thread can call this to terminate its execution
        static void exit() {pthread_exit(0);}
        
        /// Get the pthread id of this thread (if running)
        const pthread_t& get_id() const {return id;}

        /// Get index of thread in ThreadPool (0,...,nthread-1) or -1 if not in ThreadPool
        int get_pool_thread_index() const {return pool_num;}
        
        /// Cancel this thread
        int cancel() const {return pthread_cancel(get_id());}

        /// Specify the affinity pattern or how to bind threads to cpus
        static void set_affinity_pattern(const bool bind[3], const int cpu[3]) {
            memcpy(ThreadBase::bind, bind, 3*sizeof(bool));
            memcpy(ThreadBase::cpulo, cpu, 3*sizeof(int));

            int ncpu = sysconf(_SC_NPROCESSORS_CONF);
            if (ncpu <= 0) throw "ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)";

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
            return;
            if (logical_id < 0 || logical_id > 2) {
                std::cout << "ThreadBase: set_affinity: logical_id bad?" << std::endl;
                return;
            }

            if (!bind[logical_id]) return;
                
            int ncpu = sysconf(_SC_NPROCESSORS_CONF);
            if (ncpu <= 0) {
                std::cout << "ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)" << std::endl;
                return;
            }

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

            //std::cout << "BINDING THREAD: id " << logical_id << " ind " << ind << " lo " << lo << " hi " << hi << " ncpu " << ncpu << std::endl;
            
            cpu_set_t mask;
            CPU_ZERO(&mask);
            for (int i=lo; i<=hi; i++) CPU_SET(i,&mask);
            if(sched_setaffinity( 0, sizeof(mask), &mask ) == -1 ) {
                perror("system error message");
                std::cout << "ThreadBase: set_affinity: Could not set cpu Affinity" << std::endl;

            }
        }
    };
    
    
    /// Simplified thread wrapper to hide pthread complexity
    class Thread : public ThreadBase {
        void* (*f)(void *); 
        void* args;
        
        void run() {f(args);}
    public:
        /// Default constructor ... must invoke \c start() to actually begin the thread.
        Thread() : f(0), args(0) {};
        
        /// Create a thread and start it running f(args)
        Thread(void* (*f) (void *), void* args=0) 
            : f(f), args(args) 
        {
            ThreadBase::start();
        }
        
        void start(void* (*f) (void *), void* args=0) {
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
    /// \c highpriority : indicates a high priority task.
    class TaskAttributes {
        unsigned long flags;
        static const unsigned long one = 1ul;
    public:
        static const unsigned long GENERATOR    = one;
        static const unsigned long STEALABLE    = one<<1;
        static const unsigned long HIGHPRIORITY = one<<2;

        TaskAttributes(unsigned long flags = 0) : flags(flags) {}

        bool is_generator() const {return flags&GENERATOR;}

        bool is_stealable() const {return flags&STEALABLE;}

        bool is_high_priority() const {return flags&HIGHPRIORITY;}

        void set_generator (bool generator_hint) {
            if (generator_hint) flags |= GENERATOR;
            else flags &= ~GENERATOR;
        }

        void set_stealable (bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
        }

        void set_highpriority (bool hipri) {
            if (hipri) flags |= HIGHPRIORITY;
            else flags &= ~HIGHPRIORITY;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & flags;
        }

        static TaskAttributes generator() {return TaskAttributes(GENERATOR);}

        static TaskAttributes hipri() {return TaskAttributes(HIGHPRIORITY);}
    };


    class PoolTaskInterface : public TaskAttributes {
    public:
        PoolTaskInterface() {}

        explicit PoolTaskInterface(const TaskAttributes& attr)
            : TaskAttributes(attr) {}

        virtual void run() = 0;

        virtual ~PoolTaskInterface() {}
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
        MADATOMIC_INT nfinished;

        static ThreadPool* instance_ptr;

        /// The constructor is private to enforce the singleton model
        ThreadPool(int nthread=-1) : nthreads(nthread), finish(false) {
            MADATOMIC_INT_SET(&nfinished,0);
            instance_ptr = this;
            if (nthreads < 0) nthreads = default_nthread();
            //std::cout << "POOL " << nthreads << std::endl;

            try {
                if (nthreads > 0) threads = new Thread[nthreads];
            else threads = 0;
            }
            catch (...) {
                throw "memory allocation failed";
            }

            for (int i=0; i<nthreads; i++) {
                threads[i].set_pool_thread_index(i);
                threads[i].start(pool_thread_main, (void *) (threads+i));
            }
        }

        ThreadPool(const ThreadPool&);           // Verboten
        void operator=(const ThreadPool&);       // Verboten

        /// Get number of threads from the environment
        int default_nthread() {
            int nthread;
            char *cnthread = getenv("POOL_NTHREAD");

            if (cnthread) {
                if (sscanf(cnthread, "%d", &nthread) != 1) 
                    throw "POOL_NTHREAD is not an integer";
            }
            else {
                nthread = sysconf(_SC_NPROCESSORS_CONF);
                if (nthread < 2) nthread = 2;
                nthread = nthread - 1; // One less than # physical processors
            }
            return nthread;
        }

        /// Run next task ... returns true if one was run ... blocks if wait is true
        bool run_task(bool wait) {
            std::pair<PoolTaskInterface*,bool> t = queue.pop_front(wait);
            if (t.second) {
                t.first->run();          // What we are here to do
                delete t.first;
            }
            return t.second;
        }

        void thread_main(Thread* thread) {
            thread->set_affinity(2, thread->get_pool_thread_index());
            while (!finish) 
                run_task(true);
            MADATOMIC_INT_INC(&nfinished);
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

        class PoolTaskNull : public PoolTaskInterface {
            void run() {};
        };

    public:
        /// Please invoke while in single threaded environment
        static void begin(int nthread=-1) {instance(nthread);}

        static void end() {
		    if (!instance_ptr) return;
            instance()->finish = true;
            for (int i=0; i<instance()->nthreads; i++) {
                add(new PoolTaskNull);
            }
            while (MADATOMIC_INT_GET(&instance()->nfinished) != instance()->nthreads);
        }

        /// Add a new task to the pool
        static void add(PoolTaskInterface* task) {
            if (!task) throw "ThreadPool: inserting a NULL task pointer";
            if (task->is_high_priority()) {
                instance()->queue.push_front(task);
            }
            else {
                instance()->queue.push_back(task);
            }
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
            return instance()->run_task(false);
        }

        /// Returns number of threads in the pool
        static size_t size() {
            return instance()->nthreads;
        }

        /// Returns queue statistics
        static const DQStats& get_stats() {
            return instance()->queue.get_stats();
        }

        ~ThreadPool() {};
    };

}
    

#endif
