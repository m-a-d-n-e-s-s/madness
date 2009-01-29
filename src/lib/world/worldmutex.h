#ifndef WORLD_MUTEX_H
#define WORLD_MUTEX_H

#include <madness_config.h>
#ifdef ON_A_MAC 
#include <libkern/OSAtomic.h>
typedef OSSpinLock pthread_spinlock_t;

inline void pthread_spin_init(pthread_spinlock_t* p, int mode) {*p=0;}
inline int pthread_spin_trylock(pthread_spinlock_t* p) {return !OSSpinLockTry(p);}
inline int pthread_spin_lock(pthread_spinlock_t* p) {OSSpinLockLock(p); return 0;}
inline int pthread_spin_unlock(pthread_spinlock_t* p) {OSSpinLockUnlock(p); return 0;}
inline void pthread_spin_destroy(pthread_spinlock_t* p) {}
#endif

#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <utility>
#include <vector>
#include <world/madatomic.h>
#include <world/nodefaults.h>

/// \file worldmutex.h
/// \brief Implements Mutex, MutexFair, Spinlock, ConditionVariable


namespace madness {

    inline void cpu_relax(){asm volatile ( "rep;nop" : : : "memory" );}

    class MutexWaiter {
    private:
        unsigned int count;

        /// Yield for specified number of microseconds unless dedicated CPU
        void yield(int us) {
            usleep(us);
        }
        
    public:
        MutexWaiter() : count(0) {}

        void reset() {count = 0;}

        void wait() {
#ifdef HAVE_CRAYXT
            // The value of 300 below is purely empirical but apparently a substantial
            // backoff (or something that is a by-product of waiting) is necessary
            // to avoid clobbering the memory subsystem while spinning on the taskq.
            // The time is  "Time to  run 100000 chain of tasks" from running world.
            // 1000--> 2+us
            // 400 --> 1.7us
            // 300 --> 1.7us
            // 250 --> 2us
            // 200 --> 3.6us
            // 100 --> 40+us (ouch!)
            for (int i=0; i<300; i++)  cpu_relax();
#else
            const unsigned int nspin  = 1000;    // Spin for 1,000 calls
            const unsigned int nsleep = 100000;  // Sleep 10us for 100,000 calls = 1s
            if (count++ < nspin) return;
            else if (count < nsleep) yield(10);
            else yield(10000);
#endif
        }
    };


#ifdef WORLD_MUTEX_ATOMIC
    /// Mutex using spin locks and atomic operations
    
    /// Possibly much *slower* than pthread operations.  Also cannot
    /// use these Mutexes with pthread condition variables.
    class Mutex : NO_DEFAULTS {
    private:
        mutable MADATOMIC_INT flag;
    public:
        /// Make and initialize a mutex ... initial state is unlocked
        Mutex() {
            MADATOMIC_INT_SET(&flag,1L);
        }
        
        /// Try to acquire the mutex ... return true on success, false on failure
        bool try_lock() const {
            if (MADATOMIC_INT_DEC_AND_TEST(&flag)) return true;
            MADATOMIC_INT_INC(&flag);
            return false;
        }
            
        /// Acquire the mutex waiting if necessary
        void lock() const {
            MutexWaiter waiter;
            while (!try_lock()) waiter.wait();
        }
        
        /// Free a mutex owned by this thread
        void unlock() const {
            MADATOMIC_INT_INC(&flag);
        }

        //pthread_mutex_t* ptr() const {throw "there is no pthread_mutex";}
        
        virtual ~Mutex() {};
    };

#else
    
    /// Mutex using pthread mutex operations
    class Mutex {
    private:
        mutable pthread_mutex_t mutex;

        /// Copy constructor is forbidden
        Mutex(const Mutex& m) {}
        
        /// Assignment is forbidden
        void operator=(const Mutex& m) {}
        
    public:
        /// Make and initialize a mutex ... initial state is unlocked
        Mutex() {
            pthread_mutex_init(&mutex, 0);
        }

        /// Try to acquire the mutex ... return true on success, false on failure
        bool try_lock() const {
            return pthread_mutex_trylock(&mutex)==0;
        }
        
        /// Acquire the mutex waiting if necessary
        void lock() const {
            if (pthread_mutex_lock(&mutex)) throw "failed acquiring mutex";
        }
        
        /// Free a mutex owned by this thread
        void unlock() const {
            if (pthread_mutex_unlock(&mutex)) throw "failed releasing mutex";
        }
        
        /// Return a pointer to the pthread mutex for use by a condition variable
        pthread_mutex_t* ptr() const {
            return &mutex;
        }

        virtual ~Mutex() {
            pthread_mutex_destroy(&mutex);
        };
    };

 
#endif


    /// Mutex that is applied/released at start/end of a scope

    /// The mutex must provide lock and unlock methods
    template <class mutexT = Mutex>
    class ScopedMutex {
        const mutexT* m;
    public:
        ScopedMutex(const mutexT* m) : m(m) {m->lock();}
        ScopedMutex(const mutexT& m) : m(&m) {m.lock();}
        virtual ~ScopedMutex() {m->unlock();}
    };

    /// Spinlock using pthread spinlock operations
    class Spinlock {
    private:
        mutable pthread_spinlock_t spinlock;

        /// Copy constructor is forbidden
        Spinlock(const Spinlock& m) {}
        
        /// Assignment is forbidden
        void operator=(const Spinlock& m) {}
        
    public:
        /// Make and initialize a spinlock ... initial state is unlocked
        Spinlock() {
            pthread_spin_init(&spinlock, PTHREAD_PROCESS_PRIVATE);
        }

        /// Try to acquire the spinlock ... return true on success, false on failure
        bool try_lock() const {
            return pthread_spin_trylock(&spinlock)==0;
        }
        
        /// Acquire the spinlock waiting if necessary
        void lock() const {
            if (pthread_spin_lock(&spinlock)) throw "failed acquiring spinlock";
        }
        
        /// Free a spinlock owned by this thread
        void unlock() const {
            if (pthread_spin_unlock(&spinlock)) throw "failed releasing spinlock";
        }
        
        virtual ~Spinlock() {
            pthread_spin_destroy(&spinlock);
        };
    };


    class MutexReaderWriter : private Spinlock, NO_DEFAULTS {
        volatile mutable int nreader;
        volatile mutable bool writeflag;
    public:
        static const int NOLOCK=0;
        static const int READLOCK=1;
        static const int WRITELOCK=2;

        MutexReaderWriter() : nreader(0), writeflag(false) 
        {}

        bool try_read_lock() const {
            ScopedMutex<Spinlock> protect(this);
            bool gotit = !writeflag;
            if (gotit) nreader++;
            return gotit;
        }
        
        bool try_write_lock() const {
            ScopedMutex<Spinlock> protect(this);
            bool gotit = (!writeflag) && (nreader==0);
            if (gotit) writeflag = true;
            return gotit;
        }

        bool try_lock(int lockmode) const {
            if (lockmode == READLOCK) {
                return try_read_lock();
            }
            else if (lockmode == WRITELOCK) {
                return try_write_lock();
            }
            else if (lockmode == NOLOCK) {
                return true;
            }
            else {
                throw "MutexReaderWriter: try_lock: invalid lock mode";
            }
        }

        bool try_convert_read_lock_to_write_lock() const {
            ScopedMutex<Spinlock> protect(this);
            bool gotit = (!writeflag) && (nreader==1);
            if (gotit) {
                nreader = 0;
                writeflag = true;
            }
            return gotit;
        }

        void read_lock() const {
            //MutexWaiter waiter;
            while (!try_read_lock()) cpu_relax(); //waiter.wait();
        }

        void write_lock() const {
            //MutexWaiter waiter;
            while (!try_write_lock()) cpu_relax(); //waiter.wait();
        }

        void lock(int lockmode) const {
            //MutexWaiter waiter;
            while (!try_lock(lockmode)) cpu_relax(); //waiter.wait();
        }

        void read_unlock() const {
            ScopedMutex<Spinlock> protect(this);
            nreader--;
        }

        void write_unlock() const {
            // Only a single thread should be setting writeflag but
            // probably still need the mutex just to get memory fence?
            ScopedMutex<Spinlock> protect(this);
            writeflag = false;
        }

        void unlock(int lockmode) const {
            if (lockmode == READLOCK) read_unlock();
            else if (lockmode == WRITELOCK) write_unlock();
            else if (lockmode != NOLOCK) throw "MutexReaderWriter: try_lock: invalid lock mode";            
        }

        /// Converts read to write lock without releasing the read lock

        /// Note that deadlock is guaranteed if two+ threads wait to convert at the same time.
        void convert_read_lock_to_write_lock() const {
            //MutexWaiter waiter;
            while (!try_convert_read_lock_to_write_lock()) cpu_relax(); //waiter.wait();
        }

        /// Always succeeds immediately
        void convert_write_lock_to_read_lock() const {
            ScopedMutex<Spinlock> protect(this);
            nreader++;
            writeflag=false;
        }

        virtual ~MutexReaderWriter(){};
    };

    /// Scalable and fair condition variable (spins on local value)
    class ConditionVariable : public Spinlock {
    public:
        static const int MAX_NTHREAD = 64;
        mutable volatile int front;
        mutable volatile int back;
        mutable volatile bool* volatile q[MAX_NTHREAD]; // Circular buffer of flags

    public:
        ConditionVariable() : front(0), back(0) {}

        /// You should acquire the mutex before waiting
        void wait() const {
            // We put a pointer to a thread-local variable at the
            // end of the queue and wait for that value to be set,
            // thus generate no memory traffic while waiting.
            volatile bool myturn = false;
            int b = back;
            q[b] = &myturn;
            b++;
            if (b >= MAX_NTHREAD) back = 0;
            else back = b;
            
            unlock(); // Release lock before blocking
            while (!myturn) cpu_relax();
            lock();
        }
        
        /// You should acquire the mutex before signalling
        void signal() const {
            if (front != back) {
                int f = front;
                int ff = f + 1;
                if (ff >= MAX_NTHREAD) 
                    front = 0;
                else 
                    front = ff;

                *q[f] = true;
            }
        }


        /// You should acquire the mutex before broadcasting
        void broadcast() const {
            while (front != back) 
                signal();
        }

        virtual ~ConditionVariable() {}
    };

    /// A scalable and fair mutex (not recursive)

    /// Needs rewriting to use the CV above and do we really
    /// need this if using pthread_mutex .. why not pthread_cv?
    class MutexFair : private Spinlock {
    private:
        static const int MAX_NTHREAD = 64;
        mutable volatile bool* volatile q[MAX_NTHREAD]; 
        mutable volatile int n;
        mutable volatile int front;
        mutable volatile int back;

    public:
        MutexFair() : n(0), front(0), back(0) {};

        void lock() const {
            volatile bool myturn = false;
            Spinlock::lock();
            n++;
            if (n == 1) {
                myturn = true;
            }
            else {
                int b = back + 1;
                if (b >= MAX_NTHREAD) b = 0;
                q[b] = &myturn;
                back = b;
            }
            Spinlock::unlock();

            while (!myturn) cpu_relax();
        }

        void unlock() const {
            volatile bool* p = 0;
            Spinlock::lock();
            n--;
            if (n > 0) {
                int f = front + 1;
                if (f >= MAX_NTHREAD) f = 0;
                p = q[f];
                front = f;
            }
            Spinlock::unlock();
            if (p) *p = true;
        }

        bool try_lock() const {
            bool got_lock;

            Spinlock::lock();
            int nn = n;
            got_lock = (nn == 0);
            if (got_lock) n = nn + 1;
            Spinlock::unlock();

            return got_lock;
        }
             
    };


    /// Attempt to acquire two locks without blocking holding either one

    /// The code will first attempt to acquire mutex m1 and if successful
    /// will then attempt to acquire mutex m2.  
    inline bool try_two_locks(const Mutex& m1, const Mutex& m2) {
        if (!m1.try_lock()) return false;
        if (m2.try_lock()) return true;
        m1.unlock();
        return false;
    }


    /// Simple wrapper for Pthread condition variable with its own mutex

    /// Use this when you need to block without consuming cycles.
    /// Scheduling granularity is at the level of kernel ticks.
    class PthreadConditionVariable : NO_DEFAULTS {
    private:
        mutable pthread_cond_t cv;
        mutable pthread_mutex_t mutex;

    public:
        PthreadConditionVariable() {
            pthread_cond_init(&cv, NULL);
            pthread_mutex_init(&mutex, 0);
        }            

        pthread_mutex_t& get_pthread_mutex() {return mutex;}

        void lock() const {
            if (pthread_mutex_lock(&mutex)) throw "ConditionVariable: acquiring mutex";
        }

        void unlock() const {
            if (pthread_mutex_unlock(&mutex)) throw "ConditionVariable: releasing mutex";
        }

        /// You should have acquired the mutex before entering here
        void wait() const {
            pthread_cond_wait(&cv,&mutex);
        }

        void signal() const {
            if (pthread_cond_signal(&cv)) throw "ConditionalVariable: signalling failed";
        }

        void broadcast() const {
            if (pthread_cond_broadcast(&cv)) throw "ConditionalVariable: signalling failed";
        }

        virtual ~PthreadConditionVariable() {
            pthread_mutex_destroy(&mutex);
            pthread_cond_destroy(&cv);
        }
    };

#ifdef USE_SPINLOCKS
    typedef ConditionVariable CONDITION_VARIABLE_TYPE ;
    typedef Spinlock SPINLOCK_TYPE;
    typedef MutexFair SCALABLE_MUTEX_TYPE;
#else
    typedef PthreadConditionVariable CONDITION_VARIABLE_TYPE;
    typedef Mutex SPINLOCK_TYPE;
    typedef Mutex SCALABLE_MUTEX_TYPE;
#endif

}

#endif
