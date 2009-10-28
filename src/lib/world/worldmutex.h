#ifndef WORLD_MUTEX_H
#define WORLD_MUTEX_H

#include <madness_config.h>
#ifdef ON_A_MAC
#include <libkern/OSAtomic.h>
typedef OSSpinLock pthread_spinlock_t;

inline void pthread_spin_init(pthread_spinlock_t* p, int mode) {
    *p=0;
}
inline int pthread_spin_trylock(pthread_spinlock_t* p) {
    return !OSSpinLockTry(p);
}
inline int pthread_spin_lock(pthread_spinlock_t* p) {
    OSSpinLockLock(p);
    return 0;
}
inline int pthread_spin_unlock(pthread_spinlock_t* p) {
    OSSpinLockUnlock(p);
    return 0;
}
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
#include <world/nodefaults.h>
#include <world/worldtime.h>
#include <world/atomicint.h>
#include <world/worldexc.h>

/// \file worldmutex.h
/// \brief Implements Mutex, MutexFair, Spinlock, ConditionVariable


namespace madness {

    class MutexWaiter {
    private:
        unsigned int count;

        /// Yield for specified number of microseconds unless dedicated CPU
        void yield(int us) {
            usleep(us);
        }

    public:
        MutexWaiter() : count(0) {}

        void reset() {
            count = 0;
        }

        void wait() {
            //#ifdef HAVE_CRAYXT
#ifdef USE_SPINLOCKS
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
            int result = pthread_mutex_lock(&mutex);
            if (result) MADNESS_EXCEPTION("failed acquiring mutex", result);
        }

        /// Free a mutex owned by this thread
        void unlock() const {
            int result = pthread_mutex_unlock(&mutex);
            if (result) MADNESS_EXCEPTION("failed releasing mutex", result);
        }

        /// Return a pointer to the pthread mutex for use by a condition variable
        pthread_mutex_t* ptr() const {
            return &mutex;
        }

        virtual ~Mutex() {
            pthread_mutex_destroy(&mutex);
        };
    };

    /// Mutex that is applied/released at start/end of a scope

    /// The mutex must provide lock and unlock methods
    template <class mutexT = Mutex>
    class ScopedMutex {
        const mutexT* m;
    public:
        ScopedMutex(const mutexT* m) : m(m) {
            m->lock();
        }
        ScopedMutex(const mutexT& m) : m(&m) {
            m.lock();
        }
        virtual ~ScopedMutex() {
            m->unlock();
        }
    };

#ifdef NEVER_SPIN
    typedef Mutex Spinlock;
#else
    /// Spinlock using pthread spinlock operations
    class Spinlock {
    private:
        //mutable pthread_spinlock_t spinlock  __attribute__ ((aligned (64)));
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
            int result = pthread_spin_lock(&spinlock);
            if (result) MADNESS_EXCEPTION("failed acquiring spinlock", result);
        }

        /// Free a spinlock owned by this thread
        void unlock() const {
            int result = pthread_spin_unlock(&spinlock);
            if (result) MADNESS_EXCEPTION("failed releasing spinlock", result);
        }

        virtual ~Spinlock() {
            pthread_spin_destroy(&spinlock);
        };
    };
#endif


    class MutexReaderWriter : private Spinlock, NO_DEFAULTS {
        volatile mutable int nreader;
        volatile mutable bool writeflag;
    public:
        static const int NOLOCK=0;
        static const int READLOCK=1;
        static const int WRITELOCK=2;

        MutexReaderWriter() : nreader(0), writeflag(false) {}

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
               MADNESS_EXCEPTION("MutexReaderWriter: try_lock: invalid lock mode", lockmode);
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
            while (!try_read_lock()) cpu_relax();
        }

        void write_lock() const {
            while (!try_write_lock()) cpu_relax();
        }

        void lock(int lockmode) const {
            while (!try_lock(lockmode)) cpu_relax();
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
            else if (lockmode != NOLOCK) MADNESS_EXCEPTION("MutexReaderWriter: try_lock: invalid lock mode", lockmode);
        }

        /// Converts read to write lock without releasing the read lock

        /// Note that deadlock is guaranteed if two+ threads wait to convert at the same time.
        void convert_read_lock_to_write_lock() const {
            while (!try_convert_read_lock_to_write_lock()) cpu_relax();
        }

        /// Always succeeds immediately
        void convert_write_lock_to_read_lock() const {
            ScopedMutex<Spinlock> protect(this);
            nreader++;
            writeflag=false;
        }

        virtual ~MutexReaderWriter() {};
    };

    /// Scalable and fair condition variable (spins on local value)
    class ConditionVariable : public Spinlock {
    public:
        static const int MAX_NTHREAD = 64;
        mutable volatile int back;
        mutable volatile int front;
        mutable volatile bool* volatile q[MAX_NTHREAD]; // Circular buffer of flags

    public:
        ConditionVariable() : back(0), front(0) {}

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

        pthread_mutex_t& get_pthread_mutex() {
            return mutex;
        }

        void lock() const {
            int result = pthread_mutex_lock(&mutex);
            if (result) MADNESS_EXCEPTION("ConditionVariable: acquiring mutex", result);
        }

        void unlock() const {
            int result = pthread_mutex_unlock(&mutex);
            if (result) MADNESS_EXCEPTION("ConditionVariable: releasing mutex", result);
        }

        /// You should have acquired the mutex before entering here
        void wait() const {
            pthread_cond_wait(&cv,&mutex);
        }

        void signal() const {
            int result = pthread_cond_signal(&cv);
            if (result) MADNESS_EXCEPTION("ConditionalVariable: signalling failed", result);
        }

        void broadcast() const {
            int result = pthread_cond_broadcast(&cv);
            if (result) MADNESS_EXCEPTION("ConditionalVariable: signalling failed", result);
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

    class Barrier {
        const int nthread;
        volatile bool sense;
        AtomicInt nworking;
        volatile bool* pflags[64];
        
    public:
        Barrier(int nthread)
            : nthread(nthread)
            , sense(true)
        {
            nworking = nthread;
        }
        
        /// Each thread calls this once before first use
        
        /// id should be the thread id (0,..,nthread-1) and pflag a pointer to
        /// thread-local bool (probably in the thread's stack)
        void register_thread(int id, volatile bool* pflag) {
            if (id > 63) MADNESS_EXCEPTION("Barrier : hard dimension failed", id);
            pflags[id] = pflag;
            *pflag=!sense;
        }
        
        /// Each thread calls this with its id (0,..,nthread-1) to enter the barrier
        
        /// The thread last to enter the barrier returns true.  Others return false.
        /// 
        /// All calls to the barrier must use the same value of nthread.
        bool enter(const int id) {
            if (nthread <= 1) {
                return true;
            }
            else {
                if (id > 63) MADNESS_EXCEPTION("Barrier : hard dimension failed", id);
                bool lsense = sense; // Local copy of sense
                bool result = nworking.dec_and_test();
                if (result) {
                    // Reset counter and sense for next entry
                    nworking = nthread;
                    sense = !sense;
                    __asm__ __volatile__("" : : : "memory");

                    // Notify everyone including me
                    for (int i = 0; i < nthread; i++)
                        *(pflags[i]) = lsense;
                } else {
                    volatile bool* myflag = pflags[id]; // Local flag;
                    while (*myflag != lsense) {
                        cpu_relax();
                    }
                }
                return result;
            }
        }
    };
}

#endif
