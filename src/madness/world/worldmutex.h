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
#ifndef MADNESS_WORLD_WORLDMUTEX_H__INCLUDED
#define MADNESS_WORLD_WORLDMUTEX_H__INCLUDED

#include <madness/madness_config.h>
#include <pthread.h>
#include <thread>
#include <cstdio>
#ifdef ON_A_MAC
#if __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 101200

#include <os/lock.h>
#include <type_traits>

typedef std::remove_pointer<os_unfair_lock_t>::type pthread_spinlock_t;

inline void pthread_spin_init(pthread_spinlock_t* p, int /*mode*/) {
    *p = OS_UNFAIR_LOCK_INIT;
}
inline int pthread_spin_trylock(pthread_spinlock_t* p) {
    return !os_unfair_lock_trylock(p);
}
inline int pthread_spin_lock(pthread_spinlock_t* p) {
    os_unfair_lock_lock(p);
    return 0;
}
inline int pthread_spin_unlock(pthread_spinlock_t* p) {
    os_unfair_lock_unlock(p);
    return 0;
}

inline int pthread_spin_destroy(pthread_spinlock_t*) {
    return 0;
}

#else

#include <libkern/OSAtomic.h>
typedef OSSpinLock pthread_spinlock_t;

inline void pthread_spin_init(pthread_spinlock_t* p, int /*mode*/) {
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
inline void pthread_spin_destroy(pthread_spinlock_t* /*p*/) {}
#endif // __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 101200
#endif // ON_A_MAC


#include <madness/world/nodefaults.h>
#include <madness/world/timers.h>
#include <madness/world/atomicint.h>
#include <madness/world/madness_exception.h>

/// \file worldmutex.h
/// \brief Implements Mutex, MutexFair, Spinlock, ConditionVariable
/// \addtogroup mutexes
///@{



namespace madness {

    namespace detail {
        void print_mutex_error(int error_number);
    }

    class MutexWaiter {
    private:
        unsigned int count;

        /// Yield for specified number of microseconds unless dedicated CPU
        /// Blue Gene always has a dedicated hw thread
        void yield(int us) {
#if !defined(HAVE_IBMBGP) && !defined(HAVE_IBMBGQ)
            myusleep(us);
#endif
        }

    public:
        MutexWaiter() : count(0) { }

        void reset() { count = 0; }

        void wait();
    }; // class MutexWaiter


    /// Mutex using pthread mutex operations
    class Mutex {
    private:
        mutable pthread_mutex_t mutex;

        /// Copy constructor is forbidden
        Mutex(const Mutex&);

        /// Assignment is forbidden
        void operator=(const Mutex&);

    public:
        /// Make and initialize a mutex ... initial state is unlocked
        Mutex(int junk=0) // Junk so that can force initializer in mraX.cc
        {
            const int result = pthread_mutex_init(&mutex, 0);
            if (result) MADNESS_EXCEPTION("failed to initialize mutex", result);
        }

        /// Try to acquire the mutex ... return true on success, false on failure
        bool try_lock() const {
            return pthread_mutex_trylock(&mutex)==0;
        }

        /// Acquire the mutex waiting if necessary
        void lock() const {
            const int result = pthread_mutex_lock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: Mutex::lock() failed acquiring mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("Mutex::lock() failed acquiring mutex", result);
            }
        }

        /// Free a mutex owned by this thread
        void unlock() const {
            const int result = pthread_mutex_unlock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: Mutex::unlock() failed releasing mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("Mutex::unlock() failed releasing mutex", result);
            }
        }

        /// Return a pointer to the pthread mutex for use by a condition variable
        pthread_mutex_t* ptr() const {
            return &mutex;
        }

        virtual ~Mutex() {
            pthread_mutex_destroy(&mutex);
        }
    }; // class Mutex

    /// Recursive mutex using pthread mutex operations
    class RecursiveMutex {
    private:
        mutable pthread_mutex_t mutex;

        /// Copy constructor is forbidden
        RecursiveMutex(const RecursiveMutex&);

        /// Assignment is forbidden
        void operator=(const RecursiveMutex&);

    public:
        /// Make and initialize a mutex ... initial state is unlocked
        RecursiveMutex();

        /// Try to acquire the mutex ... return true on success, false on failure
        bool try_lock() const {
            return pthread_mutex_trylock(&mutex)==0;
        }

        /// Acquire the mutex waiting if necessary
        void lock() const {
            int result = pthread_mutex_lock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: RecursiveMutex::lock() failed acquiring mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("RecursiveMutex::lock() failed acquiring mutex", result);
            }
        }

        /// Free a mutex owned by this thread
        void unlock() const {
            int result = pthread_mutex_unlock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: RecursiveMutex::unlock() failed releasing mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("RecursiveMutex::unlock() failed releasing mutex", result);
            }
        }

        /// Return a pointer to the pthread mutex for use by a condition variable
        pthread_mutex_t* ptr() const {
            return &mutex;
        }

        ~RecursiveMutex() {
            pthread_mutex_destroy(&mutex);
        }
    }; // class Mutex


    /// Mutex that is applied/released at start/end of a scope

    /// The mutex must provide lock and unlock methods
    template <class mutexT = Mutex>
    class ScopedMutex {
        const mutexT* mutex;
    public:
        ScopedMutex(const mutexT* m) : mutex(m) { mutex->lock(); }

        ScopedMutex(const mutexT& m) : mutex(&m) { mutex->lock(); }

        virtual ~ScopedMutex() { mutex->unlock(); }
    }; // class ScopedMutex

#ifdef NEVER_SPIN
    typedef Mutex Spinlock;
#else
    /// Spinlock using pthread spinlock operations
    class Spinlock {
    private:
        //mutable pthread_spinlock_t spinlock  __attribute__ ((aligned (64)));
        mutable pthread_spinlock_t spinlock;

        /// Copy constructor is forbidden
        Spinlock(const Spinlock&);

        /// Assignment is forbidden
        void operator=(const Spinlock&);

    public:
        /// Make and initialize a spinlock ... initial state is unlocked
        Spinlock(int junk=0)  // Junk so that can force initializer in mraX.cc
        {
            pthread_spin_init(&spinlock, PTHREAD_PROCESS_PRIVATE);
        }

        /// Try to acquire the spinlock ... return true on success, false on failure
        bool try_lock() const {
            return pthread_spin_trylock(&spinlock)==0;
        }

        /// Acquire the spinlock waiting if necessary
        void lock() const {
            int result = pthread_spin_lock(&spinlock);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: Spinlock::lock() failed acquiring spinlock\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("Spinlock::lock() failed acquiring spinlock", result);
            }
        }

        /// Free a spinlock owned by this thread
        void unlock() const {
            int result = pthread_spin_unlock(&spinlock);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: Spinlock::unlock() failed releasing spinlock\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("Spinlock::unlock() failed releasing spinlock", result);
            }
        }

        virtual ~Spinlock() {
            pthread_spin_destroy(&spinlock);
        }
    }; // class Spinlock
#endif


#define OLDXXX
#ifdef OLDXXX
    // This version uses a spin lock
    class MutexReaderWriter : private Spinlock, private NO_DEFAULTS {
        mutable int nreader; // used to be volatile but is protected by mutex and associated barriers
        mutable bool writeflag; // ditto
    public:
        static const int NOLOCK=0;
        static const int READLOCK=1;
        static const int WRITELOCK=2;

        MutexReaderWriter() : nreader(0), writeflag(false) {}

        bool try_read_lock() const {
            ScopedMutex<Spinlock> protect(this);
            bool gotit = !writeflag;
            if (gotit) ++nreader;
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
            ++nreader;
            writeflag=false;
        }
        virtual ~MutexReaderWriter() {}
    };

#else

    // This version uses AtomicInt and CAS
    class MutexReaderWriter : private NO_DEFAULTS {
        mutable AtomicInt nreader;
        mutable AtomicInt writeflag;
        enum {UNLOCKED, LOCKED};

    public:
        enum lockT {NOLOCK, READLOCK, WRITELOCK};

        MutexReaderWriter() {nreader=0; writeflag=0;}

        bool try_read_lock() const {
            nreader++;
            if (writeflag == UNLOCKED) return true;
            nreader--;
            return false;
        }

        bool try_write_lock() const {
            return (writeflag.compare_and_swap((int) UNLOCKED, (int) LOCKED) == 0);
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
            if (!try_write_lock()) return false;
            if (nreader > 1) {
                write_unlock();
                return false;
            }
            nreader = 0;
            return true;
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
            nreader--;
        }

        void write_unlock() const {
            writeflag = UNLOCKED;
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
            nreader++;
            writeflag = UNLOCKED;
        }
    };
#endif

    /// wait policies supported by ConditionVariable/DQueue/ThreadPool
    enum class WaitPolicy {
      Busy = 1, Yield, Sleep
    };

    /// Scalable and fair condition variable (spins on local value)
    class ConditionVariable : public Spinlock {
    public:
        static const int MAX_NTHREAD = 128;
        mutable int back; // used to be volatile, but is protected by mutex and associated barriers
        mutable int front; // ditto
        mutable volatile bool* fifo[MAX_NTHREAD]; // volatile needed here; Circular buffer of flags

        void set_wait_policy(WaitPolicy p, int us = 0) {
          wait_policy_ = p;
          wait_usleep_ = std::chrono::microseconds(us);
        }

    public:
      ConditionVariable() : back(0), front(0), fifo() { }

        /// You should acquire the mutex before waiting
        void wait() const {
            // We put a pointer to a thread-local variable at the
            // end of the queue and wait for that value to be set,
            // thus generate no memory traffic while waiting.
            volatile bool myturn = false;
            int b = this->back;
            fifo[b] = &myturn;
            this->back = (b+1 < MAX_NTHREAD ? b+1 : 0);

            unlock(); // Release lock before blocking
            switch (this->wait_policy_) {
              case WaitPolicy::Yield:
                while (!myturn) std::this_thread::yield();
              case WaitPolicy::Sleep:
                while (!myturn) std::this_thread::sleep_for(this->wait_usleep_);
              default:
                while (!myturn) cpu_relax();
            }
            lock();
        }

        /// You should acquire the mutex before signalling
        void signal() const {
          int f = this->front;
          if (f == this->back) return;
          *fifo[f] = true;
          int next = (f+1 < MAX_NTHREAD ? f+1 : 0);
          this->front = next;
        }

        /// You should acquire the mutex before broadcasting
        void broadcast() const {
            while (front != back)
                signal();
        }

        virtual ~ConditionVariable() {}

      private:
        WaitPolicy wait_policy_ = WaitPolicy::Busy;
        std::chrono::microseconds wait_usleep_ = std::chrono::microseconds(0);

    };


    /// A scalable and fair mutex (not recursive)

    /// Needs rewriting to use the CV above and do we really
    /// need this if using pthread_mutex .. why not pthread_cv?
    class MutexFair : private Spinlock {
    private:
        static const int MAX_NTHREAD = 128;
        mutable volatile bool* q[MAX_NTHREAD]; // volatile needed
        mutable int n; // volatile not needed due to use of spinlock and associated barriers
        mutable int front;
        mutable int back;

    public:
        MutexFair() : n(0), front(0), back(0) {};

        void lock() const {
            volatile bool myturn = false;
            Spinlock::lock();
            ++n;
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
    class PthreadConditionVariable : private NO_DEFAULTS {
    private:
        mutable pthread_cond_t cv;
        mutable pthread_mutex_t mutex;

    public:
        PthreadConditionVariable() {
            pthread_cond_init(&cv, nullptr);
            pthread_mutex_init(&mutex, 0);
        }

        pthread_mutex_t& get_pthread_mutex() {
            return mutex;
        }

        void lock() const {
            int result = pthread_mutex_lock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: PthreadConditionVariable::lock() failed acquiring mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("PthreadConditionVariable::lock() failed acquiring mutex", result);
            }
        }

        void unlock() const {
            int result = pthread_mutex_unlock(&mutex);
            if (result) {
                fprintf(stderr, "!! MADNESS ERROR: PthreadConditionVariable::unlock() failed releasing mutex\n");
                detail::print_mutex_error(result);
                MADNESS_EXCEPTION("PthreadConditionVariable::unlock() failed releasing mutex", result);
            }
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
    }; // class PthreadConditionVariable

#ifdef USE_SPINLOCKS
    typedef ConditionVariable CONDITION_VARIABLE_TYPE ;
    typedef Spinlock SPINLOCK_TYPE;
    typedef MutexFair SCALABLE_MUTEX_TYPE;
#else
    typedef PthreadConditionVariable CONDITION_VARIABLE_TYPE;
    typedef Mutex SPINLOCK_TYPE;
    typedef Mutex SCALABLE_MUTEX_TYPE;
#endif

    // I THINK THIS IS NO LONGER USED???????????????????????????????????
    class Barrier {
        const int nthread;
        volatile bool sense;
        AtomicInt nworking;
        volatile bool* pflags[128];

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
                    for (int i = 0; i < nthread; ++i)
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
    }; // class Barrier

    namespace detail {
        extern Mutex printmutex;
    }
}

///@}


#endif // MADNESS_WORLD_WORLDMUTEX_H__INCLUDED
