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

  $Id: $
*/

#include <world/worldmutex.h>
#include <world/worldtime.h>

/// \file worldmutex.h
/// \brief Implements Mutex, MutexFair, Spinlock, ConditionVariable


namespace madness {

    // MutexWaiter

    void MutexWaiter::wait() {
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

        for (int i=0; i<300; ++i)  cpu_relax();
#else
        const unsigned int nspin  = 1000;    // Spin for 1,000 calls
        const unsigned int nsleep = 100000;  // Sleep 10us for 100,000 calls = 1s
        if (count++ < nspin) return;
        else if (count < nsleep) yield(10);
        else yield(10000);
#endif
    }

    RecursiveMutex::RecursiveMutex() {
        // Create recursive mutex attribute
        pthread_mutexattr_t attr;
        int result = pthread_mutexattr_init(&attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex attribute initialization failed.", result);
        result = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        if (result) MADNESS_EXCEPTION("RecursiveMutex attribute set type failed.", result);

        // Initialize the mutex
        result = pthread_mutex_init(&mutex, &attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex initialization failed.", result);

        // Destroy the mutex attribute
        result = pthread_mutexattr_destroy(&attr);
        if (result) MADNESS_EXCEPTION("RecursiveMutex initialization failed.", result);
    }

#define OLDXXX
#ifdef OLDXXX

    // MutexReaderWriter

    bool MutexReaderWriter::try_read_lock() const {
        ScopedMutex<Spinlock> protect(this);
        bool gotit = !writeflag;
        if (gotit) ++nreader;
        return gotit;
    }

    bool MutexReaderWriter::try_write_lock() const {
        ScopedMutex<Spinlock> protect(this);
        bool gotit = (!writeflag) && (nreader==0);
        if (gotit) writeflag = true;
        return gotit;
    }

    bool MutexReaderWriter::try_lock(int lockmode) const {
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

    bool MutexReaderWriter::try_convert_read_lock_to_write_lock() const {
        ScopedMutex<Spinlock> protect(this);
        bool gotit = (!writeflag) && (nreader==1);
        if (gotit) {
            nreader = 0;
            writeflag = true;
        }
        return gotit;
    }

    void MutexReaderWriter::read_lock() const {
        while (!try_read_lock()) cpu_relax();
    }

    void MutexReaderWriter::write_lock() const {
        while (!try_write_lock()) cpu_relax();
    }

    void MutexReaderWriter::lock(int lockmode) const {
        while (!try_lock(lockmode)) cpu_relax();
    }

    void MutexReaderWriter::read_unlock() const {
        ScopedMutex<Spinlock> protect(this);
        nreader--;
    }

    void MutexReaderWriter::write_unlock() const {
        // Only a single thread should be setting writeflag but
        // probably still need the mutex just to get memory fence?
        ScopedMutex<Spinlock> protect(this);
        writeflag = false;
    }

    void MutexReaderWriter::unlock(int lockmode) const {
        if (lockmode == READLOCK) read_unlock();
        else if (lockmode == WRITELOCK) write_unlock();
        else if (lockmode != NOLOCK) MADNESS_EXCEPTION("MutexReaderWriter: try_lock: invalid lock mode", lockmode);
    }

    void MutexReaderWriter::convert_read_lock_to_write_lock() const {
        while (!try_convert_read_lock_to_write_lock()) cpu_relax();
    }

    void MutexReaderWriter::convert_write_lock_to_read_lock() const {
        ScopedMutex<Spinlock> protect(this);
        ++nreader;
        writeflag=false;
    }

#endif

    // ConditionVariable

    void ConditionVariable::wait() const {
        // We put a pointer to a thread-local variable at the
        // end of the queue and wait for that value to be set,
        // thus generate no memory traffic while waiting.
        volatile bool myturn = false;
        int b = back;
        q[b] = &myturn;
        ++b;
        if (b >= MAX_NTHREAD) back = 0;
        else back = b;

        unlock(); // Release lock before blocking
        while (!myturn) cpu_relax();
        lock();
    }

    void ConditionVariable::signal() const {
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

    void ConditionVariable::broadcast() const {
        while (front != back)
            signal();
    }

    // MutexFair

    void MutexFair::lock() const {
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

    void MutexFair::unlock() const {
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

    bool MutexFair::try_lock() const {
        bool got_lock;

        Spinlock::lock();
        int nn = n;
        got_lock = (nn == 0);
        if (got_lock) n = nn + 1;
        Spinlock::unlock();

        return got_lock;
    }

    // Barrier

    void Barrier::register_thread(int id, volatile bool* pflag) {
        if (id > 63) MADNESS_EXCEPTION("Barrier : hard dimension failed", id);
        pflags[id] = pflag;
        *pflag=!sense;
    }

    bool Barrier::enter(const int id) {
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

} // namespace madness

