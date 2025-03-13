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


#ifndef MADNESS_WORLD_DQUEUE_H__INCLUDED
#define MADNESS_WORLD_DQUEUE_H__INCLUDED

// If defined aggregate q insertions to reduce contention on accessing the q
//#define MADNESS_DQ_USE_PREBUF // now in config file

// If defined capture stats on dqueue class --- seems to have small overhead
#define MADNESS_DQ_STATS

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <madness/config.h>
#include <madness/world/thread_info.h>
#include <madness/world/worldmutex.h>
#include <stdint.h>
#include <utility>

/// \file dqueue.h
/// \brief Implements DQueue

namespace madness {

    struct DQStats { // Dilly bar, blizzard, ...
        uint64_t npush_back;    ///< #calls to push_back
        uint64_t npush_front;   ///< #calls to push_front
        uint64_t npop_front;    ///< #calls to pop_front
        uint64_t ngrow;         ///< #calls to grow
        uint64_t nmax;          ///< Lifetime max. entries in the queue

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
        char pad[64]; ///< To put the lock and the data in separate cache lines
        
        // n, sz, buf, _front, _back used to be volatile, but not actually needed since every access
        // happens with the mutex and its implied barriers.
        size_t n __attribute__((aligned(64)));        ///< Number of elements in the buffer
        size_t sz;              ///< Current capacity
        T* buf;        ///< Actual buffer
        int _front;  ///< Index of element at front of buffer
        int _back;    ///< Index of element at back of buffer
        
        DQStats stats;

#ifdef MADNESS_DQ_USE_PREBUF
	static const size_t NPREBUF=MADNESS_DQ_PREBUF_SIZE;
	inline static thread_local T prebuf[NPREBUF] = {T{}}; // relies on this being a singleton class!!!!!!!!!!!!!!!!!!
	inline static thread_local T prebufhi[NPREBUF] = {T{}}; // relies on this being a singleton class!!!!!!!!!!!!!!!!!!
	inline static thread_local size_t ninprebuf = 0, ninprebufhi = 0;
        static auto prebuf_info() { return std::make_tuple(ninprebuf, prebuf, ninprebufhi, prebufhi); }
#endif

        void grow() {
            // ASSUME WE ALREADY HAVE THE MUTEX WHEN IN HERE
#ifdef MADNESS_DQ_STATS
	  ++(stats.ngrow);
#endif
            if (sz != n) MADNESS_EXCEPTION("assertion failure in dqueue::grow", static_cast<int>(sz));
            size_t oldsz = sz;
            if (sz < 32768)
                sz = 65536;
            else if (sz <= 1048576)
                sz *= 2;
            else
                sz += 1048576;
            //volatile T* volatile nbuf = new T[sz];
            T* nbuf = new T[sz];
            int lo = sz/2 - oldsz/2;
            for (int i=_front; i<int(oldsz); ++i,++lo) {
                nbuf[lo] = buf[i];
            }
            if (_front > 0) {
                for (int i=0; i<=_back; ++i,++lo) {
                    nbuf[lo] = buf[i];
                }
            }
            _front = sz/2 - oldsz/2;
            _back = _front + n - 1;
            delete [] buf;
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
            ++nn;
#ifdef MADNESS_DQ_STATS
            if (nn > stats.nmax) stats.nmax = nn;
#endif
            n = nn;

            int b = _back + 1;
            if (b >= int(ss)) b = 0;
            buf[b] = value;
            _back = b;
#ifdef MADNESS_DQ_STATS
            ++(stats.npush_back);
#endif

            signal();
        }

        void push_front_with_lock(const T& value) {
            //sanity_check();

            size_t nn = n;
            size_t ss = sz;
            if (nn == ss) {
                grow();
                ss = sz;
            }
            ++nn;
#ifdef MADNESS_DQ_STATS
            if (nn > stats.nmax) stats.nmax = nn;
#endif
            n = nn;

            int f = _front - 1;
            if (f < 0) f = ss - 1;
            buf[f] = value;
            _front = f;
#ifdef MADNESS_DQ_STATS
            ++(stats.npush_front);
#endif

            //sanity_check();
            signal();
            //broadcast();
        }

	void flush_prebufhi() {
#ifdef MADNESS_DQ_USE_PREBUF
          if (ninprebufhi) {
            // in reverse order of insertion
            for (size_t i=ninprebufhi; i!=0;) {
              push_front_with_lock(prebufhi[--i]);
            }
            ninprebufhi = 0;
          }
#endif
	}

        void flush_prebuflo() {
#ifdef MADNESS_DQ_USE_PREBUF
          if (ninprebuf) {
            for (size_t i=0; i<ninprebuf; i++) push_back_with_lock(prebuf[i]);
            ninprebuf = 0;
          }
#endif
        }

        void flush_prebuf() {
          flush_prebufhi();
          flush_prebuflo();
        }

    public:

        void set_wait_policy(WaitPolicy p, int us = 0) {
#ifdef USE_SPINLOCKS
          ConditionVariable::set_wait_policy(p, us);
#endif
        }

        DQueue(size_t hint=200000) // was 32768
                : n(0)
                , sz(hint>2 ? hint : 2)
                , buf(new T[sz])
                , _front(sz/2)
	        , _back(_front-1) {}

        virtual ~DQueue() {
            delete [] buf;
        }

	void lock_and_flush_prebuf();

        /// Insert value at front of queue
        void push_front(const T& value);

        /// Insert element at back of queue (default is just one copy)
        void push_back(const T& value, int ncopy=1);

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
                ++f;
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
	    flush_prebuf();

            size_t nn = n;

            if (nn==0 && wait) {
                while (n == 0) // !!! Must be n (memory) not nn (local copy)
                    CONDITION_VARIABLE_TYPE::wait(); // might release then reacquire the lock, acts as barrier

                nn = n;
            }

#ifdef MADNESS_DQ_STATS
            ++(stats.npop_front);
#endif
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
                    if (ptr == *(r-1)) {
                        break;
                    }
                    else if (ptr) { // Null pointer indicates stolen task
                        *r++ = ptr;
                        ++f;
                        if (f >= int(thesize)) f = 0;
                        ++retval;
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

        bool empty() const;

        const DQStats& get_stats() const {
            return stats;
        }
    };

    template <typename T>
    void DQueue<T>::lock_and_flush_prebuf() {
#ifdef MADNESS_DQ_USE_PREBUF
        MADNESS_ASSERT(ninprebufhi <= NPREBUF && ninprebuf <= NPREBUF);
        if (ninprebuf+ninprebufhi) {
             madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
             flush_prebuf();
        }
#endif
    }

    template <typename T>
    void DQueue<T>::push_front(const T& value) {
#ifdef MADNESS_DQ_USE_PREBUF
        if (is_madness_thread() && ninprebufhi < NPREBUF) {
             prebufhi[ninprebufhi++] = value;
             return;
        }
        MADNESS_ASSERT(ninprebufhi <= NPREBUF && ninprebuf <= NPREBUF);
#endif
        {
             madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
             push_front_with_lock(value);
             flush_prebufhi();
        }
    }

    template <typename T>
    void DQueue<T>::push_back(const T& value, int ncopy) {
#ifdef MADNESS_DQ_USE_PREBUF
        if (is_madness_thread() && ncopy==1 && ninprebuf < NPREBUF) {
             prebuf[ninprebuf++] = value;
             return;
        }
        MADNESS_ASSERT(ninprebufhi <= NPREBUF && ninprebuf <= NPREBUF);
#endif
        {
             madness::ScopedMutex<CONDITION_VARIABLE_TYPE> obolus(this);
             flush_prebuf();
             //sanity_check();
             while (ncopy--)
                 push_back_with_lock(value);
             //sanity_check();
             //broadcast();
        }
    }

    template <typename T>
    bool DQueue<T>::empty() const {
#ifdef MADNESS_DQ_USE_PREBUF
      return (ninprebuf+ninprebufhi+n)==0; // this is just from the perspective of this thread!!!!!
#else
      return (n==0);
#endif
    }

}  // namespace madness

#endif // MADNESS_WORLD_DQUEUE_H__INCLUDED
