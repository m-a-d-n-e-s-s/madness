/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef DQUEUE_H
#define DQUEUE_H

#include <vector>
#include <cstdlib>

/// \file dqueue.h
/// \brief Implements DQueue

namespace madness {
    template <typename T>

    /// A fast and simple doubled-ended queue.  You can push or pop at either end ... that's it.

    /// Since the point is speed, the implementation is a circular buffer rather than
    /// a linked list so as to avoid the new/del overhead.  It will grow as needed,
    /// but presently will not shrink.  The push/pop API is the same as the STL list.
    class DQueue {
        std::size_t sz;         ///< Current capacity
        std::vector<T> buf;     ///< Actual buffer
        int _front;              ///< Index of element at front of buffer
        int _back;               ///< Index of element at back of buffer
        std::size_t n;          ///< Number of elements in the buffer

        void grow() {
            MADNESS_ASSERT(sz == n);
            std::size_t oldsz = sz;
            if (sz < 32768) 
                sz = 65536;
            else if (sz <= 1048576) 
                sz *= 2;
            else 
                sz += 1048576;
            std::vector<T> nbuf(sz);
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
        };

        void sanity_check() const {
            int num = _back - _front + 1;
            if (num < 0) num += sz;
            if (num==int(sz) && n==0) num=0;
            if (num==0 && n==sz) num=sz;
            if (long(n) != num) print("size",sz,"front",_front,"back",_back,"n",n,"num",num);
            MADNESS_ASSERT(long(n) == num);
        };

    public:
        DQueue(std::size_t hint=32768) 
            : sz(hint)
            , buf(sz)
            , _front(sz/2)
            , _back(_front-1)
            , n(0)
        {};

        /// Insert value at front of queue
        void push_front(const T& value) {
            sanity_check();
            if (n == sz) grow();
            _front--;
            if (_front < 0) _front = sz-1;
            buf[_front] = value;
            n++;
        };

        /// Insert element at back of queue
        void push_back(const T& value) {
            sanity_check();
            if (n == sz) grow();
            _back++;
            if (_back >= int(sz)) _back = 0;
            buf[_back] = value;
            n++;
        };

        /// Return a reference to the value at the front of the queue
        T& front() {
            sanity_check();
            MADNESS_ASSERT(n);
            return buf[_front];
        };

        /// Return a const reference to the value at the front of the queue
        const T& front() const {
            sanity_check();
            MADNESS_ASSERT(n);
            return buf[_front];
        };

        /// Pop the front of queue
        void pop_front() {
            sanity_check();
            MADNESS_ASSERT(n);
            n--;
            _front++;
            if (_front >= long(sz)) _front = 0;
        };


        /// Return a reference to the value at the back of the queue
        T& back() {
            sanity_check();
            MADNESS_ASSERT(n);
            return buf[_back];
        };

        /// Return a const reference to the value at the back of the queue
        const T& back() const {
            sanity_check();
            MADNESS_ASSERT(n);
            return buf[_back];
        };

        /// Pop the back of queue
        void pop_back() {
            sanity_check();
            MADNESS_ASSERT(n);
            n--;
            _back--;
            if (_back<0) _back = sz-1;
        };

        std::size_t size() const {
            return n;
        };
        
        bool empty() const {
            return n==0;
        };

        static void self_test() {
            DQueue<int> q;

            // Alternate front/back pushing
            for (int i=0; i<100000; i++) {
                if (i&1) q.push_front(i);
                else q.push_back(i);
            };

            for (int i=99999; i>=0; i--) {
                if (i&1) {
                    MADNESS_ASSERT(i == q.front());
                    q.pop_front();
                }
                else {
                    MADNESS_ASSERT(i == q.back());
                    q.pop_back();
                }
            };
            MADNESS_ASSERT(q.size() == 0);

            // Push from front, pop from back
            for (int i=0; i<100000; i++) {
                q.push_front(i);
            };
            
            for (int i=0; i<100000; i++) {
                MADNESS_ASSERT(i == q.back());
                q.pop_back();
            };
            MADNESS_ASSERT(q.size() == 0);

            // Push from back, pop from front
            for (int i=0; i<100000; i++) {
                q.push_back(i);
            };
            
            for (int i=0; i<100000; i++) {
                MADNESS_ASSERT(i == q.front());
                q.pop_front();
            };
            MADNESS_ASSERT(q.size() == 0);


            // Random pushing and poping
            int npopfront=0, npopback=0, npushfront=0, npushback=0, maxsz=0;

            std::list<int> check;
            for (int i=0; i<2000000; i++) {
                double ran = rand()/double(RAND_MAX);
                if (maxsz < int(q.size())) maxsz = q.size();
                if (q.size() && ran<0.25) {
                    if (ran < 0.125) {
                        check.pop_front();
                        q.pop_front();
                        npopfront++;
                    }
                    else {
                        check.pop_back();
                        q.pop_back();
                        npopback++;
                    }
                }
                else if (ran < 0.5) {
                    q.push_front(i);
                    check.push_front(i);
                    npushfront++;
                }
                else {
                    q.push_back(i);
                    check.push_back(i);
                    npushback++;
                }
                if (q.size()) {
                    MADNESS_ASSERT(check.front() == q.front());
                    MADNESS_ASSERT(check.back() == q.back());
                }
            }

            print("npopfront",npopfront,"npopback",npopback,"npushfront",npushfront,"npushback",npushback,"maxsize",maxsz);

            while (q.size()) {
                MADNESS_ASSERT(check.front() == q.front());
                check.pop_front();
                q.pop_front();
            };
            

            print("DQueue<int>::self_test: OK");

        };
    };
}

#endif
