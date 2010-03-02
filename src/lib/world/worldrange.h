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
#ifndef MADNESS_WORLD_WORLDRANGE_H__INCLUDED
#define MADNESS_WORLD_WORLDRANGE_H__INCLUDED

#include <iterator>

/// \file worldrange.h
/// \brief Implement Range class for parallel iteration

namespace madness {

    /// Dummy class a la Intel TBB used to distinguish splitting constructor
    class Split {};

    /// Range vaguely a la Intel TBB encapsulates random-access STL-like start and end iterators with chunksize
    template <typename iteratorT>
    class Range {
        long n;
        iteratorT start;
        iteratorT finish;
        int chunksize;
    public:
        typedef iteratorT iterator;

        /// Makes the range [start,finish)

        /// The motivated reader should look at the Intel TBB range,
        /// partitioner, split, concepts, etc..
        ///
        /// Default chunksize is to make 10 tasks per thread to
        /// facilitate dynamic load balancing.
        Range(const iterator& start, const iterator& finish, int chunk=-1)
            : n(std::distance(start,finish))
            , start(start)
            , finish(finish)
            , chunksize(chunk) 
        {
            if (chunksize == -1) chunksize = n / (10*(ThreadPool::size()+1));
            if (chunksize < 1) chunksize = 1;
        }

        /// Copy constructor ... cost is O(1)
        Range(const Range& r)
                : n(r.n)
                , start(r.start)
                , finish(r.finish)
                , chunksize(r.chunksize) 
        {}

        /// Splits range between new and old (r) objects ... cost is O(1)
        Range(Range& left, const Split& split)
                : n(0)
                , start(left.finish)
                , finish(left.finish)
                , chunksize(left.chunksize) 
        {
            if (left.n > chunksize) {
                int nleft = (left.n+1)/2;

                start = left.start;
                std::advance(start,nleft);
                finish = left.finish;
                n = left.n - nleft;

                left.finish = start;
                left.n = nleft;
            }
        }

        /// Returns number of items in the range (cost is O(1))
        size_t size() const {
            return n;
        }

        /// Returns true if size=0
        bool empty() const {
            return n==0;
        }

        const iterator& begin() const {
            return start;
        }

        const iterator& end() const {
            return finish;
        }

        unsigned int get_chunksize() const {
            return chunksize;
        }
    };

    /// Used for iterating over an integer range
    class IntegerIterator {
        long value;

    public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef long value_type;
        typedef long difference_type;
        typedef long* pointer;
        typedef long& reference;

        IntegerIterator(long value = 0) : value(value) {}

        IntegerIterator(const IntegerIterator& other) : value(other.value) {}

        IntegerIterator& operator=(const IntegerIterator& other) {value = other.value; return *this;}

        bool operator==(const IntegerIterator& other) const {return value == other.value;}

        bool operator!=(const IntegerIterator& other) const {return value != other.value;}
        
        IntegerIterator& operator++() {value++; return *this;}

        IntegerIterator operator++(int junk) {IntegerIterator result=*this; value++; return result;}

        IntegerIterator& operator--() {value--; return *this;}

        IntegerIterator operator--(int junk) {IntegerIterator result=*this; value--; return result;}

        const long operator*() const {return value;}

        IntegerIterator operator+(int n) const {IntegerIterator result=*this; result.value += n; return result;}

        IntegerIterator operator-(int n) const {IntegerIterator result=*this; result.value -= n; return result;}

        long operator-(const IntegerIterator& other) const {return value-other.value;}

        IntegerIterator& operator+=(int n) {value += n; return *this;}

        IntegerIterator& operator-=(int n) {value -= n; return *this;}
        
        long operator[](int n) const {return value+n;}
    };
}

#endif // MADNESS_WORLD_WORLDRANGE_H__INCLUDED
