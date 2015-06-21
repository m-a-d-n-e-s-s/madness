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

#ifndef MADNESS_WORLD_RANGE_H__INCLUDED
#define MADNESS_WORLD_RANGE_H__INCLUDED

#include <type_traits>
#include <iterator>

/**
 \file range.h
 \brief Implement the \c Range class for parallel iteration.
 \ingroup parallel_runtime
*/

namespace madness {

    /// \addtogroup parallel_runtime
    /// @{

    /// Dummy class, a la Intel TBB, used to distinguish splitting constructor.
#ifdef HAVE_INTEL_TBB
    typedef tbb::split Split;
#else
    class Split {};
#endif // HAVE_INTEL_TBB

    /// \brief Range, vaguely a la Intel TBB, to encapsulate a random-access,
    ///    STL-like start and end iterator with chunksize.

    /// \tparam iteratorT The iterator type.
    template <typename iteratorT>
    class Range {
        long n; ///< Number of items to iterator over. \todo Could this be replaced by size_t?
        iteratorT start; ///< First item for iteration.
        iteratorT finish; ///< Last item for iteration (first past the end, conventionally).
        int chunksize; ///< Number of items to give to each thread/process.

    public:
        using iterator = iteratorT; ///< Alias for the iterator type.

        /// Makes the range [start, finish).

        /// The motivated reader should look at the Intel TBB range,
        /// partitioner, split, concepts, etc.
        /// \param[in] start The first item to iterate over.
        /// \param[in] finish The last item for iteration (one past the end).
        /// \param[in] chunk The number of items to give to each thread/process.
        Range(const iterator& start, const iterator& finish, int chunk=1)
            : n(distance(start,finish))
            , start(start)
            , finish(finish)
            , chunksize(chunk)
        {
            if (chunksize < 1) chunksize = 1;
        }

        /// Copy constructor. Cost is O(1).

        /// \todo Can we make this `= default`?
        /// \param[in] r The \c Range to copy.
        Range(const Range& r)
                : n(r.n)
                , start(r.start)
                , finish(r.finish)
                , chunksize(r.chunksize)
        {}

        /// Splits range between new and old (r) objects. Cost is O(1).

        /// \param[in] left The range to be split.
        Range(Range& left, const Split& /*split*/)
                : n(0)
                , start(left.finish)
                , finish(left.finish)
                , chunksize(left.chunksize)
        {
            if (left.n > chunksize) {
                int nleft = (left.n+1)/2;

                start = left.start;
                advance(start,nleft);
                finish = left.finish;
                n = left.n - nleft;

                left.finish = start;
                left.n = nleft;
            }
        }

        /// Returns the number of items in the range (cost is O(1)).

        /// \return The number of items in the range.
        size_t size() const { return n; }

        /// Returns true if `size == 0`.

        /// \return True if `size == 0`.
        bool empty() const { return n==0; }

        /// Access the beginning.

        /// \return Iterator to the first element.
        const iterator& begin() const { return start; }

        /// Access the end.

        /// \return Iterator to the last element (one past the end).
        const iterator& end() const { return finish; }

        /// \brief Return true if this iteration range can be divided; that is,
        ///    there are more items than the chunk size.

        /// \return True if this range can be divided.
        bool is_divisible() const { return n > chunksize; }

        /// Access the chunk size.

        /// \todo Should this return `long`, or `size_t`?
        /// \return The chunk size.
        unsigned int get_chunksize() const { return chunksize; }

    private:
        /// Advance by \c n elements in the range.

        /// This version is for cases where the "iterator type" is integral.
        /// \tparam integralT The integral iterator type.
        /// \tparam distanceT The distance type.
        /// \param[in,out] i The integral iterator.
        /// \param[in] n The number of elements to advance.
        template<typename integralT, typename distanceT>
        inline static typename std::enable_if<std::is_integral<integralT>::value, void>::type
        advance(integralT& i, distanceT n) { i += n; }

        /// Advance by \c n elements in the range.

        /// This version is for cases where the "iterator type" is not integral.
        /// \tparam iterT The non-integral iterator type.
        /// \tparam distanceT The distance type.
        /// \param[in,out] it The iterator.
        /// \param[in] n The number of elements to advance.
        template<typename iterT, typename distanceT>
        inline static typename std::enable_if<!std::is_integral<iterT>::value, void>::type
        advance(iterT& it, distanceT n) { std::advance(it, n); }

        /// Calculate the distance between two iterators.

        /// This version is for cases where the "iterator type" is integral.
        /// \tparam integralT The integral iterator type.
        /// \param[in] first One iterator.
        /// \param[in] last The other iterator.
        /// \return The distance between the first and last iterators.
        template<class integralT>
        inline static typename std::enable_if<std::is_integral<integralT>::value, integralT>::type
        distance(integralT first, integralT last) { return last - first; }

        /// Calculate the distance between two iterators.

        /// This version is for cases where the "iterator type" is not integral.
        /// \tparam iterT The non-integral iterator type.
        /// \param[in] first One iterator.
        /// \param[in] last The other iterator.
        /// \return The distance between the first and last iterators.
        template<class iterT>
        inline static auto
        distance(iterT first, iterT last,
                typename std::enable_if<!std::is_integral<iterT>::value>::type* = nullptr)
            -> decltype(std::distance(first, last))
        { return std::distance(first, last); }
    };

    /// @}

} // namespace madness

#endif // MADNESS_WORLD_RANGE_H__INCLUDED
