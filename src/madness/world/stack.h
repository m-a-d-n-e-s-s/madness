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

#ifndef MADNESS_WORLD_STACK_H__INCLUDED
#define MADNESS_WORLD_STACK_H__INCLUDED

/**
 \file stack.h
 \brief Implement \c Stack for a fixed-size stack container.
 \ingroup containers
*/

#include <madness/world/madness_exception.h>
#include <array>

namespace madness {

    /// A simple, fixed-size, stack.

    /// \tparam T The type of data stored in the stack.
    /// \tparam N The fixed size of the stack.
    template <typename T, std::size_t N>
    class Stack {
    private:
        std::array<T,N> t; ///< The underlying array storing the stack elements.
        std::size_t n; ///< Number of elements presently stored in the stack.

    public:
        /// Construct an empty stack.
        Stack() : n(0) {}

        /// Push a new item onto the stack.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is full.
        /// \param[in] value The item to be pushed onto the stack.
        void push(const T& value) {
            MADNESS_ASSERT(n < N);
            t[n++] = value;
        }

        /// Pop an item off of the stack.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is empty.
        /// \return The item popped from the stack.
        T& pop() {
            MADNESS_ASSERT(n > 0);
            return t[--n];
        }

        /// Look at the last item pushed onto the stack, but do not pop it off.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is empty.
        /// \return The item at the back of the stack.
        T& front() {
            MADNESS_ASSERT(n > 0);
            return t[n-1];
        }

        /// Look at the last item pushed onto the stack, but do not pop it off.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is empty.
        /// \return The item at the back of the stack.
        T& top() {
            return front();
        }

        /// Access the number of items pushed to the stack.

        /// \return The number of items pushed to the stack.
        std::size_t size() const {
            return n;
        }

        /// Determine if the stack is empty.

        /// \return True if the stack is empty; false otherwise.
        bool empty() const {
            return n==0;
        }

        /// Empty the stack.
        void clear() {
            n = 0;
        }

        /// Empty the stack.
        void reset() {
            clear();
        }

    }; // class Stack

} // namespace madness

#endif // MADNESS_WORLD_STACK_H__INCLUDED
