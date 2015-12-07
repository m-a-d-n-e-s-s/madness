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
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <new>
#include <type_traits>

namespace madness {

    namespace detail {

        /// Base class for Stack which implements basic memory operations for non-POD objects.

        /// \tparam T The data type of the stack.
        /// \tparam isPod An auxiliary template parameter to select the
        ///    POD/non-POD versions of this class.
        template <typename T, bool isPod>
        class StackBase {
        protected:

            /// Destroy a non-POD object.

            /// \param[in] ptr A pointer to the object to be destroyed.
            static void destroy(T* ptr) { ptr->~T(); }

            /// Destroy a range of non-POD objects.

            /// \param[in] first The beginning of the range to be destroyed.
            /// \param[in] last The end of the range to be destroyed.
            static void destroy(T* first, T* last) {
                while(first != last) {
                    --last;
                    destroy(last);
                }
            }

            /// Move a range of POD objects.

            /// \param[in] first The beginning of the range to be moved.
            /// \param[in] last The end of the range to be moved.
            /// \param[in] dest Pointer to the uninitialized memory range.
            static void uninitialized_move(T* first, T* last, T* dest) {
                for (; first != last; ++first, ++dest) {
                    ::new (dest) T(std::move(*first));
                    destroy(first);
                }
            }

            /// Copy a range of POD objects.

            /// \param[in] first The beginning of the range to be copied.
            /// \param[in] last The end of the range to be copied.
            /// \param[in] dest Pointer to the uninitialized memory range.
            static void uninitialized_copy(T* first, T* last, T* dest) {
                std::uninitialized_copy(first, last, dest);
            }

        }; // class StackBase


        /// Base class for `Stack` which implements basic memory operations for POD objects.

        /// \tparam T The data type of the stack.
        template <typename T>
        class StackBase<T, true> {
        protected:

            /// Destroy a POD object (no op).

            /// No need to destroy PODs.
            static void destroy(T*) { }

            /// Destroy a range of POD objects (no op).

            /// No need to destroy PODs.
            static void destroy(T*, T*) { }

            /// Move a range of POD objects to an uninitialized buffer.

            /// \param[in] first The beginning of the range to be moved.
            /// \param[in] last The end of the range to be moved.
            /// \param[in] dest Pointer to the uninitialized memory buffer.
            static void uninitialized_move(T* first, T* last, T* dest) {
                // Use simple copy for pods
                if (first != last)
                    std::memcpy(dest, first, (last - first) * sizeof(T));
            }

            /// Copy a range of POD objects to an uninitialized buffer.

            /// \param[in] first The beginning of the range to be copied.
            /// \param[in] last The end of the range to be copied.
            /// \param[in] dest Pointer to the uninitialized memory buffer.
            static void uninitialized_copy(T* first, T* last, T* dest) {
                // Use simple copy for pods
                if (first != last)
                    std::memcpy(dest, first, (last - first) * sizeof(T));
            }

        }; // class StackBase<T, true>

    } // namespace detail

    /// Dynamically sized Stack with small stack size optimization.

    /// This object is a dynamically sized stack that functions similarly to a
    /// \c std::vector. It also includes an optimization for small stack sizes
    /// that avoids memory allocations when the stack size is less than or equal
    /// to \c N.
    /// \tparam T The type of data stored in the stack.
    /// \tparam N The fixed size of the stack.
    template <typename T, unsigned int N>
    class Stack : public detail::StackBase<T, std::is_pod<T>::value> {
    public:
        typedef T value_type; ///< Type of the stack elements.
        typedef T& reference; ///< Element reference type.
        typedef const T& const_reference; ///< Element constant reference type.
        typedef unsigned int size_type; ///< An unsigned integral type.

    private:
        T* data_; ///< Pointer to the stack data.
        size_type size_; ///< Number of elements on the stack.
        size_type capacity_; ///< The maximum size, in elements, of the \c data_ buffer.
        char buffer_[sizeof(T) * N]; ///< Static buffer for storing a small number of elements.

        typedef detail::StackBase<T, std::is_pod<T>::value> StackBase_;

        using StackBase_::destroy;
        using StackBase_::uninitialized_move;
        using StackBase_::uninitialized_copy;

        /// Check if the stack is using the small buffer to store data.

        /// \return True if the small buffer is being used; false otherwise.
        bool is_small() const { return data_ == static_cast<const void*>(buffer_); }

        /// Allocate a raw buffer.

        /// Allocate an uninitialized buffer.
        /// \param[in] n The size of the new buffer.
        /// \return Pointer to the new buffer.
        T* allocate(const size_type n) {
            void* const buffer = std::malloc(n * sizeof(T));
            if(! buffer)
                throw std::bad_alloc();

            return reinterpret_cast<T*>(buffer);
        }

        /// \todo Brief description needed.

        /// \todo Parameter description needed.
        /// \param other Description needed.
        void move(Stack<T,N>& other) {
            // Move the stack data from other to this object
            if(other.is_small()) {
                // Other is using the static buffer space, so the data must
                // be moved to this object's static buffer.
                data_ = reinterpret_cast<T*>(buffer_);
                uninitialized_move(other.data_, other.data_ + other.size_, data_);
                capacity_ = N;
            } else {
                // Other is using an allocated buffer, so move the pointer
                // to the data
                data_ = other.data_;
                capacity_ = other.capacity_;

                other.data_ = reinterpret_cast<T*>(other.buffer_);
                other.capacity_ = N;
            }
            size_ = other.size_;
            other.size_ = 0u;
        }

        /// Deallocate memory.

        /// Destroy the pointer if it is a dynamically allocated buffer;
        /// otherwise do nothing.
        void deallocate() { if(! is_small()) std::free(data_); }

    public:
        /// Construct an empty stack.

        /// The capacity of the stack is \c N.
        Stack() :
            data_(reinterpret_cast<T*>(buffer_)),
            size_(0u), capacity_(N)
        { }

        /// Copy constructor.

        /// If the size of \c other is less than or equal to \c N, then this
        /// object will use the small buffer. Otherwise, it will allocate memory
        /// and copy the data of \c other. The capacity of the object will be
        /// equal to <tt>max(N, other.size())</tt>.
        /// \param[in] other The stack to be copied.
        Stack(const Stack<T,N>& other) {
            if(other.size_ > N) {
                data_ = allocate(other.size_);
                capacity_ = other.size_;
            } else {
                data_ = reinterpret_cast<T*>(buffer_);
                capacity_ = N;
            }
            size_ = other.size_;
            uninitialized_copy(other.data_, other.data_ + other.size_, data_);
        }

        /// Move constructor.

        /// Move the data from \c other to this object.
        /// \param[in] other The original stack.
        Stack(Stack<T, N>&& other) { move(other); }

        /// Assignment operator.

        /// If the size of \c other is less than or equal to \c N, then this
        /// object will use the small buffer. Otherwise, it will allocate memory
        /// and copy the data of \c other. The capacity of the object will be
        /// equal to <tt>max(N, other.size())</tt>.
        /// \param[in] other The stack to be copied.
        Stack<T,N>& operator=(const Stack<T,N>& other) {
            if(this != &other) { // avoid self assignment

                if(capacity_ < other.size_) {
                    // Allocate a larger buffer
                    T* const buffer = allocate(other.size_);
                    uninitialized_copy(other.data_, other.data_ + other.size_, buffer);

                    // Destroy the existing buffer
                    destroy(data_, data_ + size_);
                    deallocate();

                    data_ = buffer;
                    capacity_ = other.size_;
                } else {
                    destroy(data_, data_ + size_);
                    uninitialized_copy(other.data_, other.data_ + other.size_, data_);
                }

                size_ = other.size_;
            }

            return *this;
        }

        /// Move assignment operator.

        /// Move the data from \c other to this object. If \c other object is
        /// is using the static buffer, the data is moved to this object's
        /// static buffer. Otherwise, the pointer to the allocated buffer is
        /// moved.
        /// \param[in] other The other stack object to be moved
        /// \note \c other is left in a default constructed state so that it can
        /// continue to be used.
        Stack<T,N>& operator=(Stack<T,N>&& other) {
            if(this != &other) { // avoid self assignment
                destroy(data_, data_ + size_);
                deallocate();
                move(other);
            }

            return *this;
        }

        /// Destructor.
        ~Stack() {
            destroy(data_, data_ + size_);
            deallocate();
        }

        /// Push a new item onto the stack.

        /// Push an item onto the top of the stack. If the stack size is equal
        /// to the capacity, resize the stack (double).
        /// \param[in] value The item to be pushed onto the stack.
        void push(const_reference value) {
            // Grow the buffer if there is no more free space on the stack
            if(size_ == capacity_) {
                const size_type n = (size_ << 1u) + 1u;

                // Allocate new storage
                T* const new_data = allocate(n);

                // Move data to new vector
                uninitialized_move(data_, data_ + size_, new_data);

                // Deallocate the current data buffer
                deallocate();

                // Update the stack data and capacity
                data_ = new_data;
                capacity_ = n;
            }

            // Add value to the top of the stack
            ::new (data_ + size_) T(value);
            ++size_;
        }

        /// Pop an item off of the stack.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is empty.
        void pop() {
            MADNESS_ASSERT(size_);
            --size_;
            destroy(data_ + size_);
        }

        /// Get the last item pushed onto the stack.

        /// \return A reference to the top of the stack.
        /// \throw MadnessException When the stack is empty
        reference top() {
            MADNESS_ASSERT(size_);
            return data_[size_ - 1];
        }

        /// Get the last item pushed onto the stack.

        /// \return A const reference to the top of the stack.
        /// \throw MadnessException When the stack is empty.
        const_reference top() const {
            MADNESS_ASSERT(size_);
            return data_[size_ - 1];
        }

        /// Size accessor.

        /// \return The number of items pushed to the stack.
        size_type size() const { return size_; }

        /// Capacity accessor.

        /// \return The size of allocated storage capacity.
        size_type capacity() const { return capacity_; }

        /// Check if the stack is empty.

        /// \return True if the size of the stack is 0; otherwise false.
        bool empty() const { return ! size_; }

        /// Empty the stack.

        /// Destroy items on the stack (if any) and set the size to 0.
        void clear() {
            destroy(data_, data_ + size_);
            size_ = 0u;
        }

        /// Empty the stack and free memory.

        /// Destroy items on the stack (if any) and return it to the
        /// default constructed state.
        void reset() {
            destroy(data_, data_ + size_);
            deallocate();
            data_ = reinterpret_cast<T*>(buffer_);
            size_ = 0u;
            capacity_ = N;
        }

        /// Data accessor.

        /// \return A const pointer to the stack data.
        const T* data() const { return data_; }

        /// Data accessor.

        /// \return A pointer to the stack data.
        T* data() { return data_; }

    }; // class Stack

} // namespace madness

#endif // MADNESS_WORLD_STACK_H__INCLUDED
