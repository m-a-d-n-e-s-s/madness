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
#include <cstdlib>
#include <type_traits>

namespace madness {

    namespace detail {

        /// Base class for Stack which implements basic memory operations for non-pod objects

        /// \tparam T The data type of the stack
        /// \tparam isPod An auxiliary template parameter to select the
        /// pod/non-pod versions of this class
        template <typename T, bool isPod>
        class StackBase {
        protected:

            /// Destroy a non-pod object

            /// \param ptr A pointer to the object to be destroyed
            /// \throw nothing
            static void destroy(T* ptr) { ptr->~T(); }

            /// Destroy a range of non-pod objects

            /// \param first The beginning of the range to be destroyed
            /// \param last The end of the range to be destroyed
            /// \throw nothing
            static void destroy(T* first, T* last) {
                while(first != last) {
                    --last;
                    destroy(last);
                }
            }

            /// Copy a range of pod objects

            /// \param first The beginning of the range to be copied
            /// \param last The end of the range to be copied
            /// \param dest The beginning of the destination range
            static void copy(const T* first, const T* last, T* dest) {
                std::copy(first, last, dest);
            }

            /// Copy a range of pod objects

            /// \param first The beginning of the range to be moved
            /// \param last The end of the range to be moved
            /// \param dest The pointer to the uninitialized memory range
            static void uninitialized_move(T* first, T* last, T* dest) {
                for (; first != last; ++first, ++dest) {
                    ::new (dest) T(::std::move(*first));
                    destroy(first);
                }
            }

            static void uninitialized_copy(T* first, T* last, T* dest) {
                std::uninitialized_copy(first, last, dest);
            }

        }; // class StackBase


        /// Base class for Stack which implements basic memory operations for pod objects

        /// \tparam T The data type of the stack
        /// \tparam isPod An auxiliary template parameter to select the
        /// pod/non-pod versions of this class
        template <typename T>
        class StackBase<T, true> {
        protected:

            /// Destroy a pod object (no op)

            /// No need to destroy pod's
            /// \throw nothing
            static void destroy(T*) { }

            /// Destroy a range of pod objects (no op)

            /// No need to destroy pod's
            /// \throw nothing
            static void destroy(T*, T*) { }

            /// Copy a range of pod objects

            /// Pod object copy uses memcpy explicitly
            /// \param first The beginning of the range to be copied
            /// \param last The end of the range to be copied
            /// \param dest The beginning of the destination range
            static void copy(const T* first, const T* last, T* dest) {
                if (first != last)
                    memcpy(dest, first, (last - first) * sizeof(T));
            }

            /// Move a range of pod objects to an uninitialized buffer

            /// \param first The beginning of the range to be moved
            /// \param last The end of the range to be moved
            /// \param dest The pointer to the uninitialized memory buffer
            static void uninitialized_move(T* first, T* last, T* dest) {
                // Use simple copy for pods
                copy(first, last, dest);
            }

            /// Move a range of pod objects to an uninitialized buffer

            /// \param first The beginning of the range to be copied
            /// \param last The end of the range to be copied
            /// \param dest The pointer to the uninitialized memory buffer
            static void uninitialized_copy(T* first, T* last, T* dest) {
                copy(first, last, dest);
            }

        }; // class StackBase<T, true>

    } // namespace detail

    /// Dynamically sized Stack with small stack size optimization

    /// This object is a dynamically sized stack that function similarly to a
    /// \c std::vector. It also includes an optimization for small stack sizes
    /// that avoids memory allocations when the stack size is less than or equal
    /// to \c N.
    /// \tparam T The type of data stored in the stack.
    /// \tparam N The fixed size of the stack.
    template <typename T, unsigned int N>
    class Stack : public detail::StackBase<T, std::is_pod<T>::value>{
    private:
        T* data_; ///< Pointer to the stack data
        unsigned int size_; ///< Number of elements on the stack
        unsigned int capacity_; ///< The maximum size, in elements, of the \c data_ buffer
        char buffer_[sizeof(T) * N]; ///< Static buffer for storing a small number of elements

        typedef detail::StackBase<T, std::is_pod<T>::value> StackBase_;

        using StackBase_::destroy;
        using StackBase_::copy;
        using StackBase_::uninitialized_move;
        using StackBase_::uninitialized_copy;

        /// Check if the stack is using the small buffer to store data
        bool is_small() const { return data_ == static_cast<const void*>(buffer_); }

        /// Allocate a raw buffer

        /// Allocate an uninitialized buffer
        /// \param n The size of the new buffer
        T* allocate(unsigned int n) {
            void* const buffer = malloc(n * sizeof(T));
            if(! buffer)
                throw std::bad_alloc();

            return reinterpret_cast<T*>(buffer);
        }


        /// Deallocate memory

        /// Destroy the pointer if it is a dynamically allocated buffer,
        /// otherwise do nothing.
        /// \param ptr The pointer to be destroyed
        void deallocate() { if(! is_small()) free(data_); }

    public:
        /// Construct an empty stack.
        Stack() :
            data_(reinterpret_cast<T*>(buffer_)), size_(0u), capacity_(N)
        { }

        Stack(const Stack<T,N>& other) :
            data_(nullptr), size_(other.size_), capacity_(other.size_)
        {
            if(other.size() <= N) {
                data_ = reinterpret_cast<T*>(buffer_);
                capacity_ = N;
            } else {
                data_ = allocate(other.size_);
            }
            uninitialized_copy(other.data_, other.data_ + other.size_, data_);
        }

        Stack(Stack<T, N>&& other) :
            data_(other.is_small() ? reinterpret_cast<T*>(buffer_) : other.data_),
            size_(other.size_), capacity_(other.capacity_)
        {
            uninitialized_move(other.data_, other.data_ + other.size_, data_);
            other.data_ = reinterpret_cast<T*>(other.buffer_);
            other.size_ = 0u;
            other.capacity_ = N;
        }

        Stack<T,N>& operator=(const Stack<T,N>& other) {
            if(this != &other) { // avoid self assignment
                deallocate();

                // Set the data buffer and size information
                data_ = (other.is_small() ?
                        reinterpret_cast<T*>(buffer_) :
                        other.data_);
                size_ = other.size_;
                capacity_ = other.capacity_;

                // Copy the stack data from other
                copy(other.data_, other.data_ + other.size_, data_);
            }

            return *this;
        }

        /// Move assignment operator

        /// Move the data from \c other to this object. If \c other object is
        /// is using the static buffer, the data is moved to this object's
        /// static buffer. Otherwise, the pointer to the allocated buffer is
        /// moved.
        /// \param other The other stack object to be moved
        /// \note \c other is left in a default constructed state so that it can
        /// continue to be used.
        Stack<T,N>& operator=(Stack<T,N>&& other) {
            if(this != &other) { // avoid self assignment
                deallocate();

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
                }
                size_ = other.size_;


                // Reset member variables for other
                other.data_ = reinterpret_cast<T*>(other.buffer_);
                other.size_ = 0u;
                other.capacity_ = N;
            }

            return *this;
        }

        ~Stack() {
            destroy(data_, data_ + size_);
            deallocate();
        }

        /// Push a new item onto the stack.

        /// Push an item onto the top of the stack. If the stack size is equal
        /// to the capacity, resize the stack (double).
        /// \param[in] value The item to be pushed onto the stack
        void push(const T& value) {
            // Grow the buffer if there is no more free space on the stack
            if(size_ == capacity_) {
                const unsigned int n = (size_ << 1u) + 1u;

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
            new (data_ + size_) T(value);
            ++size_;
        }

        /// Pop an item off of the stack.

        /// \throw MadnessException (via MADNESS_ASSERT) if the stack is empty.
        /// \return The item popped from the stack.
        void pop() {
            MADNESS_ASSERT(size_);
            --size_;
            destroy(data_ + size_);
        }

        /// Get the last item pushed onto the stack

        /// \return A reference to the top of the stack
        /// \throw MadnessException When the stack is empty
        T& top() {
            MADNESS_ASSERT(size_);
            return data_[size_ - 1];
        }

        /// Get the last item pushed onto the stack

        /// \return A const reference to the top of the stack
        /// \throw MadnessException When the stack is empty
        const T& top() const {
            MADNESS_ASSERT(size_);
            return data_[size_ - 1];
        }

        /// Size accessor

        /// \return The number of items pushed to the stack
        unsigned int size() const {
            return size_;
        }

        /// Check if the stack is empty

        /// \return \c true if the size of the stack is 0, otherwise return
        /// \c false
        bool empty() const { return ! size_; }

        /// Empty the stack

        /// This function will destroy the items on the stack (if any), and
        /// set the size to zero.
        void clear() {
            destroy(data_, data_ + size_);
            size_ = 0u;
        }

        /// Empty the stack and free memory

        /// This function will destroy the items on the stack (if any) and
        /// return it to the default constructed state.
        void reset() {
            destroy(data_, data_ + size_);
            deallocate();
            data_ = reinterpret_cast<T*>(buffer_);
            size_ = 0u;
            capacity_ = N;
        }

    }; // class Stack

} // namespace madness

#endif // MADNESS_WORLD_STACK_H__INCLUDED
