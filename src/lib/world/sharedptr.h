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


#ifndef MADNESS_WORLD_SHAREDPTR_H__INCLUDED
#define MADNESS_WORLD_SHAREDPTR_H__INCLUDED

/// \file sharedptr.h
/// \brief Minimal, modified (and renamed) Boost-like SharedPtr & SharedArray

/// This implementation is thread safe.

#include <iostream>
#include <world/worldexc.h>
#include <world/atomicint.h>

#include <unistd.h>

namespace madness {

    /// A SharedCounter counts references to each SharedArray or SharedPtr
    class SharedCounter {
    private:
        madness::AtomicInt count;
        SharedCounter(const SharedCounter& x); // verboten
        void operator=(const SharedCounter& x); // verboten

    public:
        /// Makes a counter with initial value 1
        SharedCounter() {
            count = 1;
        }

        /// Get the count
        int get() {
            return count;
        }

        /// Increment the count
        void inc() {
            count++;
        }

        /// Decrement the count and return true if the decremented value is zero
        bool dec_and_test() {
            return count.dec_and_test();
        }
    };


    namespace detail {
        /// Function to delete arrays for shared pointers
        template <typename T>
        static void del_array(T* t) {
            delete [] t;
        }

        /// Function to delete memory using free()
        template <typename T>
        static void del_free(T* t) {
            free(t);
        }
    }

    template <typename T> class RemoteReference;

    /// A SharedPtr wraps a pointer which is deleted when the reference count goes to zero

    /// The SharedPtr works pretty much like a regular pointer except that it
    /// is reference counted so there is no need to free it.  When the last
    /// reference is destroyed the underlying pointer will be freed.
    ///
    /// By default the pointer is deleted with \c delete, but you can
    /// optionally provide an additional argument on the constructor
    /// for custom deletion.
    /// This is used by SharedArray to invoke \c delete[].
    template <typename T>
    class SharedPtr {
        friend class RemoteReference<T>;
        template <class Q> friend class SharedPtr;
    private:

    protected:
        T* p;                   ///< The pointer being wrapped
        SharedCounter *count;   ///< The counter shared by all references
        bool own;               ///< True if SharedPtr actually owns the pointer ... if not it won't be freed
        void (*deleter)(T*);    ///< Function to invoke to free memory (if null uses delete)

        /// Free the pointer if we own it and it is not null
        void free() {
            if (own && p) {
                //printf("SharedPtr free: own=%d cntptr=%p nref=%d ptr=%p\n", own, count, use_count(), p);
                MADNESS_ASSERT(use_count() == 0);
                if (deleter)
                    deleter(p);
                else
                    delete p;

                delete count;

                p = 0;
                count = 0;
                deleter = 0;
            }
        }

        /// Decrement the reference count, freeing pointer if count becomes zero
        void dec() {
            //if (own && count) print("SharedPtr  dec: own ",own, "cntptr", count, "nref", use_count(),"ptr",p);
            if (own && count && count->dec_and_test()) free();
        }

        /// Same as dec() but works on unowned pointers to avoid a race condition

        /// Race was mark_as_owned(); dec(); mark_as_unowned();
        void dec_not_owned() {
            //if (count) print("SharedPtr  dec_not_owned: own ",own, "cntptr", count, "nref", use_count(),"ptr",p);
            MADNESS_ASSERT(!own);
            if (count && count->dec_and_test()) {
                own  = true;
                free();
            }
        }

        void mark_as_unowned() {
            own = false;
        }

        void mark_as_owned() {
            own = true;
        }

    public:
        /// Default constructor makes an null pointer
        SharedPtr() : p(0), count(0), own(true), deleter(0) {
        }


        /// Wrap a pointer which may be null

        /// The explicit qualifier inhibits very dangerous automatic conversions
        explicit SharedPtr(T* ptr, void (*deleter)(T*)=0) : p(ptr), count(0), own(true), deleter(deleter) {
            if (p) count = new SharedCounter;
        }


        /// Wrap a pointer which may be null, or not owned

        /// If the pointer is not owned it will not be deleted by the destructor.
        /// If the pointer is an array, delete [] will be called.
        explicit SharedPtr(T* ptr, bool own, void (*deleter)(T*)=0) : p(ptr), count(0), own(own), deleter(deleter) {
            if (own && p) count = new SharedCounter;
        }


        /// Copy constructor generates a new reference to the same pointer
        SharedPtr(const SharedPtr<T>& s) : p(s.p), count(s.count), own(s.own), deleter(s.deleter) {
            if (own && count) count->inc();
        }


        /// Copy constructor with static type conversion generates a new reference to the same pointer

        /// !! This is a potentially unsafe conversion of the deleter.  Probably best only
        /// done if deleter=0.
        template <typename Q>
        SharedPtr(const SharedPtr<Q>& s) : p(static_cast<T*>(s.p)), count(s.count), own(s.own),
                deleter((void (*)(T*)) s.deleter) {
            if (own && count) count->inc();
        }


        /// Destructor decrements reference count freeing data only if count is zero
        virtual ~SharedPtr() {
            dec();
        }


        /// Assignment decrements reference count for current pointer and increments new count
        SharedPtr<T>& operator=(const SharedPtr<T>& s) {
            if (this != &s) {
                if (s.own && s.count) s.count->inc();
                this->dec();
                this->p = s.p;
                this->count = s.count;
                this->own = s.own;
                this->deleter = s.deleter;
            }
            return *this;
        }

        /// Returns number of references
        int use_count() const {
            if (count) return count->get();
            else return 0;
        }

        /// Returns the value of the pointer
        T* get() const {
            return p;
        }

        /// Returns true if the SharedPtr owns the pointer
        bool owned() const {
            return own;
        }

        /// Cast of SharedPtr<T> to T* returns the value of the pointer
        operator T*() const {
            return p;
        }

        /// Return pointer+offset
        T* operator+(long offset) const {
            return p+offset;
        }

        /// Return pointer-offset
        T* operator-(long offset) const {
            return p-offset;
        }

        /// Dereferencing SharedPtr<T> returns a reference to pointed value
        T& operator*() const {
            return *p;
        }

        /// Member access via pointer works as expected
        T* operator->() const {
            return p;
        }

        /// Array indexing returns reference to indexed value
        T& operator[](long index) const {
            return p[index];
        }

        /// Boolean value (test for null pointer)
        operator bool() const {
            return p;
        }

        /// Are two pointers equal?
        bool operator==(const SharedPtr<T>& other) const {
            return p == other.p;
        }

        /// Are two pointers not equal?
        bool operator!=(const SharedPtr<T>& other) const {
            return p != other.p;
        }

//         /// Steal an un-owned reference to the pointer

//         /// The returned shared pointer will contain the pointer to
//         /// the shared counter but is marked as not owned AND the
//         /// reference count is NOT incremented.  This enables an
//         /// object containing a SharedPtr to itself to be deleted and
//         /// ensures that destroying the embedded SharedPtr does not
//         /// call the object destructor again.
//         SharedPtr<T> steal() const {
//             SharedPtr<T> r(*this);
//             r.own = false;
//             r.dec();
//             return r;
//         }

        /// This just to permit use as task arguments ... throws if actually invoked
        template <typename Archive>
        void serialize(Archive& ar) {
            MADNESS_EXCEPTION("SharedPtr not serializable", 0);
        }

    };

    /// A SharedArray is just like a SharedPtr except that delete [] is used to free it
    template <class T>
    class SharedArray : public SharedPtr<T> {
    public:
        SharedArray(T* ptr = 0) : SharedPtr<T>(ptr,detail::del_array) {}
        SharedArray(const SharedArray<T>& s) : SharedPtr<T>(s) {}

        /// Assignment decrements reference count for current pointer and increments new count
        SharedArray& operator=(const SharedArray& s) {
            if (this != &s) {
                if (s.own && s.count) s.count->inc();
                this->dec();
                this->p = s.p;
                this->count = s.count;
                this->own = s.own;
                this->deleter = s.deleter;
            }
            return *this;
        }
    };
}
#endif // MADNESS_WORLD_SHAREDPTR_H__INCLUDED
