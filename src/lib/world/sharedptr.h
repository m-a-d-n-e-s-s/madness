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

  
#ifndef SHARED_PTR_H
#define SHARED_PTR_H

/// \file sharedptr.h
/// \brief Minimal, modified (and renamed) Boost-like SharedPtr & SharedArray 

/// This implementation is thread safe. However, the thread safety is
/// presently NOT turned on by default while we are implementing the
/// distributed memory code.  In addition to the usual problems with
/// building boost we we using our own shared_ptr since we need to
/// modify it to support the distributed memory interface.

#include <iostream>
#include <world/worldexc.h>

#ifdef MAD_USE_THREADS
#include "madatomic.h"
#else
#define MADATOMIC_INT long
#define MADATOMIC_INT_INC(ptr) (++(*(ptr)))
#define MADATOMIC_INT_GET(ptr) (*(ptr))
#define MADATOMIC_INT_SET(ptr,val) (*(ptr) = val)
#define MADATOMIC_INT_DEC_AND_TEST(ptr) ((--(*(ptr))) == 0)
#endif

namespace madness {
    
    /// A SharedCounter counts references to each SharedArray or SharedPtr
    class SharedCounter {
    private:
        MADATOMIC_INT count;
        SharedCounter(const SharedCounter& x); // verboten
        void operator=(const SharedCounter& x); // verboten
        
    public:
        //#ifdef IBMXLC
        //        SharedCounter() {MADATOMIC_INT_SET(&count,1);};
        //#else
        SharedCounter() : count(1) {};
        //#endif
        
        /// Get the count
        inline int get() {
            //return MADATOMIC_INT_GET(&count);
            return count;
        };
        
        /// Increment the count
        inline void inc() {
            //MADATOMIC_INT_INC(&count);
            count++;
        };
        
        /// Decrement the count and return true if the decremented value is zero
        inline bool dec_and_test() {
            //return MADATOMIC_INT_DEC_AND_TEST(&count);
            count--;
            //	    std::cout << "dec_and_test: count = " << count << std::endl;
            return count==0;
        };
    };
    
    template <class T> class SharedPtr;
    template <typename T> class RemoteReference;

    /// A SharedPtr wraps a pointer which is deleted when the reference count goes to zero
    
    /// The SharedPtr works pretty much like a regular pointer except that it
    /// is reference counted so there is no need to free it.  When the last
    /// reference is destroyed the underlying pointer will be freed.
    template <class T> 
    class SharedPtr {
        friend class RemoteReference<T>;
        template <class Q> friend class SharedPtr;
    protected:
        T* p;                   ///< The pointer being wrapped
        SharedCounter *count;   ///< The counter shared by all references
        bool isarray;           ///< If true use delete [] to free the pointer
        bool own;               ///< True if SharedPtr actually owns the pointer ... if not it won't be freed
        
        /// Free the pointer if we own it and it is not null
        void free() {
            if (own && p) {
                //print("SharedPtr free: own ",own, "cntptr", count, "nref", use_count(),"ptr",p);
                MADNESS_ASSERT(use_count() == 0);
                if (isarray) 
                    delete [] p;    
                else 
                    delete p;
                
                delete count;
                
                p = 0;      
                count = 0;
            }
        };
        
        /// Decrement the reference count, freeing pointer if count becomes zero
        void dec() {
            //print("SharedPtr  dec: own ",own, "cntptr", count, "nref", use_count(),"ptr",p);
            if (own && count && count->dec_and_test()) free();
        };
        
        void mark_as_unowned() {
            own = false;
        };
        
        void mark_as_owned() {
            own = true;
        };
        
    public:
        /// Default constructor makes an null pointer
        SharedPtr() : p(0), count(0), isarray(false), own(true) {
        };
        

        /// Wrap a pointer which may be null
        
        /// The explicit qualifier inhibits very dangerous automatic conversions
        explicit SharedPtr(T* ptr) : p(ptr), count(0), isarray(false), own(true) {
            if (p) count = new SharedCounter;
        };
        
        
        /// Wrap a pointer which may be null, or not owned, or an array
        
        /// If the pointer is not owned it will not be deleted by the destructor.
        /// If the pointer is an array, delete [] will be called.
        explicit SharedPtr(T* ptr, bool array, bool own=true) : p(ptr), count(0), isarray(array), own(own) {
            if (own && p) count = new SharedCounter;
        };
        
        
        /// Copy constructor generates a new reference to the same pointer
        SharedPtr(const SharedPtr<T>& s) : p(s.p), count(s.count), isarray(s.isarray), own(s.own) {
            if (own && count) count->inc();
        };
        

        /// Copy constructor with static type conversion generates a new reference to the same pointer
        template <typename Q>
        SharedPtr(const SharedPtr<Q>& s) : p(static_cast<T*>(s.p)), count(s.count), isarray(s.isarray), own(s.own) {
            if (own && count) count->inc();
        }
        

        /// Destructor decrements reference count freeing data only if count is zero
        virtual ~SharedPtr() {dec();};
        
        
        /// Assignment decrements reference count for current pointer and increments new count
        SharedPtr<T>& operator=(const SharedPtr<T>& s) {
            if (this != &s) {
                dec();
                p = s.p;
                count = s.count;
                isarray = s.isarray;
                own = s.own;
                if (own && count) count->inc();
            }
            return *this;
        };
        
        /// Returns number of references
        inline int use_count() const {
            if (count) return count->get();
            else return 0;
        };
        
        /// Returns the value of the pointer
        inline T* get() const {return p;};
        
        /// Returns true if the SharedPtr owns the pointer
        inline bool owned() const {return own;};
        
        /// Cast of SharedPtr<T> to T* returns the value of the pointer
        inline operator T*() const {return p;};
        
        /// Return pointer+offset
        inline T* operator+(long offset) const {
            return p+offset;
        };
        
        /// Return pointer-offset
        inline T* operator-(long offset) const {
            return p-offset;
        };
        
        /// Dereferencing SharedPtr<T> returns a reference to pointed value
        inline T& operator*() const {
            return *p;
        };
        
        /// Member access via pointer works as expected
        inline T* operator->() const {
            return p;
        };
        
        /// Array indexing returns reference to indexed value
        inline T& operator[](long index) const {
            return p[index];
        };
        
        /// Boolean value (test for null pointer)
        inline operator bool() const {
            return p;
        };
        
        /// Are two pointers equal?
        inline bool operator==(const SharedPtr<T>& other) const {
            return p == other.p;
        };
        
        /// Are two pointers not equal?
        inline bool operator!=(const SharedPtr<T>& other) const {
            return p != other.p;
        };
        
        /// Steal an un-owned reference to the pointer 
        
        /// The returned shared pointer will contain the pointer to
        /// the shared counter but is marked as not owned AND the
        /// reference count is NOT incremented.  This enables an
        /// object containing a SharedPtr to itself to be deleted and
        /// ensures that destroying the embedded SharedPtr does not
        /// call the object destructor again.
        SharedPtr<T> steal() const {
            SharedPtr<T> r(*this);
            r.dec();
            r.own = false;
            return r;
        };

    };
    
    /// A SharedArray is just like a SharedPtr except that delete [] is used to free it
    template <class T>
    class SharedArray : public SharedPtr<T> {
    public:
        SharedArray(T* ptr = 0) : SharedPtr<T>(ptr,1) {};
        SharedArray(const SharedArray<T>& s) : SharedPtr<T>(s) {};
        
        
        /// Assignment decrements reference count for current pointer and increments new count
        SharedArray& operator=(const SharedArray& s) {
            if (this != &s) {
                this->dec();
                this->p = s.p;
                this->count = s.count;
                this->isarray = s.isarray;
                this->own = s.own;
                if (this->count) this->count->inc();
            }
            return *this;
        };
    };
}
#endif
