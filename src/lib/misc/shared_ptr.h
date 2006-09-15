#ifndef SHARED_PTR_H
#define SHARED_PTR_H

/// \file shared_ptr.h
/// \brief Modified (and renamed) Boost-like SharedPtr & SharedArray 

/// This implementation is thread safe, unlike BOOST's version. However,
/// the thread safety is presently NOT turned on by default while 
/// we are implementing the distributed memory code.
/// We also eliminated automatic conversion from shared_ptr<T> to
/// T* which prompted the renaming to SharedPtr.  The logical thing
/// would be to have this wrap Boost's shared_ptr.

#include <iostream>

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
    
    /// A SharedPtr wraps a pointer which is deleted when the reference count goes to zero

    /// The SharedPtr works pretty much like a regular pointer except that it
    /// is reference counted so there is no need to free it.  When the last
    /// reference is destroyed the underlying pointer will be freed.
    template <class T> 
    class SharedPtr {
    protected:
        T* p;                   ///< The pointer being wrapped
        SharedCounter *count;  ///< The counter shared by all references
        bool isarray;           ///< If true use delete [] to free the pointer
        bool own;               ///< True if SharedPtr actually owns the pointer ... if not it won't be freed
        
        /// Decrement the reference count, freeing pointer if count becomes zero
        void dec() {
            if (count) {
                if (count->dec_and_test()) {
                    this->free();
                    delete count;
                }
            }
        };
        
        /// Free the pointer if we own it and it is not null
        void free() {
            if (own && p) {
                //std::cout << "freeing " << (void *) p << " " << isarray << std::endl;
                if (isarray) 
                    delete [] p;    
                else 
                    delete p;
            }
        };
        
    public:
        /// Default constructor makes an null pointer
        SharedPtr() : p(0), count(0), isarray(false), own(true) {
        };
        
        /// Wrap a pointer which may be null
        
        /// The explicit qualifier inhibits VERY DANGEROUS automatic conversions
        explicit SharedPtr(T* ptr) : p(ptr), count(0), isarray(false), own(true) {
            if (p) count = new SharedCounter;
        };
        
        /// Wrap a pointer which may be null, or not owned, or an array
        explicit SharedPtr(T* ptr, bool array, bool own=true) : p(ptr), count(0), isarray(array), own(own) {
            if (own && p) count = new SharedCounter;
        };
        
        /// Copy constructor generates a new reference to the same pointer
        SharedPtr(const SharedPtr<T>& s) : p(s.p), count(s.count), isarray(s.isarray) {
            if (count) {
                count->inc();
//                std::cout << "SharedPtr: copy con " << count->get() << std::endl;
//                std::cout << "copying " << (void *) p << " " << isarray << std::endl;
            }
	    if (count==0 && p!=0) throw "COPYING SHAREDPTR WITH ZERO COUNT BUT NON_ZERO POINTER\n";
        };
        
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
                if (count) count->inc();
            }
//	    std::cout << "SharedPtr: assignment operator " << count->get() << std::endl;
//            std::cout << "assigning " << (void *) p << " " << isarray << std::endl;
            return *this;
        };
        
        /// Returns number of references
        inline int use_count() const {
            if (count) return count->get();
            else return 1;
        };

        /// Returns the value of the pointer
        inline T* get() const {return p;};
       
        // Cast of SharedPtr<T> to T* returns the value of the pointer
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

	    // Less than operator (for sorting)
        // RJH ... N0!  Shared pointer should behave like a pointer.
        // If we are going to do this we should compare the value
        // of the pointers not the data that is pointed to.  If you are using
        // an STL container/algorithm for the sorting, provide your own 
        // comparison class/function at that point.  
	    //inline friend bool operator< (const SharedPtr<T>& t1, const SharedPtr<T>& t2) {
	    //   return (!(*t1 < *t2));
	    //};
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
