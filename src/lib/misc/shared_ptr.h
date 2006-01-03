#ifndef SHARED_PTR_H
#define SHARED_PTR_H

/// \file shared_ptr.h
/// \brief Minimal shared_ptr & shared_array mostly for machines on which we do not have BOOST

/// This implementation is thread safe, unlike BOOST's version. However,
/// the thread safety is presently NOT turned on by default while 
/// we are implementing the distributed memory code.

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
    
    /// A shared_counter counts references to each shared_array or shared_ptr
    class shared_counter {
    private:
        MADATOMIC_INT count;
        shared_counter(const shared_counter& x); // verboten
        void operator=(const shared_counter& x); // verboten
        
    public:
        //#ifdef IBMXLC
        //        shared_counter() {MADATOMIC_INT_SET(&count,1);};
        //#else
        shared_counter() : count(1) {};
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
            return count==0;
        };
    };
    
    /// A shared_ptr wraps a pointer which is deleted when the reference count goes to zero

    /// The shared_ptr works pretty much like a regular pointer except that it
    /// is reference counted so there is no need to free it.  When the last
    /// reference is destroyed the underlying pointer will be freed.
    template <class T> 
    class shared_ptr {
    protected:
        T* p;                   ///< The pointer being wrapped
        shared_counter *count;  ///< The counter shared by all references
        bool isarray;           ///< If true use delete [] to free the pointer
        
        /// Decrement the reference count, freeing pointer if count becomes zero
        void dec() {
            if (count) {
                if (count->dec_and_test()) {
                    this->free();
                    delete count;
                }
            }
        };
        
        /// Free the pointer if it is not null
        void free() {
            if (p) {
                //std::cout << "freeing " << (void *) p << " " << isarray << std::endl;
                if (isarray) 
                    delete [] p;    
                else 
                    delete p;
            }
        };
        
    public:
        /// Default constructor makes an null pointer
        shared_ptr() : p(0), count(0), isarray(false) {
        };
        
        /// Wrap a pointer which may be null
        shared_ptr(T* ptr) : p(ptr), count(0), isarray(false) {
            if (p) count = new shared_counter;
        };
        
        /// Wrap a pointer which may be null marking it as an array
        shared_ptr(T* ptr, bool array) : p(ptr), count(0), isarray(array) {
            if (p) count = new shared_counter;
        };
        
        /// Copy constructor generates a new reference to the same pointer
        shared_ptr(const shared_ptr& s) : p(s.p), count(s.count), isarray(s.isarray) {
            if (count) {
                count->inc();
                //std::cout << "shared_ptr: copy con " << count->get() << std::endl;
            }
        };
        
        /// Destructor decrements reference count freeing data only if count is zero
        virtual ~shared_ptr() {dec();};
        
        /// Assignment decrements reference count for current pointer and increments new count
        shared_ptr& operator=(const shared_ptr& s) {
            if (this != &s) {
                dec();
                p = s.p;
                count = s.count;
                isarray = s.isarray;
                if (count) count->inc();
            }
            return *this;
        };
        
        /// Returns number of references
        inline int use_count() const {
            if (count) return count->get();
            else return 1;
        };

        /// Returns the value of the pointer
        inline T* get() const {return p;};
        
        /// Cast of shared_ptr<T> to T* returns the value of the pointer
        inline operator T*() const {return p;};
        
        /// Return pointer+offset
        inline T* operator+(long offset) const {
            return p+offset;
        };
        
        /// Return pointer-offset
        inline T* operator-(long offset) const {
            return p-offset;
        };
        
        /// Dereferencing shared_ptr<T> returns a reference to pointed value
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
        
    };
    
    /// A shared_array is just like a shared_ptr except that delete [] is used to free it
    template <class T>
    class shared_array : public shared_ptr<T> {
    public:
        shared_array(T* ptr = 0) : shared_ptr<T>(ptr,1) {};
        shared_array(const shared_array<T>& s) : shared_ptr<T>(s) {};
        

        /// Assignment decrements reference count for current pointer and increments new count
        shared_array& operator=(const shared_array& s) {
            if (this != &s) {
                this->dec();
                this->p = s.p;
                this->count = s.count;
                this->isarray = s.isarray;
                if (this->count) this->count->inc();
            }
            return *this;
        };
    };
}
#endif
