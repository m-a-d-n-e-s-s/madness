#ifndef MADNESS_WORLD_ATOMICINT_H__INCLUDED
#define MADNESS_WORLD_ATOMICINT_H__INCLUDED

/// \file atomicint.h

/// \brief Implements AtomicInteger

namespace madness {
    
#define MADATOMIC_USE_X86_ASM
    //#define MADATOMIC_USE_GCC
    
#ifdef MADATOMIC_USE_GCC
#  ifdef GCC_ATOMICS_IN_BITS
#    include <bits/atomicity.h>
#  else
#    include <ext/atomicity.h>
#  endif
#endif
    
#ifdef MADATOMIC_USE_AIX
#  include <sys/atomic_op.h>
#endif
    
    /// An integer with atomic set, get, read+inc, read+dec, dec+test operations
    
    /// Only the default constructor is available and IT DOES NOT INITIALIZE THE VARIABLE.
    ///
    /// Conciously modeled after the TBB API to prepare for switching to it.
    class AtomicInt {
    private:
        typedef volatile int atomic_int;
        
        atomic_int value;
        int exchange_and_add(atomic_int* p, int i) {
#ifdef MADATOMIC_USE_GCC
            return __gnu_cxx::__exchange_and_add(p,i);
#elif defined(MADATOMIC_USE_X86_ASM)
            __asm__ __volatile__("lock; xaddl %0,%1" :"=r"(i) : "m"(*p), "0"(i));
            return i;
#elif defined(MADATOMIC_USE_AIX)
            return fetch_and_add(p,i);
#else 
            error ... atomic exchange_and_add operator must be implemented for this platform;
#endif
        }
        
    public:
        /// Returns the value of the counter with fence ensuring subsequent operations are not moved before the load
        operator int() const volatile {
            int result = value;
            __asm__ __volatile__ ("" : : : "memory");
            return result;
        }
        
        /// Sets the value of the counter with fence ensuring preceding operations are not moved after the store
        int operator=(int other) {
            __asm__ __volatile__ ("" : : : "memory");
            value = other;
            return other;
        }
        
        /// Sets the value of the counter with fences ensuring operations are not moved either side of the load+store
        AtomicInt& operator=(const AtomicInt& other) {
            *this = int(other);
            return *this;
        }
        
        /// Decrements the counter and returns the original value
        int operator--(int) {
            return exchange_and_add(&value, -1);
        }
        
        /// Increments the counter and returns the original value
        int operator++(int) {
            return exchange_and_add(&value,  1);
        }
        
        /// Decrements the counter and returns true if the new value is zero
        bool dec_and_test() {
            return ((*this)-- == 1);
        }
    };
    
}
#endif // MADNESS_WORLD_ATOMICINT_H__INCLUDED
