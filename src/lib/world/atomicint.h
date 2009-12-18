#ifndef MADNESS_WORLD_ATOMICINT_H__INCLUDED
#define MADNESS_WORLD_ATOMICINT_H__INCLUDED

/// \file atomicint.h

/// \brief Implements AtomicInteger

#ifdef HAVE_IBMBGP
#define MADATOMIC_USE_BGP
#elif defined(USE_X86_32_ASM) || defined(USE_X86_64_ASM)
#define MADATOMIC_USE_X86_ASM
#else
#define MADATOMIC_USE_GCC
#endif
    


#ifdef MADATOMIC_USE_GCC
#  ifdef GCC_ATOMICS_IN_BITS
#    include <bits/atomicity.h>
#  else
#    include <ext/atomicity.h>
#  endif
#endif
    
#ifdef MADATOMIC_USE_BGP
#  include <bpcore/bgp_atomic_ops.h>
#endif
    
#ifdef MADATOMIC_USE_AIX
#  include <sys/atomic_op.h>
#endif
    
namespace madness {
    
    /// An integer with atomic set, get, read+inc, read+dec, dec+test operations
    
    /// Only the default constructor is available and IT DOES NOT INITIALIZE THE VARIABLE.
    ///
    /// Conciously modeled after the TBB API to prepare for switching to it.
    class AtomicInt {
    private:

#ifdef MADATOMIC_USE_BGP
        typedef _BGP_Atomic atomic_int;
#else
        typedef volatile int atomic_int;
#endif
        atomic_int value;

        inline int exchange_and_add(int i) {
#ifdef MADATOMIC_USE_GCC
            return __gnu_cxx::__exchange_and_add(&value,i);
#elif defined(MADATOMIC_USE_X86_ASM)
            __asm__ __volatile__("lock; xaddl %0,%1" :"=r"(i) : "m"(value), "0"(i));
            return i;
#elif defined(MADATOMIC_USE_AIX)
            return fetch_and_add(&value,i);
#elif defined(MADATOMIC_USE_BGP)
            return _bgp_fetch_and_add(&value,i);
#else 
#error ... atomic exchange_and_add operator must be implemented for this platform;
#endif
        }

    public:
        /// Returns the value of the counter with fence ensuring subsequent operations are not moved before the load
        operator int() const volatile {
#if defined(MADATOMIC_USE_BGP)
            int result = value.atom;
#else
	    int result = value;
#endif
	    // BARRIER to stop instructions migrating up
            __asm__ __volatile__ ("" : : : "memory");
            return result;
        }
        
        /// Sets the value of the counter with fence ensuring preceding operations are not moved after the store
        int operator=(int other) {
	    // BARRIER to stop instructions migrating down
            __asm__ __volatile__ ("" : : : "memory");
#if defined(MADATOMIC_USE_BGP)
            value.atom = other;
#else
	    value = other;
#endif
	    return other;
        }
        
        /// Sets the value of the counter with fences ensuring operations are not moved either side of the load+store
        AtomicInt& operator=(const AtomicInt& other) {
            *this = int(other);
            return *this;
        }
        
        /// Decrements the counter and returns the original value
        int operator--(int) {
            return exchange_and_add(-1);
        }
        
        /// Increments the counter and returns the original value
        int operator++(int) {
            return exchange_and_add(1);
        }
        
        /// Decrements the counter and returns true if the new value is zero
        bool dec_and_test() {
            return ((*this)-- == 1);
        }
    };
    
}
#endif // MADNESS_WORLD_ATOMICINT_H__INCLUDED
