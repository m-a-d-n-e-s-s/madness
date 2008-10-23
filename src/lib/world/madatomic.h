#ifndef MADMADATOMIC_H
#define MADMADATOMIC_H

// A v. good source for these operations is
// http://www.cs.berkeley.edu/~bonachea/gasnet/dist/gasnet_atomicops.h

/// \file madatomic.h
/// \brief Implements atomic operations with necessary memory barriers

/// \full
/// These are defined as macros in case the inliner won't inline
/// assembler.  
///
/// You should only use atomic operations to operate on atomic variables.
///
/// \c MADATOMIC_FENCE ... memory fence that is included, if necessary, in the
///                     macros below (Cray X1 only for now)
///
/// \c MADATOMIC_INT ... the type of an atomic integer (size is platform dependent).
/// 
/// \c MADATOMIC_INT_GET(ptr) ... read an atomic integer
///
/// \c MADATOMIC_INT_SET(ptr) ... write an atomic integer
///
/// \c MADATOMIC_INT_INC(ptr) ... increment an atomic integer
///
/// \c MADATOMIC_INT_DEC_AND_TEST(ptr) ... decrement an atomic integer and return
///                                     true if the result is zero.
///
/// The unfortunate mix of macros and routines means these names
/// are sitting in the global namespace ... probably should turn all
/// of the macros into routines and rely upon inlining to give
/// the performance of macros ... not possible for INITIALIZE.


#if defined(SINGLE_THREADED)

typedef int MADATOMIC_INT;
#define MADATOMIC_FENCE 
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_INC(ptr) (++(*(ptr)))
#define MADATOMIC_INT_GET(ptr) (*(ptr))
#define MADATOMIC_INT_SET(ptr,val) (*(ptr) = val)
#define MADATOMIC_INT_DEC_AND_TEST(ptr) ((--(*(ptr))) == 0)


#elif defined(USE_GLIB_ATOMICS)

error ... not tested in a long time

#include <glib.h>

typedef gint MADATOMIC_INT;

#define MADATOMIC_FENCE 
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_INC(ptr) g_atomic_int_inc(ptr)
#define MADATOMIC_INT_GET(ptr) g_atomic_int_get(ptr)
#define MADATOMIC_INT_SET(ptr,val) g_atomic_int_set(ptr,val)
#define MADATOMIC_INT_DEC_AND_TEST(ptr) g_atomic_int_dec_and_test(ptr)

#elif defined(__GNUC__)

#define GCC_VERSION (__GNUC__*10000 + __GNUC_MINOR__*100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION < 30402
error GCC older than 3.4.3 does not seem to have working atomic operations
#endif


// version 4.* up seems to have switched to ext directory
//#include <bits/atomicity.h>
#include <ext/atomicity.h>
typedef volatile int MADATOMIC_INT;

#define MADATOMIC_FENCE 
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_INC(ptr) (__gnu_cxx::__atomic_add(ptr, 1))
#define MADATOMIC_INT_DEC(ptr) (__gnu_cxx::__atomic_add(ptr,-1))
#define MADATOMIC_INT_GET(ptr) (*(ptr))
#define MADATOMIC_INT_SET(ptr,val) (*(ptr) = val)
#define MADATOMIC_INT_DEC_AND_TEST(ptr) ((__gnu_cxx::__exchange_and_add(ptr,-1)) == 1)

#elif defined(__INTEL_COMPILER)

#if defined(__i386)

// This is just LINUX kernel source provided atomic asm operations 
// ... very generic
#include <asm/atomic.h>
typedef atomic_t MADATOMIC_INT;
#define MADATOMIC_FENCE
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_GET(ptr) atomic_read(ptr)
#define MADATOMIC_INT_SET(ptr,val) atomic_set(ptr,val)
#define MADATOMIC_INT_INC(ptr) atomic_inc(ptr)
#define MADATOMIC_INT_DEC(ptr) atomic_dec(ptr)
#define MADATOMIC_INT_DEC_AND_TEST(ptr) atomic_sub_and_test(1, ptr)

#elif defined(__ia64)
error not yet
#else
error What are we compiling on?
#endif


#elif defined(AIX)
#include <sys/atomic_op.h>
typedef struct { volatile int ctr; } MADATOMIC_INT;
#define MADATOMIC_FENCE
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_GET(ptr) ((ptr)->ctr)
#define MADATOMIC_INT_SET(ptr,val) ((ptr)->ctr = val)
#define MADATOMIC_INT_INC(ptr) (fetch_and_add((atomic_p)&((ptr)->ctr),1))
#define MADATOMIC_INT_DEC(ptr)  (fetch_and_add((atomic_p)&((ptr)->ctr),-1))
#define MADATOMIC_INT_DEC_AND_TEST(ptr)  ((fetch_and_add((atomic_p)&((ptr)->ctr),-1))==1)

#elif defined(_CRAY)

#include <intrinsics.h>
typedef volatile long MADATOMIC_INT;
#define MADATOMIC_FENCE _gsync(0x1)
#define MADATOMIC_INITIALIZE(val) (val)
static inline void MADATOMIC_INT_INC(MADATOMIC_INT *ptr) {
   MADATOMIC_FENCE; 
   _amo_aadd(ptr,1L); 
   MADATOMIC_FENCE;
}
static inline long MADATOMIC_INT_GET(MADATOMIC_INT *ptr) {
  MADATOMIC_FENCE; 
  return *ptr;
}
static inline void MADATOMIC_INT_SET(MADATOMIC_INT *ptr,long val) {
  MADATOMIC_FENCE; 
  *ptr=val;
  MADATOMIC_FENCE;
}
static inline bool MADATOMIC_INT_DEC_AND_TEST(MADATOMIC_INT *ptr) {
  MADATOMIC_FENCE; 
  bool val=(_amo_afadd(ptr,-1L) == 1);
  MADATOMIC_FENCE;
  return val;
}
static inline bool MADATOMIC_INT_COMPARE_AND_SWAP(MADATOMIC_INT *ptr, long cmpval, long swpval) {
  MADATOMIC_FENCE;
  bool val = _amo_acswap(ptr, cmpval, swpval);
  MADATOMIC_FENCE;
  return val;
}
  
#else
error need to define atomic operations or set SINGLE_THREADED

#endif

#endif
