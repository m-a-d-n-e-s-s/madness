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
/// \c void MADATOMIC_FENCE ... memory fence that is included, if necessary, in the
///                     macros below (Cray X1 only for now)
///
/// \c void MADATOMIC_INT ... the type of an atomic integer (size is platform dependent).
///
/// \c int MADATOMIC_INT_GET(ptr) ... read an atomic integer
///
/// \c void MADATOMIC_INT_SET(ptr) ... write an atomic integer
///
/// \c void MADATOMIC_INT_INC(ptr) ... increment an atomic integer
///
/// \c bool MADATOMIC_INT_DEC_AND_TEST(ptr) ... decrement an atomic integer and return
///                                             true if the result is zero.
///
/// \c int MADATOMIC_INT_READ_AND_INC(ptr) ... atomic read followed by increment
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
#define MADATOMIC_INT_READ_AND_INC(ptr) ((*(ptr))++)


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
#define MADATOMIC_INT_READ_AND_INC(ptr) g_atomic_int_exchange_and_add(ptr,1)

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
#define MADATOMIC_INT_READ_AND_INC(ptr) (__gnu_cxx::__exchange_and_add(ptr,1))

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
#define MADATOMIC_INT_READ_AND_INC(ptr) (error ... no such function in the Linux kernel API)

#elif defined(__ia64)
error not yet
#else
error What are we compiling on?
#endif


#elif defined(AIX)
#include <sys/atomic_op.h>
typedef struct {
    volatile int ctr;
} MADATOMIC_INT;
#define MADATOMIC_FENCE
#define MADATOMIC_INITIALIZE(val) (val)
#define MADATOMIC_INT_GET(ptr) ((ptr)->ctr)
#define MADATOMIC_INT_SET(ptr,val) ((ptr)->ctr = val)
#define MADATOMIC_INT_INC(ptr) (fetch_and_add((atomic_p)&((ptr)->ctr),1))
#define MADATOMIC_INT_DEC(ptr)  (fetch_and_add((atomic_p)&((ptr)->ctr),-1))
#define MADATOMIC_INT_DEC_AND_TEST(ptr)  ((fetch_and_add((atomic_p)&((ptr)->ctr),-1))==1)
#define MADATOMIC_INT_READ_AND_INC() () fetch_and_add((atomic_p)&((ptr)->ctr),1))

#else
error need to define atomic operations

#endif

#endif
