/* include/madness_config.h.  Generated from madness_config.h.in by configure.  */

#ifndef _madness_include_madnessconfig_h_
#define _madness_include_madnessconfig_h_

/* The host architecture. */
#define HOST_ARCH "x86_64-unknown-linux-gnu"

/* The target architecture. */
#define TARGET_ARCH "x86_64-unknown-linux-gnu"

/* The version number. */
#define MADNESS_VERSION "1.0.0-alpha"

/* The major version number. */
#define MADNESS_MAJOR_VERSION 1

/* The minor version number. */
#define MADNESS_MINOR_VERSION 0

/* The micro version number. */
#define MADNESS_MICRO_VERSION 0

/* have stdint.h ? */
#define HAVE_STDINT_H 1

/* is "restrict" keyword supported ? */
#define restrict __restrict

/* are unqualified static declarations considered? */
#define HAVE_UNQUALIFIED_STATIC_DECL 1

/* have nested template XLC bug? */
/* #undef HAVE_NESTED_TEMPLATE_XLC_BUG */

/* template function instantiation can be used for template argument deduction? */
#define HAVE_TFI_FOR_TEMPLATE_ARGUMENT_DEDUCTION 1

/* have std::abs(long) ? */
#define HAVE_STD_ABS_LONG 1

/* have std::labs? */
/* #undef HAVE_STD_LABS */

/* have boost/shared_array.hpp ? */
#define HAVE_BOOST_SHARED_ARRAY_HPP 1

/* can boost::shared_array<> be compiled and link? */
#define HAVE_BOOST_SHARED_ARRAY_COMPILED 1

/* have MPI library header? */
#define HAVE_MPI_H 1

/* do MPI C++ bindings work? */
#define HAVE_MPI_CXX_BINDINGS 1

/* does MPI library header break C++ compiler? */
/* #undef MPI_H_BREAKS_CXX */

/* if MPI_INIT found -- we have MPI library */
#define HAVE_MPI_LIB 1

/* can use MPI library? */
#if HAVE_MPI_H && HAVE_MPI_LIB
  #define HAVE_MPI 1
#endif

/* have Pthreads library header? */
#define HAVE_PTHREAD_H 1

/* if pthread_join() found -- we have Pthreads library? */
#define HAVE_PTHREAD_LIB 1

/* can use Pthreads library? */
#if HAVE_PTHREAD_H && HAVE_PTHREAD_LIB
  #define HAVE_PTHREADS 1
#endif

#endif 
