/*
  This file is part of MADNESS.

  Copyright (C) 2014 Virginia Tech

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

/**
 \file madness_config.h
 \brief Macros and tools pertaining to the configuration of MADNESS.
 \ingroup libraries
*/

#ifndef MADNESS_MADNESS_CONFIG_H__INCLUDED
#define MADNESS_MADNESS_CONFIG_H__INCLUDED

#include <madness/config.h>
/* undefine what every autoheader package defines to avoid clashes */
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION
#undef VERSION

/* System check */
#if defined(__CRAYXT)

#  define HAVE_CRAYXT 1
#  define MADNESS_HAVE_CRAYXT 1

#elif defined(__CRAYXE)

#  define HAVE_CRAYXE 1
#  define MADNESS_HAVE_CRAYXE 1
#  define X86_64 1
#  define MADNESS_HAVE_X86_64 1

#elif defined(__bgp__)

#  define HAVE_IBMBGP 1
#  define MADNESS_HAVE_POWERPC_32 1

#elif defined(__bgq__)

#  define HAVE_IBMBGQ 1
#  define MADNESS_HAVE_POWERPC_64 1

#endif /* System check */


/* Processor and instruction set checks */
#if defined(__x86_64__) || defined(_M_X64)
   /* x86 64-bit */
#  define X86_64 1
#  define MADNESS_HAVE_X86_64 1

/* http://lists.cs.uiuc.edu/pipermail/cfe-commits/Week-of-Mon-20130819/086386.html */
/* AVX-512 where F=foundational; ER, CD and PF extensions may also be useful some day. */
#  if defined(__AVX512F__)
#    define MADNESS_HAVE_AVX512 1
#  endif

#  if defined(__AVX2__)
#    define MADNESS_HAVE_AVX2 1
#  endif

#  if defined(__AVX__)
#    define MADNESS_HAVE_AVX 1
#  endif

#  if defined(__SSE4_2__)
#    define MADNESS_HAVE_SSE42 1
#  endif

#  if defined(__SSE4_1__)
#    define MADNESS_SSE41 1
#  endif

#  if defined(__SSSE3__)
#    define MADNESS_HAVE_SSSE3 1
#  endif

#  if defined(__SSE3__)
#    define MADNESS_HAVE_SSE3 1
#  endif

/* x86 64-bit always has SSE2 */
#  define MADNESS_HAVE_SSE2 1
#  define MADNESS_HAVE_SSE 1

#  if defined(_M_IX86_FP) /* Defined in MS compiler. 1 = SSE, 2: SSE2 */

#    if _M_IX86_FP == 2
#      define MADNESS_HAVE_SSE2 2
#    elif _M_IX86_FP == 1
#      define MADNESS_HAVE_SSE 1
#    endif

#  endif /* defined(_M_IX86_FP) */


#elif defined(__i386) || defined(_M_IX86)
   /* x86 32-bit */
#  define X86_32
#  define MADNESS_HAVE_X86_32

#  if defined(__SSE2__)
#    define MADNESS_HAVE_SSE2 2
#  endif

#  if defined(__SSE__)
#    define MADNESS_HAVE_SSE 1
#  endif

#endif /* x86 */


#if defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
  /* POWER PC */

#  if defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) || \
      defined(__64BIT__) || defined(_LP64) || defined(__LP64__)
     /* POWER PC 64-bit */
#    define MADNESS_HAVE_POWERPC_64 1

#  else
     /* POWER PC 32-bit */
#    define MADNESS_HAVE_POWERPC_32 1

#  endif

#endif /* POWERPC */

/* ----------- compiler checks -------------------*/
/*
 * - copied from https://github.com/ValeevGroup/tiledarray/blob/master/src/TiledArray/config.h.in
 * - ids taken from CMake
 * - macros are discussed at https://sourceforge.net/p/predef/wiki/Compilers/
*/
#define MADNESS_CXX_COMPILER_ID_GNU 0
#define MADNESS_CXX_COMPILER_ID_Clang 1
#define MADNESS_CXX_COMPILER_ID_AppleClang 2
#define MADNESS_CXX_COMPILER_ID_XLClang 3
#define MADNESS_CXX_COMPILER_ID_Intel 4
#if defined(__INTEL_COMPILER_BUILD_DATE)  /* macros like __ICC and even __INTEL_COMPILER can be affected by command options like -no-icc */
# define MADNESS_CXX_COMPILER_ID MADNESS_CXX_COMPILER_ID_Intel
# define MADNESS_CXX_COMPILER_IS_ICC 1
#endif
#if defined(__clang__) && !defined(MADNESS_CXX_COMPILER_IS_ICC)
# define MADNESS_CXX_COMPILER_IS_CLANG 1
# if defined(__apple_build_version__)
#  define MADNESS_CXX_COMPILER_ID MADNESS_CXX_COMPILER_ID_AppleClang
# elif defined(__ibmxl__)
#  define MADNESS_CXX_COMPILER_ID MADNESS_CXX_COMPILER_ID_XLClang
# else
#  define MADNESS_CXX_COMPILER_ID MADNESS_CXX_COMPILER_ID_Clang
# endif
#endif
#if defined(__GNUG__) && !defined(MADNESS_CXX_COMPILER_IS_ICC) && !defined(MADNESS_CXX_COMPILER_IS_CLANG)
# define MADNESS_CXX_COMPILER_ID MADNESS_CXX_COMPILER_ID_GNU
# define MADNESS_CXX_COMPILER_IS_GCC 1
#endif

/* ----------- preprocessor checks ---------------*/
#define MADNESS_PRAGMA(x) _Pragma(#x)
/* same as MADNESS_PRAGMA(x), but expands x */
#define MADNESS_XPRAGMA(x) MADNESS_PRAGMA(x)
/* "concats" a and b with a space in between */
#define MADNESS_CONCAT(a,b) a b
#ifdef MADNESS_CXX_COMPILER_IS_CLANG
#define MADNESS_PRAGMA_CLANG(x) MADNESS_XPRAGMA( MADNESS_CONCAT(clang,x) )
#else
#define MADNESS_PRAGMA_CLANG(x)
#endif
#ifdef MADNESS_CXX_COMPILER_IS_GCC
#define MADNESS_PRAGMA_GCC(x) MADNESS_XPRAGMA( MADNESS_CONCAT(GCC,x) )
#else
#define MADNESS_PRAGMA_GCC(x)
#endif

/* ----------- end of preprocessor checks ---------*/

#endif // MADNESS_MADNESS_CONFIG_H__INCLUDED
