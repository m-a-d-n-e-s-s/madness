/*
  This file is part of MADNESS.

  Copyright (C) 2019 Virginia Tech

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


#ifndef MADNESS_LINALG_CBLAS_TYPES_H__INCLUDED
#define MADNESS_LINALG_CBLAS_TYPES_H__INCLUDED

/// \file cblas_types.h
/// \brief Define types used by CBLAS API

#include <madness/madness_config.h>
#include <complex>
#if defined(HAVE_INTEL_MKL) && defined(MKL_DIRECT_CALL)
// calling MKL directly needs its complex types; via the F77 prototypes no conversion is needed
#  include <mkl_types.h>
#endif

namespace madness {
namespace cblas {

    /// The complex types of the BLAS being called. These used to be the macros
    /// blas_complex_float/double, in the LAPACKE style; blaspp declares typedefs of
    /// those names and does not tolerate a pre-existing macro, so any translation
    /// unit that included MADNESS before TiledArray failed to compile.
#if defined(HAVE_INTEL_MKL) && defined(MKL_DIRECT_CALL)
    using complex_float = MKL_Complex8;
    using complex_double = MKL_Complex16;
#else
    using complex_float = std::complex<float>;
    using complex_double = std::complex<double>;
#endif
    static_assert(sizeof(complex_float) == sizeof(std::complex<float>), "the BLAS single-precision complex type must match std::complex<float>");
    static_assert(sizeof(complex_double) == sizeof(std::complex<double>), "the BLAS double-precision complex type must match std::complex<double>");

    /// Matrix operations for BLAS function calls
    typedef enum {
      NoTrans=0,
      Trans=1,
      ConjTrans=2
    }  CBLAS_TRANSPOSE;

    /////////// legalized conversions between C++ and CBLAS types //////////
    template <typename T>
    const complex_float*
    to_cptr(const T* ptr) {
      static_assert(sizeof(T)==sizeof(complex_float), "sizes of complex_float and T given to madness::cblas::to_cptr do not match");
      return reinterpret_cast<const complex_float*>(ptr);
    }
    template <typename T>
    typename std::enable_if<!std::is_const<T>::value, complex_float*>::type
    to_cptr(T* ptr) {
      static_assert(sizeof(T)==sizeof(complex_float), "sizes of complex_float and T given to madness::cblas::to_cptr do not match");
      return reinterpret_cast<complex_float*>(ptr);
    }

    template <typename T>
    const complex_double*
    to_zptr(const T* ptr) {
      static_assert(sizeof(T)==sizeof(complex_double), "sizes of complex_double and T given to madness::cblas::to_zptr do not match");
      return reinterpret_cast<const complex_double*>(ptr);
    }
    template <typename T>
    typename std::enable_if<!std::is_const<T>::value, complex_double*>::type
    to_zptr(T* ptr) {
      static_assert(sizeof(T)==sizeof(complex_double), "sizes of complex_double and T given to madness::cblas::to_zptr do not match");
      return reinterpret_cast<complex_double*>(ptr);
    }

} // namespace cblas
} // namespace madness

#endif // MADNESS_LINALG_CBLAS_TYPES_H__INCLUDED

