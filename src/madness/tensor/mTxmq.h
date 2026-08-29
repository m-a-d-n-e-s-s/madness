/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  Copyright (C) 2026 MADNESS developers

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
*/

#ifndef MADNESS_TENSOR_MTXMQ_H__INCLUDED
#define MADNESS_TENSOR_MTXMQ_H__INCLUDED

#include <madness/madness_config.h>
#include <cstddef>
#include <type_traits>

#if !defined(MADNESS_RESTRICT)
#  if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#    define MADNESS_RESTRICT __restrict__
#  else
#    define MADNESS_RESTRICT
#  endif
#endif

namespace madness {

namespace mTxmq_detail {
    constexpr long MAX_DIMI = 400;
    constexpr long MAX_DIMJ = 24;
    constexpr long MAX_DIMK = 24;
}

/// Initialize libxsmm and pre-dispatch JIT kernels for the small matrix domain
/// Calling this is optional; initialization will happen lazily on first call if omitted.
void mTxmq_init();

/// Reference implementation for verification / fallback
template <typename aT, typename bT, typename cT>
void mTxmq_reference(long dimi, long dimj, long dimk,
                     cT* MADNESS_RESTRICT c,
                     const aT* a,
                     const bT* b,
                     long ldb = -1) {
    if (ldb == -1) ldb = dimj;
    if (dimi <= 0 || dimj <= 0) return;
    if (dimk <= 0) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = cT(0);
        return;
    }
    for (long i = 0; i < dimi; ++i, c += dimj, ++a) {
        for (long j = 0; j < dimj; ++j) c[j] = 0.0;
        const aT *aik_ptr = a;
        for (long k = 0; k < dimk; ++k, aik_ptr += dimi) {
            aT aki = *aik_ptr;
            for (long j = 0; j < dimj; ++j) {
                c[j] += aki * b[k * ldb + j];
            }
        }
    }
}

/// Base generic template function matching MADNESS signature:
/// Matrix = Matrix transpose * matrix
/// \code
///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
/// \endcode
template <typename T>
void mTxmq(long dimi, long dimj, long dimk,
           T* MADNESS_RESTRICT c,
           const T* a,
           const T* b,
           long ldb = -1);

#ifdef HAVE_MTXMQ
/// Explicit template specialization for double precision using libxsmm automatic dispatch.
/// Optimized for small matrices: dimi in [1, 400], dimj in [1, 24], dimk in [1, 24].
template <>
void mTxmq(long dimi, long dimj, long dimk,
           double* MADNESS_RESTRICT c,
           const double* a,
           const double* b,
           long ldb);
#endif

/// Overload for mixed precision types
template <typename aT, typename bT, typename cT>
void mTxmq(long dimi, long dimj, long dimk,
           cT* MADNESS_RESTRICT c,
           const aT* a,
           const bT* b,
           long ldb = -1);

/// Print call profile and shape statistics for mTxmq (automatically called at exit)
void print_mtxmq_profile();

/// Reset profiling statistics counters
void reset_mtxmq_profile();

} // namespace madness

#endif // MADNESS_TENSOR_MTXMQ_H__INCLUDED
