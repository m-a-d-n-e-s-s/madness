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
#include <madness/world/madness_exception.h>
#include <madness/world/thread_specific.h>
#include <cstddef>
#include <type_traits>
#include <complex>
#include <vector>

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

    /// Per-thread grow-on-demand scratch for the complex*real decomposition.
    ///
    /// Backed by a reclaimable thread_specific pool rather than a large stack
    /// array.  The four real buffers scale as dimi*(dimj+dimk), which for 3-D
    /// nonstandard-form blocks (dimi = (2k)^2) already exceeds any sane stack
    /// budget, and an unconditional array sized for the good case enlarges
    /// every frame -- including the calls that end up on the heap path.
    /// Buffers grow monotonically to the largest need seen on that thread and
    /// are never shared across threads; mTxmq_scratch_clear<T>() reclaims them.
    template <typename T> struct Scratch { std::vector<T> buf; };

    template <typename T>
    ::madness::detail::thread_specific<Scratch<T>>& scratch_pool() {
        static ::madness::detail::thread_specific<Scratch<T>> pool;
        return pool;
    }

    template <typename T>
    T* scratch(std::size_t need) {
        Scratch<T>& s = scratch_pool<T>().local();
        if (s.buf.size() < need) s.buf.resize(need);
        return s.buf.data();
    }
}

/// Free every thread's mTxmq scratch buffers.  Call only at a quiescent point.
template <typename T>
void mTxmq_scratch_clear() { mTxmq_detail::scratch_pool<T>().clear(); }

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
    MADNESS_ASSERT(ldb >= dimj);
    if (dimi <= 0 || dimj <= 0) return;
    if (dimk <= 0) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = cT(0);
        return;
    }
    for (long i = 0; i < dimi; ++i, c += dimj, ++a) {
        for (long j = 0; j < dimj; ++j) c[j] = cT(0);
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

/// Explicit template specializations for homogeneous types (implemented in mTxmq.cc)
template <>
void mTxmq(long dimi, long dimj, long dimk,
           double* MADNESS_RESTRICT c,
           const double* a,
           const double* b,
           long ldb);

template <>
void mTxmq(long dimi, long dimj, long dimk,
           float* MADNESS_RESTRICT c,
           const float* a,
           const float* b,
           long ldb);

template <>
void mTxmq(long dimi, long dimj, long dimk,
           std::complex<double>* MADNESS_RESTRICT c,
           const std::complex<double>* a,
           const std::complex<double>* b,
           long ldb);

template <>
void mTxmq(long dimi, long dimj, long dimk,
           std::complex<float>* MADNESS_RESTRICT c,
           const std::complex<float>* a,
           const std::complex<float>* b,
           long ldb);

/// Mixed precision: complex matrix multiplied by real matrix
/// Decomposes into 2 real mTxmq calls over per-thread scratch (see
/// mTxmq_detail::scratch), so steady-state calls allocate nothing.
template <typename T>
inline void mTxmq(long dimi, long dimj, long dimk,
                  std::complex<T>* MADNESS_RESTRICT c,
                  const std::complex<T>* a,
                  const T* b,
                  long ldb = -1) {
    if (ldb == -1) ldb = dimj;
    MADNESS_ASSERT(ldb >= dimj);
    if (dimi <= 0 || dimj <= 0) return;
    if (dimk <= 0) {
        for (long i = 0; i < dimi * dimj; ++i) c[i] = std::complex<T>(0, 0);
        return;
    }

    const long c_sz = dimi * dimj;
    const long a_sz = dimi * dimk;

    T* Ra = mTxmq_detail::scratch<T>(static_cast<std::size_t>(2 * a_sz + 2 * c_sz));
    T* Ia = Ra + a_sz;
    T* Rc = Ia + a_sz;
    T* Ic = Rc + c_sz;

    const T* a_raw = reinterpret_cast<const T*>(a);
    for (long i = 0; i < a_sz; ++i) {
        Ra[i] = a_raw[2 * i];
        Ia[i] = a_raw[2 * i + 1];
    }

    mTxmq(dimi, dimj, dimk, Rc, Ra, b, ldb);
    mTxmq(dimi, dimj, dimk, Ic, Ia, b, ldb);

    T* c_raw = reinterpret_cast<T*>(c);
    for (long i = 0; i < c_sz; ++i) {
        c_raw[2 * i]     = Rc[i];
        c_raw[2 * i + 1] = Ic[i];
    }
}

/// Generic template definitions for any other unspecialized types
template <typename T>
inline void mTxmq(long dimi, long dimj, long dimk,
                  T* MADNESS_RESTRICT c,
                  const T* a,
                  const T* b,
                  long ldb) {
    mTxmq_reference(dimi, dimj, dimk, c, a, b, ldb);
}

template <typename aT, typename bT, typename cT>
inline void mTxmq(long dimi, long dimj, long dimk,
                  cT* MADNESS_RESTRICT c,
                  const aT* a,
                  const bT* b,
                  long ldb = -1) {
    mTxmq_reference(dimi, dimj, dimk, c, a, b, ldb);
}

/// Print call profile and shape statistics for mTxmq (automatically called at exit)
void print_mtxmq_profile();

/// Reset profiling statistics counters
void reset_mtxmq_profile();

} // namespace madness

#endif // MADNESS_TENSOR_MTXMQ_H__INCLUDED
