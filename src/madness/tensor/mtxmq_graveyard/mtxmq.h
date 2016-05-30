/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

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

  $Id$
*/
#ifndef MADNESS_TENSOR_MTXMQ_H__INCLUDED
#define MADNESS_TENSOR_MTXMQ_H__INCLUDED

#include <madness/madness_config.h>

#ifdef HAVE_INTEL_MKL
#include <madness/tensor/cblas.h>
#endif

typedef std::complex<double> double_complex;

namespace madness {

#ifdef HAVE_INTEL_MKL
    /// Matrix = Matrix transpose * matrix ... MKL interface version
    
    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    ///
    /// This is the SLOW reference implementation
    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c, const aT* a, const bT* b) {
        const cT one = 1.0;  // alpha in *gemm
        const cT zero = 0.0; // beta  in *gemm
        //std::cout << "IN MKL version mTxmq " << tensor_type_names[TensorTypeData<aT>::id] << " " << tensor_type_names[TensorTypeData<bT>::id] << " " << tensor_type_names[TensorTypeData<cT>::id] << "\n";
        cblas::gemm(cblas::NoTrans,cblas::Trans,dimj,dimi,dimk,one,b,dimj,a,dimi,zero,c,dimj);
    }

#elif

    /// Matrix = Matrix transpose * matrix ... slow reference implementation
    
    /// This routine does \c C=AT*B whereas mTxm does C=C+AT*B.
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c, const aT* a, const bT* b) {
        //std::cout << "IN GENERIC mTxmq " << tensor_type_names[TensorTypeData<aT>::id] << " " << tensor_type_names[TensorTypeData<bT>::id] << " " << tensor_type_names[TensorTypeData<cT>::id] << "\n";
        for (long i=0; i<dimi; ++i,c+=dimj,++a) {
            for (long j=0; j<dimj; ++j) c[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; ++k,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; ++j) {
                    c[j] += aki*b[k*dimj+j];
                }
            }
        }
    }

#endif

    /*
     * mtxm, but with padded buffers.
     *
     * ext_b is the extent of the b array, so shrink() isn't needed.
     */
    template <typename aT, typename bT, typename cT>
    void mTxmq_padding(long dimi, long dimj, long dimk, long ext_b,
               cT* c, const aT* a, const bT* b) {
        const int alignment = 4;
        bool free_b = false;
        long effj = dimj;

        /* Setup a buffer for c if needed */
        cT* c_buf = c;
        if (dimj%alignment) {
            effj = (dimj | 3) + 1;
            c_buf = (cT*)malloc(sizeof(cT)*dimi*effj);
        }

        /* Copy b into a buffer if needed */
        if (ext_b%alignment) {
            free_b = true;
            bT* b_buf = (bT*)malloc(sizeof(bT)*dimk*effj);

            bT* bp = b_buf;
            for (long k=0; k<dimk; k++, bp += effj, b += ext_b)
                memcpy(bp, b, sizeof(bT)*dimj);

            b = b_buf;
            ext_b = effj;
        }

        cT* c_work = c_buf;
        /* mTxm */
        for (long i=0; i<dimi; ++i,c_work+=effj,++a) {
            for (long j=0; j<dimj; ++j) c_work[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; ++k,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; ++j) {
                    c_work[j] += aki*b[k*ext_b+j];
                }
            }
        }

        /* Copy c out if needed */
        if (dimj%alignment) {
            cT* ct = c_buf;
            for (long i=0; i<dimi; i++, ct += effj, c += dimj)
                memcpy(c, ct, sizeof(cT)*dimj);

            free(c_buf);
        }

        /* Free the buffer for b */
        if (free_b) free((bT*)b);
    }
#ifdef HAVE_IBMBGQ
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            double* c, const double* a, const double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const __complex__ double* a, const __complex__ double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const double* a, const __complex__ double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const __complex__ double* a, const double* b);

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                double* c, const double* a, const double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const __complex__ double* a, const __complex__ double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const double* a, const __complex__ double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const __complex__ double* a, const double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }
#elif defined(HAVE_IBMBGP)
    extern void bgpmTxmq(long ni, long nj, long nk, double* restrict c, 
                         const double* a, const double* b);
    extern void bgpmTxmq(long ni, long nj, long nk, double_complex* restrict c, 
                         const double_complex* a, const double_complex* b);
 
    template <>
    inline void mTxmq(long ni, long nj, long nk, double* restrict c, const double* a, const double* b) {
        bgpmTxmq(ni, nj, nk, c, a, b);
    }

    template <>
    inline void mTxmq(long ni, long nj, long nk, double_complex* restrict c, const double_complex* a, const double_complex* b) {
        bgpmTxmq(ni, nj, nk, c, a, b);
    }

// #elif defined(X86_64) && !defined(DISABLE_SSE3)
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double* restrict c, const double* a, const double* b);

//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double_complex* restrict c, const double_complex* a, const double_complex* b);

// #ifndef __INTEL_COMPILER
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double_complex* restrict c, const double_complex* a, const double* b);
// #endif

// #elif defined(X86_32)
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double* restrict c, const double* a, const double* b);
#endif

}

#endif // MADNESS_TENSOR_MTXMQ_H__INCLUDED
