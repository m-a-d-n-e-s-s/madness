/*
  This file is part of MADNESS.
  
  Copyright (C) 2007-10 Oak Ridge National Laboratory
  
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

  
#include <cstdio>
#include <complex>

#include <linalg/cblas.h>

typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;

extern "C" void xerbla_(char *message, integer *info, int length) {
    std::fprintf(stderr,
                 " ** On entry to  %6s, parameter number %2ld had an illegal value\n", message,
                 (long) *info);
    throw "XERBLA";
}

/// \file cblas.cc
/// \brief This file provides gemm template BLAS.

extern "C" void dgemm_(const char *opa, const char *opb, const long *m, const long *n, const long *k,
                           const double *alpha, const double *a, const long *lda, const double *b, const long *ldb,
                           const double *beta, double *c, const long *ldc);

extern "C" void sgemm_(const char *opa, const char *opb, const long *m, const long *n, const long *k,
                           const float *alpha, const float *a, const long *lda, const float *b, const long *ldb,
                           const float *beta, float *c, const long *ldc);

extern "C" void zgemm_(const char *opa, const char *opb, const long *m, const long *n, const long *k,
                           const double_complex *alpha,
                           const double_complex *a, const long *lda, const double_complex *b, const long *ldb,
                           const double_complex *beta, double_complex *c, const long *ldc);

extern "C" void cgemm_(const char *opa, const char *opb, const long *m, const long *n, const long *k,
                           const float_complex *alpha,
                           const float_complex *a, const long *lda, const float_complex *b, const long *ldb,
                           const float_complex *beta, float_complex *c, const long *ldc);

namespace madness {


    template <> void gemm<double> (bool transa, bool transb,
                                   integer m, integer n, integer k,
                                   double alpha, const double* a, integer lda,
                                   const double* b, integer ldb,
                                   double beta, double* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        dgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc);
    }

    template <> void gemm<float> (bool transa, bool transb,
                                  integer m, integer n, integer k,
                                  float alpha, const float* a, integer lda,
                                  const float* b, integer ldb,
                                  float beta, float* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        sgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc);
    }

    template <> void gemm<double_complex> (bool transa, bool transb,
                                           integer m, integer n, integer k,
                                           double_complex alpha, const double_complex* a, integer lda,
                                           const double_complex* b, integer ldb,
                                           double_complex beta, double_complex* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        zgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc);
    }

    template <> void gemm<float_complex> (bool transa, bool transb,
                                          integer m, integer n, integer k,
                                          float_complex alpha, const float_complex* a, integer lda,
                                          const float_complex* b, integer ldb,
                                          float_complex beta, float_complex* c, integer ldc) {
        const char *op[] = {"n","t"
                           };
        cgemm_(op[transa], op[transb], &m, &n, &k,
               &alpha, a, &lda, b, &ldb,
               &beta,  c, &ldc);
    }

}


