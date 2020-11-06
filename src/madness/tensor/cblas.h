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
*/


#ifndef MADNESS_LINALG_CBLAS_H__INCLUDED
#define MADNESS_LINALG_CBLAS_H__INCLUDED

/// \file cblas.h
/// \brief Define BLAS like functions


#include <madness/fortran_ctypes.h>
#include <madness/madness_config.h>
#include <madness/world/madness_exception.h>

// MKL direct macros produce a zillion warning messages about unused variables --- disable this warning just in this header
MADNESS_PRAGMA_GCC(diagnostic push)
MADNESS_PRAGMA_GCC(diagnostic ignored "-Wunused-value")
MADNESS_PRAGMA_CLANG(diagnostic push)
MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wunused-value")

#if defined(FORTRAN_LINKAGE_LC) || (defined(HAVE_INTEL_MKL) && defined(MKL_DIRECT_CALL))

#   define F77_SGEMM sgemm
#   define F77_DGEMM dgemm
#   define F77_CGEMM cgemm
#   define F77_ZGEMM zgemm
#ifdef HAVE_INTEL_MKL
#   define F77_SCGEMM scgemm
#   define F77_DZGEMM dzgemm
#endif
#   define F77_SGEMV sgemv
#   define F77_DGEMV dgemv
#   define F77_CGEMV cgemv
#   define F77_ZGEMV zgemv
#   define F77_SSCAL sscal
#   define F77_DSCAL dscal
#   define F77_CSCAL cscal
#   define F77_ZSCAL zscal
#   define F77_CSSCAL csscal
#   define F77_ZDSCAL zdscal
#   define F77_SDOT sdot
#   define F77_DDOT ddot
#   define F77_CDOTU cdotu
#   define F77_ZDOTU zdotu
#   define F77_SAXPY saxpy
#   define F77_DAXPY daxpy
#   define F77_CAXPY caxpy
#   define F77_ZAXPY zaxpy

#elif defined(FORTRAN_LINKAGE_LCU)

#   define F77_SGEMM sgemm_
#   define F77_DGEMM dgemm_
#   define F77_CGEMM cgemm_
#   define F77_ZGEMM zgemm_
#ifdef HAVE_INTEL_MKL
#  define F77_SCGEMM scgemm_
#  define F77_DZGEMM dzgemm_
#endif
#   define F77_SGEMV sgemv_
#   define F77_DGEMV dgemv_
#   define F77_CGEMV cgemv_
#   define F77_ZGEMV zgemv_
#   define F77_SSCAL sscal_
#   define F77_DSCAL dscal_
#   define F77_CSCAL cscal_
#   define F77_ZSCAL zscal_
#   define F77_CSSCAL csscal_
#   define F77_ZDSCAL zdscal_
#   define F77_SDOT sdot_
#   define F77_DDOT ddot_
#   define F77_CDOTU cdotu_
#   define F77_ZDOTU zdotu_
#   define F77_SAXPY saxpy_
#   define F77_DAXPY daxpy_
#   define F77_CAXPY caxpy_
#   define F77_ZAXPY zaxpy_

#elif defined(FORTRAN_LINKAGE_LCUU)

#   define F77_SGEMM  sgemm__
#   define F77_DGEMM  dgemm__
#   define F77_CGEMM  cgemm__
#   define F77_ZGEMM  zgemm__
#ifdef HAVE_INTEL_MKL
#   define F77_SCGEMM scgemm__
#   define F77_DZGEMM dzgemm__
#endif
#   define F77_SGEMV  sgemv__
#   define F77_DGEMV  dgemv__
#   define F77_CGEMV  cgemv__
#   define F77_ZGEMV  zgemv__
#   define F77_SSCAL  sscal__
#   define F77_DSCAL  dscal__
#   define F77_CSCAL  cscal__
#   define F77_ZSCAL  zscal__
#   define F77_CSSCAL csscal__
#   define F77_ZDSCAL zdscal__
#   define F77_SDOT   sdot__
#   define F77_DDOT   ddot__
#   define F77_CDOTU  cdotu__
#   define F77_ZDOTU  zdotu__
#   define F77_SAXPY  saxpy__
#   define F77_DAXPY  daxpy__
#   define F77_CAXPY  caxpy__
#   define F77_ZAXPY  zaxpy__

#elif defined(FORTRAN_LINKAGE_UC)

#   define F77_SGEMM  SGEMM
#   define F77_DGEMM  DGEMM
#   define F77_CGEMM  CGEMM
#   define F77_ZGEMM  ZGEMM
#ifdef HAVE_INTEL_MKL
#   define F77_SCGEMM SCGEMM
#   define F77_DZGEMM DZGEMM
#endif
#   define F77_SGEMV  SGEMV
#   define F77_DGEMV  DGEMV
#   define F77_CGEMV  CGEMV
#   define F77_ZGEMV  ZGEMV
#   define F77_SSCAL  SSCAL
#   define F77_DSCAL  DSCAL
#   define F77_CSCAL  CSCAL
#   define F77_ZSCAL  ZSCAL
#   define F77_CSSCAL CSSCAL
#   define F77_ZDSCAL ZDSCAL
#   define F77_SDOT   SDOTU
#   define F77_DDOT   DDOTU
#   define F77_CDOTU  CDOTU
#   define F77_ZDOTU  ZDOTU
#   define F77_SAXPY  SAXPY
#   define F77_DAXPY  DAXPY
#   define F77_CAXPY  CAXPY
#   define F77_ZAXPY  ZAXPY

#elif defined(FORTRAN_LINKAGE_UCU)

#   define F77_SGEMM  SGEMM_
#   define F77_DGEMM  DGEMM_
#   define F77_CGEMM  CGEMM_
#   define F77_ZGEMM  ZGEMM_
#ifdef HAVE_INTEL_MKL
#   define F77_SCGEMM SCGEMM_
#   define F77_DZGEMM DZGEMM_
#endif
#   define F77_SGEMV  SGEMV_
#   define F77_DGEMV  DGEMV_
#   define F77_CGEMV  CGEMV_
#   define F77_ZGEMV  ZGEMV_
#   define F77_SSCAL  SSCAL_
#   define F77_DSCAL  DSCAL_
#   define F77_CSCAL  CSCAL_
#   define F77_ZSCAL  ZSCAL_
#   define F77_CSSCAL CSSCAL_
#   define F77_ZDSCAL ZDSCAL_
#   define F77_SDOT   SDOT_
#   define F77_DDOT   DDOTSUB_
#   define F77_CDOTU  CDOTU_
#   define F77_ZDOTU  ZDOTU_
#   define F77_SAXPY  SAXPY_
#   define F77_DAXPY  DAXPY_
#   define F77_CAXPY  CAXPY_
#   define F77_ZAXPY  ZAXPY_

#else
// If detected another convention complain loudly.
#   error "cblas.h does not support the current Fortran symbol convention -- please, edit and check in the changes."
#endif

// process BLAS parts that are not directly callable in MKL
#if defined(FORTRAN_LINKAGE_LC)
#   define F77_SGER sger
#   define F77_DGER dger
#   define F77_CGER cger
#   define F77_ZGER zger
#elif defined(FORTRAN_LINKAGE_LCU)
#   define F77_SGER sger_
#   define F77_DGER dger_
#   define F77_CGER cger_
#   define F77_ZGER zger_
#elif defined(FORTRAN_LINKAGE_LCUU)
#   define F77_SGER   sger__
#   define F77_DGER   dger__
#   define F77_CGER   cger__
#   define F77_ZGER   zger__
#elif defined(FORTRAN_LINKAGE_UC)
#   define F77_SGER   SGER
#   define F77_DGER   DGER
#   define F77_CGER   CGER
#   define F77_ZGER   ZGER
#elif defined(FORTRAN_LINKAGE_UCU)
#   define F77_SGER   SGER_
#   define F77_DGER   DGER_
#   define F77_CGER   CGER_
#   define F77_ZGER   ZGER_
#else
// If detected another convention complain loudly.
#   error "cblas.h does not support the current Fortran symbol convention -- please, edit and check in the changes."
#endif

extern "C" {

// BLAS _GER declarations, not directly callable via MKL
void F77_SGER(const integer *, const integer *, const float *, const float *,
              const integer *, const float *, const integer *, float *,
              const integer *);
void F77_DGER(const integer *, const integer *, const double *, const double *,
              const integer *, const double *, const integer *, double *,
              const integer *);
void F77_CGER(const integer *, const integer *, const complex_real4 *,
              const complex_real4 *, const integer *, const complex_real4 *,
              const integer *, complex_real4 *, const integer *);
void F77_ZGER(const integer *, const integer *, const complex_real8 *,
              const complex_real8 *, const integer *, const complex_real8 *,
              const integer *, complex_real8 *, const integer *);
}

#ifndef MKL_DIRECT_CALL

extern "C" {

    // BLAS _GEMM declarations
    void F77_SGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const float*, const float*, const integer*,
            const float*, const integer*, const float*, float*, const integer*);
    void F77_DGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const double*, const double*, const integer*,
            const double*, const integer*, const double*, double*, const integer*);
    void F77_CGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const complex_real4*, const complex_real4*,
            const integer*, const complex_real4*, const integer*,
            const complex_real4*, complex_real4*, const integer*);
    void F77_ZGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const complex_real8*, const complex_real8*,
            const integer*, const complex_real8*, const integer*,
            const complex_real8*, complex_real8*, const integer*);

#ifdef HAVE_INTEL_MKL
    void F77_SCGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const complex_real4*, const real4*,
            const integer*, const complex_real4*, const integer*,
            const complex_real4*, complex_real4*, const integer*);
    void F77_DZGEMM(const char*, const char*, const integer*, const integer*,
            const integer*, const complex_real8*, const real8*,
            const integer*, const complex_real8*, const integer*,
            const complex_real8*, complex_real8*, const integer*);
#endif

    // BLAS _GEMV declarations
    void F77_SGEMV(const char*, const integer*, const integer*, const float*,
            const float*, const integer*, const float*, const integer*,
            const float*, float*, const integer*);
    void F77_DGEMV(const char*, const integer*, const integer*, const double*,
            const double*, const integer*, const double*, const integer*,
            const double*, double*, const integer*);
    void F77_CGEMV(const char*, const integer*, const integer*, const complex_real4*,
            const complex_real4*, const integer*, const complex_real4*,
            const integer*, const complex_real4*, complex_real4*, const integer*);
    void F77_ZGEMV(const char*, const integer*, const integer*, const complex_real8*,
            const complex_real8*, const integer*, const complex_real8*,
            const integer*, const complex_real8*, complex_real8*, const integer*);

    // BLAS _SCAL declarations
    void F77_SSCAL(const integer*, const float*, float*, const integer*);
    void F77_DSCAL(const integer*, const double*, double*, const integer*);
    void F77_CSCAL(const integer*, const complex_real4*, complex_real4*, const integer*);
    void F77_CSSCAL(const integer*, const float*, complex_real4*, const integer*);
    void F77_ZSCAL(const integer*, const complex_real8*, complex_real8*, const integer*);
    void F77_ZDSCAL(const integer*, const double*, complex_real8*, const integer*);

    // BLAS _DOT declarations
    float F77_SDOT(const integer*, const float*, const integer*, const float*,
            const integer*);
    double F77_DDOT(const integer*, const double *, const integer*,
            const double *, const integer*);
    void F77_CDOTU(complex_real4*, const integer*, const complex_real4*, const integer*,
            const complex_real4*, const integer*);
    void F77_ZDOTU(complex_real8*, const integer*, const complex_real8*, const integer*,
            const complex_real8*, const integer*);
    //
    // BLAS _AXPY declarations (INTEGER n, NUMERICAL alpha, NUMERICAL x, INTEGER incx, NUMERICAL y, INTEGER incy )
    void F77_SAXPY(const integer*, const float*, const float*, const integer*,
            float*, const integer*);
    void F77_DAXPY(const integer*, const double*, const double*, const integer*,
            double*, const integer*);
    void F77_CAXPY(const integer*, const complex_real4*, const complex_real4*,
            const integer*, complex_real4*, const integer*);
    void F77_ZAXPY(const integer*, const complex_real8*, const complex_real8*,
            const integer*, complex_real8*, const integer*);
}
#else

# include <mkl.h>

#endif // !defined(MKL_DIRECT_CALL)

// some BLAS libraries use custom complex types in their interface, so need to include their definitions here
#include <madness/tensor/cblas_types.h>

namespace madness {
namespace cblas {

    /// Multiplies a matrix by a vector

    /// \f[
    /// \mathbf{C} \leftarrow \alpha \mathbf{A}^{\mathrm{OpA}} \mathbf{B}^{\mathrm{OpB}} + \beta \mathbf{C}
    /// \f]
    /// \param OpA Operation to be applied to matrix \f$ \mathbf{A} \f$
    /// \param OpB Operation to be applied to matrix \f$ \mathbf{B} \f$
    /// \param m Rows in matrix \f$ \mathbf{C} \f$
    /// \param n Columns in matrix \f$ \mathbf{C} \f$
    /// \param k Inner dimension size for matrices \f$ \mathbf{A} \f$ and \f$ \mathbf{B} \f$
    /// \param alpha Scaling factor applied to \f$ \mathbf{A} \f$ \c * \f$ \mathbf{B} \f$
    /// \param a Pointer to matrix \f$ \mathbf{A} \f$
    /// \param lda The size of the leading-order dimension of matrix \f$ \mathbf{A} \f$
    /// \param b Pointer to matrix \f$ \mathbf{A} \f$
    /// \param ldb The size of the leading-order dimension of matrix \f$ \mathbf{B} \f$
    /// \param beta Scaling factor for matrix \f$ \mathbf{C} \f$
    /// \param c Pointer to matrix \f$ \mathbf{C} \f$
    /// \param ldc The size of the leading-order dimension of matrix \f$ \mathbf{C} \f$
    ///@{
    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
            const integer m, const integer n, const integer k, const float alpha,
            const float* a, const integer lda, const float* b, const integer ldb,
            const float beta, float* c, const integer ldc)
    {
        MADNESS_ASSERT(OpA != ConjTrans);
        MADNESS_ASSERT(OpB != ConjTrans);
        static const char *op[] = { "n","t" };
        F77_SGEMM(op[OpA], op[OpB], &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }

    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
            const integer m, const integer n, const integer k, const double alpha,
            const double* a, const integer lda, const double* b, const integer ldb,
            const double beta, double* c, const integer ldc) {
        MADNESS_ASSERT(OpA != ConjTrans);
        MADNESS_ASSERT(OpB != ConjTrans);
        static const char *op[] = { "n","t" };
        F77_DGEMM(op[OpA], op[OpB], &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }

    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
            const integer m, const integer n, const integer k,
            const complex_real4 alpha, const complex_real4* a, const integer lda,
            const complex_real4* b, const integer ldb, const complex_real4 beta,
            complex_real4* c, const integer ldc) {
      static const char *op[] = {"n", "t", "c"};
      F77_CGEMM(op[OpA], op[OpB], &m, &n, &k, cblas::to_cptr(&alpha),
                cblas::to_cptr(a), &lda, cblas::to_cptr(b), &ldb,
                cblas::to_cptr(&beta), cblas::to_cptr(c), &ldc);
    }

    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
            const integer m, const integer n, const integer k,
            const complex_real8 alpha, const complex_real8* a, const integer lda,
            const complex_real8* b, const integer ldb, const complex_real8 beta,
            complex_real8* c, const integer ldc) {
      static const char *op[] = {"n", "t", "c"};
      F77_ZGEMM(op[OpA], op[OpB], &m, &n, &k, cblas::to_zptr(&alpha),
                cblas::to_zptr(a), &lda, cblas::to_zptr(b), &ldb,
                cblas::to_zptr(&beta), cblas::to_zptr(c), &ldc);
    }

#ifdef HAVE_INTEL_MKL
    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
                     const integer m, const integer n, const integer k,
                     const complex_real4 alpha, const complex_real4* a, const integer lda,
                     const real4* b, const integer ldb, const complex_real4 beta,
                     complex_real4* c, const integer ldc) {

        //static const char *op[] = { "n","t","c" };
        //F77_CSGEMM(op[OpA], op[OpB], &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

        //Don't have CSGEMM ... only SCGEMM ... so use A*B = (BT * AT)T

      //complex_real4 ctrans[m*n]; // Here assume matrices are small and can be allocated on the stack
      complex_real4* ctrans = new complex_real4[m*n];
        static const char *opT[] = { "t","n","c" }; // Transpose of op ... conj-transpose not working yet
        MADNESS_ASSERT(OpA!=ConjTrans && OpB!=ConjTrans);
        const complex_real4 zero = 0.0;
        F77_SCGEMM(opT[OpB], opT[OpA], &n, &m, &k, cblas::to_cptr(&alpha),
                   b, &ldb, cblas::to_cptr(a), &lda,
                   cblas::to_cptr(&zero), cblas::to_cptr(ctrans), &n);

        // In fortran have CTRANS(N,M) and fortran CTRANS(i,j) maps to C ctrans[j*n+i]

        if (beta == zero) {
            for (integer i=0; i<n; i++) {
                for (integer j=0; j<m; j++) {
                    c[i*ldc+j] = ctrans[j*n+i];
                }
            }
        }
        else
            for (integer i=0; i<n; i++) {
                for (integer j=0; j<m; j++) {
                    c[i*ldc+j] = beta*c[i*ldc+j] + ctrans[j*n+i];
                }
            }
	delete [] ctrans;
    }

    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
                     const integer m, const integer n, const integer k,
                     const complex_real4 alpha, const real4* a, const integer lda,
                     const complex_real4* b, const integer ldb, const complex_real4 beta,
                     complex_real4* c, const integer ldc) {
      static const char *op[] = {"n", "t", "c"};
      F77_SCGEMM(op[OpA], op[OpB], &m, &n, &k, cblas::to_cptr(&alpha),
                 a, &lda, cblas::to_cptr(b), &ldb,
                 cblas::to_cptr(&beta), cblas::to_cptr(c), &ldc);
    }

    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
                     const integer m, const integer n, const integer k,
                     const complex_real8 alpha, const complex_real8* a, const integer lda,
                     const real8* b, const integer ldb, const complex_real8 beta,
                     complex_real8* c, const integer ldc) {
        
        //static const char *op[] = { "n","t","c" };
        //F77_ZDGEMM(op[OpA], op[OpB], &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
        
        //Don't have ZDGEMM ... only DZGEMM ... so use A*B = (BT * AT)T
        
      //complex_real8 ctrans[m*n]; // Here assume matrices are small and can be allocated on the stack
      complex_real8* ctrans = new complex_real8[m*n];
        static const char *opT[] = { "t","n","c" }; // Transpose of op ... conj-transpose not working yet
        MADNESS_ASSERT(OpA!=ConjTrans && OpB!=ConjTrans);
        const complex_real8 zero = 0.0;
        F77_DZGEMM(opT[OpB], opT[OpA], &n, &m, &k, cblas::to_zptr(&alpha),
                   b, &ldb, cblas::to_zptr(a), &lda,
                   cblas::to_zptr(&zero), cblas::to_zptr(ctrans), &n);

        // In fortran have CTRANS(N,M) and fortran CTRANS(i,j) maps to C ctrans[j*n+i]
        
        if (beta == zero) {
            for (integer i=0; i<n; i++) {
                for (integer j=0; j<m; j++) {
                    c[i*ldc+j] = ctrans[j*n+i];
                }
            }
        }
        else 
            for (integer i=0; i<n; i++) {
                for (integer j=0; j<m; j++) {
                    c[i*ldc+j] = beta*c[i*ldc+j] + ctrans[j*n+i];
                }
            }
	delete [] ctrans;
    }
    
    inline void gemm(const CBLAS_TRANSPOSE OpA, const CBLAS_TRANSPOSE OpB,
                     const integer m, const integer n, const integer k,
                     const complex_real8 alpha, const real8* a, const integer lda,
                     const complex_real8* b, const integer ldb, const complex_real8 beta,
                     complex_real8* c, const integer ldc) {
      static const char *op[] = {"n", "t", "c"};
      F77_DZGEMM(op[OpA], op[OpB], &m, &n, &k, cblas::to_zptr(&alpha), a, &lda,
                 cblas::to_zptr(b), &ldb, cblas::to_zptr(&beta),
                 cblas::to_zptr(c), &ldc);
    }

#endif


    ///@}

    /// Multiplies a matrix by a vector

    /// \f[
    /// \mathbf{y} \leftarrow  \alpha \mathbf{A}^{\mathrm{OpA}} \mathbf{x} + \beta \mathbf{y}
    /// \f]
    /// \param OpA Operation to be applied to matrix \f$ \mathbf{A} \f$
    /// \param m Rows in matrix \f$ \mathbf{A} \f$
    /// \param n Columns in matrix \f$ \mathbf{A} \f$
    /// \param alpha Scaling factor applied to \f$ \mathbf{A} \f$ \c * \f$ \mathbf{x} \f$
    /// \param A Pointer to matrix \f$ \mathbf{A} \f$
    /// \param lda The size of the leading-order dimension of matrix \f$ \mathbf{A} \f$
    /// \param x Pointer to vector \f$ \mathbf{x} \f$
    /// \param incx Stride of vector \f$ \mathbf{x} \f$
    /// \param beta Scaling factor for vector \f$ \mathbf{y} \f$
    /// \param y Pointer to vector \f$ \mathbf{y} \f$
    /// \param incy Stride of vector \f$ \mathbf{y} \f$
    ///@{
    inline void gemv(const CBLAS_TRANSPOSE OpA, const integer m, const integer n,
       const float alpha, const float *A, const integer lda, const float *x,
       const integer incx, const float beta, float *y, const integer incy)
    {
        MADNESS_ASSERT(OpA != ConjTrans);
        static const char *op[] = { "n","t" };
        F77_SGEMV(op[OpA], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    }

    inline void gemv(const CBLAS_TRANSPOSE OpA, const integer m, const integer n,
       const double alpha, const double *A, const integer lda, const double *x,
       const integer incx, const double beta, double *y, const integer incy)
    {
        MADNESS_ASSERT(OpA != ConjTrans);
        static const char *op[] = { "n","t" };
        F77_DGEMV(op[OpA], &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    }

    inline void gemv(const CBLAS_TRANSPOSE OpA, const integer m, const integer n,
       const complex_real4 alpha, const complex_real4 *A, const integer lda,
       const complex_real4 *x, const integer incx, const complex_real4 beta,
       complex_real4 *y, const integer incy) {
      static const char *op[] = {"n", "t", "c"};
      F77_CGEMV(op[OpA], &m, &n, cblas::to_cptr(&alpha), cblas::to_cptr(A),
                &lda, cblas::to_cptr(x), &incx, cblas::to_cptr(&beta),
                cblas::to_cptr(y), &incy);
    }

    inline void gemv(const CBLAS_TRANSPOSE OpA, const integer m, const integer n,
       const complex_real8 alpha, const complex_real8 *A, const integer lda,
       const complex_real8 *x, const integer incx, const complex_real8 beta,
       complex_real8 *y, const integer incy) {
      static const char *op[] = {"n", "t", "c"};
      F77_ZGEMV(op[OpA], &m, &n, cblas::to_zptr(&alpha), cblas::to_zptr(A),
                &lda, cblas::to_zptr(x), &incx, cblas::to_zptr(&beta),
                cblas::to_zptr(y), &incy);
    }
    ///@}

    /// Multiplies vector \f$ \mathbf{x} \f$ by the transform of vector \f$ \mathbf{y} \f$

    /// \f[
    /// \mathbf{A} \leftarrow  \alpha \mathbf{x} \mathbf{y}^{\mathrm{T}} + \mathbf{A}
    /// \f]
    /// \param m Rows in matrix \f$ \mathbf{A} \f$
    /// \param n Columns in matrix \f$ \mathbf{A} \f$
    /// \param alpha Scaling factor applied to \f$ \mathbf{x} \mathbf{y}^{\mathrm{T}} \f$
    /// \param x Pointer to vector \f$ \mathbf{x} \f$
    /// \param incx Stride of vector \f$ \mathbf{x} \f$
    /// \param y Pointer to vector \f$ \mathbf{y} \f$
    /// \param incy Stride of vector \f$ \mathbf{y} \f$
    /// \param A Pointer to matrix \f$ \mathbf{A} \f$
    /// \param lda The size of the leading-order dimension of matrix \f$ \mathbf{A} \f$
    ///@{
    inline void ger(const integer m, const integer n, const float alpha,
        const float *x, const integer incx, const float *y, const integer incy,
        float *A, const integer lda)
    {
        F77_SGER(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
    }

    inline void ger(const integer m, const integer n, const double alpha,
        const double *x, const integer incx, const double *y, const integer incy,
        double *A, const integer lda)
    {
        F77_DGER(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
    }

    inline void ger(const integer m, const integer n, const complex_real4 alpha,
        const complex_real4 *x, const integer incx, const complex_real4 *y,
        const integer incy, complex_real4 *A, const integer lda) {
      F77_CGER(&m, &n, &alpha, x, &incx,
               y, &incy, A, &lda);
    }

    inline void ger(const integer m, const integer n, const complex_real8 alpha,
        const complex_real8 *x, const integer incx, const complex_real8 *y,
        const integer incy, complex_real8 *A, const integer lda) {
      F77_ZGER(&m, &n, &alpha, x, &incx,
               y, &incy, A, &lda);
    }
    ///@}

    /// Compute the dot product of vectors \f$ \mathbf{x} \f$ and \f$ \mathbf{y} \f$

    /// \f[
    /// u \leftarrow  \alpha \mathbf{x} \cdot \mathbf{y}
    /// \f]
    /// \param n Size of the vectors  \f$ \mathbf{x} \f$ and \f$ \mathbf{y} \f$
    /// \param x Pointer to vector \f$ \mathbf{x} \f$
    /// \param incx Stride of vector \f$ \mathbf{x} \f$
    /// \param y Pointer to vector \f$ \mathbf{y} \f$
    /// \param incy Stride of vector \f$ \mathbf{y} \f$
    /// \return The dot product of \c x and \c y
    ///@{
    inline float dot(const integer n, const float* x, const integer incx,
        const float* y, const integer incy)
    {
        return F77_SDOT(&n, x, &incx, y, &incy);
    }

    inline double dot(const integer n, const double* x, const integer incx,
        const double* y, const integer incy)
    {
        return F77_DDOT(&n, x, &incx, y, &incy);
    }

    inline complex_real4 dot(const integer n, const complex_real4* x,
        const integer incx, const complex_real4* y, const integer incy)
    {
        complex_real4 result(0.0, 0.0);
        F77_CDOTU(cblas::to_cptr(&result), &n, cblas::to_cptr(x), &incx, cblas::to_cptr(y), &incy);
        return result;
    }

    inline complex_real8 dot(const integer n, const complex_real8* x,
        const integer incx, const complex_real8* y, const integer incy)
    {
        complex_real8 result(0.0, 0.0);
        F77_ZDOTU(cblas::to_zptr(&result), &n, cblas::to_zptr(x), &incx, cblas::to_zptr(y), &incy);
        return result;
    }
    ///@}

    /// Scale a vector

    /// \f[
    /// \mathbf{x} \leftarrow \alpha \mathbf{x}
    /// \f]
    /// \param n The size of the vector
    /// \param alpha The scaling factor for vector \f$ \mathbf{x} \f$
    /// \param x Pointer to vector \f$ \mathbf{x} \f$
    /// \param incx Stride for vector \f$ \mathbf{x} \f$
    ///@{
    inline void scal(const integer n, const float alpha, float* x, const integer incx) {
      F77_SSCAL(&n, &alpha, x, &incx);
    }

    inline void scal(const integer n, const double alpha, double* x, const integer incx) {
      F77_DSCAL(&n, &alpha, x, &incx);
    }

    inline void scal(const integer n, const complex_real4 alpha, complex_real4* x, const integer incx) {
      F77_CSCAL(&n, cblas::to_cptr(&alpha), cblas::to_cptr(x), &incx);
    }

    inline void scal(const integer n, const complex_real8 alpha, complex_real8* x, const integer incx) {
      F77_ZSCAL(&n, cblas::to_zptr(&alpha), cblas::to_zptr(x), &incx);
    }

    inline void scal(const integer n, const float alpha, complex_real4* x, const integer incx) {
      F77_CSSCAL(&n, &alpha, cblas::to_cptr(x), &incx);
    }

    inline void scal(const integer n, const double alpha, complex_real8* x, const integer incx) {
      F77_ZDSCAL(&n, &alpha, cblas::to_zptr(x), &incx);
    }
    ///@}

    /// Scale and add a vector to another

    /// \f[
    /// \mathbf{y} \leftarrow \alpha \mathbf{x} + \mathbf{y}
    /// \f]
    /// \param n The size of the vector
    /// \param alpha The scaling factor for vector \f$ \mathbf{x} \f$
    /// \param x Pointer to vector \f$ \mathbf{x} \f$
    /// \param incx Stride for vector \f$ \mathbf{x} \f$
    /// \param y Pointer to vector \f$ \mathbf{y} \f$
    /// \param incy Stride for vector \f$ \mathbf{y} \f$
    ///@{
    inline void axpy(const integer n, const float alpha, float* x, const integer incx,
                     float* y, const integer incy) {
      F77_SAXPY(&n, &alpha, x, &incx, y, &incy);
    }

    inline void axpy(const integer n, const double alpha, double* x, const integer incx,
                     double* y, const integer incy) {
      F77_DAXPY(&n, &alpha, x, &incx, y, &incy);
    }

    inline void axpy(const integer n, const complex_real4 alpha, complex_real4* x, const integer incx,
                     complex_real4* y, const integer incy) {
      F77_CAXPY(&n, cblas::to_cptr(&alpha), cblas::to_cptr(x), &incx, cblas::to_cptr(y), &incy);
    }

    inline void axpy(const integer n, const complex_real8 alpha, complex_real8* x, const integer incx,
                     complex_real8* y, const integer incy) {
      F77_ZAXPY(&n, cblas::to_zptr(&alpha), cblas::to_zptr(x), &incx, cblas::to_zptr(y), &incy);
    }
    ///@}


} // namespace cblas
} // namespace madness

MADNESS_PRAGMA_CLANG(diagnostic pop)
MADNESS_PRAGMA_GCC(diagnostic pop)

#endif // MADNESS_LINALG_CBLAS_H__INCLUDED

