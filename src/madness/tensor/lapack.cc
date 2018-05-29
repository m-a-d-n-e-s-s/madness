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


#include <madness/madness_config.h>
#include <madness/tensor/tensor.h>

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::min;
using std::max;

/// \file lapack.cc
/// \brief Partial interface from Tensor to LAPACK

// Leave the test routines in this file ... they will force instantiation
// of the necessary templates.

using madness::Tensor;

#include <madness/tensor/tensor_lapack.h>
#include <madness/tensor/clapack.h>



double tt1, ss1;

#ifdef STATIC
#  undef STATIC
#endif

#if HAVE_UNQUALIFIED_STATIC_DECL
#  define STATIC static
#else
// Cray X1 compiler won't instantiate static function templates (mxm*)
#  define STATIC
#endif

/// These oddly-named wrappers enable the generic svd iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             real4 *a, integer *lda, real4 *s, real4 *u, integer *ldu,
             real4 *vt, integer *ldvt, real4 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen) {
    //std::cout << "n " << *n << " m " << *m << " lwork " << *lwork << std::endl;
    //std::cout << " sizeof(integer) " << sizeof(integer) << std::endl;
#if MADNESS_LINALG_USE_LAPACKE
    sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, info);
#else
    sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, info, jobulen, jobvtlen);
#endif
}

STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             real8 *a, integer *lda, real8 *s, real8 *u, integer *ldu,
             real8 *vt, integer *ldvt, real8 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen){
#if MADNESS_LINALG_USE_LAPACKE
  dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
          vt, ldvt, work, lwork, info);
#else
  dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, info, jobulen, jobvtlen);
#endif
}

STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             complex_real4 *a, integer *lda, real4 *s, complex_real4 *u, integer *ldu,
             complex_real4 *vt, integer *ldvt, complex_real4 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen) {
  Tensor<float> rwork(5*min(*m,*n));
#if MADNESS_LINALG_USE_LAPACKE
    cgesvd_(jobu, jobvt, m, n, reinterpret_cast<lapack_complex_float*>(a), lda, s, reinterpret_cast<lapack_complex_float*>(u), ldu,
            reinterpret_cast<lapack_complex_float*>(vt), ldvt, reinterpret_cast<lapack_complex_float*>(work), lwork, rwork.ptr(), info);
#else
    cgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, rwork.ptr(), info, jobulen, jobvtlen);
#endif
}

STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             complex_real8 *a, integer *lda, real8 *s, complex_real8 *u, integer *ldu,
             complex_real8 *vt, integer *ldvt, complex_real8 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen) {
    Tensor<double> rwork(5*min(*m,*n));
#if MADNESS_LINALG_USE_LAPACKE
    zgesvd_(jobu, jobvt, m, n, reinterpret_cast<lapack_complex_double*>(a), lda, s, reinterpret_cast<lapack_complex_double*>(u), ldu,
            reinterpret_cast<lapack_complex_double*>(vt), ldvt, reinterpret_cast<lapack_complex_double*>(work), lwork, rwork.ptr(), info);
#else
    zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, rwork.ptr(), info, jobulen, jobvtlen);
#endif
}

/// These oddly-named wrappers enable the generic gesv iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline void dgesv_(integer* n, integer* nrhs, float* AT, integer* lda,
                          integer* piv, float* x, integer* ldx, integer* info) {
    sgesv_(n, nrhs, AT, lda, piv, x, ldx, info);
}
STATIC inline void dgesv_(integer* n, integer* nrhs, float_complex* AT, integer* lda,
                          integer* piv, float_complex* x, integer* ldx, integer* info) {
#if MADNESS_LINALG_USE_LAPACKE
    cgesv_(n, nrhs, reinterpret_cast<lapack_complex_float*>(AT), lda, piv, reinterpret_cast<lapack_complex_float*>(x), ldx, info);
#else
    cgesv_(n, nrhs, AT, lda, piv, x, ldx, info);
#endif
}
STATIC inline void dgesv_(integer* n, integer* nrhs, double_complex* AT, integer* lda,
                          integer* piv, double_complex* x, integer* ldx, integer* info) {
#if MADNESS_LINALG_USE_LAPACKE
    zgesv_(n, nrhs, reinterpret_cast<lapack_complex_double*>(AT), lda, piv, reinterpret_cast<lapack_complex_double*>(x), ldx, info);
#else
    zgesv_(n, nrhs, AT, lda, piv, x, ldx, info);
#endif
}
/// These oddly-named wrappers enable the generic gelss iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.

STATIC inline void dgelss_(integer *m, integer *n, integer *nrhs,
                           float *a, integer *lda, float *b, integer *ldb, float *sOUT,
                           float *rcondIN, integer *rankOUT, float *work,
                           integer *lwork, integer *infoOUT) {
    sgelss_(m, n, nrhs, a, lda, b, ldb, sOUT, rcondIN, rankOUT, work, lwork, infoOUT);
}

STATIC inline void dgelss_(integer *m, integer *n, integer *nrhs,
                           float_complex *a, integer *lda, float_complex *b,
                           integer *ldb, float *sOUT,
                           float *rcondIN, integer *rankOUT, float_complex *work,
                           integer *lwork, integer *infoOUT) {
    Tensor<float> rwork((5*min(*m,*n)));
#if MADNESS_LINALG_USE_LAPACKE
  cgelss_(m, n, nrhs, reinterpret_cast<lapack_complex_float*>(a), lda, reinterpret_cast<lapack_complex_float*>(b), ldb, sOUT,
          rcondIN, rankOUT, reinterpret_cast<lapack_complex_float*>(work), lwork, rwork.ptr(),infoOUT);
#else
    cgelss_(m, n, nrhs, a, lda, b, ldb, sOUT, rcondIN, rankOUT, work,
            lwork, rwork.ptr(),infoOUT);
#endif
}


STATIC inline void dgelss_(integer *m, integer *n, integer *nrhs,
                           double_complex *a, integer *lda, double_complex *b,
                           integer *ldb, double *sOUT,
                           double *rcondIN, integer *rankOUT, double_complex *work,
                           integer *lwork, integer *infoOUT) {
    Tensor<double> rwork((5*min(*m,*n)));
#if MADNESS_LINALG_USE_LAPACKE
  zgelss_(m, n, nrhs, reinterpret_cast<lapack_complex_double*>(a), lda, reinterpret_cast<lapack_complex_double*>(b), ldb, sOUT,
          rcondIN, rankOUT, reinterpret_cast<lapack_complex_double*>(work), lwork, rwork.ptr(),infoOUT);
#else
    zgelss_(m, n, nrhs, a, lda, b, ldb, sOUT, rcondIN, rankOUT, work,
            lwork, rwork.ptr(),infoOUT);
#endif
}

/// These oddly-named wrappers enable the generic sygv/hegv iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            real4 *a, integer *lda, real4 *b, integer *ldb,
            real4 *w,  real4 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
#if MADNESS_LINALG_USE_LAPACKE
    ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
#else
    ssygv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, info,
           jobzlen,uplo_len);
#endif
}

STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            real8 *a, integer *lda, real8 *b, integer *ldb,
            real8 *w,  real8 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
#if MADNESS_LINALG_USE_LAPACKE
  dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
#else
  dsygv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, info,
           jobzlen,uplo_len);
#endif
}

STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            complex_real4 *a, integer *lda, complex_real4 *b, integer *ldb,
            real4 *w,  complex_real4 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<float> rwork(max((integer) 1, (integer) (3*(*n)-2)));
#if MADNESS_LINALG_USE_LAPACKE
    chegv_(itype, jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, reinterpret_cast<lapack_complex_float*>(b),
            ldb, w, reinterpret_cast<lapack_complex_float*>(work), lwork, rwork.ptr(), info);
#else
    chegv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, rwork.ptr(), info,
           jobzlen, uplo_len);
#endif
}

STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            complex_real8 *a, integer *lda, complex_real8 *b, integer *ldb,
            real8 *w,  complex_real8 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<double> rwork(max((integer) 1, (integer) (3*(*n)-2)));
#if MADNESS_LINALG_USE_LAPACKE
    zhegv_(itype, jobz, uplo,n, reinterpret_cast<lapack_complex_double*>(a), lda, reinterpret_cast<lapack_complex_double*>(b),
            ldb, w, reinterpret_cast<lapack_complex_double*>(work), lwork, rwork.ptr(), info);
#else
    zhegv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, rwork.ptr(), info,
           jobzlen, uplo_len);
#endif
}

/// These oddly-named wrappers enable the generic syev/heev iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline void dsyev_(const char* jobz, const char* uplo, integer *n,
                          real4 *a, integer *lda, real4 *w,  real4 *work,  integer *lwork,
                          integer *info, char_len jobzlen, char_len uplo_len ) {
#if MADNESS_LINALG_USE_LAPACKE
  ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
#else
  ssyev_(jobz, uplo, n, a, lda, w,  work,  lwork, info, jobzlen, uplo_len );
#endif
}

STATIC inline void dsyev_(const char* jobz, const char* uplo, integer *n,
                          real8 *a, integer *lda, real8 *w,  real8 *work,  integer *lwork,
                          integer *info, char_len jobzlen, char_len uplo_len ) {
#if MADNESS_LINALG_USE_LAPACKE
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
#else
  dsyev_(jobz, uplo, n, a, lda, w,  work,  lwork, info, jobzlen, uplo_len );
#endif
}

STATIC void dsyev_(const char* jobz, const char* uplo, integer *n,
                   complex_real4 *a, integer *lda, real4 *w,
                   complex_real4 *work,  integer *lwork,
                   integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<float> rwork(max((integer) 1, (integer) (3* (*n)-2)));
    //std::cout << *n << " " << *lda << " " << *lwork <<std::endl;
#if MADNESS_LINALG_USE_LAPACKE
    cheev_(jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, w,
           reinterpret_cast<lapack_complex_float*>(work), lwork, rwork.ptr(), info);
#else
    cheev_(jobz, uplo, n, a, lda, w,  work,  lwork, rwork.ptr(),
           info, jobzlen, uplo_len );
#endif
}

STATIC void dsyev_(const char* jobz, const char* uplo, integer *n,
                   complex_real8 *a, integer *lda, real8 *w,
                   complex_real8 *work,  integer *lwork,
                   integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<double> rwork(max((integer) 1, (integer) (3* (*n)-2)));
#if MADNESS_LINALG_USE_LAPACKE
    zheev_(jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, w,
           reinterpret_cast<lapack_complex_double*>(work), lwork, rwork.ptr(), info);
#else
    zheev_(jobz, uplo, n, a, lda, w,  work,  lwork, rwork.ptr(),
           info, jobzlen, uplo_len );
#endif
}



/// These oddly-named wrappers enable the generic orgqr/unggr iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline void dorgqr_(integer *m, integer *n, integer *k,
        real4 *a, integer *lda, real4 *tau,
        real4 *work, integer *lwork, integer *info) {
    sorgqr_(m, n, k, a, m, tau, work, lwork, info);
}

STATIC void dorgqr_(integer *m, integer *n, integer *k,
		 complex_real4 *a, integer *lda, complex_real4 *tau,
		 complex_real4 *work, integer *lwork, integer *info) {
#if MADNESS_LINALG_USE_LAPACKE
  cungqr_(m, n, k, reinterpret_cast<lapack_complex_float*>(a), lda,
          reinterpret_cast<lapack_complex_float*>(tau),
          reinterpret_cast<lapack_complex_float*>(work), lwork, info);
#else
	cungqr_(m, n, k, a, m, tau, work, lwork, info);
#endif
}

STATIC void dorgqr_(integer *m, integer *n, integer *k,
		 complex_real8 *a, integer *lda, complex_real8 *tau,
	 	 complex_real8 *work, integer *lwork, integer *info) {
#if MADNESS_LINALG_USE_LAPACKE
  zungqr_(m,n,k, reinterpret_cast<lapack_complex_double*>(a), m,
          reinterpret_cast<lapack_complex_double*>(tau),
          reinterpret_cast<lapack_complex_double*>(work), lwork, info);
#else
	zungqr_(m, n, k, a, m, tau, work, lwork, info);
#endif
}

namespace madness {

    static void mask_info(integer& info) {
        if ( (info&0xffffffff) == 0) info = 0;
    }

    /** \brief   Compute the singluar value decomposition of an n-by-m matrix using *gesvd.

    Returns via arguments U, s, VT where

    A = U * diag(s) * VT    for A real
    A = U * diag(s) * VH    for A complex

    or

    UT * A * V = diag(s)   for A real
    UH * A * V = diag(s)   for A complex

    If A is [m,n] and r=min(m,n) then we have U[m,r], s[r], and VT[r,n]

    On failure, throws TensorException with value set to Lapack's info.
    */
    template <typename T>
    void svd(const Tensor<T>& a, Tensor<T>& U,
             Tensor< typename Tensor<T>::scalar_type >& s, Tensor<T>& VT) {
        TENSOR_ASSERT(a.ndim() == 2, "svd requires matrix",a.ndim(),&a);
        integer m = a.dim(0), n = a.dim(1), rmax = min<integer>(m,n);
        integer lwork = max<integer>(3*min(m,n)+max(m,n),5*min(m,n)-4)*32;
        integer info;
        Tensor<T> A(copy(a)), work(lwork);

        s = Tensor< typename Tensor<T>::scalar_type >(rmax);
        U = Tensor<T>(m,rmax);
        VT = Tensor<T>(rmax,n);
        dgesvd_("S","S", &n, &m, A.ptr(), &n, s.ptr(),
                VT.ptr(), &n, U.ptr(), &rmax, work.ptr(), &lwork,
                &info, (char_len) 1, (char_len) 1);

        //std::cout << "n " << n << " m " << m << " lwork " << lwork << std::endl;
	//std::cout << sizeof(long) << " " << sizeof(int) << " " << sizeof(integer) << std::endl;
	//std::cout << sizeof(real4) << " " << sizeof(float)  << std::endl;
        mask_info(info);

        TENSOR_ASSERT(info == 0, "svd: Lapack failed", info, &a);
    }

    /// same as svd, but it optimizes away the tensor construction: a = U * diag(s) * VT

    /// note that S and VT are swapped in the calling list for c/fortran consistency!
    /// @param[in,out]	a	a (m,n) matrix to be svd'ed; upon return will hold VT
    /// 					the first min(m,n) rows of VT, stored rowwise
    /// @param[in,out]	U	left singular vectors, stored columnwise
    /// @param[in,out]	s	the singular values
    /// @param[in,out]	VT	not referenced
    /// @param[in,out]	work work array; optimial length lwork:
    ///						lwork = max<integer>(3*min(m,n)+max(m,n),5*min(m,n)-4)*32;
    template <typename T>
    void svd_result(Tensor<T>& a, Tensor<T>& U,
             Tensor< typename Tensor<T>::scalar_type >& s, Tensor<T>& VT, Tensor<T>& work) {
        TENSOR_ASSERT(a.ndim() == 2, "svd requires matrix",a.ndim(),&a);

        integer m = a.dim(0), n = a.dim(1), rmax = min<integer>(m,n);
//        integer lwork = max<integer>(3*min(m,n)+max(m,n),5*min(m,n)-4)*32;
        integer lwork=work.size();
        integer info;

        // calling list is swapped
        dgesvd_("O", "S", &n, &m, a.ptr(), &n, s.ptr(),
                VT.ptr(), &n, U.ptr(), &rmax, work.ptr(), &lwork,
                &info, (char_len) 1, (char_len) 1);

        mask_info(info);

        TENSOR_ASSERT(info == 0, "svd: Lapack failed", info, &a);
    }

    /** \brief  Solve Ax = b for general A using the LAPACK *gesv routines.

    A should be a square matrix (float, double, float_complex,
    double_complex) and b should be either a vector, or a matrix with
    each vector stored in a column (i.e., b[n,nrhs]).

    If the LAPACK routine fails, it throws a TensorException with the
    LAPACK info as the value.  Otherwise, it returns the solution(s).
    The input A and b are unchanged.  There is no need to worry about
    Python/C/Fortran ordering issues.  It will solve Ax=b as written.
    */
    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x) {
        TENSOR_ASSERT(a.ndim() == 2, "gesv requires matrix",a.ndim(),&a);
        integer n = a.dim(0), m = a.dim(1), nrhs = b.dim(1);
        TENSOR_ASSERT(m == n, "gesv requires square matrix",0,&a);
        TENSOR_ASSERT(b.ndim() <= 2, "gesv require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gesv matrix and RHS must conform",b.ndim(),&b);

        // The input matrix & vectors are destroyed by gesv and we also need Fortran order
        Tensor<T> AT = transpose(a);
        if (b.ndim() == 1)
            x = copy(b);
        else
            x = transpose(b);

        Tensor<integer> piv(n);
        integer info;

        // note overriding of dgesv for other types above
        dgesv_(&n, &nrhs, AT.ptr(), &n, piv.ptr(), x.ptr(), &n, &info);
        mask_info(info);

        TENSOR_ASSERT((info == 0), "gesv failed", info, &a);

        if (b.ndim() == 2) x = transpose(x);
    }

    /// invert general square matrix A
    template<typename T>
    Tensor<T> inverse(const Tensor<T>& a_in) {
        Tensor<T> a=copy(a_in);
        TENSOR_ASSERT(a.ndim() == 2, "inverse requires matrix",a.ndim(),&a);
        TENSOR_ASSERT(a.dim(0) == a.dim(1), "inverse requires square matrix",a.ndim(),&a);
        integer n = a.dim(0);;

        Tensor<integer> ipiv(n);
        integer info;

        // DGETRF computes an LU factorization of a general M-by-N matrix A
        // using partial pivoting with row interchanges.
        dgetrf_(&n,&n,a.ptr(),&n,ipiv.ptr(),&info);

        // DGETRI computes the inverse of a matrix using the LU factorization
        // computed by DGETRF.
        integer lwork=(10*n);
        Tensor<T> work(lwork);
        dgetri_(&n,a.ptr(),&n,ipiv.ptr(),work.ptr(),&lwork,&info);

        mask_info(info);
        TENSOR_ASSERT((info == 0), "inverse failed", info, &a);
        return a;
    }

    /** \brief  Solve Ax = b for general A using the LAPACK *gelss routines.

    A should be a matrix (float, double, float_complex,
    double_complex) and b should be either a vector, or a matrix with
    each vector stored in a column (i.e., b[n,nrhs]).

    If the LAPACK routine fails, it throws a TensorException with the
    LAPACK info as the value.  Otherwise, it returns the solution(s).
    The input A and b are unchanged.  There is no need to worry about
    Python/C/Fortran ordering issues.  It will solve Ax=b as written.

    This from the LAPACK documentation
    \verbatim
    RCOND   (input) REAL
    RCOND is used to determine the effective  rank  of A.
    Singular values S(i) <= RCOND*S(1) are treated
    as zero.  If RCOND < 0, machine precision is  used
    instead.

    RANK    (output) INTEGER
    The  effective rank of A, i.e., the number of singular
    values which are greater than RCOND*S(1).
    \endverbatim

    Finally, the optional vector sumsq will store the sum-of-squares
    residual in the case of a rectangular matrix (least squares regression).
    */
    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>& x, Tensor< typename Tensor<T>::scalar_type >& s,
               long& rank, Tensor<typename Tensor<T>::scalar_type>& sumsq) {
        TENSOR_ASSERT(a.ndim() == 2, "gelss requires matrix",a.ndim(),&a);
        integer m = a.dim(0), n = a.dim(1), nrhs = b.dim(1);
        TENSOR_ASSERT(b.ndim() <= 2, "gelss require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gelss matrix and RHS must conform",b.ndim(),&b);

        // The input matrix & vectors are destroyed by gelss and we also need Fortran order
        integer maxmn = max(m, n);
        Tensor<T> AT = transpose(a);
        Tensor<T> lapack_inout;

        if (b.ndim() == 1)
            lapack_inout = copy(b);
        else {
            if (m >= n)
                lapack_inout = transpose(b);
            else {
                // dgelss_ uses the same physical array for both b (input) and
                // x (output).  for a rectangular matrix A with more columns,
                // b is a smaller array than x, and the data needs to be
                // manipulated (more than just a transpose) to make it work
                // with LAPACK
                lapack_inout = Tensor<T>(nrhs, maxmn);
                lapack_inout(Slice(0, nrhs-1), Slice(0, m-1)) =
                    transpose(b);
            }
        }

        integer lwork=(3*min(m,n)+max(max(2*min(m,n),maxmn),nrhs))*32;
        Tensor<T> work(lwork);
        typedef typename TensorTypeData<T>::scalar_type scalar_type;
        s = Tensor< scalar_type >(n);
        integer info;
        scalar_type rrcond = rcond;
        integer rrank=0;

        dgelss_(&m, &n, &nrhs, AT.ptr(), &m, lapack_inout.ptr(), &maxmn,
                s.ptr(), &rrcond, &rrank, work.ptr(), &lwork, &info);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "gelss failed", info, &a);

        rank = rrank;

        if(m > n) {
            // have a similar problem where the lapack_inout tensor is padded
            // the padding gives information on the fit

            // get the sum-of-squares for the various fits
            sumsq = Tensor<scalar_type>(nrhs);
            if(nrhs == 1) {
                sumsq[0] = lapack_inout(Slice(n, m-1)).normf();
            }
            else {
                for(integer i = 0; i < nrhs; ++i) {
                    sumsq[i] =
                        lapack_inout(Slice(i, i), Slice(n, m-1)).normf();
                }
            }

            if(b.ndim() == 1)
                x = lapack_inout(Slice(0,n-1));
            else
                x = transpose(lapack_inout(Slice(0,nrhs-1), Slice(0, n-1)));
        }
        else if(b.ndim() == 2)
            x = transpose(lapack_inout);
        else
            x = lapack_inout;
    }

    /** \brief   Real-symmetric or complex-Hermitian eigenproblem.

    A is a real symmetric or complex Hermitian matrix.  Return V and e
    where V is a matrix whose columns are the eigenvectors and e is a
    vector containing the corresponding eigenvalues.  If the LAPACK
    routine fails, it raises a TensorException with value=infor.  The
    input matrix is unchanged.  The eigenvalues are sorted into ascending
    order.  s/dsyev are used for real symmetric matrices; c/zheev are used
    for complex Hermitian.

    The results will satisfy A*V(_,i) = V(_,i)*e(i).
    */
    template <typename T>
    void syev(const Tensor<T>& A,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(A.ndim() == 2, "syev requires a matrix",A.ndim(),&A);
        TENSOR_ASSERT(A.dim(0) == A.dim(1), "syev requires square matrix",0,&A);
        integer n = A.dim(0);
        integer lwork = max(max((integer) 1,(integer) (3*n-1)),(integer) (34*n));
        integer info;
        Tensor<T> work(lwork);
        V = transpose(A);		// For Hermitian case
        e = Tensor<typename Tensor<T>::scalar_type>(n);
        dsyev_("V", "U", &n, V.ptr(), &n, e.ptr(), work.ptr(), &lwork, &info,
               (char_len) 1, (char_len) 1);

        mask_info(info);
        TENSOR_ASSERT(info == 0, "(s/d)syev/(c/z)heev failed", info, &A);
        V = transpose(V);
    }

    /** \brief  Generalized real-symmetric or complex-Hermitian eigenproblem.

    This from the LAPACK documentation

    \verbatim
    S/DSYGV computes all the eigenvalues, and optionally, the eigenvectors
    of a real generalized symmetric-definite eigenproblem, of the form
    A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x.  Here A and B
    are assumed to be symmetric and B is also positive definite.

    C/ZHEGV computes all the eigenvalues, and optionally, the eigenvectors
    of a complex generalized Hermitian-definite eigenproblem, of the form
    A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. Here A and B
    are assumed to be Hermitian and B is also positive definite.

    ITYPE   (input) INTEGER
    Specifies the problem type to be solved:
    = 1:  A*x = (lambda)*B*x
    = 2:  A*B*x = (lambda)*x
    = 3:  B*A*x = (lambda)*x
    \endverbatim

    */
    template <typename T>
    void sygv(const Tensor<T>& A, const Tensor<T>& B, int itype,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(A.ndim() == 2, "sygv requires a matrix",A.ndim(),&A);
        TENSOR_ASSERT(A.dim(0) == A.dim(1), "sygv requires square matrix",0,&A);
        TENSOR_ASSERT(B.ndim() == 2, "sygv requires a matrix",B.ndim(),&A);
        TENSOR_ASSERT(B.dim(0) == B.dim(1), "sygv requires square matrix",0,&A);
        integer ity = itype;
        integer n = A.dim(0);
        integer lwork = max((integer)1,(integer)(3*n-1))*32;
        integer info;
        Tensor<T> work(lwork);
        Tensor<T> b = transpose(B);	// For Hermitian case
        V = transpose(A);		// For Hermitian case
        e = Tensor<typename Tensor<T>::scalar_type>(n);
        dsygv_(&ity, "V", "U", &n, V.ptr(), &n, b.ptr(), &n,
               e.ptr(), work.ptr(), &lwork, &info,
               (char_len) 1, (char_len) 1);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "sygv/hegv failed", info, &A);
        V = transpose(V);
    }

    /** \brief  Compute the Cholesky factorization.

    Compute the Cholesky factorization of the symmetric positive definite matrix A

    For memory efficiency A is modified inplace.  Its upper
    triangle will hold the result and the lower trianlge will be
    zeroed such that input = inner(transpose(output),output).

    */
    template <typename T>
    void cholesky(Tensor<T>& A) {
        integer n = A.dim(0);
        integer info;

#if MADNESS_LINALG_USE_LAPACKE
        dpotrf_("L", &n, A.ptr(), &n, &info);
#else
        dpotrf_("L", &n, A.ptr(), &n, &info, 1);
#endif
        mask_info(info);
        TENSOR_ASSERT(info == 0, "cholesky: Lapack failed", info, &A);

        for (int i=0; i<n; ++i)
            for (int j=0; j<i; ++j)
                A(i,j) = 0.0;
    }

    /** \brief  Compute the QR factorization.

	Q is returned in the lapack-specific format

    */
    template<typename T>
    void geqp3(Tensor<T>& A, Tensor<T>& tau, Tensor<integer>& jpvt) {

    	TENSOR_ASSERT(A.ndim() == 2, "geqp requires a matrix",A.ndim(),&A);

    	A=transpose(A);
    	integer m=A.dim(0);
    	integer n=A.dim(1);
    	jpvt=Tensor<integer>(n);
    	tau=Tensor<T>(std::min(n,m));
    	integer lwork=2*n+(n+1)*64;
    	Tensor<T> work(lwork);
    	geqp3_result(A,tau,jpvt,work);
    	A=transpose(A);
    }

    template<typename T>
    void geqp3_result(Tensor<T>& A, Tensor<T>& tau, Tensor<integer>& jpvt,
    		Tensor<T>& work) {
    	integer m=A.dim(0);
    	integer n=A.dim(1);
    	integer lwork=work.size();
    	integer info;

//    	dgeqp3(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO );
    	std::cout << jpvt[0] << std::endl;
    	dgeqp3_(&m, &n, A.ptr(), &m, jpvt.ptr(), tau.ptr(), work.ptr(),
    			&lwork, &info);
    	std::cout << jpvt[0] << std::endl;
        mask_info(info);
        TENSOR_ASSERT(info == 0, "dgeqp3: Lapack failed", info, &A);
    }

    /// compute the QR decomposition of the matrix A

    /// @param[in,out]	A	on entry the (n,m) matrix to be decomposed
    ///						on exit the Q matrix
    /// @param[out]		R	the (n,n) matrix R (square form)
    template<typename T>
    void qr(Tensor<T>& A, Tensor<T>& R) {
    	TENSOR_ASSERT(A.ndim() == 2, "qr requires a matrix",A.ndim(),&A);
    	integer m=A.dim(0);
    	integer n=A.dim(1);
    	Tensor<T> tau(std::min(n,m));
    	integer lwork=2*n+(n+1)*64;
    	Tensor<T> work(lwork);

        integer r_rows=std::min(m,n);
        integer r_cols=n;
        R=Tensor<T>(r_rows,r_cols);

        A=transpose(A);
    	lq_result(A,R,tau,work,true);

    	A=transpose(A);
    }

    /// compute the LQ decomposition of the matrix A = L Q

    /// @param[in,out]	A	on entry the (n,m) matrix to be decomposed
    ///						on exit the Q matrix
    /// @param[out]		R	the (n,n) matrix L (square form) \todo MGR: Check this.
    template<typename T>
    void lq(Tensor<T>& A, Tensor<T>& R) {
    	TENSOR_ASSERT(A.ndim() == 2, "lq requires a matrix",A.ndim(),&A);

    	integer m=A.dim(0);
    	integer n=A.dim(1);
    	Tensor<T> tau(std::min(n,m));
    	integer lwork=2*n+(n+1)*64;
    	Tensor<T> work(lwork);

        integer r_rows=std::min(m,n);
        integer r_cols=m;
        R=Tensor<T>(r_cols,r_rows);

    	lq_result(A,R,tau,work,false);
    }

    /// compute the LQ decomposition of the matrix A = L Q

    /// @param[in,out]	A	on entry the (n,m) matrix to be decomposed
    ///						on exit the Q matrix
    /// @param[in,out]	R	the (n,n) matrix L (square form) \todo MGR: Check this.
    /// @param[in]      do_qr   do QR instead of LQ
    template<typename T>
    void lq_result(Tensor<T>& A, Tensor<T>& R, Tensor<T>& tau, Tensor<T>& work,
    		bool do_qr) {
    	TENSOR_ASSERT(A.ndim() == 2, "lq requires a matrix",A.ndim(),&A);

    	integer m=A.dim(1);		// -- use transpose(A)
    	integer n=A.dim(0);
    	integer lwork=work.size();
    	integer info;

    	dgeqrf_(&m, &n, A.ptr(), &m, tau.ptr(), work.ptr(),
    			&lwork, &info);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "dgeqrf_: Lapack failed", info, &A);

    	// reconstruction of R	 -- use transpose(A)
        integer r_rows= (m>=n) ? n : m;
        integer r_cols=n;
		if (do_qr) {
//			R=Tensor<T>(r_rows,r_cols);
	        TENSOR_ASSERT(r_rows==R.dim(0),"confused dimensions 0",r_rows,&R);
	        TENSOR_ASSERT(r_cols==R.dim(1),"confused dimensions 1",r_cols,&R);
			for (int i=0; i<r_rows; ++i) {
				for (int j=i; j<r_cols; ++j) {
					R(i,j)=A(j,i);	// <- transpose(A)
				}
			}
		} else {
			//R=Tensor<T>(r_cols,r_rows);
            TENSOR_ASSERT(r_rows==R.dim(1),"confused dimensions 1",r_rows,&R);
            TENSOR_ASSERT(r_cols==R.dim(0),"confused dimensions 0",r_cols,&R);
			for (int i=0; i<r_rows; ++i) {
				for (int j=i; j<r_cols; ++j) {
					R(j,i)=A(j,i);
				}
			}
		}

		// reconstruction of Q
    	integer k=std::min(m,n);//tau.size(); // size of tau is tau(min(m,n))
    	integer q_rows=m;
    	integer q_cols= (m>=n) ? n : m;
    	dorgqr_(&q_rows, &q_cols, &k, A.ptr(), &q_rows, const_cast<T*>(tau.ptr()),
    			work.ptr(), &lwork, &info);
    	A=A(Slice(0,q_cols-1),Slice(0,q_rows-1));	// -- use transpose(A)
        mask_info(info);
        TENSOR_ASSERT(info == 0, "dorgqr_: Lapack failed", info, &A);

    }


    /// reconstruct the orthogonal matrix Q (e.g. from QR factorization)

    /// @param[in,out]	A
    ///    On entry, the i-th column must contain the vector which
    ///    defines the elementary reflector H(i), for i = 1,2,...,k, as
    ///    returned by DGEQRF in the first k columns of its array
    ///    argument A.
    ///    On exit, the M-by-N matrix Q.
    /// @param[in]	tau
    ///	   TAU(i) must contain the scalar factor of the elementary
    ///    reflector H(i), as returned by DGEQRF.
    template<typename T>
    void orgqr(Tensor<T>& A, const Tensor<T>& tau) {
    	A=transpose(A);
    	integer n=A.dim(0);
    	integer m=A.dim(1);
    	integer k=tau.size();
    	integer lwork=64*n;
    	Tensor<T> work(lwork);
    	integer info;
    	dorgqr_(&m, &n, &k, A.ptr(), &m, const_cast<T*>(tau.ptr()),
    			work.ptr(), &lwork, &info);
        mask_info(info);
    	A=transpose(A);
        TENSOR_ASSERT(info == 0, "xorgqr: Lapack failed", info, &A);
    }


//     template <typename T>
//     void triangular_solve(const Tensor<T>& L, Tensor<T>& B, const char* side, const char* transa) {
//         integer n = L.dim(0);  // ????
//         integer m = L.dim(1);
//         double one = 1.0;
//         integer lda = n; // ???
//         integer ldb = m; //???

//         dtrsm(side, "L", transa, "N", m, n, one, L.ptr(), lda, B.ptr, ldb, 1, 1, 1, 1);
//     }


    // This stupidity since current defn of conj_tranpose() only
    // takes complex types ... was that a sensible design choice?
    STATIC Tensor<float> my_conj_transpose(Tensor<float> a) {
        return transpose(a);
    }
    STATIC Tensor<double> my_conj_transpose(Tensor<double> a) {
        return transpose(a);
    }
    STATIC Tensor<float_complex> my_conj_transpose(Tensor<float_complex> a) {
        return conj_transpose(a);
    }
    STATIC Tensor<double_complex> my_conj_transpose(Tensor<double_complex> a) {
        return conj_transpose(a);
    }



    /// Example and test code for interface to LAPACK SVD interfae
    template <typename T>
    double test_svd(int n, int m) {
        Tensor<T> a(n,m), U, VT;
        typedef typename TensorTypeData<T>::scalar_type scalar_type;
        Tensor<scalar_type> s;

        a.fillrandom();
        svd(a,U,s,VT);

        long rank = s.dim(0);
        Tensor<T> b(n,m);
        for (long i=0; i<n; ++i)
            for (long j=0; j<m; ++j)
                for (long k=0; k<rank; ++k)
                    b(i,j) += U(i,k) * T(s(k)) * VT(k,j);

        b -= a;

        return b.absmax();
    }

    /// Example and test code for interface to LAPACK SVD interfae
    template <typename T>
    double test_inverse(int n) {
        Tensor<T> A(n,n);

        A.fillrandom();
        Tensor<T> Ainv=inverse(A);
        Tensor<T> id1=inner(A,Ainv);
        Tensor<T> id2=inner(Ainv,A);

        for (int i=0; i<n; ++i) {
            id1(i,i)-=1.0;
            id2(i,i)-=1.0;
        }

        double error1=id1.normf()/id1.size();
        double error2=id2.normf()/id2.size();

        return error1+error2;
    }

    template <typename T>
    double test_gesv(int n, int nrhs) {
        Tensor<T> a(n,n), b1(n), b(n,nrhs), x1, x;

        a.fillrandom();
        b1.fillrandom();
        b.fillrandom();

        gesv(a,b,x);
        gesv(a,b1,x1);

        double err = 0;
        err = (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
        return err;
    }
    template <typename T>
    double test_syev(int n) {
        Tensor<T> a(n,n), V;
        Tensor< typename Tensor<T>::scalar_type > e;
        a.fillrandom();
        //a += madness::my_conj_transpose(a);
        a += madness::transpose(a);
        syev(a,V,e);
        double err = 0.0;
        for (int i=0; i<n; ++i) {
            err = max(err,(double) (inner(a,V(_,i)) - V(_,i)*e(i)).normf());
          //err = max(err,(double) (inner(a,V(_,i)) - V(_,i)*((T) e(i))).normf());
        }
        return err;
    }

    template <typename T>
    double test_gelss(int n, int nrhs) {
        Tensor<T> a(n,n), b1(n), b(n,nrhs), x1, x;
        Tensor< typename Tensor<T>::scalar_type > s, sumsq;
        long rank;

        a.fillrandom();
        b1.fillrandom();
        b.fillrandom();

        gelss(a,b,1e-5,x,s,rank,sumsq);
        gelss(a,b1,1e-5,x1,s,rank,sumsq);

        return (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
    }

    template <typename T>
    double test_sygv(int n) {
        Tensor<T> a(n,n), V, b(n,n);
        Tensor< typename Tensor<T>::scalar_type > e;

        a.fillrandom();
        b.fillrandom();
        a += madness::my_conj_transpose(a);
        b += madness::my_conj_transpose(b);

        for (int i=0; i<n; ++i) b(i,i) = 2*n;	// To make pos-def
        sygv(a,b,1,V,e);
        double err = 0.0;
        for (int i=0; i<n; ++i) {
            err = max(err,(double) (inner(a,V(_,i)) - inner(b,V(_,i))*(T) e(i)).normf());
        }
        return err;
    }

    template <typename T>
    double test_cholesky(int n) {
        Tensor<T> a(n,n);
        a.fillrandom();
        a += madness::my_conj_transpose(a);
        for (int i=0; i<n; ++i) a(i,i) += n;

        Tensor<T> aa = copy(a);
        cholesky(a);
        Tensor<T> LLT = inner(my_conj_transpose(a),a);
        return (LLT - aa).normf()/n;
    }

    template <typename T>
    double test_qr() {

		Tensor<T> R;
		double error=0.0;

		// test QR and LQ with random dimensions
		for (int i=0; i<20; ++i) {
			madness::Random(100);
			int n= (madness::RandomValue<long>() % 4) + 2;
			int m= (madness::RandomValue<long>() % 3) + 3;

			{	// QR
				Tensor<T> AA(n,m);
				AA.fillrandom();
				Tensor<T> A=copy(AA);
				qr(A,R);
				Tensor<T> A_reconstructed=inner(A,R);
				double err1=(AA-A_reconstructed).normf();
				std::cout << "test QR: n, m, error" << " " << n << " " << m
						 << " " << err1 << std::endl;;
				error+=err1;
			}

			{	// LQ
				Tensor<T> AA(n,m);
				AA.fillrandom();
				Tensor<T> A=copy(AA);
				lq(A,R);
				Tensor<T> A_reconstructed=inner(R,A);
				double err1=(AA-A_reconstructed).normf();
				std::cout << "test LQ: n, m, error" << " " << n << " " << m
						 << " " << err1 << std::endl;;
				error+=err1;
			}
		}
		return error;
    }


    void init_tensor_lapack() {
	char e[] = "e";
#if MADNESS_LINALG_USE_LAPACKE
    dlamch(e);
	slamch(e);
#else
	dlamch_(e,1);
	slamch_(e,1);
#endif

// 	char modes[] = "esbpnrmulo";
// 	for (int i=0; i<10; ++i) {
// 	    cout << "init_tensor_lapack: dlamch: " << modes[i] << " = " << dlamch_(modes+i,1) << endl;
// 	}
    }


    /// Test the Tensor-LAPACK interface ... currently always returns true!
    bool test_tensor_lapack() {
        try {
            cout << "error in float svd " << test_svd<float>(20,30) << endl;
            cout << "error in double svd " << test_svd<double>(30,20) << endl;
            cout << "error in float_complex svd " << test_svd<float_complex>(23,27) << endl;
            cout << "error in double_complex svd " << test_svd<double_complex>(37,19) << endl;
            cout << endl;


            cout << "error in float  gelss " << test_gelss<float>(20,30) << endl;
            cout << "error in double gelss " << test_gelss<double>(30,20) << endl;
            cout << "error in float_complex gelss " << test_gelss<float_complex>(23,27) << endl;
            cout << "error in double_complex gelss " << test_gelss<double_complex>(37,19) << endl;
            cout << endl;

            cout << "error in double syev " << test_syev<double>(21) << endl;
            cout << "error in float syev " << test_syev<float>(21) << endl;
            cout << "error in float_complex syev " << test_syev<float_complex>(21) << endl;
            cout << "error in double_complex syev " << test_syev<double_complex>(21) << endl;
            cout << endl;


            cout << "error in float sygv " << test_sygv<float>(20) << endl;
            cout << "error in double sygv " << test_sygv<double>(20) << endl;
            cout << "error in float_complex sygv " << test_sygv<float_complex>(23) << endl;
            cout << "error in double_complex sygv " << test_sygv<double_complex>(24) << endl;
            cout << endl;

            cout << "error in float gesv " << test_gesv<float>(20,30) << endl;
            cout << "error in double gesv " << test_gesv<double>(20,30) << endl;
            cout << "error in float_complex gesv " << test_gesv<float_complex>(23,27) << endl;
            cout << "error in double_complex gesv " << test_gesv<double_complex>(37,19) << endl;
            cout << endl;
            cout << "error in double cholesky " << test_cholesky<double>(22) << endl;
            cout << endl;

            cout << endl;
            cout << "error in double QR/LQ " << test_qr<double>() << endl;
            cout << endl;

            cout << "error in double inverse " << test_inverse<double>(32) << endl;
            cout << "error in double inverse " << test_inverse<double>(47) << endl;
            cout << endl;
        }

        catch (TensorException e) {
            cout << "Caught a tensor exception in test_tensor_lapack\n";
            cout << e;
            return false;
        }

        return true;			//
    }

    // int main() {
    //   test_tensor_lapack();
    //   return 0;
    // }

    // GCC 4.4.3 seems to want these explicitly instantiated whereas previous
    // versions were happy with the instantiations caused by the test code above

    template
    void svd_result(Tensor<float>& a, Tensor<float>& U,
             Tensor<Tensor<float>::scalar_type >& s, Tensor<float>& VT, Tensor<float>& work);

    template
    void orgqr(Tensor<float>& A, const Tensor<float>& tau);


    template
    void svd(const Tensor<double>& a, Tensor<double>& U,
             Tensor<Tensor<double>::scalar_type >& s, Tensor<double>& VT);

    template
    void svd_result(Tensor<double>& a, Tensor<double>& U,
             Tensor<Tensor<double>::scalar_type >& s, Tensor<double>& VT, Tensor<double>& work);

    template
    void gelss(const Tensor<double>& a, const Tensor<double>& b, double rcond,
               Tensor<double>& x, Tensor<Tensor<double>::scalar_type >& s,
               long &rank, Tensor<Tensor<double>::scalar_type>& sumsq);

    template
    void syev(const Tensor<double>& A,
              Tensor<double>& V, Tensor<Tensor<double>::scalar_type >& e);


    template
    void cholesky(Tensor<double>& A);

    template
    Tensor<double> inverse(const Tensor<double>& A);

    template
    void qr(Tensor<double>& A, Tensor<double>& R);

    template
    void lq(Tensor<double>& A, Tensor<double>& L);

    template
    void lq_result(Tensor<double>& A, Tensor<double>& R, Tensor<double>& tau, Tensor<double>& work,
    		bool do_qr);


    template
    void geqp3(Tensor<double>& A, Tensor<double>& tau, Tensor<integer>& jpvt);

//     template
//     void triangular_solve(const Tensor<double>& L, Tensor<double>& B,
//                           const char* side, const char* transa);

    template
    void orgqr(Tensor<double>& A, const Tensor<double>& tau);


    template
    void svd_result(Tensor<float_complex>& a, Tensor<float_complex>& U,
             Tensor<Tensor<float_complex>::scalar_type >& s, Tensor<float_complex>& VT,
             Tensor<float_complex>& work);


    template
    void svd(const Tensor<double_complex>& a, Tensor<double_complex>& U,
             Tensor<Tensor<double_complex>::scalar_type >& s, Tensor<double_complex>& VT);

    template
    void svd_result(Tensor<double_complex>& a, Tensor<double_complex>& U,
             Tensor<Tensor<double_complex>::scalar_type >& s, Tensor<double_complex>& VT,
             Tensor<double_complex>& work);

    template
    void gelss(const Tensor<double_complex>& a, const Tensor<double_complex>& b, double rcond,
               Tensor<double_complex>& x, Tensor<Tensor<double_complex>::scalar_type >& s,
               long &rank, Tensor<Tensor<double_complex>::scalar_type>& sumsq);

    template
    void syev(const Tensor<double_complex>& A,
              Tensor<double_complex>& V, Tensor<Tensor<double_complex>::scalar_type >& e);


//     template
//     void triangular_solve(const Tensor<double_complex>& L, Tensor<double_complex>& B,
//                           const char* side, const char* transa);
    template
    void gesv(const Tensor<double>& a, const Tensor<double>& b, Tensor<double>& x);

    template
    void gesv(const Tensor<double_complex>& a, const Tensor<double_complex>& b, Tensor<double_complex>& x);

    template
    void sygv(const Tensor<double>& A, const Tensor<double>& B, int itype,
              Tensor<double>& V, Tensor<Tensor<double>::scalar_type >& e);
    template
    void sygv(const Tensor<double_complex>& A, const Tensor<double_complex>& B, int itype,
              Tensor<double_complex>& V, Tensor<Tensor<double_complex>::scalar_type >& e);

    template
    void orgqr(Tensor<complex_real4>& A, const Tensor<complex_real4>& tau);

    template
    void orgqr(Tensor<double_complex>& A, const Tensor<double_complex>& tau);


} // namespace madness
