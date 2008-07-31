/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

#include <tensor/tensor.h>
using madness::Tensor;

#include <linalg/tensor_lapack.h>
#include <linalg/clapack.h>

#define STATIC static

#ifdef HAVE_UNQUALIFIED_STATIC_DECL
#  ifdef STATIC
#    undef STATIC
#  endif
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
    sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, info, jobulen, jobvtlen);
}

STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             complex_real4 *a, integer *lda, real4 *s, complex_real4 *u, integer *ldu,
             complex_real4 *vt, integer *ldvt, complex_real4 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen) {
    Tensor<float> rwork(5*min(*m,*n));
    cgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, rwork.ptr(), info, jobulen, jobvtlen);
}

STATIC inline
void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
             complex_real8 *a, integer *lda, real8 *s, complex_real8 *u, integer *ldu,
             complex_real8 *vt, integer *ldvt, complex_real8 *work, integer *lwork,
             integer *info, char_len jobulen, char_len jobvtlen) {
    Tensor<double> rwork(5*min(*m,*n));
    zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu,
            vt, ldvt, work, lwork, rwork.ptr(), info, jobulen, jobvtlen);
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
    cgesv_(n, nrhs, AT, lda, piv, x, ldx, info);
}
STATIC inline void dgesv_(integer* n, integer* nrhs, double_complex* AT, integer* lda,
                          integer* piv, double_complex* x, integer* ldx, integer* info) {
    zgesv_(n, nrhs, AT, lda, piv, x, ldx, info);
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
    cgelss_(m, n, nrhs, a, lda, b, ldb, sOUT, rcondIN, rankOUT, work,
            lwork, rwork.ptr(),infoOUT);
}


STATIC inline void dgelss_(integer *m, integer *n, integer *nrhs,
                           double_complex *a, integer *lda, double_complex *b,
                           integer *ldb, double *sOUT,
                           double *rcondIN, integer *rankOUT, double_complex *work,
                           integer *lwork, integer *infoOUT) {
    Tensor<double> rwork((5*min(*m,*n)));
    zgelss_(m, n, nrhs, a, lda, b, ldb, sOUT, rcondIN, rankOUT, work,
            lwork, rwork.ptr(),infoOUT);
}

/// These oddly-named wrappers enable the generic sygv/hegv iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            real4 *a, integer *lda, real4 *b, integer *ldb,
            real4 *w,  real4 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
    ssygv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, info,
           jobzlen,uplo_len);
}

STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            complex_real4 *a, integer *lda, complex_real4 *b, integer *ldb,
            real4 *w,  complex_real4 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<float> rwork(max((integer) 1, (integer) (3*(*n)-2)));
    chegv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, rwork.ptr(), info,
           jobzlen, uplo_len);
}

STATIC inline
void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
            complex_real8 *a, integer *lda, complex_real8 *b, integer *ldb,
            real8 *w,  complex_real8 *work,  integer *lwork,
            integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<double> rwork(max((integer) 1, (integer) (3*(*n)-2)));
    zhegv_(itype, jobz, uplo, n,
           a, lda, b, ldb, w,  work,  lwork, rwork.ptr(), info,
           jobzlen, uplo_len);
}



/// These oddly-named wrappers enable the generic syev/heev iterface to get
/// the correct LAPACK routine based upon the argument type.  Internal
/// use only.
STATIC inline void dsyev_(const char* jobz, const char* uplo, integer *n,
                          real4 *a, integer *lda, real4 *w,  real4 *work,  integer *lwork,
                          integer *info, char_len jobzlen, char_len uplo_len ) {
    ssyev_(jobz, uplo, n, a, lda, w,  work,  lwork, info, jobzlen, uplo_len );
}

STATIC void dsyev_(const char* jobz, const char* uplo, integer *n,
                   complex_real4 *a, integer *lda, real4 *w,
                   complex_real4 *work,  integer *lwork,
                   integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<float> rwork(max((integer) 1, (integer) (3* (*n)-2)));
    //std::cout << *n << " " << *lda << " " << *lwork <<std::endl;
    cheev_(jobz, uplo, n, a, lda, w,  work,  lwork, rwork.ptr(),
           info, jobzlen, uplo_len );
}

STATIC void dsyev_(const char* jobz, const char* uplo, integer *n,
                   complex_real8 *a, integer *lda, real8 *w,
                   complex_real8 *work,  integer *lwork,
                   integer *info, char_len jobzlen, char_len uplo_len ) {
    Tensor<double> rwork(max((integer) 1, (integer) (3* (*n)-2)));
    zheev_(jobz, uplo, n, a, lda, w,  work,  lwork, rwork.ptr(),
           info, jobzlen, uplo_len );
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
    void svd(const Tensor<T>& a, Tensor<T>* U,
             Tensor< typename Tensor<T>::scalar_type >* s, Tensor<T>* VT) {
        TENSOR_ASSERT(a.ndim == 2, "svd requires matrix",a.ndim,&a);
        integer m = a.dim[0], n = a.dim[1], rmax = min<integer>(m,n);
        integer lwork = max<integer>(3*min(m,n)+max(m,n),5*min(m,n)-4)*32;
        integer info;
        Tensor<T> A(copy(a)), work(lwork);

        *s = Tensor< typename Tensor<T>::scalar_type >(rmax);
        *U = Tensor<T>(m,rmax);
        *VT = Tensor<T>(rmax,n);

        //std::cout << "n " << n << " m " << m << " lwork " << lwork << std::endl;
        dgesvd_("S","S", &n, &m, A.ptr(), &n, s->ptr(),
                VT->ptr(), &n, U->ptr(), &rmax, work.ptr(), &lwork,
                &info, (char_len) 1, (char_len) 1);

        mask_info(info);

        TENSOR_ASSERT(info == 0, "svd: Lapack failed", info, &a);
    }

    /// Example and test code for interface to LAPACK SVD interfae
    template <typename T>
    double test_svd(int n, int m) {
        Tensor<T> a(n,m), U, VT;
        typedef typename TensorTypeData<T>::scalar_type scalar_type;
        Tensor<scalar_type> s;

        a.fillrandom();
        svd(a,&U,&s,&VT);

        long rank = s.dim[0];
        Tensor<T> b(n,m);
        for (long i=0; i<n; i++)
            for (long j=0; j<m; j++)
                for (long k=0; k<rank; k++)
                    b(i,j) += U(i,k) * T(s(k)) * VT(k,j);

        b -= a;
        return b.absmax();
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
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>* x) {
        TENSOR_ASSERT(a.ndim == 2, "gesv requires matrix",a.ndim,&a);
        integer n = a.dim[0], m = a.dim[1], nrhs = b.dim[1];
        TENSOR_ASSERT(m == n, "gesv requires square matrix",0,&a);
        TENSOR_ASSERT(b.ndim <= 2, "gesv require a vector or matrix for the RHS",b.ndim,&b);
        TENSOR_ASSERT(a.dim[0] == b.dim[0], "gesv matrix and RHS must conform",b.ndim,&b);

        // The input matrix & vectors are destroyed by gesv and we also need Fortran order
        Tensor<T> AT = transpose(a);
        if (b.ndim == 1)
            *x = copy(b);
        else
            *x = transpose(b);

        Tensor<integer> piv(n);
        integer info;

        // note overriding of dgesv for other types above
        dgesv_(&n, &nrhs, AT.ptr(), &n, piv.ptr(), x->ptr(), &n, &info);
        mask_info(info);

        TENSOR_ASSERT((info == 0), "gesv failed", info, &a);

        if (b.ndim == 2) *x = transpose(*x);
    }

    template <typename T>
    double test_gesv(int n, int nrhs) {
        Tensor<T> a(n,n), b1(n), b(n,nrhs), x1, x;

        a.fillrandom();
        b1.fillrandom();
        b.fillrandom();

        gesv(a,b,&x);
        gesv(a,b1,&x1);

        return (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
    }


    /** \brief  Solve Ax = b for general A using the LAPACK *gelss routines.
        
    A should be a square matrix (float, double, float_complex,
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
    */
    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>* x, Tensor< typename Tensor<T>::scalar_type >* s, long *rank) {
        TENSOR_ASSERT(a.ndim == 2, "gelss requires matrix",a.ndim,&a);
        integer n = a.dim[0], m = a.dim[1], nrhs = b.dim[1];
        TENSOR_ASSERT(b.ndim <= 2, "gelss require a vector or matrix for the RHS",b.ndim,&b);
        TENSOR_ASSERT(a.dim[0] == b.dim[0], "gelss matrix and RHS must conform",b.ndim,&b);

        // The input matrix & vectors are destroyed by gelss and we also need Fortran order
        Tensor<T> AT = transpose(a);
        if (b.ndim == 1)
            *x = copy(b);
        else
            *x = transpose(b);

        integer lwork=(3*min(m,n)+max(max(2*min(m,n),max(m,n)),nrhs))*32;
        Tensor<T> work(lwork);
        typedef typename TensorTypeData<T>::scalar_type scalar_type;
        *s = Tensor< scalar_type >(n);
        integer info;
        scalar_type rrcond = rcond;
        integer rrank=0;

        dgelss_(&m, &n, &nrhs, AT.ptr(), &m, x->ptr(), &m, s->ptr(),
                &rrcond, &rrank, work.ptr(), &lwork, &info);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "gelss failed", info, &a);

        *rank = rrank;

        if (b.ndim == 2) *x = transpose(*x);
    }

    template <typename T>
    double test_gelss(int n, int nrhs) {
        Tensor<T> a(n,n), b1(n), b(n,nrhs), x1, x;
        Tensor< typename Tensor<T>::scalar_type > s;
        long rank;

        a.fillrandom();
        b1.fillrandom();
        b.fillrandom();

        gelss(a,b,1e-5,&x,&s,&rank);
        gelss(a,b1,1e-5,&x1,&s,&rank);

        return (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
    }

    /** \brief   Real-symmetric or complex-Hermitian eigenproblem.
        
    A is a real symmetric or complex Hermitian matrix.  Return V and e
    where V is a matrix whose columns are the eigenvectors and e is a
    vector containing the corresponding eigenvalues.  If the LAPACK
    routine fails, it raises a TensorException with value=infor.  The
    input matrix is unchanged.  The eigenvalues are sorted into ascending
    order.  s/dsyev are used for real symmetric matrices; c/zheev are used
    for complex Hermitian.

    The reults will satisfy A*V(_,i) = V(_,i)*e(i).
    */
    template <typename T>
    void syev(const Tensor<T>& A,
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e) {
        TENSOR_ASSERT(A.ndim == 2, "syev requires a matrix",A.ndim,&A);
        TENSOR_ASSERT(A.dim[0] == A.dim[1], "syev requires square matrix",0,&A);
        integer n = A.dim[0];
        integer lwork = max(max((integer) 1,(integer) (3*n-1)),(integer) (34*n));
        integer info;
        Tensor<T> work(lwork);
        *V = transpose(A);		// For Hermitian case
        *e = Tensor<typename Tensor<T>::scalar_type>(n);
        dsyev_("V", "U", &n, V->ptr(), &n, e->ptr(), work.ptr(), &lwork, &info,
               (char_len) 1, (char_len) 1);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "(s/d)syev/(c/z)heev failed", info, &A);
        *V = transpose(*V);
    }

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

    template <typename T>
    double test_syev(int n) {
        Tensor<T> a(n,n), V;
        Tensor< typename Tensor<T>::scalar_type > e;
        a.fillrandom();
        a += madness::my_conj_transpose(a);
        syev(a,&V,&e);
        double err = 0.0;
        for (int i=0; i<n; i++) {
            err = max(err,(double) (inner(a,V(_,i)) - V(_,i)*e(i)).normf());
            //err = max(err,(double) (inner(a,V(_,i)) - V(_,i)*((T) e(i))).normf());
        }
        return err;
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
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e) {
        TENSOR_ASSERT(A.ndim == 2, "sygv requires a matrix",A.ndim,&A);
        TENSOR_ASSERT(A.dim[0] == A.dim[1], "sygv requires square matrix",0,&A);
        TENSOR_ASSERT(B.ndim == 2, "sygv requires a matrix",B.ndim,&A);
        TENSOR_ASSERT(B.dim[0] == B.dim[1], "sygv requires square matrix",0,&A);
        integer ity = itype;
        integer n = A.dim[0];
        integer lwork = max((integer)1,(integer)(3*n-1))*32;
        integer info;
        Tensor<T> work(lwork);
        Tensor<T> b = transpose(B);	// For Hermitian case
        *V = transpose(A);		// For Hermitian case
        *e = Tensor<typename Tensor<T>::scalar_type>(n);
        dsygv_(&ity, "V", "U", &n, V->ptr(), &n, b.ptr(), &n,
               e->ptr(), work.ptr(), &lwork, &info,
               (char_len) 1, (char_len) 1);
        mask_info(info);
        TENSOR_ASSERT(info == 0, "sygv/hegv failed", info, &A);
        *V = transpose(*V);
    }


    template <typename T>
    double test_sygv(int n) {
        Tensor<T> a(n,n), V, b(n,n);
        Tensor< typename Tensor<T>::scalar_type > e;
        a.fillrandom();
        b.fillrandom();
        a += madness::my_conj_transpose(a);
        b += madness::my_conj_transpose(b);
        for (int i=0; i<n; i++) b(i,i) = 2*n;	// To make pos-def
        sygv(a,b,1,&V,&e);
        double err = 0.0;
        for (int i=0; i<n; i++) {
            err = max(err,(double) (inner(a,V(_,i)) - inner(b,V(_,i))*(T) e(i)).normf());
        }
        return err;
    }


    /// Test the Tensor-LAPACK interface ... currently always returns true!
    bool test_tensor_lapack() {
        try {
            cout << "error in float svd " << test_svd<float>(20,30) << endl;
            cout << "error in double svd " << test_svd<double>(30,20) << endl;
            cout << "error in float_complex svd " << test_svd<float_complex>(23,27) << endl;
            cout << "error in double_complex svd " << test_svd<double_complex>(37,19) << endl;
            cout << endl;
            cout << "error in float gesv " << test_gesv<float>(20,30) << endl;
            cout << "error in double gesv " << test_gesv<double>(30,20) << endl;
            cout << "error in float_complex gesv " << test_gesv<float_complex>(23,27) << endl;
            cout << "error in double_complex gesv " << test_gesv<double_complex>(37,19) << endl;
            cout << endl;
            cout << "error in float gelss " << test_gelss<float>(20,30) << endl;
            cout << "error in double gelss " << test_gelss<double>(30,20) << endl;
            cout << "error in float_complex gelss " << test_gelss<float_complex>(23,27) << endl;
            cout << "error in double_complex gelss " << test_gelss<double_complex>(37,19) << endl;
            cout << endl;
            cout << "error in float syev " << test_syev<float>(21) << endl;
            cout << "error in double syev " << test_syev<double>(21) << endl;
            cout << "error in float_complex syev " << test_syev<float_complex>(21) << endl;
            cout << "error in double_complex syev " << test_syev<double_complex>(21) << endl;
            cout << endl;
            cout << "error in float sygv " << test_sygv<float>(21) << endl;
            cout << "error in double sygv " << test_sygv<double>(22) << endl;
            cout << "error in float_complex sygv " << test_sygv<float_complex>(23) << endl;
            cout << "error in double_complex sygv " << test_sygv<double_complex>(24) << endl;
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

}
