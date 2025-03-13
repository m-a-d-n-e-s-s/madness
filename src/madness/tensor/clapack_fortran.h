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

  
#ifndef MADNESS_LINALG_CLAPACK_FORTRAN_H__INCLUDED
#define MADNESS_LINALG_CLAPACK_FORTRAN_H__INCLUDED

/// \file clapack_fortran.h
/// \brief Legacy C++ prototypes for Fortran LAPACK with associated typedefs and macos

#include <madness/fortran_ctypes.h>

#ifdef FORTRAN_LINKAGE_LC
#  define sgesvd_ sgesvd
#  define dgesvd_ dgesvd
#  define cgesvd_ cgesvd
#  define zgesvd_ zgesvd

#  define sgesv_ sgesv
#  define dgesv_ dgesv
#  define cgesv_ cgesv
#  define zgesv_ zgesv

#  define sgelss_ sgelss
#  define dgelss_ dgelss
#  define cgelss_ cgelss
#  define zgelss_ zgelss

#  define sgels_ sgels
#  define dgels_ dgels
#  define cgels_ cgels
#  define zgels_ zgels

#  define ssyev_ ssyev
#  define dsyev_ dsyev
#  define cheev_ cheev
#  define zheev_ zheev

#  define sggev_ sggev
#  define dggev_ dggev
#  define cggev_ cggev
#  define zggev_ zggev

#ifndef MADNESS_HAS_ELEMENTAL
#  define ssygv_ ssygv
#  define dsygv_ dsygv
#  define chegv_ chegv
#  define zhegv_ zhegv
#endif

#  define spotrf_ spotrf
#  define cpotrf_ cpotrf
#  define dpotrf_ dpotrf
#  define zpotrf_ zpotrf

#  define sgetrf_ sgetrf
#  define cgetrf_ cgetrf
#  define dgetrf_ dgetrf
#  define zgetrf_ zgetrf

#  define sgetri_ sgetri
#  define cgetri_ cgetri
#  define dgetri_ dgetri
#  define zgetri_ zgetri

#  define strsm_ strsm
#  define ctrsm_ ctrsm
#  define dtrsm_ dtrsm
#  define ztrsm_ ztrsm

#  define dlamch_ dlamch
#  define slamch_ slamch

#  define sgeev_ sgeev
#  define cgeev_ cgeev
#  define dgeev_ dgeev
#  define zgeev_ zgeev

#else
  // only lowercase with zero and one underscores are handled -- if detected another convention complain loudly
#  ifndef FORTRAN_LINKAGE_LCU
#    error "clapack.h does not support the current Fortran symbol convention -- please, edit and check in the changes."
#  endif
#endif

// SUBROUTINE DLAMCH( CMACH, RESULT )

// PURPOSE
//     DLAMCH determines double precision machine parameters.

extern "C"
    float slamch_(const char* mode, int modelen);
extern "C"
    double dlamch_(const char* mode, int modelen);

// SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

// PURPOSE
//     DGESVD computes the singular value decomposition (SVD) of a real
//     M-by-N matrix A, optionally computing the left and/or right singular
//     vectors.

extern "C"
    void sgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
                 real4 *a, integer *lda, real4 *s, real4 *u, integer *ldu,
                 real4 *vt, integer *ldvt, real4 *work, integer *lwork,
                 integer *info, char_len jobulen, char_len jobvtlen);
extern "C"
    void dgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
                 real8 *a, integer *lda, real8 *s, real8 *u, integer *ldu,
                 real8 *vt, integer *ldvt, real8 *work, integer *lwork,
                 integer *info, char_len jobulen, char_len jobvtlen);
extern "C"
    void cgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
                 complex_real4 *a, integer *lda, real4 *s, complex_real4 *u,
                 integer *ldu, complex_real4 *vt, integer *ldvt, complex_real4 *work,
                 integer *lwork, real4 *rwork,
                 integer *info, char_len jobulen, char_len jobvtlen);
extern "C"
    void zgesvd_(const char *jobu, const char *jobvt, integer *m, integer *n,
                 complex_real8 *a, integer *lda, real8 *s, complex_real8 *u,
                 integer *ldu, complex_real8 *vt, integer *ldvt, complex_real8 *work,
                 integer *lwork, real8 *rwork,
                 integer *info, char_len jobulen, char_len jobvtlen);

// SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

// PURPOSE
//     DGESV computes the solution to a real system of linear equations
//        A * X = B,
//     where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

extern "C"
    void sgesv_(integer* n, integer* nrhs, real4* AT, integer* lda,
                integer* piv, real4* x, integer* ldx, integer* info);
extern "C"
    void dgesv_(integer* n, integer* nrhs, real8* AT, integer* lda,
                integer* piv, real8* x, integer* ldx, integer* info);
extern "C"
    void cgesv_(integer* n, integer* nrhs, complex_real4* AT, integer* lda,
                integer* piv, complex_real4* x, integer* ldx, integer* info);
extern "C"
    void zgesv_(integer* n, integer* nrhs, complex_real8* AT, integer* lda,
                integer* piv, complex_real8* x, integer* ldx, integer* info);

// SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )

// PURPOSE
//     DGELSS computes the minimum norm solution to a real linear least
//     squares problem:
//          Minimize 2-norm(| b - A*x |)
//     using the singular value decomposition (SVD) of A. A is an M-by-N
//     matrix which may be rank-deficient.

extern "C"
    void sgelss_(integer *m, integer *n, integer *nrhs,
                 real4 *a, integer *lda, real4 *b, integer *ldb, real4 *sOUT,
                 real4 *rcondIN, integer *rankOUT, real4 *work,
                 integer *lwork, integer *infoOUT);
extern "C"
    void dgelss_(integer *m, integer *n, integer *nrhs,
                 real8 *a, integer *lda, real8 *b, integer *ldb, real8 *sOUT,
                 real8 *rcondIN, integer *rankOUT, real8 *work,
                 integer *lwork, integer *infoOUT);
extern "C"
    void cgelss_(integer *m, integer *n, integer *nrhs,
                 complex_real4 *a, integer *lda, complex_real4 *b, integer *ldb,
                 real4 *sOUT,
                 real4 *rcondIN, integer *rankOUT, complex_real4 *work,
                 integer *lwork, real4 *rwork, integer *infoOUT);
extern "C"
    void zgelss_(integer *m, integer *n, integer *nrhs,
                 complex_real8 *a, integer *lda, complex_real8 *b, integer *ldb,
                 real8 *sOUT,
                 real8 *rcondIN, integer *rankOUT, complex_real8 *work,
                 integer *lwork, real8 *rwork, integer *infoOUT);

// SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

// PURPOSE
//     DGELS solves overdetermined or underdetermined real linear systems
//     involving an M-by-N matrix A, or its transpose, using a QR or LQ
//     factorization of A.  It is assumed that A has full rank.

extern "C"
    void sgels_(const char *trans, integer *m, integer *n, integer *nrhs,
            real4 *a, integer *lda, real4 *b, integer *ldb, real4 *work,
            integer *lwork, integer *infoOUT, char_len translen);
extern "C"
    void dgels_(const char *trans, integer *m, integer *n, integer *nrhs,
            real8 *a, integer *lda, real8 *b, integer *ldb, real8 *work,
            integer *lwork, integer *infoOUT, char_len translen);
extern "C"
    void cgels_(const char *trans, integer *m, integer *n, integer *nrhs,
            complex_real4 *a, integer *lda, complex_real4 *b, integer *ldb,
            complex_real4 *work,
            integer *lwork, real4 *rwork, integer *infoOUT, char_len translen);
extern "C"
    void zgels_(const char *trans, integer *m, integer *n, integer *nrhs,
            complex_real8 *a, integer *lda, complex_real8 *b, integer *ldb,
            complex_real8 *work,
            integer *lwork, real8 *rwork, integer *infoOUT, char_len translen);

// SUBROUTINE DGGEV( JOBZ, UPLO, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

// PURPOSE
//     DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
//     the generalized eigenvalues, and optionally, the left and/or right
//     generalized eigenvectors.

extern "C"
    void sggev_(const char* jobz, const char* uplo, integer *n,
                real4* a, integer* lda, real4* b, integer* ldb,
                real4* alphar, real4* alphai, real4* beta, 
                real4* vl, integer* ldvl, real4* vr, integer* ldvr,
                real4* work,  integer* lwork, integer* info,
                char_len jobzlen, char_len uplo_len);
extern "C"
     void dggev_(const char* jobl, const char* jobr, integer *n,
		 real8 *a, integer *lda, real8 *b, integer *ldb,
		 real8 *w_real, real8 *w_imag, real8 *beta,
		 real8 *vl, integer *ldvl, real8 *vr, integer *ldvr,
		 real8 *work,  integer *lwork, integer *info,
		 char_len jobzlen, char_len uplo_len);
extern "C"
    void cggev_(const char* jobz, const char* uplo, integer *n,
                complex_real4* a, integer* lda, complex_real4* b, integer* ldb,
                complex_real4* alpha, complex_real4* beta, 
                complex_real4* vl, integer* ldvl, complex_real4* vr, integer* ldvr,
                complex_real4* work,  integer* lwork, real4* rwork, integer* info,
                char_len jobzlen, char_len uplo_len);
extern "C"
    void zggev_(const char* jobz, const char* uplo, integer *n,
                complex_real8* a, integer* lda, complex_real8* b, integer* ldb,
                complex_real8* alpha, complex_real8* beta, 
                complex_real8* vl, integer* ldvl, complex_real8* vr, integer* ldvr,
                complex_real8* work,  integer* lwork, real8* rwork, integer* info,
                char_len jobzlen, char_len uplo_len);

// SUBROUTINE DGEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

// PURPOSE
//     DGEEV computes for an N-by-N real nonsymmetric matrix A, the
//     eigenvalues and, optionally, the left and/or right eigenvectors.

extern "C"
    void sgeev_(const char* jobz, const char* uplo, integer *n, real4* a, integer* lda,
                real4* w_real, real4* w_imag, real4* v, integer* ldv, real4* vr, integer* ldvr,
                real4* work,  integer* lwork, integer* info,
                char_len jobzlen, char_len uplo_len );
extern "C"
    void dgeev_(const char* jobz, const char* uplo, integer *n,
                real8* a, integer* lda, real8* w_real, real8* w_imag, real8* v, integer* ldv,
                real8* vr, integer* ldvr, real8* work,  integer* lwork, integer* info,
                char_len jobzlen, char_len uplo_len );
extern "C"
    void cgeev_(const char* jobz, const char* uplo, integer *n, complex_real4* a, integer* lda,
                complex_real4* w, complex_real4* vl, integer* ldvl, complex_real4* vr, integer* ldvr,
                complex_real4* work,  integer* lwork, real4* rwork, integer* info,
                char_len jobzlen, char_len uplo_len );
extern "C"
    void zgeev_(const char* jobz, const char* uplo, integer *n, complex_real8* a, integer* lda,
                complex_real8* w, complex_real8* vl, integer* ldvl, complex_real8* vr, integer* ldvr,
                complex_real8* work,  integer* lwork, real8* rwork, integer* info,
                char_len jobzlen, char_len uplo_len );

// SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

// PURPOSE
//     DSYEV computes all eigenvalues and, optionally, eigenvectors of a
//     real symmetric matrix A.

extern "C"
    void ssyev_(const char* jobz, const char* uplo, integer *n,
                real4 *a, integer *lda, real4 *w,  real4 *work,  integer *lwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void dsyev_(const char* jobz, const char* uplo, integer *n,
                real8 *a, integer *lda, real8 *w,  real8 *work,  integer *lwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void cheev_(const char* jobz, const char* uplo, integer *n,
                complex_real4 *a, integer *lda, real4 *w,  complex_real4 *work,
                integer *lwork, real4 *rwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void zheev_(const char* jobz, const char* uplo, integer *n,
                complex_real8 *a, integer *lda, real8 *w,  complex_real8 *work,
                integer *lwork, real8 *rwork,
                integer *info, char_len jobzlen, char_len uplo_len );

// SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )
// SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )

// PURPOSE
//     DSYGV computes all the eigenvalues, and optionally, the eigenvectors
//     of a real generalized symmetric-definite eigenproblem, of the form
//     A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.

extern "C"
    void ssygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
                real4 *a, integer *lda, real4 *b, integer *ldb,
                real4 *w,  real4 *work,  integer *lwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void dsygv_(integer *itype, const char* jobz, const char* uplo, integer *n,
                real8 *a, integer *lda, real8 *b, integer *ldb,
                real8 *w,  real8 *work,  integer *lwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void chegv_(integer *itype, const char* jobz, const char* uplo, integer *n,
                complex_real4 *a, integer *lda, complex_real4 *b, integer *ldb,
                real4 *w,  complex_real4 *work,  integer *lwork, real4 *rwork,
                integer *info, char_len jobzlen, char_len uplo_len );
extern "C"
    void zhegv_(integer *itype, const char* jobz, const char* uplo, integer *n,
                complex_real8 *a, integer *lda, complex_real8 *b, integer *ldb,
                real8 *w,  complex_real8 *work,  integer *lwork, real8 *rwork,
                integer *info, char_len jobzlen, char_len uplo_len );

// SUBROUTINE DGEQRF (M, N, A, LDA, TAU, WORK, LWORK, INFO)
//
// PURPOSE
//		DGEQRF computes a QR factorization of a real M-by-N matrix A:
//		A = Q * R.

extern "C"
	void sgeqrf_(integer *m, integer *n,
            	 real4 *a, integer *lda, real4 *tau,
            	 real4 *work, integer *lwork, integer *infoOUT);

extern "C"
	void dgeqrf_(integer *m, integer *n,
            	 real8 *a, integer *lda, real8 *tau,
            	 real8 *work, integer *lwork, integer *infoOUT);

extern "C"
	void cgeqrf_(integer *m, integer *n,
				 complex_real4 *a, integer *lda, complex_real4 *tau,
				 complex_real4 *work, integer *lwork, integer *infoOUT);

extern "C"
	void zgeqrf_(integer *m, integer *n,
				 complex_real8 *a, integer *lda, complex_real8 *tau,
				 complex_real8 *work, integer *lwork, integer *infoOUT);

// SUBROUTINE DGEQP3(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO );

// PURPOSE
//		DGEQP3 computes a QR factorization with column pivoting of a
//		matrix A:  A*P = Q*R  using Level 3 BLAS.

extern "C"
    void sgeqp3_(integer *m, integer *n,
                 real4 *a, integer *lda, integer *jpvt, real4 *tau,
                 real4 *work, integer *lwork, integer *infoOUT);

extern "C"
    void dgeqp3_(integer *m, integer *n,
            	 real8 *a, integer *lda, integer *jpvt, real8 *tau,
            	 real8 *work, integer *lwork, integer *infoOUT);

extern "C"
    void cgeqp3_(integer *m, integer *n, complex_real4 *a,
    			 integer *lda, integer *jpvt, complex_real4 *tau,
    			 complex_real4 *work, integer *lwork, real4 *rwork,
    			 integer *infoOUT);

extern "C"
    void zgeqp3_(integer *m, integer *n, complex_real8 *a,
			 	 integer *lda, integer *jpvt, complex_real8 *tau,
			 	 complex_real8 *work, integer *lwork, real8 *rwork,
			 	 integer *infoOUT);

// SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
// SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

// PURPOSE
//     DORGQR generates an M-by-N real matrix Q with orthonormal
//     columns, which is defined as the first N columns of a pro-
//     duct of K elementary reflectors of order M

extern "C"
    void sorgqr_(integer *m, integer *n, integer *k,
                 real4 *a, integer *lda, real4 *tau,
                 real4 *work, integer *lwork, integer *info);
extern "C"
    void dorgqr_(integer *m, integer *n, integer *k,
            	 real8 *a, integer *lda, real8 *tau,
            	 real8 *work, integer *lwork, integer *info);
extern "C"
    void cungqr_(integer *m, integer *n, integer *k,
    			 complex_real4 *a, integer *lda, complex_real4 *tau,
    			 complex_real4 *work, integer *lwork, integer *info);
extern "C"
    void zungqr_(integer *m, integer *n, integer *k,
    			 complex_real8 *a, integer *lda, complex_real8 *tau,
			 	 complex_real8 *work, integer *lwork, integer *info);

// SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )

// PURPOSE
//     computes the Cholesky factorization of a real symmetric
//     positive definite matrix A.

extern "C"
void spotrf_(const char *uplo, const integer* n, real4 *a, const integer *lda, integer *info, char_len uplo_len);
extern "C"
void dpotrf_(const char *uplo, const integer* n, real8 *a, const integer *lda, integer *info, char_len uplo_len);
extern "C"
void cpotrf_(const char *uplo, const integer* n, complex_real4 *a, const integer *lda, integer *info, char_len uplo_len);
extern "C"
void zpotrf_(const char *uplo, const integer* n, complex_real8 *a, const integer *lda, integer *info, char_len uplo_len);

// SUBROUTINE DPSTRF( UPLO, N, A, LDA, IPIV, RANK, TOL, WORK, INFO )

// PURPOSE
//     DPSTRF computes the Cholesky factorization with complete
//     pivoting of a real symmetric positive semidefinite matrix A.

extern "C"
void spstrf_(const char *uplo, const integer* n, real4 *a, const integer *lda, integer* ipiv, integer* rank, real4* tol,
		real4* work, integer *info);
extern "C"
void dpstrf_(const char *uplo, const integer* n, real8 *a, const integer *lda, integer* ipiv, integer* rank, real8* tol,
		real8* work, integer *info);
extern "C"
void cpstrf_(const char *uplo, const integer* n, complex_real4 *a, const integer *lda, integer* ipiv, integer* rank, real4* tol,
		complex_real4* work, integer *info);
extern "C"
void zpstrf_(const char *uplo, const integer* n, complex_real8 *a, const integer *lda, integer* ipiv, integer* rank, real8* tol,
		complex_real8* work, integer *info);

// SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

// PURPOSE
//     DGETRF computes an LU factorization of a general M-by-N matrix A
//     using partial pivoting with row interchanges.

extern "C"
void sgetrf_(const integer* m, const integer* n, real4 *a, const integer *lda,
        integer* ipiv, integer *info);
extern "C"
void dgetrf_(const integer* m, const integer* n, real8 *a, const integer *lda,
        integer* ipiv, integer *info);
extern "C"
void cgetrf_(const integer* m, const integer* n, complex_real4 *a, const integer *lda,
        integer* ipiv, integer *info);
extern "C"
void zgetrf_(const integer* m, const integer* n, complex_real8 *a, const integer *lda,
        integer* ipiv, integer *info);

// SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

// PURPOSE
//     DGETRI computes the inverse of a matrix using the LU factorization
//     computed by DGETRF.

extern "C"
void sgetri_(const integer* n, real4 *a, const integer *lda, const integer* ipiv,
        real4 *work, const integer *lwork, integer *info);
extern "C"
void dgetri_(const integer* n, real8 *a, const integer *lda, const integer* ipiv,
        real8 *work, const integer *lwork, integer *info);
extern "C"
void cgetri_(const integer* n, complex_real4 *a, const integer *lda, const integer* ipiv,
        complex_real4 *work, const integer *lwork, integer *info);
extern "C"
void zgetri_(const integer* n, complex_real8 *a, const integer *lda, const integer* ipiv,
        complex_real8 *work, const integer *lwork, integer *info);

// SUBROUTINE DTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )

// PURPOSE
//     DTRSM solves one of the matrix equations
//        op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
//     where alpha is a scalar, X and B are m by n matrices, A is a unit, or
//     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
//        op( A ) = A   or   op( A ) = A**T.

extern "C"
void strsm_(const char* side, const char* uplo, const char* transa, const char* diag,
            const integer* m, const integer* n, const real4* alpha,
            const real4* a, const integer* lda, real4* b, const integer* ldb,
            char_len sidelen, char_len uplolen, char_len transalen, char_len diaglen);
extern "C"
void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
            const integer* m, const integer* n, const real8* alpha,
            const real8* a, const integer* lda, real8* b, const integer* ldb,
            char_len sidelen, char_len uplolen, char_len transalen, char_len diaglen);
extern "C"
void ctrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
            const integer* m, const integer* n, const complex_real4* alpha,
            const complex_real4* a, const integer* lda, complex_real4* b, const integer* ldb,
            char_len sidelen, char_len uplolen, char_len transalen, char_len diaglen);
extern "C"
void ztrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
            const integer* m, const integer* n, const complex_real8* alpha,
            const complex_real8* a, const integer* lda, complex_real8* b, const integer* ldb,
            char_len sidelen, char_len uplolen, char_len transalen, char_len diaglen);

// SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )

// PURPOSE
//     DTRTRI computes the inverse of a real upper or lower triangular
//     matrix A.

extern "C"
void strtri_(const char* uplo, const char* diag, const integer* n, const real4* a,
            const integer* lda, integer *info);
extern "C"
void dtrtri_(const char* uplo, const char* diag, const integer* n, const real8* a,
            const integer* lda, integer *info);
extern "C"
void ctrtri_(const char* uplo, const char* diag, const integer* n, const complex_real4* a,
            const integer* lda, integer *info);
extern "C"
void ztrtri_(const char* uplo, const char* diag, const integer* n, const complex_real8* a,
            const integer* lda, integer *info);

#endif // MADNESS_LINALG_CLAPACK_FORTRAN_H__INCLUDED
