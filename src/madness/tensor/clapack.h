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

  
#ifndef MADNESS_LINALG_CLAPACK_H__INCLUDED
#define MADNESS_LINALG_CLAPACK_H__INCLUDED

/// \file clapack.h
/// \brief C++ prototypes for Fortran LAPACK with associated typedefs and macos

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

#  define ssyev_ ssyev
#  define dsyev_ dsyev
#  define cheev_ cheev
#  define zheev_ zheev

#ifndef MADNESS_HAS_ELEMENTAL
#  define ssygv_ ssygv
#  define dsygv_ dsygv
#  define chegv_ chegv
#  define zhegv_ zhegv
#endif

#  define dpotrf_ dpotrf
#  define dgetrf_ dgetrf
#  define dgetri_ dgetri

#  define dtrsm_ dtrsm

#  define dlamch_ dlamch
#  define slamch_ slamch
#else
  // only lowercase with zero and one underscores are handled -- if detected another convention complain loudly
#  ifndef FORTRAN_LINKAGE_LCU
#    error "clapack.h does not support the current Fortran symbol convention -- please, edit and check in the changes."
#  endif
#endif

extern "C"
    double dlamch_(const char* mode, int modelen);

extern "C"
    float slamch_(const char* mode, int modelen);


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
// dgeqrf (M, N, A, LDA, TAU, WORK, LWORK, INFO)
//
// DGEQRF computes a QR factorization of a real M-by-N matrix A:
// A = Q * R.

extern "C"
	void sgeqrf_(integer *m, integer *n,
            	 real4 *a, integer *lda, real4 *tau,
            	 real4 *work, integer *lwork, integer *infoOUT);

extern "C"
	void dgeqrf_(integer *m, integer *n,
            	 real8 *a, integer *lda, real8 *tau,
            	 real8 *work, integer *lwork, integer *infoOUT);

//    	dgeqp3(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO );

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

extern "C"
void dpotrf_(const char *uplo, const integer* n, real8 *a, const integer *lda, integer *info, char_len uplo_len);


extern "C"
void dgetrf_(const integer* m, const integer* n, real8 *a, const integer *lda,
        integer* ipiv, integer *info);

extern "C"
void dgetri_(const integer* n, real8 *a, const integer *lda, const integer* ipiv,
        real8 *work, const integer *lwork, integer *info);

extern "C"
void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
            const integer* m, const integer* n, const real8* alpha, 
            const real8* a, const integer* lda, real8* b, const integer* ldb, 
            char_len sidelen, char_len uplolen, char_len transalen, char_len diaglen);

//			SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
extern "C"
void dtrtri_(const char* uplo, const char* diag, const integer* n, const real8* a,
            const integer* lda, integer *info);

#endif // MADNESS_LINALG_CLAPACK_H__INCLUDED
