#ifndef CLAPACK_H
#define CLAPACK_H

/// \file clapack.h
/// \brief C++ prototypes for Fortran LAPACK with associated typedefs and macos

#include <fortran_ctypes.h>
#include <linalg/lapack_functions.h>

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

#endif
