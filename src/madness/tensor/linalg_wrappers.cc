#include <vector>
#include <madness/madness_config.h>
#include <madness/tensor/linalg_wrappers.h>

#ifdef MADNESS_LINALG_USE_LAPACKE
using madness::lapacke::to_cptr;
using madness::lapacke::to_zptr;
#endif


namespace madness {

    template <>
    void svd( char jobu, char jobvt, integer m, integer n, real4* A, integer lda,
              real4* S, real4* U, integer ldu, real4* VT, integer ldvt ) { 
    
    
        integer lwork = -1;
        integer info;
    
        real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        sgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, &info );
    #else
        sgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        sgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, &info );
    #else
        sgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, &info, sizeof(char), sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "SVD Failed", info);
    
    }
    
    template <>
    void svd( char jobu, char jobvt, integer m, integer n, real8* A, integer lda,
              real8* S, real8* U, integer ldu, real8* VT, integer ldvt ) { 
    
    
        integer lwork = -1;
        integer info;
    
        real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        dgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, &info );
    #else
        dgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        dgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, &info );
    #else
        dgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, &info, sizeof(char), sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "SVD Failed", info);
    }
    
    
    
    
    template <>
    void svd( char jobu, char jobvt, integer m, integer n, complex_real4* A, 
              integer lda, real4* S, complex_real4* U, integer ldu, complex_real4* VT, 
              integer ldvt ) { 
    
    
        integer lwork = -1;
        integer lrwork = 5 * std::min(m,n);
        integer info;
    
        std::vector<real4> rwork( lrwork );
    
        complex_real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        cgesvd_( &jobu, &jobvt, &m, &n, to_cptr(A), &lda, S, to_cptr(U), &ldu, to_cptr(VT), &ldvt, to_cptr(&lwork_dummy),
                 &lwork, rwork.data(), &info );
    #else
        cgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        cgesvd_( &jobu, &jobvt, &m, &n, to_cptr(A), &lda, S, to_cptr(U), &ldu, to_cptr(VT), &ldvt, to_cptr(work.data()),
                 &lwork, rwork.data(), &info );
    #else
        cgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "SVD Failed", info);
    }
    
    template <>
    void svd( char jobu, char jobvt, integer m, integer n, complex_real8* A, 
              integer lda, real8* S, complex_real8* U, integer ldu, complex_real8* VT, 
              integer ldvt ) { 
    
    
        integer lwork = -1;
        integer lrwork = 5 * std::min(m,n);
        integer info;
    
        std::vector<real8> rwork( lrwork );
    
        complex_real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        zgesvd_( &jobu, &jobvt, &m, &n, to_zptr(A), &lda, S, to_zptr(U), &ldu, to_zptr(VT), &ldvt, to_zptr(&lwork_dummy),
                 &lwork, rwork.data(), &info );
    #else
        zgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, &lwork_dummy,
                 &lwork, rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        zgesvd_( &jobu, &jobvt, &m, &n, to_zptr(A), &lda, S, to_zptr(U), &ldu, to_zptr(VT), &ldvt, to_zptr(work.data()),
                 &lwork, rwork.data(), &info );
    #else
        zgesvd_( &jobu, &jobvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work.data(),
                 &lwork, rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "SVD Failed", info);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    template <>
    void hereig( char jobz, char uplo, integer n, real4* A, integer lda, real4* W ) {
    
        integer lwork = -1;
        integer info;
    
        real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        ssyev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, &info );
    #else
        ssyev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, &info, 
                sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        ssyev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, &info );
    #else
        ssyev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, &info, 
                sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig( char jobz, char uplo, integer n, real8* A, integer lda, real8* W ) {
    
        integer lwork = -1;
        integer info;
    
        real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        dsyev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, &info );
    #else
        dsyev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, &info, 
                sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        dsyev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, &info );
    #else
        dsyev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, &info, 
                sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig( char jobz, char uplo, integer n, complex_real4* A, integer lda, 
      real4* W ) {
    
        integer lwork = -1;
        integer lrwork = std::max(integer(1), 3*n-2);
        integer info;
    
        std::vector<real4> rwork( lrwork );
    
        complex_real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        cheev_( &jobz, &uplo, &n, to_cptr(A), &lda, W, to_cptr(&lwork_dummy), &lwork, rwork.data(), &info );
    #else
        cheev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, rwork.data(), &info, 
                sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        cheev_( &jobz, &uplo, &n, to_cptr(A), &lda, W, to_cptr(work.data()), &lwork, rwork.data(), &info );
    #else
        cheev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, rwork.data(), &info, 
                sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig( char jobz, char uplo, integer n, complex_real8* A, integer lda, 
      real8* W ) {
    
        integer lwork = -1;
        integer lrwork = std::max(integer(1), 3*n-2);
        integer info;
    
        std::vector<real8> rwork( lrwork );
    
        complex_real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        zheev_( &jobz, &uplo, &n, to_zptr(A), &lda, W, to_zptr(&lwork_dummy), &lwork, rwork.data(), &info );
    #else
        zheev_( &jobz, &uplo, &n, A, &lda, W, &lwork_dummy, &lwork, rwork.data(), &info, 
                sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        zheev_( &jobz, &uplo, &n, to_zptr(A), &lda, W, to_zptr(work.data()), &lwork, rwork.data(), &info );
    #else
        zheev_( &jobz, &uplo, &n, A, &lda, W, work.data(), &lwork, rwork.data(), &info, 
                sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    template <>
    void hereig_gen( integer itype, char jobz, char uplo, integer n, real4* A, 
                     integer lda, real4* B, integer ldb, real4* W ) {
    
        integer lwork = -1;
        integer info;
    
        real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        ssygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                &info );
    #else
        ssygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        ssygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                &info );
    #else
        ssygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                &info, sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig_gen( integer itype, char jobz, char uplo, integer n, real8* A, 
                     integer lda, real8* B, integer ldb, real8* W ) {
    
        integer lwork = -1;
        integer info;
    
        real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        dsygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                &info );
    #else
        dsygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy);
    
        std::vector<real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        dsygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                &info );
    #else
        dsygv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                &info, sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig_gen( integer itype, char jobz, char uplo, integer n, complex_real4* A, 
                     integer lda, complex_real4* B, integer ldb, real4* W ) {
    
        integer lwork = -1;
        integer lrwork = std::max(integer(1), 3*n-2);
        integer info;
    
        std::vector<real4> rwork( lrwork );
    
        complex_real4 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        chegv_( &itype, &jobz, &uplo, &n, to_cptr(A), &lda, to_cptr(B), &ldb, W, to_cptr(&lwork_dummy), &lwork, 
                rwork.data(), &info );
    #else
        chegv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real4> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        chegv_( &itype, &jobz, &uplo, &n, to_cptr(A), &lda, to_cptr(B), &ldb, W, to_cptr(work.data()), &lwork, 
                rwork.data(), &info );
    #else
        chegv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    template <>
    void hereig_gen( integer itype, char jobz, char uplo, integer n, complex_real8* A, 
                     integer lda, complex_real8* B, integer ldb, real8* W ) {
    
        integer lwork = -1;
        integer lrwork = std::max(integer(1), 3*n-2);
        integer info;
      
        std::vector<real8> rwork( lrwork );
      
        complex_real8 lwork_dummy;
    
    #if MADNESS_LINALG_USE_LAPACKE
        zhegv_( &itype, &jobz, &uplo, &n, to_zptr(A), &lda, to_zptr(B), &ldb, W, to_zptr(&lwork_dummy), &lwork, 
                rwork.data(), &info );
    #else
        zhegv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, &lwork_dummy, &lwork, 
                rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
        lwork = integer(lwork_dummy.real());
    
        std::vector<complex_real8> work( lwork );
    
    #if MADNESS_LINALG_USE_LAPACKE
        zhegv_( &itype, &jobz, &uplo, &n, to_zptr(A), &lda, to_zptr(B), &ldb, W, to_zptr(work.data()), &lwork, 
                rwork.data(), &info );
    #else
        zhegv_( &itype, &jobz, &uplo, &n, A, &lda, B, &ldb, W, work.data(), &lwork, 
                rwork.data(), &info, sizeof(char), sizeof(char) );
    #endif
    
    
        LINALG_ASSERT( (info==0), "EVP Failed", info);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    template <>
    void cholesky( char uplo, integer n, real4* A, integer lda ) {
    
        integer info;
    
    #if MADNESS_LINALG_USE_LAPACKE
        spotrf_( &uplo, &n, A, &lda, &info );
    #else
        spotrf_( &uplo, &n, A, &lda, &info, sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "Cholesky Failed", info);
    }
    
    template <>
    void cholesky( char uplo, integer n, real8* A, integer lda ) {
    
        integer info;
    
    #if MADNESS_LINALG_USE_LAPACKE
        dpotrf_( &uplo, &n, A, &lda, &info );
    #else
        dpotrf_( &uplo, &n, A, &lda, &info, sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "Cholesky Failed", info);
    }
    
    template <>
    void cholesky( char uplo, integer n, complex_real4* A, integer lda ) {
    
        integer info;
    
    #if MADNESS_LINALG_USE_LAPACKE
        cpotrf_( &uplo, &n, to_cptr(A), &lda, &info );
    #else
        cpotrf_( &uplo, &n, A, &lda, &info, sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "Cholesky Failed", info);
    }
    
    template <>
    void cholesky( char uplo, integer n, complex_real8* A, integer lda ) {
    
        integer info;
    
    #if MADNESS_LINALG_USE_LAPACKE
        zpotrf_( &uplo, &n, to_zptr(A), &lda, &info );
    #else
        zpotrf_( &uplo, &n, A, &lda, &info, sizeof(char) );
    #endif
    
        LINALG_ASSERT( (info==0), "Cholesky Failed", info);
    }



}
