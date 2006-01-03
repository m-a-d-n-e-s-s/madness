#ifndef CBLAS_H
#define CBLAS_H

#include <fortran_ctypes.h>

namespace madness {
    
#ifdef USE_BLAS
    
    /// BLAS gemm exists only for float, double, float_complex, double_complex.
    
    template <typename T> void gemm(bool transa, bool transb, 
                                    integer m, integer n, integer k,
                                    T alpha, const T* a, integer lda, const T*b, integer ldb,
                                    T beta, T* c, integer ldc);
    
#endif

}

#endif

