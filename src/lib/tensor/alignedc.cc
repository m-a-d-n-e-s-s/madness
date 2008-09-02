#include <madness_config.h>

#include <tensor/aligned.h>
#include <tensor/tensor.h>

namespace madness {

#if (defined(X86_32) || defined(X86_64))
    template <> 
    void aligned_zero<double_complex>(long n, double_complex* a) {
      aligned_zero(2*n, (double *) a);
    }
#endif

    void aligned_add(long n, double_complex* restrict a, const double_complex* restrict b) {
        aligned_add(2*n, (double*) a, (const double*) b);
    }

    void aligned_sub(long n, double_complex* restrict a, const double_complex* restrict b) {
        aligned_sub(2*n, (double*) a, (const double*) b);
    }
}


