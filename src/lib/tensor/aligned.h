#ifndef MAD_ALIGNED_H
#define MAD_ALIGNED_H

#include <madness_config.h>
#include <tensor/tensor.h>

namespace madness {
    template <typename T>
    void aligned_zero(long n, T* a) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4) { 
            a[0] = 0;
            a[1] = 0;
            a[2] = 0;
            a[3] = 0;
        }
        for (long i=0; i<rem; i++) *a++ = 0;
    }


#if (defined(X86_32) || defined(X86_64))
  template <> 
  void aligned_zero<double>(long n, double* a);
  
  template <> 
  void aligned_zero<double_complex>(long n, double_complex* a);
#endif

  void aligned_add(long n, double* RESTRICT a, const double* RESTRICT b);
  void aligned_sub(long n, double* RESTRICT a, const double* RESTRICT b);
  void aligned_add(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b);
  void aligned_sub(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b);
}

#endif
