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

    template <typename T, typename Q>
    void aligned_axpy(long n, T* restrict a, const T* restrict b, Q s) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4,b+=4) {
            a[0] += s*b[0];
            a[1] += s*b[1];
            a[2] += s*b[2];
            a[3] += s*b[3];
        }
        for (long i=0; i<rem; i++) *a++ += s * *b++;
    }



#if (defined(X86_32) || defined(X86_64))
  template <>
  void aligned_zero<double>(long n, double* a);

  template <>
  void aligned_zero<double_complex>(long n, double_complex* a);

//  template <>
//  void aligned_axpy(long n, double* restrict a, const double* restrict b, double s);
#endif

  void aligned_add(long n, double* restrict a, const double* restrict b);
  void aligned_sub(long n, double* restrict a, const double* restrict b);
  void aligned_add(long n, double_complex* restrict a, const double_complex* restrict b);
  void aligned_sub(long n, double_complex* restrict a, const double_complex* restrict b);

}

#endif
