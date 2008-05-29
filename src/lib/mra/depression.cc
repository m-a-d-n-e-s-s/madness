#include <madness_config.h>
#include <tensor/tensor.h>

// <insert whining about compilers and x86 architecture here>

namespace madness {

    void aligned_add(long n, double* RESTRICT a, const double* RESTRICT b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        if (n4) {
#if (defined(X86_32) || defined(X86_64))
            // On core-2 this will give 2 cycles/element - optimal is 1.5
            __asm__ __volatile__ (
                                  ".UGHLOOPXX_99:\n"
                                  "movapd   (%0),%%xmm0; addpd   (%1),%%xmm0; movapd %%xmm0,  (%0);\n"
                                  "movapd 16(%0),%%xmm0; addpd 16(%1),%%xmm0; movapd %%xmm0,16(%0);\n"
                                  "add $32,%0; add $32,%1; sub $4,%2; jnz .UGHLOOPXX_99;\n"
                                  : 
                                  : "r"(a), "r"(b), "r"(n4)
                                  : "0","1","2", "memory");
            a+=n4; b+=n4;
#else
            for (long i=0; i<n4; i+=4,a+=4,b+=4) { 
                a[0] += b[0];
                a[1] += b[1];
                a[2] += b[2];
                a[3] += b[3];
            }
#endif    
        }
        for (long i=0; i<rem; i++) *a++ += *b++;
    }

    void aligned_sub(long n, double* RESTRICT a, const double* RESTRICT b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        if (n4) {
#if (defined(X86_32) || defined(X86_64))
            // On core-2 this will give 2 cycles/element - optimal is 1.5
            __asm__ __volatile__ (
                                  ".UGHLOOPXXX_98:\n"
                                  "movapd   (%0),%%xmm0; subpd   (%1),%%xmm0; movapd %%xmm0,  (%0);\n"
                                  "movapd 16(%0),%%xmm0; subpd 16(%1),%%xmm0; movapd %%xmm0,16(%0);\n"
                                  "add $32,%0; add $32,%1; sub $4,%2; jnz .UGHLOOPXXX_98;\n"
                                  : 
                                  : "r"(a), "r"(b), "r"(n4)
                                  : "0","1","2","memory");
            a+=n4; b+=n4;
#else
            for (long i=0; i<n4; i+=4,a+=4,b+=4) { 
                a[0] -= b[0];
                a[1] -= b[1];
                a[2] -= b[2];
                a[3] -= b[3];
            }
#endif    
        }
        for (long i=0; i<rem; i++) *a++ -= *b++;
    }


    void aligned_add(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b) {
        aligned_add(2*n, (double*) a, (const double*) b);
    }

    void aligned_sub(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b) {
        aligned_sub(2*n, (double*) a, (const double*) b);
    }


}
