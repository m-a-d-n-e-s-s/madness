#include <madness_config.h>


#include <tensor/aligned.h>
#include <tensor/tensor.h>
namespace madness {
#if (defined(X86_32) || defined(X86_64))
    template <>
    void aligned_zero<double>(long n, double* a) {
        if ((((unsigned long) a) & 0x0f)) throw "WTF";
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        if (n4) {
            //std::cout << "entering asm " << (void *) a << " " << n4 << std::endl;
            __asm__ __volatile__ (
                                  "pxor %%xmm0,%%xmm0;\n"
                                  ".UGHLOOP_47:\n"
                                  "movapd   %%xmm0,  (%0);\n"
                                  "movapd   %%xmm0,16(%0);\n"
                                  "add $32,%0; sub $4,%1; jnz .UGHLOOP_47;\n"
                                  :
                                  : "r"(a), "r"(n4)
                                  : "0","1","xmm0", "memory");
            //std::cout << "leaving asm " << (void *) a << " " << n4 << std::endl;
            a+=n4;
        }
        for (long i=0; i<rem; i++) *a++ = 0.0;
    }
#endif

    void aligned_add(long n, double* restrict a, const double* restrict b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        if (n4) {
#if ( (!defined(ON_A_MAC)) && (defined(X86_32) || defined(X86_64)) )
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


//#if (defined(X86_32) || defined(X86_64))
//    template <>
//    void aligned_axpy(long n, double* restrict a, const double* restrict b, double s) {
//        // In practice this is not all that much faster than the code generated
//        // by g++ 4.4.0 from the 4-way unrolled source code ... have not timed
//        // things to see if this is g++ doing well or me doing badly.
//        long n16 = (n>>4)<<4;
//        long rem = n-n16;
//        if (n16) {
//            __asm__ __volatile__ ("movddup (%3),%%xmm1;\n"
//                                  ".UGHLOOPRR_99:\n"
//                                  "movapd    (%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd    (%0),%%xmm0; movapd %%xmm0,   (%0);\n"
//                                  "movapd  16(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  16(%0),%%xmm0; movapd %%xmm0, 16(%0);\n"
//                                  "movapd  32(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  32(%0),%%xmm0; movapd %%xmm0, 32(%0);\n"
//                                  "movapd  48(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  48(%0),%%xmm0; movapd %%xmm0, 48(%0);\n"
//                                  "movapd  64(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  64(%0),%%xmm0; movapd %%xmm0, 64(%0);\n"
//                                  "movapd  80(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  80(%0),%%xmm0; movapd %%xmm0, 80(%0);\n"
//                                  "movapd  96(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd  96(%0),%%xmm0; movapd %%xmm0, 96(%0);\n"
//                                  "movapd 112(%1),%%xmm0; mulpd %%xmm1,%%xmm0; addpd 112(%0),%%xmm0; movapd %%xmm0,112(%0);\n"
//                                  "add $128,%0; add $128,%1; sub $16,%2; jnz .UGHLOOPRR_99;\n"
//                                  :
//                                  : "r"(a), "r"(b), "r"(n16), "r"(&s)
//                                  : "0","1","2", "memory");
//        }
//        for (long i=0; i<rem; i++) *a++ += s * *b++;
//    }
//#endif

    void aligned_sub(long n, double* restrict a, const double* restrict b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        if (n4) {
#if ( (!defined(ON_A_MAC)) && (defined(X86_32) || defined(X86_64)) )
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

}
