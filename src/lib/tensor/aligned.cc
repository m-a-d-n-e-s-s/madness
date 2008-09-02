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

    void aligned_sub(long n, double* restrict a, const double* restrict b) {
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

}
