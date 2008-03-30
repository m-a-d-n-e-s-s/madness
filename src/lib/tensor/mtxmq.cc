#include <madness_config.h>
#include <tensor/tensor.h>
#include <tensor/mtxmq.h>

// For x86-32/64 have assembly versions for double precision
// For x86-64 have assembly versions for complex double precision

#if defined(X86_32) || defined(X86_64) 

#ifdef X86_64
extern "C" void mTxm26(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm24(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm22(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm20(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm18(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm16(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm14(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm12(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;
#endif
extern "C" void mTxm10(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm8(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm6(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm4(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm2(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;


#ifdef X86_64
extern "C" void TmTxm26(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm24(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm22(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm20(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm18(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm16(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm14(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm12(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;
#endif
extern "C" void TmTxm10(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm8(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm6(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm4(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm2(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

namespace madness {
    
    template<>
    void mTxmq(const long dimi, const long dimj, const long dimk,
               double* RESTRICT c, const double* a, const double* b) {
        //std::cout << "IN DOUBLE ASM VERSION " << dimi << " " << dimj << " " << dimk << "\n";

        
#define IS_ODD(n) (n&1)
#define IS_UNALIGNED(p) (((unsigned long) p)&7)
        if (IS_ODD(dimi) || IS_ODD(dimj) || IS_ODD(dimk) ||
            IS_UNALIGNED(a) || IS_UNALIGNED(b) || IS_UNALIGNED(c)) {
            //std::cout << "slow\n";
            // CALL SLOW CODE
            for (long i=0; i<dimi; i++,c+=dimj,a++) {
                for (long j=0; j<dimj; j++) c[j] = 0.0;
                const double *ai = a;
                for (long k=0; k<dimk; k++,ai+=dimi) {
                    double aki = *ai;
                    for (long j=0; j<dimj; j++) {
                        c[j] += aki*b[k*dimj+j];
                    }
                }
            }
            return;
        }
        
        /* 
           Choice is to unroll i or j 
        */
        
#ifdef OPTERON_TUNE
        bool test = dimi <= dimj; /* Based on times from X86_64 Opteron ... an old one */
#elif defined(CORE_DUO_TUNE)
        bool test = true; /* Based on times from X86_32 Core Duo ... my old laptop */
#elif (defined(CORE2_TUNE) && defined(X86_32))
        bool test = false; /* Based on times from Core2 running in 32-bit mode ... a sad thing */
#elif (defined(CORE2_TUNE) && defined(X86_64))
        bool test = dimj > 12 || dimi <= dimj; /* Based on times from X86_64 Core2 */
#else
        bool test = dimj > 12 || dimi <= dimj; /* Based on times from X86_64 Core2 */
#endif
        if (test) {
            long nj = dimj;
            do {
#ifdef X86_64
                long numj = (nj>26) ? 26 : nj;
#else
                long numj = (nj>10) ? 10 : nj;
#endif
                
                switch (numj) {
#ifdef X86_64
                case 26:
                    TmTxm26(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 24:
                    TmTxm24(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 22:
                    TmTxm22(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 20:
                    TmTxm20(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 18:
                    TmTxm18(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 16:
                    TmTxm16(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 14:
                    TmTxm14(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 12:
                    TmTxm12(dimj, dimi, dimk, c, b, a) ;
                    break; 
#endif
                    
                case 10:
                    TmTxm10(dimj, dimi, dimk, c, b, a) ;
                    break; 
                    
                case 8:
                    TmTxm8(dimj, dimi, dimk, c, b, a) ;
                    break;
                    
                case 6:
                    TmTxm6(dimj, dimi, dimk, c, b, a) ;
                    break;
                    
                case 4:
                    TmTxm4(dimj, dimi, dimk, c, b, a) ;
                    break;
                    
                case 2:
                    TmTxm2(dimj, dimi, dimk, c, b, a) ;
                    break;
                    
                default:
                    throw "mtxmq_byj: should not be here";
                    
                }
                nj -= numj;
                c += numj;
                b += numj;
            } while (nj);
        }
        else {
            long ni = dimi;
            do {
#ifdef X86_64
                long numi = (ni>26) ? 26 : ni;
#else
                long numi = (ni>10) ? 10 : ni;
#endif
                
                switch (numi) {
#ifdef X86_64
                case 26:
                    mTxm26(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 24:
                    mTxm24(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 22:
                    mTxm22(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 20:
                    mTxm20(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 18:
                    mTxm18(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 16:
                    mTxm16(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 14:
                    mTxm14(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 12:
                    mTxm12(dimi, dimj, dimk, c, a, b) ;
                    break; 
#endif
                    
                case 10:
                    mTxm10(dimi, dimj, dimk, c, a, b) ;
                    break; 
                    
                case 8:
                    mTxm8(dimi, dimj, dimk, c, a, b) ;
                    break;
                    
                case 6:
                    mTxm6(dimi, dimj, dimk, c, a, b) ;
                    break;
                    
                case 4:
                    mTxm4(dimi, dimj, dimk, c, a, b) ;
                    break;
                    
                case 2:
                    mTxm2(dimi, dimj, dimk, c, a, b) ;
                    break;
                    
                default:
                    throw "mtxmq: should not be here!";
                }
                ni -= numi;
                c += numi*dimj;
                a += numi;
            } while (ni);
            
        }
    }
}

#else

// On other platforms used to punt to the system blas but this was so slow that
// now we use the generic code which itself could be improved
namespace madness {
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double* RESTRICT c, const double* a, const double* b) {
        for (long i=0; i<dimi; i++,c+=dimj,a++) {
            for (long j=0; j<dimj; j++) c[j] = 0.0;
            const double *ai = a;
            for (long k=0; k<dimk; k++,ai+=dimi) {
                const double aki = *ai;
                for (long j=0; j<dimj; j++) {
                    c[j] += aki*b[k*dimj+j];
                }
            }
        }
        // double zero(0.0);
        // double one(1.0);
        // gemm(false, true, dimj, dimi, dimk, one, b, dimj, a, dimi, zero, c, dimj);
    }
}
#endif

#if defined(X86_64) 
namespace madness {
    template <>
    void mTxmq(const long dimi, const long dimj, const long dimk,
               double_complex* RESTRICT c, const double_complex* a, const double_complex* b) {

        const long dimi16 = dimi<<4;
        const long dimj16 = dimj<<4;

#define ZERO(c) "pxor " #c "," #c ";\n"
#define LOADA   "movddup  (%%r9), %%xmm0; mov %%r10,%%r8; movddup 8(%%r9), %%xmm1; add %2,%%r9; add %3,%%r10; prefetcht0 (%%r9);\n"
#define ENTRY(loop) "mov %0,%%r9; mov %1, %%r10; mov %4,%%r11;.align 16;"#loop": "
#define DOIT(c) "movaps (%%r8),%%xmm2; movaps %%xmm2,%%xmm3; mulpd %%xmm0,%%xmm2; addpd %%xmm2,"#c"; shufpd $1,%%xmm3,%%xmm3; mulpd %%xmm1,%%xmm3; addsubpd %%xmm3, "#c"; \n"
#define NEXT(loop) "sub $1,%%r11; jnz "#loop";"
#define STORE(c) "movaps " #c ", (%%r8); add $16,%%r8;\n"
#define INCB    "add $16,%%r8;\n"

        const long jtile = 12;
        const double_complex* asave = a;
        for (long jlo=0; jlo<dimj; jlo+=jtile,c+=jtile,b+=jtile) {
            long nj = std::min(dimj-jlo,jtile);
            double_complex* RESTRICT ci = c;
            a = asave;
            for (long i=0; i<dimi; i++,ci+=dimj,a++) {
                const double_complex *ai = a;
                const double_complex *bk = b;
                switch (nj) {
                case 1:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 

                                          ENTRY(.KLOOP1)
                                          LOADA
                                          DOIT(%%xmm4)
                                          NEXT(.KLOOP1)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 2:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 

                                          ENTRY(.KLOOP2)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)
                                          NEXT(.KLOOP2)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 3:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 

                                          ENTRY(.KLOOP3)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)
                                          NEXT(.KLOOP3)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 4:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 

                                          ENTRY(.KLOOP4)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)
                                          NEXT(.KLOOP4)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 5:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 

                                          ENTRY(.KLOOP5)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)
                                          NEXT(.KLOOP5)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 6:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 

                                          ENTRY(.KLOOP6)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)
                                          NEXT(.KLOOP6)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );

                    break;

                case 7:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 

                                          ENTRY(.KLOOP7)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10)
                                          NEXT(.KLOOP7)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 8:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 
                                          ZERO(%%xmm11)

                                          ENTRY(.KLOOP8)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11)
                                          NEXT(.KLOOP8)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 9:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)

                                          ENTRY(.KLOOP9)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12)
                                          NEXT(.KLOOP9)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 10:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)

                                          ENTRY(.KLOOP10)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13)
                                          NEXT(.KLOOP10)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 11:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)
                                          ZERO(%%xmm14)

                                          ENTRY(.KLOOP11)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13) INCB
                                          DOIT(%%xmm14)
                                          NEXT(.KLOOP11)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          STORE(%%xmm14)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 12:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4) 
                                          ZERO(%%xmm5) 
                                          ZERO(%%xmm6) 
                                          ZERO(%%xmm7) 
                                          ZERO(%%xmm8) 
                                          ZERO(%%xmm9) 
                                          ZERO(%%xmm10) 
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)
                                          ZERO(%%xmm14)
                                          ZERO(%%xmm15)

                                          ENTRY(.KLOOP12)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13) INCB
                                          DOIT(%%xmm14) INCB
                                          DOIT(%%xmm15)
                                          NEXT(.KLOOP12)

                                          "mov %5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          STORE(%%xmm14)
                                          STORE(%%xmm15)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;
                }
            }
        }
    }
}
#else

// On other platforms used to punt to the system blas but this was so slow that
// now we use the generic code which itself could be improved
namespace madness {
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double_complex* RESTRICT c, const double_complex* a, const double_complex* b) {
        for (long i=0; i<dimi; i++,c+=dimj,a++) {
            for (long j=0; j<dimj; j++) c[j] = 0.0;
            const double_complex *ai = a;
            for (long k=0; k<dimk; k++,ai+=dimi) {
                const double_complex aki = *ai;
                for (long j=0; j<dimj; j++) {
                    c[j] += aki*b[k*dimj+j];
                }
            }
            // double_complex zero(0.0,0.0);
            // double_complex one(1.0,0.0);
            // gemm(false, true, dimj, dimi, dimk, one, b, dimj, a, dimi, zero, c, dimj);
        }
        
    }
}
#endif





