#include <tensor/mtxmq.h>
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
        
        // The caller is responsible for checking assumptions!
// #define IS_ODD(n) (n&1)
// #define IS_UNALIGNED(p) (p&7)
//         if (IS_ODD(dimi) || IS_ODD(dimj) || IS_ODD(dimk) ||
//             IS_UNALIGNED(a) || IS_UNALIGNED(b) || IS_UNALIGNED(c)) {
//             // CALL SLOW CODE
//         }

        /* 
           Choice is to unroll i or j 
        */
        
#ifdef OPTERON_TUNE
        if (dimi <= dimj) { /* Based on times from X86_64 Opteron ... an old one */
#elif defined(CORE_DUO_TUNE)
        if (1) { /* Based on times from X86_32 Core Duo ... my laptop */
#else
        if (dimj > 12 || dimi <= dimj) { /* Based on times from X86_64 Core2 */
#endif
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
