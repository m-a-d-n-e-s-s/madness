#include <builtins.h>
#include <complex>

typedef std::complex<double> double_complex;

#ifdef HAVE_IBMBGP
namespace madness {
    void bgpmTxmq(long dimi, long dimj, long dimk, double_complex *  c_x,  const double_complex *  a_x,  const double_complex *  b_x) {
        int i, j, k, ii;
        double *  c = (double*)c_x;
         double *  a = (double*)a_x;
         double *  b = (double*)b_x;
        __complex__ double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_4_0, _c_4_1, _c_4_2, _c_4_3, _c_4_4, _c_4_5, _c_4_6, _c_4_7, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7;
        double _a_0_0, _a_0_1, _a_0_2, _a_0_3, _a_0_4;
         double _ai_0_0, _ai_0_1, _ai_0_2, _ai_0_3, _ai_0_4;
        for (i=0; i+5<=dimi ; i+=5) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4*2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                _c_2_6 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                _c_3_4 = __cmplx(0.0,0.0);
                _c_3_6 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_4_2 = __cmplx(0.0,0.0);
                _c_4_4 = __cmplx(0.0,0.0);
                _c_4_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _a_0_3 = *((pa+6));
                    _ai_0_3 = *((pa+6)+1);
                    _a_0_4 = *((pa+8));
                    _ai_0_4 = *((pa+8)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_3_0 = __fxcxnpma(_c_3_0, _b_0_0, _ai_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_4_0 = __fxcxnpma(_c_4_0, _b_0_0, _ai_0_4);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_3_2 = __fxcxnpma(_c_3_2, _b_0_2, _ai_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_4_2 = __fxcxnpma(_c_4_2, _b_0_2, _ai_0_4);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _b_0_4, _a_0_3);
                    _c_3_4 = __fxcxnpma(_c_3_4, _b_0_4, _ai_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _b_0_4, _a_0_4);
                    _c_4_4 = __fxcxnpma(_c_4_4, _b_0_4, _ai_0_4);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _c_1_6 = __fxcxnpma(_c_1_6, _b_0_6, _ai_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _b_0_6, _a_0_2);
                    _c_2_6 = __fxcxnpma(_c_2_6, _b_0_6, _ai_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _b_0_6, _a_0_3);
                    _c_3_6 = __fxcxnpma(_c_3_6, _b_0_6, _ai_0_3);
                    _c_4_6 = __fxcpmadd(_c_4_6, _b_0_6, _a_0_4);
                    _c_4_6 = __fxcxnpma(_c_4_6, _b_0_6, _ai_0_4);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
                __stfpd(xc+(i+2)*dimj*2+6, _c_2_6);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj*2+2, _c_3_2);
                __stfpd(xc+(i+3)*dimj*2+4, _c_3_4);
                __stfpd(xc+(i+3)*dimj*2+6, _c_3_6);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj*2+2, _c_4_2);
                __stfpd(xc+(i+4)*dimj*2+4, _c_4_4);
                __stfpd(xc+(i+4)*dimj*2+6, _c_4_6);
            }
            if (j>3) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                _c_2_6 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                _c_3_4 = __cmplx(0.0,0.0);
                _c_3_6 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_4_2 = __cmplx(0.0,0.0);
                _c_4_4 = __cmplx(0.0,0.0);
                _c_4_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _a_0_3 = *((pa+6));
                    _ai_0_3 = *((pa+6)+1);
                    _a_0_4 = *((pa+8));
                    _ai_0_4 = *((pa+8)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_3_0 = __fxcxnpma(_c_3_0, _b_0_0, _ai_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_4_0 = __fxcxnpma(_c_4_0, _b_0_0, _ai_0_4);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_3_2 = __fxcxnpma(_c_3_2, _b_0_2, _ai_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_4_2 = __fxcxnpma(_c_4_2, _b_0_2, _ai_0_4);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _b_0_4, _a_0_3);
                    _c_3_4 = __fxcxnpma(_c_3_4, _b_0_4, _ai_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _b_0_4, _a_0_4);
                    _c_4_4 = __fxcxnpma(_c_4_4, _b_0_4, _ai_0_4);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _c_1_6 = __fxcxnpma(_c_1_6, _b_0_6, _ai_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _b_0_6, _a_0_2);
                    _c_2_6 = __fxcxnpma(_c_2_6, _b_0_6, _ai_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _b_0_6, _a_0_3);
                    _c_3_6 = __fxcxnpma(_c_3_6, _b_0_6, _ai_0_3);
                    _c_4_6 = __fxcpmadd(_c_4_6, _b_0_6, _a_0_4);
                    _c_4_6 = __fxcxnpma(_c_4_6, _b_0_6, _ai_0_4);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
                __stfpd(xc+(i+2)*dimj*2+6, _c_2_6);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj*2+2, _c_3_2);
                __stfpd(xc+(i+3)*dimj*2+4, _c_3_4);
                __stfpd(xc+(i+3)*dimj*2+6, _c_3_6);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj*2+2, _c_4_2);
                __stfpd(xc+(i+4)*dimj*2+4, _c_4_4);
                __stfpd(xc+(i+4)*dimj*2+6, _c_4_6);
            }
            else if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                _c_3_4 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_4_2 = __cmplx(0.0,0.0);
                _c_4_4 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _a_0_3 = *((pa+6));
                    _ai_0_3 = *((pa+6)+1);
                    _a_0_4 = *((pa+8));
                    _ai_0_4 = *((pa+8)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_3_0 = __fxcxnpma(_c_3_0, _b_0_0, _ai_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_4_0 = __fxcxnpma(_c_4_0, _b_0_0, _ai_0_4);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_3_2 = __fxcxnpma(_c_3_2, _b_0_2, _ai_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_4_2 = __fxcxnpma(_c_4_2, _b_0_2, _ai_0_4);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _b_0_4, _a_0_3);
                    _c_3_4 = __fxcxnpma(_c_3_4, _b_0_4, _ai_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _b_0_4, _a_0_4);
                    _c_4_4 = __fxcxnpma(_c_4_4, _b_0_4, _ai_0_4);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj*2+2, _c_3_2);
                __stfpd(xc+(i+3)*dimj*2+4, _c_3_4);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj*2+2, _c_4_2);
                __stfpd(xc+(i+4)*dimj*2+4, _c_4_4);
            }
            else if (j>1) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_4_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _a_0_3 = *((pa+6));
                    _ai_0_3 = *((pa+6)+1);
                    _a_0_4 = *((pa+8));
                    _ai_0_4 = *((pa+8)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_3_0 = __fxcxnpma(_c_3_0, _b_0_0, _ai_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_4_0 = __fxcxnpma(_c_4_0, _b_0_0, _ai_0_4);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_3_2 = __fxcxnpma(_c_3_2, _b_0_2, _ai_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_4_2 = __fxcxnpma(_c_4_2, _b_0_2, _ai_0_4);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj*2+2, _c_3_2);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj*2+2, _c_4_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _a_0_3 = *((pa+6));
                    _ai_0_3 = *((pa+6)+1);
                    _a_0_4 = *((pa+8));
                    _ai_0_4 = *((pa+8)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_3_0 = __fxcxnpma(_c_3_0, _b_0_0, _ai_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_4_0 = __fxcxnpma(_c_4_0, _b_0_0, _ai_0_4);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
            }
        }
        for (; i+3<=dimi ; i+=3) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4*2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                _c_2_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _c_1_6 = __fxcxnpma(_c_1_6, _b_0_6, _ai_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _b_0_6, _a_0_2);
                    _c_2_6 = __fxcxnpma(_c_2_6, _b_0_6, _ai_0_2);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
                __stfpd(xc+(i+2)*dimj*2+6, _c_2_6);
            }
            if (j>3) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                _c_2_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _c_1_6 = __fxcxnpma(_c_1_6, _b_0_6, _ai_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _b_0_6, _a_0_2);
                    _c_2_6 = __fxcxnpma(_c_2_6, _b_0_6, _ai_0_2);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
                __stfpd(xc+(i+2)*dimj*2+6, _c_2_6);
            }
            else if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_2_4 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _c_1_4 = __fxcxnpma(_c_1_4, _b_0_4, _ai_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _b_0_4, _a_0_2);
                    _c_2_4 = __fxcxnpma(_c_2_4, _b_0_4, _ai_0_2);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+2)*dimj*2+4, _c_2_4);
            }
            else if (j>1) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_1_2 = __fxcxnpma(_c_1_2, _b_0_2, _ai_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_2_2 = __fxcxnpma(_c_2_2, _b_0_2, _ai_0_2);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _a_0_1 = *((pa+2));
                    _ai_0_1 = *((pa+2)+1);
                    _a_0_2 = *((pa+4));
                    _ai_0_2 = *((pa+4)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_1_0 = __fxcxnpma(_c_1_0, _b_0_0, _ai_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_2_0 = __fxcxnpma(_c_2_0, _b_0_0, _ai_0_2);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
            }
        }
        for (; i+1<=dimi ; i+=1) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4*2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
            }
            if (j>3) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_0_6 = __fxcxnpma(_c_0_6, _b_0_6, _ai_0_0);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
            }
            else if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_0_4 = __fxcxnpma(_c_0_4, _b_0_4, _ai_0_0);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
            }
            else if (j>1) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_0_2 = __fxcxnpma(_c_0_2, _b_0_2, _ai_0_0);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                    _a_0_0 = *((pa+0));
                    _ai_0_0 = *((pa+0)+1);
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_0_0 = __fxcxnpma(_c_0_0, _b_0_0, _ai_0_0);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
            }
        }
    }
}
#endif
