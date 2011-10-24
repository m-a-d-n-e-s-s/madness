#include <builtins.h>
#include <complex>

#ifdef HAVE_IBMBGP
namespace madness {
    void bgpmTxmq(long dimi, long dimj, long dimk, double  *  c_x,  const double  *  a_x,  const double  *  b_x) {
        int i, j, k, ii;
        double *  c = (double*)c_x;
         double *  a = (double*)a_x;
         double *  b = (double*)b_x;
        __complex__ double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_4_0, _c_4_1, _c_4_2, _c_4_3, _c_5_0, _c_5_1, _c_5_2, _c_5_3, _c_6_0, _c_6_1, _c_6_2, _c_6_3, _c_7_0, _c_7_1, _c_7_2, _c_7_3, _b_0_0, _b_0_1, _b_0_2, _b_0_3;
        double _a_0_0, _a_0_1, _a_0_2, _a_0_3, _a_0_4, _a_0_5, _a_0_6, _a_0_7;
        for (i=0; i+8<=dimi ; i+=8) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4,xb+=4) {
                 double*  pb = xb;
                 double*  pa = a+i;
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                _c_6_0 = __cmplx(0.0,0.0);
                _c_6_2 = __cmplx(0.0,0.0);
                _c_7_0 = __cmplx(0.0,0.0);
                _c_7_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _a_0_6 = *((pa+6));
                    _a_0_7 = *((pa+7));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                    _c_6_0 = __fxcpmadd(_c_6_0, _b_0_0, _a_0_6);
                    _c_7_0 = __fxcpmadd(_c_7_0, _b_0_0, _a_0_7);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _b_0_2, _a_0_5);
                    _c_6_2 = __fxcpmadd(_c_6_2, _b_0_2, _a_0_6);
                    _c_7_2 = __fxcpmadd(_c_7_2, _b_0_2, _a_0_7);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj+2, _c_4_2);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj+2, _c_5_2);
                __stfpd(xc+(i+6)*dimj+0, _c_6_0);
                __stfpd(xc+(i+6)*dimj+2, _c_6_2);
                __stfpd(xc+(i+7)*dimj+0, _c_7_0);
                __stfpd(xc+(i+7)*dimj+2, _c_7_2);
            }
            if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i;
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                _c_6_0 = __cmplx(0.0,0.0);
                _c_6_2 = __cmplx(0.0,0.0);
                _c_7_0 = __cmplx(0.0,0.0);
                _c_7_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _a_0_6 = *((pa+6));
                    _a_0_7 = *((pa+7));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                    _c_6_0 = __fxcpmadd(_c_6_0, _b_0_0, _a_0_6);
                    _c_7_0 = __fxcpmadd(_c_7_0, _b_0_0, _a_0_7);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _b_0_2, _a_0_5);
                    _c_6_2 = __fxcpmadd(_c_6_2, _b_0_2, _a_0_6);
                    _c_7_2 = __fxcpmadd(_c_7_2, _b_0_2, _a_0_7);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj+2, _c_4_2);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj+2, _c_5_2);
                __stfpd(xc+(i+6)*dimj+0, _c_6_0);
                __stfpd(xc+(i+6)*dimj+2, _c_6_2);
                __stfpd(xc+(i+7)*dimj+0, _c_7_0);
                __stfpd(xc+(i+7)*dimj+2, _c_7_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_5_0 = __cmplx(0.0,0.0);
                _c_6_0 = __cmplx(0.0,0.0);
                _c_7_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _a_0_6 = *((pa+6));
                    _a_0_7 = *((pa+7));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                    _c_6_0 = __fxcpmadd(_c_6_0, _b_0_0, _a_0_6);
                    _c_7_0 = __fxcpmadd(_c_7_0, _b_0_0, _a_0_7);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
                __stfpd(xc+(i+6)*dimj+0, _c_6_0);
                __stfpd(xc+(i+7)*dimj+0, _c_7_0);
            }
        }
        for (; i+6<=dimi ; i+=6) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4,xb+=4) {
                 double*  pb = xb;
                 double*  pa = a+i;
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _b_0_2, _a_0_5);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj+2, _c_4_2);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj+2, _c_5_2);
            }
            if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i;
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _b_0_2, _a_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _b_0_2, _a_0_5);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+4)*dimj+2, _c_4_2);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj+2, _c_5_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_5_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _a_0_4 = *((pa+4));
                    _a_0_5 = *((pa+5));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _b_0_0, _a_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _b_0_0, _a_0_5);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+4)*dimj+0, _c_4_0);
                __stfpd(xc+(i+5)*dimj+0, _c_5_0);
            }
        }
        for (; i+4<=dimi ; i+=4) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4,xb+=4) {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
            }
            if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_2_2 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_3_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _b_0_2, _a_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _b_0_2, _a_0_3);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj+2, _c_3_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _a_0_2 = *((pa+2));
                    _a_0_3 = *((pa+3));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _b_0_0, _a_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _b_0_0, _a_0_3);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj+0, _c_3_0);
            }
        }
        for (; i+2<=dimi ; i+=2) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4,xb+=4) {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
            }
            if (j>2) {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj+2, _c_1_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                }
                __stfpd(xc+(i+0)*dimj+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj+0, _c_1_0);
            }
        }
    }
}
#endif
