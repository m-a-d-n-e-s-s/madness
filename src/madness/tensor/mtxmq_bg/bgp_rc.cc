#include <builtins.h>
#include <complex>

typedef std::complex<double> double_complex;

#ifdef HAVE_IBMBGP
namespace madness {
    void bgpmTxmq(long dimi, long dimj, long dimk, double_complex *  c_x,  const double  *  a_x, const  double_complex *  b_x) {
        int i, j, k, ii;
        double *  c = (double*)c_x;
         double *  a = (double*)a_x;
         double *  b = (double*)b_x;
        __complex__ double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_0_14, _c_0_15, _c_0_16, _c_0_17, _c_0_18, _c_0_19, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _c_1_14, _c_1_15, _c_1_16, _c_1_17, _c_1_18, _c_1_19, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13, _b_0_14, _b_0_15, _b_0_16, _b_0_17, _b_0_18, _b_0_19;
        double _a_0_0, _a_0_1;
        for (j=dimj; j>10; j-=10,c+=10*2,b+=10*2) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_0_12 = __cmplx(0.0,0.0);
                _c_0_14 = __cmplx(0.0,0.0);
                _c_0_16 = __cmplx(0.0,0.0);
                _c_0_18 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                _c_1_12 = __cmplx(0.0,0.0);
                _c_1_14 = __cmplx(0.0,0.0);
                _c_1_16 = __cmplx(0.0,0.0);
                _c_1_18 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                    _b_0_12 = __lfpd(pb+12);
                    _c_0_12 = __fxcpmadd(_c_0_12, _b_0_12, _a_0_0);
                    _c_1_12 = __fxcpmadd(_c_1_12, _b_0_12, _a_0_1);
                    _b_0_14 = __lfpd(pb+14);
                    _c_0_14 = __fxcpmadd(_c_0_14, _b_0_14, _a_0_0);
                    _c_1_14 = __fxcpmadd(_c_1_14, _b_0_14, _a_0_1);
                    _b_0_16 = __lfpd(pb+16);
                    _c_0_16 = __fxcpmadd(_c_0_16, _b_0_16, _a_0_0);
                    _c_1_16 = __fxcpmadd(_c_1_16, _b_0_16, _a_0_1);
                    _b_0_18 = __lfpd(pb+18);
                    _c_0_18 = __fxcpmadd(_c_0_18, _b_0_18, _a_0_0);
                    _c_1_18 = __fxcpmadd(_c_1_18, _b_0_18, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+0)*dimj*2+12, _c_0_12);
                __stfpd(c+(i+0)*dimj*2+14, _c_0_14);
                __stfpd(c+(i+0)*dimj*2+16, _c_0_16);
                __stfpd(c+(i+0)*dimj*2+18, _c_0_18);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
                __stfpd(c+(i+1)*dimj*2+12, _c_1_12);
                __stfpd(c+(i+1)*dimj*2+14, _c_1_14);
                __stfpd(c+(i+1)*dimj*2+16, _c_1_16);
                __stfpd(c+(i+1)*dimj*2+18, _c_1_18);
            }
        }
        if (j>9) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_0_12 = __cmplx(0.0,0.0);
                _c_0_14 = __cmplx(0.0,0.0);
                _c_0_16 = __cmplx(0.0,0.0);
                _c_0_18 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                _c_1_12 = __cmplx(0.0,0.0);
                _c_1_14 = __cmplx(0.0,0.0);
                _c_1_16 = __cmplx(0.0,0.0);
                _c_1_18 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                    _b_0_12 = __lfpd(pb+12);
                    _c_0_12 = __fxcpmadd(_c_0_12, _b_0_12, _a_0_0);
                    _c_1_12 = __fxcpmadd(_c_1_12, _b_0_12, _a_0_1);
                    _b_0_14 = __lfpd(pb+14);
                    _c_0_14 = __fxcpmadd(_c_0_14, _b_0_14, _a_0_0);
                    _c_1_14 = __fxcpmadd(_c_1_14, _b_0_14, _a_0_1);
                    _b_0_16 = __lfpd(pb+16);
                    _c_0_16 = __fxcpmadd(_c_0_16, _b_0_16, _a_0_0);
                    _c_1_16 = __fxcpmadd(_c_1_16, _b_0_16, _a_0_1);
                    _b_0_18 = __lfpd(pb+18);
                    _c_0_18 = __fxcpmadd(_c_0_18, _b_0_18, _a_0_0);
                    _c_1_18 = __fxcpmadd(_c_1_18, _b_0_18, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+0)*dimj*2+12, _c_0_12);
                __stfpd(c+(i+0)*dimj*2+14, _c_0_14);
                __stfpd(c+(i+0)*dimj*2+16, _c_0_16);
                __stfpd(c+(i+0)*dimj*2+18, _c_0_18);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
                __stfpd(c+(i+1)*dimj*2+12, _c_1_12);
                __stfpd(c+(i+1)*dimj*2+14, _c_1_14);
                __stfpd(c+(i+1)*dimj*2+16, _c_1_16);
                __stfpd(c+(i+1)*dimj*2+18, _c_1_18);
            }
        }
        else if (j>8) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_0_12 = __cmplx(0.0,0.0);
                _c_0_14 = __cmplx(0.0,0.0);
                _c_0_16 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                _c_1_12 = __cmplx(0.0,0.0);
                _c_1_14 = __cmplx(0.0,0.0);
                _c_1_16 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                    _b_0_12 = __lfpd(pb+12);
                    _c_0_12 = __fxcpmadd(_c_0_12, _b_0_12, _a_0_0);
                    _c_1_12 = __fxcpmadd(_c_1_12, _b_0_12, _a_0_1);
                    _b_0_14 = __lfpd(pb+14);
                    _c_0_14 = __fxcpmadd(_c_0_14, _b_0_14, _a_0_0);
                    _c_1_14 = __fxcpmadd(_c_1_14, _b_0_14, _a_0_1);
                    _b_0_16 = __lfpd(pb+16);
                    _c_0_16 = __fxcpmadd(_c_0_16, _b_0_16, _a_0_0);
                    _c_1_16 = __fxcpmadd(_c_1_16, _b_0_16, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+0)*dimj*2+12, _c_0_12);
                __stfpd(c+(i+0)*dimj*2+14, _c_0_14);
                __stfpd(c+(i+0)*dimj*2+16, _c_0_16);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
                __stfpd(c+(i+1)*dimj*2+12, _c_1_12);
                __stfpd(c+(i+1)*dimj*2+14, _c_1_14);
                __stfpd(c+(i+1)*dimj*2+16, _c_1_16);
            }
        }
        else if (j>7) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_0_12 = __cmplx(0.0,0.0);
                _c_0_14 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                _c_1_12 = __cmplx(0.0,0.0);
                _c_1_14 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                    _b_0_12 = __lfpd(pb+12);
                    _c_0_12 = __fxcpmadd(_c_0_12, _b_0_12, _a_0_0);
                    _c_1_12 = __fxcpmadd(_c_1_12, _b_0_12, _a_0_1);
                    _b_0_14 = __lfpd(pb+14);
                    _c_0_14 = __fxcpmadd(_c_0_14, _b_0_14, _a_0_0);
                    _c_1_14 = __fxcpmadd(_c_1_14, _b_0_14, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+0)*dimj*2+12, _c_0_12);
                __stfpd(c+(i+0)*dimj*2+14, _c_0_14);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
                __stfpd(c+(i+1)*dimj*2+12, _c_1_12);
                __stfpd(c+(i+1)*dimj*2+14, _c_1_14);
            }
        }
        else if (j>6) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_0_12 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                _c_1_12 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                    _b_0_12 = __lfpd(pb+12);
                    _c_0_12 = __fxcpmadd(_c_0_12, _b_0_12, _a_0_0);
                    _c_1_12 = __fxcpmadd(_c_1_12, _b_0_12, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+0)*dimj*2+12, _c_0_12);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
                __stfpd(c+(i+1)*dimj*2+12, _c_1_12);
            }
        }
        else if (j>5) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_0_10 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                _c_1_10 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                    _b_0_10 = __lfpd(pb+10);
                    _c_0_10 = __fxcpmadd(_c_0_10, _b_0_10, _a_0_0);
                    _c_1_10 = __fxcpmadd(_c_1_10, _b_0_10, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+0)*dimj*2+10, _c_0_10);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
                __stfpd(c+(i+1)*dimj*2+10, _c_1_10);
            }
        }
        else if (j>4) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_0_8 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                _c_1_8 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                    _b_0_8 = __lfpd(pb+8);
                    _c_0_8 = __fxcpmadd(_c_0_8, _b_0_8, _a_0_0);
                    _c_1_8 = __fxcpmadd(_c_1_8, _b_0_8, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+0)*dimj*2+8, _c_0_8);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
                __stfpd(c+(i+1)*dimj*2+8, _c_1_8);
            }
        }
        else if (j>3) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_0_6 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                _c_1_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                    _b_0_6 = __lfpd(pb+6);
                    _c_0_6 = __fxcpmadd(_c_0_6, _b_0_6, _a_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _b_0_6, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(c+(i+1)*dimj*2+6, _c_1_6);
            }
        }
        else if (j>2) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_0_4 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                _c_1_4 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                    _b_0_4 = __lfpd(pb+4);
                    _c_0_4 = __fxcpmadd(_c_0_4, _b_0_4, _a_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _b_0_4, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(c+(i+1)*dimj*2+4, _c_1_4);
            }
        }
        else if (j>1) {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                    _b_0_2 = __lfpd(pb+2);
                    _c_0_2 = __fxcpmadd(_c_0_2, _b_0_2, _a_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _b_0_2, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(c+(i+1)*dimj*2+2, _c_1_2);
            }
        }
        else {
            for (i=0; i+2<=dimi ; i+=2) {
                 double*  pb = b;
                 double*  pa = a+i;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi) {
                    _a_0_0 = *((pa+0));
                    _a_0_1 = *((pa+1));
                    _b_0_0 = __lfpd(pb+0);
                    _c_0_0 = __fxcpmadd(_c_0_0, _b_0_0, _a_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _b_0_0, _a_0_1);
                }
                __stfpd(c+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(c+(i+1)*dimj*2+0, _c_1_0);
            }
        }
    }
}
#endif
