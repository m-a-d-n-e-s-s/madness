#include <builtins.h>
#include <complex>

typedef std::complex<double> double_complex;

#ifdef HAVE_IBMBGP
namespace madness {
    void bgpmTxmq(long dimi, long dimj, long dimk, double_complex *  c_x,  const double_complex *  a_x,  const double  *  b_x) {
        int i, j, k, ii;
        double *  c = (double*)c_x;
         double *  a = (double*)a_x;
         double *  b = (double*)b_x;
        __complex__ double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_4_0, _c_4_1, _c_4_2, _c_4_3, _c_4_4, _c_4_5, _c_4_6, _c_4_7, _c_5_0, _c_5_1, _c_5_2, _c_5_3, _c_5_4, _c_5_5, _c_5_6, _c_5_7, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7;
        double _a_0_0, _a_0_1, _a_0_2, _a_0_3, _a_0_4, _a_0_5;
         __complex__ double _az_0_0, _az_0_1, _az_0_2, _az_0_3, _az_0_4, _az_0_5;
         double _bz_0_0, _bz_0_1, _bz_0_2, _bz_0_3, _bz_0_4, _bz_0_5, _bz_0_6, _bz_0_7;
        for (i=0; i+6<=dimi ; i+=6) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4) {
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                _c_5_4 = __cmplx(0.0,0.0);
                _c_5_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _az_0_4 = __lfpd((pa+8));
                    _az_0_5 = __lfpd((pa+10));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _bz_0_0, _az_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _bz_0_0, _az_0_5);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _bz_0_2, _az_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _bz_0_2, _az_0_5);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _bz_0_4, _az_0_4);
                    _c_5_4 = __fxcpmadd(_c_5_4, _bz_0_4, _az_0_5);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _bz_0_6, _az_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _bz_0_6, _az_0_3);
                    _c_4_6 = __fxcpmadd(_c_4_6, _bz_0_6, _az_0_4);
                    _c_5_6 = __fxcpmadd(_c_5_6, _bz_0_6, _az_0_5);
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
                __stfpd(xc+(i+5)*dimj*2+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj*2+2, _c_5_2);
                __stfpd(xc+(i+5)*dimj*2+4, _c_5_4);
                __stfpd(xc+(i+5)*dimj*2+6, _c_5_6);
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                _c_5_4 = __cmplx(0.0,0.0);
                _c_5_6 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _az_0_4 = __lfpd((pa+8));
                    _az_0_5 = __lfpd((pa+10));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _bz_0_0, _az_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _bz_0_0, _az_0_5);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _bz_0_2, _az_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _bz_0_2, _az_0_5);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _bz_0_4, _az_0_4);
                    _c_5_4 = __fxcpmadd(_c_5_4, _bz_0_4, _az_0_5);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _bz_0_6, _az_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _bz_0_6, _az_0_3);
                    _c_4_6 = __fxcpmadd(_c_4_6, _bz_0_6, _az_0_4);
                    _c_5_6 = __fxcpmadd(_c_5_6, _bz_0_6, _az_0_5);
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
                __stfpd(xc+(i+5)*dimj*2+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj*2+2, _c_5_2);
                __stfpd(xc+(i+5)*dimj*2+4, _c_5_4);
                __stfpd(xc+(i+5)*dimj*2+6, _c_5_6);
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                _c_5_4 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _az_0_4 = __lfpd((pa+8));
                    _az_0_5 = __lfpd((pa+10));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _bz_0_0, _az_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _bz_0_0, _az_0_5);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _bz_0_2, _az_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _bz_0_2, _az_0_5);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
                    _c_4_4 = __fxcpmadd(_c_4_4, _bz_0_4, _az_0_4);
                    _c_5_4 = __fxcpmadd(_c_5_4, _bz_0_4, _az_0_5);
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
                __stfpd(xc+(i+5)*dimj*2+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj*2+2, _c_5_2);
                __stfpd(xc+(i+5)*dimj*2+4, _c_5_4);
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
                _c_5_0 = __cmplx(0.0,0.0);
                _c_5_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _az_0_4 = __lfpd((pa+8));
                    _az_0_5 = __lfpd((pa+10));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _bz_0_0, _az_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _bz_0_0, _az_0_5);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _c_4_2 = __fxcpmadd(_c_4_2, _bz_0_2, _az_0_4);
                    _c_5_2 = __fxcpmadd(_c_5_2, _bz_0_2, _az_0_5);
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
                __stfpd(xc+(i+5)*dimj*2+0, _c_5_0);
                __stfpd(xc+(i+5)*dimj*2+2, _c_5_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                _c_4_0 = __cmplx(0.0,0.0);
                _c_5_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _az_0_4 = __lfpd((pa+8));
                    _az_0_5 = __lfpd((pa+10));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _c_4_0 = __fxcpmadd(_c_4_0, _bz_0_0, _az_0_4);
                    _c_5_0 = __fxcpmadd(_c_5_0, _bz_0_0, _az_0_5);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+4)*dimj*2+0, _c_4_0);
                __stfpd(xc+(i+5)*dimj*2+0, _c_5_0);
            }
        }
        for (; i+4<=dimi ; i+=4) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4) {
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _bz_0_6, _az_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _bz_0_6, _az_0_3);
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                    _c_2_6 = __fxcpmadd(_c_2_6, _bz_0_6, _az_0_2);
                    _c_3_6 = __fxcpmadd(_c_3_6, _bz_0_6, _az_0_3);
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _c_2_4 = __fxcpmadd(_c_2_4, _bz_0_4, _az_0_2);
                    _c_3_4 = __fxcpmadd(_c_3_4, _bz_0_4, _az_0_3);
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _c_2_2 = __fxcpmadd(_c_2_2, _bz_0_2, _az_0_2);
                    _c_3_2 = __fxcpmadd(_c_3_2, _bz_0_2, _az_0_3);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+2)*dimj*2+2, _c_2_2);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
                __stfpd(xc+(i+3)*dimj*2+2, _c_3_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_2_0 = __cmplx(0.0,0.0);
                _c_3_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _az_0_2 = __lfpd((pa+4));
                    _az_0_3 = __lfpd((pa+6));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _c_2_0 = __fxcpmadd(_c_2_0, _bz_0_0, _az_0_2);
                    _c_3_0 = __fxcpmadd(_c_3_0, _bz_0_0, _az_0_3);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+2)*dimj*2+0, _c_2_0);
                __stfpd(xc+(i+3)*dimj*2+0, _c_3_0);
            }
        }
        for (; i+2<=dimi ; i+=2) {
             double*  xb = b;
            double*  xc = c;
            for (j=dimj; j>4; j-=4,xc+=4*2,xb+=4) {
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                    _bz_0_6 = *((pb+3));
                    _c_0_6 = __fxcpmadd(_c_0_6, _bz_0_6, _az_0_0);
                    _c_1_6 = __fxcpmadd(_c_1_6, _bz_0_6, _az_0_1);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+0)*dimj*2+6, _c_0_6);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
                __stfpd(xc+(i+1)*dimj*2+6, _c_1_6);
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
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                    _bz_0_4 = *((pb+2));
                    _c_0_4 = __fxcpmadd(_c_0_4, _bz_0_4, _az_0_0);
                    _c_1_4 = __fxcpmadd(_c_1_4, _bz_0_4, _az_0_1);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+0)*dimj*2+4, _c_0_4);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
                __stfpd(xc+(i+1)*dimj*2+4, _c_1_4);
            }
            else if (j>1) {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_0_2 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                _c_1_2 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                    _bz_0_2 = *((pb+1));
                    _c_0_2 = __fxcpmadd(_c_0_2, _bz_0_2, _az_0_0);
                    _c_1_2 = __fxcpmadd(_c_1_2, _bz_0_2, _az_0_1);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+0)*dimj*2+2, _c_0_2);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
                __stfpd(xc+(i+1)*dimj*2+2, _c_1_2);
            }
            else {
                 double*  pb = xb;
                 double*  pa = a+i*2;
                _c_0_0 = __cmplx(0.0,0.0);
                _c_1_0 = __cmplx(0.0,0.0);
                for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                    _az_0_0 = __lfpd((pa+0));
                    _az_0_1 = __lfpd((pa+2));
                    _bz_0_0 = *((pb+0));
                    _c_0_0 = __fxcpmadd(_c_0_0, _bz_0_0, _az_0_0);
                    _c_1_0 = __fxcpmadd(_c_1_0, _bz_0_0, _az_0_1);
                }
                __stfpd(xc+(i+0)*dimj*2+0, _c_0_0);
                __stfpd(xc+(i+1)*dimj*2+0, _c_1_0);
            }
        }
    }
}
#endif
