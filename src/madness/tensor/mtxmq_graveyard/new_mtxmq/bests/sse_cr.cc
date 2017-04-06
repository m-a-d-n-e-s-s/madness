#include <immintrin.h>
#include <complex.h>

void mtxmq(long dimi, long dimj, long dimk, double complex * __restrict__ c_x, const double complex * __restrict__ a_x, const double  * __restrict__ b_x) {
    int i, j, k, ii;
    double * __restrict__ c = (double*)c_x;
    const double * __restrict__ a = (double*)a_x;
    const double * __restrict__ b = (double*)b_x;
    __m128d _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11;
    __m128d _a_0_0, _a_0_1;
     __m128d _az_0_0, _az_0_1;
     __m128d _bz_0_0, _bz_0_1, _bz_0_2, _bz_0_3, _bz_0_4, _bz_0_5, _bz_0_6, _bz_0_7, _bz_0_8, _bz_0_9, _bz_0_10, _bz_0_11;
    for (j=dimj; j>6; j-=6,c+=6*2,b+=6) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            _c_1_8 = _mm_setzero_pd();
            _c_1_10 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_1), _c_1_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_1), _c_1_8);
                _bz_0_10 = _mm_load1_pd((pb+5));
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_0), _c_0_10);
                _c_1_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_1), _c_1_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj*2+10, _c_0_10);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj*2+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj*2+8, _c_1_8);
            _mm_storeu_pd(c+(i+1)*dimj*2+10, _c_1_10);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _bz_0_10 = _mm_load1_pd((pb+5));
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_0), _c_0_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj*2+10, _c_0_10);
        }
    }
    if (j>5) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            _c_1_8 = _mm_setzero_pd();
            _c_1_10 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_1), _c_1_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_1), _c_1_8);
                _bz_0_10 = _mm_load1_pd((pb+5));
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_0), _c_0_10);
                _c_1_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_1), _c_1_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj*2+10, _c_0_10);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj*2+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj*2+8, _c_1_8);
            _mm_storeu_pd(c+(i+1)*dimj*2+10, _c_1_10);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _bz_0_10 = _mm_load1_pd((pb+5));
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_bz_0_10, _az_0_0), _c_0_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj*2+10, _c_0_10);
        }
    }
    else if (j>4) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            _c_1_8 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_1), _c_1_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_1), _c_1_8);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj*2+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj*2+8, _c_1_8);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _bz_0_8 = _mm_load1_pd((pb+4));
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj*2+8, _c_0_8);
        }
    }
    else if (j>3) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_1), _c_1_6);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj*2+6, _c_1_6);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_6 = _mm_load1_pd((pb+3));
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_bz_0_6, _az_0_0), _c_0_6);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj*2+6, _c_0_6);
        }
    }
    else if (j>2) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+4, _c_1_4);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _bz_0_4 = _mm_load1_pd((pb+2));
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj*2+4, _c_0_4);
        }
    }
    else if (j>1) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_1), _c_1_2);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+2, _c_1_2);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_2 = _mm_load1_pd((pb+1));
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_bz_0_2, _az_0_0), _c_0_2);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj*2+2, _c_0_2);
        }
    }
    else {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _az_0_1 = _mm_loadu_pd((pa+2));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
            _mm_storeu_pd(c+(i+1)*dimj*2+0, _c_1_0);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm_loadu_pd((pa+0));
                _bz_0_0 = _mm_load1_pd((pb+0));
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
            }
            _mm_storeu_pd(c+(i+0)*dimj*2+0, _c_0_0);
        }
    }
}
