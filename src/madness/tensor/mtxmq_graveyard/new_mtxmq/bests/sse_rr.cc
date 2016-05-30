#include <immintrin.h>
#include <complex.h>

void mtxmq(long dimi, long dimj, long dimk, double  * __restrict__ c_x, const double  * __restrict__ a_x, const double  * __restrict__ b_x) {
    int i, j, k, ii;
    double * __restrict__ c = (double*)c_x;
    const double * __restrict__ a = (double*)a_x;
    const double * __restrict__ b = (double*)b_x;
    __m128d _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13;
    __m128d _a_0_0, _a_0_1;
    for (j=dimj; j>14; j-=14,c+=14,b+=14) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_0_12 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            _c_1_8 = _mm_setzero_pd();
            _c_1_10 = _mm_setzero_pd();
            _c_1_12 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_1), _c_1_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
                _c_1_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_1), _c_1_10);
                _b_0_12 = _mm_loadu_pd(pb+12);
                _c_0_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_1_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_1), _c_1_12);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
            _mm_storeu_pd(c+(i+0)*dimj+12, _c_0_12);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj+8, _c_1_8);
            _mm_storeu_pd(c+(i+1)*dimj+10, _c_1_10);
            _mm_storeu_pd(c+(i+1)*dimj+12, _c_1_12);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_0_12 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
                _b_0_12 = _mm_loadu_pd(pb+12);
                _c_0_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_0), _c_0_12);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
            _mm_storeu_pd(c+(i+0)*dimj+12, _c_0_12);
        }
    }
    if (j>12) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_0_12 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            _c_1_8 = _mm_setzero_pd();
            _c_1_10 = _mm_setzero_pd();
            _c_1_12 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_1), _c_1_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
                _c_1_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_1), _c_1_10);
                _b_0_12 = _mm_loadu_pd(pb+12);
                _c_0_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_1_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_1), _c_1_12);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
            _mm_storeu_pd(c+(i+0)*dimj+12, _c_0_12);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj+8, _c_1_8);
            _mm_storeu_pd(c+(i+1)*dimj+10, _c_1_10);
            _mm_storeu_pd(c+(i+1)*dimj+12, _c_1_12);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            _c_0_12 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
                _b_0_12 = _mm_loadu_pd(pb+12);
                _c_0_12 = _mm_add_pd(_mm_mul_pd(_b_0_12, _a_0_0), _c_0_12);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
            _mm_storeu_pd(c+(i+0)*dimj+12, _c_0_12);
        }
    }
    else if (j>10) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
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
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_1), _c_1_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
                _c_1_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_1), _c_1_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj+8, _c_1_8);
            _mm_storeu_pd(c+(i+1)*dimj+10, _c_1_10);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            _c_0_10 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _b_0_10 = _mm_loadu_pd(pb+10);
                _c_0_10 = _mm_add_pd(_mm_mul_pd(_b_0_10, _a_0_0), _c_0_10);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+0)*dimj+10, _c_0_10);
        }
    }
    else if (j>8) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
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
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_1), _c_1_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_1), _c_1_8);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj+6, _c_1_6);
            _mm_storeu_pd(c+(i+1)*dimj+8, _c_1_8);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_0_8 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _b_0_8 = _mm_loadu_pd(pb+8);
                _c_0_8 = _mm_add_pd(_mm_mul_pd(_b_0_8, _a_0_0), _c_0_8);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+0)*dimj+8, _c_0_8);
        }
    }
    else if (j>6) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            _c_1_6 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
                _c_1_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_1), _c_1_6);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
            _mm_storeu_pd(c+(i+1)*dimj+6, _c_1_6);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_0_6 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_6 = _mm_loadu_pd(pb+6);
                _c_0_6 = _mm_add_pd(_mm_mul_pd(_b_0_6, _a_0_0), _c_0_6);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+0)*dimj+6, _c_0_6);
        }
    }
    else if (j>4) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            _c_1_4 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_1), _c_1_4);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
            _mm_storeu_pd(c+(i+1)*dimj+4, _c_1_4);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_0_4 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _b_0_4 = _mm_loadu_pd(pb+4);
                _c_0_4 = _mm_add_pd(_mm_mul_pd(_b_0_4, _a_0_0), _c_0_4);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+0)*dimj+4, _c_0_4);
        }
    }
    else if (j>2) {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            _c_1_2 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
                _c_1_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_1), _c_1_2);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
            _mm_storeu_pd(c+(i+1)*dimj+2, _c_1_2);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_0_2 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_2 = _mm_loadu_pd(pb+2);
                _c_0_2 = _mm_add_pd(_mm_mul_pd(_b_0_2, _a_0_0), _c_0_2);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+0)*dimj+2, _c_0_2);
        }
    }
    else {
        for (i=0; i+2<=dimi ; i+=2) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            _c_1_0 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _a_0_1 = _mm_load1_pd((pa+1));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_1), _c_1_0);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
            _mm_storeu_pd(c+(i+1)*dimj+0, _c_1_0);
        }
        for (; i+1<=dimi ; i+=1) {
            const double* __restrict__ pb = b;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi) {
                _a_0_0 = _mm_load1_pd((pa+0));
                _b_0_0 = _mm_loadu_pd(pb+0);
                _c_0_0 = _mm_add_pd(_mm_mul_pd(_b_0_0, _a_0_0), _c_0_0);
            }
            _mm_storeu_pd(c+(i+0)*dimj+0, _c_0_0);
        }
    }
}
