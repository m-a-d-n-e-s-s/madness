#include <immintrin.h>
#include <complex.h>

void mtxmq(long dimi, long dimj, long dimk, double complex * __restrict__ c_x, const double  * __restrict__ a_x, const double complex * __restrict__ b_x) {
    int i, j, k;
    double * __restrict__ c = (double*)c_x;
    const double * __restrict__ a = (double*)a_x;
    const double * __restrict__ b = (double*)b_x;
    long effj = dimj;
    __m256d _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_2_8, _c_2_9, _c_2_10, _c_2_11, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11;
    __m256d _a_0_0, _a_0_1, _a_0_2;
    
    __m256i mask;
    j = effj % 2;
    switch (j) {
        case 0:
            mask = _mm256_set_epi32(-1,-1,-1,-1,-1,-1,-1,-1);
            break;
        case 1:
            mask = _mm256_set_epi32( 0, 0, 0, 0,-1,-1,-1,-1);
            break;
        default:
            return;
    }
    for (i=0; i+3<=dimi; i+=3) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=effj; j>6; j-=6,xc+=6*2,xb+=6*2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            _c_2_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _a_0_1 = _mm256_broadcast_sd((pa+1));
                _a_0_2 = _mm256_broadcast_sd((pa+2));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_2), _c_2_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_2), _c_2_4);
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _c_2_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_2), _c_2_8);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*effj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*effj*2+8, _c_0_8);
            _mm256_storeu_pd(xc+(i+1)*effj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*effj*2+4, _c_1_4);
            _mm256_storeu_pd(xc+(i+1)*effj*2+8, _c_1_8);
            _mm256_storeu_pd(xc+(i+2)*effj*2+0, _c_2_0);
            _mm256_storeu_pd(xc+(i+2)*effj*2+4, _c_2_4);
            _mm256_storeu_pd(xc+(i+2)*effj*2+8, _c_2_8);
        }
        if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            _c_2_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _a_0_1 = _mm256_broadcast_sd((pa+1));
                _a_0_2 = _mm256_broadcast_sd((pa+2));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_2), _c_2_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_2), _c_2_4);
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _c_2_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_2), _c_2_8);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*effj*2+4, _c_0_4);
            _mm256_maskstore_pd(xc+(i+0)*effj*2+8, mask, _c_0_8);
            _mm256_storeu_pd(xc+(i+1)*effj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*effj*2+4, _c_1_4);
            _mm256_maskstore_pd(xc+(i+1)*effj*2+8, mask, _c_1_8);
            _mm256_storeu_pd(xc+(i+2)*effj*2+0, _c_2_0);
            _mm256_storeu_pd(xc+(i+2)*effj*2+4, _c_2_4);
            _mm256_maskstore_pd(xc+(i+2)*effj*2+8, mask, _c_2_8);
        }
        else if (j>2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _a_0_1 = _mm256_broadcast_sd((pa+1));
                _a_0_2 = _mm256_broadcast_sd((pa+2));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_2), _c_2_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_2), _c_2_4);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*effj*2+4, mask, _c_0_4);
            _mm256_storeu_pd(xc+(i+1)*effj*2+0, _c_1_0);
            _mm256_maskstore_pd(xc+(i+1)*effj*2+4, mask, _c_1_4);
            _mm256_storeu_pd(xc+(i+2)*effj*2+0, _c_2_0);
            _mm256_maskstore_pd(xc+(i+2)*effj*2+4, mask, _c_2_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _a_0_1 = _mm256_broadcast_sd((pa+1));
                _a_0_2 = _mm256_broadcast_sd((pa+2));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_2), _c_2_0);
            }
            _mm256_maskstore_pd(xc+(i+0)*effj*2+0, mask, _c_0_0);
            _mm256_maskstore_pd(xc+(i+1)*effj*2+0, mask, _c_1_0);
            _mm256_maskstore_pd(xc+(i+2)*effj*2+0, mask, _c_2_0);
        }
    }
    for (; i+1<=dimi; i+=1) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=effj; j>6; j-=6,xc+=6*2,xb+=6*2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*effj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*effj*2+8, _c_0_8);
        }
        if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*effj*2+4, _c_0_4);
            _mm256_maskstore_pd(xc+(i+0)*effj*2+8, mask, _c_0_8);
        }
        else if (j>2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
            }
            _mm256_storeu_pd(xc+(i+0)*effj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*effj*2+4, mask, _c_0_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i;
            _c_0_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=effj*2,pa+=dimi) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
            }
            _mm256_maskstore_pd(xc+(i+0)*effj*2+0, mask, _c_0_0);
        }
    }
}
