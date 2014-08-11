#include <immintrin.h>
#include <complex.h>

void mtxmq(long dimi, long dimj, long dimk, double complex * __restrict__ c_x, const double complex * __restrict__ a_x, const double complex * __restrict__ b_x) {
    int i, j, k, ii;
    double * __restrict__ c = (double*)c_x;
    const double * __restrict__ a = (double*)a_x;
    const double * __restrict__ b = (double*)b_x;
    __m256d _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_0_14, _c_0_15, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _c_1_14, _c_1_15, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13, _b_0_14, _b_0_15;
    __m256d _a_0_0, _a_0_1;
     __m256d _br_0_0, _br_0_1, _br_0_2, _br_0_3, _br_0_4, _br_0_5, _br_0_6, _br_0_7, _br_0_8, _br_0_9, _br_0_10, _br_0_11, _br_0_12, _br_0_13, _br_0_14, _br_0_15;
     __m256d _ai_0_0, _ai_0_1;
    
    __m256i mask;
    j = dimj % 2;
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
    for (i=0; i+2<=dimi ; i+=2) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=dimj; j>8; j-=8,xc+=8*2,xb+=8*2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_0_12 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_1_12 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _a_0_1 = _mm256_broadcast_sd((pa+2));
                _ai_0_1 = _mm256_broadcast_sd((pa+2)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_1_0 = _mm256_addsub_pd(_c_1_0, _mm256_mul_pd(_ai_0_1, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_1_4 = _mm256_addsub_pd(_c_1_4, _mm256_mul_pd(_ai_0_1, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _c_1_8 = _mm256_addsub_pd(_c_1_8, _mm256_mul_pd(_ai_0_1, _br_0_8));
                _b_0_12 = _mm256_loadu_pd(pb+12);
                _br_0_12 = _mm256_permute_pd(_b_0_12, 5);
                _c_0_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_0_12 = _mm256_addsub_pd(_c_0_12, _mm256_mul_pd(_ai_0_0, _br_0_12));
                _c_1_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_1), _c_1_12);
                _c_1_12 = _mm256_addsub_pd(_c_1_12, _mm256_mul_pd(_ai_0_1, _br_0_12));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+12, _c_0_12);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+4, _c_1_4);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+8, _c_1_8);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+12, _c_1_12);
        }
        if (j>6) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_0_12 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_1_12 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _a_0_1 = _mm256_broadcast_sd((pa+2));
                _ai_0_1 = _mm256_broadcast_sd((pa+2)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_1_0 = _mm256_addsub_pd(_c_1_0, _mm256_mul_pd(_ai_0_1, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_1_4 = _mm256_addsub_pd(_c_1_4, _mm256_mul_pd(_ai_0_1, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _c_1_8 = _mm256_addsub_pd(_c_1_8, _mm256_mul_pd(_ai_0_1, _br_0_8));
                _b_0_12 = _mm256_loadu_pd(pb+12);
                _br_0_12 = _mm256_permute_pd(_b_0_12, 5);
                _c_0_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_0_12 = _mm256_addsub_pd(_c_0_12, _mm256_mul_pd(_ai_0_0, _br_0_12));
                _c_1_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_1), _c_1_12);
                _c_1_12 = _mm256_addsub_pd(_c_1_12, _mm256_mul_pd(_ai_0_1, _br_0_12));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+12, mask, _c_0_12);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+4, _c_1_4);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+8, _c_1_8);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+12, mask, _c_1_12);
        }
        else if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _a_0_1 = _mm256_broadcast_sd((pa+2));
                _ai_0_1 = _mm256_broadcast_sd((pa+2)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_1_0 = _mm256_addsub_pd(_c_1_0, _mm256_mul_pd(_ai_0_1, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_1_4 = _mm256_addsub_pd(_c_1_4, _mm256_mul_pd(_ai_0_1, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_1), _c_1_8);
                _c_1_8 = _mm256_addsub_pd(_c_1_8, _mm256_mul_pd(_ai_0_1, _br_0_8));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+8, mask, _c_0_8);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+4, _c_1_4);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+8, mask, _c_1_8);
        }
        else if (j>2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _a_0_1 = _mm256_broadcast_sd((pa+2));
                _ai_0_1 = _mm256_broadcast_sd((pa+2)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_1_0 = _mm256_addsub_pd(_c_1_0, _mm256_mul_pd(_ai_0_1, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_1), _c_1_4);
                _c_1_4 = _mm256_addsub_pd(_c_1_4, _mm256_mul_pd(_ai_0_1, _br_0_4));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+4, mask, _c_0_4);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+4, mask, _c_1_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _a_0_1 = _mm256_broadcast_sd((pa+2));
                _ai_0_1 = _mm256_broadcast_sd((pa+2)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_1), _c_1_0);
                _c_1_0 = _mm256_addsub_pd(_c_1_0, _mm256_mul_pd(_ai_0_1, _br_0_0));
            }
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+0, mask, _c_0_0);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+0, mask, _c_1_0);
        }
    }
    for (; i+1<=dimi ; i+=1) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=dimj; j>8; j-=8,xc+=8*2,xb+=8*2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_0_12 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
                _b_0_12 = _mm256_loadu_pd(pb+12);
                _br_0_12 = _mm256_permute_pd(_b_0_12, 5);
                _c_0_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_0_12 = _mm256_addsub_pd(_c_0_12, _mm256_mul_pd(_ai_0_0, _br_0_12));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+12, _c_0_12);
        }
        if (j>6) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_0_12 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
                _b_0_12 = _mm256_loadu_pd(pb+12);
                _br_0_12 = _mm256_permute_pd(_b_0_12, 5);
                _c_0_12 = _mm256_add_pd(_mm256_mul_pd(_b_0_12, _a_0_0), _c_0_12);
                _c_0_12 = _mm256_addsub_pd(_c_0_12, _mm256_mul_pd(_ai_0_0, _br_0_12));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+12, mask, _c_0_12);
        }
        else if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
                _b_0_8 = _mm256_loadu_pd(pb+8);
                _br_0_8 = _mm256_permute_pd(_b_0_8, 5);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_b_0_8, _a_0_0), _c_0_8);
                _c_0_8 = _mm256_addsub_pd(_c_0_8, _mm256_mul_pd(_ai_0_0, _br_0_8));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+8, mask, _c_0_8);
        }
        else if (j>2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
                _b_0_4 = _mm256_loadu_pd(pb+4);
                _br_0_4 = _mm256_permute_pd(_b_0_4, 5);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_b_0_4, _a_0_0), _c_0_4);
                _c_0_4 = _mm256_addsub_pd(_c_0_4, _mm256_mul_pd(_ai_0_0, _br_0_4));
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+4, mask, _c_0_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj*2,pa+=dimi*2) {
                _a_0_0 = _mm256_broadcast_sd((pa+0));
                _ai_0_0 = _mm256_broadcast_sd((pa+0)+1);
                _b_0_0 = _mm256_loadu_pd(pb+0);
                _br_0_0 = _mm256_permute_pd(_b_0_0, 5);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_b_0_0, _a_0_0), _c_0_0);
                _c_0_0 = _mm256_addsub_pd(_c_0_0, _mm256_mul_pd(_ai_0_0, _br_0_0));
            }
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+0, mask, _c_0_0);
        }
    }
}
