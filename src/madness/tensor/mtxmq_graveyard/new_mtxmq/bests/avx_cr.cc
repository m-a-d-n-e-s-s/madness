#include <immintrin.h>
#include <complex.h>

void mtxmq(long dimi, long dimj, long dimk, double complex * __restrict__ c_x, const double complex * __restrict__ a_x, const double  * __restrict__ b_x) {
    int i, j, k, ii;
    double * __restrict__ c = (double*)c_x;
    const double * __restrict__ a = (double*)a_x;
    const double * __restrict__ b = (double*)b_x;
    __m256d _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_2_8, _c_2_9, _c_2_10, _c_2_11, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_3_8, _c_3_9, _c_3_10, _c_3_11, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11;
    __m256d _a_0_0, _a_0_1, _a_0_2, _a_0_3;
     __m256d _az_0_0, _az_0_1, _az_0_2, _az_0_3;
     __m256d _bz_0_0, _bz_0_1, _bz_0_2, _bz_0_3, _bz_0_4, _bz_0_5, _bz_0_6, _bz_0_7, _bz_0_8, _bz_0_9, _bz_0_10, _bz_0_11;
    
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
    for (i=0; i+4<=dimi ; i+=4) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=dimj; j>6; j-=6,xc+=6*2,xb+=6) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            _c_2_8 = _mm256_setzero_pd();
            _c_3_0 = _mm256_setzero_pd();
            _c_3_4 = _mm256_setzero_pd();
            _c_3_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _az_0_1 = _mm256_broadcast_pd((const __m128d*)(pa+2));
                _az_0_2 = _mm256_broadcast_pd((const __m128d*)(pa+4));
                _az_0_3 = _mm256_broadcast_pd((const __m128d*)(pa+6));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_2), _c_2_0);
                _c_3_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_3), _c_3_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_2), _c_2_4);
                _c_3_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_3), _c_3_4);
                _bz_0_8 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+4)),12);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_1), _c_1_8);
                _c_2_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_2), _c_2_8);
                _c_3_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_3), _c_3_8);
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+4, _c_1_4);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+8, _c_1_8);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+0, _c_2_0);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+4, _c_2_4);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+8, _c_2_8);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+0, _c_3_0);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+4, _c_3_4);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+8, _c_3_8);
        }
        if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_1_8 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            _c_2_8 = _mm256_setzero_pd();
            _c_3_0 = _mm256_setzero_pd();
            _c_3_4 = _mm256_setzero_pd();
            _c_3_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _az_0_1 = _mm256_broadcast_pd((const __m128d*)(pa+2));
                _az_0_2 = _mm256_broadcast_pd((const __m128d*)(pa+4));
                _az_0_3 = _mm256_broadcast_pd((const __m128d*)(pa+6));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_2), _c_2_0);
                _c_3_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_3), _c_3_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_2), _c_2_4);
                _c_3_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_3), _c_3_4);
                _bz_0_8 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+4)),12);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
                _c_1_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_1), _c_1_8);
                _c_2_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_2), _c_2_8);
                _c_3_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_3), _c_3_8);
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+8, mask, _c_0_8);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+4, _c_1_4);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+8, mask, _c_1_8);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+0, _c_2_0);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+4, _c_2_4);
            _mm256_maskstore_pd(xc+(i+2)*dimj*2+8, mask, _c_2_8);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+0, _c_3_0);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+4, _c_3_4);
            _mm256_maskstore_pd(xc+(i+3)*dimj*2+8, mask, _c_3_8);
        }
        else if (j>2) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_1_4 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_2_4 = _mm256_setzero_pd();
            _c_3_0 = _mm256_setzero_pd();
            _c_3_4 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _az_0_1 = _mm256_broadcast_pd((const __m128d*)(pa+2));
                _az_0_2 = _mm256_broadcast_pd((const __m128d*)(pa+4));
                _az_0_3 = _mm256_broadcast_pd((const __m128d*)(pa+6));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_2), _c_2_0);
                _c_3_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_3), _c_3_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _c_1_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_1), _c_1_4);
                _c_2_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_2), _c_2_4);
                _c_3_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_3), _c_3_4);
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+4, mask, _c_0_4);
            _mm256_storeu_pd(xc+(i+1)*dimj*2+0, _c_1_0);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+4, mask, _c_1_4);
            _mm256_storeu_pd(xc+(i+2)*dimj*2+0, _c_2_0);
            _mm256_maskstore_pd(xc+(i+2)*dimj*2+4, mask, _c_2_4);
            _mm256_storeu_pd(xc+(i+3)*dimj*2+0, _c_3_0);
            _mm256_maskstore_pd(xc+(i+3)*dimj*2+4, mask, _c_3_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_1_0 = _mm256_setzero_pd();
            _c_2_0 = _mm256_setzero_pd();
            _c_3_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _az_0_1 = _mm256_broadcast_pd((const __m128d*)(pa+2));
                _az_0_2 = _mm256_broadcast_pd((const __m128d*)(pa+4));
                _az_0_3 = _mm256_broadcast_pd((const __m128d*)(pa+6));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _c_1_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_1), _c_1_0);
                _c_2_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_2), _c_2_0);
                _c_3_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_3), _c_3_0);
            }
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+0, mask, _c_0_0);
            _mm256_maskstore_pd(xc+(i+1)*dimj*2+0, mask, _c_1_0);
            _mm256_maskstore_pd(xc+(i+2)*dimj*2+0, mask, _c_2_0);
            _mm256_maskstore_pd(xc+(i+3)*dimj*2+0, mask, _c_3_0);
        }
    }
    for (; i+1<=dimi ; i+=1) {
        const double* __restrict__ xb = b;
        double* __restrict__ xc = c;
        for (j=dimj; j>6; j-=6,xc+=6*2,xb+=6) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_8 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+4)),12);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+4, _c_0_4);
            _mm256_storeu_pd(xc+(i+0)*dimj*2+8, _c_0_8);
        }
        if (j>4) {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            _c_0_4 = _mm256_setzero_pd();
            _c_0_8 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
                _bz_0_8 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+4)),12);
                _c_0_8 = _mm256_add_pd(_mm256_mul_pd(_bz_0_8, _az_0_0), _c_0_8);
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
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
                _bz_0_4 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+2)),12);
                _c_0_4 = _mm256_add_pd(_mm256_mul_pd(_bz_0_4, _az_0_0), _c_0_4);
            }
            _mm256_storeu_pd(xc+(i+0)*dimj*2+0, _c_0_0);
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+4, mask, _c_0_4);
        }
        else {
            const double* __restrict__ pb = xb;
            const double* __restrict__ pa = a+i*2;
            _c_0_0 = _mm256_setzero_pd();
            for (k=0; k<dimk; k+=1,pb+=dimj,pa+=dimi*2) {
                _az_0_0 = _mm256_broadcast_pd((const __m128d*)(pa+0));
                _bz_0_0 = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*)(pb+0)),12);
                _c_0_0 = _mm256_add_pd(_mm256_mul_pd(_bz_0_0, _az_0_0), _c_0_0);
            }
            _mm256_maskstore_pd(xc+(i+0)*dimj*2+0, mask, _c_0_0);
        }
    }
}
