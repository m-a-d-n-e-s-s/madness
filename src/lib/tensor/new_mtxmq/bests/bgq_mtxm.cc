#include <hwi/include/bqc/A2_inlines.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#if !defined(__xlc__) && !defined(__xlC__)
#include <qpxintrin.h>
#endif

namespace madness {
    void bgq_mtxmq_padded(long dimi, long dimj, long dimk, long extb, __complex__ double *  c_x,  const __complex__ double *  a_x,  const __complex__ double *  b_x) {
        int i, j, k;
        double *  c = (double*)c_x;
        double *  a = (double*)a_x;
        double *  b = (double*)b_x;
        long effj = dimj;
        double * c_buf;
        double * b_buf;
        bool free_b = false;
        /* Setup a buffer for c if needed */
        double* c_out = c;
        if (dimj%2) {
            effj = (dimj | 1) + 1;
            posix_memalign((void **) &c, 32, dimi*effj*sizeof(double)*2);
            c_buf = c;
        }
        /* Copy b into a buffer if needed */
        if (extb%2) {
            long t_extb = (dimj | 1) + 1;
            free_b = true;
            posix_memalign((void **) &b_buf, 32, dimk*t_extb*sizeof(double)*2);
            double* bp = b_buf;
            for (k=0; k<dimk; k++, bp += t_extb*2, b += extb*2)
                memcpy(bp, b, sizeof(double)*dimj*2);
            b = b_buf;
            extb = t_extb;
        }
        vector4double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_0_14, _c_0_15, _c_0_16, _c_0_17, _c_0_18, _c_0_19, _c_0_20, _c_0_21, _c_0_22, _c_0_23, _c_0_24, _c_0_25, _c_0_26, _c_0_27, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _c_1_14, _c_1_15, _c_1_16, _c_1_17, _c_1_18, _c_1_19, _c_1_20, _c_1_21, _c_1_22, _c_1_23, _c_1_24, _c_1_25, _c_1_26, _c_1_27, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_2_8, _c_2_9, _c_2_10, _c_2_11, _c_2_12, _c_2_13, _c_2_14, _c_2_15, _c_2_16, _c_2_17, _c_2_18, _c_2_19, _c_2_20, _c_2_21, _c_2_22, _c_2_23, _c_2_24, _c_2_25, _c_2_26, _c_2_27, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13, _b_0_14, _b_0_15, _b_0_16, _b_0_17, _b_0_18, _b_0_19, _b_0_20, _b_0_21, _b_0_22, _b_0_23, _b_0_24, _b_0_25, _b_0_26, _b_0_27;
        vector4double _a_0_0, _a_0_1, _a_0_2;
        for (j=effj; j>14; j-=14,c+=14*2,b+=14*2) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                _c_0_24 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_1_20 = vec_splats(0.0);
                _c_1_24 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_2_20 = vec_splats(0.0);
                _c_2_24 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _c_1_12 = vec_xmadd(_a_0_1, _b_0_12, _c_1_12);
                    _c_1_12 = vec_xxnpmadd(_b_0_12, _a_0_1, _c_1_12);
                    _c_2_12 = vec_xmadd(_a_0_2, _b_0_12, _c_2_12);
                    _c_2_12 = vec_xxnpmadd(_b_0_12, _a_0_2, _c_2_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _c_1_16 = vec_xmadd(_a_0_1, _b_0_16, _c_1_16);
                    _c_1_16 = vec_xxnpmadd(_b_0_16, _a_0_1, _c_1_16);
                    _c_2_16 = vec_xmadd(_a_0_2, _b_0_16, _c_2_16);
                    _c_2_16 = vec_xxnpmadd(_b_0_16, _a_0_2, _c_2_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                    _c_1_20 = vec_xmadd(_a_0_1, _b_0_20, _c_1_20);
                    _c_1_20 = vec_xxnpmadd(_b_0_20, _a_0_1, _c_1_20);
                    _c_2_20 = vec_xmadd(_a_0_2, _b_0_20, _c_2_20);
                    _c_2_20 = vec_xxnpmadd(_b_0_20, _a_0_2, _c_2_20);
                    _b_0_24 = vec_ld(0, pb+24);
                    _c_0_24 = vec_xmadd(_a_0_0, _b_0_24, _c_0_24);
                    _c_0_24 = vec_xxnpmadd(_b_0_24, _a_0_0, _c_0_24);
                    _c_1_24 = vec_xmadd(_a_0_1, _b_0_24, _c_1_24);
                    _c_1_24 = vec_xxnpmadd(_b_0_24, _a_0_1, _c_1_24);
                    _c_2_24 = vec_xmadd(_a_0_2, _b_0_24, _c_2_24);
                    _c_2_24 = vec_xxnpmadd(_b_0_24, _a_0_2, _c_2_24);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
                vec_st(_c_0_24, 0, c+(i+0)*effj*2+24);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_1_20, 0, c+(i+1)*effj*2+20);
                vec_st(_c_1_24, 0, c+(i+1)*effj*2+24);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
                vec_st(_c_2_20, 0, c+(i+2)*effj*2+20);
                vec_st(_c_2_24, 0, c+(i+2)*effj*2+24);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                _c_0_24 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                    _b_0_24 = vec_ld(0, pb+24);
                    _c_0_24 = vec_xmadd(_a_0_0, _b_0_24, _c_0_24);
                    _c_0_24 = vec_xxnpmadd(_b_0_24, _a_0_0, _c_0_24);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
                vec_st(_c_0_24, 0, c+(i+0)*effj*2+24);
            }
        }
        if (j>12) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                _c_0_24 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_1_20 = vec_splats(0.0);
                _c_1_24 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_2_20 = vec_splats(0.0);
                _c_2_24 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _c_1_12 = vec_xmadd(_a_0_1, _b_0_12, _c_1_12);
                    _c_1_12 = vec_xxnpmadd(_b_0_12, _a_0_1, _c_1_12);
                    _c_2_12 = vec_xmadd(_a_0_2, _b_0_12, _c_2_12);
                    _c_2_12 = vec_xxnpmadd(_b_0_12, _a_0_2, _c_2_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _c_1_16 = vec_xmadd(_a_0_1, _b_0_16, _c_1_16);
                    _c_1_16 = vec_xxnpmadd(_b_0_16, _a_0_1, _c_1_16);
                    _c_2_16 = vec_xmadd(_a_0_2, _b_0_16, _c_2_16);
                    _c_2_16 = vec_xxnpmadd(_b_0_16, _a_0_2, _c_2_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                    _c_1_20 = vec_xmadd(_a_0_1, _b_0_20, _c_1_20);
                    _c_1_20 = vec_xxnpmadd(_b_0_20, _a_0_1, _c_1_20);
                    _c_2_20 = vec_xmadd(_a_0_2, _b_0_20, _c_2_20);
                    _c_2_20 = vec_xxnpmadd(_b_0_20, _a_0_2, _c_2_20);
                    _b_0_24 = vec_ld(0, pb+24);
                    _c_0_24 = vec_xmadd(_a_0_0, _b_0_24, _c_0_24);
                    _c_0_24 = vec_xxnpmadd(_b_0_24, _a_0_0, _c_0_24);
                    _c_1_24 = vec_xmadd(_a_0_1, _b_0_24, _c_1_24);
                    _c_1_24 = vec_xxnpmadd(_b_0_24, _a_0_1, _c_1_24);
                    _c_2_24 = vec_xmadd(_a_0_2, _b_0_24, _c_2_24);
                    _c_2_24 = vec_xxnpmadd(_b_0_24, _a_0_2, _c_2_24);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
                vec_st(_c_0_24, 0, c+(i+0)*effj*2+24);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_1_20, 0, c+(i+1)*effj*2+20);
                vec_st(_c_1_24, 0, c+(i+1)*effj*2+24);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
                vec_st(_c_2_20, 0, c+(i+2)*effj*2+20);
                vec_st(_c_2_24, 0, c+(i+2)*effj*2+24);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                _c_0_24 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                    _b_0_24 = vec_ld(0, pb+24);
                    _c_0_24 = vec_xmadd(_a_0_0, _b_0_24, _c_0_24);
                    _c_0_24 = vec_xxnpmadd(_b_0_24, _a_0_0, _c_0_24);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
                vec_st(_c_0_24, 0, c+(i+0)*effj*2+24);
            }
        }
        else if (j>10) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_1_20 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_2_20 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _c_1_12 = vec_xmadd(_a_0_1, _b_0_12, _c_1_12);
                    _c_1_12 = vec_xxnpmadd(_b_0_12, _a_0_1, _c_1_12);
                    _c_2_12 = vec_xmadd(_a_0_2, _b_0_12, _c_2_12);
                    _c_2_12 = vec_xxnpmadd(_b_0_12, _a_0_2, _c_2_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _c_1_16 = vec_xmadd(_a_0_1, _b_0_16, _c_1_16);
                    _c_1_16 = vec_xxnpmadd(_b_0_16, _a_0_1, _c_1_16);
                    _c_2_16 = vec_xmadd(_a_0_2, _b_0_16, _c_2_16);
                    _c_2_16 = vec_xxnpmadd(_b_0_16, _a_0_2, _c_2_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                    _c_1_20 = vec_xmadd(_a_0_1, _b_0_20, _c_1_20);
                    _c_1_20 = vec_xxnpmadd(_b_0_20, _a_0_1, _c_1_20);
                    _c_2_20 = vec_xmadd(_a_0_2, _b_0_20, _c_2_20);
                    _c_2_20 = vec_xxnpmadd(_b_0_20, _a_0_2, _c_2_20);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_1_20, 0, c+(i+1)*effj*2+20);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
                vec_st(_c_2_20, 0, c+(i+2)*effj*2+20);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_0_20 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _b_0_20 = vec_ld(0, pb+20);
                    _c_0_20 = vec_xmadd(_a_0_0, _b_0_20, _c_0_20);
                    _c_0_20 = vec_xxnpmadd(_b_0_20, _a_0_0, _c_0_20);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_0_20, 0, c+(i+0)*effj*2+20);
            }
        }
        else if (j>8) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _c_1_12 = vec_xmadd(_a_0_1, _b_0_12, _c_1_12);
                    _c_1_12 = vec_xxnpmadd(_b_0_12, _a_0_1, _c_1_12);
                    _c_2_12 = vec_xmadd(_a_0_2, _b_0_12, _c_2_12);
                    _c_2_12 = vec_xxnpmadd(_b_0_12, _a_0_2, _c_2_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                    _c_1_16 = vec_xmadd(_a_0_1, _b_0_16, _c_1_16);
                    _c_1_16 = vec_xxnpmadd(_b_0_16, _a_0_1, _c_1_16);
                    _c_2_16 = vec_xmadd(_a_0_2, _b_0_16, _c_2_16);
                    _c_2_16 = vec_xxnpmadd(_b_0_16, _a_0_2, _c_2_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_xmadd(_a_0_0, _b_0_16, _c_0_16);
                    _c_0_16 = vec_xxnpmadd(_b_0_16, _a_0_0, _c_0_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
            }
        }
        else if (j>6) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                    _c_1_12 = vec_xmadd(_a_0_1, _b_0_12, _c_1_12);
                    _c_1_12 = vec_xxnpmadd(_b_0_12, _a_0_1, _c_1_12);
                    _c_2_12 = vec_xmadd(_a_0_2, _b_0_12, _c_2_12);
                    _c_2_12 = vec_xxnpmadd(_b_0_12, _a_0_2, _c_2_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_xmadd(_a_0_0, _b_0_12, _c_0_12);
                    _c_0_12 = vec_xxnpmadd(_b_0_12, _a_0_0, _c_0_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
            }
        }
        else if (j>4) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                    _c_1_8 = vec_xmadd(_a_0_1, _b_0_8, _c_1_8);
                    _c_1_8 = vec_xxnpmadd(_b_0_8, _a_0_1, _c_1_8);
                    _c_2_8 = vec_xmadd(_a_0_2, _b_0_8, _c_2_8);
                    _c_2_8 = vec_xxnpmadd(_b_0_8, _a_0_2, _c_2_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_xmadd(_a_0_0, _b_0_8, _c_0_8);
                    _c_0_8 = vec_xxnpmadd(_b_0_8, _a_0_0, _c_0_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
            }
        }
        else if (j>2) {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                    _c_1_4 = vec_xmadd(_a_0_1, _b_0_4, _c_1_4);
                    _c_1_4 = vec_xxnpmadd(_b_0_4, _a_0_1, _c_1_4);
                    _c_2_4 = vec_xmadd(_a_0_2, _b_0_4, _c_2_4);
                    _c_2_4 = vec_xxnpmadd(_b_0_4, _a_0_2, _c_2_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_xmadd(_a_0_0, _b_0_4, _c_0_4);
                    _c_0_4 = vec_xxnpmadd(_b_0_4, _a_0_0, _c_0_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
            }
        }
        else {
            for (i=0; i+3<=dimi; i+=3) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _a_0_1 = vec_ld2(0, (pa+2));
                    _a_0_2 = vec_ld2(0, (pa+4));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                    _c_1_0 = vec_xmadd(_a_0_1, _b_0_0, _c_1_0);
                    _c_1_0 = vec_xxnpmadd(_b_0_0, _a_0_1, _c_1_0);
                    _c_2_0 = vec_xmadd(_a_0_2, _b_0_0, _c_2_0);
                    _c_2_0 = vec_xxnpmadd(_b_0_0, _a_0_2, _c_2_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi*2) {
                    _a_0_0 = vec_ld2(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_xmadd(_a_0_0, _b_0_0, _c_0_0);
                    _c_0_0 = vec_xxnpmadd(_b_0_0, _a_0_0, _c_0_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
            }
        }
        /* Copy c out if needed */
        if (dimj%2) {
            double* ct = c_buf;
            for (i=0; i<dimi; i++, ct += effj*2, c_out += dimj*2)
                memcpy(c_out, ct, sizeof(double)*dimj*2);
            free(c_buf);
        }
        /* Free the buffer for b */
        if (free_b) free(b_buf);
    }
    void bgq_mtxmq_padded(long dimi, long dimj, long dimk, long extb, __complex__ double *  c_x,  const __complex__ double *  a_x,  const double  *  b_x) {
        int i, j, k;
        double *  c = (double*)c_x;
        double *  a = (double*)a_x;
        double *  b = (double*)b_x;
        long effj = dimj;
        double * c_buf;
        double * b_buf;
        bool free_b = false;
        /* Setup a buffer for c if needed */
        double* c_out = c;
        if (dimj%2) {
            effj = (dimj | 1) + 1;
            posix_memalign((void **) &c, 32, dimi*effj*sizeof(double)*2);
            c_buf = c;
        }
        /* Copy b into a buffer if needed */
        if (extb%4) {
            long t_extb = (dimj | 3) + 1;
            free_b = true;
            posix_memalign((void **) &b_buf, 32, dimk*t_extb*sizeof(double));
            double* bp = b_buf;
            for (k=0; k<dimk; k++, bp += t_extb, b += extb)
                memcpy(bp, b, sizeof(double)*dimj);
            b = b_buf;
            extb = t_extb;
        }
        vector4double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_4_0, _c_4_1, _c_4_2, _c_4_3, _c_4_4, _c_4_5, _c_4_6, _c_4_7, _c_5_0, _c_5_1, _c_5_2, _c_5_3, _c_5_4, _c_5_5, _c_5_6, _c_5_7, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7;
        vector4double _a_0_0, _a_0_1, _a_0_2, _a_0_3, _a_0_4, _a_0_5;
        vector4double _az_0_0, _az_0_1, _az_0_2, _az_0_3, _az_0_4, _az_0_5;
        vector4double _bz_0_0, _bz_0_1, _bz_0_2, _bz_0_3, _bz_0_4, _bz_0_5, _bz_0_6, _bz_0_7;
        vector4double _cr_perm = vec_gpci(0x9);
        for (i=0; i+6<=dimi; i+=6) {
            double*  xb = b;
            double*  xc = c;
            for (j=effj; j>4; j-=4,xc+=4*2,xb+=4) {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_4_0 = vec_splats(0.0);
                _c_4_4 = vec_splats(0.0);
                _c_5_0 = vec_splats(0.0);
                _c_5_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _az_0_1 = vec_ld2(0, (pa+2));
                    _az_0_2 = vec_ld2(0, (pa+4));
                    _az_0_3 = vec_ld2(0, (pa+6));
                    _az_0_4 = vec_ld2(0, (pa+8));
                    _az_0_5 = vec_ld2(0, (pa+10));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_az_0_1, _bz_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_az_0_2, _bz_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_az_0_3, _bz_0_0, _c_3_0);
                    _c_4_0 = vec_madd(_az_0_4, _bz_0_0, _c_4_0);
                    _c_5_0 = vec_madd(_az_0_5, _bz_0_0, _c_5_0);
                    _bz_0_4 = vec_ld2(0, pb+2);
                    _bz_0_4 = vec_perm(_bz_0_4, _bz_0_4, _cr_perm);
                    _c_0_4 = vec_madd(_az_0_0, _bz_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_az_0_1, _bz_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_az_0_2, _bz_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_az_0_3, _bz_0_4, _c_3_4);
                    _c_4_4 = vec_madd(_az_0_4, _bz_0_4, _c_4_4);
                    _c_5_4 = vec_madd(_az_0_5, _bz_0_4, _c_5_4);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, xc+(i+0)*effj*2+4);
                vec_st(_c_1_0, 0, xc+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, xc+(i+1)*effj*2+4);
                vec_st(_c_2_0, 0, xc+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, xc+(i+2)*effj*2+4);
                vec_st(_c_3_0, 0, xc+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, xc+(i+3)*effj*2+4);
                vec_st(_c_4_0, 0, xc+(i+4)*effj*2+0);
                vec_st(_c_4_4, 0, xc+(i+4)*effj*2+4);
                vec_st(_c_5_0, 0, xc+(i+5)*effj*2+0);
                vec_st(_c_5_4, 0, xc+(i+5)*effj*2+4);
            }
            if (j>2) {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_4_0 = vec_splats(0.0);
                _c_4_4 = vec_splats(0.0);
                _c_5_0 = vec_splats(0.0);
                _c_5_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _az_0_1 = vec_ld2(0, (pa+2));
                    _az_0_2 = vec_ld2(0, (pa+4));
                    _az_0_3 = vec_ld2(0, (pa+6));
                    _az_0_4 = vec_ld2(0, (pa+8));
                    _az_0_5 = vec_ld2(0, (pa+10));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_az_0_1, _bz_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_az_0_2, _bz_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_az_0_3, _bz_0_0, _c_3_0);
                    _c_4_0 = vec_madd(_az_0_4, _bz_0_0, _c_4_0);
                    _c_5_0 = vec_madd(_az_0_5, _bz_0_0, _c_5_0);
                    _bz_0_4 = vec_ld2(0, pb+2);
                    _bz_0_4 = vec_perm(_bz_0_4, _bz_0_4, _cr_perm);
                    _c_0_4 = vec_madd(_az_0_0, _bz_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_az_0_1, _bz_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_az_0_2, _bz_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_az_0_3, _bz_0_4, _c_3_4);
                    _c_4_4 = vec_madd(_az_0_4, _bz_0_4, _c_4_4);
                    _c_5_4 = vec_madd(_az_0_5, _bz_0_4, _c_5_4);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, xc+(i+0)*effj*2+4);
                vec_st(_c_1_0, 0, xc+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, xc+(i+1)*effj*2+4);
                vec_st(_c_2_0, 0, xc+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, xc+(i+2)*effj*2+4);
                vec_st(_c_3_0, 0, xc+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, xc+(i+3)*effj*2+4);
                vec_st(_c_4_0, 0, xc+(i+4)*effj*2+0);
                vec_st(_c_4_4, 0, xc+(i+4)*effj*2+4);
                vec_st(_c_5_0, 0, xc+(i+5)*effj*2+0);
                vec_st(_c_5_4, 0, xc+(i+5)*effj*2+4);
            }
            else {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_4_0 = vec_splats(0.0);
                _c_5_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _az_0_1 = vec_ld2(0, (pa+2));
                    _az_0_2 = vec_ld2(0, (pa+4));
                    _az_0_3 = vec_ld2(0, (pa+6));
                    _az_0_4 = vec_ld2(0, (pa+8));
                    _az_0_5 = vec_ld2(0, (pa+10));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_az_0_1, _bz_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_az_0_2, _bz_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_az_0_3, _bz_0_0, _c_3_0);
                    _c_4_0 = vec_madd(_az_0_4, _bz_0_0, _c_4_0);
                    _c_5_0 = vec_madd(_az_0_5, _bz_0_0, _c_5_0);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
                vec_st(_c_1_0, 0, xc+(i+1)*effj*2+0);
                vec_st(_c_2_0, 0, xc+(i+2)*effj*2+0);
                vec_st(_c_3_0, 0, xc+(i+3)*effj*2+0);
                vec_st(_c_4_0, 0, xc+(i+4)*effj*2+0);
                vec_st(_c_5_0, 0, xc+(i+5)*effj*2+0);
            }
        }
        for (; i+1<=dimi; i+=1) {
            double*  xb = b;
            double*  xc = c;
            for (j=effj; j>4; j-=4,xc+=4*2,xb+=4) {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                    _bz_0_4 = vec_ld2(0, pb+2);
                    _bz_0_4 = vec_perm(_bz_0_4, _bz_0_4, _cr_perm);
                    _c_0_4 = vec_madd(_az_0_0, _bz_0_4, _c_0_4);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, xc+(i+0)*effj*2+4);
            }
            if (j>2) {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                    _bz_0_4 = vec_ld2(0, pb+2);
                    _bz_0_4 = vec_perm(_bz_0_4, _bz_0_4, _cr_perm);
                    _c_0_4 = vec_madd(_az_0_0, _bz_0_4, _c_0_4);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, xc+(i+0)*effj*2+4);
            }
            else {
                double*  pb = xb;
                double*  pa = a+i*2;
                _c_0_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi*2) {
                    _az_0_0 = vec_ld2(0, (pa+0));
                    _bz_0_0 = vec_ld2(0, pb+0);
                    _bz_0_0 = vec_perm(_bz_0_0, _bz_0_0, _cr_perm);
                    _c_0_0 = vec_madd(_az_0_0, _bz_0_0, _c_0_0);
                }
                vec_st(_c_0_0, 0, xc+(i+0)*effj*2+0);
            }
        }
        /* Copy c out if needed */
        if (dimj%2) {
            double* ct = c_buf;
            for (i=0; i<dimi; i++, ct += effj*2, c_out += dimj*2)
                memcpy(c_out, ct, sizeof(double)*dimj*2);
            free(c_buf);
        }
        /* Free the buffer for b */
        if (free_b) free(b_buf);
    }
    void bgq_mtxmq_padded(long dimi, long dimj, long dimk, long extb, __complex__ double *  c_x,  const double  *  a_x,  const __complex__ double *  b_x) {
        int i, j, k;
        double *  c = (double*)c_x;
        double *  a = (double*)a_x;
        double *  b = (double*)b_x;
        long effj = dimj;
        double * c_buf;
        double * b_buf;
        bool free_b = false;
        /* Setup a buffer for c if needed */
        double* c_out = c;
        if (dimj%2) {
            effj = (dimj | 1) + 1;
            posix_memalign((void **) &c, 32, dimi*effj*sizeof(double)*2);
            c_buf = c;
        }
        /* Copy b into a buffer if needed */
        if (extb%2) {
            long t_extb = (dimj | 1) + 1;
            free_b = true;
            posix_memalign((void **) &b_buf, 32, dimk*t_extb*sizeof(double)*2);
            double* bp = b_buf;
            for (k=0; k<dimk; k++, bp += t_extb*2, b += extb*2)
                memcpy(bp, b, sizeof(double)*dimj*2);
            b = b_buf;
            extb = t_extb;
        }
        vector4double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_0_14, _c_0_15, _c_0_16, _c_0_17, _c_0_18, _c_0_19, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _c_1_14, _c_1_15, _c_1_16, _c_1_17, _c_1_18, _c_1_19, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_2_8, _c_2_9, _c_2_10, _c_2_11, _c_2_12, _c_2_13, _c_2_14, _c_2_15, _c_2_16, _c_2_17, _c_2_18, _c_2_19, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_3_8, _c_3_9, _c_3_10, _c_3_11, _c_3_12, _c_3_13, _c_3_14, _c_3_15, _c_3_16, _c_3_17, _c_3_18, _c_3_19, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13, _b_0_14, _b_0_15, _b_0_16, _b_0_17, _b_0_18, _b_0_19;
        vector4double _a_0_0, _a_0_1, _a_0_2, _a_0_3;
        for (j=effj; j>10; j-=10,c+=10*2,b+=10*2) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                _c_3_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                    _c_1_16 = vec_madd(_a_0_1, _b_0_16, _c_1_16);
                    _c_2_16 = vec_madd(_a_0_2, _b_0_16, _c_2_16);
                    _c_3_16 = vec_madd(_a_0_3, _b_0_16, _c_3_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj*2+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj*2+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj*2+12);
                vec_st(_c_3_16, 0, c+(i+3)*effj*2+16);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
            }
        }
        if (j>8) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                _c_3_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                    _c_1_16 = vec_madd(_a_0_1, _b_0_16, _c_1_16);
                    _c_2_16 = vec_madd(_a_0_2, _b_0_16, _c_2_16);
                    _c_3_16 = vec_madd(_a_0_3, _b_0_16, _c_3_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj*2+16);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj*2+16);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj*2+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj*2+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj*2+12);
                vec_st(_c_3_16, 0, c+(i+3)*effj*2+16);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj*2+16);
            }
        }
        else if (j>6) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj*2+12);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj*2+12);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj*2+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj*2+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj*2+12);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj*2+12);
            }
        }
        else if (j>4) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj*2+8);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj*2+8);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj*2+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj*2+8);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj*2+8);
            }
        }
        else if (j>2) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj*2+4);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj*2+4);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj*2+4);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj*2+4);
            }
        }
        else {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
                vec_st(_c_1_0, 0, c+(i+1)*effj*2+0);
                vec_st(_c_2_0, 0, c+(i+2)*effj*2+0);
                vec_st(_c_3_0, 0, c+(i+3)*effj*2+0);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb*2,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj*2+0);
            }
        }
        /* Copy c out if needed */
        if (dimj%2) {
            double* ct = c_buf;
            for (i=0; i<dimi; i++, ct += effj*2, c_out += dimj*2)
                memcpy(c_out, ct, sizeof(double)*dimj*2);
            free(c_buf);
        }
        /* Free the buffer for b */
        if (free_b) free(b_buf);
    }
    void bgq_mtxmq_padded(long dimi, long dimj, long dimk, long extb, double  *  c_x,  const double  *  a_x,  const double  *  b_x) {
        int i, j, k;
        double *  c = (double*)c_x;
        double *  a = (double*)a_x;
        double *  b = (double*)b_x;
        long effj = dimj;
        double * c_buf;
        double * b_buf;
        bool free_b = false;
        /* Setup a buffer for c if needed */
        double* c_out = c;
        if (dimj%4) {
            effj = (dimj | 3) + 1;
            posix_memalign((void **) &c, 32, dimi*effj*sizeof(double));
            c_buf = c;
        }
        /* Copy b into a buffer if needed */
        if (extb%4) {
            long t_extb = (dimj | 3) + 1;
            free_b = true;
            posix_memalign((void **) &b_buf, 32, dimk*t_extb*sizeof(double));
            double* bp = b_buf;
            for (k=0; k<dimk; k++, bp += t_extb, b += extb)
                memcpy(bp, b, sizeof(double)*dimj);
            b = b_buf;
            extb = t_extb;
        }
        vector4double _c_0_0, _c_0_1, _c_0_2, _c_0_3, _c_0_4, _c_0_5, _c_0_6, _c_0_7, _c_0_8, _c_0_9, _c_0_10, _c_0_11, _c_0_12, _c_0_13, _c_0_14, _c_0_15, _c_0_16, _c_0_17, _c_0_18, _c_0_19, _c_1_0, _c_1_1, _c_1_2, _c_1_3, _c_1_4, _c_1_5, _c_1_6, _c_1_7, _c_1_8, _c_1_9, _c_1_10, _c_1_11, _c_1_12, _c_1_13, _c_1_14, _c_1_15, _c_1_16, _c_1_17, _c_1_18, _c_1_19, _c_2_0, _c_2_1, _c_2_2, _c_2_3, _c_2_4, _c_2_5, _c_2_6, _c_2_7, _c_2_8, _c_2_9, _c_2_10, _c_2_11, _c_2_12, _c_2_13, _c_2_14, _c_2_15, _c_2_16, _c_2_17, _c_2_18, _c_2_19, _c_3_0, _c_3_1, _c_3_2, _c_3_3, _c_3_4, _c_3_5, _c_3_6, _c_3_7, _c_3_8, _c_3_9, _c_3_10, _c_3_11, _c_3_12, _c_3_13, _c_3_14, _c_3_15, _c_3_16, _c_3_17, _c_3_18, _c_3_19, _b_0_0, _b_0_1, _b_0_2, _b_0_3, _b_0_4, _b_0_5, _b_0_6, _b_0_7, _b_0_8, _b_0_9, _b_0_10, _b_0_11, _b_0_12, _b_0_13, _b_0_14, _b_0_15, _b_0_16, _b_0_17, _b_0_18, _b_0_19;
        vector4double _a_0_0, _a_0_1, _a_0_2, _a_0_3;
        for (j=effj; j>20; j-=20,c+=20,b+=20) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                _c_3_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                    _c_1_16 = vec_madd(_a_0_1, _b_0_16, _c_1_16);
                    _c_2_16 = vec_madd(_a_0_2, _b_0_16, _c_2_16);
                    _c_3_16 = vec_madd(_a_0_3, _b_0_16, _c_3_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj+16);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj+16);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj+16);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj+12);
                vec_st(_c_3_16, 0, c+(i+3)*effj+16);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj+16);
            }
        }
        if (j>16) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_1_16 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_2_16 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                _c_3_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                    _c_1_16 = vec_madd(_a_0_1, _b_0_16, _c_1_16);
                    _c_2_16 = vec_madd(_a_0_2, _b_0_16, _c_2_16);
                    _c_3_16 = vec_madd(_a_0_3, _b_0_16, _c_3_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj+16);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj+12);
                vec_st(_c_1_16, 0, c+(i+1)*effj+16);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj+12);
                vec_st(_c_2_16, 0, c+(i+2)*effj+16);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj+12);
                vec_st(_c_3_16, 0, c+(i+3)*effj+16);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_0_16 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _b_0_16 = vec_ld(0, pb+16);
                    _c_0_16 = vec_madd(_a_0_0, _b_0_16, _c_0_16);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
                vec_st(_c_0_16, 0, c+(i+0)*effj+16);
            }
        }
        else if (j>12) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_1_12 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_2_12 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                _c_3_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                    _c_1_12 = vec_madd(_a_0_1, _b_0_12, _c_1_12);
                    _c_2_12 = vec_madd(_a_0_2, _b_0_12, _c_2_12);
                    _c_3_12 = vec_madd(_a_0_3, _b_0_12, _c_3_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj+8);
                vec_st(_c_1_12, 0, c+(i+1)*effj+12);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj+8);
                vec_st(_c_2_12, 0, c+(i+2)*effj+12);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj+8);
                vec_st(_c_3_12, 0, c+(i+3)*effj+12);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_0_12 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _b_0_12 = vec_ld(0, pb+12);
                    _c_0_12 = vec_madd(_a_0_0, _b_0_12, _c_0_12);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_0_12, 0, c+(i+0)*effj+12);
            }
        }
        else if (j>8) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_1_8 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_2_8 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                _c_3_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                    _c_1_8 = vec_madd(_a_0_1, _b_0_8, _c_1_8);
                    _c_2_8 = vec_madd(_a_0_2, _b_0_8, _c_2_8);
                    _c_3_8 = vec_madd(_a_0_3, _b_0_8, _c_3_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj+4);
                vec_st(_c_1_8, 0, c+(i+1)*effj+8);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj+4);
                vec_st(_c_2_8, 0, c+(i+2)*effj+8);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj+4);
                vec_st(_c_3_8, 0, c+(i+3)*effj+8);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_0_8 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _b_0_8 = vec_ld(0, pb+8);
                    _c_0_8 = vec_madd(_a_0_0, _b_0_8, _c_0_8);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_0_8, 0, c+(i+0)*effj+8);
            }
        }
        else if (j>4) {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_1_4 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_2_4 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                _c_3_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                    _c_1_4 = vec_madd(_a_0_1, _b_0_4, _c_1_4);
                    _c_2_4 = vec_madd(_a_0_2, _b_0_4, _c_2_4);
                    _c_3_4 = vec_madd(_a_0_3, _b_0_4, _c_3_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_1_4, 0, c+(i+1)*effj+4);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_2_4, 0, c+(i+2)*effj+4);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
                vec_st(_c_3_4, 0, c+(i+3)*effj+4);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_0_4 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _b_0_4 = vec_ld(0, pb+4);
                    _c_0_4 = vec_madd(_a_0_0, _b_0_4, _c_0_4);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_0_4, 0, c+(i+0)*effj+4);
            }
        }
        else {
            for (i=0; i+4<=dimi; i+=4) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                _c_1_0 = vec_splats(0.0);
                _c_2_0 = vec_splats(0.0);
                _c_3_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _a_0_1 = vec_lds(0, (pa+1));
                    _a_0_2 = vec_lds(0, (pa+2));
                    _a_0_3 = vec_lds(0, (pa+3));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_0_1, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_0_2, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_0_3, _b_0_0, _c_3_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
                vec_st(_c_1_0, 0, c+(i+1)*effj+0);
                vec_st(_c_2_0, 0, c+(i+2)*effj+0);
                vec_st(_c_3_0, 0, c+(i+3)*effj+0);
            }
            for (; i+1<=dimi; i+=1) {
                double*  pb = b;
                double*  pa = a+i;
                _c_0_0 = vec_splats(0.0);
                for (k=0; k<dimk; k+=1,pb+=extb,pa+=dimi) {
                    _a_0_0 = vec_lds(0, (pa+0));
                    _b_0_0 = vec_ld(0, pb+0);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                }
                vec_st(_c_0_0, 0, c+(i+0)*effj+0);
            }
        }
        /* Copy c out if needed */
        if (dimj%4) {
            double* ct = c_buf;
            for (i=0; i<dimi; i++, ct += effj, c_out += dimj)
                memcpy(c_out, ct, sizeof(double)*dimj);
            free(c_buf);
        }
        /* Free the buffer for b */
        if (free_b) free(b_buf);
    }
}
