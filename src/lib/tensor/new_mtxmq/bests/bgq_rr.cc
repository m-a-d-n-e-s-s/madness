#include <hwi/include/bqc/A2_inlines.h>
#include <stdlib.h>
#include <string.h>

namespace madness {
    void bgq_mtxm_padded(long dimi, long dimj, long dimk, long ext_b,
            double * restrict c_x, const double * restrict a_x, const double * restrict b_x) {
        bool free_b = false;
        long effj = dimj;
        long i, j, k;

        double *c = (double*)c_x;
        double *a = (double*)a_x;
        double *b = (double*)b_x;

        /* Setup a buffer for c if needed */
        double* c_out = c;
        if (dimj%4) {
            effj = (dimj | 3) + 1;
            c = (double*)malloc(sizeof(double)*dimi*effj);
        }

        /* Copy b into a buffer if needed */
        if (ext_b%4) {
            free_b = true;
            double* b_buf = (double*)malloc(sizeof(double)*dimk*effj);

            double* bp = b_buf;
            for (k=0; k<dimk; k++, bp += effj, b += ext_b)
                memcpy(bp, b, sizeof(double)*dimj);

            b = b_buf;
            ext_b = effj;
        }

        vector4double _c_0_0, _c_1_0, _c_2_0, _c_3_0, _c_4_0, _c_5_0, _c_6_0, _c_7_0;
        vector4double _b_0_0;
        vector4double _a_0_0, _a_1_0, _a_2_0, _a_3_0, _a_4_0, _a_5_0, _a_6_0, _a_7_0; 

        for (i=0; i+8<=dimi; i+= 8) {
            double *xb = b;
            double *xc = c;
            for (j=effj; j>0; j-=4,xc+=4,xb+=4) {
                double *pb = xb;
                double *pa = a+i;
                _c_0_0 = (vector4double) (0.0);
                _c_1_0 = (vector4double) (0.0);
                _c_2_0 = (vector4double) (0.0);
                _c_3_0 = (vector4double) (0.0);
                _c_4_0 = (vector4double) (0.0);
                _c_5_0 = (vector4double) (0.0);
                _c_6_0 = (vector4double) (0.0);
                _c_7_0 = (vector4double) (0.0);
                for (k=0; k<dimk; k+=1, pb+=ext_b, pa+=dimi) {
                    _a_0_0 = vec_lds(sizeof(double)*0, pa);
                    _a_1_0 = vec_lds(sizeof(double)*1, pa);
                    _a_2_0 = vec_lds(sizeof(double)*2, pa);
                    _a_3_0 = vec_lds(sizeof(double)*3, pa);
                    _a_4_0 = vec_lds(sizeof(double)*4, pa);
                    _a_5_0 = vec_lds(sizeof(double)*5, pa);
                    _a_6_0 = vec_lds(sizeof(double)*6, pa);
                    _a_7_0 = vec_lds(sizeof(double)*7, pa);
                    _b_0_0 = vec_ld(0, pb);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_1_0, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_2_0, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_3_0, _b_0_0, _c_3_0);
                    _c_4_0 = vec_madd(_a_4_0, _b_0_0, _c_4_0);
                    _c_5_0 = vec_madd(_a_5_0, _b_0_0, _c_5_0);
                    _c_6_0 = vec_madd(_a_6_0, _b_0_0, _c_6_0);
                    _c_7_0 = vec_madd(_a_7_0, _b_0_0, _c_7_0);
                }
                vec_st(_c_0_0, sizeof(double)*(i+0)*effj, xc);
                vec_st(_c_1_0, sizeof(double)*(i+1)*effj, xc);
                vec_st(_c_2_0, sizeof(double)*(i+2)*effj, xc);
                vec_st(_c_3_0, sizeof(double)*(i+3)*effj, xc);
                vec_st(_c_4_0, sizeof(double)*(i+4)*effj, xc);
                vec_st(_c_5_0, sizeof(double)*(i+5)*effj, xc);
                vec_st(_c_6_0, sizeof(double)*(i+6)*effj, xc);
                vec_st(_c_7_0, sizeof(double)*(i+7)*effj, xc);
            }
        }
        for (; i+6<=dimi; i+= 6) {
            double *xb = b;
            double *xc = c;
            for (j=effj; j>0; j-=4,xc+=4,xb+=4) {
                double *pb = xb;
                double *pa = a+i;
                _c_0_0 = (vector4double) (0.0);
                _c_1_0 = (vector4double) (0.0);
                _c_2_0 = (vector4double) (0.0);
                _c_3_0 = (vector4double) (0.0);
                _c_4_0 = (vector4double) (0.0);
                _c_5_0 = (vector4double) (0.0);
                for (k=0; k<dimk; k+=1, pb+=ext_b, pa+=dimi) {
                    _a_0_0 = vec_lds(sizeof(double)*0, pa);
                    _a_1_0 = vec_lds(sizeof(double)*1, pa);
                    _a_2_0 = vec_lds(sizeof(double)*2, pa);
                    _a_3_0 = vec_lds(sizeof(double)*3, pa);
                    _a_4_0 = vec_lds(sizeof(double)*4, pa);
                    _a_5_0 = vec_lds(sizeof(double)*5, pa);
                    _b_0_0 = vec_ld(0, pb);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_1_0, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_2_0, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_3_0, _b_0_0, _c_3_0);
                    _c_4_0 = vec_madd(_a_4_0, _b_0_0, _c_4_0);
                    _c_5_0 = vec_madd(_a_5_0, _b_0_0, _c_5_0);
                }
                vec_st(_c_0_0, sizeof(double)*(i+0)*effj, xc);
                vec_st(_c_1_0, sizeof(double)*(i+1)*effj, xc);
                vec_st(_c_2_0, sizeof(double)*(i+2)*effj, xc);
                vec_st(_c_3_0, sizeof(double)*(i+3)*effj, xc);
                vec_st(_c_4_0, sizeof(double)*(i+4)*effj, xc);
                vec_st(_c_5_0, sizeof(double)*(i+5)*effj, xc);
            }
        }

        for (; i+4<=dimi; i+= 4) {
            double *xb = b;
            double *xc = c;
            for (j=effj; j>0; j-=4,xc+=4,xb+=4) {
                double *pb = xb;
                double *pa = a+i;
                _c_0_0 = (vector4double) (0.0);
                _c_1_0 = (vector4double) (0.0);
                _c_2_0 = (vector4double) (0.0);
                _c_3_0 = (vector4double) (0.0);
                for (k=0; k<dimk; k+=1, pb+=ext_b, pa+=dimi) {
                    _a_0_0 = vec_lds(sizeof(double)*0, pa);
                    _a_1_0 = vec_lds(sizeof(double)*1, pa);
                    _a_2_0 = vec_lds(sizeof(double)*2, pa);
                    _a_3_0 = vec_lds(sizeof(double)*3, pa);
                    _b_0_0 = vec_ld(0, pb);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_1_0, _b_0_0, _c_1_0);
                    _c_2_0 = vec_madd(_a_2_0, _b_0_0, _c_2_0);
                    _c_3_0 = vec_madd(_a_3_0, _b_0_0, _c_3_0);
                }
                vec_st(_c_0_0, sizeof(double)*(i+0)*effj, xc);
                vec_st(_c_1_0, sizeof(double)*(i+1)*effj, xc);
                vec_st(_c_2_0, sizeof(double)*(i+2)*effj, xc);
                vec_st(_c_3_0, sizeof(double)*(i+3)*effj, xc);
            }
        }
        for (; i+2<=dimi; i+= 2) {
            double *xb = b;
            double *xc = c;
            for (j=effj; j>0; j-=4,xc+=4,xb+=4) {
                double *pb = xb;
                double *pa = a+i;
                _c_0_0 = (vector4double) (0.0);
                _c_1_0 = (vector4double) (0.0);
                for (k=0; k<dimk; k+=1, pb+=ext_b, pa+=dimi) {
                    _a_0_0 = vec_lds(sizeof(double)*0, pa);
                    _a_1_0 = vec_lds(sizeof(double)*1, pa);
                    _b_0_0 = vec_ld(0, pb);
                    _c_0_0 = vec_madd(_a_0_0, _b_0_0, _c_0_0);
                    _c_1_0 = vec_madd(_a_1_0, _b_0_0, _c_1_0);
                }
                vec_st(_c_0_0, sizeof(double)*(i+0)*effj, xc);
                vec_st(_c_1_0, sizeof(double)*(i+1)*effj, xc);
            }
        }

        /* Copy c out if needed */
        if (dimj%4) {
            double* ct = c;
            for (i=0; i<dimi; i++, ct += effj, c_out += dimj)
                memcpy(c_out, ct, sizeof(double)*dimj);

            free(c);
        }

        /* Free the buffer for b */
        if (free_b) free(b);
    }

}
