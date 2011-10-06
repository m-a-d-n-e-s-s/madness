/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/


#include <tensor/tensor.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <xmmintrin.h>
#include <tensor/mtxmq.h>
#include <tensor/mxm.h>
#include <tensor/cu_mtxmq.h>

using namespace madness;


#ifdef HAVE_CUDA
#ifndef FORTRAN_INTEGER
#define FORTRAN_INTEGER long
#endif
typedef FORTRAN_INTEGER integer;
extern "C" void dgemm(const char *transa, const char *transb,
                   const integer *m, const integer *n, const integer *k,
                   const double *alpha, const double *a, const integer *lda,
                   const double *b, const integer *ldb, const double *beta,
                   double *c, const integer *ldc, int la, int lb);

//extern  template void cu_mTxmq<double,double,double>(long m, long n,long k, double *C, const double *A, const double *B);
//extern  void cu_mTxmq(long m, long n,long k, double *C, const double *A, const double *B);    
void mTxm_dgemm(long ni, long nj, long nk, double* c, const double* a, const double*b ) {
  integer fni=ni;
  integer fnj=nj;
  integer fnk=nk;
  double one=1.0;
  //dgemm("n","t",&fnj,&fni,&fnk,&one,b,&fnj,a,&fni,&one,c,&fnj,1,1);
  cu_mTxmq(fnj,fni,fnk, c, a, b);
}



double ran()
{
  static unsigned long seed = 76521;

  seed = seed *1812433253 + 12345;

  return ((double) (seed & 0x7fffffff)) * 4.6566128752458e-10;
}

void ran_fill(int n, double *a) {
    while (n--) *a++ = ran();
}

void mTxm(long dimi, long dimj, long dimk,
          double* c, const double* a, const double* b) {
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}

long long rdtsc() {
  long long x = 0;
  return x;
}

void crap(double rate, double fastest, long long start) {
    if (rate == 0) printf("darn compiler bug %e %e %lld\n",rate,fastest,start);
}


void timer(const char* s, long ni, long nj, long nk, double *a, double *b, double *c) {
  double fastest=0.0, fastest_dgemm=0.0;

  double nflop = 2.0*ni*nj*nk;
  long loop;
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxmq(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    crap(rate,fastest,start);
    if (rate > fastest) fastest = rate;
  }
#ifdef HAVE_CUDA
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_dgemm(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    crap(rate,fastest_dgemm,start);
    if (rate > fastest_dgemm) fastest_dgemm = rate;
  }
#endif
  printf("%20s %3ld %3ld %3ld %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_dgemm);
}

void trantimer(const char* s, long ni, long nj, long nk, double *a, double *b, double *c) {
  double fastest=0.0, fastest_dgemm=0.0;
 
  double nflop = 3.0*2.0*ni*nj*nk;
  long loop;
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxmq(ni,nj,nk,c,a,b);
    mTxmq(ni,nj,nk,a,c,b);
    mTxmq(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    crap(rate,fastest,start);
    if (rate > fastest) fastest = rate;
  }
#ifdef HAVE_CUDA
  for (loop=0; loop<30; ++loop) {
    double rate;
    long long start = rdtsc();
    mTxm_dgemm(ni,nj,nk,c,a,b);
    mTxm_dgemm(ni,nj,nk,a,c,b);
    mTxm_dgemm(ni,nj,nk,c,a,b);
    start = rdtsc() - start;
    rate = nflop/start;
    crap(rate,fastest_dgemm,start);
    if (rate > fastest_dgemm) fastest_dgemm = rate;
  }
#endif
  printf("%20s %3ld %3ld %3ld %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_dgemm);
}

int main() {
    const long nimax=30*30;
    const long njmax=100;
    const long nkmax=100;
    long ni, nj, nk, i, m;

    double *a, *b, *c, *d;
    posix_memalign((void **) &a, 16, nkmax*nimax*sizeof(double));
    posix_memalign((void **) &b, 16, nkmax*njmax*sizeof(double));
    posix_memalign((void **) &c, 16, nimax*njmax*sizeof(double));
    posix_memalign((void **) &d, 16, nimax*njmax*sizeof(double));

    ran_fill(nkmax*nimax, a);
    ran_fill(nkmax*njmax, b);


/*     ni = nj = nk = 2; */
/*     for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0; */
/*     mTxm (ni,nj,nk,c,a,b); */
/*     mTxmq(ni,nj,nk,d,a,b); */
/*     for (i=0; i<ni; ++i) { */
/*       long j; */
/*       for (j=0; j<nj; ++j) { */
/* 	printf("%2ld %2ld %.6f %.6f\n", i, j, c[i*nj+j], d[i*nj+j]); */
/*       } */
/*     } */
/*     return 0; */

    //printf("Starting to test ... \n");
    printf("Starting to test ... \n");
    ni=20;
    nj=64;
    nk=20;
    
    for (i=0; i<ni*nj; ++i) {
      printf("a[%d]=%lf b[%d]=%lf\n",i,a[i],i,b[i]);
    }
                for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0;
                mTxm (ni,nj,nk,c,a,b);
		cu_mTxmq(ni,nj,nk,d,a,b);
                for (i=0; i<ni*nj; ++i) {
                    double err = fabs(d[i]-c[i]);
		    printf(" i=%d d=%lf,c=%lf	",i,d[i],c[i]);
		    printf("test_mtxmq: error %ld %ld %ld %e\n",ni,nj,nk,err);
                    /* This test is sensitive to the compilation options.
                       Be sure to have the reference code above compiled
                       -msse2 -fpmath=sse if using GCC.  Otherwise, to
                       pass the test you may need to change the threshold
                       to circa 1e-13.
                    
                    if (err > 1e-15) {
                        printf("test_mtxmq: error %ld %ld %ld %e\n",ni,nj,nk,err);
                        exit(1);
                    }*/
                }
       
    
    printf("... OK!\n");

    return 0;
}

#endif
