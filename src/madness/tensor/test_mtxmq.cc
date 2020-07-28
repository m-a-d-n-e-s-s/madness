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

#include <madness/madness_config.h>

bool smalltest = false;

// Disable for now to facilitate CI 
//#if !(defined(X86_32) || defined(X86_64))

//#include <iostream>
//int main() {std::cout << "x86 only\n"; return 0;}

//#else

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <xmmintrin.h>

#include <madness/world/safempi.h>
#include <madness/world/posixmem.h>
#include <madness/tensor/cblas.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/mxm.h>

using namespace madness;


#define TIME_DGEMM
#ifdef TIME_DGEMM

void mTxm_dgemm(long ni, long nj, long nk, double* c, const double* a, const double*b ) {
  double one=1.0;
  cblas::gemm(cblas::NoTrans,cblas::Trans,nj,ni,nk,one,b,nj,a,ni,one,c,nj);
}

#endif

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

void crap(double rate, double fastest, double start) {
    if (rate == 0) printf("darn compiler bug %e %e %lf\n",rate,fastest,start);
}


void timer(const char* s, long ni, long nj, long nk, double *a, double *b, double *c) {
  double fastest=0.0, fastest_dgemm=0.0;

  double nflop = 2.0*ni*nj*nk;
  long loop;
  for (int t=0; t<100; t++) {
    double rate;
    double start = SafeMPI::Wtime();
    for (loop=0; loop<100; ++loop) {
      mTxmq(ni,nj,nk,c,a,b);
    }
    start = SafeMPI::Wtime() - start;
    rate = 1.e-9*nflop/(start/100.0);
    crap(rate,fastest,start);
    if (rate > fastest) fastest = rate;
  }
#ifdef TIME_DGEMM
  for (int t=0; t<100; t++) {
    double rate;
    double start = SafeMPI::Wtime();
    for (loop=0; loop<100; ++loop) {
      mTxm_dgemm(ni,nj,nk,c,a,b);
    }
    start = SafeMPI::Wtime() - start;
    rate = 1.e-9*nflop/(start/100.0);
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
  for (int t=0; t<100; t++) {
    double rate;
    double start = SafeMPI::Wtime();
    for (loop=0; loop<100; ++loop) {
      mTxmq(ni,nj,nk,c,a,b);
      mTxmq(ni,nj,nk,a,c,b);
      mTxmq(ni,nj,nk,c,a,b);
    }
    start = SafeMPI::Wtime() - start;
    rate = 1.e-9*nflop/(start/100.0);
    crap(rate,fastest,start);
    if (rate > fastest) fastest = rate;
  }
#ifdef TIME_DGEMM
  for (int t=0; t<100; t++) {
    double rate;
    double start = SafeMPI::Wtime();
    for (loop=0; loop<100; ++loop) {
      mTxm_dgemm(ni,nj,nk,c,a,b);
      mTxm_dgemm(ni,nj,nk,a,c,b);
      mTxm_dgemm(ni,nj,nk,c,a,b);
    }
    start = SafeMPI::Wtime() - start;
    rate = 1.e-9*nflop/(start/100.0);
    crap(rate,fastest_dgemm,start);
    if (rate > fastest_dgemm) fastest_dgemm = rate;
  }
#endif
  printf("%20s %3ld %3ld %3ld %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_dgemm);
}

int main(int argc, char * argv[]) {

    if (getenv("MAD_SMALL_TESTS")) smalltest=true;
    for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
    std::cout << "small test : " << smalltest << std::endl;
    
    const long nimax=!smalltest ? 30*30 : 8*8;
    const long njmax=!smalltest ? 100 : 20;
    const long nkmax=!smalltest ? 100 : 20;
    long ni, nj, nk, i, m;
    double *a, *b, *c, *d;

    SafeMPI::Init_thread(argc, argv, MPI_THREAD_SINGLE);

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

    printf("Starting to test ... \n");
    for (ni=1; ni<std::min(60L,nimax); ni+=1) {
        for (nj=1; nj<std::min(60L,njmax); nj+=1) {
            for (nk=1; nk<std::min(60L,nkmax); nk+=1) {
                for (i=0; i<ni*nj; ++i) d[i] = c[i] = 0.0;
                mTxm (ni,nj,nk,c,a,b);
                mTxmq(ni,nj,nk,d,a,b);
                for (i=0; i<ni*nj; ++i) {
                    double err = std::abs(d[i]-c[i]);
                    /* This test is sensitive to the compilation options.
                       Be sure to have the reference code above compiled
                       -msse2 -fpmath=sse if using GCC.  Otherwise, to
                       pass the test you may need to change the threshold
                       to circa 1e-13.
                    */
                    if (err > 1e-13) {
                        printf("test_mtxmq: error %ld %ld %ld %e\n",ni,nj,nk,err);
                        exit(1);
                    }
                }
            }
        }
    }
    printf("... OK!\n");

    if (!smalltest) {
        printf("%20s %3s %3s %3s %8s %8s (GF/s)\n", "type", "M", "N", "K", "LOOP", "BLAS");
        for (ni=2; ni<60; ni+=2) timer("(m*m)T*(m*m)", ni,ni,ni,a,b,c);
        for (m=2; m<=30; m+=2) timer("(m*m,m)T*(m*m)", m*m,m,m,a,b,c);
        for (m=2; m<=30; m+=2) trantimer("tran(m,m,m)", m*m,m,m,a,b,c);
        for (m=2; m<=20; m+=2) timer("(20*20,20)T*(20,m)", 20*20,m,20,a,b,c);
    }

    SafeMPI::Finalize();

    return 0;
}

//#endif
