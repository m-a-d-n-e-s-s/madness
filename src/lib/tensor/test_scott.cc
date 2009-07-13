#include <tensor/tensor.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <xmmintrin.h>
#include <pmmintrin.h>


#include <tensor/mtxmq.h>

using namespace madness;

void mTxmSCOTT(long dimi, long dimj, long dimk,
           std::complex<double>* restrict c, const std::complex<double>* a, const double* b) {

    const long jtile = 14;
    for (long i=0; i<dimi; i++) {
    	std::complex<double>* cij = c + i*dimj;
    	for (long jlo=0; jlo<dimj; jlo+=jtile, cij+=jtile) {
    		int nj = dimj-jlo;
    		if (nj > jtile) nj = jtile;
    		switch (nj) {
    		case 14:
				{__asm__ __volatile__ (
						"pxor %xmm2,%xmm2;"
						"pxor %xmm3,%xmm3;"
						"pxor %xmm4,%xmm4;"
						"pxor %xmm5,%xmm5;"
						"pxor %xmm6,%xmm6;"
						"pxor %xmm7,%xmm7;"
						"pxor %xmm8,%xmm8;"
						"pxor %xmm9,%xmm9;"
						"pxor %xmm10,%xmm10;"
						"pxor %xmm11,%xmm11;"
						"pxor %xmm12,%xmm12;"
						"pxor %xmm13,%xmm13;"
						"pxor %xmm14,%xmm14;"
						"pxor %xmm15,%xmm15;"
				);
				const std::complex<double>* aki = a + i;
				const double* bkj = b + jlo;
				for (long k=0; k<dimk; k++,aki+=dimi,bkj+=dimj) {
					__asm__ volatile(
						"movapd (%0), %%xmm0;"

						"movddup    (%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;"
						"movddup   8(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;"
						"movddup  16(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;"
						"movddup  24(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;"
						"movddup  32(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;"
						"movddup  40(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;"
						"movddup  48(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;"
						"movddup  56(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;"
						"movddup  64(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;"
						"movddup  72(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;"
						"movddup  80(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;"
						"movddup  88(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13;"
						"movddup  96(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm14;"
						"movddup 104(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm15;"

					:
					: "r"(aki), "r"(bkj)
					:
					);
				}
				__asm__ __volatile__ (
						"movapd   %%xmm2,    (%0);"
						"movapd   %%xmm3,  16(%0);"
						"movapd   %%xmm4,  32(%0);"
						"movapd   %%xmm5,  48(%0);"
						"movapd   %%xmm6,  64(%0);"
						"movapd   %%xmm7,  80(%0);"
						"movapd   %%xmm8,  96(%0);"
						"movapd   %%xmm9, 112(%0);"
						"movapd  %%xmm10, 128(%0);"
						"movapd  %%xmm11, 144(%0);"
						"movapd  %%xmm12, 160(%0);"
						"movapd  %%xmm13, 176(%0);"
						"movapd  %%xmm14, 192(%0);"
						"movapd  %%xmm15, 208(%0);"
				:
				: "r"(cij)
				:
				);}
				break;

    		default:
    		    for (long j=jlo; j<dimj; j++) {
					std::complex<double> cij(0.0,0.0);
					for (long k=0; k<dimk; k++) {
						cij += a[k*dimi+i]*b[k*dimj+j];
					}
					c[i*dimj + j] = cij;
    		    }
    		    break;
    		}
        }
    }



}

int main(int argc, char** argv)
{
	int k = 28;

	Tensor< std::complex<double> > dc(k,k*k);
	Tensor< std::complex<double> > res(k,k*k), res2(k,k*k);
	Tensor< double > d(k,k);

	dc.fillrandom();
	d.fillrandom();

	mTxmq(k*k, k, k, res.ptr(), dc.ptr(), d.ptr());
	mTxmSCOTT(k*k, k, k, res2.ptr(), dc.ptr(), d.ptr());
	print("ERROR IS", (res-res2).normf());

	double start = cpu_time();
	for (int i = 0; i < 10000; ++i)
	{
		mTxmSCOTT(k*k, k, k, res.ptr(), dc.ptr(), d.ptr());
	}
	double tused = cpu_time() - start;
	print("Cpu time used: ", tused, "flops/s: ", 2.0*k*k*k*k*1e4 / tused);

	return 0;
}
