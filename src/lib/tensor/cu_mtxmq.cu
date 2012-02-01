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
  
  Part of the code is adopted from Nvidia CUDA sample code and NOT OPTMIZED
*/
#ifndef MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED
#define MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED


#include <madness_config.h>
#define ENABLE_CUBLAS 1 
#include <tensor/dgemm_madness_kernels.cu>
#include <tensor/sgemm_kernel_magma.cu>
#include <tensor/cu_mtxmq_kernels.cu>
#include <tensor/cu_mtxmq.h>
#include <tensor/cu_axpy.cu>
//#include <world/cuda_streams.h>
#include <string.h>
#include <cublas_v2.h>
//namespace madness {
//#include <cublas_v2.h>
#define BLOCK_SIZE_INTEGRAL 128 

    /// Matrix = Matrix transpose * matrix ... reference implementation

    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  It also supposed
    /// to be fast which it achieves thru restrictions
    ///   * All dimensions even
    ///   * All pointers aligned
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
/*
 template <typename T, typename T, typename tensorT>
 void padwrapper(long dimi, long dimj, const double* pc, T* t0, T* t1, tensorT d,unsigned int i){

	tensorT t = d;
	long nij = dimi*dimj;
	if (IS_ODD(dimi) || IS_ODD(dimj) ||
		IS_UNALIGNED(pc) || IS_UNALIGNED(t0) || IS_UNALIGNED(t1)) {
	    for (long i=0; i<nij; ++i) t0[i] = 0.0;
	    mTxm(dimi, dimj, dimj, t0, t.ptr(), pc);
	    for (int n=1; n<t.ndim(); ++n) {
		for (long i=0; i<nij; ++i) t1[i] = 0.0;
		mTxm(dimi, dimj, dimj, t1, t0, pc);
		std::swap(t0,t1);
	    }
	}
	else {
	    //mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc);
	    print("CUDA KERNEL (dim = ",dimi,",",dimj,")\n");
	    cu_mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc, i);
	    for (int n=1; n<t.ndim(); ++n) {
		//mTxmq(dimi, dimj, dimj, t1, t0, pc);
		cu_mTxmq(dimi, dimj, dimj, t1, t0, pc, i);
		std::swap(t0,t1);
	    }
	}

 }
 */
#define CUDA_CHECK_ERROR(k) \
  {\
    cudaError_t ce = cudaGetLastError();\
    if(ce != cudaSuccess) {\
      printf("%d %s\n",k, cudaGetErrorString(ce));\
      exit(EXIT_FAILURE);\
    }\
  }

 template <typename aT, typename bT, typename cT>
    void cu_mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c, const  aT* a, const bT* b,void *GPU_stream,int ndim,long tsize) {
	//printf("template code");
	const aT *h_A= a;
        const bT *h_B= b;
        cT *h_C= c;
        aT *d_A;
        aT *d_odata;
        bT *d_B, *hb;
        cT *d_C, *hc;
        aT *ha;
        dim3 threads = dim3(BLOCK_SIZE,BLOCK_SIZE);
        dim3 grid = dim3(dimj / threads.x, dimi / threads.y);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        int size_i,size_j,size_k,A_size,B_size,C_size;
        dim3 threads_rem, grid_rem;
        if (dimi%BLOCK_SIZE !=0 || dimj%BLOCK_SIZE!=0 || dimk%BLOCK_SIZE!=0){
        size_i = dimi + (BLOCK_SIZE-(dimi%BLOCK_SIZE));
        size_k = dimk + (BLOCK_SIZE-(dimk%BLOCK_SIZE));
        size_j = dimj + (BLOCK_SIZE-(dimj%BLOCK_SIZE));
        A_size = size_i*size_k*sizeof(aT);
        B_size = size_k*size_j*sizeof(bT);
        C_size = size_i*size_j*sizeof(cT);
        cudaMallocHost((void**)&ha,A_size,cudaHostAllocDefault);
        cudaMallocHost((void**)&hb,B_size,cudaHostAllocDefault);
        cudaMallocHost((void**)&hc,C_size,cudaHostAllocDefault);
        cudaMemset(ha,0,size_i*size_k*sizeof(aT));
        cudaMemset(hb,0,size_j*size_k*sizeof(bT));
        cudaMemset(hc,0,size_i*size_j*sizeof(cT));
        }
        else
        {
        size_i = dimi;
        size_k = dimk;
        size_j = dimj;
	A_size = size_i*size_k*sizeof(aT);
        B_size = size_k*size_j*sizeof(bT);
        C_size = size_i*size_j*sizeof(cT);
        cudaMallocHost((void**)&ha,A_size,cudaHostAllocDefault);
        cudaMallocHost((void**)&hb,B_size,cudaHostAllocDefault);
        cudaMallocHost((void**)&hc,C_size,cudaHostAllocDefault);

        }
        int i,j;
        if (dimi%BLOCK_SIZE !=0 || dimj%BLOCK_SIZE!=0 || dimk%BLOCK_SIZE!=0){


        for ( i=0, j=0;i<dimi*dimk; i++, j++)
        {
          if (i%dimk==0 && i!=0)
            j=j+ (BLOCK_SIZE-(dimk%BLOCK_SIZE));
          cudaMemcpy(&ha[j],&h_A[i],sizeof(aT),cudaMemcpyHostToHost);
        }


        for (  i=0, j=0;i<dimj*dimk; i++, j++)
        {
          if (i%dimj==0 && j!=0)
            j=j+ (BLOCK_SIZE-(dimj%BLOCK_SIZE));
//          hb[j]=h_B[i];
	cudaMemcpy(&hb[j],&h_B[i],sizeof(bT),cudaMemcpyHostToHost);
        }
        }


        cudaMalloc((void**)&d_A, A_size) ;
        cudaMalloc((void**)&d_B, B_size) ;
        cudaMalloc((void**)&d_C, C_size) ;
        cudaMalloc((void**)&d_odata, A_size) ;
        cudaMemcpyAsync(d_A, ha, A_size, cudaMemcpyHostToDevice,*stream) ;
        cudaMemcpyAsync(d_B, hb, B_size, cudaMemcpyHostToDevice,*stream) ;
    
 if (dimi%BLOCK_SIZE !=0 || dimj%BLOCK_SIZE!=0 || dimk%BLOCK_SIZE!=0){
	
        dim3 grid_rem = dim3(size_j / threads.x, size_i / threads.y);
        transposeNoBankConflicts<aT><<<grid_rem, threads,0,*stream>>>(d_odata,d_A, size_k, size_i, 1);
	CUDA_CHECK_ERROR(1);
        grid_rem = dim3(size_j / BLOCK_SIZE, size_k / BLOCK_SIZE);

        matrixMul_coalescing<aT,bT,cT><<< grid_rem, threads,0,*stream >>>(d_C, d_odata, d_B, size_i, size_j);
	CUDA_CHECK_ERROR(2);
        }
        else
        { transposeNoBankConflicts<aT><<<grid, threads,0,*stream>>>(d_A,d_A, dimk, dimi, 1);

        matrixMul_coalescing<aT,bT,cT><<< grid, threads,0,*stream >>>(d_C, d_A, d_B, dimi, dimj);
        }

        // copy result from device to host
        cudaMemcpyAsync((void *)hc,(void *) d_C, C_size, cudaMemcpyDeviceToHost,*stream) ;
          if (dimi%BLOCK_SIZE !=0 || dimj%BLOCK_SIZE!=0 || dimk%BLOCK_SIZE!=0){
       for (  i=0, j=0;i<dimj*dimi; i++, j++)
        {
          if (i%dimj==0 && i!=0)
            j=j+ (BLOCK_SIZE-(dimj%BLOCK_SIZE));
       //   h_C[i]=hc[j];
	cudaMemcpy(&h_C[i],&hc[j],sizeof(cT),cudaMemcpyHostToHost);
        }
          }

        cudaFreeHost(ha);
        cudaFreeHost(hb);
        cudaFreeHost(hc);
        cudaFree(d_odata);
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);

}

template <> void cu_mTxmq(long m, long n,long k, std::complex<double> *C, const std::complex<double> *A, const double *B,void *stream,int ndim,long tsize){}    
#if !ENABLE_CUBLAS 
      template void cu_mTxmq(long dimi, long dimj, long dimk, float*  c,const  float* a, const float* b,void *stream,int ndim,long tsize) ;
     
      template void cu_mTxmq(long m, long n,long k, double *C, const double *A, const double *B,void *stream,int ndim,long tsize);
    
  template <> void cu_mTxmq(long m, long n,long k, std::complex<double> *C,const  std::complex<double> *A, const std::complex<double> *B,void *stream,int ndim,long tsize){}
	 
       
#else

  template<>   void cu_mTxmq(long m, long n,long k, double *C, const double *A, const double *B,void *GPU_stream,int ndim,long tsize){

        double *ha,*hb, *hc; 
	double one=1.0;
	double zero=0.0;
	//printf(" GPU Scublas code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	//cublasStatus_t statt;
	cudaError_t stat;	
	double *devPtrA, *devPtrB, *devPtrC;
        cublasHandle_t handle;	
	cublasCreate(&handle);
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(handle, *stream);
	
	stat = cudaMallocHost ( (void**)&ha,M*K*sizeof(double) ) ;
	if (stat != cudaSuccess) {
	printf ("adevice memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&hb,K*N*sizeof(double) ) ;
	if (stat != cudaSuccess) {
	printf ("bdevice memory allocation failed");
	return ;
	}
	stat = cudaMallocHost ((void**)&hc,M*N*sizeof(double) ) ;
	if (stat != cudaSuccess) {
	printf ("cdevice memory allocation failed");
	return ;
	}
	cudaMalloc((void**)&devPtrA,M*K*sizeof(double));
	cudaMalloc((void**)&devPtrB,N*K*sizeof(double));
	cudaMalloc((void**)&devPtrC,M*N*sizeof(double));

	cudaMemcpy((void *)ha,(void *) A, M*K*sizeof(double), cudaMemcpyHostToHost) ;
	cudaMemcpy((void *)hb,(void *) B, N*K*sizeof(double), cudaMemcpyHostToHost) ;
	
	int b=cublasSetMatrixAsync (M, K, sizeof(double), (void *)ha, M, (void *)devPtrA, M,*stream);
	//int  b=cublasGetError();
//	int b=cublasSetMatrix (M, K, sizeof(double), (void *)ha, M, (void *)devPtrA, M);
//	cudaStreamSynchronize(*stream);
        if (b == CUBLAS_STATUS_INVALID_VALUE)
          printf("CUBLAS_STATUS_INVALID_VALUE");
        else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
          printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("setCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");


	b=cublasSetMatrixAsync (K, N, sizeof(double), (void *)hb, K, (void *)devPtrB, K,*stream);
//	cudaStreamSynchronize(*stream);
	//b=cublasSetMatrix (K, N, sizeof(double), (void *)hb, K, (void *)devPtrB, K);
	 // b=cublasGetError();
        if (b == CUBLAS_STATUS_INVALID_VALUE)
          printf("CUBLAS_STATUS_INVALID_VALUE");
        else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
          printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("setbCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
	//dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);


	//cublasDgemm('t','n',M,N,K,one,devPtrA,K,devPtrB,K,zero,devPtrC,M);
	if (ndim ==0)
            ndim=2;
     for (int i=0;i<ndim-1;i++){
//do{   
//	cublasSetStream(handle, *stream);
	b=cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,devPtrB,N,devPtrA,M,&zero,devPtrC,N);
  //     cudaStreamSynchronize(*stream);
//}while(b==CUBLAS_STATUS_EXECUTION_FAILED);
//	 b=cublasGetError();
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("kernelCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
	//else
	 // printf("kernel execution success");

//cudaStreamSynchronize(*stream);


	if (tsize==1){
	//	cublasSetStream(handle, *stream);
		//printf("INSIDE SWAP");
		b=cublasDswap (handle,M*N, devPtrA, 1,devPtrC ,1);
		//b=cublasGetError();
//		cudaStreamSynchronize(*stream);
		if (b == CUBLAS_STATUS_INVALID_VALUE)
		  printf("CUBLAS_STATUS_INVALID_VALUE");
		else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
		  printf("CUBLAS_STATUS_ARCH_MISMATCH");
		else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
		  printf("swapCUBLAS_STATUS_EXECUTION_FAILED");
		else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
		  printf("CUBLAS_STATUS_MAPPING_ERROR");
		else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
		  printf("CUBLAS_STATUS_ALLOC_FAILED");
		else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
		  printf("init CUBLAS_STATUS_NOT_INITIALIZED");
		else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
		  printf("CUBLAS_STATUS_INTERNAL_ERROR");
		if (i<ndim-2)
		stat=cudaMemsetAsync((void*)devPtrC,0,M*N*sizeof(double),*stream);

	}
		
}	
	//else
	  //printf("Error=%d",b);
	//cublasGetMatrix (K, K, sizeof(double), (void *)devPtrC, K, (void *)C, K);
	//dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
	//cublasSgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,N);
	if(tsize==1)
	b=cublasGetMatrixAsync (M, N, sizeof(double), (void *)devPtrA, M, (void *)hc, M,*stream);
	else
	b=cublasGetMatrixAsync (M, N, sizeof(double), (void *)devPtrC, M, (void *)hc, M,*stream);
	//b=cublasGetMatrix (M, N, sizeof(double), (void *)devPtrC, M, (void *)C, M);
	 if (b == CUBLAS_STATUS_INVALID_VALUE)
          printf("CUBLAS_STATUS_INVALID_VALUE");
        else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
          printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("getCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
//	cudaStreamSynchronize(*stream);
	stat=cudaMemcpy((void *)C,(void *) hc, M*N*sizeof(double), cudaMemcpyHostToHost) ;
	 if (stat != cudaSuccess) {
                        printf ("setdevice memory allocation failed");
                        return ;
                        }

        

	cudaFree (devPtrA);
	cudaFree (devPtrB);
	cudaFree (devPtrC);
	cudaFreeHost (ha);
	cudaFreeHost (hb);
	cudaFreeHost (hc);
	cublasDestroy(handle);

    }


/*	
//cudaSetDevice(3);
	double one=1.0;
	double zero=0.0;
	printf(" GPU cublas code execution m=%d, n%d,k=%d",m,n,k);
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	cublasStatus_t stat,b;
	double *devPtrA, *devPtrB, *devPtrC;
 		
	do{
	stat=cublasCreate(cublasHandle_t *handle);
	}while(	stat!= cudaSuccess);
    b=cublasGetError();
        if (b == CUBLAS_STATUS_INVALID_VALUE)
          printf("CUBLAS_STATUS_INVALID_VALUE");
        else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
          printf("CUBLAS_STATUS_ARCH_MISMATCH");
else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("CUBLAS_STATUS_EXECUTION_FAILED");
else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");

        else
          printf("cudaSuccess");
	
	stat = cublasAlloc (M*K, sizeof(double), (void**)&devPtrA);
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cublasAlloc (K*N, sizeof(double), (void**)&devPtrB);
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cublasAlloc (M*N*4, sizeof(double), (void**)&devPtrC);
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	cublasSetMatrixAsync (M, K, sizeof(double), (void *)A, M, (void *)devPtrA, M,stream[i]);
	cublasSetMatrixAsync (K, N, sizeof(double), (void *)B, K, (void *)devPtrB, K,stream[i]);

//do{	
//	cublasFree (devPtrC);
//	stat = cublasAlloc (M*N, sizeof(double), (void**)&devPtrC);
  //      if (stat != cudaSuccess) {
    //    printf ("device memory allocation failed");
      //  return ;
       // }

	cublasDgemm('t','n',M,N,K,one,devPtrA,M,devPtrB,K,zero,devPtrC,K);
//cudaThreadSynchronize();

	//cublasDgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,N);
	 b=cublasGetError();
//}while(b != cudaSuccess);
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("CUBLAS_STATUS_EXECUTION_FAILED");
else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("CUBLAS_STATUS_NOT_INITIALIZED");
else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");

	else
	  printf("cudaSuccess");

// make sure Dgemm is finished
        cudaError_t cudaErr = cudaThreadSynchronize();
        if( cudaErr != cudaSuccess ) {
    printf( "Dgemm failed on invocation \n" );
        }
	//cublasGetMatrix (K, K, sizeof(double), (void *)devPtrC, K, (void *)C, K);
	//dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
	//cublasDgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,N);
	cublasGetMatrixAsync (M, N, sizeof(double), (void *)devPtrC, M, (void *)C, M,stream[i]);
	cublasFree (devPtrA);
	cublasFree (devPtrB);
	cublasFree (devPtrC);
	cublasDestroy(cublasHandle_t handle);
    }
    
  */  
  template<>   void cu_mTxmq(long m, long n,long k, std::complex<double> *C, const  std::complex<double> *A, const std::complex<double> *B,void *GPU_stream,int ndim,long tsize){
	cuDoubleComplex one;
	one.x=1.0;
	one.y=0.0;
	cuDoubleComplex zero;
	zero.x=0.0;
	zero.y=0.0;
	printf(" complx GPU code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	
	cudaError_t  stat;
	cuDoubleComplex *devPtrA, *devPtrB, *devPtrC;
	cuDoubleComplex *A1=(cuDoubleComplex *)A;
	cuDoubleComplex *B1=(cuDoubleComplex *)B;
	//cuDoubleComplex *C1=(cuDoubleComplex *)C;
	cublasHandle_t handle;
	cublasCreate(&handle);
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(handle, *stream);
	stat = cudaMallocHost ( (void**)&devPtrA,M*K*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrB,K*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrC,M*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	cublasSetMatrixAsync (M, K, sizeof(cuDoubleComplex), (void *)A1, M, (void *)devPtrA, M, *stream);
	cublasSetMatrixAsync (K, N, sizeof(cuDoubleComplex), (void *)B1, K, (void *)devPtrB, K, *stream);
	
	//cublasZgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,K);
	cublasZgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,devPtrB,N,devPtrA,M,&zero,devPtrC,N);
	cublasGetMatrixAsync (N, M, sizeof(*C), (void *)devPtrC, N, (void *)C, N, *stream);
	cudaFree (devPtrA);
	cudaFree (devPtrB);
	cudaFree (devPtrC);
	cublasDestroy( handle);

    }


template<>  void cu_mTxmq(long m, long n,long k,float *C, const float *A, const float *B,void *GPU_stream,int ndim,long tsize){
	float one=1.0;
	float zero=0.0;
	printf(" GPU Scublas code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	cudaError_t stat;
	float *devPtrA, *devPtrB, *devPtrC;
	cublasHandle_t handle;
        cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        cublasSetStream(handle, *stream);	
	stat = cudaMallocHost ( (void**)&devPtrA,M*K*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrB,K*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrC,M*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	cublasSetMatrixAsync (M, K, sizeof(float), (void *)A, M, (void *)devPtrA, M,*stream);
	cublasSetMatrixAsync (K, N, sizeof(float), (void *)B, K, (void *)devPtrB, K,*stream);
	//dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
	//cublasDgemm('t','n',M,N,K,one,devPtrA,K,devPtrB,K,zero,devPtrC,M);
	int b=cublasSgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,devPtrB,N,devPtrA,M,&zero,devPtrC,N);
//	int  b=cublasGetError();
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
	else
	  printf("Error=%d",b);
	//cublasGetMatrix (K, K, sizeof(float), (void *)devPtrC, K, (void *)C, K);
	//dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
	//cublasSgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,N);
	cublasGetMatrixAsync (M, N, sizeof(float), (void *)devPtrC, M, (void *)C, M,*stream);
	cudaFree (devPtrA);
	cudaFree (devPtrB);
	cudaFree (devPtrC);
	cublasDestroy(handle);

    }
    
    
  template<>   void cu_mTxmq(long m, long n,long k, std::complex<float> *C,const std::complex<float> *A, const std::complex<float> *B,void *GPU_stream,int ndim,long tsize){
	cuComplex one;
	one.x=1.0;
	one.y=0.0;
	cuComplex zero;
	zero.x=0.0;
	zero.y=0.0;
	printf(" complx GPU code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	
	cudaError_t  stat;
	cuComplex *devPtrA, *devPtrB, *devPtrC;
	cuComplex *A1=(cuComplex *)A;
	cuComplex *B1=(cuComplex *)B;
	//cuDoubleComplex *C1=(cuDoubleComplex *)C;
	cublasHandle_t handle;
        cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        cublasSetStream(handle, *stream);
	stat = cudaMallocHost ( (void**)&devPtrA,M*K*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrB,K*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	
	stat = cudaMallocHost ((void**)&devPtrC,M*N*sizeof(double),cudaHostAllocDefault ) ;
	if (stat != cudaSuccess) {
	printf ("device memory allocation failed");
	return ;
	}
	cublasSetMatrixAsync (M, K, sizeof(cuComplex), (void *)A1, M, (void *)devPtrA, M,*stream);
	cublasSetMatrixAsync (K, N, sizeof(cuComplex), (void *)B1, K, (void *)devPtrB, K,*stream);
	
	//cublasCgemm('n','t',N,M,K,one,devPtrB,N,devPtrA,M,zero,devPtrC,K);
	cublasCgemm(handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,devPtrB,N,devPtrA,M,&zero,devPtrC,N);
	cublasGetMatrixAsync (N, M, sizeof(*C), (void *)devPtrC, N, (void *)C, N,*stream);
	cudaFree (devPtrA);
	cudaFree (devPtrB);
	cudaFree (devPtrC);
	cublasDestroy( handle);

    }
#endif

//}
template <typename T>
T* GPUtransfer_buffer(T* CPU_buf, unsigned int offset, bool copy){
	T *GPU_buf;
	cudaError_t err = cudaMalloc((void **)&GPU_buf,offset*sizeof(T));
        if (err != cudaSuccess){
          perror("Could not allocate GPU memory   ");
          exit(-1);
        }
        if (copy){
	  err = cudaMemcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T),cudaMemcpyHostToDevice);
          if (err != cudaSuccess){
            perror("Could not memcpy to just allocated GPU memory   ");
            exit(-1);
          }
        }
	return GPU_buf;
}
template double* GPUtransfer_buffer(double* CPU_buf, unsigned int offset, bool copy);
template std::complex<double>* GPUtransfer_buffer(std::complex<double>* CPU_buf, unsigned int offset, bool copy);
template float* GPUtransfer_buffer(float* CPU_buf, unsigned int offset, bool copy);
template std::complex<float>* GPUtransfer_buffer(std::complex<float>* CPU_buf, unsigned int offset, bool copy);
template long* GPUtransfer_buffer(long* CPU_buf, unsigned int offset, bool copy);
template bool* GPUtransfer_buffer(bool* CPU_buf, unsigned int offset, bool copy);

template <typename T>
T* alloc_host(T** CPU_buf, unsigned int size){
	cudaError_t err = cudaMallocHost((void **)CPU_buf,size*sizeof(T));
        if (err != cudaSuccess){
          perror("Could not allocate CPU host memory   ");
          exit(-1);
        }
	return *CPU_buf;
}
template double* alloc_host(double** CPU_buf, unsigned int size);
template std::complex<double>* alloc_host(std::complex<double>** CPU_buf, unsigned int size);
template float* alloc_host(float** CPU_buf, unsigned int size);
template std::complex<float>* alloc_host(std::complex<float>** CPU_buf, unsigned int size);
template long* alloc_host(long** CPU_buf, unsigned int size);
template bool* alloc_host(bool** CPU_buf, unsigned int size);


template <typename T>
T* GPUtransfer_buffernoalloc(T* GPU_buf, T* CPU_buf, unsigned int offset){
	cudaError_t err = cudaMemcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T),cudaMemcpyHostToDevice);
        if (err != cudaSuccess){
          perror("Could not memcpy to GPU memory   ");
          exit(-1);
        }
	return GPU_buf;
}
template double* GPUtransfer_buffernoalloc(double* GPU_buf, double* CPU_buf, unsigned int offset);
template std::complex<double>* GPUtransfer_buffernoalloc(std::complex<double>* GPU_buf, std::complex<double>* CPU_buf, unsigned int offset);
template float* GPUtransfer_buffernoalloc(float* GPU_buf, float* CPU_buf, unsigned int offset);
template std::complex<float>* GPUtransfer_buffernoalloc(std::complex<float>* GPU_buf, std::complex<float>* CPU_buf, unsigned int offset);
template long* GPUtransfer_buffernoalloc(long* GPU_buf, long* CPU_buf, unsigned int offset);


template <typename T>
void  CPUtransfer_buffer(T* CPU_buf, T *GPU_buf,unsigned int offset){
	cudaMemcpy((void*)CPU_buf,(void*)GPU_buf,offset*sizeof(T),cudaMemcpyDeviceToHost);
}
template  void  CPUtransfer_buffer(double* CPU_buf, double *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(std::complex<double>* CPU_buf, std::complex<double> *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(float* CPU_buf, float *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(std::complex<float>* CPU_buf, std::complex<float> *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(long* CPU_buf, long *GPU_buf,unsigned int offset);


template <typename W>
       void GPUdelete_buffer(W* buf){
	cudaFree(buf);
}
template   void GPUdelete_buffer(double* buf);
template   void GPUdelete_buffer(std::complex<double>* buf);
template   void GPUdelete_buffer(float* buf);
template   void GPUdelete_buffer(std::complex<float>* buf);
template   void GPUdelete_buffer(long* buf);
template   void GPUdelete_buffer(bool* buf);

template <typename W>
       void dealloc_host(W* buf){
	cudaFreeHost(buf);
}
template   void dealloc_host(double* buf);
template   void dealloc_host(std::complex<double>* buf);
template   void dealloc_host(float* buf);
template   void dealloc_host(std::complex<float>* buf);
template   void dealloc_host(long* buf);
template   void dealloc_host(bool* buf);

template <typename T>
T* GPUSimtransfer_buffer(T* CPU_buf, unsigned int offset, bool copy){
	T *GPU_buf;
	GPU_buf = new T[offset];
        if (copy)
	  memcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T));
	return GPU_buf;
}
template double* GPUSimtransfer_buffer(double* CPU_buf, unsigned int offset, bool copy);
template std::complex<double>* GPUSimtransfer_buffer(std::complex<double>* CPU_buf, unsigned int offset, bool copy);
template float* GPUSimtransfer_buffer(float* CPU_buf, unsigned int offset, bool copy);
template std::complex<float>* GPUSimtransfer_buffer(std::complex<float>* CPU_buf, unsigned int offset, bool copy);


template <typename T>
void  CPUSimtransfer_buffer(T* CPU_buf, T *GPU_buf,unsigned int offset){
	memcpy((void*)CPU_buf,(void*)GPU_buf,offset*sizeof(T));
}
template  void  CPUSimtransfer_buffer(double* CPU_buf, double *GPU_buf,unsigned int offset);
template  void  CPUSimtransfer_buffer(std::complex<double>* CPU_buf, std::complex<double> *GPU_buf,unsigned int offset);
template  void  CPUSimtransfer_buffer(float* CPU_buf, float *GPU_buf,unsigned int offset);
template  void  CPUSimtransfer_buffer(std::complex<float>* CPU_buf, std::complex<float> *GPU_buf,unsigned int offset);


template <typename W>
       void GPUSimdelete_buffer(W* buf){
       delete[] buf;
}
template   void GPUSimdelete_buffer(double* buf);
template   void GPUSimdelete_buffer(std::complex<double>* buf);
template   void GPUSimdelete_buffer(float* buf);
template   void GPUSimdelete_buffer(std::complex<float>* buf);

template <typename aT, typename bT, typename cT>
    void cu_mTxmqq(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}

  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqq(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,int ndim,long tsize, void *h){}

  void setStream(void * GPU_stream, void * h){  
        cublasHandle_t *handle=(cublasHandle_t *)h;	
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(*handle, *stream);
  }

  template<>   void cu_mTxmqq(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){

       // double *ha,*hb, *hc; 
	double one=1.0;
	double zero=0.0;
	//printf(" GPU Scublas code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	//cublasStatus_t statt;
	cudaError_t stat;	
	//double *devPtrA, *devPtrB, *devPtrC;
        cublasHandle_t *handle=(cublasHandle_t *)h;	
//	cublasCreate(&handle);
	//cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	//cublasSetStream(*handle, *stream);

        	
	if (ndim ==0)
            ndim=2;
//#pragma unroll ndim-1
     for (int i=0;i<ndim-1;i++){
//do{   
//	cublasSetStream(handle, *stream);
        int b;
        //if (i % 2 == 0)
	b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,A,M,&zero,C,N);
        //else
	//b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,C,M,&zero,A,N);
  //     cudaStreamSynchronize(*stream);
//}while(b==CUBLAS_STATUS_EXECUTION_FAILED);
//	 b=cublasGetError();
        
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("kernelCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
        
	//else
	  //printf("kernel execution success");

//cudaStreamSynchronize(*stream);


        
	if (tsize==1 && i < ndim - 2){
	//	cublasSetStream(handle, *stream);
		//printf("INSIDE SWAP");
		b=cublasDswap (*handle,M*N, A, 1,C ,1);
                //double * temp = A;
                //A = C;
                //C = temp; 
		//b=cublasGetError();
//		cudaStreamSynchronize(*stream);
                
		if (b == CUBLAS_STATUS_INVALID_VALUE)
		  printf("CUBLAS_STATUS_INVALID_VALUE");
		else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
		  printf("CUBLAS_STATUS_ARCH_MISMATCH");
		else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
		  printf("swapCUBLAS_STATUS_EXECUTION_FAILED");
		else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
		  printf("CUBLAS_STATUS_MAPPING_ERROR");
		else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
		  printf("CUBLAS_STATUS_ALLOC_FAILED");
		else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
		  printf("init CUBLAS_STATUS_NOT_INITIALIZED");
		else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
		  printf("CUBLAS_STATUS_INTERNAL_ERROR");
                
		//if (i<ndim-2)
		//stat=cudaMemsetAsync((void*)C,0,M*N*sizeof(double),*stream);

	}
        
		
}
//	cublasDestroy(handle);

    }

template <typename aT, typename bT, typename cT>
    void cu_mTxmqnew(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}

  template<>   void cu_mTxmqnew(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqnew(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqnew(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqnew(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqnew(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){

       // double *ha,*hb, *hc; 
	double one=1.0;
	double zero=0.0;
	//printf(" GPU Scublas code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	//cublasStatus_t statt;
	cudaError_t stat;	
	//double *devPtrA, *devPtrB, *devPtrC;
        cublasHandle_t *handle=(cublasHandle_t *)h;	
//	cublasCreate(&handle);
	//cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	//cublasSetStream(*handle, *stream);

        	
        int b;
	b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,A,M,&zero,C,N);
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("kernelCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
        	
}
template <typename aT, typename bT, typename cT>
    void cu_mTxmqnewstream(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}

  template<>   void cu_mTxmqnewstream(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqnewstream(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqnewstream(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqnewstream(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqnewstream(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){

       // double *ha,*hb, *hc; 
	double one=1.0;
	double zero=0.0;
	//printf(" GPU Scublas code execution");
	//sleep(100);
	int M = (int)m;
	int N = (int)n;
	int K = (int)k;
	//cublasStatus_t statt;
	cudaError_t stat;	
	//double *devPtrA, *devPtrB, *devPtrC;
        cublasHandle_t *handle=(cublasHandle_t *)h;	
//	cublasCreate(&handle);
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(*handle, *stream);

        	
        int b;
	b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,A,M,&zero,C,N);
	if (b == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("kernelCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
        	
}
template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%BLOCK_SIZE_INTEGRAL)==0)
                grid=dim3(m/BLOCK_SIZE_INTEGRAL,1,1);                                  
        else
                grid=dim3(m/BLOCK_SIZE_INTEGRAL+1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_20<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral4sep(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m, int offset)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral4sep(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m, int offset){printf("not imp");}
  template<>   void cu_mTxmq_integral4sep(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m, int offset){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral4sep(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m, int offset){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral4sep(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m, int offset){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral4sep(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi, int offset){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_20bl<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0, offset);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}
 
template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral1tb(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral1tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral1tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral1tb(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral1tb(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral1tb(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202hundredtimes<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralOneWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOneWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralOneWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralOneWrite(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOneWrite(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integralOneWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202onewrite<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralOptCWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOptCWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralOptCWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralOptCWrite(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOptCWrite(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integralOptCWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202optCwrite<<<grid,threads,5*512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralOneNoWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOneNoWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralOneNoWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralOneNoWrite(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralOneNoWrite(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integralOneNoWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202onenowrite<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralhundredOneWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralhundredOneWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralhundredOneWrite(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralhundredOneWrite(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralhundredOneWrite(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integralhundredOneWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202hundredonewrite<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral111tb(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral111tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral111tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral111tb(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral111tb(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral111tb(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( /*BLOCK_SIZE_INTEGRAL*/400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202hundrednowrite<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral11tb(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral11tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral11tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral11tb(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral11tb(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral11tb(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( 400,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(1,1,1);                                  
        else
                grid=dim3(1,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_202<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral2tb(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral2tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral2tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral2tb(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral2tb(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral2tb(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(2,1,1);                                  
        else
                grid=dim3(2,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_20hundredtimes<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integral4tb(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral4tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integral4tb(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integral4tb(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integral4tb(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m){printf("not imp");}                                    
  
  template<>   void cu_mTxmq_integral4tb(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%(BLOCK_SIZE_INTEGRAL))==0)
                grid=dim3(4,1,1);                                  
        else
                grid=dim3(4,1,1);
//      STARTt_TIMER;
        switch (k){
        case 2:
                                cu_mtxmq_integral_2<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 4:
                                cu_mtxmq_integral_4<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 6:
                                cu_mtxmq_integral_6<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 8:
                                cu_mtxmq_integral_8<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 10:
                                cu_mtxmq_integral_10<<<grid,threads,4*128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 12:
                                cu_mtxmq_integral_12<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 14:
                                cu_mtxmq_integral_14<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 16:
                                cu_mtxmq_integral_16<<<grid,threads,2*256*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 18:
                                cu_mtxmq_integral_18<<<grid,threads,2*384*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 20:
                                cu_mtxmq_integral_20hundredtimes<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 22:
                                cu_mtxmq_integral_22<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 24:
                                cu_mtxmq_integral_24<<<grid,threads,640*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 26:
                                cu_mtxmq_integral_26<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        case 28:
                                cu_mtxmq_integral_28<<<grid,threads,768*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
        break;
        default:
                                printf("Kernel does not exist for k=%d n=%d",k,n);
                                exit(-1);
        break;
        }

}
 
template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralbatch2(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, cT* restrict c1,  aT* a1,  bT* b1, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralbatch2(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B, std::complex<double> *C1, std::complex<double> *A1, std::complex<double> *B1,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralbatch2(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B, std::complex<double> *C1, std::complex<double> *A1, double *B1,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralbatch2(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B, std::complex<float> *C1, std::complex<float> *A1, std::complex<float> *B1,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralbatch2(long m, long n,long k, float *C, float *A, float *B, float *C1, float *A1, float *B1, void *GPU_stream,long prev_m){printf("not imp");}                                    
 template<>   void cu_mTxmq_integralbatch2(long m, long n,long k, double *C, double *A, double *B, double *C1, double *A1, double *B1, void *GPU_stream,long prev_dimi ){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%BLOCK_SIZE_INTEGRAL)==0)
                grid=dim3(m/BLOCK_SIZE_INTEGRAL,1,1);                                  
        else
                grid=dim3(m/BLOCK_SIZE_INTEGRAL+1,1,1);
                                
         cu_mtxmq_integralbatch2<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,A1,B1,C1,prev_dimi, 0.0,1.0);

}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralbatch3(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, cT* restrict c1,  aT* a1,  bT* b1, cT* restrict c2,  aT* a2,  bT* b2, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralbatch3(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B, std::complex<double> *C1, std::complex<double> *A1, std::complex<double> *B1,std::complex<double> *C2, std::complex<double> *A2, std::complex<double> *B2,void *GPU_stream,long prev_m){printf("not imp");}
  template<>   void cu_mTxmq_integralbatch3(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B, std::complex<double> *C1, std::complex<double> *A1, double *B1,std::complex<double> *C2, std::complex<double> *A2, double *B2,void *GPU_stream,long prev_m){printf("not imp");}     
                                                                              
  template<>   void cu_mTxmq_integralbatch3(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B, std::complex<float> *C1, std::complex<float> *A1, std::complex<float> *B1,std::complex<float> *C2, std::complex<float> *A2, std::complex<float> *B2,void *GPU_stream,long prev_m){printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_integralbatch3(long m, long n,long k, float *C, float *A, float *B, float *C1, float *A1, float *B1, float *C2, float *A2, float *B2,void *GPU_stream,long prev_m){printf("not imp");}                                    
 template<>   void cu_mTxmq_integralbatch3(long m, long n,long k, double *C, double *A, double *B, double *C1, double *A1, double *B1, double *C2, double *A2, double *B2,void *GPU_stream,long prev_dimi ){                               

/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);                                  
exit(-1);
        }*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;                       
        dim3 threads, grid;
        threads=dim3( BLOCK_SIZE_INTEGRAL,1,1 );                                       
        if ((m%BLOCK_SIZE_INTEGRAL)==0)
                grid=dim3(m/BLOCK_SIZE_INTEGRAL,1,1);                                  
        else
                grid=dim3(m/BLOCK_SIZE_INTEGRAL+1,1,1);
                                
         cu_mtxmq_integralbatch3<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,A1,B1,C1,A2,B2,C2,prev_dimi, 0.0,1.0);

}

template <typename T, typename Q>
      void cu_axpy(long n, T* a,  T*  b, Q s, void *GPU_stream, void *h) {

}
template <typename T, typename Q>
      void cu_axpystream(long n, T* a,  T*  b, Q s, void *GPU_stream, void *h) {

}

/*
template<> void cu_axpy(long n , double *a, double *b, double s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        //cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        //cublasSetStream(*handle, *stream);
        //cublasDaxpy(*handle,n,&s,b,1,a,1);
        dim3 threads(512,1,1);
        dim3 grid;
        if ( n==8000){
        grid = dim3 (4,1,1);
        cu_axpy_20<<<grid,threads,0,*stream>>>(n,a,b,s);
        //cu_axpy_20<<<grid,threads,0,0>>>(n,a,b,s);
        }
        else{
        grid = dim3 (1,1,1);
        cu_axpy_10<<<grid,threads,0,*stream>>>(n,a,b,s);
        //cu_axpy_10<<<grid,threads,0,0>>>(n,a,b,s);
        }
//cudaDeviceSynchronize();
//int  f =cudaGetLastError();
//if (f !=cudaSuccess){printf("axpy error= %d",f);
//exit(-1);
//  }
}
*/

template<> void cu_axpy(long n , double *a, double *b, double s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        //cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        //cublasSetStream(*handle, *stream);
        cublasDaxpy(*handle,n,&s,b,1,a,1);
cudaDeviceSynchronize();
int  f =cudaGetLastError();
if (f !=cudaSuccess){printf("axpy error= %d",f);
exit(-1);
  }


}
template<> void cu_axpystream(long n , double *a, double *b, double s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        cublasSetStream(*handle, *stream);
        int b1 = cublasDaxpy(*handle,n,&s,b,1,a,1);
	if (b1 == CUBLAS_STATUS_INVALID_VALUE)
	  printf("CUBLAS_STATUS_INVALID_VALUE");
	else if (b1 == CUBLAS_STATUS_ARCH_MISMATCH)
	  printf("CUBLAS_STATUS_ARCH_MISMATCH");
        else if (b1 ==CUBLAS_STATUS_EXECUTION_FAILED )
          printf("kernelCUBLAS_STATUS_EXECUTION_FAILED");
        else if (b1 ==CUBLAS_STATUS_MAPPING_ERROR )
          printf("CUBLAS_STATUS_MAPPING_ERROR");
        else if (b1 ==CUBLAS_STATUS_ALLOC_FAILED )
          printf("CUBLAS_STATUS_ALLOC_FAILED");
        else if (b1 ==CUBLAS_STATUS_NOT_INITIALIZED )
          printf("init CUBLAS_STATUS_NOT_INITIALIZED");
        else if (b1 ==CUBLAS_STATUS_INTERNAL_ERROR )
          printf("CUBLAS_STATUS_INTERNAL_ERROR");
cudaDeviceSynchronize();
int  f =cudaGetLastError();
if (f !=cudaSuccess){printf("axpy error= %d",f);
exit(-1);
  }


}

template<> void cu_axpy(long n , float *a, float *b, float s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        //cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        //cublasSetStream(*handle, *stream);
        cublasSaxpy(*handle,n,&s,b,1,a,1);

        printf("Shouldn't be here!");

}
template<> void cu_axpystream(long n , float *a, float *b, float s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        cublasSetStream(*handle, *stream);
        cublasSaxpy(*handle,n,&s,b,1,a,1);

        printf("Shouldn't be here!");

}

template<> void cu_axpy(long n , std::complex<double> *a, std::complex<double> *b, std::complex<double> s,  void *GPU_stream, void *h) {
}
template<> void cu_axpystream(long n , std::complex<double> *a, std::complex<double> *b, std::complex<double> s,  void *GPU_stream, void *h) {
}

template<> void cu_axpy(long n , std::complex<float> *a, std::complex<float> *b, std::complex<float> s,  void *GPU_stream, void *h) {
}
template<> void cu_axpystream(long n , std::complex<float> *a, std::complex<float> *b, std::complex<float> s,  void *GPU_stream, void *h) {
}

void lsk(int i){
  dim3 grid = dim3(1,1,1);
  dim3 threads=dim3(1024,1,1);
  for (int j = 0; j < i; j++){
    simple_kernel<<<grid,threads>>>();
  }
}

void lsk1(int i){
  double A[1024];
  dim3 grid = dim3(1,1,1);
  dim3 threads=dim3(1024,1,1);

  double * devA;
  cudaHostAlloc((void **)&devA, 1024*8, 0);

  cudaMemcpy(devA, A, 1024*8, cudaMemcpyHostToDevice);

  //for (int j = 0; j < i; j++){
    simple_kernel1<<<grid,threads>>>(devA,i);
  //}
}

template <typename aT, typename bT, typename cT>
    void cu_mTxmq_integralop(long dimi, long dimj, long dimk,
               cT* /*restrict*/ c,  aT* a,  bT* b, void *GPU_stream,long prev_m, bT* b1, bool *doit, bT* mufac, bT* result, int rank, long *transr, bT* bU)
{}

  template<>   void cu_mTxmq_integralop(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,long prev_m,std::complex<double>* b1, bool *doit, std::complex<double>* mufac,std::complex<double>* result, int rank, long *transr, std::complex<double> *BU){}
  template<>   void cu_mTxmq_integralop(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,long prev_m,double* b1, bool *doit,double *mufac, double *result, int rank, long *transr, double *BU ){}

  template<>   void cu_mTxmq_integralop(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,long prev_m,std::complex<float>* b1, bool *doit,  std::complex<float> *mufac,  std::complex<float> *result, int rank, long *transr, std::complex<float> *BU){}

  template<>   void cu_mTxmq_integralop(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,long prev_m,float* b1, bool *doit, float *mufac, float *result, int rank, long *transr, float *BU){}
  
  template<>   void cu_mTxmq_integralop(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi,double* BVT, bool *doit, double *mufac, double *result, int rank, long *transr, double *BU){


/*int b =cudaGetLastError();
if (b !=cudaSuccess){printf("errpr = %d",b);
exit(-1);
  	}*/
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	dim3 threads, grid;
	/*threads=dim3( BLOCK_SIZE,1,1 );
	if ((m%BLOCK_SIZE)==0)
		grid=dim3(m/BLOCK_SIZE,1,1);
	else
		grid=dim3(m/BLOCK_SIZE+1,1,1);*/
//	STARTt_TIMER;
	switch (k){
	case 10:
				threads=dim3( 128,1,1 );
				grid=dim3(1,1,1);
				//cu_mtxmq_integral_101<<<grid,threads,128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0, BVT, doit, mufac, result, rank, transr);
				cu_mtxmq_integral_110<<<grid,threads,128*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0, BVT, doit, mufac, result, rank, transr, BU);
	break;
        case 20:
				//cu_mtxmq_integral_20<<<grid,threads,512*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0);
				//threads=dim3(160,1,1 );
				//grid=dim3(3,1,1);
				//cu_mtxmq_integral_201<<<grid,threads,400*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0,BVT, doit, mufac, result, rank);
				threads=dim3(128,1,1 );
				grid=dim3(2,1,1);
				cu_mtxmq_integral_211<<<grid,threads,400*8,*stream>>>( A,m,B,n, C, k,prev_dimi, 0.0,1.0,BVT, doit, mufac, result, rank, transr, BU);
	break;
	default:
				printf("Kernel does not exist for k=%d n=%d",k,n);
				exit(-1);
	break;
	}

//cudaDeviceSynchronize();
//int  b =cudaGetLastError();
//if (b !=cudaSuccess){printf("error= %d k=%d m=%d n=%d\n",b,k,m,n);
//exit(-1);
//  }
}

template <typename T>
void fast_cutrans(long dimk, long dimi, T* mat_out, T* mat_in, void * GPU_stream){}

template <>
void fast_cutrans(long dimk, long dimi, float * mat_out, float * mat_in, void * GPU_stream){}

template <>
void fast_cutrans(long dimk, long dimi, std::complex<double> * mat_out, std::complex<double> * mat_in, void * GPU_stream){}

template <>
void fast_cutrans(long dimk, long dimi, std::complex<float> * mat_out, std::complex<float> * mat_in, void * GPU_stream){}

template <>
void fast_cutrans(long dimk, long dimi, double * mat_out, double * mat_in, void * GPU_stream){
    cudaStream_t *stream=(cudaStream_t*)GPU_stream;
    dim3 threads;
    dim3 grid;
    switch (dimk){
        case 20: 
          threads = dim3(256, 1, 1);
          grid = dim3(2, 1, 1);
          cu_transpose_21<<<grid, threads, 0, *stream>>>(mat_out, mat_in);
          break;
        case 10: 
          threads = dim3(128, 1, 1);
          grid = dim3(1, 1, 1);
          cu_transpose_11<<<grid, threads, 0, *stream>>>(mat_out, mat_in);
          break;
        default:
          printf("Only 10 and 20 are supported now for transpose \n");
          exit(-1);
          break;
    }
}

void  cu_memset(){

	int val[8];
	for (int i=0;i<8;i++)
	val[i]=0;
	cudaMemcpyToSymbol(count, val, sizeof(int));
cudaDeviceSynchronize();
int b =cudaGetLastError();
if (b !=cudaSuccess){printf("memset cpy error = %d",b);
exit(-1);
  	}
}
	
#endif // MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED

