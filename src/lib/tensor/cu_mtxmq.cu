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
#include <tensor/cu_mtxmq_kernels.cu>
#include <tensor/cu_mtxmq.h>
//#include <world/cuda_streams.h>
#include <string.h>
#include <cublas_v2.h>
//namespace madness {
//#include <cublas_v2.h>
 

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
T* GPUtransfer_buffer(T* CPU_buf, unsigned int offset){
	T *GPU_buf;
	cudaMalloc((void **)&GPU_buf,offset*sizeof(T));
	cudaMemcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T),cudaMemcpyHostToDevice);
	return GPU_buf;
}
template double* GPUtransfer_buffer(double* CPU_buf, unsigned int offset);
template std::complex<double>* GPUtransfer_buffer(std::complex<double>* CPU_buf, unsigned int offset);
template float* GPUtransfer_buffer(float* CPU_buf, unsigned int offset);
template std::complex<float>* GPUtransfer_buffer(std::complex<float>* CPU_buf, unsigned int offset);


template <typename T>
void  CPUtransfer_buffer(T* CPU_buf, T *GPU_buf,unsigned int offset){
	cudaMemcpy((void*)CPU_buf,(void*)GPU_buf,offset*sizeof(T),cudaMemcpyDeviceToHost);
}
template  void  CPUtransfer_buffer(double* CPU_buf, double *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(std::complex<double>* CPU_buf, std::complex<double> *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(float* CPU_buf, float *GPU_buf,unsigned int offset);
template  void  CPUtransfer_buffer(std::complex<float>* CPU_buf, std::complex<float> *GPU_buf,unsigned int offset);


template <typename W>
       void GPUdelete_buffer(W* buf){
	cudaFree(buf);
}
template   void GPUdelete_buffer(double* buf);
template   void GPUdelete_buffer(std::complex<double>* buf);
template   void GPUdelete_buffer(float* buf);
template   void GPUdelete_buffer(std::complex<float>* buf);

template <typename aT, typename bT, typename cT>
    void cu_mTxmqq(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}

  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, std::complex<double> *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<double> *C, std::complex<double> *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){}
  template<>   void cu_mTxmqq(long m, long n,long k, std::complex<float> *C, std::complex<float> *A, std::complex<float> *B,void *GPU_stream,int ndim,long tsize, void *h){}

  template<>   void cu_mTxmqq(long m, long n,long k, float *C, float *A, float *B,void *GPU_stream,int ndim,long tsize, void *h){}
  
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
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(*handle, *stream);
	
	if (ndim ==0)
            ndim=2;
//#pragma unroll ndim-1
     for (int i=0;i<ndim-1;i++){
//do{   
//	cublasSetStream(handle, *stream);
        int b;
        if (i % 2 == 0)
	b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,A,M,&zero,C,N);
        else
	b=cublasDgemm(*handle,CUBLAS_OP_N,CUBLAS_OP_T,N,M,K,&one,B,N,C,M,&zero,A,N);
  //     cudaStreamSynchronize(*stream);
//}while(b==CUBLAS_STATUS_EXECUTION_FAILED);
//	 b=cublasGetError();
        /*
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
        */
	//else
	 // printf("kernel execution success");

//cudaStreamSynchronize(*stream);


        
	//if (tsize==1){
	//	cublasSetStream(handle, *stream);
		//printf("INSIDE SWAP");
		//b=cublasDswap (*handle,M*N, A, 1,C ,1);
                //double * temp = A;
                //A = C;
                //C = temp; 
		//b=cublasGetError();
//		cudaStreamSynchronize(*stream);
                /*
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
                */
		//if (i<ndim-2)
		//stat=cudaMemsetAsync((void*)C,0,M*N*sizeof(double),*stream);

	//}
        
		
}	
//	cublasDestroy(handle);

    }
#endif // MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED

