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

#include <string.h>
#include <cublas_v2.h>
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

template <typename T>
T* GPUtransfer_buffer(T* CPU_buf, unsigned int offset, bool copy){
	T *GPU_buf;
	cudaMalloc((void **)&GPU_buf,offset*sizeof(T));
        if (copy)
	  cudaMemcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T),cudaMemcpyHostToDevice);
	return GPU_buf;
}
template double* GPUtransfer_buffer(double* CPU_buf, unsigned int offset, bool copy);


template <typename T>
void  CPUtransfer_buffer(T* CPU_buf, T *GPU_buf,unsigned int offset){
	cudaMemcpy((void*)CPU_buf,(void*)GPU_buf,offset*sizeof(T),cudaMemcpyDeviceToHost);
}
template  void  CPUtransfer_buffer(double* CPU_buf, double *GPU_buf,unsigned int offset);


template <typename W>
       void GPUdelete_buffer(W* buf){
	cudaFree(buf);
}
template   void GPUdelete_buffer(double* buf);

template <typename T>
T* GPUSimtransfer_buffer(T* CPU_buf, unsigned int offset, bool copy){
	T *GPU_buf;
	GPU_buf = new T[offset];
        if (copy)
	  memcpy((void*)GPU_buf,(void*)CPU_buf,offset*sizeof(T));
	return GPU_buf;
}
template double* GPUSimtransfer_buffer(double* CPU_buf, unsigned int offset, bool copy);


template <typename T>
void  CPUSimtransfer_buffer(T* CPU_buf, T *GPU_buf,unsigned int offset){
	memcpy((void*)CPU_buf,(void*)GPU_buf,offset*sizeof(T));
}
template  void  CPUSimtransfer_buffer(double* CPU_buf, double *GPU_buf,unsigned int offset);


template <typename W>
       void GPUSimdelete_buffer(W* buf){
       delete[] buf;
}
template   void GPUSimdelete_buffer(double* buf);

  void setStream(void * GPU_stream, void * h){  
        cublasHandle_t *handle=(cublasHandle_t *)h;	
	cudaStream_t *stream=(cudaStream_t*)GPU_stream;
	cublasSetStream(*handle, *stream);
  }


template <typename aT, typename bT, typename cT>
    void cu_mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}

  template<>   void cu_mTxmq(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){

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
    void cu_mTxmq_stream(long dimi, long dimj, long dimk,
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,int ndim,long tsize, void  *handle)
{}
  
  template<>   void cu_mTxmq_stream(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,int ndim,long tsize, void *h){

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
    void cu_mTxmq_quarter(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m, int offset)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_quarter(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi, int offset){                               

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
    void cu_mTxmq_Oneblock(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_Oneblock(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

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
    void cu_mTxmq_OneWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
  template<>   void cu_mTxmq_OneWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

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
    void cu_mTxmq_OptCWrite(long dimi, long dimj, long dimk,                   
               cT* restrict c,  aT* a,  bT* b, void *GPU_stream,long prev_m)
{printf("not imp");}
                                                                              
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
                                                                              
  template<>   void cu_mTxmq_OneNoWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

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
                                                                              
  template<>   void cu_mTxmq_hundredOneWrite(long m, long n,long k, double *C, double *A, double *B,void *GPU_stream,long prev_dimi){                               

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
        cublasDaxpy(*handle,n,&s,b,1,a,1);
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


}
template<> void cu_axpystream(long n , float *a, float *b, float s,  void *GPU_stream, void *h) {

        // alpha= reinterpret_cast<T> (s);
        cublasHandle_t *handle=(cublasHandle_t *)h;
//      cublasCreate(&handle);
        cudaStream_t *stream=(cudaStream_t*)GPU_stream;
        cublasSetStream(*handle, *stream);
        cublasSaxpy(*handle,n,&s,b,1,a,1);


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

#endif // MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED

