/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/* 
 * Kernel of dense matrix-matrix multiplication kernel.
 * The algorithm is based on CUDA sgemm code from Vasily Volkov
 * at UC Berkeley.
 */
#include <stdio.h>
#include <cublas_v2.h>
#define CHECK_ERROR(errorMessage) {                                    \
  cudaError_t err = cudaGetLastError();                                    \
  if( cudaSuccess != err) {                                                \
    fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
	errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
    exit(EXIT_FAILURE);                                                  \
  }                                                                        \
}
volatile __device__ int count[8];

__device__ void start_global_barrier(int x, int id){
	//This synchronization is missing in the report
	__syncthreads();
	if(threadIdx.x == 0){
        atomicAdd((int*)&count[id], 1);
		while( count[id] < x*3){
		    ;
		}
	//count=0;
    }
    __syncthreads();
    
}


__global__ void cu_mtxmq_integral_2(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<2;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<2;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_4(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[4];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<4;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<4;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<4;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_6(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[6];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<6;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<6;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<6;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_8(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[8];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<8;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<8;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<8;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_10(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[10];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<10;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<10;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_10triple(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[10];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<10;i++){
                AL1[i]=A[x+i*lda];
                AL1[i]+=A[x+i*lda];
                AL1[i]-=A[x+i*lda];
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	        Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	        Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
        }
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<10;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		Cs -=AL1[i]*Bs[j+i*prev_dimi];
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<10;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_10fivetimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[10];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<10;i++){
                AL1[i]=A[x+i*lda];
                AL1[i]+=A[x+i*lda];
                AL1[i]-=A[x+i*lda];
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
                for (int k = 0; k < 2; k++){
	            Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	            Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
                }
        }
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<10;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                for (int k = 0; k < 2; k++){
                		    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                		    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
                                 for (int k = 0; k < 2; k++){
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<10;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        for (int k = 0; k < 2; k++){
                			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                                 for (int k = 0; k < 2; k++){
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_10ninetimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[10];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<10;i++){
                AL1[i]=A[x+i*lda];
                AL1[i]+=A[x+i*lda];
                AL1[i]-=A[x+i*lda];
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
                for (int k = 0; k < 4; k++){
	            Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	            Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
                }
        }
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<10;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                for (int k = 0; k < 4; k++){
                		    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                		    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
                                 for (int k = 0; k < 4; k++){
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<10;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        for (int k = 0; k < 4; k++){
                			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                                 for (int k = 0; k < 4; k++){
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_12(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[12];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<12;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<12;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<12;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}


__global__ void cu_mtxmq_integral_14(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[14];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<14;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<14;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<14;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_16(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[16];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<16;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<16;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        		
		}		
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<16;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_18(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[18];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<18;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<18;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<18;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_20(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_20bl(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, int offset)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=offset+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (offset/blockDim.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                for (int ii = 0; ii < 50; ii++){
                		  Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		  Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - offset))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        for (int ii = 0; ii < 50; ii++){
                			  Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			  Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integralbatch2(const double *A, int lda, const double *B, int ldb, double* C, int ldc, const double *A1, const double *B1, double* C1, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
	
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A1[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B1[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C1[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C1[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
}
__global__ void cu_mtxmq_integralbatch3(const double *A, int lda, const double *B, int ldb, double* C, int ldc, const double *A1, const double *B1, double* C1, const double *A2, const double *B2, double* C2, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
	
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A1[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B1[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C1[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C1[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
	
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A2[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B2[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C2[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C2[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
}
__global__ void cu_mtxmq_integral_20twotask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
	
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
}
__global__ void cu_mtxmq_integral_20threetask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
	
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
}
__global__ void cu_mtxmq_integral_20tentask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
	
for (int ii = 0; ii < 9; ii++){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        }
}
__global__ void cu_mtxmq_integral_20hundredtask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
	
for (int ii = 0; ii < 99; ii++){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}
        __syncthreads();
        }
}
__global__ void cu_mtxmq_integral_20triple(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++){
                AL1[i]=A[x+i*lda];
                AL1[i]+=A[x+i*lda];
                AL1[i]-=A[x+i*lda];
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
	       Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	       Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
        }
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		Cs -=AL1[i]*Bs[j+i*prev_dimi];
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_20fivetimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++){
                AL1[i]=A[x+i*lda];
                for (int k = 0; k < 2; k++){
                    AL1[i]+=A[x+i*lda];
                    AL1[i]-=A[x+i*lda];
                }
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
               for (int k = 0; k < 2; k++){
	           Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	           Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
               }
        }
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                for (int k = 0; k < 2; k++){
                		    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
                                 for (int k = 0; k < 2; k++){
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        for (int k = 0; k < 2; k++){
                			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                                 for (int k = 0; k < 2; k++){
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_20ninetimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++){
		AL1[i]=A[x+i*lda];
                for (int k = 0; k < 4; k++){
                    AL1[i]+=A[x+i*lda];
                    AL1[i]-=A[x+i*lda];
                }
        }
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++){
	       Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
               for (int k = 0; k < 4; k++){
	           Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	           Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
               }
        }
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<20;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                for (int k = 0; k < 4; k++){
                		    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
                                 for (int k = 0; k < 4; k++){
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;  // single
				     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;  // single
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        for (int k = 0; k < 4; k++){
                			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                                 for (int k = 0; k < 4; k++){
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
                                 }
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
__global__ void cu_mtxmq_integral_20hundredtimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<20;i++){
		AL1[i]=A[x+i*lda];
                /*
                #pragma unroll 
                for (int k = 0; k < 50; k++){
                    AL1[i]+=A[x+i*lda];
                    AL1[i]-=A[x+i*lda];
                }
                */
        }
        //__syncthreads();
        
	//#pragma unroll 
        //for (int i =0;i<4;i++){
	       Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
               /*
	       #pragma unroll 
               for (int k = 0; k < 50; k++){
	           Bs[i*blockDim.x+threadIdx.x]+=B[i*blockDim.x+ threadIdx.x];
	           Bs[i*blockDim.x+threadIdx.x]-=B[i*blockDim.x+ threadIdx.x];
               }
              */
        //}
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		//#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        //#pragma unroll 
		        for (int  i =0;i<20;i++){ // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
	                        #pragma unroll 
	                        #pragma unroll 
                                for (int k = 0; k < 50; k++){
                		    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                		    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                }
                        }
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[x*ldb+j]=Cs;  // single
                                  
                                 /* 
	                         #pragma unroll
                                 for (int k = 0; k < 50; k++){
				     C[x*ldb+j]+=Cs;  // single
				     C[x*ldb+j]-=Cs;  // single
                                 }
                                 */
                                 
                                 
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<20;i++){ // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];
                                        #pragma unroll 
                                        for (int k = 0; k < 50; k++){
                			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
                			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
                                        }
                                }

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[x*ldb+j]=Cs;
                                 /*
	                         #pragma unroll 
                                 for (int k = 0; k < 50; k++){
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]+=Cs;
                		     C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]-=Cs;
                                 }
                                 */
				 ////C[blockIdx.x*blockDim.x+threadIdx.x+j*lda]=Cs;  // single
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_202(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	{
		#pragma unroll 
		for (int  j=0;j<ldb;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =0;i<2;i++){ // number of rows of B
				Cs +=AL1[i]*Bs[j+i*prev_dimi];
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 C[x*ldb+j]=Cs;
			//__syncthreads();

		}
	}

        for (int ik = 2; ik < 20; ik += 2){
          #pragma unroll 
          for ( int i =0;i<2;i++){
		AL1[i]=A[x+(ik+5)*lda];
                /*
                */
          }

	  if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	  {
		#pragma unroll 
		for (int  j=0;j<ldb;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 C[x*ldb+j] += Cs;
			//__syncthreads();

		}
	  }
        }
}

__global__ void cu_mtxmq_integral_202hundredtimes(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	{
		#pragma unroll 
		for (int  j=0;j<ldb;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =0;i<2;i++){ // number of rows of B
				Cs +=AL1[i]*Bs[j+i*prev_dimi];
				#pragma unroll 
				for (int k = 0; k < 50; k++){
				    Cs +=AL1[i]*Bs[j+i*prev_dimi];
				    Cs -=AL1[i]*Bs[j+i*prev_dimi];
				}
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 C[x*ldb+j]=Cs;
			//__syncthreads();

		}
	}

        for (int ik = 2; ik < 20; ik += 2){
        for (int jk = 0; jk < 20; jk += 2){
          #pragma unroll 
          for ( int i =0;i<2;i++){
		AL1[i]=A[x+(ik+i)*lda];
                /*
                */
          }

	  if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	  {
		#pragma unroll 
		for (int  j=jk;j<jk+2;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				#pragma unroll 
				for (int k = 0; k < 50; k++){
				    Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				    Cs -=AL1[i-ik]*Bs[j+i*prev_dimi];
				}
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 C[x*ldb+j] += Cs;
			//__syncthreads();

		}
	  }
        }
        }
}

__global__ void cu_mtxmq_integral_202hundrednowrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	//if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	//{
		#pragma unroll 
		for (int  j=0;j<ldb;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =0;i<2;i++){ // number of rows of B
				Cs +=AL1[i]*Bs[j+i*prev_dimi];
				#pragma unroll 
				for (int k = 0; k < 50; k++){
				    Cs +=AL1[i]*Bs[j+i*prev_dimi];
				    Cs -=AL1[i]*Bs[j+i*prev_dimi];
				}
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j]=Cs;
			//__syncthreads();

		}
	//}

        for (int ik = 2; ik < 20; ik += 2){
        for (int jk = 0; jk < 20; jk += 2){
          #pragma unroll 
          for ( int i =0;i<2;i++){
		AL1[i]=A[x+(ik+i)*lda];
                /*
                */
          }

	  //if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
	  //{
		#pragma unroll 
		for (int  j=jk;j<jk+2;j++)//number of cols of B
		{
			Cs=0.0;
			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				#pragma unroll 
				for (int k = 0; k < 50; k++){
				    Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				    Cs -=AL1[i-ik]*Bs[j+i*prev_dimi];
				}
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j] += Cs;
			//__syncthreads();

		}
	  //}
        }
        }

        for (int jk = 0; jk < 20; jk++){
            C[x*ldb + jk] = 2.325263;
        }
}

__global__ void cu_mtxmq_integral_202onewrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j++)//number of cols of B
	{
		Cs=0.0;
		#pragma unroll 
		for (int  i =0;i<2;i++){ // number of rows of B
			Cs +=AL1[i]*Bs[j+i*prev_dimi];
		}

		 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
		 //C[x*ldb+j]=Cs;
		//__syncthreads();
		#pragma unroll 
		for (int ik = 2; ik < 20; ik += 2){
			#pragma unroll 
			for ( int i =0;i<2;i++){
		 	  AL1[i]=A[x+(ik+i)*lda];
				/*
				*/
			}

			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j] += Cs;
			//__syncthreads();

		}
                C[x*ldb+j] = Cs;
	}

}

__global__ void cu_mtxmq_integral_202onenowrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j++)//number of cols of B
	{
		Cs=0.0;
		#pragma unroll 
		for (int  i =0;i<2;i++){ // number of rows of B
			Cs +=AL1[i]*Bs[j+i*prev_dimi];
		}

		 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
		 //C[x*ldb+j]=Cs;
		//__syncthreads();
		#pragma unroll 
		for (int ik = 2; ik < 20; ik += 2){
			#pragma unroll 
			for ( int i =0;i<2;i++){
		 	  AL1[i]=A[x+(ik+i)*lda];
				/*
				*/
			}

			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j] += Cs;
			//__syncthreads();

		}
               // C[x*ldb+j] = Cs;
	}

        for (int jk = 0; jk < 20; jk++){
            C[x*ldb + jk] = Cs;
        }
}

__global__ void cu_mtxmq_integral_202optCwrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        extern __shared__ double Cshare[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        //#pragma unroll 
        //for ( int i =0;i<2;i++){
	//	AL1[i]=A[x+i*lda];
                /*
                */
        //}
	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j+=4){//number of cols of B
               Cs=0.0;
	       Cshare[4*x]=0.0;
	       Cshare[4*x+1]=0.0;
	       Cshare[4*x+2]=0.0;
	       Cshare[4*x+3]=0.0;
               for (int  jk=j;jk<j+4;jk++){
                        Cs = 0.0;
			//#pragma unroll 
			//for (int  i =0;i<2;i++){ // number of rows of B
			//	Cs +=AL1[i]*Bs[j+i*prev_dimi];
			//}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j]=Cs;
			//__syncthreads();
			#pragma unroll 
			for (int ik = 0; ik < 20; ik += 2){
				#pragma unroll 
				for ( int i =0;i<2;i++){
				  AL1[i]=A[x+(ik+i)*lda];
					/*
					*/
				}

				#pragma unroll 
				for (int  i =ik;i<ik+2;i++){ // number of rows of B
					Cs +=AL1[i-ik]*Bs[jk+i*prev_dimi];
				}

				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
				 //C[x*ldb+j] += Cs;
				//__syncthreads();

			}
		        Cshare[4*x+jk-j] = Cs;
		}
                __syncthreads();
                for (int  jk=j;jk<j+4;jk++){
                       C[(16*(x/16)+x%4)*ldb+ jk]=Cshare[16*(x/16)+x%4+ 4*(jk-j)];
                }
        }

}

__global__ void cu_mtxmq_integral_202hundredonewrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

	Bs[/*i*blockDim.x+(*/threadIdx.x]=B[/*i*blockDim.x+ */threadIdx.x];
	       // Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j++)//number of cols of B
	{
		#pragma unroll 
		for ( int i =0;i<2;i++){
			AL1[i]=A[x+i*lda];
			/*
			*/
		}
		Cs=0.0;
		#pragma unroll 
		for (int  i =0;i<2;i++){ // number of rows of B
			Cs +=AL1[i]*Bs[j+i*prev_dimi];
			#pragma unroll 
			for (int k = 0; k < 50; k++){
			    Cs +=AL1[i]*Bs[j+i*prev_dimi];
			    Cs -=AL1[i]*Bs[j+i*prev_dimi];
			}
		}

		 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
		 //C[x*ldb+j]=Cs;
		//__syncthreads();
		#pragma unroll 
		for (int ik = 2; ik < 20; ik += 2){
			#pragma unroll 
			for ( int i =0;i<2;i++){
		 	  AL1[i]=A[x+(ik+i)*lda];
				/*
				*/
			}

			#pragma unroll 
			for (int  i =ik;i<ik+2;i++){ // number of rows of B
				Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				#pragma unroll 
				for (int k = 0; k < 50; k++){
				    Cs +=AL1[i-ik]*Bs[j+i*prev_dimi];
				    Cs -=AL1[i-ik]*Bs[j+i*prev_dimi];
				}
			}

			 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
			 //C[x*ldb+j] += Cs;
			//__syncthreads();

		}
                C[x*ldb+j] = Cs;
	}

}

__global__ void cu_mtxmq_integral_22(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[22];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<22;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<22;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<22;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}


__global__ void cu_mtxmq_integral_24(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[24];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<24;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<24;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<24;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}

__global__ void cu_mtxmq_integral_26(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[26];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<26;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<26;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        		
		}	
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<26;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}


__global__ void cu_mtxmq_integral_28(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[28];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<28;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<28;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<28;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb/*+g*blockDim.x*/+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
/*
__global__ void cu_mtxmq_integral(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[SIZE];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        //int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;  // batch code
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<SIZE;i++)
                AL1[i]=A[x+i*prev_dimi];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<4;i++)
	        Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
		#pragma unroll 
        	for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
		        #pragma unroll 
		        for (int  i =0;i<SIZE;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*prev_dimi];
				 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // batched code
				 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  // single
		              __syncthreads();

        	}	
			
	}
	else{
		if (threadIdx.x<(lda - blockDim.x*blockIdx.x))
		{
			#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
		        {
			        Cs=0.0;
			        #pragma unroll 
			        for (int  i =0;i<SIZE;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*prev_dimi];

		                 //C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;  //batched code
                		 C[blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

        //        __syncthreads();
	}	
}
*/

__global__ void simple_kernel(){
 int x=threadIdx.x; 
}

__global__ void simple_kernel1(double * A, int i){
  __shared__ double x; 
  x = A[threadIdx.x];
  for (int j = 1; j < i; j++){
    x = A[threadIdx.x];
  } 
}

__global__ void cu_mtxmq_integral_211( const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, const double* BVT, bool *doit, double *mufac, double *result, int rank, long *transr, const double *BU)
{
        extern __shared__ double Bs[];
        __shared__ volatile double  Cs[2560];
//	double Ctemp;
	double AL1[20];
	//double resultS[20];
	__shared__ volatile bool  doitS[128];
	__shared__ volatile double mufacS[128];
	//__shared__ volatile double  As[2560];
	__shared__  double  resultS[2560];
	__shared__ int trank[512];
//	int num_threads;

	if (threadIdx.x <rank)
	doitS[threadIdx.x]=doit[threadIdx.x];
        __syncthreads();
	if (threadIdx.x <rank)
	mufacS[threadIdx.x]=mufac[threadIdx.x];
        __syncthreads();
	
     	#pragma unroll 
	for (int i=0;i<4;i++)	{
	int index =threadIdx.x+i*128;
	if ( index <rank*3)
		trank[index]=transr[index];
	}
        __syncthreads();
//        int x=blockIdx.x*136+threadIdx.x;   // single 

/*      if (blockIdx.x == (gridDim.x-1) ) 
	num_threads=72;
      else
	num_threads=72;
  */     
     	/*#pragma unroll 
        for ( int i =0;i<20;i++)
                As[threadIdx.x*20+i]=A[blockIdx.x*200+threadIdx.x+i*lda];
        __syncthreads();*/
     	#pragma unroll 
        for (int i=0;i<20;i++)
                resultS[threadIdx.x+i*128]=0.0;
        __syncthreads();


      for (int mu=0;mu<rank;mu++){
//int mu=0;
	if (mufacS[mu]){
 
	//if (threadIdx.x==0)
//printf("rank1=%d \n",trank[mu]);	
 //       __syncthreads();
                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400){
			if (trank[mu*3]== ldb)
			Bs[index]=B[index+mu*1200];
			else
			Bs[index]=BU[index+mu*1200];
		}
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=A[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(1,prev_dimi);
     	#pragma unroll 
	for ( int i =0;i<20;i++)
	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
        __syncthreads();
        

	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	C[blockIdx.x*200*20+128*20+index]=Cs[index];
	}
        __syncthreads();
	
	start_global_barrier(2,prev_dimi);

//	if(threadIdx.x==0)
//	printf("blockIdx.x=%d count=%d ",blockIdx.x,count);
  //      __syncthreads();
        
                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400){
			if (trank[mu*3+1]== ldb)
			Bs[index]=B[index+mu*1200+400];
			else
			Bs[index]=BU[index+mu*1200+400];
		}
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(3,prev_dimi);
     	#pragma unroll 
	for ( int i =0;i<20;i++)
	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
        __syncthreads();
        
	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	C[blockIdx.x*200*20+128*20+index]=Cs[index];
	}
        __syncthreads();
	
	start_global_barrier(4,prev_dimi);
//	if(threadIdx.x==0)
//	printf("blockIdx.x=%d count=%d ",blockIdx.x,count);
  //      __syncthreads();
//
                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400){
			if (trank[mu*3+2]== ldb)
			Bs[index]=B[index+mu*1200+800];
			else
			Bs[index]=BU[index+mu*1200+800];
		}
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(5,prev_dimi);

     if (doitS[mu]==0){	
        #pragma unroll 
        for (int i=0;i<20;i++){
		int index =threadIdx.x*20+i;
                resultS[index]= resultS[index]+Cs[index]*mufacS[mu];
        }__syncthreads();
      }
      else{
     	#pragma unroll 
	for ( int i =0;i<20;i++)
	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
        __syncthreads();
      }  
	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     if (doitS[mu]){	
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	C[blockIdx.x*200*20+128*20+index]=Cs[index];
	}
        __syncthreads();
	}
	start_global_barrier(6,prev_dimi);
	

	if (threadIdx.x ==0){
	   if(count[prev_dimi] ==12)
		count[prev_dimi]=0;	
	}
        __syncthreads();

//	printf("outside doit\n");
       
     if (doitS[mu]){	
	
                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400)
			Bs[index]=BVT[index+mu*1200];
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(1,prev_dimi);
     	#pragma unroll 
	for ( int i =0;i<20;i++)
	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
        __syncthreads();
        
	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	C[blockIdx.x*200*20+128*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
	}
        __syncthreads();

	start_global_barrier(2,prev_dimi);
	
                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400)
			Bs[index]=BVT[index+mu*1200+400];
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(3,prev_dimi);
     	#pragma unroll 
	for ( int i =0;i<20;i++)
	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
        __syncthreads();
        
	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	C[blockIdx.x*200*20+128*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
	}
        __syncthreads();

	start_global_barrier(4,prev_dimi);

                        #pragma unroll 
	for (int i=0;i<4;i++){
		int index=threadIdx.x + i*128;
		if (index <400)
			Bs[index]=BVT[index+mu*1200+800];
	}
        __syncthreads();


 
	//double AL2[20];
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+threadIdx.x+i*lda];

                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;

                        }
        __syncthreads();
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=C[blockIdx.x*200+128+threadIdx.x+i*lda];
	}
        __syncthreads();
	start_global_barrier(5,prev_dimi);
        #pragma unroll 
        for (int i=0;i<20;i++){
		int index =threadIdx.x*20+i;
                resultS[index]= resultS[index]+Cs[index]*mufacS[mu];
	}
        __syncthreads();
//     	#pragma unroll 
//	for ( int i =0;i<20;i++)
//	C[blockIdx.x*200*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
  //      __syncthreads();
        
	if (threadIdx.x<72){
                        #pragma unroll 
                        for (int  j=0;j<20;j++)//number of cols of B
                        {
                                double Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<20;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
				Cs[threadIdx.x*20+j]=Ctemp;
       // __syncthreads();

                        }
	}
        __syncthreads();
     //	#pragma unroll 
//	for ( int i =0;i<12;i++){
//	int index=threadIdx.x + i*128;
//	if (index<1440)
//	C[blockIdx.x*200*20+128*20+threadIdx.x+i*128]=Cs[threadIdx.x+i*128];
//	}
  //      __syncthreads();
     
	start_global_barrier(6,prev_dimi);
	
	if (threadIdx.x ==0){
	   if(count[prev_dimi] ==12)
		count[prev_dimi]=0;	
	}
        __syncthreads();
//	if(threadIdx.x==0)
//	printf("inside doit\n");
  //      __syncthreads();
}
	if (threadIdx.x<72){
        #pragma unroll 
        for ( int i =0;i<20;i++)
                AL1[i]=result[blockIdx.x*200*20+128*20+threadIdx.x+i*72];
  	#pragma unroll 
	for (int i=0;i<20;i++)
                Cs[threadIdx.x+i*72]= AL1[i]+Cs[threadIdx.x+i*72]*mufacS[mu];
        }
                __syncthreads();
     	#pragma unroll 
	for ( int i =0;i<12;i++){
	int index=threadIdx.x + i*128;
	if (index<1440)
	result[blockIdx.x*200*20+128*20+index]=Cs[index];
	}
                __syncthreads();
}
}
        #pragma unroll 
	for (int i=0;i<20;i++)
                //result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
	result[blockIdx.x*200*20+threadIdx.x+i*128]=resultS[threadIdx.x+i*128];
//                result[(blockIdx.x*200+threadIdx.x)*20+i]= resultS[i];
        __syncthreads();

}


__global__ void cu_mtxmq_integral_110(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, const double* BVT, bool *doit, double *mufac, double *result, int rank, long *transr, const double *BU)
{
        extern __shared__ double Bs[];
        __shared__ double Cs[1024];
	 double AL1[10];
	double Ctemp;
	__shared__ double resultS[1024];
	__shared__ bool  doitS[128];
	__shared__ double mufacS[128];
	__shared__ double As[1024];
	__shared__ int trank[512];
	if (threadIdx.x <rank)
	doitS[threadIdx.x]=doit[threadIdx.x];
        __syncthreads();
	if (threadIdx.x <rank)
	mufacS[threadIdx.x]=mufac[threadIdx.x];
        __syncthreads();
	
       
     	#pragma unroll 
	for (int i=0;i<4;i++)	{
	int index =threadIdx.x+i*128;
	if ( index <rank*3)
		trank[index]=transr[index];
	}
        __syncthreads();
      //  if (threadIdx.x<100){
     	#pragma unroll 
        for ( int i =0;i<8;i++)
                As[threadIdx.x+i*128]=A[threadIdx.x+i*128];
//	}
        __syncthreads();
//	if (threadIdx.x<lda){
     	#pragma unroll 
        for (int i=0;i<8;i++)
                //resultS[threadIdx.x+lda*i]=result[threadIdx.x+lda*i];
                resultS[threadIdx.x+128*i]=0.0;
//	}
        __syncthreads();


      for (int mu=0;mu<rank;mu++){
	if (mufacS[mu]){
 
//	if (threadIdx.x==0)
//printf("rank=%d \n",trank[mu]);	
  //      __syncthreads();
	if (trank[mu*3]== ldb){
        	if (threadIdx.x<100)
			Bs[threadIdx.x]=B[threadIdx.x+mu*300];
	}
	else{
        if (threadIdx.x<100)
			Bs[threadIdx.x]=BU[threadIdx.x+mu*300];
	}
        __syncthreads();


 
        if (threadIdx.x<lda){
        for ( int i =0;i<10;i++)
                AL1[i]=As[threadIdx.x+i*lda];
        }__syncthreads();
// lda replaced by dim/trans.r
//ldb replaced by trans.r
// 10 in the loop replaced by trans.r 

        if (threadIdx.x<lda){

                        #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
                                Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;
			}
        }
        __syncthreads();

	
	if (trank[mu*3+1]== ldb){
        if (threadIdx.x<100)
	        Bs[threadIdx.x]=*(B+threadIdx.x+mu*300+100);
    }
    else{
        if (threadIdx.x<100)
			Bs[threadIdx.x]=*(BU+threadIdx.x+mu*300+100);
	}
        __syncthreads();
        
if (threadIdx.x<lda)        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*lda];
        }__syncthreads();

        if (threadIdx.x<lda){

		       #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
                                Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;	
			}
	                       
                }
        __syncthreads();


	if (trank[mu*3+2]== ldb){
        if (threadIdx.x<100)
        Bs[threadIdx.x]=*(B+threadIdx.x+mu*300+200);
	}
    else{
        if (threadIdx.x<100)
			Bs[threadIdx.x]=*(BU+threadIdx.x+mu*300+200);
	}
        __syncthreads();

        if (threadIdx.x<lda)        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*lda];
        }__syncthreads();

        if (threadIdx.x<lda){

                        #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
				Ctemp=0.0;
     				#pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;
			}
                        
               }
        __syncthreads();

        
     if (doitS[mu]){	
        if (threadIdx.x<100)
	Bs[threadIdx.x]=BVT[threadIdx.x+mu*300];
        __syncthreads();
       
        if (threadIdx.x<lda){
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*lda];
        }__syncthreads();



        if (threadIdx.x<lda){
                        #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
                                Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;
			}
                        
	}
        __syncthreads();

	
        if (threadIdx.x<100)
        Bs[threadIdx.x]=*(BVT+threadIdx.x+mu*300+100);
        __syncthreads();
        
if (threadIdx.x<lda)
        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*lda];
        }__syncthreads();


	if (threadIdx.x<lda){
		       #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
                                Ctemp=0.0;
                                #pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;	
			}
                        
                }
        __syncthreads();




        if (threadIdx.x<100)
        Bs[threadIdx.x]=*(BVT+threadIdx.x+mu*300+200);
        __syncthreads();

        if (threadIdx.x<lda)
        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*lda];
        }__syncthreads();

        if (threadIdx.x<lda){
                        #pragma unroll 
                        for (int  j=0;j<ldb;j++)//number of cols of B
                        {
				Ctemp=0.0;
     				#pragma unroll 
                                for (int  i =0;i<10;i++) // number of rows of B
                                        Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
				Cs[threadIdx.x*ldb+j]=Ctemp;

                        }
                }
        __syncthreads();

}

      //  if (threadIdx.x<lda){
  	#pragma unroll 
	for (int i=0;i<8;i++)
                resultS[threadIdx.x+i*128]= resultS[threadIdx.x+i*128]-Cs[threadIdx.x+i*128]*mufacS[mu];
                __syncthreads();
        //}
          //      __syncthreads();
}
}
        if (threadIdx.x<lda){
        #pragma unroll 
	for (int i=0;i<10;i++)
                //result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
                result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
        }
                __syncthreads();

}

__global__ void cu_transpose_11(double *odata, const double  *idata){
double AL[10];
__shared__ double  AS[1000];

if (threadIdx.x<100){
#pragma unroll 
        for ( int i =0;i<10;i++)
                AL[i]=idata[threadIdx.x+i*100];
}
        __syncthreads();

if (threadIdx.x<100){
#pragma unroll 
        for ( int i =0;i<10;i++)

AS[threadIdx.x*10+i] = AL[i];

}
        __syncthreads();

for (int i=0;i<10;i++){
int index = threadIdx.x + i*128;
if (index<1000)
odata[index]=AS[index];
}
}
__global__ void cu_transpose_21(double *odata, const double  *idata){
double AL[20];
__shared__ double  AS[4000];

if (threadIdx.x<200){
#pragma unroll 
        for ( int i =0;i<20;i++)
                AL[i]=idata[blockIdx.x*200+threadIdx.x+i*400];
}
        __syncthreads();

if (threadIdx.x<200){
#pragma unroll 
        for ( int i =0;i<20;i++)

AS[threadIdx.x*20+i] = AL[i];

}
        __syncthreads();

for (int i=0;i<20;i++){
int index = threadIdx.x + i*256;
if (index<4000)
odata[blockIdx.x*4000+index]=AS[index];
}

}





__global__ void cu_mtxmq_integral_204( const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, const double* BVT,const int8_t *doit,const double *mufacS, double *result, int rank,const long *transr, const double *BU)
{
   extern __shared__ double Bs[];
   __shared__  double  Cs[2720];
   double AL1[20];
   //	double Ctemp;
   double resultS[20];
   __shared__  int8_t  doitS[112];
   //	__shared__  double mufacS[112];
   __shared__  double  As[2720];
   __shared__ int trank[336];
   int num_threads;

   if (threadIdx.x <rank)
      doitS[threadIdx.x]=doit[threadIdx.x];
   __syncthreads();
   //	if (threadIdx.x <rank)
   //	mufacS[threadIdx.x]=mufac[threadIdx.x];
   //      __syncthreads();

#pragma unroll 
   for (int i=0;i<4;i++)   {
      int index =threadIdx.x+i*160;
      if ( index <rank*3)
	 trank[index]=transr[index];
   }
   __syncthreads();

   //        int x=blockIdx.x*136+threadIdx.x;   // single 

   if (blockIdx.x == (gridDim.x-1) ) 
      num_threads=128;
   else
      num_threads=136;

   int index= blockIdx.x*136+threadIdx.x;  
   if (threadIdx.x<num_threads){
#pragma unroll 
      for ( int i =0;i<20;i++)
	 As[threadIdx.x*20+i]=A[index+i*400];
   }
   __syncthreads();
#pragma unroll 
   for (int i=0;i<20;i++)
      resultS[i]=0.0;
   __syncthreads();


   for (int mu=0;mu<rank;mu++){
      //int mu=0;
      if (mufacS[mu]){

#pragma unroll 
	 for (int i=0;i<3;i++){
	    index=threadIdx.x + i*160;
	    if (index < 400){
	       if (trank[mu*3]== 20)
		  Bs[index]=B[index+mu*1200];
	       else
		  Bs[index]=BU[index+mu*1200];
	    }
	 }
	 __syncthreads();



	 if (threadIdx.x<num_threads){
#pragma unroll 
	    for ( int i =0;i<20;i++)
	       AL1[i]=As[i+threadIdx.x*20];
	 } // __syncthreads();
	 __syncthreads();



	 if (threadIdx.x<num_threads){
#pragma unroll 
	    for (int  j=0;j<20;j++)//number of cols of B
	    {
	       double Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<20;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
	       Cs[threadIdx.x*20+j]=Ctemp;
	       // __syncthreads();

	    }
	 }
	 __syncthreads();
	 index= blockIdx.x*2720+threadIdx.x;
	 if (blockIdx.x <2){
	    for ( int i =0;i<17;i++)
	       C[index+i*160]=Cs[threadIdx.x+i*160];
	 }
	 else{
	    for ( int i =0;i<16;i++)
	       C[index+i*160]=Cs[threadIdx.x+i*160];
	 }
	 __syncthreads();
	 //	__threadfence();


	 start_global_barrier(1,prev_dimi);

	 //	if(threadIdx.x==0)
	 //	printf("blockIdx.x=%d count=%d ",blockIdx.x,count);
	 //      __syncthreads();
	 for (int i=0;i<3;i++){
	    index=threadIdx.x + i*160;
	    if (index < 400){
	       if (trank[mu*3+1]== 20)
		  Bs[index]=B[index+mu*1200+400];
	       else
		  Bs[index]=BU[index+mu*1200+400];
	    }
	 }
	 __syncthreads();

	 index= blockIdx.x*136+threadIdx.x;
	 if (threadIdx.x<num_threads)        {
#pragma unroll 
	    for ( int i =0;i<20;i++)
	       AL1[i]=C[index+i*400];
	 }__syncthreads();
	 if (threadIdx.x<num_threads)        {


#pragma unroll 
	    for (int  j=0;j<20;j++)//number of cols of B
	    {
	       double Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<20;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
	       Cs[threadIdx.x*20+j]=Ctemp;	
	       //       __syncthreads();

	    }
	 }
	 __syncthreads();
	 index= blockIdx.x*2720+threadIdx.x;
	 if (blockIdx.x <2){
	    for ( int i =0;i<17;i++)
	       C[index+i*160]=Cs[threadIdx.x+i*160];
	 }
	 else{
	    for ( int i =0;i<16;i++)
	       C[index+i*160]=Cs[threadIdx.x+i*160];
	 }
	 //	__threadfence();
	 __syncthreads();

	 start_global_barrier(2,prev_dimi);
	 //	if(threadIdx.x==0)
	 //	printf("blockIdx.x=%d count=%d ",blockIdx.x,count);
	 //      __syncthreads();
	 //
	 for (int i=0;i<3;i++){
	    index=threadIdx.x + i*160;
	    if (index < 400){
	       if (trank[mu*3+2]== 20)
		  Bs[index]=B[index+mu*1200+800];
	       else
		  Bs[index]=BU[index+mu*1200+800];
	    }
	 }__syncthreads();

	 index= blockIdx.x*136+threadIdx.x;
	 if (threadIdx.x<num_threads)        {
#pragma unroll 
	    for ( int i =0;i<20;i++)
	       AL1[i]=C[index+i*400];
	 } __syncthreads();
	 if (threadIdx.x<num_threads)        {

#pragma unroll 
	    for (int  j=0;j<20;j++)//number of cols of B
	    {
	       double Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<20;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
	       Cs[threadIdx.x*20+j]=Ctemp;
	       //        __syncthreads();

	    }
	 }
	 __syncthreads();
	 if (doitS[mu]){	
	    index= blockIdx.x*2720+threadIdx.x;
	    if (blockIdx.x <2){
	       for ( int i =0;i<17;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    else{
	       for ( int i =0;i<16;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    //__threadfence();
	    __syncthreads();

	    start_global_barrier(3,prev_dimi);


	    if (threadIdx.x ==0){
	       if(count[prev_dimi] ==9)
		  count[prev_dimi]=0;	
	    }
	    __syncthreads();
	 }
	 else{

	    if (threadIdx.x ==0){
	       if(count[prev_dimi] ==6)
		  count[prev_dimi]=0;	
	    }
	    __syncthreads();


	 }

	 if (doitS[mu]){	

	if ( doitS[mu] & 1){
	    for (int i=0;i<3;i++){
	       index=threadIdx.x + i*160;
	       if (index <400)
		  Bs[index]=BVT[index+mu*1200];
	    }
	    __syncthreads();

	    index= blockIdx.x*136+threadIdx.x;
	    if (threadIdx.x<num_threads){
#pragma unroll 
	       for ( int i =0;i<20;i++)
		  AL1[i]=C[index+i*400];
	    }__syncthreads();



	    if (threadIdx.x<num_threads){
#pragma unroll 
	       for (int  j=0;j<20;j++)//number of cols of B
	       {
		  double Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<20;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
		  Cs[threadIdx.x*20+j]=Ctemp;

	       }
	    }
	    __syncthreads();
	    index= blockIdx.x*2720+threadIdx.x;
	    if (blockIdx.x <2){
	       for ( int i =0;i<17;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    else{
	       for ( int i =0;i<16;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    __syncthreads();
	}
	else{
	index= blockIdx.x*136+threadIdx.x;
            if (threadIdx.x<num_threads){
#pragma unroll 
               for ( int i =0;i<20;i++)
                  AL1[i]=C[index+i*400];
            }__syncthreads();

	if (threadIdx.x<num_threads){
#pragma unroll 
        for ( int i =0;i<20;i++)
		Cs[threadIdx.x*20+i] = AL1[i];

	}
        __syncthreads();
	
	index= blockIdx.x*2720+threadIdx.x;
            if (blockIdx.x <2){
#pragma unroll 
               for ( int i =0;i<17;i++)
                  C[index+i*160]=Cs[threadIdx.x+i*160];
            }
            else{
#pragma unroll 
               for ( int i =0;i<16;i++)
                  C[index+i*160]=Cs[threadIdx.x+i*160];
            }
            __syncthreads();


	}
	    start_global_barrier(1,prev_dimi);

	if ( doitS[mu] & 2){
	    for (int i=0;i<3;i++){
	       index=threadIdx.x + i*160;
	       if (index <400)
		  Bs[index]=BVT[index+mu*1200+400];
	    }
	    __syncthreads();

	    if (threadIdx.x<num_threads)
	    {
	       index= blockIdx.x*136+threadIdx.x;
#pragma unroll 
	       for ( int i =0;i<20;i++)
		  AL1[i]=C[index+i*400];
	    }__syncthreads();


	    if (threadIdx.x<num_threads){
#pragma unroll 
	       for (int  j=0;j<20;j++)//number of cols of B
	       {
		  double Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<20;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
		  Cs[threadIdx.x*20+j]=Ctemp;	

	       }
	    }
	    __syncthreads();
	    index= blockIdx.x*2720+threadIdx.x;
	    if (blockIdx.x <2){
	       for ( int i =0;i<17;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    else{
	       for ( int i =0;i<16;i++)
		  C[index+i*160]=Cs[threadIdx.x+i*160];
	    }
	    __syncthreads();
	}
	else{
	index= blockIdx.x*136+threadIdx.x;
            if (threadIdx.x<num_threads){
#pragma unroll 
               for ( int i =0;i<20;i++)
                  AL1[i]=C[index+i*400];
            }__syncthreads();

	if (threadIdx.x<num_threads){
#pragma unroll 
        for ( int i =0;i<20;i++)
		Cs[threadIdx.x*20+i] = AL1[i];

	}
        __syncthreads();
	
	index= blockIdx.x*2720+threadIdx.x;
            if (blockIdx.x <2){
#pragma unroll 
               for ( int i =0;i<17;i++)
                  C[index+i*160]=Cs[threadIdx.x+i*160];
            }
            else{
#pragma unroll 
               for ( int i =0;i<16;i++)
                  C[index+i*160]=Cs[threadIdx.x+i*160];
            }
            __syncthreads();


	}
	    start_global_barrier(2,prev_dimi);



	if ( doitS[mu] & 4){
	    for (int i=0;i<3;i++){
	       index=threadIdx.x + i*160;
	       if (index <400)
		  Bs[index]=BVT[index+mu*1200+800];
	    }__syncthreads();

	    if (threadIdx.x<num_threads)
	    {
	       index= blockIdx.x*136+threadIdx.x;
#pragma unroll 
	       for ( int i =0;i<20;i++)
		  AL1[i]=C[index+i*400];
	    }__syncthreads();

	    if (threadIdx.x<num_threads){
#pragma unroll 
	       for (int  j=0;j<20;j++)//number of cols of B
	       {
		  double Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<20;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*20];//4 way bank conflict
		  Cs[threadIdx.x*20+j]=Ctemp;

	       }
	    }
	    __syncthreads();
	}
	else{
	index= blockIdx.x*136+threadIdx.x;
            if (threadIdx.x<num_threads){
#pragma unroll 
               for ( int i =0;i<20;i++)
                  AL1[i]=C[index+i*400];
            }__syncthreads();

	if (threadIdx.x<num_threads){
#pragma unroll 
        for ( int i =0;i<20;i++)
		Cs[threadIdx.x*20+i] = AL1[i];

	}
        __syncthreads();
	}
	    /*if (blockIdx.x <2){
	      for ( int i =0;i<17;i++)
	      C[blockIdx.x*136*ldb+threadIdx.x+i*160]=Cs[threadIdx.x+i*160];
	      }
	      else{
	      for ( int i =0;i<16;i++)
	      C[blockIdx.x*136*ldb+threadIdx.x+i*160]=Cs[threadIdx.x+i*160];
	      }
	      __syncthreads();

	      start_global_barrier(3,prev_dimi);
	     */
	    if (threadIdx.x ==0){
	       if(count[prev_dimi] ==6)
		  count[prev_dimi]=0;	
	    }
	    __syncthreads();
	 }

	 if (threadIdx.x<num_threads){
#pragma unroll 
	    for (int i=0;i<20;i++)
	       resultS[i]= resultS[i]+Cs[threadIdx.x*20+i]*mufacS[mu];
	 }
	 __syncthreads();
      }
   }
   if (threadIdx.x<num_threads){
      index= blockIdx.x*136+threadIdx.x;
#pragma unroll 
      for (int i=0;i<20;i++)
	 //result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
	 result[(index)*20+i]= resultS[i];
   }
   __syncthreads();

}


__global__ void cu_mtxmq_integral_113(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, const double* BVT, int8_t *doit, double *mufac, double *result, int rank, long *transr, const double *BU)
{
   extern __shared__ double Bs[];
   __shared__ double Cs[1024];
   double AL1[10];
   double Ctemp;
   __shared__ double resultS[1024];
   __shared__ int8_t  doitS[128];
   __shared__ double mufacS[128];
   __shared__ double As[1024];
   __shared__ int trank[512];
   if (threadIdx.x <rank)
      doitS[threadIdx.x]=doit[threadIdx.x];
   __syncthreads();
   if (threadIdx.x <rank)
      mufacS[threadIdx.x]=mufac[threadIdx.x];
   __syncthreads();


#pragma unroll 
   for (int i=0;i<4;i++)	{
      int index =threadIdx.x+i*128;
      if ( index <rank*3)
	 trank[index]=transr[index];
   }
   __syncthreads();
   //  if (threadIdx.x<100){
#pragma unroll 
   for ( int i =0;i<8;i++)
      As[threadIdx.x+i*128]=A[threadIdx.x+i*128];
   //	}
   __syncthreads();
   //	if (threadIdx.x<lda){
#pragma unroll 
   for (int i=0;i<8;i++)
      //resultS[threadIdx.x+lda*i]=result[threadIdx.x+lda*i];
      resultS[threadIdx.x+128*i]=0.0;
   //	}
   __syncthreads();


   for (int mu=0;mu<rank;mu++){
      if (mufacS[mu]){

	 //	if (threadIdx.x==0)
	 //printf("rank=%d \n",trank[mu]);	
	 //      __syncthreads();
	 if (trank[mu*3]== ldb){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=B[threadIdx.x+mu*300];
	 }
	 else{
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=BU[threadIdx.x+mu*300];
	 }
	 __syncthreads();



	 if (threadIdx.x<lda){
	    for ( int i =0;i<10;i++)
	       AL1[i]=As[threadIdx.x+i*lda];
	 }__syncthreads();
	 // lda replaced by dim/trans.r
	 //ldb replaced by trans.r
	 // 10 in the loop replaced by trans.r 

	 if (threadIdx.x<lda){

#pragma unroll 
	    for (int  j=0;j<ldb;j++)//number of cols of B
	    {
	       Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<10;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
	       Cs[threadIdx.x*ldb+j]=Ctemp;
	    }
	 }
	 __syncthreads();


	 if (trank[mu*3+1]== ldb){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(B+threadIdx.x+mu*300+100);
	 }
	 else{
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(BU+threadIdx.x+mu*300+100);
	 }
	 __syncthreads();

	 if (threadIdx.x<lda)        {
#pragma unroll 
	    for ( int i =0;i<10;i++)
	       AL1[i]=Cs[threadIdx.x+i*lda];
	 }__syncthreads();

	 if (threadIdx.x<lda){

#pragma unroll 
	    for (int  j=0;j<ldb;j++)//number of cols of B
	    {
	       Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<10;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
	       Cs[threadIdx.x*ldb+j]=Ctemp;	
	    }

	 }
	 __syncthreads();


	 if (trank[mu*3+2]== ldb){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(B+threadIdx.x+mu*300+200);
	 }
	 else{
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(BU+threadIdx.x+mu*300+200);
	 }
	 __syncthreads();

	 if (threadIdx.x<lda)        {
#pragma unroll 
	    for ( int i =0;i<10;i++)
	       AL1[i]=Cs[threadIdx.x+i*lda];
	 }__syncthreads();

	 if (threadIdx.x<lda){

#pragma unroll 
	    for (int  j=0;j<ldb;j++)//number of cols of B
	    {
	       Ctemp=0.0;
#pragma unroll 
	       for (int  i =0;i<10;i++) // number of rows of B
		  Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
	       Cs[threadIdx.x*ldb+j]=Ctemp;
	    }

	 }
	 __syncthreads();


	 if (doitS[mu]){	
	
//	 	if (threadIdx.x==0)
//	 printf("doit10[%d]=%d \n",mu,doitS[mu]);	
//	       __syncthreads();
	if ( doitS[mu] & 1){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=BVT[threadIdx.x+mu*300];
	    __syncthreads();

	    if (threadIdx.x<lda){
#pragma unroll 
	       for ( int i =0;i<10;i++)
		  AL1[i]=Cs[threadIdx.x+i*lda];
	    }__syncthreads();



	    if (threadIdx.x<lda){
#pragma unroll 
	       for (int  j=0;j<ldb;j++)//number of cols of B
	       {
		  Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<10;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
		  Cs[threadIdx.x*ldb+j]=Ctemp;
	       }

	    }
	    __syncthreads();
	}
	else{

	if (threadIdx.x<100)        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*100];
        }__syncthreads();
	
	if (threadIdx.x<100){
	#pragma unroll 
        for ( int i =0;i<10;i++)
		Cs[threadIdx.x*10+i] = AL1[i];

	}
        __syncthreads();


	}

	if ( doitS[mu] & 2){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(BVT+threadIdx.x+mu*300+100);
	    __syncthreads();

	    if (threadIdx.x<lda)
	    {
#pragma unroll 
	       for ( int i =0;i<10;i++)
		  AL1[i]=Cs[threadIdx.x+i*lda];
	    }__syncthreads();


	    if (threadIdx.x<lda){
#pragma unroll 
	       for (int  j=0;j<ldb;j++)//number of cols of B
	       {
		  Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<10;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
		  Cs[threadIdx.x*ldb+j]=Ctemp;	
	       }

	    }
	    __syncthreads();
	}
	else{

	if (threadIdx.x<100)        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*100];
        }__syncthreads();
	
	if (threadIdx.x<100){
	#pragma unroll 
        for ( int i =0;i<10;i++)
		Cs[threadIdx.x*10+i] = AL1[i];

	}
        __syncthreads();


	}



	if ( doitS[mu] & 4){
	    if (threadIdx.x<100)
	       Bs[threadIdx.x]=*(BVT+threadIdx.x+mu*300+200);
	    __syncthreads();

	    if (threadIdx.x<lda)
	    {
#pragma unroll 
	       for ( int i =0;i<10;i++)
		  AL1[i]=Cs[threadIdx.x+i*lda];
	    }__syncthreads();

	    if (threadIdx.x<lda){
#pragma unroll 
	       for (int  j=0;j<ldb;j++)//number of cols of B
	       {
		  Ctemp=0.0;
#pragma unroll 
		  for (int  i =0;i<10;i++) // number of rows of B
		     Ctemp +=AL1[i]*Bs[j+i*ldb];//4 way bank conflict
		  Cs[threadIdx.x*ldb+j]=Ctemp;

	       }
	    }
	    __syncthreads();
	}
	else{

	if (threadIdx.x<100)        {
        #pragma unroll 
        for ( int i =0;i<10;i++)
                AL1[i]=Cs[threadIdx.x+i*100];
        }__syncthreads();
	
	if (threadIdx.x<100){
	#pragma unroll 
        for ( int i =0;i<10;i++)
		Cs[threadIdx.x*10+i] = AL1[i];

	}
        __syncthreads();


	}
	 }

	 //  if (threadIdx.x<lda){
#pragma unroll 
	 for (int i=0;i<8;i++)
	    resultS[threadIdx.x+i*128]= resultS[threadIdx.x+i*128]-Cs[threadIdx.x+i*128]*mufacS[mu];
	 __syncthreads();
	 //}
	 //      __syncthreads();
      }
   }
   if (threadIdx.x<lda){
#pragma unroll 
      for (int i=0;i<10;i++)
	 //result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
	 result[threadIdx.x+i*lda]= resultS[threadIdx.x+i*lda];
   }
   __syncthreads();

}
