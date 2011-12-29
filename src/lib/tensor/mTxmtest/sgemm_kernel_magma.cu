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
__global__ void cu_mtxmq_integral_10TripleWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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
__global__ void cu_mtxmq_integral_10FivetimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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
__global__ void cu_mtxmq_integral_10NinetimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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
__global__ void cu_mtxmq_integral_20PerBlock100times(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta, int offset)
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

//One kernel does two matmuls
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

//One kernel does three matmuls
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

//One kernel does the same task twice, sequentially
__global__ void cu_mtxmq_integral_20TwoSequentialTask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//One kernel does the same task three times, sequentially
__global__ void cu_mtxmq_integral_20ThreeSequentialTask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//One kernel does the same task ten times, sequentially
__global__ void cu_mtxmq_integral_20TenSeqquentialTask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//One kernel does the same task hundred times, sequentially
__global__ void cu_mtxmq_integral_20HundredSequentialTask(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//The work and access are done 3 times
__global__ void cu_mtxmq_integral_20TripleWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//The work and access are done 5 times
__global__ void cu_mtxmq_integral_20FivetimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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
//The work and access are done nine times
__global__ void cu_mtxmq_integral_20NinetimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
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

//The work and access are done hundred times
__global__ void cu_mtxmq_integral_20HundredtimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
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
	       
        Bs[threadIdx.x]=B[threadIdx.x];
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
		        C[x*ldb+j]=Cs;  // single
                                  
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
                		 C[x*ldb+j]=Cs;
		                __syncthreads();

		        }
		}

	}	
}

//The matmul is all one block, trying to access all data from registers/shared mem
//Data is loaded from input A 10% by 10%, so that everything fits in registers
//Only writes to C are not local
__global__ void cu_mtxmq_integral_20OneblockLocality(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[threadIdx.x]=B[threadIdx.x];
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

			 C[x*ldb+j]=Cs;
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
			 C[x*ldb+j] += Cs;
		}
	  }
        }
}

//The matmul is all one block, trying to access all data from registers/shared mem
//Data is loaded from input A 10% by 10%, so that everything fits in registers
//work is done 100 times, to see the limits of the kernel
//Only writes to C are not local
__global__ void cu_mtxmq_integral_20_OneblockLocality_HundredtimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[threadIdx.x]=B[threadIdx.x];
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

			 C[x*ldb+j]=Cs;

		}
	}

        for (int ik = 2; ik < 20; ik += 2){
          #pragma unroll 
          for ( int i =0;i<2;i++){
		AL1[i]=A[x+(ik+i)*lda];
                /*
                */
          }
          
          for (int jk = 0; jk < 20; jk += 2){

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

//The matmul is all one block, trying to access all data from registers/shared mem
//Data is loaded from input A 10% by 10%, so that everything fits in registers
//work is done 100 times, to see the limits of the kernel
//C accesses are replaced by writing a constant value, to emphasize the rest of the computation
__global__ void cu_mtxmq_integral_20OneblockLocality_HundredtimesWork_nowrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

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

	}

        for (int ik = 2; ik < 20; ik += 2){
          #pragma unroll 
          for ( int i =0;i<2;i++){
		AL1[i]=A[x+(ik+i)*lda];
                /*
                */
          }
          
          for (int jk = 0; jk < 20; jk += 2){

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

		}
        
          }
        }

        for (int jk = 0; jk < 20; jk++){
            C[x*ldb + jk] = 2.325263;
        }
}

__global__ void cu_mtxmq_integral_20OneblockLocality(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j++)//number of cols of B
	{
		Cs=0.0;
		#pragma unroll 
		for (int  i =0;i<2;i++){ // number of rows of B
			Cs +=AL1[i]*Bs[j+i*prev_dimi];
		}

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
		}
                C[x*ldb+j] = Cs;
	}

}

__global__ void cu_mtxmq_integral_20OneblockLocality_WriteSame(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

        #pragma unroll 
        for ( int i =0;i<2;i++){
		AL1[i]=A[x+i*lda];
                /*
                */
        }
	Bs[threadIdx.x]=B[threadIdx.x];
        __syncthreads();

	#pragma unroll 
	for (int  j=0;j<ldb;j++)//number of cols of B
	{
		Cs=0.0;
		#pragma unroll 
		for (int  i =0;i<2;i++){ // number of rows of B
			Cs +=AL1[i]*Bs[j+i*prev_dimi];
		}

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

		}
	}

        //just write some value, which is incorrect, 
        //but emphasizes the rest of the computation
        for (int jk = 0; jk < 20; jk++){
            C[x*ldb + jk] = Cs;
        }
}

//this version first writes to shared memory, then to C,
//attempting to coalesce writes to C
__global__ void cu_mtxmq_integral_202optCwrite(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        extern __shared__ double Cshare[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

	Bs[threadIdx.x]=B[threadIdx.x];
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

			}
		        Cshare[4*x+jk-j] = Cs;
		}
                __syncthreads();
                //attempting to coalesce writes
                for (int  jk=j;jk<j+4;jk++){
                       C[(16*(x/16)+x%4)*ldb+ jk]=Cshare[16*(x/16)+x%4+ 4*(jk-j)];
                }
        }

}

//400x20 * 20x20 matmul all in one block
//work done 100 times
//the inputs for work are saved in registers (A) and shared memory (B)
//input A is loaded 10% by 10%, to avoid register overflow
//J on the outside
__global__ void cu_mtxmq_integral_20OneblockLocality_HundredtimesWork(const double *A, int lda, const double *B, int ldb, double* C, int ldc, int prev_dimi, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[2];
        register double Cs;//[SIZE];
        //fetch part of A from global memory into L1 cache
        int x=blockIdx.x*blockDim.x+threadIdx.x;   // single 

	Bs[threadIdx.x]=B[threadIdx.x];
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
