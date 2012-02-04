__global__ void dgemm_madness_batched_12(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[12];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<12;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<2;i++)
        	Bs[blockIdx.y*ldb*ldb+i*blockDim.x+threadIdx.x]=B[blockIdx.y*ldb*ldb+i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<12;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<12;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_14(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[14];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<14;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<2;i++)
        	Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<14;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<14;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_16(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[16];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<16;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<2;i++)
        	Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<16;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<16;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_18(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[18];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<18;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<3;i++)
        	Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<18;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<18;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_20(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[20];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<20;i++)
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
		        for (int  i =0;i<20;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<20;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_22(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[22];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
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
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<22;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_24(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[24];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<24;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<5;i++)
        	Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<24;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<24;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

__global__ void dgemm_madness_batched_26(double *A, int lda, double *B, int ldb, double* C, int ldc, int k, double alpha, double beta)
{
        extern __shared__ double Bs[];
        register double AL1[26];
        register double Cs;
        
        int x=blockIdx.y*lda*ldb+blockIdx.x*blockDim.x+threadIdx.x;
        
	//fetch part of A from global memory into registers
	#pragma unroll 
        for ( int i =0;i<26;i++)
                AL1[i]=A[x+i*lda];
        __syncthreads();
        
	#pragma unroll 
        for (int i =0;i<6;i++)
        	Bs[i*blockDim.x+threadIdx.x]=B[i*blockDim.x+ threadIdx.x];
        __syncthreads();

        if (blockIdx.x < gridDim.x-1){
        
		#pragma unroll 
	        for (int  j=0;j<ldb;j++)//number of cols of B
        	{
        		Cs=0.0;
	 		#pragma unroll 
		        for (int  i =0;i<26;i++) // number of rows of B
                		Cs +=AL1[i]*Bs[j+i*ldb];
                 
			C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                	__syncthreads();

        	}
        }

        else{
	        if (threadIdx.x<(lda - blockDim.x*blockIdx.x)){
        		#pragma unroll 
		        for (int  j=0;j<ldb;j++)//number of cols of B
       			{
        			Cs=0.0;
       				#pragma unroll 
			        for (int  i =0;i<26;i++) // number of rows of B
                			Cs +=AL1[i]*Bs[j+i*ldb];
		                 C[blockIdx.y*lda*ldb+blockIdx.x*blockDim.x*ldb+threadIdx.x*ldb+j]=Cs;
                		__syncthreads();

        		}
        	}

        }
}

