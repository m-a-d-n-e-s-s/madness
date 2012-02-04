__global__ void cu_axpy_20(long n , double *a, double *b, double alpha) {

int index;
double x,y;
 
index =blockIdx.x*blockDim.x + threadIdx.x;
 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();

index = 2048+index;
 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();

index = 2048 + index;
 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();

 index =2048+index;
if (index < n){
 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();
}

}

__global__ void cu_axpy_10(long n , double *a, double *b, double alpha) {

int index;
double x,y;

 index = blockIdx.x*blockDim.x+threadIdx.x;

 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();
 

index = 512 + index;
if (index<n){
 x = a[index];
 y= b[index];
__syncthreads();

x=alpha*y+x;
a[index]=x;
__syncthreads();
}
}
