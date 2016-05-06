
#include <iostream>

// CUDA runtime
#include <cuda_runtime.h>

// helper functions and utilities to work with CUDA
//#include <helper_functions.h>
//#include <helper_cuda.h>

#include <cstdio>

__global__ void hello_world() {
  printf("CUDA thread [%d, %d] says \"hello, world!\"\n",\
            blockIdx.y*gridDim.x+blockIdx.x,\
            threadIdx.z*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x);
}

void __cuda_hello_world() {
  int device_count;
  cudaGetDeviceCount(&device_count);
  if (device_count > 0) {
    std::cout << "in __cuda_hello_world: device_count=" << device_count;
    if (device_count == 1) {
  	  cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, 0);
      std::cout << " device_name=\"" << deviceProp.name << "\"";
    }
    std::cout << std::endl;
    hello_world<<<1,1>>>();
    cudaDeviceSynchronize();
  }
}
