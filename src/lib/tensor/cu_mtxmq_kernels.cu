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
#ifndef MADNESS_TENSOR_CU_MTXMQ_KERNEL__INCLUDED
#define MADNESS_TENSOR_CU_MTXMQ_KERNEL__INCLUDED


//#include <madness_config.h>
#include <stdio.h>
#include <fortran_ctypes.h>
#include <iostream>
#include <cuda.h>
#include <cublas.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
//#include <tensor/cu_mtxmq.h>
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#define BLOCK_SIZE 16 
#define TILE_DIM    16
#define BLOCK_ROWS  16
// Coalesced transpose with no bank conflicts
template <typename aT>
__global__ void transposeNoBankConflicts(aT *odata, const aT *idata, int width, int height, int nreps)
{
   int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;

  int index_in  = xIndex + width * yIndex;
  int index_out = yIndex + height * xIndex;
  for (int r=0; r < nreps; r++) {
    for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
      odata[index_out+i] = idata[index_in+i*width];
    }
  }
}

template <typename aT, typename bT, typename cT>
    __global__ void
    matrixMul_coalescing_rem( cT* C, aT* A, bT* B, int wA, int wB, int hA)
    {
      //printf("cuda kernel execution");
      // Block index
  int bx = blockIdx.x;
  int by = blockIdx.y;

  // Thread index
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  // Accumulate row i of A and column j of B
  int i = by * blockDim.y + ty + (hA/BLOCK_SIZE)*BLOCK_SIZE;
  int j = bx * blockDim.x + tx + (hA/BLOCK_SIZE)*BLOCK_SIZE;

  cT accu;

  for(int k=0; k<wA; k++){
    accu = accu + A[ i * wA + k ] * B[ k * wB + j ];
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  C[ i * wB + j ] = accu;
    }
    
    template <typename aT, typename bT, typename cT>
    __global__ void
    matrixMul_coalescing( cT* C, aT* A, bT* B, int wA, int wB)
    {
  // Block index
  int bx = blockIdx.x;
  int by = blockIdx.y;

  // Thread index
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  // Accumulate row i of A and column j of B
  int i = by * blockDim.y + ty;
  int j = bx * blockDim.x + tx;

  cT accu;

  for(int k=0; k<wA; k++){
    accu = accu + A[ i * wA + k ] * B[ k * wB + j ];
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  C[ i * wB + j ] = accu;

    }
#endif // MADNESS_TENSOR_CU_MTXMQ_KERNEL__INCLUDED

