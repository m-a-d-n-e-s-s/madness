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

// extern  void cu_mTxmq(int m, int n,int k, double *C, const double *A, const double *B);

#ifdef __cplusplus

template <typename aT, typename bT, typename cT>
    void cu_mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c, const aT* a, const bT* b);
    
#endif

    /*
   extern  void cu_mTxmq(long dimi, long dimj, long dimk,
               float* restrict c, const float a, const float* b) ;
     
    extern  void cu_mTxmq(long m, long n,long k, double *C, const double *A, const double *B);
    
     extern  void cu_mTxmq(long m, long n,long k, std::complex<double> *C, const std::complex<double> *A, const std::complex<double> *B);
	
 
*/

#endif // MADNESS_TENSOR_CU_MTXMQ_H__INCLUDED

