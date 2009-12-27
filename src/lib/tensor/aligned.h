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
*/
#ifndef MADNESS_TENSOR_ALIGNED_H__INCLUDED
#define MADNESS_TENSOR_ALIGNED_H__INCLUDED

#include <madness_config.h>
#include <tensor/tensor.h>

namespace madness {
    template <typename T>
    void aligned_zero(long n, T* a) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4) {
            a[0] = 0;
            a[1] = 0;
            a[2] = 0;
            a[3] = 0;
        }
        for (long i=0; i<rem; i++) *a++ = 0;
    }

    template <typename T, typename Q>
    void aligned_axpy(long n, T* restrict a, const T* restrict b, Q s) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4,b+=4) {
            a[0] += s*b[0];
            a[1] += s*b[1];
            a[2] += s*b[2];
            a[3] += s*b[3];
        }
        for (long i=0; i<rem; i++) *a++ += s * *b++;
    }



#if (defined(X86_32) || defined(X86_64))
    template <>
    void aligned_zero<double>(long n, double* a);

    template <>
    void aligned_zero<double_complex>(long n, double_complex* a);

//  template <>
//  void aligned_axpy(long n, double* restrict a, const double* restrict b, double s);
#endif

    void aligned_add(long n, double* restrict a, const double* restrict b);
    void aligned_sub(long n, double* restrict a, const double* restrict b);
    void aligned_add(long n, double_complex* restrict a, const double_complex* restrict b);
    void aligned_sub(long n, double_complex* restrict a, const double_complex* restrict b);

}

#endif // MADNESS_TENSOR_ALIGNED_H__INCLUDED
