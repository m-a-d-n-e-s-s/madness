/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef TENSOR_LAPACK_H
#define TENSOR_LAPACK_H

#include <tensor/tensor.h>

/// \file tensor_lapack.h
/// \brief Prototypes for a partial interface from Tensor to LAPACK

namespace madness {

    template <typename T>
    void svd(const Tensor<T>& a, Tensor<T>* U,
             Tensor< typename Tensor<T>::scalar_type >* s, Tensor<T>* VT);

    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>* x);

    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>* x, Tensor< typename Tensor<T>::scalar_type >* s, long *rank);
    template <typename T>
    void syev(const Tensor<T>& A,
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e);

    template <typename T>
    void sygv(const Tensor<T>& A, const Tensor<T>& B, int itype,
              Tensor<T>* V, Tensor< typename Tensor<T>::scalar_type >* e);

    template <typename T>
    void cholesky(const Tensor<T>& A, Tensor<T>* L);

    template <typename T>
    void triangular_solve(const Tensor<T>& L, Tensor<T>& B, 
                          const char* side, const char* uplo, const char* transa);

    bool test_tensor_lapack();
}

#endif
