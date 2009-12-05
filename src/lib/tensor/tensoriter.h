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


#ifndef MADNESS_TENSOR_TENSORITER_H__INCLUDED
#define MADNESS_TENSOR_TENSORITER_H__INCLUDED

/// \file tensoriter.h
/// \brief Declares TensorIterator

namespace madness {

    template <class T> class Tensor;
    template <class T> class SliceTensor;

    static const long default_jdim = 5551212; // Never a valid dimension num.

    /// Optimized iterator for tensors supporting unary, binary and ternary operations.

    template <class T, class Q = T, class R = T> class TensorIterator {
        T* _p0_save;
        Q* _p1_save;
        R* _p2_save;
    public:
        T* _p0;
        Q* _p1;
        R* _p2;
        long ndim;
        long dimj;
        long _s0;
        long _s1;
        long _s2;
        long dim[TENSOR_MAXDIM];
        long ind[TENSOR_MAXDIM];
        long stride0[TENSOR_MAXDIM];
        long stride1[TENSOR_MAXDIM];
        long stride2[TENSOR_MAXDIM];

        TensorIterator(const Tensor<T>* t0, const Tensor<Q>* t1=0, const Tensor<R>* t2=0,
                       long iterlevel=0,
                       bool optimize=true, bool fusedim=true,
                       long jdim=default_jdim);

        TensorIterator<T,Q,R>& operator++();

        inline bool operator == (const TensorIterator<T,Q,R>& a) const {
            return _p0==a._p0;
        }

        inline bool operator != (const TensorIterator<T,Q,R>& a) const {
            return _p0!=a._p0;
        }

        inline TensorIterator<T,Q,R>* operator->() {
            return this;
        };

        inline T& operator*() const {
            return *_p0;
        };

        void reset();

        void reuse(const Tensor<T>* t0, const Tensor<Q>* t1=0, const Tensor<R>* t2=0);
    };

}
#endif // MADNESS_TENSOR_TENSORITER_H__INCLUDED
