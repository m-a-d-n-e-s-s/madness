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

  $LastChangedDate$
  $Rev$
*/

  
#include <iostream>
#include <algorithm>
#include <complex>

#include <cmath>


/// \file tensoriter.cc
/// \brief Implements TensorIterator

// These here only to instantiate specializations
typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;

#include "tensor.h"

namespace madness {

    /// General iterator to compose operations over up to three tensors
    template <class T, class Q, class R>
    TensorIterator<T,Q,R>::TensorIterator(const Tensor<T>* t0,
                                          const Tensor<Q>* t1,
                                          const Tensor<R>* t2,
                                          long iterlevel,
                                          bool optimize,
                                          bool fusedim,
                                          long jdim) {

        /// optimize -> reorder dimensions for optimal strides (only
        /// applies if iterlevel=1).  If jdim==default_jdim, all dimensions
        /// are reordered for optimal stride.  If jdim is not the default
        /// value, then dimension jdim is excluded from the set of
        /// dimensions being optimized.
        ///
        /// fusedim -> concatenate contiguous dimensions into the inner
        /// loop (only applies if iterlevel=1 and jdim=default_jdim).
        ///
        /// iterlevel -> 0=elementwise, 1=vectorwise.  Elementwise implies
        /// that the iterator returns successive elements, and explicitly
        /// in the expected order for the tensor (i.e., with the last index
        /// varying fastest).  Vectorwise implies that the user is responsible
        /// for iterating over dimension jdim.
        ///
        /// jdim -> if iterlevel=1, then jdim determines which dimension is
        /// left for iteration with an explicit for loop.  Negative values
        /// implies jdim+=ndim following the convention of Slice (and
        /// Python).  If (optimize) the default is the fastest varying
        /// dimension.  If (!optimize) the default is the last dimension
        /// (jdim=-1). If (fusedim) contiguous dimensions are fused into
        /// this inner loop.  Specifying a non-default value for jdim
        /// disables fusedim and restricts optimization to reordering only
        /// the exterior loops (so that the loop the user is iterating over
        /// corresponds exactly to those in dimension jdim).
        ///
        /// During iteration:
        ///
        /// ind[] will contain the current iteration indices .. BUT
        /// if optimize=true, they will not necessarily be in the
        /// order of those of the tensor.  if fusedim is true, then
        /// there may be fewer dimensions than the input tensor.
        ///
        /// _p0, _p1, _p2 will point to the current elements of t0,t1,t2
        /// noting that if iterlevel>0, the user must provide additional
        /// \c for loops to iterate over the additional dimensions
        ///
        /// stride0[], stride1[], stride2[] will contain the strides for each
        /// (possibly reordered) dimension.  If iterlevel=1, _s0, _s1, s2
        /// contain the stride info for the dimension that the user is
        /// responsible for iterating over.
        ///
        /// dimj -> the size of the j'th dimension


        if (!t0) {
            // Used to indicate end of iteration.
            _p0 = 0;
            return;
        }

        //std::printf("t0=%p t1=%p t2=%p optimize=%d fusedim=%d iterlevel=%ld jdim=%ld\n",
        //t0,t1,t2,optimize,fusedim,iterlevel,jdim);

        TENSOR_ASSERT(iterlevel==0 || iterlevel==1,"invalid iteration level",iterlevel,t0);

        // First copy basic info over
        ndim = t0->ndim;
        _p0_save = _p0 = t0->ptr();
        for (int i=0; i<ndim; i++) {
            dim[i] = t0->dim[i];
            stride0[i] = t0->stride[i];
        }
        if (t1) {
            TENSOR_ASSERT(t0->conforms(*t1),"first and second tensors do not conform",
                          0, t0);
            _p1_save = _p1 = t1->ptr();
            for (int i=0; i<ndim; i++) stride1[i] = t1->stride[i];
        } else {
            _p1_save = _p1 = 0;
        }
        if (t2) {
            TENSOR_ASSERT(t0->conforms(*t2),"first and third tensors do not conform",
                          0, t0);
            _p2_save = _p2 = t2->ptr();
            for (int i=0; i<ndim; i++) stride2[i] = t2->stride[i];
        } else {
            _p2_save = _p2 = 0;
        }

        if (iterlevel == 0) {
            // Iteration will include all dimensions
            fusedim = false;
            jdim = ndim;
            dimj = 0;
            _s0 = 0;
            _s1 = 0;
            _s2 = 0;
        } else if (iterlevel == 1) {
            // Apply -ve indexing convention for dimensions
            if (jdim < 0) jdim += ndim;

            // If permissible optimize the order of dimensions excluding
            // any non-default value of jdim.
            if (optimize) {
                for (int i=0; i<ndim; i++) {
                    if (i != jdim) {
                        for (int j=i; j<ndim; j++) {
                            if (j != jdim) {
                                if (std::abs(stride0[i]) < std::abs(stride0[j])) {
                                    std::swap(stride0[i],stride0[j]);
                                    if (t1) std::swap(stride1[i],stride1[j]);
                                    if (t2) std::swap(stride2[i],stride2[j]);
                                    std::swap(dim[i],dim[j]);
                                }
                            }
                        }
                    }
                    //std::cout << "stride0[" << i << "]=" << stride0[i] << std::endl;
                }
            }

            // Iterations will exclude dimension jdim, default is last one
            if (jdim == default_jdim) {
                jdim = ndim-1;
            } else {
                fusedim = false;
            }
            TENSOR_ASSERT(jdim>=0 && jdim < ndim, "invalid index for external iteration",
                          jdim, t0);
            ndim--;

            // Stride and dimension info for the excluded dimension
            _s0 = stride0[jdim];
            if (t1) {
                _s1 = stride1[jdim];
            } else {
                _s1 = 0;
            }
            if (t2) {
                _s2 = stride2[jdim];
            } else {
                _s2 = 0;
            }
            dimj = dim[jdim];

            // Collapse stride and dimension info for remaining dimensions
            for (int i=jdim+1; i<=ndim; i++) {
                dim[i-1] = dim[i];
                stride0[i-1] = stride0[i];
            }
            if (t1) {
                for (int i=jdim+1; i<=ndim; i++) stride1[i-1] = stride1[i];
            }
            if (t2) {
                for (int i=jdim+1; i<=ndim; i++) stride2[i-1] = stride2[i];
            }

            if (fusedim) {		// Only if jdim=default_jdim && iterlevel=1
                if (t2) {
                    for (int i=ndim-1; i>=0; i--) {
                        if (_s0*dimj == stride0[i] &&
                                _s1*dimj == stride1[i] &&
                                _s2*dimj == stride2[i]) {
                            dimj *= dim[i];
                            ndim--;
                        } else {
                            break;
                        }
                    }
                } else if (t1) {
                    for (int i=ndim-1; i>=0; i--) {
                        if (_s0*dimj == stride0[i] &&
                                _s1*dimj == stride1[i]) {
                            dimj *= dim[i];
                            ndim--;
                        } else {
                            break;
                        }
                    }
                } else {
                    for (int i=ndim-1; i>=0; i--) {
                        if (_s0*dimj == stride0[i]) {
                            dimj *= dim[i];
                            ndim--;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        // Initialize indices for the counter Use TENSOR_MAXDIM so reference to
        // optimized-away dimensions is vaguely meaningful.
        for (int i=0; i<TENSOR_MAXDIM; i++) ind[i] = 0;

        //   std::printf("ndim=%ld dimj=%ld _s0=%ld _s1=%ld _s2=%ld _p0=%p _p1=%p _p2=%p\n",
        // 	      ndim,dimj,_s0,_s1,_s2,_p0,_p1,_p2);
        //   for (int i=0; i<ndim; i++) {
        //     std::printf("   %d dim=%ld stride0=%ld stride1=%ld stride2=%ld\n",
        // 		i,dim[i],stride0[i],stride1[i],stride2[i]);
        //   }
    }

    template <class T, class Q, class R>
    TensorIterator<T,Q,R>& TensorIterator<T,Q,R>::operator++() {
        long d = ndim-1;
        if (d<0 || _p0==0) {
            _p0 = 0;
            return *this;
        }
        while (ind[d] >= (dim[d] - 1)) {
            _p0 -= ind[d] * stride0[d];
            if (_p1)	_p1 -= ind[d] * stride1[d];
            if (_p2)	_p2 -= ind[d] * stride2[d];
            ind[d] = 0;
            d--;
            if (d < 0) {
                _p0 = 0;
                return *this;
            }
        }
        _p0 += stride0[d];
        if (_p1) _p1 += stride1[d];
        if (_p2) _p2 += stride2[d];
        ind[d]++;
        return *this;
    }

    /// Reset the iterator back to the start ...
    template <class T, class Q, class R>
    void TensorIterator<T,Q,R>::reset() {
        _p0 = _p0_save;
        _p1 = _p1_save;
        _p2 = _p2_save;
        for (int i=0; i<TENSOR_MAXDIM; i++) ind[i] = 0;
    }


    /// Reuse this iterator to iterate over other Tensors

    /// The point of this method is to optimize away the construction of a
    /// TensorIterator when applying the same operation repeatedly to
    /// multiple tensors with identical shapes & strides.  We trust the
    /// caller to ensure all shapes & strides match between the new
    /// tensors and the ones used in the original constructor.
    template <class T, class Q, class R>
    void TensorIterator<T,Q,R>::reuse(const Tensor<T>* t0,
                                      const Tensor<Q> *t1,
                                      const Tensor<R> *t2) {
        _p0 = _p0_save = t0->ptr();
        if (t1) _p1 = _p1_save = t1->ptr();
        if (t2) _p2 = _p2_save = t2->ptr();
        for (int i=0; i<TENSOR_MAXDIM; i++) ind[i] = 0;
    }

#include "tensoriter_spec.h"

    /*
    template class TensorIterator<double,double,double>;
    template class TensorIterator<double,double,double_complex>;
    template class TensorIterator<double,double,float>;
    template class TensorIterator<double,double,float_complex>;
    template class TensorIterator<double,double,long>;
    template class TensorIterator<double,double_complex,double>;
    template class TensorIterator<double,double_complex,double_complex>;
    template class TensorIterator<double,double_complex,float>;
    template class TensorIterator<double,double_complex,float_complex>;
    template class TensorIterator<double,double_complex,long>;
    template class TensorIterator<double,float,double>;
    template class TensorIterator<double,float,double_complex>;
    template class TensorIterator<double,float,float>;
    template class TensorIterator<double,float,float_complex>;
    template class TensorIterator<double,float,long>;
    template class TensorIterator<double,float_complex,double>;
    template class TensorIterator<double,float_complex,double_complex>;
    template class TensorIterator<double,float_complex,float>;
    template class TensorIterator<double,float_complex,float_complex>;
    template class TensorIterator<double,float_complex,long>;
    template class TensorIterator<double,long,double>;
    template class TensorIterator<double,long,double_complex>;
    template class TensorIterator<double,long,float>;
    template class TensorIterator<double,long,float_complex>;
    template class TensorIterator<double,long,long>;
    template class TensorIterator<double_complex,double,double>;
    template class TensorIterator<double_complex,double,double_complex>;
    template class TensorIterator<double_complex,double,float>;
    template class TensorIterator<double_complex,double,float_complex>;
    template class TensorIterator<double_complex,double,long>;
    template class TensorIterator<double_complex,double_complex,double>;
    template class TensorIterator<double_complex,double_complex,double_complex>;
    template class TensorIterator<double_complex,double_complex,float>;
    template class TensorIterator<double_complex,double_complex,float_complex>;
    template class TensorIterator<double_complex,double_complex,long>;
    template class TensorIterator<double_complex,float,double>;
    template class TensorIterator<double_complex,float,double_complex>;
    template class TensorIterator<double_complex,float,float>;
    template class TensorIterator<double_complex,float,float_complex>;
    template class TensorIterator<double_complex,float,long>;
    template class TensorIterator<double_complex,float_complex,double>;
    template class TensorIterator<double_complex,float_complex,double_complex>;
    template class TensorIterator<double_complex,float_complex,float>;
    template class TensorIterator<double_complex,float_complex,float_complex>;
    template class TensorIterator<double_complex,float_complex,long>;
    template class TensorIterator<double_complex,long,double>;
    template class TensorIterator<double_complex,long,double_complex>;
    template class TensorIterator<double_complex,long,float>;
    template class TensorIterator<double_complex,long,float_complex>;
    template class TensorIterator<double_complex,long,long>;
    template class TensorIterator<float,double,double>;
    template class TensorIterator<float,double,double_complex>;
    template class TensorIterator<float,double,float>;
    template class TensorIterator<float,double,float_complex>;
    template class TensorIterator<float,double,long>;
    template class TensorIterator<float,double_complex,double>;
    template class TensorIterator<float,double_complex,double_complex>;
    template class TensorIterator<float,double_complex,float>;
    template class TensorIterator<float,double_complex,float_complex>;
    template class TensorIterator<float,double_complex,long>;
    template class TensorIterator<float,float,double>;
    template class TensorIterator<float,float,double_complex>;
    template class TensorIterator<float,float,float>;
    template class TensorIterator<float,float,float_complex>;
    template class TensorIterator<float,float,long>;
    template class TensorIterator<float,float_complex,double>;
    template class TensorIterator<float,float_complex,double_complex>;
    template class TensorIterator<float,float_complex,float>;
    template class TensorIterator<float,float_complex,float_complex>;
    template class TensorIterator<float,float_complex,long>;
    template class TensorIterator<float,long,double>;
    template class TensorIterator<float,long,double_complex>;
    template class TensorIterator<float,long,float>;
    template class TensorIterator<float,long,float_complex>;
    template class TensorIterator<float,long,long>;
    template class TensorIterator<float_complex,double,double>;
    template class TensorIterator<float_complex,double,double_complex>;
    template class TensorIterator<float_complex,double,float>;
    template class TensorIterator<float_complex,double,float_complex>;
    template class TensorIterator<float_complex,double,long>;
    template class TensorIterator<float_complex,double_complex,double>;
    template class TensorIterator<float_complex,double_complex,double_complex>;
    template class TensorIterator<float_complex,double_complex,float>;
    template class TensorIterator<float_complex,double_complex,float_complex>;
    template class TensorIterator<float_complex,double_complex,long>;
    template class TensorIterator<float_complex,float,double>;
    template class TensorIterator<float_complex,float,double_complex>;
    template class TensorIterator<float_complex,float,float>;
    template class TensorIterator<float_complex,float,float_complex>;
    template class TensorIterator<float_complex,float,long>;
    template class TensorIterator<float_complex,float_complex,double>;
    template class TensorIterator<float_complex,float_complex,double_complex>;
    template class TensorIterator<float_complex,float_complex,float>;
    template class TensorIterator<float_complex,float_complex,float_complex>;
    template class TensorIterator<float_complex,float_complex,long>;
    template class TensorIterator<float_complex,long,double>;
    template class TensorIterator<float_complex,long,double_complex>;
    template class TensorIterator<float_complex,long,float>;
    template class TensorIterator<float_complex,long,float_complex>;
    template class TensorIterator<float_complex,long,long>;
    template class TensorIterator<long,double,double>;
    template class TensorIterator<long,double,double_complex>;
    template class TensorIterator<long,double,float>;
    template class TensorIterator<long,double,float_complex>;
    template class TensorIterator<long,double,long>;
    template class TensorIterator<long,double_complex,double>;
    template class TensorIterator<long,double_complex,double_complex>;
    template class TensorIterator<long,double_complex,float>;
    template class TensorIterator<long,double_complex,float_complex>;
    template class TensorIterator<long,double_complex,long>;
    template class TensorIterator<long,float,double>;
    template class TensorIterator<long,float,double_complex>;
    template class TensorIterator<long,float,float>;
    template class TensorIterator<long,float,float_complex>;
    template class TensorIterator<long,float,long>;
    template class TensorIterator<long,float_complex,double>;
    template class TensorIterator<long,float_complex,double_complex>;
    template class TensorIterator<long,float_complex,float>;
    template class TensorIterator<long,float_complex,float_complex>;
    template class TensorIterator<long,float_complex,long>;
    template class TensorIterator<long,long,double>;
    template class TensorIterator<long,long,double_complex>;
    template class TensorIterator<long,long,float>;
    template class TensorIterator<long,long,float_complex>;
    template class TensorIterator<long,long,long>;
    */
}
