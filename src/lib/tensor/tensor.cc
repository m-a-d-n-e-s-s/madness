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

  
#define TENSOR_CC

/// \file tensor.cc
/// \brief Completes the implementation of Tensor and instantiates all specializations for fast compiles.


#ifdef STATIC
#  undef STATIC
#endif

#if HAVE_UNQUALIFIED_STATIC_DECL
#  define STATIC static
#else
// Cray X1 compiler won't instantiate static function templates (mxm*)
#  define STATIC
#endif

#include <algorithm>
#include <complex>
#include <cmath>
#include <cstring>
#include <iostream>

#include <madness_config.h>
#include <tensor/tensor.h>
#include <tensor/random.h>
#include <tensor/mtxmq.h>
#include <tensor/aligned.h>

#if MISSING_POSIX_MEMALIGN_PROTO
  extern "C"  int posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
#endif

#define IS_ODD(n) ((n)&0x1)
#define IS_UNALIGNED(p) (((unsigned long)(p))&0x7)

namespace madness {

#include <tensor/mxm.h>

    template <class T> class Tensor;

    template <class T> class SliceTensor;

    template <> double RandomNumber<double> () {
        //return std::rand()/(RAND_MAX+1.0);
        return genrand_res53();
    }
    template <> float RandomNumber<float> () {
        //return std::rand()/(RAND_MAX+1.0);
        return float(genrand_real2());
    }
    template <> double_complex RandomNumber<double_complex> () {
        return double_complex(RandomNumber<double>(),RandomNumber<double>());
    }
    template <> float_complex RandomNumber<float_complex> () {
        return float_complex(RandomNumber<float>(),RandomNumber<float>());
    }


    /// All new tensors are initialized by init except for the default constructor
    template <class T>
    void Tensor<T>::init(long nd, const long d[], bool dozero) {
        id = TensorTypeData<T>::id;
        TENSOR_ASSERT(nd>0 && nd <= TENSOR_MAXDIM,"invalid ndim in new tensor", nd, 0);
        // sanity check ... 2GB in doubles
        for (int i=0; i<nd; i++) {
            TENSOR_ASSERT(d[i]>=0 && d[i]<268435456, "invalid dimension size in new tensor",d[i],0);
        }
        set_dims_and_size(nd, d);
        if (size) {
            TENSOR_ASSERT(size>=0 && size<268435456, "invalid size in new tensor",size,0);
            try {

#define TENSOR_ALIGNMENT 16
#ifdef HAVE_POSIX_MEMALIGN
                T* q;
                if (posix_memalign((void **) &q, TENSOR_ALIGNMENT, sizeof(T)*size)) throw 1;
                p = SharedPtr<T>(q, ::madness::detail::del_free);
#elif defined(WORLD_GATHER_MEM_STATS)
                p = SharedPtr<T>(new T[size]);
#else
#  error Need a mechanism to allocated aligned memory
#endif

            } catch (...) {
                std::printf("new failed nd=%ld type=%ld size=%ld\n", nd, id, size);
                std::printf("  %ld %ld %ld %ld %ld %ld\n",
                            d[0], d[1], d[2], d[3], d[4], d[5]);
                TENSOR_EXCEPTION("new failed",size,this);
            }
            pointer = p.get();		// Null for ndim=-1
            //std::printf("allocated %p [%ld]  %ld\n", pointer, size, p.use_count());
            if (dozero) {
                // T zero = 0; for (long i=0; i<size; i++) pointer[i] = zero;
                // or
                // memset((void *) pointer, 0, size*sizeof(T));
                // or
                aligned_zero(size, pointer);
            }
        } else {
            pointer = 0;
        }
    }

    /// Routine for internal class use to perform shallow copy from \c t to \c *this
    template <class T>
    void Tensor<T>::internal_shallow_copy(const Tensor<T>& t) {
        ndim = t.ndim;
        size = t.size;
        pointer = t.pointer;
        id = t.id;
        p = t.p;
        for (long i=0; i<ndim; i++) {
            dim[i] = t.dim[i];
            stride[i] = t.stride[i];
        }
        for (long i=ndim; i<TENSOR_MAXDIM; i++) { // So can iterate over missing dimensions
            dim[i] = 1;
            stride[i] = 0;
        }
    }

    /// Construct a 1d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0) : p(0) {
        dim[0] = d0;
        init(1, dim);
    }

    /// Construct a 1d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(int d0) : p(0) {
        dim[0] = d0;
        init(1, dim);
    }

    /// Construct a 2d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1) : p(0){
        dim[0] = d0;
        dim[1] = d1;
        init(2, dim);
    }

    /// Construct a 3d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2) : p(0) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        init(3, dim);
    }

    /// Construct a 4d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3) : p(0) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        init(4, dim);
    }

    /// Construct a 5d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3, long d4) : p(0) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        dim[4] = d4;
        init(5, dim);
    }

    /// Construct a 6d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3, long d4, long d5) : p(0) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        dim[4] = d4;
        dim[5] = d5;
        init(6, dim);
    }

    /// Assigning a scalar to a Tensor fills the tensor with that value
    template <class T>
    Tensor<T>& Tensor<T>::operator=(T x) { // Fill with scalar ... violates shallow = ?
        UNARY_OPTIMIZED_ITERATOR(T,(*this),*_p0 = x);
        return *this;
    }

    /// Tensor + Tensor yielding new Tensor
    template <class T>
    Tensor<T> Tensor<T>::operator+(const Tensor<T>& t) const { // new = this + other
        Tensor<T> result = Tensor<T>(ndim,dim,false);
        TERNARY_OPTIMIZED_ITERATOR(T, result, T, (*this), T, t, *_p0 = *_p1 + *_p2);
        return result;
    }

    /// Tensor - Tensor yielding new Tensor
    template <class T>
    Tensor<T> Tensor<T>::operator-(const Tensor<T>& t) const { // new = this - other
        Tensor<T> result = Tensor<T>(ndim,dim,false);
        TERNARY_OPTIMIZED_ITERATOR(T, result, T, (*this), T, t, *_p0 = *_p1 - *_p2);
        return result;
    }

    /// Tensor + scalar yields new Tensor by adding scalar to each element
    template <class T>
    Tensor<T> Tensor<T>::operator+(T x) const { // new = this+scalar
        Tensor<T> result = Tensor<T>(ndim,dim,false);
        BINARY_OPTIMIZED_ITERATOR(T, result, T, (*this), *_p0 = *_p1 + x);
        return result;
    }

    /// Unary negation yields a new Tensor
    template <class T>
    Tensor<T> Tensor<T>::operator-() const { // new = -this
        Tensor<T> result = Tensor<T>(ndim,dim,false);
        BINARY_OPTIMIZED_ITERATOR(T, result, T, (*this), *(_p0) = - (*_p1));
        return result;
    }

    /// Tensor - scalar yields new Tensor by subtracting scalar from each element
    template <class T>
    Tensor<T> Tensor<T>::operator-(T x) const { // new = this - scalar
        Tensor<T> result = Tensor<T>(ndim,dim,false);
        BINARY_OPTIMIZED_ITERATOR(T, result, T, (*this), *_p0 = *_p1 - x);
        return result;
    }

    /// Inplace increment each element of Tensor by scalar
    template <class T>
    Tensor<T>& Tensor<T>::operator+=(T x) { //  this += scalar
        UNARY_OPTIMIZED_ITERATOR(T,(*this), *_p0 += x);
        return *this;
    }

    /// Inplace decrement each element of Tensor by scalar
    template <class T>
    Tensor<T>& Tensor<T>::operator-=(T x) { //  this -= scalar
        UNARY_OPTIMIZED_ITERATOR(T,(*this), *_p0 -= x);
        return *this;
    }

    /// Return a reference to the element of the Tensor indexed by the vector
    template <class T>
    T& Tensor<T>::operator()(const std::vector<long> i) const {
        TENSOR_ASSERT(i.size()==(unsigned int) ndim,"invalid number of dimensions",i.size(),this);
        long index=0;
        for (long d=0; d<ndim; d++) {
            TENSOR_ASSERT(i[d]>=0 && i[d]<dim[d],"out-of-bounds access",i[d],this);
            index += i[d]*stride[d];
        }
        return pointer[index];
    }

    /// Return a 1d SliceTensor that views the specified range of the 1d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0) const {
        TENSOR_ASSERT(this->ndim==1,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[1] = {s0};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1d SliceTensor that views the specified range of the 2d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(long i, const Slice& s1) const {
        TENSOR_ASSERT(this->ndim==2,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[2] = {Slice(i,i,0),s1};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1d SliceTensor that views the specified range of the 2d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, long j) const {
        TENSOR_ASSERT(this->ndim==2,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[2] = {s0,Slice(j,j,0)};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 2d SliceTensor that views the specified range of the 2d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1) const {
        TENSOR_ASSERT(this->ndim==2,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[2] = {s0,s1};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 3d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1, const Slice& s2) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {s0,s1,s2};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 2d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(long i, const Slice& s1, const Slice& s2) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {Slice(i,i,0),s1,s2};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 2d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, long j, const Slice& s2) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {s0,Slice(j,j,0),s2};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 2d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1, long k) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {s0,s1,Slice(k,k,0)};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(long i, long j, const Slice& s2) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {Slice(i,i,0),Slice(j,j,0),s2};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(long i, const Slice& s1, long k) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {Slice(i,i,0),s1,Slice(k,k,0)};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1d SliceTensor that views the specified range of the 3d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, long j, long k) const {
        TENSOR_ASSERT(this->ndim==3,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[3] = {s0,Slice(j,j,0),Slice(k,k,0)};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1-4d SliceTensor that views the specified range of the 4d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                         const Slice& s3) const {
        TENSOR_ASSERT(this->ndim==4,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[4] = {s0,s1,s2,s3};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1-5d SliceTensor that views the specified range of the 5d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                         const Slice& s3, const Slice& s4) const {
        TENSOR_ASSERT(this->ndim==5,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[5] = {s0,s1,s2,s3,s4};
        return SliceTensor<T>(*this,s);
    }

    /// Return a 1-6d SliceTensor that views the specified range of the 6d Tensor
    template <class T>
    SliceTensor<T> Tensor<T>::operator()(const Slice& s0, const Slice& s1, const Slice& s2,
                                         const Slice& s3, const Slice& s4, const Slice& s5) const {
        TENSOR_ASSERT(this->ndim==6,"invalid number of dimensions",
                      this->ndim,this);
        Slice s[6] = {s0,s1,s2,s3,s4,s5};
        return SliceTensor<T>(*this,s);
    }

    /// Return a new Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(const std::vector<long>& d) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(d);
        return result;
    }

    /// Return a new 1d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0));
        return result;
    }

    /// Return a new 2d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0, long dim1) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0,dim1));
        return result;
    }

    /// Return a new 3d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0, long dim1, long dim2) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0,dim1,dim2));
        return result;
    }

    /// Return a new 4d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0, long dim1, long dim2,
                                 long dim3) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0,dim1,dim2,dim3));
        return result;
    }

    /// Return a new 5d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0, long dim1, long dim2,
                                 long dim3, long dim4) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0,dim1,dim2,dim3,dim4));
        return result;
    }

    /// Return a new 6d Tensor using a different, but conforming, shape to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::reshape(long dim0, long dim1, long dim2,
                                 long dim3, long dim4, long dim5) const {
        Tensor<T> result(*this);
        result.reshape_inplace_base(vector_factory(dim0,dim1,dim2,dim3,dim4,dim5));
        return result;
    }

    /// Return a new 1d Tensor to view a contiguous tensor
    template <class T>
    Tensor<T> Tensor<T>::flat() const {
        Tensor<T> result(*this);
        result.flat_inplace_base();
        return result;
    }

    /// Return a new Tensor view by splitting dimension i. dimi0*dimi1 must equal the original size
    template <class T>
    Tensor<T> Tensor<T>::splitdim(long i, long dimi0, long dimi1) const {
        Tensor<T> result(*this);
        result.splitdim_inplace_base(i, dimi0, dimi1);
        return result;
    }

    /// Return a new Tensor view by swapping dimensions i and j
    template <class T>
    Tensor<T> Tensor<T>::swapdim(long i, long j) const {
        Tensor<T> result(*this);
        result.swapdim_inplace_base(i, j);
        return result;
    }

    /// Return a new Tensor view by fusing contiguous dimensions i and i+1
    template <class T>
    Tensor<T> Tensor<T>::fusedim(long i) const {
        Tensor<T> result(*this);
        result.fusedim_inplace_base(i);
        return result;
    }

    /// Return a new Tensor view by cyclically permutating indices (start,...,end) nshift times
    template <class T>
    Tensor<T> Tensor<T>::cycledim(long nshift, long start, long end) const {
        Tensor<T> result(*this);
        result.cycledim_inplace_base(nshift, start, end);
        return result;
    }

    /// Return a new Tensor view with a general permutaiton ... idim_new = map[idim_old]
    template <class T>
    Tensor<T> Tensor<T>::mapdim(const std::vector<long>& map) const {
        Tensor<T> result(*this);
        result.mapdim_inplace_base(map);
        return result;
    }

    /// Fill the Tensor with random values ... see RandomNumber for details
    template <class T>
    Tensor<T>& Tensor<T>::fillrandom() {
        UNARY_OPTIMIZED_ITERATOR(T,(*this), *_p0 = RandomNumber<T>());
        return *this;
    }

    /// Fill the tensor with a scalar
    template <class T>
    Tensor<T>& Tensor<T>::fill(T x) {
        UNARY_OPTIMIZED_ITERATOR(T,(*this), *_p0 = x);
        return *this;
    }

    /// Fill the tensor with the canonical index values 0,1,2,... in row major order
    template <class T>
    Tensor<T>& Tensor<T>::fillindex() {
        long count = 0;
        UNARY_UNOPTIMIZED_ITERATOR(T,(*this), *_p0 = count++); // Fusedim would be OK
        return *this;
    }

    /// Zero elements of the tensor less than x in absolute value
    template <class T>
    Tensor<T>& Tensor<T>::screen(double x) {
        T zero = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this), if (std::abs(*_p0)<x) *_p0=zero);
        return *this;
    }

    template <class T>
    T Tensor<T>::sum() const {
        T result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += *_p0);
        return result;
    }

    template <class T>
    T Tensor<T>::sumsq() const {
        T result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += (*_p0) * (*_p0));
        return result;
    }

    template <class T>
    T Tensor<T>::product() const {
        T result = 1;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result *= *_p0);
        return result;
    }

    template <class T>
    T Tensor<T>::min(long* ind) const {
        T result = *(this->pointer);
        if (ind) {
            for (long i=0; i<ndim; i++) ind[i]=0;
            long nd = ndim-1;
            UNARY_UNOPTIMIZED_ITERATOR(T,(*this),
                                       if (result > *_p0) {
                                       result = *_p0;
                                       for (long i=0; i<nd; i++) ind[i]=iter.ind[i];
                                           ind[nd] = _j;
                                       }
                                      );
        } else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::min<T>(result,*_p0));
        }
        return result;
    }

    template <class T>
    T Tensor<T>::max(long* ind) const {
        T result = *(this->pointer);
        if (ind) {
            for (long i=0; i<ndim; i++) ind[i]=0;
            long nd = ndim-1;
            UNARY_UNOPTIMIZED_ITERATOR(T,(*this),
                                       if (result < *_p0) {
                                       result = *_p0;
                                       for (long i=0; i<nd; i++) ind[i]=iter.ind[i];
                                           ind[nd] = _j;
                                       }
                                      );
        } else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::max<T>(result,*_p0));
        }
        return result;
    }

    template <class T>
    typename Tensor<T>::float_scalar_type Tensor<T>::normf() const {
        float_scalar_type result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += (*_p0) * (*_p0));
        return (float_scalar_type) std::sqrt(result);
    }

    template <class T>
    typename Tensor<T>::scalar_type Tensor<T>::absmin(long* ind) const {
        scalar_type result = std::abs(*(this->pointer));
        if (ind) {
            for (long i=0; i<ndim; i++) ind[i]=0;
            long nd = ndim-1;
            UNARY_UNOPTIMIZED_ITERATOR(T,(*this),
                                       scalar_type absval = std::abs(*_p0);
                                       if (result > absval) {
                                       result = absval;
                                       for (long i=0; i<nd; i++) ind[i]=iter.ind[i];
                                           ind[nd] = _j;
                                       }
                                      );
        } else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::min<scalar_type>(result,std::abs(*_p0)));
        }
        return result;
    }

    template <class T>
    typename Tensor<T>::scalar_type Tensor<T>::absmax(long* ind) const {
        scalar_type result = std::abs(*(this->pointer));
        if (ind) {
            for (long i=0; i<ndim; i++) ind[i]=0;
            long nd = ndim-1;
            UNARY_UNOPTIMIZED_ITERATOR(T,(*this),
                                       scalar_type absval = std::abs(*_p0);
                                       if (result < absval) {
                                       result = absval;
                                       for (long i=0; i<nd; i++) ind[i]=iter.ind[i];
                                           ind[nd] = _j;
                                       }
                                      );
        } else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::max<scalar_type>(result,std::abs(*_p0)));
        }
        return result;
    }

    template <class T>
    Tensor<T>& Tensor<T>::unaryop(T (*op) (T)) {
        UNARY_OPTIMIZED_ITERATOR(T,(*this),*_p0=op(*_p0));
        return *this;
    }

    template <class T>
    T Tensor<T>::trace(const Tensor<T>& t) const {
        T result = 0;
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t,result += (*_p0)*(*_p1));
        return result;
    }


    template <class T>
    Tensor<T>& Tensor<T>::emul(const Tensor<T>& t) {
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t,*_p0 *= *_p1);
        return *this;
    }

    template <class T>
    Tensor<T>& Tensor<T>::gaxpy(T alpha, const Tensor<T>& t, T beta) {
        //BINARYITERATOR(T,(*this),T,t, (*_p0) = alpha*(*_p0) + beta*(*_p1));
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t, (*_p0) = alpha*(*_p0) + beta*(*_p1));
        //ITERATOR((*this),(*this)(IND) = alpha*(*this)(IND) + beta*t(IND));
        return *this;
    }

    /// Construct a new SliceTensor to view a slice of a Tensor
    template <class T>
    SliceTensor<T>::SliceTensor(const Tensor<T>& t, const Slice s[]) : Tensor<T>(t) {
        // C++ standard says class derived from parameterized base class cannot
        // directly access the base class elements ... must explicitly reference.
        long nd = 0, size=1;
        for (long i=0; i<t.ndim; i++) {
            long start=s[i].start, end=s[i].end, step=s[i].step;
            //std::printf("%ld input start=%ld end=%ld step=%ld\n",
            //i, start, end, step);
            if (start < 0) start += this->dim[i];
            if (end < 0) end += this->dim[i];
            long len = end-start+1;
            if (step) len /= step;	// Rounds len towards zero

            // if input length is not exact multiple of step, round end towards start
            // for the same behaviour of for (i=start; i<=end; i+=step);
            end = start + (len-1)*step;

            //std::printf("%ld munged start=%ld end=%ld step=%ld len=%ld dim=%ld\n",
            //		i, start, end, step, len, this->dim[i]);


            TENSOR_ASSERT(start>=0 && start<this->dim[i],"slice start invalid",start,this);
            TENSOR_ASSERT(end>=0 && end<this->dim[i],"slice end invalid",end,this);
            TENSOR_ASSERT(len>0,"slice length must be non-zero",len,this);

            this->pointer += start * t.stride[i];

            if (step) {
                size *= len;
                this->dim[nd] = len;
                this->stride[nd] = step * t.stride[i];
                nd++;
            }
        }
        //For Python interface need to be able to return a scalar inside a tensor with nd=0
        //TENSOR_ASSERT(nd>0,"slicing produced a scalar, but cannot return one",nd,this);
        for (long i=nd; i<TENSOR_MAXDIM; i++) { // So can iterate over missing dimensions
            this->dim[i] = 1;
            this->stride[i] = 0;
        }
        this->ndim = nd;
        this->size = size;
    }

    /// Assigning to a SliceTensor causes a deep copy ... the data is copied
    template <class T>
    SliceTensor<T>& SliceTensor<T>::operator=(const SliceTensor<T>& t) {
        BINARY_OPTIMIZED_ITERATOR(T, (*this), T, t, *_p0 = (*_p1));
        return *this;
    }

    /// Assigning to a SliceTensor causes a deep copy ... the data is copied
    template <class T>
    SliceTensor<T>& SliceTensor<T>::operator=(const Tensor<T>& t) {
        BINARY_OPTIMIZED_ITERATOR(T, (*this), T, t, *_p0 = (*_p1));
        return *this;
    }

    /// Assigning a scalar to a SliceTensor fills the tensor with that value
    template <class T>
    SliceTensor<T>& SliceTensor<T>::operator=(T x) {
        UNARY_OPTIMIZED_ITERATOR(T, (*this), *_p0 = x);
        return *this;
    }

    /// Print (for human consumption) a tensor to the stream
    template <class T>
    std::ostream& operator << (std::ostream& s, const Tensor<T>& t) {
        using namespace std;

        long maxdim = 0;
        long index_width = 0;
        for (int i = 0; i<(t.ndim-1); i++) {
            if (maxdim < t.dim[i]) maxdim = t.dim[i];
        }
        if (maxdim < 10)
            index_width = 1;
        else if (maxdim < 100)
            index_width = 2;
        else if (maxdim < 1000)
            index_width = 3;
        else if (maxdim < 10000)
            index_width = 4;
        else
            index_width = 6;

        ios::fmtflags oldflags = s.setf(ios::scientific);
        long oldprec = s.precision();
        long oldwidth = s.width();

        // C++ formatted IO is worse than Fortran !!
        for (TensorIterator<T> iter=t.unary_iterator(1,false,false); iter!=t.end(); ++iter) {
            const T* p = iter._p0;
            long inc = iter._s0;
            long dimj = iter.dimj;
            s.unsetf(ios::scientific);
            s << '[';
            for (long i=0; i<iter.ndim; i++) {
                s.width(index_width);
                s << iter.ind[i];
                if (i != iter.ndim) s << ",";
            }
            s << "*]";
            s.setf(ios::scientific);
            for (long j=0; j<dimj; j++, p+=inc) {
                s.precision(4);
                s.width(12);
                s << *p;
            }
            s.unsetf(ios::scientific);
            s << endl;
        }
        s.setf(oldflags);
        s.precision(oldprec);
        s.width(oldwidth);

        return s;
    }

    /// Print a TensorException to the stream (for human consumption)
    std::ostream& operator <<(std::ostream& out, const TensorException& e) {
        out << "TensorException: msg='";
        if (e.msg) out << e.msg;
        out << "'\n";
        if (e.assertion) out << "                 failed assertion='" <<
            e.assertion << "'\n";
        out << "                 value=" << e.value << "\n";
        if (e.line) out << "                 line=" << e.line << "\n";
        if (e.function) out << "                 function='" <<
            e.function << "'\n";
        if (e.filename) out << "                 filename='" <<
            e.filename << "'\n";
        if (e.t) {
            out << "                 tensor=Tensor<";
            if (e.t->id>=0 && e.t->id<=TENSOR_MAX_TYPE_ID) {
                out << tensor_type_names[e.t->id] << ">(";
            } else {
                out << "invalid_type_id>(";
            }
            if (e.t->ndim>=0 && e.t->ndim<TENSOR_MAXDIM) {
                for (int i=0; i<e.t->ndim; i++) {
                    out << e.t->dim[i];
                    if (i != (e.t->ndim-1)) out << ",";
                }
                out << ")";
            } else {
                out << "invalid_dimensions)";
            }
            out << " at 0x" << (void *) (e.t) << "\n";
        }

        return out;
    }


    /// Specializations of normf, min, max, for complex types ...

    /// Errors will be resported at link time whereas if we inlined the
    /// code they would be detected at compile time ... but then you have
    /// to processes this stuff in every file.
    template <> Tensor<float_complex>::scalar_type Tensor<float_complex>::normf() const {
        double result = 0.0;
        UNARY_OPTIMIZED_ITERATOR(float_complex,(*this),result += std::norm(*_p0));
        return (scalar_type) std::sqrt(result);
    }

    template <> Tensor<double_complex>::scalar_type Tensor<double_complex>::normf() const {
        double result = 0.0;
        UNARY_OPTIMIZED_ITERATOR(double_complex,(*this),result += std::norm(*_p0));
        return std::sqrt(result);
    }
    template<> float_complex Tensor<float_complex>::min(long* ind) const {
        TENSOR_EXCEPTION("cannot perform min on complex types",0,this);
        return 0;
    }
    template<> double_complex Tensor<double_complex>::min(long* ind) const {
        TENSOR_EXCEPTION("cannot perform min on complex types",0,this);
        return 0;
    }
    template<> float_complex Tensor<float_complex>::max(long* ind) const {
        TENSOR_EXCEPTION("cannot perform max on complex types",0,this);
        return 0;
    }
    template<> double_complex Tensor<double_complex>::max(long* ind) const {
        TENSOR_EXCEPTION("cannot perform max on complex types",0,this);
        return 0;
    }

    // Functions operating upon tensors ... below, we only instantiate
    // the type appropriate functions.

    /// Make a deep copy of the tensor (i.e., make an independent and contiguous copy of the data)
    template <class T>
    Tensor<T> copy(const Tensor<T>& t) {	// Deep, contiguous copy of t
        if (t.size) {
            Tensor<T> result = Tensor<T>(t.ndim,t.dim,false);
            BINARY_OPTIMIZED_ITERATOR(T, result, T, t, *_p0 = *_p1);
            return result;
        }
        else {
            return Tensor<T>();
        }
    }

    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)
    template <class T>
    Tensor<T> outer(const Tensor<T>& left, const Tensor<T>& right) {
        long nd = left.ndim + right.ndim;
        TENSOR_ASSERT(nd <= TENSOR_MAXDIM,"too many dimensions in result",
                      nd,0);
        long d[TENSOR_MAXDIM];
        for (long i=0; i<left.ndim; i++) d[i] = left.dim[i];
        for (long i=0; i<right.ndim; i++) d[i+left.ndim] = right.dim[i];
        Tensor<T> result(nd,d,false);
        T* ptr = result.ptr();

        TensorIterator<T> iter=right.unary_iterator(1,false,true);
        for (TensorIterator<T> p=left.unary_iterator(); p!=left.end(); ++p) {
            T val1 = *p;
            // Cannot reorder dimensions, but can fuse contiguous dimensions
            for (iter.reset(); iter._p0; ++iter) {
                long dimj = iter.dimj;
                T* _p0 = iter._p0;
                long Tstride = iter._s0;
                for (long _j=0; _j<dimj; _j++, _p0+=Tstride) {
                    *ptr++ = val1 * (*_p0);
                }
            }
        }

        return result;
    }


    /// Inner product ... result(i,j,...,p,q,...) = sum(z) left(i,j,...,z)*right(z,p,q,...)

    /// By default it contracts the last dimension of the left tensor and
    /// the first dimension of the right tensor.  These defaults can be
    /// changed by specifying \c k0 and \c k1 , the index to contract in
    /// the left and right side tensors, respectively.  The defaults
    /// correspond to (\c k0=-1 and \c k1=0 ).
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> inner(const Tensor<T>& left, const Tensor<Q>& right,
                    long k0, long k1) {
        if (k0 < 0) k0 += left.ndim;
        if (k1 < 0) k1 += right.ndim;
        long nd = left.ndim + right.ndim - 2;
        TENSOR_ASSERT(nd!=0, "result is a scalar but cannot return one ... use dot",
                      nd, &left);
        TENSOR_ASSERT(left.dim[k0] == right.dim[k1],"common index must be same length",
                      right.dim[k1], &left);

        TENSOR_ASSERT(nd > 0 && nd <= TENSOR_MAXDIM,
                      "invalid number of dimensions in the result", nd,0);

        long d[TENSOR_MAXDIM];

        long base=0;
        for (long i=0; i<k0; i++) d[i] = left.dim[i];
        for (long i=k0+1; i<left.ndim; i++) d[i-1] = left.dim[i];
        base = left.ndim-1;
        for (long i=0; i<k1; i++) d[i+base] = right.dim[i];
        base--;
        for (long i=k1+1; i<right.ndim; i++) d[i+base] = right.dim[i];

        Tensor<TENSOR_RESULT_TYPE(T,Q)> result(nd,d);

        inner_result(left,right,k0,k1,result);

        return result;
    }

    /// Accumulate inner product into user provided, contiguous, correctly sized result tensor

    /// This routine may be used to optimize away the tensor constructor
    /// of the result tensor in inner loops when the result tensor may be
    /// reused or accumulated into.  If the user calls this routine
    /// directly very little checking is done since it is intended as an
    /// optmization for small tensors.  As far as the result goes, the
    /// caller is completely responsible for providing a contiguous tensor
    /// that has the correct dimensions and is appropriately initialized.
    /// The inner product is accumulated into result.
    template <class T, class Q>
    void inner_result(const Tensor<T>& left, const Tensor<Q>& right,
                      long k0, long k1, Tensor< TENSOR_RESULT_TYPE(T,Q) >& result) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        // Need to include explicit optimizations for common special cases
        // E.g., contiguous, matrix-matrix, and 3d-tensor*matrix

        resultT* ptr = result.ptr();

        if (k0 < 0) k0 += left.ndim;
        if (k1 < 0) k1 += right.ndim;

        if (left.iscontiguous() && right.iscontiguous()) {
            if (k0==0 && k1==0) {
                // c[i,j] = a[k,i]*b[k,j] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimj = right.stride[0];
                long dimi = left.stride[0];
                mTxm(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            } else if (k0==(left.ndim-1) && k1==(right.ndim-1)) {
                // c[i,j] = a[i,k]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.size/dimk;
                long dimj = right.size/dimk;
                mxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            } else if (k0==0 && k1==(right.ndim-1)) {
                // c[i,j] = a[k,i]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.stride[0];
                long dimj = right.size/dimk;
                mTxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            } else if (k0==(left.ndim-1) && k1==0) {
                // c[i,j] = a[i,k]*b[k,j] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.size/dimk;
                long dimj = right.stride[0];
                mxm(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
        }

        long dimj = left.dim[k0];
        TensorIterator<Q> iter1=right.unary_iterator(1,false,false,k1);

        for (TensorIterator<T> iter0=left.unary_iterator(1,false,false,k0);
                iter0._p0; ++iter0) {
            T* xp0 = iter0._p0;
            long s0 = iter0._s0;
            for (iter1.reset(); iter1._p0; ++iter1) {
                T* p0 = xp0;
                Q* p1 = iter1._p0;
                long s1 = iter1._s0;
                resultT sum = 0;
                for (long j=0; j<dimj; j++,p0+=s0,p1+=s1) {
                    sum += (*p0) * (*p1);
                }
                *ptr++ += sum;
            }
        }
    }

    /// Transform all dimensions of the tensor t by the matrix c

    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
    /// The input dimensions of \c t must all be the same and agree with
    /// the first dimension of \c c .  The dimensions of \c may differ in
    /// size.  If the dimensions of \c c are the same, and the operation
    /// is being performed repeatedly, then you might consider calling \c
    /// fast_transform instead which enables additional optimizations and
    /// can eliminate all constructor overhead and improve cache locality.
    ///
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> transform(const Tensor<T>& t, const Tensor<Q>& c) {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        TENSOR_ASSERT(c.ndim == 2,"second argument must be a matrix",c.ndim,&c);
        if (c.dim[0]==c.dim[1] && t.iscontiguous() && c.iscontiguous()) {
            Tensor<resultT> result(t.ndim,t.dim,false);
            Tensor<resultT> work(t.ndim,t.dim,false);
            return fast_transform(t, c, result, work);
        }
        else {
            Tensor<resultT> result = t;
            for (long i=0; i<t.ndim; i++) {
                result = inner(result,c,0,0);
            }
            return result;
        }
    }

    /// Transform all dimensions of the tensor t by distinct matrices c

    /// Similar to transform but each dimension is transformed with a 
    /// distinct matrix.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c[0](i',i) c[1](j',j) c[2](k',k) ...
    /// \endcode
    /// The first dimension of the matrices c must match the corresponding
    /// dimension of t.
    template <class T, class Q>
    Tensor<TENSOR_RESULT_TYPE(T,Q)> general_transform(const Tensor<T>& t, const Tensor<Q> c[]) {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        Tensor<resultT> result = t;
        for (long i=0; i<t.ndim; i++) {
            result = inner(result,c[i],0,0);
        }
        return result;
    }

    /// Restricted but heavily optimized form of transform()

    /// Both dimensions of \c c must be the same and match all dimensions
    /// of the input tensor \c t.  All tensors must be contiguous.
    ///
    /// Performs the same operation as \c transform but it requires
    /// that the caller pass in workspace and a preallocated result,
    /// hoping that that both can be reused.  If the result and
    /// workspace are reused between calls, then no tensor
    /// constructors need be called and cache locality should be
    /// improved.  By passing in the workspace, this routine is kept
    /// thread safe.
    ///
    /// The input, result and workspace tensors must be distinct.
    ///
    /// All input tensors must be contiguous and fastest execution
    /// will result if all dimensions are even and data is aligned on
    /// 16-byte boundaries.  The workspace and the result must be of
    /// the same size as the input \c t .  The result tensor need not
    /// be initialized before calling fast_transform.
    ///
    /// \code 
    ///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...  
    /// \endcode
    ///
    /// The input dimensions of \c t must all be the same .
    template <class T, class Q>
    Tensor< TENSOR_RESULT_TYPE(T,Q) >& fast_transform(const Tensor<T>& t, const Tensor<Q>& c,  Tensor< TENSOR_RESULT_TYPE(T,Q) >& result,
                              Tensor< TENSOR_RESULT_TYPE(T,Q) >& workspace) {
        typedef  TENSOR_RESULT_TYPE(T,Q) resultT;
        const Q *pc=c.ptr();
        resultT *t0=workspace.ptr(), *t1=result.ptr();
        if (t.ndim&1) {
            t0 = result.ptr();
            t1 = workspace.ptr();
        }

        long dimj = c.dim[1];
        long dimi = 1;
        for (int n=1; n<t.ndim; n++) dimi *= dimj;
        long nij = dimi*dimj;


        if (IS_ODD(dimi) || IS_ODD(dimj) || 
            IS_UNALIGNED(pc) || IS_UNALIGNED(t0) || IS_UNALIGNED(t1)) {
            for (long i=0; i<nij; i++) t0[i] = 0.0;
            mTxm(dimi, dimj, dimj, t0, t.ptr(), pc);
            for (int n=1; n<t.ndim; n++) {
                for (long i=0; i<nij; i++) t1[i] = 0.0;
                mTxm(dimi, dimj, dimj, t1, t0, pc);
                std::swap(t0,t1);
            }
        }
        else {
            mTxmq(dimi, dimj, dimj, t0, t.ptr(), pc); 
            for (int n=1; n<t.ndim; n++) {
                mTxmq(dimi, dimj, dimj, t1, t0, pc);
                std::swap(t0,t1);
            }
        }

        return result;
    }

    /// Return a new tensor holding the absolute value of each element of t
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > abs(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,T,t,*_p0 = std::abs(*_p1));
        return result;
    }

    /// Return a new tensor holding the argument of each element of t (complex types only)
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > arg(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,T,t,*_p0 = std::arg(*_p1));
        return result;
    }

    /// Return a new tensor holding the real part of each element of t (complex types only)
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > real(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,T,t,*_p0 = std::real(*_p1));
        return result;
    }

    /// Return a new tensor holding the imaginary part of each element of t (complex types only)
    template <class T>
    Tensor< typename Tensor<T>::scalar_type > imag(const Tensor<T>& t) {
        typedef typename Tensor<T>::scalar_type scalar_type;
        Tensor<scalar_type> result(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(scalar_type,result,T,t,*_p0 = std::imag(*_p1));
        return result;
    }

    /// Returns a new deep copy of the complex conjugate of the input tensor (complex types only)
    template <class T>
    Tensor<T> conj(const Tensor<T>& t) {
        Tensor<T> result(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(T,result,T,t,*_p0 = std::conj(*_p1));
        return result;
    }

    std::ostream& operator<<(std::ostream& stream, const Slice& s) {
        stream << "Slice(" << s.start << "," << s.end << "," << s.step << ")";
        return stream;
    }


#include "tensor_spec.h"
}
