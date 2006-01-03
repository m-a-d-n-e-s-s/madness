#define TENSOR_CC

/// \file tensor.cc
/// \brief Completes the implementation of Tensor and instantiates all specializations for fast compiles.

#ifdef _CRAY
// Cray won't instantiate static function templates (mxm*)
#define STATIC 
#elif defined(IBMXLC)
// Compiler will only use fully resolved static names?
#define STATIC 
#else
#define STATIC static
//#define STATIC 
#endif

//#include <cstdio>
//#include <cstdlib>
#include <algorithm>

#include <complex>

#include <cmath>

#include <iostream>

#include "tensor.h"

#include "mtrand.h"

namespace madness {

    
    template <class T>
    T RandomNumber() {
        //return (T) std::rand();
        return (T) genrand_int31();
    }
    template <> double RandomNumber<double> () {
        //return std::rand()/(RAND_MAX+1.0);
        return genrand_res53();
    }
    template <> float RandomNumber<float> () {
        //return std::rand()/(RAND_MAX+1.0);
        return genrand_real2();
    }
    template <> double_complex RandomNumber<double_complex> () {
        return double_complex(RandomNumber<double>(),RandomNumber<double>());
    }
    template <> float_complex RandomNumber<float_complex> () {
        return float_complex(RandomNumber<float>(),RandomNumber<float>());
    }
    
    template <class T> class Tensor;
    
    template <class T> class SliceTensor;
    
    /// All new tensors are initialized by init except for the default constructor
    template <class T>
    void Tensor<T>::init(long nd, const long d[], bool dozero) {
        id = TensorTypeData<T>::id;
        TENSOR_ASSERT(nd>0 && nd <= TENSOR_MAXDIM,"invalid ndim in new tensor", nd, 0);
        for (int i=0; i<nd; i++) {
           TENSOR_ASSERT(d[i]>=0 && d[i]<200000, "invalid dimension size in new tensor",i,0);
        }
        set_dims_and_size(nd, d);
        if (size) {
            try {
                p = shared_array<T>(new T[size]);
            }
            catch (...) {
                std::printf("new failed nd=%ld type=%ld size=%ld\n", nd, id, size);
                std::printf("  %ld %ld %ld %ld %ld %ld\n", 
                            d[0], d[1], d[2], d[3], d[4], d[5]);
                TENSOR_EXCEPTION("new failed",size,this);
            }
            pointer = p.get();		// Null for ndim=-1
            //std::printf("allocated %p [%ld]  %ld\n", pointer, size, p.use_count());
            if (dozero) {
                T zero = 0;
                for (long i=0; i<size; i++) pointer[i] = zero;
            }
        }
        else {
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
    Tensor<T>::Tensor(long d0) {
        dim[0] = d0;
        init(1, dim);
    }
    
    /// Construct a 1d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(int d0) {
        dim[0] = d0;
        init(1, dim);
    }
    
    /// Construct a 2d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1) {
        dim[0] = d0;
        dim[1] = d1;
        init(2, dim);
    }
    
    /// Construct a 3d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        init(3, dim);
    }
    
    /// Construct a 4d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        init(4, dim);
    }
    
    /// Construct a 5d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3, long d4) {
        dim[0] = d0;
        dim[1] = d1;
        dim[2] = d2;
        dim[3] = d3;
        dim[4] = d4;
        init(5, dim);
    }
    
    /// Construct a 6d tensor initialized to zero
    template <class T>
    Tensor<T>::Tensor(long d0, long d1, long d2, long d3, long d4, long d5) {
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
    
    /// Return the sum of all elements of the tensor
    template <class T>
    T Tensor<T>::sum() const {
        T result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += *_p0);
        return result;
    }
    
    /// Return the sum of the squares of all elements of the tensor
    template <class T>
    T Tensor<T>::sumsq() const {
        T result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += (*_p0) * (*_p0));
        return result;
    }
    
    /// Return the product of all elements of the tensor
    template <class T>
    T Tensor<T>::product() const {
        T result = 1;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result *= *_p0);
        return result;
    }
    
    /// Return the minimum value (and if ind is non-null, its index) in the Tensor
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
                                       });
        }
        else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::min<T>(result,*_p0));
        }
        return result;
    }
    
    /// Return the maximum value (and if ind is non-null, its index) in the Tensor
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
                                       });
        }
        else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::max<T>(result,*_p0));
        }
        return result;
    }
    
    /// Return the Frobenius norm of the tensor (needs correcting for complex??)
    template <class T>
    typename Tensor<T>::float_scalar_type Tensor<T>::normf() const {
        float_scalar_type result = 0;
        UNARY_OPTIMIZED_ITERATOR(T,(*this),result += (*_p0) * (*_p0));
        return (float_scalar_type) std::sqrt(result);
    }
    
    /// Return the absolute minimum value (and if ind is non-null, its index) in the Tensor
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
                                       });
        }
        else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::min<scalar_type>(result,std::abs(*_p0)));
        }
        return result;
    }
    
    /// Return the absolute maximum value (and if ind is non-null, its index) in the Tensor
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
                                       });
        }
        else {
            UNARY_OPTIMIZED_ITERATOR(T,(*this),result=std::max<scalar_type>(result,std::abs(*_p0)));
        }
        return result;
    }
    
    /// Inplace apply a unary function to each element of the tensor
    template <class T>
    Tensor<T>& Tensor<T>::unaryop(T (*op) (T)) {
        UNARY_OPTIMIZED_ITERATOR(T,(*this),*_p0=op(*_p0));
        return *this;
    }
    
    /// Return the trace of two tensors (no complex conjugate invoked)
    template <class T>
    T Tensor<T>::trace(const Tensor<T>& t) const {
        T result = 0;
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t,result += (*_p0)*(*_p1));
        return result;
    }
    
    /// Inplace multiply by corresponding elements of argument Tensor
    template <class T>
    Tensor<T>& Tensor<T>::emul(const Tensor<T>& t) {
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t,*_p0 *= *_p1);
        return *this;
    }
    
    /// Inplace generalized saxpy ... this = this*alpha + other*beta
    template <class T>
    Tensor<T>& Tensor<T>::gaxpy(T alpha, const Tensor<T>& t, T beta) {
        BINARY_OPTIMIZED_ITERATOR(T,(*this),T,t, (*_p0) = alpha*(*_p0) + beta*(*_p1));
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
            }
            else {
                out << "invalid_type_id>(";
            }
            if (e.t->ndim>=0 && e.t->ndim<TENSOR_MAXDIM) {
                for (int i=0; i<e.t->ndim; i++) {
                    out << e.t->dim[i];
                    if (i != (e.t->ndim-1)) out << ",";
                }
                out << ")";
            }
            else {
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
        Tensor<T> result = Tensor<T>(t.ndim,t.dim,false);
        BINARY_OPTIMIZED_ITERATOR(T, result, T, t, *_p0 = *_p1);
        return result;
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
    
    /// Matrix transpose * matrix ... reference implementation (slow but correct)
    template <typename T> 
    STATIC inline void mTxm_reference(long dimi, long dimj, long dimk, 
                                      T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticpated application
        */
        
        for (long k=0; k<dimk; k++) {
            for (long i=0; i<dimi; i++) {
                for (long j=0; j<dimj; j++) {
                    c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }

    /// Matrix transpose * matrix ... double_complex * double

    /// This version is not optimized ... it can be greatly accelerated
    void mTxm(long dimi, long dimj, long dimk, 
              double_complex* RESTRICT c, 
              const double_complex* RESTRICT a, const double* RESTRICT b) 
    {
        for (long k=0; k<dimk; k++) {
            for (long i=0; i<dimi; i++) {
                for (long j=0; j<dimj; j++) {
                    c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }

    
#ifdef _CRAY
    /// Matrix transpose * matrix (Cray version)
    template <typename T> 
    STATIC void mTxm(long dimi, long dimj, long dimk, 
                     T* RESTRICT c, const T* RESTRICT a, const T* RESTRICT b) 
    {
        for (long j=0; j<dimj; j++) {
            for (long k=0; k<dimk; k++) {
                for (long i=0; i<dimi; i++) {
                    c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }
#else
    /// Matrix transpose * matrix (hand unrolled version)
    template <typename T> 
    STATIC inline void mTxm(long dimi, long dimj, long dimk, 
                            T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticpated application
          
          4-way unrolled k loop ... empirically fastest on PIII 
          compared to 2/3 way unrolling (though not by much).
        */
        
        long dimk4 = (dimk/4)*4;
        for (long i=0; i<dimi; i++,c+=dimj) {
            const T* ai = a+i;
            const T* p = b;
            for (long k=0; k<dimk4; k+=4,ai+=4*dimi,p+=4*dimj) {
                T ak0i = ai[0   ];
                T ak1i = ai[dimi];
                T ak2i = ai[dimi+dimi];
                T ak3i = ai[dimi+dimi+dimi];
                const T* bk0 = p;
                const T* bk1 = p+dimj;
                const T* bk2 = p+dimj+dimj;
                const T* bk3 = p+dimj+dimj+dimj;
                for (long j=0; j<dimj; j++) {
                    c[j] += ak0i*bk0[j] + ak1i*bk1[j] + ak2i*bk2[j] + ak3i*bk3[j];
                }
            }
            for (long k=dimk4; k<dimk; k++) {
                T aki = a[k*dimi+i];
                const T* bk = b+k*dimj;
                for (long j=0; j<dimj; j++) {
                    c[j] += aki*bk[j];
                }
            }
        }
    }
#endif
    
    /// Matrix * matrix transpose ... reference implementation (slow but correct)
    template <typename T> 
    STATIC inline void mxmT_reference(long dimi, long dimj, long dimk, 
                                      T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticpated application
        */
        
        for (long i=0; i<dimi; i++) {
            for (long j=0; j<dimj; j++) {
                T sum = 0;
                for (long k=0; k<dimk; k++) {
                    sum += a[i*dimk+k]*b[j*dimk+k];
                }
                c[i*dimj+j] += sum;
            }
        }
    }
    
    /// Matrix * matrix transpose (hand unrolled version)
    template <typename T> 
    STATIC inline void mxmT(long dimi, long dimj, long dimk, 
                            T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          j loop might be long in anticpated application
          
          Unrolled i loop.  Empirically fastest on PIII compared
          to unrolling j, or both i&j.  
        */
        
        long dimi2 = (dimi/2)*2;
        for (long i=0; i<dimi2; i+=2) {
            const T* ai0 = a+i*dimk;
            const T* ai1 = a+i*dimk+dimk;
            T* ci0 = c+i*dimj;
            T* ci1 = c+i*dimj+dimj;
            for (long j=0; j<dimj; j++) {
                T sum0 = 0;
                T sum1 = 0;
                const T* bj = b + j*dimk;
                for (long k=0; k<dimk; k++) {
                    sum0 += ai0[k]*bj[k];
                    sum1 += ai1[k]*bj[k];
                }
                ci0[j] += sum0;
                ci1[j] += sum1;
            }
        }
        for (long i=dimi2; i<dimi; i++) {
            const T* ai = a+i*dimk;
            T* ci = c+i*dimj;
            for (long j=0; j<dimj; j++) {
                T sum = 0;
                const T* bj = b+j*dimk;
                for (long k=0; k<dimk; k++) {
                    sum += ai[k]*bj[k];
                }
                ci[j] += sum;
            }
        }
    }
    
    /// Matrix * matrix reference implementation (slow but correct)
    template <typename T> 
    STATIC inline void mxm_reference(long dimi, long dimj, long dimk, 
                                     T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
        */
        
        for (long i=0; i<dimi; i++) {
            for (long k=0; k<dimk; k++) {
                for (long j=0; j<dimj; j++) {
                    c[i*dimj+j] += a[i*dimk+k]*b[k*dimj+j];
                }
            }
        }
    }
    
    /// Matrix * matrix (hand unrolled version)
    template <typename T> 
    STATIC inline void mxm(long dimi, long dimj, long dimk, 
                           T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          4-way unrolled k loop ... empirically fastest on PIII 
          compared to 2/3 way unrolling (though not by much).
        */
        
        long dimk4 = (dimk/4)*4;
        for (long i=0; i<dimi; i++, c+=dimj,a+=dimk) {
            const T* p = b;
            for (long k=0; k<dimk4; k+=4,p+=4*dimj) {
                T aik0 = a[k  ];
                T aik1 = a[k+1];
                T aik2 = a[k+2];
                T aik3 = a[k+3];
                const T* bk0 = p;
                const T* bk1 = bk0+dimj;
                const T* bk2 = bk1+dimj;
                const T* bk3 = bk2+dimj;
                for (long j=0; j<dimj; j++) {
                    c[j] += aik0*bk0[j] + aik1*bk1[j] + aik2*bk2[j] + aik3*bk3[j];
                }
            }
            for (long k=dimk4; k<dimk; k++) {
                T aik = a[k];
                for (long j=0; j<dimj; j++) {
                    c[j] += aik*b[k*dimj+j];
                }
            }
        }
    }
    
    /// Matrix transpose * matrix transpose reference implementation (slow but correct)
    template <typename T> 
    STATIC inline void mTxmT_reference(long dimi, long dimj, long dimk, 
                                       T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
        */
        
        for (long i=0; i<dimi; i++) {
            for (long j=0; j<dimj; j++) {
                for (long k=0; k<dimk; k++) {
                    c[i*dimj+j] += a[k*dimi+i]*b[j*dimk+k];
                }
            }
        }
    }
    
    /// Matrix transpose * matrix transpose (hand tiled and unrolled)
    template <typename T> 
    STATIC inline void mTxmT(long dimi, long dimj, long dimk, 
                             T* c, const T*a, const T*b) 
    {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          Tiled k, copy row of a into tempory, and uroll j once.
        */
        
        const int ktile=32;
        T ai[ktile];
        long dimj2 = (dimj/2)*2;
        
        T* csave = c;
        const T* asave = a;
        for (long klo=0; klo<dimk; klo+=ktile, asave+=ktile*dimi, b+=ktile) {
            long khi = klo+ktile;
            if (khi > dimk) khi = dimk;
            long nk = khi-klo;
            
            a = asave;
            c = csave;
            for (long i=0; i<dimi; i++,c+=dimj,a++) {
                const T* q = a;
                for (long k=0; k<nk; k++,q+=dimi) ai[k] = *q;
                
                const T* bj0 = b;
                for (long j=0; j<dimj2; j+=2,bj0+=2*dimk) {
                    const T* bj1 = bj0+dimk;
                    T sum0 = 0;
                    T sum1 = 0;
                    for (long k=0; k<nk; k++) {
                        sum0 += ai[k]*bj0[k];
                        sum1 += ai[k]*bj1[k];
                    }
                    c[j  ] += sum0;
                    c[j+1] += sum1;
                }
                
                for (long j=dimj2; j<dimj; j++,bj0+=dimk) {
                    T sum = 0;
                    for (long k=0; k<nk; k++) {
                        sum += ai[k]*bj0[k];
                    }
                    c[j] += sum;
                }
            }
        }
    }
    
    /// Inner product ... result(i,j,...,p,q,...) = sum(z) left(i,j,...,z)*right(z,p,q,...)
    
    /// By default it contracts the last dimension of the left tensor and
    /// the first dimension of the right tensor.  These defaults can be
    /// changed by specifying \c k0 and \c k1 , the index to contract in
    /// the left and right side tensors, respectively.  The defaults
    /// correspond to (\c k0=-1 and \c k1=0 ).
    template <class T>
    Tensor<T> inner(const Tensor<T>& left, const Tensor<T>& right,
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
        
        Tensor<T> result(nd,d);
        
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
    template <class T>
    void inner_result(const Tensor<T>& left, const Tensor<T>& right,
                      long k0, long k1, Tensor<T>& result) {
        
        
        // Need to include explicit optimizations for common special cases
        // E.g., contiguous, matrix-matrix, and 3d-tensor*matrix
        
        T* ptr = result.ptr();
        
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
            }
            else if (k0==(left.ndim-1) && k1==(right.ndim-1)) {
                // c[i,j] = a[i,k]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.size/dimk;
                long dimj = right.size/dimk;
                mxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
            else if (k0==0 && k1==(right.ndim-1)) {
                // c[i,j] = a[k,i]*b[j,k] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.stride[0];
                long dimj = right.size/dimk;
                mTxmT(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
            else if (k0==(left.ndim-1) && k1==0) {
                // c[i,j] = a[i,k]*b[k,j] ... collapsing extra indices to i & j
                long dimk = left.dim[k0];
                long dimi = left.size/dimk;
                long dimj = right.stride[0];
                mxm(dimi,dimj,dimk,ptr,left.ptr(),right.ptr());
                return;
            }
        }
        
        long dimj = left.dim[k0];
        TensorIterator<T> iter1=right.unary_iterator(1,false,false,k1);
        
        for (TensorIterator<T> iter0=left.unary_iterator(1,false,false,k0); 
             iter0._p0; ++iter0) {
            T* xp0 = iter0._p0;
            long s0 = iter0._s0;
            for (iter1.reset(); iter1._p0; ++iter1) {
                T* p0 = xp0;
                T* p1 = iter1._p0;
                long s1 = iter1._s0;
                T sum = 0;
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
    template <class T> 
    Tensor<T> transform(const Tensor<T>& t, const Tensor<T>& c) {
        TENSOR_ASSERT(c.ndim == 2,"second argument must be a matrix",c.ndim,&c);
        Tensor<T> result = t;
        for (long i=0; i<t.ndim; i++) {
            result = inner(result,c,0,0);
        }
        return result;
    }
    
    /// Restricted but heavily optimized form of transform()
    
    /// Performs the same operation as \c transform but with the
    /// restriction that the dimensions of \c c are identical, and it
    /// requires that the caller pass in workspace and a preallocated
    /// result, hoping that that both can be reused.  If the result &
    /// workspace are reused between calls, then no tensor constructors
    /// need be called and cache locality should be improved.  
    /// By passing in the workspace, this routine is
    /// kept thread safe.  All input tensors must be contiguous.  The
    /// workspace and the result must be of the same size as the input \c
    /// t .  The result tensor need not be initialized before calling
    /// fast_transform.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode 
    /// The input dimensions of \c t must all be the same and agree with the
    /// first dimension of \c c.
    template <class T> 
    void fast_transform(const Tensor<T>& t, const Tensor<T>& c, Tensor<T>& result, 
                        Tensor<T>& workspace) {
        TENSOR_ASSERT(c.ndim == 2,"second argument must be a matrix",c.ndim,&c);
        TENSOR_ASSERT(c.dim[0] == c.dim[1],"dimensions of c must match",0,&c);
        
        Tensor<T> *t0, *t1;
        int odd = t.ndim&1;
        if (odd) {
            t0 = &result;
            t1 = &workspace;
        }
        else {
            t0 = &workspace;
            t1 = &result;
        }
        t0->fill(0);
        inner_result(t,c,0,0,*t0);
        for (int n=1; n<t.ndim; n++) {
            t1->fill(0);
            inner_result(*t0,c,0,0,*t1);
            std::swap(t0,t1);
        }
        return;
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
    
#include "tensor_spec.h"
/*
    // Instantiations for double
    template class Tensor<double>;
    template class SliceTensor<double>;
    template std::ostream& operator << (std::ostream& s, const Tensor<double>& t);
    template Tensor<double> copy(const Tensor<double>& t);
    template Tensor<double> outer(const Tensor<double>& left, const Tensor<double>& right);
    template void inner_result(const Tensor<double>& left, const Tensor<double>& right,
                               long k0, long k1, Tensor<double>& result);
    template Tensor<double> inner(const Tensor<double>& left, const Tensor<double>& right,
                                  long k0, long k1);
    template Tensor<double> transform(const Tensor<double>& t, const Tensor<double>& c);
    template void fast_transform(const Tensor<double>& t, const Tensor<double>& c, Tensor<double>& result, Tensor<double>& workspace);
    template Tensor< Tensor<double>::scalar_type > abs(const Tensor<double>& t);
    template Tensor<double> transpose(const Tensor<double>& t);
    
    // Instantiations for float
    template class Tensor<float>;
    template class SliceTensor<float>;
    template std::ostream& operator << (std::ostream& s, const Tensor<float>& t);
    template Tensor<float> copy(const Tensor<float>& t);
    template Tensor<float> outer(const Tensor<float>& left, const Tensor<float>& right);
    template void inner_result(const Tensor<float>& left, const Tensor<float>& right,
                               long k0, long k1, Tensor<float>& result);
    template Tensor<float> inner(const Tensor<float>& left, const Tensor<float>& right,
                                 long k0, long k1);
    template Tensor<float> transform(const Tensor<float>& t, const Tensor<float>& c);
    template void fast_transform(const Tensor<float>& t, const Tensor<float>& c, Tensor<float>& result, Tensor<float>& workspace);
    template Tensor< Tensor<float>::scalar_type > abs(const Tensor<float>& t);
    template Tensor<float> transpose(const Tensor<float>& t);
    
    // Instantiations for long
    template class Tensor<long>;
    template class SliceTensor<long>;
    template std::ostream& operator << (std::ostream& s, const Tensor<long>& t);
    template Tensor<long> copy(const Tensor<long>& t);
    template Tensor<long> outer(const Tensor<long>& left, const Tensor<long>& right);
    template void inner_result(const Tensor<long>& left, const Tensor<long>& right,
                               long k0, long k1, Tensor<long>& result);
    template Tensor<long> inner(const Tensor<long>& left, const Tensor<long>& right,
                                long k0, long k1);
    template Tensor<long> transform(const Tensor<long>& t, const Tensor<long>& c);
    template void fast_transform(const Tensor<long>& t, const Tensor<long>& c, Tensor<long>& result, Tensor<long>& workspace);
    template Tensor< Tensor<long>::scalar_type > abs(const Tensor<long>& t);
    template Tensor<long> transpose(const Tensor<long>& t);
    
    // Instantiations for double_complex
    template class Tensor<double_complex>;
    template class SliceTensor<double_complex>;
    template std::ostream& operator << (std::ostream& s, const Tensor<double_complex>& t);
    template Tensor<double_complex> copy(const Tensor<double_complex>& t);
    template Tensor<double_complex> outer(const Tensor<double_complex>& left, const Tensor<double_complex>& right);
    template void inner_result(const Tensor<double_complex>& left, const Tensor<double_complex>& right,
                               long k0, long k1, Tensor<double_complex>& result);
    template Tensor<double_complex> inner(const Tensor<double_complex>& left, const Tensor<double_complex>& right,
                                          long k0, long k1);
    template Tensor<double_complex> transform(const Tensor<double_complex>& t, const Tensor<double_complex>& c);
    template void fast_transform(const Tensor<double_complex>& t, const Tensor<double_complex>& c, Tensor<double_complex>& result, Tensor<double_complex>& workspace);
    template Tensor< Tensor<double_complex>::scalar_type > abs(const Tensor<double_complex>& t);
    template Tensor<double_complex> transpose(const Tensor<double_complex>& t);
    
    // Instantiations for float_complex
    template class Tensor<float_complex>;
    template class SliceTensor<float_complex>;
    template std::ostream& operator << (std::ostream& s, const Tensor<float_complex>& t);
    template Tensor<float_complex> copy(const Tensor<float_complex>& t);
    template Tensor<float_complex> outer(const Tensor<float_complex>& left, const Tensor<float_complex>& right);
    template void inner_result(const Tensor<float_complex>& left, const Tensor<float_complex>& right,
                               long k0, long k1, Tensor<float_complex>& result);
    template Tensor<float_complex> inner(const Tensor<float_complex>& left, const Tensor<float_complex>& right,
                                         long k0, long k1);
    template Tensor<float_complex> transform(const Tensor<float_complex>& t, const Tensor<float_complex>& c);
    template void fast_transform(const Tensor<float_complex>& t, const Tensor<float_complex>& c, Tensor<float_complex>& result, Tensor<float_complex>& workspace);
    template Tensor< Tensor<float_complex>::scalar_type > abs(const Tensor<float_complex>& t);
    template Tensor<float_complex> transpose(const Tensor<float_complex>& t);
    
    // Instantiations only for complex types
    
    // Instantiations for double_complex
    template Tensor< Tensor<double_complex>::scalar_type > arg(const Tensor<double_complex>& t);
    template Tensor< Tensor<double_complex>::scalar_type > real(const Tensor<double_complex>& t);
    template Tensor< Tensor<double_complex>::scalar_type > imag(const Tensor<double_complex>& t);
    template Tensor<double_complex> conj(const Tensor<double_complex>& t);
    template Tensor<double_complex> conj_transpose(const Tensor<double_complex>& t);
    
    // Instantiations for float_complex
    template Tensor< Tensor<float_complex>::scalar_type > arg(const Tensor<float_complex>& t);
    template Tensor< Tensor<float_complex>::scalar_type > real(const Tensor<float_complex>& t);
    template Tensor< Tensor<float_complex>::scalar_type > imag(const Tensor<float_complex>& t);
    template Tensor<float_complex> conj(const Tensor<float_complex>& t);
    template Tensor<float_complex> conj_transpose(const Tensor<float_complex>& t);
*/
    
    template void mTxm(long dimi, long dimj, long dimk, double* c, const double* a, const double* b) ;
   
}
