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

#ifndef MADNESS_TENSOR_LOWRANKTENSOR_H_
#define MADNESS_TENSOR_LOWRANKTENSOR_H_

#include <memory>
#include <vector>

#include "../world/madness_exception.h"
#include "../world/print.h"
#include "gentensor.h"
#include "slice.h"
#include "srconf.h"
#include "tensor.h"
#include "tensortrain.h"
#include "type_data.h"

namespace madness {

template<typename T>
class SVDTensor : public SRConf<T> {
public:
    explicit SVDTensor() : SRConf<T>() {    }
    explicit SVDTensor(const Tensor<T>& rhs) : SRConf<T>(rhs) {    }
    explicit SVDTensor(const SVDTensor<T>& rhs) : SRConf<T>(static_cast<SRConf<T> >(rhs)) {    }
    explicit SVDTensor(const SRConf<T>& rhs) : SRConf<T>(SRConf<T>(rhs)) {    }
    explicit SVDTensor(const std::vector<long>& dims) : SRConf<T>(dims.size(),dims[0],TT_2D) {    }
    explicit SVDTensor(const Tensor<double>& weights, const Tensor<T>& vector1,
            const Tensor<T>& vector2, const unsigned int& dim,
            const unsigned int maxk) : SRConf<T>(weights,
                    vector1, vector2, dim, maxk ) {}

};

// forward declaration
template <class T> class SliceLowRankTensor;


template<typename T>
class LowRankTensor {

public:

    /// C++ typename of the real type associated with a complex type.
    typedef typename TensorTypeData<T>::scalar_type scalar_type;

    /// C++ typename of the floating point type associated with scalar real type
    typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

      struct implT {
        std::shared_ptr<Tensor<T> > full;
        std::shared_ptr<TensorTrain<T> > tt;
        std::shared_ptr<SVDTensor<T> > svd;

        /// ctor resets all shared_ptr
        implT() : full(), tt(), svd() {}

        /// make sure only a single pointer is set
        void check_unique() const {
            MADNESS_ASSERT(full.get() xor tt.get() xor svd.get());
        }
    };

    TensorType type;
    implT impl;

    /// empty ctor
    LowRankTensor() : type(TT_NONE) {}

    /// copy ctor, shallow
    LowRankTensor(const LowRankTensor<T>& other) : type(other.type),
            impl(other.impl) {
    }

    /// ctor with dimensions; constructs tensor filled with zeros
    LowRankTensor(const std::vector<long>& dim, const TensorType& tt)
        : type(tt) {

        if (type==TT_FULL) impl.full=std::shared_ptr<Tensor<T> >(new Tensor<T>(dim));
        if (type==TT_2D) impl.svd=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(dim));
        if (type==TT_TENSORTRAIN) impl.tt=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(dim));
    }

    /// ctor with dimensions; constructs tensor filled with zeros
    LowRankTensor(const std::vector<long>& dim, const TensorArgs& targs) :
        LowRankTensor(dim, targs.tt) {
    }

    /// ctor with dimensions; all dims have the same value k
    LowRankTensor(const TensorType& tt, const long k, const long ndim) :
        LowRankTensor(std::vector<long>(ndim,k), tt) {
    }

    /// ctor with a regular Tensor and arguments, deep
    LowRankTensor(const Tensor<T>& rhs, const double& thresh, const TensorType& tt) :
        LowRankTensor(rhs,TensorArgs(thresh,tt)) {
    }

    /// ctor with a regular Tensor and arguments, deep
    LowRankTensor(const Tensor<T>& rhs, const TensorArgs& targs)
        : type(targs.tt) {

        // deep copy, either explicitly (full) or implicitly (SVD, TensorTrain)
        if (type==TT_FULL) {
            impl.full=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(rhs)));

        } else if ((type==TT_2D) or (type==TT_TENSORTRAIN)) {

            // construct tensortrain first, convert into SVD format afterwards
            std::shared_ptr<TensorTrain<T> > tt(new TensorTrain<T>(rhs,targs.thresh*facReduce()));

            if (type==TT_2D) {
                Tensor<T> U,VT;
                Tensor< typename Tensor<T>::scalar_type > s;
                tt->two_mode_representation(U,VT,s);
                const long r=VT.dim(0);
                const long nd=VT.ndim();
                MADNESS_ASSERT(U.dim(nd-1)==r);

                // empty SVDTensor due to zero rank
                if (r==0) {
                    impl.svd.reset(new SVDTensor<T>(tt->dims()));

                } else {
                    Tensor<T> UU=U.reshape(U.size()/r,r);
                    SVDTensor<T> sr(s, copy(transpose(UU)), VT.reshape(r,VT.size()/r),
                                          rhs.ndim(), rhs.dim(0));
                    sr.normalize();
                    impl.svd.reset(new SVDTensor<T>(sr));
                }
            } else if (type==TT_TENSORTRAIN) {
                impl.tt=tt;
            }
        }
    }

    /// ctor with a regular Tensor, deep
    LowRankTensor(const Tensor<T>& other) : type(TT_FULL), impl() {
        impl.full=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(other)));
    }

    /// ctor with a TensorTrain as argument, shallow
    LowRankTensor(const TensorTrain<T>& tt) : type(TT_TENSORTRAIN), impl() {
        impl.tt=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(tt));
    }

    /// ctor with a SVDTensor as argument, shallow
    LowRankTensor(const SVDTensor<T>& t) : type(TT_2D), impl() {
        impl.svd=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(t));
    }

    /// ctor with a SliceLowRankTensor as argument, deep
    LowRankTensor(const SliceLowRankTensor<T>& other) : type(other.lrt.type), impl() {
        *this=other;
    }

    /// shallow assignment operator
    LowRankTensor& operator=(const LowRankTensor<T>& other) {
        if (this!=&other) {
            type=other.type;
            impl=other.impl;
        }
        return *this;
    }

    /// deep assignment with slices: g0 = g1(s)
    LowRankTensor& operator=(const SliceLowRankTensor<T>& other) {
        type=other.lrt.type;
        const std::vector<Slice>& s=other.thisslice;
        if (this->type==TT_FULL)
            impl.full.reset(new Tensor<T>(copy((*other.lrt.impl.full)(s))));
        else if (this->type==TT_2D)
            impl.svd.reset(new SVDTensor<T>(other.lrt.impl.svd->copy_slice(s)));
        else if (this->type==TT_TENSORTRAIN)
            impl.tt.reset(new TensorTrain<T>(copy(*(other.lrt.impl.tt),s)));
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    /// Type conversion makes a deep copy
    template <class Q> operator LowRankTensor<Q>() const { // type conv => deep copy

        LowRankTensor<Q> result;
        if (this->type==TT_FULL) {
            result.impl.full.reset(new Tensor<Q>(*impl.full));
        } else if (this->type==TT_2D) {
            MADNESS_EXCEPTION("no type conversion for TT_2D yes=t",1);
        } else if (this->type==TT_TENSORTRAIN) {
            result.impl.tt.reset(new TensorTrain<Q>(*impl.tt));
        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return result;
    }


    /// general slicing, shallow; for temporary use only!
    SliceLowRankTensor<T> operator()(const std::vector<Slice>& s) {
        return SliceLowRankTensor<T>(*this,s);
    }

    /// general slicing, shallow; for temporary use only!
    const SliceLowRankTensor<T> operator()(const std::vector<Slice>& s) const {
        return SliceLowRankTensor<T>(*this,s);
    }



//    /// shallow assignment operator
//    LowRankTensor& operator=(const double& fac) {
//        MADNESS_EXCEPTION("no assignment of numbers in LowRankTensor",1);
//        return *this;
//    }

    /// deep copy
    friend LowRankTensor copy(const LowRankTensor& other) {
        LowRankTensor result;
        result.type=other.type;
        if (other.type==TT_FULL)
            result.impl.full=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(*other.impl.full)));
        if (other.type==TT_2D)
            result.impl.svd=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(copy(*other.impl.svd)));
        if (other.type==TT_TENSORTRAIN)
            result.impl.tt=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(copy(*other.impl.tt)));
        return result;
    }

    /// return the tensor type
    TensorType tensor_type() const {return type;}

    long ndim() const {
        if (type==TT_NONE) return -1;
        else if (type==TT_FULL) return impl.full->ndim();
        else if (type==TT_2D) return impl.svd->dim();
        else if (type==TT_TENSORTRAIN) return impl.tt->ndim();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return -1;
    }

    /// return the number of entries in dimension i
    long dim(const int i) const {
        if (type==TT_NONE) return -1;
        else if (type==TT_FULL) return impl.full->dim(i);
        else if (type==TT_2D) return impl.svd->get_k();
        else if (type==TT_TENSORTRAIN) return impl.tt->dim(i);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return -1;
    }


    void normalize() {
        if (type==TT_2D) impl.svd->normalize();
    }

    float_scalar_type normf() const {
        if (type==TT_FULL) return impl.full->normf();
        else if (type==TT_2D) return impl.svd->normf();
        else if (type==TT_TENSORTRAIN) return impl.tt->normf();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return 0;
    }

    float_scalar_type svd_normf() const {
        if (type==TT_FULL) return impl.full->normf();
        else if (type==TT_2D) return impl.svd->svd_normf();
        else if (type==TT_TENSORTRAIN) return impl.tt->normf();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return 0;
    }


    /// Inplace multiplication by scalar of supported type (legacy name)

    /// @param[in] x Scalar value
    /// @return %Reference to this tensor
    template <typename Q>
    typename IsSupported<TensorTypeData<Q>,LowRankTensor<T>&>::type
    scale(Q fac) {
        if (type==TT_NONE) return *this;
        else if (type==TT_FULL) impl.full->scale(T(fac));
        else if (type==TT_2D) impl.svd->scale(T(fac));
        else if (type==TT_TENSORTRAIN) impl.tt->scale(T(fac));
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    Tensor<T> full_tensor_copy() const {
        if (type==TT_NONE) return Tensor<T>();
        else if (type==TT_FULL) return copy(*impl.full);
        else if (type==TT_2D) return impl.svd->reconstruct();
        else if (type==TT_TENSORTRAIN) return impl.tt->reconstruct();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return Tensor<T>();
    }

    const Tensor<T>& full_tensor() const {
        MADNESS_ASSERT(type==TT_FULL);
        return *impl.full;
    }

    Tensor<T>& full_tensor() {
        MADNESS_ASSERT(type==TT_FULL);
        return *impl.full;
    }


    /// reconstruct this to return a full tensor
    Tensor<T> reconstruct_tensor() const {

        if (type==TT_FULL) return full_tensor();
        else if (type==TT_2D or type==TT_TENSORTRAIN) return full_tensor_copy();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return Tensor<T>();
    }


    static double facReduce() {return 1.e-3;}
    static double fac_reduce() {return 1.e-3;}

    long rank() const {
        if (type==TT_NONE) return 0;
        else if (type==TT_FULL) return -1;
        else if (type==TT_2D) return impl.svd->rank();
        else if (type==TT_TENSORTRAIN) return -2;   // not defined;
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return 0l;
    }

    bool has_data() const {return (type != TT_NONE);}

    bool has_no_data() const {return (not has_data());}

    long size() const {
        if (type==TT_NONE) return 0l;
        else if (type==TT_FULL) return impl.full->size();
        else if (type==TT_2D) return impl.svd->nCoeff();
        else if (type==TT_TENSORTRAIN) return impl.tt->size();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return false;
    }

    long real_size() const {
        if (type==TT_NONE) return 0l;
        else if (type==TT_FULL) return impl.full->size();
        else if (type==TT_2D) return impl.svd->real_size();
        else if (type==TT_TENSORTRAIN) return impl.tt->real_size();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return false;
    }

    /// returns the trace of <this|rhs>
    template<typename Q>
    TENSOR_RESULT_TYPE(T,Q) trace_conj(const LowRankTensor<Q>& rhs) const {

        if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in LowRankTensor, sorry",1);
        if (TensorTypeData<Q>::iscomplex) MADNESS_EXCEPTION("no complex trace in LowRankTensor, sorry",1);

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        // fast return if possible
        if ((this->rank()==0) or (rhs.rank()==0)) return resultT(0.0);

        MADNESS_ASSERT(this->type==rhs.type);

        if (type==TT_NONE) return resultT(0.0);
        else if (type==TT_FULL) return impl.full->trace_conj(*(rhs.impl.full));
        else if (type==TT_2D) return overlap(*impl.svd,*(rhs.impl.svd));
        else if (type==TT_TENSORTRAIN) return impl.tt->trace(*rhs.impl.tt);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return TENSOR_RESULT_TYPE(T,Q)(0);
    }

    LowRankTensor& operator+=(const LowRankTensor& other) {
        if (type==TT_FULL) impl.full->operator+=(*other.impl.full);
        else if (type==TT_2D) impl.svd->append((*other.impl.svd),1.0);
        else if (type==TT_TENSORTRAIN) impl.tt->operator+=(*other.impl.tt);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    LowRankTensor& operator-=(const LowRankTensor& other) {
        if (type==TT_FULL) impl.full->operator-=(*other.impl.full);
        else if (type==TT_2D) impl.svd->append((*other.impl.svd),-1.0);
        else if (type==TT_TENSORTRAIN) impl.tt->operator-=(*other.impl.tt);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    LowRankTensor& gaxpy(const T alpha, const LowRankTensor& other, const T beta) {
        MADNESS_ASSERT(this->type==other.type);

        if (type==TT_FULL) {
            impl.full->gaxpy(alpha,*other.impl.full,beta);
        } else if (type==TT_2D) {
            if (not (alpha==1.0)) this->scale(alpha);
            impl.svd->append(*other.impl.svd,beta);
        } else if (type==TT_TENSORTRAIN) {
            impl.tt->gaxpy(alpha,*other.impl.tt,beta);
        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    void add_SVD(const LowRankTensor& other, const double& thresh) {
        if (type==TT_FULL) impl.full->operator+=(*other.impl.full);
        else if (type==TT_2D) impl.svd->add_SVD((*other.impl.svd),thresh*facReduce());
        else if (type==TT_TENSORTRAIN) impl.tt->operator+=(*other.impl.tt);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
    }

    /// multiply with a number
    template<typename Q>
    LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> operator*(const Q& x) const {
        LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> result(copy(*this));
        result.scale(x);
        return result;
    }

    void reduce_rank(const double& thresh) {
        if ((type==TT_FULL) or (type==TT_NONE)) return;
        else if (type==TT_2D) impl.svd->divide_and_conquer_reduce(thresh*facReduce());
        else if (type==TT_TENSORTRAIN) impl.tt->truncate(thresh*facReduce());
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
    }

    /// Returns a pointer to the internal data

    /// @param[in]  ivec    index of core vector to which the return values points
    T* ptr(const int ivec=0) {
        if (type==TT_NONE) return 0;
        else if (type==TT_FULL) return impl.full->ptr();
        else if (type==TT_2D) return impl.svd->ref_vector(ivec).ptr();
        else if (type==TT_TENSORTRAIN) return impl.tt->ptr(ivec);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return 0;
    }

    /// Returns a pointer to the internal data

    /// @param[in]  ivec    index of core vector to which the return values points
    const T* ptr(const int ivec=0) const {
        if (type==TT_NONE) return 0;
        else if (type==TT_FULL) return impl.full->ptr();
        else if (type==TT_2D) return impl.svd->ref_vector(ivec).ptr();
        else if (type==TT_TENSORTRAIN) return impl.tt->ptr(ivec);
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return 0;
    }


    SRConf<T>& config() {
        MADNESS_EXCEPTION("implement config in LowRankTensor",1);
    }

    void method1() {
        if (type == TT_FULL) impl.full->method1();
        if (type == TT_2D) impl.svd->method1();
        if (type == TT_TENSORTRAIN) impl.tt->method1();
        // etc.
    }


    /// Transform all dimensions of the tensor t by the matrix c

    /// \ingroup tensor
    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
    /// The input dimensions of \c t must all be the same and agree with
    /// the first dimension of \c c .  The dimensions of \c c may differ in
    /// size.
    template <typename Q>
    LowRankTensor< TENSOR_RESULT_TYPE(T,Q) > transform(const Tensor<Q>& c) const {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        if (type==TT_NONE) return LowRankTensor<resultT>();
        else if (type==TT_FULL) {
            Tensor<resultT> result=madness::transform((*impl.full),c);
            return LowRankTensor<resultT>(result,0.0,TT_FULL);
        } else if (type==TT_2D) {
            SVDTensor<resultT> result(impl.svd->transform(c));
            return LowRankTensor<resultT>(result);
        } else if (type==TT_TENSORTRAIN) {
            TensorTrain<resultT> tt=madness::transform(*impl.tt,c);
            return LowRankTensor(tt);
        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        if (has_no_data()) return LowRankTensor<resultT>();
        return LowRankTensor<resultT>();
    }


    /// Transform all dimensions of the tensor t by distinct matrices c

    /// \ingroup tensor
    /// Similar to transform but each dimension is transformed with a
    /// distinct matrix.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c[0](i',i) c[1](j',j) c[2](k',k) ...
    /// \endcode
    /// The first dimension of the matrices c must match the corresponding
    /// dimension of t.
    template <typename Q>
    LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> general_transform(const Tensor<Q> c[]) const {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        if (type==TT_NONE) return LowRankTensor<resultT>();
        else if (type==TT_FULL) {
            Tensor<resultT> result=madness::general_transform((*impl.full),c);
            return LowRankTensor<resultT>(result,0.0,TT_FULL);
        } else if (type==TT_2D) {
            SVDTensor<resultT> result(impl.svd->general_transform(c));
            return LowRankTensor<resultT>(result);
        } else if (type==TT_TENSORTRAIN) {
            TensorTrain<resultT> tt=madness::general_transform((*impl.tt),c);
            return LowRankTensor(tt);
        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        if (has_no_data()) return LowRankTensor<resultT>();
        return LowRankTensor<resultT>();
    }


    /// Transforms one dimension of the tensor t by the matrix c, returns new contiguous tensor

    /// \ingroup tensor
    /// \code
    /// transform_dir(t,c,1) = r(i,j,k,...) = sum(j') t(i,j',k,...) * c(j',j)
    /// \endcode
    /// @param[in] t Tensor to transform (size of dimension to be transformed must match size of first dimension of \c c )
    /// @param[in] c Matrix used for the transformation
    /// @param[in] axis Dimension (or axis) to be transformed
    /// @result Returns a new, contiguous tensor
    template <typename Q>
    LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> transform_dir(const Tensor<Q>& c, int axis) const {
        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        if (type==TT_NONE) return LowRankTensor<resultT>();
        else if (type==TT_FULL) {
            Tensor<resultT> result=madness::transform_dir((*impl.full),c,axis);
            return LowRankTensor<resultT>(result,0.0,TT_FULL);
        } else if (type==TT_2D) {
            SVDTensor<resultT> result(impl.svd->transform_dir(c,axis));
            return LowRankTensor<resultT>(result);
        } else if (type==TT_TENSORTRAIN) {
            TensorTrain<resultT> tt=madness::transform_dir((*impl.tt),c,axis);
            return LowRankTensor(tt);
        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        if (has_no_data()) return LowRankTensor<resultT>();
        return LowRankTensor<resultT>();
    }

};


template <class T, class Q>
LowRankTensor< TENSOR_RESULT_TYPE(T,Q) > transform(const LowRankTensor<Q>& t, const Tensor<T>& c) {
    return t.transform(c);
}

template <class T, class Q>
LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> general_transform(const LowRankTensor<T>& t, const Tensor<Q> c[]) {
    return t.general_transform(c);
}

/// Transforms one dimension of the tensor t by the matrix c, returns new contiguous tensor

/// \ingroup tensor
/// \code
/// transform_dir(t,c,1) = r(i,j,k,...) = sum(j') t(i,j',k,...) * c(j',j)
/// \endcode
/// @param[in] t Tensor to transform (size of dim to be transformed must match size of first dim of \c c )
/// @param[in] c Matrix used for the transformation
/// @param[in] axis Dimension (or axis) to be transformed
/// @result Returns a new, contiguous tensor
template <class T, class Q>
LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> transform_dir(const LowRankTensor<Q>& t,
        const Tensor<T>& c, int axis) {
    return t.transform_dir(c,axis);
}

/// outer product of two Tensors, yielding a low rank tensor

/// do the outer product of two tensors; distinguish these tensortype cases by
/// the use of final_tensor_type
///  - full x full -> full
///  - full x full -> SVD                           ( default )
///  - TensorTrain x TensorTrain -> TensorTrain
/// all other combinations are currently invalid.
template <class T, class Q>
LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const LowRankTensor<T>& t1,
        const LowRankTensor<Q>& t2, const TensorType final_tensor_type=TT_2D) {

    typedef TENSOR_RESULT_TYPE(T,Q) resultT;


    MADNESS_ASSERT(t1.tensor_type()==t2.tensor_type());

    if (final_tensor_type==TT_FULL) {
        MADNESS_ASSERT(t1.tensor_type()==TT_FULL);
        Tensor<resultT> t(outer(*t1.impl.full,*t2.impl.full));
        return LowRankTensor<resultT>(t);

    } else if (final_tensor_type==TT_2D) {
        MADNESS_ASSERT(t1.tensor_type()==TT_FULL);

        // srconf is shallow, do deep copy here
        const Tensor<T> lhs=t1.full_tensor_copy();
        const Tensor<Q> rhs=t2.full_tensor_copy();

        const long k=lhs.dim(0);
        const long dim=lhs.ndim()+rhs.ndim();
        long size=1;
        for (int i=0; i<lhs.ndim(); ++i) size*=k;
        MADNESS_ASSERT(size==lhs.size());
        MADNESS_ASSERT(size==rhs.size());
        MADNESS_ASSERT(lhs.size()==rhs.size());

        Tensor<double> weights(1);
        weights=1.0;

        SRConf<resultT> srconf(weights,lhs.reshape(1,lhs.size()),rhs.reshape(1,rhs.size()),dim,k);
//        srconf.normalize();
        return LowRankTensor<resultT>(SVDTensor<resultT>(srconf));

    } else if (final_tensor_type==TT_TENSORTRAIN) {
        MADNESS_EXCEPTION("implement outer for TensorTrain",1);
    } else {
        MADNESS_EXCEPTION("you should not be here",1);
    }
    return LowRankTensor<TENSOR_RESULT_TYPE(T,Q)>();

}


namespace archive {
    /// Serialize a tensor
    template <class Archive, typename T>
    struct ArchiveStoreImpl< Archive, LowRankTensor<T> > {

        friend class LowRankTensor<T>;
        /// Stores the GenTensor to an archive
        static void store(const Archive& ar, const LowRankTensor<T>& t) {
            bool exist=t.has_data();
            int i=int(t.type);
            ar & exist & i;
            if (exist) {
                if (t.impl.svd) ar & *t.impl.svd.get();
                if (t.impl.full) ar & *t.impl.full.get();
                if (t.impl.tt) ar & *t.impl.tt.get();
            }
        };
    };


    /// Deserialize a tensor ... existing tensor is replaced
    template <class Archive, typename T>
    struct ArchiveLoadImpl< Archive, LowRankTensor<T> > {

        friend class GenTensor<T>;
        /// Replaces this GenTensor with one loaded from an archive
        static void load(const Archive& ar, LowRankTensor<T>& t) {
            // check for pointer existence
            bool exist=false;
            int i=-1;
            ar & exist & i;
            t.type=TensorType(i);

            if (exist) {
                if (t.type==TT_2D) {
                    SVDTensor<T> svd;
                    ar & svd;
                    t.impl.svd.reset(new SVDTensor<T>(svd));
                } else if (t.type==TT_FULL) {
                    Tensor<T> full;
                    ar & full;
                    t.impl.full.reset(new Tensor<T>(full));
                } else if (t.type==TT_TENSORTRAIN) {
                    TensorTrain<T> tt;
                    ar & tt;
                    t.impl.tt.reset(new TensorTrain<T>(tt));
                }

            }
        };
    };
};


/// The class defines tensor op scalar ... here define scalar op tensor.
template <typename T, typename Q>
typename IsSupported < TensorTypeData<Q>, LowRankTensor<T> >::type
operator*(const Q& x, const LowRankTensor<T>& t) {
    return t*x;
}

/// add all the GenTensors of a given list

 /// If there are many tensors to add it's beneficial to do a sorted addition and start with
 /// those tensors with low ranks
 /// @param[in]  addends     a list with gentensors of same dimensions; will be destroyed upon return
 /// @param[in]  eps         the accuracy threshold
 /// @param[in]  are_optimal flag if the GenTensors in the list are already in SVD format (if TT_2D)
 /// @return     the sum GenTensor of the input GenTensors
 template<typename T>
 LowRankTensor<T> reduce(std::list<LowRankTensor<T> >& addends, double eps, bool are_optimal=false) {
     typedef typename std::list<LowRankTensor<T> >::iterator iterT;
     LowRankTensor<T> result=copy(addends.front());
      for (iterT it=++addends.begin(); it!=addends.end(); ++it) {
          result+=*it;
      }
      result.reduce_rank(eps);
      return result;

 }


/// implements a temporary(!) slice of a LowRankTensor
template<typename T>
class SliceLowRankTensor {
public:

    LowRankTensor<T> lrt;
    std::vector<Slice> thisslice;

    // all ctors are private, only accessible by GenTensor

    /// default ctor
    SliceLowRankTensor<T> () {}

    /// ctor with a GenTensor; shallow
    SliceLowRankTensor<T> (const LowRankTensor<T>& gt, const std::vector<Slice>& s)
            : lrt(gt), thisslice(s) {}

public:

    /// assignment as in g(s) = g1;
    SliceLowRankTensor<T>& operator=(const LowRankTensor<T>& rhs) {
        print("You don't want to assign to a SliceLowRankTensor; use operator+= instead");
        MADNESS_ASSERT(0);
        return *this;
    };

    /// assignment as in g(s) = g1(s);
    SliceLowRankTensor<T>& operator=(const SliceLowRankTensor<T>& rhs) {
        print("You don't want to assign to a SliceLowRankTensor; use operator+= instead");
        MADNESS_ASSERT(0);
        return *this;
    };

    /// inplace addition as in g(s)+=g1
    SliceLowRankTensor<T>& operator+=(const LowRankTensor<T>& rhs) {

        // fast return if possible
        if (rhs.has_no_data() or rhs.rank()==0) return *this;

        if (lrt.has_data()) MADNESS_ASSERT(lrt.tensor_type()==rhs.tensor_type());

        // no fast return possible!!!
        //          if (this->rank()==0) {
        //              // this is a deep copy
        //              *this=rhs(rhs_s);
        //              this->scale(beta);
        //              return;
        //          }

        std::vector<Slice> rhs_slice(rhs.ndim(),Slice(_));

        if (lrt.type==TT_FULL) {
            (*lrt.impl.full)(thisslice).gaxpy(1.0,(*rhs.impl.full)(rhs_slice),1.0);

        } else if (lrt.type==TT_2D) {
            MADNESS_ASSERT(lrt.impl.svd->has_structure());
            MADNESS_ASSERT(rhs.impl.svd->has_structure());
            lrt.impl.svd->inplace_add(*rhs.impl.svd,thisslice,rhs_slice, 1.0, 1.0);

        } else if (lrt.type==TT_TENSORTRAIN) {
            lrt.impl.tt->gaxpy(thisslice,*rhs.impl.tt,1.0,rhs_slice);
        }
        return *this;
    }

    /// inplace subtraction as in g(s)-=g1
    SliceLowRankTensor<T>& operator-=(const LowRankTensor<T>& rhs) {

        // fast return if possible
        if (rhs.has_no_data() or rhs.rank()==0) return *this;

        if (lrt.has_data()) MADNESS_ASSERT(lrt.tensor_type()==rhs.tensor_type());

        // no fast return possible!!!
        //          if (lrt.rank()==0) {
        //              // this is a deep copy
        //              *this=rhs(rhs_s);
        //              lrt.scale(beta);
        //              return;
        //          }

        std::vector<Slice> rhs_slice(rhs.ndim(),Slice(_));

        if (lrt.type==TT_FULL) {
            (*lrt.impl.full)(thisslice).gaxpy(1.0,(*rhs.impl.full)(rhs_slice),-1.0);

        } else if (lrt.type==TT_2D) {
            MADNESS_ASSERT(lrt.impl.svd->has_structure());
            MADNESS_ASSERT(rhs.impl.svd->has_structure());
            lrt.impl.svd->inplace_add(*rhs.impl.svd,thisslice,rhs_slice, 1.0, -1.0);

        } else if (lrt.type==TT_TENSORTRAIN) {
            lrt.impl.tt->gaxpy(thisslice,*rhs.impl.tt,-1.0,rhs_slice);
        }
        return *this;
    }

    /// inplace addition as in g(s)+=g1(s)
    SliceLowRankTensor<T>& operator+=(const SliceLowRankTensor<T>& rhs) {
        // fast return if possible
        if (rhs.lrt.has_no_data() or rhs.lrt.rank()==0) return *this;

        if (lrt.has_data()) MADNESS_ASSERT(lrt.tensor_type()==rhs.lrt.tensor_type());

        // no fast return possible!!!
        //          if (lrt.rank()==0) {
        //              // this is a deep copy
        //              *this=rhs(rhs_s);
        //              lrt.scale(beta);
        //              return;
        //          }

        std::vector<Slice> rhs_slice=rhs.thisslice;

        if (lrt.type==TT_FULL) {
            (*lrt.impl.full)(thisslice).gaxpy(1.0,(*rhs.lrt.impl.full)(rhs_slice),1.0);

        } else if (lrt.type==TT_2D) {
            MADNESS_ASSERT(lrt.impl.svd->has_structure());
            MADNESS_ASSERT(rhs.lrt.impl.svd->has_structure());
            lrt.impl.svd->inplace_add(*rhs.lrt.impl.svd,thisslice,rhs_slice, 1.0, 1.0);

        } else if (lrt.type==TT_TENSORTRAIN) {
            lrt.impl.tt->gaxpy(thisslice,*rhs.lrt.impl.tt,1.0,rhs_slice);
        }
        return *this;

    }

    /// inplace zero-ing as in g(s)=0.0
    SliceLowRankTensor<T>& operator=(const T& number) {
        MADNESS_ASSERT(number==T(0.0));

        if (lrt.type==TT_FULL) {
            (*lrt.impl.full)(thisslice)=0.0;

        } else if (lrt.type==TT_2D) {
            MADNESS_ASSERT(lrt.impl.svd->has_structure());
            LowRankTensor<T> tmp(lrt);
            lrt.impl.svd->inplace_add(*tmp.impl.svd,thisslice,thisslice, 1.0, -1.0);

        } else if (lrt.type==TT_TENSORTRAIN) {
            lrt.impl.tt->gaxpy(thisslice,*lrt.impl.tt,-1.0,thisslice);
        }
        return *this;
    }

    friend LowRankTensor<T> copy(const SliceLowRankTensor<T>& other) {
        LowRankTensor<T> result;
        const std::vector<Slice> s=other.thisslice;
        result.type=other.lrt.type;
        if (result.type==TT_FULL)
            result.impl.full.reset(new Tensor<T>(copy((*other.lrt.impl.full)(s))));
        else if (result.type==TT_2D)
            result.impl.svd.reset(new SVDTensor<T>(other.lrt.impl.svd->copy_slice(s)));
        else if (result.type==TT_TENSORTRAIN)
            result.impl.tt.reset(new TensorTrain<T>(copy(*(other.lrt.impl.tt),s)));
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return result;

    }


};



} // namespace madness

#endif /* MADNESS_TENSOR_LOWRANKTENSOR_H_ */
