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
#include <variant>

#include <madness/world/madness_exception.h>
#include <madness/world/print.h>
#include <madness/tensor/slice.h>
#include <madness/tensor/SVDTensor.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensortrain.h>
#include "type_data.h"
#include <madness/tensor/RandomizedMatrixDecomposition.h>


namespace madness {


// forward declaration
template <class T> class SliceLowRankTensor;


template<typename T>
class GenTensor {

public:

	friend class SliceLowRankTensor<T>;

    /// C++ typename of the real type associated with a complex type.
    typedef typename TensorTypeData<T>::scalar_type scalar_type;

    /// C++ typename of the floating point type associated with scalar real type
    typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

    /// empty ctor
	GenTensor() = default;

    /// copy ctor, shallow
	GenTensor(const GenTensor<T>& other) = default;

	GenTensor(const long ndim, const long* dims, const TensorType& tt) {
		if (tt==TT_FULL) tensor=std::shared_ptr<Tensor<T> >(new Tensor<T>(ndim, dims));
		if (tt==TT_2D) tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(ndim, dims));
		if (tt==TT_TENSORTRAIN) tensor=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(ndim, dims));
    }

    /// ctor with dimensions; constructs tensor filled with zeros
	GenTensor(const std::vector<long>& dim, const TensorType& tt) :
		GenTensor(long(dim.size()),&dim.front(),tt) {
    }

    /// ctor with dimensions; constructs tensor filled with zeros
	GenTensor(const std::vector<long>& dim, const TensorArgs& targs) :
		GenTensor(dim, targs.tt) {
    }

    /// ctor with dimensions; all dims have the same value k
	GenTensor(const TensorType& tt, const long k, const long ndim) :
		GenTensor(std::vector<long>(ndim,k), tt) {
    }

    /// ctor with a regular Tensor and arguments, deep
	GenTensor(const Tensor<T>& rhs, const double& thresh, const TensorType& tt) :
		GenTensor(rhs,TensorArgs(thresh,tt)) {
    }

    /// ctor with a regular Tensor and arguments, deep
	GenTensor(const Tensor<T>& rhs, const TensorArgs& targs) {
		if (targs.tt==TT_FULL) *this=copy(rhs);
		else if (targs.tt==TT_2D) {
			if (rhs.size()==0) {
				tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(rhs,targs.thresh*facReduce()));
        } else {
				TensorTrain<T> tt(rhs,targs.thresh*facReduce());
				GenTensor<T> tmp=tt;
				*this=tmp.convert(targs);
			}
//		} else if (targs.tt==TT_DYNAMIC) {
//			if (rhs.size()==0) {
//				tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(rhs,targs.thresh*facReduce()));
//			} else {
//
//				long maxrank=std::max(50.0,floor(0.3*sqrt(rhs.size())));
//				RandomizedMatrixDecomposition<T> rmd=RMDFactory().maxrank(maxrank);
//				Tensor<T> Q=rmd.compute_range(rhs,targs.thresh*facReduce()*0.1,{0,0});
//				if (Q.size()==0) {
//					*this=SVDTensor<T>(rhs.ndim(),rhs.dims());
//				} else if (not rmd.exceeds_maxrank()) {
//					SVDTensor<T> result(rhs.ndim(),rhs.dims());
//					result=SVDTensor<T>::compute_svd_from_range(Q,rhs);
//					*this=result;
//				} else {
//					*this=copy(rhs);
////					TensorTrain<T> tt(rhs,targs.thresh*facReduce());
////					GenTensor<T> tmp=tt;
////					*this=tmp.convert(targs);
//				}
////				tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(rhs,targs.thresh*facReduce()));
//			}
		} else if (targs.tt==TT_TENSORTRAIN) {
			tensor=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(rhs,targs.thresh*facReduce()));
                    } else {
			MADNESS_EXCEPTION("unknown tensor type in LowRankTensor constructor",1);
        }
    }

    /// ctor with a regular Tensor, deep
	GenTensor(const Tensor<T>& other)  {
		tensor=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(other)));
    }

    /// ctor with a TensorTrain as argument, shallow
	GenTensor(const TensorTrain<T>& other) {
		tensor=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(copy(other))) ;
    }

    /// ctor with a SVDTensor as argument, shallow
	GenTensor(const SVDTensor<T>& other) {
		tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(copy(other)));
    }

    /// ctor with a SliceLowRankTensor as argument, deep
	GenTensor(const SliceLowRankTensor<T>& other) {
        *this=other;
    }

    /// shallow assignment operator
	GenTensor& operator=(const GenTensor<T>& other) {
		if (this!=&other) tensor=other.tensor;
		return *this;
	}

	/// deep assignment operator
	GenTensor& operator=(const Tensor<T>& other) {
		tensor=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(other)));
		return *this;
        }

	/// deep assignment operator
	GenTensor& operator=(const SVDTensor<T>& other) {
		tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(copy(other)));
		return *this;
	}

	/// deep assignment operator
	GenTensor& operator=(const TensorTrain<T>& other) {
		tensor=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(copy(other)));
        return *this;
    }

    /// deep assignment with slices: g0 = g1(s)
	GenTensor& operator=(const SliceLowRankTensor<T>& other) {
		const std::array<Slice,TENSOR_MAXDIM>& s=other.thisslice;
		MADNESS_ASSERT(other.is_assigned());
		if (other.is_full_tensor())
			tensor=std::shared_ptr<Tensor<T> >(new Tensor<T>(copy(other.get_tensor()(s))));
		else if (other.is_svd_tensor())
			tensor=std::shared_ptr<SVDTensor<T> >(new SVDTensor<T>(other.get_svdtensor().copy_slice(s)));
		else if (other.is_tensortrain())
			tensor=std::shared_ptr<TensorTrain<T> >(new TensorTrain<T>(copy(other.get_tensortrain(),s)));
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

    /// Type conversion makes a deep copy
	template <class Q> operator GenTensor<Q>() const { // type conv => deep copy

		GenTensor<Q> result;
		if (is_full_tensor()) {
			result=Tensor<Q>(get_tensor());
		} else if (is_svd_tensor()) {
            MADNESS_EXCEPTION("no type conversion for TT_2D yes=t",1);
		} else if (is_tensortrain()) {
			MADNESS_EXCEPTION("no type conversion for TT_2D yes=t",1);
        }
        return result;
    }

	SVDTensor<T>& get_svdtensor() {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<1>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return SVDTensor from LowRankTensor",1);
		}
	}

	const SVDTensor<T>& get_svdtensor() const {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<1>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return SVDTensor from LowRankTensor",1);
		}
	}

	Tensor<T>& get_tensor() {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<0>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return Tensor from LowRankTensor",1);
		}
	}

	const Tensor<T>& get_tensor() const {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<0>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return Tensor from LowRankTensor",1);
		}
	}

	TensorTrain<T>& get_tensortrain() {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<2>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return TensorTrain from LowRankTensor",1);
		}
	}

	const TensorTrain<T>& get_tensortrain() const {
		MADNESS_ASSERT(is_assigned());
		try {
			return *(std::get<2>(tensor).get());
		} catch (...) {
			MADNESS_EXCEPTION("failure to return TensorTrain from LowRankTensor",1);
		}
	}

    /// general slicing, shallow; for temporary use only!
    SliceLowRankTensor<T> operator()(const std::vector<Slice>& s) {
        return SliceLowRankTensor<T>(*this,s);
    }

    /// general slicing, shallow; for temporary use only!
    const SliceLowRankTensor<T> operator()(const std::vector<Slice>& s) const {
        return SliceLowRankTensor<T>(*this,s);
    }


    /// deep copy
	friend GenTensor copy(const GenTensor& other) {
		GenTensor<T> result;
		if (other.is_assigned()) std::visit([&result](auto& obj) {result=copy(*obj);}, other.tensor);
        return result;
    }

    /// return the tensor type
	TensorType tensor_type() const {
		if (index()==0) return TT_FULL;
		if (index()==1) return TT_2D;
		if (index()==2) return TT_TENSORTRAIN;
		MADNESS_EXCEPTION("confused tensor types ",1);
	}

    constexpr bool is_full_tensor() const {
        return (index()==0);
    }

	constexpr bool is_svd_tensor() const {
        return (index()==1);
	}

	constexpr bool is_tensortrain() const {
        return (index()==2);
	}

	bool is_of_tensortype(const TensorType& tt) const {
		if ((index()==0) and (tt==TT_FULL)) return true;
		if ((index()==1) and (tt==TT_2D)) return true;
		if ((index()==2) and (tt==TT_TENSORTRAIN)) return true;
		return false;
	}

	template<typename Q, typename R>
	friend bool is_same_tensor_type(const GenTensor<R>& rhs, const GenTensor<Q>& lhs);

	int index() const {
		return is_assigned() ? tensor.index() : -1;
                }

	GenTensor& convert_inplace(const TensorArgs& targs) {

        // fast return
        if (not is_assigned()) return *this;
        if (is_of_tensortype(targs.tt)) return *this;
//		if (targs.tt==TT_DYNAMIC) if (is_svd_tensor()) return *this;

        // target is full tensor
        if (targs.tt == TT_FULL) {
            *this = this->full_tensor_copy();
        }

        // source is full tensor: construct the corresponding representation
        else if (is_full_tensor()) {
            *this = GenTensor<T>(get_tensor(), targs);
        }

        // TT_TENSORTRAIN TO TT_2D
        else if ((is_tensortrain()) and (targs.tt == TT_2D)) {
            Tensor<T> U, VT;
            Tensor<typename Tensor<T>::scalar_type> s;
            get_tensortrain().two_mode_representation(U, VT, s);
            long rank = s.size();
            if (rank == 0) {
                *this = SVDTensor<T>(get_tensortrain().ndim(), get_tensortrain().dims(), ndim() / 2);
                return *this;
            }

            long n = 1, m = 1;
            for (int i = 0; i < U.ndim() - 1; ++i) n *= U.dim(i);
            for (int i = 1; i < VT.ndim(); ++i) m *= VT.dim(i);
            MADNESS_ASSERT(rank * n == U.size());
            MADNESS_ASSERT(rank * m == VT.size());
            U = copy(transpose(U.reshape(n, rank)));   // make it contiguous
            VT = VT.reshape(rank, m);
            SVDTensor<T> svdtensor(s, U, VT, ndim(), dims());
            svdtensor.normalize();
            *this = svdtensor;
        }
        else if ((is_svd_tensor()) and (targs.tt == TT_TENSORTRAIN)) {
            TensorTrain<T> tt(this->full_tensor_copy(),targs.thresh);
            *this=tt;
        } else {
            print("conversion from type ", index(), "to type", targs.tt, "not supported");
            MADNESS_EXCEPTION("type conversion not supported in LowRankTensor::convert ", 1);
        }
        return *this;
    }

	/// convert this to a new LowRankTensor of given tensor type
	GenTensor convert(const TensorArgs& targs) const {

		// deep copy for same type
		if (is_of_tensortype(targs.tt)) return copy(*this);

		// new LRT will be newly constructed anyways
		if (is_full_tensor()) return GenTensor<T>(get_tensor(),targs);

		GenTensor<T> result(*this);	// shallow
		result.convert_inplace(targs);
        return result;
    }

    long ndim() const {
		return (is_assigned()) ? ptr()->ndim() : -1;
    }

    /// return the number of entries in dimension i
    long dim(const int i) const {
		MADNESS_ASSERT(is_assigned());
		return ptr()->dim(i);
    }

	/// return the number of entries in dimension i
	const long* dims() const {
		MADNESS_ASSERT(is_assigned());
		return ptr()->dims();
	}

    void normalize() {
		if (is_svd_tensor()) get_svdtensor().normalize();
    }

    float_scalar_type normf() const {
		float_scalar_type norm;
		std::visit([&norm](auto& obj) {norm=obj->normf();}, tensor);
		return norm;
    }

    float_scalar_type svd_normf() const {
		float_scalar_type norm;
		if (is_svd_tensor()) return get_svdtensor().svd_normf();
		std::visit([&norm](auto& obj) {norm=obj->normf();}, tensor);
		return norm;
    }


    /// Inplace multiplication by scalar of supported type (legacy name)

    /// @param[in] x Scalar value
    /// @return %Reference to this tensor
    template <typename Q>
	typename IsSupported<TensorTypeData<Q>,GenTensor<T>&>::type
    scale(Q fac) {
		if (not is_assigned()) return *this;
		std::visit([&fac](auto& obj) {obj->scale(T(fac));}, tensor);
        return *this;
    }

    Tensor<T> full_tensor_copy() const {
		if (not is_assigned()) return Tensor<T>();
		else if (is_full_tensor()) return copy(get_tensor());
		else if (is_svd_tensor()) return get_svdtensor().reconstruct();
		else if (is_tensortrain()) return get_tensortrain().reconstruct();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return Tensor<T>();
    }

	/// return a view (shallow copy) of the full tensor, if possible,
	/// otherwise reconstruct to full tensor and return a deep copy
    Tensor<T> full_tensor() const {
		if (is_full_tensor()) return get_tensor();
		return full_tensor_copy();
    }

	/// return a reference of the full tensor, presumably for in-place modification
	/// WARNING: works only if this is indeed a full tensor
	Tensor<T>& full_tensor() {
		MADNESS_ASSERT(is_full_tensor());
		return get_tensor();
    }


    /// reconstruct this to return a full tensor
    Tensor<T> reconstruct_tensor() const {

		if (is_full_tensor()) return full_tensor();
		if (is_svd_tensor() or is_tensortrain()) return full_tensor_copy();
        return Tensor<T>();
    }


    static double facReduce() {return 1.e-3;}
    static double fac_reduce() {return 1.e-3;}

    long rank() const {
		if (is_full_tensor()) return -1;
		else if (is_svd_tensor()) return get_svdtensor().rank();
		else if (is_tensortrain()) {
			std::vector<long> r=get_tensortrain().ranks();
            return *(std::max_element(r.begin(), r.end()));
        }
        return 0l;
    }

	bool is_assigned() const {
		return ptr() ? true : false;
	}

	bool has_data() const {return size()>0;}

    bool has_no_data() const {return (not has_data());}

    long size() const {
		return (is_assigned()) ? ptr()->size() : 0;
	}

	long nCoeff() const {
		if (is_full_tensor()) return get_tensor().size();
		else if (is_svd_tensor()) return get_svdtensor().nCoeff();
		else if (is_tensortrain()) return get_tensortrain().real_size();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return false;
    }

    long real_size() const {
		if (is_full_tensor()) return get_tensor().size();
		else if (is_svd_tensor()) return get_svdtensor().real_size();
		else if (is_tensortrain()) return get_tensortrain().real_size();
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return false;
    }

    /// returns the trace of <this|rhs>
    template<typename Q>
	TENSOR_RESULT_TYPE(T,Q) trace_conj(const GenTensor<Q>& rhs) const {

        if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in LowRankTensor, sorry",1);
        if (TensorTypeData<Q>::iscomplex) MADNESS_EXCEPTION("no complex trace in LowRankTensor, sorry",1);

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        // fast return if possible
        if ((this->rank()==0) or (rhs.rank()==0)) return resultT(0.0);

		MADNESS_ASSERT(is_same_tensor_type(*this,rhs));

		if (is_full_tensor()) return get_tensor().trace_conj(rhs.get_tensor());
		else if (is_svd_tensor()) return trace(get_svdtensor(),rhs.get_svdtensor());
		else if (is_tensortrain()) return get_tensortrain().trace(rhs.get_tensortrain());
        else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
        return TENSOR_RESULT_TYPE(T,Q)(0);
    }

    /// multiply with a number
    template<typename Q>
	GenTensor<TENSOR_RESULT_TYPE(T,Q)> operator*(const Q& x) const {
		GenTensor<TENSOR_RESULT_TYPE(T,Q)> result(copy(*this));
        result.scale(x);
        return result;
    }

	GenTensor operator+(const GenTensor& other) {
		GenTensor<T> result=copy(*this);
		result.gaxpy(1.0,other,1.0);
		return result;
	}

	GenTensor operator+(const SliceLowRankTensor<T>& other) {
		GenTensor<T> result=copy(*this);
		std::array<Slice,TENSOR_MAXDIM> s0;
		s0.fill(_);
		result.gaxpy(1.0,s0,other,1.0,other.thisslice);
		return result;
	}

	GenTensor& operator+=(const GenTensor& other) {
        gaxpy(1.0,other,1.0);
        return *this;
    }

	GenTensor& operator+=(const SliceLowRankTensor<T>& other) {
		std::array<Slice,TENSOR_MAXDIM> s0;
		s0.fill(_);
		this->gaxpy(1.0,s0,other,1.0,other.thisslice);
		return *this;
	}

	GenTensor operator-(const GenTensor& other) {
		GenTensor<T> result=copy(*this);
		result.gaxpy(1.0,other,-1.0);
		return result;
	}

	GenTensor& operator-=(const GenTensor& other) {
        gaxpy(1.0,other,-1.0);
        return *this;
    }

	GenTensor& operator-=(const SliceLowRankTensor<T>& other) {
		std::array<Slice,TENSOR_MAXDIM> s0;
		s0.fill(_);
		this->gaxpy(1.0,s0,other,-1.0,other.thisslice);
		return *this;
	}

	GenTensor& gaxpy(const T alpha, const GenTensor& other, const T beta) {

		// deliberately excluding gaxpys for different tensors due to efficiency considerations!
		MADNESS_ASSERT(is_same_tensor_type(*this,other));
		MADNESS_ASSERT(is_assigned());
		if (is_full_tensor()) get_tensor().gaxpy(alpha,other.get_tensor(),beta);
		else if (is_svd_tensor()) get_svdtensor().gaxpy(alpha,other.get_svdtensor(),beta);
		else if (is_tensortrain()) get_tensortrain().gaxpy(alpha,other.get_tensortrain(),beta);
		else {
			MADNESS_EXCEPTION("unknown tensor type in LowRankTensor::gaxpy",1);
		}
		return *this;
	}

	GenTensor& gaxpy(const T alpha, std::array<Slice,TENSOR_MAXDIM> s0,
			const GenTensor& other, const T beta, std::array<Slice,TENSOR_MAXDIM> s1) {

		// deliberately excluding gaxpys for different tensors due to efficiency considerations!
		MADNESS_ASSERT(is_same_tensor_type(*this,other));
		MADNESS_ASSERT(is_assigned());

		if (is_full_tensor()) {
			get_tensor()(s0).gaxpy(alpha,other.get_tensor()(s1),beta);
		} else if (is_svd_tensor()) {
			get_svdtensor().inplace_add(other.get_svdtensor(),s0,s1,alpha,beta);
		} else if (is_tensortrain()) {
			MADNESS_ASSERT(alpha==1.0);
			get_tensortrain().gaxpy(s0, other.get_tensortrain(), beta, s1);
        } else {
			MADNESS_EXCEPTION("unknown tensor type in LowRankTensor::gaxpy",1);
        }
        return *this;
    }

    /// assign a number to this tensor
	GenTensor& operator=(const T& number) {
		std::visit([&number](auto& obj) {*obj=number;}, tensor);
        return *this;

    }

	void add_SVD(const GenTensor& other, const double& thresh) {
		if (is_full_tensor()) get_tensor()+=other.get_tensor();
		else if (is_svd_tensor()) get_svdtensor().add_SVD(other.get_svdtensor(),thresh*facReduce());
		else if (is_tensortrain()) get_tensortrain()+=(other.get_tensortrain());
        else {
			MADNESS_EXCEPTION("unknown tensor type in LowRankTensor::add_SVD",1);
        }
    }

    /// Inplace multiply by corresponding elements of argument Tensor
	GenTensor<T>& emul(const GenTensor<T>& other) {

		// deliberately excluding emuls for different tensors due to efficiency considerations!
		MADNESS_ASSERT(is_same_tensor_type(*this,other));

		// binary operation with the visitor pattern
		//        std::visit([&other](auto& obj) {obj.emul(other.tensor);}, tensor);
		if (is_full_tensor()) get_tensor().emul(other.get_tensor());
		else if (is_svd_tensor()) get_svdtensor().emul(other.get_svdtensor());
		else if (is_tensortrain()) get_tensortrain().emul(other.get_tensortrain());
		else {
			MADNESS_EXCEPTION("unknown tensor type in LowRankTensor::gaxpy",1);
        }
        return *this;

	}

    void reduce_rank(const double& thresh) {
		if (is_svd_tensor()) get_svdtensor().divide_and_conquer_reduce(thresh*facReduce());
		if (is_tensortrain()) get_tensortrain().truncate(thresh*facReduce());
    }


public:

    /// Transform all dimensions of the tensor t by the matrix c

    /// \ingroup tensor
    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
    /// The input dimensions of \c t must all be the same and agree with
    /// the first dimension of \c c .  The dimensions of \c c may differ in
    /// size.
	template <typename R, typename Q>
	friend GenTensor<TENSOR_RESULT_TYPE(R,Q)> transform(
			const GenTensor<R>& t, const Tensor<Q>& c);

    /// Transform all dimensions of the tensor t by distinct matrices c

    /// \ingroup tensor
    /// Similar to transform but each dimension is transformed with a
    /// distinct matrix.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c[0](i',i) c[1](j',j) c[2](k',k) ...
    /// \endcode
    /// The first dimension of the matrices c must match the corresponding
	/// dimension of t.    template <typename R, typename Q>
	template <typename R, typename Q>
	friend GenTensor<TENSOR_RESULT_TYPE(R,Q)> general_transform(
			const GenTensor<R>& t, const Tensor<Q> c[]);

    /// Transforms one dimension of the tensor t by the matrix c, returns new contiguous tensor

    /// \ingroup tensor
    /// \code
    /// transform_dir(t,c,1) = r(i,j,k,...) = sum(j') t(i,j',k,...) * c(j',j)
    /// \endcode
	/// @param[in] t Tensor to transform (size of dim to be transformed must match size of first dim of \c c )
    /// @param[in] c Matrix used for the transformation
    /// @param[in] axis Dimension (or axis) to be transformed
	/// @result Returns a new, contiguous tensor    template <typename R, typename Q>
	template <typename R, typename Q>
	friend GenTensor<TENSOR_RESULT_TYPE(R,Q)> transform_dir(
			const GenTensor<R>& t, const Tensor<Q>& c, const int axis);


	std::string what_am_i() const {
		TensorType tt;
		if (this->is_full_tensor()) tt=TT_FULL;
		if (this->is_svd_tensor()) tt=TT_2D;
		if (this->is_tensortrain()) tt=TT_TENSORTRAIN;
		return TensorArgs::what_am_i(tt);
	};


	/// might return a NULL pointer!
	const BaseTensor* ptr() const {
		const BaseTensor* p;
		std::visit([&p](auto& obj) {p=dynamic_cast<const BaseTensor*>(obj.get());}, tensor);
		return p;
        }

private:

	/// holding the implementation of the low rank tensor representations
	//	std::variant<Tensor<T>, SVDTensor<T>, TensorTrain<T> > tensor;
	std::variant<std::shared_ptr<Tensor<T> >,
	std::shared_ptr<SVDTensor<T> >,
	std::shared_ptr<TensorTrain<T> > > tensor;

};



namespace archive {
/// Serialize a tensor
template <class Archive, typename T>
struct ArchiveStoreImpl< Archive, GenTensor<T> > {

	friend class GenTensor<T>;
	/// Stores the GenTensor to an archive
	static void store(const Archive& ar, const GenTensor<T>& t) {
		int index1=t.index();
		ar & index1;
		if (index1==0) {
			const Tensor<T>& tt=t.get_tensor();
			ar & tt;
		} else if (index1==1) {
			const SVDTensor<T>& tt=t.get_svdtensor();
			ar & tt;
		} else if (index1==2) {
			const TensorTrain<T>& tt=t.get_tensortrain();
			ar & tt;
		}
	};
};


/// Deserialize a tensor ... existing tensor is replaced
template <class Archive, typename T>
struct ArchiveLoadImpl< Archive, GenTensor<T> > {

	friend class GenTensor<T>;
	/// Replaces this GenTensor with one loaded from an archive
	static void load(const Archive& ar, GenTensor<T>& tensor) {
		int index=-2;
		ar & index;
		if (index==0) {
			Tensor<T> tt;
			ar & tt;
			tensor=tt;
		} else if (index==1) {
			SVDTensor<T> tt;
			ar & tt;
			tensor=tt;
		} else if (index==2) {
			TensorTrain<T> tt;
			ar & tt;
			tensor=tt;
		} else if (index==-1) {	 // defined value: empty tensor
			;
		} else {
			MADNESS_EXCEPTION("unknow tensor type",1);
		}


	};
};
};

/// type conversion implies a deep copy

/// @result Returns a new tensor that is a deep copy of the input
template <class Q, class T>
GenTensor<Q> convert(const GenTensor<T>& other) {

	// simple return
	if (std::is_same<Q, T>::value) return copy(other);

	GenTensor<Q> result;
	if (other.is_full_tensor())
		result=Tensor<Q>(convert<Q,T>(other.get_tensor()));
	if (other.is_svd_tensor())
        MADNESS_EXCEPTION("no type conversion for SVDTensors",1);
	if (other.is_tensortrain())
        MADNESS_EXCEPTION("no type conversion for TensorTrain",1);
    return result;
}


/// change representation to targ.tt
template<typename T>
void change_tensor_type(GenTensor<T>& t, const TensorArgs& targs) {
	t.convert_inplace(targs);
}

/// outer product of two Tensors, yielding a low rank tensor

/// do the outer product of two tensors; distinguish these tensortype cases by
/// the use of final_tensor_type
///  - full x full -> full
///  - full x full -> SVD                           ( default )
///  - TensorTrain x TensorTrain -> TensorTrain
/// all other combinations are currently invalid.
template <class T, class Q>
GenTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const GenTensor<T>& t1,
		const GenTensor<Q>& t2, const TensorArgs final_tensor_args=TensorArgs(-1.0,TT_2D)) {

    typedef TENSOR_RESULT_TYPE(T,Q) resultT;


	MADNESS_ASSERT(is_same_tensor_type(t1,t2));

    if (final_tensor_args.tt==TT_FULL) {
		MADNESS_ASSERT(t1.is_full_tensor());
		Tensor<resultT> t(outer(t1.get_tensor(),t2.get_tensor()));
		return GenTensor<resultT>(t);

    } else if (final_tensor_args.tt==TT_2D) {
		MADNESS_ASSERT(t1.is_full_tensor());

        // srconf is shallow, do deep copy here
        const Tensor<T> lhs=t1.full_tensor_copy();
        const Tensor<Q> rhs=t2.full_tensor_copy();

        const long k=lhs.dim(0);
		const long ndim=lhs.ndim()+rhs.ndim();
        long size=1;
        for (int i=0; i<lhs.ndim(); ++i) size*=k;
        MADNESS_ASSERT(size==lhs.size());
        MADNESS_ASSERT(size==rhs.size());
        MADNESS_ASSERT(lhs.size()==rhs.size());

        Tensor<double> weights(1);
        weights=1.0;

		std::array<long,TENSOR_MAXDIM> dims;
		for (int i=0; i<t1.ndim(); ++i) dims[i]=t1.dim(i);
		for (int i=0; i<t2.ndim(); ++i) dims[i+t1.ndim()]=t2.dim(i);

		SRConf<resultT> srconf(weights,lhs.reshape(1,lhs.size()),rhs.reshape(1,rhs.size()),ndim,dims.data(),t1.ndim());
//        srconf.normalize();
		return GenTensor<resultT>(SVDTensor<resultT>(srconf));

    } else if (final_tensor_args.tt==TT_TENSORTRAIN) {
		MADNESS_ASSERT(t1.is_tensortrain());
		MADNESS_ASSERT(t2.is_tensortrain());
		return outer(t1.get_tensortrain(),t2.get_tensortrain());
    } else {
        MADNESS_EXCEPTION("you should not be here",1);
    }
	return GenTensor<TENSOR_RESULT_TYPE(T,Q)>();

            }


/// outer product of two Tensors, yielding a low rank tensor
template <class T, class Q>
GenTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const Tensor<T>& lhs2,
		const Tensor<Q>& rhs2, const TensorArgs final_tensor_args) {

	typedef TENSOR_RESULT_TYPE(T,Q) resultT;

	// prepare lo-dim tensors for the outer product
	TensorArgs targs;
	targs.thresh=final_tensor_args.thresh;
	if (final_tensor_args.tt==TT_FULL) targs.tt=TT_FULL;
	else if (final_tensor_args.tt==TT_2D) targs.tt=TT_FULL;
	else if (final_tensor_args.tt==TT_TENSORTRAIN) targs.tt=TT_TENSORTRAIN;
	else {
		MADNESS_EXCEPTION("confused tensor args in outer_low_rank",1);
                }

	GenTensor<T> lhs(lhs2,targs);
	GenTensor<Q> rhs(rhs2,targs);
	GenTensor<resultT> result=outer(lhs,rhs,final_tensor_args);
	return result;
            }


/// The class defines tensor op scalar ... here define scalar op tensor.
template <typename T, typename Q>
typename IsSupported < TensorTypeData<Q>, GenTensor<T> >::type
operator*(const Q& x, const GenTensor<T>& t) {
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
GenTensor<T> reduce(std::list<GenTensor<T> >& addends, double eps, bool are_optimal=false) {

	// fast return
	addends.remove_if([](auto element) {return not element.is_assigned();});
	addends.remove_if([](auto element) {return element.rank()==0;});
	if (addends.size()==0) return GenTensor<T>();


	if (addends.front().is_svd_tensor()) {
		std::list<SVDTensor<T> > addends1;
		for (auto a : addends) addends1.push_back(a.get_svdtensor());
		return reduce(addends1,eps*GenTensor<T>::facReduce());
      }
	// make error relative
	eps=eps/addends.size();

	// if the addends are not in SVD format do that now so that we can call add_svd later
	if (not are_optimal) {
		for (auto element : addends) element.reduce_rank(eps);
	}

	// remove zero ranks and sort the list according to the gentensor's ranks
	addends.remove_if([](auto element) {return element.rank()==0;});
	if (addends.size()==0) return GenTensor<T>();
	addends.sort([](auto element1, auto element2) {return element1.rank()<element2.rank();});

	// do the additions
	GenTensor<T> result=copy(addends.front());
	addends.pop_front();
	for (auto element : addends) result.add_SVD(element,eps);
	addends.clear();

      return result;
}




/// implements a temporary(!) slice of a LowRankTensor
template<typename T>
class SliceLowRankTensor : public GenTensor<T> {
	//class SliceLowRankTensor {
public:

	std::array<Slice,TENSOR_MAXDIM> thisslice;
	//    GenTensor<T>* lrt;

    // all ctors are private, only accessible by GenTensor

    /// default ctor
    SliceLowRankTensor<T> () {}

    /// ctor with a GenTensor; shallow
	SliceLowRankTensor<T> (const GenTensor<T>& gt, const std::vector<Slice>& s)
    				: GenTensor<T>(const_cast<GenTensor<T>& > (gt)) {
		//        : Tensor<T>(const_cast<Tensor<T>&>(t)) //!!!!!!!!!!!
		for (int i=0; i<s.size(); ++i) thisslice[i]=s[i];
	}

	/// ctor with a GenTensor; shallow
	SliceLowRankTensor<T> (const GenTensor<T>& gt, const std::array<Slice,TENSOR_MAXDIM>& s)
    				: GenTensor<T>(&gt), thisslice(s) {}

public:

    /// assignment as in g(s) = g1;
	SliceLowRankTensor<T>& operator=(const GenTensor<T>& rhs) {
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
	SliceLowRankTensor<T>& operator+=(const GenTensor<T>& rhs) {
		std::array<Slice,TENSOR_MAXDIM> rhs_slice;
		rhs_slice.fill(_);
		gaxpy(thisslice,rhs,rhs_slice,1.0);
        return *this;
    }

    /// inplace subtraction as in g(s)-=g1
	SliceLowRankTensor<T>& operator-=(const GenTensor<T>& rhs) {
		std::array<Slice,TENSOR_MAXDIM> rhs_slice;
		rhs_slice.fill(_);
		gaxpy(thisslice,rhs,rhs_slice,-1.0);
        return *this;
    }

    /// inplace addition as in g(s)+=g1(s)
    SliceLowRankTensor<T>& operator+=(const SliceLowRankTensor<T>& rhs) {
		gaxpy(thisslice,rhs,rhs.thisslice,1.0);
		return *this;
	}

	/// inplace addition as in g(s)-=g1(s)
	SliceLowRankTensor<T>& operator-=(const SliceLowRankTensor<T>& rhs) {
		gaxpy(thisslice,rhs,rhs.thisslice,-1.0);
		return *this;
	}

	/// *this = *this(s) + beta * rhs
	void gaxpy(const std::array<Slice,TENSOR_MAXDIM>& lslice, const GenTensor<T>& rhs,
			const std::array<Slice,TENSOR_MAXDIM>& rslice, const double& beta) {

        // fast return if possible
		if (rhs.has_no_data() or rhs.rank()==0) return;

		if (this->has_data()) MADNESS_ASSERT(is_same_tensor_type(*this,rhs));

		if (this->is_full_tensor()) {
			this->get_tensor()(thisslice).gaxpy(1.0,rhs.get_tensor()(rslice),beta);

		} else if (this->is_svd_tensor()) {
			this->get_svdtensor().inplace_add(rhs.get_svdtensor(),thisslice,rslice, 1.0, beta);

		} else if (this->is_tensortrain()) {
			this->get_tensortrain().gaxpy(thisslice,rhs.get_tensortrain(),beta,rslice);
        }
		return ;
    }

    /// inplace zero-ing as in g(s)=0.0
    SliceLowRankTensor<T>& operator=(const T& number) {
        MADNESS_ASSERT(number==T(0.0));

		if (this->is_full_tensor()) {
			this->get_tensor()(thisslice)=0.0;

		} else if (this->is_svd_tensor()) {
			MADNESS_ASSERT(this->get_svdtensor().has_structure());
			SliceLowRankTensor<T> tmp(*this);
			this->get_svdtensor().inplace_add(tmp.get_svdtensor(),thisslice,thisslice, 1.0, -1.0);

		} else if (this->is_tensortrain()) {
			this->get_tensortrain().gaxpy(thisslice,this->get_tensortrain(),-1.0,thisslice);
		} else {
			MADNESS_EXCEPTION("you should not be here",1);
        }
        return *this;
    }

	friend GenTensor<T> copy(const SliceLowRankTensor<T>& other) {
		GenTensor<T> result;
		const std::array<Slice,TENSOR_MAXDIM> s=other.thisslice;
		if (other.is_full_tensor())
			result=Tensor<T>(copy(other.get_tensor()(s)));
		else if (other.is_svd_tensor())
			result=SVDTensor<T>(other.get_svdtensor().copy_slice(s));
		else if (other.is_tensortrain())
			result=TensorTrain<T>(copy(other.get_tensortrain(),s));
        else {
		}
		return result;
	}


};


template<typename Q, typename R>
bool is_same_tensor_type(const GenTensor<R>& rhs, const GenTensor<Q>& lhs) {
	return (rhs.tensor.index()==lhs.tensor.index());
}

template <typename R, typename Q>
GenTensor<TENSOR_RESULT_TYPE(R,Q)> transform(
		const GenTensor<R>& t, const Tensor<Q>& c) {
	typedef TENSOR_RESULT_TYPE(R,Q) resultT;
	GenTensor<resultT> result;
	std::visit([&result, &c](auto& obj) {result=transform(*obj,c);}, t.tensor);
	return result;
        }

template <typename R, typename Q>
GenTensor<TENSOR_RESULT_TYPE(R,Q)> general_transform(
		const GenTensor<R>& t, const Tensor<Q> c[]) {
	typedef TENSOR_RESULT_TYPE(R,Q) resultT;
	GenTensor<resultT> result;
	std::visit([&result, &c](auto& obj) {result=general_transform(*obj,c);}, t.tensor);
	return result;
}

template <typename R, typename Q>
GenTensor<TENSOR_RESULT_TYPE(R,Q)> transform_dir(
		const GenTensor<R>& t, const Tensor<Q>& c, const int axis) {
	GenTensor<TENSOR_RESULT_TYPE(R,Q)> result;
	std::visit([&result, &c, &axis](auto& obj) {result=transform_dir(*obj,c,axis);}, t.tensor);
        return result;

    }



} // namespace madness

#endif /* MADNESS_TENSOR_LOWRANKTENSOR_H_ */
