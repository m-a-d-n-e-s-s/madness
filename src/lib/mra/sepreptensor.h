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

#ifndef SEPREPTENSOR_H_
#define SEPREPTENSOR_H_
/*!
 *
 * \file SepRepTensor.h
 * \brief Provides a tensor with taking advantage of possibly low rank
 *
 * a SepRepTensor is a generalized form of a Tensor
 *
 * for now only little functionality shall be implemented; feel free to extend
 *
 * a consequence is that we (usually) can't directly access individual
 * matrix elements
 *
 * LowRankTensors' ctors will most likely not be shallow
 *
 * note that a LowRankTensor might have zero rank, but is still a valid
 * tensor, and should therefore return a FullTensor with zeros in it.
 *
 * \par Slicing in LowRankTensors:
 *
 * LowRankTensors differ from FullTensors in that we can't directly access
 * individual matrix elements, and thus can't directly assign or manipulate slices
 * as lvalues. For rvalues we simply provide a slices of the constituent vectors
 * in SRConf, which are valid LowRankTensors by themselves
 * \code
 *  lhs = rhs(s)
 * \endcode
 *
 * However, for manipulations like
 * \code
 *  lhs(s) = rhs
 * \endcode
 * we introduce the SliceLowRankTensor, which has no further effect than to
 * indicate that the operation is done on a slice. Effectively we turn the rhs
 * into a slice (which is larger than the original LRT), and assign
 * that to lhs. Actually, we wouldn't need to do that, but if we don't we
 * render it impossible for the user to check if the shapes conform.
 *
 * Manipulations of slices of LowRankTensors are heavily restricted, but should
 * cover the most important cases:
 * - assignment to a slice that was zero before, as in (t being some Tensor of
 *   dimensions 2,2,2,2, and s being the (2,2,2,2) slice); this is performed
 *   by inplace addition
 * 	  \code
 *    LowRankTensor lrt1(3,3,3,3), lrt2(t);
 *    lrt1(s) = lrt2;
 * 	  \endcode
 * - assignment to zero; done by inplace subtraction of the slice
 *    \code
 *    LowRankTensor lrt1(t);
 *    lrt1(s) = 0.0;
 *    \endcode
 * - in-place addition
 *    \code
 *    LowRankTensor lrt1(3,3,3,3), lrt2(t);
 *    lrt1(s) += lrt2;
 *    \endcode
 *
 * Note that *all* of these operation increase the rank of lhs
 */


#include "tensor/tensor.h"
#include <mra/funcdefaults.h>
#include "mra/seprep.h"

namespace madness {


    /// abstract base class for the SepRepTensor
    template<typename T>
    class SepRepTensor {

    public:

    	typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

    	/// default ctor
    	SepRepTensor() {};

    	/// implement the virtual constructor idiom
    	virtual SepRepTensor* clone() const =0;

    	/// add-to-this operator; in LowRank requires reallocation and is expensive!
//    	virtual SepRepTensor& operator+=(const SepRepTensor& rhs)=0;

    	/// virtual destructor
    	virtual ~SepRepTensor() {};

    	/// return the type of the derived class
        virtual TensorType type() const =0;

        /// return the type of the derived class for me
        virtual std::string what_am_i() const =0;

        /// Returns if this Tensor exists
        virtual bool has_data() const =0;

    	virtual long size() const =0;

    	virtual long dim(const int) const =0;

    	virtual long ndim() const =0;

    	virtual long rank() const =0;

    	virtual Tensor<T> reconstructTensor() const =0;
    	virtual Tensor<T>& fullTensor() =0;
    	virtual const Tensor<T>& fullTensor() const =0;

    	virtual T operator()(long i, long j, long k) const=0;

    	/// compute the Frobenius norm
    	virtual float_scalar_type normf() const =0;

    	/// scale by a number
    	virtual void scale(T a) {throw;};

    };

    template<typename T>
    class LowRankTensor : public SepRepTensor<T> {

    private:
    	typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

    public:
    	SepRep<T> data;		///< the tensor data


       	/// default constructor (tested)
    	LowRankTensor<T>() : data(SepRep<T>()) {
    	};

    	/// copy constructor (tested)
    	LowRankTensor<T>(const LowRankTensor<T>& rhs) {
    		data=rhs.data;
    	};

    	/// ctor with a SepRep
    	LowRankTensor<T>(const SepRep<T>& rhs) : data(rhs) {
    	}

    	/// ctor with regular Tensor (tested)
    	LowRankTensor<T>(const Tensor<T>& t, const double& eps, const TensorType& tt) {
			data=SepRep<T>(t,eps,tt);
    	}

    	/// ctor with ndim tensors for each dimension
    	LowRankTensor<T>(const Tensor<double>& weights,
    			const std::vector<Tensor<T> >& t, const TensorType& tt)
    			: data(weights,t,tt) {
    	}

    	/// dtor
    	~LowRankTensor<T>() {}

    	/// assignment with rhs=LowRankTensor
    	LowRankTensor<T>& operator=(const LowRankTensor<T>& rhs) {
    		MADNESS_ASSERT(0);
    	};

    	/// assignment to a number
    	LowRankTensor<T>& operator=(const T& a) {
    		MADNESS_ASSERT(0);
    	}

    	/// slicing to assign lhs=rhs(s) (tested)
    	LowRankTensor<T> operator()(const std::vector<Slice>& s) const {

    		LowRankTensor<T> result(this->data(s));
    		return result;
    	}


#if 0
    	/// enlarge the dimensions of this, while keeping the information
    	void project(const std::vector<Slice>& s) {
    		// symbolically: MADNESS_ASSERT( s > dim );
    		MADNESS_ASSERT(0);
    	}

    	/// set a slices of this to zero; will double the rank iff not the
    	/// whole tensor is set to zero
    	/// *this(s) = 0.0;
    	void zero_out(const std::vector<Slice>& s) {
    		MADNESS_ASSERT(0);
    	}
#endif

    	/// implement the virtual constructor idiom
    	LowRankTensor<T>* clone() const {
    		return new LowRankTensor<T>(*this);
    	}

       	/// return this tensor type
        TensorType type() const {return data.tensor_type();};

        /// return the type of the derived class for me
        std::string what_am_i() const {return "LowRank";};

        /// return the rank of this
        long rank() const {return data.rank();};

        /// Returns if this Tensor exists
        bool has_data() const {return data.is_valid();};

        /// Return the number of coefficients (valid for all SPR)
        long size() const {
        	if (!this->has_data()) return 0;
        	return data.nCoeff();
        };

    	virtual long dim(const int) const {return data.get_k();};

		virtual long ndim() const {return data.dim();};

		/// reconstruct this to a full tensor
		Tensor<T> reconstructTensor() const {
			MADNESS_ASSERT(data.is_valid());
			return data.reconstructTensor();

		};

		/// (Don't) return a full tensor representation
		Tensor<T>& fullTensor() {
			MADNESS_EXCEPTION("No fullTensor in LowRankTensor; reconstruct first",0);
		};

		/// (Don't) return a full tensor representation
		Tensor<T>& fullTensor() const {
			MADNESS_EXCEPTION("No fullTensor in LowRankTensor; reconstruct first",0);
		};

    	T operator()(long i, long j, long k) const {throw;};

    	/// compute the Frobenius norm
    	float_scalar_type normf() const {throw;};

    	/// scale by a number
    	void scale(T a) {throw;};

    };

    /// derived from LowRankTensor, intended only for assignment

    /// we need this for all assignments of the form (s being Slices):
    /// \code
    ///  lhs(s) = rhs
    /// \endcode
    /// for assignments like
    /// \code
    /// lhs = rhs(s)
    /// use the LowRankTensor class directly
    template<typename T>
    class SliceLowRankTensor : private LowRankTensor<T> {

    private:
    	SliceLowRankTensor<T>() {};

    public:

    	/// ctor
    	SliceLowRankTensor<T>(const LowRankTensor<T>& lrt, const std::vector<Slice>& s) {
    	}

    	/// assignment
    	SliceLowRankTensor<T>& operator=(const LowRankTensor<T>& rhs) {
    	}

    	SliceLowRankTensor<T>& operator=(const SliceLowRankTensor<T>& rhs) {
    	}


    };

    /// the case of a full rank tensor, as tensor.h

    /**
     * this class merely wraps the tensor class, but we need to do so for
     * we want it to be a derived class of SepRepTensor;
     */
    template<typename T>
    class FullTensor : public SepRepTensor<T> {

    private:

    	typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

    	Tensor<T> data;	///< the tensor data

    public:

    	/// default constructor constructs an empty tensor
    	FullTensor<T>() {
    		data=Tensor<T>() ;
    	};

    	/// copy constructor (tested)
		FullTensor<T>(const FullTensor<T>& rhs) {
			*this=rhs;
		};

		/// ctor with a regular Tensor (tested)
		FullTensor<T>(const Tensor<T>& t) {
			data=t;
		}

    	/// assignment with rhs=FullTensor (tested)
    	FullTensor<T>& operator=(const FullTensor<T>& rhs) {
    		if (this!=&rhs) {
    			data=rhs.data;
    		}
    		return *this;
    	}

    	/// assignment with rhs=Tensor (tested)
    	FullTensor<T>& operator=(const Tensor<T>& rhs) {
    		if (&this->data!=&rhs) {
    			data=rhs;
    		}
    		return *this;
    	}

    	/// assignment to a number (tested)
    	FullTensor<T>& operator=(const T& a) {
    		data=a;
    		return *this;
    	}

    	/// general slicing (tested)
    	SliceTensor<T> operator()(const std::vector<Slice>& s) {
//    		return SliceTensor<T> (data,&(s[0]));
    		return data(s);
    	}

    	/// add-to-this operator
    	FullTensor<T>& operator+=(const FullTensor<T>& rhs) {
    		MADNESS_ASSERT(0);
    	}

    	/// add-to-this operator
    	FullTensor<T>& operator+=(const Tensor<T>& rhs) {
    		MADNESS_ASSERT(0);
    	}

        /// Inplace generalized saxpy ... this = this*alpha + other*beta
    	FullTensor<T>& gaxpy(T alpha, const FullTensor<T>& t, T beta) {
    		MADNESS_ASSERT(0);
    	}

        /// Inplace generalized saxpy ... this = this*alpha + other*beta
     	FullTensor<T>& gaxpy(T alpha, const Tensor<T>& t, T beta) {
     		MADNESS_ASSERT(0);
     	}



    	/// implement the virtual constructor idiom
    	FullTensor<T>* clone() const {
    		return new FullTensor<T>(*this);
    	}

    	/// return this tensor type
    	TensorType type() const {return TT_FULL;};

        /// return the type of the derived class for me
        std::string what_am_i() const {return "FullRank";};

        /// return the rank of this
        long rank() const {return -1;};

        /// Returns if this Tensor exists
        bool has_data() const {return size()!=0;};

    	/// the number of elements of this tensor
    	long size() const {return data.size();};

    	/// the number of dimensions (number of indices)
    	long ndim() const {return data.ndim();};

    	/// the length of each dimension (range for each index)
    	long dim(const int i) const {return data.dim(i);};

    	/// return a tensor of this
    	Tensor<T>& fullTensor() {return data;};

    	/// return a const tensor of this
    	const Tensor<T>& fullTensor() const {return data;};

    	Tensor<T> reconstructTensor() const {throw;};

    	/// compute the Frobenius norm
       	float_scalar_type normf() const {return data.normf();};

    	/// scale by a number
    	void scale(T a) {data*=a;};


    	T operator()(long i, long j, long k) const {return data(i,j,k);};
    };

    /// transform the argument SepRepTensor to FullTensor form
    template <typename T>
    void to_full_rank(SepRepTensor<T>*& arg) {

//    	print("to_full_rank");
    	if (arg->type()==TT_FULL) {
    		;
    	} else if (arg->type()==TT_3D) {
    		Tensor<T> t;
    		if (arg->has_data()) {
    			t=arg->reconstructTensor();
        		delete arg;
        		arg=new FullTensor<T>(t);
    		} else {
        		delete arg;
        		arg=new FullTensor<T>();
    		}
    	} else {
    		throw std::runtime_error("unknown TensorType in to_full_tensor");
    	}
    }

    /// transform the argument SepRepTensor to LowRankTensor form
    template <typename T>
    void to_low_rank(SepRepTensor<T>*& arg, const double& eps, const TensorType& target_type) {

//    	print("to_low_rank");

    	if (arg->type()==TT_FULL) {
			if (arg->has_data()) {
				const Tensor<T> t1=arg->fullTensor();
				delete arg;
				arg=new LowRankTensor<T>(t1,eps,target_type);
			} else {
				delete arg;
				arg=new LowRankTensor<T>();
			}
     	} else if (arg->type()==TT_3D) {
     		;
     	} else {
     		throw std::runtime_error("unknown TensorType in to_full_tensor");
     	}
     }


    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
	/// cf tensor/tensor.h
    template <class T, class Q>
    FullTensor<TENSOR_RESULT_TYPE(T,Q)> transform(const Tensor<T>& t, const FullTensor<Q>& c) {

    	FullTensor<TENSOR_RESULT_TYPE(T,Q)> result=FullTensor<TENSOR_RESULT_TYPE(T,Q)>(transform(t.data,c));
    	return result;
    }

    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
	/// cf tensor/tensor.h
    template <class T, class Q>
    LowRankTensor<TENSOR_RESULT_TYPE(T,Q)> transform(const Tensor<T>& t, const LowRankTensor<Q>& c) {

    	MADNESS_ASSERT(c.dim(0)==c.dim(1) && c.iscontiguous() and t.dim(0)==c.dim(0));
    	typedef TENSOR_RESULT_TYPE(T,Q) TQ;

    	Tensor<TQ> result(c.dim(0),c.dim(1),false);
    	for (long idim=0; idim<c.ndim(); idim++) {
    	}

    }

}

#endif /* SEPREPTENSOR_H_ */
