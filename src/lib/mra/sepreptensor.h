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

/// \file SepRepTensor.h
/// \brief Provides a tensor with taking advantage of possibly low rank

/**
 * note that a LowRankTensor might have zero rank, but is still a valid
 * tensor, and should therefore return a FullTensor with zeros in it.
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

    	SepRep<T> data;		///< the tensor data

    public:

       	/// default constructor
    	LowRankTensor<T>() : data(SepRep<T>()) {
    	};

    	/// copy constructor
    	LowRankTensor<T>(const LowRankTensor<T>& rhs) {
    		data=rhs.data;
    	};

    	/// ctor with regular Tensor
    	LowRankTensor<T>(const Tensor<T>& t, const double& eps, const TensorType& tt) {
			data=SepRep<T>(t,eps,tt);
    	}

    	/// dtor
    	~LowRankTensor<T>() {}

    	/// virtual assignment
    	LowRankTensor<T>& operator=(const LowRankTensor<T>& rhs) {};

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

    	/// copy constructor
		FullTensor<T>(const FullTensor<T>& rhs) {
			*this=rhs;
		};

		/// ctor with a regular Tensor
		FullTensor<T>(const Tensor<T>& t) {
			data=t;
		}

    	/// assignment operator
    	FullTensor<T>& operator=(const FullTensor<T>& rhs) {
    		if (this!=&rhs) {
    			data=rhs.data;
    		}
    		return *this;
    	}

    	/// assignment operator
    	FullTensor<T>& operator=(const Tensor<T>& rhs) {
    		if (&this->data!=&rhs) {
    			data=rhs;
    		}
    		return *this;
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

    	/// reimplement some basic inquiries about the shape

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

}

#endif /* SEPREPTENSOR_H_ */
