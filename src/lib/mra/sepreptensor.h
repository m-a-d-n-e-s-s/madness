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
 ****************************************************************
 * MAIN DIFFERENCES (Tensor t; GenTensor g)
 *  t=t1(s) is shallow
 *  g=g1(s) is deep
 ****************************************************************
 *
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

	template <class T> class SliceLowRankTensor;
	template <class T> class SepRepTensor;
	template <class T> class FullTensor;
	template <class T> class LowRankTensor;
	template <class T> class SliceGenTensor;

	TensorType default_tt(const int& ndim) {
		TensorType tt=TT_NONE;
		if (ndim==1) tt=FunctionDefaults<1>::get_tensor_type();
		if (ndim==2) tt=FunctionDefaults<2>::get_tensor_type();
		if (ndim==3) tt=FunctionDefaults<3>::get_tensor_type();
		if (ndim==4) tt=FunctionDefaults<4>::get_tensor_type();
		if (ndim==5) tt=FunctionDefaults<5>::get_tensor_type();
		if (ndim==6) tt=FunctionDefaults<6>::get_tensor_type();
		MADNESS_ASSERT(tt!=TT_NONE);
		return tt;
	}


	/// A GenTensor is an interface to a Tensor or a LowRankTensor

	/// This class provides a wrapper for an abstract base class SepRepTensor,
	/// which implements FullTensor (aka Tensor), and LowRankTensor (aka SepRep).
	/// Since not all operations are possible (or sensible) for LowRankTensors,
	/// only those that are are provided. There are no assignment for slices,
	/// neither for numbers nor for other GenTensors, use inplace-addition instead.
	///
	/// Behavior (in particular shallow and deep construction/assignment) should
	/// resemble the one of Tensor as far as possible:
	/// Assignments to/from other GenTensors are shallow
	/// Assignments from Tensor is deep
	/// Assignments to/from Slices are deep
	template<typename T>
	class GenTensor {

	private:

		friend class SliceGenTensor<T>;

		/// pointer to the abstract base class; implements FullTensor and LowRankTensor
		SepRepTensor<T>* _ptr;

	public:

		/// default ctor
		GenTensor() : _ptr(0) {}

		/// ctor with a SepRepTensor, shallow
		GenTensor(SepRepTensor<T>* sr) : _ptr(sr) {}

		/// copy ctor, shallow
		GenTensor(const GenTensor<T>& rhs) : _ptr(rhs._ptr) {};

		/// ctor with dimensions
		GenTensor(const std::vector<long>& dim) {
			const long ndim=dim.size();
			TensorType tt=default_tt(ndim);

			if (tt==TT_FULL) _ptr=new FullTensor<T>(dim);
			else _ptr=new LowRankTensor<T>(dim);
		}

		/// ctor with a regular Tensor, deep
		GenTensor(const Tensor<T>& rhs, double eps=0.0, TensorType tt=TT_NONE) {

			MADNESS_ASSERT(rhs.iscontiguous());
			if (eps==0.0) tt=TT_FULL;
			if (tt==TT_NONE) tt=TT_FULL;

			if (tt==TT_FULL) _ptr=new FullTensor<T> (copy(rhs));
			else if (tt==TT_3D or tt==TT_2D) _ptr=new LowRankTensor<T> (rhs,eps,tt);

		}

		/// ctor with a SliceGenTensor, deep
		GenTensor(const SliceGenTensor<T>& rhs) {
			*this=rhs;
		}

		/// dtor
		~GenTensor() {
			this->clear();
			_ptr=0;
		}

		/// release memory
		void clear() {if (not _ptr) delete _ptr; _ptr=0;};

		/// shallow assignment
		GenTensor<T>& operator=(const GenTensor<T>& rhs) {
			this->clear();
			_ptr=rhs._ptr;
			return *this;
		}

		/// deep copy
		friend GenTensor<T> copy(const GenTensor<T>& rhs) {
			return GenTensor<T> (rhs._ptr->copy_this());
		}

		/// general slicing
		SliceGenTensor<T> operator()(const std::vector<Slice>& s) {
			return SliceGenTensor<T>(*this,s);
		}

		/// general slicing, for g0 = g1(s)
		const SliceGenTensor<T> operator()(const std::vector<Slice>& s) const {
			return SliceGenTensor<T>(*this,s);
		}

		/// assignment: g0 = g1(s)
		GenTensor<T>& operator=(const SliceGenTensor<T>& rhs) {
			this->clear();
			_ptr=rhs._refGT._ptr->clone(rhs._s);
			return *this;
		}

		/// assign a number
		GenTensor<T>& operator=(const double& fac) {
			MADNESS_ASSERT(0);
		}

		/// inplace addition
		GenTensor<T>& operator+=(const GenTensor<T>& rhs) {

			MADNESS_ASSERT(this->_ptr->type()==rhs._ptr->type());
			static const std::vector<Slice> s(this->ndim(),Slice(_));
			this->_ptr->inplace_add(rhs._ptr,s,s);
			return *this;

		}

		/// returns if this GenTensor exists
		bool exists() const {
			if (_ptr) return false;
			return true;
		}

		/// returns the number of coefficients (might return zero, although tensor exists)
		size_t size() const {return _ptr->size();};

		/// returns the TensorType of this
		TensorType type() const {return _ptr->type();};

        /// return the type of the derived class for me
        std::string what_am_i() const {return _ptr->what_am_i();};

		/// returns the rank of this; if FullTensor returns -1
		long rank() const {return _ptr->rank();};

		/// returns the dimensions
		long dim(const int& i) const {return _ptr->dim(i);};

		/// returns the number of dimensions
		long ndim() const {return _ptr->ndim();};

		/// returns the Frobenius norm
		double normf() const {return _ptr->normf();};

		/// returns the trace of <this|rhs>
		T trace_conj(const GenTensor<T>& rhs) {
			MADNESS_ASSERT(0);
		}

		/// scale this by a number
		void scale(const T& fac) {_ptr->scale(fac);};

		/// returns a reference to FullTensor of this; no reconstruction
		Tensor<T>& full_tensor() const {
			return _ptr->fullTensor();
		}

		/// returns a FullTensor of this; reconstructs a LowRankTensor
		Tensor<T> reconstruct_tensor() const {
			return _ptr->reconstructTensor();
		}

		/// return a Tensor, no matter what
		Tensor<T> full_tensor_copy() const {
			if (_ptr->type()==TT_FULL) return _ptr->fullTensor();
			else if (_ptr->type()==TT_NONE) return Tensor<T>();
			else if (_ptr->type()==TT_3D or _ptr->type()==TT_2D) return _ptr->reconstructTensor();
			else {
				print(_ptr->type(),"unknown tensor type");
				MADNESS_ASSERT(0);
			}
		}

		template <typename Archive>
        void serialize(Archive& ar) {}

	private:

		/// inplace add rhs to this, provided slices: *this(s1)+=rhs(s2)
		GenTensor<T>& inplace_add(const GenTensor<T>& rhs, const std::vector<Slice>& lhs_s,
				const std::vector<Slice>& rhs_s) {

			this->_ptr->inplace_add(rhs._ptr,lhs_s,rhs_s);
			return *this;
		}


	};


	/// implements a slice of a GenTensor
	template <typename T>
	class SliceGenTensor {

	private:
		friend class GenTensor<T>;

		GenTensor<T>& _refGT;
		std::vector<Slice> _s;

		// all ctors are private, only accessible by GenTensor

		/// default ctor
		SliceGenTensor<T> () {}

		/// ctor with a GenTensor;
		SliceGenTensor<T> (const GenTensor<T>& gt, const std::vector<Slice>& s)
				: _refGT(const_cast<GenTensor<T>& >(gt))
				, _s(s) {}

	public:

		/// assignment as in g(s) = g1;
		SliceGenTensor<T>& operator=(const GenTensor<T>& rhs) {
//			_refGT._ptr->inplace_add()
			MADNESS_ASSERT(0);
			return *this;
		};

		/// assignment as in g(s) = g1(s);
		SliceGenTensor<T>& operator=(const SliceGenTensor<T>& rhs) {
			MADNESS_ASSERT(0);
			return *this;
		};

		/// inplace addition
		SliceGenTensor<T>& operator+=(const GenTensor<T>& rhs) {
			static std::vector<Slice> s(this->_refGT.ndim(),Slice(_));
			_refGT.inplace_add(rhs,this->_s,s);
			return *this;
		}

		/// inplace addition
		SliceGenTensor<T>& operator+=(const SliceGenTensor<T>& rhs) {
			_refGT.inplace_add(rhs._refGT,this->_s,rhs._s);
			return *this;
		}

		/// inplace zero-ing
		SliceGenTensor<T>& operator=(const double& number) {
			MADNESS_ASSERT(number==0.0);
			return *this;
		}

	};

    /// abstract base class for the SepRepTensor
    template<typename T>
    class SepRepTensor {

    public:

    	typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

    	/// default ctor
    	SepRepTensor() {};

    	/// implement the virtual constructor idiom; "default constructor"
    	virtual SepRepTensor* create() const =0;

    	/// implement the virtual constructor idiom; "shaped constructor"
    	virtual SepRepTensor* create(const std::vector<long>& s) const =0;

    	/// implement the virtual constructor idiom; "sliced copy constructor"
    	virtual SepRepTensor* clone(const std::vector<Slice>& s) const =0;

//    	/// implement the virtual constructor idiom; "copy constructor"
//    	virtual SepRepTensor* clone() const =0;

    	/// implement the virtual constructor idiom, deep copy
    	virtual SepRepTensor* copy_this() const =0;

    	/// add-to-this operator; in LowRank requires reallocation and is expensive!
//    	virtual SepRepTensor& operator+=(const SepRepTensor<T>* rhs)=0;

    	/// sliced assignment
//    	virtual SepRepTensor<T>* assign(const SepRepTensor<T>* rhs, const std::vector<Slice>& s) =0;

    	/// virtual destructor
    	virtual ~SepRepTensor() {};

       	/// inplace add
    	virtual SepRepTensor<T>* inplace_add(const SepRepTensor<T>* rhs,
    			const std::vector<Slice>& lhs_s, const std::vector<Slice>& rhs_s) =0;

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

    	virtual T trace_conj(const SepRepTensor<T>* rhs) const =0;

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
    	SepRep<T> _data;		///< the tensor data


       	/// default constructor (tested)
    	LowRankTensor<T>() : _data(SepRep<T>()) {
    	};

    	/// construct am empty tensor
    	LowRankTensor<T> (const std::vector<long>& s) {
    		// check consistency
    		const long ndim=s.size();
    		const long maxk=s[0];
    		for (long idim=0; idim<ndim; idim++) {
    			MADNESS_ASSERT(maxk==s[0]);
    		}
    		_data=SepRep<T>(this->type(),maxk,ndim);
    	}

    	/// copy constructor, shallow (tested)
    	LowRankTensor<T>(const LowRankTensor<T>& rhs) : _data(rhs._data) {
    	};

    	/// constructor w/ slice on rhs, deep
    	LowRankTensor<T>(const SliceLowRankTensor<T>& rhs) : _data(rhs.ref_data()(rhs._s)) {
    	};

    	/// ctor with a SepRep (shallow)
    	LowRankTensor<T>(const SepRep<T>& rhs) : _data(rhs) {
    	}

    	/// ctor with regular Tensor (tested)
    	LowRankTensor<T>(const Tensor<T>& t, const double& eps, const TensorType& tt) {
			_data=SepRep<T>(t,eps,tt);
    	}

    	/// dtor
    	~LowRankTensor<T>() {}

    	/// shallow assignment with rhs=LowRankTensor (tested)
    	LowRankTensor<T>& operator=(const LowRankTensor<T>& rhs) {
    		_data=rhs._data;
    		return *this;
    	};

    	/// assignment with rhs=SliceLowRankTensor, deep
    	LowRankTensor<T>& operator=(const SliceLowRankTensor<T>& rhs) {
    		_data=rhs.ref_data()(rhs._s);
    		return *this;
    	};

    	/// assignment to a number
    	LowRankTensor<T>& operator=(const T& a) {
    		MADNESS_ASSERT(0);
    	}

    	/// implement the virtual constructor idiom; "default constructor"
    	LowRankTensor<T>* create() const {
    		return new LowRankTensor<T>();
    	}

    	/// implement the virtual constructor idiom; "shaped constructor"
    	LowRankTensor<T>* create(const std::vector<long>& s) const {
    		return new LowRankTensor(s);
    	}

    	/// implement the virtual constructor idiom; "sliced copy constructor"
    	LowRankTensor<T>* clone(const std::vector<Slice>& s) const {
    		return new LowRankTensor<T>((this->_data)(s));
    	}


//    	/// implement the virtual constructor idiom, shallow
//    	LowRankTensor<T>* clone() const {
//    		return new LowRankTensor<T>(*this);
//    	}

    	/// implement the virtual constructor idiom, deep
    	LowRankTensor* copy_this() const {
    		return new LowRankTensor<T>(copy(*this));
    	}

    	/// deep copy
    	friend LowRankTensor<T> copy(const LowRankTensor<T>& rhs) {
    		return LowRankTensor<T>(copy(rhs._data));
    	}

       	/// return this tensor type
        TensorType type() const {return _data.tensor_type();};

        /// return the type of the derived class for me
        std::string what_am_i() const {
        	if (this->type()==TT_2D) return "LowRank-2D";
        	else if (this->type()==TT_3D) return "LowRank-3D";
        	else {
        		print("unknown tensor type",this->type());
        		MADNESS_ASSERT(0);
        	}
        };

        /// return the rank of this
        long rank() const {return _data.rank();};

        /// Returns if this Tensor exists
        bool has_data() const {return _data.is_valid();};

        /// Return the number of coefficients (valid for all SPR)
        long size() const {
        	if (!this->has_data()) return 0;
        	return _data.nCoeff();
        };

    	virtual long dim(const int) const {return _data.get_k();};

		virtual long ndim() const {return _data.dim();};

		/// reconstruct this to a full tensor
		Tensor<T> reconstructTensor() const {
			MADNESS_ASSERT(_data.is_valid());
			return _data.reconstructTensor();

		};

		/// (Don't) return a full tensor representation
		Tensor<T>& fullTensor() {
			MADNESS_EXCEPTION("No fullTensor in LowRankTensor; reconstruct first",0);
		};

		/// (Don't) return a full tensor representation
		Tensor<T>& fullTensor() const {
			MADNESS_EXCEPTION("No fullTensor in LowRankTensor; reconstruct first",0);
		};

       	/// compute the inner product
    	T trace_conj(const SepRepTensor<T>* rhs) const {
    		const LowRankTensor<T>* other=dynamic_cast<const LowRankTensor<T>* >(rhs);
    		T result=overlap(this->_data,other->_data);
    		return result;
    	}


    	T operator()(long i, long j, long k) const {throw;};

    	/// compute the Frobenius norm
    	float_scalar_type normf() const {throw;};

    	/// scale by a number
    	void scale(T a) {_data.scale(a);};

    private:

    	/// inplace add: this(lhs_s) += rhs(rhs_s)
    	LowRankTensor<T>* inplace_add(const SepRepTensor<T>* rhs,
    			const std::vector<Slice>& lhs_s, const std::vector<Slice>& rhs_s) {

    		// some checks
    		MADNESS_ASSERT(this->type()==rhs->type());

    		const LowRankTensor<T>* rhs2=dynamic_cast<const LowRankTensor<T>* > (rhs);
    		this->_data.inplace_add(rhs2->_data,lhs_s,rhs_s);
    		return this;
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

		/// ctor with an empty Tensor
		FullTensor<T>(const std::vector<long>& s) {
			data=Tensor<T>(s);
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
    	FullTensor<T> operator()(const std::vector<Slice>& s) const {
//    		return SliceTensor<T> (data,&(s[0]));
    		return data(s);
    	}

    	/// sliced assignment
    	FullTensor<T>* assign(const SepRepTensor<T>* rhs, const std::vector<Slice>& s) {
    		const FullTensor<T>* rhs2=dynamic_cast<const FullTensor<T>* > (rhs);
    		*this=(*rhs2)(s);
    		return this;
    	}

        /// Inplace generalized saxpy ... this = this*alpha + other*beta
    	FullTensor<T>& gaxpy(T alpha, const FullTensor<T>& t, T beta) {
    		MADNESS_ASSERT(0);
    	}

        /// Inplace generalized saxpy ... this = this*alpha + other*beta
     	FullTensor<T>& gaxpy(T alpha, const Tensor<T>& t, T beta) {
     		MADNESS_ASSERT(0);
     	}


       	/// implement the virtual constructor idiom; "default constructor"
       	virtual FullTensor<T>* create() const {
       		return new FullTensor<T>();
       	}

       	/// implement the virtual constructor idiom; "shaped constructor"
       	virtual FullTensor<T>* create(const std::vector<long>& s) const {
       		return new FullTensor(s);
       	}

    	/// implement the virtual constructor idiom
    	FullTensor<T>* clone() const {
    		return new FullTensor<T>(*this);
    	}

    	/// implement the virtual constructor idiom, "sliced copy ctor"
    	FullTensor<T>* clone(const std::vector<Slice>& s) const {
    		return new FullTensor<T>((*this)(s));
    	}

    	/// implement the virtual constructor idiom, deep
    	FullTensor<T>* copy_this() const {
    		return new FullTensor<T>(copy(*this));
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

       	/// compute the inner product
    	T trace_conj(const SepRepTensor<T>* rhs) const {
    		const FullTensor<T>* other=dynamic_cast<const FullTensor<T>* >(rhs);
    		return this->data.trace_conj(other->data);
    	}

    	/// scale by a number
    	void scale(T a) {data*=a;};

        template <typename Archive>
        void serialize(Archive& ar) {}

    	T operator()(long i, long j, long k) const {return data(i,j,k);};

    	friend FullTensor<T> copy(const FullTensor<T>& rhs) {
    		return FullTensor<T>(copy(rhs.data));
    	}

    private:

    	/// inplace add: this(lhs_s) += rhs(rhs_s)
    	FullTensor<T>* inplace_add(const SepRepTensor<T>* rhs,
    			const std::vector<Slice>& lhs_s, const std::vector<Slice>& rhs_s) {

    		// some checks
    		MADNESS_ASSERT(this->type()==rhs->type());

    		const FullTensor<T>* rhs2=dynamic_cast<const FullTensor<T>* >(rhs);
    		this->data(lhs_s)+=rhs2->data(rhs_s);
    		return this;
    	}

    };

    /// transform the argument SepRepTensor to FullTensor form
    template <typename T>
    void to_full_rank(GenTensor<T>& arg) {

//    	print("to_full_rank");
    	if (arg.type()==TT_FULL) {
    		;
    	} else if (arg.type()==TT_3D) {
    		Tensor<T> t;
    		if (arg.exists()) {
    			t=arg.reconstruct_tensor();
        		arg.clear();
        		arg=new FullTensor<T>(t);
    		} else {
        		arg.clear();
        		arg=new FullTensor<T>();
    		}
    	} else {
    		throw std::runtime_error("unknown TensorType in to_full_tensor");
    	}
    }

    /// transform the argument SepRepTensor to LowRankTensor form
    template <typename T>
    void to_low_rank(GenTensor<T>& arg, const double& eps, const TensorType& target_type) {

//    	print("to_low_rank");

    	if (arg.type()==TT_FULL) {
			if (arg.exists()) {
				const Tensor<T> t1=arg.full_tensor();
				arg.clear();
				arg=(new LowRankTensor<T>(t1,eps,target_type));
			} else {
				arg=new LowRankTensor<T>();
			}
     	} else if (arg.type()==TT_3D) {
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
