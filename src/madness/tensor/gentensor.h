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

#ifndef GENTENSOR_H_
#define GENTENSOR_H_


/*!
 *
 * \file GenTensor.h
 * \brief Provides a tensor with taking advantage of possibly low rank
 *
 ****************************************************************
 * MAIN DIFFERENCES (Tensor t; GenTensor g)
 *  t=t1(s) is shallow
 *  g=g1(s) is deep
 ****************************************************************
 *
 * a GenTensor is a generalized form of a Tensor
 * for now only little functionality shall be implemented; feel free to extend
 * a consequence is that we (usually) can't directly access individual
 * matrix elements
 * note that a GenTensor might have zero rank, but is still a valid
 * tensor, and should therefore return a FullTensor with zeros in it.
 *
 * \par Slicing in GenTensors:
 *
 * A GenTensor differs from a Tensor in that we can't directly access
 * individual matrix elements, and thus can't directly assign or manipulate slices
 * as lvalues. For rvalues we simply provide a slices of the constituent vectors
 * in SRConf, which is a valid GenTensor by itself
 * \code
 *  lhs = rhs(s)
 * \endcode
 *
 * Manipulations of slices of a GenTensor are heavily restricted, but should
 * cover the most important cases:
 * - assignment to zero; done by inplace subtraction of the slice
 *    \code
 *    GenTensor lrt1(t);
 *    lrt1(s) = 0.0;
 *    \endcode
 * - in-place addition
 *    \code
 *    GenTensor lrt1(3,3,3,3), lrt2(t);
 *    lrt1(s) += lrt2;
 *    \endcode
 *
 * Note that *all* of these operation increase the rank of lhs
 *
 * \par Addition in GenTensors
 *
 * Addition in a GenTensor is a fairly complicated issue, and there are several
 * algorithms you can use, each with certain pros and cons
 *
 * - append()
 *   plain concatenating of the configurations will increase the rank and will
 *   introduce linear dependencies and redundancies. Also, you might have to
 *   reallocate memory and copy a lot of data around. However, no accuracy is lost
 *   and it is fairly fast. Final rank reduction might be expensive, since you
 *   might have accumulated a huge amount of terms.
 *
 * - low_rank_add_sequential()
 *   only for SVD (2-way decomposition)
 *   take the 2nd vector as a basis of a vector space, project the overlap of the
 *   old and new terms on the lhs basis, perform modified Gram-Schmidt
 *   orthogonalization on the rhs basis and increase the basis if necessary.
 *   This will require the left hand side to be right-orthonormal, and the caller is
 *   responsible for doing that. If a new GenTensor is created and only additions
 *   using this method are performed the tensor will automatically be orthonormal.
 *   Cost depends on the number of terms of the rhs and the lhs
 *
 * - addition of slices
 *   as of now we can add slices only using append()
 *
 * - addition in full rank form
 *   at times it might be sensible to accumulate your result in a full rank tensor,
 *   and after collecting all the contribution you can transform it into a low
 *   rank tensor. This will also maintain "all" digits, and cost depends not on
 *   the number of terms on the lhs.
 */

#include <madness/tensor/tensor.h>
#include <madness/tensor/srconf.h>
#include <madness/tensor/tensortrain.h>
#include <stdexcept>


// you can use low-rank tensors only when you use gentensor
#if HAVE_GENTENSOR
#define USE_LRT
#else
#undef USE_LRT
#endif

namespace madness {

	// a forward declaration
	template <class T> class SliceGenTensor;
	template <class T> class GenTensor;


	/// TensorArgs holds the arguments for creating a LowRankTensor
	struct TensorArgs {
		double thresh;
		TensorType tt;
        TensorArgs() : thresh(-1.0), tt(TT_NONE) {}
		TensorArgs(const double& thresh1, const TensorType& tt1)
			: thresh(thresh1)
			, tt(tt1) {
		}
		static std::string what_am_i(const TensorType& tt) {
			if (tt==TT_2D) return "TT_2D";
			if (tt==TT_FULL) return "TT_FULL";
            if (tt==TT_TENSORTRAIN) return "TT_TENSORTRAIN";
			return "unknown tensor type";
		}
		template <typename Archive>
		void serialize(const Archive& ar) {
		    int i=int(tt);
		    ar & thresh & i;
		    tt=TensorType(i);
		}
	};


    /// true if GenTensor has zero rank or no data
    template<typename T>
    bool has_zero_rank(const GenTensor<T>& g) {
    	return ((g.rank()==0) or (g.tensor_type()==TT_NONE));
    }


#if !HAVE_GENTENSOR

	template <typename T>
	class GenTensor : public Tensor<T> {

	public:

		GenTensor<T>() : Tensor<T>() {}

		GenTensor<T>(const Tensor<T>& t1) : Tensor<T>(t1) {}
		GenTensor<T>(const Tensor<T>& t1, const TensorArgs& targs) : Tensor<T>(t1) {}
		GenTensor<T>(const Tensor<T>& t1, double eps, const TensorType tt) : Tensor<T>(t1) {}
		GenTensor<T>(const TensorType tt): Tensor<T>() {}
		GenTensor<T>(std::vector<long> v, const TensorType& tt) : Tensor<T>(v) {}
		GenTensor<T>(std::vector<long> v, const TensorArgs& targs) : Tensor<T>(v) {}
		GenTensor<T>(const SRConf<T>& sr1) : Tensor<T>() {MADNESS_EXCEPTION("no ctor with SRConf: use HAVE_GENTENSOR",1);}

        /// Type conversion makes a deep copy
        template <class Q> operator GenTensor<Q>() const { // type conv => deep copy
            Tensor<Q> result = Tensor<Q>(this->_ndim,this->_dim,false);
            BINARY_OPTIMIZED_ITERATOR(Q, result, const T, (*this), *_p0 = (Q)(*_p1));
            return result;
        }

        GenTensor convert(const TensorArgs& targs) const {return copy(*this);}

		GenTensor<T> reconstruct_tensor() const {return *this;}
		GenTensor<T> full_tensor() const {return *this;}
		GenTensor<T>& full_tensor() {return *this;}
        GenTensor<T> full_tensor_copy() const {return *this;}
        GenTensor<T> full_tensor_copy() {return *this;}

        bool has_data() const {return this->size()>0;};
        bool has_no_data() const {return not has_data();};
		long rank() const {return -1;}
		double svd_normf() const {return this->normf();}
		size_t real_size() const {return this->size();}

        void reduce_rank(const double& eps) {return;};
        void normalize() {return;}

        std::string what_am_i() const {return "GenTensor, aliased to Tensor";};
		TensorType tensor_type() const {return TT_FULL;}

		void add_SVD(const GenTensor<T>& rhs, const double& eps) {*this+=rhs;}

		SRConf<T> config() const {MADNESS_EXCEPTION("no SRConf in complex GenTensor",1);}
        SRConf<T> get_configs(const int& start, const int& end) const {MADNESS_EXCEPTION("no SRConf in complex GenTensor",1);}

        template<typename Q>
        GenTensor<T> general_transform(const Tensor<Q> c[]) const {

            return madness::general_transform(*static_cast<const Tensor<T>*>(this),c);
        }


		/// return the additional safety for rank reduction
		static double fac_reduce() {return -1.0;};

	};

    template <class T>
	GenTensor<T> reduce(std::list<GenTensor<T> >& addends, double eps, bool are_optimal=false) {
    	typedef typename std::list<GenTensor<T> >::iterator iterT;
    	GenTensor<T> result=copy(addends.front());
    	for (iterT it=++addends.begin(); it!=addends.end(); ++it) {
    		result+=*it;
    	}
	    return result;
    }

    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)

    /// \ingroup tensor
    template <class T>
    GenTensor<T> outer(const GenTensor<T>& left, const GenTensor<T>& right,
            const TensorArgs final_tensor_args) {
        return madness::outer(static_cast<Tensor<T> >(left),static_cast<Tensor<T> >(right));
    }

    /// Outer product ... result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)

     /// \ingroup tensor
     template <class T>
     GenTensor<T> outer(const Tensor<T>& left, const Tensor<T>& right,
             const TensorArgs final_tensor_args) {
         return madness::outer(left,right);
     }


	namespace archive {
	/// Serialize a tensor
	template <class Archive, typename T>
	struct ArchiveStoreImpl< Archive, GenTensor<T> > {
		static void store(const Archive& s, const GenTensor<T>& t) {
			if (t.iscontiguous()) {
				s & t.size() & t.id();
				if (t.size()) s & t.ndim() & wrap(t.dims(),TENSOR_MAXDIM) & wrap(t.ptr(),t.size());
			}
			else {
				s & copy(t);
			}
		};
	};


	/// Deserialize a tensor ... existing tensor is replaced
	template <class Archive, typename T>
	struct ArchiveLoadImpl< Archive, GenTensor<T> > {
		static void load(const Archive& s, GenTensor<T>& t) {
			long sz = 0l, id =0l;
			s & sz & id;
			if (id != t.id()) throw "type mismatch deserializing a tensor";
			if (sz) {
				long _ndim=0l, _dim[TENSOR_MAXDIM];
				s & _ndim & wrap(_dim,TENSOR_MAXDIM);
				t = Tensor<T>(_ndim, _dim, false);
				if (sz != t.size()) throw "size mismatch deserializing a tensor";
				s & wrap(t.ptr(), t.size());
			}
			else {
				t = Tensor<T>();
			}
		};
	};

	}

#else

#ifndef USE_LRT

	/// A GenTensor is a generalized tensor, possibly in a low rank representation

	/// Since not all operations are possible (or sensible) for low rank tensors,
	/// only those that are are provided. There are no assignment for slices,
	/// neither for numbers nor for other GenTensors, use inplace-addition instead.
	///
	/// Behavior (in particular shallow and deep construction/assignment) should
	/// resemble the one of Tensor as far as possible:
	/// assignments/construction to/from other GenTensors are shallow
	/// assignments/construction from Tensor is deep
	/// assignments/construction to/from Slices are deep
	template<typename T>
	class GenTensor {

	private:

		friend class SliceGenTensor<T>;

        /// C++ typename of the real type associated with a complex type.
        typedef typename TensorTypeData<T>::scalar_type scalar_type;

        /// C++ typename of the floating point type associated with scalar real type
        typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;


		typedef SRConf<T> configT;
		typedef Tensor<T> tensorT;
		typedef GenTensor<T> gentensorT;
		typedef std::shared_ptr<configT> sr_ptr;

		/// pointer to the low rank tensor
		sr_ptr _ptr;

		/// the machine precision
		static double machinePrecision() {return 1.e-14;}

		/// safety for rank reduction
		static double facReduce() {return 1.e-3;}

	public:

		/// empty ctor
		GenTensor() : _ptr() {
		}

		/// copy ctor, shallow
//		GenTensor(const GenTensor<T>& rhs) : _ptr(rhs._ptr) { // DON'T DO THIS: USE_COUNT BLOWS UP
		GenTensor(const GenTensor<T>& rhs) : _ptr() {
			if (rhs.has_data()) _ptr=rhs._ptr;
		};

		/// ctor with dimensions
		GenTensor(const std::vector<long>& dim, const TensorType tt) : _ptr() {

    		// check consistency
    		const long ndim=dim.size();
    		const long maxk=dim[0];
    		for (long idim=0; idim<ndim; idim++) {
    			MADNESS_ASSERT(maxk==dim[0]);
    		}

   			_ptr=sr_ptr(new configT(dim.size(),dim[0],tt));

		}

		/// ctor with dimensions
		GenTensor(const std::vector<long>& dim, const TensorArgs& targs) : _ptr() {

			// check consistency
    		const long ndim=dim.size();
    		const long maxk=dim[0];
    		for (long idim=0; idim<ndim; idim++) {
    			MADNESS_ASSERT(maxk==dim[0]);
    		}

   			_ptr=sr_ptr(new configT(dim.size(),dim[0],targs.tt));

		}

		/// ctor with dimensions
		GenTensor(const TensorType& tt, const unsigned int& k, const unsigned int& dim) {
   			_ptr=sr_ptr(new configT(dim,k,tt));
		}

		/// ctor with a regular Tensor and arguments, deep
		GenTensor(const Tensor<T>& rhs, const TensorArgs& targs) {

			// fast return if possible
			if (not rhs.has_data()) {
				_ptr.reset();
				return;
			}

			MADNESS_ASSERT(rhs.ndim()>0);
			for (long idim=0; idim<rhs.ndim(); idim++) {
			    MADNESS_ASSERT(rhs.dim(0)==rhs.dim(idim));
			}

			_ptr=sr_ptr(new configT(rhs.ndim(),rhs.dim(0),targs.tt));

			// direct reduction on the polynomial values on the Tensor
			TensorType ttype=tensor_type();
			if (ttype==TT_2D) {

#if 1
				Tensor<T> U,VT;
				Tensor< typename Tensor<T>::scalar_type > s;

				// add some extra dimensions if this is supposedly NS form:
				// separate sum and wavelet coeffs
				if (rhs.dim(0)%2==0) {
					std::vector<long> dims(rhs.ndim()*2);
					for (int i=0; i<rhs.ndim(); ++i) {
						int k=rhs.dim(i);
						dims[2*i]=k/2;
						dims[2*i+1]=2;
					}
					TensorTrain<T> tt(rhs,targs.thresh*facReduce(),dims);
					// fuse sum and wavelet coeffs back together
					for (int i=0; i<rhs.ndim(); ++i) tt.fusedim(i);
					// fuse dimensions into particles 1 and 2
					tt.two_mode_representation(U,VT,s);

				} else {
					TensorTrain<T> tt(rhs,targs.thresh*facReduce());
					tt.two_mode_representation(U,VT,s);
				}

				const long r=VT.dim(0);
				const long nd=VT.ndim();
				if (r==0) {
                    _ptr=sr_ptr(new configT(dim(),get_k(),ttype));
				} else {
					MADNESS_ASSERT(U.dim(nd-1)==r);
					Tensor<T> UU=U.reshape(U.size()/r,r);
					_ptr=sr_ptr(new configT(s, copy(transpose(UU)), VT.reshape(r,VT.size()/r),
				                          dim(), get_k()));
					this->normalize();
				}
#else
				// adapt form of values
				MADNESS_ASSERT(rhs.iscontiguous());
				std::vector<long> d(_ptr->dim_eff(),_ptr->kVec());
				Tensor<T> values_eff=rhs.reshape(d);

				this->computeSVD(targs.thresh,values_eff);
#endif
			} else if (ttype==TT_FULL){
				_ptr.reset(new configT(copy(rhs)));
			} else {
				MADNESS_EXCEPTION("unknown tensor type in GenTensor(tensorT,targs)",0);
			}
		}

		/// ctor with a regular Tensor and arguments, deep
		GenTensor(const Tensor<T>& rhs, const double& thresh, const TensorType& tt) {
			*this=gentensorT(rhs,TensorArgs(thresh,tt));
		}

		/// ctor with a SliceGenTensor, deep
		GenTensor(const SliceGenTensor<T>& rhs) : _ptr() {
			*this=rhs;
		}

		/// dtor
		~GenTensor() {
			this->clear();
		}

		/// shallow assignment operator: g0 = g1
		gentensorT& operator=(const gentensorT& rhs) {
			if (this != &rhs) _ptr=rhs._ptr;
			return *this;
		}

		/// deep assignment with slices: g0 = g1(s)
		GenTensor& operator=(const SliceGenTensor<T>& rhs) {
			*this=rhs._refGT.copy_slice(rhs._s);
			return *this;
		}

		/// deep copy of rhs by deep copying rhs.configs
		friend gentensorT copy(const gentensorT& rhs) {
			if (rhs._ptr) return gentensorT(copy(*rhs._ptr));
			return gentensorT();
		}


        /// Type conversion makes a deep copy
        template <class Q> operator GenTensor<Q>() const { // type conv => deep copy
            GenTensor<Q> result;
            MADNESS_EXCEPTION("no type conversion in GenTensor -- use NO_GENTENSOR",1);
            return result;
        }

		/// return some of the terms of the SRConf (start,..,end), inclusively
		/// shallow copy
		const GenTensor get_configs(const int& start, const int& end) const {
			return gentensorT(config().get_configs(start,end));
		}

		/// general slicing, shallow; for temporary use only!
		SliceGenTensor<T> operator()(const std::vector<Slice>& s) {
			return SliceGenTensor<T>(*this,s);
		}

		/// general slicing, for g0 = g1(s), shallow, for temporary use only!
		const SliceGenTensor<T> operator()(const std::vector<Slice>& s) const {
			return SliceGenTensor<T>(*this,s);
		}

		/// assign a number
		GenTensor& operator=(const T& fac) {
            if (this->tensor_type()==TT_FULL) full_tensor()=fac;
            else config()=fac;
            return *this;
		}

	private:
		/// disable direct construction with a tensor since it's not shallow
		GenTensor(const Tensor<T>& t) {
			print("in wrong constructor");
		}
	public:

		/// ctor w/ configs, shallow (indirectly, via vector_)
		explicit GenTensor(const SRConf<T>& config) : _ptr(new configT(config)) {
		}

	private:
		/// return a slice of this (deep copy)
		gentensorT copy_slice(const std::vector<Slice>& s) const {

            // fast return if possible
            if (this->has_no_data()) {
                int k_new=s[0].end-s[0].start+1;
                return gentensorT (this->tensor_type(),k_new,this->dim());
            }

			// consistency check
			MADNESS_ASSERT(s.size()==this->dim());
			MADNESS_ASSERT(s[0].step==1);

			// fast return for full rank tensors
			if (tensor_type()==TT_FULL) {
				tensorT a=copy(full_tensor()(s));
				return gentensorT(configT(a));
			} else if (tensor_type()==TT_2D) {
			    return gentensorT(config().copy_slice(s));
			} else {
			    MADNESS_EXCEPTION("fault in GenTensor::copy_slice",1);
			    return gentensorT();
			}
		}

	public:

		/// access the tensor values, iff this is full rank representation
		const tensorT& full_tensor() const {
			MADNESS_ASSERT(tensor_type()==TT_FULL);
			return _ptr->ref_vector(0);
		}

		/// access the tensor values, iff this is full rank representation
		tensorT& full_tensor() {
			MADNESS_ASSERT(tensor_type()==TT_FULL);
			return _ptr->ref_vector(0);
		}

		/// add another SepRep to this one
		gentensorT& operator+=(const gentensorT& rhs) {

		    if (rhs.has_no_data()) return *this;
		    if (this->has_no_data()) {
		    	*this=copy(rhs);
		    	return *this;
		    }
		    this->gaxpy(1.0,rhs,1.0);
			return *this;
		}

		/// inplace subtraction
		gentensorT& operator-=(const gentensorT& rhs) {

		    if (rhs.has_no_data()) return *this;
		    if (this->has_no_data()) {
		    	*this=copy(rhs);
		    	return *this;
		    }
		    this->gaxpy(1.0,rhs,-1.0);
			return *this;
		}

		/// inplace addition
		gentensorT& operator+=(const SliceGenTensor<T>& rhs) {
			const std::vector<Slice> s(this->ndim(),Slice(0,get_k()-1,1));
			this->_ptr->inplace_add(*rhs._refGT._ptr,s,rhs._s,1.0,1.0);
			return *this;
		}

		/// inplace subtraction
		gentensorT& operator-=(const SliceGenTensor<T>& rhs) {
			const std::vector<Slice> s(this->ndim(),Slice(0,get_k()-1,1));
			this->_ptr->inplace_add(*rhs._refGT._ptr,s,rhs._s,1.0,-1.0);
			return *this;
		}

		/// multiply with a number
		template<typename Q>
		GenTensor<TENSOR_RESULT_TYPE(T,Q)> operator*(const Q& x) const {
			GenTensor<TENSOR_RESULT_TYPE(T,Q)> result(copy(*this));
        	result.scale(x);
        	return result;
		}

	    /// Inplace generalized saxpy ... this = this*alpha + other*beta
	    gentensorT& gaxpy(const T alpha, const gentensorT& rhs, const T beta) {
			MADNESS_ASSERT(this->tensor_type()==rhs.tensor_type());

	    	if (tensor_type()==TT_FULL) {
                full_tensor().gaxpy(alpha,rhs.full_tensor(),beta);
                return *this;
            }
	    	if (not (alpha==1.0)) this->scale(alpha);
	    	rhs.append(*this,beta);
	    	return *this;
	    }

		/// multiply with a scalar
	    template<typename Q>
	    GenTensor<TENSOR_RESULT_TYPE(T,Q)>& scale(const Q& dfac) {
			if (!_ptr) return *this;
			if (tensor_type()==TT_FULL) {
				full_tensor().scale(dfac);
			} else {
				_ptr->scale(dfac);
			}
			return *this;
		};

		/// normalize
		void normalize() {
			if (tensor_type()==TT_2D) _ptr->normalize();
		};

		void fillrandom(const int r=1) {
			if (tensor_type()==TT_FULL) full_tensor().fillrandom();
			else _ptr->fillWithRandom(r);
		}

		/// do we have data? note difference to SRConf::has_data() !
		bool has_data() const {
		    if (_ptr) return true;
		    return false;
		}

		/// do we have data?
		bool has_no_data() const {return (!has_data());}

		/// return the separation rank
		long rank() const {
			if (_ptr) return _ptr->rank();
			else return 0;
		};

		/// return the dimension
		unsigned int dim() const {return _ptr->dim();};

		/// returns the dimensions
		long dim(const int& i) const {return _ptr->get_k();};

		/// returns the number of dimensions
		long ndim() const {
			if (_ptr) return _ptr->dim();
			return -1;
		};

		/// return the polynomial order
		unsigned int get_k() const {return _ptr->get_k();};

		/// returns the TensorType of this
		TensorType tensor_type() const {
			if (_ptr) return _ptr->type();
			return TT_NONE;
		};

        /// return the type of the derived class for me
        std::string what_am_i() const {return TensorArgs::what_am_i(_ptr->type());};

		/// returns the number of coefficients (might return zero, although tensor exists)
		size_t size() const {
			if (_ptr) return _ptr->nCoeff();
			return 0;
		};

		/// returns the number of coefficients (might return zero, although tensor exists)
		size_t real_size() const {
			if (_ptr) return _ptr->real_size()+sizeof(*this);
			return 0;
		};

		/// returns the Frobenius norm
		float_scalar_type normf() const {
			if (has_no_data()) return 0.0;
			return config().normf();
		};

        /// returns the Frobenius norm; expects full rank or SVD!
        double svd_normf() const {
            if (has_no_data()) return 0.0;
            if (tensor_type()==TT_2D) return config().svd_normf();
            return config().normf();
        };

        /// returns the trace of <this|rhs>
		T trace(const GenTensor<T>& rhs) const {
			return this->trace_conj(rhs);
		}

        /// returns the trace of <this|rhs>
		template<typename Q>
		TENSOR_RESULT_TYPE(T,Q) trace_conj(const GenTensor<Q>& rhs) const {
            if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in GenTensor, sorry",1);
            if (TensorTypeData<Q>::iscomplex) MADNESS_EXCEPTION("no complex trace in GenTensor, sorry",1);

            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
			// fast return if possible
			if ((this->rank()==0) or (rhs.rank()==0)) return resultT(0.0);

			MADNESS_ASSERT(compatible(*this,rhs));
			MADNESS_ASSERT(this->tensor_type()==rhs.tensor_type());

			return overlap(*(this->_ptr),*rhs._ptr);
		}

        /// returns the trace of <this|rhs>
		template<typename Q>
		TENSOR_RESULT_TYPE(T,Q) trace_conj(const Tensor<Q>& rhs) const {
            if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in GenTensor, sorry",1);
            if (TensorTypeData<Q>::iscomplex) MADNESS_EXCEPTION("no complex trace in GenTensor, sorry",1);

            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
			// fast return if possible
			if ((this->rank()==0)) return resultT(0.0);

			// fast return if this is a full tensor
			// otherwise reconstruct this to a full tensor, since it's presumably
			// faster than turning the full tensor rhs into a low-rank approximation
			if (this->tensor_type()==TT_FULL) return full_tensor().trace_conj(rhs);
			else return full_tensor_copy().trace_conj(rhs);

		}

        /// Inplace multiply by corresponding elements of argument Tensor
        GenTensor<T>& emul(const GenTensor<T>& t) {
            MADNESS_ASSERT(t.tensor_type()==this->tensor_type());

            if (this->tensor_type()==TT_FULL) full_tensor().emul(t.full_tensor());
            else config().emul(t.config());
            return *this;
        }

        GenTensor convert(const TensorArgs& targs) const {
            GenTensor<T> result=copy(*this);
            change_tensor_type(result,targs);
            return result;
        }


        /// return a Tensor, no matter what
		Tensor<T> full_tensor_copy() const {
			const TensorType tt=tensor_type();
			if (tt==TT_NONE) return Tensor<T>();
			else if (tt==TT_2D) return this->reconstruct_tensor();
			else if (tt==TT_FULL) {
				return copy(full_tensor());
			} else {
				print(TensorArgs::what_am_i(tt),"unknown tensor type");
				MADNESS_ASSERT(0);
			}
		}

		/// reduce the rank of this
		void reduce_rank(const double& eps) {

			if (rank()==0) return;

			// direct reduction on the polynomial values on the Tensor
			if (tensor_type()==TT_FULL or tensor_type()==TT_NONE) {
				return;
			} else if (this->tensor_type()==TT_2D) {
				config().divide_and_conquer_reduce(eps*facReduce());
			} else {
				MADNESS_EXCEPTION("unknown tensor type in GenTensor::reduceRank()",0);
			}
			MADNESS_ASSERT(this->_ptr->has_structure() or this->rank()==0);
		}

		/// print this' coefficients
		void printCoeff(const std::string title) const {
			print("printing SepRep",title);
			print(_ptr->weights_);
			for (unsigned int idim=0; idim<this->_ptr->dim_eff(); idim++) {
				print("coefficients for dimension",idim);
				print(_ptr->vector_[idim]);
			}
		}

		/// reconstruct this to return a full tensor
		Tensor<T> reconstruct_tensor() const {

			if (tensor_type()==TT_FULL) return full_tensor();
			else if (tensor_type()==TT_2D) return config().reconstruct();
			else {
			    MADNESS_EXCEPTION("you should not be here",1);
			}
            return Tensor<T>();
		}

		/// append this to rhs, shape must conform
		void append(gentensorT& rhs, const T fac=1.0) const {
			rhs.config().append(*this->_ptr,fac);
		}

		/// add SVD
		void add_SVD(const gentensorT& rhs, const double& thresh) {
			if (rhs.has_no_data()) return;
            if (has_no_data()) {
                *this=rhs;
                return;
            }
			if (tensor_type()==TT_FULL or tensor_type()==TT_NONE) {
				this->full_tensor()+=rhs.full_tensor();
				return;
			}
			config().add_SVD(rhs.config(),thresh*facReduce());
		}

	    /// check compatibility
		friend bool compatible(const gentensorT& rhs, const gentensorT& lhs) {
			return ((rhs.tensor_type()==lhs.tensor_type()) and (rhs.get_k()==lhs.get_k())
					and (rhs.dim()==lhs.dim()));
		};

		/// transform the Legendre coefficients with the tensor
		gentensorT transform(const Tensor<T> c) const {
//			_ptr->make_structure();
		    if (has_no_data()) return gentensorT();
			MADNESS_ASSERT(_ptr->has_structure());
			return gentensorT (this->_ptr->transform(c));
		}

		/// transform the Legendre coefficients with the tensor
		template<typename Q>
		gentensorT general_transform(const Tensor<Q> c[]) const {
//		    this->_ptr->make_structure();
		    if (has_no_data()) return gentensorT();
            MADNESS_ASSERT(_ptr->has_structure());
			return gentensorT (this->config().general_transform(c));
		}

		/// inner product
		gentensorT transform_dir(const Tensor<T>& c, const int& axis) const {
//            this->_ptr->make_structure();
            MADNESS_ASSERT(_ptr->has_structure());
            return GenTensor<T>(this->_ptr->transform_dir(c,axis));
		}

		/// return a reference to the SRConf
		const SRConf<T>& config() const {return *_ptr;}

		/// return a reference to the SRConf
		SRConf<T>& config() {return *_ptr;}

		/// return the additional safety for rank reduction
		static double fac_reduce() {return facReduce();};

	private:

		/// release memory
		void clear() {_ptr.reset();};

		/// same as operator+=, but handles non-conforming vectors (i.e. slices)
		void inplace_add(const gentensorT& rhs, const std::vector<Slice>& lhs_s,
				const std::vector<Slice>& rhs_s, const double alpha, const double beta) {

			// fast return if possible
			if (rhs.has_no_data() or rhs.rank()==0) return;

			if (this->has_data()) MADNESS_ASSERT(this->tensor_type()==rhs.tensor_type());

			// no fast return possible!!!
			//			if (this->rank()==0) {
			//				// this is a deep copy
			//				*this=rhs(rhs_s);
			//				this->scale(beta);
			//				return;
			//			}

			if (tensor_type()==TT_FULL) {
				full_tensor()(lhs_s).gaxpy(alpha,rhs.full_tensor()(rhs_s),beta);
			} else {
//				rhs._ptr->make_structure();
//				_ptr->make_structure();
				MADNESS_ASSERT(_ptr->has_structure());
                MADNESS_ASSERT(rhs._ptr->has_structure());
				this->_ptr->inplace_add(*rhs._ptr,lhs_s,rhs_s, alpha, beta);
			}
		}

		/// reduce the rank using SVD
		void computeSVD(const double& eps,const Tensor<T>& values_eff) {

			// SVD works only with matrices (2D)
			MADNESS_ASSERT(values_eff.ndim()==2);
			MADNESS_ASSERT(this->tensor_type()==TT_2D);

			// fast return if possible
			if (values_eff.normf()<eps*facReduce()) {
				_ptr=sr_ptr(new configT(_ptr->dim(),_ptr->get_k(),tensor_type()));
				return;
			}

			// output from svd
			Tensor<T> U;
			Tensor<T> VT;
			Tensor< typename Tensor<T>::scalar_type > s;

			svd(values_eff,U,s,VT);

			// find the maximal singular value that's supposed to contribute
			// singular values are ordered (largest first)
			const double thresh=eps*facReduce();
			long i=SRConf<T>::max_sigma(thresh,s.dim(0),s);

			// convert SVD output to our convention
			if (i>=0) {
			    // copy to have contiguous and tailored singular vectors
				_ptr=sr_ptr(new configT(copy(s(Slice(0,i))), copy(transpose(U(_,Slice(0,i)))),
				        copy(VT(Slice(0,i),_)), dim(), get_k()));
                MADNESS_ASSERT(this->_ptr->get_k()==this->_ptr->vector_[0].dim(1));
                MADNESS_ASSERT(this->_ptr->rank()==this->_ptr->vector_[0].dim(0));
                MADNESS_ASSERT(this->_ptr->rank()==this->_ptr->weights_.dim(0));

			} else {
				_ptr=sr_ptr(new configT(dim(),get_k(),tensor_type()));
			}
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
			print("You don't want to assign to a SliceGenTensor; use operator+= instead");
			MADNESS_ASSERT(0);
			return *this;
		};

		/// assignment as in g(s) = g1(s);
		SliceGenTensor<T>& operator=(const SliceGenTensor<T>& rhs) {
			print("You don't want to assign to a SliceGenTensor; use operator+= instead");
			MADNESS_ASSERT(0);
			return *this;
		};

		/// inplace addition
		SliceGenTensor<T>& operator+=(const GenTensor<T>& rhs) {
			std::vector<Slice> s(this->_refGT.ndim(),Slice(_));
			_refGT.inplace_add(rhs,this->_s,s,1.0,1.0);
			return *this;
		}

		/// inplace subtraction
		SliceGenTensor<T>& operator-=(const GenTensor<T>& rhs) {
			std::vector<Slice> s(this->_refGT.ndim(),Slice(_));
			_refGT.inplace_add(rhs,this->_s,s,1.0,-1.0);
			return *this;
		}

		/// inplace addition
		SliceGenTensor<T>& operator+=(const SliceGenTensor<T>& rhs) {
			_refGT.inplace_add(GenTensor<T>(*rhs._refGT._ptr),this->_s,rhs._s,1.0,1.0);
			return *this;
		}

		/// inplace zero-ing
		SliceGenTensor<T>& operator=(const T& number) {
			MADNESS_ASSERT(number==T(0.0));
			if (this->_refGT.tensor_type()==TT_FULL) {
				_refGT.full_tensor()(_s)=0.0;
			} else {
				const GenTensor<T>& tmp=(this->_refGT);
				_refGT.inplace_add(tmp,_s,_s,1.0,-1.0);
			}
			return *this;
		}

		/// for compatibility with tensor
		friend GenTensor<T> copy(const SliceGenTensor<T>& rhs) {
		    if (rhs._refGT.has_data()) return GenTensor<T>(rhs);
		    return GenTensor<T>();
		}

	};


    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
	/// cf tensor/tensor.h
    template <class T, class Q>
    GenTensor< TENSOR_RESULT_TYPE(T,Q) > transform(const GenTensor<Q>& t, const Tensor<T>& c) {
    	return t.transform(c);
    }

    /// helper struct for reduce
    struct RankReduceWrapper {
    	double eps;
    	RankReduceWrapper(const double& e) : eps(e) {}

    	template<typename T>
    	void operator()(GenTensor<T>& g) const {
    		g.reduce_rank(eps);
    	}
    };

    /// compare the rank of two GenTensors for sorting
    template<typename T>
    bool compare_rank(const GenTensor<T>& g1, const GenTensor<T>& g2) {
    	return g1.rank()<g2.rank();
    }


    /// add all the GenTensors of a given list

    /// If there are many tensors to add it's beneficial to do a sorted addition and start with
    /// those tensors with low ranks
    /// @param[in]	addends		a list with gentensors of same dimensions; will be destroyed upon return
    /// @param[in]	eps			the accuracy threshold
    /// @param[in]	are_optimal	flag if the GenTensors in the list are already in SVD format (if TT_2D)
    ///	@return		the sum GenTensor of the input GenTensors
    template<typename T>
	GenTensor<T> reduce(std::list<GenTensor<T> >& addends, double eps, bool are_optimal=false) {

    	// fast return
    	if (addends.size()==0) return GenTensor<T>();

    	typedef typename std::list<GenTensor<T> >::iterator iterT;

    	// make error relative
    	eps=eps/addends.size();

    	// if the addends are not in SVD format do that now so that we can call add_svd later
    	if (not are_optimal) {
    		RankReduceWrapper f(eps);
    		std::for_each(addends.begin(),addends.end(),f);
    	}

    	// remove zero ranks and sort the list according to the gentensor's ranks
		addends.remove_if(has_zero_rank<T>);
    	if (addends.size()==0) return GenTensor<T>();
		addends.sort(compare_rank<T>);

    	// do the additions
    	GenTensor<T> result=copy(addends.front());
    	for (iterT it=++addends.begin(); it!=addends.end(); ++it) {
    		GenTensor<T>& rhs=*it;
    		MADNESS_ASSERT(&rhs!=&result);
    		result.add_SVD(rhs,eps);
    	}
    	addends.clear();

    	return result;
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
    GenTensor<TENSOR_RESULT_TYPE(T,Q)> transform_dir(const GenTensor<Q>& t, const Tensor<T>& c,
                                          int axis) {

    	return t.transform_dir(c,axis);
    }

    template<typename T>
    static inline
    std::ostream& operator<<(std::ostream& s, const GenTensor<T>& g) {
        std::string str="GenTensor has no data";
        if (g.has_no_data()) {
            std::string str="GenTensor has no data";
            s << str.c_str() ;
        } else {
            str="GenTensor has data";
            s << str.c_str() << g.config();
        }
        return s;
    }

    namespace archive {
        /// Serialize a tensor
        template <class Archive, typename T>
		struct ArchiveStoreImpl< Archive, GenTensor<T> > {

			friend class GenTensor<T>;
			/// Stores the GenTensor to an archive
			static void store(const Archive& ar, const GenTensor<T>& t) {
				bool exist=t.has_data();
				ar & exist;
				if (exist) ar & t.config();
			};
		};


		/// Deserialize a tensor ... existing tensor is replaced
		template <class Archive, typename T>
		struct ArchiveLoadImpl< Archive, GenTensor<T> > {

			friend class GenTensor<T>;
			/// Replaces this GenTensor with one loaded from an archive
			static void load(const Archive& ar, GenTensor<T>& t) {
				// check for pointer existence
				bool exist=false;
				ar & exist;
				if (exist) {
					SRConf<T> conf;
					ar & conf;
					//t.config()=conf;
					t=GenTensor<T>(conf);
				}
			};
		};
    };

    /// outer product of two Tensors, yielding a low rank tensor
     template <class T, class Q>
     GenTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const GenTensor<T>& lhs2,
             const GenTensor<Q>& rhs2, const TensorArgs final_tensor_args) {
     	return outer(lhs2.full_tensor(),rhs2.full_tensor(), final_tensor_args);
     }

     /// outer product of two Tensors, yielding a low rank tensor
     template <class T, class Q>
     GenTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const Tensor<T>& lhs2,
             const Tensor<Q>& rhs2, const TensorArgs final_tensor_args) {

         MADNESS_ASSERT(final_tensor_args.tt==TT_2D);
     	typedef TENSOR_RESULT_TYPE(T,Q) resultT;

     	// srconf is shallow, do deep copy here
     	const Tensor<T> lhs=copy(lhs2);
     	const Tensor<Q> rhs=copy(rhs2);

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
 		GenTensor<resultT> coeff(srconf);
 		coeff.normalize();
 		return coeff;
     }

#endif /* not USE_LRT */

    #endif /* HAVE_GENTENSOR */

    /// change representation to targ.tt
    template<typename T>
    void change_tensor_type(GenTensor<T>& t, const TensorArgs& targs) {

        // fast return if possible
        const TensorType current_type=t.tensor_type();
        if (current_type==targs.tt) return;
        if (t.has_no_data()) return;

        // for now
        MADNESS_ASSERT(targs.tt==TT_FULL or targs.tt==TT_2D);
        MADNESS_ASSERT(current_type==TT_FULL or current_type==TT_2D);

        GenTensor<T> result;
        if (targs.tt==TT_FULL) {
            result=GenTensor<T>(t.full_tensor_copy(),targs);
        } else if (targs.tt==TT_2D) {
            MADNESS_ASSERT(current_type==TT_FULL);
            result=GenTensor<T>(t.full_tensor(),targs);
        }

        t=result;

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
    GenTensor<TENSOR_RESULT_TYPE(T,Q)> general_transform(const GenTensor<T>& t, const Tensor<Q> c[]) {
        return t.general_transform(c);
    }

    template <class T>
    GenTensor<T> general_transform(const GenTensor<T>& t, const Tensor<T> c[]) {
        return t.general_transform(c);
    }

    /// The class defines tensor op scalar ... here define scalar op tensor.
    template <typename T, typename Q>
    typename IsSupported < TensorTypeData<Q>, GenTensor<T> >::type
    operator*(const Q& x, const GenTensor<T>& t) {
        return t*x;
    }


}   // namespace madness

#if HAVE_GENTENSOR
#ifdef USE_LRT
#include <madness/tensor/lowranktensor.h>

namespace madness {

template<typename T>
class GenTensor : public LowRankTensor<T> {

public:

    using LowRankTensor<T>::LowRankTensor;

    GenTensor<T>() : LowRankTensor<T>() {}
    GenTensor<T>(const GenTensor<T>& g) : LowRankTensor<T>(g) {}
    GenTensor<T>(const LowRankTensor<T>& g) : LowRankTensor<T>(g) {}
    GenTensor<T>(const SRConf<T>& sr1) : LowRankTensor<T>(SVDTensor<T>(sr1)) {
    }

    operator LowRankTensor<T>() const {return *this;}
    operator LowRankTensor<T>() {return *this;}

    /// general slicing, shallow; for temporary use only!
    SliceGenTensor<T> operator()(const std::vector<Slice>& s) {
        return SliceGenTensor<T>(*this,s);
    }

    /// general slicing, shallow; for temporary use only!
    const SliceGenTensor<T> operator()(const std::vector<Slice>& s) const {
        return SliceGenTensor<T>(*this,s);
    }

    /// assign a number to this tensor
    GenTensor<T>& operator=(const T& number) {
        LowRankTensor<T>& base=*this;
        base=(number);
        return *this;
    }


    std::string what_am_i() const {return TensorArgs::what_am_i(this->tensor_type());};

    SRConf<T>& config() const {
        MADNESS_ASSERT(this->type==TT_2D and (this->impl.svd));
        return *this->impl.svd.get();
    }
    GenTensor<T> get_configs(const int& start, const int& end) const {
        MADNESS_ASSERT(this->type==TT_2D and (this->impl.svd));
        return GenTensor<T>(config().get_configs(start,end));
    }

    /// deep copy of rhs by deep copying rhs.configs
    friend GenTensor<T> copy(const GenTensor<T>& rhs) {
        return GenTensor<T>(copy(LowRankTensor<T>(rhs)));
    }

};

/// implements a slice of a GenTensor
template <typename T>
class SliceGenTensor : public SliceLowRankTensor<T> {
public:
    using SliceLowRankTensor<T>::SliceLowRankTensor;

    SliceGenTensor<T>(const SliceGenTensor<T>& g) : SliceLowRankTensor<T>(g) {}
    SliceGenTensor<T>(const SliceLowRankTensor<T>& g) : SliceLowRankTensor<T>(g) {}

    operator SliceLowRankTensor<T>() const {return *this;}
    operator SliceLowRankTensor<T>() {return *this;}

    /// inplace zero-ing as in g(s)=0.0
    SliceGenTensor<T>& operator=(const T& number) {
        SliceLowRankTensor<T>& base=*this;
        base=number;
        return *this;
    }

};

/// add all the GenTensors of a given list

 /// If there are many tensors to add it's beneficial to do a sorted addition and start with
 /// those tensors with low ranks
 /// @param[in]  addends     a list with gentensors of same dimensions; will be destroyed upon return
/// @param[in]  eps         the accuracy threshold
/// @param[in]  are_optimal flag if the GenTensors in the list are already in SVD format (if TT_2D)
/// @return     the sum GenTensor of the input GenTensors
template<typename T>
GenTensor<T> reduce(std::list<GenTensor<T> >& addends, double eps, bool are_optimal=false) {
    typedef typename std::list<GenTensor<T> >::iterator iterT;
    GenTensor<T> result=copy(addends.front());
    for (iterT it=++addends.begin(); it!=addends.end(); ++it) {
        result+=*it;
    }
    result.reduce_rank(eps);
    return result;

}

/// outer product of two Tensors, yielding a low rank tensor
 template <class T, class Q>
 GenTensor<TENSOR_RESULT_TYPE(T,Q)> outer(const GenTensor<T>& lhs2,
         const GenTensor<Q>& rhs2, const TensorArgs final_tensor_args) {
//    return outer_low_rank(lhs2.full_tensor(),rhs2.full_tensor(), final_tensor_args);
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

     LowRankTensor<T> lhs=lhs2.convert(targs);
     LowRankTensor<Q> rhs=rhs2.convert(targs);

     LowRankTensor<resultT> result=outer(lhs,rhs,final_tensor_args);
     return result;


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

     LowRankTensor<T> lhs(lhs2,targs);
     LowRankTensor<Q> rhs(rhs2,targs);
     LowRankTensor<resultT> result=outer(lhs,rhs,final_tensor_args);
     return result;
 }



namespace archive {
/// Serialize a tensor
template <class Archive, typename T>
struct ArchiveStoreImpl< Archive, GenTensor<T> > {

    friend class GenTensor<T>;
    /// Stores the GenTensor to an archive
    static void store(const Archive& ar, const GenTensor<T>& t) {
        LowRankTensor<T> tt(t);
        ar & tt;
    };
};


/// Deserialize a tensor ... existing tensor is replaced
template <class Archive, typename T>
struct ArchiveLoadImpl< Archive, GenTensor<T> > {

    friend class GenTensor<T>;
    /// Replaces this GenTensor with one loaded from an archive
    static void load(const Archive& ar, GenTensor<T>& t) {
        LowRankTensor<T> tt;
        ar & tt;
        t=tt;
    };
};
};

}
#endif /* USE_LRT */
#endif /* HAVE_GENTENSOR */
#endif /* GENTENSOR_H_ */
