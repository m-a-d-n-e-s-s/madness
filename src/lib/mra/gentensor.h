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


#include "tensor/tensor.h"
#include "mra/srconf.h"
#include <stdexcept>

namespace madness {

	template <class T> class SliceGenTensor;


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
			if (tt==TT_3D) return "TT_3D";
			if (tt==TT_FULL) return "TT_FULL";
			return "unknown tensor type";
		}
		template <typename Archive>
		void serialize(const Archive& ar) {
		    int i=int(tt);
		    ar & thresh & i;
		    tt=TensorType(i);
		}
	};


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


		GenTensor<T> reconstruct_tensor() const {return *this;}
		GenTensor<T> full_tensor() const {return *this;}
		GenTensor<T>& full_tensor() {return *this;}
        GenTensor<T> full_tensor_copy() const {return *this;}
        GenTensor<T> full_tensor_copy() {return *this;}

        bool has_data() const {return this->size()>0;};
		long rank() const {return -1;}

        void reduceRank(const double& eps) {return;};
        void right_orthonormalize(const double& eps) {return;};
        void reconstruct_and_decompose(const double& eps) {return;}

        std::string what_am_i() const {return "GenTensor, aliased to Tensor";};
		TensorType tensor_type() const {return TT_FULL;}

		void accumulate_into(Tensor<T>& t, const double& eps, const double& fac) const {t+=*this *fac;}
		void accumulate_into(Tensor<T>& t, const std::complex<double>& fac) const {MADNESS_EXCEPTION("",0);}
		void accumulate_into(GenTensor<T>& t, const double& eps, const double& fac) const {t+=*this*fac;}
		void accumulate_into(GenTensor<T>& t, const std::complex<double>& fac) const {t+=*this*fac;}
		void add_SVD(const GenTensor<T>& rhs, const double& eps) {*this+=rhs;}

		/// return the additional safety for rank reduction
		static double fac_reduce() {return -1.0;};

	};

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
			long sz, id;
			s & sz & id;
			if (id != t.id()) throw "type mismatch deserializing a tensor";
			if (sz) {
				long _ndim, _dim[TENSOR_MAXDIM];
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

		typedef SRConf<T> configT;
		typedef Tensor<T> tensorT;
		typedef GenTensor<T> gentensorT;
		typedef std::shared_ptr<configT> sr_ptr;

		/// pointer to the low rank tensor
		sr_ptr _ptr;

		/// the machine precision
		static const double machinePrecision=1.e-14;

		/// sqrt of the machine precision
		static const double sqrtMachinePrecision=1.e-7;

		/// safety for rank reduction
		static const double facReduce=1.e-3;

	public:

		/// empty ctor
		GenTensor() : _ptr() {
		}

		/// copy ctor, shallow
//		GenTensor(const GenTensor<T>& rhs) : _ptr(rhs._ptr) { // DON'T DO THIS: USE_COUNT BLOWS UP
		GenTensor(const GenTensor<T>& rhs) {
			_ptr=rhs._ptr;
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
		GenTensor(const Tensor<T>& rhs, const TensorArgs& targs)
			: _ptr(new configT(rhs.ndim(),rhs.dim(0),targs.tt)) {

			MADNESS_ASSERT(rhs.ndim()>0);
			MADNESS_ASSERT(rhs.iscontiguous());
    		for (long idim=0; idim<rhs.ndim(); idim++) {
    			MADNESS_ASSERT(rhs.dim(0)==rhs.dim(idim));
    		}

			// adapt form of values
			std::vector<long> d(_ptr->dim_eff(),_ptr->kVec());
			Tensor<T> values_eff=rhs.reshape(d);

			// direct reduction on the polynomial values on the Tensor
			TensorType tt=tensor_type();
			if (tt==TT_3D) {
				this->doReduceRank(targs.thresh,values_eff);
			} else if (tt==TT_2D) {
				this->computeSVD(targs.thresh,values_eff);
			} else if (tt==TT_FULL){
				_ptr.reset(new configT(copy(rhs)));
			} else {
				MADNESS_EXCEPTION("unknown tensor type in GenTensor(tensorT,targs)",0);
			}
		}

		/// ctor with a regular Tensor and arguments, deep
		GenTensor(const Tensor<T>& rhs, const double& thresh, const TensorType& tt) {
			gentensorT tmp(rhs,TensorArgs(thresh,tt));
			*this=tmp;
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

//        /// Replaces this GenTensor with one loaded from an archive
//        template <typename Archive>
//        void load(World& world, Archive& ar) {
//            bool exist;
//            ar & exist;
//            if (exist) {
//               _ptr.reset(new configT());
//               _ptr->load(ar);
//            } else {
//               _ptr.reset();
//            }
//        }
//
//
//        /// Stores the GenTensor to an archive
//        template <typename Archive>
//        void store(Archive& ar) const {
//            bool exist=(_ptr);
//            ar & exist;
//            if (exist) _ptr->store(ar);
//        }


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
			MADNESS_EXCEPTION("no assignment with a number for GenTensor",0);
		}


		/// ctor w/ configs, shallow (indirectly, via vector_)
		GenTensor(const SRConf<T>& config) : _ptr(new configT(config)) {
		}

	private:
		/// return a slice of this (deep copy)
		gentensorT copy_slice(const std::vector<Slice>& s) const {

			// consistency check
			MADNESS_ASSERT(s.size()==this->dim());
			MADNESS_ASSERT(s[0].step==1);

			// fast return if possible
			if (this->has_no_data()) {
				int k_new=s[0].end-s[0].start+1;
				return gentensorT (this->tensor_type(),k_new,this->dim());
			}

			// fast return for full rank tensors
			if (tensor_type()==TT_FULL) {
				tensorT a=copy(full_tensor())(s);
				return gentensorT(configT(a));
			}

			_ptr->make_structure();

			// get dimensions
			const TensorType tt=this->tensor_type();
			const int merged_dim=this->_ptr->dim_per_vector();
			const int dim_eff=this->_ptr->dim_eff();
			const int rank=this->rank();
			int k_new=s[0].end-s[0].start+1;
			if (s[0].end<0) k_new+=this->get_k();

			// get and reshape the vectors, slice and re-reshape again;
			// this is shallow
			const gentensorT& sr=*this;

			std::vector<Tensor<T> > vectors(dim_eff,Tensor<T>());

			for (int idim=0; idim<dim_eff; idim++) {

				// assignment from/to slice is deep-copy
				if (merged_dim==1) {
					if (rank>0) {
						vectors[idim]=copy(sr._ptr->ref_vector(idim)(Slice(0,rank-1),s[idim]));
					} else {
						vectors[idim]=Tensor<T>(0,s[idim].end-s[idim].start+1);
					}
				} else if (merged_dim==2) {
					if (rank>0) {
						vectors[idim]=copy(sr._ptr->ref_vector(idim)(Slice(0,rank-1),s[2*idim],s[2*idim+1]));
					} else {
						vectors[idim]=tensorT(0,s[2*idim].end-s[2*idim].start+1,
												s[2*idim+1].end-s[2*idim+1].start+1);
					}
				} else if (merged_dim==3) {
					if (rank>0) {
						vectors[idim]=copy(sr._ptr->ref_vector(idim)(Slice(0,rank-1),
								s[3*idim],s[3*idim+1],s[3*idim+2]));
					} else {
						vectors[idim]=tensorT(0,s[3*idim].end-s[3*idim].start+1,
								s[3*idim+1].end-s[3*idim+1].start+1,
								s[3*idim+2].end-s[3*idim+2].start+1);

					}
				} else MADNESS_EXCEPTION("unknown number of dimensions in GenTensor::copy_slice()",0);
			}

			// work-around for rank==0
			Tensor<double> weights;
			if (rank>0) {
				weights=copy(this->_ptr->weights_(Slice(0,rank-1)));
			} else {
				weights=Tensor<double>(int(0));
			}
			const configT conf(weights,vectors,this->dim(),k_new,tt);

			return gentensorT(conf);

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
//			rhs.accumulate_into(*this,1.0);
			if (tensor_type()==TT_FULL) {
				if (not rhs.tensor_type()==TT_FULL)
					MADNESS_EXCEPTION("gentensorT::operator+= failure; call accumulate_into",0);
				this->full_tensor()+=rhs.full_tensor();
				return *this;
			}
			rhs.append(*this,1.0);
			return *this;
		}

		/// inplace addition
		gentensorT& operator+=(const SliceGenTensor<T>& rhs) {

			MADNESS_ASSERT(this->_ptr->type()==rhs._refGT._ptr->type());
			const std::vector<Slice> s(this->ndim(),Slice(0,get_k()-1,1));
			this->_ptr->inplace_add(*rhs._refGT._ptr,s,rhs._s,1.0,1.0);
			return *this;

		}

		/// multiply with a number
		gentensorT operator*(const double& x) const {
			gentensorT result(copy(*this));
        	result.scale(x);
        	return result;
		}

		/// multiply with a number
		gentensorT operator*(const std::complex<double>& x) const {
			MADNESS_EXCEPTION("only double factors in GenTensor::operator*",0);
		}

	    /// Inplace generalized saxpy ... this = this*alpha + other*beta
	    gentensorT& gaxpy(std::complex<double> alpha, const gentensorT& rhs, std::complex<double> beta) {
	    	MADNESS_EXCEPTION("GenTensor::gaxpy only with double factors",0);
	    }

	    /// Inplace generalized saxpy ... this = this*alpha + other*beta
	    gentensorT& gaxpy(const double alpha, const gentensorT& rhs, const double beta) {
	    	if (not alpha==1.0) this->scale(alpha);
//	    	rhs.accumulate_into(*this,beta);
	    	rhs.append(*this,beta);
	    	return *this;
	    }

		/// multiply with a scalar
		gentensorT& scale(const double& dfac) {
			if (!_ptr) return *this;
			if (tensor_type()==TT_FULL) {
				full_tensor().scale(dfac);
			} else {
				_ptr->scale(dfac);
			}
			return *this;
		};

		/// multiply with a scalar
		gentensorT& scale(const std::complex<double>& dfac) {
			MADNESS_EXCEPTION("only double factors in GenTensor::scale",0);
			return *this;
		};

		/// orthonormalize the right subspace and shift the weights to the left one
		void right_orthonormalize(const double& thresh) {
			if (has_no_data()) return;
			if (tensor_type()==TT_2D) {
				_ptr->undo_structure();
				if (rank()==1) _ptr->normalize_and_shift_weights_to_x();
				else _ptr->right_orthonormalize(thresh*facReduce);
			}
		}

		/// normalize
		void normalize() {
			TensorType tt=tensor_type();
			if ((tt==TT_2D) or (tt==TT_3D)) _ptr->normalize();
		};

		void fillrandom(const int r=1) {
			if (tensor_type()==TT_FULL) full_tensor().fillrandom();
			else _ptr->fillWithRandom(r);
		}

		/// do we have data? note difference to SRConf::has_data() !
		bool has_data() const {return (_ptr);}

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

		/// return the number length of the underlying vectors
		unsigned int kVec() const {return _ptr->kVec();};

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

		/// returns the Frobenius norm
		double normf() const {
			if (has_no_data()) return 0.0;
			_ptr->undo_structure();
			return config().normf();
		};

		/// returns the trace of <this|rhs>
		T trace_conj(const GenTensor<T>& rhs) const {
			return overlap(*_ptr,*rhs._ptr);
		}

		/// overlap between two SepReps (Frobenius inner product)
		friend T overlap(const gentensorT& rhs, const gentensorT& lhs)  {

			// fast return if possible
			if ((lhs.rank()==0) or (rhs.rank()==0)) return 0.0;

			MADNESS_ASSERT(compatible(lhs,rhs));
			MADNESS_ASSERT(lhs.tensor_type()==rhs.tensor_type());

			const T ovlp=overlap(*lhs._ptr,*rhs._ptr);
			return ovlp;
		}

        /// Inplace multiply by corresponding elements of argument Tensor
        GenTensor<T>& emul(const GenTensor<T>& t) {
        	print("no GenTensor<T>::emul yet");
        	MADNESS_ASSERT(0);
        }

		/// return a Tensor, no matter what
		Tensor<T> full_tensor_copy() const {
			const TensorType tt=tensor_type();
			if (tt==TT_NONE) return Tensor<T>();
			else if (tt==TT_3D or tt==TT_2D) return this->reconstruct_tensor();
			else if (tt==TT_FULL) {
				return copy(full_tensor());
			} else {
				print(TensorArgs::what_am_i(tt),"unknown tensor type");
				MADNESS_ASSERT(0);
			}
		}

		// reduce the rank of this
		void reduceRank(const double& eps) {

			if (rank()==0) return;

			// direct reduction on the polynomial values on the Tensor
			if (tensor_type()==TT_FULL or tensor_type()==TT_NONE) {
				return;
			} else if (this->tensor_type()==TT_3D) {
				this->doReduceRank(eps,Tensor<T>());
				this->_ptr->make_structure();
			} else if (this->tensor_type()==TT_2D) {

				// determine what is faster: reconstruction or direct rank reduction
				const double ratio=(_ptr->kVec())/(2.0*rank());
				if (ratio>1.0) {

					config().undo_structure();
					if (OrthoMethod::om==ortho3_ or OrthoMethod::om==ortho6_) {
						config().divide_and_conquer_reduce(eps*facReduce);

					} else if (OrthoMethod::om==sequential_) {
						config().sequential_orthogonalization(eps*facReduce);

					} else if (OrthoMethod::om==reconstruct_) {
						reconstruct_and_decompose(eps);

					} else {
						MADNESS_EXCEPTION("confused about orthogonalization method??",0);
					}
					config().make_structure();
				} else {
					reconstruct_and_decompose(eps);
				}
			} else {
				MADNESS_EXCEPTION("unknown tensor type in GenTensor::reduceRank()",0);
			}
			MADNESS_ASSERT(this->_ptr->has_structure() or this->rank()==0);
		}

		/// reduce rank by reconstruction of the full tensor and subsequent SVD decomposition
		void reconstruct_and_decompose(const double& eps) {
			Tensor<T> values=this->reconstruct_tensor();
			std::vector<long> d(_ptr->dim_eff(),_ptr->kVec());
			Tensor<T> values_eff=values.reshape(d);
			this->computeSVD(eps,values_eff);
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

			/*
			 * reconstruct the tensor first to the configurational dimension,
			 * then to the real dimension
			 */

			// fast return for full rank tensors
			if (tensor_type()==TT_FULL) return full_tensor();

			// for convenience
			const unsigned int conf_dim=this->_ptr->dim_eff();
			const unsigned int conf_k=this->kVec();			// possibly k,k*k,..
			const unsigned int rank=this->rank();
			long d[TENSOR_MAXDIM];

			// fast return if possible
			if (rank==0) {
				// reshape the tensor to the really required one
				const long k=this->get_k();
				const long dim=this->dim();
				for (long i=0; i<dim; i++) d[i] = k;

				return Tensor<T> (dim,d,true);
			}


			// set up result Tensor (in configurational dimensions)
			for (long i=0; i<conf_dim; i++) d[i] = conf_k;
			tensorT s(conf_dim,d,true);

			// flatten this
			gentensorT sr=*this;
			sr._ptr->undo_structure();
			//			sr._ptr->semi_flatten();

			// and a scratch Tensor
			Tensor<T>  scr(rank);
			Tensor<T>  scr1(rank);
			Tensor<T>  scr2(rank);

			if (conf_dim==1) {

				for (unsigned int i0=0; i0<conf_k; i0++) {
					scr=sr._ptr->weights_;
					//					scr.emul(F[0][i0]);
					T buffer=scr.sum();
					s(i0)=buffer;
				}

			} else if (conf_dim==2) {

				//				tensorT weight_matrix(rank,rank);
				//				for (unsigned int r=0; r<rank; r++) {
				//					weight_matrix(r,r)=this->weight(r);
				//				}
				//				s=inner(weight_matrix,sr._ptr->refVector(0));
				//				s=inner(s,sr._ptr->refVector(1),0,0);
				tensorT sscr=copy(sr._ptr->ref_vector(0)(sr._ptr->c0()));
				for (unsigned int r=0; r<rank; r++) {
					const double w=_ptr->weights(r);
					for (unsigned int k=0; k<conf_k; k++) {
						sscr(r,k)*=w;
					}
				}
				inner_result(sscr,sr._ptr->ref_vector(1)(sr._ptr->c0()),0,0,s);


			} else if (conf_dim==3) {

				for (unsigned int i0=0; i0<conf_k; i0++) {
					scr=copy(sr._ptr->weights_(Slice(0,sr.rank()-1)));
					scr.emul(sr._ptr->ref_vector(0)(Slice(0,rank-1),i0));
					for (unsigned int i1=0; i1<conf_k; i1++) {
						scr1=copy(scr);
						scr1.emul(sr._ptr->ref_vector(1)(Slice(0,rank-1),i1));
						for (unsigned int i2=0; i2<conf_k; i2++) {
							scr2=copy(scr1);
							scr2.emul(sr._ptr->ref_vector(2)(Slice(0,rank-1),i2));
							s(i0,i1,i2)=scr2.sum();
						}
					}
				}
			} else {
				print("only config_dim=1,2,3 in GenTensor::reconstructTensor");
				MADNESS_ASSERT(0);
			}


			// reshape the tensor to the really required one
			const long k=this->get_k();
			const long dim=this->dim();
			for (long i=0; i<dim; i++) d[i] = k;

			Tensor<T> s2=s.reshape(dim,d);
			return s2;
		}

	    /// accumulate this into t
	    void accumulate_into(GenTensor<T>& t, const double& thresh, const std::complex<double>& fac) const {
	    	MADNESS_EXCEPTION("no GenTensor::accumulate_into with complex fac",0);
	    }

	    /// accumulate this into t
	    void accumulate_into(GenTensor<T>& t, const double& thresh, const double& fac) const {

	    	if (has_no_data()) return;

	    	// need case-by-case decision
	    	if (tensor_type()==TT_FULL) {
	    		if (t.has_no_data()) {
	    			t=copy(*this);
	    			t.full_tensor()*=fac;
	    		} else {
	    			t.full_tensor()+=this->full_tensor()*fac;
	    		}
	    	} else if (tensor_type()==TT_2D) {
	    		if (t.tensor_type()==TT_FULL) {
	    			accumulate_into(t.full_tensor(),fac);
	    		} else {

		    		t._ptr->undo_structure();
		    		_ptr->undo_structure();
		    		t.config().low_rank_add_sequential(*_ptr,thresh*facReduce,fac);
	    		}
	    	} else if (tensor_type()==TT_3D) {
	    		t._ptr->undo_structure();
	    		_ptr->undo_structure();
	    		t._ptr->append(*this->_ptr,fac);
	    	} else {
	    		MADNESS_EXCEPTION("unknown tensor type in GenTensor::accumulate_into",0);
	    	}
	    }

	    /// accumulate this into t
	    void accumulate_into(Tensor<T>& t, const std::complex<double>& fac) const {
	    	MADNESS_EXCEPTION("no GenTensor::accumulate_into with complex fac",0);
	    }

		/// reconstruct this to full rank, and accumulate into t
		void accumulate_into(Tensor<T>& t, const double fac) const {

			// fast return if possible
			if (this->has_no_data()) return;	// no SRConf at all
			if (_ptr->has_no_data()) return;	// rank==0

			// fast return for full rank tensors
			if (tensor_type()==TT_FULL) {
				t+=this->full_tensor()*fac;
				return;
			}

			MADNESS_ASSERT(t.iscontiguous());

			// for convenience
			const unsigned int conf_dim=this->_ptr->dim_eff();
			const unsigned int conf_k=this->kVec();			// possibly k,k*k,..
			const unsigned int rank=this->rank();
			long d[TENSOR_MAXDIM];
			MADNESS_ASSERT(conf_dim==2 or conf_dim==3);

			// set up result Tensor (in configurational dimensions)
			for (long i=0; i<conf_dim; i++) d[i] = conf_k;
			tensorT s=t.reshape(conf_dim,d);

			// flatten this
			_ptr->undo_structure();

			if (conf_dim==2) {

				tensorT sscr=copy(_ptr->ref_vector(0)(_ptr->c0()));
				if (fac!=1.0) sscr.scale(fac);
				for (unsigned int r=0; r<rank; r++) {
					const double w=_ptr->weights(r);
					for (unsigned int k=0; k<conf_k; k++) {
						sscr(r,k)*=w;
					}
				}
				inner_result(sscr,_ptr->ref_vector(1)(_ptr->c0()),0,0,s);

			} else {

				// include weights in first vector
				tensorT scr1=copy(_ptr->ref_vector(0)(_ptr->c0()));
				if (fac!=1.0) scr1.scale(fac);
				for (unsigned int r=0; r<rank; r++) {
					scr1(r,Slice(_))*=_ptr->weights(r);
				}

				// merge second and third vector to G(r,k1,k2)
				tensorT scr2(rank,conf_k,conf_k);
				for (unsigned int r=0; r<rank; r++) {
					scr2(r,Slice(_),Slice(_))=outer(_ptr->ref_vector(1)(r,Slice(_)),
													_ptr->ref_vector(2)(r,Slice(_)));
				}

				inner_result(scr1,scr2,0,0,s);

			}

		}

		/// append this to rhs, shape must conform
		void append( gentensorT& rhs, const double fac=1.0) const {
			rhs._ptr->undo_structure();
			_ptr->undo_structure();
			rhs.config().append(*this->_ptr,fac);
		}

		/// add SVD
		void add_SVD(const gentensorT& rhs, const double& thresh) {
			if (rhs.has_no_data()) return;
			if (tensor_type()==TT_FULL or tensor_type()==TT_NONE) {
				this->full_tensor()+=rhs.full_tensor();
				return;
			}
			if (has_no_data()) {
				*this=rhs;
				return;
			}
			config().undo_structure();
			rhs._ptr->undo_structure();
			config().add_SVD(rhs.config(),thresh);
		}

	    /// check compatibility
		friend bool compatible(const gentensorT& rhs, const gentensorT& lhs) {
			return ((rhs.tensor_type()==lhs.tensor_type()) and (rhs.get_k()==lhs.get_k())
					and (rhs.dim()==lhs.dim()));
		};

		/// compute a best one-term approximation wrt to this
		gentensorT oneTermApprox(const double& eps, std::vector<tensorT>& B1) const {

			/*
			 * return a SepRep that represents this with only one term
			 */

			// new random SepRep of rank 1 and optimize wrt this
			configT residual(dim(),get_k(),tensor_type());
			configT dummy(residual);

			residual.fillWithRandom();

			Tensor<T> t1;
			Tensor<T> B2;
			//			std::vector<Tensor<T> > B2(B1);

			// optimize the new SepRep wrt the residual; *this is the residual
			// maxloop can be large because we have an additional stopping criterion
			const unsigned int maxloop=50;
			bool successful=residual.optimize(*this,1.0,dummy,0.0,t1,0.0,eps,maxloop,B1,B2);
			if (not successful) {
				std::cout << "NaNs in oneTermApprox " << std::endl;
				MADNESS_ASSERT(0);
			}
			return residual;
		}

		/// transform the Legendre coefficients with the tensor
		gentensorT transform(const Tensor<T> c) const {
			_ptr->make_structure();
			return gentensorT (this->_ptr->transform(c));
		}

		/// transform the Legendre coefficients with the tensor
		template<typename Q>
		gentensorT general_transform(const Tensor<Q> c[]) const {
			return gentensorT (this->config().general_transform(c));
		}

		/// inner product
		gentensorT transform_dir(const Tensor<T>& c, const int& axis) const {
			return this->_ptr->transform_dir(c,axis);
		}

		/// return a reference to the SRConf
		const SRConf<T>& config() const {return *_ptr;}

		/// return a reference to the SRConf
		SRConf<T>& config() {return *_ptr;}

		/// return the additional safety for rank reduction
		static double fac_reduce() {return facReduce;};

	private:

		/// optimize this wrt reference, and return the error norm
		bool optimize(const gentensorT& ref1, const double& fac1,
				const gentensorT& ref2, const double fac2,
				const tensorT& ref3, const double fac3,
				const double& eps, const unsigned int& maxloop,
				std::vector<tensorT>& B1, std::vector<tensorT>& B2) {


			// for convenience
			const unsigned int config_dim=this->_ptr->dim_eff();
			const unsigned int rF=this->rank();
			const unsigned int rG1=ref1.rank();
			const unsigned int rG2=ref2.rank();

			// some checks
			if (fac1!=0) MADNESS_ASSERT(compatible(*this,ref1));
			if (fac2!=0) MADNESS_ASSERT(compatible(ref1,ref2));

			double oldnorm=1.0;

			// reshape the scratch Tensor
			for (unsigned int idim=0; idim<config_dim; idim++) {
				//				B1[idim].reshape(rF,rG1);
				//				B2[idim].reshape(rF,rG2);
				B1[idim]=Tensor<T>(rF,rG1);
				B2[idim]=Tensor<T>(rF,rG2);				// 0.2 sec
			}

			// keep optimizing until either the norm doesn't change
			// or we hit some max number of runs
			unsigned int iloop=0;
			for (iloop=0; iloop<maxloop; iloop++) {

				// optimize once for all dimensions
				try {
					this->generalizedALS(ref1,fac1,ref2,fac2,ref3,fac3,B1,B2);
				} catch (std::runtime_error) {
					throw std::runtime_error("rank reduction failed");
					return false;
				}

				/*
				 * for residuals: also exit if norm vanishes, or if
				 * norm doesn't change any more
				 */
				//		const double norm=fabs(this->_ptr->weights(this->rank()-1));
				const double norm=this->normf();
				if (iloop>1) {
					const double ratio=oldnorm/norm;
					//			std::cout << "  ratio " << ratio << " norm " << norm << std::endl;
					if (fabs(ratio-1.0)<0.003) break;
				}
				oldnorm=norm;

			}
			return true;
		}

		/// perform the alternating least squares algorithm directly on function values,
		/// minus the difference SR
		void generalizedALS(const gentensorT& ref1, const double& fac1,
				const gentensorT& ref2, const double& fac2,
				const Tensor<T> & ref3, const double& fac3,
				std::vector<Tensor<T> >& B1, std::vector<Tensor<T> >& B2) {

			// for convenience
			const unsigned int dim=this->_ptr->dim_eff();
			const unsigned int kvec=this->kVec();
			gentensorT& trial=*this;

			const bool have1=(fac1!=0.0 and ref1.rank()>0);
			const bool have2=(fac2!=0.0 and ref2.rank()>0);
			const bool have3=(fac3!=0.0);

			/*
			 * rF is the rank of the trial SepRep
			 * rG1 is the rank of the SepRep ref1
			 * rG2 is the rank of the SepRep ref2
			 */
			const unsigned int rF=trial.rank();
			const unsigned int rG1=ref1.rank();
			const unsigned int rG2=ref2.rank();

			// some checks
			MADNESS_ASSERT(ref1.tensor_type()==this->tensor_type());
			if (have2) MADNESS_ASSERT(ref2.tensor_type()==this->tensor_type());
			if (have2) MADNESS_ASSERT(compatible(ref1,ref2));
			MADNESS_ASSERT(dim==B1.size());
			MADNESS_ASSERT(dim==B2.size());
			MADNESS_ASSERT(rG1>=0);
			MADNESS_ASSERT(rG2>=0);
			MADNESS_ASSERT(rF>0);
			MADNESS_ASSERT(B1[0].dim(0)==rF);
			MADNESS_ASSERT(B1[0].dim(1)==rG1);
			MADNESS_ASSERT(B2[0].dim(0)==rF);
			MADNESS_ASSERT(B2[0].dim(1)==rG2);

			// for controlling the condition number, sec. (3.2) of BM2005
			const double alpha=machinePrecision;
			Tensor<T> unity(trial.rank(),trial.rank());
			for (unsigned int i=0; i<trial.rank(); i++) {
				unity(i,i)=alpha;
			}

			// some scratch Tensors
			Tensor<T>  fvec(kvec);
			Tensor<T>  vecb(rF,kvec);
			//			Tensor<T>  vecb(kvec,rF);
			// no copy ctor here for B, since it is shallow!
			//			std::vector<Tensor<T> > B(dim,Tensor<T> (rF,rF));
			std::vector<Tensor<T> > B(dim);
			for (unsigned int idim=0; idim<dim; idim++) B[idim]=Tensor<T> (rF,rF);


			/*
			 * first make all factors of the two B matrices of eq. (3.3) BM2005,
			 * 	- B is the overlap <F | F> for each dimension
			 * 	- B_GF is the overlap <F | G> for each dimension
			 * 	- as B[idim] and B_GF[idim] is not required, use it as scratch for constructing
			 * 		the product \Prod_{dim\idim} <F | G>
			 * 	- include the weights in B_GF[idim], and construct the vectors b as
			 * 		b(l',k) = B_GF(l',l) * F(l,k)
			 */

			// leave out idim=0, since it not required in the first alteration
			for (unsigned int idim=1; idim<dim; idim++) {
				if (have1) makeB(B1[idim],idim,*trial._ptr,*ref1._ptr);
				if (have2) makeB(B2[idim],idim,*trial._ptr,*ref2._ptr);
				makeB(B[idim],idim,*trial._ptr,*trial._ptr);
				B[idim]+=unity;
				//				print("B[idim]");
				//				print(B[idim]);
			}

			// next loop over all dimensions
			for (unsigned int idim=0; idim<dim; idim++) {

				// reconstruct B and B_GF for the dimension that has been
				// altered before, include the unit matrix
				if (idim>0) {
					if (have1) makeB(B1[idim-1],idim-1,*trial._ptr,*ref1._ptr);
					if (have2) makeB(B2[idim-1],idim-1,*trial._ptr,*ref2._ptr);
					makeB(B[idim-1],idim-1,*trial._ptr,*trial._ptr);
					B[idim-1]+=unity;

				}

				// construct the products of the B's and B_GF's
				B1[idim]=1.0;
				B2[idim]=1.0;
				B[idim]=1.0;
				for (unsigned int jdim=0; jdim<dim; jdim++) {
					if (jdim!=idim) {
						if (have1) B1[idim].emul(B1[jdim]);
						if (have2) B2[idim].emul(B2[jdim]);
						B[idim].emul(B[jdim]);
					}
				}


				/*
				 * now construct the b vector of eq (3.4) BM2005
				 */
				vecb.fill(0.0);

				// bring the quantity \prod_(i/=k) < G|F > in some efficient form
				// it is independent of jk, and include the weights
				if (have1) {
					for (unsigned int l=0; l<rG1; l++) {
						const double w=ref1._ptr->weights(l);
						for (unsigned int l1=0; l1<rF; l1++) {
							B1[idim](l1,l)*=w*fac1;
						}
					}
					Tensor<T> tmp=ref1._ptr->ref_vector(idim)(ref1._ptr->c0());
					vecb+=madness::inner(B1[idim],tmp);
				}
				if (have2) {
					for (unsigned int l=0; l<rG2; l++) {
						const double w=ref2._ptr->weights(l);
						for (unsigned int l1=0; l1<rF; l1++) {
							B2[idim](l1,l)*=w*fac2;
						}
					}
					Tensor<T>  tmp=ref2._ptr->ref_vector(idim)(ref2._ptr->c0());
					vecb+=madness::inner(B2[idim],tmp);
				}

				/*
				 * now construct the b(r,k) vector for the Tensor values
				 */
				if (have3) {

					MADNESS_ASSERT((dim==3) or (dim==2));
					if (dim==3) {

						// b[jk][rF] += inner( s[jk][k1,k2] , [rF][\prod[k1,k2])

						// reorder s[k1,k2,k3] to s[jk][k1,k2]
						// need to cycle backwards..
						const Tensor<T> weights=copy(ref3.cycledim(dim-idim,0,dim-1));
						const Tensor<T> w=weights.reshape(kvec,kvec*kvec);

						const unsigned int idim0=(idim+1)%dim;
						const unsigned int idim1=(idim+2)%dim;

						// set up \Prod_{i\neq k} <G | F>
						Tensor<T> prod(rF,kvec,kvec);
						for (unsigned int r=0; r<rF; r++) {
							for (unsigned int i0=0; i0<kvec; i0++) {
								const T F1=trial._ptr->vector_[idim0](r,i0);
								for (unsigned int i1=0; i1<kvec; i1++) {
									const T F2=trial._ptr->vector_[idim1](r,i1);
									prod(r,i0,i1)=F1*F2;
								}
							}
						}
						const Tensor<T> p=prod.reshape(rF,kvec*kvec);

						// compute the contrib to b
						vecb+=madness::inner(p,w,-1,-1);

					} else if (dim==2) {

						// b[jk][rF] += inner( s[jk][k1,k2] , [rF][\prod[k1,k2])

						// reorder s[k1,k2,k3] to s[jk][k1,k2]
						// need to cycle backwards..
						const Tensor<T> weights=copy(ref3.cycledim(dim-idim,0,dim-1));
						const Tensor<T> w=weights.reshape(kvec,kvec);

						const unsigned int idim0=(idim+1)%dim;

						// set up \Prod_{i\neq k} <G | F>
						Tensor<T> prod(rF,kvec);
						for (unsigned int r=0; r<rF; r++) {
							for (unsigned int i0=0; i0<kvec; i0++) {
								const T F1=trial._ptr->vector_[idim0](r,i0);
								prod(r,i0)=F1;
							}
						}
						const Tensor<T> p=prod.reshape(rF,kvec);

						// compute the contrib to b
						vecb+=madness::inner(p,w,-1,-1);
					}

				}


				// solve the linear system
				// note that gesv requires: vecb(kvec,rF) -> vecb(rF,kvec)p;
				// x can be empty for now
				Tensor<T> x;

				try {
#if !bench
					gesv(B[idim],vecb,x);
#else
					x=vecb;
					x=1.0;
#endif

				} catch (std::exception) {
					print("gesv failed..");
					print(B);
					print(vecb);
					MADNESS_ASSERT(0);
				}

				vecb=x;

				for (unsigned int l=0; l<rF; l++) {

					// calculate the new weights s_
					typename madness::Tensor<T>::float_scalar_type norm=0.0;
					for (unsigned int jk=0; jk<kvec; jk++) {
						const T val=vecb(l,jk);
						norm+= madness::detail::mynorm(val);
						//						print("norm in ALS", norm);
					}

					/*
					 * check for NaNs somewhere
					 */
					if (not (norm==norm)) {
						std::cout << "NaNs in ALS" << std::endl;
						print("B[idim]",B[idim]);
						//				vecbb.print("old vecb");
						std::cout << "idim " << idim << std::endl;
						std::cout << "weight_l in NaN; norm: " << norm << std::endl;
						MADNESS_EXCEPTION("NaNs in ALS",0);
					}
					MADNESS_ASSERT(norm>=0.0);

					/*
					 * if trial is a residual the weights might be zero;
					 * fill the vectors F with random numbers, not with zeros,
					 * otherwise they will screw up the ALS algorithm in the next
					 * iteration!
					 */
					double weight_l=sqrt(norm);
					//					print("weight in ALS", weight_l);
					if (norm==0.0) {
						weight_l=0.0;
						fvec.fillrandom();
						fvec=1.0;
						fvec.scale(0.01);
					} else {
						for (unsigned int jk=0; jk<kvec; jk++) {
							const T c_jk_l=vecb(l,jk);
							fvec(jk)=c_jk_l/weight_l;
						}
					}
					//					print("fvec in ALS", fvec);

					// use this->kVec(), as this might be an SVR vector or
					// an SPR operator
					trial._ptr->reassign(idim,l,weight_l,fvec,ref1.kVec());
				}
			}
		}

		/// release memory
		void clear() {_ptr.reset();};

		/// same as operator+=, but handles non-conforming vectors (i.e. slices)
		void inplace_add(const gentensorT& rhs, const std::vector<Slice>& lhs_s,
				const std::vector<Slice>& rhs_s, const double alpha, const double beta) {

			// fast return if possible
			if (rhs.rank()==0) return;

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
				rhs._ptr->make_structure();
				_ptr->make_structure();
				this->_ptr->inplace_add(*rhs._ptr,lhs_s,rhs_s, alpha, beta);
			}
		}

		/// inplace add
		void update_by(const gentensorT& rhs2) {
			if (this->rank()==0) {
				*this+=rhs2;
				return;
			}

			if (rhs2.rank()==1) {
				_ptr->rank_n_update_sequential(*rhs2._ptr);

			} else {
				gentensorT rhs3=copy(rhs2);
				rhs3._ptr->undo_structure();										// 0.3 s

				const long chunk_size=8;

				for (unsigned int i=0; i<rhs2.rank(); i+=chunk_size) {
					int begin=i;
					int end=std::min(rhs2.rank()-1,i+chunk_size-1);

					// very hard-wired
					configT rhs(rhs2.dim(),rhs2.get_k(),rhs2.tensor_type());
					rhs.weights_=rhs3._ptr->weights_(Slice(begin,end));
					rhs.ref_vector(0)=rhs3._ptr->ref_vector(0)(Slice(begin,end),Slice(_));
					rhs.ref_vector(1)=rhs3._ptr->ref_vector(1)(Slice(begin,end),Slice(_));
					rhs.rank_=end-begin+1;
					rhs.make_slices();

					rhs.orthonormalize();

					const tensorT a=(rhs.ref_vector(0)(rhs.c0()));
					const tensorT b=(rhs.ref_vector(1)(rhs.c0()));
					_ptr->rank_n_update_chunkwise(a,b,rhs.weights_(Slice(0,rhs.rank()-1)));
				}
			}


		}

		void finalize_accumulate() {_ptr->finalize_accumulate();}

		/// reduce the separation rank of this to a near optimal value
		/// follow section 3 in BM2005
		void doReduceRank(const double& eps, const Tensor<T>& values=Tensor<T>()){//,
			//				const SepRep& trial2=SepRep()) {
			/*
			 * basic idea is to use the residual Frobenius norm to check
			 * convergence. Don't know if this is rigorous, probably not..
			 *
			 * convergence criterion: optimize the trial SepRep until the
			 * residual doesn't change any more
			 */

			//			madness::print(values);

			/*
			 * figure out what to do:
			 * 	1. this exists and is to be reduced
			 * 	2. this doesn't exist, but values are provided
			 */
			gentensorT& reference=*this;
			const bool haveSR=(reference.rank()!=0);
			const bool haveVal=(values.ndim()!=-1);

			// set factors
			double facSR=0.0;
			if (haveSR) facSR=1.0;

			double facVal=0.0;
			if (haveVal) facVal=1.0;

			// fast return if possible
			if ((not haveSR) and (not haveVal)) return;


			//			reference.configs_=this->_ptr->semi_flatten();
			reference._ptr->undo_structure();

			//			const bool useTrial=(not (trial2.tensor_type()==TT_NONE));

			//			timeReduce_.start();

			/*
			 * some constants
			 */

			// the threshold
			const double threshold=eps*facReduce;

			// what we expect the trial rank might be (engineering problem)
			const unsigned int maxTrialRank=300;
			const unsigned int maxloop=300;

			const bool print=false;
			double norm=1.0;

			const unsigned int config_dim=this->_ptr->dim_eff();
			const unsigned int rG1=reference.rank();

			// scratch Tensor for this and the reference
			std::vector<Tensor<T> > B1(config_dim);
			for (unsigned int idim=0; idim<config_dim; idim++) B1[idim]=Tensor<T> (maxTrialRank,rG1);
			// scratch Tensor for this and the residual
			std::vector<Tensor<T> > B2(config_dim);
			for (unsigned int idim=0; idim<config_dim; idim++) B2[idim]=Tensor<T> (maxTrialRank,1);

			// set up a trial function
			gentensorT trial(reference.tensor_type(),reference.get_k(),reference.dim());
			trial._ptr->undo_structure();
			trial._ptr->reserve(10);
			//			trial._ptr->semi_flatten();
			//			if (useTrial) trial=trial2;
			//			else trial._ptr->ensureSpace(maxTrialRank);

			// and the residual
			gentensorT residual(trial.tensor_type(),trial.get_k(),trial.dim());

			// loop while || F-G || > epsilon
			for (unsigned int iloop=0; iloop<maxloop; iloop++) {

				// compute the residual wrt reference minus trial
				residual._ptr->fillWithRandom(1);
				residual._ptr->undo_structure();

				residual.optimize(reference,facSR,trial,-1.0,values,facVal,threshold,50,B1,B2);

				// exit if residual is supposedly small
				norm=residual.normf();
				if (print) printf("trial norm in reduceRank %d %12.8f\n", int(trial.rank()), norm);
#if bench
				if (iloop>5) break;
#else
				if (norm<threshold) break;
#endif
				// otherwise add residual to the trial function ..
				//				trial._ptr->unflatten();
				//				residual._ptr->unflatten();
				trial._ptr->make_structure();
				residual._ptr->make_structure();
				trial+=residual;
				trial._ptr->undo_structure();
				residual._ptr->undo_structure();
				//				trial._ptr->semi_flatten();
				//				residual._ptr->semi_flatten();

				// .. and optimize trial wrt the reference
				bool successful=trial.optimize(reference,facSR,residual,0.0,
						values,facVal,threshold,10,B1,B2);

				MADNESS_ASSERT(successful);

			}
			if (print) std::cout << "final trial norm in reduceRank " << trial.rank() << " " << norm  << std::endl;
			if (print) std::cout << "threshold " << threshold << std::endl;

#if !bench
			// check actual convergence
			if (norm>threshold) {
				//		trial.printCoeff("trial");
				//		this->printCoeff("failed SepRep");
				std::cout << "failed to reduce rank in SepRep::reduceRank() " << std::endl;
				printf("initial rank %d, trial rank %d\n", int(this->rank()), int(trial.rank()));
				printf("residual's norm         %24.16f\n", norm);
				printf("norm(this) %12.8f\n", this->normf());
				printf("no convergence in SepRep::reduceRank()");
				MADNESS_ASSERT(0);
			}
#endif

			// copy shrinks the matrices
			*this=copy(*trial._ptr);
			this->_ptr->make_structure();
			//			timeReduce_.end();

		}

		/// reduce the rank using SVD
		void computeSVD(const double& eps,const Tensor<T>& values_eff) {

			// SVD works only with matrices (2D)
			MADNESS_ASSERT(values_eff.ndim()==2);
			MADNESS_ASSERT(this->tensor_type()==TT_2D);

			// fast return if possible
			if (values_eff.normf()<eps*facReduce) {
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
			const double threshold=eps*eps*facReduce*facReduce;
			double residual=0.0;
			long i;
			for (i=s.dim(0)-1; i>=0; i--) {
				residual+=s(i)*s(i);
				if (residual>threshold) break;
			}

			// convert SVD output to our convention
			if (i>=0) {
				this->_ptr->weights_=s(Slice(0,i));
				this->_ptr->vector_[0]=transpose(U(Slice(_),Slice(0,i)));
				this->_ptr->vector_[1]=(VT(Slice(0,i),Slice(_)));
				this->_ptr->rank_=i+1;
				MADNESS_ASSERT(this->_ptr->kVec()==this->_ptr->vector_[0].dim(1));
				MADNESS_ASSERT(this->_ptr->rank()==this->_ptr->vector_[0].dim(0));
				MADNESS_ASSERT(this->_ptr->rank()==this->_ptr->weights_.dim(0));
			} else {
				_ptr=sr_ptr(new configT(dim(),get_k(),tensor_type()));
			}
			this->_ptr->make_structure();
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

		/// inplace addition
		SliceGenTensor<T>& operator+=(const SliceGenTensor<T>& rhs) {
			_refGT.inplace_add(*rhs._refGT._ptr,this->_s,rhs._s,1.0,1.0);
			return *this;
		}

		/// inplace zero-ing
		SliceGenTensor<T>& operator=(const double& number) {
			MADNESS_ASSERT(number==0.0);
			GenTensor<T> tmp=*this;
			if (this->_refGT.tensor_type()==TT_FULL) tmp=copy(tmp);
			tmp.scale(-1.0);
			_refGT.inplace_add(tmp,_s,_s,1.0,1.0);
			return *this;
		}

		/// for compatibility with tensor
		friend GenTensor<T> copy(const SliceGenTensor<T>& rhs) {
			return GenTensor<T>(rhs);
		}

	};


    /// Often used to transform all dimensions from one basis to another
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// \endcode
	/// cf tensor/tensor.h
    template <class T, class Q>
    GenTensor< TENSOR_RESULT_TYPE(T,Q) > transform(const GenTensor<Q>& t, const Tensor<T>& c) {
    	MADNESS_ASSERT(0);
//    	return t.transform(c);
    }

    template <class T>
    GenTensor<T> transform(const GenTensor<T>& t, const Tensor<T>& c) {
    	return t.transform(c);
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
				bool exist;
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



    #endif /* HAVE_GENTENSOR */


    /// transform the argument SepRepTensor to FullTensor form
    template <typename T>
    void to_full_rank(GenTensor<T>& arg) {

    	if (arg.has_data()) {
        	if (arg.tensor_type()==TT_FULL) {
        		;
        	} else if (arg.tensor_type()==TT_3D or arg.tensor_type()==TT_2D) {
        		Tensor<T> t=arg.reconstruct_tensor();
            	arg=GenTensor<T>(t,0.0,TT_FULL);
        	} else {
        		throw std::runtime_error("unknown TensorType in to_full_tensor");
        	}
    	}
    }

    /// transform the argument SepRepTensor to LowRankTensor form
    template <typename T>
    void to_low_rank(GenTensor<T>& arg, const double& eps, const TensorType& target_type) {

    	if (arg.has_data()) {
        	if (arg.tensor_type()==TT_FULL) {
				const Tensor<T> t1=arg.reconstruct_tensor();
				arg=(GenTensor<T>(t1,eps,target_type));
	     	} else if (arg.tensor_type()==TT_2D or arg.tensor_type()==TT_NONE) {
         		;
         	} else {
         		throw std::runtime_error("unknown TensorType in to_full_tensor");
         	}
    	}
    }







}

#endif /* GENTENSOR_H_ */
