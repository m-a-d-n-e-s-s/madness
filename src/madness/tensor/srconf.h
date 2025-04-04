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

/// \file SRConf.h
/// \brief handles the low-level details of a separated representation tensor

#ifndef SRCONF_H_
#define SRCONF_H_

//#define BENCH 0

#include <madness/world/print.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/clapack.h>
#include <madness/tensor/tensor_lapack.h>
#include <list>
#include<array>

namespace madness {



#ifdef BENCH
	static double time_[30];
#endif


	template <class T> class GenTensor;
    template <class T> class GenTensor;
    template <class T> class SliceLowRankTensor;

	/**
	 * A SRConf handles all the configurations in a Separated Representation.
	 */

	template <typename T>
	class SRConf : public BaseTensor {
		friend class GenTensor<T>;
        friend class SliceLowRankTensor<T>;

		/// the scalar type of T
		typedef typename Tensor<T>::scalar_type scalar_type;
		typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;
	public:

#ifdef BENCH
		static double& time(int i) {return time_[i];}
#endif

		typedef Tensor<T> tensorT;

		/// check orthonormality at low rank additions
		static const bool check_orthonormality=false;

		/// for each configuration the weight; length should be r
		Tensor< typename Tensor<T>::scalar_type >  weights_;

		/// for each (physical) dimension one Tensor of (logical) dimension (r,k)
		/// for vectors or (r,kprime,k) for operators
		std::array<Tensor<T>,2> vector_;

		/// separation dimensions: A(n,m) -> A(r,n) B(r,m), with n={k1,k2},m={k3,k4,k5..) multi-indices
		long nci_left=-1;

		/// Slice containing the actual data in each vector, ignoring "empty" configurations;
		/// will maintain contiguity of the data.
		std::vector<Slice> s0,s1;
		
	public:

    	/// return the index of the last singular vector/value to meet the threshold
		///        (returns -1 if all meet threshold, i.e. || A ||_2 < threshold)
    	/// given a matrix A in SVD form, truncate the singular values such that the
    	/// accuracy threshold is still met.
    	/// @param[in]	thresh	the threshold eps: || A - A(truncated) || < eps
    	/// @param[in] 	rank	the number of singular values in w
    	/// @param[in]	w		the weights/singular values of A
    	/// @return		i		the index of s_max to contribute: w(Slice(0,i)); i.e. inclusive!
    	static int max_sigma(const double& thresh, const long& rank, const Tensor<double>& w) {

    	    if (thresh<0.0) return rank-1;
    		// find the maximal singular value that's supposed to contribute
    		// singular values are ordered (largest first)
    		double residual=0.0;
    		long i;
    		for (i=rank-1; i>=0; i--) {
    			residual+=w(i)*w(i);
    			if (residual>thresh*thresh) break;
    		}
    		return i;
    	}

		/// default ctor
		SRConf() : nci_left(-1) {};

		SRConf(const long& ndim, const long* dimensions,
				const long nci) : nci_left(nci) {
    		BaseTensor::set_dims_and_size(ndim,dimensions);
			if (nci_left<0) nci_left=ndim/2;	// integer division
			make_structure();
		}

		SRConf(const long& ndim, const std::array<long,TENSOR_MAXDIM>& dimensions,
				const long nci) : SRConf(ndim,dimensions.data(),nci) {
		}

		/// copy ctor (tested); shallow copy
		SRConf(const SRConf& rhs) = default;

		/// ctor with provided weights and effective vectors; shallow copy
		SRConf(const Tensor<double>& weights, const std::vector<Tensor<T> >& vectors,
				const long& ndim, const long& dims, const long nci)
			: SRConf(ndim,dims,nci) {
			MADNESS_ASSERT(vectors.size()==2);
			set_vectors_and_weights(weights,vectors[0],vectors[1]);
			make_structure();
            MADNESS_ASSERT(has_structure());
			MADNESS_ASSERT(check_dimensions());
		}

		/// explicit ctor with two vectors (aka SVD), shallow
		SRConf(const Tensor<double>& weights, const tensorT& vector1, const tensorT& vector2,
				const long& ndim, const long* dims, const long nci)
			: SRConf(ndim,dims,nci) {
			set_vectors_and_weights(weights,vector1,vector2);
			make_structure();
			MADNESS_ASSERT(check_dimensions());
		}

		/// assignment operator (tested), shallow copy of vectors
		SRConf& operator=(const SRConf& rhs)  {

			// check for self-assignment
			if (&rhs==this) return *this;
			if (rhs.has_no_data()) {
				clear();
				return *this;
			}

			// these always hold
			nci_left=rhs.nci_left;
    		BaseTensor::set_dims_and_size(rhs.ndim(),rhs.dims());
			s0=rhs.s0;
			s1=rhs.s1;

			if (rhs.rank()==0) {
				// construct empty vector
				make_empty_vectors_and_weights(0);
				make_structure();

			} else {
				// assign vectors; shallow copy
				for (unsigned int i=0; i<rhs.vector_.size(); i++) {
					vector_[i]=rhs.vector_[i];
				}

				// shallow copy
				weights_=(rhs.weights_);

			}
            MADNESS_ASSERT(has_structure());
			return *this;
		}

		/// assign a number to this;
		SRConf& operator=(const T& number) {

			// rank will be one
			this->make_empty_vectors_and_weights(1);
			vector_[0]=number;
			vector_[1]=1.0;
			weights_(0l)=1.0;
			make_structure();
			normalize();
		    return *this;
		}


    	void set_size_and_dim(long ndim, long k) {
    		std::array<long,TENSOR_MAXDIM> dims;
    		dims.fill(k);
    		BaseTensor::set_dims_and_size(ndim,dims.data());
    	}

    	/// deduce the dimensions of the left and right singular vectors from the tensor dimensions
    	std::array<std::array<long,TENSOR_MAXDIM>, 2> make_vector_dimensions(const long rank) const {
    		std::array<std::array<long,TENSOR_MAXDIM>, 2> dimensions;
			for (int i=0; i<nci_left; ++i) dimensions[0][i+1]=this->dim(i);
			for (int i=nci_left; i<ndim(); ++i) dimensions[1][i+1-nci_left]=this->dim(i);

			dimensions[0][0]=rank;
			dimensions[1][0]=rank;
			return dimensions;
    	}

    	void make_empty_vectors_and_weights(const long rank) {
    		auto dimensions=make_vector_dimensions(rank);
    		weights_=Tensor<float_scalar_type>(rank);
    		vector_[0]=Tensor<T>(nci_left+1,dimensions[0].data());
    		vector_[1]=Tensor<T>(ndim()-nci_left+1,dimensions[1].data());
    		make_structure();
    	}

    	void set_vectors_and_weights(const Tensor< typename Tensor<T>::scalar_type >&  weights,
    			const Tensor<T>& vector1, const Tensor<T>& vector2) {
    		weights_=weights;
    		vector_={vector1,vector2};
    		make_structure();
    	}

    	void clear() {
    		weights_.clear();
    		vector_[0].clear();
    		vector_[1].clear();
            _size = 0;
            _ndim = -1;
    	}

		/// return some of the terms of the SRConf (start,..,end), inclusively
		/// shallow copy
		const SRConf get_configs(const int& start, const int& end) const {

			MADNESS_ASSERT((start>=0) and (end<=rank()));

			Slice s(start,end);
			tensorT v0,v1;

			// slice vectors
			v0=flat_vector(0)(s,_);
			v1=flat_vector(1)(s,_);

			SRConf<T> result(ndim(),dims(),nci_left);
			result.set_vectors_and_weights(weights_(s),v0,v1);
            MADNESS_ASSERT(result.has_structure());
			return result;
		}

		/// dtor
		~SRConf() {}

        template <typename Archive>
        void serialize(Archive& ar) {
              	ar & weights_ & vector_[0] & vector_[1] & nci_left & _ndim & _size
					& _id &  archive::wrap(_dim,TENSOR_MAXDIM) & s0 & s1;
//              	make_slices();
                MADNESS_ASSERT(has_structure());
        }

		/// does this have any data?
		bool has_data() const {
			return (size()!=0);
		}

		/// does this have any data?
		bool has_no_data() const {return !has_data();}

		/// return a Slice that corresponds the that part of vector_ that holds coefficients
		const std::vector<Slice>& c0(const int idim) const {
			if (idim==0) return s0;
			else if (idim==1) return s1;
			else {
				MADNESS_EXCEPTION("invalid idim in SRConf::idim",1);
			}
		}

		/// append configurations of rhs to this

		/// simplified version of inplace_add for flattened configurations
		/// *this += fac*rhs
		void append(const SRConf<T>& rhs, const double fac=1.0) {

			// fast return if possible
			if (rhs.has_no_data() or rhs.rank()==0) return;
			if (this->has_no_data() or rank()==0) {
				*this=copy(rhs);
				this->scale(fac);
				return;
			}

    		auto dimensions=make_vector_dimensions(rank()+rhs.rank());
			Tensor<float_scalar_type> weights(rank()+rhs.rank());
			Tensor<T> vector0(nci_left+1,dimensions[0].data());
			Tensor<T> vector1(ndim()-nci_left+1,dimensions[1].data());

			// assign weights
			weights(Slice(0,rank()-1))=weights_(Slice(0,rank()-1));
			weights(Slice(rank(),rank()+rhs.rank()-1))=rhs.weights_(Slice(0,rhs.rank()-1))*fac;

			vector0(c0(0))=vector_[0](c0(0));
			vector1(c0(1))=vector_[1](c0(1));

			auto s00=s0;
			auto s10=s1;
			s00[0]=Slice(rank(),rank()+rhs.rank()-1);
			s10[0]=Slice(rank(),rank()+rhs.rank()-1);

			vector0(s00)=rhs.vector_[0](rhs.c0(0));
			vector1(s10)=rhs.vector_[1](rhs.c0(1));

			std::swap(weights,weights_);
			std::swap(vector0,vector_[0]);
			std::swap(vector1,vector_[1]);

			make_slices();
            MADNESS_ASSERT(has_structure());

		}

		void append(const SRConf<T>& rhs, const double_complex fac=1.0) {
			MADNESS_EXCEPTION("no complex in SRConf",1);
		}

	public:
		/// add two orthonormal configurations, yielding an optimal SVD decomposition
		void add_SVD(const SRConf<T>& rhs, const double& thresh) {
#ifdef BENCH
			double cpu0=wall_time();
#endif
			if (rhs.has_no_data() or rhs.rank()==0) return;
			if (has_no_data() or rank()==0) {
				*this=rhs;
				return;
			}

			if (check_orthonormality) check_right_orthonormality();
            if (check_orthonormality) rhs.check_right_orthonormality();

//            Tensor<T> x1=flat_vector(0);
//            Tensor<T> x2=flat_vector(1);
            Tensor<T> x1=vector_[0].reshape(rank(),kVec(0));
            Tensor<T> x2=vector_[1].reshape(rank(),kVec(1));
            ortho5(x1,x2,weights_,
					rhs.flat_vector(0),rhs.flat_vector(1),rhs.weights_,thresh);
            std::swap(x1,vector_[0]);
            std::swap(x2,vector_[1]);

			make_structure();
			make_slices();
            MADNESS_ASSERT(has_structure());
#ifdef BENCH
			double cpu1=wall_time();
			time(25)+=cpu1-cpu0;
#endif
		}

	protected:
		/// alpha * this(lhs_s) + beta * rhs(rhs_s)

		/// bounds checking should have been performed by caller
		/// s denotes where in lhs the new contribution from rhs will be inserted
		void inplace_add(const SRConf<T>& rhs, std::array<Slice,TENSOR_MAXDIM> lhs_s,
				std::array<Slice,TENSOR_MAXDIM> rhs_s, const double alpha, const double beta) {

			// cannot scale a slice only...
			MADNESS_ASSERT(alpha==1.0);

			// fast return if possible; no fast return for this.rank()==0
			// since we might work with slices!
			if (rhs.has_no_data() or rhs.rank()==0) return;

			// prepare the vectors
			SRConf<T> result(ndim(),dims(),nci_left);
			result.make_empty_vectors_and_weights(rank()+rhs.rank());

			// insert lhs into result
			if (rank()>0) {
				result.vector_[0](s0)=vector_[0];
				result.vector_[1](s1)=vector_[1];
				result.weights_(Slice(0,rank()-1))=weights_;
			}
			// insert rhs into result
			{
	            auto [sr0,sr1]=rhs.make_slices(rhs_s);
	            auto [sl0,sl1]=make_slices(lhs_s);
	            sl0[0]=Slice(rank(),result.rank()-1);
	            sl1[0]=Slice(rank(),result.rank()-1);
	            result.vector_[0](sl0)=rhs.vector_[0](sr0);
	            result.vector_[1](sl1)=rhs.vector_[1](sr1);
	            result.weights_(Slice(rank(),result.rank()-1))=rhs.weights_*beta;

			}
			std::swap(*this,result);
            MADNESS_ASSERT(has_structure());
		}

		/// deep copy of rhs, shrink
		friend SRConf<T> copy(const SRConf<T>& rhs) {

			if (rhs.has_no_data()) return SRConf<T>();

			SRConf<T> result(rhs.ndim(),rhs.dims(),rhs.nci_left);

			// if rhs is non-existent simply construct a new SRConf
			if (rhs.has_data() and rhs.rank()>0) {
				result.set_vectors_and_weights(copy(rhs.weights_(Slice(0,rhs.rank()-1))),
						copy(rhs.vector_[0](rhs.c0(0))),copy(rhs.vector_[1](rhs.c0(1))));
			}

			return result;
		}

public:
        /// return a slice of this (deep copy)
        SRConf<T> copy_slice(const std::array<Slice,TENSOR_MAXDIM>& s) const {

        	std::array<long,TENSOR_MAXDIM> k;
        	for (int i=0; i<s.size(); ++i) {
        		if (s[i].end==-1) k[i]=dim(i);
        		else k[i]= s[i].end-s[i].start+1;
        	}
            SRConf<T> result(ndim(),k,nci_left);

            // fast return if possible
            if (this->has_no_data() or rank()==0) return result;

            auto [s00,s11]=make_slices(s);
            Tensor<T> vector0=copy(vector_[0](s00));
            Tensor<T> vector1=copy(vector_[1](s11));

            Tensor<double> weights=copy(this->weights_(Slice(0,rank()-1)));
            result.set_vectors_and_weights(weights,vector0,vector1);

            return result;
        }

        /// perform elementwise Hadamard product
        SRConf<T>& emul(const SRConf<T>& other) {
            // consistency check
            MADNESS_ASSERT(compatible(*this,other));

            long finalrank=this->rank()*other.rank();

            SRConf<T> result(ndim(),dims(),nci_left);	// empty tensor

            if ((this->rank()==0) or (other.rank()==0)) {
                ;   // pass
            } else {

            	result.make_empty_vectors_and_weights(finalrank);
                result.weights_=outer(weights_,other.weights_).flat();

                // left vector
                for (int i=0; i<2; ++i) {
                    Tensor<T> a1=flat_vector(i);
                    Tensor<T> b1=other.flat_vector(i);
                	Tensor<T> r1=result.vector_[i].reshape(finalrank,kVec(i));
                    Tensor<T> tmp(finalrank,1);
                    for (int k=0; k<a1.dim(1); ++k) {
//                    	r1(_,Slice(k,k))=outer(a1(_,k),b1(_,k)).reshape(finalrank,1);
                        outer_result(a1(_,k),b1(_,k),tmp);
                        r1(_,Slice(k,k))=tmp;
                    }
                }
            }
            result.make_structure();
            result.normalize();

            *this=result;
            return *this;
        }

protected:
		/// redo the Slices for getting direct access to the configurations
		void make_slices() {
			s0=std::vector<Slice>(nci_left+1,_);				// first dim is the rank
			s1=std::vector<Slice>(ndim()-nci_left+1,_);				// first dim is the rank

			s0[0]=Slice(0,rank()-1);
			s1[0]=Slice(0,rank()-1);
		}

		std::array<std::vector<Slice>,2> make_slices(const std::array<Slice,TENSOR_MAXDIM>& s) const {
			std::vector<Slice> s00(nci_left+1,_);				// first dim is the rank
			std::vector<Slice> s10(ndim()-nci_left+1,_);				// first dim is the rank

			s00[0]=Slice(0,rank()-1);
			s10[0]=Slice(0,rank()-1);

			for (int idim=0; idim<ndim(); ++idim) {
				if (idim<nci_left) s00[idim+1]=s[idim];
				if (idim>=nci_left) s10[idim-nci_left+1]=s[idim];
			}

			std::array<std::vector<Slice>,2> result={s00,s10};
			return result;
		}

		void make_structure(bool force=false) {

			auto dimensions=make_vector_dimensions(rank());
			vector_[0]=vector_[0].reshape(nci_left+1,&dimensions[0][0]);
			vector_[1]=vector_[1].reshape(ndim()-nci_left+1,&dimensions[1][0]);

			this->make_slices();

		}

	public:
		/// return reference to one of the vectors F
		Tensor<T>& ref_vector(const unsigned int& idim) {
			return vector_[idim];
		}

		/// return reference to one of the vectors F
		const Tensor<T>& ref_vector(const unsigned int& idim) const {
			return vector_[idim];
		}


		long kVec(const int idim) const {
			return vector_[idim].size()/vector_[idim].dim(0);
		}

		/// return shallow copy of a slice of one of the vectors, flattened to (r,kVec)
		const Tensor<T> flat_vector(const unsigned int& idim) const {
		    MADNESS_ASSERT(rank()>0);
//		    const Tensor<T> result=vector_[idim](c0(idim)).reshape(rank(),kVec(idim));
		    const Tensor<T> result=vector_[idim].reshape(rank(),kVec(idim));
		    return result;
		}

		/// return shallow copy of a slice of one of the vectors, flattened to (r,kVec)
		Tensor<T> flat_vector(const unsigned int& idim) {
		    MADNESS_ASSERT(rank()>0);
//		    Tensor<T> result=vector_[idim](c0(idim)).reshape(rank(),kVec(idim));
		    Tensor<T> result=vector_[idim].reshape(rank(),kVec(idim));
		    return result;
		}

	protected:
//		/// fill this SRConf with 1 flattened random configurations (tested)
//		void fillWithRandom(const long& rank=1) {
//
//
//			// assign; note that Slice(0,_) is inclusive
//			weights_=Tensor<double>(rank);
//			weights_=1.0;
//
//			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
//				vector_[idim]=Tensor<T>(rank,this->kVec());
//				vector_[idim].fillrandom();
//			}
//
//			this->normalize();
//			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
//				vector_[idim].scale(madness::RandomValue<T>()*scalar_type(10.0));
//			}
//			weights_(Slice(0,this->rank()-1)).fillrandom().scale(10.0);
//			make_slices();
//            MADNESS_ASSERT(has_structure());
//		}

		/// normalize the vectors (tested)
		void normalize() {

			if (rank()==0) return;
            MADNESS_ASSERT(has_structure());

	        // for convenience
	        const unsigned int rank=this->rank();
//	        std::vector<Slice> s(2,_);

	        // we calculate the norm sum_i < F^r_i | F^r_i > for each dimension for each r

	        // loop over all configurations
	        for (unsigned int r=0; r<rank; r++) {
	        	// loop over all dimensions
	        	for (unsigned int idim=0; idim<2; idim++) {
	    	        std::vector<Slice> s(dim_per_vector(idim)+1,_);
		            s[0]=Slice(r,r);

//	        		const double norm=this->ref_vector(idim)(s).normf();
	        		const double norm=this->vector_[idim](s).normf();
	        		const double fac=norm;
	        		double oofac=1.0/fac;
	        		if (fac<1.e-13) oofac=0.0;
	        		weights_(r)*=fac;
	        		vector_[idim](s).scale(oofac);
	        	}
	        }
            MADNESS_ASSERT(has_structure());
		}

		bool check_dimensions() const {
			bool correct=true;
			if (vector_[0].dim(0)!=rank()) return false;
			if (vector_[0].ndim()+vector_[1].ndim()!=ndim()+2) return false;
			if (vector_[0].ndim()!=nci_left+1) return false;
			if (vector_[1].ndim()!=ndim()-nci_left+1) return false;

			return correct;
		}

		/// check if the terms are orthogonal
		bool check_right_orthonormality() const {

			// fast return if possible
			if (rank()==0) return true;

			const tensorT t1=ref_vector(1)(c0(1)).reshape(rank(),kVec(1));
			tensorT S=inner(t1,t1,1,1);
			for (int i=0; i<S.dim(0); i++) S(i,i)-=1.0;

			// error per matrix element
			double norm=S.normf();
			double small=sqrt(norm*norm/S.size());
			return (small<1.e-13);
		}

		/// return if this has only one additional dimension (apart from rank)
		bool is_flat() const {
			return (vector_[0].ndim()==2);
		}

	public:
		/// return if this has a tensor structure (has not been flattened)
		bool has_structure() const {
			if (vector_.size()==2) {
				auto vector_dimensions=make_vector_dimensions(vector_[0].dim(0));
				for (int i=0; i<vector_[0].ndim(); ++i) {
					if (vector_dimensions[0][i] != vector_[0].dim(i)) return false;
				}
				for (int i=0; i<vector_[1].ndim(); ++i) {
					if (vector_dimensions[1][i] != vector_[1].dim(i)) return false;
				}
			}
			return true;
		}

		/// return the logicalrank
		long rank() const {return weights_.size();};

	public:
		/// return the number of physical dimensions
		int dim_per_vector(int idim) const {
                    MADNESS_ASSERT(vector_.size()>size_t(idim));
			return vector_[idim].ndim()-1;		// remove dimension for the rank
		}

		/// return the weight
		double weights(const unsigned int& i) const {return weights_(i);};

        /// reconstruct this to return a full tensor
        Tensor<T> reconstruct() const {

            // fast return if possible
            if (rank()==0) return Tensor<T> (ndim(),dims(),true);

            // include weights in left vector
            Tensor<T> scr=make_left_vector_with_weights();

			Tensor<T> result=inner(conj(scr),flat_vector(1),0,0);
            return result.reshape(ndim(),dims());
        }

	protected:

        Tensor<T> make_left_vector_with_weights() const {
            return make_vector_with_weights(0);
        }

    public:
        Tensor<T> make_vector_with_weights(const int dim) const {
            Tensor<T> v=copy(vector_[dim].reshape(rank(),vector_[dim].size()/rank()));
            for (unsigned int r=0; r<rank(); r++) v(r,_)*=weights(r);
            v=v.reshape(vector_[dim].ndim(),vector_[dim].dims());
            return v;
        }

        /// return flat (r,i) view of the tensor with the weights multiplied in

        /// return a(r,i) = vec(dim)(r,i) * w(r)
        Tensor<T> flat_vector_with_weights(const int dim) const {
            return make_vector_with_weights(dim).reshape(rank(),vector_[dim].size()/rank());
        }

	protected:
		/// return the number of coefficients
		unsigned int nCoeff() const {
			return vector_[0].size()+vector_[1].size()+weights_.size();
		};

		/// return the real size of this
		size_t real_size() const {
			size_t n=0;
			for (size_t i=0; i<vector_.size(); ++i) {
				n+=vector_[i].size()*sizeof(T) + sizeof(tensorT);
			}
			n+=weights_.size()*sizeof(double) + sizeof(Tensor<double>);
			n+=sizeof(*this);
			n+=2*s1.size()*sizeof(Slice);
			return n;
		}


		template<typename Q>
	    typename std::enable_if<(TensorTypeData<T>::iscomplex or TensorTypeData<Q>::iscomplex), TENSOR_RESULT_TYPE(T,Q)>::type
		friend  trace(const SRConf<T>& rhs, const SRConf<Q>& lhs) {
			MADNESS_EXCEPTION("no complex trace in srconf.h",1);
			return T(0.0);
		}

		/// calculate the Frobenius inner product (tested)
		template<typename Q>
	    typename std::enable_if<!(TensorTypeData<T>::iscomplex or TensorTypeData<Q>::iscomplex) , TENSOR_RESULT_TYPE(T,Q)>::type
		friend  trace(const SRConf<T>& rhs, const SRConf<Q>& lhs) {

			// fast return if either rank is 0
			if ((lhs.has_no_data()) or (rhs.has_no_data())) return 0.0;
			if ((lhs.rank()==0) or (rhs.rank()==0)) return 0.0;

			/*
			 * the structure of an SRConf is (r,k) or (r,k',k), with
			 * r the slowest index; the overlap is therefore simply calculated
			 * as the matrix multiplication rhs*lhs^T
			 */

			// some checks
			MADNESS_ASSERT(rhs.ndim()==lhs.ndim());
			MADNESS_ASSERT(rhs.ndim()>0);

			typedef TENSOR_RESULT_TYPE(T,Q) resultT;

			const unsigned int dim_eff=2;

			// get the weight matrix
			Tensor<resultT> weightMatrix=outer(lhs.weights_(Slice(0,lhs.rank()-1)),
					rhs.weights_(Slice(0,rhs.rank()-1)));

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<dim_eff; idim++) {
				const Tensor<T> lhs2=lhs.flat_vector(idim);
				const Tensor<Q> rhs2=rhs.flat_vector(idim);
				Tensor<resultT> ovlp(lhs.rank(),rhs.rank());
				inner_result(lhs2,rhs2,-1,-1,ovlp);

			    // multiply all overlap matrices with the weight matrix
				weightMatrix.emul(ovlp);
			}

			//	return weightMatrix;
			const TENSOR_RESULT_TYPE(T,Q) overlap=weightMatrix.sum();
			return overlap;
		}

		/// calculate the Frobenius norm, if this is in SVD form
        typename TensorTypeData<T>::float_scalar_type svd_normf() const {
            if (has_no_data() or rank()==0) return 0.0;
            return weights_(Slice(0,rank()-1)).normf();
        }

		/// calculate the Frobenius norm
		typename TensorTypeData<T>::float_scalar_type normf() const {

			// fast return if possible
			if (has_no_data() or rank()==0) return 0.0;

			// some checks
			MADNESS_ASSERT(ndim()>0);
			MADNESS_ASSERT(not TensorTypeData<T>::iscomplex);

			// get the weight matrix
			Tensor<T> weightMatrix=outer(weights_(Slice(0,rank()-1)),weights_(Slice(0,rank()-1)));

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<2; idim++) {
				const Tensor<T> vec=flat_vector(idim);
				Tensor<T> ovlp(rank(),rank());
				inner_result(vec,vec,-1,-1,ovlp);

			    // multiply all overlap matrices with the weight matrix
				weightMatrix.emul(ovlp);
			}

			typedef typename TensorTypeData<T>::float_scalar_type resultT;
			const resultT overlap=std::abs(weightMatrix.sum());
			return sqrt(overlap);
		}

	protected:
		/// scale this by a number
		void scale(const double& fac) {weights_.scale(fac);};

		void scale(const double_complex& fac) {
			MADNESS_EXCEPTION("no complex scaling in SRConf",1);
		}

		/// check compatibility
		friend bool compatible(const SRConf& lhs, const SRConf& rhs) {
			return (lhs.conforms(&rhs) and (lhs.dim_per_vector(0)==rhs.dim_per_vector(0)));
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		SRConf<T> transform(const Tensor<T>& c) const {

			// fast return if possible
			if (this->has_no_data() or rank()==0) return SRConf<T>(ndim(),dims(),nci_left);

			// transpose both singular vectors from U(r,i,j,k) to U(i,j,k,r)
			// and run the contraction as in tensor: U(i,j,k,r) t(i,i') = U(j,k,r,i')
			// and so on, yielding U(r,i',j',k')
			Tensor<T> left=copy(vector_[0].cycledim(-1,0,-1));
			Tensor<T> right=copy(vector_[1].cycledim(-1,0,-1));

	        for (long i=0; i<nci_left; ++i) left = inner(left,c,0,0);
	        for (long i=nci_left; i<ndim(); ++i) right = inner(right,c,0,0);

	        SRConf<T> result(ndim(),dims(),this->nci_left);
	        result.set_vectors_and_weights(copy(weights_),left,right);

            MADNESS_ASSERT(result.has_structure());
			return result;
		}
public:
	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		template<typename Q>
		SRConf<TENSOR_RESULT_TYPE(T,Q) > general_transform(const Tensor<Q> c[]) const {

			// fast return if possible
			if (this->has_no_data() or rank()==0) return SRConf<T>(ndim(),dims(),nci_left);

			// copying shrinks the vectors to (r,k,k,..)
			MADNESS_ASSERT(this->has_structure());

			// transpose both singular vectors from U(r,i,j,k) to U(i,j,k,r)
			// and run the contraction as in tensor: U(i,j,k,r) t(i,i') = U(j,k,r,i')
			// and so on, yielding U(r,i',j',k')
			Tensor<T> left=copy(vector_[0].cycledim(-1,0,-1));
			Tensor<T> right=copy(vector_[1].cycledim(-1,0,-1));

	        for (long i=0; i<nci_left; ++i) left = inner(left,c[i],0,0);
	        for (long i=nci_left; i<ndim(); ++i) right = inner(right,c[i],0,0);

	        SRConf<T> result(ndim(),dims(),this->nci_left);
	        result.set_vectors_and_weights(copy(weights_),left,right);

            MADNESS_ASSERT(result.has_structure());
			return result;
		}

		SRConf<T> transform_dir(const Tensor<T>& c, const int& axis) const {

			if (this->has_no_data() or rank()==0) return SRConf<T>(ndim(),dims(),nci_left);

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);

			// only a matrix is allowed for c
			MADNESS_ASSERT(c.ndim()==2);

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			// compute idim for accessing the vector_, and the dimension inside vector_
			// the +1 on jdim for the rank
			const long idim=axis/this->dim_per_vector(0);
			const long jdim=axis%this->dim_per_vector(0)+1;

			// note: no slicing necessary, since we've copied this to result (incl shrinking)
			result.ref_vector(idim)=madness::transform_dir(this->ref_vector(idim),c,jdim);
            MADNESS_ASSERT(result.has_structure());

			return result;
		}

	};

	/// sophisticated version of ortho2

	/// after calling this we will have an optimally rank-reduced representation
	/// with the left and right subspaces being bi-orthogonal and normalized;
	/// outline of the algorithm:
	///  - canonical orthogonalization of the subspaces (screen for small eigenvalues)
	///  - SVD of the modified overlap (incorporates the roots of eigenvalues)
	/// operation count is O(kr^2 + r^3)
	///
	/// @param[in,out]	x normalized left subspace
	/// @param[in,out]	y normalize right subspace
	/// @param[in,out]	weights weights
	/// @param[in]		thresh	truncation threshold
	template<typename T>
	void ortho3(Tensor<T>& x, Tensor<T>& y, Tensor<typename Tensor<T>::scalar_type>& weights, const double& thresh) {

#ifdef BENCH
		double cpu0=wall_time();
#endif
		typedef Tensor<T> tensorT;
		typedef typename Tensor<T>::scalar_type scalar_type;

		const long rank=x.dim(0);
		const double w_max=weights.absmax()*rank;		// max Frobenius norm


		// overlap of 1 and 2
		tensorT S1=inner(x,x,1,1);
		tensorT S2=inner(y,y,1,1);	// 0.5 / 2.1
#ifdef BENCH
		double cpu1=wall_time();
		SRConf<T>::time(1)+=(cpu1-cpu0);
#endif

//	    print("norm(S1)",S1.normf());
//	    print("norm(S2)",S2.normf());

		// diagonalize
		tensorT U1, U2;
		Tensor<scalar_type> e1, e2;
	    syev(S1,U1,e1);
	    syev(S2,U2,e2);										// 2.3 / 4.0
#ifdef BENCH
		double cpu3=wall_time();
		SRConf<T>::time(3)+=cpu3-cpu1;
#endif

	    const scalar_type e1_max=e1.absmax();
	    const scalar_type e2_max=e2.absmax();

		// fast return if possible
		if ((e1_max*w_max<thresh) or (e2_max*w_max<thresh)) {
			x.clear();
			y.clear();
			weights.clear();
			return;
		}

	    // remove small negative eigenvalues
	    e1.screen(1.e-13);
	    e2.screen(1.e-13);
	    Tensor<scalar_type> sqrt_e1(rank), sqrt_e2(rank);


	    // shrink U1, U2
	    int lo1=0;
	    int lo2=0;
	    for (unsigned int r=0; r<rank; r++) {
	    	if (e1(r)*w_max<thresh) lo1=r+1;
	    	if (e2(r)*w_max<thresh) lo2=r+1;
	    	sqrt_e1(r)=sqrt(std::abs(e1(r)));
	    	sqrt_e2(r)=sqrt(std::abs(e2(r)));
	    }

	    U1=U1(Slice(_),Slice(lo1,-1));
	    U2=U2(Slice(_),Slice(lo2,-1));
	    sqrt_e1=sqrt_e1(Slice(lo1,-1));
	    sqrt_e2=sqrt_e2(Slice(lo2,-1));
	    unsigned int rank1=rank-lo1;
	    unsigned int rank2=rank-lo2;						// 0.0 / 0.0


	    MADNESS_ASSERT(sqrt_e1.size()==rank1);
	    MADNESS_ASSERT(sqrt_e2.size()==rank2);
#ifdef BENCH
		double cpu4=wall_time();
		SRConf<T>::time(4)+=cpu4-cpu3;
#endif


	    // set up overlap M; include X+
	    tensorT M(rank1,rank2);
	    for (unsigned int i=0; i<rank1; i++) {
	    	for (unsigned int j=0; j<rank2; j++) {
	    		for (unsigned int r=0; r<rank; r++) {
		    		M(i,j)+=U1(r,i)*sqrt_e1(i)*weights(r)*U2(r,j) * sqrt_e2(j);
//			    		M(i,j)+=U1(r,i)*weights(r)*U2(r,j);
	    		}
	    	}
	    }


	    // include X-
    	for (unsigned int r=0; r<rank1; r++) {
    		scalar_type fac=1.0/sqrt_e1(r);
    		for (unsigned int t=0; t<rank; t++) {
	    		U1(t,r)*=fac;
//	    		if (sqrt_e1(r)<thresh) throw;
    		}
    	}

	   	for (unsigned int r=0; r<rank2; r++) {
	   		scalar_type fac=1.0/sqrt_e2(r);
    		for (unsigned int t=0; t<rank; t++) {
	    		U2(t,r)*=fac;
//	    		if (sqrt_e2(r)<thresh) throw;
	    	}
	    }													// 0.2 / 0.6
#ifdef BENCH
		double cpu5=wall_time();
		SRConf<T>::time(5)+=cpu5-cpu4;
#endif

	    // decompose M
		tensorT Up,VTp;
		Tensor<scalar_type> Sp;
		svd(M,Up,Sp,VTp);									// 1.5 / 3.0
#ifdef BENCH
		double cpu6=wall_time();
		SRConf<T>::time(6)+=cpu6-cpu5;
#endif

		// make transformation matrices
		Up=inner(Up,U1,0,1);
		VTp=inner(VTp,U2,1,1);



//		// find the maximal singular value that's supposed to contribute
//		// singular values are ordered (largest first)
//		double residual=0.0;
//		long i;
//		for (i=Sp.dim(0)-1; i>=0; i--) {
//			residual+=Sp(i)*Sp(i);
//			if (residual>thresh*thresh) break;
//		}
		long i=SRConf<T>::max_sigma(thresh,Sp.dim(0),Sp);

#ifdef BENCH
		double cpu7=wall_time();
		SRConf<T>::time(7)+=cpu7-cpu6;
#endif

//		i=std::min(i,long(0));

	    Up=Up(Slice(0,i),Slice(_));
	    VTp=VTp(Slice(0,i),Slice(_));


		// convert SVD output to our convention
		if (i>=0) {

			// transform 1 and 2
		    x=inner(Up,x,1,0);
		    y=inner(VTp,y,1,0);				// 0.5 / 2.5
		    weights=Sp(Slice(0,i));

		} else {
			x.clear();
			y.clear();
			weights.clear();
		}
#ifdef BENCH
		double cpu8=wall_time();
		SRConf<T>::time(8)+=cpu8-cpu7;
		SRConf<T>::time(0)+=cpu8-cpu0;
#endif
		return;
	}

	/// specialized version of ortho3

	/// does the same as ortho3, but takes two bi-orthonormal configs as input
	/// and saves on the inner product. Result will be written onto the first config
	///
	/// @param[in,out]	x1	left subspace, will hold the result on exit
	/// @param[in,out]	y1	right subspace, will hold the result on exit
	/// @param[in]		x2	left subspace, will be accumulated onto x1
	/// @param[in]		y2	right subspace, will be accumulated onto y1
	template<typename T>
	void ortho5(Tensor<T>& x1, Tensor<T>& y1, Tensor<typename Tensor<T>::scalar_type>& w1,
				const Tensor<T>& x2, const Tensor<T>& y2, const Tensor<typename Tensor<T>::scalar_type>& w2,
				const double& thresh) {

#ifdef BENCH
		double cpu0=wall_time();
#endif
		typedef Tensor<T> tensorT;

		const long rank1=x1.dim(0);
		const long rank2=x2.dim(0);
		const long rank=rank1+rank2;

		// for convenience: blocks of the matrices
		const Slice s0(0,rank1-1), s1(rank1,rank-1);

		const double w_max=std::max(w1.absmax(),w2.absmax());
		const double norm_max=w_max*rank;		// max Frobenius norm

		// the overlap between 1 and 2;
		// the overlap of 1 and 1, and 2 and 2 is assumed to be the identity matrix
		tensorT Sx12=inner(x1,x2,1,1);
		tensorT Sy12=inner(y1,y2,1,1);
#ifdef BENCH
		double cpu1=wall_time();
		SRConf<T>::time(11)+=cpu1-cpu0;
#endif

		tensorT Sx(rank,rank);
		tensorT Sy(rank,rank);

		// the identity matrix (half of it)
		for (long i=0; i<rank; i++) {
			Sx(i,i)=0.5;
			Sy(i,i)=0.5;
		}
		Sx(s0,s1)=Sx12;
		Sy(s0,s1)=Sy12;
		Sx+=transpose(Sx);
		Sy+=transpose(Sy);

		// overlap of 1 and 2
//		tensorT S1=inner(x,x,1,1);
//		tensorT S2=inner(y,y,1,1);	// 0.5 / 2.1
#ifdef BENCH
		double cpu2=wall_time();
		SRConf<T>::time(12)+=cpu2-cpu1;
#endif

		// diagonalize
		tensorT U1, U2;
		Tensor<typename Tensor<T>::scalar_type> e1, e2;
	    syev(Sx,U1,e1);
	    syev(Sy,U2,e2);										// 2.3 / 4.0
#ifdef BENCH
		double cpu3=wall_time();
		SRConf<T>::time(13)+=cpu3-cpu2;
#endif

//	    print("norm(Sx)",Sx.normf());
//	    print("norm(Sy)",Sy.normf());

	    const double e1_max=e1.absmax();
	    const double e2_max=e2.absmax();

		// fast return if possible
		if ((e1_max*norm_max<thresh) or (e2_max*norm_max<thresh)) {
			x1.clear();
			y1.clear();
			w1.clear();
			return;
		}

	    // remove small negative eigenvalues
	    e1.screen(1.e-13);
	    e2.screen(1.e-13);
	    Tensor<double> sqrt_e1(rank), sqrt_e2(rank);


	    // shrink U1, U2
	    int lo1=0;
	    int lo2=0;
	    for (unsigned int r=0; r<rank; r++) {
	    	if (e1(r)<thresh/norm_max) lo1=r+1;
	    	else sqrt_e1(r)=sqrt(std::abs(e1(r)));
	    	if (e2(r)<thresh/norm_max) lo2=r+1;
	    	else sqrt_e2(r)=sqrt(std::abs(e2(r)));
	    }

	    U1=U1(Slice(_),Slice(lo1,-1));
	    U2=U2(Slice(_),Slice(lo2,-1));
	    sqrt_e1=sqrt_e1(Slice(lo1,-1));
	    sqrt_e2=sqrt_e2(Slice(lo2,-1));
	    unsigned int rank_x=rank-lo1;
	    unsigned int rank_y=rank-lo2;						// 0.0 / 0.0


//	    MADNESS_ASSERT(sqrt_e1.size()==rank_x);
//	    MADNESS_ASSERT(sqrt_e2.size()==rank_y);

	    // set up overlap M; include X+

//	    for (unsigned int i=0; i<rank_x; ++i) U(i,_)*=sqrt_e1(i);
//	    for (unsigned int i=0; i<rank_y; ++i) U(i,_)*=sqrt_e2(i);

	    tensorT UU1=copy(U1);
	    for (unsigned int i=0; i<rank1; ++i) UU1(i,_)*=w1(i);
	    for (unsigned int i=rank1; i<rank; ++i) UU1(i,_)*=w2(i-rank1);

	    tensorT M=inner(UU1,U2,0,0);
	    tensorT ee=outer(sqrt_e1,sqrt_e2);
	    M.emul(ee);

	    // include X-
    	for (unsigned int r=0; r<rank_x; r++) {
    		double fac=1.0/sqrt_e1(r);
    		U1(_,r)*=fac;
    	}

	   	for (unsigned int r=0; r<rank_y; r++) {
    		double fac=1.0/sqrt_e2(r);
    		U2(_,r)*=fac;
	    }													// 0.2 / 0.6
#ifdef BENCH
		double cpu4=wall_time();
		SRConf<T>::time(14)+=cpu4-cpu3;
#endif


	    // decompose M
		tensorT Up,VTp;
		Tensor<typename Tensor<T>::scalar_type> Sp;
		svd(M,Up,Sp,VTp);									// 1.5 / 3.0
#ifdef BENCH
		double cpu5=wall_time();
		SRConf<T>::time(15)+=cpu5-cpu4;
#endif

		// make transformation matrices
		Up=inner(Up,U1,0,1);
		VTp=inner(VTp,U2,1,1);
#ifdef BENCH
		double cpu6=wall_time();
		SRConf<T>::time(16)+=cpu6-cpu5;
#endif

		// find the maximal singular value that's supposed to contribute
		// singular values are ordered (largest first)
		double residual=0.0;
		long i;
		for (i=Sp.dim(0)-1; i>=0; i--) {
			residual+=Sp(i)*Sp(i);
			if (residual>thresh*thresh) break;
		}

		// convert SVD output to our convention
		if (i>=0) {

			// make it contiguous
		    tensorT Up1=transpose(Up(Slice(0,i),s0));
		    tensorT Up2=transpose(Up(Slice(0,i),s1));
		    tensorT VTp1=transpose(VTp(Slice(0,i),s0));
		    tensorT VTp2=transpose(VTp(Slice(0,i),s1));

			// transform 1 and 2
		    x1=inner(Up1,x1,0,0);
		    inner_result(Up2,x2,0,0,x1);
		    y1=inner(VTp1,y1,0,0);
		    inner_result(VTp2,y2,0,0,y1);
		    w1=Sp(Slice(0,i));

		} else {
			x1.clear();
			y1.clear();
			w1.clear();
		}
#ifdef BENCH
		double cpu7=wall_time();
		SRConf<T>::time(17)+=cpu7-cpu6;
		SRConf<T>::time(10)+=cpu7-cpu0;
#endif
		return;
	}

	template<typename T>
	static inline
	std::ostream& operator<<(std::ostream& s, const SRConf<T>& sr) {

		s << "dim_          " << sr.ndim() << "\n";;
		s << "rank_         " << sr.rank_ << "\n";;
		s << "vector_.size()" << sr.vector_.size() << "\n";
		s << "has_data()    " << sr.has_data() << "\n";
		s << "TensorType    " << sr.type() << "\n\n";
		return s;
	}
}

#endif /* SRCONF_H_ */
