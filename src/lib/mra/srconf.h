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

#include "tensor/tensor.h"
#include "mra/funcdefaults.h"

namespace madness {

	/// return the number of vectors (i.e. dim_eff) according to the TensorType
	static unsigned int compute_nvec(const TensorType& tt) {
		if (tt==TT_2D) return 2;
		if (tt==TT_3D) return 3;
		print("unknown TensorType",tt);
		MADNESS_ASSERT(0);
	}



	/*!
	 * A SRConf handles all the configurations in a Separated Representation.
	 *
	 *
	 */

	template <typename T>
	class SRConf {

	public:

		/// the number of dimensions (the order of the tensor)
		unsigned int dim_;

		/// for each configuration the weight; length should be r
		Tensor<double>  weights_;

		/// for each (physical) dimension one Tensor of (logical) dimension (r,k)
		/// for vectors or (r,kprime,k) for operators
		std::vector<Tensor<T> > vector_;

		/// what is the rank of this
		long rank_;

		/// the number of underlying basis functions
		/// the dimensions of vector_ will be
		/// vector_(rank,maxk),
		/// vector_(rank,maxk,maxk), etc
		unsigned int maxk_;

		/// how will this be represented
		TensorType tensortype_;

	public:

		/// default ctor
		SRConf() : dim_(0), rank_(0), maxk_(0), tensortype_(TT_NONE) {
		};

		/// default ctor
		SRConf(const TensorType& tt) : dim_(0), rank_(0), maxk_(0), tensortype_(tt) {
		};


		/// ctor with dimensions for a vector configuration (tested)
		SRConf(const unsigned int& dim, const unsigned int& k, const TensorType& tt)
			: dim_(dim)
			, rank_(0)
			, maxk_(k)
			, tensortype_(tt) {

			// make sure dim is integer multiple of requested TT
			const unsigned int nvec=compute_nvec(tt);
			MADNESS_ASSERT(dim%nvec==0);

			// construct empty vector
			weights_=Tensor<double>(int(0));
			vector_=std::vector<Tensor<T> > (nvec);
			for (unsigned int idim=0; idim<nvec; idim++) vector_[idim]=Tensor<T>(0,this->kVec());

		}

		/// copy ctor (tested); shallow copy
		SRConf(const SRConf& rhs)  {
			*this=rhs;
		}


		/// ctor with provided weights and effective vectors; shallow copy
		SRConf(const Tensor<double>& weights, const std::vector<Tensor<T> >& vectors,
				const unsigned int& dim, const unsigned int maxk, const TensorType& tt)
			: dim_(dim)
			, maxk_(maxk)
			, tensortype_(tt) {

			// consistency check
			MADNESS_ASSERT(vectors.size()>0);
			MADNESS_ASSERT(weights.ndim()==1 and weights.dim(0)==vectors[0].dim(0));

			// compute dimension
			unsigned int nvec=compute_nvec(tt);
			MADNESS_ASSERT(vectors.size()==nvec);
			MADNESS_ASSERT(dim%nvec==0);

			rank_=weights.dim(0);
			weights_=weights;
			vector_=std::vector<Tensor<T> > (vectors.size());
			for (unsigned int idim=0; idim<vectors.size(); idim++) {
				vector_[idim]=vectors[idim];
			}
		}

		/// assignment operator (tested), shallow copy of vectors
		SRConf& operator=(const SRConf& rhs)  {

			// check for self-assignment
			if (&rhs==this) return *this;

			// these always hold
			dim_=rhs.dim_;
			tensortype_=rhs.tensortype_;
			maxk_=rhs.maxk_;

			// assign vectors; shallow copy
			vector_.resize(rhs.vector_.size());
			for (unsigned int i=0; i<rhs.vector_.size(); i++) {
				vector_[i]=(rhs.refVector(i));
			}

			// shallow copy
			weights_=(rhs.weights_);

			rank_=rhs.rank();

			// consistency check
			for (unsigned int idim=0; idim<dim_eff(); idim++) {
				assert(weights_.dim(0)==vector_[idim].dim(0));
			}

			return *this;
		}

		/// dtor
		virtual ~SRConf() {
			vector_.clear();
		}

		/// append an SRConf to this
		SRConf& operator+=(const SRConf& rhs) {

			// fast return if possible
			if (this->rank()==0) {
				*this=copy(rhs);
				return *this;
			} else if (rhs.rank()==0) {
				return *this;
			}

			assert(compatible(*this,rhs));

			// pass some dummy slice
			std::vector<Slice> s(this->dim(),Slice(_));
			this->inplace_add(rhs,s,s);
			return *this;

		}

		/// same as operator+=, but handles non-conforming vectors (i.e. slices)

		/// bounds checking should have been performed by caller
		/// s denotes where in lhs the new contribution from rhs will be inserted
		void inplace_add(const SRConf<T>& rhs2, std::vector<Slice> lhs_s,
				std::vector<Slice> rhs_s) {

			// fast return if possible; no fast return for this.rank()==0
			// since we might work with slices!
			if (rhs2.rank()==0) {
				return;
			}

			// unflatten this and rhs; shallow wrt vector_
			SRConf<T> lhs=this->unflatten();
			const SRConf<T> rhs=rhs2.unflatten();

			// for convenience
			const long lhsRank=lhs.rank();
			const long rhsRank=rhs.rank();
			const long newRank=lhs.rank()+rhs.rank();

			const long rhs_k=rhs.get_k();
			const long lhs_k=lhs.get_k();

			const long dim_pv=lhs.dim_per_vector();

			// adapt slices for use
			for (unsigned int idim=0; idim<lhs.dim(); idim++) {
				if (lhs_s[idim].end<0) lhs_s[idim].end+=lhs_k;
				if (rhs_s[idim].end<0) rhs_s[idim].end+=rhs_k;
				// make sure slices conform
				MADNESS_ASSERT((lhs_s[idim].end-lhs_s[idim].start) == (rhs_s[idim].end-rhs_s[idim].start));
				// make sure lhs can actually hold rhs(s)
				MADNESS_ASSERT(lhs_k>=(rhs_s[idim].end-rhs_s[idim].start+1));
			}

			// assign weights
			Tensor<double> newWeights(newRank);
			if (lhsRank>0) newWeights(Slice(0,lhsRank-1))=lhs.weights_(Slice(0,lhsRank-1));
			newWeights(Slice(lhsRank,newRank-1))=rhs.weights_(Slice(0,rhsRank-1));
			std::swap(lhs.weights_,newWeights);


			// assign vectors
			for (unsigned int idim=0; idim<lhs.dim_eff(); idim++) {

				if (dim_pv==1) {
					Tensor<T> newVector(newRank,lhs_k);

					// lhs unchanged
					if (lhsRank>0) newVector(Slice(0,lhsRank-1),Slice(_))=
							lhs.refVector(idim)(Slice(0,lhsRank-1),Slice(_));

					// insert rhs at the right place
					newVector(Slice(lhsRank,newRank-1),lhs_s[idim])=
							rhs.refVector(idim)(Slice(0,rhsRank-1),rhs_s[idim]);

					std::swap(lhs.refVector(idim),newVector);

				} else if (dim_pv==2) {

					Tensor<T> newVector(newRank,lhs_k,lhs_k);

					// lhs unchanged
					if (lhsRank>0) newVector(Slice(0,lhsRank-1),Slice(_),Slice(_))=
							lhs.refVector(idim)(Slice(0,lhsRank-1),Slice(_),Slice(_));

					// insert rhs at the right place
					newVector(Slice(lhsRank,newRank-1),lhs_s[2*idim],lhs_s[2*idim+1])=
							rhs.refVector(idim)(Slice(0,rhsRank-1),rhs_s[2*idim],rhs_s[2*idim+1]);

					std::swap(lhs.refVector(idim),newVector);

				} else if (dim_pv==3) {

					Tensor<T> newVector(newRank,lhs_k,lhs_k,lhs_k);

					// lhs unchanged
					if (lhsRank>0) newVector(Slice(0,lhsRank-1),Slice(_),Slice(_),Slice(_))=
							lhs.refVector(idim)(Slice(0,lhsRank-1),Slice(_),Slice(_),Slice(_));

					// insert rhs at the right place
					newVector(Slice(lhsRank,newRank-1),lhs_s[3*idim],lhs_s[3*idim+1],lhs_s[3*idim+2])=
							rhs.refVector(idim)(Slice(0,rhsRank-1),rhs_s[3*idim],rhs_s[3*idim+1],rhs_s[3*idim+2]);

					std::swap(lhs.refVector(idim),newVector);
				} else {
					print("extend dim_pv in srconf::inplace_add");
					MADNESS_ASSERT(0);
				}
			}

			lhs.rank_=newRank;
			*this=lhs.semi_flatten();
		}

		/// append an SRConf to this
//		SRConf& operator-=(const SRConf& rhs);

		/// deep copy of rhs
		friend SRConf<T> copy(const SRConf<T>& rhs) {

			// if rhs is non-existent simply construct a new SRConf
			if (rhs.rank()==0) return SRConf<T>(rhs.dim(),rhs.get_k(),rhs.type());

			// pass a copy of the weights and vectors of rhs to ctor
			std::vector<Tensor<T> > vector(rhs.dim_eff());
			for (unsigned int idim=0; idim<rhs.dim_eff(); idim++)
				vector[idim]=copy(rhs.refVector(idim));

			return SRConf<T>(copy(rhs.weights_),vector,rhs.dim(),rhs.get_k(),rhs.type());
		}

		/// reassign weight and data for a specific SRConf only
		void reassign(const unsigned int& idim, const unsigned int& r,
				const double& weight, const Tensor<T> & data, const unsigned int& maxk) {

			// some checks
	        assert(idim<this->dim_eff());
	        assert(r<this->rank());
	        assert(this->is_flat());

	        // assign weight
	        weights_(r)=weight;

	        // assign data
	        assert(data.size()==maxk);
	        for (unsigned int k=0; k<maxk; k++) {
	        	refVector(idim)(r,k)=data(k);
	        }

	        // make sure the ranks comply for all dimensions
	        for (unsigned int idim=0; idim<dim_eff(); idim++) {
                assert(weights_.dim(0)==vector_[idim].dim(0));
	        }

		}

		/// cast this into a SRConf with different maxk
//		SRConf cast(const unsigned int maxk) const;

		/// perform pointwise multiplication
//		SRConf times(const SRConf& rhs, const unsigned int& maxk) const;

		/// transform this using the matrices and accumulate on target
		/// note that the number of matrices must equal the physical dimension
		/// of this, not necessarily the number of vector_.

//		void fast_transform(SRConf& target,const Tensor<Q>& matrix, const double& dfac) const {
//
//			typedef TENSOR_RESULT_TYPE(T,Q) TQ;
//			typedef resultT Tensor<TQ>;
//
//			// some checks
//			MADNESS_ASSERT(target.dim_eff()==this->dim_eff());
//			MADNESS_ASSERT(matrices.size()==this->dim());
//
//			// workspace
//            Tensor<resultT> work(this->refVector(0).ndim(),this->refVector(0).dims(),false);
//
//
//			for (unsigned int idim=0; idim<this->ndim(); idim++) {
//				target.refVector(idim)
//			}
//		}

		/// transform this using the matrices
//		void selfTransform(std::vector<const Tensor<T> *>& matrices, const double& dfac);

		/// return const reference to one of the vectors F
		const Tensor<T>& refVector(const unsigned int& idim) const {
//				if (this->rank()==0) return Tensor<T>(0,vector_[0].dim(1));
				MADNESS_ASSERT(this->rank()==vector_[idim].dim(0));
				return vector_[idim];
//				return vector_[idim](Slice(0,this->rank()-1),Slice(_));
		}


		/// return reference to one of the vectors F
		Tensor<T>& refVector(const unsigned int& idim) {
//			if (this->rank()==0) return Tensor<T>(0,vector_[0].dim(1));
			MADNESS_ASSERT(this->rank()==vector_[idim].dim(0));
			return vector_[idim];
//			return vector_[idim](Slice(0,this->rank()-1),Slice(_));
		}

		/// zero out
		void zeroOut() {MADNESS_ASSERT(0); this->chopoff(0);};

		/// fill this SRConf with 1 random configurations (tested)
		void fillWithRandom(const unsigned int& rank=1) {

			rank_=rank;

			// assign; note that Slice(0,_) is inclusive
			weights_=Tensor<double>(rank);
			weights_=1.0;

			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				vector_[idim]=Tensor<T>(rank_,this->kVec());
				vector_[idim].fillrandom();
			}

			this->normalize();
			weights_(Slice(0,this->rank()-1))=1.0;
		}

		/// normalize the vectors (tested)
		void normalize() {

	        // for convenience
	        const unsigned int rank=this->rank();

	        // we calculate the norm sum_i < F^r_i | F^r_i > for each dimension for each r

	        // loop over all configurations
	        for (unsigned int r=0; r<rank; r++) {

	        	// loop over all dimensions
	        	for (unsigned int idim=0; idim<dim_eff(); idim++) {

	        		Tensor<T> config=this->refVector(idim)(Slice(r,r),Slice(_));

	        		const double norm=config.normf();
	        		const double fac=norm;
	        		double oofac=1.0/fac;
	        		if (fac==0.0) oofac=0.0;

	        		weights_(r)*=fac;
	        		config.scale(oofac);
	        	}
	        }
		}

		/// remove all configurations whose weight is smaller than the threshold
//		void truncate(const double& thresh);

		/// return if this has only one additional dimension (apart from rank)
		bool is_flat() const {
			return (vector_[0].ndim()==2);
		}

		/// return if this has a tensor structure (has not been flattened)
		bool has_structure() const {
			return (vector_[0].dim(1)==this->get_k());
		}

		/// return the dimension of this
		unsigned int dim() const {return dim_;}

		/// return the number of vectors
		unsigned int dim_eff() const {return vector_.size();}

		/// return the logicalrank
		unsigned int rank() const {return rank_;};

		/// return the number of physical matrix elements per dimension
		unsigned int get_k() const {return maxk_;};

		/// return the length of the vector (dim_pv*maxk)
		unsigned int kVec() const {return pow(this->get_k(),this->dim_per_vector());};

		/// return the tensor type
		TensorType type() const {return tensortype_;};

		/// return the number of physical dimensions
		int dim_per_vector() const {
			const int nvec=vector_.size();
			const int dim=this->dim();
			MADNESS_ASSERT(dim%nvec==0);
			return dim/nvec;
		}

		/// shallow flatten the vectors (reverse of unflatten)
		/// vector_[i](rank,maxk,maxk) -> vector_[i](rank,maxk*maxk)
		SRConf<T> semi_flatten() const {

			// shallow ctor
			SRConf<T> result(*this);
			const int dim_pv=this->dim_per_vector();
			const int maxk=this->get_k();

			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				if (dim_pv==1) ;
				else if (dim_pv==2)
					result.refVector(idim)=this->refVector(idim).reshape(this->rank(),maxk*maxk);
				else if (dim_pv==3)
					result.refVector(idim)=this->refVector(idim).reshape(this->rank(),maxk*maxk*maxk);
			}
			return result;
		}

		/// shallow unflatten the vectors (reverse of semi_flatten)
		/// vector_[i](rank,maxk*maxk) -> vector_[i](rank,maxk,maxk)
		SRConf unflatten() const {

			// shallow assignment
			SRConf<T> result(*this);

			const int dim_pv=this->dim_per_vector();
			const int maxk=this->get_k();

			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				if (dim_pv==1) ;
				else if (dim_pv==2)
					result.refVector(idim)=this->refVector(idim).reshape(this->rank(),maxk,maxk);
				else if (dim_pv==3)
					result.refVector(idim)=this->refVector(idim).reshape(this->rank(),maxk,maxk,maxk);
			}
			return result;
		}


		/// return the number of coefficients
		unsigned int nCoeff() const {
			return this->dim()*maxk_*rank_;
		};

		/// does what you think it does (tested)
	//	unsigned int get_kprime() const;

		/// reserve enough space to hold maxRank configurations
//		void reserve(const unsigned int& maxRank);


//		/// make sure there is at least maxRank space in this
//		void ensureSpace(const unsigned int& minRank) {
//			MADNESS_ASSERT(0);
//
//			// if the first (slow) dimension of the vectors is too small, reallocate
//			// keep the second dimension
//			if (weights_.dim(0)<minRank) {
//
//				// for convenience
//				const long k=vector_[0].dim(1);
//				const long r=this->rank();
//
//				// reallocate the weigths
//				Tensor<double> tmp(minRank);
//				if (r>0) tmp(Slice(0,r-1))=weights_(Slice(0,r-1));
//				std::swap(weights_,tmp);
//
//				// reallocate all vectors
//				for (unsigned int idim=0; idim<this->dim(); idim++) {
//					Tensor<T> t2(minRank,k);
//					if (r>0) t2(Slice(0,r-1),Slice(_))=this->refVector(idim);
//					std::swap(t2,vector_[idim]);
//				}
//			}
//		}

		/// chop off all configurations after r
//		void chopoff(const unsigned int& r);

		/// physically chop off all configurations after rank
//		void pchopoff();

		/// return the value at point r
	//	double value(const std::vector<double>& r, const double& pow2n,
	//			const double& vol, const MultiIndex& index, const unsigned int& maxk) const;

		/// calculate the Frobenius inner product (tested)
		friend T overlap(const SRConf<T>& rhs, const SRConf<T>& lhs) {

			/*
			 * the structure of an SRConf is (r,k) or (r,k',k), with
			 * r the slowest index; the overlap is therefore simply calculated
			 * as the matrix multiplication rhs*lhs^T
			 */

			// some checks
			assert(rhs.is_flat()==lhs.is_flat());
			assert(rhs.dim()==lhs.dim());
			assert(rhs.dim()>0);

			SRConf<T> rhs3=rhs.semi_flatten();
			SRConf<T> lhs3=lhs.semi_flatten();

			// fast return if either rank is 0
			if ((lhs.rank()==0) or (rhs.rank()==0)) return 0.0;

			const unsigned int dim_eff=rhs.dim_eff();
//			const unsigned int k=rhs.vector_[0].dim(1);

			// get the weight matrix
			Tensor<T> weightMatrix=outer(lhs.weights_,rhs.weights_);

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<dim_eff; idim++) {
//				dgemm_NT(lhs.refVector(idim),rhs.refVector(idim),ovlp,
//						lhs.rank(),rhs.rank());
				const Tensor<T>& lhs2=lhs3.refVector(idim);
				const Tensor<T>& rhs2=rhs3.refVector(idim);
				Tensor< TENSOR_RESULT_TYPE(T,T) > ovlp(lhs.rank(),rhs.rank());
//				inner_result(lhs.refVector(idim),rhs,refVector(idim),-1,-1,ovlp);
				inner_result(lhs2,rhs2,-1,-1,ovlp);

			    // multiply all overlap matrices with the weight matrix
				weightMatrix.emul(ovlp);
			}

		//	return weightMatrix;
			const T overlap=weightMatrix.sum();
			return overlap;
		}

		/// scale this by a number
		void scale(const double& fac) {weights_.scale(fac);};

		/// return the weight
		double weights(const unsigned int& i) const {return weights_(i);};

		/// return the maximum weight
		double maxWeight() const {return weights_(Slice(0,this->rank()-1)).max();};

//		/// return a specific elements of this
//		T getElement(const unsigned int& r, const unsigned int& idim, const unsigned int& k) const
//		{assert(this->isVector()); assert(r<rank()); assert(idim<dim()); return refVector(idim)(r,k);};

		/// the weights and vectors
//		void print(const std::string& title="") const;

		/// check compatibility
		friend bool compatible(const SRConf& lhs, const SRConf& rhs) {
			return ((lhs.dim()==rhs.dim()) and (lhs.vector_[0].dim(1)==rhs.vector_[0].dim(1))
					and (lhs.is_flat()==rhs.is_flat()));
		}

		/// make one of the terms in the B matrix (BM2005)
		void friend makeB(Tensor<T> & B, const unsigned int& idim, const SRConf<T>& lhs, const SRConf<T>& rhs) {
			// some checks
			assert(compatible(rhs,lhs));
			assert(lhs.rank()==B.dim(0));
			assert(rhs.rank()==B.dim(1));
			assert(idim<rhs.dim_eff());

//			dgemm_NT(lhs.refVector(idim),rhs.refVector(idim),B,lhs.rank(),rhs.rank());
			Tensor<T> lhs2=lhs.refVector(idim)(Slice(0,lhs.rank()-1),Slice(_));
			Tensor<T> rhs2=rhs.refVector(idim)(Slice(0,rhs.rank()-1),Slice(_));
			B=0.0;
			inner_result(lhs2,rhs2,-1,-1,B);
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		SRConf<T> transform(const Tensor<T>& c) const {

			SRConf<T> result=copy(*this);

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			// these two loops go over all physical dimensions (dim = dim_eff * merged_dim)
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				for (unsigned int jdim=1; jdim<this->refVector(idim).ndim(); jdim++) {

					// note tricky ordering (jdim is missing): this is actually correct!
					result.refVector(idim)=madness::inner(result.refVector(idim),c,1,0);

				}
			}
			return result;
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		SRConf<T> general_transform(const Tensor<T> c[]) const {

			SRConf<T> result=copy(*this);

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			long i=0;
			// these two loops go over all physical dimensions (dim = dim_eff * merged_dim)
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				for (unsigned int jdim=1; jdim<this->refVector(idim).ndim(); jdim++) {

					// note tricky ordering (jdim is missing): this is actually correct!
					result.refVector(idim)=madness::inner(result.refVector(idim),c[i],1,0);
					i++;

				}
			}
			return result;
		}


		SRConf<T> transform_dir(const Tensor<T>& c, const int& axis) const {

			// only a matrix is allowed for c
			MADNESS_ASSERT(c.ndim()==2);

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			// compute idim for accessing the vector_, and the dimension inside vector_
			// the +1 on jdim for the rank
			const long idim=axis/this->dim_per_vector();
			const long jdim=axis%this->dim_per_vector()+1;

			SRConf<T> result=copy(*this);
			result.refVector(idim)=madness::transform_dir(this->refVector(idim),c,jdim);

			return result;
		}


	public:


	private:

//		/// return the physical rank (the size of the Tensor)
//		unsigned int physicalRank() const {return weights_.pSize();}
//
//		/// make a matrix of the weights (lhs.r, rhs.r)
//		void friend makeWeightMatrix(Tensor<T> & matrix, const SRConf& lhs, const SRConf& rhs) {
//			// some checks
//			assert(compatible(rhs,lhs));
//
//			matrix=outer(lhs.weights_,rhs.weights_);
//		}

	};


}

#endif /* SRCONF_H_ */
