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

	/*!
	 * A SRConf handles all the configurations in a Separated Representation.
	 *
	 *
	 */

	template <typename T>
	class SRConf {
	public:

		/// default ctor
		SRConf() : rank_(0), maxk_(0) {
		};


		/// ctor with dimensions for a vector configuration (tested)
		SRConf(const unsigned int& dim, const unsigned int& k)
			: rank_(0)
			, maxk_(k) {

			// construct empty vector
			vector_=std::vector<Tensor<T> > (dim);

		}

		/// copy ctor (tested); deep copy
		SRConf(const SRConf& rhs)  {
			*this=rhs;
		}


		/// ctor with provided weights and effective vectors; deep copy
		SRConf(const Tensor<double>& weights, const std::vector<Tensor<T> >& vectors) {

			// consistency check
			MADNESS_ASSERT(vectors.size()>0);
			MADNESS_ASSERT(weights.ndim()==1 and weights.dim(0)==vectors[0].dim(0));

			rank_=weights.dim(0);
			maxk_=vectors[0].dim(1);
			weights_=copy(weights);
			vector_=std::vector<Tensor<T> > (vectors.size());
			for (unsigned int idim=0; idim<vectors.size(); idim++) {
				vector_[idim]=copy(vectors[idim]);
			}

		}

		/// assignment operator (tested), deep copy of vectors
		SRConf& operator=(const SRConf& rhs)  {

			// check for self-assignment
			if (&rhs==this) return *this;

			// fast return iff rhs is invalid
			if (rhs.vector_.size()==0) {
				weights_=Tensor<double>();
				vector_.clear();
				maxk_=rhs.maxk_;
				rank_=0;
				return *this;
			}

			// fast return if possible; but construct a valid configuration
			if (rhs.rank()==0) {

				maxk_=rhs.maxk_;
				rank_=0;

				const unsigned int dim=rhs.dim();
				weights_=Tensor<double>();
				vector_=std::vector<Tensor<T> > (dim);

				return *this;
			}

			// assign vectors; deep copy
			vector_.resize(rhs.vector_.size());
			for (unsigned int i=0; i<rhs.vector_.size(); i++) {
				vector_[i]=copy(rhs.refVector(i));
			}

			// deep copy; note that rank is >=1 here
//			weights_=copy(rhs.weights_(Slice(0,rhs.rank()-1)));
			weights_=copy(rhs.weights_);

			rank_=rhs.rank();
			maxk_=rhs.maxk_;

			// consistency check
			for (unsigned int idim=0; idim<dim(); idim++) {
				assert(weights_.dim(0)==vector_[idim].dim(0));
				assert(vector_[idim].dim(1)==maxk_);
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
				*this=rhs;
				return *this;
			} else if (rhs.rank()==0) {
				return *this;
			}

			assert(compatible(*this,rhs));

			// for convenience
			const long oldRank=this->rank();
			const long rhsRank=rhs.rank();
			const long k=rhs.vector_[0].dim(1);

			// make sure we have enough space
			const long newRank=this->rank()+rhs.rank();
#if 1
			// assign weights
			Tensor<double> newWeights(newRank);
			newWeights(Slice(0,oldRank-1))=this->weights_(Slice(0,oldRank-1));
			newWeights(Slice(oldRank,newRank-1))=rhs.weights_(Slice(0,rhsRank-1));
			weights_=newWeights;

			// assign vectors
			for (unsigned int idim=0; idim<this->dim(); idim++) {
				Tensor<T> newVector(newRank,k);

				newVector(Slice(0,oldRank-1),Slice(_))=this->refVector(idim);
				newVector(Slice(oldRank,newRank-1),Slice(_))=rhs.refVector(idim);
				vector_[idim]=copy(newVector);

//				vector_[idim](Slice(oldRank,newRank-1),Slice(_))=rhs.refVector(idim);
			}


#else
			this->ensureSpace(newRank);


			// assign weights
			weights_(Slice(oldRank,newRank-1))=rhs.weights_(Slice(0,rhsRank-1));

			// assign vectors
			for (unsigned int idim=0; idim<this->dim(); idim++) {
				vector_[idim](Slice(oldRank,newRank-1),Slice(_))=rhs.refVector(idim);
			}
#endif
			rank_=newRank;
			return *this;

		}

		/// append an SRConf to this
//		SRConf& operator-=(const SRConf& rhs);

		/// reassign weight and data for a specific SRConf only
		void reassign(const unsigned int& idim, const unsigned int& r,
				const double& weight, const Tensor<T> & data, const unsigned int& maxk) {

			// some checks
	        assert(idim<this->dim());
	        assert(r<this->rank());
	        assert(this->isVector());

	        // assign weight
	        weights_(r)=weight;

	        // assign data
	        assert(data.size()==maxk);
	        for (unsigned int k=0; k<maxk; k++) {
	        	refVector(idim)(r,k)=data(k);
	        }

	        // make sure the ranks comply for all dimensions
	        for (unsigned int idim=0; idim<dim(); idim++) {
                assert(weights_.dim(0)==vector_[idim].dim(0));
	        }

		}

		/// cast this into a SRConf with different maxk
//		SRConf cast(const unsigned int maxk) const;

		/// perform pointwise multiplication
//		SRConf times(const SRConf& rhs, const unsigned int& maxk) const;

		/// transform this using the matrices and accumulate on target
//		void transform(SRConf& target, std::vector<const Tensor<T> *>& matrices, const double& dfac) const;

		/// transform this using the matrices
//		void selfTransform(std::vector<const Tensor<T> *>& matrices, const double& dfac);

		/// return const reference to one of the vectors F
		const Tensor<T>  refVector(const unsigned int& idim) const {
//				if (this->rank()==0) return Tensor<T>(0,vector_[0].dim(1));
				MADNESS_ASSERT(this->rank()==vector_[idim].dim(0));
				return vector_[idim];
//				return vector_[idim](Slice(0,this->rank()-1),Slice(_));
		}


		/// return reference to one of the vectors F
		Tensor<T>  refVector(const unsigned int& idim) {
			if (this->rank()==0) return Tensor<T>(0,vector_[0].dim(1));
			MADNESS_ASSERT(this->rank()==vector_[idim].dim(0));
			return vector_[idim];
//			return vector_[idim](Slice(0,this->rank()-1),Slice(_));
		}

		/// zero out
		void zeroOut() {MADNESS_ASSERT(0); this->chopoff(0);};

		/// fill this SRConf with 1 random configurations (tested)
		void fillWithRandom(const unsigned int& rank=1) {

//			this->ensureSpace(rank);
			rank_=rank;

			// assign; note that Slice(0,_) is inclusive
//			weights_(Slice(0,this->rank()-1))=1.0;
			weights_=Tensor<double>(rank);
			weights_=1.0;

			for (unsigned int idim=0; idim<this->dim(); idim++) {
//				vector_[idim](Slice(0,this->rank()-1),Slice(_)).fillrandom();
//				this->refVector(idim).fillrandom();
				vector_[idim]=Tensor<T>(rank_,maxk_);
				vector_[idim].fillrandom();
			}

			this->normalize();
			weights_(Slice(0,this->rank()-1))=1.0;
		}

		/// normalize the vectors (tested)
		void normalize() {

	        // for convenience
	        const unsigned int dim=this->dim();
	        const unsigned int rank=this->rank();

	        // some checks
	        assert(dim!=0);

	        // we calculate the norm sum_i < F^r_i | F^r_i > for each dimension for each r

	        // loop over all configurations
	        for (unsigned int r=0; r<rank; r++) {

	        	// loop over all dimensions
	        	for (unsigned int idim=0; idim<dim; idim++) {

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

		/// return if this is a vector or an operator (tested)
		bool isVector() const {
			return vector_[0].ndim()==2;
		}

		/// return the dimension of this
		unsigned int dim() const {return vector_.size();}

		/// return the logicalrank
		unsigned int rank() const {return rank_;};

		/// return the number of vector elements (tested)
		unsigned int get_k() const {return maxk_;};

		/// return the number of coefficients
		unsigned int nCoeff() const {
			return this->dim()*maxk_*rank_;
		};

		/// does what you think it does (tested)
	//	unsigned int get_kprime() const;

		/// reserve enough space to hold maxRank configurations
//		void reserve(const unsigned int& maxRank);


		/// make sure there is at least maxRank space in this
		void ensureSpace(const unsigned int& minRank) {
			MADNESS_ASSERT(0);

			// if the first (slow) dimension of the vectors is too small, reallocate
			// keep the second dimension
			if (weights_.dim(0)<minRank) {

				// for convenience
				const long k=vector_[0].dim(1);
				const long r=this->rank();

				// reallocate the weigths
				Tensor<double> tmp(minRank);
				if (r>0) tmp(Slice(0,r-1))=weights_(Slice(0,r-1));
				std::swap(weights_,tmp);

				// reallocate all vectors
				for (unsigned int idim=0; idim<this->dim(); idim++) {
					Tensor<T> t2(minRank,k);
					if (r>0) t2(Slice(0,r-1),Slice(_))=this->refVector(idim);
					std::swap(t2,vector_[idim]);
				}
			}
		}

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
			assert(rhs.isVector()==lhs.isVector());
			assert(rhs.dim()==lhs.dim());
			assert(rhs.dim()>0);

			// for now
			assert(rhs.isVector());

			// fast return if either rank is 0
			if ((lhs.rank()==0) or (rhs.rank()==0)) return 0.0;

			const unsigned int dim=rhs.dim();
//			const unsigned int k=rhs.vector_[0].dim(1);

			// get the weight matrix
			Tensor<T> weightMatrix=outer(lhs.weights_,rhs.weights_);

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<dim; idim++) {
//				dgemm_NT(lhs.refVector(idim),rhs.refVector(idim),ovlp,
//						lhs.rank(),rhs.rank());
				const Tensor<T>& lhs2=lhs.refVector(idim);
				const Tensor<T>& rhs2=rhs.refVector(idim);
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

		/// return a specific elements of this
		T getElement(const unsigned int& r, const unsigned int& idim, const unsigned int& k) const
		{assert(this->isVector()); assert(r<rank()); assert(idim<dim()); return refVector(idim)(r,k);};

		/// the weights and vectors
//		void print(const std::string& title="") const;

		/// check compatibility
		friend bool compatible(const SRConf& lhs, const SRConf& rhs) {
			return ((lhs.dim()==rhs.dim()) and (lhs.vector_[0].dim(1)==rhs.vector_[0].dim(1))
					and (lhs.isVector()==rhs.isVector()));
		}

		/// make one of the terms in the B matrix (BM2005)
		void friend makeB(Tensor<T> & B, const unsigned int& idim, const SRConf<T>& lhs, const SRConf<T>& rhs) {
			// some checks
			assert(compatible(rhs,lhs));
			assert(lhs.rank()==B.dim(0));
			assert(rhs.rank()==B.dim(1));
			assert(idim<rhs.dim());

//			dgemm_NT(lhs.refVector(idim),rhs.refVector(idim),B,lhs.rank(),rhs.rank());
			Tensor<T> lhs2=lhs.refVector(idim)(Slice(0,lhs.rank()-1),Slice(_));
			Tensor<T> rhs2=rhs.refVector(idim)(Slice(0,rhs.rank()-1),Slice(_));
			B=0.0;
			inner_result(lhs2,rhs2,-1,-1,B);
		}

	public:

		/// for each configuration the weight; length should be r
		Tensor<double>  weights_;

	public:
		/// for each (physical) dimension one Tensor of (logical) dimension (r,k)
		/// for vectors or (r,kprime,k) for operators
		std::vector<Tensor<T> > vector_;

		/// what is the rank of this
		long rank_;

		/// the number of underlying basis functions
		unsigned int maxk_;


	private:

		/// return the physical rank (the size of the Tensor)
//		unsigned int physicalRank() const {return weights_.pSize();}

		/// make a matrix of the weights (lhs.r, rhs.r)
		void friend makeWeightMatrix(Tensor<T> & matrix, const SRConf& lhs, const SRConf& rhs) {
			// some checks
			assert(compatible(rhs,lhs));

			matrix=outer(lhs.weights_,rhs.weights_);
		}

	};
}

#endif /* SRCONF_H_ */
