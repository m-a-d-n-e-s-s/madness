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

		typedef Tensor<T> tensorT;

		/// the number of dimensions (the order of the tensor)
		unsigned int dim_;

		/// for each configuration the weight; length should be r
		Tensor<double>  weights_;

		/// for each (physical) dimension one Tensor of (logical) dimension (r,k)
		/// for vectors or (r,kprime,k) for operators
		std::vector<tensorT> vector_;

		/// what is the rank of this
		long rank_;

		/// the number of underlying basis functions
		/// the dimensions of vector_ will be
		/// vector_(rank,maxk),
		/// vector_(rank,maxk,maxk), etc
		unsigned int maxk_;

		/// Slice containing the actual data, ignoring "empty" configurations
		std::vector<Slice> s_;

		/// how will this be represented
		TensorType tensortype_;

	public:

		/// default ctor
		SRConf() : dim_(0), rank_(0), maxk_(0), s_(), tensortype_(TT_NONE) {
		};

		/// default ctor
		SRConf(const TensorType& tt) : dim_(0), rank_(0), maxk_(0), s_(), tensortype_(tt) {
		};


		/// ctor with dimensions for a vector configuration (tested)
		SRConf(const unsigned int& dim, const unsigned int& k, const TensorType& tt)
			: dim_(dim)
			, rank_(0)
			, maxk_(k)
			, s_()
			, tensortype_(tt) {

			// make sure dim is integer multiple of requested TT
			const unsigned int nvec=compute_nvec(tt);
			MADNESS_ASSERT(dim%nvec==0);

			// construct empty vector
			weights_=Tensor<double>(int(0));
			vector_=std::vector<Tensor<T> > (nvec);
			for (unsigned int idim=0; idim<nvec; idim++) vector_[idim]=Tensor<T>(0,this->kVec());
			make_structure();
		}


		/// ctor with dimensions for a vector configuration (tested)
		SRConf(const unsigned int& dim, const unsigned int& k, const unsigned int& rank, const TensorType& tt)
			: dim_(dim)
			, rank_(0)
			, maxk_(k)
			, tensortype_(tt) {

			// make sure dim is integer multiple of requested TT
			const unsigned int nvec=compute_nvec(tt);
			MADNESS_ASSERT(dim%nvec==0);

			// construct empty vector
			weights_=Tensor<double>(rank);
			vector_=std::vector<Tensor<T> > (nvec);
			for (unsigned int idim=0; idim<nvec; idim++) vector_[idim]=Tensor<T>(rank,this->kVec());
			make_structure();
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
			make_slices();
		}

		/// assignment operator (tested), shallow copy of vectors
		SRConf& operator=(const SRConf& rhs)  {

			// check for self-assignment
			if (&rhs==this) return *this;

			// these always hold
			dim_=rhs.dim_;
			tensortype_=rhs.tensortype_;
			maxk_=rhs.maxk_;
			s_=rhs.s_;

			// assign vectors; shallow copy
			vector_.resize(rhs.vector_.size());
			for (unsigned int i=0; i<rhs.vector_.size(); i++) {
//				vector_[i]=(rhs.refVector_struct(i));
				vector_[i]=rhs.vector_[i];
			}

			// shallow copy
			weights_=(rhs.weights_);

			rank_=rhs.rank();

			// consistency check
			for (unsigned int idim=0; idim<dim_eff(); idim++) {
				MADNESS_ASSERT(weights_.dim(0)==vector_[idim].dim(0));
			}

			return *this;
		}

		/// dtor
		virtual ~SRConf() {
			vector_.clear();
			s_.clear();
		}

		/// reserve enough space to hold at least r configurations
		void reserve(const long r) {

			// this should at least hold the current information
			MADNESS_ASSERT(r>=this->rank());
			MADNESS_ASSERT(this->rank()==0 or vector_.size()>0);

			// fast return if possible
			// nothing to be done
			if (r==0) return;
			// already large enuff?
			if (this->vector_[0].dim(0)>=r) return;

			// for convenience
			const long rank=this->rank();
			const long kvec=this->kVec();
			this->undo_structure();

			// transfer weights
			Tensor<double> newWeights(r);
			if (rank>0) newWeights(Slice(0,rank-1))=weights_(Slice(0,rank-1));
			std::swap(weights_,newWeights);

			// transfer vectors
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {

				tensorT newVector(r,kvec);
				if (rank>0) newVector(this->c0())=vector_[idim](this->c0());
				std::swap(vector_[idim],newVector);

			}
			MADNESS_ASSERT(weights_.dim(0)==vector_[0].dim(0));
			this->make_structure();

		}

		/// return a Slice that corresponds the that part of vector_ that holds coefficients
		const std::vector<Slice>& c0() const {
			MADNESS_ASSERT(s_.size()>0);
			return s_;
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

			MADNESS_ASSERT(compatible(*this,rhs));

			// pass some dummy slice
			std::vector<Slice> s(this->dim(),Slice(_));
			this->inplace_add(rhs,s,s,1.0,1.0);
			return *this;

		}


		/// rank-1 update of this as in:	 *this += alpha * rhs

		/// M. Brand, Linear Algebra Appl 2006 vol. 415 (1) pp. 20-30
		void rank1_update_slow(const SRConf<T>& rhs2, const double& alpha) {

			// works only for SVD
			MADNESS_ASSERT(this->dim_eff()==2);
			if (rhs2.rank()!=1) print("rhs rank != 1");
			MADNESS_ASSERT(rhs2.rank()==1);

			SRConf<T> rhs=rhs2;
			rhs.undo_structure();
			this->undo_structure();

			const long rank=this->rank();
			const long kvec=this->kVec();

			// use the language of the article
			const tensorT a=(rhs.ref_vector(0)(0,Slice(_)))*rhs.weights(0);
			const tensorT bT=(rhs.ref_vector(1)).reshape(1,kvec);
			const tensorT b=bT.reshape(kvec);
			const tensorT UT=this->ref_vector(0)(this->c0());
			const tensorT VT=this->ref_vector(1)(this->c0());

			const tensorT U=transpose(UT);
			const tensorT V=transpose(VT);

			// eq (6)
			const tensorT m=inner(UT,a);
			const tensorT p=a-inner(U,m);
			const double Ra=p.normf();
			const tensorT P=p*(1.0/Ra);

			// eq (7)
			const tensorT n=inner(VT,b);
			const tensorT q=b-inner(V,n);
			const double Rb=q.normf();
			const tensorT Q=q*(1.0/Rb);

//			print("m,n");
//			print(m);
//			print(n);
//
//			print("a,b");
//			print(a);
//			print(b);
//
//			print("p,q");
//			print(p);
//			print(q);
//
//			print("P,Q");
//			print(P);
//			print(Q);
//
//			print("Ra,Rb",Ra,Rb);

			// eq (8)
			tensorT mp(rank+1);
			mp(Slice(0,rank-1))=m;
			mp(rank)=Ra;
			tensorT nq(rank+1);
			nq(Slice(0,rank-1))=n;
			nq(rank)=Rb;
			tensorT K=outer(mp,nq);
			for (long i=0; i<rank; i++) {
				K(i,i)+=this->weights(i);
			}

			// diagonalize K
			tensorT Up,Sp,VTp;
			svd(K,Up,Sp,VTp);
			tensorT UTp=transpose(Up);
			tensorT Vp=transpose(VTp);



			// rotate U, VT
			tensorT UP(kvec,rank+1);
			tensorT VQ(kvec,rank+1);
			UP(Slice(_),Slice(0,rank-1))=U;
			UP(Slice(_),rank)=P;
			VQ(Slice(_),Slice(0,rank-1))=V;
			VQ(Slice(_),rank)=Q;

			tensorT U_new=inner(UP,Up);
			tensorT V_new=inner(VQ,Vp);

			// insert new vectors
//			print("weights");
//			print(weights_);
//			print("Sp");
//			print(Sp);
			weights_=Sp;
			vector_[0]=transpose(U_new);
			vector_[1]=transpose(V_new);
			rank_+=1;

		}


		/// alpha * this(lhs_s) + beta * rhs(rhs_s)

		/// bounds checking should have been performed by caller
		/// s denotes where in lhs the new contribution from rhs will be inserted
		void inplace_add(const SRConf<T>& rhs2, std::vector<Slice> lhs_s,
				std::vector<Slice> rhs_s, const double alpha, const double beta) {

			// fast return if possible; no fast return for this.rank()==0
			// since we might work with slices!
			if (rhs2.rank()==0) return;


			// unflatten this and rhs; shallow wrt vector_
			SRConf<T>& lhs=*this;
			const SRConf<T>& rhs=rhs2;
			MADNESS_ASSERT(lhs.has_structure() or (lhs.rank()==0));
			MADNESS_ASSERT(rhs.has_structure());

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

			lhs.reserve(newRank);

			// assign weights, and include factors alpha and beta
			if (alpha!=1.0) lhs.scale(alpha);
			lhs.weights_(Slice(lhsRank,newRank-1))=rhs.weights_(Slice(0,rhsRank-1))*beta;


			// assign vectors
			for (unsigned int idim=0; idim<lhs.dim_eff(); idim++) {

				// insert rhs at the right place
				if (dim_pv==1) {
					lhs.ref_vector(idim)(Slice(lhsRank,newRank-1),lhs_s[idim])=
							rhs.ref_vector(idim)(Slice(0,rhsRank-1),rhs_s[idim]);

				} else if (dim_pv==2) {
					lhs.ref_vector(idim)(Slice(lhsRank,newRank-1),lhs_s[2*idim],lhs_s[2*idim+1])=
							rhs.ref_vector(idim)(Slice(0,rhsRank-1),rhs_s[2*idim],rhs_s[2*idim+1]);

				} else if (dim_pv==3) {
					lhs.ref_vector(idim)(Slice(lhsRank,newRank-1),lhs_s[3*idim],lhs_s[3*idim+1],lhs_s[3*idim+2])=
							rhs.ref_vector(idim)(Slice(0,rhsRank-1),rhs_s[3*idim],rhs_s[3*idim+1],rhs_s[3*idim+2]);

				} else {
					MADNESS_EXCEPTION("extend dim_pv in srconf::inplace_add",0);
				}
			}

			lhs.rank_=newRank;
			lhs.make_slices();
		}



		/// deep copy of rhs, shrink
		friend SRConf<T> copy(const SRConf<T>& rhs) {

			// if rhs is non-existent simply construct a new SRConf
			if (rhs.rank()==0) return SRConf<T>(rhs.dim(),rhs.get_k(),rhs.type());

			// pass a copy of the weights and vectors of rhs to ctor
			std::vector<tensorT> vector(rhs.dim_eff());
			for (unsigned int idim=0; idim<rhs.dim_eff(); idim++)
				vector[idim]=copy(rhs.ref_vector(idim)(rhs.c0()));

			return SRConf<T>(copy(rhs.weights_(Slice(0,rhs.rank()-1))),vector,rhs.dim(),rhs.get_k(),rhs.type());
		}

		/// reassign weight and data for a specific SRConf only
		void reassign(const unsigned int& idim, const unsigned int& r,
				const double& weight, const Tensor<T> & data, const unsigned int& maxk) {

			// some checks
	        MADNESS_ASSERT(idim<this->dim_eff());
	        MADNESS_ASSERT(r<this->rank());
	        MADNESS_ASSERT(this->is_flat());

	        // assign weight
	        weights_(r)=weight;

	        // assign data
	        MADNESS_ASSERT(data.size()==maxk);
	        for (unsigned int k=0; k<maxk; k++) {
	        	ref_vector(idim)(r,k)=data(k);
	        }

	        // make sure the ranks comply for all dimensions
	        for (unsigned int idim=0; idim<dim_eff(); idim++) {
                MADNESS_ASSERT(weights_.dim(0)==vector_[idim].dim(0));
	        }
		}

		/// redo the Slices for getting direct access to the configurations
		void make_slices() {
			if (this->rank()==0) {
				s_.clear();
			} else {
				// first dim is the rank
				s_.resize(vector_[0].ndim());
				s_[0]=Slice(0,this->rank()-1);
				for (int i=1; i<vector_[0].ndim(); i++) {
					s_[i] = Slice(_);
				}
			}
		}

		void make_structure() {

			// fast return if rank is zero
			if (this->rank()==0) return;

			const int dim_pv=this->dim_per_vector();
			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			const int rr=weights_.dim(0);	// not the rank!
			const int k=this->get_k();

			// reshape the vectors and adapt the Slices
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				if (dim_pv==2) vector_[idim]=vector_[idim].reshape(rr,k,k);
				if (dim_pv==3) vector_[idim]=vector_[idim].reshape(rr,k,k,k);
			}

			make_slices();

		}

		void undo_structure() {

			// fast return if rank is zero
			if (this->rank()==0) return;

			const int dim_pv=this->dim_per_vector();
			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			const int rr=weights_.dim(0);	// not the rank!
			const int kvec=this->kVec();

			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				vector_[idim]=vector_[idim].reshape(rr,kvec);
			}

			make_slices();
		}

		/// return reference to one of the vectors F
		Tensor<T>& ref_vector(const unsigned int& idim) {
			return vector_[idim];
		}

		/// return reference to one of the vectors F
		const Tensor<T>& ref_vector(const unsigned int& idim) const {
			return vector_[idim];
		}

		/// zero out
		void zeroOut() {MADNESS_ASSERT(0); this->chopoff(0);};

		/// fill this SRConf with 1 flattened random configurations (tested)
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
	        MADNESS_ASSERT(this->is_flat());
	        // we calculate the norm sum_i < F^r_i | F^r_i > for each dimension for each r

	        // loop over all configurations
	        for (unsigned int r=0; r<rank; r++) {

	        	// loop over all dimensions
	        	for (unsigned int idim=0; idim<dim_eff(); idim++) {

	        		Tensor<T> config=this->ref_vector(idim)(Slice(r,r),Slice(_));

	        		const double norm=config.normf();
	        		const double fac=norm;
	        		double oofac=1.0/fac;
	        		if (fac==0.0) oofac=0.0;

	        		weights_(r)*=fac;
	        		config.scale(oofac);
	        	}
	        }
		}

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

		/// return the number of coefficients
		unsigned int nCoeff() const {
//			return this->dim()*maxk_*rank_;
			return this->dim_eff()*this->kVec()*this->rank();
		};

		/// calculate the Frobenius inner product (tested)
		friend T overlap(const SRConf<T>& rhs, const SRConf<T>& lhs) {

			/*
			 * the structure of an SRConf is (r,k) or (r,k',k), with
			 * r the slowest index; the overlap is therefore simply calculated
			 * as the matrix multiplication rhs*lhs^T
			 */

			// some checks
			MADNESS_ASSERT(rhs.is_flat()==lhs.is_flat());
			MADNESS_ASSERT(rhs.dim()==lhs.dim());
			MADNESS_ASSERT(rhs.dim()>0);

			// fast return if either rank is 0
			if ((lhs.rank()==0) or (rhs.rank()==0)) return 0.0;

			const unsigned int dim_eff=rhs.dim_eff();

			// get the weight matrix
			Tensor<T> weightMatrix=outer(lhs.weights_,rhs.weights_);
			SRConf<T> lhs3(lhs);
			lhs3.undo_structure();
			SRConf<T> rhs3(rhs);
			rhs3.undo_structure();

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<dim_eff; idim++) {
				const Tensor<T> lhs2=lhs3.ref_vector(idim)(lhs.c0());
				const Tensor<T> rhs2=rhs3.ref_vector(idim)(rhs.c0());
				Tensor< TENSOR_RESULT_TYPE(T,T) > ovlp(lhs.rank(),rhs.rank());
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

		/// check compatibility
		friend bool compatible(const SRConf& lhs, const SRConf& rhs) {
			return ((lhs.dim()==rhs.dim()) and (lhs.dim_per_vector()==rhs.dim_per_vector()));
		}

		/// make one of the terms in the B matrix (BM2005)
		void friend makeB(Tensor<T> & B, const unsigned int& idim, const SRConf<T>& lhs, const SRConf<T>& rhs) {
			// some checks
			MADNESS_ASSERT(compatible(rhs,lhs));
			MADNESS_ASSERT(lhs.rank()==B.dim(0));
			MADNESS_ASSERT(rhs.rank()==B.dim(1));
			MADNESS_ASSERT(idim<rhs.dim_eff());

//			dgemm_NT(lhs.refVector(idim),rhs.refVector(idim),B,lhs.rank(),rhs.rank());
			Tensor<T> lhs2=lhs.ref_vector(idim)(lhs.c0());
			Tensor<T> rhs2=rhs.ref_vector(idim)(rhs.c0());
			B=0.0;
			inner_result(lhs2,rhs2,-1,-1,B);
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		SRConf<T> transform(const Tensor<T>& c) const {

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);
			if (this->rank()==0) return result;

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			// these two loops go over all physical dimensions (dim = dim_eff * merged_dim)
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				for (unsigned int jdim=1; jdim<this->ref_vector(idim).ndim(); jdim++) {

					// note: tricky ordering (jdim is missing): this is actually correct!
					// note: no slicing necessary, since we've copied this to result (incl shrinking)
//					result.refVector_struct(idim)=madness::inner(result.refVector_struct(idim),c,1,0);
					result.ref_vector(idim)=madness::inner(result.ref_vector(idim),c,1,0);

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

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);
			if (this->rank()==0) return result;

			// make sure this is not flattened
			if (not this->has_structure()) {
				print("no structure!");
			}
			MADNESS_ASSERT(this->has_structure());

			long i=0;
			// these two loops go over all physical dimensions (dim = dim_eff * merged_dim)
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				for (unsigned int jdim=1; jdim<this->ref_vector(idim).ndim(); jdim++) {

					// note tricky ordering (jdim is missing): this is actually correct!
					// note: no slicing necessary, since we've copied this to result (incl shrinking)
					result.ref_vector(idim)=madness::inner(result.ref_vector(idim),c[i],1,0);
					i++;

				}
			}
			return result;
		}


		SRConf<T> transform_dir(const Tensor<T>& c, const int& axis) const {

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);
			if (this->rank()==0) return result;

			// only a matrix is allowed for c
			MADNESS_ASSERT(c.ndim()==2);

			// make sure this is not flattened
			MADNESS_ASSERT(this->has_structure());

			// compute idim for accessing the vector_, and the dimension inside vector_
			// the +1 on jdim for the rank
			const long idim=axis/this->dim_per_vector();
			const long jdim=axis%this->dim_per_vector()+1;

			// note: no slicing necessary, since we've copied this to result (incl shrinking)
			result.ref_vector(idim)=madness::transform_dir(this->ref_vector(idim),c,jdim);

			return result;
		}

	};

}

#endif /* SRCONF_H_ */
