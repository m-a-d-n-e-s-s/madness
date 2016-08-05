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

namespace madness {



#ifdef BENCH
	static double time_[30];
#endif


	template <class T> class GenTensor;
    template <class T> class LowRankTensor;
    template <class T> class SliceLowRankTensor;

	/**
	 * A SRConf handles all the configurations in a Separated Representation.
	 */

	template <typename T>
	class SRConf {
		friend class GenTensor<T>;
		friend class LowRankTensor<T>;
        friend class SliceLowRankTensor<T>;

		/// return the number of vectors (i.e. dim_eff) according to the TensorType
		static unsigned int compute_nvec(const TensorType& tt) {
			if (tt==TT_FULL) return 1;
			if (tt==TT_2D) return 2;
			print("unknown TensorType",tt);
			MADNESS_ASSERT(0);
		}

		/// the scalar type of T
		typedef typename Tensor<T>::scalar_type scalar_type;
	public:

#ifdef BENCH
		static double& time(int i) {return time_[i];}
#endif

		typedef Tensor<T> tensorT;

		/// check orthonormality at low rank additions
		static const bool check_orthonormality=false;

		/// the number of dimensions (the order of the tensor)
		unsigned int dim_;

		/// for each configuration the weight; length should be r
		Tensor< typename Tensor<T>::scalar_type >  weights_;

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

		/// Slice containing the actual data in each vector, ignoring "empty" configurations;
		/// will maintain contiguity of the data.
		std::vector<Slice> s_;

		/// how will this be represented
		TensorType tensortype_;

		
	public:

    	/// return the index of the last singular vector/value to meet the threshold
		///        (returns -1 if all meet threshold, i.e. || A ||_2 < threshold)
    	/// given a matrix A in SVD form, truncate the singular values such that the
    	/// accuracy threshold is still met.
    	/// @param[in]	thresh	the threshold eps: || A - A(truncated) || < eps
    	/// @param[in] 	rank	the number of singular values in w
    	/// @param[in]	w		the weights/singular values of A
    	/// @return		i		the index of s_max to contribute: w(Slice(0,i)); i.e. inclusive!
    	static int max_sigma(const double& thresh, const int& rank, const Tensor<double>& w) {

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
		SRConf() : dim_(0), rank_(0), maxk_(0), s_(), tensortype_(TT_NONE) {
		};

		/// ctor with dimensions for a vector configuration (tested)
		SRConf(const unsigned int& dim, const unsigned int& k, const TensorType& tt)
			: dim_(dim)
			, rank_(0)
			, maxk_(k)
			, s_()
			, tensortype_(tt) {

			// make sure dim is integer multiple of requested TT
			const long nvec=compute_nvec(tt);
			MADNESS_ASSERT(dim%nvec==0);

			// construct empty vector
			weights_=Tensor<double>(int(0));
			vector_=std::vector<Tensor<T> > (nvec);

			if (tt==TT_FULL) {
				vector_[0]=tensorT(std::vector<long>(dim,k));
				rank_=-1;
			} else {
				for (unsigned int idim=0; idim<nvec; idim++) vector_[idim]=Tensor<T>(0,kVec());
			}
			make_structure();
			MADNESS_ASSERT(has_structure());
		}

		/// copy ctor (tested); shallow copy
		SRConf(const SRConf& rhs)  {
			*this=rhs;
            MADNESS_ASSERT(has_structure());
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
			vector_=std::vector<Tensor<T> > (long(vectors.size()));
			for (unsigned int idim=0; idim<vectors.size(); idim++) {
				vector_[idim]=vectors[idim];
			}
			make_slices();
            MADNESS_ASSERT(has_structure());
		}

		/// explicit ctor with one vector (aka full representation), shallow
		SRConf(const tensorT& vector1)
			: dim_(vector1.ndim())
			, weights_(Tensor<double>())
			, rank_(-1)
			, maxk_(vector1.dim(0))
			, tensortype_(TT_FULL) {

			vector_.resize(1);
			vector_[0]=vector1;
			MADNESS_ASSERT(has_structure());
		}

		/// explicit ctor with two vectors (aka SVD), shallow
		SRConf(const Tensor<double>& weights, const tensorT& vector1, const tensorT& vector2,
				const unsigned int& dim, const unsigned int maxk)
			: dim_(dim)
			, maxk_(maxk)
			, tensortype_(TT_2D) {

			MADNESS_ASSERT(weights.ndim()==1);
			MADNESS_ASSERT(vector1.ndim()==2);
			MADNESS_ASSERT(vector2.ndim()==2);
			MADNESS_ASSERT(weights.dim(0)==vector1.dim(0));
			MADNESS_ASSERT(vector2.dim(0)==vector1.dim(0));
			vector_.resize(2);
			vector_[0]=vector1;
			vector_[1]=vector2;
			weights_=weights;
			rank_=weights.dim(0);
			make_structure();
			make_slices();
            MADNESS_ASSERT(has_structure());
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

			if (rhs.has_no_data()) {
				// construct empty vector
				weights_=Tensor<double>(0);
				vector_=std::vector<Tensor<T> > (rhs.dim_eff());
				rank_=0;
				for (unsigned int idim=0; idim<dim_eff(); idim++) vector_[idim]=Tensor<T>(0,long(this->kVec()));
				make_structure();


			} else if (rhs.type()==TT_FULL) {
				weights_=Tensor<double>();
				rank_=-1;
				vector_.resize(1);
				vector_[0]=rhs.ref_vector(0);

			} else {
				// assign vectors; shallow copy
				vector_.resize(rhs.vector_.size());
				for (unsigned int i=0; i<rhs.vector_.size(); i++) {
					vector_[i]=rhs.vector_[i];
				}

				// shallow copy
				weights_=(rhs.weights_);
				rank_=rhs.rank();

				// consistency check
				for (unsigned int idim=0; idim<dim_eff(); idim++) {
					MADNESS_ASSERT(weights_.dim(0)==vector_[idim].dim(0));
				}
			}
            MADNESS_ASSERT(has_structure());
			return *this;
		}

		/// assign a number to this;
		SRConf& operator=(const T& number) {

		    if (type()==TT_2D) {
                // rank will be one
                rank_=1;
                vector_[0]=Tensor<T>(1,kVec());
                vector_[1]=Tensor<T>(1,kVec());
                vector_[0]=number;
                vector_[1]=1.0;
                weights_=Tensor< typename Tensor<T>::scalar_type > (1);
                weights_(0l)=1.0;
                make_structure();
                normalize();
		    } else if (type()==TT_FULL) {
		        vector_[0]=number;

		    }
		    return *this;
		}

		/// return some of the terms of the SRConf (start,..,end), inclusively
		/// shallow copy
		const SRConf get_configs(const int& start, const int& end) const {

			MADNESS_ASSERT((start>=0) and (end<=rank()));
			MADNESS_ASSERT(s_.size()>1);
			const long nvec=dim_eff();
			const long dim_pv_eff=s_.size()-1;	// #dim per vector minus rank-dim

			Slice s(start,end);
			std::vector<tensorT> v(nvec);

			// slice vectors
			if (dim_pv_eff==1) {
				for (long i=0; i<nvec; i++) v[i]=ref_vector(i)(s,_);
			} else if (dim_pv_eff==2) {
				for (long i=0; i<nvec; i++) v[i]=ref_vector(i)(s,_,_);
			} else if (dim_pv_eff==3) {
				for (long i=0; i<nvec; i++) v[i]=ref_vector(i)(s,_,_,_);
			} else {
				MADNESS_EXCEPTION("faulty dim_pv in SRConf::get_configs",0);
			}

			SRConf<T> result(weights_(s),v,dim(),get_k(),type());
            MADNESS_ASSERT(result.has_structure());
			return result;
		}

		/// dtor
		~SRConf() {}

        template <typename Archive>
        void serialize(Archive& ar) {
              	int i=int(tensortype_);
              	ar & dim_ & weights_ & vector_ & rank_ & maxk_ & i;
              	tensortype_=TensorType(i);
              	make_slices();
                MADNESS_ASSERT(has_structure());
        }

		/// return the tensor type
		TensorType type() const {return tensortype_;};

		/// does this have any data?
		bool has_data() const {
			if (tensortype_==TT_FULL) return (vector_.size()>0 and vector_[0].has_data());
			return rank()>0;
		}

	private:

		/// does this have any data?
		bool has_no_data() const {return !has_data();}

		/// reserve enough space to hold at least r configurations
		void reserve(long r) {

			// this should at least hold the current information
			MADNESS_ASSERT(r>=this->rank());
			MADNESS_ASSERT(has_data() or vector_.size()>0);

			// fast return if possible
			// nothing to be done
			if (r==0) return;
			// already large enuff?
			if (this->vector_[0].dim(0)>=r) return;

			// to avoid incremental increase of the rank
			r+=3;

			// for convenience
			const long rank=this->rank();
			const long kvec=this->kVec();
			const bool had_structure=this->has_structure();
			if (had_structure) this->undo_structure();

			// transfer weights
			Tensor<scalar_type> newWeights(r);
			if (rank>0) newWeights(Slice(0,rank-1))=weights_(Slice(0,rank-1));
			std::swap(weights_,newWeights);

			// transfer vectors
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {

				tensorT newVector(r,kvec);
				if (rank>0) newVector(this->c0())=vector_[idim](this->c0());
				std::swap(vector_[idim],newVector);

			}
			MADNESS_ASSERT(weights_.dim(0)==vector_[0].dim(0));
			if (had_structure) this->make_structure(true);
            MADNESS_ASSERT(has_structure());

		}

		/// return a Slice that corresponds the that part of vector_ that holds coefficients
		const std::vector<Slice>& c0() const {
			MADNESS_ASSERT(s_.size()>0);
			return s_;
		}

		/// reduce the rank using a divide-and-conquer approach
		void divide_and_conquer_reduce(const double& thresh) {

			if (has_no_data()) return;
			if (rank()==1) {
				normalize();
				return;
			}

			// divide the SRConf into two
			const long chunksize=8;
			if (rank()>chunksize) {
        		SRConf<T> chunk1=this->get_configs(0,rank()/2);
        		SRConf<T> chunk2=this->get_configs(rank()/2+1,rank()-1);
        		chunk1.divide_and_conquer_reduce(thresh*0.5);
        		chunk2.divide_and_conquer_reduce(thresh*0.5);

        		// collect the two SRConfs
        		*this=chunk1;
        		this->add_SVD(chunk2,thresh);

			} else {

				// and reduce the rank
				this->orthonormalize(thresh);
			}
            MADNESS_ASSERT(has_structure());
		}

	public:
		/// orthonormalize this
		void orthonormalize(const double& thresh) {

			if (type()==TT_FULL) return;
			if (has_no_data()) return;
			if (rank()==1) {
				normalize();
				return;
			}
#ifdef BENCH
			double cpu0=wall_time();
#endif
//			vector_[0]=copy(vector_[0](c0()));
//			vector_[1]=copy(vector_[1](c0()));
//			weights_=weights_(Slice(0,rank()-1));
            normalize();
#ifdef BENCH
			double cpu1=wall_time();
#endif
//            this->undo_structure();
			weights_=weights_(Slice(0,rank()-1));
//			ortho3(vector_[0],vector_[1],weights_,thresh);
            tensorT v0=flat_vector(0);
            tensorT v1=flat_vector(1);
#ifdef BENCH
			double cpu2=wall_time();
#endif
			ortho3(v0,v1,weights_,thresh);
            std::swap(vector_[0],v0);
            std::swap(vector_[1],v1);
#ifdef BENCH
			double cpu3=wall_time();
#endif
			rank_=weights_.size();
			MADNESS_ASSERT(rank_>=0);
			this->make_structure();
			make_slices();
            MADNESS_ASSERT(has_structure());
#ifdef BENCH
			double cpu4=wall_time();
			SRConf<T>::time(21)+=cpu1-cpu0;
			SRConf<T>::time(22)+=cpu2-cpu1;
			SRConf<T>::time(23)+=cpu3-cpu2;
			SRConf<T>::time(24)+=cpu4-cpu3;
			SRConf<T>::time(20)+=cpu4-cpu0;
#endif

		}

	private:
		/// append configurations of rhs to this

		/// simplified version of inplace_add for flattened configurations
		/// *this += fac*rhs
		void append(const SRConf<T>& rhs, const double fac=1.0) {

			// fast return if possible
			if (rhs.has_no_data()) return;
			if (this->has_no_data()) {
				*this=copy(rhs);
				this->scale(fac);
				return;
			}

			const long newRank=this->rank()+rhs.rank();
			const long lhsRank=this->rank();
			const long rhsRank=rhs.rank();
			reserve(newRank);

			// assign weights
			this->weights_(Slice(lhsRank,newRank-1))=rhs.weights_(Slice(0,rhsRank-1))*fac;
			std::vector<Slice> s(dim_per_vector()+1,_);
			s[0]=Slice(lhsRank,newRank-1);

			// assign vectors
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
//              vector_[idim](Slice(lhsRank,newRank-1),_)=rhs.vector_[idim](rhs.c0());
				vector_[idim](s)=rhs.vector_[idim](rhs.c0());
			}

			rank_=newRank;
			make_slices();
            MADNESS_ASSERT(has_structure());

		}
		void append(const SRConf<T>& rhs, const double_complex fac=1.0) {
			MADNESS_EXCEPTION("no complex in SRConf",1);
		}

		/// add two orthonormal configurations, yielding an optimal SVD decomposition
		void add_SVD(const SRConf<T>& rhs, const double& thresh) {
#ifdef BENCH
			double cpu0=wall_time();
#endif
			if (rhs.has_no_data()) return;
			if (has_no_data()) {
				*this=rhs;
				return;
			}

			if (check_orthonormality) check_right_orthonormality();
            if (check_orthonormality) rhs.check_right_orthonormality();

            this->undo_structure();
            ortho5(ref_vector(0),ref_vector(1),weights_,
					rhs.flat_vector(0),rhs.flat_vector(1),rhs.weights_,thresh);
			rank_=weights_.size();
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
		void inplace_add(const SRConf<T>& rhs2, std::vector<Slice> lhs_s,
				std::vector<Slice> rhs_s, const double alpha, const double beta) {

			// fast return if possible; no fast return for this.rank()==0
			// since we might work with slices!
			if (rhs2.has_no_data()) return;

			// fast return for full rank tensors
			if (type()==TT_FULL) {
				vector_[0](lhs_s)+=rhs2.vector_[0](rhs_s);
				return;
			}

			// unflatten this and rhs; shallow wrt vector_
			SRConf<T>& lhs=*this;
			const SRConf<T>& rhs=rhs2;
			if (lhs.has_no_data()) lhs.make_structure(true);
			MADNESS_ASSERT(lhs.has_structure() or (lhs.has_no_data()));
			MADNESS_ASSERT(rhs.has_structure());

			// conflicts with lhs_s ??
			MADNESS_ASSERT(alpha==1.0);

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
            MADNESS_ASSERT(has_structure());
		}

		/// deep copy of rhs, shrink
		friend SRConf<T> copy(const SRConf<T>& rhs) {

			// if rhs is non-existent simply construct a new SRConf
			if (rhs.has_no_data()) return SRConf<T>(rhs.dim(),rhs.get_k(),rhs.type());

			if (rhs.type()==TT_FULL) return SRConf<T>(copy(rhs.ref_vector(0)));

			// pass a copy of the weights and vectors of rhs to ctor
			std::vector<tensorT> vector(rhs.dim_eff());
			for (unsigned int idim=0; idim<rhs.dim_eff(); idim++)
				vector[idim]=copy(rhs.ref_vector(idim)(rhs.c0()));

			return SRConf<T>(copy(rhs.weights_(Slice(0,rhs.rank()-1))),vector,rhs.dim(),rhs.get_k(),rhs.type());
		}

public:
        /// return a slice of this (deep copy)
        SRConf<T> copy_slice(const std::vector<Slice>& s) const {

            // fast return if possible
            if (this->has_no_data()) {
                int k_new=s[0].end-s[0].start+1;
                return SRConf<T>(dim(),k_new,this->type());
            }

            // consistency check
            MADNESS_ASSERT(s.size()==this->dim());
            MADNESS_ASSERT(s[0].step==1);

            // fast return for full rank tensors
            if (type()==TT_FULL) {
                tensorT a=copy(ref_vector(0)(s));
                return SRConf<T>(a);
            }


            MADNESS_ASSERT(has_structure());
//          _ptr->make_structure();

            // get dimensions
            const TensorType tt=this->type();
            const int merged_dim=this->dim_per_vector();
            const int dim_eff=this->dim_eff();
            const int rank=this->rank();
            int k_new=s[0].end-s[0].start+1;
            if (s[0].end<0) k_new+=this->get_k();

            // get and reshape the vectors, slice and re-reshape again;
            // this is shallow
            const SRConf<T>& sr=*this;

            std::vector<Tensor<T> > vectors(dim_eff,Tensor<T>());

            for (int idim=0; idim<dim_eff; idim++) {

                // assignment from/to slice is deep-copy
                if (merged_dim==1) {
                    if (rank>0) {
                        vectors[idim]=copy(sr.ref_vector(idim)(Slice(0,rank-1),s[idim]));
                    } else {
                        vectors[idim]=Tensor<T>(0,s[idim].end-s[idim].start+1);
                    }
                } else if (merged_dim==2) {
                    if (rank>0) {
                        vectors[idim]=copy(sr.ref_vector(idim)(Slice(0,rank-1),s[2*idim],s[2*idim+1]));
                    } else {
                        vectors[idim]=tensorT(0,s[2*idim].end-s[2*idim].start+1,
                                                s[2*idim+1].end-s[2*idim+1].start+1);
                    }
                } else if (merged_dim==3) {
                    if (rank>0) {
                        vectors[idim]=copy(sr.ref_vector(idim)(Slice(0,rank-1),
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
                weights=copy(this->weights_(Slice(0,rank-1)));
            } else {
                weights=Tensor<double>(int(0));
            }
            return SRConf<T>(weights,vectors,this->dim(),k_new,tt);
        }

        /// perform elementwise Hadamard product
        SRConf<T>& emul(const SRConf<T>& other) {
            // consistency check
            MADNESS_ASSERT(this->dim()==other.dim());
            MADNESS_ASSERT(this->get_k()==other.get_k());

            long finalrank=this->rank()*other.rank();
            SRConf<T> result(dim(),get_k(),TT_2D);
            if ((this->rank()==0) or (other.rank()==0)) {
                ;   // pass
            } else {

                result.vector_[0]=Tensor<T>(finalrank,kVec());
                result.vector_[1]=Tensor<T>(finalrank,kVec());
                result.weights_=outer(weights_,other.weights_).flat();
                result.rank_=finalrank;

                for (int k=0; k<kVec(); ++k) {
                    Tensor<T> a1=flat_vector(0)(_,Slice(k,k));   // (1,k)->(k)
                    Tensor<T> a2=flat_vector(1)(_,Slice(k,k));
                    Tensor<T> b1=other.flat_vector(0)(_,Slice(k,k));
                    Tensor<T> b2=other.flat_vector(1)(_,Slice(k,k));

                    result.vector_[0](_,Slice(k,k))=outer(a1,b1).reshape(finalrank,1);
                    result.vector_[1](_,Slice(k,k))=outer(a2,b2).reshape(finalrank,1);
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
			if (type()==TT_FULL) return;
			if (this->has_no_data()) {
				s_.clear();
			} else {
				// first dim is the rank
				if (vector_[0].ndim()>TENSOR_MAXDIM) {
					print(*this);
					MADNESS_EXCEPTION("serializing failed",0);
				}
				s_.resize(vector_[0].ndim());
				s_[0]=Slice(0,this->rank()-1);
				for (int i=1; i<vector_[0].ndim(); i++) {
					s_[i] = Slice(_);
				}
			}
		}


		void make_structure(bool force=false) {

			// fast return if rank is zero
			if ((not force) and this->has_no_data()) return;
			if (type()==TT_FULL) return;

			const int dim_pv=this->dim_per_vector();
			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			int rr=weights_.dim(0);	// not the rank!
			if (weights_.size()==0) rr=0;
			const int k=this->get_k();

			// reshape the vectors and adapt the Slices
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				if (dim_pv==2) this->vector_[idim]=vector_[idim].reshape(rr,k,k);
				if (dim_pv==3) this->vector_[idim]=vector_[idim].reshape(rr,k,k,k);
			}

			this->make_slices();

		}

		void undo_structure(bool force=false) {

			// fast return if rank is zero
			if ((not force) and this->has_no_data()) return;
			if (type()==TT_FULL) return;

			const int dim_pv=this->dim_per_vector();
			MADNESS_ASSERT(dim_pv>0 and dim_pv<=3);
			int rr=weights_.dim(0);	// not the rank!
			if (weights_.size()==0) rr=0;
			const int kvec=this->kVec();

			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				this->vector_[idim]=this->vector_[idim].reshape(rr,kvec);
			}

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

	private:
		/// return shallow copy of a slice of one of the vectors, flattened to (r,kVec)
		const Tensor<T> flat_vector(const unsigned int& idim) const {
		    MADNESS_ASSERT(rank()>0);
		    return vector_[idim](c0()).reshape(rank(),kVec());
		}

		/// return shallow copy of a slice of one of the vectors, flattened to (r,kVec)
		Tensor<T> flat_vector(const unsigned int& idim) {
		    MADNESS_ASSERT(rank()>0);
		    return vector_[idim](c0()).reshape(rank(),kVec());
		}

		/// fill this SRConf with 1 flattened random configurations (tested)
		void fillWithRandom(const long& rank=1) {

			rank_=rank;

			// assign; note that Slice(0,_) is inclusive
			weights_=Tensor<double>(rank);
			weights_=1.0;

			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				vector_[idim]=Tensor<T>(rank_,this->kVec());
				vector_[idim].fillrandom();
			}

			this->normalize();
			for (unsigned int idim=0; idim<this->dim_eff(); idim++) {
				vector_[idim].scale(madness::RandomValue<T>()*scalar_type(10.0));
			}
			weights_(Slice(0,this->rank()-1)).fillrandom().scale(10.0);
			make_slices();
            MADNESS_ASSERT(has_structure());
		}

		/// normalize the vectors (tested)
		void normalize() {

			if (type()==TT_FULL) return;
			if (rank()==0) return;
            MADNESS_ASSERT(has_structure());

	        // for convenience
	        const unsigned int rank=this->rank();
	        std::vector<Slice> s(dim_per_vector()+1,_);
//	        std::vector<Slice> s(2,_);

	        // we calculate the norm sum_i < F^r_i | F^r_i > for each dimension for each r

	        // loop over all configurations
	        for (unsigned int r=0; r<rank; r++) {
	            s[0]=Slice(r,r);
	        	// loop over all dimensions
	        	for (unsigned int idim=0; idim<dim_eff(); idim++) {

//	        		Tensor<T> config=this->ref_vector(idim)(s);
//	        		const double norm=config.normf();
//	        		const double fac=norm;
//	        		double oofac=1.0/fac;
//	        		if (fac<1.e-13) oofac=0.0;
//	        		weights_(r)*=fac;
//	        		config.scale(oofac);

//	        		const double norm=this->ref_vector(idim)(s).normf();
	        		const double norm=this->vector_[idim](s).normf();
	        		const double fac=norm;
	        		double oofac=1.0/fac;
	        		if (fac<1.e-13) oofac=0.0;
	        		weights_(r)*=fac;
//	        		this->ref_vector(idim)(s).scale(oofac);
//	        		this->flat_vector(idim)(s).scale(oofac);
	        		vector_[idim](s).scale(oofac);

	        	}
	        }
            MADNESS_ASSERT(has_structure());
		}

		/// check if the terms are orthogonal
		bool check_right_orthonormality() const {

			// fast return if possible
			if (rank()==0) return true;

			MADNESS_ASSERT(type()==TT_2D);

			const tensorT t1=ref_vector(1)(c0()).reshape(rank(),kVec());
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
            return (type()==TT_FULL or has_no_data() or vector_[0].dim(1)==this->get_k());
		}

	private:
		/// return the dimension of this
		unsigned int dim() const {return dim_;}

		/// return the number of vectors
		unsigned int dim_eff() const {return vector_.size();}

		/// return the logicalrank
		long rank() const {return rank_;};

		/// return the number of physical matrix elements per dimension
		unsigned int get_k() const {return maxk_;};

		/// return the length of the vector (dim_pv*maxk)
		long kVec() const {
			const int dimpv=this->dim_per_vector();
			int kv=1;
			for (int i=0; i<dimpv; ++i) kv*=this->get_k();
//			const int kv1= pow(this->get_k(),this->dim_per_vector());
//			MADNESS_ASSERT(kv==kv1);
			return kv;
		}

	public:
		/// return the number of physical dimensions
		int dim_per_vector() const {
			const int nvec=vector_.size();
			const int dim=this->dim();
			MADNESS_ASSERT(dim%nvec==0);
			return dim/nvec;
		}

		/// return the weight
		double weights(const unsigned int& i) const {return weights_(i);};

        /// reconstruct this to return a full tensor
        Tensor<T> reconstruct() const {

            if (type()==TT_FULL) return ref_vector(0);

            /*
             * reconstruct the tensor first to the configurational dimension,
             * then to the real dimension
             */

            // for convenience
            const unsigned int conf_dim=this->dim_eff();
            const unsigned int conf_k=this->kVec();           // possibly k,k*k,..
            const long rank=this->rank();
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
            SRConf<T> sr=*this;

            // and a scratch Tensor
            Tensor<T>  scr(rank);
            Tensor<T>  scr1(rank);
            Tensor<T>  scr2(rank);

            if (conf_dim==1) {

                for (unsigned int i0=0; i0<conf_k; i0++) {
                    scr=sr.weights_;
                    //                  scr.emul(F[0][i0]);
                    T buffer=scr.sum();
                    s(i0)=buffer;
                }

            } else if (conf_dim==2) {


                //              tensorT weight_matrix(rank,rank);
                //              for (unsigned int r=0; r<rank; r++) {
                //                  weight_matrix(r,r)=this->weight(r);
                //              }
                //              s=inner(weight_matrix,sr._ptr->refVector(0));
                //              s=inner(s,sr._ptr->refVector(1),0,0);
//                tensorT sscr=copy(sr._ptr->ref_vector(0)(sr._ptr->c0()));
                tensorT sscr=copy(sr.flat_vector(0));
                for (unsigned int r=0; r<rank; r++) {
                    const double w=weights(r);
                    for (unsigned int k=0; k<conf_k; k++) {
                        sscr(r,k)*=w;
                    }
                }
                inner_result(sscr,sr.flat_vector(1),0,0,s);


            } else {
                print("only config_dim=1,2 in SRConf::reconstruct");
                MADNESS_ASSERT(0);
            }


            // reshape the tensor to the really required one
            const long k=this->get_k();
            const long dim=this->dim();
            for (long i=0; i<dim; i++) d[i] = k;

            Tensor<T> s2=s.reshape(dim,d);
            return s2;
        }


	private:
		/// return the number of coefficients
		unsigned int nCoeff() const {
			if (type()==TT_FULL) return ref_vector(0).size();
			return this->dim_eff()*this->kVec()*this->rank();
		};

		/// return the real size of this
		size_t real_size() const {
			size_t n=0;
			for (size_t i=0; i<vector_.size(); ++i) {
				n+=vector_[i].size()*sizeof(T) + sizeof(tensorT);
			}
			n+=weights_.size()*sizeof(double) + sizeof(Tensor<double>);
			n+=sizeof(*this);
			n+=s_.size()*sizeof(Slice);
			return n;
		}

		/// calculate the Frobenius inner product (tested)
		template<typename Q>
		friend TENSOR_RESULT_TYPE(T,Q) overlap(const SRConf<T>& rhs, const SRConf<Q>& lhs) {

			// fast return if either rank is 0
			if ((lhs.has_no_data()) or (rhs.has_no_data())) return 0.0;

			/*
			 * the structure of an SRConf is (r,k) or (r,k',k), with
			 * r the slowest index; the overlap is therefore simply calculated
			 * as the matrix multiplication rhs*lhs^T
			 */

			// some checks
			MADNESS_ASSERT(rhs.dim()==lhs.dim());
			MADNESS_ASSERT(rhs.dim()>0);

			typedef TENSOR_RESULT_TYPE(T,Q) resultT;

			if (rhs.type()==TT_FULL) {
				return rhs.ref_vector(0).trace(lhs.ref_vector(0));
			}

			const unsigned int dim_eff=rhs.dim_eff();

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
            if (has_no_data()) return 0.0;
            MADNESS_ASSERT(type()==TT_2D);
            return weights_(Slice(0,rank()-1)).normf();
        }

		/// calculate the Frobenius norm
		typename TensorTypeData<T>::float_scalar_type normf() const {

			// fast return if possible
			if (has_no_data()) return 0.0;
			if (type()==TT_FULL) return ref_vector(0).normf();

			// some checks
			MADNESS_ASSERT(dim()>0);
			MADNESS_ASSERT(not TensorTypeData<T>::iscomplex);

			// get the weight matrix
			Tensor<T> weightMatrix=outer(weights_(Slice(0,rank()-1)),weights_(Slice(0,rank()-1)));

			// calculate the overlap matrices for each dimension at a time
			for (unsigned int idim=0; idim<dim_eff(); idim++) {
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

		/// scale this by a number
		void scale(const double& fac) {weights_.scale(fac);};

		void scale(const double_complex& fac) {
			MADNESS_EXCEPTION("no complex scaling in SRConf",1);
		}

		/// check compatibility
		friend bool compatible(const SRConf& lhs, const SRConf& rhs) {
			return ((lhs.dim()==rhs.dim()) and (lhs.dim_per_vector()==rhs.dim_per_vector()));
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		SRConf<T> transform(const Tensor<T>& c) const {

			// fast return if possible
			if (this->has_no_data()) {
				return copy(*this);
			}

			// fast return for full rank tensor
			if (type()==TT_FULL) {
				return SRConf<T> (madness::transform(this->vector_[0],c));
			}

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);

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
            MADNESS_ASSERT(result.has_structure());
			return result;
		}

	    /// \code
		///     result(i,j,k,...) <-- sum(i',j', k',...) t(i',j',k',...)  c(i',i) c(j',j) c(k',k) ...
		/// \endcode
		///
		/// The input dimensions of \c t must all be the same .
		template<typename Q>
		SRConf<TENSOR_RESULT_TYPE(T,Q) > general_transform(const Tensor<Q> c[]) const {

			// fast return if possible
			if (this->has_no_data()) return SRConf<T>(copy(*this));
			if (type()==TT_FULL) {
				return SRConf<T> (madness::general_transform(this->vector_[0],c));
			}

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);

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
            MADNESS_ASSERT(result.has_structure());
			return result;
		}

		SRConf<T> transform_dir(const Tensor<T>& c, const int& axis) const {

			if (this->has_no_data()) {
				return SRConf<T>(copy(*this));
			}

			// fast return for full rank tensor
			if (type()==TT_FULL) {
				return SRConf<T> (madness::transform_dir(this->vector_[0],c,axis));
			}

			// copying shrinks the vectors to (r,k,k,..)
			SRConf<T> result=copy(*this);

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
	void ortho3(Tensor<T>& x, Tensor<T>& y, Tensor<double>& weights, const double& thresh) {

#ifdef BENCH
		double cpu0=wall_time();
#endif
		typedef Tensor<T> tensorT;

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
		Tensor<double> e1, e2;
	    syev(S1,U1,e1);
	    syev(S2,U2,e2);										// 2.3 / 4.0
#ifdef BENCH
		double cpu3=wall_time();
		SRConf<T>::time(3)+=cpu3-cpu1;
#endif

	    const double e1_max=e1.absmax();
	    const double e2_max=e2.absmax();

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
	    Tensor<double> sqrt_e1(rank), sqrt_e2(rank);


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
    		double fac=1.0/sqrt_e1(r);
    		for (unsigned int t=0; t<rank; t++) {
	    		U1(t,r)*=fac;
//	    		if (sqrt_e1(r)<thresh) throw;
    		}
    	}

	   	for (unsigned int r=0; r<rank2; r++) {
    		double fac=1.0/sqrt_e2(r);
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
		Tensor<double> Sp;
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
	void ortho5(Tensor<T>& x1, Tensor<T>& y1, Tensor<double>& w1,
				const Tensor<T>& x2, const Tensor<T>& y2, const Tensor<double>& w2,
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
		Tensor<double> e1, e2;
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
//    		for (unsigned int t=0; t<rank; t++) {
//	    		U1(t,r)*=fac;
////	    		if (sqrt_e1(r)<thresh) throw;
//    		}
    	}

	   	for (unsigned int r=0; r<rank_y; r++) {
    		double fac=1.0/sqrt_e2(r);
    		U2(_,r)*=fac;
//    		for (unsigned int t=0; t<rank; t++) {
//	    		U2(t,r)*=fac;
////	    		if (sqrt_e2(r)<thresh) throw;
//	    	}
	    }													// 0.2 / 0.6
#ifdef BENCH
		double cpu4=wall_time();
		SRConf<T>::time(14)+=cpu4-cpu3;
#endif


	    // decompose M
		tensorT Up,VTp;
		Tensor<double> Sp;
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

		s << "dim_          " << sr.dim_ << "\n";;
		s << "rank_         " << sr.rank_ << "\n";;
		s << "maxk_         " << sr.maxk_ << "\n";;
		s << "vector_.size()" << sr.vector_.size() << "\n";
		s << "has_data()    " << sr.has_data() << "\n";
		s << "TensorType    " << sr.type() << "\n\n";
		return s;
	}
}

#endif /* SRCONF_H_ */
