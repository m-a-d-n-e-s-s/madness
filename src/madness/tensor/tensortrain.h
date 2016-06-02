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

  $Id: srconf.h 2816 2012-03-23 14:59:52Z 3ru6ruWu $
*/

#ifndef TENSORTRAIN_H_
#define TENSORTRAIN_H_

#include <madness/tensor/tensor.h>
#include <madness/tensor/srconf.h>
#include <madness/tensor/clapack.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/fortran_ctypes.h>


/**
  \file tensortrain.h
  \brief Defines and implements the tensor train decomposition as described in
         I.V. Oseledets, Siam J. Sci. Comput. 33, 2295 (2011).
  \ingroup tensor
  \addtogroup tensor

*/


namespace madness {

	/**
	 * A tensor train is a multi-modal representation of a tensor t
	 * \code
	 *  t(i,j,k,l) = \sum G^1_{a1,i,a2} G^2_{a2,j,a3} G^3_{a3,k,a4} G^4_{a4,l,a5}
	 * \endcode
	 * The "core" tensors G are connected via a linear index network, where the
	 * first index a1 and the last index a5 are boundary indices and are set to 1.
	 *
	 * The tensor train representation is suited for any number of dimensions and
	 * in general at least as fast as the 2-way decomposition SVD. If the tensor
	 * has full rank it will need about twice the storage space of the full tensor
	 */
	template<typename T>
	class TensorTrain {

	    // make all types of TensorTrain friends of each other (for type conversion)
	    template<typename Q>
	    friend class TensorTrain;

        /// C++ typename of this tensor.
        typedef T type;

        /// C++ typename of the real type associated with a complex type.
        typedef typename TensorTypeData<T>::scalar_type scalar_type;

        /// C++ typename of the floating point type associated with scalar real type
        typedef typename TensorTypeData<T>::float_scalar_type float_scalar_type;

		/// holding the core tensors of a tensor train
		/// the tensors have the shape (k,r0) (r0,k,r1) (r1,k,r2) .. (rn-1,k)
		std::vector<Tensor<T> > core;
		/// true if rank is zero
		bool zero_rank;

	public:

		/// empty constructor
		TensorTrain() : core(), zero_rank(true) {
		}

		/// ctor for a TensorTrain, with the tolerance eps

		/// The tensor train will represent the input tensor with
		/// accuracy || t - this ||_2 < eps
		///
		/// Note that we rely on specific layout of the memory in the tensors, e.g.
		/// we pass SliceTensors on to lapack. This will only work if the slices are
		/// contiguous.
		///
		/// @param[in]	t	full representation of a tensor
		/// @param[in]	eps	the accuracy threshold
		TensorTrain(const Tensor<T>& t, double eps)
			: core(), zero_rank(false) {

		    MADNESS_ASSERT(t.size() != 0);
            MADNESS_ASSERT(t.ndim() != 0);

            std::vector<long> dims(t.ndim());
            for (int d=0; d<t.ndim(); ++d) dims[d]=t.dim(d);
            decompose(t.flat(),eps,dims);

		}

		/// ctor for a TensorTrain, with the tolerance eps

		/// The tensor train will represent the input tensor with
		/// accuracy || t - this ||_2 < eps
		///
		/// Note that we rely on specific layout of the memory in the tensors, e.g.
		/// we pass SliceTensors on to lapack. This will only work if the slices are
		/// contiguous.
		///
		/// @param[in]	t		full representation of a tensor
		/// @param[in]	eps		the accuracy threshold
		/// @param[in]	dims	the tt structure
		TensorTrain(const Tensor<T>& t, double eps, const std::vector<long> dims)
			: core(), zero_rank(false) {

		    MADNESS_ASSERT(t.size() != 0);
            MADNESS_ASSERT(t.ndim() != 0);
            decompose(t,eps,dims);
		}

		/// ctor for a TensorTrain, set up only the dimensions, no data
		TensorTrain(const std::vector<long>& dims) {
            zero_rank = true;

            core.resize(dims.size());
            // first and last core tensor
            core[0] = Tensor<T>(dims[0],long(0));
            core[dims.size()-1] = Tensor<T>(long(0),dims[dims.size()-1]);

            // iterate through the rest -- fast forward
            for (int d=1; d<dims.size()-1; ++d) {
                core[d] = Tensor<T>(long(0),dims[d],long(0));
		    }

		}

		/// copy constructor, shallow
		TensorTrain(const TensorTrain& other) : core(other.core),
		        zero_rank(other.zero_rank) {
		}

        /// Type conversion makes a deep copy
        template <class Q> operator TensorTrain<Q>() const { // type conv => deep copy

            TensorTrain<Q> result(this->dims());
            result.zero_rank=zero_rank;
            for (const Tensor<T>& c : core) {
                result.core.push_back(Tensor<Q>(c));
            }
            return result;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
                ar & core & zero_rank;
        }


		/// deep copy of the whole tensor

		/// if argument has zero rank return a zero-rank tensor of the same dimensions
		friend TensorTrain copy(const TensorTrain& other) {

		    // fast return
		    if (other.zero_rank) return TensorTrain(other.dims());

            TensorTrain result;
            for (const Tensor<T>& t: other.core) {
                if (t.size()==0) return TensorTrain(other.dims());  // also zero
                result.core.push_back(madness::copy(t));
            }
            result.zero_rank=other.zero_rank;
            return result;
		}

		/// deep copy of a slice of the tensor

		/// this operation does not change the ranks, i.e. the resulting
		/// tensor is most likely not in an optimal compression state
		/// @param[in]  other   tensor to be sliced
		/// @param[in]  s       vector of slices
		friend TensorTrain copy(const TensorTrain& other, const std::vector<Slice>& s) {

		    MADNESS_ASSERT(other.ndim()==s.size());
		    if (other.zero_rank) return TensorTrain(other.dims());

		    TensorTrain result;
		    const long nd=other.ndim();
		    result.zero_rank=other.zero_rank;
		    result.core.resize(nd);

		    // special treatment for first and last core tensor
		    // slice dim only, keep ranks
		    result.core[0]=copy(other.core[0](s[0],_));
		    for (long i=1; i<nd-1; ++i) {
		        result.core[i]=copy(other.core[i](_,s[i],_));
		    }

		    if (other.core.size()>1) result.core[nd-1]=copy(other.core[nd-1](_,s[nd-1]));
		    return result;
		}

		/// decompose the input tensor into a TT representation

		/// @param[in]	t		tensor in full rank
		/// @param[in]	eps		the precision threshold
		/// @param[in]	dims	the tt structure
		void decompose(const Tensor<T>& t, double eps,
				const std::vector<long>& dims) {

			core.resize(dims.size());
			eps=eps/sqrt(dims.size()-1);	// error is relative

			// the maximum rank is the smaller one of the ranks unfolded from the
			// left and the right.
			int rmax1=1;
			int rmax2=1;
			long n1=1, n2=1;
			long lwork=0;
			for (int i=0, j=dims.size()-1; i<=j; ++i, --j) {
				rmax1*=dims[i];
				rmax2*=dims[j];

				// work array for dgesvd
				n1*=dims[i];
				long m1=t.size()/dims[i];
				n2*=dims[j];
				long m2=t.size()/dims[j];
				long max1=std::max(n1,m1);
				long max2=std::max(n2,m2);
				long min1=std::min(n1,m1);
				long min2=std::min(n2,m2);

				lwork=std::max(lwork,std::max(3*min1+max1,5*min1));
				lwork=std::max(lwork,std::max(3*min2+max2,5*min2));
			}
			const int rmax=std::min(rmax1,rmax2);

			// these are max dimensions, so we can avoid frequent reallocation
			Tensor<T> u(rmax1*rmax2);
			Tensor<T> dummy;
			Tensor< typename Tensor<T>::scalar_type > s(rmax);

			// the dimension of the remainder tensor; will be cut down in each iteration
			long vtdim=t.size();

			// c will be destroyed, and assignment is only shallow, so need to deep copy
			Tensor<T> c=madness::copy(t);

			// work array for dgesvd
			Tensor<T> work(lwork);


			// this keeps track of the ranks
			std::vector<long> r(dims.size()+1,0l);
			r[0] = r[dims.size()] = 1;

			for (std::size_t d=1; d<dims.size(); ++d) {

				// the core tensors will have dimensions c(rank_left,k,rank_right)
				// or for lapack's purposes c(d1,rank_right)
				const long k=dims[d-1];
				const long d1=r[d-1]*k;
				c=c.reshape(d1,c.size()/d1);
				const long rmax=std::min(c.dim(0),c.dim(1));
				vtdim=vtdim/k;

				// testing
#if 0
				// c will be destroyed upon return
				Tensor<T> aa=copy(c);
#endif
				// The svd routine assumes lda=a etc. Pass in a flat tensor and reshape
				// and slice it after processing.
				u=u.flat();
				svd_result(c,u,s,dummy,work);

				// this is rank_right
				r[d]=SRConf<T>::max_sigma(eps,rmax,s)+1;
				const long rank=r[d];

				// this is for testing
#if 0
				if (0) {
					Tensor<T> uu=(u(Slice(0,c.dim(0)*rmax-1))).reshape(c.dim(0),rmax);
					Tensor<T> vtt=copy(c(Slice(0,rmax-1),_));
					Tensor<T> b(d1,vtdim);
					for (long i=0; i<c.dim(0); ++i)
						for (long j=0; j<c.dim(1); ++j)
							for (long k=0; k<rmax; ++k)
								b(i,j) += uu(i,k) * T(s(k)) * vtt(k,j);
					b -= aa;
					print("b.conforms c",b.conforms(c));
					print("early error:",b.absmax());
				}
#endif

				//        U = Tensor<T>(m,rmax);
				//        VT = Tensor<T>(rmax,n);

				// handle rank=0 explicitly
				if (r[d]) {

					// done with this dimension -- slice and deep-copy
					core[d-1]=madness::copy((u(Slice(0,c.dim(0)*rmax-1)))
							.reshape(c.dim(0),rmax)(_,Slice(0,rank-1)));
					core[d-1]=core[d-1].reshape(r[d-1],k,r[d]);

					// continue with the next dimension
					c=c(Slice(0,rank-1),_);

					for (int i=0; i<rank; ++i) {
						for (int j=0; j<c.dim(1); ++j) {
							c(i,j)*=s(i);
						}
					}

					if (d == dims.size()-1) core[d]=c;
				}
				else {
					zero_rank = true;
					core[d-1] = Tensor<T>(r[d-1],k,long(0));
					// iterate through the rest -- fast forward
					for(++d; d<dims.size(); ++d) {
						const long k=dims[d-1];
						core[d-1] = Tensor<T>(long(0),k,long(0));
					}
					core[dims.size()-1] = Tensor<T>(long(0),dims[dims.size()-1]);
				}
			}
			core[0]=core[0].fusedim(0);


		}


		/// turn this into an empty tensor with all cores properly shaped
		void zero_me() {
		    *this=TensorTrain<T>(this->dims());
		}

		/// return this multiplied by a scalar

		/// @return new tensor
		TensorTrain<T> operator*(const T& factor) const {
		    TensorTrain result=copy(*this);
		    result.scale(factor);
		    return result;
		}

		/// inplace addition of two Tensortrains; will increase ranks of this

		/// inefficient if many additions are performed, since it requires
		/// many calls of new.
		/// @param[in]	rhs	a TensorTrain to be added
		TensorTrain<T>& operator+=(const TensorTrain<T>& rhs) {
		    gaxpy(1.0,rhs,1.0);
		    return *this;
		}

        /// inplace subtraction of two Tensortrains; will increase ranks of this

        /// inefficient if many subtractions are performed, since it requires
        /// many calls of new.
        /// @param[in]  rhs a TensorTrain to be added
        TensorTrain<T>& operator-=(const TensorTrain<T>& rhs) {
            gaxpy(1.0,rhs,-1.0);
            return *this;
        }

		/// Inplace generalized saxpy ... this = this*alpha + other*beta
        TensorTrain<T>& gaxpy(T alpha, const TensorTrain<T>& rhs, T beta) {

			// make sure dimensions conform
			MADNESS_ASSERT(this->ndim()==rhs.ndim());

			if (this->zero_rank or (alpha==0.0)) {
				*this=rhs*beta;
			} else if (rhs.zero_rank or (beta==0.0)) {
				scale(alpha);
			} else {

				// special treatment for first border cores (k,r1)

			    // alpha and beta are only included in the first core(!)
				{
					long k=core[0].dim(0);
					long r1_this=core[0].dim(1);
					long r1_rhs=rhs.core[0].dim(1);
					Tensor<T> core_new(k,r1_this+r1_rhs);
					core_new(_,Slice(0,r1_this-1))=alpha*core[0];
					core_new(_,Slice(r1_this,r1_this+r1_rhs-1))=beta*rhs.core[0];
					core[0]=core_new;
				}

				// interior cores (r0,k,r1)
				for (std::size_t i=1; i<core.size()-1; ++i) {
					MADNESS_ASSERT(core[i].ndim()==3 or i==core.size()-1);
					long r0_this=core[i].dim(0);
					long r0_rhs=rhs.core[i].dim(0);
					long k=core[i].dim(1);
					long r1_this=core[i].dim(2);
					long r1_rhs=rhs.core[i].dim(2);
					Tensor<T> core_new(r0_this+r0_rhs,k,r1_this+r1_rhs);
					core_new(Slice(0,r0_this-1),_,Slice(0,r1_this-1))=core[i];
					core_new(Slice(r0_this,r0_this+r0_rhs-1),_,Slice(r1_this,r1_this+r1_rhs-1))=rhs.core[i];
					core[i]=core_new;
				}

				// special treatment for last border core (r0,k)
				{
					std::size_t d=core.size()-1;
					long r0_this=core[d].dim(0);
					long r0_rhs=rhs.core[d].dim(0);
					long k=core[d].dim(1);
					Tensor<T> core_new(r0_this+r0_rhs,k);
					core_new(Slice(0,r0_this-1),_)=core[d];
					core_new(Slice(r0_this,r0_this+r0_rhs-1),_)=rhs.core[d];
					core[d]=core_new;
				}
			}
			if (not verify()) MADNESS_EXCEPTION("ranks in TensorTrain inconsistent",1);
			return *this;
		}

        /// Inplace generalized saxpy with slices and without alpha

        /// return this = this(s1) + other(s2) * beta
        TensorTrain<T>& gaxpy(const std::vector<Slice>& s1,
                const TensorTrain<T>& rhs, T beta, const std::vector<Slice>& s2) {

            // make sure dimensions conform
            MADNESS_ASSERT(this->ndim()==rhs.ndim());

//            if (this->zero_rank) {
//                *this=rhs*beta;

//            } else {
            if (true) {
                // special treatment for first border cores (k,r1)

                // alpha and beta are only included in the first core(!)
                {
                    long k=core[0].dim(0);  // this is max dimension: k>=slice1; k>=slice2
                    long r1_this=core[0].dim(1);
                    long r1_rhs=rhs.core[0].dim(1);

                    Tensor<T> core_new(k,r1_this+r1_rhs);
                    if (r1_this>0) core_new(_,Slice(0,r1_this-1))=core[0];
                    if (r1_rhs>0)  core_new(s1[0],Slice(r1_this,r1_this+r1_rhs-1))=beta*rhs.core[0](s2[0],_);
                    core[0]=core_new;
                }

                // interior cores (r0,k,r1)
                for (std::size_t i=1; i<core.size()-1; ++i) {
                    MADNESS_ASSERT(core[i].ndim()==3 or i==core.size()-1);
                    long r0_this=core[i].dim(0);
                    long r0_rhs=rhs.core[i].dim(0);
                    long k=core[i].dim(1);
                    long r1_this=core[i].dim(2);
                    long r1_rhs=rhs.core[i].dim(2);
                    Tensor<T> core_new(r0_this+r0_rhs,k,r1_this+r1_rhs);
                    if (r1_this>0) core_new(Slice(0,r0_this-1),_,Slice(0,r1_this-1))=core[i];
                    if (r1_rhs>0)  core_new(Slice(r0_this,r0_this+r0_rhs-1),s1[i],Slice(r1_this,r1_this+r1_rhs-1))=rhs.core[i](_,s2[i],_);
                    core[i]=core_new;
                }

                // special treatment for last border core (r0,k)
                {
                    std::size_t d=core.size()-1;
                    long r0_this=core[d].dim(0);
                    long r0_rhs=rhs.core[d].dim(0);
                    long k=core[d].dim(1);
                    Tensor<T> core_new(r0_this+r0_rhs,k);
                    if (r0_this>0) core_new(Slice(0,r0_this-1),_)=core[d];
                    if (r0_rhs>0)  core_new(Slice(r0_this,r0_this+r0_rhs-1),s1[d])=rhs.core[d](_,s2[d]);
                    core[d]=core_new;
                }
            }
            if (not rhs.zero_rank) zero_rank=false;
            if (not verify()) MADNESS_EXCEPTION("ranks in TensorTrain inconsistent",1);
            return *this;
        }


		/// merge two dimensions into one

		/// merge dimension i and i+1 into new dimension i
		/// @param[in]	i	the first dimension
		void fusedim(const long i) {
			// core_new = left * right
			// (r1, k1*k2, r3) = sum_r2 (r1, k1, r2) * (r2, k2, r3)

			// determine index
			const int index=core[i].ndim()-2;	// (r-1, k, k, .. , k, r1)

			if (not zero_rank) core[i]=inner(core[i],core[i+1]);
			core[i]=core[i].fusedim(index);

			// shift all subsequent cores and remove the last one
			for (std::size_t d=i+1; d<core.size()-1; ++d) core[d]=core[d+1];
			core.pop_back();

		}


		/// reconstruct this to a full representation

		/// @param[in]	flat	return this in flat representation
		/// @return	this in full rank representation
		Tensor<T> reconstruct(const bool flat=false) const {

			if (zero_rank) {
				if (flat) {
					long size=1;
					for (int i=1; i<this->ndim(); ++i) size*=core[i].dim(1);
					return Tensor<T>(size);
				} else {
					std::vector<long> d(this->ndim());
					d[0]=core[0].dim(0);	// first core tensor has shape (k,r1)
					for (int i=1; i<this->ndim(); ++i) d[i]=core[i].dim(1);
					return Tensor<T>(d);
				}
			}

			Tensor<T> result=core.front();
			typename std::vector<Tensor<T> >::const_iterator it;

			for (it=++core.begin(); it!=core.end(); ++it) {
				result=inner(result,*it);
				if (flat) result=result.fusedim(0);
			}
			return result;
		}

		/// construct a two-mode representation (aka unnormalized SVD)

		/// @param[out] U The left singular vectors, ({i},rank)
		/// @param[out] VT The right singular vectors, (rank,{i})
		/// @param[out]	s Vector holding 1's, (rank)
		void two_mode_representation(Tensor<T>& U, Tensor<T>& VT,
				Tensor< typename Tensor<T>::scalar_type >& s) {

			// number of dimensions needs to be even
			MADNESS_ASSERT(ndim()%2==0);

			if (not zero_rank) {
			  typename std::vector<Tensor<T> >::const_iterator it1, it2;
			  U=core.front();
			  VT=core.back();
			  for (it1=++core.begin(), it2=--(--core.end()); it1<it2; ++it1, --it2) {
			    U=inner(U,*it1);
				VT=inner(*it2,VT);
			  }
			  s=Tensor< typename Tensor<T>::scalar_type >(VT.dim(0));
			  s=1.0;
			}
			else {
			  long dim1 = core.front().dim(0);
			  long dim2 = core.back().dim(1);
              for (int d1=1, d2=core.size()-2; d1<d2; ++d1, --d2) {
                dim1 *= core[d1].dim(1);
                dim2 *= core[d2].dim(1);
              }
              U = Tensor<T>(dim1,long(0));
              VT = Tensor<T>(long(0),dim2);
              s = Tensor< typename Tensor<T>::scalar_type >(VT.dim(0));
			}
		}

        /// recompress and truncate this TT representation

        /// this in recompressed TT form with optimal rank
        /// @param[in]  eps the truncation threshold
		template<typename R=T>
        typename std::enable_if<!std::is_arithmetic<R>::value, void>::type
        truncate(double eps) {
            MADNESS_EXCEPTION("no complex truncate in TensorTrain",1);
        }


		/// recompress and truncate this TT representation

		/// this in recompressed TT form with optimal rank
		/// @param[in]	eps	the truncation threshold
        template<typename R=T>
		typename std::enable_if<std::is_arithmetic<R>::value, void>::type
		truncate(double eps) {

		    // fast return
		    if (zero_rank) return;

		    for (long i=0; i<core.size(); ++i) if (ranks(i)==0) zero_rank=true;

		    if (zero_rank) {
		        zero_me();
		        return;
		    }


			eps=eps/sqrt(this->ndim());
            if (not verify()) MADNESS_EXCEPTION("ranks in TensorTrain inconsistent",1);

			// right-to-left orthogonalization (line 4)
			Tensor<T> L;
			long dims[TENSOR_MAXDIM];
			for (std::size_t d=core.size()-1; d>0; --d) {

				// save tensor structure
				const long ndim=core[d].ndim();
				for (int i=0; i<ndim; ++i) dims[i]=core[d].dim(i);

				// G(r0, k*r1)
				const long r0=core[d].dim(0);
				core[d]=core[d].reshape(r0,core[d].size()/r0);

				// decompose the core tensor (line 5)
				lq(core[d],L);  // might shrink the core
				dims[0]=std::min(r0,core[d].dim(0));
				core[d]=core[d].reshape(ndim,dims);

				// multiply to the left (line 6)
				core[d-1]=inner(core[d-1],L);
			}

			// left-to-right SVD (line 9)
			for (std::size_t d=0; d<core.size()-1; ++d) {

				// save tensor structure
				const long ndim=core[d].ndim();
				for (int i=0; i<ndim; ++i) dims[i]=core[d].dim(i);

				// reshape the core tensor (r0*k, r1)
                long r1=core[d].dim(core[d].ndim()-1);
//				long r1=core[d].dim(1);
				core[d]=core[d].reshape(core[d].size()/r1,r1);

				Tensor<T> U,VT;
				long rmax=std::min(core[d].dim(0),core[d].dim(1));
				Tensor< typename Tensor<T>::scalar_type > s(rmax);

				// decompose (line 10)
				svd(core[d],U,s,VT);

				// truncate the SVD
				int r_truncate=SRConf<T>::max_sigma(eps,rmax,s)+1;
				if (r_truncate==0) {
				    zero_me();
				    return;
				}
				U=madness::copy(U(_,Slice(0,r_truncate-1)));
				VT=madness::copy(VT(Slice(0,r_truncate-1),_));

				dims[ndim-1]=r_truncate;
				core[d]=U.reshape(ndim,dims);


				for (int i=0; i<VT.dim(0); ++i) {
					for (int j=0; j<VT.dim(1); ++j) {
						VT(i,j)*=s(i);
					}
				}

				// multiply to the right (line 11)
				core[d+1]=inner(VT,core[d+1]);

			}

            if (not verify()) MADNESS_EXCEPTION("ranks in TensorTrain inconsistent",1);
		}

		/// return the number of dimensions
		long ndim() const {return core.size();}

		/// return the number of coefficients in all core tensors
		long size() const {
			if (zero_rank) return 0;
			long n=0;
			typename std::vector<Tensor<T> >::const_iterator it;
			for (it=core.begin(); it!=core.end(); ++it) n+=it->size();
			return n;
		}

		/// return the size of this instance, including static memory for vectors and such
		long real_size() const {
			long n=this->size()*sizeof(T);
			n+=core.size() * sizeof(Tensor<T>);
			n+=sizeof(*this);
			return n;
		}

		/// return the number of entries in dimension i
		long dim(const int i) const {
			if (i==0) return core[0].dim(0);
			return core[i].dim(1);
		}

		/// return the dimensions of this tensor
		std::vector<long> dims() const {
		    std::vector<long> d(ndim());
            d[0]=core[0].dim(0);    // dim,rank
		    for (long i=1; i<ndim(); ++i) d[i]=core[i].dim(1);  // rank,dim,rank
		    return d;
		}

		/// check that the ranks of all core tensors are consistent
		bool verify() const {
		    if (core[0].dim(1)!=core[1].dim(0)) return false;
		    for (int d=2; d<ndim(); ++d) {
		        if (core[d-1].dim(2)!=core[d].dim(0)) return false;
		    }
		    for (const Tensor<T>& c : core) {
		        int size=1;
		        for (int i=0; i<c.ndim(); ++i) size*=c.dim(i);
		        if (size!=c.size()) return false;
		        if (not c.iscontiguous()) return false;
		    }
		    return true;
		}

		/// if rank is zero
		bool is_zero_rank() const {return zero_rank;}

		/// return the TT ranks
		std::vector<long> ranks() const {
			if (zero_rank) return std::vector<long>(core.size()-1,0);
			std::vector<long> r(core.size()-1);
			for (std::size_t i=0; i<r.size(); ++i) r[i]=core[i+1].dim(0);
			return r;
		}

        /// return the TT ranks for dimension i (to i+1)
       long ranks(const int i) const {
            if (zero_rank) return 0;
            if (i==0) {
                return core[0].dim(1);
            } else if (i<core.size()) {
                return core[i].dim(2);
            } else {
                print("ndim ",ndim());
                print("i    ",i);
                MADNESS_EXCEPTION("requested invalid rank in TensorTrain",1);
                return 0;
            }
        }

        /// returns the Frobenius norm
        float_scalar_type normf() const {
            return sqrt(float_scalar_type(std::abs(trace(*this))));
        };

        /// scale this by a number

        /// @param[in]  fac the factor to multiply
        /// @return *this * fac
        void scale(T fac) {
            if (not zero_rank and (core.size()>0)) core.front().scale(fac);
        }

        /// Returns a pointer to the internal data

        /// @param[in]  ivec    index of core vector to which the return values points
        T* ptr(const int ivec=0) {
            if (core.size()) return core[ivec].ptr();
            return 0;
        }

        /// Returns a pointer to the internal data

        /// @param[in]  ivec    index of core vector to which the return values points
        const T* ptr(const int ivec=0) const {
            if (core.size()) return core[ivec].ptr();
            return 0;
        }


        /// Return the trace of two tensors, no complex conjugate involved

        /// @return <this | B>
        template <class Q>
        TENSOR_RESULT_TYPE(T,Q) trace(const TensorTrain<Q>& B) const {
            if (TensorTypeData<T>::iscomplex) MADNESS_EXCEPTION("no complex trace in TensorTrain, sorry",1);
            if (TensorTypeData<Q>::iscomplex) MADNESS_EXCEPTION("no complex trace in TensorTrain, sorry",1);

            typedef TENSOR_RESULT_TYPE(T,Q) resultT;

            // alias
            const TensorTrain<T>& A=*this;

            MADNESS_ASSERT(A.ndim()==B.ndim());   // number of dimensions

            // fast return
            if (A.zero_rank or B.zero_rank) return resultT(0.0);

            // set up temporary tensors for intermediates
            long size1=A.ranks(0)*B.ranks(0);
            long size2=A.ranks(0)*A.dim(0);

            for (int d=1; d<A.ndim(); ++d) {
                size1=std::max(size1,A.ranks(d)*B.ranks(d));
                size2=std::max(size2,A.ranks(d)*B.ranks(d-1)*A.dim(d));
            }
            Tensor<resultT> tmp1(size1), tmp2(size2);       // scratch
            Tensor<resultT> Aprime, AB;                     // for flat views of tmp

            // loop over all dimensions but the last one
            for (int d=0; d<A.ndim()-1; ++d) {

                // contract cores to matrix AB(rank(A), rank(B))
                // AB(ra1,rb1) = sum_(r0,i0) A(ra0,i0,ra1) B(rb0,i0,rb1)
                // index dimension is the first one: core((r * i), r)

                // calculate dimensions
                long rA= (d==0) ? A.core[d].dim(1) : Aprime.dim(2);     //  r_d (A)
                long rB= (d==0) ? B.core[d].dim(1) : B.core[d].dim(2);                           //  r_d (B)
                MADNESS_ASSERT(rA*rB<=size1);
                if (d>0) tmp1(Slice(0,rA*rB-1))=0.0;         // zero out old stuff

//                Tensor<resultT> AB;
                if (d==0) {
//                    AB=inner(A.core[d],B.core[d],0,0);
                    inner_result(A.core[d],B.core[d],0,0,tmp1);
                } else {
//                    AB=inner(Aprime.fusedim(0),B.core[d].fusedim(0),0,0);
                    inner_result(Aprime.fusedim(0),B.core[d].fusedim(0),0,0,tmp1);
                }
                AB=tmp1(Slice(0,rA*rB-1)).reshape(rA,rB);

                // contract temp matrix AB into core of A of dimension i+1
                // into temp core A_i+1
                // Atmp(rank(B1), i(2)) = sum_ rank(A1) ABi(rank(A1),rank(B1)) A.core[1](rank(A1),i(2),rank(A2))

                // calculate dimensions
                long d1=AB.dim(1);              //  r_d
                long d2=A.core[d+1].dim(1);     //  i_d
                long d3=A.core[d+1].dim(2);     //  r_{d+1}
                MADNESS_ASSERT(d1*d2*d3<=size2);

                // repeated zero-ing probably much faster than reallocation
                //                Aprime=inner(AB,A.core[d+1],0,0);
                if (d>0) tmp2(Slice(0,d1*d2*d3-1))=0.0;
                inner_result(AB,A.core[d+1],0,0,tmp2);
                Aprime=tmp2(Slice(0,d1*d2*d3-1)).reshape(d1,d2,d3);

            }

            // special treatment for the last dimension
            resultT result=Aprime.trace(B.core[ndim()-1]);
            return result;
        }


        template <typename R, typename Q>
        friend TensorTrain<TENSOR_RESULT_TYPE(R,Q)> transform(
                const TensorTrain<R>& t, const Tensor<Q>& c);

        template <typename R, typename Q>
        friend TensorTrain<TENSOR_RESULT_TYPE(R,Q)> general_transform(
                const TensorTrain<R>& t, const Tensor<Q> c[]);

        template <typename R, typename Q>
        friend TensorTrain<TENSOR_RESULT_TYPE(R,Q)> transform_dir(
                const TensorTrain<R>& t, const Tensor<Q>& c, const int axis);

	};


	/// transform each dimension with the same operator matrix

    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c(i',i) c(j',j) c(k',k) ...
    /// TODO: merge this with general_transform
	template <class T, class Q>
    TensorTrain<TENSOR_RESULT_TYPE(T,Q)> transform(const TensorTrain<T>& t,
            const Tensor<Q>& c) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        // fast return if possible
        if (t.zero_rank or (t.ndim()==0)) return TensorTrain<resultT>(t.dims());

        const long ndim=t.ndim();

        TensorTrain<resultT> result;
        result.zero_rank=false;
        result.core.resize(ndim);
        // special treatment for first core(i1,r1) and last core (rd-1, id)
        result.core[0]=inner(c,t.core[0],0,0);
        if (ndim>1) result.core[ndim-1]=inner(t.core[ndim-1],c,1,0);

        // other cores have dimensions core(r1,i2,r2);

        // set up scratch tensor
        long size=0;
        for (int d=1; d<ndim-1; ++d) size=std::max(size,t.core[d].size());
        Tensor<resultT> tmp(size);

        for (int d=1; d<ndim-1; ++d) {
            long r1=t.core[d].dim(0);
            long i2=t.core[d].dim(1);
            long r2=t.core[d].dim(2);

            // zero out old stuff from the scratch tensor
            if (d>1) tmp(Slice(0,r1*i2*r2-1))=0.0;
            inner_result(t.core[d],c,1,0,tmp);
            result.core[d]=copy(tmp(Slice(0,r1*i2*r2-1)).reshape(r1,r2,i2).swapdim(1,2));
        }
        return result;
    }


    /// Transform all dimensions of the tensor t by distinct matrices c

    /// \ingroup tensor
    /// Similar to transform but each dimension is transformed with a
    /// distinct matrix.
    /// \code
    /// result(i,j,k...) <-- sum(i',j', k',...) t(i',j',k',...) c[0](i',i) c[1](j',j) c[2](k',k) ...
    /// \endcode
    /// The first dimension of the matrices c must match the corresponding
    /// dimension of t.
    template <class T, class Q>
    TensorTrain<TENSOR_RESULT_TYPE(T,Q)> general_transform(const TensorTrain<T>& t,
            const Tensor<Q> c[]) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        // fast return if possible
        if (t.zero_rank or (t.ndim()==0)) return TensorTrain<resultT>(t.dims());

        const long ndim=t.ndim();

        TensorTrain<resultT> result;
        result.zero_rank=false;
        result.core.resize(ndim);
        // special treatment for first core(i1,r1) and last core (rd-1, id)
        result.core[0]=inner(c[0],t.core[0],0,0);
        if (ndim>1) result.core[ndim-1]=inner(t.core[ndim-1],c[ndim-1],1,0);

        // other cores have dimensions core(r1,i2,r2);

        // set up scratch tensor
        long size=0;
        for (int d=1; d<ndim-1; ++d) size=std::max(size,t.core[d].size());
        Tensor<resultT> tmp(size);

        for (int d=1; d<ndim-1; ++d) {
            long r1=t.core[d].dim(0);
            long i2=t.core[d].dim(1);
            long r2=t.core[d].dim(2);

            // zero out old stuff from the scratch tensor
            if (d>1) tmp(Slice(0,r1*i2*r2-1))=0.0;
            inner_result(t.core[d],c[d],1,0,tmp);
            result.core[d]=copy(tmp(Slice(0,r1*i2*r2-1)).reshape(r1,r2,i2).swapdim(1,2));
        }
        return result;
    }

    /// Transforms one dimension of the tensor t by the matrix c, returns new contiguous tensor

    /// \ingroup tensor
    /// \code
    /// transform_dir(t,c,1) = r(i,j,k,...) = sum(j') t(i,j',k,...) * c(j',j)
    /// \endcode
    /// @param[in] t Tensor to transform (size of dimension to be transformed must match size of first dimension of \c c )
    /// @param[in] c Matrix used for the transformation
    /// @param[in] axis Dimension (or axis) to be transformed
    /// @result Returns a new tensor train
    template <class T, class Q>
    TensorTrain<TENSOR_RESULT_TYPE(T,Q)> transform_dir(const TensorTrain<T>& t,
            const Tensor<Q>& c, const int axis) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;

        // fast return if possible
        if (t.zero_rank or (t.ndim()==0)) return TensorTrain<resultT>(t.dims());

        const long ndim=t.ndim();
        MADNESS_ASSERT(axis<ndim and axis>=0);
        MADNESS_ASSERT(c.ndim()==2);

        TensorTrain<resultT> result=copy(t);

        if (axis==0) {
            result.core[0]=inner(c,t.core[0],0,0);
        } else if (axis==ndim-1) {
            result.core[ndim-1]=inner(t.core[ndim-1],c,1,0);
        } else {
            Tensor<resultT> tmp=inner(t.core[axis],c,1,0);  // G~(r1,r2,i')
            result.core[axis]=copy(tmp.swapdim(1,2));
        }
        return result;

    }



}

#endif /* TENSORTRAIN_H_ */
