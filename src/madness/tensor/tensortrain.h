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

		/// holding the core tensors of a tensor train
		/// the tensors have the shape (k,r0) (r0,k,r1) (r1,k,r2) .. (rn-1,k)
		std::vector<Tensor<T> > core;
		/// true if rank is zero
		bool zero_rank;

	public:

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
			Tensor<T> c=copy(t);

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
					core[d-1]=copy((u(Slice(0,c.dim(0)*rmax-1)))
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

		/// inplace addition of two Tensortrains; will increase ranks of this

		/// inefficient if many additions are performed, since it requires
		/// many calls of new.
		/// @param[in]	rhs	a TensorTrain to be added
		TensorTrain<T>& operator+=(const TensorTrain<T>& rhs) {

			// make sure dimensions conform
			MADNESS_ASSERT(this->ndim()==rhs.ndim());

			if (this->zero_rank) {
				*this=rhs;
			} else if (rhs.zero_rank) {
				;
			} else {

				// special treatment for first border cores (k,r1)
				{
					long k=core[0].dim(0);
					long r1_this=core[0].dim(1);
					long r1_rhs=rhs.core[0].dim(1);
					Tensor<T> core_new(k,r1_this+r1_rhs);
					core_new(_,Slice(0,r1_this-1))=core[0];
					core_new(_,Slice(r1_this,r1_this+r1_rhs-1))=rhs.core[0];
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
		/// @param[in]	eps	the truncation threshold
		void truncate(double eps) {
			eps=eps/sqrt(this->ndim());

			// right-to-left orthogonalization (line 4)
			Tensor<T> R;
			long dims[TENSOR_MAXDIM];
			for (std::size_t d=core.size()-1; d>0; --d) {

				// save tensor structure
				const long ndim=core[d].ndim();
				for (int i=0; i<ndim; ++i) dims[i]=core[d].dim(i);

				// G(r0, k*r1)
				const long r0=core[d].dim(0);
				core[d]=core[d].reshape(r0,core[d].size()/r0);

				// decompose the core tensor (line 5)
				lq(core[d],R);
				core[d]=core[d].reshape(ndim,dims);

				// multiply to the left (line 6)
				core[d-1]=inner(core[d-1],R);
			}

			// left-to-right SVD (line 9)
			for (std::size_t d=0; d<core.size()-1; ++d) {

				// save tensor structure
				const long ndim=core[d].ndim();
				for (int i=0; i<ndim; ++i) dims[i]=core[d].dim(i);

				// reshape the core tensor (r0*k, r1)
				long r1=core[d].dim(core[d].ndim()-1);
				core[d]=core[d].reshape(core[d].size()/r1,r1);

				Tensor<T> U,VT;
				long rmax=std::min(core[d].dim(0),core[d].dim(1));
				Tensor< typename Tensor<T>::scalar_type > s(rmax);

				// decompose (line 10)
				svd(core[d],U,s,VT);

				// truncate the SVD
				int r_truncate=SRConf<T>::max_sigma(eps,rmax,s)+1;
				U=copy(U(_,Slice(0,r_truncate-1)));
				VT=copy(VT(Slice(0,r_truncate-1),_));

				dims[ndim-1]=r_truncate;
				core[d]=U.reshape(ndim,dims);


				for (int i=0; i<VT.dim(0); ++i) {
					for (int j=0; j<VT.dim(1); ++j) {
						VT(i,j)*=s(j);
					}
				}

				// multiply to the right (line 11)
				core[d+1]=inner(VT,core[d+1]);

			}

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

		/// if rank is zero
		bool is_zero_rank() const {return zero_rank;}

		/// return the TT ranks
		std::vector<long> ranks() const {
			if (zero_rank) return std::vector<long>(0,core.size()-1);
			std::vector<long> r(core.size()-1);
			for (std::size_t i=0; i<r.size(); ++i) r[i]=core[i+1].dim(0);
			return r;
		}

	};


}

#endif /* TENSORTRAIN_H_ */
