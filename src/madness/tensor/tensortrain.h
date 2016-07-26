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

    /// decompose the input tensor A into U and V, skipping singular vectors of
    /// small singular values. One deep copy/allocation for U is performed

    /// @param[inout]   A the input matrix, on exit the matrix VT (right sing. vectors)
    /// @param[out]     U contiguous new tensor holding the left sing. vectors
    /// @param[in]      thresh threshold for truncation of the singular vectors
    /// @param[in]      s   scratch tensor for the singular values, dimension min(n,m)
    /// @param[in]      scr scratch tensor
    /// the dimension of the scratch tensor may be computed as
    /// long lwork=std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n)) + n*m;
    template<typename T>
    long rank_revealing_decompose(Tensor<T>& A, Tensor<T>& U,
            const double thresh, Tensor< typename Tensor<T>::scalar_type > & s,
            Tensor<T>& scr) {

        MADNESS_ASSERT(A.ndim()==2);    // must be a matrix
        const long n=A.dim(0);
        const long m=A.dim(1);
        const long rmax=std::min(n,m);
        long lwork=std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n));

        // set up contiguous scratch arrays
        MADNESS_ASSERT(scr.size()>lwork+n*m);
        scr=scr.flat();
        Tensor<T> work=scr(Slice(0,lwork-1));
        Tensor<T> utmp=scr(Slice(lwork,lwork+n*m-1));
        Tensor<T> dummy;    // real dummy

        svd_result(A,utmp,s,dummy,work);

        // this is rank_right
        const long R1=SRConf<T>::max_sigma(thresh,rmax,s)+1;

        // skip if rank=0
        if (R1>0) {

            U=madness::copy((utmp(Slice(0,n*rmax-1)))
                    .reshape(n,rmax)(_,Slice(0,R1-1)));

            A=A(Slice(0,R1-1),_);

            // continue with the next dimension
            for (int j=0; j<m; ++j) {
                for (int i=0; i<R1; ++i) {
                    A(i,j)*=s(i);
                }
            }
        } else {
            U=Tensor<T>(n,0l);
        }
        return R1;
    }



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


        /// ctor for a TensorTrain, with core tensors explicitly given

        /// the core tensors must have the shape (k1,r1) (r1,k2,r2) .. (r2,k3)
		/// @param[in]  core    vector of core tensors, properly shaped
        TensorTrain(const std::vector<Tensor<T> >& c) : core(c) {
            zero_rank=false;

            // check for zero ranks
            for (int d=1; d<core.size(); ++d) if (core[d].dim(0)==0) zero_me();
        }


		/// copy constructor, shallow
		TensorTrain(const TensorTrain& other) : core(other.core),
		        zero_rank(other.zero_rank) {
		}

		/// assigment operator
		TensorTrain& operator=(const TensorTrain& other) {
		    if (this!=&other) {
                zero_rank=other.zero_rank;
                core=other.core;
		    }
            return *this;
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

        /// serialize this
        template <typename Archive>
        void serialize(Archive& ar) {
            long dim=ndim();
            ar & zero_rank & dim;

            // no empty tensor
            if (dim>0) {

                ar & core;

                // need this because tensor serialization does not preserve the
                // dimensions if the tensor is empty, but we need to!

                // existing tensor filled with zeros
                if (zero_rank) {
                    std::vector<long> d=this->dims();
                    ar & d;
                    zero_me(d);
                }
            }
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
		    if (other.zero_rank) {
		        std::vector<long> dims(s.size());
		        for (int i=0; i<dims.size(); ++i) dims[i]=s[i].end-s[i].start+1;
		        return TensorTrain(dims);
		    }

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

        /// turn this into an empty tensor with all cores properly shaped
        void zero_me(const std::vector<long>& dim) {
            *this=TensorTrain<T>(dim);
        }

        /// assign a number to this tensor
        TensorTrain<T>& operator=(const T& number) {

            // tensor will have ranks = 1 all over
            core[0]=Tensor<T>(dim(0),1);
            for (int i=1; i<ndim()-1; ++i) core[i]=Tensor<T>(1,dim(i),1);
            if (ndim()>1) core[ndim()-1]=Tensor<T>(1,dim(ndim()-1));

            core[0]=number;
            for (int i=1; i<ndim(); ++i) core[i]=1.0;

            zero_rank=false;
            return *this;
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
			} else if (ndim()==1) {     // simple algorithm for ndim=1
			    core[0].gaxpy(alpha,rhs.core[0],beta);
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

        /// compute the Hadamard product of two TensorTrains

        /// see Sec. 4.2 of the TT paper
        /// @return *this with the ranks of the input tensors multiplied!
        TensorTrain<T>& emul(const TensorTrain<T>& other) {

            // consistency checks
            MADNESS_ASSERT(ndim()==other.ndim());
            for (int i=0; i<ndim(); ++i) MADNESS_ASSERT(dim(i)==other.dim(i));

            // set up the result tensor
            std::vector<long> cranks(ndim()-1);       // ranks of result C
            for (int i=0; i<ndim()-1; ++i) cranks[i]=ranks(i)*other.ranks(i);

            // check for large sizes
            for (int i=0; i<ndim()-2; ++i) {
                long size=cranks[i]*cranks[i+1]*dim(i);
                if (size>10000000)
                    MADNESS_EXCEPTION("emul in TensorTrain too large -- use full rank tenspr",1);
            }
            TensorTrain result(dims());

            // fast return for zero ranks
            if (zero_rank or other.zero_rank) {
                ;
            } else if (ndim()==1) {
                // fast return for one core only (all ranks equal 1)
                result.core[0]=copy(core[0]);
                result.core[0].emul(other.core[0]);

            } else {
                result.zero_rank=false;
                // compute outer product for each core for each k

                // special treatment for first core
                result.core[0]=Tensor<T>(dim(0),cranks[0]);
                for (int k=0; k<dim(0); ++k) {
                    Tensor<T> a1=core[0](Slice(k,k),_);
                    Tensor<T> b1=other.core[0](Slice(k,k),_);
                    result.core[0](Slice(k,k),_)=outer(a1,b1).reshape(1,cranks[0]);
                    // before reshape rhs has dimensions (1, r1, 1, r1')
                    // after reshape rhs has dimensions (1, r1*r1')
                }

                for (int i=1; i<ndim()-1; ++i) {
                    result.core[i]=Tensor<T>(cranks[i-1],dim(i),cranks[i]);
                    for (int k=0; k<dim(i); ++k) {
                        Tensor<T> a1=core[i](_,Slice(k,k),_).fusedim(1);    // (r1,r2)
                        Tensor<T> b1=other.core[i](_,Slice(k,k),_).fusedim(1);  // (r1',r2')
                        Tensor<T> tmp=copy(outer(a1,b1).swapdim(1,2));  // make contiguous
                        result.core[i](_,Slice(k,k),_)=tmp.reshape(cranks[i-1],1,cranks[i]);
                        // before swap/fuse/splitdim rhs has dimensions (r1, r2, r1', r2')
                        // after fusedim/splitdim rhs has dimensions (r1*r1', 1, r2*r2')
                    }
                }

                // special treatment for last core
                const long n=ndim()-1;
                result.core[n]=Tensor<T>(cranks[n-1],dim(n));
                for (int k=0; k<dim(n); ++k) {
                    Tensor<T> a1=core[n](_,Slice(k,k));
                    Tensor<T> b1=other.core[n](_,Slice(k,k));
                    result.core[n](_,Slice(k,k))=outer(a1,b1).reshape(cranks[n-1],1);
                    // before reshape rhs has dimensions (r1,1, r1',1)
                    // after reshape rhs has dimensions (r1*r1', 1)

                }
            }
            *this=result;
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

			if (not zero_rank) {
			    core[i]=inner(core[i],core[i+1]);
	            core[i]=core[i].fusedim(index);
			} else {
			    if (i==0) { // (k1*k2, r2=0)
	                core[i]=Tensor<T>(core[i].dim(0)*core[i+1].dim(1),0l);
			    } else {                /// (r1=0, k1*k2, r2=0)
                    core[i]=Tensor<T>(0l,core[i].dim(0)*core[i+1].dim(1),0l);
			    }
			}

			// shift all subsequent cores and remove the last one
			for (std::size_t d=i+1; d<core.size()-1; ++d) core[d]=core[d+1];
			core.pop_back();

		}

        /// Returns new view/tensor splitting dimension \c i as \c dimi0*dimi1
		/// to produce conforming d+1 dimension tensor

		/// @param[in]  idim    the dimension to be split
		/// @param[in]  k1      new first dimension of idim
		/// @param[in]  k2      new second dimension of idim
		/// @param[in]  eps     threshold for SVD (choose negative to keep all terms)
        /// @return new deep copy of this with split dimensions
        TensorTrain<T> splitdim(long idim, long k1, long k2, const double eps) const {
            // core_new = left * right
            // (r1, k1*k2, r3) = sum_r2 (r1, k1, r2) * (r2, k2, r3)

            // check for consistency
            MADNESS_ASSERT(k1*k2==dim(idim));

            if (zero_rank) {
                std::vector<long> newdims(this->ndim()+1);
                for (long i=0; i<idim; ++i) newdims[i]=this->dim(i);
                newdims[idim]=k1;
                newdims[idim+1]=k2;
                for (long i=idim+1; i<ndim(); ++i) newdims[i+1]=dim(i);
                return TensorTrain(newdims);
            }

            TensorTrain<T> result;

            long r1= (idim==0) ? 1 : ranks(idim-1);       // left-side rank
            long r2= (idim==ndim()-1) ? 1 : ranks(idim);  // right-side rank
            long k12=dim(idim);

            Tensor<T> A=core[idim].reshape(r1*k1,r2*k2);
            Tensor<T> U,VT;
            Tensor< typename Tensor<T>::scalar_type > s(k12);
            svd(A,U,s,VT);

            // this is the new interior rank
            long r=SRConf<T>::max_sigma(eps,std::min(A.dim(0),A.dim(1)),s)+1;

            // handle rank=0 explicitly
            if (r==0) {
                std::vector<long> newdims(this->ndim()+1);
                for (long i=0; i<idim; ++i) newdims[i]=this->dim(i);
                newdims[idim]=k1;
                newdims[idim+1]=k2;
                for (long i=idim+1; i<ndim(); ++i) newdims[i+1]=dim(i);
                return TensorTrain(newdims);
            } else {

                // convolve the singular values into the right singular vectors
                for (int ii=0; ii<r; ++ii) {
                    for (int j=0; j<VT.dim(1); ++j) {
                        VT(ii,j)*=s(ii);
                    }
                }

                for (long ii=0; ii<idim; ++ii) result.core.push_back(copy(core[ii]));
                result.core.push_back(copy(U(_,Slice(0,r-1))).reshape(r1,k1,r));
                result.core.push_back(copy(VT(Slice(0,r-1),_)).reshape(r,k2,r2));
                for (long ii=idim+1; ii<ndim(); ++ii) result.core.push_back(core[ii]);

                // post-processing
                if (result.core.front().ndim()==3) result.core.front()=result.core.front().fusedim(0);
                if (result.core.back().ndim()==3) result.core.back()=result.core.back().fusedim(1);
                result.zero_rank=false;
            }
            return result;

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
		        if (ndim()==2) {
                    U = Tensor<T>(dim(0),long(0));
                    VT = Tensor<T>(long(0),dim(1));
		        } else if (ndim()==4) {
                    U = Tensor<T>(dim(0),dim(1),long(0));
                    VT = Tensor<T>(long(0),dim(2),dim(3));
		        } else if (ndim()==6) {
                    U = Tensor<T>(dim(0),dim(1),dim(2),long(0));
                    VT = Tensor<T>(long(0),dim(3),dim(4),dim(5));
		        } else {
		            MADNESS_EXCEPTION("ndim>6 in tensortrain::two_mode_representation",1);
		        }
		        s = Tensor< typename Tensor<T>::scalar_type >(0l);
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

		    for (long i=0; i<core.size()-1; ++i) if (ranks(i)==0) zero_rank=true;

		    if (zero_rank) {
		        zero_me();
		        return;
		    }

		    std::vector<long> tt_dims=this->dims();

		    eps=eps/sqrt(this->ndim());
		    if (not verify()) MADNESS_EXCEPTION("ranks in TensorTrain inconsistent",1);

		    // right-to-left orthogonalization (line 4)
		    // get maximum rank and maximum k-value
		    // cores are (k0,r0), (r0,k1,r1), (r1,k2,r2), ... (rd-1,kd)
		    long rmax = core[0].dim(1);
		    long kmax = core[0].dim(0);
		    for(size_t i=1;i<core.size();i++){
		        rmax = std::max(rmax,core[i].dim(0));
		        kmax = std::max(kmax,core[i].dim(1));
		    }

		    Tensor<T> L;//L_buffer(rmax,rmax*kmax);
		    Tensor<T> lq_tau(rmax);
		    long max_rk = std::max(rmax,kmax);
		    long lq_work_dim = 2*max_rk+(max_rk+1)*64;
		    Tensor<T> lq_work(lq_work_dim);
		    Tensor<T> L_buffer(max_rk,max_rk);
		    long dimensions[TENSOR_MAXDIM];
		    // last tensor differs in dimension (first tensor also, but we dont reach it in the loop)
		    if(core.size()>1){
		        const long n_dim = core.back().ndim();
		        for (int i=0; i<n_dim; ++i) dimensions[i]=core.back().dim(i);

		        const long r0 = core.back().dim(0);
		        const long r1 = core.back().size()/r0;
		        core.back()=core.back().reshape(r0,r1);

		        // assignement of L with the L_buffer tensor
		        // works only if the bool for lq_result is set to false
		        {
		            long r_rows= (core.back().dim(1)>=core.back().dim(0)) ? core.back().dim(0) : core.back().dim(1);
		            long r_cols=core.back().dim(0);
		            L = L_buffer(Slice(0,r_cols-1),Slice(0,r_rows-1));
		            L = 0.0;
		        }
		        lq_result(core.back(),L,lq_tau,lq_work,false);
		        //Tensor<T> L = L_buffer(Slice(0,r0-1),Slice(0,r1-1));

		        dimensions[0]=std::min(dimensions[0],core.back().dim(0));
		        core.back()=core.back().reshape(n_dim,dimensions);
		        // multiply to the left (line 6)
		        core[core.size()-2]=inner(core[core.size()-2],L);

		    }


		    for (std::size_t d=core.size()-2; d>0; --d) {

		        // save tensor structure
		        const long ndim=core[d].ndim();
		        for (int i=0; i<ndim; ++i) dimensions[i]=core[d].dim(i);

		        // G(r0, k*r1)
		        const long r0=core[d].dim(0);
		        const long r1=core[d].size()/r0;
		        core[d]=core[d].reshape(r0,r1);

		        // decompose the core tensor (line 5)
		        //lq(core[d],L_buffer);  // might shrink the core

		        // assignement of L with the inner_buffer tensor
		        // works only if the bool for lq_result is set to false
		        {
		            long r_rows= (core[d].dim(1)>=core[d].dim(0)) ? core[d].dim(0) : core[d].dim(1);
		            long r_cols=core[d].dim(0);
		            L = L_buffer(Slice(0,r_cols-1),Slice(0,r_rows-1));
		            L = 0.0;
		        }

		        // workaround for LQ decomposition to avoid reallocations
		        //L_buffer = 0.0;
		        lq_tau = 0.0;
		        lq_work = 0.0;
		        lq_result(core[d],L,lq_tau,lq_work,false);
		        // slice L to the right size
		        //Tensor<T> L = L_buffer(Slice(0,r0-1),Slice(0,r1-1));

		        dimensions[0]=std::min(r0,core[d].dim(0));
		        core[d]=core[d].reshape(ndim,dimensions);

		        // multiply to the left (line 6)
		        core[d-1]=inner(core[d-1],L);
		    }

		    // svd buffer tensor (see svd_results in lapack.cc)
		    long m =rmax*kmax;
		    long n =rmax;
		    long k =std::min<long>(m,n);
		    long svd_buffer_size = std::max<long>(3*std::min<long>(m,n)+std::max<long>(m,n),5*std::min<long>(m,n)-4)*32;
		    Tensor<T> svd_buffer(svd_buffer_size);
		    Tensor<T> U_buffer(m,k);
		    Tensor<T> dummy;
		    Tensor< typename Tensor<T>::scalar_type > s_buffer(k);

		    // left-to-right SVD (line 9)
		    for (std::size_t d=0; d<core.size()-1; ++d) {

		        // save tensor structure
		        const long ndim=core[d].ndim();
		        for (int i=0; i<ndim; ++i) dimensions[i]=core[d].dim(i);

		        // reshape the core tensor (r0*k, r1)
		        long r1=core[d].dim(core[d].ndim()-1);
		        //				long r1=core[d].dim(1);
		        core[d]=core[d].reshape(core[d].size()/r1,r1);

		        Tensor<T> U,VT;
		        long ds=std::min(core[d].dim(0),core[d].dim(1));
		        s_buffer = 0.0;
		        //Tensor< typename Tensor<T>::scalar_type > s(ds)
		        Tensor< typename Tensor<T>::scalar_type > s = s_buffer(Slice(0,ds-1));

		        // decompose (line 10)
		        //svd(core[d],U,s,VT);
		        // get the dimensions of U and V
		        long du = core[d].dim(0);
		        long dv = core[d].dim(1);
		        U_buffer = 0.0;
		        svd_buffer = 0.0;
		        U_buffer = U_buffer.flat();
		        // VT is written on core[d] input
		        svd_result(core[d],U_buffer,s,dummy,svd_buffer);

		        // truncate the SVD
		        int r_truncate=SRConf<T>::max_sigma(eps,ds,s)+1;
		        if (r_truncate==0) {
		            zero_me(tt_dims);
		            return;
		        }

		        // make tensors contiguous
		        U=madness::copy(U_buffer(Slice(0,(du*ds)-1)).reshape(du,ds)(_,Slice(0,r_truncate-1)));
		        VT=madness::copy(core[d](Slice(0,r_truncate-1),Slice(0,dv-1)));

		        dimensions[ndim-1]=r_truncate;
		        core[d]=U.reshape(ndim,dimensions);

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
		    if (core.size()<2) return true; // no ranks
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
			MADNESS_ASSERT(is_tensor());
			std::vector<long> r(core.size()-1);
			for (std::size_t i=0; i<r.size(); ++i) r[i]=core[i+1].dim(0);
			return r;
		}

        /// return the TT ranks for dimension i (to i+1)
       long ranks(const int i) const {
            if (zero_rank) return 0;
            if (i<core.size()-1) {
                return core[i].dim(core[i].ndim()-1);
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

        /// check if this is a tensor (r,k,r)
        bool is_tensor() const {
            if (ndim()>0) return (core[0].ndim()==2);
            return false;
        }

        /// check if this is an operator (r,k',k,r)
        bool is_operator() const {
            return (!is_tensor());
        }

        /// convert this into a tensor representation (r,k,r)
        TensorTrain<T>& make_tensor() {
            if (is_tensor()) return *this;

            long nd=this->ndim();
            core[0]=core[0].fusedim(0);         // (k,k,r) -> (k,r)
            core[nd-1]=core[nd-1].fusedim(1);   // (r,k,k) -> (r,k)
            for (int i=1; i<nd-1; ++i) core[i]=core[i].fusedim(1);  // (r,k',k,r) -> (r,k,r)
            return *this;

        }

        /// convert this into an operator representation (r,k',k,r)
        TensorTrain<T>& make_operator() {
            if (is_operator()) return *this;

            long nd=this->ndim();

            long k2=core[0].dim(0);
            long k0=sqrt(k2);
            MADNESS_ASSERT(k0*k0==k2);
            MADNESS_ASSERT(core[0].ndim()==2);
            core[0]=core[0].splitdim(0,k0,k0);         // (k*k,r) -> (k,k,r)

            k2=core[nd-1].dim(1);
            k0=sqrt(k2);
            MADNESS_ASSERT(k0*k0==k2);
            MADNESS_ASSERT(core[nd-1].ndim()==2);
            core[nd-1]=core[nd-1].splitdim(1,k0,k0);   // (r,k*k) -> (r,k,k)

            for (int i=1; i<nd-1; ++i) {
                k2=core[i].dim(1);
                k0=sqrt(k2);
                MADNESS_ASSERT(k0*k0==k2);
                MADNESS_ASSERT(core[i].ndim()==3);
                core[i]=core[i].splitdim(1,k0,k0);  // (r,k*k,r) -> (r,k',k,r)
            }
            return *this;
        }


        /// reference to the internal core
        Tensor<T>& get_core(const int i) {
            MADNESS_ASSERT(i<ndim());
            return core[i];
        }

        /// const reference to the internal core
        const Tensor<T>& get_core(const int i) const {
            MADNESS_ASSERT(i<ndim());
            return core[i];
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
            if (A.ndim()==1) {
                return A.core[0].trace(B.core[0]);
            }

            // set up temporary tensors for intermediates
            long size1=A.ranks(0)*B.ranks(0);                   // first contraction
            long size2=B.ranks(ndim()-2)*A.dim(ndim()-1);       // last contraction

            for (int d=1; d<A.ndim()-1; ++d) {
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

        template <typename R, typename Q>
        friend TensorTrain<TENSOR_RESULT_TYPE(R,Q)> outer(
                const TensorTrain<R>& t1, const TensorTrain<Q>& t2);

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


    /// apply an operator in TT format on a tensor in TT format

    /// @param[in]  op  operator in TT format ..(r_1,k',k,r_2)..
    /// @param[in]  t   tensor in TT format  ..(r_1,k',r_2)..
    /// the result tensor will be
    /// .. (r_1,k,r_2) = \sum_k' ..(r_1,k',k,r_2)..  ..(r_1,k',r_2)..
    /// during the apply a rank reduction will be performed
    /// 2*ndim allocates are needed
    template <class T, class Q>
    TensorTrain<TENSOR_RESULT_TYPE(T,Q)> apply(const TensorTrain<T>& op,
            const TensorTrain<Q>& t, const double thresh) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;
        MADNESS_ASSERT(op.ndim()==t.ndim());

        const long nd=t.ndim();

        // need to add special cases for low dimensions
        MADNESS_ASSERT(nd>2);

        std::vector<Tensor<resultT> > B(t.ndim());  // will be the result cores

        // set up scratch tensors
        long maxk=0;           // max dimension of the tensor t
        long maxr_t=1;         // max rank of the input tensor t
        long maxr_op=1;        // max rank of the operator op
        for (int i=0; i<nd-1; ++i) {
            maxk=std::max(maxk,t.dim(i));
            maxr_t=std::max(t.ranks(i),maxr_t);
            maxr_op=std::max(op.ranks(i),maxr_op);
        }

        long maxr_r=1;           // max rank of the result tensor
        for (int i=0, j=nd-1; i<j; ++i, --j) {
            maxr_r*=t.dim(i);
        }

        long maxn=0, maxm=0;
//        {
            long R11=t.dim(0);         // max final ranks
            long rR=R11*t.ranks(1);      // R1 r2
            long maxR=R11, maxr=t.ranks(1);            //
            for (int i=1; i<nd-2; ++i) {
                long k=t.dim(i);
                R11*=k;     // max final ranks
                R11=std::min(R11,maxr_r);
                long r2=t.ranks(i+1);   // initial ranks
                if (rR<R11*r2) {
                    maxR=std::max(maxR,R11);
                    maxr=std::max(maxr,r2);
                    rR=R11*r2;
                }
            }
            // max matrix dimensions to be svd'ed
            maxn=maxR*maxk;
            maxm=maxr_op*maxr;
//        }


        if (maxm*maxn>5e7) {
            print("huge scratch spaces!! ",maxn*maxm/1024/1024,"MByte");
        }
        long lscr=std::max(3*std::min(maxm,maxn)+std::max(maxm,maxn),
                5*std::min(maxm,maxn)) + maxn*maxm;
        Tensor<resultT> scr(lscr);   // scratch
        Tensor< typename Tensor<T>::scalar_type > s(std::min(maxn,maxm));

        // scratch space for contractions
        Tensor<resultT> scr3(2*maxn*maxm);
        Tensor<resultT> scr1=scr3(Slice(0,maxn*maxm-1));
        Tensor<resultT> scr2=scr3(Slice(maxn*maxm,-1));


        // contract first core
        const long r0=t.ranks(0l);
        const long q0=op.get_core(0).dim(2);
        const long k0=t.dim(0);
        inner_result(op.get_core(0),t.get_core(0),0,0,scr2);
        Tensor<resultT> AC=scr2(Slice(0,r0*q0*k0-1)).reshape(k0,r0*q0);

        // SVD on first core, skip small singular values
        long R=rank_revealing_decompose(AC,B[0],thresh,s,scr1);
        if (R==0) return TensorTrain<resultT>(t.dims());    // fast return for zero ranks
        B[0]=B[0].reshape(k0,R);

        // AC has dimensions R1,(q1,r1)
        Tensor<resultT> VT=AC.reshape(R,q0,r0);

        // loop over all dimensions 1,..,nd-1
        for (int d=1; d<nd; ++d) {

            // alias
            Tensor<T> C=t.get_core(d);
            Tensor<T> A=op.get_core(d);

            // alias dimensions
            long R1=VT.dim(0);        // left rank of the result tensor
            long q1=A.dim(0);         // left rank of the operator
            long k=C.dim(1);          // true for d>0
            long r2=1;                // right rank of the input tensor, true for d=nd-1
            long q2=1;                // right rank of the operator, true for d=nd-1
            if (d<nd-1) {
                r2=C.dim(2);          // right rank of the input tensor, true for d<nd-1
                q2=A.dim(3);          // right rank of the operator
            }

            // contract VT into the next core
            Tensor<resultT> VC=scr1(Slice(0,R1*q1*k*r2-1));
            VC(Slice(0,R1*q1*k*r2-1))=resultT(0.0);         // zero out result tensor
            inner_result(VT,C,-1,0,VC);          // VT(R1,q1,r1) * C(r1,k',r2)

            // contract A into VC (R,q,k,r2)
            VC=VC.reshape(R1,q1*k,r2);           // VC(R,(q,k),r2)
            A=A.fusedim(0);                      // A((q1,k),k,q2);
            Tensor<resultT> AVC=scr2(Slice(0,R1*q2*k*r2-1));
            AVC(Slice(0,R1*q2*k*r2-1))=resultT(0.0);        // zero out result tensor
            inner_result(VC,A,1,0,AVC);          // AVC(R1,r2,k,q2)
            AVC=AVC.reshape(R1,r2,k,q2);

            // AVC is the final core if we have reached the end of the tensor train
            if (d==nd-1) {
                B[d]=copy(AVC.reshape(R1,k));
                break;
            }

            // SVD on current core
            AVC=copy(AVC.cycledim(2,1,3));       // AVC(R,k,q2,r2); deep copy necessary

            MADNESS_ASSERT(AVC.dim(0)==R1);
            MADNESS_ASSERT(AVC.dim(1)==k);
            MADNESS_ASSERT(AVC.dim(2)==q2);

            AVC=AVC.reshape(R1*k,q2*r2);
            long R2=rank_revealing_decompose(AVC,B[d],thresh,s,scr1);
            if (R2==0) return TensorTrain<resultT>(t.dims());    // fast return for zero ranks
            B[d]=B[d].reshape(R1,k,R2);
            VT=AVC.reshape(R2,q2,r2);
        }

        TensorTrain<T> result(B);
        return result;

    }


    /// compute the n-D identity operator with k elements per dimension
    template<typename T>
    TensorTrain<T> tt_identity(const long ndim, const long k) {
        Tensor<T> id(k,k);
        for (int i=0; i<k; ++i) id(i,i)=1.0;
        id=id.reshape(1,k,k,1);
        std::vector<Tensor<T> > cores(ndim,id);
        TensorTrain<T> result(cores);
        if (ndim>1) {
            result.get_core(0)=result.get_core(0).reshape(k,k,1);
            result.get_core(ndim-1)=result.get_core(ndim-1).reshape(1,k,k);
        } else {
            result.get_core(0)=result.get_core(0).reshape(k,k);
        }
        return result;
    }


    /// computes the outer product of two tensors

    /// result(i,j,...,p,q,...) = left(i,k,...)*right(p,q,...)
    /// @result Returns a new tensor train
    template <class T, class Q>
    TensorTrain<TENSOR_RESULT_TYPE(T,Q)> outer(const TensorTrain<T>& t1,
            const TensorTrain<Q>& t2) {

        typedef TENSOR_RESULT_TYPE(T,Q) resultT;


        // fast return if possible
        if (t1.zero_rank or t2.zero_rank) {
            // compute new dimensions
            std::vector<long> dims(t1.ndim()+t2.ndim());
            for (int i=0; i<t1.ndim(); ++i) dims[i]=t1.dim(i);
            for (int i=0; i<t2.ndim(); ++i) dims[t1.ndim()+i]=t2.dim(i);

            return TensorTrain<resultT>(dims);
        }

        TensorTrain<resultT> result;
        for (int i=0; i<t1.ndim(); ++i) result.core.push_back(copy(t1.core[i]));
        for (int i=0; i<t2.ndim(); ++i) result.core.push_back(copy(t2.core[i]));

        // reshape the new interior cores
        long core_dim=t1.core.back().ndim();       // 2 for tensors, 3 for operators
        long k1=t1.core.back().dim(core_dim-1);    // (r,k) for tensors, (r,k',k) for operators
        long k2=t2.core.front().dim(0);
        result.core[t1.ndim()-1]=result.core[t1.ndim()-1].splitdim(core_dim-1,k1,1);
        result.core[t1.ndim()]=result.core[t1.ndim()].splitdim(0,1,k2);
        result.zero_rank=false;

        return result;

    }

}

#endif /* TENSORTRAIN_H_ */
