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
#ifndef MADNESS_MRA_VMRA_H__INCLUDED
#define MADNESS_MRA_VMRA_H__INCLUDED

/*!
	\file vmra.h
	\brief Defines operations on vectors of Functions
	\ingroup mra

	This file defines a number of operations on vectors of functions.
	Assume v is a vector of NDIM-D functions of a certain type.


	Operations on array of functions

	*) copying: deep copying of vectors of functions to vector of functions
	\code
	vector2 = copy(world, vector1,fence);
	\endcode

	*) compress: convert multiwavelet representation to legendre representation
	\code
	compress(world, vector, fence);
	\endcode

	*) reconstruct: convert representation to multiwavelets
	\code
	reconstruct(world, vector, fence);
	\endcode

	*) nonstandard: convert to non-standard form
	\code
	nonstandard(world, v, fence);
	\endcode

	*) standard: convert to standard form
	\code
	standard(world, v, fence);
	\endcode

	*) truncate: truncating vectors of functions to desired precision
	\code
	truncate(world, v, tolerance, fence);
	\endcode


	*) zero function: create a vector of zero functions of length n
	\code
	v=zero(world, n);
	\endcode

	*) transform: transform a representation from one basis to another
	\code
	transform(world, vector, tensor, tolerance, fence )
	\endcode

	Setting thresh-hold for precision

	*) set_thresh: setting a finite thresh-hold for a vector of functions
	\code
	void set_thresh(World& world, std::vector< Function<T,NDIM> >& v, double thresh, bool fence=true);
	\endcode

	Arithmetic Operations on arrays of functions

	*) conjugation: conjugate a vector of complex functions

	*) add
	*) sub
	*) mul
	   - mul_sparse
	*) square
	*) gaxpy
	*) apply

	Norms, inner-products, blas-1 like operations on vectors of functions

	*) inner
	*) matrix_inner
	*) norm_tree
	*) normalize
	*) norm2
	    - norm2s
	*) scale(world, v, alpha);




*/

#include <madness/mra/mra.h>
#include <madness/mra/derivative.h>
#include <madness/tensor/distributed_matrix.h>
#include <cstdio>

namespace madness {



    /// Compress a vector of functions
    template <typename T, std::size_t NDIM>
    void compress(World& world,
                  const std::vector< Function<T,NDIM> >& v,
                  bool fence=true) {

        PROFILE_BLOCK(Vcompress);
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); ++i) {
            if (!v[i].is_compressed()) {
                v[i].compress(false);
                must_fence = true;
            }
        }

        if (fence && must_fence) world.gop.fence();
    }


    /// Reconstruct a vector of functions
    template <typename T, std::size_t NDIM>
    void reconstruct(World& world,
                     const std::vector< Function<T,NDIM> >& v,
                     bool fence=true) {
        PROFILE_BLOCK(Vreconstruct);
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); ++i) {
            if (v[i].is_compressed()) {
                v[i].reconstruct(false);
                must_fence = true;
            }
        }

        if (fence && must_fence) world.gop.fence();
    }

    /// refine the functions according to the autorefine criteria
    template <typename T, std::size_t NDIM>
    void refine(World& world, const std::vector<Function<T,NDIM> >& vf,
            bool fence=true) {
        for (const auto& f : vf) f.refine(false);
        if (fence) world.gop.fence();
    }

    /// refine all functions to a common (finest) level

    /// if functions are not initialized (impl==NULL) they are ignored
    template <typename T, std::size_t NDIM>
    void refine_to_common_level(World& world, std::vector<Function<T,NDIM> >& vf,
            bool fence=true) {

        reconstruct(world,vf);
        Key<NDIM> key0(0, Vector<Translation, NDIM> (0));
        std::vector<FunctionImpl<T,NDIM>*> v_ptr;

        // push initialized function pointers into the vector v_ptr
        for (unsigned int i=0; i<vf.size(); ++i) {
            if (vf[i].is_initialized()) v_ptr.push_back(vf[i].get_impl().get());
        }

        // sort and remove duplicates to not confuse the refining function
        std::sort(v_ptr.begin(),v_ptr.end());
        typename std::vector<FunctionImpl<T, NDIM>*>::iterator it;
        it = std::unique(v_ptr.begin(), v_ptr.end());
        v_ptr.resize( std::distance(v_ptr.begin(),it) );

        std::vector< Tensor<T> > c(v_ptr.size());
        v_ptr[0]->refine_to_common_level(v_ptr, c, key0);
        if (fence) v_ptr[0]->world.gop.fence();
        if (VERIFY_TREE)
            for (unsigned int i=0; i<vf.size(); i++) vf[i].verify_tree();
    }

    /// Generates non-standard form of a vector of functions
    template <typename T, std::size_t NDIM>
    void nonstandard(World& world,
                     std::vector< Function<T,NDIM> >& v,
                     bool fence=true) {
        PROFILE_BLOCK(Vnonstandard);
        reconstruct(world, v);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].nonstandard(false,false);
        }
        if (fence) world.gop.fence();
    }


    /// Generates standard form of a vector of functions
    template <typename T, std::size_t NDIM>
    void standard(World& world,
                  std::vector< Function<T,NDIM> >& v,
                  bool fence=true) {
        PROFILE_BLOCK(Vstandard);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].standard(false);
        }
        if (fence) world.gop.fence();
    }


    /// Truncates a vector of functions
    template <typename T, std::size_t NDIM>
    void truncate(World& world,
                  std::vector< Function<T,NDIM> >& v,
                  double tol=0.0,
                  bool fence=true) {
        PROFILE_BLOCK(Vtruncate);

        compress(world, v);

        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].truncate(tol, false);
        }

        if (fence) world.gop.fence();
    }

    /// Truncates a vector of functions

    /// @return the truncated vector for chaining
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> > truncate(std::vector< Function<T,NDIM> > v,
                  double tol=0.0, bool fence=true) {
        if (v.size()>0) truncate(v[0].world(),v,tol,fence);
        return v;
    }


    /// Applies a derivative operator to a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    apply(World& world,
          const Derivative<T,NDIM>& D,
          const std::vector< Function<T,NDIM> >& v,
          bool fence=true)
    {
        reconstruct(world, v);
        std::vector< Function<T,NDIM> > df(v.size());
        for (unsigned int i=0; i<v.size(); ++i) {
            df[i] = D(v[i],false);
        }
        if (fence) world.gop.fence();
        return df;
    }

    /// Generates a vector of zero functions (reconstructed)
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    zero_functions(World& world, int n, bool fence=true) {
        PROFILE_BLOCK(Vzero_functions);
        std::vector< Function<T,NDIM> > r(n);
        for (int i=0; i<n; ++i)
  	    r[i] = Function<T,NDIM>(FunctionFactory<T,NDIM>(world).fence(false));

	if (n && fence) world.gop.fence();

        return r;
    }

    /// Generates a vector of zero functions (compressed)
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    zero_functions_compressed(World& world, int n, bool fence=true) {
        PROFILE_BLOCK(Vzero_functions);
        std::vector< Function<T,NDIM> > r(n);
        for (int i=0; i<n; ++i)
  	    r[i] = Function<T,NDIM>(FunctionFactory<T,NDIM>(world).fence(false).compressed(true).initial_level(1));

	if (n && fence) world.gop.fence();

        return r;
    }

    /// symmetric orthonormalization (see e.g. Szabo/Ostlund)
    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_symmetric(
    		const std::vector<Function<T,NDIM> >& v,
			  const Tensor<T>& ovlp) {
    	if(v.empty()) return v;

    	Tensor<T> U;
    	Tensor< typename Tensor<T>::scalar_type > s;
    	syev(ovlp,U,s);

    	// transform s to s^{-1}
    	for(size_t i=0;i<v.size();++i) s(i)=1.0/(sqrt(s(i)));

    	// save Ut before U gets modified with s^{-1}
    	const Tensor<T> Ut=transpose(U);
    	for(size_t i=0;i<v.size();++i){
    		for(size_t j=0;j<v.size();++j){
    			U(i,j)=U(i,j)*s(j);
    		}
    	}

    	Tensor<T> X=inner(U,Ut,1,0);

    	World& world=v.front().world();
    	return transform(world,v,X);

    }
    /// convenience routine for symmetric orthonormalization (see e.g. Szabo/Ostlund)
    /// overlap matrix is calculated
    /// @param[in] the vector to orthonormalize
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_symmetric(const std::vector<Function<T,NDIM> >& v){
    	if(v.empty()) return v;


    	World& world=v.front().world();
    	Tensor<T> ovlp = matrix_inner(world, v, v);

    	return orthonormalize_symmetric(v,ovlp);
    }

    /// canonical orthonormalization (see e.g. Szabo/Ostlund)
    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix
    /// @param[in]	lindep	linear dependency threshold relative to largest eigenvalue
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_canonical(
    		const std::vector<Function<T,NDIM> >& v,
			const Tensor<T>& ovlp,
			double lindep) {

    	if(v.empty()) return v;

    	Tensor<T> U;
    	Tensor< typename Tensor<T>::scalar_type > s;
    	syev(ovlp,U,s);
    	lindep*=s(s.size()-1);	// eigenvalues are in ascending order

    	// transform s to s^{-1}
    	int rank=0,lo=0;
    	Tensor< typename Tensor<T>::scalar_type > sqrts(v.size());
    	for(size_t i=0;i<v.size();++i) {
    		if (s(i)>lindep) {
    			sqrts(i)=1.0/(sqrt(s(i)));
        		rank++;
    		} else {
    			sqrts(i)=0.0;
    			lo++;
    		}
    	}
    	MADNESS_ASSERT(lo+rank==v.size());

    	for(size_t i=0;i<v.size();++i){
    		for(size_t j=0;j<v.size();++j){
    			U(i,j)=U(i,j)*(sqrts(j));
    		}
    	}
    	Tensor<double> X=U(_,Slice(lo,-1));

    	World& world=v.front().world();
    	return transform(world,v,X);

    }

    /// convenience routine for canonical routine for symmetric orthonormalization (see e.g. Szabo/Ostlund)
    /// overlap matrix is calculated
    /// @param[in] the vector to orthonormalize
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_canonical(const std::vector<Function<T,NDIM> >& v,
    		const double lindep){
    	if(v.empty()) return v;

    	World& world=v.front().world();
    	Tensor<T> ovlp = matrix_inner(world, v, v);

    	return orthonormalize_canonical(v,ovlp,lindep);
    }

    /// cholesky orthonormalization without pivoting
    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix, destroyed on return!
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_cd(
    		const std::vector<Function<T,NDIM> >& v,
			Tensor<T>& ovlp) {

    	if (v.empty()) return v;

    	cholesky(ovlp); // destroys ovlp and gives back Upper ∆ Matrix from CD

    	Tensor<T> L = transpose(ovlp);
    	Tensor<T> Linv = inverse(L);
    	Tensor<T> U = transpose(Linv);

    	World& world=v.front().world();
    	return transform(world, v, U);

    }

    /// convenience routine for cholesky orthonormalization without pivoting
    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_cd(const std::vector<Function<T,NDIM> >& v){
    	if(v.empty()) return v;

    	World& world=v.front().world();
    	Tensor<T> ovlp = matrix_inner(world, v, v);

    	return orthonormalize_cd(v,ovlp);
    }

    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix, will be destroyed on return!
    /// @param[in] tolerance for numerical rank reduction
    /// @param[out] pivoting vector, no allocation on input needed
    /// @param[out] rank
    /// @return orthonrormalized vector (may or may not be truncated)
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_rrcd(
    		const std::vector<Function<T,NDIM> >& v,
			Tensor<T>& ovlp,
			const double tol,
			Tensor<integer>& piv,
			int& rank) {

    	if (v.empty()) {
    		return v;
    	}

    	rr_cholesky(ovlp,tol,piv,rank); // destroys ovlp and gives back Upper ∆ Matrix from CCD

    	// rearrange and truncate the functions according to the pivoting of the rr_cholesky
    	std::vector<Function<T,NDIM> > pv(rank);
    	for(integer i=0;i<rank;++i){
    		pv[i]=v[piv[i]];
    	}
    	ovlp=ovlp(Slice(0,rank-1),Slice(0,rank-1));

    	Tensor<T> L = transpose(ovlp);
    	Tensor<T> Linv = inverse(L);
    	Tensor<T> U = transpose(Linv);

    	World& world=v.front().world();
    	return transform(world, pv, U);
    }

    /// convenience routine for orthonromalize_cholesky: orthonromalize_cholesky without information on pivoting and rank
    /// @param[in] the vector to orthonormalize
    /// @param[in] overlap matrix
    /// @param[in] tolerance for numerical rank reduction
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_rrcd(const std::vector<Function<T,NDIM> >& v, Tensor<T> ovlp , const double tol) {
    	Tensor<integer> piv;
    	int rank;
    	return orthonormalize_rrcd(v,ovlp,tol,piv,rank);
    }

    /// convenience routine for orthonromalize_cholesky: computes the overlap matrix and then calls orthonromalize_cholesky
    /// @param[in] the vector to orthonormalize
    /// @param[in] tolerance for numerical rank reduction
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > orthonormalize_rrcd(const std::vector<Function<T,NDIM> >& v, const double tol) {
    	if (v.empty()) {
    		return v;
    	}
    	// compute overlap
    	World& world=v.front().world();
    	Tensor<T> ovlp = matrix_inner(world, v, v);
    	return orthonormalize_rrcd(v,ovlp,tol);
    }

    /// combine two vectors
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > append(const std::vector<Function<T,NDIM> > & lhs, const std::vector<Function<T,NDIM> > & rhs){
    	std::vector<Function<T,NDIM> >  v=lhs;
    	for (std::size_t i = 0; i < rhs.size(); ++i) v.push_back(rhs[i]);
    	return v;
    }

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > append(const std::vector< std::vector<Function<T,NDIM> > > vv){
    	std::vector<Function<T,NDIM> >result;
    	for(const auto& x:vv) result=append(result,x);
    	return result;
    }

    /// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]

    /// Uses sparsity in the transformation matrix --- set small elements to
    /// zero to take advantage of this.
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> >
    transform(World& world,
              const std::vector< Function<T,NDIM> >& v,
              const Tensor<R>& c,
              bool fence=true) {

        PROFILE_BLOCK(Vtransformsp);
        typedef TENSOR_RESULT_TYPE(T,R) resultT;
        int n = v.size();  // n is the old dimension
        int m = c.dim(1);  // m is the new dimension
        MADNESS_ASSERT(n==c.dim(0));

        std::vector< Function<resultT,NDIM> > vc = zero_functions_compressed<resultT,NDIM>(world, m);
        compress(world, v);

        for (int i=0; i<m; ++i) {
            for (int j=0; j<n; ++j) {
                if (c(j,i) != R(0.0)) vc[i].gaxpy(1.0,v[j],c(j,i),false);
            }
        }

        if (fence) world.gop.fence();
        return vc;
    }


    template <typename L, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(L,R),NDIM> >
    transform(World& world,  const std::vector< Function<L,NDIM> >& v,
            const Tensor<R>& c, double tol, bool fence) {
        PROFILE_BLOCK(Vtransform);
        MADNESS_ASSERT(v.size() == (unsigned int)(c.dim(0)));

        std::vector< Function<TENSOR_RESULT_TYPE(L,R),NDIM> > vresult
            = zero_functions_compressed<TENSOR_RESULT_TYPE(L,R),NDIM>(world, c.dim(1));

        compress(world, v, true);
        vresult[0].vtransform(v, c, vresult, tol, fence);
        return vresult;
    }

    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> >
    transform(World& world,
              const std::vector< Function<T,NDIM> >& v,
              const DistributedMatrix<R>& c,
              bool fence=true) {
        PROFILE_FUNC;

        typedef TENSOR_RESULT_TYPE(T,R) resultT;
        long n = v.size();    // n is the old dimension
        long m = c.rowdim();  // m is the new dimension
        MADNESS_ASSERT(n==c.coldim());

        // new(i) = sum(j) old(j) c(j,i)

        Tensor<T> tmp(n,m);
        c.copy_to_replicated(tmp); // for debugging
        tmp = transpose(tmp);

        std::vector< Function<resultT,NDIM> > vc = zero_functions_compressed<resultT,NDIM>(world, m);
        compress(world, v);

        for (int i=0; i<m; ++i) {
            for (int j=0; j<n; ++j) {
                if (tmp(j,i) != R(0.0)) vc[i].gaxpy(1.0,v[j],tmp(j,i),false);
            }
        }

        if (fence) world.gop.fence();
        return vc;
    }


    /// Scales inplace a vector of functions by distinct values
    template <typename T, typename Q, std::size_t NDIM>
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v,
               const std::vector<Q>& factors,
               bool fence=true) {
        PROFILE_BLOCK(Vscale);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(factors[i],false);
        if (fence) world.gop.fence();
    }

    /// Scales inplace a vector of functions by the same
    template <typename T, typename Q, std::size_t NDIM>
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v,
	       const Q factor,
               bool fence=true) {
        PROFILE_BLOCK(Vscale);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(factor,false);
        if (fence) world.gop.fence();
    }

    /// Computes the 2-norms of a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector<double> norm2s(World& world,
                              const std::vector< Function<T,NDIM> >& v) {
        PROFILE_BLOCK(Vnorm2);
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = sqrt(norms[i]);
        world.gop.fence();
        return norms;
    }

    /// Computes the 2-norm of a vector of functions
    template <typename T, std::size_t NDIM>
    double norm2(World& world,
                              const std::vector< Function<T,NDIM> >& v) {
        PROFILE_BLOCK(Vnorm2);
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=1; i<v.size(); ++i) norms[0] += norms[i];
        world.gop.fence();
        return sqrt(norms[0]);
    }

    inline double conj(double x) {
        return x;
    }

    inline double conj(float x) {
        return x;
    }

// !!! FIXME: this task is broken because FunctionImpl::inner_local forces a
// future on return from WorldTaskQueue::reduce, which will causes a deadlock if
// run inside a task. This behavior must be changed before this task can be used
// again.
//
//    template <typename T, typename R, std::size_t NDIM>
//    struct MatrixInnerTask : public TaskInterface {
//        Tensor<TENSOR_RESULT_TYPE(T,R)> result; // Must be a copy
//        const Function<T,NDIM>& f;
//        const std::vector< Function<R,NDIM> >& g;
//        long jtop;
//
//        MatrixInnerTask(const Tensor<TENSOR_RESULT_TYPE(T,R)>& result,
//                        const Function<T,NDIM>& f,
//                        const std::vector< Function<R,NDIM> >& g,
//                        long jtop)
//                : result(result), f(f), g(g), jtop(jtop) {}
//
//        void run(World& world) {
//            for (long j=0; j<jtop; ++j) {
//                result(j) = f.inner_local(g[j]);
//            }
//        }
//
//    private:
//        /// Get the task id
//
//        /// \param id The id to set for this task
//        virtual void get_id(std::pair<void*,unsigned short>& id) const {
//            PoolTaskInterface::make_id(id, *this);
//        }
//    }; // struct MatrixInnerTask



    template <typename T, std::size_t NDIM>
    DistributedMatrix<T> matrix_inner(const DistributedMatrixDistribution& d,
                                      const std::vector< Function<T,NDIM> >& f,
                                      const std::vector< Function<T,NDIM> >& g,
                                      bool sym=false)
    {
        PROFILE_FUNC;
        DistributedMatrix<T> A(d);
        const int64_t n = A.coldim();
        const int64_t m = A.rowdim();
        MADNESS_ASSERT(int64_t(f.size()) == n && int64_t(g.size()) == m);

        // Assume we can always create an ichunk*jchunk matrix locally
        const int ichunk = 1000;
        const int jchunk = 1000; // 1000*1000*8 = 8 MBytes
        for (int64_t ilo=0; ilo<n; ilo+=ichunk) {
            int64_t ihi = std::min(ilo + ichunk, n);
            std::vector< Function<T,NDIM> > ivec(f.begin()+ilo, f.begin()+ihi);
            for (int64_t jlo=0; jlo<m; jlo+=jchunk) {
                int64_t jhi = std::min(jlo + jchunk, m);
                std::vector< Function<T,NDIM> > jvec(g.begin()+jlo, g.begin()+jhi);

                Tensor<T> P = matrix_inner(A.get_world(),ivec,jvec);
                A.copy_from_replicated_patch(ilo, ihi-1, jlo, jhi-1, P);
            }
        }
        return A;
    }

    /// Computes the matrix inner product of two function vectors - q(i,j) = inner(f[i],g[j])

    /// For complex types symmetric is interpreted as Hermitian.
    ///
    /// The current parallel loop is non-optimal but functional.
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > matrix_inner(World& world,
                                                   const std::vector< Function<T,NDIM> >& f,
                                                   const std::vector< Function<R,NDIM> >& g,
                                                   bool sym=false) 
    {
        world.gop.fence();
        compress(world, f);
        if ((void*)(&f) != (void*)(&g)) compress(world, g);

        std::vector<const FunctionImpl<T,NDIM>*> left(f.size());
        std::vector<const FunctionImpl<R,NDIM>*> right(g.size());
        for (unsigned int i=0; i<f.size(); i++) left[i] = f[i].get_impl().get();
        for (unsigned int i=0; i<g.size(); i++) right[i]= g[i].get_impl().get();

        Tensor< TENSOR_RESULT_TYPE(T,R) > r= FunctionImpl<T,NDIM>::inner_local(left, right, sym);

        world.gop.fence();
        world.gop.sum(r.ptr(),f.size()*g.size());

        return r;
    }

    /// Computes the matrix inner product of two function vectors - q(i,j) = inner(f[i],g[j])

    /// For complex types symmetric is interpreted as Hermitian.
    ///
    /// The current parallel loop is non-optimal but functional.
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > matrix_inner_old(World& world,
            const std::vector< Function<T,NDIM> >& f,
            const std::vector< Function<R,NDIM> >& g,
            bool sym=false) {
        PROFILE_BLOCK(Vmatrix_inner);
        long n=f.size(), m=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n,m);
        if (sym) MADNESS_ASSERT(n==m);

        world.gop.fence();
        compress(world, f);
        if ((void*)(&f) != (void*)(&g)) compress(world, g);

        for (long i=0; i<n; ++i) {
            long jtop = m;
            if (sym) jtop = i+1;
            for (long j=0; j<jtop; ++j) {
                r(i,j) = f[i].inner_local(g[j]);
                if (sym) r(j,i) = conj(r(i,j));
            }
         }
        
//        for (long i=n-1; i>=0; --i) {
//            long jtop = m;
//            if (sym) jtop = i+1;
//            world.taskq.add(new MatrixInnerTask<T,R,NDIM>(r(i,_), f[i], g, jtop));
//        }
        world.gop.fence();
        world.gop.sum(r.ptr(),n*m);

//        if (sym) {
//            for (int i=0; i<n; ++i) {
//                for (int j=0; j<i; ++j) {
//                    r(j,i) = conj(r(i,j));
//                }
//            }
//        }
        return r;
    }

    /// Computes the element-wise inner product of two function vectors - q(i) = inner(f[i],g[i])
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const std::vector< Function<T,NDIM> >& f,
                                            const std::vector< Function<R,NDIM> >& g) {
        PROFILE_BLOCK(Vinnervv);
        long n=f.size(), m=g.size();
        MADNESS_ASSERT(n==m);
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n);

        compress(world, f);
        compress(world, g);

        for (long i=0; i<n; ++i) {
            r(i) = f[i].inner_local(g[i]);
        }

        world.taskq.fence();
        world.gop.sum(r.ptr(),n);
        world.gop.fence();
        return r;
    }


    /// Computes the inner product of a function with a function vector - q(i) = inner(f,g[i])
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const Function<T,NDIM>& f,
                                            const std::vector< Function<R,NDIM> >& g) {
        PROFILE_BLOCK(Vinner);
        long n=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n);

        f.compress();
        compress(world, g);

        for (long i=0; i<n; ++i) {
            r(i) = f.inner_local(g[i]);
        }

        world.taskq.fence();
        world.gop.sum(r.ptr(),n);
        world.gop.fence();
        return r;
    }

    /// inner function with right signature for the nonlinear sovler
    /// this is needed for the KAIN solvers and other functions
    template <typename T, typename R, std::size_t NDIM>
    double inner( const std::vector< Function<T,NDIM> >& f,
	                                            const std::vector< Function<R,NDIM> >& g){
      MADNESS_ASSERT(f.size()==g.size());
      if(f.empty()) return 0.0;
      else return inner(f[0].world(),f,g).sum();
    }


    /// Multiplies a function against a vector of functions --- q[i] = a * v[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul(World& world,
        const Function<T,NDIM>& a,
        const std::vector< Function<R,NDIM> >& v,
        bool fence=true) {
        PROFILE_BLOCK(Vmul);
        a.reconstruct(false);
        reconstruct(world, v, false);
        world.gop.fence();
        return vmulXX(a, v, 0.0, fence);
    }

    /// Multiplies a function against a vector of functions using sparsity of a and v[i] --- q[i] = a * v[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul_sparse(World& world,
               const Function<T,NDIM>& a,
               const std::vector< Function<R,NDIM> >& v,
               double tol,
               bool fence=true) {
        PROFILE_BLOCK(Vmulsp);
        a.reconstruct(false);
        reconstruct(world, v, false);
        world.gop.fence();
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].norm_tree(false);
        }
        a.norm_tree();
        return vmulXX(a, v, tol, fence);
    }

    /// Makes the norm tree for all functions in a vector
    template <typename T, std::size_t NDIM>
    void norm_tree(World& world,
                   const std::vector< Function<T,NDIM> >& v,
                   bool fence=true)
    {
        PROFILE_BLOCK(Vnorm_tree);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].norm_tree(false);
        }
        if (fence) world.gop.fence();
    }

    /// Multiplies two vectors of functions q[i] = a[i] * b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vmulvv);
        reconstruct(world, a, true);
        if (&a != &b) reconstruct(world, b, true);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > q(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            q[i] = mul(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return q;
    }


    /// Computes the square of a vector of functions --- q[i] = v[i]**2
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    square(World& world,
           const std::vector< Function<T,NDIM> >& v,
           bool fence=true) {
        return mul<T,T,NDIM>(world, v, v, fence);
//         std::vector< Function<T,NDIM> > vsq(v.size());
//         for (unsigned int i=0; i<v.size(); ++i) {
//             vsq[i] = square(v[i], false);
//         }
//         if (fence) world.gop.fence();
//         return vsq;
    }


    /// Sets the threshold in a vector of functions
    template <typename T, std::size_t NDIM>
    void set_thresh(World& world, std::vector< Function<T,NDIM> >& v, double thresh, bool fence=true) {
        for (unsigned int j=0; j<v.size(); ++j) {
            v[j].set_thresh(thresh,false);
        }
        if (fence) world.gop.fence();
    }

    /// Returns the complex conjugate of the vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    conj(World& world,
         const std::vector< Function<T,NDIM> >& v,
         bool fence=true) {
        PROFILE_BLOCK(Vconj);
        std::vector< Function<T,NDIM> > r = copy(world, v); // Currently don't have oop conj
        for (unsigned int i=0; i<v.size(); ++i) {
            r[i].conj(false);
        }
        if (fence) world.gop.fence();
        return r;
    }


    /// Returns a deep copy of a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    copy(World& world,
         const std::vector< Function<T,NDIM> >& v,
         bool fence=true) {
        PROFILE_BLOCK(Vcopy);
        std::vector< Function<T,NDIM> > r(v.size());
        for (unsigned int i=0; i<v.size(); ++i) {
            r[i] = copy(v[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

   /// Returns a vector of deep copies of of a function
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    copy(World& world,
         const Function<T,NDIM>& v,
         const unsigned int n,
         bool fence=true) {
        PROFILE_BLOCK(Vcopy1);
        std::vector< Function<T,NDIM> > r(n);
        for (unsigned int i=0; i<n; ++i) {
            r[i] = copy(v, false);
        }
        if (fence) world.gop.fence();
        return r;
    }

    /// Returns new vector of functions --- q[i] = a[i] + b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vadd);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            r[i] = add(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

     /// Returns new vector of functions --- q[i] = a + b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const Function<T,NDIM> & a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vadd1);
        a.compress();
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(b.size());
        for (unsigned int i=0; i<b.size(); ++i) {
            r[i] = add(a, b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }
    template <typename T, typename R, std::size_t NDIM>
    inline std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const std::vector< Function<R,NDIM> >& b,
        const Function<T,NDIM> & a,
        bool fence=true) {
        return add(world, a, b, fence);
    }

    /// Returns new vector of functions --- q[i] = a[i] - b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    sub(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vsub);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            r[i] = sub(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

    /// Returns new function --- q = sum_i f[i]
    template <typename T, std::size_t NDIM>
    Function<T, NDIM> sum(World& world, const std::vector<Function<T,NDIM> >& f,
        bool fence=true) {

        compress(world, f);
        Function<T,NDIM> r=real_factory_3d(world).compressed();

        for (unsigned int i=0; i<f.size(); ++i) r.gaxpy(1.0,f[i],1.0,false);
        if (fence) world.gop.fence();
        return r;
    }


    /// Multiplies and sums two vectors of functions r = \sum_i a[i] * b[i]
    template <typename T, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(T,R), NDIM>
    dot(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        return sum(world,mul(world,a,b,true),fence);
    }



    /// Generalized A*X+Y for vectors of functions ---- a[i] = alpha*a[i] + beta*b[i]
    template <typename T, typename Q, typename R, std::size_t NDIM>
    void gaxpy(World& world,
               Q alpha,
               std::vector< Function<T,NDIM> >& a,
               Q beta,
               const std::vector< Function<R,NDIM> >& b,
               bool fence=true) {
        PROFILE_BLOCK(Vgaxpy);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        for (unsigned int i=0; i<a.size(); ++i) {
            a[i].gaxpy(alpha, b[i], beta, false);
        }
        if (fence) world.gop.fence();
    }


    /// Applies a vector of operators to a vector of functions --- q[i] = apply(op[i],f[i])
    template <typename opT, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> >
    apply(World& world,
          const std::vector< std::shared_ptr<opT> >& op,
          const std::vector< Function<R,NDIM> > f) {

        PROFILE_BLOCK(Vapplyv);
        MADNESS_ASSERT(f.size()==op.size());

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); ++i) {
            MADNESS_ASSERT(not op[i]->is_slaterf12);
            result[i] = apply_only(*op[i], f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
        world.gop.fence();
        reconstruct(world, result);

        return result;
    }


    /// Applies an operator to a vector of functions --- q[i] = apply(op,f[i])
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    apply(World& world,
          const SeparatedConvolution<T,NDIM>& op,
          const std::vector< Function<R,NDIM> > f) {
        PROFILE_BLOCK(Vapply);

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); ++i) {
            result[i] = apply_only(op, f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
        reconstruct(world, result);

        if (op.is_slaterf12) {
        	MADNESS_ASSERT(not op.destructive());
            for (unsigned int i=0; i<f.size(); ++i) {
            	double trace=f[i].trace();
                result[i]=(result[i]-trace).scale(-0.5/op.mu());
            }
        }

        return result;
    }

    /// Normalizes a vector of functions --- v[i] = v[i].scale(1.0/v[i].norm2())
    template <typename T, std::size_t NDIM>
    void normalize(World& world, std::vector< Function<T,NDIM> >& v, bool fence=true) {
        PROFILE_BLOCK(Vnormalize);
        std::vector<double> nn = norm2s(world, v);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(1.0/nn[i],false);
        if (fence) world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    void print_size(World &world, const std::vector<Function<T,NDIM> > &v, const std::string &msg = "vectorfunction" ){
    	if(v.empty()){
    		if(world.rank()==0) std::cout << "print_size: " << msg << " is empty" << std::endl;
    	}else if(v.size()==1){
    		v.front().print_size(msg);
    	}else{
    		for(auto x:v){
    			x.print_size(msg);
    		}
    	}
    }

    // gives back the size in GB
    template <typename T, std::size_t NDIM>
    double get_size(World& world, const std::vector< Function<T,NDIM> >& v){

    	if (v.empty()) return 0.0;

    	const double d=sizeof(T);
        const double fac=1024*1024*1024;

        double size=0.0;
        for(unsigned int i=0;i<v.size();i++){
            if (v[i].is_initialized()) size+=v[i].size();
        }

        return size/fac*d;

    }

    // gives back the size in GB
    template <typename T, std::size_t NDIM>
    double get_size(const Function<T,NDIM> & f){
    	const double d=sizeof(T);
        const double fac=1024*1024*1024;
        double size=f.size();
        return size/fac*d;
    }

    /// apply op on the input vector yielding an output vector of functions

    /// @param[in]  op   the operator working on vin
    /// @param[in]  vin  vector of input Functions; needs to be refined to common level!
    /// @return vector of output Functions vout = op(vin)
    template <typename T, typename opT, std::size_t NDIM>
    std::vector<Function<T,NDIM> > multi_to_multi_op_values(const opT& op,
            const std::vector< Function<T,NDIM> >& vin,
            const bool fence=true) {
        MADNESS_ASSERT(vin.size()>0);
        MADNESS_ASSERT(vin[0].is_initialized()); // might be changed
        World& world=vin[0].world();
        Function<T,NDIM> dummy;
        dummy.set_impl(vin[0], false);
        std::vector<Function<T,NDIM> > vout=zero_functions<T,NDIM>(world, op.get_result_size());
        for (auto& out : vout) out.set_impl(vin[0],false);
        dummy.multi_to_multi_op_values(op, vin, vout, fence);
        return vout;
    }




    // convenience operators

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator+(const std::vector<Function<T,NDIM> >& lhs,
            const std::vector<Function<T,NDIM>>& rhs) {
        return add(lhs[0].world(),lhs,rhs);
    }

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator-(const std::vector<Function<T,NDIM> >& lhs,
            const std::vector<Function<T,NDIM> >& rhs) {
        return sub(lhs[0].world(),lhs,rhs);
    }

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator*(const double fac,
            const std::vector<Function<T,NDIM> >& rhs) {
        std::vector<Function<T,NDIM> > tmp=copy(rhs[0].world(),rhs);
        scale(tmp[0].world(),tmp,fac);
        return tmp;
    }

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator*(const std::vector<Function<T,NDIM> >& rhs,
            const double fac) {
        std::vector<Function<T,NDIM> > tmp=copy(rhs[0].world(),rhs);
        scale(tmp[0].world(),tmp,fac);
        return tmp;
    }

    /// multiply a vector of functions with a function: r[i] = v[i] * a
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator*(const Function<T,NDIM>& a,
            const std::vector<Function<T,NDIM> >& v) {
        return mul(v[0].world(),a,v,true);
    }


    /// multiply a vector of functions with a function: r[i] = a * v[i]
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator*(const std::vector<Function<T,NDIM> >& v,
            const Function<T,NDIM>& a) {
        return mul(v[0].world(),a,v,true);
    }


    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator+=(std::vector<Function<T,NDIM> >& rhs,
            const std::vector<Function<T,NDIM> >& lhs) {
        rhs=add(rhs[0].world(),rhs,lhs);
        return rhs;
    }

    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > operator-=(std::vector<Function<T,NDIM> >& rhs,
            const std::vector<Function<T,NDIM> >& lhs) {
        rhs=sub(rhs[0].world(),rhs,lhs);
        return rhs;
    }

    /// shorthand gradient operator

    /// returns the differentiated function f in all NDIM directions
    /// @param[in]  f       the function on which the grad operator works on
    /// @param[in]  refine  refinement before diff'ing makes the result more accurate
    /// @param[in]  fence   fence after completion; if reconstruction is needed always fence
    /// @return     the vector \frac{\partial}{\partial x_i} f
    template <typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > grad(const Function<T,NDIM>& f,
            bool refine=false, bool fence=true) {

        World& world=f.world();
        f.reconstruct();
        if (refine) f.refine();      // refine to make result more precise

        std::vector< std::shared_ptr< Derivative<T,NDIM> > > grad=
                gradient_operator<T,NDIM>(world);

        std::vector<Function<T,NDIM> > result(NDIM);
        for (int i=0; i<NDIM; ++i) result[i]=apply(*(grad[i]),f,false);
        if (fence) world.gop.fence();
        return result;
    }

    /// shorthand div operator

    /// returns the dot product of nabla with a vector f
    /// @param[in]  f       the vector of functions on which the div operator works on
    /// @param[in]  refine  refinement before diff'ing makes the result more accurate
    /// @param[in]  fence   fence after completion; currently always fences
    /// @return     the vector \frac{\partial}{\partial x_i} f
    /// TODO: add this to operator fusion
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> div(const std::vector<Function<T,NDIM> >& v,
            bool do_refine=false, bool fence=true) {

        MADNESS_ASSERT(v.size()>0);
        World& world=v[0].world();
        reconstruct(world,v);
        if (do_refine) refine(world,v);      // refine to make result more precise

        std::vector< std::shared_ptr< Derivative<T,NDIM> > > grad=
                gradient_operator<T,NDIM>(world);

        std::vector<Function<T,NDIM> > result(NDIM);
        for (int i=0; i<NDIM; ++i) result[i]=apply(*(grad[i]),v[i],false);
        world.gop.fence();
        return sum(world,result,fence);
    }

    /// load a vector of functions
    template<typename T, size_t NDIM>
    void load_function(World& world, std::vector<Function<T,NDIM> >& f,
            const std::string name) {
        if (world.rank()==0) print("loading vector of functions",name);
        archive::ParallelInputArchive ar(world, name.c_str(), 1);
        std::size_t fsize=0;
        ar & fsize;
        f.resize(fsize);
        for (std::size_t i=0; i<fsize; ++i) ar & f[i];
    }

    /// save a vector of functions
    template<typename T, size_t NDIM>
    void save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) {
        if (f.size()>0) {
            World& world=f.front().world();
            if (world.rank()==0) print("saving vector of functions",name);
            archive::ParallelOutputArchive ar(world, name.c_str(), 1);
            std::size_t fsize=f.size();
            ar & fsize;
            for (std::size_t i=0; i<fsize; ++i) ar & f[i];
        }
    }


}
#endif // MADNESS_MRA_VMRA_H__INCLUDED
