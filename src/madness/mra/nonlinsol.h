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

#ifndef MADNESS_EXAMPLES_NONLINSOL_H__INCLUDED
#define MADNESS_EXAMPLES_NONLINSOL_H__INCLUDED

/*!
  \file nonlinsol.h
  \brief Implementation of Krylov-subspace nonlinear equation solver
  \defgroup nonlinearsolve Simple Krylov-subspace nonlinear equation solver 
  \ingroup mra

  This class implements the solver described in 
  \verbatim
   R. J. Harrison, Krylov subspace accelerated inexact newton method for linear
   and nonlinear equations, J. Comput. Chem. 25 (2004), no. 3, 328-334.
  \endverbatim
 */

#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>

namespace madness {


	/// check for subspace linear dependency

	/// @param[in]     Q	the input matrix for KAIN
	/// @param[in,out] c	the coefficients for constructing the new solution
    /// @param[in]     rcondtol rcond less than this will cause the subspace to be shrunk due to linear dependence
    /// @param[in]     cabsmax  maximum element of c greater than this will cause the subspace to be shrunk due to linear dependence
	template<typename C>
	void check_linear_dependence(const Tensor<C>& Q, Tensor<C>& c, const double rcondtol, const double cabsmax,
			bool do_print=true) {
		double rcond = 1e-12;
		int m = c.dim(0);

		while(1){
			c = KAIN(Q, rcond);
			//if (world.rank() == 0) print("kain c:", c);
//			if(std::abs(c[m - 1]) < 3.0){
			if (c.absmax()<cabsmax) {
				break;
			} else  if(rcond < rcondtol){
				if (do_print) print("Increasing subspace singular value threshold ", c[m - 1], rcond);
				rcond *= 100;
			} else {
				if (do_print) print("Forcing full step due to subspace malfunction");
				c = 0.0;
				c[m - 1] = 1.0;
				break;
			}
		}
	}

	/// A simple Krylov-subspace nonlinear equation solver

    /// \ingroup nonlinearsolve
	template<size_t NDIM>
	class NonlinearSolverND {
		unsigned int maxsub; ///< Maximum size of subspace dimension

		std::vector<Function<double,NDIM> > ulist, rlist;
		real_tensor Q;

	public:
		void set_maxsub(const unsigned int &new_maxsub){
			maxsub = new_maxsub;
		}
		unsigned int get_maxsub()const{
			const unsigned int tmp = maxsub;
			return tmp;
		}
		bool do_print;

		NonlinearSolverND(unsigned int maxsub = 10) : maxsub(maxsub), do_print(false) {}

		/// Computes next trial solution vector

		/// You are responsible for performing step restriction or line search
		/// (not necessary for linear problems).
		///
		/// @param u Current solution vector
		/// @param r Corresponding residual
		/// @return Next trial solution vector
		/// @param[in]          rcondtol rcond less than this will cause the subspace to be shrunk due to linear dependence
		/// @param[in]          cabsmax  maximum element of c greater than this will cause the subspace to be shrunk due to linear dependence
		Function<double,NDIM> update(const Function<double,NDIM>& u, const Function<double,NDIM>& r,
				const double rcondtol=1e-8, const double cabsmax=1000.0) {
			if (maxsub==1) return u-r;
			int iter = ulist.size();
			ulist.push_back(u);
			rlist.push_back(r);

			// Solve subspace equations
			real_tensor Qnew(iter+1,iter+1);
			if (iter>0) Qnew(Slice(0,-2),Slice(0,-2)) = Q;
			for (int i=0; i<=iter; i++) {
				Qnew(i,iter) = inner(ulist[i],rlist[iter]);
				Qnew(iter,i) = inner(ulist[iter],rlist[i]);
			}
			Q = Qnew;
			real_tensor c = KAIN(Q);
			check_linear_dependence(Q,c,rcondtol,cabsmax);
			if (do_print) print("subspace solution",c);

			// Form new solution in u
			Function<double,NDIM> unew = FunctionFactory<double,NDIM>(u.world());
			if (ulist[0].is_compressed()) unew.compress();
			for (int i=0; i<=iter; i++) {
				unew.gaxpy(1.0,ulist[i], c[i]);
				unew.gaxpy(1.0,rlist[i],-c[i]);
			}
			unew.truncate();

			if (ulist.size() == maxsub) {
				ulist.erase(ulist.begin());
				rlist.erase(rlist.begin());
				Q = copy(Q(Slice(1,-1),Slice(1,-1)));
			}
			return unew;
		}
	};

	// backwards compatibility
	typedef NonlinearSolverND<3> NonlinearSolver;


	template <class T>
	struct default_allocator {
        T operator()() {return T();}
    };


    // The default constructor for functions does not initialize
    // them to any value, but the solver needs functions initialized
    // to zero for which we also need the world object.
    template<typename T, std::size_t NDIM>
    struct vector_function_allocator {
        World& world;
        const int n=-1;

        /// @param[in]	world	the world
        /// @param[in]	nn		the number of functions in a given vector
        vector_function_allocator(World& world, const int nn) : world(world), n(nn) {}

        /// allocate a vector of n empty functions
        std::vector<Function<T, NDIM> > operator()() {
            return zero_functions<T, NDIM>(world, n);
        }
    };



    /// Generalized version of NonlinearSolver not limited to a single madness function

    /// \ingroup nonlinearsolve 
    ///
    /// This solves the equation \f$r(u) = 0\f$ where u and r are both
    /// of type \c T and inner products between two items of type \c T
    /// produce a number of type \c C (defaulting to double).  The type \c T
    /// must support storage in an STL vector, scaling by a constant
    /// of type \c C, inplace addition (+=), subtraction, allocation with
    /// value zero, and inner products computed with the interface \c
    /// inner(a,b).  Have a look in examples/testsolver.cc for a
    /// simple but complete example, and in examples/h2dynamic.cc for a 
    /// more complex example.
    ///
    /// I've not yet tested with anything except \c C=double and I think
    /// that the KAIN routine will need extending for anything else.
    template <class T, class C = double, class Alloc = default_allocator<T> >
    class XNonlinearSolver {
        unsigned int maxsub; ///< Maximum size of subspace dimension
        Alloc alloc;
        std::vector<T> ulist, rlist; ///< Subspace information
        Tensor<C> Q;
        Tensor<C> c;		///< coefficients for linear combination
    public:
        bool do_print;

	XNonlinearSolver(const Alloc& alloc = Alloc(),bool print=false)
            : maxsub(10)
            , alloc(alloc)
    		, do_print(print)
        {}

	XNonlinearSolver(const XNonlinearSolver& other)
            : maxsub(other.maxsub)
            , alloc(other.alloc)
			, do_print(other.do_print)
        {}


	std::vector<T>& get_ulist() {return ulist;}
	std::vector<T>& get_rlist() {return rlist;}

	void set_maxsub(int maxsub) {this->maxsub = maxsub;}
	Tensor<C> get_c() const {return c;}

	void clear_subspace() {
		ulist.clear();
		rlist.clear();
		Q=Tensor<C>();
	}

	/// Computes next trial solution vector

	/// You are responsible for performing step restriction or line search
	/// (not necessary for linear problems).
	///
	/// @param u Current solution vector
	/// @param r Corresponding residual
	/// @return Next trial solution vector
        /// @param[in]          rcondtol rcond less than this will cause the subspace to be shrunk due to linear dependence
        /// @param[in]          cabsmax  maximum element of c greater than this will cause the subspace to be shrunk due to li
	T update(const T& u, const T& r, const double rcondtol=1e-8, const double cabsmax=1000.0) {
		if (maxsub==1) return u-r;
		int iter = ulist.size();
		ulist.push_back(u);
		rlist.push_back(r);

		// Solve subspace equations
		Tensor<C> Qnew(iter+1,iter+1);
		if (iter>0) Qnew(Slice(0,-2),Slice(0,-2)) = Q;
		for (int i=0; i<=iter; i++) {
			Qnew(i,iter) = inner(ulist[i],rlist[iter]);
			Qnew(iter,i) = inner(ulist[iter],rlist[i]);
		}
		Q = Qnew;
		c = KAIN(Q);

		check_linear_dependence(Q,c,rcondtol,cabsmax,do_print);
		if (do_print) print("subspace solution",c);

		// Form new solution in u
		T unew = alloc();
		for (int i=0; i<=iter; i++) {
			unew += (ulist[i] - rlist[i])*c[i];
		}

		if (ulist.size() == maxsub) {
			ulist.erase(ulist.begin());
			rlist.erase(rlist.begin());
			Q = copy(Q(Slice(1,-1),Slice(1,-1)));
		}
		return unew;
	}

    };


    template<typename T, std::size_t NDIM>
    static inline XNonlinearSolver<std::vector<Function<T,NDIM>>,T,vector_function_allocator<T,NDIM>>
    nonlinear_vector_solver(World& world, const long nvec) {
        auto alloc=vector_function_allocator<T,NDIM>(world,nvec);
        return XNonlinearSolver<std::vector<Function<T,NDIM>>,T,vector_function_allocator<T,NDIM>>(alloc);
    };


    typedef XNonlinearSolver<std::vector<Function<double,3>>,double,vector_function_allocator<double,3>> NonlinearVectorSolver_3d;
    typedef XNonlinearSolver<std::vector<Function<double,6>>,double,vector_function_allocator<double,6>> NonlinearVectorSolver_6d;


}
#endif
