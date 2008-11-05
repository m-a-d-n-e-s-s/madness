/*
  This file is part of MADNESS.

  Copyright (C) <2007> <Oak Ridge National Laboratory>

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


#ifndef MAD_MRA_H
#define MAD_MRA_H

#include <world/world.h>
#include <misc/misc.h>
#include <tensor/tensor.h>

#define FUNCTION_INSTANTIATE_1
#define FUNCTION_INSTANTIATE_2
#define FUNCTION_INSTANTIATE_3
#define FUNCTION_INSTANTIATE_4

static const bool VERIFY_TREE = false; //true;


namespace madness {
    void startup(World& world, int argc, char** argv);
}

#include <mra/key.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <mra/indexit.h>
#include <mra/funcimpl.h>
#include <mra/operator.h>
#include <mra/loadbal.h>

namespace madness {

    // Forward declarations of friends demanded by PGI

    template <typename T, int NDIM>
    Function<T,NDIM> square(const Function<T,NDIM>& f, bool fence = true);

    template <typename T, int NDIM>
    Function<T,NDIM>
    diff(const Function<T,NDIM>& f, int axis, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R),D>
    mul(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

    template <typename Q, typename R, int D>
    Function<TENSOR_RESULT_TYPE(Q,R),D>
    mul(const Q alpha, const Function<R,D>& f, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R),D>
    mul_sparse(const Function<L,D>& left, const Function<R,D>& right, double tol, bool fence=true);

    template <typename Q, typename R, int D>
    Function<TENSOR_RESULT_TYPE(Q,R),D>
    mul_sparse(const Q alpha, const Function<R,D>& f, double tol, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R),D>
    gaxpy_oop(TENSOR_RESULT_TYPE(L,R) alpha, const Function<L,D>& left,
              TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,D>& right, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R),D>
    add(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R),D>
    sub(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R), D>
    operator+(const Function<L,D>& left, const Function<R,D>& right);

    template <typename L, typename R, int D>
    Function<TENSOR_RESULT_TYPE(L,R), D>
    operator-(const Function<L,D>& left, const Function<R,D>& right);

    template <typename opT, typename R, int D>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), D>
    apply_only(const opT& op, const Function<R,D>& f, bool fence=true);

    template <typename opT, typename R, int D>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), D>
    apply(const opT& op, const Function<R,D>& f, bool fence=true);


    template <typename T, int NDIM>
    Function<T,NDIM>
    project(const Function<T,NDIM>& other,
            int k=FunctionDefaults<NDIM>::get_k(),
            double thresh=FunctionDefaults<NDIM>::get_thresh(),
            bool fence=true);

    template <typename T, int NDIM>
    class Function : public ParallelSerializableObject {

        template<typename Q, int D>
        friend
        class Function;

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R),D>
        madness::mul(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

        template <typename Q, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(Q,R),D>
        madness::mul(const Q alpha, const Function<R,D>& f, bool fence=true);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R),D>
        madness::mul_sparse(const Function<L,D>& left, const Function<R,D>& right, double tol, bool fence=true);

        template <typename Q, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(Q,R),D>
        madness::mul_sparse(const Q alpha, const Function<R,D>& f, double tol, bool fence=true);

        template <typename Q, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(Q,R),D>
        madness::operator*(const Function<R,D>& f, const Q alpha);

        template <typename Q, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(Q,R),D>
        madness::operator*(const Q alpha, const Function<R,D>& f);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R),D>
        madness::gaxpy_oop(TENSOR_RESULT_TYPE(L,R) alpha, const Function<L,D>& left,
                           TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,D>& right, bool fence=true);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R),D>
        madness::add(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R),D>
        madness::sub(const Function<L,D>& left, const Function<R,D>& right, bool fence=true);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R), D>
        madness::operator+(const Function<L,D>& left, const Function<R,D>& right);

        template <typename L, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(L,R), D>
        madness::operator-(const Function<L,D>& left, const Function<R,D>& right);

        template <typename opT, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), D>
        madness::apply_only(const opT& op, const Function<R,D>& f, bool fence=true);

        template <typename opT, typename R, int D>
        friend
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), D>
        madness::apply(const opT& op, const Function<R,D>& f, bool fence=true);

        friend Function<T,NDIM> square<T,NDIM>(const Function<T,NDIM>&, bool);

        friend Function<T,NDIM> project<T,NDIM>(const Function<T,NDIM>&, int, double, bool);

        friend Function<T,NDIM> diff<T,NDIM>(const Function<T,NDIM>& f, int axis, bool fence);

    private:
        SharedPtr< FunctionImpl<T,NDIM> > impl;


    public:
        /// Asserts that the function is initialized
        inline void verify() const {
            MADNESS_ASSERT(impl);
        }

        typedef FunctionImpl<T,NDIM> implT;
        typedef FunctionFactory<T,NDIM> factoryT;
        typedef typename implT::coordT coordT; ///< Type of vector holding coordinates

        /// Default constructor makes uninitialized function.  No communication.

        /// An unitialized function can only be assigned to.  Any other operation will throw.
        Function()
            : impl(0)
        {}


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const factoryT& factory)
            : impl(new FunctionImpl<T,NDIM>(factory))
        {
            PROFILE_MEMBER_FUNC(Function);
        }


        /// Copy constructor is \em shallow.  No communication, works in either basis.
        Function(const Function<T,NDIM>& f)
            : impl(f.impl)
        {}


        /// Assignment is \em shallow.  No communication, works in either basis.
        Function<T,NDIM>& operator=(const Function<T,NDIM>& f) {
            if (this != &f) impl = f.impl;
            return *this;
        }

        /// Destruction of any underlying implementation is deferred to next global fence.
        ~Function(){}

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        Future<T> eval(const coordT& xuser) const {
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(!is_compressed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (int d=0; d<NDIM; d++) {
                if (xsim[d] < -eps) {
                    MADNESS_EXCEPTION("eval: coordinate lower-bound error in dimension", d);
                }
                else if (xsim[d] < eps) {
                    xsim[d] = eps;
                }

                if (xsim[d] > 1.0+eps) {
                    MADNESS_EXCEPTION("eval: coordinate upper-bound error in dimension", d);
                }
                else if (xsim[d] > 1.0-eps) {
                    xsim[d] = 1.0-eps;
                }
            }

            Future<T> result;
            impl->eval(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        }

        /// Evaluates a cube/slice of points (probably for plotting) ... collective but no fence necessary

        /// All processes recieve the entire result (which is a rather severe limit
        /// on the size of the cube that is possible).
        Tensor<T> eval_cube(const Tensor<double>& cell, const vector<long>& npt) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-14;
            verify();
            reconstruct();
            coordT simlo, simhi;
            for (int d=0; d<NDIM; d++) {
                simlo[d] = cell(d,0);
                simhi[d] = cell(d,1);
            }
            user_to_sim(simlo, simlo);
            user_to_sim(simhi, simhi);

            // Move the bounding box infintesimally inside dyadic
            // points so that the evaluation logic does not fail
            for (int d=0; d<NDIM; d++) {
                MADNESS_ASSERT(simhi[d] >= simlo[d]);
                MADNESS_ASSERT(simlo[d] >= 0.0);
                MADNESS_ASSERT(simhi[d] <= 1.0);

                double delta = eps*simhi[d];
                simlo[d] += delta;
                simhi[d] -= 2*delta;  // deliberate asymmetry
            }
            //madness::print("plotbox in sim", simlo, simhi);
            Tensor<T> r = impl->eval_plot_cube(simlo, simhi, npt);
            impl->world.gop.sum(r.ptr(), r.size);
            impl->world.gop.fence();
            return r;
        }


        /// Evaluates the function at a point in user coordinates.  Collective operation.

        /// *CHANGE* All processess receive the result.
        ///
        /// Throws if function is not initialized.
        ///
        /// This function calls eval, blocks until the result is
        /// available and then broadcasts the result to everyone.
        /// Therefore, if you are evaluating many points in parallel
        /// it is \em vastly less efficient than calling eval
        /// directly, saving the futures, and then forcing all of the
        /// results.
        T operator()(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            T result;
            if (impl->world.rank() == 0) result = eval(xuser).get();
            impl->world.gop.broadcast(result);
            impl->world.gop.fence();
            return result;
        }

        /// Returns an estimate of the difference ||this-func||^2 from local data

        /// No communication is performed.  If the function is not
        /// reconstructed, it throws an exception.  To get the global
        /// value either do a global sum of the local values or call
        /// errsq
        template <typename funcT>
        double errsq_local(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) MADNESS_EXCEPTION("Function:errsq_local:not reconstructed",0);
            return impl->errsq_local(func);
        }


        /// Returns an estimate of the difference ||this-func|| ... global sum performed

        /// If the function is compressed, it is reconstructed first.  For efficient use
        /// especially with many functions, reconstruct them all first, and use errsq_local
        /// instead so you can perform a global sum on all at the same time.
        template <typename funcT>
        double err(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            double local = impl->errsq_local(func);
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return sqrt(local);
        }

        /// Verifies the tree data structure ... global sync implied
        void verify_tree() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->verify_tree();
        }


        /// Returns true if compressed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_compressed() const {
            if (impl)
                return impl->is_compressed();
            else
                return false;
        }


        /// Returns the number of nodes in the function tree ... collective global sum
	std::size_t tree_size() const {
            PROFILE_MEMBER_FUNC(Function);
	    if (!impl) return 0;
	    return impl->tree_size();
	}


	/// Returns the maximum depth of the function tree
	std::size_t max_depth() const {
            PROFILE_MEMBER_FUNC(Function);
	    if (!impl) return 0;
	    return impl->max_depth();
	}


        /// Returns the max number of nodes on a processor
	std::size_t max_nodes() const {
            PROFILE_MEMBER_FUNC(Function);
	    if (!impl) return 0;
	    return impl->max_nodes();
	}

        /// Returns the min number of nodes on a processor
	std::size_t min_nodes() const {
            PROFILE_MEMBER_FUNC(Function);
	    if (!impl) return 0;
	    return impl->min_nodes();
	}


        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t size() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->size();
        }


        /// Returns value of autorefine flag.  No communication.
        bool autorefine() const {
            if (!impl) return true;
            return impl->autorefine;
        }


        /// Sets the value of the autorefine flag.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_autorefine(bool value, bool fence = true) {
            verify();
            impl->autorefine = value;
            if (fence) impl->world.gop.fence();
        }


        /// Returns value of truncation threshold.  No communication.
        double thresh() const {
            if (!impl) return 0.0;
            return impl->thresh;
        }


        /// Sets the vaule of the truncation threshold.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_thresh(double value, bool fence = true) {
            verify();
            impl->thresh = value;
            if (fence) impl->world.gop.fence();
        }


        /// Returns the number of multiwavelets (k).  No communication.
        int k() const {
            verify();
            return impl->k;
        }


        /// Truncate the function with optional fence.  Compresses with fence if not compressed.

        /// If the truncation threshold is less than or equal to zero the default value
        /// specified when the function was created is used.
        /// If the function is not initialized, it just returns.
        ///
        /// Returns this for chaining.
        Function<T,NDIM>& truncate(double tol = 0.0, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return *this;
            verify();
            if (!is_compressed()) compress();
            impl->truncate(tol,fence);
            if (VERIFY_TREE) verify_tree();
            return *this;
        }



	/// Returns a shared-pointer to the implementation
	const SharedPtr< FunctionImpl<T,NDIM> >& get_impl() const {
	    verify();
	    return impl;
	}

	/// Returns the world
	World& world() const {
	  verify();
	  return  impl->world;
	}


        /// Returns a shared pointer to the process map
        const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
            verify();
            return impl->get_pmap();
        }



        /// Returns the square of the norm of the local function ... no communication

        /// Works in either basis
        double norm2sq_local() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->norm2sq_local();
        }



        /// Returns the 2-norm of the function ... global sum ... works in either basis

        /// See comments for err() w.r.t. applying to many functions.
        double norm2() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            double local = impl->norm2sq_local();

            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return sqrt(local);
        }


        /// Initializes information about the function norm at all length scales
        void norm_tree(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            const_cast<Function<T,NDIM>*>(this)->impl->norm_tree(fence);
        }


        /// Compresses the function, transforming into wavelet basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        ///
        /// Since reconstruction/compression do not discard information we define them
        /// as const ... "logical constness" not "bitwise constness".
        const Function<T,NDIM>& compress(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl || is_compressed()) return *this;
            if (VERIFY_TREE) verify_tree();
            const_cast<Function<T,NDIM>*>(this)->impl->compress(false, false, fence);
            return *this;
        }


        /// Compresses the function retaining scaling function coeffs.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        void nonstandard(bool keepleaves, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
	    if (impl->nonstandard) return;
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            impl->compress(true, keepleaves, fence);
        }


        void standard(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            MADNESS_ASSERT(is_compressed());
            impl->standard(fence);
            if (fence && VERIFY_TREE) verify_tree();
        }

        /// Reconstructs the function, transforming into scaling function basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already reconstructed or if not initialized.
        ///
        /// Since reconstruction/compression do not discard information we define them
        /// as const ... "logical constness" not "bitwise constness".
        void reconstruct(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl || !is_compressed()) return;
            const_cast<Function<T,NDIM>*>(this)->impl->reconstruct(fence);
            if (fence && VERIFY_TREE) verify_tree(); // Must be after in case nonstandard
        }


        /// Inplace autorefines the function using the same test as for squaring.  Possible non-blocking comm.

        /// This needs generalizing to a user-defined threshold and criterion.
        void refine(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) reconstruct();
            impl->refine(fence);
        }

	/// Get the scaling function coeffs at level n starting from NS form
	Tensor<T> coeffs_for_jun(Level n, long mode=0) {
	    nonstandard(true,true);
	    return impl->coeffs_for_jun(n,mode);
	    //return impl->coeffs_for_jun(n);
	}

	  ///Change bv on the fly. Temporary workaround until better bc handling is introduced.
		void set_bc(const Tensor<int>& value) { impl->set_bc(value);}

	/// Clears the function as if constructed uninitialized.  Optional fence.

        /// Any underlying data will not be freed until the next global fence.
        void clear(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) {
                World& world = impl->world;
                impl = SharedPtr< FunctionImpl<T,NDIM> >(0);
                if (fence) world.gop.fence();
            }
        }

        /// Process 0 prints a summary of all nodes in the tree (collective)
        void print_tree() const {
            if (impl) impl->print_tree();
        }

        /// Print a summary of the load balancing info

        /// This is serial and VERY expensive
        void print_info() const {
            if (impl) impl->print_info();
        }

        /// Type conversion implies a deep copy.  No communication except for optional fence.

        /// Works in either basis but any loss of precision may result in different errors
        /// in applied in a different basis.
        ///
        /// The new function is formed with the options from the default constructor.
        ///
        /// There is no automatic type conversion since this is generally a rather dangerous
        /// thing and because there would be no way to make the fence optional.
        template <typename Q>
        Function<Q,NDIM> convert(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            Function<Q,NDIM> result;
	    result.impl = SharedPtr< FunctionImpl<Q,NDIM> >(new FunctionImpl<Q,NDIM>(*impl));
	    result.impl->copy_coeffs(*impl, fence);
	    return result;
	}

        /// Inplace unary operation on function values with optional autorefining and fence
        template <typename opT>
        void unaryop(const opT& op, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (autorefine()) {
                impl->unary_op_value_inplace(&implT::autorefine_square_test, op, fence);
            }
            else {
                impl->unary_op_value_inplace(op, fence);
            }
        }


        /// Unary operation applied inplace to the coefficients with optional autorefining and fence
        template <typename opT>
        void unaryop_coeff(const opT& op,
                           bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (autorefine()) {
                impl->unary_op_coeff_inplace(&implT::autorefine_square_test, op, fence);
            }
            else {
                impl->unary_op_coeff_inplace(op, fence);
            }
        }


    private:
        static void doconj(const Key<NDIM>, Tensor<T>& t) {
            t.conj();
        }
    public:

        /// Inplace complex conjugate.  No communication except for optional fence.

        /// Returns this for chaining.  Works in either basis.
        Function<T,NDIM> conj(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            unaryop_coeff(&Function<T,NDIM>::doconj, fence);
            return *this;
        }

        /// Deep copy generating a new function (same distribution).  No communication except optional fence.

        /// Works in either basis.
	Function<T,NDIM> copy(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return this->copy(get_pmap(), fence);
	}


        /// Deep copy generating a new function with change of process map and optional fence

        /// Works in either basis.  Different distributions imply
        /// asynchronous communication and the optional fence is
        /// collective.
	Function<T,NDIM> copy(const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap, bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
	    Function<T,NDIM> result;
	    result.impl = SharedPtr<implT>(new implT(*impl, pmap, false));
	    result.impl->copy_coeffs(*impl, fence);
            if (VERIFY_TREE) result.verify_tree();
	    return result;
	}


        /// Inplace, scale the function by a constant.  No communication except for optional fence.

        /// Works in either basis.  Returns reference to this for chaining.
        template <typename Q>
        Function<T,NDIM>& scale(const Q q, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            impl->scale_inplace(q,fence);
            return *this;
        }

        /// Inplace add scalar.  No communication except for optional fence.
        Function<T,NDIM>& add_scalar(T t, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            impl->add_scalar_inplace(t,fence);
            return *this;
        }

        /// Inplace, general bi-linear operation in wavelet basis.  No communication except for optional fence.

        /// If the functions are not in the wavelet basis an exception is thrown since this routine
        /// is intended to be fast and unexpected compression is assumed to be a performance bug.
        ///
        /// Returns this for chaining.
        ///
        /// this <-- this*alpha + other*beta
        template <typename Q, typename R>
        Function<T,NDIM>& gaxpy(const T& alpha,
                                const Function<Q,NDIM>& other, const R& beta, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            other.verify();
            MADNESS_ASSERT(is_compressed() && other.is_compressed());
            impl->gaxpy_inplace(alpha,*other.impl,beta,fence);
            return *this;
        }


        /// Inplace addition of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation.
        /// Functions are compressed if not already so.
        template <typename Q>
        Function<T,NDIM>& operator+=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_compressed()) compress();
            if (!other.is_compressed()) other.compress();
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(1.0), true);
        }


        /// Inplace subtraction of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        Function<T,NDIM>& operator-=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_compressed()) compress();
            if (!other.is_compressed()) other.compress();
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(-1.0), true);
        }


        /// Inplace scaling by a constant

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        Function<T,NDIM>& operator*=(const Q q) {
            PROFILE_MEMBER_FUNC(Function);
            scale(q,true);
            return *this;
        }


        /// Inplace squaring of function ... global comm only if not reconstructed

        /// Returns *this for chaining.
        Function<T,NDIM>& square(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->square_inplace(fence);
            return *this;
        }

        /// Returns local contribution to \c int(f(x),x) ... no communication

        /// In the wavelet basis this is just the coefficient of the first scaling
        /// function which is a constant.  In the scaling function basis we
        /// must add up contributions from each box.
        T trace_local() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            if (VERIFY_TREE) verify_tree();
            return impl->trace_local();
        }


        /// Returns global value of \c int(f(x),x) ... global comm required
        T trace() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            T sum = impl->trace_local();
            impl->world.gop.sum(sum);
            impl->world.gop.fence();
            return sum;
        }


        /// Returns local part of inner product ... throws if both not compressed
	template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);
            MADNESS_ASSERT(is_compressed());
            MADNESS_ASSERT(g.is_compressed());
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();
            return impl->inner_local(*(g.impl));
	}


	/// Returns the inner product

	/// Not efficient for computing multiple inner products
	template <typename R>
	  TENSOR_RESULT_TYPE(T,R) inner(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);
	  if (!is_compressed()) compress();
	  if (!g.is_compressed()) g.compress();
          if (VERIFY_TREE) verify_tree();
          if (VERIFY_TREE) g.verify_tree();

	  TENSOR_RESULT_TYPE(T,R) local = impl->inner_local(*(g.impl));
	  impl->world.gop.sum(local);
          world.gop.fence();
	  return local;
	}


        /// Replaces this function with one loaded from an archive using the default processor map

        /// Archive can be sequential or parallel.
        ///
        /// The & operator for serializing will only work with parallel archives.
        template <typename Archive>
        void load(World& world, Archive& ar) {
            PROFILE_MEMBER_FUNC(Function);
            // Type checking since we are probably circumventing the archive's own type checking
            long magic, id, ndim, k;
            ar & magic & id & ndim & k;
            MADNESS_ASSERT(magic == 7776768); // Mellow Mushroom Pizza tel.# in Knoxville
            MADNESS_ASSERT(id == TensorTypeData<T>::id);
            MADNESS_ASSERT(ndim == NDIM);

            impl = SharedPtr<implT>(new implT(FunctionFactory<T,NDIM>(world).k(k).empty()));

            impl->load(ar);
        }


        /// Stores the function to an archive

        /// Archive can be sequential or parallel.
        ///
        /// The & operator for serializing will only work with parallel archives.
        template <typename Archive>
        void store(Archive& ar) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            // For type checking, etc.
            ar & long(7776768) & long(TensorTypeData<T>::id) & long(NDIM) & long(k());

            impl->store(ar);
        }


      void set_apply_time_ptr(SharedPtr<ApplyTime<NDIM> > ptr) {
	impl->set_apply_time_ptr(ptr);
      }

    private:

        /// Projects inplace function to new order basis ... private
        Function<T,NDIM>& project(const Function<T,NDIM>& other, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            impl->project(*other.impl, fence);
            return *this;
        }

        /// This is replaced with left*right ...  private
        template <typename L, typename R>
        Function<T,NDIM>& mul(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            left.verify();
            right.verify();
            MADNESS_ASSERT(!(left.is_compressed() || right.is_compressed()));
            if (VERIFY_TREE) left.verify_tree();
            if (VERIFY_TREE) right.verify_tree();
            impl = SharedPtr<implT>(new implT(*left.impl, left.get_pmap(), false));
            impl->mul(*left.impl,*right.impl,fence);
            return *this;
        }

        /// This is replaced with left*right using sparsity ...  private
        template <typename L, typename R>
        Function<T,NDIM>& mul_sparse(const Function<L,NDIM>& left, const Function<R,NDIM>& right, double tol, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            left.verify();
            right.verify();
            MADNESS_ASSERT(!(left.is_compressed() || right.is_compressed()));
            if (VERIFY_TREE) left.verify_tree();
            if (VERIFY_TREE) right.verify_tree();
            impl = SharedPtr<implT>(new implT(*left.impl, left.get_pmap(), false));
            impl->mul_sparse(*left.impl,*right.impl,tol,fence);
            return *this;
        }

        /// This is replaced with alpha*left + beta*right ...  private
        template <typename L, typename R>
        Function<T,NDIM>& gaxpy_oop(T alpha, const Function<L,NDIM>& left,
                                    T beta,  const Function<R,NDIM>& right, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            left.verify();
            right.verify();
            MADNESS_ASSERT(left.is_compressed() && right.is_compressed());
            if (VERIFY_TREE) left.verify_tree();
            if (VERIFY_TREE) right.verify_tree();
            impl = SharedPtr<implT>(new implT(*left.impl, left.get_pmap(), false));
            impl->gaxpy(alpha,*left.impl,beta,*right.impl,fence);
            return *this;
        }

        /// This is replaced with alpha*f ...  private
        template <typename Q, typename L>
        Function<T,NDIM>& scale_oop(const Q alpha, const Function<L,NDIM>& f, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            impl = SharedPtr<implT>(new implT(*f.impl, f.get_pmap(), false));
            impl->scale_oop(alpha,*f.impl,fence);
            return *this;
        }

        /// This is replaced with df/dx ...  private.
        Function<T,NDIM>& diff(const Function<T,NDIM>& f, int axis, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            impl = SharedPtr<implT>(new implT(*f.impl, f.get_pmap(), false));
            impl->diff(*f.impl,axis,fence);
            return *this;
        }

        /// This is replaced with op(f) in NS form ...  private.
        template <typename opT, typename R>
        Function<T,NDIM>& apply(opT& op, const Function<R,NDIM>& f, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            impl = SharedPtr<implT>(new implT(*f.impl, f.get_pmap(), true));
            impl->apply(op, *f.impl, fence);
            return *this;
        }

        /// This is replaced with mapdim(f) ...  private
        Function<T,NDIM>& mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            for (int i=0; i<NDIM; i++) MADNESS_ASSERT(map[i]>=0 && map[i]<NDIM);
            impl = SharedPtr<implT>(new implT(*f.impl, f.get_pmap(), false));
            impl->mapdim(*f.impl,map,fence);
            return *this;
        }

    };


    /// Returns new function equal to alpha*f(x) with optional fence
    template <typename Q, typename T, int NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Q alpha, const Function<T,NDIM>& f, bool fence) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(Q,T),NDIM> result;
        return result.scale_oop(alpha, f, fence);
    }


    /// Returns new function equal to f(x)*alpha with optional fence
    template <typename Q, typename T, int NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Function<T,NDIM>& f, const Q alpha, bool fence) {
        PROFILE_FUNC;
        return mul(alpha,f,fence);
    }


    /// Returns new function equal to f(x)*alpha

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, int NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Function<T,NDIM>& f, const Q alpha) {
        return mul(alpha, f, true);
    }

    /// Returns new function equal to alpha*f(x)

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, int NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Q alpha, const Function<T,NDIM>& f) {
        return mul(alpha, f, true);
    }


    /// Same as \c operator* but with optional fence and no automatic reconstruction
    template <typename L, typename R,int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence) {
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        return result.mul(left,right,fence);
    }


    /// Sparse multiplication --- right *must* have tree of norms already created
    template <typename L, typename R,int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul_sparse(const Function<L,NDIM>& left, const Function<R,NDIM>& right, double tol, bool fence) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        return result.mul_sparse(left,right,tol,fence);
    }



    /// Multiplies two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation but also
    /// enables us to automatically reconstruct the input functions as required.
    template <typename L, typename R, int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator*(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (left.is_compressed())  left.reconstruct();
        if (right.is_compressed()) right.reconstruct();
        return mul(left,right,true);
    }


    /// Returns new function alpha*left + beta*right optional fence and no automatic compression
    template <typename L, typename R,int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    gaxpy_oop(TENSOR_RESULT_TYPE(L,R) alpha, const Function<L,NDIM>& left,
          TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,NDIM>& right, bool fence) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        return result.gaxpy_oop(alpha, left, beta, right, fence);
    }

    /// Same as \c operator+ but with optional fence and no automatic compression
    template <typename L, typename R,int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    add(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)(1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(1.0), right, fence);
    }


    /// Adds two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator+(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (VERIFY_TREE) left.verify_tree();
        if (VERIFY_TREE) right.verify_tree();
        if (!left.is_compressed())  left.compress();
        if (!right.is_compressed()) right.compress();
        return add(left,right,true);
    }

    /// Same as \c operator- but with optional fence and no automatic compression
    template <typename L, typename R,int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    sub(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)( 1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(-1.0), right, fence);
    }


    /// Subtracts two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, int NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator-(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        PROFILE_FUNC;
        if (!left.is_compressed())  left.compress();
        if (!right.is_compressed()) right.compress();
        return sub(left,right,true);
    }


    /// Create a new function that is the square of f - global comm only if not reconstructed
    template <typename T, int NDIM>
    Function<T,NDIM> square(const Function<T,NDIM>& f, bool fence) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.square(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }


    /// Create a new copy of the function with different distribution and optional fence
    template <typename T, int NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f,
			  const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
                          bool fence = true) {
        PROFILE_FUNC;
	return f.copy(pmap,fence);
    }

    /// Create a new copy of the function with the same distribution and optional fence
    template <typename T, int NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
	return f.copy(fence);
    }

    /// Return the complex conjugate of the input function with the same distribution and optional fence
    template <typename T, int NDIM>
    Function<T,NDIM> conj(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);
	return result.conj(fence);
    }

    /// Differentiate w.r.t. given coordinate (x=0, y=1, ...) with optional fence

    /// Returns a new function with the same distribution
    template <typename T, int NDIM>
    Function<T,NDIM>
    diff(const Function<T,NDIM>& f, int axis, bool fence) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        if (f.is_compressed()) {
            if (fence) {
                f.reconstruct();
            }
            else {
                MADNESS_EXCEPTION("diff: trying to diff a compressed function without fencing",0);
            }
        }

        return result.diff(f,axis, fence);
    }


    /// Apply operator in non-standard form

    /// Returns a new function with the same distribution
    ///
    /// !!! For the moment does NOT respect fence option ... always fences
    template <typename opT, typename R, int NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply(const opT& op, const Function<R,NDIM>& f, bool fence) {
        PROFILE_FUNC;
	Function<R,NDIM>& ff = const_cast< Function<R,NDIM>& >(f);
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;
        if (VERIFY_TREE) ff.verify_tree();
	ff.reconstruct();
	ff.nonstandard(op.doleaves, true);
	result.apply(op, f, true);
	result.reconstruct();
        if (VERIFY_TREE) result.verify_tree();
	ff.standard();
        return result;
    }

    /// Apply operator ONLY in non-standard form - required other steps missing !!
    template <typename opT, typename R, int NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply_only(const opT& op, const Function<R,NDIM>& f, bool fence) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;
	return result.apply(op, f, fence);
    }


    /// Generate a new function by reordering dimensions ... optional fence

    /// You provide an array of dimension NDIM that maps old to new dimensions
    /// according to
    /// \code
    ///    newdim = mapdim[olddim]
    /// \endcode
    /// Otherwise the process map of the input function is used.
    ///
    /// Works in either scaling function or wavelet basis.
    ///
    /// Would be easy to modify this to also change the procmap here
    /// if desired but presently it uses the same procmap as f.
    template <typename T, int NDIM>
    Function<T,NDIM>
    mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        return result.mapdim(f,map,fence);
    }

    template <typename T, int NDIM>
    Function<T,NDIM>
    project(const Function<T,NDIM>& other,
            int k,
            double thresh,
            bool fence)
    {
        PROFILE_FUNC;
        Function<T,NDIM> r = FunctionFactory<T,NDIM>(other.impl->world).k(k).thresh(thresh).empty();
        other.reconstruct();
        r.project(other,fence);
        return r;
    }


    /// Computes the scalar/inner product between two functions

    /// In Maple this would be \c int(conjugate(f(x))*g(x),x=-infinity..infinity)
    template <typename T, typename R, int NDIM>
      TENSOR_RESULT_TYPE(T,R) inner(const Function<T,NDIM>& f, const Function<R,NDIM>& g) {
        PROFILE_FUNC;
      return f.inner(g);
    }


    /// Writes an OpenDX format file with a cube/slice of points on a uniform grid

    /// Collective operation but only process 0 writes the file.  By convention OpenDX
    /// files end in ".dx" but this choice is up to the user.  The binary format is
    /// more compact and vastly faster to both write and load but is not as portable.
    ///
    /// Now follow some brief tips about how to look at files inside OpenDX.
    ///
    /// To view a 1D function \c file-selector-->import-->plot-->image.
    ///
    /// To view a 2D function as a colored plane \c file-selector-->import-->autocolor-->image.
    ///
    /// To view a 2D function as a 3D surface \c file-selector-->import-->rubbersheet-->image.
    ///
    /// To view a 3D function as an isosurface \c file-selector-->import-->isosurface-->image.
    ///
    /// To select the real/imaginary/absolute value of a complex number insert a compute
    /// element after the import.
    template <typename T, int NDIM>
    void plotdx(const Function<T,NDIM>& f,
                const char* filename,
                const Tensor<double>& cell,
                const std::vector<long>& npt,
                bool binary=true);


    static inline void plot_line_print_value(FILE* f, double_complex v) {
        fprintf(f, "    %.6e %.6e   ", real(v), imag(v));
    }

    static inline void plot_line_print_value(FILE* f, double v) {
        fprintf(f, " %.6e", v);
    }


    /// Generates ASCI file tabulating f(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, int NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (int i=0; i<NDIM; i++) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; i++) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.6e ", i*sum);
                plot_line_print_value(file, f(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    /// Generates ASCI file tabulating f(r) and g(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, int NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (int i=0; i<NDIM; i++) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; i++) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.6e ", i*sum);
                plot_line_print_value(file, f(r));
                plot_line_print_value(file, g(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }


    /// Generates ASCI file tabulating f(r), g(r), and a(r) at npoints along line r=lo,...,hi

    /// The ordinate is distance from lo
    template <typename T, typename U, typename V, int NDIM>
    void plot_line(const char* filename, int npt, const Vector<double,NDIM>& lo, const Vector<double,NDIM>& hi,
                   const Function<T,NDIM>& f, const Function<U,NDIM>& g, const Function<V,NDIM>& a) {
        typedef Vector<double,NDIM> coordT;
        coordT h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (int i=0; i<NDIM; i++) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f.world();
        f.reconstruct();
        g.reconstruct();
        a.reconstruct();
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; i++) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.6e ", i*sum);
                plot_line_print_value(file, f(r));
                plot_line_print_value(file, g(r));
                plot_line_print_value(file, a(r));
                fprintf(file,"\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }


    namespace archive {
        template <class T, int NDIM>
        struct ArchiveLoadImpl< ParallelInputArchive, Function<T,NDIM> > {
            static inline void load(const ParallelInputArchive& ar, Function<T,NDIM>& f) {
                f.load(*ar.get_world(), ar);
            }
        };

        template <class T, int NDIM>
        struct ArchiveStoreImpl< ParallelOutputArchive, Function<T,NDIM> > {
            static inline void store(const ParallelOutputArchive& ar, const Function<T,NDIM>& f) {
                f.store(ar);
            }
        };
    }


}

#include <mra/vmra.h>


#endif
