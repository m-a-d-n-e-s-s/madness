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
#include <tensor/mtrand.h>
#include <tensor/tensor.h>

#define FUNCTION_INSTANTIATE_1
#define FUNCTION_INSTANTIATE_2
#define FUNCTION_INSTANTIATE_3

namespace madness {
    void startup(World& world, int argc, char** argv);

    /// Translation in 1d ... more than 31 levels of refinement will require wide integers
    typedef unsigned long Translation;

    /// Level
    typedef long Level;
}

#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <mra/key.h>
#include <mra/funcimpl.h>
#include <mra/loadbal.h>

namespace madness {

    template <typename T, int NDIM>
    class Function {
    private:
        SharedPtr< FunctionImpl<T,NDIM> > impl;

        inline void verify() const {
            MADNESS_ASSERT(impl);
        };

    public:
        typedef FunctionImpl<T,NDIM> implT;
        typedef FunctionFactory<T,NDIM> factoryT;
        typedef typename implT::coordT coordT; ///< Type of vector holding coordinates 

        /// Default constructor makes uninitialized function.  No communication.

        /// An unitialized function can only be assigned to.  Any other operation will throw.
        Function()
            : impl(0)
        {};


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const factoryT& factory)
            : impl(new FunctionImpl<T,NDIM>(factory))
        {};


        /// Copy constructor is \em shallow.  No communication, works in either basis.
        Function(const Function<T,NDIM>& f)
            : impl(f.impl)
        {};


        /// Assignment is \em shallow.  No communication, works in either basis.
        Function<T,NDIM>& operator=(const Function<T,NDIM>& f) {
            if (this != &f) impl = f.impl;
            return *this;
        };

        /// Destruction of any underlying implementation is deferred to next global fence.
        ~Function(){};

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        ///
        /// Needs a lot of optimization for efficient parallel execution.
        Future<T> eval(const coordT& xuser) {
            verify();
            coordT xsim;
            impl->user_to_sim(xuser,xsim);
            Future<T> result;
            impl->eval(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        };

        /// Evaluates the function at a point in user coordinates.  Possible \em blocking comm.

        /// Only the invoking process will receive the result.
        ///
        /// Throws if function is not initialized.
        ///
        /// This function calls eval and blocks until the result is available.  Therefore,
        /// if you are evaluating many points in parallel it is \em vastly less efficient than
        /// calling eval directly.
        T operator()(const coordT& xuser) {
            return eval(xuser).get();
        };


        /// Returns an estimate of the difference ||this-func||^2 from local data

        /// No communication is performed.  If the function is not
        /// reconstructed, it throws an exception.  To get the global
        /// value either do a global sum of the local values or call
        /// errsq
        template <typename funcT>
        double errsq_local(const funcT& func) const {
            verify();
            if (is_compressed()) MADNESS_EXCEPTION("Function:errsq_local:not reconstructed",0);
            return impl->errsq_local(func);
        };


        /// Returns an estimate of the difference ||this-func|| ... global sum performed

        /// If the function is compressed, it is reconstructed first.  For efficient use
        /// especially with many functions, reconstruct them all first, and use errsq_local
        /// instead so you can perform a global sum on all at the same time.
        template <typename funcT>
        double err(const funcT& func) const {
            verify();
            if (is_compressed()) const_cast<Function<T,NDIM>*>(this)->reconstruct();
            double local = impl->errsq_local(func);
            impl->world.gop.sum(local);
            return sqrt(local);
        };


        /// Returns true if compressed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_compressed() const {
            if (impl) 
                return impl->is_compressed();
            else
                return false;
        };

        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t size() const {
            if (!impl) return 0;
            return impl->size();
        };

        /// Returns the number of multiwavelets (k).  No communication.
        bool k() const {
            verify();
            return impl->get_k();
        };

        /// Truncate the function with optional fence.  Compresses with fence if not compressed.

        /// If the truncation threshold is less than or equal to zero the default value
        /// specified when the function was created is used.
        /// If the function is not initialized, it just returns.
        void truncate(double tol = 0.0, bool fence = true) {
            if (!impl) return;
            verify();
            impl->truncate(tol,fence);
        };

	/// Returns a shared-pointer to the implementation
	const SharedPtr< FunctionImpl<T,NDIM> >& get_impl() const {
	    verify();
	    return impl;
	};


        /// Returns a shared pointer to the process map
        const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
            verify();
            return impl->get_pmap();
        };

        
        /// Returns the square of the norm of the local function ... no communication
        
        /// Works in either basis
        double norm2sq_local() const {
            verify();
            return impl->norm2sq_local();
        };


        /// Returns the 2-norm of the function ... global sum ... works in either basis
        
        /// See comments for err() w.r.t. applying to many functions.
        double norm2() const {
            verify();
            double local = impl->norm2sq_local();
            impl->world.gop.sum(local);
            return sqrt(local);
        };


        /// Compresses the function, transforming into wavelet basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        void compress(bool fence = true) {
            if (!impl || is_compressed()) return;
            impl->compress(fence);
        };
        
        /// Reconstructs the function, transforming into scaling function basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already reconstructed or if not initialized.
        void reconstruct(bool fence = true) {
            if (!impl || !is_compressed()) return;
            impl->reconstruct(fence);
        };


        /// Process 0 prints a summary of all nodes in the tree (collective)
        void print_tree() const {if (impl) impl->print_tree();};


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
            verify();
            Function<Q,NDIM> result;
	    result.impl = SharedPtr< FunctionImpl<Q,NDIM> >(new typename FunctionImpl<Q,NDIM>::implT(*impl));
	    result.impl->copy_coeffs(*impl, fence);
	    return result;
	};


        /// Deep copy generating a new function (with same distribution).  No communication except due to optional fence.

        /// Works in either basis.  
	Function<T,NDIM> copy(bool fence = true) const {
            verify();
            return this->copy(get_pmap(), fence);
	};


        /// Deep copy generating a new function with change of process map and optional fence

        /// Works in either basis.  Different distributions imply
        /// asynchronous communication and the optional fence is
        /// collective.
	Function<T,NDIM> copy(const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap, bool fence = true) const {
            verify();
	    Function<T,NDIM> result;
	    result.impl = SharedPtr<implT>(new implT(*impl, pmap));
	    result.impl->copy_coeffs(*impl, fence);
	    return result;
	};


        /// Inplace, scale the function by a constant.  No communication except for optional fence.

        /// Works in either basis.  Returns reference to this for chaining.
        template <typename Q>
        Function<T,NDIM>& scale_inplace(const Q q, bool fence=true) {
            verify();
            impl->scale_inplace(q,fence);
            return *this;
        };

        /// Inplace, general bi-linear operation in wavelet basis.  No communication except for optional fence.

        /// If the functions are not in the wavelet basis they are compressed with implied communication
        /// and a forced global fence.  Returns this for chaining.
        ///
        /// this <-- this*alpha + other*beta
        template <typename Q, typename R>
        Function<T,NDIM>& gaxpy_inplace(const T& alpha, 
                                        const Function<Q,NDIM>& other, const R& beta, bool fence=true) {
            verify();
            other.verify();
            if (!is_compressed()) compress();
            if (!other.is_compressed()) const_cast<Function<Q,NDIM>*>(&other)->compress();
            impl->gaxpy_inplace(alpha,*other.impl,beta,fence);
            return *this;
        };
    private:

    };

    /// Create a new copy of the function with different distribution and optional fence
    template <typename T, int NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f, 
			  const SharedPtr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
                          bool fence = true) {
	return f.copy(pmap,fence);
    };
	    
    /// Create a new copy of the function with the same distribution and optional fence
    template <typename T, int NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f, bool fence = true) {
	return f.copy(fence);
    };
	    




}



#endif
