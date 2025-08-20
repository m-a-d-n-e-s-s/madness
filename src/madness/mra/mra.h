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
*/

#ifndef MADNESS_MRA_MRA_H__INCLUDED
#define MADNESS_MRA_MRA_H__INCLUDED

/*!
  \file mra/mra.h
  \brief Main include file for MADNESS and defines \c Function interface

  \addtogroup mra

*/


#include <madness/world/MADworld.h>
#include <madness/misc/misc.h>
#include <madness/tensor/tensor.h>

#define FUNCTION_INSTANTIATE_1
#define FUNCTION_INSTANTIATE_2
#define FUNCTION_INSTANTIATE_3
#if !defined(HAVE_IBMBGP) || !defined(HAVE_IBMBGQ)
#define FUNCTION_INSTANTIATE_4
#define FUNCTION_INSTANTIATE_5
#define FUNCTION_INSTANTIATE_6
#endif

static const bool VERIFY_TREE = false; //true


namespace madness {
    /// @brief initialize the internal state of the MADmra library
    ///
    /// Reads in (and broadcasts across \p world) the twoscale and autocorrelation coefficients,
    /// Gauss-Legendre quadrature roots/weights, function defaults and operator displacement lists.
    /// \warning By default this generates operator displacement lists (see Displacements) for up to 6-d free
    ///          and 3-d periodic boundary conditions. For optimal support for mixed boundary conditions
    ///          (periodic along some axes only) assign the desired boundary conditions
    ///          as default (e.g. `FunctionDefaults<3>::set_bc(BoundaryConditions<3>({BC_FREE, BC_FREE, BC_FREE, BC_FREE, BC_PERIODIC, BC_PERIODIC})`)
    ///          prior to calling this. This will make operator application with such boundary conditions
    ///          as efficient as possible, but will not allow the use of operators with
    ///          other boundary conditions that include periodic axes until Displacements::reset_periodic_axes is invoked.
    ///          By default efficiency is sacrificed for generality.
    /// \param world broadcast data across this World
    /// \param argc command-line parameter count
    /// \param argv command-line parameters array
    /// \param doprint if true, will log status to std::cout on rank 0 [default=false]
    /// \param make_stdcout_nice_to_reals if true, will configure std::cout to print reals prettily, according to the MADNESS convention [default=true]
    void startup(World& world, int argc, char** argv, bool doprint=false, bool make_stdcout_nice_to_reals = true);
    std::string get_mra_data_dir();
}

#include <madness/mra/key.h>
#include <madness/mra/twoscale.h>
#include <madness/mra/legendre.h>
#include <madness/mra/indexit.h>
#include <madness/world/parallel_archive.h>
#include <madness/world/worlddc.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/function_factory.h>
#include <madness/mra/lbdeux.h>
#include <madness/mra/funcimpl.h>

// some forward declarations
namespace madness {

    template<typename T, std::size_t NDIM>
    class FunctionImpl;

    template<typename T, std::size_t NDIM>
    class Function;

    template<typename T, std::size_t NDIM>
    class FunctionNode;

    template<typename T, std::size_t NDIM>
    class FunctionFactory;

    template<typename T, std::size_t NDIM>
    class FunctionFunctorInterface;

    template<typename T, std::size_t NDIM>
    struct leaf_op;

    template<typename T, std::size_t NDIM>
    struct mul_leaf_op;

    template<typename T, std::size_t NDIM>
    struct hartree_leaf_op;

    template<typename T, std::size_t NDIM, std::size_t LDIM, typename opT>
    struct hartree_convolute_leaf_op;

    template<typename T, std::size_t NDIM, typename opT>
    struct op_leaf_op;

    template<typename T, std::size_t NDIM>
    struct error_leaf_op;

}


namespace madness {

    /// \ingroup mra
    /// \addtogroup function

    /// A multiresolution adaptive numerical function
    template <typename T, std::size_t NDIM>
    class Function : public archive::ParallelSerializableObject {
        // We make all of the content of Function and FunctionImpl
        // public with the intent of avoiding the cumbersome forward
        // and friend declarations.  However, this open access should
        // not be abused.

    private:
        std::shared_ptr< FunctionImpl<T,NDIM> > impl;

    public:
        bool impl_initialized()const{
        	if(impl==NULL) return false;
        	else return true;
        }
        typedef FunctionImpl<T,NDIM> implT;
        typedef FunctionNode<T,NDIM> nodeT;
        typedef FunctionFactory<T,NDIM> factoryT;
        typedef Vector<double,NDIM> coordT; ///< Type of vector holding coordinates
        typedef T typeT;
        static constexpr std::size_t dimT=NDIM;


        /// Asserts that the function is initialized
        inline void verify() const {
            MADNESS_ASSERT(impl);
        }

        /// Returns true if the function is initialized
        bool is_initialized() const {
            return impl.get();
        }

        /// Default constructor makes uninitialized function.  No communication.

        /// An uninitialized function can only be assigned to.  Any other operation will throw.
        Function() : impl() {}


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const factoryT& factory)
                : impl(new FunctionImpl<T,NDIM>(factory)) {
            PROFILE_MEMBER_FUNC(Function);
        }


        /// Copy constructor is \em shallow.  No communication, works in either basis.
        Function(const Function<T,NDIM>& f)
                : impl(f.impl) {
        }


        /// Assignment is \em shallow.  No communication, works in either basis.
        Function<T,NDIM>& operator=(const Function<T,NDIM>& f) {
            PROFILE_MEMBER_FUNC(Function);
            if (this != &f) impl = f.impl;
            return *this;
        }

        /// Destruction of any underlying implementation is deferred to next global fence.
        ~Function() {}

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        Future<T> eval(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(is_reconstructed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
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

        /// Evaluate function only if point is local returning (true,value); otherwise return (false,0.0)

        /// maxlevel is the maximum depth to search down to --- the max local depth can be
        /// computed with max_local_depth();
        std::pair<bool,T> eval_local_only(const Vector<double,NDIM>& xuser, Level maxlevel) const {
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(is_reconstructed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
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
            return impl->eval_local_only(xsim,maxlevel);
        }

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        ///
        /// This function is a minimally-modified version of eval()
        Future<Level> evaldepthpt(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(is_reconstructed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
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

            Future<Level> result;
            impl->evaldepthpt(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        }


        /// Evaluates the function rank at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        Future<long> evalR(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(is_reconstructed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
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

            Future<long> result;
            impl->evalR(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        }

        /// Evaluates a cube/slice of points (probably for plotting) ... collective but no fence necessary

        /// All processes recieve the entire result (which is a rather severe limit
        /// on the size of the cube that is possible).

        /// Set eval_refine=true to return the refinment levels of
        /// the given function.

        /// @param[in] cell A Tensor describe the cube where the function to be evaluated in
        /// @param[in] npt How many points to evaluate in each dimension
        /// @param[in] eval_refine Wether to return the refinment levels of the given function
        Tensor<T> eval_cube(const Tensor<double>& cell,
                            const std::vector<long>& npt,
                            bool eval_refine = false) const {
            MADNESS_ASSERT(static_cast<std::size_t>(cell.dim(0))>=NDIM && cell.dim(1)==2 && npt.size()>=NDIM);
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-14;
            verify();
            reconstruct();
            coordT simlo, simhi;
            for (std::size_t d=0; d<NDIM; ++d) {
                simlo[d] = cell(d,0);
                simhi[d] = cell(d,1);
            }
            user_to_sim(simlo, simlo);
            user_to_sim(simhi, simhi);

            // Move the bounding box infintesimally inside dyadic
            // points so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
                MADNESS_ASSERT(simhi[d] >= simlo[d]);
                MADNESS_ASSERT(simlo[d] >= 0.0);
                MADNESS_ASSERT(simhi[d] <= 1.0);

                double delta = eps*(simhi[d]-simlo[d]);
                simlo[d] += delta;
                simhi[d] -= 2*delta;  // deliberate asymmetry
            }
            return impl->eval_plot_cube(simlo, simhi, npt, eval_refine);
        }


        /// Evaluates the function at a point in user coordinates.  Collective operation.

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
            verify();
            if (!is_reconstructed()) reconstruct();
            T result;
            if (impl->world.rank() == 0) result = eval(xuser).get();
            impl->world.gop.broadcast(result);
            //impl->world.gop.fence();
            return result;
        }

        /// Evaluates the function at a point in user coordinates.  Collective operation.

        /// See "operator()(const coordT& xuser)" for more info
        T operator()(double x, double y=0, double z=0, double xx=0, double yy=0, double zz=0) const {
            coordT r;
            r[0] = x;
            if (NDIM>=2) r[1] = y;
            if (NDIM>=3) r[2] = z;
            if (NDIM>=4) r[3] = xx;
            if (NDIM>=5) r[4] = yy;
            if (NDIM>=6) r[5] = zz;
            return (*this)(r);
        }

        /// Throws if function is not initialized.
        ///
        /// This function mimics operator() by going through the
        /// tree looking for the depth of the tree at the point.
        /// It blocks until the result is
        /// available and then broadcasts the result to everyone.
        /// Therefore, if you are evaluating many points in parallel
        /// it is \em vastly less efficient than calling evaldepthpt
        /// directly, saving the futures, and then forcing all of the
        /// results.
        Level depthpt(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (!is_reconstructed()) reconstruct();
            Level result;
            if (impl->world.rank() == 0) result = evaldepthpt(xuser).get();
            impl->world.gop.broadcast(result);
            //impl->world.gop.fence();
            return result;
        }

        /// Returns an estimate of the difference ||this-func||^2 from local data

        /// No communication is performed.  If the function is not
        /// reconstructed, it throws an exception.  To get the global
        /// value either do a global sum of the local values or call
        /// errsq
        /// @param[in] func Templated interface to the a user specified function
        template <typename funcT>
        double errsq_local(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (!is_reconstructed()) MADNESS_EXCEPTION("Function:errsq_local:not reconstructed",0);
            return impl->errsq_local(func);
        }


        /// Returns an estimate of the difference ||this-func|| ... global sum performed

        /// If the function is compressed, it is reconstructed first.  For efficient use
        /// especially with many functions, reconstruct them all first, and use errsq_local
        /// instead so you can perform a global sum on all at the same time.
        /// @param[in] func Templated interface to the a user specified function
        template <typename funcT>
        double err(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            if (!is_reconstructed()) reconstruct();
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
            PROFILE_MEMBER_FUNC(Function);
            if (impl)
                return impl->is_compressed();
            else
                return false;
        }

        /// Returns true if reconstructed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_reconstructed() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl)
                return impl->is_reconstructed();
            else
                return false;
        }

        /// Returns true if nonstandard-compressed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_nonstandard() const {
            PROFILE_MEMBER_FUNC(Function);
            return impl ? impl->is_nonstandard() : false;
        }

        /// Returns true if redundant, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_redundant() const {
            PROFILE_MEMBER_FUNC(Function);
            return impl ? impl->is_redundant() : false;
        }

        /// Returns true if redundant_after_merge, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_redundant_after_merge() const {
            PROFILE_MEMBER_FUNC(Function);
            return impl ? impl->is_redundant_after_merge() : false;
        }

        /// Returns the number of nodes in the function tree ... collective global sum
        std::size_t tree_size() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->tree_size();
        }

        /// print some info about this
        void print_size(const std::string name) const {
            if (!impl) {
                print("function",name,"not assigned yet");
            } else {
                impl->print_size(name);
            }
        }

        /// Returns the maximum depth of the function tree ... collective global sum
        std::size_t max_depth() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->max_depth();
        }


        /// Returns the maximum local depth of the function tree ... no communications
        std::size_t max_local_depth() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->max_local_depth();
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

        /// Return the number of coefficients in the function on this processor
        std::size_t size_local() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->size_local();
        }


        /// Returns value of autorefine flag.  No communication.
        bool autorefine() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return true;
            return impl->get_autorefine();
        }


        /// Sets the value of the autorefine flag.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_autorefine(bool value, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->set_autorefine(value);
            if (fence) impl->world.gop.fence();
        }


        /// Returns value of truncation threshold.  No communication.
        double thresh() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            return impl->get_thresh();
        }


        /// Sets the value of the truncation threshold.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_thresh(double value, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->set_thresh(value);
            if (fence) impl->world.gop.fence();
        }


        /// Returns the number of multiwavelets (k).  No communication.
        int k() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->get_k();
        }


        /// Truncate the function with optional fence.  Compresses with fence if not compressed.

        /// If the truncation threshold is less than or equal to zero the default value
        /// specified when the function was created is used.
        /// If the function is not initialized, it just returns.
        ///
        /// Returns this for chaining.
		/// @param[in] tol Tolerance for truncating the coefficients. Default 0.0 means use the implementation's member value \c thresh instead.
		/// @param[in] fence Do fence
        Function<T,NDIM>& truncate(double tol = 0.0, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return *this;
            verify();
//            if (!is_compressed()) compress();
            impl->truncate(tol,fence);
            if (VERIFY_TREE) verify_tree();
            return *this;
        }


        /// Returns a shared-pointer to the implementation
        const std::shared_ptr< FunctionImpl<T,NDIM> >& get_impl() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl;
        }

        /// Replace current FunctionImpl with provided new one
        void set_impl(const std::shared_ptr< FunctionImpl<T,NDIM> >& impl) {
            PROFILE_MEMBER_FUNC(Function);
            this->impl = impl;
        }


        /// Replace the current functor with the provided new one

        /// presumably the new functor will be a CompositeFunctor, which will
        /// change the behavior of the function: multiply the functor with the function
        void set_functor(const std::shared_ptr<FunctionFunctorInterface<T, NDIM> > functor) {
        	this->impl->set_functor(functor);
        	print("set functor in mra.h");
        }

        bool is_on_demand() const {return this->impl->is_on_demand();}

        /// Replace current FunctionImpl with a new one using the same parameters & map as f

        /// If zero is true the function is initialized to zero, otherwise it is empty
        template <typename R>
        void set_impl(const Function<R,NDIM>& f, bool zero = true) {
            impl = std::shared_ptr<implT>(new implT(*f.get_impl(), f.get_pmap(), zero));
	    if (zero) world().gop.fence();
        }

        /// Returns the world
        World& world() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return  impl->world;
        }


        /// Returns a shared pointer to the process map
        const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->get_pmap();
        }

        /// replicate this function, generating a unique pmap
        void replicate(bool fence=true) const {
            verify();
            impl->replicate(fence);
        }

        /// distribute this function according to newmap
        void distribute(std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > newmap) const {
            verify();
            impl->distribute(newmap);
        }


        /// Returns the square of the norm of the local function ... no communication

        /// Works in either basis
        double norm2sq_local() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            MADNESS_CHECK_THROW(is_compressed() or is_reconstructed(),
                "function must be compressed or reconstructed for norm2sq_local");
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
            if (!is_reconstructed()) reconstruct();
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
            return change_tree_state(compressed,fence);
        }


        /// Compresses the function retaining scaling function coeffs.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        void make_nonstandard(bool keepleaves, bool fence=true) const {
            TreeState newstate=TreeState::nonstandard;
            if (keepleaves) newstate=nonstandard_with_leaves;
            change_tree_state(newstate,fence);
        }

        /// Converts the function standard compressed form.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Must be already compressed.
        void standard(bool fence = true) {
            change_tree_state(compressed,fence);
        }

        /// Converts the function to redundant form, i.e. sum coefficients on all levels

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Must be already compressed.
        void make_redundant(bool fence = true) {
            change_tree_state(redundant, fence);
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
        const Function<T,NDIM>& reconstruct(bool fence = true) const {
            return change_tree_state(reconstructed, fence);
        }

        /// changes tree state to given state

        /// Since reconstruction/compression do not discard information we define them
        /// as const ... "logical constness" not "bitwise constness".
        /// @param[in]  finalstate  The final state of the tree
        /// @param[in]  fence       Fence after the operation (might not be respected!!!)
        const Function<T,NDIM>& change_tree_state(const TreeState finalstate, bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (not impl) return *this;
            TreeState current_state = impl->get_tree_state();
            if (finalstate == current_state) return *this;
            MADNESS_CHECK_THROW(current_state != TreeState::unknown, "unknown tree state");

            impl->change_tree_state(finalstate, fence);
            if (fence && VERIFY_TREE) verify_tree();
            return *this;
        }

        /// Sums scaling coeffs down tree restoring state with coeffs only at leaves.  Optional fence.  Possible non-blocking comm.
        void sum_down(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            MADNESS_CHECK_THROW(impl->get_tree_state()==redundant_after_merge, "sum_down requires a redundant_after_merge state");
            const_cast<Function<T,NDIM>*>(this)->impl->sum_down(fence);
            const_cast<Function<T,NDIM>*>(this)->impl->set_tree_state(reconstructed);

            if (fence && VERIFY_TREE) verify_tree(); // Must be after in case nonstandard
        }


        /// Inplace autorefines the function.  Optional fence. Possible non-blocking comm.
        template <typename opT>
        void refine_general(const opT& op, bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (!is_reconstructed()) reconstruct();
            impl->refine(op, fence);
        }


        struct autorefine_square_op {
            bool operator()(implT* impl, const Key<NDIM>& key, const nodeT& t) const {
                return impl->autorefine_square_test(key, t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        /// Inplace autorefines the function using same test as for squaring.

        /// return this for chaining
        const Function<T,NDIM>& refine(bool fence = true) const {
            refine_general(autorefine_square_op(), fence);
            return *this;
        }

        /// Inplace broadens support in scaling function basis
        void broaden(const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                     bool fence = true) const {
            verify();
            reconstruct();
            impl->broaden(bc.is_periodic(), fence);
        }


        /// Clears the function as if constructed uninitialized.  Optional fence.

        /// Any underlying data will not be freed until the next global fence.
        void clear(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) {
                World& world = impl->world;
                impl.reset();
                if (fence) world.gop.fence();
            }
        }

        /// Process 0 prints a summary of all nodes in the tree (collective)
        void print_tree(std::ostream& os = std::cout) const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->print_tree(os);
        }

        /// same as print_tree() but produces JSON-formatted string
        /// @warning enclose the result in braces to make it a valid JSON object
        void print_tree_json(std::ostream& os = std::cout) const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->print_tree_json(os);
        }

        /// Process 0 prints a graphviz-formatted output of all nodes in the tree (collective)
        void print_tree_graphviz(std::ostream& os = std::cout) const {
            PROFILE_MEMBER_FUNC(Function);
            os << "digraph G {" << std::endl;
            if (impl) impl->print_tree_graphviz(os);
            os << "}" << std::endl;
        }

        /// Print a summary of the load balancing info

        /// This is serial and VERY expensive
        void print_info() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->print_info();
        }

        struct SimpleUnaryOpWrapper {
            T (*f)(T);
            SimpleUnaryOpWrapper(T (*f)(T)) : f(f) {}
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = f(*_p0));
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };

        /// Inplace unary operation on function values
        void unaryop(T (*f)(T)) {
            // Must fence here due to temporary object on stack
            // stopping us returning before complete
            this->unaryop(SimpleUnaryOpWrapper(f));
        }


        /// Inplace unary operation on function values
        template <typename opT>
        void unaryop(const opT& op, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            reconstruct();
            impl->unary_op_value_inplace(op, fence);
        }


        /// Unary operation applied inplace to the coefficients
        template <typename opT>
        void unaryop_coeff(const opT& op,
                           bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->unary_op_coeff_inplace(op, fence);
        }


        /// Unary operation applied inplace to the nodes
        template <typename opT>
        void unaryop_node(const opT& op,
                          bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->unary_op_node_inplace(op, fence);
        }




        static void doconj(const Key<NDIM>, Tensor<T>& t) {
            PROFILE_MEMBER_FUNC(Function);
            t.conj();
        }

        /// Inplace complex conjugate.  No communication except for optional fence.

        /// Returns this for chaining.  Works in either basis.
        Function<T,NDIM> conj(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            unaryop_coeff(&Function<T,NDIM>::doconj, fence);
            return *this;
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
        /// Returns this for chaining, can be in states compressed of redundant_after_merge.
        ///
        /// this and other may have different distributions and may even live in different worlds
        ///
        /// this <-- this*alpha + other*beta
        template <typename Q, typename R>
        Function<T,NDIM>& gaxpy(const T& alpha,
                                const Function<Q,NDIM>& other, const R& beta, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            other.verify();

            // operation is done either in compressed or reconstructed state
            TreeState operating_state=this->get_impl()->get_tensor_type()==TT_FULL ? compressed : reconstructed;

            TreeState thisstate=impl->get_tree_state();
            TreeState otherstate=other.get_impl()->get_tree_state();

            if (operating_state==compressed) {
                MADNESS_CHECK_THROW(thisstate==compressed, "gaxpy: this must be compressed");
                MADNESS_CHECK_THROW(otherstate==compressed, "gaxpy: other must be compressed");
                impl->gaxpy_inplace(alpha, *other.get_impl(), beta, fence);

            } else if (operating_state==reconstructed) {
                // this works both in reconstructed and redundant_after_merge states
                MADNESS_CHECK_THROW(thisstate==reconstructed or thisstate==redundant_after_merge,
                    "gaxpy: this must be reconstructed or redundant_after_merge");
                MADNESS_CHECK_THROW(otherstate==reconstructed or otherstate==redundant_after_merge,
                    "gaxpy: other must be reconstructed or redundant_after_merge");

                impl->gaxpy_inplace_reconstructed(alpha,*other.get_impl(),beta,fence);
            } else {
                MADNESS_EXCEPTION("unknown tree state",1);
            }
            return *this;
        }


        /// Inplace addition of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation.
        /// Functions don't need to be compressed, it's the caller's responsibility
        /// to choose an appropriate state with performance, usually compressed for 3d,
        /// reconstructed for 6d)
        template <typename Q>
        Function<T,NDIM>& operator+=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);

            // do this in reconstructed or compressed form
            TreeState operating_state=get_impl()->get_tensor_type()==TT_FULL ? compressed : reconstructed;
            this->change_tree_state(operating_state);
            other.change_tree_state(operating_state);

            MADNESS_ASSERT(impl->get_tree_state() == other.get_impl()->get_tree_state());
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(1.0), true);
        }


        /// Inplace subtraction of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        Function<T,NDIM>& operator-=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);
            if (NDIM<=3) {
                compress();
                other.compress();
            } else {
                reconstruct();
                other.reconstruct();
            }
            MADNESS_ASSERT(impl->get_tree_state() == other.get_impl()->get_tree_state());
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(-1.0), true);
        }


        /// Inplace scaling by a constant

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Function<T,NDIM> >::type &
        operator*=(const Q q) {
            PROFILE_MEMBER_FUNC(Function);
            scale(q,true);
            return *this;
        }


        /// Inplace squaring of function ... global comm only if not reconstructed

        /// Returns *this for chaining.
        Function<T,NDIM>& square(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_reconstructed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->square_inplace(fence);
            return *this;
        }

        /// Returns *this for chaining.
        Function<T,NDIM>& abs(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_reconstructed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->abs_inplace(fence);
            return *this;
        }

        /// Returns *this for chaining.
        Function<T,NDIM>& abs_square(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_reconstructed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->abs_square_inplace(fence);
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
            bool compressed=is_compressed() and g.is_compressed();
            bool redundant=is_redundant() and g.is_redundant();
            MADNESS_CHECK_THROW(compressed or redundant,"functions must be compressed or redundant in inner");
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();
            return impl->inner_local(*(g.get_impl()));
        }

        /// Returns local part of dot product ... throws if both not compressed
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) dot_local(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);
            MADNESS_ASSERT(is_compressed());
            MADNESS_ASSERT(g.is_compressed());
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();
            return impl->dot_local(*(g.get_impl()));
        }


        /// With this being an on-demand function, fill the MRA tree according to different criteria

        /// @param[in]  g   the function after which the MRA structure is modeled (any basis works)
        template<typename R>
        Function<T,NDIM>& fill_tree(const Function<R,NDIM>& g, bool fence=true) {
          MADNESS_ASSERT(g.is_initialized());
          MADNESS_ASSERT(is_on_demand());

          // clear what we have
          impl->get_coeffs().clear();

          //leaf_op<T,NDIM> gnode_is_leaf(g.get_impl().get());
          Leaf_op_other<T,NDIM> gnode_is_leaf(g.get_impl().get());
          impl->make_Vphi(gnode_is_leaf,fence);
          return *this;

        }

        /// With this being an on-demand function, fill the MRA tree according to different criteria

        /// @param[in]  op  the convolution operator for screening
        template<typename opT>
        Function<T,NDIM>& fill_tree(const opT& op, bool fence=true) {
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          Specialbox_op<T,NDIM> sbox;
          Leaf_op<T,NDIM,opT,Specialbox_op<T,NDIM> > leaf_op(this->get_impl().get(),&op,sbox);
          impl ->make_Vphi(leaf_op,fence);
          return *this;
        }

        /// With this being an on-demand function, fill the MRA tree according to different criteria
        Function<T,NDIM>& fill_tree(bool fence=true) {
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          Leaf_op<T,NDIM,SeparatedConvolution<double,NDIM>,Specialbox_op<T,NDIM> > leaf_op(this->get_impl().get());
          impl->make_Vphi(leaf_op,fence);
          return *this;
        }

        /// Special refinement on 6D boxes where the electrons come close (meet)
        /// @param[in]  op  the convolution operator for screening
        template<typename opT>
        Function<T,NDIM>& fill_cuspy_tree(const opT& op,const bool fence=true){
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          ElectronCuspyBox_op<T,NDIM> sbox;

          Leaf_op<T,NDIM,opT,ElectronCuspyBox_op<T,NDIM> > leaf_op(this->get_impl().get(),&op,sbox);
          impl ->make_Vphi(leaf_op,fence);

          return *this;
        }

        /// Special refinement on 6D boxes where the electrons come close (meet)
        Function<T,NDIM>& fill_cuspy_tree(const bool fence=true){
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          ElectronCuspyBox_op<T,NDIM> sbox;

          Leaf_op<T,NDIM,SeparatedConvolution<double,NDIM>,ElectronCuspyBox_op<T,NDIM> > leaf_op(this->get_impl().get(),sbox);
          impl ->make_Vphi(leaf_op,fence);

          return *this;
        }

        /// Special refinement on 6D boxes for the nuclear potentials (regularized with cusp, non-regularized with singularity)
        /// @param[in]  op  the convolution operator for screening
        template<typename opT>
        Function<T,NDIM>& fill_nuclear_cuspy_tree(const opT& op,const size_t particle,const bool fence=true){
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          NuclearCuspyBox_op<T,NDIM> sbox(particle);

          Leaf_op<T,NDIM,opT,NuclearCuspyBox_op<T,NDIM> > leaf_op(this->get_impl().get(),&op,sbox);
          impl ->make_Vphi(leaf_op,fence);

          return *this;
        }

        /// Special refinement on 6D boxes for the nuclear potentials (regularized with cusp, non-regularized with singularity)
        Function<T,NDIM>& fill_nuclear_cuspy_tree(const size_t particle,const bool fence=true){
          MADNESS_ASSERT(is_on_demand());
          // clear what we have
          impl->get_coeffs().clear();
          NuclearCuspyBox_op<T,NDIM> sbox(particle);

          Leaf_op<T,NDIM,SeparatedConvolution<double,NDIM>,NuclearCuspyBox_op<T,NDIM> > leaf_op(this->get_impl().get(),sbox);
          impl ->make_Vphi(leaf_op,fence);

          return *this;
        }

        /// perform the hartree product of f*g, invoked by result
        template<size_t LDIM, size_t KDIM, typename opT>
        void do_hartree_product(const std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> left,
                                const std::vector<std::shared_ptr<FunctionImpl<T,KDIM>>> right,
                                const opT* op) {

            // get the right leaf operator
            hartree_convolute_leaf_op<T,KDIM+LDIM,LDIM,opT> leaf_op(impl.get(),left,op);
            impl->hartree_product(left,right,leaf_op,true);
            impl->finalize_sum();
//            this->truncate();

        }

        /// perform the hartree product of f*g, invoked by result
        template<size_t LDIM, size_t KDIM>
        void do_hartree_product(const std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> left,
                                const std::vector<std::shared_ptr<FunctionImpl<T,KDIM>>> right) {

//            hartree_leaf_op<T,KDIM+LDIM> leaf_op(impl.get(),cdata.s0);
            hartree_leaf_op<T,KDIM+LDIM> leaf_op(impl.get(),k());
            impl->hartree_product(left,right,leaf_op,true);
            impl->finalize_sum();
//            this->truncate();

        }

        /// Returns the inner product

        /// Not efficient for computing multiple inner products
        /// @param[in]  g   Function, optionally on-demand
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);

            // fast return if possible
            if (not this->is_initialized()) return 0.0;
            if (not g.is_initialized()) return 0.0;

            // if this and g are the same, use norm2()
            if constexpr (std::is_same_v<T,R>) {
              if (this->get_impl() == g.get_impl()) {
                TreeState state = this->get_impl()->get_tree_state();
                if (not(state == reconstructed or state == compressed))
                  change_tree_state(reconstructed);
                double norm = this->norm2();
                return norm * norm;
              }
            }

            // do it case-by-case
            if constexpr (std::is_same_v<R,T>) {
              if (this->is_on_demand())
                return g.inner_on_demand(*this);
              if (g.is_on_demand())
                return this->inner_on_demand(g);
            }

            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();

            // compute in compressed form if compression is fast, otherwise in redundant form
            TreeState operating_state=get_impl()->get_tensor_type()==TT_FULL ? compressed : redundant;

            change_tree_state(operating_state,false);
            g.change_tree_state(operating_state,false);
            impl->world.gop.fence();

            TENSOR_RESULT_TYPE(T,R) local = impl->inner_local(*g.get_impl());
            impl->world.gop.sum(local);
            impl->world.gop.fence();

            // restore state -- no need for this
            // change_tree_state(state,false);
            // g.change_tree_state(gstate,false);
            // impl->world.gop.fence();

            return local;
        }

        /// Return the local part of inner product with external function ... no communication.
        /// If you are going to be doing a bunch of inner_ext calls, set
        /// keep_redundant to true and then manually undo_redundant when you
        /// are finished.
        /// @param[in] f Pointer to function of type T that take coordT arguments. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @param[in] keep_redundant boolean switch to turn on/off undo_redundant
        /// @return Returns local part of the inner product, i.e. over the domain of all function nodes on this compute node.
        T inner_ext_local(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f, const bool leaf_refine=true, const bool keep_redundant=false) const {
            PROFILE_MEMBER_FUNC(Function);
            change_tree_state(redundant);
            T local = impl->inner_ext_local(f, leaf_refine);
            if (not keep_redundant) change_tree_state(reconstructed);
            return local;
        }

        /// Return the inner product with external function ... requires communication.
        /// If you are going to be doing a bunch of inner_ext calls, set
        /// keep_redundant to true and then manually undo_redundant when you
        /// are finished.
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @param[in] keep_redundant boolean switch to turn on/off undo_redundant
        /// @return Returns the inner product
        T inner_ext(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f, const bool leaf_refine=true, const bool keep_redundant=false) const {
            PROFILE_MEMBER_FUNC(Function);
            change_tree_state(redundant);
            T local = impl->inner_ext_local(f, leaf_refine);
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            if (not keep_redundant) change_tree_state(reconstructed);
            return local;
        }

        /// Return the inner product with external function ... requires communication.
        /// If you are going to be doing a bunch of inner_ext calls, set
        /// keep_redundant to true and then manually undo_redundant when you
        /// are finished.
        /// @param[in] f Reference to FunctionFunctorInterface. This is the externally provided function
        /// @param[in] leaf_refine boolean switch to turn on/off refinement past leaf nodes
        /// @return Returns the inner product
        T inner_adaptive(const std::shared_ptr< FunctionFunctorInterface<T,NDIM> > f,
                const bool leaf_refine=true) const {
            PROFILE_MEMBER_FUNC(Function);
            reconstruct();
            T local = impl->inner_adaptive_local(f, leaf_refine);
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return local;
        }

        /// Return the local part of gaxpy with external function, this*alpha + f*beta ... no communication.
        /// @param[in] alpha prefactor for this Function
        /// @param[in] f Pointer to function of type T that take coordT arguments. This is the externally provided function
        /// @param[in] beta prefactor for f
        template <typename L>
        void gaxpy_ext(const Function<L,NDIM>& left, T (*f)(const coordT&), T alpha, T beta, double tol, bool fence=true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!left.is_reconstructed()) left.reconstruct();
            impl->gaxpy_ext(left.get_impl().get(), f, alpha, beta, tol, fence);
        }

        /// Returns the inner product for one on-demand function

        /// It does work, but it might not give you the precision you expect.
        /// The assumption is that the function g returns proper sum
        /// coefficients on the MRA tree of this. This might not be the case if
        /// g is constructed with an implicit multiplication, e.g.
        ///  result = <this|g>,   with g = 1/r12 | gg>
        /// @param[in]  g	on-demand function
        template<typename R>
        TENSOR_RESULT_TYPE(T, R) inner_on_demand(const Function<R, NDIM>& g) const {
            MADNESS_ASSERT(g.is_on_demand() and (not this->is_on_demand()));

            constexpr std::size_t LDIM=std::max(NDIM/2,std::size_t(1));
            auto func=dynamic_cast<CompositeFunctorInterface<T,NDIM,LDIM>* >(g.get_impl()->get_functor().get());
            MADNESS_ASSERT(func);
            func->make_redundant(true);
            func->replicate_low_dim_functions(true);
            this->reconstruct();        // if this == &g we don't need g to be redundant

            if (VERIFY_TREE) verify_tree();

            TENSOR_RESULT_TYPE(T, R) local = impl->inner_local_on_demand(*g.get_impl());
            impl->world.gop.sum(local);
            impl->world.gop.fence();

            return local;
        }

        /// project this on the low-dim function g: h(x) = <f(x,y) | g(y)>

        /// @param[in]  g   low-dim function
        /// @param[in]  dim over which dimensions to be integrated: 0..LDIM-1 or LDIM..NDIM-1
        /// @return     new function of dimension NDIM-LDIM
        template <typename R, size_t LDIM>
        Function<TENSOR_RESULT_TYPE(T,R),NDIM-LDIM> project_out(const Function<R,LDIM>& g, const int dim) const {
            if (NDIM<=LDIM) MADNESS_EXCEPTION("confused dimensions in project_out?",1);
            MADNESS_CHECK_THROW(dim==0 or dim==1,"dim must be 0 or 1 in project_out");
            verify();
            typedef TENSOR_RESULT_TYPE(T,R) resultT;
            static const size_t KDIM=NDIM-LDIM;

            FunctionFactory<resultT,KDIM> factory=FunctionFactory<resultT,KDIM>(world())
                    .k(g.k()).thresh(g.thresh());
            Function<resultT,KDIM> result=factory;      // no empty() here!

            change_tree_state(reconstructed,false);
            g.change_tree_state(redundant,false);
            world().gop.fence();
            this->get_impl()->project_out(result.get_impl().get(),g.get_impl().get(),dim,true);
//            result.get_impl()->project_out2(this->get_impl().get(),gimpl,dim);
            result.world().gop.fence();
            g.change_tree_state(reconstructed,false);
            result.get_impl()->trickle_down(false);
            result.get_impl()->set_tree_state(reconstructed);
            result.world().gop.fence();
            return result;
        }

        Function<T,NDIM/2> dirac_convolution(const bool fence=true) const {
            constexpr std::size_t LDIM=NDIM/2;
            MADNESS_CHECK_THROW(NDIM==2*LDIM,"NDIM must be even");
//        	// this will be the result function
        	FunctionFactory<T,LDIM> factory=FunctionFactory<T,LDIM>(world()).k(this->k());
        	Function<T,LDIM> f = factory;
        	if(!is_reconstructed()) this->reconstruct();
        	this->get_impl()->do_dirac_convolution(f.get_impl().get(),fence);
        	return f;
        }

        /// Replaces this function with one loaded from an archive using the default processor map

        /// Archive can be sequential or parallel.
        ///
        /// The & operator for serializing will only work with parallel archives.
        template <typename Archive>
        void load(World& world, Archive& ar) {
            PROFILE_MEMBER_FUNC(Function);
            // Type checking since we are probably circumventing the archive's own type checking
            long magic = 0l, id = 0l, ndim = 0l, k = 0l;
            ar & magic & id & ndim & k;
            MADNESS_ASSERT(magic == 7776768); // Mellow Mushroom Pizza tel.# in Knoxville
            MADNESS_ASSERT(id == TensorTypeData<T>::id);
            MADNESS_ASSERT(ndim == NDIM);

            impl.reset(new implT(FunctionFactory<T,NDIM>(world).k(k).empty()));

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

        /// change the tensor type of the coefficients in the FunctionNode

        /// @param[in]  targs   target tensor arguments (threshold and full/low rank)
        void change_tensor_type(const TensorArgs& targs, bool fence=true) {
            if (not impl) return;
            impl->change_tensor_type1(targs,fence);
        }


        /// This is replaced with left*right ...  private
        template <typename Q, typename opT>
        Function<typename opT::resultT,NDIM>& unary_op_coeffs(const Function<Q,NDIM>& func,
                const opT& op, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            func.verify();
            MADNESS_ASSERT(func.is_reconstructed());
            if (VERIFY_TREE) func.verify_tree();
            impl.reset(new implT(*func.get_impl(), func.get_pmap(), false));
            impl->unaryXX(func.get_impl().get(), op, fence);
            return *this;
        }

        /// Returns vector of FunctionImpl pointers corresponding to vector of functions
        template <typename Q, std::size_t D>
        static std::vector< std::shared_ptr< FunctionImpl<Q,D> > > vimpl(const std::vector< Function<Q,D> >& v) {
            PROFILE_MEMBER_FUNC(Function);
            std::vector< std::shared_ptr< FunctionImpl<Q,D> > > r(v.size());
            for (unsigned int i=0; i<v.size(); ++i) r[i] = v[i].get_impl();
            return r;
        }

        /// This is replaced with op(vector of functions) ... private
        template <typename opT>
        Function<T,NDIM>& multiop_values(const opT& op, const std::vector< Function<T,NDIM> >& vf) {
            std::vector<implT*> v(vf.size(),NULL);
            for (unsigned int i=0; i<v.size(); ++i) {
                if (vf[i].is_initialized()) v[i] = vf[i].get_impl().get();
            }
            impl->multiop_values(op, v);
            world().gop.fence();
            if (VERIFY_TREE) verify_tree();

            return *this;
        }

        /// apply op on the input vector yielding an output vector of functions

        /// (*this) is just a dummy Function to be able to call internal methods in FuncImpl
        /// @param[in]  op   the operator working on vin
        /// @param[in]  vin  vector of input Functions
        /// @param[out] vout vector of output Functions vout = op(vin)
        template <typename opT>
        void multi_to_multi_op_values(const opT& op,
                const std::vector< Function<T,NDIM> >& vin,
                std::vector< Function<T,NDIM> >& vout,
                const bool fence=true) {
            std::vector<implT*> vimplin(vin.size(),NULL);
            for (unsigned int i=0; i<vin.size(); ++i) {
                if (vin[i].is_initialized()) vimplin[i] = vin[i].get_impl().get();
            }
            std::vector<implT*> vimplout(vout.size(),NULL);
            for (unsigned int i=0; i<vout.size(); ++i) {
                if (vout[i].is_initialized()) vimplout[i] = vout[i].get_impl().get();
            }

            impl->multi_to_multi_op_values(op, vimplin, vimplout, fence);
            if (VERIFY_TREE) verify_tree();

        }


        /// Multiplication of function * vector of functions using recursive algorithm of mulxx
        template <typename L, typename R>
        void vmulXX(const Function<L,NDIM>& left,
                    const std::vector< Function<R,NDIM> >& right,
                    std::vector< Function<T,NDIM> >& result,
                    double tol,
                    bool fence) {
            PROFILE_MEMBER_FUNC(Function);

            std::vector<FunctionImpl<T,NDIM>*> vresult(right.size());
            std::vector<const FunctionImpl<R,NDIM>*> vright(right.size());
            for (unsigned int i=0; i<right.size(); ++i) {
                result[i].set_impl(left,false);
                vresult[i] = result[i].impl.get();
                vright[i] = right[i].get_impl().get();
            }

            left.world().gop.fence(); // Is this still essential?  Yes.
            vresult[0]->mulXXvec(left.get_impl().get(), vright, vresult, tol, fence);
        }

        /// Same as \c operator* but with optional fence and no automatic reconstruction

        /// f or g are on-demand functions
        template<typename L, typename R>
        void mul_on_demand(const Function<L,NDIM>& f, const Function<R,NDIM>& g, bool fence=true) {
            const FunctionImpl<L,NDIM>* fimpl=f.get_impl().get();
            const FunctionImpl<R,NDIM>* gimpl=g.get_impl().get();
            if (fimpl->is_on_demand() and gimpl->is_on_demand()) {
                MADNESS_EXCEPTION("can't multiply two on-demand functions",1);
            }

            if (fimpl->is_on_demand()) {
                leaf_op<R,NDIM> leaf_op1(gimpl);
                impl->multiply(leaf_op1,gimpl,fimpl,fence);
            } else {
                leaf_op<L,NDIM> leaf_op1(fimpl);
                impl->multiply(leaf_op1,fimpl,gimpl,fence);
            }
        }

        /// sparse transformation of a vector of functions ... private
        template <typename R, typename Q>
        void vtransform(const std::vector< Function<R,NDIM> >& v,
                        const Tensor<Q>& c,
                        std::vector< Function<T,NDIM> >& vresult,
                        double tol,
                        bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            vresult[0].impl->vtransform(vimpl(v), c, vimpl(vresult), tol, fence);
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
            impl.reset(new implT(*left.get_impl(), left.get_pmap(), false));
            impl->gaxpy(alpha,*left.get_impl(),beta,*right.get_impl(),fence);
            return *this;
        }

        /// This is replaced with mapdim(f) ...  private
        Function<T,NDIM>& mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            for (std::size_t i=0; i<NDIM; ++i) MADNESS_ASSERT(map[i]>=0 && static_cast<std::size_t>(map[i])<NDIM);
            impl.reset(new implT(*f.impl, f.get_pmap(), false));
            impl->mapdim(*f.impl,map,fence);
            return *this;
        }

        /// This is replaced with mirror(f) ...  private

        /// similar to mapdim, but maps from x to -x, y to -y, and so on
        /// Example: mirror a 3d function on the xy plane: mirror={1,1,-1}
        /// @param[in]	mirror	array of -1 and 1, corresponding to mirror or not
        Function<T,NDIM>& mirror(const Function<T,NDIM>& f, const std::vector<long>& mirrormap, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            for (std::size_t i=0; i<NDIM; ++i) MADNESS_ASSERT((mirrormap[i]==1) or (mirrormap[i]==-1));
            impl.reset(new implT(*f.impl, f.get_pmap(), false));
            impl->mirror(*f.impl,mirrormap,fence);
            return *this;
        }

        /// This is replaced with mirror(map(f)) ...  private

        /// first map then mirror!
        /// mirror is similar to mapdim, but maps from x to -x, y to -y, and so on
        /// Example: mirror a 3d function on the xy plane: mirror={1,1,-1}
        /// Example: c4 rotation of a 3d function around the z axis:
        /// 	x->y, y->-x, z->z: map(1,0,2); mirror(-1,1,1)
        /// @param[in]	map		array holding dimensions
        /// @param[in]	mirror	array of -1 and 1, corresponding to mirror or not
        Function<T,NDIM>& map_and_mirror(const Function<T,NDIM>& f,
        		const std::vector<long>& map, const std::vector<long>& mirror,
				bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            for (std::size_t i=0; i<mirror.size(); ++i) MADNESS_ASSERT((mirror[i]==1) or (mirror[i]==-1));
            for (std::size_t i=0; i<map.size(); ++i) MADNESS_ASSERT(map[i]>=0 && static_cast<std::size_t>(map[i])<NDIM);

            impl.reset(new implT(*f.impl, f.get_pmap(), false));
            impl->map_and_mirror(*f.impl,map,mirror,fence);
            return *this;
        }


        /// check symmetry of a function by computing the 2nd derivative
        double check_symmetry() const {

            change_tree_state(redundant);
            if (VERIFY_TREE) verify_tree();
            double local = impl->check_symmetry_local();
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            double asy=sqrt(local);
            if (this->world().rank()==0) print("asymmetry wrt particle",asy);
            change_tree_state(reconstructed);
            return asy;
        }

        /// reduce the rank of the coefficient tensors
        Function<T,NDIM>& reduce_rank(const double thresh=0.0, const bool fence=true) {
            verify();
            double thresh1= (thresh==0.0) ? impl->get_tensor_args().thresh : thresh;
            impl->reduce_rank(thresh1,fence);
            return *this;
        }

        /// remove all nodes with level higher than n
        Function<T,NDIM>& chop_at_level(const int n, const bool fence=true) {
            verify();
            change_tree_state(redundant);
            impl->chop_at_level(n,true);
            change_tree_state(reconstructed);
            return *this;
        }
    };

//    template <typename T, typename opT, std::size_t NDIM>
    template <typename T, typename opT, std::size_t NDIM>
    Function<T,NDIM> multiop_values(const opT& op, const std::vector< Function<T,NDIM> >& vf) {
        Function<T,NDIM> r;
        r.set_impl(vf[0], false);
        r.multiop_values(op, vf);
        return r;
    }

    /// Returns new function equal to alpha*f(x) with optional fence
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Q alpha, const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        f.verify();
        if (VERIFY_TREE) f.verify_tree();
        Function<TENSOR_RESULT_TYPE(Q,T),NDIM> result;
        result.set_impl(f, false);
        result.get_impl()->scale_oop(alpha,*f.get_impl(),fence);
        return result;
    }


    /// Returns new function equal to f(x)*alpha with optional fence
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Function<T,NDIM>& f, const Q alpha, bool fence=true) {
        PROFILE_FUNC;
        return mul(alpha,f,fence);
    }


    /// Returns new function equal to f(x)*alpha

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Function<T,NDIM>& f, const Q alpha) {
        return mul(alpha, f, true);
    }

    /// Returns new function equal to alpha*f(x)

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Q alpha, const Function<T,NDIM>& f) {
        return mul(alpha, f, true);
    }

    /// Sparse multiplication --- left and right must be reconstructed and if tol!=0 have tree of norms already created
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul_sparse(const Function<L,NDIM>& left, const Function<R,NDIM>& right, double tol, bool fence=true) {
        PROFILE_FUNC;
        left.verify();
        right.verify();
        MADNESS_ASSERT(left.is_reconstructed() and right.is_reconstructed());
        if (VERIFY_TREE) left.verify_tree();
        if (VERIFY_TREE) right.verify_tree();

        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        result.set_impl(left, false);
        result.get_impl()->mulXX(left.get_impl().get(), right.get_impl().get(), tol, fence);
        return result;
    }

    /// Same as \c operator* but with optional fence and no automatic reconstruction
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return mul_sparse(left,right,0.0,fence);
    }

    /// Generate new function = op(left,right) where op acts on the function values
    template <typename L, typename R, typename opT, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    binary_op(const Function<L,NDIM>& left, const Function<R,NDIM>& right, const opT& op, bool fence=true) {
        PROFILE_FUNC;
        if (!left.is_reconstructed()) left.reconstruct();
        if (!right.is_reconstructed()) right.reconstruct();

        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        result.set_impl(left, false);
        result.get_impl()->binaryXX(left.get_impl().get(), right.get_impl().get(), op, fence);
        return result;
    }

    /// Out of place application of unary operation to function values with optional fence
    template <typename Q, typename opT, std::size_t NDIM>
    Function<typename opT::resultT, NDIM>
    unary_op(const Function<Q,NDIM>& func, const opT& op, bool fence=true) {
        if (!func.is_reconstructed()) func.reconstruct();
        Function<typename opT::resultT, NDIM> result;
        if (VERIFY_TREE) func.verify_tree();
        result.set_impl(func, false);
        result.get_impl()->unaryXXvalues(func.get_impl().get(), op, fence);
        return result;
    }


    /// Out of place application of unary operation to scaling function coefficients with optional fence
    template <typename Q, typename opT, std::size_t NDIM>
    Function<typename opT::resultT, NDIM>
    unary_op_coeffs(const Function<Q,NDIM>& func, const opT& op, bool fence=true) {
        if (!func.is_reconstructed()) func.reconstruct();
        Function<typename opT::resultT, NDIM> result;
        return result.unary_op_coeffs(func,op,fence);
    }

    /// Use the vmra/mul(...) interface instead

    /// This so that we don't have to have friend functions in a different header.
    ///
    /// If using sparsity (tol != 0) you must have created the tree of norms
    /// already for both left and right.
    template <typename L, typename R, std::size_t D>
    std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> >
    vmulXX(const Function<L,D>& left, const std::vector< Function<R,D> >& vright, double tol, bool fence=true) {
        if (vright.size() == 0) return std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> >();
        std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> > vresult(vright.size());
        vresult[0].vmulXX(left, vright, vresult, tol, fence);
        return vresult;
    }

    /// Multiplies two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation but also
    /// enables us to automatically reconstruct the input functions as required.
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator*(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (!left.is_reconstructed())  left.reconstruct();
        if (!right.is_reconstructed()) right.reconstruct();
        MADNESS_ASSERT(not (left.is_on_demand() or right.is_on_demand()));
        return mul(left,right,true);
    }

    /// Performs a Hartree/outer product on the two given low-dimensional function vectors

    /// @return   result(x,y) = \sum_i f_i(x) g_i(y)
    template<typename T, std::size_t KDIM, std::size_t LDIM>
    Function<T,KDIM+LDIM>
    hartree_product(const std::vector<Function<T,KDIM>>& left, const std::vector<Function<T,LDIM>>& right) {

        MADNESS_CHECK_THROW(left.size()==right.size(), "hartree_product: left and right must have same size");
        if (left.size()==0) return Function<T,KDIM+LDIM>();

        const double thresh=FunctionDefaults<KDIM+LDIM>::get_thresh();

        FunctionFactory<T,KDIM+LDIM> factory=FunctionFactory<T,KDIM+LDIM>(left.front().world())
                .k(left.front().k()).thresh(thresh);
        Function<T,KDIM+LDIM> result=factory.empty();

        // some prep work
        change_tree_state(left,nonstandard_with_leaves);
        change_tree_state(right,nonstandard_with_leaves);
        std::vector<std::shared_ptr<FunctionImpl<T,KDIM>>> vleft=get_impl(left);
        std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> vright=get_impl(right);

        result.do_hartree_product(vleft,vright);

        return result;

    }

    /// Performs a Hartree product on the two given low-dimensional functions
    template<typename T, std::size_t KDIM, std::size_t LDIM>
    Function<T,KDIM+LDIM>
    hartree_product(const Function<T,KDIM>& left2, const Function<T,LDIM>& right2) {
        typedef std::vector<Function<T,KDIM>> vector;
        return hartree_product(vector({left2}),vector({right2}));
    }

    /// Performs a Hartree product on the two given low-dimensional functions
    template<typename T, std::size_t KDIM, std::size_t LDIM, typename opT>
    Function<T,KDIM+LDIM>
    hartree_product(const Function<T,KDIM>& left2, const Function<T,LDIM>& right2,
            const opT& op) {

        // we need both sum and difference coeffs for error estimation
    	Function<T,KDIM>& left = const_cast< Function<T,KDIM>& >(left2);
    	Function<T,LDIM>& right = const_cast< Function<T,LDIM>& >(right2);

    	const double thresh=FunctionDefaults<KDIM+LDIM>::get_thresh();

    	FunctionFactory<T,KDIM+LDIM> factory=FunctionFactory<T,KDIM+LDIM>(left.world())
    			.k(left.k()).thresh(thresh);
        Function<T,KDIM+LDIM> result=factory.empty();

    	if (result.world().rank()==0) {
            print("incomplete FunctionFactory in Function::hartree_product");
            print("thresh: ", thresh);
        }
        bool same=(left2.get_impl()==right2.get_impl());

        // some prep work
        left.make_nonstandard(true, true);
        right.make_nonstandard(true, true);

        std::vector<std::shared_ptr<FunctionImpl<T,KDIM>>> vleft;
        std::vector<std::shared_ptr<FunctionImpl<T,LDIM>>> vright;
        vleft.push_back(left.get_impl());
        vright.push_back(right.get_impl());
        result.do_hartree_product(vleft,right,&op);

        left.standard(false);
        if (not same) right.standard(false);
        left2.world().gop.fence();

        return result;
    }

    /// adds beta*right only left:  alpha*left + beta*right optional fence and no automatic compression

    /// left and right might live in different worlds, the accumulation is non-blocking
    template <typename L, typename R,std::size_t NDIM>
    void
    gaxpy(TENSOR_RESULT_TYPE(L,R) alpha, Function<L,NDIM>& left,
              TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,NDIM>& right, bool fence=true) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        left.gaxpy(alpha, right, beta, fence);
    }

    /// Returns new function alpha*left + beta*right optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    gaxpy_oop(TENSOR_RESULT_TYPE(L,R) alpha, const Function<L,NDIM>& left,
              TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,NDIM>& right, bool fence=true) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        return result.gaxpy_oop(alpha, left, beta, right, fence);
    }

    /// Same as \c operator+ but with optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    add(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)(1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(1.0), right, fence);
    }


    /// Returns new function alpha*left + beta*right optional fence, having both addends reconstructed
    template<typename T, std::size_t NDIM>
    Function<T,NDIM> gaxpy_oop_reconstructed(const double alpha, const Function<T,NDIM>& left,
            const double beta, const Function<T,NDIM>& right, const bool fence=true) {
        Function<T,NDIM> result;
        result.set_impl(right,false);

        MADNESS_ASSERT(left.is_reconstructed());
        MADNESS_ASSERT(right.is_reconstructed());
        result.get_impl()->gaxpy_oop_reconstructed(alpha,*left.get_impl(),beta,*right.get_impl(),fence);
        return result;

    }

    /// Adds two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator+(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (VERIFY_TREE) left.verify_tree();
        if (VERIFY_TREE) right.verify_tree();

        // no compression for high-dimensional functions
        if (NDIM==6) {
            left.reconstruct();
            right.reconstruct();
            return gaxpy_oop_reconstructed(1.0,left,1.0,right,true);
        } else {
            if (!left.is_compressed())  left.compress();
            if (!right.is_compressed()) right.compress();
            return add(left,right,true);
        }
    }

    /// Same as \c operator- but with optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    sub(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)(1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(-1.0), right, fence);
    }


    /// Subtracts two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator-(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        PROFILE_FUNC;
        // no compression for high-dimensional functions
        if (NDIM==6) {
            left.reconstruct();
            right.reconstruct();
            return gaxpy_oop_reconstructed(1.0,left,-1.0,right,true);
        } else {
            if (!left.is_compressed())  left.compress();
            if (!right.is_compressed()) right.compress();
            return sub(left,right,true);
        }
    }

    /// Create a new copy of the function with different distribution and optional fence

    /// Works in either basis.  Different distributions imply
    /// asynchronous communication and the optional fence is
    /// collective.
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f,
                          const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
                          bool fence = true) {
        PROFILE_FUNC;
        f.verify();
        Function<T,NDIM> result;
        typedef FunctionImpl<T,NDIM> implT;
        result.set_impl(std::shared_ptr<implT>(new implT(*f.get_impl(), pmap, false)));
        result.get_impl()->copy_coeffs(*f.get_impl(), fence);
        if (VERIFY_TREE) result.verify_tree();
        return result;
    }

    /// Create a new copy of the function with the same distribution and optional fence
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        return copy(f, f.get_pmap(), fence);
    }

    /// Create a new copy of function f living in world (might differ from f.world)

    /// uses the default processor map of world
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> copy(World& world, const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        typedef FunctionImpl<T,NDIM> implT;
        auto pmap=FunctionDefaults<NDIM>::get_pmap();

        // create a new function with pmap distribution, same parameters as f, but no coeffs
        Function<T,NDIM> result;
        result.set_impl(std::make_shared<implT>(world,*f.get_impl(), pmap, false));
        // copy f's coefficients to result
        result.get_impl()->copy_coeffs(*f.get_impl(), fence);
        return result;
    }

    /// Type conversion implies a deep copy.  No communication except for optional fence.

    /// Works in either basis but any loss of precision may result in different errors
    /// in applied in a different basis.
    ///
    /// The new function is formed with the options from the default constructor.
    ///
    /// There is no automatic type conversion since this is generally a rather dangerous
    /// thing and because there would be no way to make the fence optional.
    template <typename T, typename Q, std::size_t NDIM>
    Function<Q,NDIM> convert(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        f.verify();
        Function<Q,NDIM> result;
        result.set_impl(f, false);
        result.get_impl()->copy_coeffs(*f.get_impl(), fence);
        return result;
    }


    /// Return the complex conjugate of the input function with the same distribution and optional fence

    /// !!! The fence is actually not optional in the current implementation !!!
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> conj(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);
        return result.conj(fence);
    }

    /// Apply operator on a hartree product of two low-dimensional functions

    /// Supposed to be something like result= G( f(1)*f(2))
    /// the hartree product is never constructed explicitly, but its coeffs are
    /// constructed on the fly and processed immediately.
    /// @param[in]	op	the operator
    /// @param[in]	f1	function of particle 1
    /// @param[in]	f2	function of particle 2
    /// @param[in]	fence if we shall fence
    /// @return		a function of dimension NDIM=LDIM+LDIM
    template <typename opT, typename T, std::size_t LDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,T), LDIM+LDIM>
    apply(const opT& op, const std::vector<Function<T,LDIM>>& f1, const std::vector<Function<T,LDIM>>& f2, bool fence=true) {

        World& world=f1.front().world();

        typedef TENSOR_RESULT_TYPE(T,typename opT::opT) resultT;
        typedef std::vector<Function<T,LDIM>> vecfuncL;

    	vecfuncL& ff1 = const_cast< vecfuncL& >(f1);
    	vecfuncL& ff2 = const_cast< vecfuncL& >(f2);

    	bool same=(ff1[0].get_impl()==ff2[0].get_impl());

        reconstruct(world,f1,false);
        reconstruct(world,f2,false);
        world.gop.fence();
        // keep the leaves! They are assumed to be there later
        // even for modified op we need NS form for the hartree_leaf_op
        for (auto& f : f1) f.make_nonstandard(true,false);
        for (auto& f : f2) f.make_nonstandard(true,false);
        world.gop.fence();


        FunctionFactory<T,LDIM+LDIM> factory=FunctionFactory<resultT,LDIM+LDIM>(world)
                .k(f1.front().k()).thresh(FunctionDefaults<LDIM+LDIM>::get_thresh());
    	Function<resultT,LDIM+LDIM> result=factory.empty().fence();

    	result.get_impl()->reset_timer();
    	op.reset_timer();

		// will fence here
        for (size_t i=0; i<f1.size(); ++i)
            result.get_impl()->recursive_apply(op, f1[i].get_impl().get(),f2[i].get_impl().get(),false);
        world.gop.fence();

        if (op.print_timings) {
            result.get_impl()->print_timer();
            op.print_timer();
        }

		result.get_impl()->finalize_apply();	// need fence before reconstruct

        if (op.modified()) {
            result.get_impl()->trickle_down(true);
        } else {
    		result.get_impl()->reconstruct(true);
        }
    	standard(world,ff1,false);
    	if (not same) standard(world,ff2,false);

		return result;
    }


    /// Apply operator ONLY in non-standard form - required other steps missing !!
    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply_only(const opT& op, const Function<R,NDIM>& f, bool fence=true) {
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;

        constexpr std::size_t OPDIM=opT::opdim;
        constexpr bool low_dim=(OPDIM*2==NDIM);     // apply on some dimensions only

        // specialized version for 3D
        if (NDIM <= 3 and (not low_dim)) {
            result.set_impl(f, false);
            result.get_impl()->apply(op, *f.get_impl(), fence);

        } else {        // general version for higher dimension
	  //bool print_timings=false;
            Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> r1;

            result.set_impl(f, false);
            r1.set_impl(f, false);

            result.get_impl()->reset_timer();
            op.reset_timer();

            result.get_impl()->apply_source_driven(op, *f.get_impl(), fence);

            // recursive_apply is about 20% faster than apply_source_driven
            //result.get_impl()->recursive_apply(op, f.get_impl().get(),
            //        r1.get_impl().get(),true);          // will fence here

        }

        return result;
    }

    /// Apply operator in non-standard form

    /// Returns a new function with the same distribution
    ///
    /// !!! For the moment does NOT respect fence option ... always fences
    /// if the operator acts on one particle only the result will be sorted as
    /// g.particle=1:     g(f) = \int g(x,x') f(x',y) dx' = result(x,y)
    /// g.particle=2:     g(f) = \int g(y,y') f(x,y') dy' = result(x,y)
    /// for the second case it will notably *not* be as it is implemented in the partial inner product!
    /// g.particle=2                          g(f) = result(x,y)
    ///                 inner(g(y,y'),f(x,y'),1,1) = result(y,x)
    /// also note the confusion with the counting of the particles/integration variables
    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply(const opT& op, const Function<R,NDIM>& f, bool fence=true) {

    	typedef TENSOR_RESULT_TYPE(typename opT::opT,R) resultT;
    	Function<R,NDIM>& ff = const_cast< Function<R,NDIM>& >(f);
    	Function<resultT, NDIM> result;

		MADNESS_ASSERT(not f.is_on_demand());
		bool print_timings=op.print_timings;

    	if (VERIFY_TREE) ff.verify_tree();
    	ff.reconstruct();
        if (print_timings) ff.print_size("ff in apply after reconstruct");

    	if (op.modified()) {

            ff.change_tree_state(redundant);
//    	    ff.get_impl()->make_redundant(true);
            result = apply_only(op, ff, fence);
            ff.get_impl()->undo_redundant(false);
            result.get_impl()->trickle_down(true);

    	} else {

    		// saves the standard() step, which is very expensive in 6D
//    		Function<R,NDIM> fff=copy(ff);
    		Function<R,NDIM> fff=(ff);
            fff.make_nonstandard(op.doleaves, true);
            if (print_timings) fff.print_size("ff in apply after make_nonstandard");
            if ((print_timings) and (f.world().rank()==0)) {
                fff.get_impl()->timer_filter.print("filter");
                fff.get_impl()->timer_compress_svd.print("compress_svd");
            }
            result = apply_only(op, fff, fence);
            result.get_impl()->set_tree_state(nonstandard_after_apply);
        	ff.world().gop.fence();
            if (print_timings) result.print_size("result after apply_only");

        	// svd-tensors need some post-processing
        	if (result.get_impl()->get_tensor_type()==TT_2D) {
            	double elapsed=result.get_impl()->finalize_apply();
                if (print_timings) printf("time in finalize_apply        %8.2f\n",elapsed);
			}
			if (print_timings) {
				result.get_impl()->print_timer();
				op.print_timer();
			}

            result.get_impl()->reconstruct(true);

//            fff.clear();
            if (op.destructive()) {
            	ff.world().gop.fence();
            	ff.clear();
            } else {
            	// ff.standard();
            	ff.reconstruct();
            }

    	}
        if (print_timings) result.print_size("result after reconstruction");
        return result;
    }


    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply_1d_realspace_push(const opT& op, const Function<R,NDIM>& f, int axis, bool fence=true) {
        PROFILE_FUNC;
        Function<R,NDIM>& ff = const_cast< Function<R,NDIM>& >(f);
        if (VERIFY_TREE) ff.verify_tree();
        ff.reconstruct();

        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;

        result.set_impl(ff, false);
        result.get_impl()->apply_1d_realspace_push(op, ff.get_impl().get(), axis, fence);
        result.get_impl()->set_tree_state(redundant_after_merge);
        return result;
    }


    /// Generate a new function by reordering dimensions ... optional fence

    /// You provide an array of dimension NDIM that maps old to new dimensions
    /// according to
    /// \code
    ///    newdim = mapdim[olddim]
    /// \endcode
    /// Works in either scaling function or wavelet basis.
    ///
    /// Would be easy to modify this to also change the procmap here
    /// if desired but presently it uses the same procmap as f.
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        return result.mapdim(f,map,fence);
    }

    /// Generate a new function by mirroring within the dimensions .. optional fence

    /// similar to mapdim
    /// @param[in]	mirror	array with -1 and 1, corresponding to mirror this dimension or not
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    mirror(const Function<T,NDIM>& f, const std::vector<long>& mirrormap, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        return result.mirror(f,mirrormap,fence);
    }

    /// This is replaced with mirror(map(f)), optional fence

    /// first map then mirror!
    /// mirror is similar to mapdim, but maps from x to -x, y to -y, and so on
    /// Example: mirror a 3d function on the xy plane: mirror={1,1,-1}
    /// Example: c4 rotation of a 3d function around the z axis:
    /// 	x->y, y->-x, z->z: map(1,0,2); mirror(-1,1,1)
    /// @param[in]	map		array holding dimensions
    /// @param[in]	mirror	array of -1 and 1, corresponding to mirror or not
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    map_and_mirror(const Function<T,NDIM>& f, const std::vector<long>& map,
    		const std::vector<long>& mirror, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        return result.map_and_mirror(f,map,mirror,fence);
    }


    /// swap particles 1 and 2

    /// param[in]	f	a function of 2 particles f(1,2)
    /// return	the input function with particles swapped g(1,2) = f(2,1)
    template <typename T, std::size_t NDIM>
    typename std::enable_if_t<NDIM%2==0, Function<T,NDIM>>
    swap_particles(const Function<T,NDIM> & f){
      // this could be done more efficiently for SVD, but it works decently
      std::vector<long> map(NDIM);
      constexpr std::size_t LDIM=NDIM/2;
      static_assert(LDIM*2==NDIM);
      for (auto d=0; d<LDIM; ++d) {
          map[d]=d+LDIM;
          map[d+LDIM]=d;
      }
//      map[0]=3;
//      map[1]=4;
//      map[2]=5;     // 2 -> 1
//      map[3]=0;
//      map[4]=1;
//      map[5]=2;     // 1 -> 2
      return mapdim(f,map);
    }

    /// symmetrize a function

    /// @param[in]  symmetry possibilities are:
    ///                 (anti-) symmetric particle permutation ("sy_particle", "antisy_particle")
    ///                 symmetric mirror plane ("xy", "xz", "yz")
    /// @return     a new function symmetrized according to the input parameter
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    symmetrize(const Function<T,NDIM>& f, const std::string symmetry, bool fence=true) {
        Function<T,NDIM> result;

        MADNESS_ASSERT(NDIM==6);            // works only for pair functions
        std::vector<long> map(NDIM);

        // symmetric particle permutation
        if (symmetry=="sy_particle") {
            map[0]=3; map[1]=4; map[2]=5;
            map[3]=0; map[4]=1; map[5]=2;
        } else if (symmetry=="cx") {
            map[0]=0; map[1]=2; map[2]=1;
            map[3]=3; map[4]=5; map[5]=4;

        } else if (symmetry=="cy") {
            map[0]=2; map[1]=1; map[2]=0;
            map[3]=5; map[4]=4; map[5]=3;

        } else if (symmetry=="cz") {
            map[0]=1; map[1]=0; map[2]=2;
            map[3]=4; map[4]=3; map[5]=5;

        } else {
            if (f.world().rank()==0) {
                print("unknown parameter in symmetrize:",symmetry);
            }
            MADNESS_EXCEPTION("unknown parameter in symmetrize",1);
        }

        result.mapdim(f,map,true);  // need to fence here
        result.get_impl()->average(*f.get_impl());

        return result;
    }



    /// multiply a high-dimensional function with a low-dimensional function

    /// @param[in]  f   NDIM function of 2 particles: f=f(1,2)
    /// @param[in]  g   LDIM function of 1 particle: g=g(1) or g=g(2)
    /// @param[in]  particle    if g=g(1) or g=g(2)
    /// @return     h(1,2) = f(1,2) * g(p)
    template<typename T, std::size_t NDIM, std::size_t LDIM>
    Function<T,NDIM> multiply(const Function<T,NDIM> f, const Function<T,LDIM> g, const int particle, const bool fence=true) {

        static_assert(LDIM+LDIM==NDIM);
        MADNESS_ASSERT(particle==1 or particle==2);

        Function<T,NDIM> result;
        result.set_impl(f, false);

//        Function<T,NDIM>& ff = const_cast< Function<T,NDIM>& >(f);
//        Function<T,LDIM>& gg = const_cast< Function<T,LDIM>& >(g);

        f.change_tree_state(redundant,false);
        g.change_tree_state(redundant);
		FunctionImpl<T,NDIM>* fimpl=f.get_impl().get();
		FunctionImpl<T,LDIM>* gimpl=g.get_impl().get();

        result.get_impl()->multiply(fimpl,gimpl,particle);
        result.world().gop.fence();

        f.change_tree_state(reconstructed,false);
        g.change_tree_state(reconstructed);
        return result;
    }


    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    project(const Function<T,NDIM>& other,
            int k=FunctionDefaults<NDIM>::get_k(),
            double thresh=FunctionDefaults<NDIM>::get_thresh(),
            bool fence=true)
    {
        PROFILE_FUNC;
        Function<T,NDIM> result = FunctionFactory<T,NDIM>(other.world()).k(k).thresh(thresh).empty();
        other.reconstruct();
        result.get_impl()->project(*other.get_impl(),fence);
        return result;
    }


    /// Computes the scalar/inner product between two functions

    /// In Maple this would be \c int(conjugate(f(x))*g(x),x=-infinity..infinity)
    template <typename T, typename R, std::size_t NDIM>
    TENSOR_RESULT_TYPE(T,R) inner(const Function<T,NDIM>& f, const Function<R,NDIM>& g) {
        PROFILE_FUNC;
        return f.inner(g);
    }


    /// Computes the partial scalar/inner product between two functions, returns a low-dim function

    /// syntax similar to the inner product in tensor.h
    /// e.g result=inner<3>(f,g),{0},{1}) : r(x,y) = int f(x1,x) g(y,x1) dx1
    /// @param[in]  task    0: everything, 1; prepare only (fence), 2: work only (no fence), 3: finalize only (fence)
    template<std::size_t NDIM, typename T, std::size_t LDIM, typename R, std::size_t KDIM,
            std::size_t CDIM = (KDIM + LDIM - NDIM) / 2>
    std::vector<Function<TENSOR_RESULT_TYPE(T, R), NDIM>>
    innerXX(const Function<T, LDIM>& f, const std::vector<Function<R, KDIM>>& vg, const std::array<int, CDIM> v1,
           const std::array<int, CDIM> v2, int task=0) {
        bool prepare = ((task==0) or (task==1));
        bool work = ((task==0) or (task==2));
        bool finish = ((task==0) or (task==3));

        static_assert((KDIM + LDIM - NDIM) % 2 == 0, "faulty dimensions in inner (partial version)");
        static_assert(KDIM + LDIM - 2 * CDIM == NDIM, "faulty dimensions in inner (partial version)");

        // contraction indices must be contiguous and either in the beginning or at the end
        for (int i=0; i<CDIM-1; ++i) MADNESS_CHECK((v1[i]+1)==v1[i+1]);
        MADNESS_CHECK((v1[0]==0) or (v1[CDIM-1]==LDIM-1));

        for (int i=0; i<CDIM-1; ++i) MADNESS_CHECK((v2[i]+1)==v2[i+1]);
        MADNESS_CHECK((v2[0]==0) or (v2[CDIM-1]==KDIM-1));

        MADNESS_CHECK(f.is_initialized());
        MADNESS_CHECK(vg[0].is_initialized());
        MADNESS_CHECK(f.world().id() == vg[0].world().id());
        // this needs to be run in a single world, so that all coefficients are local.
        // Use macrotasks if run on multiple processes.
        World& world=f.world();
        MADNESS_CHECK(world.size() == 1);

        if (prepare) {
            f.change_tree_state(nonstandard);
            change_tree_state(vg,nonstandard);
            world.gop.fence();
            f.get_impl()->compute_snorm_and_dnorm(false);
            for (auto& g : vg) g.get_impl()->compute_snorm_and_dnorm(false);
            world.gop.fence();
        }

        typedef TENSOR_RESULT_TYPE(T, R) resultT;
        std::vector<Function<resultT,NDIM>> result(vg.size());
        if (work) {
            world.gop.set_forbid_fence(true);
            for (int i=0; i<vg.size(); ++i)  {
                result[i]=FunctionFactory<resultT,NDIM>(world)
                        .k(f.k()).thresh(f.thresh()).empty().nofence();
                result[i].get_impl()->partial_inner(*f.get_impl(),*(vg[i]).get_impl(),v1,v2);
                result[i].get_impl()->set_tree_state(nonstandard_after_apply);
            }
            world.gop.set_forbid_fence(false);
        }

        if (finish) {

            world.gop.fence();
//            result.get_impl()->reconstruct(true);

            change_tree_state(result,reconstructed);
//            result.reconstruct();
            // restore initial state of g and h
            auto erase_list = [] (const auto& funcimpl) {
                typedef typename std::decay_t<decltype(funcimpl)>::keyT keyTT;
                std::list<keyTT> to_be_erased;
                for (auto it=funcimpl.get_coeffs().begin(); it!=funcimpl.get_coeffs().end(); ++it) {
                    const auto& key=it->first;
                    const auto& node=it->second;
                    if (not node.has_children()) to_be_erased.push_back(key);
                }
                return to_be_erased;
            };

            FunctionImpl<T,LDIM>& f_nc=const_cast<FunctionImpl<T,LDIM>&>(*f.get_impl());
            for (auto& key : erase_list(f_nc)) f_nc.get_coeffs().erase(key);
            for (auto& g : vg) {
                FunctionImpl<R,KDIM>& g_nc=const_cast<FunctionImpl<R,KDIM>&>(*g.get_impl());
                for (auto& key : erase_list(g_nc)) g_nc.get_coeffs().erase(key);
            }
            world.gop.fence();
            change_tree_state(vg,reconstructed);
            f_nc.reconstruct(false);
            world.gop.fence();

        }

        return result;
    }


    /// Computes the partial scalar/inner product between two functions, returns a low-dim function

    /// syntax similar to the inner product in tensor.h
    /// e.g result=inner<3>(f,g),{0},{1}) : r(x,y) = int f(x1,x) g(y,x1) dx1
    /// @param[in]  task    0: everything, 1; prepare only (fence), 2: work only (no fence), 3: finalize only (fence)
    template<std::size_t NDIM, typename T, std::size_t LDIM, typename R, std::size_t KDIM,
            std::size_t CDIM = (KDIM + LDIM - NDIM) / 2>
    Function<TENSOR_RESULT_TYPE(T, R), NDIM>
    innerXX(const Function<T, LDIM>& f, const Function<R, KDIM>& g, const std::array<int, CDIM> v1,
            const std::array<int, CDIM> v2, int task=0) {
        return innerXX<NDIM,T,LDIM,R,KDIM>(f,std::vector<Function<R,KDIM>>({g}),v1,v2,task)[0];
    }

    /// Computes the partial scalar/inner product between two functions, returns a low-dim function

    /// syntax similar to the inner product in tensor.h
    /// e.g result=inner<3>(f,g),{0},{1}) : r(x,y) = int f(x1,x) g(y,x1) dx1
    template <typename T, std::size_t LDIM, typename R, std::size_t KDIM>
    Function<TENSOR_RESULT_TYPE(T,R),KDIM+LDIM-2>
    inner(const Function<T,LDIM>& f, const Function<R,KDIM>& g, const std::tuple<int> v1, const std::tuple<int> v2) {
        return innerXX<KDIM+LDIM-2>(f,g,
                     std::array<int,1>({std::get<0>(v1)}),
                     std::array<int,1>({std::get<0>(v2)}));
    }

    /// Computes the partial scalar/inner product between two functions, returns a low-dim function

    /// syntax similar to the inner product in tensor.h
    /// e.g result=inner<3>(f,g),{0,1},{1,2}) : r(y) = int f(x1,x2) g(y,x1,x2) dx1 dx2
    template <typename T, std::size_t LDIM, typename R, std::size_t KDIM>
    Function<TENSOR_RESULT_TYPE(T,R),KDIM+LDIM-4>
    inner(const Function<T,LDIM>& f, const Function<R,KDIM>& g, const std::tuple<int,int> v1, const std::tuple<int,int> v2) {
        return innerXX<KDIM+LDIM-4>(f,g,
                                  std::array<int,2>({std::get<0>(v1),std::get<1>(v1)}),
                                  std::array<int,2>({std::get<0>(v2),std::get<1>(v2)}));
    }

    /// Computes the partial scalar/inner product between two functions, returns a low-dim function

    /// syntax similar to the inner product in tensor.h
    /// e.g result=inner<3>(f,g),{1},{2}) : r(x,y,z) = int f(x,x1) g(y,z,x1) dx1
    template <typename T, std::size_t LDIM, typename R, std::size_t KDIM>
    Function<TENSOR_RESULT_TYPE(T,R),KDIM+LDIM-6>
    inner(const Function<T,LDIM>& f, const Function<R,KDIM>& g, const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
        return innerXX<KDIM+LDIM-6>(f,g,
                                  std::array<int,3>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)}),
                                  std::array<int,3>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)}));
    }



    /// Computes the scalar/inner product between an MRA function and an external functor

    /// Currently this defaults to inner_adaptive, which might be more expensive
    /// than inner_ext since it loops over all leaf nodes. If you feel inner_ext
    /// is more efficient you need to call it directly
    /// @param[in]  f   MRA function
    /// @param[in]  g   functor
    /// @result     inner(f,g)
    template <typename T, typename opT, std::size_t NDIM>
    TENSOR_RESULT_TYPE(T,typename opT::value_type) inner(const Function<T,NDIM>& f, const opT& g) {
        PROFILE_FUNC;
        std::shared_ptr< FunctionFunctorInterface<double,3> > func(new opT(g));
        return f.inner_adaptive(func);
    }

    /// Computes the scalar/inner product between an MRA function and an external functor

    /// Currently this defaults to inner_adaptive, which might be more expensive
    /// than inner_ext since it loops over all leaf nodes. If you feel inner_ext
    /// is more efficient you need to call it directly
    /// @param[in]  g   functor
    /// @param[in]  f   MRA function
    /// @result     inner(f,g)
    template <typename T, typename opT, std::size_t NDIM>
    TENSOR_RESULT_TYPE(T,typename opT::value_type) inner(const opT& g, const Function<T,NDIM>& f) {
        return inner(f,g);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator+(const Function<T,NDIM>& f, R r) {
        return (f*R(1.0)).add_scalar(r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator+(R r, const Function<T,NDIM>& f) {
        return (f*R(1.0)).add_scalar(r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator-(const Function<T,NDIM>& f, R r) {
        return (f*R(1.0)).add_scalar(-r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator-(R r, const Function<T,NDIM>& f) {
        return (f*R(-1.0)).add_scalar(r);
    }

    namespace detail {
        template <std::size_t NDIM>
        struct realop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                return real(t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        template <std::size_t NDIM>
        struct imagop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                return imag(t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        template <std::size_t NDIM>
        struct abssqop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                Tensor<double> r = abs(t);
                return r.emul(r);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        template <std::size_t NDIM>
        struct absop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                Tensor<double> r = abs(t);
                return r;
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

    }

    /// Returns a new function that is the real part of the input
    template <std::size_t NDIM>
    Function<double,NDIM> real(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op_coeffs(z, detail::realop<NDIM>(), fence);
    }

    /// Returns a new function that is the real part of the input
    template <std::size_t NDIM>
    Function<double,NDIM> real(const Function<double,NDIM>& z, bool fence=true) {
    	return copy(z);
    }

    /// Returns a new function that is the imaginary part of the input
    template <std::size_t NDIM>
    Function<double,NDIM> imag(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op_coeffs(z, detail::imagop<NDIM>(), fence);
    }


    /// Create a new function that is the square of f - global comm only if not reconstructed
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> square(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.square(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /// Create a new function that is the abs of f - global comm only if not reconstructed
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> abs(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.abs(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /// Create a new function that is the abs_square of f - global comm only if not reconstructed
    template <typename T, std::size_t NDIM>
    typename std::enable_if<!TensorTypeData<T>::iscomplex, Function<T,NDIM> >::type
    abs_square(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.abs_square(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /// Create a new function that is the abs_square of f - global comm only if not reconstructed
    template <typename T, std::size_t NDIM>
    typename std::enable_if<TensorTypeData<T>::iscomplex, Function<typename Tensor<T>::scalar_type,NDIM> >::type
   	abs_square(const Function<T,NDIM>& f, bool fence=true) {
        return unary_op(f, detail::abssqop<NDIM>(), fence);
    }

    /// Returns a new function that is the square of the absolute value of the input
    template <std::size_t NDIM>
    Function<double,NDIM> abssq(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op(z, detail::abssqop<NDIM>(), fence);
    }

    /// Returns a new function that is the absolute value of the input
    template <std::size_t NDIM>
    Function<double,NDIM> abs(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op(z, detail::absop<NDIM>(), fence);
    }

    /// get tree state of a function

    /// there is a corresponding function in vmra.h
    /// @param[in]  f   function
    /// @return TreeState::unknown if the function is not initialized
    template <typename T, std::size_t NDIM>
    TreeState get_tree_state(const Function<T,NDIM>& f) {
        if (f.is_initialized()) return f.get_impl()->get_tree_state();
        return TreeState::unknown;
    }

    /// change tree state of a function

    /// there is a corresponding function in vmra.h
    /// return this for chaining
    /// @param[in]  f   function
    /// @param[in]  finalstate  the new state
    /// @return this in the requested state
    template <typename T, std::size_t NDIM>
    const Function<T,NDIM>& change_tree_state(const Function<T,NDIM>& f,
            const TreeState finalstate, bool fence=true) {
        return f.change_tree_state(finalstate,fence);
    }


}

#include <madness/mra/funcplot.h>

namespace madness {
    namespace archive {
        template <class archiveT, class T, std::size_t NDIM>
        struct ArchiveLoadImpl< ParallelInputArchive<archiveT>, Function<T,NDIM> > {
            static inline void load(const ParallelInputArchive<archiveT>& ar, Function<T,NDIM>& f) {
                f.load(*ar.get_world(), ar);
            }
        };

        template <class archiveT, class T, std::size_t NDIM>
        struct ArchiveStoreImpl< ParallelOutputArchive<archiveT>, Function<T,NDIM> > {
            static inline void store(const ParallelOutputArchive<archiveT>& ar, const Function<T,NDIM>& f) {
                f.store(ar);
            }
        };
    }

    template <class T, std::size_t NDIM>
    void save(const Function<T,NDIM>& f, const std::string name) {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar2(f.world(), name.c_str(), 1);
        ar2 & f;
    }

    template <class T, std::size_t NDIM>
    void load(Function<T,NDIM>& f, const std::string name) {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar2(f.world(), name.c_str(), 1);
        ar2 & f;
    }

}

namespace madness {
    // type traits to check if a template parameter is a Function
    template<typename>
    struct is_madness_function : std::false_type {};

    template<typename T, std::size_t NDIM>
    struct is_madness_function<madness::Function<T, NDIM>> : std::true_type {};

    template<typename>
    struct is_madness_function_vector : std::false_type {
    };

    template<typename T, std::size_t NDIM>
    struct is_madness_function_vector<std::vector<typename madness::Function<T, NDIM>>> : std::true_type {
};

}


/* @} */

#include <madness/mra/derivative.h>
#include <madness/mra/operator.h>
#include <madness/mra/functypedefs.h>
#include <madness/mra/vmra.h>
// #include <madness/mra/mraimpl.h> !!!!!!!!!!!!! NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO  !!!!!!!!!!!!!!!!!!

#endif // MADNESS_MRA_MRA_H__INCLUDED
