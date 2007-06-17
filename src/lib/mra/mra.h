#ifndef MAD_MRA_H
#define MAD_MRA_H

#include <world/world.h>
#include <misc/misc.h>
#include <tensor/mtrand.h>
#include <tensor/tensor.h>

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

namespace madness {

    //template <typename T, int NDIM, typename Pmap> class LoadBalImpl;

    template <typename T, int NDIM>
    class Function {
	//friend class LoadBalImpl<T,NDIM,Pmap>;
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

        ~Function(){};

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        Future<T> eval(const coordT& xuser) {
            verify();
            coordT xsim;
            impl->user_to_sim(xuser,xsim);
            Future<T> result;
            impl->eval(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        };

        /// Returns true if compressed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_compressed() const {
            if (impl) 
                return impl->is_compressed();
            else
                return false;
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
        
//         /// Reconstructs the function, transforming into scaling function basis.  Possible non-blocking comm.

//         /// By default fence=true meaning that this operation completes before returning,
//         /// otherwise if fence=false it returns without fencing and the user must invoke
//         /// world.gop.fence() to assure global completion before using the function
//         /// for other purposes.
//         ///
//         /// Noop if already reconstructed or if not initialized.
//         void compress(bool fence = true) {
//             if (!impl || !is_compressed()) return;
//             impl->compress(fence);
//         };


    private:

    };
}
#endif
