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

    template <typename T, int NDIM, typename Pmap> class LoadBalImpl;

    template <typename T, int NDIM, typename Pmap=DCDefaultProcmap<Key<NDIM> > >
    class Function {
	friend class LoadBalImpl<T,NDIM,Pmap>;
    private:
        SharedPtr< FunctionImpl<T,NDIM,Pmap> > impl;
    public:
        typedef FunctionImpl<T,NDIM,Pmap> implT;
        typedef FunctionFactory<T,NDIM,Pmap> factoryT;
        typedef typename implT::coordT coordT; ///< Type of vector holding coordinates 

        /// Default constructor makes uninitialized function.  No communication.

        /// An unitialized function can only be assigned to.  Any other operation will throw.
        Function()
            : impl(0)
        {};


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const factoryT& factory)
            : impl(new FunctionImpl<T,NDIM,Pmap>(factory))
        {};


        /// Copy constructor is \em shallow.  No communication, works in either basis.
        Function(const Function<T,NDIM,Pmap>& f)
            : impl(f.impl)
        {};


        /// Assignment is \em shallow.  No communication, works in either basis.
        Function<T,NDIM,Pmap>& operator=(const Function<T,NDIM,Pmap>& f) {
            if (this != &f) impl = f.impl;
            return *this;
        };

        ~Function(){};

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        Future<T> eval(const coordT& xuser) {
            coordT xsim;
            impl->user_to_sim(xuser,xsim);
            Future<T> result;
            impl->eval(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        };


    private:

    };
}
#endif
