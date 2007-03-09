#ifndef MAD_MRA_H
#define MAD_MRA_H

#include <world/world.h>
#include <misc/misc.h>
#include <tensor/mtrand.h>
#include <tensor/tensor.h>

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

    template <typename T, int NDIM>
    class Function {
    private:
        SharedPtr< FunctionImpl<T,NDIM> > impl;

    public:

        /// Default constructor makes uninitialized function.  No communication.

        /// An unitialized function can only be assigned to.  Any other operation
        /// will throw an exception.

        Function()
            : impl(0)
        {};


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const FunctionFactory<T,NDIM>& factory)
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

    private:

    };
}
#endif
