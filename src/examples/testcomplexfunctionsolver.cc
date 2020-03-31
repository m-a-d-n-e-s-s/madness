#include <iostream>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>
#include <complex>
#include <cmath>

using namespace madness;

const size_t N = 1;
using dataT = std::complex<double>;
using functionT = Function<dataT,N>;
using factoryT = FunctionFactory<dataT,N>;

class F;
F allocator();

class F {
    functionT x;

public:
    F(const functionT& x) : x(copy(x)) {}

    F() : x(allocator().x) {}      // Default constructor necessary for storage in vector

    F(const F& a) : x(copy(a.x)) {} // Copy constructor necessary

    F operator=(const F& f) { // Assignment required for storage in vector
        if (this != &f) {
            x = copy(f.x);
        }
        return *this;
    }

    F operator-(const F& b) const { // Operator- necessary
        return F(x-b.x);
    }

    F& operator+=(const F& b) { // Operator+= necessary
        x += b.x;
        return *this;
    }

    F operator*(const dataT& a) const { // Scale by a constant necessary
        return F(x*a);
    }

    const functionT& get() const {return x;}
};

// If the default constructor does not make a zero value need an
// allocator.  It can be a function or a class.
F allocator() {
    return F(functionT(factoryT(World::get_default())));
}

// This interface is necessary to compute inner products
dataT inner(const F& a, const F& b) {
    return conj(a.get()).inner(b.get());
}

// void residualx(const Key<N>& key, Tensor<dataT>& t) {
//     UNARY_OPTIMIZED_ITERATOR(dataT, t, *_p0 = std::exp(-*_p0) - *_p0;);
// }

void residualx(const Key<N>& key, Tensor<dataT>& t) {
    UNARY_OPTIMIZED_ITERATOR(dataT, t, *_p0 = (*_p0)*(*_p0)+1.0;);
}

F residual(const F& f) {
    functionT r = copy(f.get());
    r.unaryop(&residualx);
    return F(r);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world,argc,argv);
        XNonlinearSolver<F,dataT,F(*)()> solver(allocator);
        solver.set_maxsub(5);
        
        // This line should compile but won't work because the
        // default constructor F() sets x=99 not zero
        //XNonlinearSolver<F,double> solver;
        
        functionT guess = factoryT(world);
        guess.add_scalar(std::complex<double>(0.5,0.5));
        F x(guess);
        for (int iter=0; iter<20; iter++) {
            auto r = residual(x);
            std::cout << iter << " " << residual(x).get().norm2() << " " << x.get()(coord_1d({0.5})) << std::endl;
            x = solver.update(x, residual(x));
        }
        world.gop.fence();
        world.gop.fence();
    }
    finalize();

    return 0;
}

