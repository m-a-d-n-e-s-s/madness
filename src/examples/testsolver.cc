#include <iostream>
#include <madness/mra/nonlinsol.h>
#include <cmath>

using namespace madness;

class F {
    double x;

public:

    F(double x) : x(x) {}

    F() : x(99) {}      // Default constructor necessary for storage in vector

    F(const F& a) : x(a.x) {} // Copy constructor necessary

    F operator=(const F& f) { // Assignment required for storage in vector
        if (this != &f) {
            x = f.x;
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

    F operator*(double a) const { // Scale by a constant necessary
        return F(x*a);
    }

    double get() const {return x;}
};

// This interface is necessary to compute inner products
double inner(const F& a, const F& b) {
    return a.get()*b.get();
}

// If the default constructor does not make a zero value need an
// allocator.  It can be a function or a class.
F allocator() {
    return F(0.0);
}

// The test code solves r(x) = exp(-x) - x = 0
F residual(const F& f) {
    double x = f.get();
    return F(std::exp(-x)-x);
}

int main() {
    XNonlinearSolver<F,double,F(*)()> solver(allocator);

    // This line should compile but won't work because the
    // default constructor F() sets x=99 not zero
    //XNonlinearSolver<F,double> solver;

    std::cout << std::setprecision(10);
    F x = 0.5;
    for (int iter=0; iter<8; iter++) {
        std::cout << iter << " " << x.get() << std::endl;
        x = solver.update(x, residual(x));
    }

    // again without computing the residual
    x = 0.5;
    XNonlinearSolver<F,double,F(*)()> solver1(allocator);
    solver1.initialize(F(0.5));
    for (int iter=0; iter<8; iter++) {
        std::cout << iter << " " << x.get() << std::endl;
        auto xpreliminary=x-residual(x);
        x = solver1.update(xpreliminary);
    }

    return 0;
}

