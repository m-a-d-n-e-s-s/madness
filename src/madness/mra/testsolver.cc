#include <iostream>
#include <madness/mra/nonlinsol.h>
#include <cmath>
#include<madness/world/test_utilities.h>

using namespace madness;

/// simple class for testing the solver
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

F exact_solution() {
    return F(0.5671432904097839);
}


// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
template<typename T, std::size_t NDIM>
struct function_allocator {
    World& world;
    const int n=-1;

    /// @param[in]	world	the world
    function_allocator(World& world) : world(world) {}

    /// allocate a vector of n empty functions
    Function<T, NDIM> operator()() {
        auto f=FunctionFactory<T,NDIM>(world);
        return f;
    }
};

int test_simple(World& world) {
    test_output t1("testing test_simple with class F");
    {
        double cpu0=cpu_time();
        XNonlinearSolver<F,double,F(*)()> solver(allocator);
        F x = 0.5;
        for (int iter=0; iter<8; iter++) {
            std::cout << iter << " " << x.get() << std::endl;
            x = solver.update(x, residual(x));
        }
        double cpu1=cpu_time();
        double error=fabs(x.get()-exact_solution().get());
        t1.checkpoint(error,1.e-6,"error after 8 iterations using residual",cpu1-cpu0);
    }
    {
        double cpu0=cpu_time();
        XNonlinearSolver<F,double,F(*)()> solver(allocator);
        F x = 0.5;
        solver.initialize(F(0.5));
        for (int iter=0; iter<8; iter++) {
            std::cout << iter << " " << x.get() << std::endl;
            auto xpreliminary=x-residual(x);    // looks stupid, but it is what it is here
            x = solver.update(xpreliminary);
        }
        double cpu1=cpu_time();
        double error=fabs(x.get()-exact_solution().get());
        t1.checkpoint(error,1.e-6,"error after 8 iterations using update",cpu1-cpu0);
    }


    return t1.end();
}


template<typename T, std::size_t NDIM>
Function<T,NDIM> compute_update(World& world, const Function<T,NDIM>& f) {
    auto op=BSHOperator<NDIM>(world,-1.0,1.e-4,FunctionDefaults<NDIM>::get_thresh());
    auto pot_functor=[](const Vector<double,NDIM>& r) {
        return inner(r,r)-5.0;
    };
    Function<T,NDIM> V=FunctionFactory<T,NDIM>(world).f(pot_functor);
    return op(-2.0*V*f);
}

template<typename T, std::size_t NDIM>
int test_with_function(World& world) {
    test_output t1("test with function");
    Function<T,NDIM> guess=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r){return exp(-inner(r,r));});
    double n=guess.norm2();
    guess.scale(1.0/n);
    {
        double cpu0=cpu_time();
        auto f=copy(guess);
        double error1=1.0;
        for (int i=0; i<10; ++i) {
            auto fnew=compute_update(world,f);
            double n2=fnew.norm2();
            fnew.scale(1.0/n2);
            error1=(f-fnew).norm2();
            print("error in iteration",i,error1);
            f=fnew;
        }
        double cpu1=cpu_time();
        t1.checkpoint(error1,1.e-6,"without kain ",cpu1-cpu0);
    }
    {
        double cpu0=cpu_time();
        auto alloc=function_allocator<T,NDIM>(world);
        XNonlinearSolver<Function<T,NDIM>,double,function_allocator<T,NDIM>> solver(alloc);

        auto f=copy(guess);
        double error1=1.0;
        for (int i=0; i<10; ++i) {
            auto up=compute_update(world,f);
            double n2=up.norm2();
            up.scale(1.0/n2);
            auto residual = f-up;
            auto fnew=solver.update(f,residual);
            error1=(f-fnew).norm2();
            print("error in iteration",i,error1);
            f=fnew;
        }
        double cpu1=cpu_time();
        t1.checkpoint(error1,1.e-7,"with KAIN/residual",cpu1-cpu0);
    }
    {
        double cpu0=cpu_time();
        auto alloc=function_allocator<T,NDIM>(world);
        XNonlinearSolver<Function<T,NDIM>,double,function_allocator<T,NDIM>> solver(alloc);
        guess.reconstruct();
        solver.initialize(guess);

        auto f=copy(guess);
        double error1=1.0;
        for (int i=0; i<10; ++i) {
            auto fnew=compute_update(world,f);
            fnew.reconstruct();
            double n2=fnew.norm2();
            fnew.scale(1.0/n2);
            fnew=solver.update(fnew);
            error1=(f-fnew).norm2();
            print("error in iteration",i,error1);
            f=fnew;
        }
        double cpu1=cpu_time();
        t1.checkpoint(error1,1.e-7,"with KAIN/update",cpu1-cpu0);
    }


    return t1.end();
}


/// test solver with a vector of functions

/// keep it simple: both functions converge to the same state, no orthogonalization needed
template<typename T, std::size_t NDIM>
int test_with_function_vector(World& world) {
    test_output t1("test with function vector with dim "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    Function<T,NDIM> guess1=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r){return exp(-inner(r,r));});
    Function<T,NDIM> guess2=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r){return exp(-2.0*inner(r,r));});
    std::vector<Function<T,NDIM>> guess({guess1,guess2});
    normalize(world,guess);
    reconstruct(guess);
    {
        double cpu0=cpu_time();
        auto solver=nonlinear_vector_solver<T,NDIM>(world,2);
        // solver.do_print=true;

        auto f=copy(guess);
        double error1=1.0;
        for (int i=0; i<10; ++i) {
            auto f1=compute_update(world,f[0]);
            auto f2=compute_update(world,f[1]);
            auto fnew=std::vector<Function<T,NDIM>>({f1,f2});
            reconstruct(fnew);
            normalize(world,fnew);
            auto residual=f-fnew;
            error1=norm2(world,residual);
            print("error in iteration",i,error1);
            f=solver.update(f,residual);
        }
        double cpu1=cpu_time();
        t1.checkpoint(error1,5.e-6,"with KAIN/residual",cpu1-cpu0);
    }
    {
        double cpu0=cpu_time();
        auto solver=nonlinear_vector_solver<T,NDIM>(world,2);
        auto f=copy(guess);
        solver.initialize(guess);
        // solver.do_print=true;

        double error1=1.0;
        for (int i=0; i<10; ++i) {
            auto f1=compute_update(world,f[0]);
            auto f2=compute_update(world,f[1]);
            auto fnew=std::vector<Function<T,NDIM>>({f1,f2});
            reconstruct(fnew);
            normalize(world,fnew);
            error1=norm2(world,f-fnew);
            print("error in iteration",i,error1);
            f=solver.update(fnew);
        }
        double cpu1=cpu_time();
        t1.checkpoint(error1,5.e-6,"with KAIN/update",cpu1-cpu0);
    }


    return t1.end();
}

int main(int argc, char** argv) {

    // This line should compile but won't work because the
    // default constructor F() sets x=99 not zero
    //XNonlinearSolver<F,double> solver;

    World& world=initialize(argc,argv);
    world.gop.fence();
    startup(world,argc,argv);
    FunctionDefaults<1>::set_cubic_cell(-2.0,2.0);
    FunctionDefaults<1>::set_thresh(1.e-6);
    FunctionDefaults<1>::set_k(9);
    FunctionDefaults<2>::set_cubic_cell(-2.0,2.0);
    FunctionDefaults<2>::set_thresh(1.e-7);
    FunctionDefaults<2>::set_k(9);
    FunctionDefaults<2>::set_tensor_type(TT_2D);




    std::cout << std::setprecision(10);
    int ierr=0;

    ierr+=test_simple(world);
    ierr+=test_with_function<double,1>(world);
    ierr+=test_with_function_vector<double,1>(world);
    ierr+=test_with_function_vector<double,2>(world);



    madness::finalize();
    return ierr;
}

