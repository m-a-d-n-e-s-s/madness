
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <complex>

using namespace madness;

typedef Vector<double,1> coordT;
typedef SharedPtr< FunctionFunctorInterface<std::complex<double>,1> > functorT;
typedef Function<std::complex<double>,1> cfunctionT;
typedef FunctionFactory<std::complex<double>,1> factoryT;
typedef SeparatedConvolution<std::complex<double>,1> operatorT;

static const double R = 1.4;    // bond length
static const double L = 32.0*R; // box size
static const long k = 5;        // wavelet order
static const double thresh = 1e-3; // precision
static const double thresh1 = thresh*0.1;

static std::complex<double> f(const coordT& r)
{
    const double x=r[0];
    return std::complex<double>(0.0,2.0);
}

template <typename T>
inline static void csqrt_op(const Key<1>& key, Tensor< T >& t) {
    UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = sqrt(*_p0));
}



int main(int argc, char** argv)
{
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_refine(true);
    FunctionDefaults<1>::set_initial_level(2);
    FunctionDefaults<1>::set_truncate_mode(0);
    FunctionDefaults<1>::set_cubic_cell(-L/2, L/2);

    cfunctionT boring = factoryT(world).f(f);
    boring.truncate();

    cfunctionT sqrt_of_boring = copy(boring);
    sqrt_of_boring.unaryop(&csqrt_op<double_complex>);
    double error = (square(sqrt_of_boring) - boring).norm2();
    print("error is ", error);
    world.gop.fence();

    finalize();
    return 0;
}
