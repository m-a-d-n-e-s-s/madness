#include <complex>
#include <memory>
#include <madness/mra/mra.h>
#include <madness/world/worldmutex.h>

using namespace madness;

static const size_t D = 2;
typedef Vector<double,D> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<std::complex<double>,D> > functorT;
typedef Function<std::complex<double>,D> cfunctionT;
typedef FunctionFactory<std::complex<double>,D> factoryT;
typedef SeparatedConvolution<std::complex<double>,D> operatorT;

static const double R = 1.4;    // bond length
static const double L = 32.0*R; // box size
static const long k = 3;        // wavelet order
static const double thresh = 1e-3; // precision

static std::complex<double> f(const coordT& r)
{
    return std::complex<double>(0.0,2.0);
}

template <typename T, size_t D>
class WriteCoeffImpl : public Mutex {
    const int k;
    std::ostream& out;
public:
    WriteCoeffImpl(const int k, std::ostream& out)
        : k(k)
        , out(out)
    {}

    void operator()(const Key<D>& key, const Tensor< T >& t) const {
        ScopedMutex obolus(*this);
        out << key << " " << t << std::endl;
    }
};

template <typename T, size_t D>
class WriteCoeff {
    std::shared_ptr<WriteCoeffImpl<T,D>> impl;

public:
    WriteCoeff(const int k, std::ostream& out)
        : impl(new WriteCoeffImpl<T,D>(k, out))
    {}

    void operator()(const Key<D>& key, const Tensor< T >& t) const {
        (*impl)(key, t);
    }
};



int main(int argc, char** argv)
{
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<D>::set_k(k);
    FunctionDefaults<D>::set_thresh(thresh);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_truncate_mode(0);
    FunctionDefaults<D>::set_cubic_cell(-L/2, L/2);

    cfunctionT fun = factoryT(world).f(f);
    fun.truncate();

    cfunctionT sqrt_of_fun = copy(fun);
    auto op = WriteCoeff<std::complex<double>,D>(k, std::cout);
    fun.unaryop(op);
    world.gop.fence();

    finalize();
    return 0;
}
