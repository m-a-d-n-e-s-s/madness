#include <complex>
#include <memory>
#include <madness/mra/mra.h>

using namespace madness;

static const size_t D = 2;
typedef Vector<double,D> coordT;
typedef Key<D> keyT;
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

template <typename T, std::size_t NDIM>
void write_function_coeffs(const typename FunctionImpl<T,NDIM>::dcT& coeffs, std::ostream& out, const Key<NDIM>& key) {
    auto it = coeffs.find(key).get();
    if (it == coeffs.end()) {
        for (int i=0; i<key.level(); ++i) out << "  ";
        out << key << "  missing --> " << coeffs.owner(key) << "\n";
    }
    else {
        const auto& node = it->second;
        for (int i=0; i<key.level(); ++i) out << "  ";
        out << key << "  " << node << " --> " << "\n";
        if (node.has_children()) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                write_function_coeffs<T,NDIM>(coeffs, out, kit.key());
            }
        }
    }
}

template <typename T, std::size_t NDIM>
void write_function(const Function<T,NDIM>& f, std::ostream& out) {
    if (f.get_impl()->world.rank() == 0) {
        out << NDIM << " " << f.k() << " " << FunctionDefaults<NDIM>::get_cell() << std::endl;
        
        write_function_coeffs<T,NDIM>(f.get_impl()->get_coeffs(), out, Key<NDIM>(0));
    }
    f.get_impl()->world.gop.fence();
}

void test(World& world) {
    cfunctionT fun = factoryT(world).f(f);
    fun.truncate();
    write_function(fun,std::cout);    
}

int main(int argc, char** argv)
{
    World& world = initialize(argc, argv);
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<D>::set_k(k);
    FunctionDefaults<D>::set_thresh(thresh);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_truncate_mode(0);
    FunctionDefaults<D>::set_cubic_cell(-L/2, L/2);

    test(world);

    world.gop.fence();
    finalize();
    return 0;
}
