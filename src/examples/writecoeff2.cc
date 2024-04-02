#include <complex>
#include <memory>
#include <madness/mra/mra.h>

using namespace madness;

static const size_t D = 2;
typedef Vector<double,D> coordT;
typedef Key<D> keyT;
typedef double dataT; // was std::complex<double>
typedef std::shared_ptr< FunctionFunctorInterface<dataT,D> > functorT;
typedef Function<dataT,D> cfunctionT;
typedef FunctionFactory<dataT,D> factoryT;
typedef SeparatedConvolution<dataT,D> operatorT;

static const double L = 4.0;
static const long k = 5;        // wavelet order
static const double thresh = 1e-3; // precision

static dataT f(const coordT& r)
{
    double R = r.normf();
    return std::exp(-R*R);
}

template <typename T, std::size_t NDIM>
void write_function_coeffs(const Function<T,NDIM>& f, std::ostream& out, const Key<NDIM>& key) {
    const auto& coeffs = f.get_impl()->get_coeffs();
    auto it = coeffs.find(key).get();
    if (it == coeffs.end()) {
        for (int i=0; i<key.level(); ++i) out << "  ";
        out << key << "  missing --> " << coeffs.owner(key) << "\n";
    }
    else {
        const auto& node = it->second;
        if (node.has_coeff()) {
            auto values = f.get_impl()->coeffs2values(key, node.coeff());
            for (int i=0; i<key.level(); ++i) out << "  ";
            out << key << std::endl;
            for (size_t i=0; i< values.size(); i++) out << values.ptr()[i] << " ";
            out << std::endl;
        }
        if (node.has_children()) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                write_function_coeffs<T,NDIM>(f, out, kit.key());
            }
        }
    }
}

template <typename T, std::size_t NDIM>
void write_function(const Function<T,NDIM>& f, std::ostream& out) {
    if (f.get_impl()->world.rank() == 0) {
        out << NDIM << " " << f.k() << " " << FunctionDefaults<NDIM>::get_cell() << std::endl;
        
        write_function_coeffs(f, out, Key<NDIM>(0));
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
