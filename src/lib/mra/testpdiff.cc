#include <mra/mra.h>
#include <constants.h>

using namespace std;
using namespace madness;

static const int NDIM=3;
static const double L=2.0; // Simulate in [0,L]^NDIM

typedef Vector<double,NDIM> coordT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<double,NDIM> factoryT;

double func(const coordT& r) {
    // Test is product of cos(2*Pi*(d+1)/L + 0.2)
    double result = 1.0;
    for (int d=0; d<NDIM; d++) {
        double fac = 2.0*constants::pi/L;
        result *= cos(fac*r[d] + 0.2);
    }
    return result;
}

int axis; // Axis we are diffing

double dfunc(const coordT& r) {
    // Test is product of cos(2*Pi*(d+1)/L + 0.2)
    double result = 1.0;
    for (int d=0; d<NDIM; d++) {
        double fac = 2.0*constants::pi/L;
        if (d == axis) {
            result *= -fac*sin(fac*r[d] + 0.2);
        }
        else {
            result *= cos(fac*r[d] + 0.2);
        }
    }
    return result;
}

struct FunctorInterfaceWrapper : public FunctionFunctorInterface<double,NDIM> {
    double (*f)(const coordT&);
    
    FunctorInterfaceWrapper(double (*f)(const coordT&)) : f(f) {}
    
    double operator()(const coordT& x) const {return f(x);}
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    FunctionDefaults<NDIM>::set_cubic_cell(0.0,L);

    Tensor<int> bc_periodic(NDIM,2);
    bc_periodic.fill(1);

    Tensor<int> bc_zero(NDIM,2);
    bc_zero.fill(0);

    functionT f = factoryT(world).f(func);
    f.set_bc(bc_periodic);

    for (axis=0; axis<NDIM; axis++) {
        functionT df = diff(f,axis);
        print(axis,"error",df.err(FunctorInterfaceWrapper(dfunc)));
    }

    finalize();
    return 0;
}
