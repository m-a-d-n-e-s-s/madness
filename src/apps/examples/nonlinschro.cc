/// \file nonlinschro.cc
/// \brief Solves 1D nonlinear Schrodinger equation


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <algorithm>

using namespace madness;

typedef Vector<double,1> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,1> > functorT;
typedef Function<double,1> functionT;
typedef FunctionFactory<double,1> factoryT;
typedef SeparatedConvolution<double,1> operatorT;

static const double L = 10.0;           // box size
static const long k = 10;                // wavelet order
static const double thresh = 1e-8;      // precision
static const double ntol = thresh*10.0; // Cutoff for numerically small density

// V := -a*exp(-b*x^2)+c/(n(x)+eta)^(1/3)-d*n(x)^(5/3) + vshift
// asymptotic value is c/eta^1/3 + vshift

static const double a = 0.125;
static const double b = 1.0;
static const double c = 1.0;
static const double d = 1.0;
static const double eta_end = 1e-5; // Target value for eta
static double eta = eta_end*1024;;  // Initial value foe eta

static const double vasymp = c/std::pow(eta_end,1.0/3.0);


static const double vshift = -3.0; // Only impacts BSH operator

// For x<xmax, smoothly restricts x to be greater than or equal to xmin>0
double munge(double x, double xmin, double xmax) {
    // A quartic polyn that smoothly interpolates between
    //
    // x=0    where it has value xmin and zero slope, and
    //
    // x=xmax where it has value xmax, unit slope, and zero second derivative
    
    if (x > xmax) {
        return x;
    }
    else if (x <= 0.0) {
        return xmin;
    }
    else {
        double x1 = x/xmax;
        double x2 = x1*x1;
        double x3 = x1*x2;
        double x4 = x2*x2;
        return xmin+(3.0*xmax-6.0*xmin)*x2 - (3*xmax-8*xmin)*x3 + (-3*xmin+xmax)*x4;
    }    
}

/// This invoked to compute the part of the potential that is a function of the density
template <typename T>
static void Vdynamic(const Key<1> & key, Tensor<T>& t)
{
    static const double onethird = 1.0/3.0;
    static const double fivethird = 5.0/3.0;
    UNARY_OPTIMIZED_ITERATOR(T, t, 
                             double n=munge(*_p0, ntol, 10.0*ntol); 
                             *_p0 = c/pow(n+eta,onethird) - d*pow(n+eta,fivethird));
}

/// The part of the potential that does not depend on the density
static double Vstatic(const coordT& r) {
    const double x=r[0];
    return -a*exp(-b*x*x) + vshift;
}

static double guess(const coordT& r) {
    const double x=r[0];
    return exp(-0.2*x*x);
}

functionT make_potential(World& world, const functionT& rho) {
    functionT vstatic = factoryT(world).f(Vstatic);
    functionT vdynamic = copy(rho);
    vdynamic.unaryop(&Vdynamic<double>);
    return vstatic + vdynamic;
}

void iterate(World& world, functionT& psi) {
    // Compute density and full potential
    functionT rho = psi*psi;
    functionT v = make_potential(world, rho);

    // Compute energy components and print
    functionT dpsi = diff(psi,0);
    double kinetic_energy = inner(dpsi,dpsi);
    double potential_energy = inner(rho, v);
    double energy = potential_energy + kinetic_energy;

    // Update the wave function
    functionT Vpsi = v*psi;
    Vpsi.scale(-1.0).truncate();
    operatorT op = BSHOperator<double,1>(world, sqrt(-energy), k, 0.001, 1e-6);
    functionT tmp = apply(op,Vpsi).truncate();
    functionT r = tmp-psi;
    double rnorm = r.norm2();
    //double step = std::min(1.0,0.05/rnorm);
    double step = 0.025;
    if (rnorm > 10) step = 0.1/rnorm;
    psi = psi + step*r;
    psi.scale(1.0/psi.norm2());

    print("KE =", kinetic_energy,"  PE =",potential_energy,"  E =",energy, "  err(psi) =", rnorm, "  step =", step);
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_truncate_mode(1);  
    FunctionDefaults<1>::set_cubic_cell(-L,L);

    functionT psi = factoryT(world).f(guess);
    print(psi.norm2());
    psi.scale(1.0/psi.norm2());
    print(psi.norm2());

    while (1) {
        print("\nSOLVING WITH ETA", eta);
        for (int iter=0; iter<150; iter++)
            iterate(world, psi);
        if (eta <= eta_end) break;
        eta *= 0.5;
    }

    functionT rho = psi*psi;
    functionT v = make_potential(world, rho);

    coordT lo(-L), hi(L);
    double scale = vasymp/psi(coordT(0.0));
    plot_line("psi.txt", 201, lo, hi, psi*scale, v);

    world.gop.fence();

    finalize();
    return 0;
}
