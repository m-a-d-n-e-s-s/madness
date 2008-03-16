/// \file tdse.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;

// Convenient but sleazy use of global variables to define simulation parameters
static const double L = 64.0;
static const long k = 6;        // wavelet order
static const double thresh = 1e-4; // precision

/// Regularized 1/r potential.

/// Invoke as \c u(r/c)/c where \c c is the radius of the
/// smoothed volume.
static double smoothed_potential(double r) {
    const double PI = 3.1415926535897932384;
    const double THREE_SQRTPI = 5.31736155271654808184;
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-8){
        pot = erf(r)/r + (exp(-r2) + 16.0*exp(-4.0*r2))/(THREE_SQRTPI);
    } else{
        pot = (2.0 + 17.0/3.0)/sqrt(PI);
    }
    
    return pot;
}

/// Smoothed 1/r nuclear potential
static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double rr = sqrt(x*x+y*y+z*z);
    const double cut = 0.1;
    return -smoothed_potential(rr/cut)/cut;
}

/// Initial guess wave function
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-1.5*sqrt(x*x+y*y+z*z+0.001)); // Deliberately bad guess
}

/// Given psi and V evaluate the energy
double energy(World& world, const functionT& psi, const functionT& potn) {
    double S = psi.inner(psi);
    double PE = psi.inner(psi*potn);
    double KE = 0.0;
    for (int axis=0; axis<3; axis++) {
        functionT dpsi = diff(psi,axis);
        KE += inner(dpsi,dpsi)*0.5;
    }
    double E = (KE+PE)/S;

    if (world.rank() == 0) {
        print("the overlap integral is",S);
        print("the kinetic energy integral",KE);
        print("the potential energy integral",PE);
        print("the total energy",E);
    }
    return E;
}


struct realexp {
    void operator()(const Key<3>& key, Tensor<double>& t) const {
        UNARY_OPTIMIZED_ITERATOR(double, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


/// Evolve the wave function in imaginary time to converge to ground state
void converge(World& world, functionT& potn, functionT& psi, double& eps) {
    // The imaginary time Green's function is just the heat kernel.
    // In 1D it is 1/sqrt(2*pi*t) * exp(-x^2/2t)

    // We will evolve with the simple Trotter form exp(-T*t/2) exp(-V*t) exp(-T*t/2)
    // so to keep things vaguely stable we want max(exp(-V*t)) (which is at the origin)
    // to be circa 10.0.  Hence our largest time step is -ln(10)/V(0).  As we 
    // approach convergence we shrink the time step to the desired precision.

    double tmax = -2.3/V(coordT(0.0));
    if (world.rank() == 0) print("tmax", tmax);

    functionT expV = (-tmax)*potn;
    expV.unaryop(realexp());

    Tensor<double> coeff(1); coeff[0] = 1.0/sqrt(2*constants::pi*tmax);
    Tensor<double> expnt(1); expnt[0] = 0.5/tmax;
    operatorT* op = new operatorT(world,k,coeff,expnt);

    for (int iter=0; iter<200; iter++) {
        if (world.rank() == 0) print("ITER",iter);
        psi = apply(*op, expV*psi);
        psi.truncate();
        double norm = psi.norm2();
        if (world.rank() == 0) print("new norm", norm);
        psi.scale(1.0/norm);
        eps = energy(world, psi,potn);
        if (((iter+1)%50)==0 && tmax>=0.001) {
            tmax *= 0.5;
            expV = (-tmax)*potn;
            expV.unaryop(realexp());

            delete op;
            coeff[0] = 1.0/sqrt(2*constants::pi*tmax);
            expnt[0] = 0.5/tmax;
            op = new operatorT(world,k,coeff,expnt);
        }        
    }
    delete op;
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);

    FunctionDefaults<3>::k = k;                 // Wavelet order
    FunctionDefaults<3>::thresh = thresh;         // Accuracy
    FunctionDefaults<3>::refine = true;         // Enable adaptive refinement
    FunctionDefaults<3>::initial_level = 2;     // Initial projection level
    for (int i=0; i<3; i++) {
        FunctionDefaults<3>::cell(i,0) = -L;   // User-defined volume
        FunctionDefaults<3>::cell(i,1) =  L;
    }

    if (world.rank() == 0) print("V(0)", V(coordT(0)));

    functionT psi = factoryT(world).f(guess);
    functionT potn = factoryT(world).f(V);
    psi.truncate();

    energy(world, psi, potn);

    double eps = -0.5;
    converge(world, potn, psi, eps);

    world.gop.fence();

    MPI::Finalize();
    return 0;
}
