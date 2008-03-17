/// \file tdse.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

// Convenient but sleazy use of global variables to define simulation parameters
static const double L = 23.0;
static const long k = 12;        // wavelet order
static const double thresh = 1e-6; // precision
static const double cut = 0.1;  // smoothing parameter for 1/r
static const double F = 1.0;  // Laser field strength
static const double omega = 1.0; // Laser frequency

// typedefs to make life less verbose
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;

typedef SharedPtr< FunctionFunctorInterface<double_complex,3> > complex_functorT;
typedef Function<double_complex,3> complex_functionT;
typedef FunctionFactory<double_complex,3> complex_factoryT;
typedef SeparatedConvolution<double_complex,3> complex_operatorT;

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
    return -smoothed_potential(rr/cut)/cut;
}

/// Initial guess wave function
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-1.0*sqrt(x*x+y*y+z*z+cut*cut)); // Change 1.0 to 0.6 to make bad guess
}

/// z-dipole
double zdipole(const coordT& r) {
    return r[2];
}

/// Strength of the laser field at time t ... 1 full cycle
double laser(double t) {
    double omegat = omega*t;
    if (omegat < 0.0 || omegat > 1) return 0.0;
    else return F*sin(2*constants::pi*omegat);
}

/// Given psi and V evaluate the energy
template <typename T>
double energy(World& world, const Function<T,3>& psi, const functionT& potn) {
    T S = psi.inner(psi);
    T PE = psi.inner(psi*potn);
    T KE = 0.0;
    for (int axis=0; axis<3; axis++) {
        Function<T,3> dpsi = diff(psi,axis);
        KE += inner(dpsi,dpsi)*0.5;
    }
    T E = (KE+PE)/S;

    if (world.rank() == 0) {
        print("the overlap integral is",S);
        print("the kinetic energy integral",KE);
        print("the potential energy integral",PE);
        print("the total energy",E);
    }
    return -std::abs(E); // ugh
}

template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

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
    expV.unaryop(unaryexp<double,3>());

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
            expV.unaryop(unaryexp<double,3>());

            delete op;
            coeff[0] = 1.0/sqrt(2*constants::pi*tmax);
            expnt[0] = 0.5/tmax;
            op = new operatorT(world,k,coeff,expnt);
        }        
    }
    delete op;
}

/// Evolve the wave function in real time
void propagate(World& world, functionT& potn, functionT& psi0, double& eps) {
    // In the absense of a time-dependent potential we should just have the
    // rotating phase of the ground state wave function

    double ctarget = constants::pi/(0.5*cut);
    double c = 1.86*ctarget; //1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);
    double tstep = 10.0*tcrit;
    int nstep = 100.0/tstep;
    tstep = 0.005;

    if (world.rank() == 0) {
        print("bandlimit",ctarget,"effband",c,"tcrit",tcrit,"tstep",tstep,"nstep",nstep);
    }

    complex_functionT psi = double_complex(sqrt(0.5),sqrt(0.5))*psi0;
    SeparatedConvolution<double_complex,3> G = qm_free_particle_propagator<3>(world, k, c, 0.5*tstep, 2*L);

    functionT zdip = factoryT(world).f(zdipole);

    if (world.rank() == 0) print("truncating");
    psi.truncate();
    if (world.rank() == 0) print("initial normalize");    
    psi.scale(1.0/psi.norm2());
    for (int step=-4.0; step<nstep; step++) {
        if (world.rank() == 0) print("\nStarting time step",step,tstep*step);

        energy(world, psi, potn+laser(step*tstep)*zdip);  // Use field at current time to evaluate energy

        double_complex dipole = inner(psi,zdip*psi);
        if (world.rank() == 0) print("THE DIPOLE IS", dipole);

        if (world.rank() == 0) print("load balancing");
        LoadBalImpl<3> lb(psi, lbcost<double_complex,3>);
        lb.load_balance();
        FunctionDefaults<3>::pmap = lb.load_balance();
        world.gop.fence();
        psi = copy(psi, FunctionDefaults<3>::pmap, false);
        potn = copy(potn, FunctionDefaults<3>::pmap, true);
        zdip = copy(zdip, FunctionDefaults<3>::pmap, true);

        if (world.rank() == 0) print("making expV");
        functionT totalpotn = potn + laser((step+0.5)*tstep)*zdip; // Use midpoint field to advance in time
        totalpotn.refine();
        complex_functionT expV = double_complex(0.0,-tstep)*totalpotn;
        expV.unaryop(unaryexp<double_complex,3>());
        expV.truncate();
        expV.refine();
        double expVnorm = expV.norm2();
        if (world.rank() == 0) print("expVnorm", expVnorm);

        long sz = psi.size();
        psi.refine();
        if (world.rank() == 0) print("applying operator 1",sz);
        psi = apply(G, psi);
        psi.truncate();
        psi.refine();
        sz = psi.size();
        if (world.rank() == 0) print("multipling by expV",sz);
        psi = expV*psi;
        psi.truncate();
        psi.refine();
        sz = psi.size();
        if (world.rank() == 0) print("applying operator 2",sz);
        psi = apply(G, psi);
        psi.truncate();
        sz = psi.size();
        double norm = psi.norm2();
        if (world.rank() == 0) print(step, step*tstep, norm,sz);
    }
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
    //converge(world, potn, psi, eps);
    propagate(world, potn, psi, eps);

    world.gop.fence();

    MPI::Finalize();
    return 0;
}
