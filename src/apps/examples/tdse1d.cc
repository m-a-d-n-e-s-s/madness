/// \file tdse1d.cc
/// \brief Example propagation of TDSE (translating atom) using various propagators


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

// typedefs to make life less verbose
typedef Vector<double,1> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,1> > functorT;
typedef Function<double,1> functionT;
typedef FunctionFactory<double,1> factoryT;
typedef SeparatedConvolution<double,1> operatorT;

typedef SharedPtr< FunctionFunctorInterface<double_complex,1> > complex_functorT;
typedef Function<double_complex,1> complex_functionT;
typedef FunctionFactory<double_complex,1> complex_factoryT;
typedef SeparatedConvolution<double_complex,1> complex_operatorT;

// Simulation parameters
static const double L = 100.0; // Simulation in [-L,L]
static const double x0 = -L + 10.0; // Initial position of the atom
static const double energy_exact = -6.188788775728796797594788; // From Maple
static const long k = 12;        // wavelet order
static const double thresh = 1e-10; // precision
static const double velocity = 3.0;
//static const double eshift = energy_exact - 0.5*velocity*velocity; // Use this value to remove rotating phase
static const double eshift = 0.0;

static double current_time = 0.0; // Lazy but easier than making functors for everything

// Position of atom at current time
double atom_position() {
    return x0 + velocity*current_time;
}

// Exact solution ... (H-E)psi is accurate to 2e-7 or better inside Maple
double_complex psi_exact(const coordT& r) {
    const double x = r[0] - atom_position();

    if (fabs(x) > 9.0) return 0.0;

    const double xsq = x*x;
    
    // Horner's form for stability ... yes, it is 70-order polyn ... don't panic ... all is OK
    const double psi = exp(-1.30*xsq)*(-1.02151632756275513018719494826+(.522210612113707231536059971069+(-.378478352719362210706386739834+(.183732263756009855463432582593+(-0.866826311335724362186706464428e-1+(0.364601910940641762284283418688e-1+(-0.144289291226899801775738667691e-1+(0.536464813679927807734520149659e-2+(-0.188945345474975346004237994967e-2+(0.628725522158030920219217207400e-3+(-0.195986657875763402765072721284e-3+(0.563993909330309881048156461300e-4+(-0.147273758530730072646826354339e-4+(0.343202525037692780236348165792e-5+(-7.03765391498970506694108123682e-7+(1.25577395245191089173671652503e-7+(-1.93270666918809750376161513191e-8+(2.54624395753990033043923369519e-9+(-2.84968491109847694778452937732e-10+(2.68398879018928855922002181223e-11+(-2.09811331054703124038096336404e-12+(1.32869596552058891546970129944e-13+(-6.47453843054578193348503832223e-15+(2.08146181239264250219157355910e-16+(-8.27692933283146365185698371000e-19+(-4.21076100567125673420604632290e-19+(3.34873683682880953223636900553e-20+(-1.62840449033428572820374528252e-21+(5.80781234060753332169593869292e-23+(-1.58754665363960354625758075480e-24+(3.34792415609741263184593450566e-26+(-5.37727687523701580992550859153e-28+(6.37706272777548158355242766766e-30+(-5.27154378837014394526580093435e-32+(2.71430533655911926094030668033e-34-6.55694230766452931026127498540e-37*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq);

    // Galilean translation factor
    const double arg = x*velocity - (energy_exact - velocity*velocity*0.5 - eshift)*current_time;
    const double_complex tranfac = exp(double_complex(0,arg));
    return psi*tranfac;
}

// Time-dependent potential for translating atom
double V(const coordT& r) {
    const double x = r[0] - atom_position();
    if (fabs(x) > 5.5) return 0.0;

    return -8.0*exp(-x*x) - eshift;
}

template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Evolve forward one time step using Trotter ... G = G0(tstep/2)
complex_functionT trotter(World& world, const complex_operatorT& G, const complex_functionT& psi0, const double tstep) {
    complex_functionT psi = apply(G, psi0);
    psi.truncate();
    current_time += 0.5*tstep;

    functionT potn = factoryT(world).f(V).truncate_on_project();
    potn.truncate();
    complex_functionT expV = double_complex(0.0,-tstep)*potn;
    expV.unaryop(unaryexp<double_complex,1>());

    psi = expV*psi; psi.truncate();
    psi = apply(G,psi);
    psi.truncate();
    
    current_time += 0.5*tstep;

    return psi;
}

void print_info(World& world, const complex_functionT& psi, int step) {
    functionT potn = factoryT(world).f(V).truncate_on_project();
    complex_functionT dpsi = diff(psi,0);
    double ke = inner(dpsi,dpsi).real() * 0.5;
    double pe = psi.inner(psi*potn).real();
    double norm = psi.norm2();

    complex_functionT psiX = complex_factoryT(world).f(psi_exact);
    double err = (psiX - psi).norm2();
    //abs(inner(psiX,psi));//

    if ((step%40) == 0) {
        printf("\n");
        printf(" step    time      atom x        norm        kinetic    potential      energy      err norm\n");
        printf("------ -------- ------------ ------------ ------------ ------------ ------------ ------------\n");
    }
    printf("%6d %8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", 
           step, current_time, atom_position(), norm, ke, pe, ke+pe, err);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    cout.precision(8);

    FunctionDefaults<1>::set_k(k);                 // Wavelet order
    FunctionDefaults<1>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<1>::set_refine(false);        // Disable adaptive refinement
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_initial_level(8);     // Initial projection level
    FunctionDefaults<1>::set_truncate_mode(1);

    // Estimate the bandwidth and largest practical time step using G0
    double ctarget = 20.0; // Estimated from FT of potential
    double c = 1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);

    print("Critical time step is", tcrit);

    complex_functionT psi = complex_factoryT(world).f(psi_exact);
    psi.truncate();

    double tstep = tcrit*3; // This choice for Trotter

    int nstep = velocity==0 ? 100 : (L - 10 - x0)/velocity/tstep;

    print("No. of time steps is", nstep);

    // This section does Trotter
    SeparatedConvolution<double_complex,1> G0 = qm_free_particle_propagator<1>(world, k, c, 0.5*tstep, 2*L);
    for (int step=0; step<nstep; step++) {
        print_info(world, psi,step);
        psi = trotter(world, G0, psi, tstep);
    }

    world.gop.fence();

    finalize();
    return 0;
}

    
