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
static const double energy_exact = -6.1887887757287948; // From Maple, probably OK to 1e-14
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision
static const double velocity = 0.0;
//static const double eshift = energy_exact + 0.5*velocity*velocity; // Use this value to remove rotating phase
static const double eshift = 0.0;

static double current_time = 0.0; // Lazy but easier than making functors for everything

// Position of atom at current time
double atom_position() {
    return x0 + velocity*current_time;
}

// Exact solution ... (H-E)psi is accurate to 2e-7 or better inside Maple
double_complex psi_exact(const coordT& r) {
    const double x = r[0] - atom_position();

    if (fabs(x) > 8.0) return 0.0;

    const double xsq = x*x;
    
    // Horner's form for stability

    const double psi = exp(-1.30*xsq)*(-1.02151632756275513018719494826+(.522210612113707231536059971069+(-.378478352719362210706386739834+(.183732263756009855463432582593+(-0.866826311335724362186706464428e-1+(0.364601910940641762284283418688e-1+(-0.144289291226899801775738667691e-1+(0.536464813679927807734520149659e-2+(-0.188945345474975346004237994967e-2+(0.628725522158030920219217207400e-3+(-0.195986657875763402765072721284e-3+(0.563993909330309881048156461300e-4+(-0.147273758530730072646826354339e-4+(0.343202525037692780236348165792e-5+(-7.03765391498970506694108123682e-7+(1.25577395245191089173671652503e-7+(-1.93270666918809750376161513191e-8+(2.54624395753990033043923369519e-9+(-2.84968491109847694778452937732e-10+(2.68398879018928855922002181223e-11+(-2.09811331054703124038096336404e-12+(1.32869596552058891546970129944e-13+(-6.47453843054578193348503832223e-15+(2.08146181239264250219157355910e-16+(-8.27692933283146365185698371000e-19+(-4.21076100567125673420604632290e-19+(3.34873683682880953223636900553e-20+(-1.62840449033428572820374528252e-21+(5.80781234060753332169593869292e-23+(-1.58754665363960354625758075480e-24+(3.34792415609741263184593450566e-26+(-5.37727687523701580992550859153e-28+(6.37706272777548158355242766766e-30+(-5.27154378837014394526580093435e-32+(2.71430533655911926094030668033e-34-6.55694230766452931026127498540e-37*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq);

    // Galilean translation factor
    const double arg = x*velocity - (energy_exact + velocity*velocity*0.5 - eshift)*current_time;
    const double_complex tranfac = exp(double_complex(0,arg));
    return psi*tranfac;
}

// Time-dependent potential for translating atom
double V(const coordT& r) {
    const double x = r[0] - atom_position();
    return -8.0*exp(-x*x) - eshift;
}


// // Evolve the wave function one step in real time using Trotter.
// // G0 = G0(tstep/2)
// complex_functionT propagate(World& world, const complex_operatorT& G0, const complex_functionT& psi0, const double tstep) {
//     functionT zdip = factoryT(world).f(zdipole);

//     if (world.rank() == 0) print("truncating");
//     psi.truncate();
//     if (world.rank() == 0) print("initial normalize");    
//     psi.scale(1.0/psi.norm2());
//     int steplo = -1;
//     for (int step=steplo; step<nstep; step++) {
//         if (world.rank() == 0) print("\nStarting time step",step,tstep*step);

//         energy(world, psi, potn+laser((step)*tstep)*zdip);  // Use field at current time to evaluate energy

//         double_complex dipole = inner(psi,zdip*psi);
//         if (world.rank() == 0) print("THE DIPOLE IS", dipole);

//         if (step > steplo) {
//             if (world.rank() == 0) print("load balancing");
//             LoadBalImpl<3> lb(psi, lbcost<double_complex,3>);
//             lb.load_balance();
//             FunctionDefaults<3>::set_pmap(lb.load_balance());
//             world.gop.fence();
//             psi = copy(psi, FunctionDefaults<3>::get_pmap(), false);
//             potn = copy(potn, FunctionDefaults<3>::get_pmap(), true);
//             zdip = copy(zdip, FunctionDefaults<3>::get_pmap(), true);
//         }

//         if (world.rank() == 0) print("making expV");
//         functionT totalpotn = potn + laser((step+0.5)*tstep)*zdip; // Use midpoint field to advance in time
//         //totalpotn.refine();
//         totalpotn.reconstruct();
//         complex_functionT expV = double_complex(0.0,-tstep)*totalpotn;
//         expV.unaryop(unaryexp<double_complex,3>());
//         expV.truncate();
//         //expV.refine();
//         expV.reconstruct();
//         double expVnorm = expV.norm2();
//         if (world.rank() == 0) print("expVnorm", expVnorm);

//         long sz = psi.size();
//         //psi.refine();
//         if (world.rank() == 0) print("applying operator 1",sz);
//         psi = apply(G, psi);
//         psi.truncate();
//         //psi.refine();

//         sz = psi.size();
//         if (world.rank() == 0) print("multipling by expV",sz);
//         psi = expV*psi;
//         psi.truncate();
//         //psi.refine();
//         sz = psi.size();
//         if (world.rank() == 0) print("applying operator 2",sz);
//         psi = apply(G, psi);
//         psi.truncate();
//         sz = psi.size();
//         double norm = psi.norm2();
//         if (world.rank() == 0) print(step, step*tstep, norm,sz);
//     }
// }


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
    double ctarget = 10.0; // Estimated from FT of potential
    double c = 1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);

    print("Critical time step is", tcrit);

    complex_functionT psi = complex_factoryT(world).f(psi_exact);
    functionT potn = factoryT(world).f(V);

    complex_functionT dpsi = diff(psi,0);
    double ke = inner(dpsi,dpsi).real() * 0.5;
    double pe = psi.inner(psi*potn).real();
    print(ke, pe, pe+ke);

    print("initial norm", psi.norm2());

    double tstep = 10.0*tcrit;

    print(tstep/tcrit);

    // This section does Trotter
    //SeparatedConvolution<double_complex,1> G0 = qm_free_particle_propagator<3>(world, k, c, tstep, 2*L);

    world.gop.fence();

    finalize();
    return 0;
}

    
