
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/mra/mra.h>
#include <iostream>
#include <vector>

// These forward declarations are here for SpectralPropagator::step() in spectralprop.h.
class Fred;
namespace madness {
    double distance(const Fred& a, const Fred& b);
    double distance(madness::Function<std::complex<double>, 1ul>& a, madness::Function<std::complex<double>, 1ul>& b);
}

#include "spectralprop.h"

using namespace madness;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// test0 is a simple scalar problem
//
// -u + u^2 = du/dt    u(0)=-0.1
//
// test1 is the same differential equation
// but simultaneously propagating 3 different
// initial conditions (-0.1,-0.2,-0.3) to
// illustrate use of a custom class.
//
// test2 is the same differential equation
// as the heat2.cc example
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/// Wrapper around vector demonstrating interface necessary
class Fred {
    std::vector<double> v;

    // Default constructor not needed
    Fred();

public:
    // Specific to this example
    Fred(double a, double b, double c)
        : v(3) {v[0] = a; v[1] = b; v[2] = c;}
    void set(int i, double a) {v[i] = a;}
    double get(int i) const {return v[i];}

    // The following interfaces are ncessary

    // Copy constructor
    Fred(const Fred& f) : v(f.v) {}

    // Assignment
    Fred& operator=(const Fred& f) {
        if (this != &f) v = f.v;
        return *this;
    }

    // Inplace addition
    Fred& operator+=(const Fred& f) {
        for (int i=0; i<3; i++) v[i] += f.v[i];
        return *this;
    }

    // Scale from right by double
    Fred operator*(const double& d) const {
        Fred f(*this);
        for (int i=0; i<3; i++) f.v[i] *= d;
        return f;
    }
};

namespace madness {
    // distance is in madness namespace to disambiguate it from std::distance.
    double distance(const Fred& a, const Fred& b)
    {
        double xsq = 0.0;
        for (int i=0; i<3; i++) {
            double xx = a.get(i) - b.get(i);
            xsq += xx*xx;
        }
        return std::sqrt(xsq);
    }
} // namespace madness


// Not required
std::ostream& operator<<(std::ostream& s, const Fred& f) {
    s << "[" << f.get(0) << ", " << f.get(1) << ", " << f.get(2) << "]";
    return s;
}

double expL_double(double dt, double u) {
    return exp(-dt)*u;
}

double N_double(double t, const double u) {
    return u*u;
}

Fred expL(double dt, const Fred& f) {
    Fred r(f);
    for (int i=0; i<3; i++) r.set(i,exp(-dt)*f.get(i));
    return r;
}

Fred N(double t, const Fred& f) {
    Fred r(f);
    for (int i=0; i<3; i++) r.set(i,f.get(i)*f.get(i));
    return r;
}

double exact(double t, double u0) {
    return -u0/(-u0+std::exp(t)*u0-std::exp(t));
}

void test0(World& world) {
    if (world.rank() == 0) {
        print("Testing Gauss Legendre --- double");
        const double t = 0.0;
        const double dt = 0.1;
        double u0 = -0.1;

        // Solutions from maple
        double M[4] = {
             0.0,
            -0.089631631149550335419,
            -0.089630792118101842108,
            -0.089630792044157819486
        };

        for (int NPT=1; NPT<=3; NPT++) {
            SpectralPropagator P(NPT);
            double u = P.step(t, dt, u0, expL_double, N_double);

            // Error relative to numerical solution computed by Maple
            // ... it should be near machine precision
            double err = std::fabs(u - M[NPT]);

            // Error relative to exact solution ... should converge
            // rapidly as order is increased.
            double errexact = std::fabs(u - exact(t+dt, u0));

            print("NPT =", NPT, "err(Maple)", err, "err(exact)", errexact);
        }
    }
}

void test0GaussLobatto(World& world) {
    if (world.rank() == 0) {
        print("Testing Gauss Lobatto --- double");
        const double dt = 0.1;
        double u0 = -0.1;

        for (int NPT=2; NPT<=6; NPT++) {
            double t = 0.0;
            SpectralPropagatorGaussLobatto P(NPT);
            double u = u0;
            for (int step=0; step<11; step++) {
                u = P.step(t, dt, u, expL_double, N_double, 1e-12, false);
                t += dt;
            }
            double errexact = u - exact(t,u0);
            print("NPT =", NPT, "err(exact t=11dt)", errexact);
        }
    }
}

void test1(World& world) {
    if (world.rank() == 0) {
        print("Testing Gauss Legendre --- user-defined type");
        const double t = 0.0;
        const double dt = 0.1;
        Fred u0(-0.1,-0.2,-0.3);

        // Solutions from maple
        Fred M[4] = {
            Fred(0.0,0.0,0.0),
            Fred(-0.089631631149550335419,-0.17759183895709322265,-0.26392868193278156650),
            Fred(-0.089630792118101842108,-0.17758754617993929487,-0.26391672762980049068),
            Fred(-0.089630792044157819486,-0.17758754573188727124,-0.26391672628171338046)
        };

        for (int NPT=1; NPT<=3; NPT++) {
            SpectralPropagator P(NPT);
            Fred u = P.step(t, dt, u0, expL, N);

            // Error relative to numerical solution computed by Maple
            // ... it should be near machine precision
            double err = 0.0;

            // Error relative to exact solution ... should converge
            // rapidly as order is increased.
            double errexact = 0.0;
            for (int i=0; i<3; i++) {
                err += std::fabs(u.get(i) - M[NPT].get(i));
                errexact += std::fabs(u.get(i) - exact(t+dt, u0.get(i)));
            }

            print("NPT =", NPT, "err(Maple)", err, "err(exact)", errexact);
        }

        for (int NPT=1; NPT<=3; NPT++) {
            double t = 0.0;
            SpectralPropagator P(NPT);
            Fred u = u0;
            for (int step=0; step<11; step++) {
                u = P.step(t, dt, u, expL, N);
                t += dt;
            }
            double errexact = 0.0;
            for (int i=0; i<3; i++) {
                errexact += std::fabs(u.get(i) - exact(t, u0.get(i)));
            }

            print("NPT =", NPT, "err(exact t=11dt)", errexact);

        }
    }
}

///////////////////////////////////////////////////
// All this crap for the non-linear TDSE problem //
///////////////////////////////////////////////////

#include <madness/mra/mra.h>
#include <madness/mra/qmprop.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>

namespace madness {

double distance(madness::Function<std::complex<double>, 1ul>& a, madness::Function<std::complex<double>, 1ul>& b) {
    return (a-b).norm2();
}
}



static const double L = 100.0; // Simulation in [-L,L]
static const double x0 = -L + 10.0; // Initial position of the atom
static const double energy_exact = -6.188788775728796797594788; // From Maple
static const double velocity = 3.0;
static const long k = 12;        // wavelet order
static const double thresh = 1e-8; // precision
static const double eshift = 0.0;
//static const double eshift = energy_exact - 0.5*velocity*velocity; // Use this value to remove rotating phase

// Estimate the bandwidth and largest practical time step using G0
static double ctarget = 20.0; // Estimated from FT of exact solution ... was 20
static double c = 1.86*ctarget;
static double tcrit = 2*constants::pi/(c*c);

typedef Convolution1D<double_complex> complex_operatorT;
typedef std::shared_ptr<complex_operatorT> pcomplex_operatorT;

// Position of atom at given time
double atom_position(double current_time) {
    return x0 + velocity*current_time;
}

class PsiExact : public FunctionFunctorInterface<double_complex,1>
{
    const double current_time;
    const double atom_x;
public:
    PsiExact(double current_time)
        : current_time(current_time)
        , atom_x(atom_position(current_time))
    {}

    // Exact solution ... (H-E)psi is accurate to 2e-7 or better inside Maple
    double_complex operator()(const coord_1d& r) const {
        const double x = r[0] - atom_x;

        if (fabs(x) > 9.0) return 0.0;

        const double xsq = x*x;

        // Horner's form for stability ... yes, it is 70-order polyn ... don't panic ... all is OK

        const double psi = exp(-1.30*xsq)*(-1.02151632756275513018719494826+(.522210612113707231536059971069+(-.378478352719362210706386739834+(.183732263756009855463432582593+(-0.866826311335724362186706464428e-1+(0.364601910940641762284283418688e-1+(-0.144289291226899801775738667691e-1+(0.536464813679927807734520149659e-2+(-0.188945345474975346004237994967e-2+(0.628725522158030920219217207400e-3+(-0.195986657875763402765072721284e-3+(0.563993909330309881048156461300e-4+(-0.147273758530730072646826354339e-4+(0.343202525037692780236348165792e-5+(-7.03765391498970506694108123682e-7+(1.25577395245191089173671652503e-7+(-1.93270666918809750376161513191e-8+(2.54624395753990033043923369519e-9+(-2.84968491109847694778452937732e-10+(2.68398879018928855922002181223e-11+(-2.09811331054703124038096336404e-12+(1.32869596552058891546970129944e-13+(-6.47453843054578193348503832223e-15+(2.08146181239264250219157355910e-16+(-8.27692933283146365185698371000e-19+(-4.21076100567125673420604632290e-19+(3.34873683682880953223636900553e-20+(-1.62840449033428572820374528252e-21+(5.80781234060753332169593869292e-23+(-1.58754665363960354625758075480e-24+(3.34792415609741263184593450566e-26+(-5.37727687523701580992550859153e-28+(6.37706272777548158355242766766e-30+(-5.27154378837014394526580093435e-32+(2.71430533655911926094030668033e-34-6.55694230766452931026127498540e-37*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq);

        // Galilean translation factor
        const double arg = x*velocity - (energy_exact - velocity*velocity*0.5 - eshift)*current_time;
        const double_complex tranfac = exp(double_complex(0,arg));
        return psi*tranfac;
    }
};

class Vnuclear : public FunctionFunctorInterface<double,1>
{
    const double atom_x;
public:
    Vnuclear(double current_time) : atom_x(atom_position(current_time)) {}

    double operator()(const coord_1d& r) const {
        // Time-dependent potential for translating atom
        const double x = r[0] - atom_x;
        if (fabs(x) > 6.2) return 0.0;

        return -8.0*exp(-x*x) - eshift;
    }
};

complex_function_1d APPLY(complex_operatorT* q1d, const complex_function_1d& psi) {
    psi.reconstruct();

    psi.broaden();
    psi.broaden();
    psi.broaden();
    psi.broaden();
    psi.broaden();

    complex_function_1d r = apply_1d_realspace_push(*q1d, psi, 0);
    r.sum_down();
    return r;
}

complex_operatorT* MAKE_PROPAGATOR(World& world, double t)
{
    return qm_1d_free_particle_propagator(k, c, t, 2.0*L);
}

template <typename T>
std::pair<bool, T> find_fuzzy(double t, const std::vector< std::pair<double,T> >& cache) {
    for (unsigned int i=0; i<cache.size(); i++) {
        if (std::fabs(cache[i].first - t) < 1e-12*t) {
            return std::make_pair(true, cache[i].second);
        }
    }
    return std::make_pair(false, T());
}

static std::vector< std::pair<double, complex_operatorT*> > G_cache;
complex_operatorT* G(World& world, double dt) {
    std::pair<bool, complex_operatorT*> r = find_fuzzy(dt, G_cache);
    complex_operatorT* p;
    if (r.first) {
        p = r.second;
    }
    else {
        //print("making operator for dt =", dt);
        G_cache.push_back(std::make_pair(dt, p = MAKE_PROPAGATOR(world, dt)));
    }
    return p;
}

std::vector< std::pair<double, real_function_1d> > vnuc_cache;
real_function_1d vnuc(World& world, double t) {
    std::pair<bool,real_function_1d> r = find_fuzzy(t, vnuc_cache);
    real_function_1d vnuc;
    if (r.first) {
        vnuc = r.second;
    }
    else {
        //print("making vnuc for t =", t);
        vnuc_cache.push_back(std::make_pair(t, vnuc = real_factory_1d(world).functor(real_functor_1d(new Vnuclear(t)))));
    }
    return vnuc;
}

complex_function_1d applyexpLt(double dt, const complex_function_1d& u) {
    if (dt < 1e-12) return copy(u);
    return APPLY(G(u.world(), dt), u).truncate();
}

complex_function_1d applyN(double t, const complex_function_1d& u) {
    return (u*(vnuc(u.world(), t)*double_complex(0.0,-1.0))).truncate(); ///  + conj(u)*u
}

void print_info(World& world, double current_time, const complex_function_1d& psi, int step) {
    real_function_1d potn = vnuc(world, current_time);
    complex_derivative_1d D = free_space_derivative<double_complex,1>(world,0);
    complex_function_1d dpsi = D(psi);
    double ke = inner(dpsi,dpsi).real() * 0.5;
    double pe = psi.inner(psi*(double_complex(1,0)*potn)).real(); // +psi*conj(psi)
    double norm = psi.norm2();
    double err = psi.err(PsiExact(current_time));
    ke /= norm*norm;
    pe /= norm*norm;

    if (world.rank() > 0) return;
    if ((step%40) == 0) {
        printf("\n");
        printf(" step    time      atom x              norm               kinetic              potential             energy              err norm     depth   size  \n");
        printf("------ -------- ------------     ----------------     ----------------     ----------------     ----------------     ---------------- ----  --------\n");
    }
    printf("%6d %8.4f %12.8f %20.13f %20.13f %20.13f %20.13f %20.13f %4d %9ld\n",
           step, current_time, atom_position(current_time), norm, ke, pe, ke+pe, err, int(psi.max_depth()), long(psi.size()));
}


void test2(World& world) {
    if (world.rank() == 0) print("Testing Gauss Legendre --- TDSE 1D example");
    std::cout.precision(8);
    FunctionDefaults<1>::set_k(k);                 // Wavelet order
    FunctionDefaults<1>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<1>::set_autorefine(false);
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_initial_level(8);     // Initial projection level
    FunctionDefaults<1>::set_truncate_mode(1);

    complex_function_1d psi0 = complex_factory_1d(world).functor(complex_functor_1d(new PsiExact(0.0)));
    psi0.truncate();
    complex_function_1d psi = copy(psi0);

    // These should be input
    double tstep = 10.0*tcrit;
    int NPT = 3;

    world.gop.broadcast(NPT);
    world.gop.broadcast(tstep);

    int nstep = (velocity==0) ? 100 : int((L - 10 - x0)/velocity/tstep);

    if (world.rank() == 0) {
        print("");
        print(" Wavelet order", k);
        print("     Threshold", thresh);
        print("      Velocity", velocity);
        print("     Time step", tstep);
        print("      No.steps", nstep);
        print("           NPT", NPT);
    }

    SpectralPropagator P(NPT);
    double t = 0.0;
    for (int step=0; step<nstep; step++) {
        //if ((step%10) == 0)
        print_info(world, t, psi, step);
        if (velocity==0) print("PHASE", inner(psi,psi0), exp(double_complex(0.0,-energy_exact*t)));
        psi = P.step(t, tstep, psi, applyexpLt, applyN, thresh*10.0); //, true, false);
        psi.truncate();

        t += tstep;
        vnuc_cache.clear();
    }
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    test0(world);
    print("GL");
    test0GaussLobatto(world);
    test1(world);
    test2(world);

    finalize();
    return 0;
}




