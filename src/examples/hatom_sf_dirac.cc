#include <fstream>
#include <iostream>
#include <algorithm>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/lbdeux.h>
#include <madness/mra/nonlinsol.h>
#include <madness/chem/atomutil.h>

using namespace madness; 

static const double Z = 80.0;      // nuclear charge
static const double L = 40.0/Z;    // L=40/Z [-L,L] box size so exp(-Zr)=1e-16 padded a bit since we are masking
static const long k = 8;           // wavelet order
static const double thresh = 1e-6; // precision
static const double c = 137.035999679; // speed of light
static const int Vtruncmode = (Z<=10 && thresh>1e-10) ? 1 : 0; // Point nucleus potential woes
static const double eprec = 0.0; //std::min(1e-4,thresh*Z*Z); // Note this is for NR wavefun, and for now we just want relative not absolute error
static const double rcut = 1e-4; // typical size of heavy nucleus //smoothing_parameter(Z, eprec);
static const double v = std::sqrt(1.0 - Z*Z/(c*c)) - 1.0;

inline double mask1(double x) {
    x = (x * x * (3. - 2. * x));
    return x;
}

static double mask3(const coord_3d& ruser) {
    coord_3d rsim;
    user_to_sim(ruser, rsim);
    double x = rsim[0], y = rsim[1], z = rsim[2];
    double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
    double rlo = 1.0 / lo;

    if (x < lo)
        result *= mask1(mask1(x * rlo));
    else if (x > hi)
        result *= mask1(mask1((1.0 - x) * rlo));
    if (y < lo)
        result *= mask1(mask1(y * rlo));
    else if (y > hi)
        result *= mask1(mask1((1.0 - y) * rlo));
    if (z < lo)
        result *= mask1(mask1(z * rlo));
    else if (z > hi)
        result *= mask1(mask1((1.0 - z) * rlo));

    return result;
}


struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0)
        : leaf_value(leaf_value)
        , parent_value(parent_value)
    {}

    double operator()(const Key<3>& key, const FunctionNode<double,3>& node) const {
        if (node.is_leaf()) return leaf_value;
        else return parent_value;
    }
};

// This class is used to store information for the non-linear solver
struct F {
    real_function_3d psi, phi;

    F(const real_function_3d& psi, const real_function_3d& phi) 
        : psi(psi), phi(phi)
    {}

    F operator-(const F& b) const {
        return F(psi-b.psi, phi-b.phi);
    }

    F operator+=(const F& b) { // Operator+= necessarphi
        psi += b.psi; phi += b.phi;
        return *this;
    }

    F operator*(double a) const { // Scale bphi a constant necessarphi
        return F(psi*a,phi*a);
    }
};

double inner(const F& a, const F& b) {
    return inner(a.psi,b.psi) + inner(a.phi,b.phi);
}

struct allocator {
    World& world;

    allocator(World& world) : world(world) {}

    F operator()() {
        return F(real_function_3d(world),real_function_3d(world));
    }
};

double exact_energy(double Z) {
    if (Z<=10) {
        double a2 = Z*Z/(c*c);
        double a4 = a2*a2;
        double a6 = a4*a2;
        double a8 = a4*a4;
        return Z*Z*(-1.0/2.0 - 1.0/8.0*a2 - 1/16*a4 - 5.0/128.0*a6 - 7.0/256.0*a8);
    }
    double Za = Z/c;
    double s = Za / sqrt(1.0 - Za*Za);
    return c*c/std::sqrt(1.0 + s*s) - c*c;
}

static double psi_guess(const coord_3d& r) {
    double R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    const double expnt = std::sqrt(-2*exact_energy(Z)*(1+exact_energy(Z)/(2*c*c)));
    return std::pow(R+rcut,v)*exp(-expnt*R)/sqrt(constants::pi/(expnt*expnt*expnt)); // note +rcut in singular part to smooth on nuclear scale
    //const double expnt = Z*0.5;
    //return exp(-expnt*R)/sqrt(constants::pi/(expnt*expnt*expnt)); // note +rcut in singular part to smooth on nuclear scale
}

static double phi_guess(const coord_3d& r) {
    double R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    const double expnt = std::sqrt(-2*exact_energy(Z)*(1+exact_energy(Z)/(2*c*c)));
    return exp(-expnt*R)/sqrt(constants::pi/(expnt*expnt*expnt));
    // const double expnt = Z;
    // return exp(-expnt*R)/sqrt(constants::pi/(expnt*expnt*expnt)); // note +rcut in singular part to smooth on nuclear scale
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    //return -Z/(sqrt(x*x+y*y+z*z));
    return (-Z/rcut)*smoothed_potential(sqrt(x*x+y*y+z*z)/rcut);
}

class Vderiv : public FunctionFunctorInterface<double,3> {
    const int axis;
public:
    Vderiv(int axis) : axis(axis) {}

    double operator()(const coord_3d& r) const {
        double R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double rc = 1.0/rcut;
        return - Z * (r[axis] / R) * dsmoothed_potential(R * rc) * (rc * rc);
    }
};

// Applies pV.p/4m^2c^2 phi = sum_q (-dV/dq dphi/dq - V d^2phi/dq^2)
real_function_3d pVp(World& world, const real_function_3d& V, const real_function_3d gradV[3], const real_function_3d& phi) {
    real_function_3d Wphi = real_function_3d(world);
    for (int axis=0; axis<3; axis++) {
        {
            real_derivative_3d D = free_space_derivative<double,3>(world, axis);
            D.set_bspline1();
            Wphi -= D(phi)*gradV[axis];
        }
        {
            real_derivative_3d D = free_space_derivative<double,3>(world, axis);
            D.set_bspline2();
            Wphi -= D(phi)*V;
        }
    }
    Wphi *= 1.0/(4.0*c*c);
    return Wphi;
}

// Applies pV.p/4m^2c^2
real_function_3d pVp(World& world, const real_function_3d& V, const real_function_3d& phi) {
    real_function_3d Wphi = real_function_3d(world);
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        D.set_bspline1(); // try bspline smoother derivative
        Wphi -= D(D(phi)*V);
    }
    Wphi *= 1.0/(4.0*c*c);
    return Wphi;
}

// E<psi|psi> = <psi|T|phi> + <psi|V|psi>
std::tuple<double,double,double,double> compute_energy(World& world,
                                                       real_function_3d& V,
                                                       real_function_3d& psi,
                                                       real_function_3d& phi)
{
    double overlap = inner(psi,psi);
    double potential_energy = inner(psi,V*psi) / overlap;
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis); // more accurate to use ABGV derivative for KE
        real_function_3d dpsi = D(psi);
        real_function_3d dphi = D(phi);
        kinetic_energy += 0.5*inner(dpsi,dphi);
    }
    kinetic_energy /= overlap;
    double total_energy = kinetic_energy + potential_energy;
    return {total_energy, kinetic_energy, potential_energy, overlap};
}

// (T-E) phi = -V psi + E(psi-phi)
// and with s=1+E/2mc^2
// T (psi-s*phi) = - pVp phi/4mc^2
// BAD!! (T-E) (psi-s*phi) = - pVp phi - E(psi-s*phi)
template <typename solverT>
std::tuple<real_function_3d, real_function_3d, double> iterate(World& world,
                                                               const real_function_3d& V,
                                                               const real_function_3d gradV[3],
                                                               const real_function_3d& mask,
                                                               const real_function_3d& psi,
                                                               const real_function_3d& phi,
                                                               const double energy,
                                                               const int iter,
                                                               solverT& solver)
{
    real_convolution_3d coulop = CoulombOperator(world, rcut*0.1, thresh);
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*energy), rcut*0.1, thresh);
    real_function_3d rhs = -2*(psi*V - energy*(psi-phi));
    rhs.truncate();
    real_function_3d phi_new = apply(op,rhs);

    double s = 1.0 + energy/(2*c*c);
    //rhs = -2*(pVp(world,V,gradV,phi) + energy*(psi-s*phi));
    rhs = -(2/(4*constants::pi))*pVp(world,V,phi);
    rhs.truncate();
    real_function_3d psimphi_new = apply(coulop,rhs);
    real_function_3d psi_new = s*phi_new + psimphi_new;

    psi_new.truncate();
    phi_new.truncate();

    psi_new = psi_new * mask;
    phi_new = phi_new * mask;
    
    double rnormL = (psi_new-psi).norm2();
    double rnormP = (phi_new-phi).norm2();
    if (world.rank() == 0) print(rnormL, rnormP);
    double rnorm = (psi_new-psi).norm2() + (phi_new-phi).norm2();

    // Comment out next 5 lines to not use the solver (i.e., to just do fixed-point iteration)
    if (iter>1) {
        F f = solver.update(F(psi,phi), F(psi_new-psi, phi_new-phi));
        psi_new = f.psi;
        phi_new = f.phi;
        psi_new.truncate();
        phi_new.truncate();
    }

    double stepnorm = (psi_new-psi).norm2() + (phi_new-phi).norm2();
    double maxstep = 1e-1;
    if (stepnorm > maxstep) {
        if (world.rank() == 0) print("step restriction", stepnorm);
        double step = maxstep/stepnorm;
        psi_new = psi + step*(psi_new-psi);
        phi_new = phi + step*(phi_new-phi);
    }

    double norm = inner(psi_new,psi_new);
    psi_new *= 1.0/sqrt(norm);
    phi_new *= 1.0/sqrt(norm);
    
    return {psi_new, phi_new, rnorm};
}

void run(World& world) {
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    if (world.rank() == 0) {
        print("            Z:", Z);
        print("            L:", L);
        print("            k:", k);
        print("            v:", v);
        print("            c:", c);
        print("         rcut:", rcut);
        print("        eprec:", eprec);
        print("       thresh:", thresh);
        print("       Vtmode:", Vtruncmode);
        print("     exact DC:", exact_energy(Z));
        print("");
    }

    {
        real_function_3d Vnuc = real_factory_3d(world).f(V).truncate_mode(Vtruncmode);
        LoadBalanceDeux<3> lb(world);
        lb.add_tree(Vnuc, LBCost(1.0,8.0));
        FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,false));
    }
    real_function_3d Vnuc = real_factory_3d(world).f(V).truncate_mode(Vtruncmode);

    real_function_3d gradV[3];
    for (int axis=0; axis<3; axis++) {
        std::shared_ptr<FunctionFunctorInterface<double,3>> p(new Vderiv(axis));
        gradV[axis] = real_factory_3d(world).functor(p).truncate_mode(0);
    }


    coord_3d lo = {0.0,0.0,-L}, hi = {0.0,0.0,L};
    const int npt = 1001;
    
    real_function_3d mask  = real_factory_3d(world).f(mask3);
    plot_line("mask.dat", npt, lo, hi, mask);
    
    real_function_3d psi  = real_factory_3d(world).f(psi_guess);
    real_function_3d phi  = real_factory_3d(world).f(phi_guess);
    psi = psi*mask;
    phi = phi*mask;
    {
        double norm;
        norm = psi.norm2();  psi.scale(1.0/norm);
        norm = phi.norm2();  phi.scale(1.0/norm);
    }

    XNonlinearSolver<F,double,allocator> solver = XNonlinearSolver<F,double,allocator>(allocator(world));
    solver.set_maxsub(10);
    
    double energy = -0.5*Z*Z; // NR energy guess
    for (int iter=0; iter<1000; iter++) {
        char fname[256];
        snprintf(fname,256,"psi-phi-%3.3d.dat", iter);
        plot_line(fname, npt, lo, hi, psi, phi);
        
        auto [psi_new, phi_new, rnorm] = iterate(world, Vnuc, gradV, mask, psi, phi, energy, iter, solver);
        psi = psi_new;
        phi = phi_new;

        auto psisize = psi.size();
        auto phisize = phi.size();
        auto phinorm = phi.norm2();
        auto [total_energy, kinetic_energy, potential_energy, overlap] = compute_energy(world, Vnuc, psi, phi);
        
        if (world.rank() == 0) {
            print("\niteration ", iter, " rnorm ", rnorm, " psi size ", psisize, " phi size ", phisize, " phi norm ", phinorm);
            print("            Kinetic energy  ", kinetic_energy);
            print(" Nuclear attraction energy ", potential_energy);
            print("              Total energy ", total_energy, total_energy - energy);
        }

        energy = total_energy;

        if (iter>10 and rnorm<10.0*thresh) break;
    }
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world,argc,argv);
    std::cout.precision(10);

    run(world); // Put all functions in separate scope to force destruction before finalize

    world.gop.fence();

    finalize();
    return 0;
}
