#include <fstream>
#include <iostream>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;

static const double Z = 10.0;      // nuclear charge
static const double L = 23.0/Z;    // [-L,L] box size
static const long k = 8;           // wavelet order
static const double thresh = 1e-6; // precision
static const double c = 137.035999679; // speed of light

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-Z*sqrt(x*x+y*y+z*z))/sqrt(constants::pi/(Z*Z*Z));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -Z/(sqrt(x*x+y*y+z*z));
}

// Applies pV.p/4m^2c^2
real_function_3d pVp(World& world, real_function_3d& V, real_function_3d& phi) {
    real_function_3d Wphi = real_function_3d(world);
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        D.set_bspline1(); // try bspline smoother derivative
        Wphi -= D(V*D(phi));
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
// (T-E) (psi-s*phi) = - pVp phi - E(psi-s*phi)
std::tuple<real_function_3d, real_function_3d, double> iterate(World& world,
                                                               real_function_3d& V,
                                                               real_function_3d& psi,
                                                               real_function_3d& phi,
                                                               double energy)
{
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*energy), 0.001, 1e-6);
    
    real_function_3d rhs = -2*(V*psi - energy*(psi-phi));
    rhs.truncate();
    real_function_3d phi_new = apply(op,rhs);

    double s = 1.0 + energy/(2*c*c);
    rhs = -2*(pVp(world,V,phi) + energy*(psi-s*phi));
    rhs.truncate();
    real_function_3d psimphi_new = apply(op,rhs);
    real_function_3d psi_new = s*phi_new + psimphi_new;
    psi_new.truncate();
    phi_new.truncate();

    double rnorm = (psi_new-psi).norm2();
    double norm = inner(psi_new,psi_new);
    psi_new *= 1.0/sqrt(norm);
    phi_new *= 1.0/sqrt(norm);
    
    return {psi_new, phi_new, rnorm};
}

double exact_energy(double Z) {
    if (Z<10) {
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

void run(World& world) {
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    real_function_3d Vnuc = real_factory_3d(world).f(V);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    real_function_3d phi  = real_factory_3d(world).f(guess);
    print("Initial should norm be 1:", psi.norm2());
    print("Exact Dirac-Coulomb energy:", exact_energy(Z));

    double energy = -0.5*Z*Z; // NR energy guess
    for (int iter=0; iter<20; iter++) {
        auto [psi_new, phi_new, rnorm] = iterate(world, Vnuc, psi, phi, energy);
        psi = psi_new;
        phi = phi_new;
        auto psisize = psi.size();
        auto phisize = phi.size();
        auto phinorm = phi.norm2();
        if (world.rank() == 0) {
            print("\niteration ", iter, " rnorm ", rnorm, " psi size ", psisize, " phi size ", phisize, " phi norm ", phinorm);
        }

        auto [total_energy, kinetic_energy, potential_energy, overlap] = compute_energy(world, Vnuc, psi, phi);
        
        if (world.rank() == 0) {
            print("            Kinetic energy  ", kinetic_energy);
            print(" Nuclear attraction energy ", potential_energy);
            print("              Total energy ", total_energy, total_energy - energy);
        }

        energy = total_energy;

        if (iter>2 and rnorm<1e-4) break;
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
