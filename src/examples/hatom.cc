#include <fstream>
#include <iostream>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;

static const double L = 32.0;   // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-7; // precision
static const double F=0.0; // electric field

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    // this is the exact solution
    //return exp(-sqrt(x*x+y*y+z*z))/sqrt(constants::pi);
    //
    // this is too diffuse
    return exp(-sqrt(x*x+y*y+z*z)/2.0)/sqrt(8.0*constants::pi);
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/(sqrt(x*x+y*y+z*z)) + F*z;
}

inline double mask1(double x) {
    /* Iterated first beta function to switch smoothly
           from 0->1 in [0,1].  n iterations produce 2*n-1
           zero derivatives at the end points. Order of polyn
           is 3^n.

           Currently use one iteration so that first deriv.
           is zero at interior boundary and is exactly representable
           by low order multiwavelet without refinement */

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
        result *= mask1(x * rlo);
    else if (x > hi)
        result *= mask1((1.0 - x) * rlo);
    if (y < lo)
        result *= mask1(y * rlo);
    else if (y > hi)
        result *= mask1((1.0 - y) * rlo);
    if (z < lo)
        result *= mask1(z * rlo);
    else if (z > hi)
        result *= mask1((1.0 - z) * rlo);

    return result;
}

void iterate(World& world, real_function_3d& V, real_function_3d& mask, real_function_3d& psi, double& eps) {
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d tmp = apply(op,Vpsi).truncate();
    tmp = tmp*mask;
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new;
    if (rnorm > 0.2) {
        r *= 0.2/rnorm;
        print("step restriction");
        eps_new = eps;
    }
    else {
        // Only update energy once step restriction is lifted since this is only locally convergent
        eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    }
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi += r;
    psi.scale(1.0/psi.norm2());
    eps = eps_new;
}

std::tuple<double,double,double> compute_energy(World& world, real_function_3d& Vnuc, real_function_3d& psi) {
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += 0.5*inner(dpsi,dpsi);
    }
    real_function_3d rho = square(psi).truncate();
    double nuclear_attraction_energy = inner(Vnuc*psi,psi);
    double total_energy = kinetic_energy + nuclear_attraction_energy;
    return {kinetic_energy, nuclear_attraction_energy, total_energy};
}

void run(World& world) {
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

    real_function_3d Vnuc = real_factory_3d(world).f(V).truncate_mode(1);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    real_function_3d mask = real_factory_3d(world).f(mask3);
    psi = psi*mask;
    print("initial", psi.norm2());
    psi.scale(1.0/psi.norm2());

    double eps = -0.5;
    for (int iter=0; iter<15; iter++) {
        real_function_3d rho = square(psi).truncate();
        iterate(world, Vnuc, mask, psi, eps);
    }

    // Manually tabluate the orbital along a line along the z axis ... probably easier to use the lineplot routine
    coord_3d r(0.0);
    psi.reconstruct();
    for (int i=0; i<201; i++) {
        r[2] = -L/2 + L*i/200.0;
        print(r[2], psi(r));
    }
  
    auto [kinetic_energy, nuclear_attraction_energy, total_energy] = compute_energy(world, Vnuc, psi);
    if (world.rank() == 0) {
        print("            Electric Field ", F);
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("              Total energy ", total_energy);
    }
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world,argc,argv);
    std::cout.precision(6);

    run(world);
    
    world.gop.fence();
    finalize();
    return 0;
}
