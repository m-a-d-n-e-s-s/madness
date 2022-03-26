#include <madness/mra/mra.h>
#include <madness/constants.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

#include <cmath>

using namespace madness;

static const double L = 60.0;
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

inline double mask1(double x) {
    x = (x*x*(3.-2.*x));
    return x;
}

static double mask3(const coord_3d& ruser) {
    coord_3d rsim;
    user_to_sim(ruser, rsim);
    double x= rsim[0], y=rsim[1], z=rsim[2];
    double lo = 0.0625, hi = 1.0-lo, result = 1.0;
    double rlo = 1.0/lo;
    
    if (x<lo)
        result *= mask1(x*rlo);
    else if (x>hi)
        result *= mask1((1.0-x)*rlo);
    if (y<lo)
        result *= mask1(y*rlo);
    else if (y>hi)
        result *= mask1((1.0-y)*rlo);
    if (z<lo)
        result *= mask1(z*rlo);
    else if (z>hi)
        result *= mask1((1.0-z)*rlo);
    
    return result;
}

static double psi1s(const coord_3d& r) {
    const double fac = 1.0/std::sqrt(constants::pi);
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x+y*y+z*z);
    return fac*std::exp(-R);
}

static double psi2s(const coord_3d& r) {
    const double fac = 1.0/(4.0*std::sqrt(2*constants::pi));
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x+y*y+z*z);
    return fac*(2.0 - R)*std::exp(-0.5*R);
}

static double psi2p(const coord_3d& r) {
    const double fac = 1.0/(4.0*std::sqrt(2*constants::pi));
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x+y*y+z*z);
    return fac*z*std::exp(-0.5*R);
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/std::sqrt(x*x+y*y+z*z);
}

static double zdipole(const coord_3d& r) {
    return r[2];
}

double iterate(World& world, NonlinearSolver& solver, const real_function_3d& V, const real_function_3d& mask, const real_function_3d& PSI1S, real_function_3d& psi, double& eps) {
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d tmp = apply(op,Vpsi);
    tmp -= PSI1S*inner(PSI1S,tmp); // Constrain new solution to be orthogonal to 1s
    tmp.truncate();
    tmp = tmp*mask;
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    //psinew=psi+delta...delta defined in KAIN paper.  Uses residual and current
    real_function_3d psinew = solver.update(psi, r, 1e-8, 1e4);
    double step = (psinew - psi).norm2();
    if (step > 0.1) {
        print("  step restriction", step);
        psinew = psi + (0.1/step)*(psinew-psi);
    }
    psi = psinew;
    psi.truncate();
    psi.scale(1.0/psi.norm2());
    //psi = tmp.scale(1.0/norm);
    eps = eps_new;
    return rnorm;
}

void solve(World& world,
           const real_function_3d& V,
           const real_function_3d& mask,
           const real_function_3d& PSI1S,
           const real_function_3d& PSI2S,
           const real_function_3d& PSI2P,
           const real_function_3d& Z,
           real_function_3d& psi,
           double& eps,
           double& energy)
{
    NonlinearSolver solver(10);
    for (int iter=0; iter<1000; iter++) {
        double err = iterate(world, solver, V, mask, PSI1S, psi, eps);
        double proj1s = inner(psi,PSI1S);
        double proj2s = inner(psi,PSI2S);
        double proj2p = inner(psi,PSI2P);
        double dipole = inner(psi*psi,Z);
        print("<1s|psi>", proj1s, "<2s|psi>", proj2s, "<2p|psi>", proj2p, "<psi|z|psi>", dipole);
        if (iter>=3 && err<1e-4) break;
    }
    
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }
    kinetic_energy *= 0.5;
    
    double potential_energy = inner(V,psi*psi);
    double total_energy = kinetic_energy + potential_energy;
    
    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print("          Potential energy ", potential_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", potential_energy / kinetic_energy);
    }

    energy = total_energy;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(8);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    real_function_3d Vnuc = real_factory_3d(world).f(V).truncate_mode(1);
    real_function_3d PSI1S = real_factory_3d(world).f(psi1s);
    real_function_3d psi  = real_factory_3d(world).f(psi2s);
    real_function_3d mask  = real_factory_3d(world).f(mask3);
    real_function_3d PSI2S = copy(psi);
    psi = psi*mask;
    real_function_3d PSI2P = real_factory_3d(world).f(psi2p);
    real_function_3d Z = real_factory_3d(world).f(zdipole);
    psi.scale(1.0/psi.norm2());

    Z = Z*mask;

    double eps = -1.0/8.0;
    double energy;

    double field_increment = 0.0001;

    for (int f=0; f<20; f++) {
        print("Solving with field", field_increment*f);
        print("");
        solve(world, Vnuc, mask, PSI1S, PSI2S, PSI2P, Z, psi, eps, energy);
        Vnuc += field_increment * Z;
    }
    
    world.gop.fence();

    finalize();
    return 0;
}
