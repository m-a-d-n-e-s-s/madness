#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <linalg/solvers.h>
#include <examples/molecularmask.h>
#include <examples/nonlinsol.h>
#include <constants.h>
#include <vector>

using namespace madness;
using namespace std;

const int k = 6; // wavelet order
const double thresh = 1e-4; // truncation threshold
const double L = 5; // box is [-L,L]
const double sigma = 0.2; // Surface width

const double epsilon_0 = 1.0; // Interior dielectric
const double epsilon_1 =10.0; // Exterior dielectric
const double R = 2.0; // Radius of cavity

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \


template <typename T, int NDIM>
struct Reciprocal {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = 1.0/(*_p0));
    }
    template <typename Archive> void serialize(Archive& ar) {}
};

double charge_function(const coord_3d& r) {
    const double expnt = 100.0;
    const double coeff = pow(1.0/constants::pi*expnt,0.5*3);
    return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double exact_function(const coord_3d& x) {
    const double expnt = 100.0;
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    if (r > R) {
        return 1.0/(epsilon_1*r);
    }
    else {
        return erf(sqrt(expnt)*r)/(epsilon_0*r)
            + (1.0/epsilon_1 - 1.0/epsilon_0)/R;
    }
}

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    coord_3d lo(0.0), hi(0.0); // Range for line plotting
    lo[0] =-5.0;
    hi[0] = 5.0;

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_truncate_on_project(true);
    //FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_bc(BC_FREE);

    // The Coulomb operator (this is just 1/r ... whereas the notes are -1/4pir)
    real_convolution_3d op = CoulombOperator(world, sigma*0.001, thresh*0.1);

    // Derivative operators
    real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);

    // We will have one sphere of radius R centered at the origin
    vector<double> atomic_radii(1,R);
    vector<coord_3d> atomic_coords(1,coord_3d(0.0));
    print("k     ", k);
    print("thresh", thresh);
    print("L     ", L);
    print("sigma ", sigma);
    print("eps0  ", epsilon_0, "  eps1  ", epsilon_1);
    print("radii ", atomic_radii);
    print("coords", atomic_coords);

    // Functors for mask related quantities
    real_functor_3d volume_functor(new MolecularVolumeMask(sigma, atomic_radii, atomic_coords));
    real_functor_3d gradx_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 0));
    real_functor_3d grady_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 1));
    real_functor_3d gradz_functor(new MolecularVolumeMaskGrad(sigma, atomic_radii, atomic_coords, 2));
    real_functor_3d surface_functor(new MolecularSurface(sigma, atomic_radii, atomic_coords));

    // Make the actual functions
    TIME("make volume ", real_function_3d volume = real_factory_3d(world).functor(volume_functor));
    TIME("make gradx  ", real_function_3d gradx = real_factory_3d(world).functor(gradx_functor));
    TIME("make grady  ", real_function_3d grady = real_factory_3d(world).functor(grady_functor));
    TIME("make gradz  ", real_function_3d gradz = real_factory_3d(world).functor(gradz_functor));
    TIME("make surface", real_function_3d surface = real_factory_3d(world).functor(surface_functor));
    TIME("make charge ", real_function_3d charge = real_factory_3d(world).f(charge_function));
    TIME("make exact  ", real_function_3d exact = real_factory_3d(world).f(exact_function));
    
    // Reciprocal of the dielectric function
    real_function_3d rdielectric = epsilon_0*volume + epsilon_1*(1.0-volume);
    rdielectric.unaryop(Reciprocal<double,3>());

    // Gradient of the dielectric function
    real_function_3d di_gradx = (epsilon_0-epsilon_1)*gradx;
    real_function_3d di_grady = (epsilon_0-epsilon_1)*grady;
    real_function_3d di_gradz = (epsilon_0-epsilon_1)*gradz;

    // Print some values for sanity checking
    print("the volume is", volume.trace());
    print("the area   is", surface.trace());
    print("the charge is", charge.trace());

    // Free up stuff we are not using any more to save memory
    volume.clear();
    surface.clear();
    gradx.clear();
    grady.clear();
    gradz.clear();

    const double rfourpi = 1.0/(4.0*constants::pi);
    charge = (1.0/epsilon_0)*charge;

    // Initial guess is constant dielectric
    real_function_3d u = op(charge).truncate();
    double unorm = u.norm2();

    const bool USE_SOLVER = true;

    if (USE_SOLVER) {
        // This section employs a non-linear equation solver from solvers.h
        // and nonlinsol.h as described in
        // http://onlinelibrary.wiley.com/doi/10.1002/jcc.10108/abstract
        NonlinearSolver solver;
        for (int iter=0; iter<20; iter++) {
            real_function_3d surface_charge = 
                rfourpi*rdielectric*(di_gradx*Dx(u) + di_grady*Dy(u) + di_gradz*Dz(u));
            real_function_3d r = (u - op(charge + surface_charge)).truncate();

            real_function_3d unew = solver.update(u, r);

            double change = (unew-u).norm2();
            double err = (u-exact).norm2();
            print("iter", iter, "change", change, "err", err, 
                  "exact(3.0)", exact(coord_3d(3.0)), "soln(3.0)", u(coord_3d(3.0)),
                  "surface charge", surface_charge.trace());

            if (change > 0.3*unorm) 
                u = 0.5*unew + 0.5*u;
            else 
                u = unew;

            if (change < 10.0*thresh) break;
        }
    }
    else {
        // This section employs a simple iteration with damping (step restriction)
        for (int iter=0; iter<20; iter++) {
            real_function_3d u_prev = u;
            real_function_3d surface_charge = 
                rfourpi*rdielectric*(di_gradx*Dx(u) + di_grady*Dy(u) + di_gradz*Dz(u));
            u = op(charge + surface_charge).truncate();
            
            double change = (u-u_prev).norm2();
            double err = (u-exact).norm2();
            print("iteration", iter, change, err, exact(coord_3d(3.0)), u(coord_3d(3.0)), surface_charge.trace());
            if (change > 0.3*unorm) u = 0.5*u + 0.5*u_prev;
            if (change < 10.0*thresh) break;
        }
    }


    plot_line("testpot.dat", 301, lo, hi, u, exact);
    real_tensor cell(3,2);
    cell(_,0) = -4.0;
    cell(_,1) =  4.0;

    plotdx(u, "testpot.dx", cell);
    plotdx(exact, "exact.dx", cell);
    plotdx(u-exact, "err.dx", cell);

    finalize();
    return 0;
}
