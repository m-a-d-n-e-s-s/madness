#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <linalg/solvers.h>
#include <constants.h>
#include <vector>

using namespace madness;
using namespace std;

const int k = 6; // wavelet order
const double thresh = 1e-4; // truncation threshold
const double L = 5; // box is [-L,L]
const double sigma = 0.1; // Surface width

const double epsilon_0 = 1.0; // Interior dielectric
const double epsilon_1 =10.0; // Exterior dielectric
const double R = 2.0; // Radius of cavity

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \

// Distance between two points in 3D
double distance(const coord_3d& a, const coord_3d& b) {
    double x = a[0] - b[0];
    double y = a[1] - b[1];
    double z = a[2] - b[2];
    return sqrt(x*x + y*y + z*z);
}

// Basic functionality for the mask
class MolecularMaskBase {
protected:
    const double sigma;
    const vector<double> atomic_radii;
    const vector<coord_3d> atomic_coords;
    const int natom;

    // signed distance function for point r relative to sphere radius R at center
    double sdf(const coord_3d& r, const coord_3d& center, double R) const {
        return distance(r,center) - R;
    }

    // gradient of the signed distance function
    coord_3d grad_sdf(const coord_3d& r, const coord_3d& center) const {
        return (r - center)*(1.0/distance(r,center));
    }

    // Mask or characteristic function (argument s is the signed distance)
    double mask(double s) const {
        if (s > 6.0) return 0.0;
        else if (s < -6.0) return 1.0;
        else return 0.5*erfc(s);
    }

    // Complement of the mask or characteristic function (argument s is the signed distance)
    double cmask(double s) const {
        return mask(-s);
    }

    // Derivative of the mask w.r.t. s
    double dmask(double s) const {
        const double fac = 1.0/sqrt(constants::pi);
        if (fabs(s) > 6.0) return 0.0;
        return -exp(-s*s)*fac;
    }

    // Mask or characteristic function for atom i
    double atomic_mask(const coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return mask(s/sigma);
    }

    // Complement of the mask or characteristic function for atom i
    // (we use this directly to avoid numerical cancellation)
    double atomic_cmask(const coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return cmask(s/sigma);
    }

    // Gradient of the atomic mask
    coord_3d grad_atomic_mask(const coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return grad_sdf(r,atomic_coords[i])*(dmask(s/sigma)/sigma);
    }

    // Gradient of the molecular mask
    coord_3d gradient(const coord_3d& r) const {
        // precompute the atomic masks
        vector<double> m(natom);
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            m[i] = atomic_cmask(r,i);
            value *= m[i];
        }
        
        // return 0.0 if not in the surface
        if (value<1e-12 || (1.0-value)<1e-12) return 0.0;
        
        coord_3d grad(0.0);
        for (int i=0; i<natom; i++) {
            if (m[i] > 1e-12) 
                grad += grad_atomic_mask(r,i)*(value/m[i]);
        }
        return grad;
    }

public:
    MolecularMaskBase(double sigma, 
                      const vector<double> atomic_radii,
                      const vector<coord_3d> atomic_coords) 
        : sigma(sigma)
        , atomic_radii(atomic_radii)
        , atomic_coords(atomic_coords)
        , natom(atomic_coords.size())
    {
        MADNESS_ASSERT(atomic_radii.size() == atomic_coords.size());
    }
};

// This functor is 1 inside the molecule, 1/2 on the surface, and zero
// exterior to the molecule.
class MolecularVolumeMask : private MolecularMaskBase
                          , public FunctionFunctorInterface<double,3> {
public:
    MolecularVolumeMask(double sigma, 
                        const vector<double> atomic_radii,
                        const vector<coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const coord_3d& r) const {
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            value *= atomic_cmask(r,i);
        }
        return 1.0 - value;
    }
};

// This functor is a shell that limits to a delta function in the
// molecular surface and integrates to the molecular surface area.
class MolecularSurface : private MolecularMaskBase
                       , public FunctionFunctorInterface<double,3> {
public:
    MolecularSurface(double sigma, 
                     const vector<double> atomic_radii,
                     const vector<coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const coord_3d& r) const {
        coord_3d grad = gradient(r);
        return sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
    }
};

// Evaluates component i (0=x, 1=y, 2=z) of the gradient of the mask
class MolecularVolumeMaskGrad : private MolecularMaskBase
                              , public FunctionFunctorInterface<double,3> {
    const int i;
public:
    MolecularVolumeMaskGrad(double sigma, 
                     const vector<double> atomic_radii,
                     const vector<coord_3d> atomic_coords,
                     int i) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
        , i(i)
    {}

    virtual double operator()(const coord_3d& r) const {
        coord_3d grad = gradient(r);
        return grad[i];
    }
};

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

    real_tensor Q;
    vector_real_function_3d uvec, rvec;
    for (int iter=0; iter<20; iter++) {
        uvec.push_back(u);
        rvec.push_back(u - op(charge + (rfourpi)*rdielectric*(di_gradx*Dx(u) + di_grady*Dy(u) + di_gradz*Dz(u))).truncate());
  
        real_tensor newQ(iter+1,iter+1);
        if (iter>0) newQ(Slice(0,-2),Slice(0,-2)) = Q;
        Q = newQ;

        for (int jter=0; jter<=iter; jter++) {
            Q(jter,iter) = inner(uvec[jter],rvec[iter]);
            if (iter != jter) Q(iter,jter) = inner(uvec[iter],rvec[jter]);
        }
        
        real_tensor c = KAIN(Q);
        //print(Q);
        //print("KAIN", c);

        u = real_factory_3d(world); 
        for (int i=0; i<=iter; i++) u += c[i]*(uvec[i] - rvec[i]);

        real_function_3d u_prev = uvec[iter];

        double change = (u-u_prev).norm2();
        double err = (u-exact).norm2();
        print("iteration", iter, change, err, exact(coord_3d(3.0)), u(coord_3d(3.0)));

        if (change > 0.3*unorm) u = 0.5*u + 0.5*u_prev;

        if (change < 10.0*thresh) break;
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
