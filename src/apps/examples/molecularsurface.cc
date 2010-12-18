#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>
#include <mra/funcplot.h>
#include <constants.h>
#include <vector>

using namespace madness;
using namespace std;

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \

// Area of two intersecting spheres separated by d
double area_two_spheres(double r1, double r2, double d) {
    if (d > r1+r2) d = r1 + r2;
    return constants::pi*(2*d*r1*r1+2*d*r2*r2+r1*r1*r1+r1*d*d-r1*r2*r2+r2*d*d-r2*r1*r1+r2*r2*r2)/d;
}

// Volume of two intersecting spheres separated by d
double volume_two_spheres(double r1, double r2, double d) {
    if (d > r1+r2) return 0.0;
    double overlap = constants::pi*(r1+r2-d)*(r1+r2-d)*(d*d + 2*d*r1 + 2*d*r2 - 3*r1*r1 - 3*r2*r2 + 6*r1*r2)/(12*d);
    return 4.0*constants::pi*(r1*r1*r1 + r2*r2*r2)/3.0 - overlap;
}

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

    // Gradient of the moleculalr mask
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
    MolecularSurface(double sigma, 
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

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    const int k = 6; // wavelet order
    const double thresh = 1e-4; // truncation threshold
    const double L = 5; // box is [-L,L]
    const int natom = 2; // number of atoms
    const double sigma = 0.1; // Surface width

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(2);

    // Set up atomic coordinates and radii
    vector<double> atomic_radii(natom);
    vector<coord_3d> atomic_coords(natom);
    for (int i=0; i<natom; i++) {
        atomic_radii[i] = 1.0 + i*0.5;
        atomic_coords[i][0] = 0.0;
        atomic_coords[i][1] = 0.0;
        atomic_coords[i][2] = i*1.5;
        print("atom",i,atomic_coords[i],atomic_radii[i]);
    }

    real_functor_3d volume_functor(new MolecularVolumeMask(sigma, atomic_radii, atomic_coords));
    real_functor_3d surface_functor(new MolecularSurface(sigma, atomic_radii, atomic_coords));

    TIME("make volume ", real_function_3d volume = real_factory_3d(world).functor(volume_functor));
    TIME("make surface", real_function_3d surface = real_factory_3d(world).functor(surface_functor).truncate_on_project());

    print("the volume is", volume.trace());
    print("the area   is", surface.trace());

    if (natom == 2) {
        double d = distance(atomic_coords[0],atomic_coords[1]);
        double r1 = atomic_radii[0];
        double r2 = atomic_radii[1];
        print("d",d,"r1",r1,"r2",r2,"r1+r2-d",r1+r2-d);
        print("analytic volume intersecting spheres", volume_two_spheres(r1, r2, d));
        print("analytic area   intersecting spheres", area_two_spheres(r1, r2, d));
    }

    std::vector<long> npt(3,401);
    Tensor<double> cell(3,2);
    cell(_,0) = -5;
    cell(_,1) =  5;

    TIME("plot surface",plotdx(surface, "surface.dx"));
    TIME("plot volume ",plotdx(volume, "volume.dx", cell, npt));
    TIME("plot povray ",plotpovray(volume, "volume.df3", cell, npt));

    finalize();
    return 0;
}
