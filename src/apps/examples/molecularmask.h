#ifndef MOLECULAR_MASK_H
#define MOLECULAR_MASK_H

#include <mra/mra.h>
#include <constants.h>
#include <cmath>
#include <vector>

// Distance between two points in 3D
inline double distance(const madness::coord_3d& a, const madness::coord_3d& b) {
    double x = a[0] - b[0];
    double y = a[1] - b[1];
    double z = a[2] - b[2];
    return sqrt(x*x + y*y + z*z);
}

// Basic functionality for the mask
class MolecularMaskBase {
protected:
    const double sigma;
    const std::vector<double> atomic_radii;
    const std::vector<madness::coord_3d> atomic_coords;
    const int natom;

    // signed distance function for point r relative to sphere radius R at center
    double sdf(const madness::coord_3d& r, const madness::coord_3d& center, double R) const {
        return distance(r,center) - R;
    }

    // gradient of the signed distance function
    madness::coord_3d grad_sdf(const madness::coord_3d& r, const madness::coord_3d& center) const {
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
        const double fac = 1.0/sqrt(madness::constants::pi);
        if (fabs(s) > 6.0) return 0.0;
        return -exp(-s*s)*fac;
    }

    // Mask or characteristic function for atom i
    double atomic_mask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return mask(s/sigma);
    }

    // Complement of the mask or characteristic function for atom i
    // (we use this directly to avoid numerical cancellation)
    double atomic_cmask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return cmask(s/sigma);
    }

    // Gradient of the atomic mask
    madness::coord_3d grad_atomic_mask(const madness::coord_3d& r, unsigned int i) const {
        double s = sdf(r, atomic_coords[i], atomic_radii[i]);
        return grad_sdf(r,atomic_coords[i])*(dmask(s/sigma)/sigma);
    }

    // Gradient of the molecular mask
    madness::coord_3d gradient(const madness::coord_3d& r) const {
        // precompute the atomic masks
        std::vector<double> m(natom);
        double value = 1.0;
        for (int i=0; i<natom; i++) {
            m[i] = atomic_cmask(r,i);
            value *= m[i];
        }
        
        // return 0.0 if not in the surface
        if (value<1e-12 || (1.0-value)<1e-12) return 0.0;
        
        madness::coord_3d grad(0.0);
        for (int i=0; i<natom; i++) {
            if (m[i] > 1e-12) 
                grad += grad_atomic_mask(r,i)*(value/m[i]);
        }
        return grad;
    }

public:
    MolecularMaskBase(double sigma, 
                      const std::vector<double> atomic_radii,
                      const std::vector<madness::coord_3d> atomic_coords) 
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
                          , public madness::FunctionFunctorInterface<double,3> {
public:
    MolecularVolumeMask(double sigma, 
                        const std::vector<double> atomic_radii,
                        const std::vector<madness::coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
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
                       , public madness::FunctionFunctorInterface<double,3> {
public:
    MolecularSurface(double sigma, 
                     const std::vector<double> atomic_radii,
                     const std::vector<madness::coord_3d> atomic_coords) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
        madness::coord_3d grad = gradient(r);
        return sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
    }
};

// Evaluates component i (0=x, 1=y, 2=z) of the gradient of the mask
class MolecularVolumeMaskGrad : private MolecularMaskBase
                              , public madness::FunctionFunctorInterface<double,3> {
    const int i;
public:
    MolecularVolumeMaskGrad(double sigma, 
                     const std::vector<double> atomic_radii,
                     const std::vector<madness::coord_3d> atomic_coords,
                     int i) 
        : MolecularMaskBase(sigma, atomic_radii, atomic_coords)
        , i(i)
    {}

    virtual double operator()(const madness::coord_3d& r) const {
        madness::coord_3d grad = gradient(r);
        return grad[i];
    }
};

#endif
