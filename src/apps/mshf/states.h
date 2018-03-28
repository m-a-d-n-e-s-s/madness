#ifndef STATES_H
#define STATES_H

#include "input.h"

// Harmonic Oscillator wavefunction
struct HO : FunctionFunctorInterface<double_complex,3>
{   
    const int nx, ny, nz;
    const double d;
    HO(int nx, int ny, int nz, double d) : nx(nx), ny(ny), nz(nz), d(d) {}
    double_complex operator()(const coordT& r) const {
        double rsq = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        double psi = exp( (-1.0/(2.0 * d * d)) * (rsq * rsq) ); 
        if (nx) for(int inx=0; inx<nx; inx++) {psi *= (1.0/d) * r[0];}
        if (ny) for(int iny=0; iny<ny; iny++) {psi *= (1.0/d) * r[1];}
        if (nz) for(int inz=0; inz<nz; inz++) {psi *= (1.0/d) * r[2];}
        return double_complex(psi,0.0);
     }
};


// Modified Harmonic Oscillator wavefunction
// Similar to HO but with different gaussian widths in x,y,z direction 
struct HOm : FunctionFunctorInterface<double_complex,3>
{   
    const int nx, ny, nz;
    const double d; 
    HOm(int nx, int ny, int nz, double d) : nx(nx), ny(ny), nz(nz), d(d) {}
    double_complex operator()(const coordT& r) const {
    const double dx = d;    // width in x
    const double dy = dx;   // width in y
    const double dz = dx;   // width in z
    double rsq = std::sqrt( (r[0] * r[0]/(dx * dx)) + (r[1] * r[1]/(dy * dy)) + (r[2] * r[2]/(dz * dz)) );
    double psi = std::exp((-1.0/2.0) * (rsq * rsq));
    
    if (nx) for(int inx=0; inx<nx; inx++) {psi *= r[0];}
    if (ny) for(int iny=0; iny<ny; iny++) {psi *= r[1];}
    if (nz) for(int inz=0; inz<nz; inz++) {psi *= r[2];}
    return double_complex(psi, 0.0);
    }
};


// Molecular Dynamics input: Transforms nucleon coordinates from a simulation space from (0 ... 2L) 
// to a simulation space from (-L...L). Creates 27 nucleon mirror images for periodic boundary conditions 
// and folds Gaussian around the nucleon corrdinates
struct MD : FunctionFunctorInterface<double_complex,3>
{
    const double nx, ny, nz;
    const double L;
    MD(double nx, double ny, double nz, double L) : nx(nx), ny(ny), nz(nz), L(L) {}
    double_complex operator()(const coordT& r) const {
        const double dmd = 3.0; // width of Gaussian
        double psi = 0.0;

        for (int k=-1; k < 2; k++) {
            for (int j=-1; j < 2; j++) {
                for (int i=-1; i < 2; i++) {
                    const double xs = (r[0] - (nx + (2.0 * L * i - 0.0)))/dmd;
                    const double ys = (r[1] - (ny + (2.0 * L * j - 0.0)))/dmd;
                    const double zs = (r[2] - (nz + (2.0 * L * k - 0.0)))/dmd;
                    psi += exp((-1.0/2.0) * (xs*xs + ys*ys + zs*zs));
                }
            }
        }
        return double_complex(psi, 0.0);
     }
};


// Plane wave wavefunction
struct Fermi : FunctionFunctorInterface<double_complex,3>
{
    //const double_complex I = double_complex(0.0, 1.0);
    const int nx, ny, nz, px, py, pz;
    const double L;
    Fermi(int nx, int ny, int nz, int px, int py, int pz, double L) : nx(nx), ny(ny), nz(nz), px(px), py(py), pz(pz), L(L) {}
    double_complex operator()(const coordT& r) const {
        double kx = px * nx * 2.0 * M_PI / (2.0 * L);
        double_complex psi  = exp(I * kx * r[0]);
        double ky = py * ny * 2.0 * M_PI / (2.0 * L);
                       psi *= exp(I * ky * r[1]);
        double kz = pz * nz * 2.0 * M_PI / (2.0 * L);
                       psi *= exp(I * kz * r[2]);
        return psi;
    }
};


// Make vector of gaussians around nucleon coordinates
extern void make_MD(World& world, comp_vecfuncT& psi_n, 
                                  comp_vecfuncT& psi_p,
                                  const double A, 
                                  const double Z, 
                                  const double L, 
                                  const double prec);

// Make vector of modified HO wavefunctions
extern void make_HO(World& world, comp_vecfuncT& u, const double A);


// Make vector of plane waves wavefunctions
extern void make_Fermi(World& world, comp_vecfuncT& u, 
                                     const double A, 
                                     const double Z, 
                                     const double L, 
                                     const double prec);


// Normalizes two vectors of complex functions
extern void normalize_2v(World& world, comp_vecfuncT& psi_n, comp_vecfuncT& psi_p);

// Normalizes one vector of complex functions with spin up and down
extern void normalize_ud(World& world, comp_vecfuncT& psi_qu, comp_vecfuncT& psi_qd);

// Truncates two vectors of complex functions
extern void truncate2(World& world, comp_vecfuncT& psi_q1,
                                    comp_vecfuncT& psi_q2,
                                    const double prec);

// Truncates one vector of complex functions
extern void truncate1(World& world, comp_vecfuncT& psi_q, const double prec);

// Splits one vector of complex functions into two for spin up and 
// spin down
extern void spin_split(World& world, const comp_vecfuncT& psi_q,
                                           comp_vecfuncT& psi_qu,
                                           comp_vecfuncT& psi_qd);

#endif
