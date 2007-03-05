#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <madness_config.h>
#include <world/world.h>

namespace madness {
    extern void load_quadrature(World& world);
    extern void legendre_polynomials(double x, long order, double *p);
    extern void legendre_scaling_functions(double x, long k, double *p);
    extern bool gauss_legendre(int n, double xlo, double xhi, double *x, double *w);
    extern bool gauss_legendre_numeric(int n, double xlo, double xhi, double *x, double *w);
    extern bool gauss_legendre_test(bool print=false);
}

#endif

