// We will adopt the intel MKL interface as the standard
// for vector math routines.  

// Must provide a compatibility interface to ACML and also
// to no underlying math library

#include <complex>
#include <cmath>

typedef std::complex<double> double_complex;

#include <madness_config.h>
#ifdef HAVE_MKL
#include <mkl.h>

#elif defined(HAVE_ACML)
#include <acml_mv.h>

void vdSinCos(int n, const double* x, double* sinx, double* cosx) {
    vrda_sincos(n, const_cast<double*>(x), sinx, cosx);
}

void vdExp(int n, const double *x, double *y) {
    vrda_exp(n, const_cast<double*>(x), y);
}

void vzExp(int n, const double_complex* x, double_complex* y) {
    if (n <= 0) return;
    double* a = new double[n];
    double* b = new double[n];
    double* expa = new double[n];
    double* sinb = new double[n];
    double* cosb = new double[n];

    for (int i=0; i<n; i++) {
        a[i] = x[i].real();
        b[i] = x[i].imag();
    }
    vdExp(n, a, expa);
    vdSinCos(n, b, sinb, cosb);
    for (int i=0; i<n; i++) {
        y[i] = double_complex(a[i]*cosb[i],a[i]*sinb[i]);
    }
    delete cosb;
    delete sinb;
    delete expa;
    delete b;
    delete a;
}

#else

void vzExp(int n, const double_complex* x, double_complex* y) {
    for (int i=0; i<n; i++) y[i] = exp(x[i]);
}

#endif
