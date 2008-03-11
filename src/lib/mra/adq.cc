#include <mra/legendre.h>
#include <cmath>
#include <iostream>
using namespace std;

/// \file adq.cc
/// \brief Adaptive 1d quadrature using Gauss-Legendre on intervals

namespace madness {
    template <typename funcT>
    static double do_adq(double lo, double hi, const funcT& func, 
                  int n, const double* x, const double* w) {
        // x and w provide the Gauss-Legendre quadrature rule of order n on [0,1]
        double range = (hi-lo);
        double sum = 0.0;
        for (int i=0; i<n; i++) sum += func(lo + range*x[i])*w[i];
        return sum*range;
    }


    template <typename funcT>
    double adq1(double lo, double hi, const funcT& func, double thresh, 
                int n, const double* x, const double* w, int level) {
        static int MAX_LEVEL=10;
        if (level > MAX_LEVEL) throw "Adaptive quadrature: failed : runaway refinement?";
        double full = do_adq(lo, hi, func, n, x, w);
        double d = (hi-lo)/2;
        double half = do_adq(lo, lo+d, func, n, x, w) + do_adq(lo+d, hi, func, n, x, w);

        for (int i=0; i<level; i++) cout << "  ";
        cout << half << " " << full-half << endl;


        if (abs(full-half) > thresh) {
            return adq1(lo, lo+d, func, thresh*0.5, n, x, w, level+1) + 
                adq1(lo+d, hi, func, thresh*0.5, n, x, w, level+1);
        }
        else {
            return half;
        }
    }

    template <typename funcT>
    double adq(double lo, double hi, const funcT& func, double thresh) {
        const int n = 20;
        double x[n], y[n];
        gauss_legendre(n, 0.0, 1.0, x, y);

        return adq1(lo, hi, func, thresh, n, x, y, 0);
    }

    static double testfunc(double x) {
        // int(exp(-x^2)*cos(a*x),x=-inf..inf) = sqrt(pi)*exp(-a^2/4)
        return exp(-x*x)*cos(3*x);
    }

    bool testadq() {
        const double pi = 3.1415926535897932384;
        double test = adq(-6.0,6.0,testfunc,1e-14);
        double exact = sqrt(pi)*exp(-9.0/4.0);
        return std::fabs(test-exact)<1e-14;
    }
}
