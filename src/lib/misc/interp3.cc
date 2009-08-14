#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

template <typename T>
class CubicInterpolationTable {
    const double lo;            //< Interpolation is in range [lo,hi]
    const double hi;            //< Interpolation is in range [lo,hi]
    const double h;             //< Grid spacing
    const double rh;            //< 1/h
    const int npt;              //< No. of grid points
    std::vector<T> a;           //< (1+4)*npt vector of x and polynomial coefficients

    // Cubic interp thru 4 points ... not good for noisy data
    void cubic_fit(const double* x, const double* f, double* a) {
        double x0_2 = x[0] * x[0], x1_2 = x[1] * x[1], x2_2 = x[2] * x[2], x3_2 = x[3] * x[3];
        double x0_3 = x[0] * x[0] * x[0], x1_3 = x[1] * x[1] * x[1], x2_3 = x[2] * x[2] * x[2], x3_3 = x[3] * x[3] * x[3];
        
        a[0] = -(-x0_3 * x2_2 * x[3] * f[1] + x0_3 * x[2] * x3_2 * f[1] - x0_3 * f[3] * x[2] * x1_2 + x0_3 * x[3] * f[2] * x1_2 + x0_3 * f[3] * x2_2 * x[1] - x0_3 * x3_2 * f[2] * x[1] + x0_2 * x1_3 * f[3] * x[2] - x0_2 * x1_3 * f[2] * x[3] + x0_2 * x3_3 * f[2] * x[1] + x0_2 * f[1] * x2_3 * x[3] - x0_2 * f[1] * x3_3 * x[2] - x0_2 * f[3] * x2_3 * x[1] + x[0] * x3_2 * f[2] * x1_3 - x[0] * f[3] * x2_2 * x1_3 + x[0] * x1_2 * f[3] * x2_3 - x[0] * x1_2 * f[2] * x3_3 - x[0] * f[1] * x3_2 * x2_3 + x[0] * f[1] * x2_2 * x3_3 - f[0] * x2_3 * x1_2 * x[3] + f[0] * x2_2 * x1_3 * x[3] + f[0] * x3_2 * x2_3 * x[1] - f[0] * x3_3 * x2_2 * x[1] + f[0] * x3_3 * x1_2 * x[2] - f[0] * x3_2 * x1_3 * x[2]) / (-x2_2 * x[0] * x3_3 + x2_2 * x[0] * x1_3 - x0_2 * x[3] * x2_3 + x0_2 * x3_3 * x[2] - x0_2 * x[1] * x3_3 + x0_2 * x[1] * x2_3 + x0_2 * x1_3 * x[3] - x0_2 * x1_3 * x[2] + x3_2 * x[0] * x2_3 - x3_2 * x[0] * x1_3 + x[3] * x2_3 * x1_2 - x3_2 * x2_3 * x[1] + x3_3 * x2_2 * x[1] - x3_3 * x[2] * x1_2 + x[0] * x3_3 * x1_2 - x[0] * x2_3 * x1_2 - x0_3 * x3_2 * x[2] - x0_3 * x[3] * x1_2 + x0_3 * x3_2 * x[1] + x[2] * x3_2 * x1_3 - x2_2 * x[3] * x1_3 + x0_3 * x2_2 * x[3] - x0_3 * x2_2 * x[1] + x0_3 * x[2] * x1_2);
        a[1] = (-x2_3 * x1_2 * f[0] + x3_2 * x2_3 * f[0] + x2_3 * x0_2 * f[1] + x1_2 * f[3] * x2_3 - x2_3 * x0_2 * f[3] - f[1] * x3_2 * x2_3 - f[3] * x2_2 * x1_3 - x3_3 * x2_2 * f[0] + f[1] * x2_2 * x3_3 + x2_2 * x1_3 * f[0] - f[1] * x2_2 * x0_3 + f[3] * x2_2 * x0_3 - x1_3 * x3_2 * f[0] - x0_2 * x1_3 * f[2] - f[3] * x0_3 * x1_2 + f[1] * x3_2 * x0_3 + x1_2 * f[2] * x0_3 + x3_3 * f[0] * x1_2 - x3_2 * f[2] * x0_3 - f[1] * x3_3 * x0_2 + x0_2 * x3_3 * f[2] - x1_2 * f[2] * x3_3 + x3_2 * f[2] * x1_3 + x0_2 * x1_3 * f[3]) / (-x[3] + x[2]) / (-x2_2 * x0_2 * x[3] - x2_2 * x[1] * x3_2 + x2_2 * x1_2 * x[3] + x2_2 * x[0] * x3_2 - x2_2 * x[0] * x1_2 + x2_2 * x0_2 * x[1] + x[2] * x[0] * x1_3 + x[2] * x0_3 * x[3] - x[2] * x0_3 * x[1] - x[2] * x1_3 * x[3] - x[2] * x0_2 * x3_2 + x[2] * x1_2 * x3_2 - x[2] * x[3] * x[0] * x1_2 + x[2] * x[3] * x0_2 * x[1] + x0_3 * x1_2 - x0_2 * x1_3 + x[3] * x[0] * x1_3 - x[3] * x0_3 * x[1] - x3_2 * x[0] * x1_2 + x3_2 * x0_2 * x[1]);
        a[2] = -(-x1_3 * f[3] * x[2] + x1_3 * f[2] * x[3] + x1_3 * x[0] * f[3] + x1_3 * f[0] * x[2] - x1_3 * x[0] * f[2] - x1_3 * f[0] * x[3] + f[3] * x2_3 * x[1] - f[3] * x0_3 * x[1] - x[1] * x2_3 * f[0] + x[1] * f[2] * x0_3 + x3_3 * f[0] * x[1] - x3_3 * f[2] * x[1] + f[1] * x[3] * x0_3 - f[1] * x[0] * x3_3 - x3_3 * f[0] * x[2] + x3_3 * x[0] * f[2] - f[2] * x0_3 * x[3] + x2_3 * f[0] * x[3] + f[1] * x3_3 * x[2] - f[1] * x2_3 * x[3] + x2_3 * x[0] * f[1] - x2_3 * x[0] * f[3] + f[3] * x0_3 * x[2] - x[2] * f[1] * x0_3) / (x[3] * x2_2 - x3_2 * x[2] + x[1] * x3_2 - x[1] * x2_2 - x1_2 * x[3] + x1_2 * x[2]) / (-x[2] * x[1] * x[3] + x[1] * x[2] * x[0] - x0_2 * x[1] + x[1] * x[3] * x[0] + x[2] * x[3] * x[0] + x0_3 - x0_2 * x[2] - x0_2 * x[3]);
        a[3] = (x[0] * f[3] * x1_2 - x0_2 * x[3] * f[2] + x2_2 * x[0] * f[1] + x0_2 * f[3] * x[2] - x3_2 * f[2] * x[1] - f[0] * x3_2 * x[2] - f[3] * x[2] * x1_2 - x2_2 * x[0] * f[3] - f[0] * x2_2 * x[1] + f[3] * x2_2 * x[1] + x0_2 * f[1] * x[3] + x[2] * x3_2 * f[1] - x0_2 * f[1] * x[2] + f[0] * x[2] * x1_2 + x[3] * f[2] * x1_2 + f[0] * x3_2 * x[1] + x3_2 * x[0] * f[2] - x[0] * f[2] * x1_2 - f[0] * x[3] * x1_2 - x0_2 * x[1] * f[3] + x0_2 * x[1] * f[2] + f[0] * x2_2 * x[3] - x2_2 * x[3] * f[1] - x3_2 * x[0] * f[1]) / (-x2_2 * x[0] * x3_3 + x2_2 * x[0] * x1_3 - x0_2 * x[3] * x2_3 + x0_2 * x3_3 * x[2] - x0_2 * x[1] * x3_3 + x0_2 * x[1] * x2_3 + x0_2 * x1_3 * x[3] - x0_2 * x1_3 * x[2] + x3_2 * x[0] * x2_3 - x3_2 * x[0] * x1_3 + x[3] * x2_3 * x1_2 - x3_2 * x2_3 * x[1] + x3_3 * x2_2 * x[1] - x3_3 * x[2] * x1_2 + x[0] * x3_3 * x1_2 - x[0] * x2_3 * x1_2 - x0_3 * x3_2 * x[2] - x0_3 * x[3] * x1_2 + x0_3 * x3_2 * x[1] + x[2] * x3_2 * x1_3 - x2_2 * x[3] * x1_3 + x0_3 * x2_2 * x[3] - x0_3 * x2_2 * x[1] + x0_3 * x[2] * x1_2);
    }

public:
    template <typename functionT>
    CubicInterpolationTable(double lo, double hi, int npt, const functionT& f) 
        : lo(lo)
        , hi(hi)
        , h((hi-lo)/(npt-1))
        , rh(1.0/h)
        , npt(npt)
        , a(npt*5)
    {
        // Evaluate the function to be interpolated
        std::vector<T> p(npt), x(npt);
        for (int i=0; i<npt; i++) {
            x[i] = lo + i*h;
            p[i] = f(x[i]);
        }

        // Generate interior polynomial coeffs
        for (int i=1; i<=npt-3; i++) {
            double mid = (x[i] + x[i+1])*0.5;
            double y[4] = {x[i-1]-mid,x[i]-mid,x[i+1]-mid,x[i+2]-mid};
            a[i*5] = mid;
            cubic_fit(y, &p[i-1], &a[i*5+1]);
        }

        // Fixup end points
        for (int j=0; j<5; j++) {
            a[j] = a[5+j];
            a[5*npt-5+j] = a[5*npt-10+j] = a[5*npt-15+j];
        }
    }

    T operator()(double y) const {
        int i = (y-lo)*rh;
        if (i<0 || i>=npt) throw "Out of range point";
        i *= 5;
        y -= a[i];
        double yy = y*y;
        return (a[i+1] + y*a[i+2]) + yy*(a[i+3] + y*a[i+4]);
    }

    template <typename functionT>
    double err(const functionT& f) const {
        double maxabserr = 0.0;
        double h7 = h/7.0;
        for (int i=0; i<7*npt; i++) {
            double x = lo + h7*i;
            double fit = (*this)(x);
            double exact = f(x);
            maxabserr = max(fabs(fit-exact),maxabserr);
        }
        return maxabserr;
    }
};

double func(double x) {
    return sin(x);
}

int main() {
    cout.precision(12);

    // Uniform mesh for sin(x)
    CubicInterpolationTable<double> fit(-10.0,30.0,100000,func);
    cout << "maxerr " << fit.err(func) << endl;

    return 0;
}

