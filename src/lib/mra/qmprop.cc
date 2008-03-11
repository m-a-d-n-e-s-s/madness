#include <cmath>
#include <complex>
#include <constants.h>
#include <iostream>
using namespace std;


/// Class to evaluate the filtered Schrodinger free-particle propagator in real space

/// Follows the corresponding Maple worksheet
class BandlimitedPropagator {
private:
    const double c;
    const double t;
    const double ctop;
    const double L;
    const double h;
    const int n;
    const double dc;
    
    std::complex<double> ff(double k) {
        if (k>2.54*c) return std::complex<double>(0.0,0.0);
        const std::complex<double> arg(0,-k*k*t*0.5);
        return std::exp(arg)/(1.0+std::pow(k/c, 30.0));
    }

public:
    BandlimitedPropagator(double c, double t)
        : c(c)
        , t(t)
        , ctop(2.15*c)
        , L(1.2*0.5435*(3.0*pow(c,5.0/3.0)*pow(t,0.75)+400.0)/c) // 1.2 --> 2.0 ??
        , h(3.14/ctop)
        , n(2.0*L/h+1)
        , dc(2*ctop/(n-1))
    {
//         cout << " c " << c << endl;
//         cout << " t " << t << endl;
//         cout << " ctop " << ctop << endl;
//         cout << " L " << L << endl;
//         cout << " h " << h << endl;
//         cout << " n " << n << endl;
//         cout << " dc " << dc << endl;
    }
    
    std::complex<double> operator()(double x) {
        std::complex<double> base = exp(std::complex<double>(0.0,-x*ctop));
        std::complex<double>  fac = exp(std::complex<double>(0.0,x*dc));
        std::complex<double> sum(0.0,0.0);
        for (int i=0; i<n; i++) {
            double W = -ctop + i*dc;
            //std::complex<double> arg(0.0,x*W);
            //std::complex<double> base = std::exp(arg);
            sum += ff(W)*base;
            base *= fac;
        }
        return sum*dc*0.5/madness::constants::pi;
    }
};
    
int main() {
    BandlimitedPropagator bp(31.4, 0.07);
    cout << bp(0.1) << endl;
}

