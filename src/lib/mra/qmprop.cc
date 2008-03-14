#include <cmath>
#include <complex>
#include <constants.h>
#include <iostream>

/// Class to evaluate the filtered Schrodinger free-particle propagator in real space

/// Follows the corresponding Maple worksheet and the implementation notes.
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
    typedef double_complex returnT;

    BandlimitedPropagator(double c, double t)
        : c(c)
        , t(t)
        , ctop(2.15*c)
        , L(1.2*0.5435*(3.0*pow(c,5.0/3.0)*pow(t,0.75)+400.0)/c) // 1.2 --> 2.0 ??
        , h(3.14/ctop)
        , n(2.0*L/h+1)
        , dc(2*ctop/(n-1))
    {
//         std::cout << " c " << c << std::endl;
//         std::cout << " t " << t << std::endl;
//         std::cout << " ctop " << ctop << std::endl;
//         std::cout << " L " << L << std::endl;
//         std::cout << " h " << h << std::endl;
//         std::cout << " n " << n << std::endl;
//         std::cout << " dc " << dc << std::endl;
    }
    
    std::complex<double> operator()(double x) {
        std::complex<double> base = exp(std::complex<double>(0.0,-x*ctop));
        std::complex<double>  fac = exp(std::complex<double>(0.0,x*dc));
        std::complex<double> sum(0.0,0.0);
        double W = -ctop;
        for (int i=0; i<n; i++,W+=dc) {
            //std::complex<double> arg(0.0,x*W);
            //std::complex<double> base = std::exp(arg);
            sum += ff(W)*base;
            base *= fac;
        }
        return sum*dc*0.5/madness::constants::pi;
    }

    static void test() {
        std::complex<double> maple(1.13851441120840,-.986104972277800);
        BandlimitedPropagator bp(31.4, 0.07);
        if (std::abs(bp(0.1)-maple) > 1e-12) throw "BandlimitedPropagator: failed test";
        return;
    }
};



    
// int main() {
//     BandlimitedPropagator::test();
//     return 0;
// }
