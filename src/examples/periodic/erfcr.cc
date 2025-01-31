#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

#include "/home/rjh/Devel/madtran/madness/src/madness/misc/gnuplot.h"

// Gaussian expansion of erfc(a*r)/r accurate to epsilon over [rlo,inf]
// Returns pair [coeffs,expnts]
std::pair<std::vector<double>,std::vector<double>>
make_fit(double a, double epsilon, double rlo) {
    const double pi = 3.14159265358979323846264338328;
    
    // Alpert end-point correction points and weights
    const std::vector<double> x = {0.05899550614325259, 0.3082757062227814, 0.7463707253079130, 1.355993726494664, 2.112943217346336, 2.987241496545946, 3.944798920961176, 4.950269202842798, 5.972123043117706, 6.989783558137742, 7.997673019512965, 8.999694932747039, 9.999979225211805, 10.99999938266130, 11.99999999462073, 13.00000000000000};
    const std::vector<double> omega = {0.1511076023874179, 0.3459395921169090, 0.5273502805146873, 0.6878444094543021, 0.8210319140034114, 0.9218382875515803, 0.9873027487553060, 1.018251913441155, 1.021933430349293, 1.012567983413513, 1.004052289554521, 1.000713413344501, 1.000063618302950, 1.000002486385216, 1.000000030404477, 1.000000000020760};
    const size_t A = 14;

    // The kernel in s
    auto K = [=](double s){return 2*std::exp(-rlo*rlo*std::exp(2*s)+s)/std::sqrt(pi);};
    const double slo = std::log(a);
    double shi = slo;
    while (K(shi) > 1e-100) { // Estimate upper bound 
        shi += 0.01;
    }
    const double T = shi - slo;
    const double hs = pi*pi / (std::log(16*pi) - 2*std::log(epsilon) + 1);
    double ht = hs/T;
    const size_t n = std::ceil(1.0/ht) - A;
    ht = 1.0/(n+A-1);

    std::cout << "n   " << n << std::endl;
    std::cout << "slo " << slo << std::endl;
    std::cout << "shi " << shi << std::endl;
    std::cout << "T   " << T << std::endl;
    std::cout << "hs  " << hs << std::endl;
    std::cout << "ht  " << ht << std::endl;

    auto CC = [&](double t){return 2*a*T*std::exp(T*t)/std::sqrt(pi);};
    auto ZZ = [&](double t){return a*a*std::exp(2*T*t);};
        
    std::vector<double> coeffs,expnts;

    for (size_t i=0; i<x.size(); i++) {
        double t = ht*x[i];
        coeffs.push_back(ht*omega[i]*CC(t));
        expnts.push_back(ZZ(t));
    }
    
    for (size_t i=0; i<n; i++) {
        double t = ht*(A+i);
        coeffs.push_back(ht*CC(t));
        expnts.push_back(ZZ(t));
    }

    // for (size_t i=0; i<coeffs.size(); i++) {
    //     std::cout << i << " " << coeffs[i] << " " << expnts[i] << std::endl;
    // }

    return {coeffs,expnts};
    
}

int main() {
    const double a = 0.25;
    const double epsilon = 1e-6;
    const double rlo = 1e-8;
    auto [c,t] = make_fit(a,epsilon,rlo);

    auto fit = [&](double r) {
        double sum = 0;
        for (size_t i=0; i<c.size(); i++) {
            sum += c[i]*std::exp(-t[i]*r*r);
        }
        return sum;
    };

    // Fill vector R and E with overall error data for plotting
    std::vector<double> R,E;
    double scale = std::pow(10.0,1.0/(5*c.size()));
    for (double r=1e-8; r<30.0; r*=scale) {
        double exact = std::erfc(0.25*r)/r;
        double approx = fit(r);
        double err = std::abs(approx-exact);
        double relerr = err/exact;
        double overallerr = std::min(err,relerr);
        R.push_back(r);
        E.push_back(overallerr);
    }
    madness::Gnuplot gp("set logscale; set style data lines; set xrange [1e-8:30]; set yrange [1e-10:1e-6]; set xlabel 'r'; set ylabel 'min(err,relerr)'; set format x '%.0e'");
    gp.plot(R,E);

    // double scale = std::pow(10.0,1.0/3.0);
    // for (double r=1e-8; r<30.0; r*=scale) {
    //     double exact = std::erfc(0.25*r)/r;
    //     double approx = fit(r);
    //     double err = std::abs(approx-exact);
    //     double relerr = err/exact;
    //     double overallerr = std::min(err,relerr);
    //     std::cout << r << " " << fit(r) << " " << overallerr << std::endl;
    // }
    
    return 0;
}
