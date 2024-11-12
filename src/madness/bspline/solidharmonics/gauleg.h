#ifndef GAULEG__H
#define GAULEG__H

#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <myqd.h>

template <typename T>
class GaussLegendre {
    size_t npts;
    T a, b;
    std::vector<T> x, w;

    static T mypow(const T& x, int degree) {
        // dd/qd cannot handle degree=0
        if (degree == 0) {
            return T(1);
        }
        else {
            return std::pow(x,degree);
        }
    }

public:
    GaussLegendre(size_t npt, T a=-1, T b=1)
        : npts(npt), a(a), b(b), x(npts), w(npts)
    {
        if (npts > 120) throw "max order is 120";
        std::ifstream f("gauleg.220bit");

        // scaling from x in [-1,1] to z in [a,b]
        // z = alpha + beta*x   and   dz = beta*dx
        bool doscale = (a!=T(-1)) or (b!=T(1));
        T alpha = T(0.5)*(b+a);
        T  beta = T(0.5)*(b-a);

        size_t mpts = 0;
        while (mpts != npts) {
            f >> mpts;
            //std::cout << "trying " << mpts << std::endl;
            for (size_t i=0; i<mpts; i++) {
                size_t junk;
                f >> junk >> x[i] >> w[i];
                if (doscale) {
                    x[i] = beta*x[i] + alpha;
                    w[i] = beta*w[i];
                }
            }
        }

        if (!test()) throw "test failed: did you compile dd/qd code with fast math options?";
    }

    const std::vector<T> pts() const {return x;}

    const std::vector<T> wts() const {return w;}

    size_t npt() const {return npts;}

    template <typename functionT>
    T integrate(functionT& f) const {
        T result = 0;
        for (size_t i=0; i<npts; i++) {
            result += w[i] * f(x[i]);
        }
        return result;
    }

    bool test() const {
        // check powers of x are exactly integrated up to 2*npts-1
        // must scale x so that high powers are representable
        // int((x/s)^k,x=a..b) = (b(b/s)^k - a(a/s)^k)/(k+1)

        T s = std::max(std::abs(a),std::abs(b));
        T maxneps = 0;
        for (int degree=0; size_t(degree)<2*npts-1; degree++) {
            auto f = [&](T x){return mypow(x/s,degree);};
            T exact = (b*mypow(b/s,degree) - a*mypow(a/s,degree))/T(degree+1);
            if (degree == 0) exact = (b-a);
            T value = integrate(f);
            T abserr = std::abs(exact-value);
            T relerr = std::abs(exact) > 0 ? abserr/std::abs(exact) : abserr;
            T neps = relerr/std::numeric_limits<T>::epsilon();
            //std::cout << npts << " " << degree << " " << value << " " << exact << " " << abserr << " " << relerr << " " << neps << std::endl;
            maxneps = std::max(neps, maxneps);
        }
        bool status = maxneps < (3.0*npts); // dd_real seems to be the problem; double and qd_real both pass with 1.0? 
        if (!status) {
            std::cout << "GaussLegendre failed:    type " << typeid(a).name() << "   npts " << npts << "   maxrelerr/epsilon " << maxneps << std::endl;
        }
        return status;
    }

    static bool testall() {
        T lo = -2.0;
        T hi = 5.0;
        for (size_t n=1; n<=120; n++) {
            if (!GaussLegendre<T>(n, lo, hi).test()) return false;
        }
        return true;
    }            
};

template <typename T>
class GaussLegendreSphere {
    size_t order; // Maximum angular momentum exactly integrated
    size_t nphi;  // Number of phi points is just 2*order
    size_t ntheta;// Number of theta points is ceil((order+1)/2)
    size_t npts;  // Total number of points in the rule = nphi*ntheta
    std::vector<T> thetas;
    std::vector<T> sinthetas;
    std::vector<T> costhetas;
    std::vector<T> phis;
    std::vector<T> sinphis;
    std::vector<T> cosphis;
    std::vector<T> weights; // GL weights * 2*Pi/nphi

    // int(f(theta) * sin(theta) dtheta, theta = 0.. pi) = int(f(acos(x)) dx, x=-1..1) ... use GL for this

    // int(f(phi) dphi, phi=0..2pi) for periodic f is accurate with trapezoid rule
    // if f(phi) = exp(I m phi) for integer m, then to integrate need up to 2m points with equal qeights

    // This rule has about 2/3 efficiency (Lebedev is close to 1 but clusters near the poles), Beylkin's rule fixes that issue.
    // But advantage of this GL+trap rule is we can form the quadrature points+weights in arbitrary precision

public:

    GaussLegendreSphere(size_t lmax)
        : order(lmax)
        , nphi(std::max(2*lmax,size_t(1)))
        , ntheta(lmax/2+1)
        , npts(nphi*ntheta)
        , thetas(ntheta)
        , sinthetas(ntheta)
        , costhetas()
        , phis(nphi)
        , sinphis(nphi)
        , cosphis(nphi)
        , weights()
    {
        T pi = from_str<T>("3.141592653589793238462643383279502884197169399375105820974944592307816");
        GaussLegendre<T> g(ntheta);
        costhetas = g.pts();
        weights = g.wts();
        for (size_t i=0; i<ntheta; i++) {
            thetas[i] = std::acos(costhetas[i]);
            //std::cout << i << " " << thetas[i] << std::endl;
            sinthetas[i] = std::sin(thetas[i]);
            weights[i] = weights[i]*2*pi/nphi;
        }
        for (size_t i=0; i<nphi; i++) {
            phis[i] = (i*2)*pi/nphi;
            sinphis[i] = std::sin(phis[i]);
            cosphis[i] = std::cos(phis[i]);
        }

        T sum = 0;
        for (size_t j=0; j<ntheta; j++) {
            sum += weights[j] * nphi;
        }
    }

    template <typename functionT>
    auto cartesian_integral(functionT f, T r=1) const {
        auto value = f(r,T(0),T(0)); // just to deduce the return type using auto
        value = T(0);
        for (size_t j=0; j<ntheta; j++) {
            T weight = weights[j];
            for (size_t i=0; i<nphi; i++) {
                T x = r*sinthetas[j]*cosphis[i];
                T y = r*sinthetas[j]*sinphis[i];
                T z = r*costhetas[j];
                value += f(x,y,z) * weight;
            }
        }
        return value * r * r;
    }

    // Cartesian points on unit sphere
    std::tuple<std::vector<T>,std::vector<T>,std::vector<T>>
    cartesian_pts() {
        std::vector<T> x(nphi*ntheta);
        std::vector<T> y(nphi*ntheta);
        std::vector<T> z(nphi*ntheta);
        for (size_t j=0; j<ntheta; j++) {
            for (size_t i=0; i<nphi; i++) {
                size_t ji = j*nphi + i;
                x[ji+i] = sinthetas[j]*cosphis[i];
                y[ji+i] = sinthetas[j]*sinphis[i];
                z[ji+i] = costhetas[j];
            }
        }
        return {x,y,z};
    }
};

#endif
