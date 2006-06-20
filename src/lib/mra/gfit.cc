#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::max;
using std::min;

#include <cmath>

#include <tensor/tensor.h>

namespace madness {
    using std::abs;

    static const double pi = 3.14159265358979323846264338328;

    /// Form separated expansion for BSH (and Coulomb, mu=0) via trapezoidal quadrature.

    /// For Coulomb, prune, assuming hi = 1.0 (was 2sqrt(3))
    /// See Harrison et al in LNCS for details
    void bsh_fit(double mu, double lo, double eps,
                 Tensor<double> *pcoeff, Tensor<double> *pexpnt,
                 bool prnt=false) {
        double T;
        double slo, hi, shi;

        if (eps >= 1e-2) T = 5;
        else if (eps >= 1e-4) T = 10;
        else if (eps >= 1e-6) T = 14;
        else if (eps >= 1e-8) T = 18;
        else T = 25;

        hi = 1.0;			// was 2.0*sqrt(3.0)
        if (mu > 0) {
            slo = -0.5*log(4.0*T/(mu*mu));
        } else {
            slo = log(eps/hi) - 1.0;
        }

        shi = 0.5*log(T/(lo*lo));
        double h = 1.0/(.2-.50*log10(eps)); // was 0.47

        long npt = long((shi-slo)/h);
        h = (shi-slo)/(npt+1.0);

        if (prnt)
            cout << "slo " << slo << "shi " << shi << "npt " << npt << endl;

        Tensor<double> coeff(npt), expnt(npt);

        for (int i=0; i<npt; i++) {
            double s = slo + h*(npt-i);	// i+1
            coeff[i] = h*2.0/sqrt(pi)*exp(-mu*mu*exp(-2.0*s)/4.0)*exp(s);
            coeff[i] = coeff[i]/(4.0*pi);
            expnt[i] = exp(2.0*s);
            if (prnt) cout << i << " " << coeff[i] << " " << expnt[i] << endl;
        }

        // Prune small exponents from Coulomb fit.  Evaluate a gaussian at
        // the range midpoint, and replace it there with the next most
        // diffuse gaussian.  Then examine the resulting error at the two
        // end points ... if this error is less than the desired
        // precision, can discard the diffuse gaussian.

        if (mu == 0.0) {
            double mid = sqrt(0.5*(hi*hi + lo*lo));
            long i;
            for (i=npt-1; i>0; i--) {
                double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
                double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
                               cnew*exp(-expnt[i-1]*lo*lo);
                double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
                               cnew*exp(-expnt[i-1]*hi*hi);
                if (max(abs(errlo),abs(errhi)) > 0.03*eps) break;
                npt--;
                coeff[i-1] = coeff[i-1] + cnew;
            }
            coeff = coeff(Slice(0,npt-1));
            expnt = expnt(Slice(0,npt-1));
        }

        if (prnt) {
            long npt = 300;
            double hi = 1.0;
            if (mu) hi = min(1.0,30.0/mu);
            cout << "       x         value   abserr   relerr" << endl;
            cout << "  ------------  ------- -------- -------- " << endl;
            double step = exp(log(hi/lo)/(npt+1));
            for (int i=0; i<=npt; i++) {
                double r = lo*(pow(step,i+0.5));
                double exact = exp(-mu*r)/r/4.0/pi;
                double test = 0.0;
                for (int j=0; j<coeff.dim[0]; j++)
                    test += coeff[j]*exp(-r*r*expnt[j]);
                double err = 0.0;
                if (exact) err = (exact-test)/exact;
                printf("  %.6e %8.1e %8.1e %8.1e\n",r, exact, exact-test, err);
            }
        }
        *pcoeff = coeff;
        *pexpnt = expnt;
    }

    // int main() {
    //   Tensor<double> coeff, expnt;
    //   bsh_fit(0.0,1e-8,1e-6,&coeff,&expnt,true);
    //   return 0;
    // }

}
