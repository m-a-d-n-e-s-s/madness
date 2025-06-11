#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h> //for erfc

#include "../../../../mpfrc++-3/mpreal.h"

static const double m = 1.0;
static const double c = 137.03599917697017; // from CODATA 2022
static const double mc2 = m*c*c;
static const double PI = 3.14159265358979323846264338328;

/// Exponents in Fourier space quadrature over t
double q(double t) {
    return std::exp(-t)/(m*m*c*c);
}

// Convert Fourier gaussian exponent a into real space exponent
double Q(double a) {
    return 0.25/a;
}

/// Additional factor when transforming gaussian exponent a from fourier to real space
double C(double a) {
    return 0.125/std::pow(a*PI,1.5);
}


/// This fast summation of an oscillating series is based upon Pade.

/// http://people.mpim-bonn.mpg.de/zagier/files/exp-math-9/fulltext.pdf
///
/// For Tbar I think it works for the right reasons for t>0 but for large negative t I think it gets lucky.
/// C(k) is a function/functor that returns the k'th term in the series for k=0,...,nmax
template <typename T>
double fastsum_oscillating(const T& C, int nmax) {
    double d, b, c, s, k;

    d = std::pow(3.0+sqrt(8.0),nmax);
    d = (d+1.0/d)*0.5;
    b = -1.0;
    c = -d;
    s = 0.0;
    for (k=0; k<nmax; k++) {
        c = b-c;
        s = s+c*C(k);
        b = (k+nmax)*(k-nmax)*b/((k+0.5)*(k+1.0));
    }
    return s/d;
}

/// Computes term k at quadrature point t in infinte Fourier-space sum defining Tbar
class Tbar_omega {
    const double t;
public:
    Tbar_omega(double t) : t(t) {}

    double operator()(int k) const {
        return 2.0*std::exp(-std::exp(-t)-0.5*(k+1)*t - std::lgamma(0.5*(k+1)));
    }
};


/// Computes term k at quadrature point t in infinte Fourier-space sum defining A
class A_omega {
    const double t;
public:
    A_omega(double t) : t(t) {}

    double operator()(int k) const {
        return 0.25*exp(-exp(-t) - 0.5*(k+1)*t + log(2.0)*(0.5-2*k) + lgamma(1.0+2*k) - 2.0*lgamma(k+1.0) - lgamma(0.5*(k+1))) / ((k+1.0));
    }
};


/// Computes term k at quadrature point t in infinte Fourier-space sum defining APbar
class APbar_omega {
    const double t;
public:
    APbar_omega(double t) : t(t) {}

    double operator()(int k) const {
        return sqrt(0.5)*exp( -exp(-t) - log(4.0)*k - 0.5*(k+1)*t + lgamma(1.0+2*k) - 2.0*lgamma(k+1.0) - lgamma(0.5*(k+1))) / (c*m);
    }
};


/// Computes term k at quadrature point t in infinte Fourier-space sum defining Pbar

/// Use inifinte sum for Pbar even though sum can be eliminated only to keep the same construction as used for A, AP, Tbar etc.
class Pbar_omega {
    const double t;
public:
    Pbar_omega(double t) : t(t) {}

    double operator()(int k) const {
        return exp(-exp(-t) - 0.5*(k+1)*t - lgamma(0.5*(k+1))) / (c*m);
    }
};


/// Computes at quadrature point t the infinte Fourier-space sum defining Tbar
double tbar_OMEGA(double t) {
    if (t < -4.0 || t > 100.0) {
        return 0.0;
    }
    else {
        double a = fastsum_oscillating(Tbar_omega(t), 35);
        //return (2.0*c/sqrt(PI)*exp(t-c*c*exp(-t))-2.0*c*c*exp(t/2.0)*erfc(c*exp(-t/2.0)))/pow(2*PI,3.0); //doesn't work
        double b = 2.0*c/sqrt(PI)*exp(-t/2.0)*exp(-c*c*exp(-t))-2.0*c*c*exp(-t)*erfc(c*exp(-t/2.0));
        //printf("t: %f \nold: %f \nnew: %f\n\n",t,a,b);
        return b;
    }
}


/// Computes at quadrature point t the infinte Fourier-space sum defining A
double a_OMEGA(double t) {
    if (t < -4.0 || t > 100.0) {
        return 0.0;
    }
    else {
        return fastsum_oscillating(A_omega(t), 35);
    }
}

/// Computes at quadrature point t the infinte Fourier-space sum defining Pbar

/// Computes at quadrature point t the infinte Fourier-space sum defining Pbar
double pbar_OMEGA(double t) {
    if (t < -4.0 || t > 100.0) {
        return 0.0;
    }
    else {
        //return fastsum_oscillating(Pbar_omega(t), 35);
        double b = 1.0/sqrt(PI)*exp(-t/2.0)*exp(-c*c*exp(-t))-c*exp(-t)*erfc(c*exp(-t/2.0));
        return b;
    }
}

/// Computes at quadrature point t the infinte Fourier-space sum defining APbar
double apbar_OMEGA(double t) {
    if (t < -4.0 || t > 100.0) {
        return 0.0;
    }
    else {
        return fastsum_oscillating(APbar_omega(t), 35);
    }
}


// adjust quadrature points onto common mesh for all operators
void  munge_quadrature_points(int& npt, double& tlo, double& thi, double& h) {
    //std::cout << "munge before " << npt << " " << tlo << " " << thi << " " << h << std::endl;
    h = floor(64.0*h)/64.0;

    // Round thi/lo up/down to an integral multiple of quadrature
    // points by using a common origin below that used by any operator

    const double base = -1000.0;
    tlo = base + floor((tlo-base)/h)*h;
    thi = base +  ceil((thi-base)/h)*h;
    npt = ceil(thi-tlo)/h+1;
    //std::cout << "munge after  " << npt << " " << tlo << " " << thi << " " << h << std::endl;
}

double relops_h(const double quadacc) {
    return -5.0/std::log10(std::min(1e-6,quadacc));
} 
    

void tbar_fit(const double dx, const double thresh, const double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta) {
    const double expnt_max = 1.0/(dx*dx);
    
    double h = relops_h(quadacc); //CHANGE THIS TO REFLECT OPERATOR USED IN MAPLE CODE
    double tlo = -4.0;
    double thi = 90.0;
    int npt;

    munge_quadrature_points(npt, tlo, thi, h);

    coeffs.clear();
    expnts.clear();

    double sum = 0.0;
    double sum_trunc = 0.0;

    Cdelta = 0.0;

    for (int i=0; i<npt; i++) {
        double t = i*h+tlo;
        double qt = q(t);
        //double coeff = h * C(qt) * tbar_OMEGA(t);  //Robert's
        //double expnt = Q(qt);                      //Robert's
        double coeff = h * C(exp(-t)) * tbar_OMEGA(t);        //Joel's
        double expnt = std::exp(t)/4.0;          //Joel's
        double norm = std::pow(PI/expnt,1.5);

        sum += coeff*norm;

        if (expnt > expnt_max) {
            Cdelta += coeff*norm; 
        }
        else if (coeff > 0.1*thresh) {
            coeffs.push_back(coeff);
            expnts.push_back(expnt);

            sum_trunc += coeff*norm;
            
            //std::cout << coeff << " " << expnt << std::endl;
        }
    }

    double trace = 1.0;
    
    // std::cout << "Tbar: " << npt << " " << coeffs.size() << std::endl;
    // std::cout << "Tbar: " << sum << " " << std::abs(trace-sum) << std::endl;
    // std::cout << "Tbar: " << sum_trunc << " " << std::abs(trace-sum_trunc) << std::endl;
    // std::cout << "Tbar: " << sum_trunc+Cdelta << " " << std::abs(trace-sum_trunc-Cdelta)<< std::endl;
    // std::cout << std::endl;
}

void pbar_fit(const double dx, const double thresh, const double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta) {
    const double expnt_max = 1.0/(dx*dx);
    
    double h = relops_h(quadacc);
    double tlo = -4.0;
    double thi = 90.0;
    int npt;

    munge_quadrature_points(npt, tlo, thi, h);

    coeffs.clear();
    expnts.clear();

    double sum = 0.0;
    double sum_trunc = 0.0;

    Cdelta = 0.0;

    for (int i=0; i<npt; i++) {
        double t = i*h+tlo;
        double qt = q(t);
        //double coeff = h * C(qt) * pbar_OMEGA(t);
        //double expnt = Q(qt);
        double coeff = h * C(exp(-t)) * pbar_OMEGA(t);        //Joel's
        double expnt = std::exp(t)/4.0;          //Joel's
        double norm = std::pow(PI/expnt,1.5);

        sum += coeff*norm;

        if (expnt > expnt_max) {
            Cdelta += coeff*norm;
        }
        else if (coeff > 0.1*thresh) {
            coeffs.push_back(coeff);
            expnts.push_back(expnt);

            sum_trunc += coeff*norm;
            
            //std::cout << coeff << " " << expnt << std::endl;
        }
    }

    double trace = 0.5/(c*m);
    
    // std::cout << "Pbar: " << npt << " " << coeffs.size() << std::endl;
    // std::cout << "Pbar: " << sum << " " << std::abs(trace-sum) << std::endl;
    // std::cout << "Pbar: " << sum_trunc << " " << std::abs(trace-sum_trunc) << std::endl;
    // std::cout << "Pbar: " << sum_trunc+Cdelta << " " << std::abs(trace-sum_trunc-Cdelta)<< std::endl;
    // std::cout << std::endl;
}

void apbar_fit(const double dx, const double thresh, const double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta) {
    const double expnt_max = 1.0/(dx*dx);
    
    double h = relops_h(quadacc);
    double tlo = -4.0;
    double thi = 90.0;
    int npt;

    munge_quadrature_points(npt, tlo, thi, h);

    coeffs.clear();
    expnts.clear();

    double sum = 0.0;
    double sum_trunc = 0.0;

    Cdelta = 0.0;

    for (int i=0; i<npt; i++) {
        double t = i*h+tlo;
        double qt = q(t);
        double coeff = h * C(qt) * apbar_OMEGA(t);
        double expnt = Q(qt);
        double norm = std::pow(PI/expnt,1.5);

        sum += coeff*norm;

        if (expnt > expnt_max) {
            Cdelta += coeff*norm;
        }
        else if (coeff > 0.1*thresh) {
            coeffs.push_back(coeff);
            expnts.push_back(expnt);

            sum_trunc += coeff*norm;
            
            //std::cout << coeff << " " << expnt << std::endl;
        }
    }

    double trace = 0.5/(c*m);
    
    // std::cout << "Apbar: " << npt << " " << coeffs.size() << std::endl;
    // std::cout << "Apbar: " << sum << " " << std::abs(trace-sum) << std::endl;
    // std::cout << "Apbar: " << sum_trunc << " " << std::abs(trace-sum_trunc) << std::endl;
    // std::cout << "Apbar: " << sum_trunc+Cdelta << " " << std::abs(trace-sum_trunc-Cdelta)<< std::endl;
    // std::cout << std::endl;
}


void a_fit(const double dx, const double thresh, const double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta) {
    const double expnt_max = 1.0/(dx*dx);
    
    double h = relops_h(quadacc);
    double tlo = -4.0;
    double thi = 90.0;
    int npt;

    munge_quadrature_points(npt, tlo, thi, h);

    coeffs.clear();
    expnts.clear();

    double sum = 0.0;
    double sum_trunc = 0.0;

    Cdelta = 0.0;

    for (int i=0; i<npt; i++) {
        double t = i*h+tlo;
        double qt = q(t);
        double coeff = h * C(qt) * a_OMEGA(t);
        double expnt = Q(qt);
        double norm = std::pow(PI/expnt,1.5);

        sum += coeff*norm;

        if (expnt > expnt_max) {
            Cdelta += coeff*norm;
        }
        else if (coeff > 0.1*thresh) {
            coeffs.push_back(coeff);
            expnts.push_back(expnt);

            sum_trunc += coeff*norm;
            
            //std::cout << coeff << " " << expnt << std::endl;
        }
    }

    Cdelta += sqrt(0.5);

    double trace = 1.0;
    // std::cout << "A: " << npt << " " << coeffs.size() << std::endl;
    // std::cout << "A: " << sum_trunc+Cdelta << " " << std::abs(trace-sum_trunc-Cdelta)<< std::endl;
    // std::cout << std::endl;
}


double bshrel_omega(double t, int k, double epsilon) {
    // return
    //     std::pow(mc2, -(k+1.0)) *
    //     std::pow(mc2 + epsilon, (double) k) *
    //     std::exp(-std::exp(-t) - 0.5*(k+1.0)*t - lgamma(0.5*(k+1)));

    // std::cout << std::endl;
    // std::cout << t << " " << k << " " << epsilon << std::endl;
    // std::cout << std::exp(-t) << std::endl;
    // std::cout << 0.5*(k+1.0)*t << std::endl;
    // std::cout << lgamma(0.5*(k+1)) << std::endl;
    // std::cout << std::endl;

    // Note, cancellation for large k and negative t so must use extended precision
    if (t>0.0) {
        // return
        //     std::pow(1 + epsilon/mc2, (double) k) *
        //     std::exp(-std::exp(-t) - 0.5*(k+1.0)*t - lgamma(0.5*(k+1))) / mc2;
        return
            std::exp(-std::exp(-t) - 0.5*(k+1.0)*t - std::lgamma(0.5*(k+1)) + k*std::log(1.0 + epsilon/mc2)) / mc2;
    }
    else {
        mpfr::mpreal mt=t, mk=k, half=0.5, one=1.0, meps=epsilon, mmc2=mc2;
        mpfr::mpreal arg = -mpfr::exp(-mt) - half*(mk+one)*mt - mpfr::lngamma(half*(mk+one)) + mk*mpfr::log(1.0 + meps/mmc2);
        mpfr::mpreal val = mpfr::exp(arg)/mmc2;
        return double(val);
    }
}


/*
double bshrel_OMEGA(double t, double epsilon) {
    const double thresh = 1e-20;
    const int kk = std::ceil(2*std::exp(-t));
    double value=0.0;

    int k = kk;
    while (true) {
        double term = bshrel_omega(t, k, epsilon);
        value += term;
        k++;
        if (term < thresh) break;
    }

    k = kk-1;
    while (true) {
        if (k<0) break;
        double term = bshrel_omega(t, k, epsilon);
        value += term;
        k--;
        if (term < thresh) break;
    }

    return value;
}
*/


double bshrel_OMEGA(double t, double epsilon) {
     //printf("bshrel_OMEGA called with t = %f\n",t);
     double result = 1/(c*sqrt(8*PI))*exp(t-c*c*exp(-t))-1/sqrt(8)*(1+epsilon/(c*c))*exp((2*epsilon+epsilon*epsilon/(c*c))*exp(-t)+t/2.0)*(erfc((c+epsilon/c)*exp(-t/2.0))-2);
     result = result/pow(2*PI,3.0/2.0);
     return result;
}


void bshrel_fit(double epsilon, const double dx, const double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta) {
    const double tlo = -5.0;
    const double thi = 80.0;
    //const int npt = 200; // about 4 digits in 2-norm of operator
    //const int npt = 300; // about 7 digits
    // const int npt = 400; // about 9 digits // was using this 
    //const int npt = 500; // about 11 digits
    //const int npt = 600; // about 12 digits
    //const int npt = 800;
     const int npt = 350;

    const double h = (thi-tlo)/(npt-1);

    const double expnt_max = 1.0/(dx*dx);

    //std::cout << "h " << h << std::endl;

    coeffs.clear();
    expnts.clear();

    double sum = 0.0;

    double sum_trunc = 0.0;

    Cdelta = 0.0;
    
    for (int i=0; i<npt; i++) {
        double t = i*h+tlo;
        double qt = q(t);
        //double coeff = h * C(qt) * bshrel_OMEGA(t,epsilon);
        double coeff = h * bshrel_OMEGA(t,epsilon);
        //double expnt = Q(qt);
        double expnt = std::exp(t)/4.0;
        double norm = std::pow(PI/expnt,1.5);

        sum += coeff*norm;

        if (expnt > expnt_max) {
            Cdelta += coeff*norm;
        }
        else if (coeff > 0.1*thresh) {
            coeffs.push_back(coeff);
            expnts.push_back(expnt);

            sum_trunc += coeff*norm;
            
            //std::cout << coeff << " " << expnt << std::endl;
        }
    }

    //std::cout << sum << " " << -1.0/epsilon << " " << std::abs((sum+1.0/epsilon)*epsilon) << std::endl;
    //std::cout << sum_trunc << " " << -1.0/epsilon << " " << std::abs((sum_trunc+1.0/epsilon)*epsilon) << std::endl;
    //std::cout << sum_trunc << " " << -1.0/epsilon << " " << std::abs((sum_trunc+1.0/epsilon)*epsilon) << std::endl;
}

bool check(double value, double correct, double thresh, const char* msg) {
    double err = std::abs(value-correct);
    if (err<=thresh) {
        std::cout << msg << " : passed " << std::endl;
        return true;
    }
    else {
        std::cout << msg << " : failed |" << value << "-" << correct << "| = " << err << " > " << thresh << std::endl;
        return false;
    }
}

// int main() {
//     mpfr::mpreal::set_default_rnd(GMP_RNDN);
//     mpfr::mpreal::set_default_prec(256);

// //     // std::cout.precision(16);
// //     // check(bshrel_omega(0.1,3,-0.01),0.0000176402701584786376862159169871,1e-20,"omega(.1, 3, -0.1e-1)");
// //     // check(bshrel_OMEGA(1.0,-0.01),0.000044133601851009832964,1e-19,"OMEGA(1.0,-0.01)");
// //     // check(bshrel_OMEGA(0.1,-0.01),0.000099327541717740664618,1e-19,"OMEGA(0.1,-0.01)");
// //     // check(bshrel_OMEGA(-1.0,-0.01),0.000289917644952772489991,1e-19,"OMEGA(-1.0,-0.01)");
// //     // check(bshrel_OMEGA(-10.0,-0.01),2.29148640330432807356574,1e-19,"OMEGA(-10.0,-0.01)");

//     std::vector<double> coeffs, expnts;
//     double Cdelta;
// //     // bshrel_fit(-6500.0,coeffs,expnts,Cdelta);
// //     // bshrel_fit(-0.1,coeffs,expnts,Cdelta);

//     check(A_omega(0.1)(3),0.00914998206504088871,1e-14,"Aomega(0.1)(3)");
//     check(Pbar_omega(0.1)(3),0.0024173557380369478775,1e-14,"Pbaromega(0.1)(3)");
//     check(APbar_omega(0.1)(3),0.000534165198408167875,1e-14,"APbaromega(0.1)(3)");
//     check(apbar_OMEGA(0.1),0.0006128434915269844710,1e-14,"apbarOMEGA(0.1)");
        
//     // tbar_fit(1e-6, 1e-6, 1e-6, coeffs,expnts,Cdelta);
//     // tbar_fit(1e-6, 1e-6, 1e-8, coeffs,expnts,Cdelta);
//     // tbar_fit(1e-6, 1e-6, 1e-10, coeffs,expnts,Cdelta);
//     // tbar_fit(1e-6, 1e-6, 1e-12, coeffs,expnts,Cdelta);
//     // tbar_fit(1e-6, 1e-6, 1e-14, coeffs,expnts,Cdelta);
//     // tbar_fit(1e-6, 1e-6, 1e-16, coeffs,expnts,Cdelta);
 
//     // a_fit(1e-6, 1e-6, 1e-6, coeffs,expnts,Cdelta);
//     // a_fit(1e-6, 1e-6, 1e-8, coeffs,expnts,Cdelta);
//     // a_fit(1e-6, 1e-6, 1e-10, coeffs,expnts,Cdelta);
//     // a_fit(1e-6, 1e-6, 1e-12, coeffs,expnts,Cdelta);
//     // a_fit(1e-6, 1e-6, 1e-14, coeffs,expnts,Cdelta);
//     // a_fit(1e-6, 1e-6, 1e-16, coeffs,expnts,Cdelta);

//     // pbar_fit(1e-6, 1e-6, 1e-6, coeffs,expnts,Cdelta);
//     // pbar_fit(1e-6, 1e-6, 1e-8, coeffs,expnts,Cdelta);
//     // pbar_fit(1e-6, 1e-6, 1e-10, coeffs,expnts,Cdelta);
//     // pbar_fit(1e-6, 1e-6, 1e-12, coeffs,expnts,Cdelta);
//     // pbar_fit(1e-6, 1e-6, 1e-14, coeffs,expnts,Cdelta);
//     // pbar_fit(1e-6, 1e-6, 1e-16, coeffs,expnts,Cdelta);
 
//     apbar_fit(1e-6, 1e-6, 1e-6, coeffs,expnts,Cdelta);
//     apbar_fit(1e-6, 1e-6, 1e-8, coeffs,expnts,Cdelta);
//     apbar_fit(1e-6, 1e-6, 1e-10, coeffs,expnts,Cdelta);
//     apbar_fit(1e-6, 1e-6, 1e-12, coeffs,expnts,Cdelta);
//     apbar_fit(1e-6, 1e-6, 1e-14, coeffs,expnts,Cdelta);
//     apbar_fit(1e-6, 1e-6, 1e-16, coeffs,expnts,Cdelta);
 
//     return 0;
// }

