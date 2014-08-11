#include <madness/tensor/tensor_lapack.h>
#include <madness/tensor/solvers.h>
#include <madness/world/print.h>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace madness;
using namespace std;

double myexp(double x) {
    if (x>0.0) return 1.0;
    else return exp(x);
}

Tensor<double> optimize_coeffs(const Tensor<double>& f,
                               const Tensor<double>& x,
                               const Tensor<double>& w,
                               const Tensor<double>& expnt,
                               const double mu,
                               const double p,
                               const bool prnt=0) {
    const long nx = x.size();
    const long ng = expnt.size();

    Tensor<double> A(nx+ng,ng), b(nx+ng);
    for (int i=0; i<nx; i++) {
        const double sqrtw = sqrt(w[i]);
        for (int j=0; j<ng; j++) {
            A(i,j) = myexp(-expnt(j)*x(i)*x(i)) * sqrtw;
        }
        b[i] = f[i] * sqrtw;
    }

    const double smu = sqrt(mu);
    for (int j=0; j<ng; j++) {
        double e = expnt(j);
        if (e<0.0) e = 1e-2;
        A(j+nx,j) = smu * pow(e, p*0.5);
    }

    const double rcond = 1e-14;
    Tensor<double> c, s, sumsq;
    long rank;
    gelss(A, b, rcond, c, s, rank, sumsq);

    ////print("A\n",A);

    // print(" solution ", c);
    // print("        s ", s);
    // print("     rank ", rank, A.dim(0), A.dim(1));
    // print("    sumsq ", sumsq);

    return c;
}

double reciprocal(double x) {return 1.0/x;}

double square(double x) {return x*x;}

Tensor<double> geometric_series(long n, double a, double r) {
    Tensor<double> q(n);
    for (int i=0; i<n; i++, a*=r) q[i] = a;
    return q;
}

Tensor<double> arithmetic_series(long n, double a, double h) {
    Tensor<double> q(n);
    for (int i=0; i<n; i++, a+=h) q[i] = a;
    return q;
}

Tensor<double> map_tensor(double (*f)(double), const Tensor<double>& x) {
    Tensor<double> y = copy(x);
    y.unaryop(f);
    return y;
}

double fit(double x, const Tensor<double>& c, const Tensor<double>& expnt) {
    const int n = expnt.size();
    double sum = 0.0;
    for (int i=0; i<n; i++) sum += c[i] * myexp(-expnt[i]*x*x);
    return sum;
}

double errsq(const Tensor<double> f, 
             const Tensor<double> x, 
             const Tensor<double> w,
             const Tensor<double> c, 
             const Tensor<double> expnt,
             const double mu,
             const double p) {
    const int n = x.size();
    double sum = 0.0;
    for (int i=0; i<n; i++) {
        double err = f[i] - fit(x[i], c, expnt);
        sum += w[i]*err*err;
    }

    for (int j=0; j<expnt.size(); j++) sum += mu*c[j]*c[j]*pow(expnt[j],p);

    return sum;
}

class Fred : public OptimizationTargetInterface {
    const Tensor<double> f, x, w;
    const double mu;
    const double p;
    const double alpha;
    const long nx;

    /// Makes the matrix of Gaussians g(i,j) = myexp(-expnt[j]*x[i]x[i])
    Tensor<double> make_g(const Tensor<double>& expnt) {
        const long ng = expnt.size();
        Tensor<double> g(nx,ng);
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ng; j++) {
                g(i,j) = myexp(-expnt[j]*x[i]*x[i]);
            }
        }
        return g;
    }

    double penalty(double expnt) const {
        if (expnt < 0) return alpha*expnt*expnt;
        else return 0.0;
    }

    double dpenalty(double expnt) const {
        if (expnt < 0) return 2.0*alpha*expnt;
        else return 0.0;
    }

public:

    Fred(const Tensor<double> f, 
         const Tensor<double> x, 
         const Tensor<double> w,
         const double mu,
         const double p)
        : f(f), x(x), w(w), mu(mu), p(p), alpha(1e12), nx(x.size())
    {}

    bool provides_gradient() const {return true;}

    void value_and_gradient(const Tensor<double>& expnt,
                            double& value,
                            Tensor<double>& gradient) 
    {
        const long ng = expnt.size();
        Tensor<double> c = optimize_coeffs(f, x, w, expnt, mu, p);

        Tensor<double> g = make_g(expnt);

        gradient = Tensor<double>(ng);

        double errsq = 0.0;
        for (int i=0; i<nx; i++) {
            double v = 0.0;
            for (int j=0; j<ng; j++) {
                v += g(i,j)*c[j];
            }
            double eps = f[i] - v;
            errsq += w[i]*eps*eps;

            eps *= 2.0* w[i] * x[i]*x[i];
            for (int j=0; j<ng; j++) {
                gradient[j] += eps * c[j] * g(i,j);
            }
        }

        for (int j=0; j<ng; j++) {
            double e = expnt(j);
            if (e<0.0) e = 1e-2;

            errsq += mu*c[j]*c[j]*pow(e,p) + penalty(expnt[j]);
            gradient[j] += p*mu*c[j]*c[j]*pow(e,p-1.0) + dpenalty(expnt[j]);
        }

        value = errsq;
    }

    double value(const Tensor<double>& expnt) {
        Tensor<double> gradient;
        double errsq;

        value_and_gradient(expnt, errsq, gradient);
        return errsq;
    }

    Tensor<double> gradient(const Tensor<double>& expnt) {
        double errsq;
        Tensor<double> g;
        value_and_gradient(expnt, errsq, g);
        return g;
    }
};

const double m = 1;
const double c = 137.0359895;
const double mc2 = m*c*c;

double E0(double p) {
    return sqrt(p*p*c*c + mc2*mc2);
}

double Akernel(double p) {
    const double e0 = E0(p);
    const double rsqrt2 = 1.0/sqrt(2.0);
    return sqrt((e0+mc2)/(2.0*e0)) - rsqrt2;
}

int main() {
    // This for fitting 1/x over large dynamic range [xlo,xhi]
    // const long ng = 32;
    // const double mu = 1e-14; // penalty term is sum(i) mu * c[i]^2 * expnt[i]^p ... 1e-12 OK
    // const double p = -0.5;
    // const double xlo = 1e-6;
    // const double xhi = 1e0;
    // const long nx = max(10*long(log10(xhi/xlo)), 10*ng);
    // const double r = pow(xhi/xlo,1.0/(nx-1));
    // Tensor<double> x = geometric_series(nx,xlo,r);
    // Tensor<double> f = map_tensor(reciprocal, x);
    // Tensor<double> w = map_tensor(square, x);   
    // const double ahi = -log(0.1)/(xlo*xlo); // 0.1 = myexp(-ahi*xlo^2)
    // const double alo = -log(0.9)/(xhi*xhi); // 0.9 = myexp(-alo*xhi^2)
    // const double ra = pow(ahi/alo,1.0/(ng-1));
    // Tensor<double> expnt = geometric_series(ng,alo,ra);

    // This for fitting A(p) over [0,10c]
    const long ng = 50;
    const double mu = 1e-12;
    const double p = 0.0; // was 0

    const long nx = 2000;

    //const double xlo = 0.0;
    //const double xhi = 10000.0*c;
    //const double h = (xhi-xlo)/(nx-1);
    //Tensor<double> x = arithmetic_series(nx,xlo,h);

    const double xlo = 1e-6;
    const double xhi = 10000.0*c;
    const double h = pow(xhi/xlo,1.0/(nx-1));
    Tensor<double> x = geometric_series(nx,xlo,h);
    Tensor<double> f = map_tensor(Akernel, x);

    //Tensor<double> w(nx);
    //w.fill(1.0);

    Tensor<double> w = map_tensor(sqrt,x);
    w = map_tensor(reciprocal,w);

    const double ahi = 1e-3; //1.0/(mc2); // was mc2**2
    const double alo = log(2.0)/xhi/xhi; //1e-8; //1e-4*ahi;
    const double ra = pow(ahi/alo,1.0/(ng-1));
    Tensor<double> expnt = geometric_series(ng,alo,ra);

    print(" x ", x);
    print(" w ", w);
    print(" f ", f);
    print(" expnt ", expnt);

    Tensor<double> c = optimize_coeffs(f, x, w, expnt, mu, p);
    print(" solution ", c);

    double epssq = errsq(f, x, w, c, expnt, mu, p);
    print(" errsq", epssq,"\n\n");

    for (int i=0; i<nx; i++) {
        print(x[i], f[i], fit(x[i],c,expnt), f[i] - fit(x[i],c,expnt), (f[i] - fit(x[i],c,expnt))/f[i]);
    }

    shared_ptr<Fred> fptr(new Fred(f, x, w, mu, p)); 
    print("fred ", fptr->value(expnt));
    fptr->test_gradient(expnt, 1e-12*epssq, true);
    //return 0;

    // SteepestDescent sopt(fptr, 1e-12, 1e-20, 1e-20); 
    // sopt.optimize(expnt);
    // sopt.optimize(expnt);
    // sopt.optimize(expnt);
    // sopt.optimize(expnt);
    // sopt.optimize(expnt);

    QuasiNewton opt(fptr, 20*ng, 1e-15, 1e-24, 1e-24); opt.set_update("BFGS");
    
    for (int iter=0; iter<4; iter++) {
        if (opt.optimize(expnt)) break; 
        opt.reset_hessian();
        sort(expnt.ptr(),expnt.ptr()+expnt.size());
        print("\n optimized exponents ", expnt);
        c = optimize_coeffs(f, x, w, expnt, mu, p);
        print("        coefficients ", c);
    }

    print("COEFF WEIGHTS");
    printf("   ");
    for (double pp=-2; pp<3.01; pp += 0.5) printf("%.1e ",pp);
    printf("\n");
    for (int i=0; i<c.size(); i++) {
        printf("%2d ", i);
        for (double pp=-2; pp<3.01; pp += 0.5) {
            printf("%.1e ", c[i]*c[i]*pow(expnt[i],pp));
        }
        printf("\n");
    }
    
    double maxerr = 0.0;
    for (int i=0; i<nx; i++) {
        double err = f[i] - fit(x[i],c,expnt);
        print(x[i], f[i], fit(x[i],c,expnt), err , err/f[i]);
        maxerr = max(abs(err),maxerr);
    }
    print("maxerr", maxerr);

    return 0;
}


    
    
