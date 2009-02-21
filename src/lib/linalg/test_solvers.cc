#include <linalg/solvers.h>

using namespace madness;
using namespace std;

struct Test : public OptimizationTargetInterface {
    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        return 0.5*3.0*x.sumsq();
    }

    Tensor<double> gradient(const Tensor<double>& x) {
        return x*3.0;
    }
};


struct Test2 : public OptimizationTargetInterface {
    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        double v = 1.0;
        for (int i=0; i<x.dim[0]; i++) {
            v *= cos((i+1)*x[i]);
        }
        return v;
    }

    Tensor<double> gradient(const Tensor<double>& x) {
        double v = value(x);
        Tensor<double> g(x.dim[0]);
        for (int i=0; i<x.dim[0]; i++) {
            g[i]= -v*(i+1)*sin((i+1)*x[i])/cos((i+1)*x[i]);
        }
        return g;
    }
};



Tensor<double> op(const Tensor<double>& x) {
    const long n = x.dim[0];
    Tensor<double> f(n);
    for (long i=0; i<n; i++) {
        f(i) = (i + 1)*x[i]; // + 0.01*i*x[i]*x[i]*x[i];
        for (long j=0; j<n; j++) 
            f(i) += 0.0001*i*j*x[i]*x[i]*x[j]*x[j]/((i+1)*(j+1));
    }
    return f;
}

double dot_product(const Tensor<double>& a, const Tensor<double>& b) {
    return a.trace(b);
}

int main() {

    Tensor<double> x(5);
    x.fillrandom();
    QuasiNewton solver(SharedPtr<OptimizationTargetInterface>(new Test2));
    solver.set_update("SR1");
    solver.optimize(x);
    return 0;


//     int n= 40;
//     int maxiter = 100;
//     int maxnvec = 20;
//     Tensor<double> f(maxnvec,n), x(maxnvec,n), Q(maxnvec,maxnvec);

//     int m = 0;
//     x(0,_).fillrandom();
//     for (int iter=0; iter<maxiter; iter++) {
//         print("\nITERATION", iter, m);
//         f(m,_) = op(x(m,_));
//         print("x");
//         print(x(m,_));
//         print(f(m,_));

//         for (int j=0; j<=m; j++) {
//             Q(j,m) = dot_product(x(j,_),f(m,_));
//             Q(m,j) = dot_product(x(m,_),f(j,_));
//         }

//         Tensor<double> c = KAIN(Q(Slice(0,m),Slice(0,m)));
//         print("THIS IS THE C I GOT");
//         print(c);
        
//         {
//             m++;

//             Tensor<double> xnew(n);
//             for (int j=0; j<m; j++) {
//                 xnew += c(j)*(x(j,_) - f(j,_));
//             }

//             double steplen = (xnew-x(m-1,_)).normf();
//             double xnorm = xnew.normf();
//             if (steplen > 0.5) {
//                 double damp = 0.5/steplen;
//                 if (damp < 0.1) damp = 0.1;
//                 print("restrictING", steplen, xnorm, damp);
//                 xnew = damp*xnew + (1.0-damp)*x(m-1,_);
//             }
            
//             if (m == maxnvec) {
//                 for (int i=1; i<m; i++) {
//                     f(i-1,_) = f(i,_);
//                     x(i-1,_) = f(i,_);
//                 }
//                 Q(Slice(0,-2),Slice(0,-2)) = copy(Q(Slice(1,-1),Slice(1,-1)));
                
//                 m--;
//             }
//             x(m,_) = xnew;
//         }
//     }
//     return 0;

}
