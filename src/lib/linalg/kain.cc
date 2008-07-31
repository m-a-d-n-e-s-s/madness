#include <iostream>
#include <tensor/tensor.h>
#include <linalg/tensor_lapack.h>
#include <world/print.h>

using namespace madness;
using namespace std;



/// Solves the KAIN equations for coefficients to compute the next vector

/// \verbatim
///   Q(i,j) = <xi|fj>
///   A(i,j) = <xi-xm | fj-fm> = Q(i,j) - Q(m,j) - Q(i,m) + Q(m,m)
///   b(i,j) =-<xi-xm | fm> = -Q(i,m) + Q(m,m)
///   A c = b
///
///   . Correction to vector m
///   .   interior = sum(i<m)[ c(i)*(x(i)-x(m)) ] = sum(i<m)[c(i)*x(i)] - x[m]*sum(i<m)[c(i)]
///   .   exterior = -f(m) - sum(i<m)[ c(i)*(f(i)-f(m)) ] = -f(m) - sum(i<m)[c(i)*f(i)] + f(m)*sum(i<m)[c(i)]
///   . New vector
///   .   define C = sum(i<m)(c(i))  (note less than)
///   .   define c(m) = 1.0 - C
///   .   xnew = sum(i<=m) [ c(i)*(x(i) - f(i)) ]
/// \endverbatim
template <typename T>
Tensor<T> KAIN(const Tensor<T>& Q) {
    const int nvec = Q.dim[0];
    const int m = nvec-1;

    if (nvec == 1) {
        Tensor<T> c(1);
        c(0L) = 1.0;
        return c;
    }

    Tensor<T> A(m,m);
    Tensor<T> b(m);
    for (long i=0; i<m; i++) {
        b(i) = Q(m,m) - Q(i,m);
        for (long j=0; j<m; j++) {
            A(i,j) = Q(i,j) - Q(m,j) - Q(i,m) + Q(m,m);
        }
    }

    print("Q");
    print(Q);
    print("A");
    print(A);
    print("b");
    print(b);

    double rcond = 1e-12;
    Tensor<T> x;
    Tensor<double> s;
    long rank;
    gelss(A, b, rcond, &x, &s, &rank);
    print("singular values", s);
    print("rank", rank);
    print("solution", x);

    Tensor<T> c(nvec);
    T sumC = 0.0;
    for (long i=0; i<m; i++) sumC += x(i);
    c(Slice(0,m-1)) = x;
    print("SUMC", nvec, m, sumC);
    c(m) = 1.0 - sumC;

    print("returned C", c);

    return c;
}

Tensor<double> op(const Tensor<double>& x) {
    const long n = x.dim[0];
    Tensor<double> f(n);
    for (long i=0; i<n; i++) {
        f(i) = (i*1)*x[i] + 0.01*i*x[i]*x[i]*x[i];
        for (long j=0; j<n; j++) 
            f(i) += 0.0001*i*j*x[i]*x[i]*x[j]*x[j]/((i+1)*(j+1));
    }
    return f;
}

double dot_product(const Tensor<double>& a, const Tensor<double>& b) {
    double sum = 0.0;
    ITERATOR(a, sum += a(IND)*b(IND));
    return sum;
}

int main() {

    int n= 40;
    int maxiter = 100;
    int maxnvec = 10;
    Tensor<double> f(maxnvec,n), x(maxnvec,n), Q(maxnvec,maxnvec);

    int m = 0;
    x(0,_).fillrandom();
    for (int iter=0; iter<maxiter; iter++) {
        print("\nITERATION", iter, m);
        f(m,_) = op(x(m,_));
        print("x");
        print(x(m,_));
        print(f(m,_));

        for (int j=0; j<=m; j++) {
            Q(j,m) = dot_product(x(j,_),f(m,_));
            Q(m,j) = dot_product(x(m,_),f(j,_));
        }

        Tensor<double> c = KAIN(Q(Slice(0,m),Slice(0,m)));
        print("THIS IS THE C I GOT");
        print(c);
        
        {
            m++;

            Tensor<double> xnew(n);
            for (int j=0; j<m; j++) {
                xnew(_) += c(j)*(x(j,_) - f(j,_));
            }

//             double steplen = (xnew-x(m-1,_)).normf();
//             double xnorm = xnew.normf();
//             if (steplen > 0.5*xnorm) {
//                 double damp = 0.3*xnorm/steplen;
//                 if (damp < 0.1) damp = 0.1;
//                 print("RESTRICTING", steplen, xnorm, damp);
//                 xnew = damp*xnew + (1.0-damp)*x(m-1,_);
//             }
            
            if (m == maxnvec) {
                for (int i=1; i<m; i++) {
                    f(i-1,_) = f(i,_);
                    x(i-1,_) = f(i,_);
                }
                Q(Slice(0,-2),Slice(0,-2)) = copy(Q(Slice(1,-1),Slice(1,-1)));
                
                m--;
            }
            x(m,_) = xnew;
        }
    }
    return 0;

}
