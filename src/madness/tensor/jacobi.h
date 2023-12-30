#ifndef MADNESS_JACOBI_H
#define MADNESS_JACOBI_H

#include <madness/world/madness_exception.h>
#include <madness/tensor/tensor.h>
#include <madness/world/print.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

namespace  madness {
void jacobi(Tensor<double>& A, Tensor<double>& V, const std::vector<int>& set) {
    int n = A.dim(0);
    MADNESS_CHECK(A.ndim() == 2);
    MADNESS_CHECK(A.dim(0) == A.dim(1));
    MADNESS_CHECK(V.ndim() == 2);
    MADNESS_CHECK(n == V.dim(0) && n == V.dim(1));
    MADNESS_CHECK(n == long(set.size()));

    int i, j, k, iter, nrot;
    double tolmin = std::sqrt(n) * 5.0e-16, tol = 0.05;
    double maxd;

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            A(i,j) = A(j,i) = 0.5*(A(i,j) + A(j,i)); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            V(i,j) = 0.0;
        }
        V(i,i) = 1.0;
    }

    maxd = 0.;
    for (i=0; i<n; i++) {
        double daii = std::fabs(A(i,i));
        maxd = std::max(maxd,daii);
    }

#define ROT(A,i,j)                \
    { \
      double * __restrict__ ai = &A(i,0); \
      double * __restrict__ aj = &A(j,0); \
      for (k=0; k<n; k++) { \
          double t = ai[k]; \
          double u = aj[k]; \
          ai[k] = c*t - s*u; \
          aj[k] = s*t + c*u; \
      } \
    }

#define ROTT(A,i,j)                \
    { \
      double * __restrict__ ai = &A(0,i);      \
      double * __restrict__ aj = &A(0,j);       \
      for (k=0; k<n*n; k+=n) { \
          double t = ai[k]; \
          double u = aj[k]; \
          ai[k] = c*t - s*u; \
          aj[k] = s*t + c*u; \
      } \
    }

    for (iter=0; iter<5000; iter++) {
        double maxdaij = 0.0;
        nrot = 0;
        for (i=0; i<n; i++) {
            for (j=i+1; j<n; j++) {
                if (set[i] != set[j]) {
                    double aii = A(i,i);
                    double ajj = A(j,j);
                    double aij = A(i,j);
                    double daij = std::fabs(aij);

                    maxdaij = std::max(maxdaij,daij/maxd);

                    if (daij > tol*maxd) {
                        double s = ajj - aii;
                        double ds = std::fabs(s);
                        if (daij > (tolmin*ds)) {
                            double c;
                            nrot++;
                            if ((tolmin*daij) > ds) {
                                c = s = 1.0/std::sqrt(2.);
                            }
                            else {
                                double t = aij/s;
                                double u = 0.25/std::sqrt(0.25+t*t);
                                c = std::sqrt(0.5+u);
                                s = 2.*t*u/c;
                            }

                            //ROTT(A,i,j);
                            for (k=0; k<n; k++) {
                                double t = A(k,i);
                                double u = A(k,j);
                                A(k,i) = c*t - s*u;
                                A(k,j) = s*t + c*u;
                            }

                            ROT(A,i,j);
                            // for (k=0; k<n; k++) {
                            //     double t = A(i,k);
                            //     double u = A(j,k);
                            //     A(i,k) = c*t - s*u;
                            //     A(j,k) = s*t + c*u;
                            // }

                            ROT(V,i,j);
                            // for (k=0; k<n; k++) {
                            //     double t = V(i,k);
                            //     double u = V(j,k);
                            //     V(i,k) = c*t - s*u;
                            //     V(j,k) = s*t + c*u;
                            // }

                            A(j,i) = A(i,j) = 0.0;

                            maxd = std::max(maxd,std::fabs(A(i,i)));
                            maxd = std::max(maxd,std::fabs(A(j,j)));
                        }
                    }
                }
            }
        }
        //printf("iter=%d nrot=%d err=%e\n", iter, nrot, maxdaij);
        if (nrot == 0 && tol <= tolmin) break;
        tol = std::min(tol,maxdaij*1e-1);
        //tol = std::min(tol,maxdaij*maxdaij); // is not quadratic if only block diagonalizing?
        tol = std::max(tol,tolmin);
    }
    if (nrot != 0)
        printf("Jacobi iteration did not converge in 5000 passes\n");
}
} // namespace madness

#endif //MADNESS_JACOBI_H
