#ifndef MAD_GMRES_H
#define MAD_GMRES_H

#include <tensor/tensor.h>
#include <world/print.h>
#include <iostream>
#include <linalg/tensor_lapack.h>

/// \file gmres.h
/// \brief Defines a general operator interface and a templated GMRES solver
/// for solving linear equations.

namespace madness {

    /// A generic operator: takes in one T and produces another T.
    /// Over-ride the protected action() function to implement an Operator.
    template <typename T>
    class Operator {
        protected:
            /// The action of the operator
            virtual void action(const T &in, T &out) const = 0;

        public:
            /// public access to the operator's action, returns out for convenience
            T & applyOp(const T &in, T &out) const {
                action(in, out);
                return out;
            }
    };

    /// A GMRES solver routine for linear systems, Ax == b.
    /// Requires the operator, A, the inhomogeneity b, an initial guess x,
    /// the maximum number of iterations, the convergence threshold, and a
    /// flag for producing output.
    ///
    /// The output, if desired, shows iteration number, the residual, and the
    /// effective rank of the GMRES matrix at that iteration.
    ///
    /// Returns the number of iterations taken.
    ///
    /// NOTE: This function assumes T has gaxpy, norm2, and scale member
    /// functions.  Specialization may be needed if this is not the case.
    template <typename T, typename fp>
    int GMRES(const Operator<T> &op, const T &bin, T &x, const int maxiters,
        const fp thresh, const bool outp = true) {

        int iter, i;
        long rank;
        std::vector<T> V;
        T r;
        Tensor<fp> H(maxiters+1, maxiters);
        Tensor<fp> betae(maxiters+1);
        Tensor<fp> y;
        Tensor<typename Tensor<fp>::scalar_type> s;
        Tensor<typename Tensor<fp>::scalar_type> sumsq;
        fp resid;

        // initialize
        H = 0.0;
        betae = 0.0;

        // construct the first subspace basis vector
        iter = 0;
        op.applyOp(x, r);
        r.gaxpy(-1.0, bin, 1.0);
        betae[0] = r.norm2();
        if(outp)
            printf("itr rnk resid\n%.3d N/A %.6e\n", iter, betae[0]);
        if(betae[0] < thresh)
            return 0;
        r.scale(1.0 / betae[0]);
        V.push_back(r);
        ++iter;

        do {
            // compute the new vector
            op.applyOp(V[iter - 1], r);

            // orthogonalize the new vector
            for(i = 0; i < iter; ++i) {
                H(i, iter-1) = inner(r, V[i]);
                r.gaxpy(1.0, V[i], -H(i, iter-1));
            }
            H(iter, iter-1) = r.norm2();

            // solve Hy == betae for y
            gelss(H(Slice(0, iter), Slice(0, iter-1)), betae(Slice(0, iter)),
                1.0e-12, &y, &s, &rank, &sumsq);

            resid = sumsq[0];
            if(outp) {
                printf("%.3d %.3ld %.6e", iter, rank, resid);
                if(iter != rank)
                    printf(" ** Questionable Progress **");
                printf("\n");
            }

            r.scale(1.0 / H(iter, iter-1));
            V.push_back(r);
            ++iter;
       } while(iter < maxiters && resid > thresh);

       // build the solution and destroy the basis vectors
       for(i = 0; i < iter; ++i) {
           x.gaxpy(1.0, V[i], y[i]);
           V[i].clear();
       }

       return iter;
    }

}

#endif
