#ifndef MADNESS_LINALG_SOLVERS_H__INCLUDED
#define MADNESS_LINALG_SOLVERS_H__INCLUDED

#include <tensor/tensor.h>
#include <world/print.h>
#include <iostream>
#include <linalg/tensor_lapack.h>

/*!
  \file solvers.h
  \brief Defines interfaces for optimization and non-linear equation solvers
  \ingroup solvers
*/

namespace madness {

    /*!
      Solves the KAIN equations for coefficients to compute the next vector
      
      \ingroup solvers

      \verbatim
      Wish to solve f(x)=0 with x and f(x) vectors 1..n
      
      Define (with i,j referring to the Krylov subspace and m the current iteration):
      
      Q(i,j) = <xi|fj>
      A(i,j) = <xi-xm | fj-fm> = Q(i,j) - Q(m,j) - Q(i,m) + Q(m,m)
      b(i) =-<xi-xm | fm> = -Q(i,m) + Q(m,m)
      
      We solve the equation
      
      A c = b
      
      and the correction to vector m
      
      . C
      .   interior = sum(i<m)[ c(i)*(x(i)-x(m)) ] = sum(i<m)[c(i)*x(i)] - x[m]*sum(i<m)[c(i)]
      .   exterior = -f(m) - sum(i<m)[ c(i)*(f(i)-f(m)) ] = -f(m) - sum(i<m)[c(i)*f(i)] + f(m)*sum(i<m)[c(i)]
      . New vector
      .   define C = sum(i<m)(c(i))  (note less than)
      .   define c(m) = 1.0 - C
      .   xnew = sum(i<=m) [ c(i)*(x(i) - f(i)) ]
      \endverbatim
    */
    template <typename T>
    Tensor<T> KAIN(const Tensor<T>& Q, double rcond=1e-12) {
        const int nvec = Q.dim(0);
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

    //     print("Q");
    //     print(Q);
    //     print("A");
    //     print(A);
    //     print("b");
    //     print(b);

        Tensor<T> x;
        Tensor<double> s, sumsq;
        long rank;
        gelss(A, b, rcond, x, s, rank, sumsq);
//         print("singular values", s);
//         print("rank", rank);
//         print("solution", x);

        Tensor<T> c(nvec);
        T sumC = 0.0;
        for (long i=0; i<m; i++) sumC += x(i);
        c(Slice(0,m-1)) = x;
//         print("SUMC", nvec, m, sumC);
        c(m) = 1.0 - sumC;

//         print("returned C", c);

        return c;
    }


    /// The interface to be provided by targets for non-linear equation solver

    /// \ingroup solvers
    struct SolverTargetInterface {
        /// Should return true if the Jacobian is implemented
        virtual bool provides_jacobian() const = 0;

        /// Should return the resdiual (vector F(x)) 
        virtual Tensor<double> residual(const Tensor<double>& x) = 0;

        /// Some solvers require the jacobian or are faster if an analytic form is available
        
        /// J(i,j) = partial F[i] over partial x[j] where F(x) is the vector valued residual
        virtual Tensor<double> jacobian(const Tensor<double>& x) {
            throw "not implemented";
        }

        /// Implement this if advantageous to compute residual and jacobian simultaneously
        virtual void residual_and_jacobian(const Tensor<double>& x, 
                                           Tensor<double>& residual, Tensor<double>& jacobian) {
            residual = this->residual(x);
            jacobian = this->jacobian(x);
        }

        virtual ~SolverTargetInterface() {}
    };



    /// The interface to be provided by functions to be optimized

    /// \ingroup solvers
    struct OptimizationTargetInterface {

        /// Should return true if the gradient is implemented
        virtual bool provides_gradient() const = 0;

        /// Should return the value of the objective function
        virtual double value(const Tensor<double>& x) = 0;

        /// Should return the derivative of the function
        virtual Tensor<double> gradient(const Tensor<double>& x) {
            throw "not implemented";
        }

        /// Reimplement if more efficient to evaluate both value and gradient in one call
        virtual void value_and_gradient(const Tensor<double>& x,
                                        double& value, 
                                        Tensor<double>& gradient) {
            value = this->value(x);
            gradient = this->gradient(x);
        }


        /// Numerical test of the derivative ... optionally prints to stdout, returns max abs error
        double test_gradient(Tensor<double>& x, double value_precision, bool doprint=true);

	virtual ~OptimizationTargetInterface(){}
    };


    /// The interface to be provided by solvers

    /// \ingroup solvers
    struct SolverInterface {
        virtual bool solve(Tensor<double>& x) = 0;
        virtual bool converged() const = 0;
        virtual double residual_norm() const = 0;
	virtual ~SolverInterface() {}
    };

    /// The interface to be provided by optimizers

    /// \ingroup solvers
    struct OptimizerInterface {
        virtual bool optimize(Tensor<double>& x) = 0;
        virtual bool converged() const = 0;
        virtual double value() const = 0;
        virtual double gradient_norm() const = 0;
	virtual ~OptimizerInterface(){}
    };


    /// Optimization via steepest descent

    /// \ingroup solvers
    class SteepestDescent : public OptimizerInterface {
        SharedPtr<OptimizationTargetInterface> target;
        const double tol;
        const double value_precision;  // Numerical precision of value
        const double gradient_precision; // Numerical precision of each element of residual
        double f;
        double gnorm;

    public:
        SteepestDescent(const SharedPtr<OptimizationTargetInterface>& target,
                        double tol = 1e-6,
                        double value_precision = 1e-12,
                        double gradient_precision = 1e-12);

        bool optimize(Tensor<double>& x);

        bool converged() const;

        double gradient_norm() const;

        double value() const;

	virtual ~SteepestDescent(){}
    };


    /// Optimization via quasi-Newton (BFGS or SR1 update)

    /// \ingroup solvers
    /// This is presently not 
    class QuasiNewton : public OptimizerInterface {
    private:
        std::string update;              // One of BFGS or SR1
        SharedPtr<OptimizationTargetInterface> target;
        const double tol;
        const double value_precision;  // Numerical precision of value
        const double gradient_precision; // Numerical precision of each element of residual
        double f;
        double gnorm;
        Tensor<double> h;
        int n;

        double line_search(double a1, double f0, double dxgrad, const Tensor<double>& x, const Tensor<double>& dx);

        void hessian_update_sr1(const Tensor<double>& s, const Tensor<double>& y);

        void hessian_update_bfgs(const Tensor<double>& dx, 
                                 const Tensor<double>& dg);

        Tensor<double> new_search_direction(const Tensor<double>& g);

    public:
        QuasiNewton(const SharedPtr<OptimizationTargetInterface>& target,
                    double tol = 1e-6,
                    double value_precision = 1e-12,
                    double gradient_precision = 1e-12);
        
        void set_update(const std::string& method);

        bool optimize(Tensor<double>& x);

        bool converged() const;

        double value() const;

        double gradient_norm() const;

	virtual ~QuasiNewton() {}
    };

}

#endif // MADNESS_LINALG_SOLVERS_H__INCLUDED
