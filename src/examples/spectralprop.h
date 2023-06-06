/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
/*!
  \file spectralprop.h
  \brief Spectral propagator in time using semigroup approach
  \defgroup spectralprop Spectral propagator in time using semigroup approach
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/spectralprop.h>here</a>.

  \par Points of interest
  - high-order integration of general time-dependent problems
  - semigroup approach

  \par Background

  Given \f$ u(0) \f$,  we seek to solve the PDE
  \f[
     \frac{du}{dt} = \hat{L} u + N(u)
  \f]
  for the function at some future time \f$ t \f$ (i.e., for \f$ u(t) \f$).
  \f$ \hat{L} \f$ is a linear operator that we are
  able to exponentiate, and \f$ N \f$ is everything else
  including linear and non-linear parts.

  In the semigroup approach the formal solution to the PDE
  is written
  \f[
     u(t) = \exp(t \hat{L} ) u(0) + \int_0^t \exp((t-\tau)\hat{L}) N(u(\tau)) d\tau
  \f]
  Numerical quadrature of the integral using Gauss-Legendre
  quadrature points is used resulting in a set of equations
  that are iteratively solved (presently using simple fixed
  point iteration from a first-order explicit rule).

  The user provides
  - functors to apply the exponential and non-linear parts, and
  - if necessary a user-defined data type that supports
    a copy constructor, assignment, inplace addition, multiplication
    from the right by a double, and computation of the distance
    between two solutions \f$ a \f$ and \f$ b \f$ with the api
    \verb+double distance(a,b)+

  Have a look in testspectralprop.cc for example use.

  With \f$ n \f$ quadrature points, the error is \f$ O\left(t^{2n+1}\right) \f$
  and the number of applications of the exponential operator per
  time step is \f$ 1+(n_{it}+1)n +n_{it}n^2 \f$ where \f$ n_{it} \f$
  is the number of iterations necessary to solve the equations
  (typically about 5 but this is problem dependent).

*/


#include <madness/mra/legendre.h>
#include <madness/world/madness_exception.h>

#include <algorithm>
#include <cmath>
#include <complex>

namespace madness {

    /// Default function for computing the distance between two doubles
    inline double distance(double a, double b)
    {
        return std::sqrt((a-b)*(a-b));
    }

    /// Default function for computing the distance between two complex numbers
    template <typename T>
    inline double distance(std::complex<T>& a, std::complex<T>& b)
    {
        return std::abs(a-b);
    }


    /// Spectral propagtor in time.  Refer to documentation of file spectralprop.h for math detail.
    class SpectralPropagator {
        const int NPT;          ///< Number of quadrature points
        std::vector<double> x;  ///< Quadrature points on [0,1] with value 0 prepended
        std::vector<double> w;  ///< Quadrature weights on [0,1] with value 0 prepended
        SpectralPropagator* q;
        
        /// Private: Makes interpolating polyn p[i](t), i=0..NPT
        double p(int i, double t) {
            double top=1.0, bot=1.0;
            for (int j=0; j<i; j++) {
                top *= (t-x[j]);
                bot *= (x[i]-x[j]);
            }
            for (int j=i+1; j<=NPT; j++) {
                top *= (t-x[j]);
                bot *= (x[i]-x[j]);
            }
            return top/bot;
        }

        /// Private: Computes interpolated solution
        template <typename uT>
        uT u(double dt, const std::vector<uT>& v) {
            uT U = v[0]*p(0,dt);
            for (int i=1; i<=NPT; i++) {
                U += v[i]*p(i,dt);
            }
            return U;
        }

        // Separated this out to enable recursive use of solver ... but this 
        // turned out to be largely not useful.
        template <typename uT, typename expLT, typename NT>
        std::vector<uT> solve(double t, double Delta, const uT& u0, const expLT& expL, const NT& N, const double eps, bool doprint, bool recurinit, int& napp) {
            if (doprint) madness::print("solving with NPT", NPT);
            std::vector<uT> v(NPT+1,u0);

            if (!recurinit || NPT==1) {
                // Initialize using explicit first-order accurate rule
                for (int i=1; i<=NPT; i++) {
                    double dt = x[i] - x[i-1];
                    v[i] = expL(dt*Delta,v[i-1]);
                    v[i] += expL(Delta*dt,N(t+Delta*x[i-1],v[i-1]))*(dt*Delta);
                }
            }
            else {
                // Recur down to lower-order quadrature rules for initial guess
                std::vector<uT> vx = q->solve(t, Delta, u0, expL, N, eps*10000.0, doprint, recurinit, napp);
                for (int i=1; i<=NPT; i++) {
                    v[i] = q->u(x[i],vx);
                }
            }

            // Precompute G(x[1]-x[0])*v[0]
            uT Gv0 = expL((x[1] - x[0])*Delta,v[0]); napp++;

            // Crude fixed point iteration
            uT vold = v[NPT]; // Track convergence thru change at last quadrature point
            bool converged = false;
            for (int iter=0; iter<20; iter++) {
                for (int i=1; i<=NPT; i++) {
                    double dt = x[i] - x[i-1];
                    uT vinew = (i==0) ? Gv0*1.0 : expL(dt*Delta,v[i-1]); // *1.0 in case copy is shallow
                    if (i != 0) napp++;
                    for (int k=1; k<=NPT; k++) {
                        double ddt = x[i-1] + dt*x[k];
                        vinew += expL(dt*Delta*(1.0-x[k]), N(t + Delta*ddt, u(ddt, v)))*(dt*Delta*w[k]); napp++;
                    }
                    v[i] = vinew;
                }
                double err = ::madness::distance(v[NPT],vold);
                vold = v[NPT];
                if (doprint) print("spectral",iter,err);
                converged = (err < eps);
                if (converged) break;
            }
            if (!converged) throw "spectral propagator failed to converge";

            return v;
        }

    public:
        /// Construct propagator using \c NPT points
	SpectralPropagator(int NPT)
            : NPT(NPT)
            , x(NPT+1)
            , w(NPT+1)
            , q(0)
        {
            x[0] = w[0] = 0.0;
            MADNESS_ASSERT(gauss_legendre(NPT, 0.0, 1.0, &x[1], &w[1]));

            // The iterative solution requires quadrature points in increasing
            // order corresponding to shooting forward.
            std::reverse(x.begin()+1, x.end());
            std::reverse(w.begin()+1, w.end());

            if (NPT > 1) q = new SpectralPropagator(NPT-1);
        }

        ~SpectralPropagator() 
        {
            delete q;
        }

        /// Step forward in time from \f$ t \f$ to \f$ t+\Delta \f$

        /// The template types should be automatically inferred from
        /// the invocation.  \c uT is the C++ type for the solution.
        ///
        /// @param[in] t The current time
        /// @param[in] Delta The time step
        /// @param[in] u0 The solution at the current time
        /// @param[in] expL A function or functor to compute \f$ \exp{\tau \hat{L}} u \f$  for \f$ \tau \in (0,\Delta) \f$
        /// @param[in] N A function to compute the non-linear part
        /// @param[in] eps Threshold for convergence of iteration on norm of change in solution at last quadrature point.
        /// @param[in] doprint If true will print some info on convergence
        /// @param[in] recurinit If true initial guess is recursion from lower order quadrature instead of default explicit rule
        /// @returns The solution at time \f$ t+\Delta \f$
        ///
        /// The user provided operators are invoked as
        /// \code
        /// uT expL(double tau, const uT& u)
        /// \endcode
        /// where \f$ \tau \in (0,\Delta) \f$
        /// and
        /// \code
        /// uT N(double T, const uT& u)
        /// \endcode
        /// where \f$ T \in (t,t+\Delta) \f$ and \f$ u \f$ is a trial solution at time \f$ T \f$.
        template <typename uT, typename expLT, typename NT>
        uT step(double t, double Delta, const uT& u0, const expLT& expL, const NT& N, const double eps=1e-12, bool doprint=false, bool recurinit=true) {
            int napp = 0;
            
            // Compute solution at quadrature points
            std::vector<uT> v = solve(t, Delta, u0, expL, N, eps, doprint, recurinit, napp);

            // Integrate forward from last quadrature point to the end point
            double dt = 1.0 - x[NPT];
            uT vinew = expL(dt*Delta,v[NPT]);
            for (int k=1; k<=NPT; k++) {
                double ddt = x[NPT] + dt*x[k];
                vinew += expL(dt*Delta*(1.0-x[k]), N(t + Delta*ddt, u(ddt, v)))*(dt*Delta*w[k]); napp++;
            }

            if (doprint) print("number of operator applications", napp);
            return vinew;
	}
    };

    class SpectralPropagatorGaussLobatto {
        const int NPT;
        std::vector<double> x;
        std::vector<double> w;

        // Makes interpolating polyn p[i](t), i=0..NPT-1
        double p(int i, double t) {
            double top=1.0, bot=1.0;
            for (int j=0; j<i; j++) {
                top *= (t-x[j]);
                bot *= (x[i]-x[j]);
            }
            for (int j=i+1; j<NPT; j++) {
                top *= (t-x[j]);
                bot *= (x[i]-x[j]);
            }
            return top/bot;
        }

        template <typename uT>
        uT u(double dt, const std::vector<uT>& v) {
            uT U = v[0]*p(0,dt);
            for (int i=1; i<NPT; i++) {
                U += v[i]*p(i,dt);
            }
            return U;
        }

    public:
        // Constructor
        // ... input ... npt,
        // ... makes ... x, w, p
        //
        // Apply
        // ... input ... u(0)
        // ... makes guess v(i), i=0..npt using explicit 1st order rule.

	SpectralPropagatorGaussLobatto(int NPT)
            : NPT(NPT)
            , x(NPT)
            , w(NPT)
        {
            MADNESS_ASSERT(NPT>=2 && NPT<=6);
            // Gauss Lobatto on [-1,1]
            x[0] = -1.0; x[NPT-1] = 1.0; w[0] = w[NPT-1] = 2.0/(NPT*(NPT-1));
            if (NPT == 2) {
            }
            else if (NPT == 3) {
                x[1] = 0.0; w[1] = 4.0/3.0;
            }
            else if (NPT == 4) {
                x[1] = -1.0/std::sqrt(5.0);  w[1] = 5.0/6.0;
                x[2] = -x[1];                w[2] = w[1];
            }
            else if (NPT == 5) {
                x[1] = -std::sqrt(3.0/7.0); w[1] = 49.0/90.0;
                x[2] = 0.0;                  w[2] = 32.0/45.0;
                x[3] = -x[1];                w[3] = w[1];
            }
            else if (NPT == 6) {
                x[1] = -std::sqrt((7.0+2.0*std::sqrt(7.0))/21.0); w[1] = (14.0-std::sqrt(7.0))/30.0;
                x[2] = -std::sqrt((7.0-2.0*std::sqrt(7.0))/21.0); w[2] = (14.0+std::sqrt(7.0))/30.0;
                x[3] = -x[2];                                     w[3] = w[2];
                x[4] = -x[1];                                     w[4] = w[1];
            }
            else {
                throw "bad NPT";
            }

            // Map from [-1,1] onto [0,1]
            for (int i=0; i<NPT; i++) {
                x[i] = 0.5*(1.0 + x[i]);
                w[i] = 0.5*w[i];
            }
        }

        template <typename uT, typename expLT, typename NT>
        uT step(double t, double Delta, const uT& u0, const expLT& expL, const NT& N, const double eps=1e-12, bool doprint=false) {
            std::vector<uT> v(NPT,u0);

            // Initialize using explicit first-order accurate rule
            for (int i=1; i<NPT; i++) {
                double dt = x[i] - x[i-1];
                v[i] = expL(dt*Delta,v[i-1]);
                v[i] += expL(Delta*dt,N(t+Delta*x[i-1],v[i-1]))*(dt*Delta);
            }

            // Crude fixed point iteration
            uT vold = v[NPT-1];
            for (int iter=0; iter<12; iter++) {
                for (int i=1; i<NPT; i++) {
                    double dt = x[i] - x[i-1];
                    uT vinew = expL(dt*Delta,v[i-1]);
                    for (int k=0; k<NPT; k++) {
                        double ddt = x[i-1] + dt*x[k];
                        vinew += expL(dt*Delta*(1.0-x[k]), N(t + Delta*ddt, u(ddt, v)))*(dt*Delta*w[k]);
                    }
                    v[i] = vinew;
                }
                double err = ::madness::distance(vold, v[NPT-1]);
                if (doprint) print("hello",iter,err);
                vold = v[NPT-1];
                if (err < eps) break;
            }

            return vold;
	}
    };

}
