/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <madness_config.h>
#include <world/world.h>

namespace madness {
    extern void load_quadrature(World& world, const char* dir);
    extern void legendre_polynomials(double x, long order, double *p);
    extern void legendre_scaling_functions(double x, long k, double *p);
    extern bool gauss_legendre(int n, double xlo, double xhi, double *x, double *w);
    extern bool gauss_legendre_numeric(int n, double xlo, double xhi, double *x, double *w);
    extern bool gauss_legendre_test(bool print=false);


    namespace detail {
        template <typename T>
        double norm(const T& t) {return std::abs(t);}

        template <typename T>
        double norm(const Tensor<T>& t) {return t.normf();}
    }

    template <typename funcT>
    typename funcT::returnT do_adq(double lo, double hi, const funcT& func, 
                                   int n, const double* x, const double* w) {
        // x and w provide the Gauss-Legendre quadrature rule of order n on [0,1]
        double range = (hi-lo);
        typename funcT::returnT sum = func(lo + range*x[0])*w[0];
        for (int i=1; i<n; i++) sum += func(lo + range*x[i])*w[i];
        return sum*range;
    }


    template <typename funcT>
    typename funcT::returnT adq1(double lo, double hi, const funcT& func, double thresh, 
                                 int n, const double* x, const double* w, int level) {
        static int MAX_LEVEL=10;
        if (level > MAX_LEVEL) throw "Adaptive quadrature: failed : runaway refinement?";
        double d = (hi-lo)/2;
        // Twoscale by any other name would smell just as sweet.
        typename funcT::returnT full = do_adq(lo, hi, func, n, x, w);
        typename funcT::returnT half = do_adq(lo, lo+d, func, n, x, w) + do_adq(lo+d, hi, func, n, x, w);

        double err = madness::detail::norm(full-half);
        //for (int i=0; i<level; i++) std::cout << "! ";
        //std::cout << err << std::endl;

        if (err > thresh) {
            return adq1(lo, lo+d, func, thresh*0.5, n, x, w, level+1) + 
                adq1(lo+d, hi, func, thresh*0.5, n, x, w, level+1);
        }
        else {
            return half;
        }
    }

    template <typename funcT>
    typename funcT::returnT adq(double lo, double hi, const funcT& func, double thresh) {
        const int n = 20;
        double x[n], y[n];
        gauss_legendre(n, 0.0, 1.0, x, y);

        return adq1(lo, hi, func, thresh, n, x, y, 0);
    }

    namespace detail {
        struct adqtest {
            typedef double returnT;
            double operator() (double x) const {
                // int(exp(-x^2)*cos(a*x),x=-inf..inf) = sqrt(pi)*exp(-a^2/4)
                return exp(-x*x)*cos(3*x);
            }

            static double exact() {
                const double pi = 3.1415926535897932384;
                return sqrt(pi)*exp(-9.0/4.0);
            }

            static bool runtest() {
                double test = madness::adq(-6.0,6.0,adqtest(),1e-14);
                return std::abs(test-exact())<1e-14;
            }

        };
    }
}

#endif

