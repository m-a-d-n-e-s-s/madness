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
*/

/**
 \file gaussian.cc
 \brief A simple example of projecting and using MADNESS functions.
 \ingroup gettingstarted

 This example projects a simple Gaussian function,
 \f[ g(\vec{x}) = e^{-|\vec{x}|^2}, \f]
 into the MADNESS adaptive basis and computes several integrals:
 \f[
    \int \mathrm{d}^3 \vec{x} \; g(\vec{x}) = \pi^{3/2} \doteq 5.5683279,
 \f]
 \f[
    \left( \int \mathrm{d}^3 \vec{x} \; g(\vec{x})^2 \right)^{1/2} = (\pi/2)^{3/4} \doteq 1.403104,
 \f]
 \f[
    \int \mathrm{d}^3 \vec{x} \; g(\vec{x}) \int \mathrm{d}^3 \vec{y} \; |\vec{x}-\vec{y}|^{-1} g(\vec{y}) = \pi^{5/2} 2^{1/2} \doteq 24.739429.
 \f]
*/

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
using namespace madness;

/// C++ function for projecting the Gaussian test function.

/// \param[in] r The coordinate at which to evaluate the function.
/// \return The value of the Gaussian at the specified coordinate.
double gaussian(const coord_3d& r) {
    double x=r[0], y=r[1], z=r[2];
    return exp(-(x*x + y*y + z*z));
}

/// Main function. Project the function and compute the integrals.

/// \param[in] argc The number of command-line arguments.
/// \param[in] argv The command-line arguments.
/// \return Exit status.
int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    FunctionDefaults<3>::set_cubic_cell(-6.0,6.0);
    real_convolution_3d op = CoulombOperator(world, 1e-4, 1e-6);
    real_function_3d g = real_factory_3d(world).f(gaussian);

    g.truncate();

    print(g.trace()); // exact trace is Pi^(3/2)=5.56832799683175
    print(g.norm2()); // exact norm2 is (Pi/2)^3/4=1.40310414553422
    print(inner(g,op(g))); // exact answer is Pi^(5/2)*sqrt(2) = 24.7394294511936

    finalize();
    return 0;
}
