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

/*!
  \file examples/sininteg.cc
  \brief Compute the integral sin(x) x=0..10
  \warning If you change this example, please update the "MADNESS Basics" module
  in the documentation (doc/getting_started/gstart_basics.dox). \defgroup
  sininteg First example from getting started guide \ingroup examples

  Computes the integral
  \f[
     \int_0^{10} sin(x) dx
  \f]
  by projecting \f$ sin(x) \f$ into the discontinuous spectral element
  basis and using the \c trace() method.
 */

// #warning If you changed this example, please update the "Getting started with
// MADNESS" documentation module.

#include <madness/mra/mra.h>

using namespace madness;

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

using namespace madness;

static const double L = 32.0;      // box size
static const long k = 8;           // wavelet order
static const double thresh = 1e-6; // precision

static const size_t nfunc = 64; // number of functions

// A class that behaves like a function to compute a Gaussian of given origin
// and exponent
class Gaussian : public FunctionFunctorInterface<double, 3> {
public:
  const coord_3d center;
  const double exponent;
  const double coefficient;
  std::vector<coord_3d> specialpt;

  Gaussian(const coord_3d &center, double exponent, double coefficient)
      : center(center), exponent(exponent), coefficient(coefficient),
        specialpt(1) {
    specialpt[0][0] = center[0];
    specialpt[0][1] = center[1];
    specialpt[0][2] = center[2];
  }

  // MADNESS will call this interface
  double operator()(const coord_3d &x) const {
    double sum = 0.0;
    for (int i = 0; i < 3; i++) {
      double xx = center[i] - x[i];
      sum += xx * xx;
    };
    return coefficient * exp(-exponent * sum);
  }

  // By default, adaptive projection into the spectral element basis
  // starts uniformly distributed at the initial level.  However, if
  // a function is "spiky" it may be necessary to project at a finer
  // level but doing this uniformly is expensive.  This method
  // enables us to tell MADNESS about points/areas needing deep
  // refinement (the default is no special points).
  std::vector<coord_3d> special_points() const { return specialpt; }
};

// Makes a new square-normalized Gaussian functor with random origin and
// exponent
real_functor_3d random_gaussian() {
  const double expntmin = 1e-1;
  const double expntmax = 1e4;
  const real_tensor &cell = FunctionDefaults<3>::get_cell();
  coord_3d origin;
  for (int i = 0; i < 3; i++) {
    origin[i] = RandomValue<double>() * (cell(i, 1) - cell(i, 0)) + cell(i, 0);
  }
  double lo = log(expntmin);
  double hi = log(expntmax);
  double expnt = exp(RandomValue<double>() * (hi - lo) + lo);
  print("expnt", expnt, origin);
  double coeff = pow(2.0 * expnt / constants::pi, 0.75);
  return real_functor_3d(new Gaussian(origin, expnt, coeff));
}

// Makes a vector of new square-normalized Gaussian functions with random origin
// and exponent
std::vector<real_function_3d> random_gaussians(size_t n, World &world) {
  std::vector<real_function_3d> result(n);
  for (size_t i = 0; i < n; i++) {
    result[i] = FunctionFactory<double, 3>(world).functor(random_gaussian());
  }
  return result;
}

void test(World &world) {

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L / 2, L / 2);

  default_random_generator.setstate(
      99); // Ensure all processes have the same state

  // Create a vector of random Gaussian functions
  std::vector<real_function_3d> a = random_gaussians(nfunc, world);
  truncate(world, a);
  compress(world, a, true);
  auto b = copy(world, a);
  reconstruct(world, b, true);
  compress(world, b, true);

  auto diff = sub(world, a, b, true);
  // compute the norm of the errors for each component
  for (size_t i = 0; i < nfunc; i++) {
    if (world.rank() == 0) {
      print("error", i, diff[i].norm2());
    }
  }
}

int main(int argc, char **argv) {
  initialize(argc, argv);
  World world(SafeMPI::COMM_WORLD);

  startup(world, argc, argv);

  if (world.rank() == 0)
    FunctionDefaults<3>::print();

  test(world);

  finalize();
  return 0;
}
