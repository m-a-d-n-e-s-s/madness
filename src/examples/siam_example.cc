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
 \file siam_example.cc
 \brief A simple example of solving PDEs with MADNESS.
 \ingroup gettingstarted

 This example solves Poisson's equation for a system with an inhomogeneous dielectric constant,
 \f[ \nabla \cdot{} \left[ \epsilon(\vec{r}) \nabla u(\vec{r}) \right] = - 4 \pi \rho(\vec{r}). \f]

 \f$ \epsilon(\vec{r}) = \epsilon_\mathrm{int}\f$ for \f$ |\vec{r}| < R\f$ and
 \f$ \epsilon(\vec{r}) = \epsilon_\mathrm{ext}\f$ for \f$ |\vec{r}| > R.\f$
 The PDE is converted to a fixed point iteration using the Green's function for
 the Laplacian (with free-space boundary conditions) and then solved using a
 fixed point iteration.

 Other objectives are to demonstrate the use of functors when projecting MADNESS
 functions.
*/

#include <madness/mra/mra.h>
#include <madness/constants.h>
#include "nonlinsol.h"

using namespace madness;
using namespace std;

const double eps_int = 1.0; ///< Interior dielectric.
const double eps_ext =10.0; ///< Exterior dielectric.
const double R = 2.0; ///< Radius of sphere.
const double xi = 100.0; ///< Exponent for delta function approx.
const double sigma = 0.1; ///< Surface "width".

/// Functor representing the log-derivative of \f$\epsilon\f$ along an axis.

/// A "mask" function is used to smoothly transition from inside the sphere
/// to outside the sphere. The width of this transition region is
/// \f$\mathcal{O}(\sigma)\f$.
class DSphere : public FunctionFunctorInterface<double, 3> {
  const double fac; ///< The factor multiplying the log-derivative.
  const int axis; ///< The component of the log-derivative we want.

protected:
  /// Derivative of the mask with respect to \f$s=sdf/\sigma\f$.

  /// \param[in] s The signed distance (in units of \f$\sigma\f$) from the
  ///    surface of the sphere.
  /// \return The derivative of the mask.
  double dmask(double s) const {
    if (fabs(s) > 6.0) return 0.0;
    return -exp(-s*s) / sqrt(constants::pi);
  }
  
public:
  /// Constructor.

  /// \param[in] Vint The dielectric inside the sphere.
  /// \param[in] Vext The dielectric outside the sphere.
  /// \param[in] axis The component of the log-derivative to project with this functor.
  DSphere(double Vint, double Vext, int axis) 
    : fac(log(Vint/Vext)), axis(axis)  {}
  
  /// Function used by MADNESS to project a function.

  /// \param[in] r The coordinate at which to evaluate the function.
  /// \return The value of the component of the log-derivative at the point.
  double operator()(const coord_3d &r) const {
    double d = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double sdf = (d-R) / sigma;
    return fac*r[axis]*dmask(sdf)/(sigma*d);
  }
};

/// Function that provides the charge density, \f$ \rho \f$.

/// \param[in] r The coordinate at which to evaluate the function.
/// \return The value of \f$ \rho \f$ at the point.
double rho_function(const coord_3d &r) {
  const double coeff = pow(xi / constants::pi, 1.5);
  return coeff*exp(-xi * (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

/// Main function. Project the functions and solve the PDE.

/// \param[in] argc The number of command-line arguments.
/// \param[in] argv The command-line arguments.
/// \return Exit status.
int main(int argc, char **argv) {
  const double rfourpi = 1.0/(4.0*constants::pi);
  initialize(argc, argv);
  World world(SafeMPI::COMM_WORLD);
  startup(world, argc, argv);

  // Function defaults
  FunctionDefaults<3>::set_k(6);
  FunctionDefaults<3>::set_thresh(1e-4);
  FunctionDefaults<3>::set_cubic_cell(-5, 5);
  FunctionDefaults<3>::set_initial_level(3);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // Make integral and derivative operators
  real_convolution_3d op = CoulombOperator(world, 1e-4, FunctionDefaults<3>::get_thresh());
  real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
  real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
  real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);

  // Make functors for dielectric related quantities
  real_functor_3d epsx_functor(new DSphere(eps_int, eps_ext, 0));
  real_functor_3d epsy_functor(new DSphere(eps_int, eps_ext, 1));
  real_functor_3d epsz_functor(new DSphere(eps_int, eps_ext, 2));

  // Make the actual numerical functions
  real_function_3d rho = real_factory_3d(world).f(rho_function);
  real_function_3d eps_x = real_factory_3d(world).functor(epsx_functor);
  real_function_3d eps_y = real_factory_3d(world).functor(epsy_functor);
  real_function_3d eps_z = real_factory_3d(world).functor(epsz_functor);

  real_function_3d u = op(rho).truncate(); // Initial guess
  NonlinearSolver solver;
  for (int iter=0; iter<20; iter++) {
    real_function_3d surf_rho = rfourpi*(eps_x*Dx(u) + eps_y*Dy(u) + eps_z*Dz(u));
    real_function_3d r = (u - op(rho + surf_rho)).truncate(); // residual
    real_function_3d unew = solver.update(u, r);
    double change = (unew - u).norm2();
    u = unew;

    print("iter", iter, "change", change, "surf charge", surf_rho.trace());
    
    if (change < 10.*FunctionDefaults<3>::get_thresh()) break;
  }
  
  coord_3d lo{0.0,0.0,0.0}, hi{5.0,0.0,0.0}; // Range for line plotting
  plot_line("testpot.dat", 501, lo, hi, u);
  
  finalize();
  return 0;
}
