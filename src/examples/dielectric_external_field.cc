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
  \file dielectric_external_field.cc
  \brief Example solution of Laplace's equations for dielectric sphere in an external field
  \defgroup exampledielectricfield Laplace's equations for dielectric sphere in an external field
  \ingroup examples

  The source is <a href=https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/dielectric_external_field.cc>here</a>.

  \par Points of interest
  - use of iterative equation solver
  - convolution with the Green's function
  - use of diffuse domain approximation to represent interior surfaces
  - line and volume plots
  - imposition of exterior Neumann boundary conditions

  \par Background

  We wish to solve Laplaces's equation within a non-homogeneous medium, i.e.,
  in which the permittivity is not constant, with an external electric field.
  \f[
  \nabla . \left( \epsilon(r) \nabla v(r)  \right) = 0
  \f]
  Expanding and rearranging yields,
  \f[
     \nabla^2 v(r) = - \frac{1}{\epsilon} \nabla \epsilon(r) .  \nabla v(r)
  \f]
  We can interpret \f$\nabla \epsilon . \nabla u / \epsilon\f$ as the
  induced surface charge density.

  MADNESS's default Coulomb Green's function (\f$-1 / |r-s|\f$) corresponds to the free-space
  Dirichlet boundary condition \f$v(\infty)=0\f$.  But we wish to
  impose the Neumann condition \f$- \partial v / \partial z(\infty) = E\f$.
  Thus, we write our solution as the sum of the asymptotic solution
  (\f$v = - E z\f$) and a component that goes to zero at infinity so that
  our MADNESS Green's function can be used.  Substituting
  \f[
  v(r) = - E z + u(r)
  \f]
  our equation becomes
  \f[
     \nabla^2 u(r) = \frac{1}{\epsilon} \nabla \epsilon(r) .  \left(E - \nabla u(r) \right)
  \f]
  The quantity in parentheses is the total electric field --- the sum of the applied
  field and that due to the induced surface charge.

  The exact solution is readily derived by expanding in spherical harmonics
  and matching the potential and electric displacement at the surface, or by
  looking inside Jackson, or by a quick google to find
  <A HREF="http://homer.phys.uh.edu/~hunger/class/fall_2010/lectures/lecture_10.pdf">
  this page.</A>

  (Note that while the Green's function to Poisson's equation is \f$-1 / 4 \pi |r-s|\f$,
  MADNESS provides convolution with \f$G(r,s) = 1/|r-s|\f$,
  so you have to keep track of the \f$-1/4\pi\f$ yourself.)

  (Note that to obtain an accurate solution with \f$\sigma<0.2\f$ you will need
  to increase the order of the wavelet and decrease the threshold, e.g.,
  \f$k=8, \mbox{thresh}=1e-5\f$.
*/


#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/tensor/solvers.h>
#include "molecularmask.h"
#include <madness/mra/nonlinsol.h>
#include <madness/constants.h>
#include <vector>

using namespace madness;
using namespace std;

const int k = 6; // wavelet order
const double thresh = 1e-4; // truncation threshold
const double L = 50; // box is [-L,L]
const double sigma = 0.3; // Surface width
const double Ez = 0.05; // External electric field

const double epsilon_0 = 100.0; // Interior dielectric
const double epsilon_1 =   1.0; // Exterior dielectric
const double R = 20.0; // Radius of cavity

// Adjustment to radius of sphere to enhance accuracy of calculation
// at finite value of sigma
const double delta = 0.65*sigma;

double exact_function(const coord_3d& x) {
    const double epsilon_r = epsilon_1 / epsilon_0;
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

    if (r < R) {
        return Ez*x[2] * (1 - 3 / (epsilon_r + 2));
    }
    else {
        return Ez*x[2] * (epsilon_r - 1)/(epsilon_r+2)*std::pow(R/r,3.0); // -1
    }
}

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    coord_3d lo(0.0), hi(0.0); // Range for line plotting
    lo[2] =-L;
    hi[2] = L;

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_truncate_on_project(true);
    FunctionDefaults<3>::set_bc(BC_FREE);

    // The Coulomb operator (this is just 1/r ... whereas the notes are -1/4pir)
    auto op = CoulombOperator(world, sigma*0.001, thresh*0.1);

    // Derivative operators
    auto Dx = free_space_derivative<double,3>(world, 0);
    auto Dy = free_space_derivative<double,3>(world, 1);
    auto Dz = free_space_derivative<double,3>(world, 2);

    // We will have one sphere of radius R centered at the origin
    vector<double> atomic_radii(1,R-delta);
    vector<coord_3d> atomic_coords(1,coord_3d(0.0));
    print("k     ", k);
    print("thresh", thresh);
    print("L     ", L);
    print("sigma ", sigma);
    print("delta ", delta);
    print("eps0  ", epsilon_0, "  eps1  ", epsilon_1);
    print("radii ", atomic_radii);
    print("coords", atomic_coords);

    print(MolecularVolumeExponentialSwitchLogGrad(sigma,epsilon_1,epsilon_0, atomic_radii, atomic_coords,0).special_points());

    // Log derivative of the dielectric function
    std::vector<real_function_3d> logd(3);
    for (int i = 0; i < 3; i++)
      logd.emplace_back(real_factory_3d(world).functor(MolecularVolumeExponentialSwitchLogGrad(sigma,epsilon_1,epsilon_0, atomic_radii, atomic_coords,i)));

    //double area = 4*madness::constants::pi*R*R;
    //double simulation_volume = 8*L*L*L;

    const double rfourpi = 1.0/(4.0*constants::pi);

    // Initial guess is zero
    real_function_3d u(world);
    double unorm = 0.0;

    NonlinearSolver solver(20);
    real_function_3d surface_charge, old_surface_charge(world);
    for (int iter=0; iter<20; iter++) {
        // Scale with 1/4pi AFTER applying operator to get one more digit of accuracy
        surface_charge = (logd[0]*Dx(u) + logd[1]*Dy(u) + logd[2]*(-Ez+Dz(u))).truncate();
        real_function_3d r = (u - op(surface_charge).scale(rfourpi)).truncate(thresh*0.032);
        surface_charge.scale(rfourpi);

        real_function_3d unew = solver.update(u, r);

        double change = (unew-u).norm2();
        double schange = (old_surface_charge-surface_charge).norm2();
        double err = u.err(exact_function);
        print("iter", iter, "change", change, "err", err,
              "surface charge", surface_charge.trace(), "schange", schange);

        old_surface_charge = surface_charge;

        if (change > 0.3*unorm)
            u = 0.5*unew + 0.5*u;
        else
            u = unew;

        if (schange < 10.0*thresh) break;
    }

    plot_line("testpot.dat", 1001, lo, hi, u, surface_charge);

    real_tensor cell(3,2);
    cell(_,0) = -L;
    cell(_,1) =  L;

    plotdx(u, "testpot.dx", cell);

    finalize();
    return 0;
}
