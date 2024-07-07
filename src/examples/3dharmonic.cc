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
  \file examples/3dharmonic.cc
  \brief Solves for the ground state of the quantum 3d harmonic oscillator
  \defgroup example3dharm Solves the 3D harmonic oscillator
  \ingroup examples

  The source is <a href=https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/3dharmonic.cc>here</a>.

  \par Points of interest
  - convolution with the Green's function
  - need to adjust the zero of energy to use the bound-state Green's function
  - failure of simple fixed-point iteration
  - use of simple non-linear equation solver
  - plotting 3D function along a line

  \par Background

  We seek the ground state of the 3D Schrodinger equation
  \f[
  \left( -\frac{1}{2} \nabla^2 + V(r) \right) \psi(r) = E \psi(r)
  \f]
  with
  \f[
  V(r) = \frac{1}{2} |r|^2
  \f]
  As usual, we rewrite the differential equation into integral form
  \f[
  \psi(r) = \left( -\frac{1}{2} \nabla^2 - E \right)^{-1} V(r) \psi(r)
  \f]
  but unfortunately we are left with two problems.

  First, recall that application
  of the inverse of the differential operator corresponds to convolution
  with the Green's function to the Helmholtz equation that satisfies
  \f[
  \left(-\nabla^2 + \mu^2 \right) G(r,r'; \mu) = \delta(r-r')
  \f]
  In 3D, we have
  \f[
  G(r,r'; \mu) = \frac{e^{-\mu |r-r'|}}{4 \pi |r-r|}
  \f]
  that MADNESS can currently only apply efficiently for real \f$\mu\f$ and since
  \f$\mu = \sqrt{-2 E}\f$ only for negative energies (hence bound
  states).  But for the harmonic oscillator there are no free particle
  states and the zero of energy is not chosen to describe the lowest
  energy of a free particle but simply as the zero of potential
  energy.  To solve this problem we can shift the zero of energy down
  by subtracting a constant (\f$\Delta\f$) from both sides of the
  equation, hence making the effective ground state energy negative.
  \f[
  \psi(r) = \left( -\frac{1}{2} \nabla^2 - E + \Delta \right)^{-1} \left( V(r) -\Delta\right)  \psi(r)
  \f]

  How negative do we need to make the energy?  To answer this we need
  to discuss the second problem.  The fixed-point iteration described
  by the integral equation only reliably converges to the ground state
  if the potential is negative everywhere the wave function is
  significant.  The exact solution is
  \f$\psi(r)=\pi^{-1/4}\exp(-r^2 / 2)\f$ (with $E=$1.5) that
  becomes 1e-6 (but how small is small enough?) at \f$r=5.3\f$ where
  \f$V\f$ is 14.0. So let's take this as the value of \f$\Delta\f$ and
  try the fixed point iteration.  Bad news.  It starts converging
  (slowly) to the right answer but then diverges and even damping
  (step restriction) does not solve the problem.  We have to make the
  shift large enough to make the potential negative in the entire
  volume to avoid the divergence, but this makes the convergence
  really slow.

  The fix is to not rely upon the simple fixed point iteration but to
  use an equation solver to force convergence.  This also enables us
  to choose the size of the shift to optimize the rate of convergence
  (empirically \f$\Delta=7\f$ is best) rather than being forced to
  pick a large value.  We use the very easy to use solver in
  mra/nonlinsol.h .

  [Aside.  It is possible to apply the operator for positive energies,
  but efficient application requires separate treatment of the
  singular and the long-range oscillatory terms, and the latter is
  presently not a production capability of MADNESS.  If you need this,
  let us know.]
*/

#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;

const double L = 7.0;
//const double DELTA = 3*L*L/2; // Use this to make fixed-point iteration converge
const double DELTA = 7.0;

// The initial guess wave function
double guess(const coord_3d& r) {
  return exp(-(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/3.0);
}

// The shifted potential
double potential(const coord_3d& r) {
  return 0.5*(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]) - DELTA;
}

// Convenience routine for plotting
void plot(const char* filename, const real_function_3d& f) {
  coord_3d lo(0.0), hi(0.0);
  lo[2] = -L; hi[2] = L;
  plot_line(filename,401,lo,hi,f);
}

double energy(World& world, const real_function_3d& phi, const real_function_3d& V) {
  double potential_energy = inner(phi,V*phi); // <phi|Vphi> = <phi|V|phi>
  double kinetic_energy = 0.0;
  for (int axis=0; axis<3; axis++) {
    real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    real_function_3d dphi = D(phi);
    kinetic_energy += 0.5*inner(dphi,dphi);  // (1/2) <dphi/dx | dphi/dx>
  }
  double energy = kinetic_energy + potential_energy;
  //print("kinetic",kinetic_energy,"potential", potential_energy, "total", energy);
  return energy;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());

    const double thresh = 1e-5;
    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    real_function_3d phi = real_factory_3d(world).f(guess);
    real_function_3d V = real_factory_3d(world).f(potential);
    plot("potential.dat", V);

    phi.scale(1.0/phi.norm2());  // phi *= 1.0/norm

    double E = energy(world,phi,V);

    NonlinearSolver solver;

    for (int iter=0; iter<100; iter++) {
      char filename[256];
      snprintf(filename, 256, "phi-%3.3d.dat", iter);
      plot(filename,phi);

      real_function_3d Vphi = V*phi;
      Vphi.truncate();
      real_convolution_3d op = BSHOperator3D(world, sqrt(-2*E), 0.01, thresh);

      real_function_3d r = phi + 2.0 * op(Vphi); // the residual
      double err = r.norm2();

      phi = solver.update(phi, r);
      //phi = phi-r;  // Replace the above line with this to use fixed-point iteration

      double norm = phi.norm2();
      phi.scale(1.0/norm);  // phi *= 1.0/norm
      E = energy(world,phi,V);

      if (world.rank() == 0)
          print("iteration", iter, "energy", E, "norm", norm, "error",err);

      if (err < 5e-4) break;
    }

    print("Final energy without shift", E+DELTA);

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
