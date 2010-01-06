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
#define WORLD_INSTANTIATE_STATIC_TEMPLATES  

/*!
  \file examples/laplace_sphere.cc
  \brief Solves Laplace's equation on the interior and exterior of a sphere
  \defgroup laplace_sphere Use of interior boundary conditions to solve Laplace's equation
  \ingroup examples

  \par Points of interest
  - Use of interior boundary conditions
  - Use of KAIN solver

  \par Background
  This example solves Laplace's equation on the interior and exterior of a sphere with
  Dirichlet boundary conditions.  Specifically, we solve 
  \f{eqnarray*}{
     \nabla^2 u(x) & = & 0 \\
     u(x) & = & \cos \theta \ \ |x|=1 
  \f}
  These simple boundary conditions are chosen to explore the accuracy
  of the solution since the exact solution is given by 
  \f{eqnarray*}{
  u(x) = \left\{  
     \begin{array}{1 1}
         |x| \cos \theta & \quad |x| \le 1 \\
         |x|^{-2} \cos \theta & \quad \mbox{otherwise}
     \end{array}
  \right.
  \f}

  For potential problems there are several ways to proceed, but we
  follow the more generally applicable diffuse domain approach of
  X. Li, J. Lowengrub, A. R&auml;tz, and A. Voight, "Solving PDEs in
  Complex Geometries: A Diffuse Domain Approach," Commun. Math. Sci.,
  7, p81-107, 2009 using approximation 2 (equation 2.22 in the paper).
  The surface is represented using a characteristic function \f$
  \phi(x) \f$ of half-width \f$ \epsilon \f$ computed here using the
  shape function library in mra/sdf_shape_3D.h .  A penalty-like term
  is introduced into the equation
  \f[
     \nabla^2 u - \epsilon^{-3} B(\phi) \left( u - g \right) = 0
  \f]
  where \f$ g(x) \f$ is the desired value of solution on the boundary
  (extended away from the boundary as a constant in the direction of
  the normal) and \f$ B(\phi) = 36 \phi^2 \left( 1 - \phi^2 \right) \f$.
  Employing the known free-space Green's function (\f$ G(x) = 1/4 \pi |x| \f$) 
  yields the working equation and expression for the residual
  \f[
      r(x) = u - G * \left( \epsilon^{-3} B(\phi) \left( u - g \right) \right) = 0
  \f]

  The initial guess is 
  

  \par Implementation

*/

#include <mra/mra.h>
#include <constants.h>
#include <mra/sdf_shape_3D.h>
#include <linalg/solvers.h>
using namespace madness;


// Returns a new functor combining two functors via multiplication.
// Look in mra/testsuite.cc for a more general version (BinaryOp)
class Product : public FunctionFunctorInterface<double,3> {
    real_functor_3d left, right;

public:
    Product(const real_functor_3d& left, const real_functor_3d& right)
        : left(left), right(right)
    {}

    double operator()(const coord_3d& x) const {
        return (*left)(x) * (*right)(x);
    }
};

// Computes cos(theta) (easier to combine with phi if we 
// use a functor rather than a function)
class CosTheta : public FunctionFunctorInterface<double,3> {
public:
    double operator()(const coord_3d& x) const {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        return x[2]/r;
    }
};

class Exact :  public FunctionFunctorInterface<double,3> {
public:
    double operator()(const coord_3d& x) const {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        double c = x[2]/r;
        
        if (r < 1.0) return c*r;
        return c/(r*r);
    }
};


int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);

  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  std::cout.precision(8);

  // Use defaults for numerical functions except for user simulation volume
  FunctionDefaults<3>::set_cubic_cell(-3,3);
  FunctionDefaults<3>::set_thresh(1e-4);
  FunctionDefaults<3>::set_k(6);

  coord_3d center;
  double epsilon = 0.05;

  // Make the characteristic function
  print("Making phi");
  real_functor_3d phi_functor(new SDF_Sphere<double>(epsilon, 1.0, center));
  real_function_3d phi = real_factory_3d(world).functor(phi_functor);  
  real_function_3d phic = 1.0 - phi;

  // Make B * epsilon^-3
  print("Making B");
  real_function_3d B = phi*phi*phic*phic*(36.0/(epsilon*epsilon*epsilon));
  B.truncate();
  print("B trace", B.trace());
  //plotdx(B, "B.dx");
  //plot_line("B.dat", 10001, coord_3d(-1.5), coord_3d(+1.5), B);

  // Make phi*g function and hence B*g*epsilon^-3
  print("Making phi * cos theta");
  real_functor_3d cos_functor(new CosTheta);
  real_functor_3d cos_phi_functor(new Product(cos_functor, phi_functor));
  real_function_3d Bg = real_factory_3d(world).functor(cos_phi_functor);
  plot_line("Bg.dat", 10001, coord_3d(-1.5), coord_3d(+1.5), Bg);
  print("Making B g");
  Bg = Bg*phi*phic*phic*(36.0/(epsilon*epsilon*epsilon));
  Bg.truncate();

  phi.clear(); /// Don't need these anymore
  phic.clear();

  // Make the Coulomb Green's function
  real_convolution_3d G = CoulombOperator<double>(world, FunctionDefaults<3>::get_k(), 
                                                  0.1*epsilon, FunctionDefaults<3>::get_thresh());
  // Initial guess for u is zero
  real_function_3d u = real_factory_3d(world);

  //u = apply(G,Bg);
  //u.scale(0.25/constants::pi/30);
  //u.truncate();

  // Iterate
  real_tensor Q;
  vector_real_function_3d ulist, rlist;
  for (int iter=0; iter<10; iter++) {
      // Compute the residual
      real_function_3d rhs = B*u - Bg;
      rhs.scale(-0.25/constants::pi);
      rhs.truncate();
      real_function_3d r = apply(G,rhs);
      r.truncate();
      r = r-u;
      ulist.push_back(u);
      rlist.push_back(r);

      // Solve subspace equations
      real_tensor Qnew(iter+1,iter+1);
      if (iter>0) Qnew(Slice(0,-2),Slice(0,-2)) = Q;
      for (int i=0; i<=iter; i++) {
          Qnew(i,iter) = inner(ulist[i],rlist[iter]);
          Qnew(iter,i) = inner(ulist[iter],rlist[i]);
      }
      Q = Qnew;
      real_tensor c = KAIN(Q);
      print("SOLUTION VECTOR", c);

      // Form new solution
      r = real_factory_3d(world);
      r.compress();
      for (int i=0; i<=iter; i++) {
          r.gaxpy(1.0,ulist[i], c[i]); 
          r.gaxpy(1.0,rlist[i],-c[i]); 
      }
      r.truncate();
      
      // Print/plot info
      print(iter, u.norm2(), r.norm2(), (u-r).norm2(), r.err(Exact()));
      u = r;
  }


  plotdx(u, "u.dx");
  plot_line("u.dat", 10001, coord_3d(-1.5), coord_3d(+1.5), u);

  finalize();

  return 0;
}

