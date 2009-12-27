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
// Solves the Navier-Stokes equations using a Chorin-Temam 
// pressure-correction scheme:
//   Step 1. A viscous prediction step
//        (1/dt) (u~(k+1) - u(k)) - nu Laplace(u(k+1) = f (k+1)
//
//   Step 2.  Projection
//        (1/dt) (u(k+1) - u~(k+1) ) + grad(phi(k+1)) = 0
//        div u(k+1) = 0


// Use solve_mode to select between operator- and iterative-based solvers.

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/mraimpl.h>
#include "poperator.h"
#include "mybicgstab.h"
#include "solutions.h"

using namespace madness;

const double deltaT = 0.01; // Size of time step
const int Nts = 100; // Number of time steps
const int k =  8; // Wavelet order (usually precision + 2)
const double pthresh = 1.e-6 ; // Precision
const double pthresh1 = 0.1*pthresh ;// * pthresh;
const double uthresh = pthresh; // Precision
const double uthresh1 = pthresh1;

const double L = 2 * WST_PI;
const double N = 8.0;

const int MAX_ITER = 1000;

// pick your preferred solution mode
const int SOLVE_MODE_OPER = 0;
const int SOLVE_MODE_ITER = 1;
const int solve_mode = SOLVE_MODE_OPER;

// shall we refine?
const bool REFINE_MODE = false;
const int  MAX_REF_LEVEL = 7;

const bool PLOT_OUTPUT = true;
//----------------------------------------------------------------------------

Tensor<int> bc(3, 2), bc0(3, 2);
	
World *pworld;
#define myfun vector< Function<T,NDIM> >

template<typename T, int NDIM> inline myfun operator-(const myfun& l, const myfun& r) { return sub(*pworld, l, r); }

// A simple laplacian operator
template<typename T, int n>
class LaplaceBiCG: public KryFun<T> {
  inline void operator() (const T& x, T& Ax) {
    FunctionDefaults<3>::set_refine(false);

    Ax = diff( diff(x,0), 0)
       + diff( diff(x,1), 1)
       + diff( diff(x,2), 2);

    Ax.truncate(pthresh1) ;
    return;
  }
};

template<typename T, int n>
class ChorinBiCG: public KryFun<T> {
  inline void operator() (const T& x, T& Ax) {
    FunctionDefaults<3>::set_refine(false);


    Ax[0] = 1./deltaT/mu*x[0];
    Ax[0] -= diff( diff(x[0],0), 0)
           + diff( diff(x[0],1), 1)
           + diff( diff(x[0],2), 2);

    Ax[1] = 1./deltaT/mu*x[1];
    Ax[1] -= diff( diff(x[1],0), 0)
           + diff( diff(x[1],1), 1)
           + diff( diff(x[1],2), 2);

    Ax[2] = 1./deltaT/mu*x[2];
    Ax[2] -= diff( diff(x[2],0), 0)
           + diff( diff(x[2],1), 1)
           + diff( diff(x[2],2), 2);

    Ax[0].truncate(pthresh1) ;
    Ax[1].truncate(pthresh1) ;
    Ax[2].truncate(pthresh1) ;
    return; 
  }
};
//--------------------------------------------------------------------

void testNavierStokes(World& world) {
  pworld = &world;

  // Function defaults
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(0.0, Lx);
  FunctionDefaults<3>::set_thresh(pthresh);

  FunctionDefaults<3>::set_autorefine(REFINE_MODE);
  FunctionDefaults<3>::set_refine(REFINE_MODE);
  FunctionDefaults<3>::set_max_refine_level(MAX_REF_LEVEL);

  bc = 1;  // periodic
  bc0 = 0;

  // create necessary operators
  FunctionDefaults<3>::set_bc(bc0);
  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
  SeparatedConvolution<double, 3> op = PeriodicCoulombOp<double, 3> (world,
 		k, pthresh1, pthresh1, cellsize);

  double const dum = 1 / deltaT / mu;
  SeparatedConvolution<double, 3> op1 = PeriodicBSHOp<double, 3> (world,
                sqrt(dum), k, uthresh1, uthresh1, cellsize);
  FunctionDefaults<3>::set_bc(bc);

  // Initialize t=0 solution
  mytime = 0.0;

  functT u(3);
  functT rhs(3);
  functT f(3);
  functT ue(3);
	
  u[0] = FunctionFactory<double, 3> (world).f(uxexact);
  u[1] = FunctionFactory<double, 3> (world).f(uyexact);
  u[2] = FunctionFactory<double, 3> (world).f(uzexact);

  ue[0] = FunctionFactory<double, 3> (world).f(uxexact);
  ue[1] = FunctionFactory<double, 3> (world).f(uyexact);
  ue[2] = FunctionFactory<double, 3> (world).f(uzexact);

  Function<double, 3> divf = FunctionFactory<double, 3> (world).f(init_zero);
  Function<double, 3> p    = FunctionFactory<double, 3> (world).f(pexact);

  BiCG0<functionT, double> BiSolver;
  BiCG1<functT, double>    BiSolver3;
  BiSolver3.setn(3);


  for (int t =  0; t < Nts; t++) {
    mytime = deltaT*(t+1);
    // Step 1.  Calculate the pressure at time t+1 explicitly 
    //          using uu (current u + last delta).
    //          Laplace p = div (f - uu \cdot grad(uu) ) 
    f[0] = FunctionFactory<double, 3> (world).f(fxexact);
    f[1] = FunctionFactory<double, 3> (world).f(fyexact);
    f[2] = FunctionFactory<double, 3> (world).f(fzexact);

    for (int i=0; i < 3; ++i) 
      rhs[i] = u[0]*diff(u[i],0) + u[1]*diff(u[i],1) + u[2]*diff(u[i],2);
	
    divf = div(f-rhs);

    divf.truncate();
    if (solve_mode == SOLVE_MODE_OPER) {
      FunctionDefaults<3>::set_bc(bc0);
      divf.set_bc(bc0);
      p.set_bc(bc0);
      p = apply(op, divf);
      p.scale(-1. / (4. * WST_PI)).set_bc(bc);
      divf.set_bc(bc);
      p.set_bc(bc);
      FunctionDefaults<3>::set_bc(bc);
    }
    else {
      int max_iter = MAX_ITER;
      double tol = pthresh;
      LaplaceBiCG<functionT,1> LapBiCG;
      int flag = BiSolver(LapBiCG, p, divf, max_iter, tol);
      if (world.rank() == 0)
        printf("  out flag = %d, error = %e, steps = %d \n", 
               flag, tol, max_iter);
    }
    divf.truncate();
    p.truncate();

    // Step 2.  Calculate the velocity at time t+1.
    //           (1/(deltaT mu) - Laplace) u_t+1 = 
    //             (f - grad p)/mu + u_t/(deltaT mu)

    rhs[0] = (f[0] - diff(p, 0) -rhs[0])*(1. / mu) + u[0]*dum;
    rhs[1] = (f[1] - diff(p, 1) -rhs[1])*(1. / mu) + u[1]*dum;
    rhs[2] = (f[2] - diff(p, 2) -rhs[2])*(1. / mu) + u[2]*dum;

    for (int i = 0; i < 3; ++i) {
      rhs[i].truncate();
      ue[i].truncate();
    }

    if (solve_mode == SOLVE_MODE_OPER) {
      FunctionDefaults<3>::set_bc(bc0);
      for (int i = 0; i < 3; ++i) {
        rhs[i].set_bc(bc0);
        ue[i].set_bc(bc0);
        ue[i] = apply(op1, rhs[i]);
	ue[i].set_bc(bc);
	rhs[i].set_bc(bc);
      }
      FunctionDefaults<3>::set_bc(bc);
    }
    else {
      int max_iter = MAX_ITER;
      double tol = pthresh;
      ChorinBiCG<functT,1> ChorBiCG;
      int flag = BiSolver3(ChorBiCG, ue, rhs, max_iter, tol);
      if (world.rank() == 0)
        printf("  out flag = %d, error = %e, steps = %d \n", 
               flag, tol, max_iter);
    }

    for (int i = 1; i < 3; ++i) {
      rhs[i].truncate();
      ue[i].truncate();
    }
    
    // update velocity
    for (int i=0; i < 3; ++i) { 
      u[i] = 2.0*ue[i] - u[i];
      u[i].truncate();
      ue[i].truncate();
    }
    ++t; 
    mytime += deltaT;

    // print norms, etc
    Function<double, 3> du = FunctionFactory<double, 3> (world).f(uxexact);
    Function<double, 3> dv = FunctionFactory<double, 3> (world).f(uyexact);
    Function<double, 3> dw = FunctionFactory<double, 3> (world).f(uzexact); 
    Function<double, 3> dp = FunctionFactory<double, 3> (world).f(pexact); 
    Function<double, 3> pe = FunctionFactory<double, 3> (world).f(pexact); 
    du -= u[0];
    dv -= u[1];
    dw -= u[2];
    dp -= p;
    pe.truncate();
    double a=div(u).norm2(), b=du.norm2(), c=dv.norm2(), 
           d=dw.norm2(),e=dp.norm2();
    if (world.rank()==0)  print("Norm", t+1, mytime, a,b,c,d,e);
			
    // Output the solution
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> numpt;
    for(int i = 0; i < 3; ++i) {
      plotlo[i] = 0;
      plothi[i] = L;
      numpt[i] = 20;
    }
    bool binary = false;

    char filename[100];
    if (solve_mode == SOLVE_MODE_OPER)  {
      sprintf(filename, "data-op%02d.vts", t);
    }
    else {
      sprintf(filename, "data-it%02d.vts", t);
    }

    if (PLOT_OUTPUT == true) {		
      plotvtk_begin(world, filename, plotlo, plothi, numpt,binary);
      plotvtk_data(u[0], "u", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(u[1], "v", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(u[2], "w", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(ue[0], "ue", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(ue[0], "uer", world, filename, plotlo, plothi, 
                   numpt, binary, true);
      plotvtk_data(ue[1], "ve", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(ue[2], "we", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(p, "p", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_data(pe, "pe", world, filename, plotlo, plothi, numpt, binary);
      plotvtk_end<3>(world, filename,binary);
    }
  } // end evolution loop
} // end testNavierStokes
//*****************************************************************************

//*****************************************************************************
int main(int argc, char**argv) {
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);
  ThreadPool::begin(); 

  startup(world,argc,argv);

  try {
    testNavierStokes(world);
  } catch (const MPI::Exception& e) {
      //print(e); std::cout.flush();
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
        print(e); std::cout.flush();
        error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
        print(e); std::cout.flush();
        error("caught a Tensor exception");
  } catch (const char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
  } catch (char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
  } catch (const std::string& s) {
        print(s); std::cout.flush();
        error("caught a string (class) exception");
  } catch (const std::exception& e) {
        print(e.what()); std::cout.flush();
        error("caught an STL exception");
  } catch (...) {
        error("caught unhandled exception");
  }


  world.gop.fence();

  ThreadPool::end();
  print_stats(world);
  finalize();
  return 0;
}
//*****************************************************************************

