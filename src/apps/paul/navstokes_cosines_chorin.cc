// Solves the Navier-Stokes equations using a Chorin pressure-correction scheme
// Evolves two solutions: one using an operator method, and another using 
//   an iterative method

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "poperator.h"
#include <linalg/gmres.h>
#include "mybicgstab.h"

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Vector<double,1> coordT1d;
typedef Function<double,3> functionT;
typedef std::vector<functionT> functT ;

// useful constants
const int XDIR = 0;
const int YDIR = 1;
const int ZDIR = 2;
const int MAX_ITER = 10000; // for iterative solvers

// problem and resolution definitions
const double mu     = 1.     ; // Effective Viscosity
const double deltaT = 0.01   ; // Size of time step
const double rho    = 1.     ; // (constant) density
const    int Nts    = 100     ; // Number of time steps
const    int k      = 8      ; // Wavelet order (usually precision + 2)
const double thresh = 1e-6   ; // Precision
const double thresh1= 0.1*thresh ;

const double L = 2*WST_PI;
const double Lx = 2. * WST_PI;
const double Ly = 1. * WST_PI;
const double Lz = 1. * WST_PI;
const double N = 8.0;

double mytime = 0.0  ; // Global variable for the current time 
                       // This should be passed in thru the class or app context
const double cc = 2 * WST_PI/deltaT/Nts;

//---------------------------------------------------------------------------**
static double init_zero(const coordT3d& r)
{
  return 0.0 ;
}


// the following define the exact body forces, velocity fields, and pressure
// for a known test problem


//---------------------------------------------------------------------------**
static double uxexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return  cos(t)*sin(x)*sin(x)* (sin(2.*y)*sin(z)*sin(z) - sin(y)*sin(y)*sin(2.*z)) ;
}
//---------------------------------------------------------------------------**
static double uyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return  cos(t)*sin(y)*sin(y)* (sin(2.*z)*sin(x)*sin(x) - sin(z)*sin(z)*sin(2.*x)) ;
}
//---------------------------------------------------------------------------**
static double uzexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return  cos(t)*sin(z)*sin(z)* (sin(2.*x)*sin(y)*sin(y) - sin(x)*sin(x)*sin(2.*y)) ;
}

//---------------------------------------------------------------------------**
static double fxexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double value = -(cos(t)*sin(x)*(cos(z) + 4*cos(y - 2*z)*sin(x))*sin(y)) +
            2*(2*cos(t)*cos(2*y - z)*pow(sin(x),2) +
            (2*cos(t)*cos(2*x) + sin(t)*pow(sin(x),2))*sin(y)*sin(y - z))*
            sin(z) - 2*pow(cos(t),2)*
            (-2*cos(x) + cos(x - 2*y) + cos(x - 2*z))*pow(sin(x),3)*
            pow(sin(y),2)*pow(sin(z),2);
  return  value ;
}
//---------------------------------------------------------------------------**
static double fyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double value = -2*sin(t)*sin(x)*pow(sin(y),2)*sin(x - z)*sin(z) -
          2*pow(cos(t),2)*(cos(2*x - y) - 2*cos(y) + cos(y - 2*z))*
          pow(sin(x),2)*pow(sin(y),3)*pow(sin(z),2) +
          cos(t)*(cos(x)*cos(y)*cos(z) +
          4*cos(x - 2*z)*sin(x)*pow(sin(y),2) -
          4*(cos(2*x - z)*pow(sin(y),2) + cos(2*y)*sin(x)*sin(x - z))*
          sin(z));
  return  value ;
}
//---------------------------------------------------------------------------**
static double fzexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double value = 2*sin(t)*sin(x)*sin(x - y)*sin(y)*pow(sin(z),2) -
          2*pow(cos(t),2)*(cos(2*x - z) + cos(2*y - z) - 2*cos(z))*
          pow(sin(x),2)*pow(sin(y),2)*pow(sin(z),3) +
          cos(t)*(4*cos(2*z)*sin(x)*sin(x - y)*sin(y) +
          sin(z)*(-(cos(x)*sin(y)) +
          2*(sin(2*x) - 2*sin(2*(x - y)) - sin(2*y))*sin(z)));
  return  value ;
}

//---------------------------------------------------------------------------**
static double pexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return cos(t)*cos(x)*sin(y)*cos(z) ;
}


//---------------------------------------------------------------------------**
// A simple Laplacian operator
class Laplace: public Operator<functionT> {
  protected:
    void action (const functionT& x, functionT& Ax) const {
      FunctionDefaults<3>::set_refine(false);

      Ax = diff( diff(x,0), 0)
         + diff( diff(x,1), 1)
         + diff( diff(x,2), 2);

      Ax.truncate(thresh1);
      return;
    }
};

template<typename T, int n>
class LaplaceBiCG: public KryFun<T> {
  inline void operator() (const T& x, T& Ax) {
    FunctionDefaults<3>::set_refine(false);

    Ax = diff( diff(x,0), 0)
       + diff( diff(x,1), 1)
       + diff( diff(x,2), 2) ;

    Ax.truncate(thresh1) ;
    return;
  }
};

/*
//----------------------------------------------------------------------------
// two interacting vortices
static double uxexact(const coordT3d& r) {
  const double x = r[0], y = r[1], z = r[2];
  double t = mytime;

  return cos(t) * sin(x) * sin(x) * (sin(2. * y) * sin(z) * sin(z) - sin(y)
      * sin(y) * sin(2. * z));
  //return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double uyexact(const coordT3d& r) {
  const double x = r[0], y = r[1], z = r[2];
  double t = mytime;

  return cos(t) * sin(y) * sin(y) * (sin(2. * z) * sin(x) * sin(x) - sin(z)
      * sin(z) * sin(2. * x));
  //return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double uzexact(const coordT3d& r) {
  const double x = r[0], y = r[1], z = r[2];
  double t = mytime;

  return cos(t) * sin(z) * sin(z) * (sin(2. * x) * sin(y) * sin(y) - sin(x)
      * sin(x) * sin(2. * y));
  //return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double fxexact(const coordT3d& r)
{
  return 0. ;
}
static double fyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double x0 = Lx / 2. ;
  double y0 = Ly / 2. ;
  double Lparam = 0.25 ;

  // Creates two interacting vortices
     double myval =  18.* exp(- ( (x-x0)*(x-x0) + (y-y0)*(y-y0) ) / Lparam)  ;
  return myval ;
}
static double fzexact(const coordT3d& r)
{
  return 0. ;
}

static double pexact(const coordT3d& r) {
  const double x = r[0], y = r[1], z = r[2];
  double t = mytime;

        return cos(t) * cos(x) * sin(y) * cos(z);
}
*/


//---------------------------------------------------------------------------**
   
inline functionT div(const functT& uint) {
	return diff(uint[0],0) + diff(uint[1],1) + diff(uint[2],2);
}

//---------------------------------------------------------------------------**

void testNavierStokes(World& world)
{
  // Function defaults
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_refine(false);
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(0,L) ;
  FunctionDefaults<3>::set_thresh(thresh);
  Tensor<int> bc(3,2), bc0(3,2);
  bc=1; // periodic
  bc0=0;
  FunctionDefaults<3>::set_bc(bc) ;
    
  // Initialize initial solution 
  functT u(3);     // velocity field
  functT us(3);    // half-step velocity field
  functionT p;     // pressure
  functionT b;     // RHS in poisson equation
 
  functT u0(3);    // t=0 exact velocity
  functionT p0;    // t=0 exact pressure

  functT ul(3);     // linear solver velocity field
  functT usl(3);    // linear solver half-step velocity field
  functionT pl;     // linear solver pressure
  functionT bl;     // linear solver RHS
 
  functT f(3);     // external body force
  functT ue(3);    // exact velocity field
  functionT pe;    // exact pressure

  // solver for the pressure
  Laplace A;
  LaplaceBiCG<functionT,1> ABiCG;
  BiCG0<functionT, double> BiSolver;

  mytime = 0.;

  u[0] = FunctionFactory<double,3>(world).f(uxexact);
  u[1] = FunctionFactory<double,3>(world).f(uyexact);
  u[2] = FunctionFactory<double,3>(world).f(uzexact);
  p    = FunctionFactory<double,3>(world).f(pexact);
  b    = FunctionFactory<double,3>(world).f(init_zero);

  ul[0] = FunctionFactory<double,3>(world).f(uxexact);
  ul[1] = FunctionFactory<double,3>(world).f(uyexact);
  ul[2] = FunctionFactory<double,3>(world).f(uzexact);
  pl    = FunctionFactory<double,3>(world).f(pexact);
  bl    = FunctionFactory<double,3>(world).f(init_zero);

  u0[0] = FunctionFactory<double,3>(world).f(uxexact);
  u0[1] = FunctionFactory<double,3>(world).f(uyexact);
  u0[2] = FunctionFactory<double,3>(world).f(uzexact);
  p0    = FunctionFactory<double,3>(world).f(pexact);

  // begin time evolution  
  for (int t = 1 ; t <= Nts; t++) { 
    if (world.rank() == 0) {
      printf("Timestep %d of %d, time %f\n",t,Nts,mytime);
    }
  
    f[0] = FunctionFactory<double,3>(world).f(fxexact);
    f[1] = FunctionFactory<double,3>(world).f(fyexact);
    f[2] = FunctionFactory<double,3>(world).f(fzexact);

    // construct u*
    for (int dir=XDIR; dir<=ZDIR; dir++) {
      us[dir] = u[dir] - deltaT*(u[XDIR]*diff(u[dir],XDIR) +
                                 u[YDIR]*diff(u[dir],YDIR) +
                                 u[ZDIR]*diff(u[dir],ZDIR) ) 
                       + deltaT*mu* ( diff(diff(u[dir],XDIR),XDIR) +
                                      diff(diff(u[dir],YDIR),YDIR) +
                                      diff(diff(u[dir],ZDIR),ZDIR) )
                       + deltaT*f[dir];
     
      usl[dir] = ul[dir] - deltaT*(ul[XDIR]*diff(ul[dir],XDIR) +
                                 ul[YDIR]*diff(ul[dir],YDIR) +
                                 ul[ZDIR]*diff(ul[dir],ZDIR) ) 
                       + deltaT*mu* ( diff(diff(ul[dir],XDIR),XDIR) +
                                      diff(diff(ul[dir],YDIR),YDIR) +
                                      diff(diff(ul[dir],ZDIR),ZDIR) )
                       + deltaT*f[dir];
   }

    world.gop.fence();

    
    // solve for the new pressure
    b = rho/deltaT*div(us);
    double tol = thresh;
    //int flag = GMRES<functionT, double>(A, b, p, MAX_ITER, tol, true);
    int max_iter = MAX_ITER;
    int flag = BiSolver(ABiCG, p, b, max_iter, tol);
    if (world.rank() == 0) {
      printf("  out flag = %d, error = %e, steps = %d \n", flag, tol, max_iter);
    }

    if (world.rank() == 0) {
      printf("  Applying linear operator...\n");
    }
    bl = rho/deltaT*div(usl);
    FunctionDefaults<3>::set_bc(bc0);
    Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
    SeparatedConvolution<double, 3> op = PeriodicCoulombOp<double, 3> (world,
       k, thresh1, thresh1, cellsize);
    bl.set_bc(bc0);
    bl.truncate();
    p.truncate();
    pl = apply(op, bl);
    pl.scale(-1. / (4. * WST_PI)).set_bc(bc);
    FunctionDefaults<3>::set_bc(bc);
    bl.set_bc(bc);
    bl.truncate();
    p.truncate();
    
    world.gop.fence();

    // update u 
    for (int dir=XDIR; dir<=ZDIR; dir++) {
      u[dir] = us[dir] - deltaT/rho*diff(p,dir);
      ul[dir] = usl[dir] - deltaT/rho*diff(pl,dir);
    }


    world.gop.fence();

    // update time and exact solutions
    mytime += deltaT;
    ue[0] = FunctionFactory<double,3>(world).f(uxexact);
    ue[1] = FunctionFactory<double,3>(world).f(uyexact);
    ue[2] = FunctionFactory<double,3>(world).f(uzexact);
    pe    = FunctionFactory<double,3>(world).f(pexact);
    world.gop.fence();

    u[0].reconstruct();
    u[1].reconstruct();
    u[2].reconstruct();
    ul[0].reconstruct();
    ul[1].reconstruct();
    ul[2].reconstruct();
    ue[0].reconstruct();
    ue[1].reconstruct();
    ue[2].reconstruct();
    p.reconstruct();
    pl.reconstruct();
    pe.reconstruct();
    
    world.gop.fence();

/*
    double bstep = L / 15.0;
    for (int i=0; i<=15; i++)
    {
      //coordT3d pt(i*bstep);
      coordT3d pt;
      pt[0] = (i*bstep);
      pt[1] = (i*bstep/2.);
      pt[2] = (i*bstep/4.);
      double uxerror = fabs(ue[0](pt) - ul[0](pt));
      //double uyerror = fabs(ue[1](pt) - u[1](pt));
      //double uzerror = fabs(ue[2](pt) - u[2](pt));
      if (world.rank() == 0) {
        printf("%.2f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", pt[0], ue[0](pt), ul[0](pt), uxerror, uxerror / ue[0](pt), pe(pt), pl(pt),fabs(pe(pt)-pl(pt))/pe(pt));
      }
    }
    world.gop.fence();
*/

    // collect norms
    functionT uxdiff, uydiff, uzdiff, pdiff;
    // norms for iterative solver
    uxdiff = ue[0] - u[0];
    uydiff = ue[1] - u[1];
    uzdiff = ue[2] - u[2];
    pdiff  = pe    - p;
    double uxnorm = uxdiff.norm2() ;
    double uynorm = uydiff.norm2() ;
    double uznorm = uzdiff.norm2() ;
    double pnorm  = pdiff.norm2() ;

    // norms for linear solver
    uxdiff = ue[0] - ul[0];
    uydiff = ue[1] - ul[1];
    uzdiff = ue[2] - ul[2];
    pdiff  = pe    - pl;
    double ulxnorm = uxdiff.norm2() ;
    double ulynorm = uydiff.norm2() ;
    double ulznorm = uzdiff.norm2() ;
    double plnorm  = pdiff.norm2() ;

    // norm of current solution versus t=0 solution
    uxdiff = ue[0] - u0[0];
    uydiff = ue[1] - u0[1];
    uzdiff = ue[2] - u0[2];
    pdiff  = pe    - p0;
    double udxnorm = uxdiff.norm2() ;
    double udynorm = uydiff.norm2() ;
    double udznorm = uzdiff.norm2() ;
    double pdnorm  = pdiff.norm2() ;

    if (world.rank() == 0) {
      printf(" %d Iterative Norms: %e %e %e %e\n", t, uxnorm, uynorm, uznorm, pnorm);
      printf(" %d Linear Norms: %e %e %e %e\n", t, ulxnorm, ulynorm, ulznorm, plnorm);
      printf(" %d Diff Norms: %e %e %e %e\n", t, udxnorm, udynorm, udznorm, pdnorm);

      // plot out data
      Vector<double, 3> plotlo, plothi;
      Vector<long, 3> numpt;
      for(int i = 0; i < 3; ++i) {
        plotlo[i] = 0;
        plothi[i] = L;
        numpt[i] = 50;
      }
      bool binary = false;

      char output_filename[100];
      sprintf(output_filename, "data-%03d.vts", t);
      plotvtk_begin(world, output_filename, plotlo, plothi, numpt, binary);
      plotvtk_data(u[0], "u", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_data(ul[0], "ul", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_data(ue[0], "ue", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_data(p, "p", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_data(pl, "pl", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_data(pe, "pe", world, output_filename, plotlo, plothi, 
                   numpt, binary);
      plotvtk_end<3>(world, output_filename, binary);

    }
    world.gop.fence();

  }
} // end testNavierStokes
//---------------------------------------------------------------------------**

//---------------------------------------------------------------------------**
int main(int argc, char**argv)
{
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);

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
//---------------------------------------------------------------------------**

