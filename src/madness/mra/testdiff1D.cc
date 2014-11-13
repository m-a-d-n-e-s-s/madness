/// \file testdiff1D.cc
/// \brief testing for diff() in 1D

#include <madness/mra/mra.h>

using namespace madness;

typedef Vector<double,1> coordT;
typedef Function<double,1> functionT;
typedef FunctionFactory<double,1> factoryT;
typedef Tensor<double> tensorT;

static const int k = 9 ; // Wavelet order (usually precision + 2)
static const double thresh = 1.e-7 ; // Precision
static const int init_lev = 2;
static const int test_axis = 0;
static const double Length = 4.;

int compare(World& world, functionT test, functionT exact, const char *str)
{
   double error = (exact - test).norm2() ;
   int success;

   if (error < thresh) {
       success=0;
   }
   else {
       success=1;
   }
   /* There should be a reduction over success here... */
   if (world.rank() == 0) {
       std::cerr << "Error in " << str << ": " << error ;
       std::cerr << ((success==0) ? " PASSED " : " FAILED ") << std::endl ;
   }
   return success;
}


// Testing the derivatives of a function 1+x
// for a variety of boundary conditions
static double u_exact(const coordT &pt) {
  double x = pt[0];
  return (1.0+x) ;
}

static double du_exact(const coordT &pt) {
  return (1.0) ;
}

static double left_dirichlet(const coordT &pt) {
  return (1.0) ;
}
static double right_dirichlet(const coordT &pt) {
  double x = Length;
  return (1.0+x) ;
}
static double left_neumann  (const coordT &pt) {
  return (1.0) ;
}
static double right_neumann  (const coordT &pt) {
  return (1.0) ;
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    int success=0;

        std::cout.precision(6);

       // Function defaults
        FunctionDefaults<1>::set_k(k);
        FunctionDefaults<1>::set_thresh(thresh);
        FunctionDefaults<1>::set_refine(true );
        FunctionDefaults<1>::set_autorefine(true );
        FunctionDefaults<1>::set_initial_level(init_lev);
        FunctionDefaults<1>::set_cubic_cell( 0. , Length);

        BoundaryConditions<1> bc;

        functionT  u      = factoryT(world).f( u_exact );
        functionT due     = factoryT(world).f(du_exact );

        functionT  left_d = factoryT(world).f( left_dirichlet) ;
        functionT right_d = factoryT(world).f(right_dirichlet) ;
        functionT  left_n = factoryT(world).f( left_neumann  ) ;
        functionT right_n = factoryT(world).f(right_neumann  ) ;

        // Right B.C.: Dirichlet
        // Left  B.C.: Free
        bc(0,0) = BC_DIRICHLET ;
        bc(0,1) = BC_FREE ;

        Derivative<double,1> dx1(world, test_axis, bc, left_d, right_d, k) ;
        functionT du1 = dx1(u) ;
        success+=compare(world, du1, due, "du1") ;

        // Right B.C.: Free
        // Left  B.C.: Dirichlet
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_DIRICHLET ;

        Derivative<double,1> dx2(world, test_axis, bc, left_d, right_d, k) ;
        functionT du2 = dx2(u) ;
        success+=compare(world, du2, due, "du2") ;

        // Right B.C.: Neumann
        // Left  B.C.: Free
        bc(0,0) = BC_NEUMANN ;
        bc(0,1) = BC_FREE ;

        Derivative<double,1> dx3(world, test_axis, bc, left_n, right_n, k) ;
        functionT du3 = dx3(u) ;
        success+=compare(world, du3, due, "du3") ;

        // Right B.C.: Free
        // Left  B.C.: Neumann
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_NEUMANN ;

        Derivative<double,1> dx4(world, test_axis, bc, left_n, right_n, k) ;
        functionT du4 = dx4(u) ;
        success+=compare(world, du4, due, "du4") ;

         world.gop.fence();

    finalize();

    return success;
}



