/// \file M3Ddifftesting.cc
/// \brief testing for diff() in 3D

/// vectorg only takes in 1 function, for now. the workaround is to call it twice.
/// i.e. set bc(0,0) to a dummy value - not 3, nor 5 (vectorg does nothing for other values of bc)
/// and bc(0,1) to whatever bc you're actually using (3, or 5)
/// then call vectorg, then switch bc(0,1) to a dummy value, and bc(0,0) to 3 or 5,
/// and call vectorg again

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <string>
#include <mra/mra.h>
#include <mra/mraimpl.h>

using namespace madness;

typedef Vector<double,3> coordT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef Tensor<double> tensorT;

static const double PI = 3.1415926535897932384;
static const int k = 9 ; // Wavelet order (usually precision + 2)
static const double thresh = 1.e-7 ; // Precision
static const int init_lev = 2; 
static int test_axis = 0; 

void compare(functionT test, functionT exact, const char *str)
{
   double error = (exact - test).norm2() ;

   std::cerr << "Error in " << str << ": " << error ;

   if (error < thresh)
     std::cerr << " PASSED " << std::endl ;
   else
     std::cerr << " FAILED " << std::endl ;

   return ;
}



// Testing the derivatives of a function 1+x
// for a variety of boundary conditions
static double u_exact(const coordT &pt) {
  return (1.0+pt[0]*pt[1]*pt[2]) ;
}

static double dudx_exact(const coordT &pt) {
  return (pt[1]*pt[2]) ;
}

static double dudy_exact(const coordT &pt) {
  return (pt[0]*pt[2]) ;
}

static double dudz_exact(const coordT &pt) {
  return (pt[0]*pt[1]) ;
}

static double xleft_dirichlet(const coordT &pt) {
  return (1.) ;
}
static double xright_dirichlet(const coordT &pt) {
  return (1.+pt[1]*pt[2]) ;
}
static double xleft_neumann  (const coordT &pt) {
  return (pt[1]*pt[2]) ;
}
static double xright_neumann  (const coordT &pt) {
  return (pt[1]*pt[2]) ;
}

static double yleft_dirichlet(const coordT &pt) {
  return (1.) ;
}
static double yright_dirichlet(const coordT &pt) {
  return (1.+pt[0]*pt[2]) ;
}
static double yleft_neumann  (const coordT &pt) {
  return (pt[0]*pt[2]) ;
}
static double yright_neumann  (const coordT &pt) {
  return (pt[0]*pt[2]) ;
}

static double zleft_dirichlet(const coordT &pt) {
  return (1.) ;
}
static double zright_dirichlet(const coordT &pt) {
  return (1.+pt[0]*pt[1]) ;
}
static double zleft_neumann  (const coordT &pt) {
  return (pt[0]*pt[1]) ;
}
static double zright_neumann  (const coordT &pt) {
  return (pt[0]*pt[1]) ;
}


int main(int argc, char** argv) {
    
	MPI::Init(argc, argv);
	World world(MPI::COMM_WORLD);
        startup(world,argc,argv);

        std::cout.precision(6);

       // Function defaults
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_refine(true );
        FunctionDefaults<3>::set_autorefine(true );
        FunctionDefaults<3>::set_initial_level(init_lev);
        FunctionDefaults<3>::set_cubic_cell( 0. , 1.);

        Tensor<int> bctensor(3,2) ;
        BoundaryConds<3> bdry_cs(bctensor);

        functionT  u      = factoryT(world).f(   u_exact );
        functionT dudxe   = factoryT(world).f(dudx_exact );
        functionT dudye   = factoryT(world).f(dudy_exact );
        functionT dudze   = factoryT(world).f(dudz_exact );

        functionT  xleft_d = factoryT(world).f( xleft_dirichlet) ;
        functionT xright_d = factoryT(world).f(xright_dirichlet) ;
        functionT  xleft_n = factoryT(world).f( xleft_neumann  ) ;
        functionT xright_n = factoryT(world).f(xright_neumann  ) ;

        functionT  yleft_d = factoryT(world).f( yleft_dirichlet) ;
        functionT yright_d = factoryT(world).f(yright_dirichlet) ;
        functionT  yleft_n = factoryT(world).f( yleft_neumann  ) ;
        functionT yright_n = factoryT(world).f(yright_neumann  ) ;

        functionT  zleft_d = factoryT(world).f( zleft_dirichlet) ;
        functionT zright_d = factoryT(world).f(zright_dirichlet) ;
        functionT  zleft_n = factoryT(world).f( zleft_neumann  ) ;
        functionT zright_n = factoryT(world).f(zright_neumann  ) ;

        // Derivative in the x-direction
        test_axis = 0;

        // X Right B.C.: Dirichlet
        // X Left  B.C.: Free
        bctensor(0,0) = 3 ;
        bctensor(0,1) = 2 ;
        bctensor(1,0) = 2 ;
        bctensor(1,1) = 2 ;
        bctensor(2,0) = 2 ;
        bctensor(2,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dx1(world, test_axis, bdry_cs, xleft_d, xright_d, k) ;
        functionT dudx1 = dx1(u) ;
        compare(dudx1, dudxe, "dudx1") ;

        // X Right B.C.: Free
        // X Left  B.C.: Dirichlet
        bctensor(0,0) = 2 ;
        bctensor(0,1) = 3 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dx2(world, test_axis, bdry_cs, xleft_d, xright_d, k) ;
        functionT dudx2 = dx2(u) ;
        compare(dudx2, dudxe, "dudx2") ;

        // X Right B.C.: Neumann
        // X Left  B.C.: Free
        bctensor(0,0) = 5 ;
        bctensor(0,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dx3(world, test_axis, bdry_cs, xleft_n, xright_n, k) ;
        functionT dudx3 = dx3(u) ;
        compare(dudx3, dudxe, "dudx3") ;

        // X Right B.C.: Free
        // X Left  B.C.: Neumann
        bctensor(0,0) = 2 ;
        bctensor(0,1) = 5 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dx4(world, test_axis, bdry_cs, xleft_n, xright_n, k) ;
        functionT dudx4 = dx4(u) ;
        compare(dudx4, dudxe, "dudx4") ;


        // Derivative in the y-direction
        test_axis = 1;

        // Y Right B.C.: Dirichlet
        // Y Left  B.C.: Free
        bctensor(0,0) = 2 ;
        bctensor(0,1) = 2 ;
        bctensor(1,0) = 3 ;
        bctensor(1,1) = 2 ;
        bctensor(2,0) = 2 ;
        bctensor(2,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dy1(world, test_axis, bdry_cs, yleft_d, yright_d, k) ;
        functionT dudy1 = dy1(u) ;
        compare(dudy1, dudye, "dudy1") ;

        // Y Right B.C.: Free
        // Y Left  B.C.: Dirichlet
        bctensor(1,0) = 2 ;
        bctensor(1,1) = 3 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dy2(world, test_axis, bdry_cs, yleft_d, yright_d, k) ;
        functionT dudy2 = dy2(u) ;
        compare(dudy2, dudye, "dudy2") ;

        // Y Right B.C.: Neumann
        // Y Left  B.C.: Free
        bctensor(1,0) = 5 ;
        bctensor(1,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dy3(world, test_axis, bdry_cs, yleft_n, yright_n, k) ;
        functionT dudy3 = dy3(u) ;
        compare(dudy3, dudye, "dudy3") ;

        // Y Right B.C.: Free
        // Y Left  B.C.: Neumann
        bctensor(1,0) = 2 ;
        bctensor(1,1) = 5 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dy4(world, test_axis, bdry_cs, yleft_n, yright_n, k) ;
        functionT dudy4 = dy4(u) ;
        compare(dudy4, dudye, "dudy4") ;


        // Derivative in the z-direction
        test_axis = 2;

        // Z Right B.C.: Dirichlet
        // Z Left  B.C.: Free
        bctensor(0,0) = 2 ;
        bctensor(0,1) = 2 ;
        bctensor(1,0) = 2 ;
        bctensor(1,1) = 2 ;
        bctensor(2,0) = 3 ;
        bctensor(2,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dz1(world, test_axis, bdry_cs, zleft_d, zright_d, k) ;
        functionT dudz1 = dz1(u) ;
        compare(dudz1, dudze, "dudz1") ;

        // Z Right B.C.: Free
        // Z Left  B.C.: Dirichlet
        bctensor(2,0) = 2 ;
        bctensor(2,1) = 3 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dz2(world, test_axis, bdry_cs, zleft_d, zright_d, k) ;
        functionT dudz2 = dz2(u) ;
        compare(dudz2, dudze, "dudz2") ;

        // Z Right B.C.: Neumann
        // Z Left  B.C.: Free
        bctensor(2,0) = 5 ;
        bctensor(2,1) = 2 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dz3(world, test_axis, bdry_cs, zleft_n, zright_n, k) ;
        functionT dudz3 = dz3(u) ;
        compare(dudz3, dudze, "dudz3") ;

        // Z Right B.C.: Free
        // Z Left  B.C.: Neumann
        bctensor(2,0) = 2 ;
        bctensor(2,1) = 5 ;
        bdry_cs = bctensor ;

        Derivative<double,3> dz4(world, test_axis, bdry_cs, zleft_n, zright_n, k) ;
        functionT dudz4 = dz4(u) ;
        compare(dudz4, dudze, "dudz4") ;

         world.gop.fence();

    finalize();
    
    return 0;
}

        
        
