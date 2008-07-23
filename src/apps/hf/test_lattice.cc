#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "poperator.h"

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Vector<double,1> coordT1d;

//*****************************************************************************
static double rho_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  return  cos(twopi*x) * cos(twopi*y) * cos(twopi*z);
}
//*****************************************************************************

//*****************************************************************************
static double phi_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  double threepi = 3 * WST_PI;
  return  (1/threepi) * cos(twopi*x) * cos(twopi*y) * cos(twopi*z);
}
//*****************************************************************************

//*****************************************************************************
static double rho_func1d(const coordT1d& r)
{
  const double x=r[0];
  double twopi = 2 * WST_PI;
  return  cos(twopi*x);
}
//*****************************************************************************

//*****************************************************************************
static double phi_func1d(const coordT1d& r)
{
  const double x=r[0];
  double twopi = 2 * WST_PI;
  return  (1/WST_PI) * cos(twopi*x);
}
//*****************************************************************************

//*****************************************************************************
template <typename Q, int NDIM>
Q laplacian(const Q& f) {
        Q lapf = diff(diff(f,0),0);
        for (int i=1; i<NDIM; ++i) lapf += diff(diff(f,i),i);
        return lapf;
};
//*****************************************************************************

//*****************************************************************************
template <typename Q> class wstFunctor
{
public:
  int kmax;
  double coeff, expnt;
  wstFunctor(int kmax, double coeff, double expnt) :
    kmax(kmax), coeff(coeff), expnt(expnt)
  {
  }
  Q operator()(double x) const
  { //x \in [-1,1] as checked
    Q sx0 = 0.0;
    for (int ki = -kmax; ki <= kmax; ki++)
    {
      double k = (double) ki;
      sx0 += exp(-expnt*(x-k)*(x-k));
    }
    return sx0*coeff;
  }
};
//*****************************************************************************

//*****************************************************************************
struct PeriodicConditionalRefineTest
{
  bool operator() (const Key<3>& key, const Tensor<double>& t) const
  {
    // Is the box above the level from where we want to refine?
    int n = key.level();
    if (n >= 6) return false;
    // Are we on the boundary?
    Translation l1 = key.translation()[0];
    Translation l2 = key.translation()[1];
    Translation l3 = key.translation()[2];
    Translation maxl = (1ul<<n) - 1;
    if ((l1 == 0 || l1 == maxl) || (l2 == 0 || l2 == maxl) || (l3 == 0 || l3 == maxl))
    {
      print("Refining ...\n", key);
      return true;
    }
    return false;
  }

  template <typename Archive>
  void serialize(const Archive& arch) {}
};
//*****************************************************************************

//*****************************************************************************
void test3d(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  double bsize = 0.5;
  int funck = 8;
  double thresh = 1e-8;
  double eps = 1e-12;
  Tensor<int> bc(3,2);
  bc(_) = 0;
  FunctionDefaults<3>::set_bc(bc);
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);

  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_func3d);
  Function<double,3> phi_exact = FunctionFactory<double,3>(world).f(phi_func3d);

  SeparatedConvolution<double,3> op = PBSHOperator<double,3>(world, 0.0, funck, 1e-8, eps);
  Function<double,3> phi_test = apply(op, rho);
  //phi_test.scale(-1.0);
  Function<double,3> rho_test = laplacian<Function<double,3>, 3>(phi_test);
  rho_test.scale(1.0/4.0/WST_PI);

//  coordT3d point1;
//  point1[0] = 0.254; point1[1] = 0.114; point1[2] = 0.854;
//  coordT3d point2;
//  point2[0] = 0.1054; point2[1] = 0.930114; point2[2] = 0.00054;
//
//  double pept1 = phi_exact(point1);
//  double pept2 = phi_exact(point2);
//  double ptpt1 = phi_test(point1);
//  double ptpt2 = phi_test(point2);
//  if (world.rank() == 0) printf("pept1 = %.8f\t\tptpt1 = %.8f\n\n", pept1, ptpt1);
//  if (world.rank() == 0) printf("pept2 = %.8f\t\tptpt2 = %.8f\n\n", pept2, ptpt2);
//  if (world.rank() == 0) printf("ptpt1 / pept1 = %.8f\n\n", ptpt1 / pept1);
//  if (world.rank() == 0) printf("ptpt2 / pept2 = %.8f\n\n", ptpt2 / pept2);

  vector<long> npt(3,101);
  Function<double,3> phi_diff = phi_exact - phi_test;
  Function<double,3> rho_diff = rho - rho_test;
  plotdx(phi_test, "phitest.dx", FunctionDefaults<3>::get_cell(), npt);
  plotdx(phi_exact, "phiexact.dx", FunctionDefaults<3>::get_cell(), npt);
  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);
  plotdx(phi_diff, "rhodiff.dx", FunctionDefaults<3>::get_cell(), npt);

  coordT3d point1;
  point1[0] = 0.49; point1[1] = 0.49; point1[2] = 0.49;
  coordT3d point2;
  point2[0] = 0.1; point2[1] = 0.1; point2[2] = 0.1;

  phi_exact.reconstruct(true);
  double pept1 = phi_exact(point1);
  double pept2 = phi_exact(point2);
  double ptpt1 = phi_test(point1);
  double ptpt2 = phi_test(point2);
  if (world.rank() == 0) printf("pept1 = %.8f\t\tptpt1 = %.8f\n\n", pept1, ptpt1);
  if (world.rank() == 0) printf("pept2 = %.8f\t\tptpt2 = %.8f\n\n", pept2, ptpt2);
  if (world.rank() == 0) printf("ptpt1 / pept1 = %.8f\n\n", ptpt1 / pept1);
  if (world.rank() == 0) printf("ptpt2 / pept2 = %.8f\n\n", ptpt2 / pept2);
  if (world.rank() == 0) printf("rept1 = %.8f\t\trtpt1 = %.8f\n\n", pept1, ptpt1);
  if (world.rank() == 0) printf("rept2 = %.8f\t\trtpt2 = %.8f\n\n", pept2, ptpt2);
  if (world.rank() == 0) printf("rtpt1 / rept1 = %.8f\n\n", ptpt1 / pept1);
  if (world.rank() == 0) printf("rtpt2 / rept2 = %.8f\n\n", ptpt2 / pept2);

  double error = (phi_exact - phi_test).norm2();
  if (world.rank() == 0) printf("Error is %.8f\n\n", error);

  double phitest_norm = phi_test.norm2();
  double phiexact_norm = phi_exact.norm2();
  double rho_norm = rho.norm2();
  if (world.rank() == 0) printf("Norm of phi_test = %.8f\n\n", phitest_norm);
  if (world.rank() == 0) printf("Norm of phi_exact = %.8f\n\n", phiexact_norm);
  if (world.rank() == 0) printf("Norm of rho = %.8f\n\n", rho_norm);

  MPI::Finalize();
}
//*****************************************************************************

////*****************************************************************************
//void test3djun(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  double bsize = 0.5;
//  int funck = 8;
//  double thresh = 1e-6;
//  FunctionDefaults<3>::bc(___) = 0;  //1 for periodic boundary conditions
//  FunctionDefaults<3>::set_k(funck);
//  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(2);
//  FunctionDefaults<3>::set_truncate_mode(1);
//  FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);
//
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_func3d);
//  SeparatedConvolution<double,3> op = PBSHOperator<double,3>(world, 0.0, funck, 1e-8, 1e-8);
//  Function<double, 3> u = apply(op, rho);
//  double ratio = 1.0/(3.0*WST_PI); print(ratio);
//  ratio = u.norm2()/rho.norm2(); printf("%5e\n", ratio);
//  print("PBSH a-error, norm2: ", (rho*ratio-u).norm2(),u.norm2());
//  FunctionDefaults<3>::bc(_,_) = 1;
//  Function<double,3> p=laplacian<Function<double,3>, 3>(rho);
//  print((rho+p.scale(1.0/(3*WST_PI))).norm2());
//
//  MPI::Finalize();
//}
////*****************************************************************************

//*****************************************************************************
void test1d(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  double bsize = 1.0;
  int funck = 8;
  double thresh = 1e-6;
  Tensor<int> bc(1,2);
  bc(0,0) = 1; bc(0,1) = 1;
  FunctionDefaults<1>::set_bc(bc);
  FunctionDefaults<1>::set_k(funck);
  FunctionDefaults<1>::set_thresh(thresh);
  FunctionDefaults<1>::set_refine(true);
  FunctionDefaults<1>::set_initial_level(2);
  FunctionDefaults<1>::set_truncate_mode(1);
  FunctionDefaults<1>::set_cubic_cell(0.0, bsize);

  Function<double,1> rho = FunctionFactory<double,1>(world).f(rho_func1d);
  Function<double,1> phi_exact = FunctionFactory<double,1>(world).f(phi_func1d);

  SeparatedConvolution<double,1> op = PBSHOperator<double,1>(world, 0.0, funck, 1e-8, 1e-8);
  Function<double,1> phi_test = apply(op, rho);

  coordT1d point1;
  point1[0] = 0.254;
  coordT1d point2;
  point2[0] = 0.1054;

  double pept1 = phi_exact(point1);
  double pept2 = phi_exact(point2);
  double ptpt1 = phi_test(point1);
  double ptpt2 = phi_test(point2);
  if (world.rank() == 0) printf("pept1 = %.8f\t\tptpt1 = %.8f\n\n", pept1, ptpt1);
  if (world.rank() == 0) printf("pept2 = %.8f\t\tptpt2 = %.8f\n\n", pept2, ptpt2);
  if (world.rank() == 0) printf("ptpt1 / pept1 = %.8f\n\n", ptpt1 / pept1);
  if (world.rank() == 0) printf("ptpt2 / pept2 = %.8f\n\n", ptpt2 / pept2);

  double error = (phi_exact - phi_test).norm2();
  if (world.rank() == 0) printf("Error is %.8f\n\n", error);

  double phitest_norm = phi_test.norm2();
  double phiexact_norm = phi_exact.norm2();
  double rho_norm = rho.norm2();
  if (world.rank() == 0) printf("Norm of phi_test = %.8f\n\n", phitest_norm);
  if (world.rank() == 0) printf("Norm of phi_exact = %.8f\n\n", phiexact_norm);
  if (world.rank() == 0) printf("Norm of rho = %.8f\n\n", rho_norm);

  MPI::Finalize();
}
//*****************************************************************************

//*****************************************************************************
void testscott(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  double bsize = 0.5;
  int funck = 16;
  double thresh = 1e-8;

  // BSH fit parameters
  double eps = 1e-6;
  double hi = 10.0;
  double lo = 1e-4;

  Tensor<int> bc(3,2);
  bc(___) = 1;
  FunctionDefaults<3>::set_bc(bc);
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(false);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);

  // Test function
  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_func3d);
//  print("before", rho.size());
//  rho.conditional_refine(PeriodicConditionalRefineTest());
//  print("after", rho.size());

  // Get exponents and coefficients
//  Tensor<double> coeff, expnt;
//  bsh_fit(0.0, lo, hi, eps, &coeff, &expnt, true);
//  // Needed to add the 1/4pi to the coefficients
//  coeff.scale(4.0*WST_PI);

  Tensor<double> coeff(3), expnt(3);
  coeff(1) = 27.0;
  coeff(2) = 64.0;
  coeff(3) = 1000000.0;
  expnt(1) = 10.0;
  expnt(2) = 1000.0;
  expnt(3) = 55000.0;

  printf("coeff(1) = %.5f\texptnt(1) = %.3f\n", coeff[1], expnt[1]);
  printf("coeff(2) = %.5f\texptnt(2) = %.3f\n", coeff[2], expnt[2]);
  printf("coeff(3) = %.5f\texptnt(3) = %.3f\n", coeff[3], expnt[3]);
  printf("coeff.dim[0] = %d", coeff.dim[0]);


  // Just looking to see of the 1/r is truly "1/r"
  double pt1 = 1e-1;
  double pt2 = 1e-2;
  double pt3 = 1e-3;
  double pt4 = 2;
  double f1 = foo(coeff, expnt, pt1);
  double f2 = foo(coeff, expnt, pt2);
  double f3 = foo(coeff, expnt, pt3);
  double f4 = foo(coeff, expnt, pt4);
  if (world.rank() == 0) printf("pt1 = %.5e\t\tf1 = %.5e\n", pt1, f1);
  if (world.rank() == 0) printf("pt2 = %.5e\t\tf2 = %.5e\n", pt2, f2);
  if (world.rank() == 0) printf("pt3 = %.5e\t\tf3 = %.5e\n", pt3, f3);
  if (world.rank() == 0) printf("pt4 = %.5e\t\tf4 = %.5e\n", pt4, f4);

  // Rescaling the coefficients for 3 dimensions
  for (int j = 1; j <= coeff.dim[0]; j++)
  {
    coeff[j] = (coeff[j] < 0.0) ? -pow(-coeff[j], 1.0/3.0) : pow(coeff[j], 1.0/3.0);
    if (world.rank() == 0) printf("expt[%d] = %.8f\t\tcoeff[%d] = %.8f\n", j, expnt[j], j, coeff[j]);
  }

  // The exponent being tested.
  double w = 80000.0;
  double c = 100.0;

  // One gaussian only (no lattice sum)
  std::vector< SharedPtr< Convolution1D<double> > > ops(1);
  ops[0]
        = SharedPtr< Convolution1D<double> >(new GenericConvolution1D< double,wstFunctor<double> >(funck,wstFunctor<double>(0, c, w)));
  SeparatedConvolution<double,3> op(world, funck, ops);
  Function<double,3> phi_test1 = apply(op, rho);

  // One gaussian only (with lattice)
  std::vector< SharedPtr< Convolution1D<double> > > opsumcoll(1);
  opsumcoll[0]
        = SharedPtr< Convolution1D<double> >(new GenericConvolution1D< double,wstFunctor<double> >(funck,wstFunctor<double>(6, c, w)));
  SeparatedConvolution<double,3> opsum(world, funck, opsumcoll);
  Function<double,3> phi_test2 = apply(opsum, rho);

  // One gaussian only (using different operator)
  Tensor<double> xcoeff(1), xexpnt(1);
  xcoeff[0] = c*c*c;
  xexpnt[0] = w;
  SeparatedConvolution<double,3> xop(world, funck, xcoeff, xexpnt);
  Function<double,3> phi_testx = apply(xop, rho);

   // Full set of gaussians (with lattice sum)
   std::vector< SharedPtr< Convolution1D<double> > > opfulllattice;
   for (int i = 1; i <= coeff.dim[0]; i++)
   {
     if (expnt[i] > 0.24)
     {
       // How many lattice spaces do I need to sum in real space?
       double dum = log(2.0/eps);
       int kmax = ceil(sqrt(dum/expnt[i]));
       printf("kmax[%d] = %d\t\tcoeff[%d] = %.8f\t\texpnt[%d] = %.8f\n", i, kmax, i, coeff[i], i, expnt[i]);
       // Add exponential with kmax, coeff, and expnt to list
       opfulllattice.push_back(SharedPtr< Convolution1D<double> >(new GenericConvolution1D< double,wstFunctor<double> >(funck,wstFunctor<double>(kmax, coeff[i], expnt[i]))));
     }
   }
   printf("creating opsumfull ...\n\n");
   SeparatedConvolution<double,3> opsumfull(world, funck, opfulllattice);
   printf("applying opsumfull ...\n\n");
   Function<double,3> phi_test3 = apply(opsumfull, rho);
   printf("done applying opsumfull ...\n\n");

  /// Point to be tested
  coordT3d point1(0.49);

   double ptpt1 = phi_test1(point1);
   double ptpt2 = phi_test2(point1);
   double ptptx = phi_testx(point1);
   double ptpt3 = phi_test3(point1);
  for (int i=0; i<101; i++) {
    coordT3d p(-0.5 + i*0.01);
    printf("%.2f  %.8f\n", p[0], phi_testx(p));
  }
   if (world.rank() == 0) printf("\n\n");
   if (world.rank() == 0) printf("ptpt1 = %.8e\n\n", ptpt1);
   if (world.rank() == 0) printf("ptpt2 = %.8e\n\n", ptpt2);
   if (world.rank() == 0) printf("ptptx = %.8e\n\n", ptptx);
//   if (world.rank() == 0) printf("ptpt3 = %.8e\n\n", ptpt3);

  MPI::Finalize();
}
//*****************************************************************************

////*****************************************************************************
//void testbob(int argc, char**argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Function defaults
//  double bsize = 0.5;
//  int funck = 8;
//  double thresh = 1e-8;
//  FunctionDefaults<3>::bc(___) = 1;  //1 for periodic boundary conditions
//  FunctionDefaults<3>::set_k(funck);
//  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(2);
//  FunctionDefaults<3>::set_truncate_mode(1);
//  FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);
//
//  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_func3d);
//
//  Function<double,3> rho_test = laplacian<Function<double,3>, 3>(rho);
//  rho_test.scale(-1.0/12.0/WST_PI/WST_PI);
//
//  Function<double,3> rho_diff = rho - rho_test;
//  if (world.rank() == 0) printf("Norm of rho_diff = %.8f\n\n", rho_diff.norm2());
//
//  vector<long> npt(3,101);
//  plotdx(rho, "rho.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(rho_test, "rho_test.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(rho_diff, "rho_diff.dx", FunctionDefaults<3>::get_cell(), npt);
//
//  MPI::Finalize();
//}
////*****************************************************************************

//*****************************************************************************
int main(int argc, char**argv)
{
  testscott(argc, argv);
  return 0;
}
//*****************************************************************************

