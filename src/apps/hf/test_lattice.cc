#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "poperator.h"

using namespace madness;

typedef Vector<double,3> coordT3d;
typedef Vector<double,1> coordT1d;

const double L = 7.0;
const double N = 8.0;

//*****************************************************************************
static double rho_bsh_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  return  cos(twopi*x/L) * cos(twopi*y/L) * cos(twopi*z/L);
}
//*****************************************************************************

//*****************************************************************************
static double rho_coulomb_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double npi = N * WST_PI;
  return  cos(npi*x/L) * cos(npi*y/L) * cos(npi*z/L);
}
//*****************************************************************************

//*****************************************************************************
static double phi_coulomb_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double npi = N * WST_PI;
  double threepi = 3 * WST_PI;
  return  (4.0*L*L/(threepi*N*N)) * cos(npi*x/L) * cos(npi*y/L) * cos(npi*z/L);
}
//*****************************************************************************

////*****************************************************************************
//static double phi_func3d(const coordT3d& r)
//{
//  const double x=r[0], y=r[1], z=r[2];
//  double twopi = 2 * WST_PI;
//  double threepi = 3 * WST_PI;
//  return  (1/threepi) * cos(twopi*x) * cos(twopi*y) * cos(twopi*z);
//}
////*****************************************************************************

//*****************************************************************************
static double phi_bsh_func3d(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double twopi = 2 * WST_PI;
  double sixteenpisquared = 16.0 * WST_PI * WST_PI;
  return  (L*L/sixteenpisquared) * cos(twopi*x/L) * cos(twopi*y/L) * cos(twopi*z/L);
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
void testPeriodicGaussian(World& world, double coeff, double expnt, int lmax, int k, double thresh, double* data)
{
  // Function defaults
  double bsize = 0.5;
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(-bsize,bsize);
  FunctionDefaults<3>::set_thresh(thresh);

  // Test function
  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_coulomb_func3d);

  // Create operator
  std::vector< SharedPtr< Convolution1D<double> > > ops(1);
  ops[0] = SharedPtr< Convolution1D<double> >(new PeriodicGaussianConvolution1D<double>(k, 2, coeff, expnt));
  SeparatedConvolution<double,3> op(world, k, ops);

  // Apply operator
  Function<double,3> phi = apply(op, rho);

  for (int i=0; i<21; i++)
  {
    coordT3d p(-0.5 + i*0.05);
    printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], phi(p), data[i], fabs(phi(p) - data[i]));
  }
}
//*****************************************************************************

//*****************************************************************************
void testSinglePeriodicGaussians(int argc, char** argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  double maple_data_2500[21] =
    {
      -44.02214685,
      -37.86955441,
      -23.31010082,
      -8.939789115,
      -1.299027398,
      0.0,
      1.299027398,
      8.939789115,
      23.31010082,
      37.86955441,
      44.02214685,
      37.86955441,
      23.31010082,
      8.939789115,
      1.299027398,
      0.0,
      -1.299027398,
      -8.939789115,
      -23.31010082,
      -37.86955441,
      -44.02214685
    };

  double maple_data_25000[21] =
    {
      -1.407020542,
      -1.210373524,
      -0.7450293332,
      -0.2857304296,
      -0.04151906172,
      0.0,
      0.04151906172,
      0.2857304296,
      0.7450293332,
      1.210373524,
      1.407020542,
      1.210373524,
      0.7450293332,
      0.2857304296,
      0.04151906172,
      0.0,
      -0.04151906172,
      -0.2857304296,
      -0.7450293332,
      -1.210373524,
      -1.407020542
    };

  double maple_data_55000[21] =
    {
      -0.4314663948,
      -0.3711640906,
      -0.2284651223,
      -0.08761995618,
      -0.01273192489,
      0.0,
      0.01273192489,
      0.08761995618,
      0.2284651223,
      0.3711640906,
      0.4314663948,
      0.3711640906,
      0.2284651223,
      0.08761995618,
      0.01273192489,
      0.0,
      -0.01273192489,
      -0.08761995618,
      -0.2284651223,
      -0.3711640906,
      -0.4314663948
    };

  double maple_data_95000[21] =
    {
      -0.1901095982,
      -0.1635396336,
      -0.1006646476,
      -0.03860647055,
      -0.005609848544,
      0.0,
      0.005609848544,
      0.03860647055,
      0.1006646476,
      0.1635396336,
      0.1901095982,
      0.1635396336,
      0.1006646476,
      0.03860647055,
      0.005609848544,
      0.0,
      -0.005609848544,
      -0.03860647055,
      -0.1006646476,
      -0.1635396336,
      -0.1901095982
    };

  int k = 8;
  double thresh = 1e-6;
  printf("\nTesting with exponent = 2500\n\n");
  testPeriodicGaussian(world, 100, 2500, 16, k, thresh, &maple_data_2500[0]);
  printf("\nTesting with exponent = 25000\n\n");
  testPeriodicGaussian(world, 100, 25000, 16, k, thresh, &maple_data_25000[0]);
  printf("\nTesting with exponent = 55000\n\n");
  testPeriodicGaussian(world, 100, 55000, 16, k, thresh, &maple_data_55000[0]);
  printf("\nTesting with exponent = 95000\n\n");
  testPeriodicGaussian(world, 100, 95000, 16, k, thresh, &maple_data_95000[0]);
  MPI::Finalize();
}
//*****************************************************************************

//*****************************************************************************
void testPeriodicCoulomb3d(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  int k = 12;
  double thresh = 1e-10;
  double eps = 1e-10;
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
  FunctionDefaults<3>::set_thresh(thresh);

  // Create test charge density and the exact solution to Poisson's equation
  // with said charge density
  printf("building rho ...\n\n");
  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_coulomb_func3d);
  printf("building phi_exact ...\n\n");
  Function<double,3> phi_exact = FunctionFactory<double,3>(world).f(phi_coulomb_func3d);

  // Create operator and apply
  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
  SeparatedConvolution<double,3> op = PeriodicCoulombOp<double,3>(world, k, 1e-6, eps, cellsize);
  printf("applying operator ...\n\n");
  Function<double,3> phi_test = apply(op, rho);

  double bstep = L / 100.0;
  for (int i=0; i<101; i++)
  {
    coordT3d p(-L/2 + i*bstep);
    double error = fabs(phi_exact(p) - phi_test(p));
    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], phi_exact(p), phi_test(p), error, error / phi_exact(p));
  }

  // Plot to OpenDX
//  vector<long> npt(3,101);
//  Function<double,3> phi_diff = phi_exact - phi_test;
//  plotdx(phi_test, "phitest.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(phi_exact, "phiexact.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);

  MPI::Finalize();
}
//*****************************************************************************

//*****************************************************************************
void testPeriodicBSH3d(int argc, char**argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  int k = 8;
  double thresh = 1e-6;
  double eps = 1e-6;
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
  FunctionDefaults<3>::set_thresh(thresh);

  // Create test charge density and the exact solution to Poisson's equation
  // with said charge density
  Function<double,3> rho = FunctionFactory<double,3>(world).f(rho_bsh_func3d);
  Function<double,3> phi_exact = FunctionFactory<double,3>(world).f(phi_bsh_func3d);

  // Create operator and apply
  double twopi = 2 * WST_PI;
  Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
  SeparatedConvolution<double,3> op = PeriodicBSHOp<double,3>(world, twopi/L, k, 1e-8, eps, cellsize);
  Function<double,3> phi_test = apply(op, rho);

  double bstep = L / 100.0;
  for (int i=0; i<101; i++)
  {
    coordT3d p(-L/2 + i*bstep);
    printf("%.2f\t\t%.8f\t%.8f\t%.8f\t%.8f\n", p[0], phi_exact(p), phi_test(p), phi_exact(p)/phi_test(p), fabs(phi_exact(p)-phi_test(p)));
  }

  // Plot to OpenDX
//  printf("plotting to openDX ...\n\n");
//  vector<long> npt(3,101);
//  Function<double,3> phi_diff = phi_exact - phi_test;
//  plotdx(phi_test, "phitest.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(phi_exact, "phiexact.dx", FunctionDefaults<3>::get_cell(), npt);
//  plotdx(phi_diff, "phidiff.dx", FunctionDefaults<3>::get_cell(), npt);
//  printf("finished plotting to openDX ...\n\n");

  MPI::Finalize();
  printf("done!\n\n");
}
//*****************************************************************************

//*****************************************************************************
int main(int argc, char**argv)
{
  testPeriodicCoulomb3d(argc, argv);
//  testPeriodicBSH3d(argc, argv);
  return 0;
}
//*****************************************************************************

