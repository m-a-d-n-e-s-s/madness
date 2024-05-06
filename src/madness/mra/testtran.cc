#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

using namespace madness;

static const double L = 32.0;   // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-7; // precision

double psi(const Vector<double,3>& r) {
  return exp(-sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6));
}

double V(const Vector<double,3>& r) {
  return -1.0/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6);
}

int main(int argc, char**argv) {
  initialize(argc,argv);
  World world(SafeMPI::COMM_WORLD);

  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  std::cout.precision(8);

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

  print(psi({0.0,0.0,0.0}), psi({0.0,0.0,1.0}));

  Function<double,3> u = FunctionFactory<double,3>(world).f(psi);
  Function<double,3> v = FunctionFactory<double,3>(world).f(V);

  u.print_tree();
  v.print_tree();
  

  std::vector<Function<double,3>> w(2);
  w[0] = u;
  w[1] = v;

  Tensor<double> c(2,2);
  c(0,0) = 9.0;
  c(1,1) = 4.0;
  c(0,1) = 0.5;
  c(1,0) =-0.5;

  auto r = transform(world, w, c, 0.0, true);

  auto r2 = transform(world, w, c);

  print("err ", norm2s(world, r-r2));

  finalize();

  return 0;
}

