#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

using namespace madness;

static const double L = 32.0;   // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-4; // precision

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

  Function<double,3> u = FunctionFactory<double,3>(world).f(psi);
  Function<double,3> v = FunctionFactory<double,3>(world).f(V);

  u.truncate();
  v.truncate();
  u.compress();

  //print("u");
  //u.print_tree();
  // v.print_tree();

  std::vector<Function<double,3>> w(2);
  w[0] = u;
  w[1] = v;

  Tensor<double> c(2,2);
  c(0,0) = 2.0;
  c(1,1) = 3.0;
  c(0,1) = 5.0;
  c(1,0) = 7.0;

  auto newu = c(0,0)*w[0] + c(0,1)*w[1];
  auto newv = c(1,0)*w[0] + c(1,1)*w[1];

  auto r = transform(world, w, c, 0.0, true);
  r[0].verify_tree();
  r[1].verify_tree();

  //print("new");
  //r[0].print_tree();
  //print("old");
  //newu.print_tree();

  print("newu", (newu-r[0]).norm2());
  print("newv", (newv-r[1]).norm2());

  finalize();

  return 0;
}

