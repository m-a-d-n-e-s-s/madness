#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>

using namespace madness;

static const double L = 32.0;   // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-4; // precision

static const size_t nfunc = 4; // number of functions

// A class that behaves like a function to compute a Gaussian of given origin and exponent
class Gaussian : public FunctionFunctorInterface<double,3> {
public:
    const coord_3d center;
    const double exponent;
    const double coefficient;
    std::vector<coord_3d> specialpt;

    Gaussian(const coord_3d& center, double exponent, double coefficient)
        : center(center), exponent(exponent), coefficient(coefficient), specialpt(1)
    {
        specialpt[0][0] = center[0];
        specialpt[0][1] = center[1];
        specialpt[0][2] = center[2];
    }

    // MADNESS will call this interface
    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        for (int i=0; i<3; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    }

    // By default, adaptive projection into the spectral element basis
    // starts uniformly distributed at the initial level.  However, if
    // a function is "spiky" it may be necessary to project at a finer
    // level but doing this uniformly is expensive.  This method
    // enables us to tell MADNESS about points/areas needing deep
    // refinement (the default is no special points).
    std::vector<coord_3d> special_points() const {
        return specialpt;
    }
};

// Makes a new square-normalized Gaussian functor with random origin and exponent
real_functor_3d random_gaussian() {
    const double expntmin=1e-1;
    const double expntmax=1e4;
    const real_tensor& cell = FunctionDefaults<3>::get_cell();
    coord_3d origin;
    for (int i=0; i<3; i++) {
        origin[i] = RandomValue<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(expntmin);
    double hi = log(expntmax);
    double expnt = exp(RandomValue<double>()*(hi-lo) + lo);
    double coeff = pow(2.0*expnt/constants::pi,0.75);
    return real_functor_3d(new Gaussian(origin,expnt,coeff));
}

// Makes a vector of new square-normalized Gaussian functions with random origin and exponent
std::vector<real_function_3d> random_gaussians(size_t n, World& world) {
    std::vector<real_function_3d> result(n);
    for (size_t i=0; i<n; i++) {
        result[i] = FunctionFactory<double,3>(world).functor(random_gaussian());
    }
    return result;
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

  default_random_generator.setstate(99); // Ensure all processes have the same state

  // Create a vector of random Gaussian functions
  std::vector<real_function_3d> v = random_gaussians(nfunc, world);

  truncate(world, v);

  // Generate a random tensor of coefficients
  Tensor<double> c(nfunc,nfunc);
  for (size_t i=0; i<nfunc; i++) {
      for (size_t j=0; j<nfunc; j++) {
          c(i,j) = RandomValue<double>()-0.5;
      }
  }

  print("c");
  print(c);

  auto rold = transform(world, v, c); // old transform uses gaxpy and only zero-based sparsity

  // manual loop no sparsity or asynchronous execution
  // std::vector<real_function_3d> rold = zero_functions_compressed<double,3>(world, v.size());
  // for (size_t i=0; i<nfunc; i++) {
  //     for (size_t j=0; j<nfunc; j++) {
  //         rold[i] += v[j] * c(j,i);
  //     }
  // }
  
  auto r = transform(world, v, c, 0.0, true); // new transform uses mxm and thresholding
  
  for (const auto& f : r) {
      f.verify_tree();
  }

  // compute the norm of the errors for each component
  for (size_t i=0; i<nfunc; i++) {
      print("error", i, (rold[i]-r[i]).norm2());
  }

  finalize();

  return 0;
}

