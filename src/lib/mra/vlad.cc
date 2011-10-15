#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/vmra.h>
#include <mra/lbdeux.h>
#include <constants.h>
using namespace madness;

static const int NFUNC = 10;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

// A class that behaves like a function to compute a Gaussian of given origin and exponent
class Gaussian : public FunctionFunctorInterface<double,3> {
public:
    const coord_3d center;
    const double exponent;
    const double coefficient;
    std::vector<coord_3d> specialpt;

    Gaussian(const coord_3d& center, double exponent, double coefficient)
        : center(center), exponent(exponent), coefficient(coefficient), specialpt(1,center)
    {}

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

// This structure is used to estimate the cost of computing on a block of coefficients
// The constructor saves an estimate of the relative cost of computing on
// leaf (scaling function) or interior (wavelet) coefficients.  Unless you 
// are doing nearly exclusively just one type of operation the speed is
// not that sensitive to the choice, so the default is usually OK.
//
// The operator() method is invoked for each block of coefficients
// (function node) to estimate the cost. Since pretty much everything
// involves work at the top of the tree we give levels 0 and 1 a 100x
// multiplier to avoid them being a bottleneck.
struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0) 
        : leaf_value(leaf_value)
        , parent_value(parent_value) 
    {}

    double operator()(const Key<3>& key, const FunctionNode<double,3>& node) const {
        if (key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

void test(World& world, bool doloadbal=false) {
    double start;
    vector_real_function_3d f(NFUNC);

    default_random_generator.setstate(99); // Ensure all processes have the same state

    // By default a sychronization (fence) is performed after projecting a new function
    // but if you are projecting many functions this is inefficient.  Here we turn
    // off the fence when projecting, and do it manually just once.  
    start = wall_time();
    for (int i=0; i<NFUNC; i++) 
        f[i] = real_factory_3d(world).functor(random_gaussian()).nofence();
    world.gop.fence();
    double projection = wall_time() - start;

    start = wall_time();
    compress(world, f);
    double compression = wall_time() - start;

    start = wall_time();
    reconstruct(world, f);
    double reconstruction = wall_time() - start;
    if (world.rank() == 0) printf("project %.2f compress %.2f reconstruct %.2f\n",
                                  projection, compression, reconstruction);
}

int main(int argc, char** argv) {
  // Initialize the parallel programming environment
  initialize(argc,argv);
  World world(MPI::COMM_WORLD);

  // Load info for MADNESS numerical routines
  startup(world,argc,argv);

  // Set/override defaults
  FunctionDefaults<3>::set_cubic_cell(-20,20);
  FunctionDefaults<3>::set_apply_randomize(false);
  FunctionDefaults<3>::set_project_randomize(false);
  FunctionDefaults<3>::set_truncate_on_project(false); // 8x more work but good for Vlad
  FunctionDefaults<3>::set_k(8);
  FunctionDefaults<3>::set_thresh(1e-6);

  START_TIMER;

  test(world);
  test(world);
  test(world);

  END_TIMER("vlad_tests");

  finalize();

  return 0;
}
