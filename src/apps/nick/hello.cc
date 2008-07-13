/// hello.cc
/// compute the 3D integral over r

#include <mra/mra.h>
#include <complex>
#include <string>
#include <fstream>
using std::ofstream;
#include <nick/wavef.h>

using namespace madness;

const string prefix = "hello";
typedef std::complex<double> complexd;
typedef Function<complexd,NDIM> complex_functionT;

void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);
const char* wave_function_filename(int step);

void wave_function_store(World& world, int step, const complex_functionT& psi) {
	int nIOProcessors =1;
    ParallelOutputArchive ar(world, wave_function_filename(step), nIOProcessors);
    ar & psi;
}

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
    return fname;
}

complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  
  // Setup defaults for numerical functions
  FunctionDefaults<NDIM>::set_k(12);             // Wavelet order
  FunctionDefaults<NDIM>::set_thresh(1e-2);       // Accuracy
  FunctionDefaults<NDIM>::set_cubic_cell(-20.0, 20.0);

  // Testing our capacity to load a wave function
  Function<complexd,NDIM> psi200 = FunctionFactory<complexd,NDIM>(world).
	  functor(functorT( new BoundWF(1.0, 1.0, 2,0,0)));
  wave_function_store(world, 1, psi200);
  Function<complexd, NDIM> loadedFunc = wave_function_load(world, 1);
  cout << "We are dotting the wave function with itself " 
	  << psi200.inner(loadedFunc) << endl;
  
    MPI::Finalize();				//FLAG
  return 0;
}
