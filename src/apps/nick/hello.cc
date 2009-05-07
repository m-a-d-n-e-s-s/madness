/// hello.cc
/// We are projecting a time evolved wave function onto some bound states

#include <mra/mra.h>
#include <complex>
#include <string>
#include <fstream>
using std::ofstream;
#include <nick/wavef.h>

#define PRINT(str) if(world.rank()==0) cout << str
#define PRINTLINE(str) if(world.rank()==0) cout << str << endl

using namespace madness;

const int nIOProcessors =1;
const string prefix = "data";
typedef std::complex<double> complexd;
typedef Function<complexd,NDIM> complex_functionT;

const char* wave_function_filename(int step);
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);


const char* wave_function_filename(int step) {
  static char fname[1024];
  sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
  return fname;
}

bool wave_function_exists(World& world, int step) {
  return ParallelInputArchive::exists(world, wave_function_filename(step));
}

void wave_function_store(World& world, int step, const complex_functionT& psi) {
  ParallelOutputArchive ar(world, wave_function_filename(step), nIOProcessors);
  ar & psi;
}

complex_functionT wave_function_load(World& world, int step) {
  complex_functionT psi;
  ParallelInputArchive ar(world, wave_function_filename(step));
  ar & psi;
  return psi;
}


void doWork(World& world) {
  PRINTLINE("Creating three basis functions");
  Function<complexd,NDIM> psi100 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 1,0,0)));
  Function<complexd,NDIM> psi200 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 2,0,0)));
  Function<complexd,NDIM> psi210 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 2,1,0)));

  int step = 0;
  PRINTLINE("Testing our capacity to load a wave function from disk");
  if(wave_function_exists(world,step)) {
    PRINTLINE("wave_function_exists = true");
    Function<complexd, NDIM> loadedFunc = wave_function_load(world, step);
    PRINT("<data|100> =  ") << loadedFunc.inner(psi100) << endl;
    PRINT("<data|200> =  ") << loadedFunc.inner(psi200) << endl;
    PRINT("<data|210> =  ") << loadedFunc.inner(psi210) << endl;
  } else PRINTLINE("LoadedFunc doesn't exist");
}


int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  // Setup defaults for numerical functions
  FunctionDefaults<NDIM>::set_k(8);             // Wavelet order
  FunctionDefaults<NDIM>::set_thresh(1e-3);       // Accuracy
  FunctionDefaults<NDIM>::set_cubic_cell(-20.0, 20.0);

  try {
    doWork(world);
  } catch (const MPI::Exception& e) {
    //print(e);
    error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
    print(e);
    error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
    print(e);
    error("caught a Tensor exception");
  } catch (const char* s) {
    print(s);
    error("caught a c-string exception");
  } catch (char* s) {
    print(s);
    error("caught a c-string exception");
  } catch (const std::string& s) {
    print(s);
    error("caught a string (class) exception");
  } catch (const std::exception& e) {
    print(e.what());
    error("caught an STL exception");
  } catch (...) {
    error("caught unhandled exception");
  }

  MPI::Finalize();				//FLAG
  return 0;
}
