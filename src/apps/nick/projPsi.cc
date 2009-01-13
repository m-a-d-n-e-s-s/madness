//\file projPsi.cc
//\brief Projects a time evolved wave function onto some bound states
/******************************************
 * By: Nick Vence
 * This code must handled with care for two reasons:
 * First it uses the Gnu Scientific Library
 * http://www.gnu.org/software/gsl/
 * Second it is loading a function from disk?  NO
 * If YES ensure the following
 * 1) k (The wavelet order) is the same as the projected functions: see main()
 * 2) functions to be loaded are in the directory
 *****************************************/
#include <mra/mra.h>
#include <complex>
#include <string>
#include <fstream>
using std::ofstream;
#include "wavef.h"
#define PRINTLINE(str) if(world.rank()==0) cout << str << endl

using namespace madness;

const int nIOProcessors =1;
const string prefix = "data";
typedef std::complex<double> complexd;
typedef Vector<double,NDIM> vector3D;
typedef Function<complexd,NDIM> complex_functionT;
typedef FunctionFactory<complexd,NDIM> complex_factoryT;

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

static complexd gaussian(const vector3D& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-x*x -y*y -z*z);
}

static complexd zgaussian(const vector3D& r) {
    const double x=r[0], y=r[1], z=r[2];
    return x * exp(-x*x -y*y -z*z);
}

void doWork(World& world) {

  PRINTLINE("Creating gaussians");
  complex_functionT gauss  = complex_factoryT(world).f(gaussian);
  complex_functionT zgauss = complex_factoryT(world).f(zgaussian);

  PRINTLINE("Creating 1s 2s and three 2p basis functions");
  complex_functionT psi100 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 1.0, 1,0,0)));
  complex_functionT psi200 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 1.0, 2,0,0)));
  complex_functionT psi210 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 1.0, 2,1,0)));

  //The complex numbers are giving us an error in a gsl file
  complex_functionT psi211 = FunctionFactory<complexd,NDIM>(world).
    functor(functorT( new BoundWF(1.0, 1.0, 2,1,1)));
  complex_functionT psi21_1 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new BoundWF(1.0, 1.0, 2,1,-1)));


  PRINTLINE("Creating |k= 0.52>");
  double k[] = {0.52, 0, 0};
  Vector<double,NDIM> k52(k);
  Function<complexd,NDIM> psi_52 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new ScatteringWF(1.0, 1.0, k52)));

  PRINTLINE("Projecting the gaussians onto the basis functions");
  PRINTLINE("<gauss |100>   = " << gauss.inner(psi100) );
  PRINTLINE("<gauss |200>   = " << gauss.inner(psi200) );
  PRINTLINE("<gauss |210>   = " << gauss.inner(psi210) );
  PRINTLINE("<gauss |k=1.0> = " << gauss.inner(psi_52) );
  PRINTLINE("<zgauss|100>   = " << zgauss.inner(psi100) );
  PRINTLINE("<zgauss|200>   = " << zgauss.inner(psi200) );
  PRINTLINE("<zgauss|210>   = " << zgauss.inner(psi210) );
  PRINTLINE("<zgauss|k=1.0> = " << zgauss.inner(psi_52) );


  /*
  PRINTLINE("Testing our capacity to load a wave function from disk");
  int step = 0;
  if(wave_function_exists(world,step)) {
    PRINTLINE("wave_function_exists = true");
    Function<complexd, NDIM> loadedFunc = wave_function_load(world, step);
    PRINTLINE("<data|100> =  " << loadedFunc.inner(psi100) );
    PRINTLINE("<data|200> =  " << loadedFunc.inner(psi200) );
    PRINTLINE("<data|210> =  " << loadedFunc.inner(psi210) );
  } else PRINTLINE("LoadedFunc doesn't exist");
  
  k = {0.50, 0, 0};
  Vector<double,NDIM> k50(ka);
  k = {0.54, 0, 0};
  Vector<double,NDIM> k54(ka);
  k = {0.48, 0, 0};
  Vector<double,NDIM> k48(ka);
  k = {0.56, 0, 0};
  Vector<double,NDIM> k56(ka);

  PRINTLINE("Creating |k= 0.50>");
  Function<complexd,NDIM> psi_50 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new ScatteringWF(1.0, 1.0, k50)));
  PRINTLINE("Creating |k= 0.54>");
  Function<complexd,NDIM> psi_54 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new ScatteringWF(1.0, 1.0, k54)));
  PRINTLINE("Creating |k= 0.48>");
  Function<complexd,NDIM> psi_48 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new ScatteringWF(1.0, 1.0, k48)));
  PRINTLINE("Creating |k= 0.56>");
  Function<complexd,NDIM> psi_56 = FunctionFactory<complexd,NDIM>(world).
      functor(functorT( new ScatteringWF(1.0, 1.0, k56)));
  */

  /* Load psi(+) for each of the cut parameters
   */

  ///* COMPLEX FUNCTIONS GIVE AN ERROR  
  PRINTLINE("<gauss|211> = " << gauss.inner(psi211) );
  PRINTLINE("<gauss|21_1> = " << gauss.inner(psi21_1) );
  //*/
}


int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  // Setup defaults for numerical functions
  FunctionDefaults<NDIM>::set_k(12);              // Wavelet order
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
