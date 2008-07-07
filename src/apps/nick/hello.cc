/// volumeINT.cc
/// compute the 1D integral over r

#include <mra/mra.h>
#include <complex>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_hyperg.h>
#include <fstream>
using std::ofstream;
#include "nick/wavef.h"

using namespace madness;

typedef std::complex<double> complexd;


int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  
  // Setup defaults for numerical functions
  FunctionDefaults<NDIM>::set_k(11);                // Wavelet order
  FunctionDefaults<NDIM>::set_thresh(1e-7);       // Accuracy
  FunctionDefaults<NDIM>::set_cubic_cell(-30.0, 30.0);

  f11Tester(world);
  MPI::Finalize();				//FLAG
  return 0;
}
