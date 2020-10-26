
//
// Created by adrian on 4/19/20.
//
/*
 *    Written by: bsundahl
 *    Date: A long time ago...
 *
 */

#include <stdlib.h>

#include "TDDFT.h"  // All response functions/objects enter through this
#include "adrian/density.h"
#include "adrian/property_functions.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && \
    defined(HAVE_UNISTD_H)
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <unistd.h>
static inline int file_exists(const char* inpname) {
  struct stat buffer;
  int rc = stat(inpname, &buffer);
  return (rc == 0);
}

#endif

int main(int argc, char** argv) {
  // Initialize MADNESS mpi
  initialize(argc, argv);
  World world(SafeMPI::COMM_WORLD);
  startup(world, argc, argv);

  // This makes a default input file name of 'input'
  const char* input = "input";
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      input = argv[i];
      break;
    }
  }
  if (!file_exists(input)) throw "input file not found";

  // Create the TDHF object

  FirstOrderDensity densityTest(world, input);
  //
  densityTest.PlotResponseDensity(world);
  densityTest.PrintDensityInformation();

  ResponseParameters Rparams = densityTest.GetResponseParameters();
  if (Rparams.property) {  //
    Property nuclear_operator(world, "nuclear", densityTest.GetMolecule());
    Property dipole_operator(world, "dipole");

    Tensor<double> alpha = densityTest.ComputeSecondOrderPropertyTensor(world);
    PrintSecondOrderAnalysis(world, alpha, densityTest.GetFrequencyOmega());
  }

  world.gop.fence();
  finalize();

  return 0;
}
