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
  TDHF my_calc(world, input);

  // Check if calculating a property
  if (my_calc.Rparams.property) {
    // compute frequency density response will compute given
    // the property asked for
    // first computes the rhs for
    // dipole
    // nuclear
    // order2
    // order3
    my_calc.compute_freq_density(world);
    // Future properties go here
    //
    // If not a property, just calculate response states
  } else {
    my_calc.solve(world);
  }

  world.gop.fence();
  finalize();

  return 0;
}
