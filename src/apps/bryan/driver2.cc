#include "TDA2.h"
#include <stdlib.h>

int main(int argc, char** argv)
{
   // Initialize MADNESS mpi 
   initialize(argc, argv);
   World world(SafeMPI::COMM_WORLD);
   startup(world,argc,argv);
   

   // Create the TDA object
   TDA my_calc(world,           // Communicator
               argv[1],         // Contains restart file
               atoi(argv[2]),   // Number of requested states
               atof(argv[3]),   // Energy range of acceptable orbitals to excite from (energy from HOMO)
               atoi(argv[4]));  // Print level

   // Have it iterate to convergence
   my_calc.solve(world);

   world.gop.fence();
   finalize();

   return 0; 
}
