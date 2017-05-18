#include "TDA.h"
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
               atoi(argv[3]));  // Print level

   // Have it iterate to convergence
   my_calc.solve(world);

   finalize();
   return 0; 
}
