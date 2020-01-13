/*
 *    Written by: bsundahl and jscanderson
 *    Date: A long time ago...
 *
 */ 


#include "DF.h"    // All response functions/objects enter through this
#include <stdlib.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
static inline int file_exists(const char * inpname)
{
     struct stat buffer;
     int rc = stat(inpname, &buffer);
     return (rc==0);
}
#endif

int main(int argc, char** argv)
{
     // Initialize MADNESS mpi 
     initialize(argc, argv);
     World world(SafeMPI::COMM_WORLD);
     startup(world,argc,argv);
     std::cout.precision(10);

     // This makes a default input file name of 'input'   
     const char * input = "input";
     for (int i=1; i<argc; i++){
          if(argv[i][0] != '-'){
               input = argv[i];
               break;
          }
     }
     if (!file_exists(input)) throw "input file not found";

     // Create the TDHF object
     DF my_calc(world, input);       
                 
     // Have it iterate to convergence
     my_calc.solve(world);
     //my_calc.virtuals(world);

     world.gop.fence();
     world.gop.fence();
     finalize();

     return 0; 
}

//kthxbye
