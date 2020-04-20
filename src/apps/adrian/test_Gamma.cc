
#include "ResponseFunction.h"
#include "DerivCoulomb.h"
#include "DerivExchange.h"
#include <madness/mra/mra.h>
//#include <madness/mra/operator.h>
//#include <madness/constants.h>
//#include <madness/mra/nonlinsol.h>  // The kain solver
#include <vector>
#include <memory>
//#include <math.h>


using namespace madness;

int main(int argc, char** argv)
{
   // Initialize MADNESS mpi
   initialize(argc, argv);
   World world(SafeMPI::COMM_WORLD);
   startup(world, argc, argv);

   auto orbitals = zero_functions_compressed<double,3>(world, 1);
   auto orbitals_ptr = std::make_shared<std::vector<Function<double,3>>>(orbitals);

   // Testing J constructor 
   DerivCoulomb<double,3> dJ(world, 1e-6, 1e-6, orbitals_ptr);

   // Testing K constructor
   DerivExchange<double,3> dK(world, 1e-6, 1e-6, orbitals_ptr, true);

   world.gop.fence();
   finalize();

   return 0;
}
