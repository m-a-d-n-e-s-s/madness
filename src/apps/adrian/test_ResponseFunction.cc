
#include "ResponseFunction.h"
#include <memory>
#include <madness/mra/mra.h>
//#include <madness/mra/operator.h>
//#include <madness/constants.h>
//#include <madness/mra/nonlinsol.h>  // The kain solver
//#include <vector>
//#include <math.h>


using namespace madness;


template<std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double,NDIM> {
    typedef Vector<double,NDIM> coordT;

public:
    GaussianGuess(const coordT& origin, const double alpha,
            const std::vector<int> ijk=std::vector<int>(NDIM))
            : origin(origin), exponent(alpha), ijk(ijk) {
    }

    coordT origin;
    double exponent;        ///< exponent of the guess
    std::vector<int> ijk;   ///< cartesian exponents

    double operator()(const coordT& xyz) const {
        double arg=0.0, prefac=1.0;
        for (std::size_t i=0; i<NDIM;++i) {
            arg+=(xyz[i]-origin[i])*(xyz[i]-origin[i]);
            prefac*=pow(xyz[i],ijk[i]);
        }
        const double e=exponent*arg;
        return prefac*exp(-e);
    }
};


int main(int argc, char** argv)
{
   // Initialize MADNESS mpi
   initialize(argc, argv);
   World world(SafeMPI::COMM_WORLD);
   startup(world, argc, argv);

   // Testing constructor 
   ResponseFunction<double,3> x(world, 1,1);
   ResponseFunction<double,3> y(world, 2,1);
   ResponseFunction<double,3> z(world, 2,2);

   // Testing member variables
   if(world.rank() == 0) print("\n size of x:", x.r_states, "by", x.g_states);   
   if(world.rank() == 0) print(" size of y:", y.r_states, "by", y.g_states);   
   if(world.rank() == 0) print(" size of z:", z.r_states, "by", z.g_states);   
   if(world.rank() == 0) print("");

    // Testing accessor
   if(world.rank() == 0) print(" Testing accessor.");

   // Give each ResponseFunction a gaussian
   x(0,0) = FunctionFactory<double,3>(world).functor(std::shared_ptr<FunctionFunctorInterface<double,3>>(new GaussianGuess<3>(Vector<double,3>{0,0,0}, 1.0, std::vector<int>{0,0,0}))); 
   y(1,0) = FunctionFactory<double,3>(world).functor(std::shared_ptr<FunctionFunctorInterface<double,3>>(new GaussianGuess<3>(Vector<double,3>{0,0,0}, 1.0, std::vector<int>{0,0,0}))); 
   z(1,1) = FunctionFactory<double,3>(world).functor(std::shared_ptr<FunctionFunctorInterface<double,3>>(new GaussianGuess<3>(Vector<double,3>{0,0,0}, 1.0, std::vector<int>{0,0,0}))); 

   // Testing print_norms()
   if(world.rank() == 0) print(" norms of x:");
   x.print_norms();
   if(world.rank() == 0) print(" norms of y:");
   y.print_norms();
   if(world.rank() == 0) print(" norms of z:");
   z.print_norms();
   if(world.rank() == 0) print("");

   // Testing same_size()
   if(world.rank() == 0) print(" x same size as x?", x.same_size(x));
   if(world.rank() == 0) print(" x same size as y?", x.same_size(y));
   if(world.rank() == 0) print(" y same size as z?", y.same_size(z));
   if(world.rank() == 0) print("");

   // Testing operator+
   try {
      x = x + x;
      if(world.rank() == 0) print(" Successfully added x + x");
      if(world.rank() == 0) print(" New norms:");
      x.print_norms();
   }
   catch (...) {
      if(world.rank() == 0) print(" Error: should be able to add same size ResponseFunctions");
   }
   try {
      x = x + y;
   }
   catch (...) {
      if(world.rank() == 0) print(" Can't add things not of the same dimensions.");
      if(world.rank() == 0) print("    Tried adding x + y.");
      if(world.rank() == 0) print("    Size of x:", x.r_states, "by", x.g_states);
      if(world.rank() == 0) print("    Size of y:", y.r_states, "by", y.g_states);
   }
   if(world.rank() == 0) print("");

   // Testing operator-
   try {
      x = x - x;
      if(world.rank() == 0) print(" Successfully subtracted x - x");
      if(world.rank() == 0) print(" New norms:");
      x.print_norms();
   }
   catch (...) {
      if(world.rank() == 0) print(" Error: should be able to subtract same size ResponseFunctions");
   }
   try {
      x = x - y;
   }
   catch (...) {
      if(world.rank() == 0) print(" Can't subtract things not of the same dimensions.");
      if(world.rank() == 0) print("    Tried subtacting x - y.");
      if(world.rank() == 0) print("    Size of x:", x.r_states, "by", x.g_states);
      if(world.rank() == 0) print("    Size of y:", y.r_states, "by", y.g_states);
   }
   if(world.rank() == 0) print("");

   world.gop.fence();
   finalize();

   return 0;
}
