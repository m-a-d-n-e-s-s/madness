#define CATCH_CONFIG_RUNNER
#include "apps/external_headers/catch.hpp"
#include "apps/chem/SCF.h"
#include "madness/world/worldmem.h"

#include "TDDFT.h"  // All response functions/objects enter through this
#include "molresponse/density.h"


int main(int argc, char *argv[]) {
  // initialize(argc, argv);// initializes a world argument with argc and argv
  // World world(SafeMPI::COMM_WORLD);
  // startup(world, argc, argv, true);
  try {
    int result = Catch::Session().run(argc, argv);
    return result;
  } catch (const SafeMPI::Exception &e) {
    print(e);
    error("caught an MPI exception");
  } catch (const madness::MadnessException &e) {
    print(e);
    error("caught a MADNESS exception");
  } catch (const madness::TensorException &e) {
    print(e);
    error("caught a Tensor exception");
  } catch (const char *s) {
    print(s);
    error("caught a string exception");
  } catch (const std::string &s) {
    print(s);
    error("caught a string (class) exception");
  } catch (const std::exception &e) {
    print(e.what());
    error("caught an STL exception");
  } catch (...) {
    error("caught unhandled exception");
  }

  // print_meminfo(world.rank(), "startup");
  // std::cout.precision(6);
  // print_stats(world);
}
