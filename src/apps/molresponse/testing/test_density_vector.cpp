//
// Created by adrianhurtado on 1/1/22.
//

#include "molresponse/ResponseExceptions.hpp"
#include "molresponse/densityVector.hpp"
#include "molresponse/global_functions.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <unistd.h>

static inline int file_exists(const char *input_name) {
  struct stat buffer {};
  size_t rc = stat(input_name, &buffer);
  return (rc == 0);
}

#endif


TEST_CASE("Test Density Vectors", "Testing Basic Functionality of Response Functions") {
  int argc = 1;
  char **argv = nullptr;
  initialize(argc, argv);  // initializes a world argument with argc and argv
  {
    World world(SafeMPI::COMM_WORLD);
    startup(world, 1, nullptr,
            true);  // TODO: ask Robert about proper startup and implement
    std::cout.precision(6);

    const char *input_file = "response.in";
    if (world.rank() == 0) {
      print("input filename: ", input_file);
    }
    if (!file_exists(input_file)) {
      throw Input_Error{};
    }

    auto calc_params = initialize_calc_params(world, std::string(input_file));
    auto &[ground_calculation, molecule, response_parameters] = calc_params;
    vecfuncT ground_orbitals = ground_calculation.orbitals();

    response_parameters.print("ResponseParameters", "Not sure footer");




        // Everything goes in here

        finalize();
  }
}
