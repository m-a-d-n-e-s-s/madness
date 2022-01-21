//
// Created by adrianhurtado on 1/1/22.
//

#include "apps/external_headers/tensor_json.hpp"
#include "response_functions.h"
#include "timer.h"
#include "string"
#include "x_space.h"

TEST_CASE("Response Parameters Test ", "Testing parameters to_json") {

  int argc = 1;
  char** argv = nullptr;
  initialize(argc, argv);  // initializes a world argument with argc and argv
  {
    World world(SafeMPI::COMM_WORLD);
    startup(world,
            1,
            nullptr,
            true);  // TODO: ask Robert about proper startup and implement
    molresponse::start_timer(world);
    // Everything goes in here
    ResponseParameters r_params;
    json j;
    r_params.to_json(j);
    world.gop.fence();
    std::cout<<j<<endl;
    molresponse::end_timer(world, "basic Response parameters testing");
    finalize();
  }


}
