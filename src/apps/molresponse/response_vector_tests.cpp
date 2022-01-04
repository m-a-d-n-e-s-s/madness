//
// Created by adrianhurtado on 1/1/22.
//
#include "catch.hpp"
#include "response_functions.h"
#include "timer.h"
#include "x_space.h"

unsigned int Factorial(unsigned int number) { return number <= 1 ? number : Factorial(number - 1) * number; }

TEST_CASE("Factorials are computed", "[factorial]") {
  REQUIRE(Factorial(1) == 1);
  REQUIRE(Factorial(2) == 2);
  REQUIRE(Factorial(3) == 6);
  REQUIRE(Factorial(10) == 3628800);
}

using namespace madness;

TEST_CASE("X_space", "[copy]") {
  unsigned int m = 5;
  unsigned int n = 4;
  int argc = 1;
  char** argv = nullptr;
  initialize(argc, argv);  // initializes a world argument with argc and argv
  {
    World world(SafeMPI::COMM_WORLD);
    startup(world, 1, nullptr, true);  // TODO: ask Robert about proper startup and implement
    molresponse::start_timer(world);
    X_space v{world, m, n};
    REQUIRE(v.num_states() == m);
    // world.gop.fence();
    REQUIRE(v.num_orbitals() == n);
    molresponse::end_timer(world, "basic x space test");
    finalize();
  }
}
