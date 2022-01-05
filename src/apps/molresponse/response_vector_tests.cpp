//
// Created by adrianhurtado on 1/1/22.
//
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/json.hpp"
#include "response_functions.h"
#include "timer.h"
#include "x_space.h"

unsigned int Factorial(unsigned int number) { return number <= 1 ? number : Factorial(number - 1) * number; }


using namespace madness;

using json=nlohmann::json;

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


TEST_CASE("Json Testing", "Simple JSON") {

  // create an empty structure (null)
  json j;
  // add a number that is stored as double (note the implicit conversion of j to an object)
  j["pi"] = 3.141;
  // add a Boolean that is stored as bool
  j["happy"] = true;
  // add a string that is stored as std::string
  j["name"] = "Niels";
  // add another null object by passing nullptr
  j["nothing"] = nullptr;
  // add an object inside the object
  j["answer"]["everything"] = 42;
  // add an array that is stored as std::vector (using an initializer list)
  j["list"] = { 1, 0, 2 };
  // add another object (using an initializer list of pairs)
  j["object"] = { {"currency", "USD"}, {"value", 42.99} };
  // instead, you could also write (which looks very similar to the JSON above)
  json j2 = {
      {"pi", 3.141},
      {"happy", true},
      {"name", "Niels"},
      {"nothing", nullptr},
      {"answer", {
                     {"everything", 42}
                 }},
      {"list", {1, 0, 2}},
      {"object", {
                     {"currency", "USD"},
                     {"value", 42.99}
                 }}
  };
  REQUIRE(j==j2);
  std::cout<<j<<std::endl;
}
TEST_CASE("Json Testing 2", "Serialization") {




}


