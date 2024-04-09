//
// Created by adrianhurtado on 1/1/22.
//

#include "madness/external/catch/catch.hpp"
#include "madness/tensor/tensor_json.hpp"
#include "response_functions.h"
#include "response_parameters.h"
#include "string"
#include "timer.h"
#include "x_space.h"

unsigned int Factorial(unsigned int number) {
  return number <= 1 ? number : Factorial(number - 1) * number;
}

TEST_CASE("X_space", "[copy]") {
  unsigned int m = 5;
  unsigned int n = 4;
  int argc = 1;
  char** argv = nullptr;
  madness::initialize(argc, argv);  // initializes a world argument with argc and argv
  {
    madness::World world(SafeMPI::COMM_WORLD);
    startup(world,
            1,
            nullptr,
            true);  // TODO: ask Robert about proper startup and implement
    molresponse::start_timer(world);
    madness::X_space v{world, m, n};
    REQUIRE(v.num_states() == m);
    // world.gop.fence();
    REQUIRE(v.num_orbitals() == n);
    molresponse::end_timer(world, "basic x space test");
    madness::finalize();
  }
}

TEST_CASE("Json Testing", "Simple JSON") {
  // create an empty structure (null)
  nlohmann::json j;
  // add a number that is stored as double (note the implicit conversion of j to
  // an object)
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
  j["list"] = {1, 0, 2};
  // add another object (using an initializer list of pairs)
  j["object"] = {{"currency", "USD"}, {"value", 42.99}};
  // instead, you could also write (which looks very similar to the JSON above)
  nlohmann::json j2 = {{"pi", 3.141},
             {"happy", true},
             {"name", "Niels"},
             {"nothing", nullptr},
             {"answer", {{"everything", 42}}},
             {"list", {1, 0, 2}},
             {"object", {{"currency", "USD"}, {"value", 42.99}}}};
  REQUIRE(j == j2);
  std::ofstream ofs("j_example.json");
  ofs << j << std::endl;
}

template <typename T>
bool operator==(const madness::Tensor<T>& a, const madness::Tensor<T>& b) {
  // check the size
  if (a.size() != b.size()) {
    return false;
  }  // check the number of dimensions
  if (a.ndim() != b.ndim()) {
    return false;
  }
  // check if the dimension sizes are the same
  auto a_dims = a.dims();
  auto b_dims = b.dims();
  bool dims_equal = std::equal(a_dims, a_dims + a.ndim(), b_dims);
  auto a_flat = a.flat();
  auto b_flat = b.flat();
  bool val_equal = std::equal(&a_flat[0], &a_flat[0] + a.size(), &b_flat[0]);

  return dims_equal && val_equal;
}

TEST_CASE("Json Testing 2", "Serialization") {
  madness::Tensor<double> a(3, 3);
  madness::Tensor<double> b(3, 2, 4);
  madness::Tensor<double> c(1, 3);
  a.fillrandom();
  b.fillindex();
  c.fill(3);

  nlohmann::json j_a = tensor_to_json(a);
  nlohmann::json j_b = tensor_to_json(b);
  nlohmann::json j_c = tensor_to_json(c);

  // How can I make this an automatic template function?
  auto a_copy = madness::tensor_from_json<double>(j_a);
  auto b_copy = madness::tensor_from_json<double>(j_b);
  auto c_copy = madness::tensor_from_json<double>(j_c);

  REQUIRE(a == a_copy);
  REQUIRE(b == b_copy);
  REQUIRE(c == c_copy);
}

TEST_CASE("Json Testing 3", "Json Tensor Indexing") {
  madness::Tensor<double> a(3, 3);
  madness::Tensor<double> b(3, 2, 4);
  madness::Tensor<double> c(1, 3);
  a.fillrandom();
  b.fillindex();
  c.fill(3);
  nlohmann::json j = {};

  nlohmann::json j_a = tensor_to_json(a);
  nlohmann::json j_b = tensor_to_json(b);
  nlohmann::json j_c = tensor_to_json(c);

  for (int i = 0; i < 3; i++) {
    nlohmann::json j_iter = {};
    j_iter["iter"] = i;
    j_iter["a"] = j_a;
    j_iter["b"] = j_b;
    j_iter["c"] = j_c;
    j.push_back(j_iter);
  }

  std::ofstream ofs("j_iters.json");
  ofs << j << std::endl;
  std::cout << j << std::endl;
  // How can I make this an automatic template function?
}

TEST_CASE("print_QCSchema Test ", "Json Tensor Indexing") {

  madness::vec_pair_ints int_vals;
  madness::vec_pair_T<double> double_vals;
  madness::vec_pair_tensor_T<double> double_tensor_vals;

  nlohmann::json j={};

  int_vals.push_back({"a", 4});
  int_vals.push_back({"b", 5});
  int_vals.push_back({"c", 4});
  int_vals.push_back({"d", 6});

  madness::to_json(j,int_vals);

  double_vals.push_back({"aa", 4});
  double_vals.push_back({"bb", 5});
  double_vals.push_back({"cc", 4});
  double_vals.push_back({"dd", 6});

  madness::to_json(j,double_vals);

  madness::Tensor<double> t_a(3, 3);
  madness::Tensor<double> t_b(3, 2, 4);
  madness::Tensor<double> t_c(1, 3);

  double_tensor_vals.push_back({"t_a",t_a});
  double_tensor_vals.push_back({"t_b",t_b});
  double_tensor_vals.push_back({"t_c",t_c});

  to_json(j,double_tensor_vals);

  madness::output_schema( "test_schema", j);

}
TEST_CASE("Response Parameters Test ", "Testing parameters to_json") {

  int argc = 1;
  char** argv = nullptr;
  madness::initialize(argc, argv);  // initializes a world argument with argc and argv
  {
    madness::World world(SafeMPI::COMM_WORLD);
    startup(world,
            1,
            nullptr,
            true);  // TODO: ask Robert about proper startup and implement
    molresponse::start_timer(world);
    // Everything goes in here
    madness::ResponseParameters r_params;
    nlohmann::json j;
    r_params.to_json(j);
    world.gop.fence();
    std::cout<<j<<std::endl;
    molresponse::end_timer(world, "basic Response parameters testing");
  }


}
