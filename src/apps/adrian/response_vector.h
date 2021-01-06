
#ifndef SRC_APPS_ADRIAN_RESPONSE_VECTOR_H_
#define SRC_APPS_ADRIAN_RESPONSE_VECTOR_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

namespace madness {

template <typename T> struct response_vector {
  size_t num_orbitals;
  std::vector<Function<T, 3>> x, y;
  typedef Function<T, 3> function_T;
  typedef T value_type;

  // copy constructor
  explicit response_vector(const response_vector& a)
      : num_orbitals(a.num_orbitals), x(madness::copy(a.x[0].world(), a.x)),
        y(madness::copy(a.y[0].world(), a.y)) {}
  // conversion from 2 vector functions to response_vector
  explicit response_vector(const function_T x_in, const function_T y_in) {
    MADNESS_ASSERT(x_in.size() == y_in.size());
    num_orbitals = x_in.size();
    x = x_in;
    y = y_in;
  }
  // default constructor
  response_vector() {}
  // destructor
  ~response_vector() {}
  // zero_constructor
  response_vector(World& world, size_t num_orbs) : num_orbitals(num_orbs) {
    x = zero_functions<T>(world, num_orbitals);
    y = zero_functions<T>(world, num_orbitals);
  }
  // assignment operator
  response_vector& operator=(const response_vector& a) {
    num_orbitals = a.num_orbitals;
    x = a.x;
    y = y.x;
    return *this;
  }
  // same_size predicate
  bool same_size(const response_vector& a) {
    return (num_orbitals == a.num_orbitals);
  }
  // same_size predicate
  friend bool same_size(const response_vector& a, const response_vector& b) {
    return (a.num_orbitals == b.num_orbitals);
  }
  // operators for KAIN...expressive basis
  // addition +=
  response_vector operator+=(const response_vector& b) {
    MADNESS_ASSERT(same_size(b));
    World world = x[0].world();
    x = add(world, x, b.x);
    y = add(world, y, b.y);
    return *this;
  }
  // addition +  one sided
  response_vector operator+(const response_vector& b) {
    MADNESS_ASSERT(same_size(b));
    World world = x[0].world();
    response_vector result(world, num_orbitals);
    x = add(world, x, b.x);
    y = add(world, y, b.y);
    return result;
  }
  // addition  symmetric addition a+b
  friend response_vector operator+(const response_vector& a,
                                   const response_vector& b) {
    MADNESS_ASSERT(a.same_size(b));
    World world = a.x[0].world();
    response_vector result(world, a.num_orbitals);
    result.x = add(world, a.x, b.x);
    result.y = add(world, a.y, b.y);
    return result;
  }
  // Subtraction
  // sub -=
  response_vector operator-=(const response_vector& b) {
    MADNESS_ASSERT(same_size(b));
    World world = x[0].world();
    x = sub(world, x, b.x);
    y = sub(world, y, b.y);
    return *this;
  }
  // sub -  one sided
  response_vector operator-(const response_vector& b) {
    MADNESS_ASSERT(same_size(b));
    World world = x[0].world();
    // return rV(x-bx,y-by)
    return response_vector(sub(world, x, b.x), sub(world, y, b.y));
  }
  // sub  symmetric sub a-b
  friend response_vector operator-(const response_vector& a,
                                   const response_vector& b) {
    MADNESS_ASSERT(a.same_size(b));
    World world = a.x[0].world();
    return response_vector(sub(world, a.x, b.x), sub(world, b.y, b.y));
  }
  // multiplication
  response_vector operator*(double num_a) {
    World world = x[0].world();
    return response_vector(madness::scale(world, x, num_a, false),
                           madness::scale(world, y, num_a, false));
  }
  // multiplication *= a number
  response_vector operator*=(double num_a) {
    World world = x[0].world();
    madness::scale(world, x, num_a, false);
    madness::scale(world, y, num_a, false);
    return *this;
  }
  // multiplication *= by funcion f
  response_vector operator*=(function_T f) {
    World world = x[0].world();
    x = mul(world, f, x, false);
    y = mul(world, f, y, false);
    return *this;
  }
  // multiplication * by function f
  response_vector operator*(function_T f) {
    World world = x[0].world();
    return response_vector(mul(world, f, x, false), mul(world, f, y, false));
  }
  // madness procedures on functions
  // zero the functions
  void zero() {
    x = zero_functions<value_type, 3>(x[0].world(), num_orbitals);
    y = zero_functions<value_type, 3>(y[0].world(), num_orbitals);
  }
  void compress() {
    compress(x[0].world(), x, true);
    compress(y[0].world(), y, true);
  }
  void reconstruct() {
    reconstruct(x[0].world(), x, true);
    reconstruct(y[0].world(), y, true);
  }
  void truncate() {
    truncate(x[0].world(), x, true);
    truncate(y[0].world(), y, true);
  }
  // norm of vector
  double norm2() { return sqrt(inner(x, x) + inner(y, y)); }
  // inner of our Hilbert space of functions
  friend double inner(const response_vector& a, const response_vector& b) {
    MADNESS_ASSERT(a.num_orbitals > 0);
    MADNESS_ASSERT(a.same_size(b));
    return inner(a.x, b.x) + inner(a.y, b.y);
  }
};
// response allocator
template <typename T> struct response_vector_allocator {
  typedef Function<T, 3> function_T;
  World& world;
  const size_t num_orbitals;
  response_vector_allocator(World& world, const size_t num_orbitals)
      : world(world), num_orbitals(num_orbitals) {}
  // overloading ()
  response_vector<T> operator()() {
    return response_vector<T>(world, num_orbitals);
  }
  // Copy constructor
  response_vector_allocator operator=(const response_vector_allocator& other) {
    return response_vector_allocator(world, other.num_orbitals);
  }
};

} // namespace madness
#endif // SRC_APPS_ADRIAN_RESPONSE_VECTOR_H_
